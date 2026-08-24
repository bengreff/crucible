"""OFFL-3 §3.3 — the ignition/flame closure surface (plan S6).

The SOLV-4 §3.6 burn-progress rate law reads two closures over the local
**unburnt** state (p, T_u, Z): the laminar flame speed ``S_L`` (the
propagation term ``ρ_u·S_T·|∇c|``) and the induction time ``τ_ign`` (the
auto-ignition term ``ρ(1−c)/τ_ign``). Both are finite-rate chemistry done at
its correct scale **offline** on the cited ``h2-kinetics-mech`` (Cantera's
bundled ``h2o2.yaml`` — GRI-Mech 3.0 H₂/O₂ subset) and homogenised into one
(p, T_u, Z) surface with two columns; the runtime only interpolates (no
runtime kinetics network).

**Extinction is the columns' own values, so a runtime cell never branches
(Rule 12):**

* ``laminar_flame_speed`` is **linear-valued** and **→ 0** where the 1-D
  flame does not light (a non-converging ``FreeFlame`` ⇒ ``S_L = 0``, the
  flammability-limit value) — a value 0 is representable, a log column could
  not carry it.
* ``ignition_delay`` is **log-valued** and **capped at a declared
  ``TAU_MAX``** — a reactor that does not ignite within the horizon returns
  ``TAU_MAX`` (never ``∞``), so the auto-ignition source ``ρ(1−c)/τ_ign``
  vanishes smoothly at the quench/low-T corner.

Determinism (META-1 §2.1, Tier-2): every solve uses declared, fixed solver
settings (flame width + refine criteria; reactor horizon + trigger), a fresh
``Solution`` per call (no carried state), and no RNG — verified
bit-reproducible cross-process by the regen probe. The Cantera pin is part
of the record: the provenance deck stamps ``cantera {version}`` + the
mechanism (the OFFL-5 transport-feed pattern).
"""

from __future__ import annotations

import sys

from dataclasses import asdict, dataclass

import numpy as np

import cantera as ct

import cea

from .chemistry import LOX_LH2, Propellant, gasify
from .surfaces import (
    SAFETY,
    SCHEMA_VERSION,
    _edge_points,
    _holdout_points,
    _multilinear,
    _provenance,
)
from .tables import Axis, TableValue, WriteSpec, write_table

#: The cited H₂/O₂ mechanism (META-3 `h2-kinetics-mech`): GRI-Mech 3.0 H₂/O₂
#: submechanism as bundled with Cantera. A pin bump is a `data_version` event.
CANTERA_MECH = "h2o2.yaml"
#: Mixture-averaged transport is required for the FreeFlame diffusion terms;
#: asserted at construction so a silent model change cannot slip in.
TRANSPORT_MODEL = "mixture-averaged"

#: Declared cap on the induction time [s]: a reactor that has not ignited by
#: `REACTOR_HORIZON` returns this (never `∞`). Chosen a decade above the
#: horizon so "capped" and "slow but real" never collide in the log column;
#: the runtime source `ρ(1−c)/τ_ign` at `τ_max = 10 s` is ~10⁻¹ s⁻¹ per unit
#: (1−c) — negligible against any resolved flow timescale, i.e. quenched.
TAU_MAX = 10.0
#: The 0-D const-P reactor integration horizon [s]. Past this, no ignition.
REACTOR_HORIZON = 1.0
#: Ignition = the instant of steepest temperature rise (the canonical
#: shock-tube inflection definition); a run whose total rise never exceeds
#: this threshold [K] never ignited (returns TAU_MAX). 50 K is far above the
#: reactor's numerical drift and far below any real exotherm.
IGN_MIN_RISE_K = 50.0

#: FreeFlame domain half-width [m] and the declared fixed refine criteria.
#: These are the determinism contract, not tunables — a change is a
#: `data_version` event. Verified bit-reproducible cross-process (regen probe).
FLAME_WIDTH_M = 0.05
FLAME_RATIO = 3.0
FLAME_SLOPE = 0.08
FLAME_CURVE = 0.08
#: Burning velocities below this [m/s] are read as non-propagating (0): a
#: FreeFlame can converge to a numerically tiny residual velocity at a
#: marginal mixture; that is the flammability limit, i.e. extinction.
SL_FLOOR_M_PER_S = 1.0e-3
#: Burning velocities above this [m/s] are the **auto-ignition regime**, not a
#: real laminar flame: past the T_u where the unburnt mixture auto-ignites
#: faster than a flame can transit its own preheat zone, `FreeFlame` no longer
#: describes a propagating front and its `velocity[0]` runs away (10³–10⁵ m/s —
#: the whole domain igniting). Real H₂/O₂ laminar S_L stays well under 300 m/s
#: (peak ~150 at 1100 K/30 bar), so a solve above the ceiling is read as
#: S_L = 0 — the propagation term switches off and the auto-ignition term
#: (`τ_ign`, tiny at those T_u) governs, exactly where it should (Rule 12: the
#: surface's own value, not a runtime branch).
SL_CEILING_M_PER_S = 300.0


#: The two closure columns and their FND-5 (interp_rule, units). S_L is
#: linear-valued (carries 0 at the limit); τ_ign is log-valued (spans
#: orders, capped finite). Axes are (p log, T_u lin, Z lin).
_IGNITION_COLUMNS: dict[str, tuple[str, str]] = {
    "laminar_flame_speed": ("log-lin-lin-lin", "m/s"),
    "ignition_delay": ("log-lin-lin-log", "s"),
}
_IGNITION_RULES = {name: rule for name, (rule, _units) in _IGNITION_COLUMNS.items()}


def _mass_fractions(propellant: Propellant, z: float) -> str:
    """The unburnt reactant mixture at mixture fraction Z as a Cantera
    mass-fraction string. Z = fuel mass fraction (Z = 1/(1+MR)); the streams
    are gasified (H₂(L)→H₂) exactly as the frozen-reactant surface (OFFL-3
    §3.3), so this surface and the unburnt-thermo surface describe the SAME
    gas."""
    fuel = [gasify(s) for s in propellant.fuel_species]
    ox = [gasify(s) for s in propellant.oxidizer_species]
    if len(fuel) != 1 or len(ox) != 1:
        raise ValueError("ignition surface supports one fuel + one oxidizer species")
    return f"{fuel[0]}:{z}, {ox[0]}:{1.0 - z}"


class IgnitionEngine:
    """One Cantera ``Solution`` on the cited mechanism; evaluates ``S_L`` (1-D
    freely-propagating flame) and ``τ_ign`` (0-D const-P reactor) on the
    unburnt reactant mixture. All boundaries SI (Pa, K, m/s, s)."""

    def __init__(self, propellant: Propellant = LOX_LH2) -> None:
        self.propellant = propellant
        # Construct once to validate the mechanism + species mapping; each
        # solve builds its own fresh Solution for state isolation.
        gas = ct.Solution(CANTERA_MECH, transport_model=TRANSPORT_MODEL)
        if gas.transport_model != TRANSPORT_MODEL:
            raise RuntimeError(
                f"mechanism {CANTERA_MECH} did not take transport model "
                f"{TRANSPORT_MODEL!r} (got {gas.transport_model!r})"
            )
        self._species = set(gas.species_names)
        for s in (*propellant.fuel_species, *propellant.oxidizer_species):
            if gasify(s) not in self._species:
                raise ValueError(
                    f"reactant species {gasify(s)!r} not in mechanism {CANTERA_MECH}"
                )

    def laminar_flame_speed(self, p: float, t_u: float, z: float) -> float:
        """S_L [m/s] of the freely-propagating premixed flame; a
        non-converging solve is the flammability limit ⇒ 0."""
        gas = ct.Solution(CANTERA_MECH, transport_model=TRANSPORT_MODEL)
        gas.TPY = t_u, p, _mass_fractions(self.propellant, z)
        flame = ct.FreeFlame(gas, width=FLAME_WIDTH_M)
        flame.set_refine_criteria(
            ratio=FLAME_RATIO, slope=FLAME_SLOPE, curve=FLAME_CURVE
        )
        flame.transport_model = TRANSPORT_MODEL
        try:
            flame.solve(loglevel=0, auto=True)
        except Exception:
            return 0.0
        v = float(flame.velocity[0])
        if not np.isfinite(v) or v < SL_FLOOR_M_PER_S or v > SL_CEILING_M_PER_S:
            return 0.0
        return v

    def ignition_delay(self, p: float, t_u: float, z: float) -> float:
        """τ_ign [s]: the time of steepest temperature rise in a 0-D const-P
        reactor; a run that never rises by IGN_MIN_RISE_K returns TAU_MAX."""
        gas = ct.Solution(CANTERA_MECH, transport_model=TRANSPORT_MODEL)
        gas.TPY = t_u, p, _mass_fractions(self.propellant, z)
        reactor = ct.IdealGasConstPressureReactor(gas)
        net = ct.ReactorNet([reactor])
        t_prev, temp_prev = 0.0, reactor.T
        best_rate, best_t = -1.0, TAU_MAX
        max_temp = reactor.T
        while net.time < REACTOR_HORIZON:
            net.step()
            temp = reactor.T
            dt = net.time - t_prev
            if dt > 0.0:
                rate = (temp - temp_prev) / dt
                if rate > best_rate:
                    best_rate, best_t = rate, net.time
            max_temp = max(max_temp, temp)
            t_prev, temp_prev = net.time, temp
        if max_temp - t_u < IGN_MIN_RISE_K:
            return TAU_MAX
        return min(best_t, TAU_MAX)

    def columns(self, p: float, t_u: float, z: float) -> dict[str, float]:
        return {
            "laminar_flame_speed": self.laminar_flame_speed(p, t_u, z),
            "ignition_delay": self.ignition_delay(p, t_u, z),
        }


@dataclass(frozen=True)
class IgnitionGrid:
    """Axis grids + declared envelopes for the (p, T_u, Z) closure surface (SI)."""

    p_points: tuple[float, ...]  # Pa, log-spaced
    tu_points: tuple[float, ...]  # K, unburnt temperature
    z_points: tuple[float, ...]  # fuel mass fraction Z = 1/(1+MR)
    p_envelope: tuple[float, float]
    tu_envelope: tuple[float, float]
    z_envelope: tuple[float, float]


def ignition_grid(
    n_p: int = 6,
    n_tu: int = 12,
    n_z: int = 6,
    p_lo: float = 3.0e3,
    p_hi: float = 1.0e7,
    tu_lo: float = 150.0,
    tu_hi: float = 1900.0,
    z_lo: float = 0.06,
    z_hi: float = 0.30,
    n_tu_ext: int = 7,
) -> IgnitionGrid:
    """A (p, T_u, Z) grid spanning the chemical-slice startup envelope: p from
    fill to chamber, T_u from cryo-adjacent to auto-ignitive, Z from lean to
    rich of the H₂/O₂ stoichiometric Z = 1/9. The grid overhangs the declared
    envelope by a small margin so envelope-edge queries interpolate.

    `n_tu_ext` (0.2.0, plan S6 close) appends nodes ABOVE the 0.1.0 T_u axis
    at the same spacing — a **strict extension** (the base `n_tu` nodes are
    bit-exact, so every in-old-envelope value is unchanged) — lifting the
    T_u ceiling ~1900 → ~3000 K to honor the two-surface envelope-consistency
    contract (OFFL-3 0.6.2): the rate law queries S_L/τ_ign at
    T_u = T_unburnt(p, h, Z), so this surface's T_u envelope must cover the
    unburnt surface's own h-ceiling (t_ceil_ext = 2900 K), or a hot
    mid-transition cell refuses on the closure query. Physically the new rows
    are the auto-ignitive regime: `FreeFlame` no longer describes a
    propagating front there (S_L reads 0 via the ceiling guard) and τ_ign is
    sub-µs — the surface's own values say "burns essentially instantly",
    which is the self-consistent statement about superheated reactants."""
    p = np.geomspace(p_lo, p_hi, n_p)
    tu_base = np.linspace(tu_lo, tu_hi, n_tu)  # the 0.1.0 axis, bit-exact
    dtu = (tu_base[1] - tu_base[0]) if n_tu > 1 else 0.0
    tu = np.concatenate([tu_base, tu_hi + dtu * np.arange(1, n_tu_ext + 1)])
    z = np.linspace(z_lo, z_hi, n_z)
    # Overhang the rectangular envelope inside the grid extremes so an
    # envelope-edge runtime query is bracketed (never an extrapolation).
    dp = (np.log(p[1]) - np.log(p[0])) if n_p > 1 else 0.0
    dz = (z[1] - z[0]) if n_z > 1 else 0.0
    return IgnitionGrid(
        p_points=tuple(p),
        tu_points=tuple(tu),
        z_points=tuple(z),
        p_envelope=(float(np.exp(np.log(p_lo) + 0.5 * dp)), float(np.exp(np.log(p_hi) - 0.5 * dp))),
        tu_envelope=(float(tu_lo + 0.5 * dtu), float(tu[-1] - 0.5 * dtu)),
        z_envelope=(float(z_lo + 0.5 * dz), float(z_hi - 0.5 * dz)),
    )


def build_ignition_surface(
    engine: IgnitionEngine,
    grid: IgnitionGrid,
    data_version: str,
    generator_commit: str,
    holdout_stride: int = 2,
) -> tuple[WriteSpec, dict[str, float]]:
    """Solve every grid node (flame + reactor), measure holdout errors in the
    columns' own rule space, return the WriteSpec + measured per-column
    bounds (also stamped)."""
    p_ax, tu_ax, z_ax = (np.asarray(a) for a in (grid.p_points, grid.tu_points, grid.z_points))
    shape = (len(p_ax), len(tu_ax), len(z_ax))
    grids = {name: np.empty(shape) for name in _IGNITION_RULES}
    for i, p in enumerate(p_ax):
        for j, t_u in enumerate(tu_ax):
            # Progress per (p, T_u) row (flame solves dominate; a silent
            # 30-minute generator is undiagnosable when a solve dies).
            print(
                f"  ignition row p={p:.3e} Pa ({i + 1}/{len(p_ax)}) "
                f"T_u={t_u:.0f} K ({j + 1}/{len(tu_ax)})",
                file=sys.stderr,
                flush=True,
            )
            for k, z in enumerate(z_ax):
                cols = engine.columns(float(p), float(t_u), float(z))
                for name in _IGNITION_RULES:
                    grids[name][i, j, k] = cols[name]

    # Holdout: fresh points (midpoint/quarter, subsampled) + envelope edges,
    # each disjoint from the tests' ⅜-offset gate. Measure in rule space.
    hold_p = np.concatenate(
        [_holdout_points(p_ax, True)[::holdout_stride], _edge_points(p_ax, True, grid.p_envelope)]
    )
    hold_tu = np.concatenate(
        [_holdout_points(tu_ax, False)[::holdout_stride], _edge_points(tu_ax, False, grid.tu_envelope)]
    )
    hold_z = np.concatenate(
        [_holdout_points(z_ax, False)[::holdout_stride], _edge_points(z_ax, False, grid.z_envelope)]
    )
    axes = [p_ax, tu_ax, z_ax]
    (plo, phi), (tlo, thi), (zlo, zhi) = grid.p_envelope, grid.tu_envelope, grid.z_envelope
    bounds = {name: 0.0 for name in _IGNITION_RULES}
    bounds_log = {name: 0.0 for name in _IGNITION_RULES}
    # Per-column ACTIVE-REGION holdout (the declared-sharp-feature pattern).
    # Each column carries a sharp EXTINCTION transition that no interpolant can
    # represent — S_L steps to 0 at the flame->auto-ignition crossover, tau_ign
    # steps to the TAU_MAX cap at the ignitable->quenched crossover — but each
    # transition sits exactly where the OTHER term dominates the SOLV-4 §3.6
    # rate law (S_L->0 where tau_ign is tiny; tau_ign->cap where a flame
    # propagates), so a column only governs the physics in its ACTIVE region.
    # The stamped bound therefore measures each column only in grid cells whose
    # EIGHT corners are all active (the smooth interior); a cell straddling the
    # extinction boundary is the declared model feature, not an interpolation
    # the bound must cover. Fully-active is tested by interpolating a 0/1 active
    # mask (lin) and checking it equals 1.0 at the holdout point. `n_active`
    # per column must be > 0 or the grid does not cover that region.
    def _active(name: str, value: float) -> bool:
        if name == "laminar_flame_speed":
            return value > SL_FLOOR_M_PER_S  # flammable (nonzero S_L)
        return value < TAU_MAX * (1.0 - 1.0e-9)  # ignited within the horizon

    masks = {
        name: np.vectorize(lambda v, n=name: 1.0 if _active(n, v) else 0.0)(grids[name])
        for name in _IGNITION_RULES
    }
    mask_rules = {name: rule.rsplit("-", 1)[0] + "-lin" for name, rule in _IGNITION_RULES.items()}

    n_holdout = 0
    n_active = {name: 0 for name in _IGNITION_RULES}
    for p in hold_p:
        if not (plo <= p <= phi):
            continue
        for t_u in hold_tu:
            if not (tlo <= t_u <= thi):
                continue
            for z in hold_z:
                if not (zlo <= z <= zhi):
                    continue
                n_holdout += 1
                q = (float(p), float(t_u), float(z))
                truth = engine.columns(*q)
                for name, rule in _IGNITION_RULES.items():
                    # Fully-active cell? (all 8 corners active ⇒ mask interp 1.0)
                    if _multilinear(mask_rules[name], axes, masks, q, name) < 1.0 - 1.0e-12:
                        continue  # extinction transition — declared, not bounded
                    est = _multilinear(rule, axes, grids, q, name)
                    t = truth[name]
                    n_active[name] += 1
                    bounds[name] = max(bounds[name], abs(est - t))
                    if rule.endswith("log") and est > 0.0 and t > 0.0:
                        bounds_log[name] = max(bounds_log[name], abs(np.log(est) - np.log(t)))
    if n_holdout <= 0:
        raise RuntimeError("ignition surface holdout set is empty — grid/envelope mismatch")
    for name, n in n_active.items():
        if n <= 0:
            raise RuntimeError(
                f"ignition column {name!r} has no in-active-region holdout points — "
                "the grid does not resolve its flammable/ignitable region"
            )
    bounds = {name: SAFETY * b for name, b in bounds.items()}
    bounds_log = {name: SAFETY * b for name, b in bounds_log.items()}

    deck = {
        "table": "ignition",
        "propellant": asdict(engine.propellant),
        "grid": asdict(grid),
        "columns": sorted(_IGNITION_COLUMNS),
        "holdout_stride": holdout_stride,
        "kinetics_engine": f"cantera {ct.__version__}",
        "kinetics_mech": CANTERA_MECH,
        "transport_model": TRANSPORT_MODEL,
        "tau_max_s": TAU_MAX,
        "reactor_horizon_s": REACTOR_HORIZON,
        "ign_min_rise_k": IGN_MIN_RISE_K,
        "flame": {
            "width_m": FLAME_WIDTH_M,
            "ratio": FLAME_RATIO,
            "slope": FLAME_SLOPE,
            "curve": FLAME_CURVE,
            "sl_floor_m_per_s": SL_FLOOR_M_PER_S,
        },
        "products": "finite-rate H2/O2 (offline); runtime interpolation only (SOLV-4 §3.6)",
    }
    base = _provenance(generator_commit, deck)
    spec = WriteSpec(
        kind="regular",
        schema_version=SCHEMA_VERSION,
        data_version=data_version,
        interp_method="multilinear",
        provenance=type(base)(
            producer="crucible-offl/cantera",
            producer_version=f"cantera {ct.__version__}",
            input_deck_hash=base.input_deck_hash,
            source_library=(
                f"Cantera {CANTERA_MECH} — GRI-Mech 3.0 H2/O2 submechanism "
                "(META-3 h2-kinetics-mech); 1-D freely-propagating flame (S_L) "
                "+ 0-D const-P reactor (tau_ign)"
            ),
            generator_commit=generator_commit,
            rng_seed=None,
        ),
        axes=(
            Axis("p", tuple(p_ax), *grid.p_envelope),
            Axis("T_u", tuple(tu_ax), *grid.tu_envelope),
            Axis("Z", tuple(z_ax), *grid.z_envelope),
        ),
        values=tuple(
            TableValue(
                name,
                tuple(grids[name].reshape(-1)),
                units=_IGNITION_COLUMNS[name][1],
                interp_rule=_IGNITION_COLUMNS[name][0],
                interp_error_bound=bounds[name],
                interp_error_bound_log=bounds_log.get(name) or None,
            )
            for name in sorted(_IGNITION_COLUMNS)
        ),
    )
    return spec, bounds


def write_ignition_table(
    path: str,
    data_version: str,
    generator_commit: str,
    engine: IgnitionEngine | None = None,
    grid: IgnitionGrid | None = None,
    holdout_stride: int = 2,
) -> dict[str, str]:
    """Write the (p, T_u, Z) ignition closure surface; return {group: digest}."""
    eng = engine if engine is not None else IgnitionEngine()
    grd = grid if grid is not None else ignition_grid()
    spec, _bounds = build_ignition_surface(eng, grd, data_version, generator_commit, holdout_stride)
    group = f"/chem/{eng.propellant.name}/ignition"
    return {group: write_table(path, group, spec)}
