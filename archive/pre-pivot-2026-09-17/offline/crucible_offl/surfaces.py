"""OFFL-3 §3.3 — the tabulated products: the (p, h, Z) equilibrium surface
and the (p_c, MR) performance reference, written to the FND-5 schema.

Both tables carry **measured** interpolation-error bounds (FND-5 §3.4): the
generator interpolates every table column at holdout points (cell midpoints
*and* ¼-offset points — for smooth columns the midpoint dominates, but a
minor species near its exponential onset peaks off-center, which a
midpoint-only sweep understates; the §6-4 gate caught exactly that) in the
same rule space the runtime uses, solves CEA directly there, and stamps
``SAFETY ×`` the observed max abs error — the margin covers the finite
sampling of a continuous sup and is declared here, not hidden. The bound is
table metadata for COUP-5's epistemic interval — never folded into a value.

Equilibrium-surface honesty note: the surface is *the equilibrium*,
including condensed H₂O where (p, h, Z) demands it (deep-cold fuel-rich
corners of the rectangular grid). The `condensed_fraction` column reports
the condensed mole fraction so a consumer can see where the surface leaves
gas-only territory; the certificate declares the gas-only region. Refusing
to tabulate those corners would leave invalid nodes inside the rectangle —
worse than reporting the truth (META-1 P6: never guess).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass

import numpy as np

import cea

from .chemistry import (
    EquilibriumEngine,
    EqState,
    FrozenReactantEngine,
    FrozenReactantState,
    Propellant,
)
from .tables import Axis, Provenance, TableValue, WriteSpec, write_table

SCHEMA_VERSION = "1.0"
GAS_SPECIES = ("H2O", "H2", "OH", "H", "O2", "O")
#: Declared sampling margin on the measured holdout max (see module doc).
#: 1.25 was observed insufficient on badly under-resolved grids (a fresh
#: ⅜-offset sweep found 1.28× the measured max on a 15×11×5 subgrid);
#: 1.5 + the §6-4 fresh-holdout CI gate on the shipped artifact is the
#: enforcement loop — a violation there means "refine the grid", never
#: "relax the gate".
SAFETY = 1.5


def _holdout_points(ax: np.ndarray, log: bool) -> np.ndarray:
    """Cell midpoints + ¼-offsets in the axis's rule space."""
    a = np.log(ax) if log else ax
    pts = np.concatenate([0.5 * (a[:-1] + a[1:]), 0.75 * a[:-1] + 0.25 * a[1:]])
    return np.exp(pts) if log else pts


def _edge_points(ax: np.ndarray, log: bool, env: tuple[float, float]) -> np.ndarray:
    """Envelope-edge coverage (session-12 review): a declared envelope
    bound can cut through an axis cell, leaving the in-envelope slab of
    that cell structurally unreachable by midpoint/quarter offsets (the
    envelope filter then discards them all) — the stamped bound silently
    excluded that slab. Sample the envelope boundary itself plus the
    midpoint between each bound and its nearest interior grid node, in
    rule space."""
    lo, hi = env
    a = np.log(ax) if log else ax
    el, eh = (np.log(lo), np.log(hi)) if log else (lo, hi)
    pts = [el, eh]
    inner = a[(a > el) & (a < eh)]
    if inner.size:
        pts.append(0.5 * (el + inner[0]))
        pts.append(0.5 * (eh + inner[-1]))
    out = np.array(pts)
    return np.exp(out) if log else out


@dataclass(frozen=True)
class EquilibriumGrid:
    """Axis grids + declared envelopes for the (p, h, Z) surface (SI)."""

    p_points: tuple[float, ...]  # Pa, log-spaced, strictly increasing
    h_points: tuple[float, ...]  # J/kg
    z_points: tuple[float, ...]
    p_envelope: tuple[float, float]
    h_envelope: tuple[float, float]
    z_envelope: tuple[float, float]


@dataclass(frozen=True)
class PerformanceGrid:
    """Axis grids + envelopes for the (p_c, MR) performance reference."""

    pc_points: tuple[float, ...]  # Pa
    mr_points: tuple[float, ...]
    pc_envelope: tuple[float, float]
    mr_envelope: tuple[float, float]


def design_window_grid() -> EquilibriumGrid:
    """The W2 design-window grid (OFFL-3 R2: dense design core; a
    resolved-tier version later widens Z as a new data_version)."""
    return EquilibriumGrid(
        p_points=tuple(np.geomspace(1.5e3, 8.0e6, 41)),
        h_points=tuple(np.linspace(-1.18e7, -1.0e5, 31)),
        z_points=tuple(np.linspace(0.10, 0.26, 13)),
        p_envelope=(2.0e3, 7.0e6),
        h_envelope=(-1.15e7, -2.0e5),
        z_envelope=(1.0 / 9.0, 0.25),  # MR 8 … 3
    )


def station5_envelope_grid() -> EquilibriumGrid:
    """The station-5 full-engine grid (data_version 0.3.x–0.4.x): the design-window
    density with the pressure axis widened DOWN to the vacuum-plume fringe
    (the epsilon = 61 exit runs ~3 kPa static; lip transients dip further)
    and — 0.3.0, session 12 — the COLD fringe covered honestly:

    - The settled coarse RL10 field held 35 gas cells pinned at the 0.2.x
      enthalpy floor (h = -1.15e7 ⇒ T ≈ 698 K at fringe pressures) and the
      dial-8 establishment transient stepped through it: the floor was
      binding physics, not margin.
    - The surface is generated GAS-ONLY (metastable equilibrium, stamped in
      the deck): the declared model for a rapidly-expanding plume — real
      plumes supersaturate (condensation kinetics slow vs flow time) — and
      the only CEA-convergent branch in the deep-cold corners. The
      `condensed_fraction` column is identically 0 here; fringe validity is
      this declared model form, not that column.
    - The Z axis is NARROWED to the premixed operating class (the prior-tier
      field holds Z = Z_inj everywhere — element advection is source-free —
      so station-5 never leaves the injector's Z): that is what buys the
      cold floor on a rectangular grid, whose fuel-rich edge binds CEA's
      low-T convergence. The MR 3–8 width remains the design-window /
      resolved-tier table's property.

    Floor placement: gas-only CEA converges at h = -1.25e7 across the full
    p × Z grid (T down to ~210 K at the rich edge, ~300 K at Z = 1/6);
    -1.26e7 fails at Z ≥ 0.195. Envelope floor -1.23e7 (T ≈ 410 K at
    Z = 1/6) with the declared ~1 mbar altitude-cell ambient floor keeps
    the settled fringe (T ≈ 460 K class) strictly inside — zero pinned
    cells, checkable in the fields artifact. Envelope widening is a
    table-version setting of the same pipeline, never a new pipeline
    (OFFL-3 §3.3 R2 doctrine)."""
    return EquilibriumGrid(
        p_points=tuple(np.geomspace(5.0, 8.0e6, 69)),  # ~11 pts/decade, as v0.1
        # 0.3.x ceiling: establishment transients overshoot the injection
        # enthalpy — the dial-16 piston start reached h ≈ -1e5 (0.3.1
        # raised to +1.3e6) and the dial-12 overexpanded-bell backflow
        # recompression then reached +1.33e6 (measured at each halt's
        # crash artifact). 0.3.2 set the ceiling with transient margin:
        # +4.0e6 grid / +3.8e6 envelope (T ~ 4100 K class).
        # 0.4.0 (plan S6, OFFL-3 0.6.2 — the ignition-headroom ruling):
        # the SOLV-4 §3.6 blend interrogates BOTH branches at the same
        # cell enthalpy, and a spark-heated igniting kernel (T_u ~
        # 1000-1300 K => h ~ 2.5-4e6 J/kg) rode the TOP of this surface
        # while still valid cold gas — the burnt ceiling sat BELOW the
        # unburnt surface's (+5.0e6), an envelope inversion. Standing
        # rule: burnt ceiling >> unburnt ceiling, so a burning cell can
        # never refuse where the same cold gas was fine. STRICT
        # EXTENSION discipline: the 0.3.2 43-node axis is reproduced
        # bit-exact and 21 nodes appended above at the SAME spacing
        # (ceiling +1.225e7; T ~ 5000-6000 K class, trivially
        # CEA-convergent, still chemistry) — every in-old-envelope
        # interpolation, i.e. all five station certificates, is
        # byte-identical; only pins/digests/bounds metadata move.
        h_points=tuple(
            np.concatenate(
                [
                    np.linspace(-1.25e7, 4.0e6, 43),  # the 0.3.2 axis, bit-exact
                    4.0e6 + (1.65e7 / 42.0) * np.arange(1, 22),  # same ~3.9e5 spacing
                ]
            )
        ),
        z_points=tuple(np.linspace(0.145, 0.195, 11)),
        p_envelope=(1.0e1, 7.0e6),
        h_envelope=(-1.23e7, 1.2e7),
        z_envelope=(0.155, 0.185),  # MR 5.45 … 4.41 (design 5.0 = Z 1/6 mid)
    )


def design_window_performance_grid() -> PerformanceGrid:
    return PerformanceGrid(
        pc_points=tuple(np.linspace(1.0e6, 6.0e6, 11)),
        mr_points=tuple(np.linspace(3.0, 8.0, 11)),
        pc_envelope=(1.2e6, 5.5e6),
        mr_envelope=(3.5, 7.5),
    )


def _condensed_fraction(state: EqState) -> float:
    return sum(v for k, v in state.mole_fractions.items() if "(" in k)


def _eq_columns(state: EqState) -> dict[str, float]:
    cols = {
        "temperature": state.T,
        "density": state.rho,
        "gamma_eff": state.gamma,
        "sound_speed": state.a,
        "mbar": state.mbar,
        "condensed_fraction": _condensed_fraction(state),
    }
    for sp in GAS_SPECIES:
        cols[f"X_{sp}"] = state.mole_fractions.get(sp, 0.0)
    return cols


#: (interp_rule, units) per equilibrium column — ONE table, so a new
#: column cannot silently ship with a fallback unit (review finding: a
#: parallel units dict with .get(name, "1") was the silent-default hazard
#: the Rust expect_units bind gate would then bake in). Axes (p, h, Z)
#: then the value; p is a log axis everywhere; density is log-valued
#: (near-linear in log p – log rho); everything else linear (mole
#: fractions can be 0).
_EQ_COLUMNS = {
    "temperature": ("log-lin-lin-lin", "K"),
    "density": ("log-lin-lin-log", "kg/m^3"),
    "gamma_eff": ("log-lin-lin-lin", "1"),
    "sound_speed": ("log-lin-lin-lin", "m/s"),
    "mbar": ("log-lin-lin-lin", "kg/kmol"),
    "condensed_fraction": ("log-lin-lin-lin", "1"),
    **{f"X_{sp}": ("log-lin-lin-lin", "1") for sp in GAS_SPECIES},
}
#: Rule-only view (the certificate + tests key interpolation on it).
_EQ_RULES = {name: rule for name, (rule, _units) in _EQ_COLUMNS.items()}


def _multilinear(rule: str, axes: list[np.ndarray], grids: dict[str, np.ndarray], q: tuple[float, ...], name: str) -> float:
    """Reference multilinear evaluation in the rule's linearizing space —
    mirrors the runtime (FND-5 §3.3) for holdout error measurement."""
    tokens = rule.split("-")
    ax_scaled = []
    q_scaled = []
    for t, ax, qi in zip(tokens[:-1], axes, q):
        if t == "log":
            ax_scaled.append(np.log(ax))
            q_scaled.append(np.log(qi))
        else:
            ax_scaled.append(ax)
            q_scaled.append(qi)
    data = grids[name]
    if tokens[-1] == "log":
        data = np.log(data)
    # Locate + fold one axis at a time.
    weights = []
    idx = []
    for ax, qi in zip(ax_scaled, q_scaled):
        i = int(np.clip(np.searchsorted(ax, qi) - 1, 0, len(ax) - 2))
        weights.append((ax[i + 1] - qi) / (ax[i + 1] - ax[i]))
        idx.append(i)
    v = 0.0
    for corner in range(2 ** len(idx)):
        w = 1.0
        pos = []
        for d in range(len(idx)):
            hi = (corner >> d) & 1
            w *= (1.0 - weights[d]) if hi else weights[d]
            pos.append(idx[d] + hi)
        v += w * data[tuple(pos)]
    return float(np.exp(v)) if tokens[-1] == "log" else float(v)


def _provenance(generator_commit: str, deck: dict) -> Provenance:
    deck_hash = hashlib.sha256(
        json.dumps(deck, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return Provenance(
        producer="crucible-offl/cea",
        producer_version=f"cea {cea.__version__} (libcea {cea.lib_version()})",
        input_deck_hash=f"sha256:{deck_hash}",
        source_library="NASA Glenn thermo database (bundled with cea, META-3 nasa-cea)",
        generator_commit=generator_commit,
        rng_seed=None,
    )


def build_equilibrium_surface(
    engine: EquilibriumEngine,
    grid: EquilibriumGrid,
    data_version: str,
    generator_commit: str,
    holdout_stride: int = 2,
) -> tuple[WriteSpec, dict[str, float]]:
    """Solve every grid node, measure holdout midpoint errors, return the
    WriteSpec + the measured per-column bounds (also stamped in the spec)."""
    p_ax, h_ax, z_ax = (np.asarray(a) for a in (grid.p_points, grid.h_points, grid.z_points))
    shape = (len(p_ax), len(h_ax), len(z_ax))
    grids = {name: np.empty(shape) for name in _EQ_RULES}
    for i, p in enumerate(p_ax):
        for j, h in enumerate(h_ax):
            for k, z in enumerate(z_ax):
                for name, v in _eq_columns(engine.state_php(p, h, z)).items():
                    grids[name][i, j, k] = v

    # Holdout: midpoints + ¼-offsets (each axis's rule space), subsampled
    # by stride, restricted to the declared envelope; SAFETY margin on the
    # observed max (module doc).
    axes_list = [p_ax, h_ax, z_ax]
    hold_p = np.concatenate(
        [_holdout_points(p_ax, log=True)[::holdout_stride], _edge_points(p_ax, True, grid.p_envelope)]
    )
    hold_h = np.concatenate(
        [_holdout_points(h_ax, log=False)[::holdout_stride], _edge_points(h_ax, False, grid.h_envelope)]
    )
    hold_z = np.concatenate(
        [_holdout_points(z_ax, log=False)[::holdout_stride], _edge_points(z_ax, False, grid.z_envelope)]
    )
    env = (grid.p_envelope, grid.h_envelope, grid.z_envelope)
    bounds = {name: 0.0 for name in _EQ_RULES}
    bounds_log = {name: 0.0 for name in _EQ_RULES if _EQ_RULES[name].endswith("log")}
    n_holdout = 0
    for p in hold_p:
        for h in hold_h:
            for z in hold_z:
                if not all(lo <= q <= hi for (lo, hi), q in zip(env, (p, h, z))):
                    continue
                n_holdout += 1
                truth = _eq_columns(engine.state_php(p, h, z))
                for name, t in truth.items():
                    est = _multilinear(_EQ_RULES[name], axes_list, grids, (p, h, z), name)
                    bounds[name] = max(bounds[name], abs(est - t))
                    if name in bounds_log:
                        # Rule-space (|Δ ln|, ≈ relative) error for
                        # log-valued columns — the scale-honest bound the
                        # runtime acceptance uses (session-12 review).
                        bounds_log[name] = max(bounds_log[name], abs(np.log(est) - np.log(t)))
    if n_holdout <= 0:
        raise RuntimeError("holdout set must not be empty")  # survives -O
    bounds = {name: SAFETY * b for name, b in bounds.items()}
    bounds_log = {name: SAFETY * b for name, b in bounds_log.items()}

    deck = {
        "table": "equilibrium_surface",
        "propellant": asdict(engine.propellant),
        "grid": asdict(grid),
        "columns": sorted(_EQ_RULES),
        "holdout_stride": holdout_stride,
        "engine": f"cea {cea.__version__}",
        "products": "gas-only-metastable" if engine.gas_only else "full-condensed",
    }
    spec = WriteSpec(
        kind="regular",
        schema_version=SCHEMA_VERSION,
        data_version=data_version,
        interp_method="multilinear",
        provenance=_provenance(generator_commit, deck),
        axes=(
            Axis("p", tuple(p_ax), *grid.p_envelope),
            Axis("h", tuple(h_ax), *grid.h_envelope),
            Axis("Z", tuple(z_ax), *grid.z_envelope),
        ),
        values=tuple(
            TableValue(
                name,
                tuple(grids[name].reshape(-1)),
                units=_EQ_COLUMNS[name][1],
                interp_rule=_EQ_COLUMNS[name][0],
                interp_error_bound=bounds[name],
                interp_error_bound_log=bounds_log.get(name),
            )
            for name in sorted(_EQ_COLUMNS)
        ),
    )
    return spec, bounds


def build_performance_reference(
    engine: EquilibriumEngine,
    grid: PerformanceGrid,
    data_version: str,
    generator_commit: str,
) -> tuple[WriteSpec, dict[str, float]]:
    """The (p_c, MR) chamber performance functional (OFFL-3 §2): c*_ideal,
    T_c, γ, M̄ — SOLV-7's anchor and SOLV-1 §3.4's knockdown reference.
    The provenance deck stamps the product model form (session-12 review:
    v0.2.0/v0.3.0 shipped identical input_deck_hash across different
    generation physics)."""
    pc_ax, mr_ax = np.asarray(grid.pc_points), np.asarray(grid.mr_points)
    names = ("c_star_ideal", "T_c", "gamma", "mbar")
    rules = {n: "lin-lin-lin" for n in names}
    grids = {n: np.empty((len(pc_ax), len(mr_ax))) for n in names}
    for i, pc in enumerate(pc_ax):
        for j, mr in enumerate(mr_ax):
            perf = engine.performance(pc, mr)
            grids["c_star_ideal"][i, j] = perf.c_star
            grids["T_c"][i, j] = perf.T_c
            grids["gamma"][i, j] = perf.gamma_c
            grids["mbar"][i, j] = perf.mbar_c

    axes_list = [pc_ax, mr_ax]
    env = (grid.pc_envelope, grid.mr_envelope)
    bounds = {n: 0.0 for n in names}
    n_holdout = 0
    for pc in _holdout_points(pc_ax, log=False):
        for mr in _holdout_points(mr_ax, log=False):
            if not all(lo <= q <= hi for (lo, hi), q in zip(env, (pc, mr))):
                continue
            n_holdout += 1
            perf = engine.performance(pc, mr)
            truth = {
                "c_star_ideal": perf.c_star,
                "T_c": perf.T_c,
                "gamma": perf.gamma_c,
                "mbar": perf.mbar_c,
            }
            for n, t in truth.items():
                est = _multilinear(rules[n], axes_list, grids, (pc, mr), n)
                bounds[n] = max(bounds[n], abs(est - t))
    if n_holdout <= 0:
        raise RuntimeError("holdout set must not be empty")  # survives -O
    bounds = {n: SAFETY * b for n, b in bounds.items()}

    deck = {
        "table": "performance_reference",
        "propellant": asdict(engine.propellant),
        "grid": asdict(grid),
        "columns": list(names),
        "engine": f"cea {cea.__version__}",
        "products": "gas-only-metastable" if engine.gas_only else "full-condensed",
    }
    spec = WriteSpec(
        kind="regular",
        schema_version=SCHEMA_VERSION,
        data_version=data_version,
        interp_method="multilinear",
        provenance=_provenance(generator_commit, deck),
        axes=(
            Axis("p_c", tuple(pc_ax), *grid.pc_envelope),
            Axis("MR", tuple(mr_ax), *grid.mr_envelope),
        ),
        values=tuple(
            TableValue(
                n,
                tuple(grids[n].reshape(-1)),
                units={"c_star_ideal": "m/s", "T_c": "K", "gamma": "1", "mbar": "kg/kmol"}[n],
                interp_rule=rules[n],
                interp_error_bound=bounds[n],
            )
            for n in sorted(names)
        ),
    )
    return spec, bounds


#: (interp_rule, units) per unburnt-reactant column — the SAME schema and
#: rule strings as the equilibrium surface's shared columns, so `TableEos`
#: binds either with no new occupant (SOLV-1 §3.4). Axes (p, h, Z): p log,
#: h/Z linear; density is log-valued (ideal gas ρ ∝ p spans orders),
#: everything else linear. No X_species columns — the frozen reactant
#: composition is fixed by Z, so tabulating it would be redundant bytes.
_UNBURNT_COLUMNS = {
    "temperature": ("log-lin-lin-lin", "K"),
    "density": ("log-lin-lin-log", "kg/m^3"),
    "gamma_eff": ("log-lin-lin-lin", "1"),
    "sound_speed": ("log-lin-lin-lin", "m/s"),
    "mbar": ("log-lin-lin-lin", "kg/kmol"),
}
_UNBURNT_RULES = {name: rule for name, (rule, _u) in _UNBURNT_COLUMNS.items()}


def unburnt_reactant_grid(
    engine: FrozenReactantEngine,
    n_p: int = 9,
    n_h: int = 121,
    n_z: int = 11,
    t_floor: float = 100.0,
    t_ceil: float = 2200.0,
    cold_face: bool = True,
) -> EquilibriumGrid:
    """The (p, h, Z) grid + envelope for the unburnt-reactant surface
    (OFFL-3 §3.3, plan S5), DERIVED from the engine so the rectangular
    h-envelope is one every Z can convergently populate.

    The frozen reactant enthalpy is Z-dependent (H₂-rich mixtures store more
    enthalpy per K), so a fixed h maps to a different T at each Z. A
    rectangular envelope valid across the whole Z band therefore runs its
    **hot floor at the H₂-poor edge and its cold floor at the H₂-rich edge**:
    `h_env_lo = max_Z h(t_floor, Z)` guarantees every Z is at least `t_floor`
    at the envelope's cold edge, and `h_env_hi = min_Z h(t_ceil, Z)` bounds
    every Z below `t_ceil` at the hot edge. The grid overhangs the envelope
    by a few percent of the h-span (the equilibrium grid's convention), so
    the envelope boundary interpolates from bracketing nodes.

    `t_floor = 100 K` is the declared cold floor of the **gas-phase** branch:
    below the liquefaction line the ideal-gas frozen mixture is a declared
    metastable model, and the real two-phase state is the SOLV-1 W4 drift-
    flux extension (plan S15). `t_ceil = 2200 K` was the 0.1.0 ceiling
    ("pre-ignition compression with headroom"); the S6 igniter mini-sims
    measured that headroom as insufficient — blast-focused compression in a
    confined ignition drives mid-transition (0 < c < 1) cells past it while
    they still read the unburnt branch (OFFL-3 0.6.2). `t_ceil_ext` (0.2.0)
    therefore appends nodes ABOVE the 0.1.0 axis — a **strict extension**
    (the base `n_h` nodes are reproduced bit-exact; same spacing above), so
    every in-old-envelope value is unchanged — up to the metastable-reactant
    ceiling the ignition surface's T_u envelope covers (the two surfaces'
    envelope-consistency contract: T_u(h_env_hi) must be inside the ignition
    surface's T_u envelope, or a hot transition cell refuses on the closure
    query instead). Reactants at ~2900 K "should have reacted" — which is
    exactly what τ_ign(T_u) says (sub-µs there): the metastable branch is
    self-consistently transient, never an equilibrium claim."""
    # 0.3.0 (plan S7): the Z band NARROWS to the burnt surface's own grid
    # band (0.145–0.195, envelope 0.155–0.185 — the premixed design-window
    # class the blend intersects to anyway): the rectangular h-envelope's
    # floor binds at the O2-rich Z edge, and the wide S5 band put the
    # mid-Z effective floor ~60 K above the declared t_floor. The narrow
    # band makes the rectangle nearly tight in Z, so the cold face reaches
    # the wall-cooled startup states (120 K coolant class) the S7 march
    # actually holds. A DECLARED re-gridded variant (OFFL-3 §3.3 R2: the
    # envelope is a table-version setting), not a strict extension — no
    # certificate consumes this surface's numbers; the wide-Z 0.2.0 grid
    # remains on disk as the S5b species-work base.
    zs = np.linspace(0.145, 0.195, n_z)
    t_ceil_ext = 2900.0
    # 0.3.0 (plan S7) — the COLD face. The rectangular h-envelope's floor
    # binds at the H2-POOR Z edge, which put the mid-Z design line's
    # effective floor at ~145 K while the S7 startup march holds REAL
    # gas-phase states there: fill gas cooled by the 120 K coolant wall at
    # 1–20 kPa sits 45+ K above its own O2 saturation (T_sat ≈ 61–82 K at
    # the band's O2 partial pressures). `t_floor_ext = 75 K` re-derives the
    # box floor so the mid-Z line is valid to ~108 K; at the binding
    # H2-poor edge 75 K is honestly gas-phase for p ≲ 2 kPa and a DECLARED
    # metastable (supersaturated-vapor) overhang toward the ~50 kPa+
    # corner — the same declared-metastable discipline as the sub-floor
    # grid overhang; real condensation is the S15 drift-flux wave. The
    # cold nodes are PREPENDED at the same spacing (the 0.1.0 base axis is
    # anchored to the old floor arithmetic and stays bit-exact — strict
    # extension).
    t_floor_ext = 75.0
    h_env_lo = max(engine.enthalpy_at(t_floor, float(z)) for z in zs)  # the wide-Z anchor arithmetic
    h_env_hi = min(engine.enthalpy_at(t_ceil, float(z)) for z in zs)
    if not h_env_lo < h_env_hi:
        raise RuntimeError(
            f"unburnt grid: empty h-envelope [{h_env_lo}, {h_env_hi}] — "
            "t_floor/t_ceil cross once the Z-dependence is folded in"
        )
    span = h_env_hi - h_env_lo
    margin = 0.03 * span
    base = np.linspace(h_env_lo - margin, h_env_hi + margin, n_h)  # the 0.1.0 axis, bit-exact
    delta = (base[-1] - base[0]) / (n_h - 1)
    h_env_ext = min(engine.enthalpy_at(t_ceil_ext, float(z)) for z in zs)
    n_ext = int(np.ceil((h_env_ext + margin - base[-1]) / delta))
    # Bracket the cold envelope edge: prepend cells only if the base axis's
    # own 3% margin does not already cover it (for the SHIPPED narrow-Z
    # grid the margin suffices — n_cold = 0, deepest node ~40 K equivalent,
    # grid-min γ 1.151; the prepend machinery serves other
    # parameterizations). `cold_face = False` (probe grids) keeps the base
    # envelope: a probe-coarse delta would put an overhang cell below the
    # engine's solvable bracket.
    if cold_face:
        h_env_lo_ext = max(engine.enthalpy_at(t_floor_ext, float(z)) for z in zs)
        n_cold = max(0, int(np.ceil((base[0] - h_env_lo_ext) / delta)))
    else:
        h_env_lo_ext = h_env_lo
        n_cold = 0
    h_points = np.concatenate(
        [
            base[0] - delta * np.arange(n_cold, 0, -1),
            base,
            base[-1] + delta * np.arange(1, n_ext + 1),
        ]
    )
    return EquilibriumGrid(
        p_points=tuple(np.geomspace(5.0, 8.0e6, n_p)),
        h_points=tuple(h_points),
        z_points=tuple(zs),
        p_envelope=(10.0, 7.0e6),
        h_envelope=(h_env_lo_ext, h_env_ext),
        z_envelope=(0.155, 0.185),  # the burnt surface's own band (MR 4.41–5.45)
    )


def _unburnt_columns(state: FrozenReactantState) -> dict[str, float]:
    return {
        "temperature": state.T,
        "density": state.rho,
        "gamma_eff": state.gamma,
        "sound_speed": state.a,
        "mbar": state.mbar,
    }


def build_unburnt_surface(
    engine: FrozenReactantEngine,
    grid: EquilibriumGrid,
    data_version: str,
    generator_commit: str,
    holdout_stride: int = 2,
) -> tuple[WriteSpec, dict[str, float]]:
    """Solve every grid node of the gas-phase frozen reactant mixture,
    measure holdout errors exactly as `build_equilibrium_surface` does
    (midpoints + ¼-offsets + envelope edges, `SAFETY` margin, absolute +
    rule-space), return the WriteSpec + per-column bounds."""
    p_ax, h_ax, z_ax = (
        np.asarray(a) for a in (grid.p_points, grid.h_points, grid.z_points)
    )
    shape = (len(p_ax), len(h_ax), len(z_ax))
    grids = {name: np.empty(shape) for name in _UNBURNT_RULES}
    # T, γ, a, M̄ are pressure-independent (ideal gas) — cache them per
    # (h, Z) so the surface is built with one T-inversion per column, not
    # one per (p, h, Z) node.
    for j, h in enumerate(h_ax):
        for k, z in enumerate(z_ax):
            base = engine.state_php(float(p_ax[0]), float(h), float(z))
            for i, p in enumerate(p_ax):
                # Ideal gas: ρ ∝ p at fixed (T, M̄), so scale the base state
                # rather than re-derive — exact in real arithmetic (~1 ULP in
                # f64, far below the stamped bound) and it keeps the sole ρ
                # formula in `state_php` (one owner).
                grids["density"][i, j, k] = base.rho * (float(p) / float(p_ax[0]))
                grids["temperature"][i, j, k] = base.T
                grids["gamma_eff"][i, j, k] = base.gamma
                grids["sound_speed"][i, j, k] = base.a
                grids["mbar"][i, j, k] = base.mbar

    axes_list = [p_ax, h_ax, z_ax]
    hold_p = np.concatenate(
        [_holdout_points(p_ax, log=True)[::holdout_stride], _edge_points(p_ax, True, grid.p_envelope)]
    )
    hold_h = np.concatenate(
        [_holdout_points(h_ax, log=False)[::holdout_stride], _edge_points(h_ax, False, grid.h_envelope)]
    )
    hold_z = np.concatenate(
        [_holdout_points(z_ax, log=False)[::holdout_stride], _edge_points(z_ax, False, grid.z_envelope)]
    )
    env = (grid.p_envelope, grid.h_envelope, grid.z_envelope)
    bounds = {name: 0.0 for name in _UNBURNT_RULES}
    bounds_log = {name: 0.0 for name in _UNBURNT_RULES if _UNBURNT_RULES[name].endswith("log")}
    n_holdout = 0
    for p in hold_p:
        for h in hold_h:
            for z in hold_z:
                if not all(lo <= q <= hi for (lo, hi), q in zip(env, (p, h, z))):
                    continue
                n_holdout += 1
                truth = _unburnt_columns(engine.state_php(float(p), float(h), float(z)))
                for name, t in truth.items():
                    est = _multilinear(_UNBURNT_RULES[name], axes_list, grids, (p, h, z), name)
                    bounds[name] = max(bounds[name], abs(est - t))
                    if name in bounds_log:
                        bounds_log[name] = max(bounds_log[name], abs(np.log(est) - np.log(t)))
    if n_holdout <= 0:
        raise RuntimeError("holdout set must not be empty")  # survives -O
    bounds = {name: SAFETY * b for name, b in bounds.items()}
    bounds_log = {name: SAFETY * b for name, b in bounds_log.items()}

    deck = {
        "table": "unburnt_reactant_surface",
        "propellant": asdict(engine.propellant),
        "grid": asdict(grid),
        "columns": sorted(_UNBURNT_RULES),
        "holdout_stride": holdout_stride,
        "engine": f"cea {cea.__version__}",
        "products": "gas-phase-frozen-reactant-mixture",  # the declared model form
    }
    spec = WriteSpec(
        kind="regular",
        schema_version=SCHEMA_VERSION,
        data_version=data_version,
        interp_method="multilinear",
        provenance=_provenance(generator_commit, deck),
        axes=(
            Axis("p", tuple(p_ax), *grid.p_envelope),
            Axis("h", tuple(h_ax), *grid.h_envelope),
            Axis("Z", tuple(z_ax), *grid.z_envelope),
        ),
        values=tuple(
            TableValue(
                name,
                tuple(grids[name].reshape(-1)),
                units=_UNBURNT_COLUMNS[name][1],
                interp_rule=_UNBURNT_COLUMNS[name][0],
                interp_error_bound=bounds[name],
                interp_error_bound_log=bounds_log.get(name),
            )
            for name in sorted(_UNBURNT_COLUMNS)
        ),
    )
    return spec, bounds


def write_unburnt_table(
    path: str,
    data_version: str,
    generator_commit: str,
    engine: FrozenReactantEngine | None = None,
    grid: EquilibriumGrid | None = None,
    holdout_stride: int = 2,
) -> dict[str, str]:
    """Write the unburnt-reactant surface (the c = 0 branch); return
    {group_path: digest} for pinning."""
    eng = engine if engine is not None else FrozenReactantEngine()
    spec, _ = build_unburnt_surface(
        eng, grid or unburnt_reactant_grid(eng), data_version, generator_commit, holdout_stride
    )
    group = f"/chem/{eng.propellant.name}/unburnt"
    return {group: write_table(path, group, spec)}


def write_station3_tables(
    path: str,
    data_version: str,
    generator_commit: str,
    propellant: Propellant | None = None,
    eq_grid: EquilibriumGrid | None = None,
    perf_grid: PerformanceGrid | None = None,
    gas_only: bool = False,
) -> dict[str, str]:
    """Generate both station-3 tables into one HDF5 file; return
    {group_path: content_digest} for pinning. `gas_only` selects the
    metastable product set (see EquilibriumEngine; deck-stamped)."""
    engine = (
        EquilibriumEngine(propellant, gas_only=gas_only)
        if propellant
        else EquilibriumEngine(gas_only=gas_only)
    )
    eq_spec, _ = build_equilibrium_surface(
        engine, eq_grid or design_window_grid(), data_version, generator_commit
    )
    perf_spec, _ = build_performance_reference(
        engine, perf_grid or design_window_performance_grid(), data_version, generator_commit
    )
    name = engine.propellant.name
    digests = {}
    for group, spec in (
        (f"/chem/{name}/equilibrium", eq_spec),
        (f"/chem/{name}/performance", perf_spec),
    ):
        digests[group] = write_table(path, group, spec)
    return digests
