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

from .chemistry import EquilibriumEngine, EqState, Propellant
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
    """The station-5 full-engine grid (data_version 0.3.x): the design-window
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
        h_points=tuple(np.linspace(-1.25e7, -1.0e5, 33)),  # ~3.9e5 J/kg spacing kept
        z_points=tuple(np.linspace(0.145, 0.195, 11)),
        p_envelope=(1.0e1, 7.0e6),
        h_envelope=(-1.23e7, -2.0e5),
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
    hold_p = _holdout_points(p_ax, log=True)[::holdout_stride]
    hold_h = _holdout_points(h_ax, log=False)[::holdout_stride]
    hold_z = _holdout_points(z_ax, log=False)[::holdout_stride]
    env = (grid.p_envelope, grid.h_envelope, grid.z_envelope)
    bounds = {name: 0.0 for name in _EQ_RULES}
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
    assert n_holdout > 0, "holdout set must not be empty"
    bounds = {name: SAFETY * b for name, b in bounds.items()}

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
    T_c, γ, M̄ — SOLV-7's anchor and SOLV-1 §3.4's knockdown reference."""
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
    assert n_holdout > 0, "holdout set must not be empty"
    bounds = {n: SAFETY * b for n, b in bounds.items()}

    deck = {
        "table": "performance_reference",
        "propellant": asdict(engine.propellant),
        "grid": asdict(grid),
        "columns": list(names),
        "engine": f"cea {cea.__version__}",
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
