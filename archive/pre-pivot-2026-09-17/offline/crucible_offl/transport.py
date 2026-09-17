"""OFFL-5 §3.1a — the constitutive spine's **chemical-regime transport
surface**, written to the FND-5 schema.

This is the product that retires the three declared transport constants in
`[mechanisms.wall]` (the degenerate FND-7 spine occupant) and gives the
resolved `F_visc` operator per-cell properties. Two engines, one state:

- **CEA** (the pinned equilibrium engine, `chemistry.py`) solves each node,
  exactly as the equilibrium surface does — same solver, same `gas_only`
  mode, same (p, h, Z) coordinate — and supplies the **caloric** columns.
  Using a second thermochemistry here would put a seam between the EOS the
  runtime projects onto and the heat capacity it linearizes with.
- **Cantera** evaluates **mixture-averaged Chapman-Enskog transport** on
  *that* composition at *that* (T, p): Lennard-Jones pure-species μ_k, k_k
  with the Wilke / Mathur-Tondon-Saxena mixture rules. It is never asked to
  equilibrate anything. [META-3: `chapman-enskog-mixavg`,
  `h2o2-transport-data`]

**Coordinate (FND-7 §3.7):** the runtime local state **(p, h, Z)** — the
same axes as the equilibrium surface, so a cell interrogates both surfaces
with one coordinate and one interpolation. The grid **spans exactly the
equilibrium surface's declared envelope**, and re-declares that envelope,
so the two surfaces refuse and accept on precisely the same set of states
(a gap between two independently declared tables is invisible to the
COUP-8 §3.3(2) coverage check, which sees each table alone).

**Why six columns and not three.** The runtime needs the *equilibrium*
conductivity — a dissociated gas moves far more energy down a temperature
gradient by recombination than by collisions — but folding that into a
tabulated `k_eq` would be a second model of a flux the operator already
carries. So the table ships the **molecular** pieces and the caloric
pieces, and the runtime composes them with the one declared Schmidt
closure (`ρD = μ/Sc`, META-3 `schmidt-combustion-gas`):

    Σ_k h_k j_k = −ρD [ (c_p,eq − c_p,fr)·∇T + (∂h/∂Z)|_{p,T}·∇Z ]

(exact on the equilibrium manifold Y_k(T, p, Z), by the chain rule). The
first term folds into the Fourier flux as `k_eff = k_fr + ρD·(c_p,eq −
c_p,fr)`, which is precisely the classical equilibrium conductivity; the
second is the resolved species-enthalpy flux. **One flux, decomposed —
never counted twice.**

The ∇p limb is neglected, and not for being small — measured, it is 22–45%
of the retained ∇T limb. It is neglected because in a boundary layer
∂p/∂n ≈ 0, so it vanishes where the diffusive flux matters, while
streamwise Pe ≫ 1 makes the whole diffusive flux negligible. See
`crates/solvers/src/transport.rs` for the full statement (one owner).

Bounds are measured exactly as `surfaces.py` measures them (holdout
midpoints + ¼-offsets + envelope edges, `SAFETY` margin, absolute and
rule-space) — the machinery is imported, not restated.

**Declared band: 10–20%** (FND-7 §3.3/§5) on μ and k — mixture-averaged vs
multicomponent, plus the pure-species LJ fits' own spread at combustion
temperatures. It is a `UncertainInput` declaration for COUP-5, never
folded into a value (FND-5 §2).

**Recorded model-form limits (owners named).** (1) Cantera's species
thermo fits stop at 3500 K while the hot corner reaches ~3900 K; the
`_CP_CROSSCHECK` gate below measures the consequence against CEA's Glenn
fits at *every* node and refuses if it exceeds the declared tolerance —
the extrapolation is bounded by test, not by hope. (2) Species outside
Cantera's set (O₃, at ≤ 1e-8 mole fraction here) are refused above
`_X_UNMAPPED_MAX` rather than silently dropped. (3) Soret/Dufour and Stefan-Maxwell
baro-diffusion are not represented (plan S5's species-vector state owns
them); the equilibrium composition's own pressure response is a separate
term, neglected for the boundary-layer reason above.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

import cantera as ct
import cea

from .chemistry import EquilibriumEngine
from .surfaces import (
    SAFETY,
    SCHEMA_VERSION,
    _edge_points,
    _holdout_points,
    _multilinear,
    _provenance,
)
from .tables import Axis, TableValue, WriteSpec, write_table

#: The Cantera deck: species set + the Lennard-Jones transport data. The
#: same file the OFFL-3 §6-2 equilibrium cross-check already pins.
CANTERA_MECH = "h2o2.yaml"
#: Chapman-Enskog mixture rules (Wilke / Mathur-Tondon-Saxena). The
#: `multicomponent` model is the band's other end, not a free choice — it
#: would be a new `data_version` with its own measured bounds.
TRANSPORT_MODEL = "mixture-averaged"

#: Central-difference step in Z for `∂h/∂Z|_{p,T}` (≈ 0.06% of the axis
#: values here — far above CEA's convergence noise, far below the axis's
#: own curvature scale). Fixed, declared: a derivative computed with an
#: adaptive step would not be bit-reproducible (META-1 §2.1).
DZ_FD = 1.0e-4

#: Refusal thresholds (META-1 P6 — fail loud at generation, never ship a
#: silently-degraded node).
#: Largest mole fraction a CEA product species absent from the Cantera set
#: may carry before the mapping is refused as lossy.
_X_UNMAPPED_MAX = 1.0e-6
#: Largest relative disagreement tolerated between Cantera's and CEA's
#: frozen c_p on the *same* composition and state, at or below `_CP_HOT_T`
#: — where both fit families (NASA-7 in `h2o2.yaml`, NASA-9 Glenn) are
#: inside their calibrated ranges: 2% is ~5× the observed agreement in the
#: chamber/throat class and still an order below the declared transport
#: band. It also fires on a composition mis-mapping (a wrong species moves
#: c_p at the ~10%+ scale), which is why it stays tight here; the
#: above-`_CP_HOT_T` extrapolation is bounded by `_CP_CROSSCHECK_HOT`.
_CP_CROSSCHECK = 2.0e-2
#: Above `_CP_HOT_T` the two fit families diverge by construction —
#: `h2o2.yaml`'s NASA-7 coefficients extrapolate past their fit ceiling
#: while Glenn NASA-9 remains calibrated — so the allowance there is the
#: MEASURED divergence with margin, not the mapping-bug gate: the 0.2.0
#: ignition-headroom envelope (OFFL-3 0.6.2; equilibrium h to +1.2e7,
#: T ≤ ~4300 K — dissociation buffers the hot corner) sweeps a smooth,
#: monotone-in-T drift peaking at 2.44% (p = 7e6 Pa, h = +1.2e7), so 4%
#: is measured-max ×1.6 — still 2.5–5× below the declared 10–20% band,
#: and far below the mapping-bug scale the check must keep catching.
_CP_CROSSCHECK_HOT = 4.0e-2
_CP_HOT_T = 3500.0


@dataclass(frozen=True)
class TransportGrid:
    """Axis grids + declared envelopes for the (p, h, Z) transport surface
    (SI). Built by `matched_transport_grid` from the equilibrium surface it
    must cover — never hand-written, so the two cannot drift apart."""

    p_points: tuple[float, ...]  # Pa, log-spaced, strictly increasing
    h_points: tuple[float, ...]  # J/kg
    z_points: tuple[float, ...]
    p_envelope: tuple[float, float]
    h_envelope: tuple[float, float]
    z_envelope: tuple[float, float]


def matched_transport_grid(
    eq_grid, n_p: int = 57, n_h: int = 37, n_z: int = 7
) -> TransportGrid:
    """The transport grid derived from an `EquilibriumGrid`: its axes span
    **exactly** that surface's declared envelope, and it re-declares the
    same envelope.

    Spanning the envelope rather than the equilibrium surface's wider grid
    is deliberate. The equilibrium grid's outer shell exists to bracket the
    envelope for interpolation; those nodes are states no run may legally
    interrogate, and two of the coldest low-pressure ones sit below the
    Cantera species-thermo floor. Generating them would mean either a
    fabricated value or a hole. Making the envelope the domain means every
    node is a state the runtime can actually reach — and a query at the
    envelope boundary is still in-domain (FND-5 §3.5 refuses *outside* the
    grid, and the boundary is not outside).

    Node counts default to roughly the equilibrium surface's density over
    the narrower span; they are a measured choice — the holdout bounds in
    `build_transport_surface` are the acceptance, and a poor bound means
    refine here."""
    for n, name in ((n_p, "n_p"), (n_h, "n_h"), (n_z, "n_z")):
        if n < 2:
            raise ValueError(f"{name}={n}: an axis needs at least 2 points")
    return TransportGrid(
        p_points=tuple(np.geomspace(*eq_grid.p_envelope, n_p)),
        h_points=tuple(np.linspace(*eq_grid.h_envelope, n_h)),
        z_points=tuple(np.linspace(*eq_grid.z_envelope, n_z)),
        p_envelope=tuple(eq_grid.p_envelope),
        h_envelope=tuple(eq_grid.h_envelope),
        z_envelope=tuple(eq_grid.z_envelope),
    )


#: (interp_rule, units) per transport column — ONE dict, so a new column
#: cannot ship with a fallback unit (the `surfaces.py` review finding).
#: Axes are (p, h, Z) with p logarithmic, then the value's own rule.
#: μ, k, and the heat capacities are strictly positive and span orders of
#: magnitude (c_p,eq runs 2.7e3 → 6e4 across the plume fringe), so they are
#: log-valued and carry the scale-honest rule-space bound. `dh_dz` is
#: linear: it is a signed quantity by nature (nothing forbids ∂h/∂Z < 0 for
#: another propellant) and a log rule would be a load-time refusal the day
#: it happened.
_TR_COLUMNS = {
    "viscosity": ("log-lin-lin-log", "Pa*s"),
    "conductivity_frozen": ("log-lin-lin-log", "W/(m*K)"),
    "cp_frozen": ("log-lin-lin-log", "J/(kg*K)"),
    "cp_equilibrium": ("log-lin-lin-log", "J/(kg*K)"),
    "cv_equilibrium": ("log-lin-lin-log", "J/(kg*K)"),
    "dh_dz": ("log-lin-lin-lin", "J/kg"),
}
#: Rule-only view (tests + the holdout gate key interpolation on it).
_TR_RULES = {name: rule for name, (rule, _units) in _TR_COLUMNS.items()}


class TransportEvaluator:
    """One CEA engine + one Cantera `Solution`, reused across many nodes.

    Deliberately *not* a subclass or a method of `EquilibriumEngine`: the
    thermochemistry engine has no business knowing about collision
    integrals, and this object has no business solving equilibria."""

    def __init__(self, engine: EquilibriumEngine):
        self.engine = engine
        self.gas = ct.Solution(CANTERA_MECH, transport_model=TRANSPORT_MODEL)
        if self.gas.transport_model != TRANSPORT_MODEL:
            raise RuntimeError(
                f"Cantera selected transport model {self.gas.transport_model!r}, "
                f"not the declared {TRANSPORT_MODEL!r}"
            )
        self._species = set(self.gas.species_names)

    def _cantera_transport(self, state) -> tuple[float, float]:
        """μ, k_frozen at the CEA state's own (T, p, composition)."""
        x, dropped = {}, 0.0
        for name, v in state.mole_fractions.items():
            if v <= 0.0:
                continue
            if name in self._species:
                x[name] = v
            else:
                dropped = max(dropped, v)
        if dropped > _X_UNMAPPED_MAX:
            raise RuntimeError(
                f"CEA product species absent from {CANTERA_MECH} carries mole "
                f"fraction {dropped:.3e} > {_X_UNMAPPED_MAX:.0e} at "
                f"p={state.p} Pa, T={state.T} K — the transport mixture would "
                "silently omit it; extend the mechanism, do not drop it"
            )
        if not x:
            raise RuntimeError(f"no mappable species at p={state.p}, T={state.T}")
        self.gas.TPX = state.T, state.p, x
        # The two-fit cross-check (module doc): Cantera's NASA-7 frozen c_p
        # against CEA's NASA-9 on the identical composition and state. It
        # is what bounds the >3500 K extrapolation, and it fires on the
        # composition mapping too — a mis-mapped species moves c_p first.
        rel = abs(self.gas.cp_mass - state.cp_fr) / state.cp_fr
        tol = _CP_CROSSCHECK if state.T <= _CP_HOT_T else _CP_CROSSCHECK_HOT
        if rel > tol:
            raise RuntimeError(
                f"Cantera/CEA frozen c_p disagree by {rel:.2%} (> "
                f"{tol:.0%} at T {'<=' if state.T <= _CP_HOT_T else '>'} "
                f"{_CP_HOT_T:.0f} K) at p={state.p} Pa, T={state.T} K: "
                f"{self.gas.cp_mass:.1f} vs {state.cp_fr:.1f} J/(kg·K) — the "
                "transport evaluation is not on the state CEA solved"
            )
        return float(self.gas.viscosity), float(self.gas.thermal_conductivity)

    def columns(self, p: float, h: float, z: float) -> dict[str, float]:
        """The six transport-surface columns at one (p, h, Z) node."""
        state = self.engine.state_php(p, h, z)
        mu, k_fr = self._cantera_transport(state)
        # ∂h/∂Z at FIXED (p, T) — two TP solves either side of the node.
        # Not (∂h/∂Z)|_{p,h}, which is identically zero, and not a
        # difference along the h axis, which would hold h fixed too: the
        # species-enthalpy flux coefficient is a constant-temperature
        # derivative because j_Z is driven by composition gradients the
        # temperature field does not already carry.
        hi = self.engine.state_tp(p, state.T, z + DZ_FD)
        lo = self.engine.state_tp(p, state.T, z - DZ_FD)
        return {
            "viscosity": mu,
            "conductivity_frozen": k_fr,
            "cp_frozen": state.cp_fr,
            "cp_equilibrium": state.cp_eq,
            "cv_equilibrium": state.cv_eq,
            "dh_dz": (hi.h - lo.h) / (2.0 * DZ_FD),
        }


def build_transport_surface(
    engine: EquilibriumEngine,
    grid: TransportGrid,
    data_version: str,
    generator_commit: str,
    holdout_stride: int = 3,
) -> tuple[WriteSpec, dict[str, float]]:
    """Evaluate every grid node, measure holdout errors, return the
    WriteSpec + the measured per-column absolute bounds (also stamped)."""
    ev = TransportEvaluator(engine)
    p_ax, h_ax, z_ax = (
        np.asarray(a) for a in (grid.p_points, grid.h_points, grid.z_points)
    )
    shape = (len(p_ax), len(h_ax), len(z_ax))
    grids = {name: np.empty(shape) for name in _TR_RULES}
    for i, p in enumerate(p_ax):
        for j, h in enumerate(h_ax):
            for k, z in enumerate(z_ax):
                for name, v in ev.columns(p, h, z).items():
                    grids[name][i, j, k] = v
    for name, rule in _TR_RULES.items():
        if rule.endswith("log") and not np.all(grids[name] > 0.0):
            bad = np.argwhere(grids[name] <= 0.0)[0]
            raise RuntimeError(
                f"column {name!r} is log-valued but node {tuple(bad)} holds "
                f"{grids[name][tuple(bad)]!r} — FND-5 would refuse the table at "
                "load; the rule or the physics is wrong"
            )

    axes_list = [p_ax, h_ax, z_ax]
    hold_p = np.concatenate(
        [
            _holdout_points(p_ax, log=True)[::holdout_stride],
            _edge_points(p_ax, True, grid.p_envelope),
        ]
    )
    hold_h = np.concatenate(
        [
            _holdout_points(h_ax, log=False)[::holdout_stride],
            _edge_points(h_ax, False, grid.h_envelope),
        ]
    )
    hold_z = np.concatenate(
        [
            _holdout_points(z_ax, log=False)[::holdout_stride],
            _edge_points(z_ax, False, grid.z_envelope),
        ]
    )
    env = (grid.p_envelope, grid.h_envelope, grid.z_envelope)
    bounds = {name: 0.0 for name in _TR_RULES}
    bounds_log = {name: 0.0 for name in _TR_RULES if _TR_RULES[name].endswith("log")}
    n_holdout = 0
    for p in hold_p:
        for h in hold_h:
            for z in hold_z:
                if not all(lo <= q <= hi for (lo, hi), q in zip(env, (p, h, z))):
                    continue
                n_holdout += 1
                truth = ev.columns(p, h, z)
                for name, t in truth.items():
                    est = _multilinear(_TR_RULES[name], axes_list, grids, (p, h, z), name)
                    bounds[name] = max(bounds[name], abs(est - t))
                    if name in bounds_log:
                        bounds_log[name] = max(
                            bounds_log[name], abs(np.log(est) - np.log(t))
                        )
    if n_holdout <= 0:
        raise RuntimeError("holdout set must not be empty")  # survives -O
    bounds = {name: SAFETY * b for name, b in bounds.items()}
    bounds_log = {name: SAFETY * b for name, b in bounds_log.items()}

    deck = {
        "table": "spine_transport_surface",
        "propellant": asdict(engine.propellant),
        "grid": asdict(grid),
        "columns": sorted(_TR_RULES),
        "holdout_stride": holdout_stride,
        "engine": f"cea {cea.__version__}",
        "products": "gas-only-metastable" if engine.gas_only else "full-condensed",
        # The model form is part of the deck hash: two different physics
        # choices must never share one (the v0.2.0/v0.3.0 review finding).
        "transport_engine": f"cantera {ct.__version__}",
        "transport_mech": CANTERA_MECH,
        "transport_model": TRANSPORT_MODEL,
        "dz_fd": DZ_FD,
    }
    prov = _provenance(generator_commit, deck)
    spec = WriteSpec(
        kind="regular",
        schema_version=SCHEMA_VERSION,
        data_version=data_version,
        interp_method="multilinear",
        provenance=type(prov)(
            producer="crucible-offl/cea+cantera",
            producer_version=(
                f"cea {cea.__version__} (libcea {cea.lib_version()}) + "
                f"cantera {ct.__version__}"
            ),
            input_deck_hash=prov.input_deck_hash,
            source_library=(
                "NASA Glenn thermo database (META-3 nasa-cea) + Cantera "
                f"{CANTERA_MECH} Lennard-Jones transport data (META-3 "
                "h2o2-transport-data), mixture-averaged Chapman-Enskog "
                "(META-3 chapman-enskog-mixavg)"
            ),
            generator_commit=generator_commit,
            rng_seed=None,
        ),
        axes=(
            Axis("p", tuple(p_ax), *grid.p_envelope),
            Axis("h", tuple(h_ax), *grid.h_envelope),
            Axis("Z", tuple(z_ax), *grid.z_envelope),
        ),
        values=tuple(
            TableValue(
                name,
                tuple(grids[name].reshape(-1)),
                units=_TR_COLUMNS[name][1],
                interp_rule=_TR_COLUMNS[name][0],
                interp_error_bound=bounds[name],
                interp_error_bound_log=bounds_log.get(name),
            )
            for name in sorted(_TR_COLUMNS)
        ),
    )
    return spec, bounds


def write_transport_table(
    path: str,
    data_version: str,
    generator_commit: str,
    grid: TransportGrid,
    engine: EquilibriumEngine | None = None,
    holdout_stride: int = 3,
) -> dict[str, str]:
    """Write the spine transport surface; return {group_path: digest}."""
    eng = engine if engine is not None else EquilibriumEngine(gas_only=True)
    spec, _bounds = build_transport_surface(
        eng, grid, data_version, generator_commit, holdout_stride
    )
    group = f"/spine/{eng.propellant.name}/transport"
    return {group: write_table(path, group, spec)}
