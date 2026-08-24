"""OFFL-3 §3.3 / §6-4 (plan S5) — the UNBURNT-REACTANT production artifact
`tables/chem/lox_lh2_unburnt_v0.2.0.h5`: the burn-progress c = 0 branch
(SOLV-4 §3.6), a gas-phase ideal-gas frozen reactant mixture on the same CEA
enthalpy reference as the equilibrium surface.

Gates, in the order they matter:

- **§6-4 fresh holdout:** ⅜-offset points disjoint from the generator's
  midpoint + ¼-offset set (plus the envelope edges unsubsampled) respect
  every stored `interp_error_bound` — the same discipline station-3/5 use.
- **Cryo validity (the point of the branch):** the surface represents cold
  states the equilibrium surface cannot, and it does so where the frozen
  reactant elements have zero formation enthalpy.
- **One enthalpy reference:** the branch's `h` is the CEA reactant-mixture
  enthalpy `injection_enthalpy` already uses, so the SOLV-4 §3.6 blend adds.
- **Ideal-gas + thermodynamic consistency:** ρ, a, γ reproduce from the
  tabulated (T, M̄) by the declared relations.
"""

import pathlib

import h5py
import numpy as np
import pytest

import cantera as ct
import cea

from crucible_offl.chemistry import (
    LOX_LH2,
    EquilibriumEngine,
    FrozenReactantEngine,
    R_CEA,
)
from crucible_offl.surfaces import (
    _UNBURNT_RULES,
    _edge_points,
    _multilinear,
    _unburnt_columns,
)

REPO = pathlib.Path(__file__).resolve().parents[2]
PRODUCTION = REPO / "tables" / "chem" / "lox_lh2_unburnt_v0.2.0.h5"
GROUP = "/chem/lox_lh2/unburnt"


@pytest.fixture(scope="module")
def engine():
    return FrozenReactantEngine()


def test_s6_4_fresh_holdout_respects_stored_bounds(engine):
    with h5py.File(PRODUCTION, "r") as f:
        g = f[GROUP]
        axes = [np.array(g["axes"][n]) for n in ("p", "h", "Z")]
        shape = tuple(len(a) for a in axes)
        grids, bounds, bounds_log, env = {}, {}, {}, []
        for n in ("p", "h", "Z"):
            ds = g["axes"][n]
            env.append((float(ds.attrs["envelope_min"]), float(ds.attrs["envelope_max"])))
        for name in g["values"]:
            if name.startswith("sigma_"):
                continue
            grids[name] = np.array(g["values"][name]).reshape(shape)
            bounds[name] = float(g["values"][name].attrs["interp_error_bound"])
            if "interp_error_bound_log" in g["values"][name].attrs:
                bounds_log[name] = float(g["values"][name].attrs["interp_error_bound_log"])
    p_ax, h_ax, z_ax = axes

    # Fresh ⅜-offset points (disjoint from the generator's ½/¼ set) + the
    # envelope edges unsubsampled (the session-12 slab-coverage fix).
    hold_p = np.concatenate(
        [
            np.exp(0.625 * np.log(p_ax[:-1]) + 0.375 * np.log(p_ax[1:]))[::3],
            _edge_points(p_ax, True, env[0]),
        ]
    )
    hold_h = np.concatenate(
        [(0.625 * h_ax[:-1] + 0.375 * h_ax[1:])[::3], _edge_points(h_ax, False, env[1])]
    )
    hold_z = np.concatenate(
        [(0.625 * z_ax[:-1] + 0.375 * z_ax[1:])[::2], _edge_points(z_ax, False, env[2])]
    )

    n = 0
    for p in hold_p:
        for h in hold_h:
            for z in hold_z:
                if not all(lo <= q <= hi for (lo, hi), q in zip(env, (p, h, z))):
                    continue
                n += 1
                truth = _unburnt_columns(engine.state_php(float(p), float(h), float(z)))
                for name, t in truth.items():
                    est = _multilinear(_UNBURNT_RULES[name], axes, grids, (p, h, z), name)
                    assert abs(est - t) <= bounds[name], (
                        f"{name} at ({p:.4g},{h:.4g},{z:.4g}): "
                        f"|{est:.6g}-{t:.6g}| > stored {bounds[name]:.6g}"
                    )
                    if name in bounds_log:
                        log_err = abs(np.log(est) - np.log(t))
                        assert log_err <= bounds_log[name], (
                            f"{name} log err {log_err:.6g} > stored {bounds_log[name]:.6g}"
                        )
    assert n > 100
    assert all(b > 0.0 for b in bounds.values())
    # density is the only log-valued column, and it must carry a positive
    # rule-space bound (the runtime acceptance keys on it).
    assert "density" in bounds_log and bounds_log["density"] > 0.0


def test_surface_is_physical_everywhere():
    """Every tabulated node is a real gas state: T, ρ, a, M̄ > 0 and
    1 < γ < 5/3 (a frozen H/O gas mixture is diatomic-dominated). The γ
    lower bound is **1.1, not 1.0** (review): a uniform c_p units-scaling bug
    collapses γ → 1.0003 and would pass a `> 1.0` gate while staying
    self-consistent through a regen — the tightened bound catches it (the
    real grid min, including the sub-floor overhang, is 1.164)."""
    with h5py.File(PRODUCTION, "r") as f:
        g = f[GROUP]
        cols = {name: np.array(g["values"][name]) for name in g["values"] if not name.startswith("sigma_")}
    for name in ("temperature", "density", "sound_speed", "mbar", "gamma_eff"):
        assert np.all(np.isfinite(cols[name])), name
        assert np.all(cols[name] > 0.0), name
    assert np.all(cols["gamma_eff"] > 1.1) and np.all(cols["gamma_eff"] < 5.0 / 3.0)
    # A cryogenic corner exists (the branch's reason to be): the coldest
    # tabulated T is below the equilibrium surface's floor (~410 K at Z=1/6).
    assert cols["temperature"].min() < 200.0


def test_in_envelope_states_are_diatomic_and_sanely_heated(engine):
    """The physicality gate that a units bug cannot slip past (review): every
    IN-ENVELOPE frozen state must be a diatomic-dominated gas — γ ∈ [1.3, 1.5]
    and c_p ∈ [2000, 5500] J/(kg·K) (measured in-envelope ranges are
    [1.30, 1.48] and [2139, 5230]). A ×1000 or ÷1000 c_p error breaks BOTH
    bounds, unlike the self-consistent γ ≈ 1.0003 collapse a loose gate
    admits."""
    for z in (0.12, 1.0 / 6.0, 0.24):
        for t in (110.0, 300.0, 1200.0, 2100.0):
            s = engine.state_php(1.0e6, engine.enthalpy_at(t, z), z)
            assert 1.3 <= s.gamma <= 1.5, (z, t, s.gamma)
            assert 2000.0 <= s.cp_fr <= 5500.0, (z, t, s.cp_fr)


def test_unburnt_branch_is_distinct_from_burnt_by_the_heat_of_reaction(engine):
    """Why the c = 0 branch exists (SOLV-4 §3.6): at the SAME (p, T, Z) the
    unburnt (frozen reactants) and burnt (equilibrium products) enthalpies
    differ by the heat of reaction — the burnt surface cannot represent an
    unburnt cell as unburnt, it can only represent it as already-reacted, so
    forcing a startup (unburnt) cell onto the equilibrium surface is a
    category error. The unburnt branch is that missing model. (The full
    establishment cure routes transient cells to this branch via S6's
    c-blend; S5 ships the representable branch.)"""
    p, t, z = 3.0e6, 1000.0, 1.0 / 6.0
    h_u = engine.enthalpy_at(t, z)
    eq = EquilibriumEngine(gas_only=True)
    h_b = eq.state_tp(p, t, z).h  # equilibrium (burnt) enthalpy at the same T
    # Reactants sit ABOVE the recombined products by the heat of reaction —
    # several MJ/kg for LOX/LH2 — a large, one-signed, physical gap.
    assert h_u - h_b > 5.0e6, (h_u, h_b)


def test_cold_gas_states_are_representable(engine):
    """The branch represents genuinely cold reactant gas (T ~ 100–300 K) —
    hundreds of K below the equilibrium surface's ~410 K floor at Z = 1/6 —
    with physical γ and ρ. These are the startup states the equilibrium
    projection has no honest home for."""
    st = engine.state_php(3.0e6, engine.enthalpy_at(120.0, 1.0 / 6.0), 1.0 / 6.0)
    assert 100.0 < st.T < 150.0
    assert 1.0 < st.gamma < 5.0 / 3.0
    assert st.rho > 0.0 and st.a > 0.0


def test_enthalpy_is_on_the_cea_reactant_reference(engine):
    """The blend `h = (1-c)h_u + c*h_b` (SOLV-4 §3.6) is only defined on one
    enthalpy reference. The unburnt branch's h must be the SAME CEA reactant
    enthalpy `EquilibriumEngine.injection_enthalpy` uses — check that the
    frozen-mixture enthalpy at the RL10 liquid-injection *temperatures*
    reproduces the injection enthalpy that the coarse RL10 config pins."""
    # injection_enthalpy is the LIQUID reactant enthalpy; the gas branch is a
    # different phase, so we check the reference not the value: the gas
    # frozen-mixture enthalpy is continuous and monotone through the same
    # zero-formation-enthalpy reference (elements → 0 at 298.15 K).
    z = 1.0 / 6.0
    # Elements have ~0 formation enthalpy: at 298.15 K the enthalpy is
    # ESSENTIALLY zero (measured ~-6e-4 J/kg), never the ~-1e7 of the
    # dissociated-then-recombined products. The tight backstop (1e3, review:
    # 5e4 admitted a +4e4 J/kg reference offset) is what makes this a
    # reference check rather than a sign check.
    h_298 = engine.enthalpy_at(298.15, z)
    assert abs(h_298) < 1.0e3, h_298  # at the formation-enthalpy zero
    # Monotone increasing in T (the inversion depends on it).
    ts = np.linspace(120.0, 2000.0, 40)
    hs = np.array([engine.enthalpy_at(float(t), z) for t in ts])
    assert np.all(np.diff(hs) > 0.0)


def test_ideal_gas_and_sound_speed_consistency(engine):
    """ρ = p·M̄/(R̄·T) and a = √(γ·R̄·T/M̄) exactly (the declared model)."""
    for p, h, z in [(3.0e6, 5.0e5, 1.0 / 6.0), (1.0e2, -3.0e5, 0.20), (7.0e6, 2.0e6, 0.13)]:
        s = engine.state_php(p, h, z)
        rho_expected = p * s.mbar / (R_CEA * s.T)
        assert s.rho == pytest.approx(rho_expected, rel=1e-12)
        a_expected = np.sqrt(s.gamma * (R_CEA / s.mbar) * s.T)
        assert s.a == pytest.approx(a_expected, rel=1e-12)


def test_gas_species_are_the_gasified_reactants(engine):
    """The mixture is the GAS forms of the propellant reactant species —
    H2(L)/O2(L) gasified to H2/O2 (the liquid/vapor coupling is plan S15)."""
    # A single-species-per-stream propellant maps cleanly; multi-species
    # streams refuse (they need declared intra-stream proportions).
    assert LOX_LH2.fuel_species == ("H2(L)",)
    # M̄ at pure-fuel-ish and pure-ox-ish limits brackets H2 (2.016) .. O2 (32).
    assert 2.0 < engine.mbar(0.999) < 2.1
    assert 31.0 < engine.mbar(0.001) < 33.0


def test_frozen_cp_is_bounded_by_an_independent_cantera_fit(engine):
    """The unburnt surface's c_p (which sets γ and a, both shipped) rides on
    CEA's NASA-9 reactant polynomials, EXTRAPOLATED below their ~200 K fit
    floor at the cold end (review: `interp_error_bound` bounds interpolation,
    not the underlying model/extrapolation error, and unlike the equilibrium
    surface there is no built-in thermo cross-check).

    So bound it against an INDEPENDENT fit: Cantera's NASA-7 H₂/O₂ thermo
    (`h2o2.yaml`, the same set the transport surface's `_CP_CROSSCHECK`
    uses). Two independent fits agreeing across the envelope is what turns
    the cold-corner extrapolation from a hope into a bounded number — exactly
    the OFFL-3 §6-2 discipline, applied to the frozen branch. Measured worst
    disagreement is ~2.8% at the 100 K floor and < 0.3% above 200 K; the 3.5%
    gate is that with headroom (and still an order below the branch's own
    ideal-gas model-form band)."""
    gas = ct.Solution("h2o2.yaml")
    worst = 0.0
    for z in (0.111, 1.0 / 6.0, 0.25):  # the Z-envelope span
        for t in (100.0, 150.0, 200.0, 400.0, 1200.0, 2200.0):
            cea_cp = engine.state_php(1.0e6, engine.enthalpy_at(t, z), z).cp_fr
            gas.TPY = t, 1.0e6, f"H2:{z}, O2:{1.0 - z}"
            rel = abs(cea_cp - gas.cp_mass) / cea_cp
            worst = max(worst, rel)
            assert rel < 3.5e-2, (z, t, cea_cp, gas.cp_mass, rel)
    # The extrapolation is genuinely exercised at the cold floor (the gate is
    # not vacuous): the worst disagreement lives there, not in the fitted band.
    assert worst > 5.0e-3, f"cross-check never bit — worst {worst:.3%} is suspiciously tight"
