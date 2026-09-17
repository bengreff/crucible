"""OFFL-3 §3.3 / VAL-2 §3.3 (plan S6) — the IGNITION closure production
artifact `tables/chem/lox_lh2_ignition_v0.3.0.h5`: the SOLV-4 §3.6 rate law's
`S_L(p, T_u, Z)` and `τ_ign(p, T_u, Z)`, from Cantera 1-D flames + 0-D
reactors on `h2o2.yaml`.

Gates, in the order they matter:

- **§6-4 active-region fresh holdout:** ⅜-offset points disjoint from the
  generator's set respect the stored `interp_error_bound` **in each column's
  active region** (`S_L > 0`, `τ_ign < τ_max`) — the extinction transitions
  are the declared sharp feature (each mechanism governs where the other is
  extinct), not an interpolation the bound must cover.
- **Physicality:** `S_L ∈ [0, 300]` m/s (0 = extinction; the ceiling is the
  auto-ignition regime), `τ_ign ∈ (0, τ_max]` s, all finite.
- **VAL-2 `h2-flame-speed` anchor:** the mechanism reproduces the measured
  stoichiometric + near-peak H₂/air laminar flame speeds within a declared
  model-form band (the surface is on H₂/O₂, so the anchor is a *mechanism*
  check via a fresh H₂/air solve).
- **VAL-2 `h2-ignition-delay` anchor:** the mechanism's shock-tube-class τ_ign
  has the right magnitude (10⁻⁵–10⁻³ s near 1000–1300 K) and the correct
  Arrhenius sign (τ falls as T rises).
"""

import pathlib

import numpy as np
import pytest

import cantera as ct

from crucible_offl.ignition import (
    SL_CEILING_M_PER_S,
    SL_FLOOR_M_PER_S,
    TAU_MAX,
    IgnitionEngine,
    _IGNITION_RULES,
)
from crucible_offl.surfaces import _edge_points, _multilinear

REPO = pathlib.Path(__file__).resolve().parents[2]
PRODUCTION = REPO / "tables" / "chem" / "lox_lh2_ignition_v0.3.0.h5"
GROUP = "/chem/lox_lh2/ignition"


@pytest.fixture(scope="module")
def engine():
    return IgnitionEngine()


def _active(name: str, value: float) -> bool:
    if name == "laminar_flame_speed":
        return value > SL_FLOOR_M_PER_S
    return value < TAU_MAX * (1.0 - 1.0e-9)


def _load():
    import h5py

    with h5py.File(PRODUCTION, "r") as f:
        g = f[GROUP]
        axes = [np.array(g["axes"][n]) for n in ("p", "T_u", "Z")]
        env = [
            (float(g["axes"][n].attrs["envelope_min"]), float(g["axes"][n].attrs["envelope_max"]))
            for n in ("p", "T_u", "Z")
        ]
        shape = tuple(len(a) for a in axes)
        grids, bounds = {}, {}
        for name in _IGNITION_RULES:
            ds = g["values"][name]
            grids[name] = np.array(ds).reshape(shape)
            bounds[name] = float(ds.attrs["interp_error_bound"])
    return axes, env, grids, bounds


def test_s6_4_active_region_fresh_holdout(engine):
    axes, env, grids, bounds = _load()
    (plo, phi), (tlo, thi), (zlo, zhi) = env
    # Fresh ⅜-offset points (disjoint from the generator's ½/¼ set), heavily
    # subsampled (each is a flame + reactor solve), + the envelope edges.
    hp = np.concatenate(
        [np.exp(0.625 * np.log(axes[0][:-1]) + 0.375 * np.log(axes[0][1:]))[::2],
         _edge_points(axes[0], True, env[0])]
    )
    ht = np.concatenate([(0.625 * axes[1][:-1] + 0.375 * axes[1][1:])[::3], _edge_points(axes[1], False, env[1])])
    hz = np.concatenate([(0.625 * axes[2][:-1] + 0.375 * axes[2][1:])[::2], _edge_points(axes[2], False, env[2])])
    # The stored bound's DECLARED DOMAIN is fully-active cells only (the
    # builder's active-region rule: a cell straddling an extinction
    # transition — an S_L→0 crossover or a τ_ign cap corner — is the
    # declared sharp model feature, not an interpolation the bound covers).
    # This test must therefore apply the builder's own 8-corner criterion,
    # not a point-activity filter: a ⅜-point in a straddling cell can have
    # active truth AND active estimate while the estimate is polluted by the
    # capped corner (measured on the 0.2.0 extended surface: τ est 3.1e-3 vs
    # truth 8.1e-6 at (6.75 kPa, 2894 K) — a cap-straddling cell).
    masks = {
        name: np.vectorize(lambda v, n=name: 1.0 if _active(n, v) else 0.0)(grids[name])
        for name in _IGNITION_RULES
    }
    mask_rules = {name: rule.rsplit("-", 1)[0] + "-lin" for name, rule in _IGNITION_RULES.items()}
    n = {name: 0 for name in _IGNITION_RULES}
    for p in hp:
        if not (plo <= p <= phi):
            continue
        for t_u in ht:
            if not (tlo <= t_u <= thi):
                continue
            for z in hz:
                if not (zlo <= z <= zhi):
                    continue
                q = (float(p), float(t_u), float(z))
                truth = engine.columns(*q)
                for name, rule in _IGNITION_RULES.items():
                    if _multilinear(mask_rules[name], axes, masks, q, name) < 1.0 - 1.0e-12:
                        continue  # straddling cell — outside the bound's domain
                    est = _multilinear(rule, axes, grids, q, name)
                    t = truth[name]
                    n[name] += 1
                    assert abs(est - t) <= bounds[name], (
                        f"{name} holdout error {abs(est - t):.4g} exceeds stored bound "
                        f"{bounds[name]:.4g} at (p={p:.3g}, T_u={t_u:.1f}, Z={z:.3f})"
                    )
    for name, count in n.items():
        assert count > 5, f"{name}: only {count} in-active-region holdout points"


def test_surface_is_physical_everywhere():
    _axes, _env, grids, _bounds = _load()
    s_l = grids["laminar_flame_speed"]
    tau = grids["ignition_delay"]
    assert np.all(np.isfinite(s_l)) and np.all(np.isfinite(tau))
    assert np.all(s_l >= 0.0) and np.all(s_l <= SL_CEILING_M_PER_S)
    assert np.all(tau > 0.0) and np.all(tau <= TAU_MAX * (1.0 + 1.0e-12))
    # There is a genuinely flammable region (some S_L well above the floor)
    # and a genuinely ignitable region (some τ well below the cap).
    assert np.any(s_l > 1.0), "no flammable region"
    assert np.any(tau < 1.0e-3), "no fast-ignition region"


def test_h2_air_flame_speed_anchor():
    """VAL-2 §3.3 `h2-flame-speed`: the mechanism reproduces the measured
    H₂/air laminar flame speed. Cached reference bands (META-3 §6.9): the
    stoichiometric value ≈ 2.1 m/s and the near-peak value ≈ 3.0 m/s
    (298 K, 1 atm), robust across the standard compilations (Law; Egolfopoulos
    & Law; Verhelst & Wallner). The declared model-form band is generous — the
    claim is that the mechanism is a real H₂ mechanism, not that it is tuned."""
    def sl_h2_air(phi: float) -> float:
        gas = ct.Solution("h2o2.yaml", transport_model="mixture-averaged")
        # 2 H2 + O2 stoichiometric; air = O2 + 3.76 N2. φ scales the fuel.
        gas.set_equivalence_ratio(phi, "H2", {"O2": 1.0, "N2": 3.76})
        gas.TP = 298.0, ct.one_atm
        f = ct.FreeFlame(gas, width=0.03)
        f.set_refine_criteria(ratio=3.0, slope=0.06, curve=0.06)
        f.transport_model = "mixture-averaged"
        f.solve(loglevel=0, auto=True)
        return float(f.velocity[0])

    s_stoich = sl_h2_air(1.0)
    s_peak = sl_h2_air(1.7)
    # Measured: φ=1 ≈ 2.1 (1.9–2.4 across studies); φ≈1.7 ≈ 3.0 (2.9–3.5).
    # Declared model-form bands (±~25%):
    assert 1.6 <= s_stoich <= 2.7, f"H2/air stoich S_L = {s_stoich:.3f} m/s off the anchor band"
    assert 2.4 <= s_peak <= 3.8, f"H2/air peak S_L = {s_peak:.3f} m/s off the anchor band"
    assert s_peak > s_stoich, "the flame speed must peak fuel-rich (φ>1)"


def test_ignition_delay_anchor(engine):
    """VAL-2 §3.3 `h2-ignition-delay`: shock-tube-class magnitude + Arrhenius
    trend for stoichiometric H₂/O₂ (Z = 1/9). τ_ign ~ 10⁻⁵–10⁻³ s near
    1000–1300 K and falls monotonically as T rises."""
    z = 1.0 / 9.0  # stoichiometric H2/O2
    p = 1.0e5
    taus = [engine.ignition_delay(p, t, z) for t in (1000.0, 1100.0, 1200.0, 1300.0)]
    # Right magnitude at the hot end:
    assert 1.0e-6 <= taus[-1] <= 1.0e-3, f"τ_ign(1300 K) = {taus[-1]:.3e} s off the shock-tube band"
    # Strictly Arrhenius (falls with T):
    assert all(taus[i] > taus[i + 1] for i in range(len(taus) - 1)), (
        f"τ_ign must fall as T rises (got {['%.2e' % t for t in taus]})"
    )
