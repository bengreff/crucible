"""OFFL-5 §6 — the spine's chemical-regime transport surface
(`tables/spine/lox_lh2_transport_v0.2.0.h5`), the artifact the RL10 configs
pin for per-cell μ/k/c_p/c_v/∂h∂Z.

Same fresh-holdout discipline as the station-5 equilibrium gate: ⅜-offset
points disjoint from the generator's own midpoint + ¼-offset measurement
set, plus the declared envelope-EDGE points unsubsampled (the session-12
finding — subsampling those away reopens the hole they closed), scored
against the stored `interp_error_bound` and, for the five log-valued
columns, the stored rule-space `interp_error_bound_log`.

Two invariants are specific to this surface and are the point of the file:

- **Envelope identity with the equilibrium surface.** OFFL-5 §3.1a requires
  the two runtime surfaces a cell interrogates to refuse and accept on the
  same set of states. COUP-8 §3.3(2) checks each table against its
  consumers, one at a time — it cannot see a gap *between* two independently
  declared tables. So the gap is closed here instead.
- **The interpolation bound sits well inside the declared physical band.**
  FND-5 §3.4's grid-sizing rule is only met if the measured bound is small
  against the 10–20% band FND-7 §3.3 declares on this stage; a table whose
  interpolation error rivalled its physics band would be a false economy.
"""

import pathlib

import h5py
import numpy as np
import pytest

from crucible_offl.chemistry import EquilibriumEngine
from crucible_offl.surfaces import _edge_points, _multilinear, station5_envelope_grid
from crucible_offl.transport import _TR_COLUMNS, _TR_RULES, TransportEvaluator

REPO = pathlib.Path(__file__).resolve().parents[2]
PRODUCTION = REPO / "tables" / "spine" / "lox_lh2_transport_v0.2.0.h5"
GROUP = "/spine/lox_lh2/transport"
#: FND-7 §3.3's declared band on this stage; the interpolation bound must
#: stay a factor below it (module doc).
DECLARED_BAND = 0.10


@pytest.fixture(scope="module")
def evaluator():
    # gas-only metastable, exactly as the surface was generated — a
    # full-condensed truth engine would score against different physics.
    return TransportEvaluator(EquilibriumEngine(gas_only=True))


@pytest.fixture(scope="module")
def surface():
    with h5py.File(PRODUCTION, "r") as f:
        g = f[GROUP]
        axes = [np.array(g["axes"][n]) for n in ("p", "h", "Z")]
        env = [
            (float(g["axes"][n].attrs["envelope_min"]), float(g["axes"][n].attrs["envelope_max"]))
            for n in ("p", "h", "Z")
        ]
        shape = tuple(len(a) for a in axes)
        grids, bounds, bounds_log, units, rules = {}, {}, {}, {}, {}
        for name in g["values"]:
            if name.startswith("sigma_"):
                continue
            ds = g["values"][name]
            grids[name] = np.array(ds).reshape(shape)
            bounds[name] = float(ds.attrs["interp_error_bound"])
            units[name] = str(ds.attrs["units"])
            rules[name] = str(ds.attrs["interp_rule"])
            if "interp_error_bound_log" in ds.attrs:
                bounds_log[name] = float(ds.attrs["interp_error_bound_log"])
    return axes, env, grids, bounds, bounds_log, units, rules


def test_columns_units_and_rules_are_exactly_the_declared_set(surface):
    _axes, _env, grids, _b, _bl, units, rules = surface
    assert set(grids) == set(_TR_COLUMNS), "shipped columns must be the declared set"
    for name, (rule, unit) in _TR_COLUMNS.items():
        # The units string crosses the language seam as data and is gated
        # once by Rust's `Table::bind` — a drift here is a load refusal.
        assert units[name] == unit
        assert rules[name] == rule


def test_envelope_is_identical_to_the_equilibrium_surface(surface):
    """OFFL-5 §3.1a: one coordinate, one pair of surfaces, one envelope."""
    _axes, env, *_ = surface
    eq = station5_envelope_grid()
    for got, want, name in zip(
        env, (eq.p_envelope, eq.h_envelope, eq.z_envelope), ("p", "h", "Z")
    ):
        assert got == pytest.approx(want, rel=0.0, abs=0.0), (
            f"{name} envelope {got} != equilibrium surface's {want} — the two "
            "runtime surfaces would refuse on different state sets"
        )


def test_axes_span_exactly_the_envelope(surface):
    """The grid IS the envelope (transport.py `matched_transport_grid`), so
    no node is a state the runtime may not legally interrogate — and an
    envelope-boundary query is still in-domain."""
    axes, env, *_ = surface
    for ax, (lo, hi) in zip(axes, env):
        assert ax[0] == pytest.approx(lo, rel=1e-12)
        assert ax[-1] == pytest.approx(hi, rel=1e-12)


def test_offl5_6_fresh_holdout_respects_stored_bounds(evaluator, surface):
    axes, env, grids, bounds, bounds_log, _u, _r = surface
    p_ax, h_ax, z_ax = axes
    hold_p = np.concatenate(
        [
            np.exp(0.625 * np.log(p_ax[:-1]) + 0.375 * np.log(p_ax[1:]))[::4],
            _edge_points(p_ax, True, env[0]),
        ]
    )
    hold_h = np.concatenate(
        [(0.625 * h_ax[:-1] + 0.375 * h_ax[1:])[::4], _edge_points(h_ax, False, env[1])]
    )
    hold_z = np.concatenate(
        [(0.625 * z_ax[:-1] + 0.375 * z_ax[1:])[::2], _edge_points(z_ax, False, env[2])]
    )

    n = 0
    max_log = {name: 0.0 for name in bounds_log}
    for p in hold_p:
        for h in hold_h:
            for z in hold_z:
                if not all(lo <= q <= hi for (lo, hi), q in zip(env, (p, h, z))):
                    continue
                n += 1
                truth = evaluator.columns(p, h, z)
                for name, t in truth.items():
                    est = _multilinear(_TR_RULES[name], axes, grids, (p, h, z), name)
                    err = abs(est - t)
                    assert err <= bounds[name], (
                        f"{name} at ({p:.4g}, {h:.4g}, {z:.4g}): "
                        f"|{est:.6g} - {t:.6g}| > stored bound {bounds[name]:.6g}"
                    )
                    if name in bounds_log:
                        log_err = abs(np.log(est) - np.log(t))
                        max_log[name] = max(max_log[name], log_err)
                        assert log_err <= bounds_log[name], (
                            f"{name} at ({p:.4g}, {h:.4g}, {z:.4g}): log error "
                            f"{log_err:.6g} > stored log bound {bounds_log[name]:.6g}"
                        )
    assert n > 100
    assert all(b > 0.0 for b in bounds.values()), "every column must carry a measured bound"
    # The five positive columns are log-valued and must have been scored in
    # rule space; dh_dz is signed by nature and is linear-only.
    assert set(bounds_log) == {n for n, r in _TR_RULES.items() if r.endswith("log")}
    assert all(m > 0.0 for m in max_log.values())


def test_interpolation_bound_sits_inside_the_declared_physical_band(surface):
    _axes, _env, grids, bounds, bounds_log, _u, _r = surface
    for name, b in bounds.items():
        scale = float(np.abs(grids[name]).max())
        assert b / scale < DECLARED_BAND, (
            f"{name}: absolute bound {b:.4g} is {b / scale:.1%} of the column's "
            f"peak — not comfortably inside the declared {DECLARED_BAND:.0%} "
            "physics band; refine the grid (never relax this gate)"
        )
    for name, b in bounds_log.items():
        assert b < DECLARED_BAND, (
            f"{name}: rule-space bound {b:.2%} rivals the declared "
            f"{DECLARED_BAND:.0%} physics band; refine the grid"
        )


def test_cantera_cea_frozen_cp_crosscheck_is_armed(evaluator):
    """The generator's per-node two-fit gate (`_CP_CROSSCHECK`) is what
    bounds Cantera's above-3500 K thermo extrapolation and catches a
    mis-mapped species. Prove it actually fires rather than sitting dead:
    a state whose composition has been corrupted must be refused."""
    from crucible_offl.transport import _CP_CROSSCHECK

    state = evaluator.engine.state_php(3.0e6, -1.083e6, 1.0 / 6.0)
    evaluator._cantera_transport(state)  # the honest state passes
    poisoned = type(state)(
        **{**state.__dict__, "cp_fr": state.cp_fr * (1.0 + 3.0 * _CP_CROSSCHECK)}
    )
    with pytest.raises(RuntimeError, match="frozen c_p disagree"):
        evaluator._cantera_transport(poisoned)
