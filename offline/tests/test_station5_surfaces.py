"""OFFL-3 §6-4 — the STATION-5 production artifact (session-12 review
finding): `tables/chem/lox_lh2_v0.3.2.h5` is the artifact the RL10 engine
configs actually pin (`configs/rl10_coarse.toml` et al., group
`/chem/lox_lh2/equilibrium`), but only the station-3 design-window artifact
(`lox_lh2_v0.1.0.h5`, see test_station3_surfaces.py) carried a fresh-holdout
gate. This module closes that gap on the artifact that matters for the
blind RL10.

Fresh ⅜-offset holdout points — disjoint from the generator's midpoint + ¼-
offset measurement set — must respect the stored `interp_error_bound` of
every column, exactly as station-3's §6-4 gate requires. Two things are
specific to station-5:

- The surface was generated GAS-ONLY (`station5_envelope_grid`, metastable
  equilibrium, see scripts/make_station5_tables.py): the truth engine here
  must be built `EquilibriumEngine(gas_only=True)` or every comparison is
  scored against the wrong physics.
- The density column carries a stored `interp_error_bound_log` (it is
  log-valued, `interp_rule` ending in "log"): fresh holdout must also
  respect that bound in rule (log) space, not just the linear one.

The holdout set additionally folds in the declared envelope-EDGE points
(`crucible_offl.surfaces._edge_points`, the session-12 fix) unsubsampled —
those are exactly the previously-unreachable-by-midpoint/quarter-offset
slabs the review finding called out, so subsampling them away here would
reopen the same hole.
"""

import pathlib

import h5py
import numpy as np
import pytest

from crucible_offl.chemistry import EquilibriumEngine
from crucible_offl.surfaces import _EQ_RULES, _edge_points, _eq_columns, _multilinear

REPO = pathlib.Path(__file__).resolve().parents[2]
PRODUCTION = REPO / "tables" / "chem" / "lox_lh2_v0.3.2.h5"


@pytest.fixture(scope="module")
def engine():
    # station5_envelope_grid tables are generated gas-only (metastable
    # equilibrium) — match the generation physics or the "truth" solves
    # are for a different product set entirely.
    return EquilibriumEngine(gas_only=True)


def test_s6_4_fresh_holdout_respects_stored_bounds_on_station5(engine):
    with h5py.File(PRODUCTION, "r") as f:
        g = f["/chem/lox_lh2/equilibrium"]
        axes = [np.array(g["axes"][n]) for n in ("p", "h", "Z")]
        shape = tuple(len(a) for a in axes)
        grids = {}
        bounds = {}
        bounds_log = {}
        env = []
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

    # Fresh ⅜-offset points (station-3's stride convention: p//3, h//3,
    # z//2), PLUS the envelope-edge points unsubsampled.
    hold_p = np.concatenate(
        [
            np.exp(0.625 * np.log(p_ax[:-1]) + 0.375 * np.log(p_ax[1:]))[::3],
            _edge_points(p_ax, True, env[0]),
        ]
    )
    hold_h = np.concatenate(
        [
            (0.625 * h_ax[:-1] + 0.375 * h_ax[1:])[::3],
            _edge_points(h_ax, False, env[1]),
        ]
    )
    hold_z = np.concatenate(
        [
            (0.625 * z_ax[:-1] + 0.375 * z_ax[1:])[::2],
            _edge_points(z_ax, False, env[2]),
        ]
    )

    n = 0
    max_abs = {name: 0.0 for name in bounds}
    max_log = {name: 0.0 for name in bounds_log}
    for p in hold_p:
        for h in hold_h:
            for z in hold_z:
                if not all(lo <= q <= hi for (lo, hi), q in zip(env, (p, h, z))):
                    continue
                n += 1
                truth = _eq_columns(engine.state_php(p, h, z))
                for name, t in truth.items():
                    est = _multilinear(_EQ_RULES[name], axes, grids, (p, h, z), name)
                    err = abs(est - t)
                    max_abs[name] = max(max_abs[name], err)
                    assert err <= bounds[name], (
                        f"{name} at ({p:.4g}, {h:.4g}, {z:.4g}): "
                        f"|{est:.6g} - {t:.6g}| > stored bound {bounds[name]:.6g}"
                    )
                    if name in bounds_log:
                        log_err = abs(np.log(est) - np.log(t))
                        max_log[name] = max(max_log[name], log_err)
                        assert log_err <= bounds_log[name], (
                            f"{name} at ({p:.4g}, {h:.4g}, {z:.4g}): "
                            f"log error |ln {est:.6g} - ln {t:.6g}| = "
                            f"{log_err:.6g} > stored log bound {bounds_log[name]:.6g}"
                        )
    assert n > 100
    # condensed_fraction is identically 0 on this gas-only surface (declared
    # model form, see station5_envelope_grid's docstring) so its measured —
    # and stored — bound is legitimately 0.0; every other column must carry
    # a strictly positive measured bound.
    assert all(b > 0.0 for name, b in bounds.items() if name != "condensed_fraction")
    assert bounds["condensed_fraction"] == 0.0
    # density must have been checked in log space at least once.
    assert bounds_log
    assert all(m > 0.0 for m in max_log.values())


def test_density_carries_finite_positive_log_bound():
    with h5py.File(PRODUCTION, "r") as f:
        ds = f["/chem/lox_lh2/equilibrium/values/density"]
        assert "interp_error_bound_log" in ds.attrs
        bound_log = float(ds.attrs["interp_error_bound_log"])
    assert np.isfinite(bound_log)
    assert bound_log > 0.0
