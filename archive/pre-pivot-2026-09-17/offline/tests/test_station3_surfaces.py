"""OFFL-3 §6-4/§6-5/§6-6 — the tabulated products themselves.

§6-4 runs against the *committed production artifact*: fresh holdout points
(⅜-offsets, disjoint from the generator's midpoint + ¼-offset measurement
set) must respect the stored `interp_error_bound` of every column. A
violation means the production grid needs refining — the gate is never
relaxed (surfaces.SAFETY doc).

§6-5: regeneration determinism — two builds → identical content digests
(Tier-2, META-1 §2.1) — and provenance completeness on the *committed*
production artifact.

§6-6: coordinate consistency (S22) on the committed artifact — the (p,h,Z)
surface and the (p_c, MR) performance reference describe the same chamber
along the design line h = h_inj(Z): their temperatures must agree within
the two stored bounds, and both must sit on the direct CEA solution in the
gas region.
"""

import pathlib

import h5py
import numpy as np
import pytest

from crucible_offl.chemistry import EquilibriumEngine, mr_to_z
from crucible_offl.surfaces import (
    _EQ_RULES,
    EquilibriumGrid,
    _eq_columns,
    _multilinear,
    build_equilibrium_surface,
)
from crucible_offl.tables import spec_digest

REPO = pathlib.Path(__file__).resolve().parents[2]
PRODUCTION = REPO / "tables" / "chem" / "lox_lh2_v0.1.0.h5"

CI_GRID = EquilibriumGrid(
    p_points=tuple(np.geomspace(1.5e3, 8.0e6, 15)),
    h_points=tuple(np.linspace(-1.18e7, -1.0e5, 11)),
    z_points=tuple(np.linspace(0.10, 0.26, 5)),
    p_envelope=(2.0e3, 7.0e6),
    h_envelope=(-1.15e7, -2.0e5),
    z_envelope=(1.0 / 9.0, 0.25),
)


@pytest.fixture(scope="module")
def engine():
    return EquilibriumEngine()


@pytest.fixture(scope="module")
def ci_surface(engine):
    return build_equilibrium_surface(engine, CI_GRID, "ci", "ci")


def test_s6_4_fresh_holdout_respects_stored_bounds_on_production(engine):
    with h5py.File(PRODUCTION, "r") as f:
        g = f["/chem/lox_lh2/equilibrium"]
        axes = [np.array(g["axes"][n]) for n in ("p", "h", "Z")]
        shape = tuple(len(a) for a in axes)
        grids = {}
        bounds = {}
        env = []
        for n in ("p", "h", "Z"):
            ds = g["axes"][n]
            env.append((float(ds.attrs["envelope_min"]), float(ds.attrs["envelope_max"])))
        for name in g["values"]:
            if name.startswith("sigma_"):
                continue
            grids[name] = np.array(g["values"][name]).reshape(shape)
            bounds[name] = float(g["values"][name].attrs["interp_error_bound"])
    p_ax, h_ax, z_ax = axes
    n = 0
    for p in np.exp(0.625 * np.log(p_ax[:-1]) + 0.375 * np.log(p_ax[1:]))[::3]:
        for h in (0.625 * h_ax[:-1] + 0.375 * h_ax[1:])[::3]:
            for z in (0.625 * z_ax[:-1] + 0.375 * z_ax[1:])[::2]:
                if not all(lo <= q <= hi for (lo, hi), q in zip(env, (p, h, z))):
                    continue
                n += 1
                truth = _eq_columns(engine.state_php(p, h, z))
                for name, t in truth.items():
                    est = _multilinear(_EQ_RULES[name], axes, grids, (p, h, z), name)
                    assert abs(est - t) <= bounds[name], (
                        f"{name} at ({p:.4g}, {h:.4g}, {z:.4g}): "
                        f"|{est:.6g} - {t:.6g}| > stored bound {bounds[name]:.6g}"
                    )
    assert n > 100
    assert all(b > 0.0 for b in bounds.values())


def test_s6_5_regeneration_is_digest_identical(engine, ci_surface):
    spec, _ = ci_surface
    spec2, _ = build_equilibrium_surface(engine, CI_GRID, "ci", "ci")
    assert spec_digest(spec) == spec_digest(spec2)


def test_s6_5_committed_artifact_carries_full_provenance():
    with h5py.File(PRODUCTION, "r") as f:
        for group in ("/chem/lox_lh2/equilibrium", "/chem/lox_lh2/performance"):
            g = f[group]
            for attr in (
                "kind",
                "schema_version",
                "data_version",
                "interp_method",
                "producer",
                "producer_version",
                "input_deck_hash",
                "source_library",
                "generator_commit",
                "axis_order",
            ):
                assert attr in g.attrs, f"{group} missing {attr}"
            assert g.attrs["input_deck_hash"].startswith("sha256:")
            commit = g.attrs["generator_commit"].removesuffix("-dirty")
            assert len(commit) == 40 and all(c in "0123456789abcdef" for c in commit)


def test_s6_6_coordinate_consistency_on_the_design_line(engine):
    with h5py.File(PRODUCTION, "r") as f:
        eq = f["/chem/lox_lh2/equilibrium"]
        perf = f["/chem/lox_lh2/performance"]
        eq_axes = [np.array(eq["axes"][n]) for n in ("p", "h", "Z")]
        eq_shape = tuple(len(a) for a in eq_axes)
        t_grid = {"temperature": np.array(eq["values/temperature"]).reshape(eq_shape)}
        t_bound = float(eq["values/temperature"].attrs["interp_error_bound"])
        p_axes = [np.array(perf["axes"][n]) for n in ("p_c", "MR")]
        p_shape = tuple(len(a) for a in p_axes)
        tc_grid = {"T_c": np.array(perf["values/T_c"]).reshape(p_shape)}
        tc_bound = float(perf["values/T_c"].attrs["interp_error_bound"])

    for pc, mr in ((32.75e5, 5.0), (20e5, 4.5), (45e5, 6.0), (15e5, 3.8), (50e5, 7.0)):
        z = mr_to_z(mr)
        h = engine.injection_enthalpy(z)
        t_surf = _multilinear(
            _EQ_RULES["temperature"], eq_axes, t_grid, (pc, h, z), "temperature"
        )
        t_perf = _multilinear("lin-lin-lin", p_axes, tc_grid, (pc, mr), "T_c")
        t_true = engine.chamber(pc, mr).T
        # The two parameterizations agree within their stored bounds…
        assert abs(t_surf - t_perf) <= t_bound + tc_bound
        # …and in the gas region both sit on the direct CEA solution far
        # inside the (kink-dominated) stored bound: measured ≤ 2.6 K.
        assert abs(t_surf - t_true) <= 4.0
        assert abs(t_perf - t_true) <= 7.0
