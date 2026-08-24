"""META-1 §2.1 Tier-2 / FND-4 §6-4 — cross-process regeneration determinism.

The [tables] pin chain (config → sidecar digest → FND-5 open gate) is only
sound if regenerating a table from its pinned pipeline reproduces the SAME
canonical content digest in a fresh process: cold-start CEA + numpy, no
shared interpreter state, no wall-clock, no ambient environment influence.

Session-10's review exercised this once by hand on the production artifact;
this test pins it in CI on a small probe surface (3×3×3 — 27 HP solves per
process, seconds) driven through the identical build path
(`build_equilibrium_surface` → `spec_digest`, digest v3 — THE cross-language
contract). Two subprocesses must print byte-identical digests.
"""

from __future__ import annotations

import subprocess
import sys

PROBE = r"""
import numpy as np
from crucible_offl.chemistry import EquilibriumEngine
from crucible_offl.surfaces import EquilibriumGrid, build_equilibrium_surface
from crucible_offl.tables import spec_digest

grid = EquilibriumGrid(
    p_points=tuple(np.geomspace(2.0e6, 4.0e6, 3)),
    h_points=(-2.0e6, -1.5e6, -1.0e6),
    z_points=(0.15, 1.0 / 6.0, 0.18),
    p_envelope=(2.2e6, 3.8e6),
    h_envelope=(-1.9e6, -1.1e6),
    z_envelope=(0.155, 0.175),
)
spec, _bounds = build_equilibrium_surface(
    EquilibriumEngine(), grid, data_version="probe-0.0.0", generator_commit="regen-probe"
)
print(spec_digest(spec))
"""


# The unburnt-reactant build path (plan S5) is a DIFFERENT pipeline —
# `FrozenReactantEngine` with a fixed-count bisection T-inversion and a
# CEA reactant `calc_property` — so it earns its own cross-process probe:
# the bisection and the derivative step are fixed/declared, but only a
# two-subprocess digest match proves the whole path is bit-reproducible.
UNBURNT_PROBE = r"""
from crucible_offl.chemistry import FrozenReactantEngine
from crucible_offl.surfaces import unburnt_reactant_grid, build_unburnt_surface
from crucible_offl.tables import spec_digest

eng = FrozenReactantEngine()
grid = unburnt_reactant_grid(eng, n_p=3, n_h=5, n_z=3)
spec, _bounds = build_unburnt_surface(
    eng, grid, data_version="probe-0.0.0", generator_commit="regen-probe"
)
print(spec_digest(spec))
"""


# The ignition build path (plan S6) is a THIRD pipeline — Cantera FreeFlame
# (S_L) + IdealGasConstPressureReactor (τ_ign) on h2o2.yaml — with its own
# fixed solver settings; a two-subprocess digest match proves the flame +
# reactor solves are bit-reproducible cross-process (a tiny 2×2×2 probe — 8
# flame+reactor grid solves plus its holdouts per process, ~seconds;
# n_tu_ext=0 keeps the probe's T_u axis 2 nodes, not the production grid's
# appended auto-ignitive rows). The probe's (p, T_u, Z) box sits inside the
# MEASURED dual-active band — every corner both flammable (S_L > floor) and
# auto-ignitive within the reactor horizon (tau < cap): 50-100 bar x
# 1000-1100 K x Z 0.10-0.20, all 8 corners measured OK — because the
# builder's active-region-holdout refusal (production-grid QA) requires each
# column a fully-active cell, and the dual-active region is a thin diagonal
# band in (p, T_u) a naive probe box misses (at 1 atm these temperatures
# already run away past the flame ceiling; at 30 mbar the flame solve dies;
# at 150-1900 K one corner is always dead).
IGNITION_PROBE = r"""
from crucible_offl.ignition import IgnitionEngine, ignition_grid, build_ignition_surface
from crucible_offl.tables import spec_digest

eng = IgnitionEngine()
grid = ignition_grid(
    n_p=2, n_tu=2, n_z=2,
    p_lo=5.0e6, p_hi=1.0e7, tu_lo=1000.0, tu_hi=1100.0, z_lo=0.10, z_hi=0.20,
    n_tu_ext=0,
)
spec, _bounds = build_ignition_surface(
    eng, grid, data_version="probe-0.0.0", generator_commit="regen-probe"
)
print(spec_digest(spec))
"""


def _run(probe: str) -> str:
    out = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=True,
    )
    digest = out.stdout.strip()
    assert digest.startswith("sha256:"), f"probe emitted {digest!r}"
    return digest


def test_fresh_process_regeneration_reproduces_the_digest() -> None:
    first = _run(PROBE)
    second = _run(PROBE)
    assert first == second, (
        "two cold-start generations of the same probe surface disagree — "
        "the regeneration key (config + pinned pipeline) cannot replay "
        f"({first} vs {second})"
    )


def test_unburnt_build_path_regenerates_deterministically() -> None:
    first = _run(UNBURNT_PROBE)
    second = _run(UNBURNT_PROBE)
    assert first == second, (
        "two cold-start generations of the unburnt probe surface disagree — "
        f"the frozen-reactant build path is not bit-reproducible ({first} vs {second})"
    )


def test_ignition_build_path_regenerates_deterministically() -> None:
    first = _run(IGNITION_PROBE)
    second = _run(IGNITION_PROBE)
    assert first == second, (
        "two cold-start generations of the ignition probe surface disagree — "
        f"the Cantera flame+reactor build path is not bit-reproducible ({first} vs {second})"
    )
