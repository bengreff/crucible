"""OFFL-3 §6-1 — RP-1311 reproduction through *our* pipeline primitives.

Reference values are the published output of RP-1311 example 8 (LOX/LH2 IAC
rocket, p_c = 53.3172 bar, o/f = 5.55157; McBride & Gordon, NASA RP-1311,
NTRS 19960044559) as printed in the NASA CEA release docs
(nasa.github.io/cea, examples/rocket/example8) — NASA's own regression
anchor for the reimplementation. Tolerances are the printed precision of
the reference, not a fit margin: the pipeline must drive CEA (weights,
injection enthalpy, HP/rocket problems) exactly as the manual case does.

§6-3 (frozen/shifting bracket) rides here too: the same case's shifting vs
frozen-from-chamber vacuum Isp must strictly bracket, and the JANNAF
kinetic-efficiency band (~0.8–1% of shifting, META-3 `jannaf-eff`) must
sit inside the measured gap — the gap *is* the model-form error bar (S3).
"""

import pytest

from crucible_offl.chemistry import EquilibriumEngine, mr_to_z

PC = 53.3172e5  # Pa
MR = 5.55157

# Published example-8 output, chamber column.
CHAMBER = {
    "T": 3383.845,
    "rho": 2.410,
    "h": -1026.05e3,
    "s": 18.659e3,
    "mbar": 12.716,
    "gamma": 1.145,
    "a": 1591.47,
}
CHAMBER_X = {
    "H2O": 0.63456,
    "H2": 0.29479,
    "H": 0.033498,
    "OH": 0.033341,
    "O": 0.0020678,
    "O2": 0.0017218,
}
C_STAR = 2332.34
T_THROAT = 3185.673
ISP_VAC_25 = 4348.510  # m/s at Ae/At = 25, shifting


@pytest.fixture(scope="module")
def engine():
    return EquilibriumEngine()


@pytest.fixture(scope="module")
def chamber(engine):
    return engine.chamber(PC, MR)


def test_chamber_state_reproduces_manual(chamber):
    assert chamber.T == pytest.approx(CHAMBER["T"], abs=5e-3)
    assert chamber.rho == pytest.approx(CHAMBER["rho"], abs=5e-4)
    assert chamber.h == pytest.approx(CHAMBER["h"], abs=5.0)
    assert chamber.s == pytest.approx(CHAMBER["s"], abs=0.5)
    assert chamber.mbar == pytest.approx(CHAMBER["mbar"], abs=5e-4)
    assert chamber.gamma == pytest.approx(CHAMBER["gamma"], abs=5e-4)
    assert chamber.a == pytest.approx(CHAMBER["a"], abs=5e-3)


def test_chamber_composition_reproduces_manual(chamber):
    for sp, x_ref in CHAMBER_X.items():
        assert chamber.mole_fractions[sp] == pytest.approx(x_ref, rel=1e-3), sp


def test_rocket_performance_reproduces_manual(engine):
    perf = engine.performance(PC, MR)
    assert perf.c_star == pytest.approx(C_STAR, abs=5e-2)
    assert perf.T_c == pytest.approx(CHAMBER["T"], abs=5e-3)
    assert perf.T_throat == pytest.approx(T_THROAT, abs=5e-3)


def test_frozen_shifting_bracket(engine):
    shifting, frozen = engine.frozen_shifting_gap(PC, MR, supar=25.0)
    assert shifting == pytest.approx(ISP_VAC_25, abs=0.5)
    # Strict bracket: frozen < shifting, always.
    assert frozen < shifting
    gap = (shifting - frozen) / shifting
    # Frozen-from-chamber at Ae/At = 25 measures a few percent for LOX/LH2;
    # sanity-band the measurement (a collapse to ~0 or a blowup would mean
    # the frozen mode is not doing what RP-1311 defines).
    assert 0.01 < gap < 0.08, gap
    # The JANNAF kinetic-efficiency knockdown (~0.8–1% of shifting Isp,
    # META-3 `jannaf-eff`) must sit inside the bracket the two tables span.
    assert gap > 0.01, "bracket must contain the JANNAF kinetic band"
