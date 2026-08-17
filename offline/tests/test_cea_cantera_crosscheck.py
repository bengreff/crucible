"""OFFL-3 §6-2 — CEA↔Cantera independent-solver agreement on chamber
states. Two codes, no shared implementation, different thermo fits (Glenn
NASA9 vs h2o2.yaml NASA7); measured agreement is ~0.13% on T and ~0.05% on
M̄ — the gates below are ~2× the measured band, so a real regression in
either wrapper trips them while data-fit noise does not."""

import pytest

from crucible_offl.chemistry import EquilibriumEngine, mr_to_z
from crucible_offl.crosscheck import cantera_state_php

STATES = [
    (32.75e5, 5.0),  # RL10-class chamber
    (53.3172e5, 5.55157),  # RP-1311 example 8
    (10.0e5, 4.0),  # low-p fuel-rich corner of the design window
]


@pytest.fixture(scope="module")
def engine():
    return EquilibriumEngine()


@pytest.mark.parametrize("p,mr", STATES)
def test_two_solvers_agree_on_chamber_state(engine, p, mr):
    z = mr_to_z(mr)
    h = engine.injection_enthalpy(z)
    a = engine.state_php(p, h, z)
    b = cantera_state_php(p, h, z)
    assert b.T == pytest.approx(a.T, rel=3e-3)
    assert b.mbar == pytest.approx(a.mbar, rel=1.5e-3)
    assert b.cp_eq == pytest.approx(a.cp_eq, rel=1e-2)
    for sp in ("H2O", "H2", "OH"):
        assert b.mole_fractions[sp] == pytest.approx(
            a.mole_fractions[sp], abs=5e-3
        ), sp
