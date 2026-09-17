"""OFFL-3 §6-2 — Cantera as the independent equilibrium cross-check.

Cantera never generates a shipped table (CEA is the engine, §3.1); it
exists to catch a wrong answer neither solver alone would reveal. The two
codes share no implementation and use different thermo fits (NASA Glenn
NASA9 vs Cantera's NASA7 ``h2o2.yaml``), so agreement within the measured
data-fit band (~0.15% on T) is a real two-solver check, and the residual
disagreement is part of the thermodynamic-data uncertainty story (§3.4).

The (H, P) dance: the local-state coordinate h is only reachable with
reacted composition (it carries H₂O formation enthalpy), so we equilibrate
hot at TP first, then pin (H, P) and re-equilibrate.
"""

from __future__ import annotations

from dataclasses import dataclass

import cantera as ct

CANTERA_MECH = "h2o2.yaml"


@dataclass(frozen=True)
class CanteraState:
    T: float  # K
    mbar: float  # kg/kmol
    cp_eq: float  # J/(kg·K), finite-difference shifting cp
    mole_fractions: dict[str, float]


def cantera_state_php(p: float, h: float, z: float, dh: float = 2.0e4) -> CanteraState:
    """Independent equilibrium at the same (p [Pa], h [J/kg], Z) coordinate
    the CEA surface uses. Same enthalpy datum (formation enthalpies over
    elements at 298.15 K) — the h axis crosses solvers unchanged."""
    gas = ct.Solution(CANTERA_MECH)
    gas.TPY = 3000.0, p, {"H2": z, "O2": 1.0 - z}
    gas.equilibrate("TP")
    gas.HP = h, p
    gas.equilibrate("HP")
    T = float(gas.T)
    mbar = float(gas.mean_molecular_weight)
    x = {sp: float(gas[sp].X[0]) for sp in gas.species_names}
    # Shifting cp = (∂h/∂T)_p along the equilibrium path, by FD in h.
    gas.HP = h + dh, p
    gas.equilibrate("HP")
    cp_eq = dh / (float(gas.T) - T)
    return CanteraState(T=T, mbar=mbar, cp_eq=cp_eq, mole_fractions=x)
