"""OFFL-3 §3.1/§3.3 — the equilibrium engine: NASA CEA behind SI boundaries.

One engine, two products (OFFL-3 §2): **local-state equilibrium** at
(p, h, Z) — the HP problem at the local elemental composition — and the
**performance reference** (c*_ideal, T_c, γ, M̄) vs (p_c, MR) via the IAC
rocket problem. Cantera is the independent §6-2 cross-check
(``crosscheck.py``), never the engine.

Units at every public boundary are SI (META-2): Pa, J/kg, K, kg/m³, m/s.
CEA internally wants bar and h/R; the conversions live here and nowhere
else (empirically pinned: ``solution.enthalpy`` [kJ/kg] = h0 · R / 1000
round-trips exactly through an HP solve).

Mixture fraction: Z = fuel-stream mass fraction, Z = 1/(1 + MR),
MR = ṁ_ox/ṁ_fuel. The advected-element coordinate of SOLV-1's shifting
mode (S19) — for a two-stream propellant the elemental composition is
fully determined by Z.

Deferred to later waves (loud, with owners — OFFL-3 §2 contract rows not
yet emitted): the frozen-path surface vs (p, h, {X_k}) (its axis set is the
SOLV-1 frozen-advection consumer's, arrives with that wave); quasi-1-D
expansion oracles (SOLV-7/VAL-2 C_F cross-check wave — station 5); B′
ablation tables (SOLV-8 wave); the Cantera transport feed (OFFL-5 spine
wave, S23). The frozen↔shifting *bracket* itself is validated now
(§6-3, ``frozen_shifting_gap``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import cea

# CEA's molar gas constant, J/(kmol·K). The h/R nondimensionalization and
# the a = sqrt(γ_s·R·T/M) sound speed both use CEA's own constant so we
# reproduce the library's algebra exactly (not CODATA-2022; the difference
# is 1.7e-6 relative and belongs to CEA's pinned data, META-3 `nasa-cea`).
R_CEA = float(cea.R)

BAR = 1.0e5  # Pa


def z_to_mr(z: float) -> float:
    """Mixture fraction (fuel-stream mass fraction) → mixture ratio o/f."""
    if not 0.0 < z < 1.0:
        raise ValueError(f"mixture fraction Z={z} outside (0, 1)")
    return (1.0 - z) / z


def mr_to_z(mr: float) -> float:
    if mr <= 0.0:
        raise ValueError(f"mixture ratio MR={mr} must be positive")
    return 1.0 / (1.0 + mr)


class ConvergenceError(RuntimeError):
    """A CEA solve did not converge — refuse, never guess (META-1 P6)."""


@dataclass(frozen=True)
class Propellant:
    """A two-stream propellant at fixed injection temperatures."""

    name: str
    fuel_species: tuple[str, ...]
    oxidizer_species: tuple[str, ...]
    injection_temps: tuple[float, ...]  # K, per species (fuel first)


#: RL10-class LOX/LH2 at the CEA standard liquid injection states
#: (RP-1311 example 8: H2(L) 20.27 K, O2(L) 90.17 K).
LOX_LH2 = Propellant(
    name="lox_lh2",
    fuel_species=("H2(L)",),
    oxidizer_species=("O2(L)",),
    injection_temps=(20.27, 90.17),
)


@dataclass(frozen=True)
class EqState:
    """One equilibrium state (SI). gamma is γ_s, the shifting isentropic
    exponent; a = sqrt(γ_s·R·T/M) — CEA's own sound-speed algebra."""

    p: float  # Pa
    T: float  # K
    rho: float  # kg/m^3
    h: float  # J/kg
    s: float  # J/(kg·K)
    gamma: float  # 1
    a: float  # m/s
    mbar: float  # kg/kmol
    cp_eq: float  # J/(kg·K)
    mole_fractions: dict[str, float]


@dataclass(frozen=True)
class RocketPerformance:
    """IAC rocket solution, chamber + throat (SI). c* is the ideal
    characteristic velocity — SOLV-7's anchor quantity."""

    p_c: float  # Pa
    mr: float
    c_star: float  # m/s
    T_c: float  # K
    gamma_c: float  # 1
    mbar_c: float  # kg/kmol
    T_throat: float  # K
    isp_vac_throat: float  # m/s (N·s/kg)


class EquilibriumEngine:
    """CEA solvers for one propellant, reused across many states."""

    def __init__(self, propellant: Propellant = LOX_LH2):
        self.propellant = propellant
        species = list(propellant.fuel_species) + list(propellant.oxidizer_species)
        # Fail loud at the boundary where the mistake is made (META-1 P6):
        # incoherent propellants otherwise surface as NaN weights deep in
        # the first solve (review finding).
        if len(propellant.injection_temps) != len(species):
            raise ValueError(
                f"propellant {propellant.name!r}: {len(propellant.injection_temps)} "
                f"injection temperatures for {len(species)} species"
            )
        overlap = set(propellant.fuel_species) & set(propellant.oxidizer_species)
        if overlap:
            raise ValueError(
                f"propellant {propellant.name!r}: species {sorted(overlap)} appear in both "
                "streams — stream weights would be degenerate (NaN); split streams must be "
                "disjoint"
            )
        if not propellant.fuel_species or not propellant.oxidizer_species:
            raise ValueError(
                f"propellant {propellant.name!r}: both streams must be non-empty"
            )
        self._reac = cea.Mixture(species)
        self._prod = cea.Mixture(species, products_from_reactants=True)
        self._fuel_w = np.array(
            [1.0 if s in propellant.fuel_species else 0.0 for s in species]
        )
        self._ox_w = 1.0 - self._fuel_w
        self._temps = np.array(propellant.injection_temps)
        self._eq_solver = cea.EqSolver(self._prod, reactants=self._reac)
        self._rocket_solver = cea.RocketSolver(self._prod, reactants=self._reac)

    def _weights(self, z: float) -> np.ndarray:
        return self._reac.of_ratio_to_weights(self._ox_w, self._fuel_w, z_to_mr(z))

    def injection_enthalpy(self, z: float) -> float:
        """h of the unreacted streams at the injection temperatures, J/kg.
        Conserved by the HP chamber solve — the chamber's h coordinate."""
        w = self._weights(z)
        return float(self._reac.calc_property(cea.ENTHALPY, w, self._temps))

    def state_php(self, p: float, h: float, z: float) -> EqState:
        """Equilibrium state at local (p [Pa], h [J/kg], Z) — the S22
        runtime-surface coordinate, one HP solve."""
        sol = cea.EqSolution(self._eq_solver)
        self._eq_solver.solve(sol, cea.HP, h / R_CEA, p / BAR, self._weights(z))
        if not sol.converged:
            raise ConvergenceError(
                f"CEA HP solve failed at p={p} Pa, h={h} J/kg, Z={z}"
            )
        T = float(sol.T)
        mbar = float(sol.M)
        gamma = float(sol.gamma_s)
        return EqState(
            p=p,
            T=T,
            rho=float(sol.density),
            h=float(sol.enthalpy) * 1.0e3,
            s=float(sol.entropy) * 1.0e3,
            gamma=gamma,
            a=float(np.sqrt(gamma * R_CEA / mbar * T)),
            mbar=mbar,
            cp_eq=float(sol.cp_eq) * 1.0e3,
            mole_fractions={k: float(v) for k, v in sol.mole_fractions.items()},
        )

    def chamber(self, p_c: float, mr: float) -> EqState:
        """Chamber stagnation state at (p_c [Pa], MR): HP at the injection
        enthalpy — the design line the §6-6 coordinate check walks."""
        z = mr_to_z(mr)
        return self.state_php(p_c, self.injection_enthalpy(z), z)

    def _rocket(
        self,
        p_c: float,
        mr: float,
        n_frz: int | None,
        supar: list[float] | None = None,
    ) -> cea.RocketSolution:
        z = mr_to_z(mr)
        w = self._weights(z)
        hc = self.injection_enthalpy(z) / R_CEA
        sol = cea.RocketSolution(self._rocket_solver)
        self._rocket_solver.solve(
            sol, w, p_c / BAR, supar=supar, hc=hc, iac=True, n_frz=n_frz
        )
        if not sol.converged:
            raise ConvergenceError(
                f"CEA rocket solve failed at p_c={p_c} Pa, MR={mr}, n_frz={n_frz}"
            )
        return sol

    def performance(self, p_c: float, mr: float) -> RocketPerformance:
        """The (p_c, MR) performance reference row (shifting equilibrium)."""
        sol = self._rocket(p_c, mr, n_frz=None)
        return RocketPerformance(
            p_c=p_c,
            mr=mr,
            c_star=float(sol.c_star[1]),
            T_c=float(sol.T[0]),
            gamma_c=float(sol.gamma_s[0]),
            mbar_c=float(sol.M[0]),
            T_throat=float(sol.T[1]),
            isp_vac_throat=float(sol.Isp_vacuum[1]),
        )

    def frozen_shifting_gap(self, p_c: float, mr: float, supar: float) -> tuple[float, float]:
        """§6-3: vacuum Isp at one area ratio, (shifting, frozen-from-chamber).
        The pair brackets delivered performance; the gap is the JANNAF
        kinetic-efficiency band (META-3 `jannaf-eff`), an epistemic UQ
        dimension — never averaged away."""
        isp = []
        for n_frz in (None, 1):
            sol = self._rocket(p_c, mr, n_frz, supar=[supar])
            isp.append(float(sol.Isp_vacuum[-1]))
        return isp[0], isp[1]
