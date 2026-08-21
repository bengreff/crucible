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

This module is also the thermochemical half of the OFFL-3 §2 **transport
feed** (v0.5, plan S4): every state it returns carries the caloric
companions (``cp_eq``, ``cp_fr``, ``cv_eq``) that ``transport.py``
assembles with Cantera's molecular μ, k into the spine's chemical-regime
surface (OFFL-5 §3.1a). One engine states the thermochemistry; Cantera
states only the collision physics.

The **unburnt-reactant surface** (the burn-progress c = 0 branch, SOLV-4
§3.6) ships here as of plan S5: ``FrozenReactantEngine`` is the gas-phase
ideal-gas frozen reactant mixture, on the SAME CEA enthalpy reference as the
equilibrium surface (``surfaces.build_unburnt_surface`` writes it). The
frozen↔shifting **bracket** is likewise promoted from a validation check to
a shipped declared band (``frozen_shifting_band``, OFFL-3 §3.2).

Deferred to later waves (loud, with owners — OFFL-3 §2 contract rows not
yet emitted): the frozen-**advection** surface vs (p, h, {X_k}) (its axis
set is the SOLV-1 frozen-advection consumer's, rides the species-vector
state — plan **S5b**); quasi-1-D expansion oracles (SOLV-7/VAL-2 C_F
cross-check wave — station 5); B′ ablation tables (SOLV-8 wave).
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
    #: The caloric companions of the OFFL-3 §2 transport feed (v0.5), from
    #: the SAME solve: frozen-composition c_p (which pairs with a frozen
    #: conductivity to form a true molecular Prandtl number) and the
    #: EQUILIBRIUM c_v (the class-D temperature solve's linearization slope
    #: ∂e/∂T|_ρ on the very surface the runtime interpolates — FND-7 §3.3).
    #: They differ by up to ~15× where dissociation is strong, so which one
    #: a consumer wants is a physical question, never a rounding one.
    cp_fr: float  # J/(kg·K)
    cv_eq: float  # J/(kg·K)
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
    """CEA solvers for one propellant, reused across many states.

    `gas_only=True` restricts the product set to gas-phase species —
    the METASTABLE (supersaturated) equilibrium: the declared model for
    rapidly-expanding plume states where condensation kinetics are slow
    against the flow time, and the only convergent branch in the deep-cold
    corners where CEA's condensed-species iteration fails (session 12).
    The mode is stamped into the table provenance deck.
    """

    def __init__(self, propellant: Propellant = LOX_LH2, gas_only: bool = False):
        self.propellant = propellant
        self.gas_only = gas_only
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
        if gas_only:
            # Condensed CEA species carry a TRAILING parenthesized phase
            # tag — (L), (cr), (a), (b), (s), (I…) — filter exactly those.
            # A bare `"(" in name` test would also drop legitimate gas
            # species with interior parentheses (e.g. NASA-Glenn organics
            # like HO(CO)2OH) for carbon-bearing propellants — refuse to
            # guess: anything filtered must match the condensed-suffix
            # form, else this propellant needs an explicit species list
            # (session-12 review).
            import re

            cond = re.compile(r"\((L|cr|s|a|b|I{1,3}|IV|V|VI)['\d]*\)$")
            names = list(self._prod.species_names)
            gas_names = [n for n in names if not cond.search(n)]
            dropped = [n for n in names if cond.search(n)]
            odd = [n for n in gas_names if "(" in n]
            if odd:
                raise ValueError(
                    f"gas_only filter: species {odd} carry parentheses but no "
                    "recognized condensed-phase suffix — extend the filter or "
                    "supply an explicit gas species list; refusing to guess"
                )
            if not dropped:
                raise ValueError(
                    "gas_only requested but no condensed species were present "
                    "to drop — the flag would be a silent no-op; remove it or "
                    "check the product set"
                )
            self._prod = cea.Mixture(gas_names)
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

    def _state_from(self, sol: "cea.EqSolution", p: float) -> EqState:
        """Shared SI unpacking of a converged CEA equilibrium solution —
        one conversion site (the module doc's rule)."""
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
            cp_fr=float(sol.cp_fr) * 1.0e3,
            cv_eq=float(sol.cv_eq) * 1.0e3,
            mole_fractions={k: float(v) for k, v in sol.mole_fractions.items()},
        )

    def state_php(self, p: float, h: float, z: float) -> EqState:
        """Equilibrium state at local (p [Pa], h [J/kg], Z) — the S22
        runtime-surface coordinate, one HP solve."""
        sol = cea.EqSolution(self._eq_solver)
        self._eq_solver.solve(sol, cea.HP, h / R_CEA, p / BAR, self._weights(z))
        if not sol.converged:
            raise ConvergenceError(
                f"CEA HP solve failed at p={p} Pa, h={h} J/kg, Z={z}"
            )
        return self._state_from(sol, p)

    def state_tp(self, p: float, t: float, z: float) -> EqState:
        """Equilibrium state at (p [Pa], T [K], Z) — the TP problem. Used
        by the OFFL-5 §3.1a transport surface for `∂h/∂Z|_{p,T}`, which is
        a derivative **at fixed temperature** and so cannot be taken on the
        (p, h, Z) coordinate without unwinding the enthalpy change first."""
        sol = cea.EqSolution(self._eq_solver)
        self._eq_solver.solve(sol, cea.TP, t, p / BAR, self._weights(z))
        if not sol.converged:
            raise ConvergenceError(
                f"CEA TP solve failed at p={p} Pa, T={t} K, Z={z}"
            )
        return self._state_from(sol, p)

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

    def frozen_shifting_band(
        self, p_c: float, mr: float, supar: float
    ) -> "KineticEfficiencyBand":
        """OFFL-3 §3.2 (plan S5) — the frozen↔shifting model-form band as a
        **shipped declared datum**, not just a validation check. The pair
        `frozen_shifting_gap` returns brackets delivered performance; this
        wraps it as the labeled `jannaf-eff` band a run records in its
        pedigree so the report is an interval, not a point. The band is the
        relative gap `(Isp_shift − Isp_frozen)/Isp_shift` (frozen always the
        lower bound — kinetics never *help*)."""
        isp_shift, isp_frozen = self.frozen_shifting_gap(p_c, mr, supar)
        if not (isp_shift > 0.0 and isp_frozen > 0.0):
            raise ConvergenceError(
                f"frozen↔shifting band: non-physical Isp pair "
                f"({isp_shift}, {isp_frozen}) at p_c={p_c}, MR={mr}, ε-ratio={supar}"
            )
        if isp_frozen > isp_shift:
            # Shifting is the equilibrium (upper) limit by construction;
            # frozen-from-chamber cannot exceed it. If it does, the rocket
            # solve is not the pair we think it is — refuse (META-1 P6).
            raise ConvergenceError(
                f"frozen↔shifting band: frozen Isp {isp_frozen} exceeds shifting "
                f"{isp_shift} — the bracket is inverted; the solve is wrong"
            )
        return KineticEfficiencyBand(
            p_c=p_c,
            mr=mr,
            area_ratio=supar,
            isp_shifting=isp_shift,
            isp_frozen=isp_frozen,
            relative_gap=(isp_shift - isp_frozen) / isp_shift,
        )


@dataclass(frozen=True)
class KineticEfficiencyBand:
    """OFFL-3 §3.2 — the shipped frozen↔shifting model-form band (plan S5).
    A run records this so its performance is reported as the `[frozen,
    shifting]` interval with the `jannaf-eff` anchor, never a hidden point
    choice.

    `relative_gap` is the **full raw bracket width** `(shift−frozen)/shift` —
    a few percent for LOX/LH₂, and it widens with the area ratio (more
    recombination energy the frozen limb leaves on the table): ~3.7–4.8% at
    ε = 61 over MR 5.0–5.5 (≈ 4.3% at the RP-1311 example-8 anchor). The
    **JANNAF kinetic-efficiency knockdown** (~0.8–1% of shifting
    Isp, META-3 `jannaf-eff`) is the *data-anchored delivered estimate that
    sits inside* this bracket — H/O kinetics are fast, so the delivered value
    hugs the shifting (equilibrium) end. The bracket is the model-form
    interval; the JANNAF factor says where in it the real engine lands."""

    p_c: float  # Pa
    mr: float
    area_ratio: float  # A_exit / A_throat the band is evaluated at
    isp_shifting: float  # m/s (N·s/kg), the equilibrium upper limit
    isp_frozen: float  # m/s, the frozen-from-chamber lower limit
    relative_gap: float  # (shift − frozen)/shift ≥ 0


#: The gas-phase suffix stripped to turn a condensed reactant species into
#: its gas form for the frozen unburnt-mixture branch (H2(L) → H2). The
#: same condensed-phase-suffix grammar the `gas_only` product filter uses
#: (chemistry.py EquilibriumEngine) — one place, one regex.
_CONDENSED_SUFFIX = __import__("re").compile(r"\((L|cr|s|a|b|I{1,3}|IV|V|VI)['\d]*\)$")


def gasify(species: str) -> str:
    """Strip a condensed-phase suffix so a liquid injection species names
    its gas-phase form (H2(L) → H2, O2(L) → O2); a bare gas name is
    unchanged. The unburnt branch is the GAS reactant mixture (SOLV-1 W4's
    two-phase drift-flux owns liquid/vapor — plan S15)."""
    return _CONDENSED_SUFFIX.sub("", species)


@dataclass(frozen=True)
class FrozenReactantState:
    """One state of the gas-phase **frozen reactant mixture** (SI) — the
    burn-progress blend's c = 0 branch (SOLV-4 §3.6). `gamma` is the frozen
    isentropic exponent c_p,fr/c_v,fr; `a = √(γ_fr·R̄·T/M̄)` — CEA's own
    sound-speed algebra, identical to `EqState`'s."""

    p: float  # Pa
    T: float  # K
    rho: float  # kg/m^3 (ideal gas: p·M̄/(R̄·T))
    h: float  # J/kg (CEA reference — the SAME as EqState.h, so the blend adds)
    gamma: float  # 1 (frozen c_p/c_v)
    a: float  # m/s
    mbar: float  # kg/kmol
    cp_fr: float  # J/(kg·K)


class FrozenReactantEngine:
    """The gas-phase ideal-gas **frozen reactant mixture** thermo for one
    propellant — the OFFL-3 §3.3 unburnt-reactant surface's engine (plan S5).

    Deliberately reuses the **same CEA reactant `Mixture` machinery** as
    `EquilibriumEngine`: `calc_property(ENTHALPY, …)` is the very call
    `injection_enthalpy` uses, so the unburnt branch's enthalpy is on the
    **same reference** as the burnt equilibrium surface. The SOLV-4 §3.6
    blend `h = (1−c)·h_u + c·h_b` is a category error on two references, so
    this shared reference is a correctness requirement, not tidiness.

    The mixture is FROZEN (no equilibration) and IDEAL (ρ = p·M̄/(R̄·T)) —
    the same ideal-gas assumption CEA makes for the equilibrium products, so
    the two branches differ only in composition.

    **Validity floor (review-clarified).** The reactant elements have zero
    formation enthalpy, so — unlike the equilibrium surface's condensing
    products — CEA *converges* down to ~35 K here (which is what lets the
    shipped grid overhang the envelope for interpolation). But convergence is
    not validity: the NASA polynomials are extrapolated below their ~200 K fit
    floor, and the extrapolated γ is still physical (~1.47) at the **100 K
    envelope floor** but degrades below ~60 K (γ → 1.14 at 40 K, `c_p` blowing
    up) and `c_p ≤ 0` below ~35 K, where `state_php` refuses. The shipped
    surface's **envelope** floor is 100 K; the sub-floor grid nodes are the
    declared metastable-model overhang, gated off at runtime. The real cold
    two-phase state is SOLV-1's W4 drift-flux extension (plan S15).
    """

    #: Deterministic T-inversion bracket [K]. Wide enough to bracket every
    #: envelope the surface declares; fixed (never data-dependent) so the
    #: bisection is bit-reproducible (META-1 §2.1). The reactant thermo is
    #: a declared metastable ideal-gas model below the liquefaction line
    #: (two-phase = plan S15); it stays convergent over this whole bracket.
    T_SOLVE_LO: float = 20.0
    T_SOLVE_HI: float = 6000.0
    #: Fixed bisection count: (hi−lo)/2^60 ≈ 5e-15 K — machine-exact for a
    #: monotone smooth h(T). Not a tunable.
    N_T_BISECT: int = 60
    #: Central-difference half-step for c_p = ∂h/∂T [K]. Fixed/declared: an
    #: adaptive step would not be bit-reproducible (META-1 §2.1).
    DT_CP: float = 0.5

    def __init__(self, propellant: Propellant = LOX_LH2):
        self.propellant = propellant
        fuel = [gasify(s) for s in propellant.fuel_species]
        ox = [gasify(s) for s in propellant.oxidizer_species]
        species = fuel + ox
        dup = {s for s in species if species.count(s) > 1}
        if dup:
            raise ValueError(
                f"propellant {propellant.name!r}: gas-phase reactant species "
                f"{sorted(dup)} collide after gasify — a frozen mixture needs "
                "distinct gas species per stream; refusing to guess"
            )
        self._mix = cea.Mixture(species)
        self._fuel_w = np.array([1.0 if s in fuel else 0.0 for s in species])
        self._ox_w = 1.0 - self._fuel_w

    def _weights(self, z: float) -> np.ndarray:
        """Mass weights of the gas reactant mixture at mixture fraction Z —
        CEA's own of_ratio split (identical routine to `EquilibriumEngine`)."""
        return self._mix.of_ratio_to_weights(self._ox_w, self._fuel_w, z_to_mr(z))

    def mbar(self, z: float) -> float:
        """Mean molar mass [kg/kmol] of the gas reactant mixture at Z.
        `of_ratio_to_weights` returns UN-normalized mass weights (CEA scales
        them internally), so M̄ = Σ(mass)/Σ(moles) — normalization-
        independent — never `1/Σ(moles)`."""
        w = self._weights(z)
        moles = self._mix.weights_to_moles(w)
        return float(np.sum(w)) / float(np.sum(moles))  # g/mol = kg/kmol

    def _h_of_t(self, t: float, w: np.ndarray) -> float:
        """Frozen-mixture specific enthalpy [J/kg] at uniform temperature t.
        Pressure-independent (ideal gas) — the surface's p-dependence is in
        ρ alone."""
        return float(self._mix.calc_property(cea.ENTHALPY, w, np.full(w.shape, t)))

    def t_of_h(self, h: float, z: float) -> float:
        """Temperature [K] of the frozen mixture at (h, Z) — the monotone
        inverse of `_h_of_t`, by deterministic fixed-count bisection."""
        w = self._weights(z)
        lo, hi = self.T_SOLVE_LO, self.T_SOLVE_HI
        h_lo, h_hi = self._h_of_t(lo, w), self._h_of_t(hi, w)
        if not (h_lo <= h <= h_hi):
            raise ConvergenceError(
                f"frozen-reactant enthalpy h={h} J/kg at Z={z} is outside the "
                f"solvable T bracket [{lo}, {hi}] K (h ∈ [{h_lo:.1f}, {h_hi:.1f}])"
            )
        for _ in range(self.N_T_BISECT):
            mid = 0.5 * (lo + hi)
            if self._h_of_t(mid, w) < h:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    def state_php(self, p: float, h: float, z: float) -> FrozenReactantState:
        """Frozen gas reactant state at local (p [Pa], h [J/kg], Z)."""
        w = self._weights(z)
        t = self.t_of_h(h, z)
        # M̄ = Σ(mass)/Σ(moles); of_ratio weights are un-normalized (see `mbar`).
        mbar = float(np.sum(w)) / float(np.sum(self._mix.weights_to_moles(w)))
        rho = p * mbar / (R_CEA * t)  # ideal gas; M̄ [kg/kmol], R_CEA [J/(kmol·K)]
        # c_p = ∂h/∂T (frozen composition); central difference of the same
        # monotone enthalpy the inversion uses.
        cp_fr = (
            self._h_of_t(t + self.DT_CP, w) - self._h_of_t(t - self.DT_CP, w)
        ) / (2.0 * self.DT_CP)
        # `c_p > 0` is dh/dT > 0 — the inversion's monotonicity assumption,
        # checked AT THE ROOT (META-1 P6, review-hardened). The bisection is
        # valid only on an increasing branch; a `c_p ≤ 0` here means the
        # NASA-polynomial extrapolation folded (it does below ~40 K, or for
        # an extreme ox-rich Z outside any shipped envelope), so the located
        # root is off the physical branch — refuse rather than return it. On
        # every shipped-envelope state this holds with wide margin.
        if cp_fr <= 0.0:
            raise ConvergenceError(
                f"frozen-reactant c_p = {cp_fr} ≤ 0 at (p={p}, h={h}, Z={z}) — the "
                "enthalpy is non-monotone in T here (the extrapolated thermo folded); "
                "the bisection root is off the physical branch, refusing"
            )
        r_specific = R_CEA / mbar
        cv_fr = cp_fr - r_specific
        if cv_fr <= 0.0:
            raise ConvergenceError(
                f"frozen-reactant c_v = {cv_fr} ≤ 0 at (p={p}, h={h}, Z={z}) — "
                "the ideal-gas relation c_p − R/M̄ degenerated; refusing"
            )
        gamma = cp_fr / cv_fr
        a = float(np.sqrt(gamma * r_specific * t))
        return FrozenReactantState(
            p=p, T=t, rho=rho, h=h, gamma=gamma, a=a, mbar=mbar, cp_fr=cp_fr
        )

    def enthalpy_at(self, t: float, z: float) -> float:
        """Frozen-mixture enthalpy [J/kg] at (T, Z) — grid/envelope planning
        (`surfaces.unburnt_reactant_grid` maps a T-range to the rectangular
        h-envelope every Z can convergently populate)."""
        return self._h_of_t(t, self._weights(z))
