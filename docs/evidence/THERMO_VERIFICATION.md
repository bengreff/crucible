# Thermochemistry verification against NASA CEA (2026-09-30)

What: `adapters/thermo.{hpp,cpp}` (namespace `crucible::thermo`) wraps Cantera 3.2.0, built from
source as a static library in `~/src/cantera`. Each `Mixture` owns its own Cantera Solution,
so there is one context per thread. Cantera types never cross the header.

`idealRocket()` solves the same problem as the CEA rocket case:
- **Chamber:** adiabatic equilibrium (HP) at Pc, with an infinite-area chamber.
- **Expansion:** isentropic at s0. In equilibrium mode the composition is re-equilibrated (SP); in frozen mode it is fixed at the chamber composition, which matches CEA frozen with freezing at the chamber.
- **Throat:** the maximum of rho*u along the isentrope, with c* = Pc / (rho*u)_throat.
- **Exit:** the supersonic point where (rho*u)*eps = (rho*u)_throat.
- **Isp_vac:** (u_e + p_e/(rho_e u_e))/g0.

Reference: NASA CEA (Glenn Fortran) through RocketCEA 1.2.3. The cases are in `tools/cea_reference.py`
and the output in `tools/cea_reference.csv`. Reactants at 298.15 K, except LOX/LH2, which uses the
CEA thermo.lib liquid enthalpies H2(L) -9012 J/mol at 20.27 K and O2(L) -12979 J/mol at 90.17 K;
the adapter is given the same values. Mechanisms: h2o2.yaml (H2/O2), gri30.yaml (CH4/O2).

Tolerance, stated before the run: Tc, c* and Isp_vac each within 0.5% of CEA, in both modes.

Result (**measured**, `build/crucible_thermo_tests`, ctest `thermo_cea_verification`): 18 rows, 0 failures,
worst |delta| 0.131%.

| Case | Tc delta | c* delta | Isp_vac delta (eq / frozen) |
|---|---|---|---|
| H2/O2 20 bar eps 40, O/F 4/6/8 | +0.07 to +0.13% | +0.01 to +0.04% | +0.002 to +0.031% / +0.03 to +0.06% |
| CH4/O2 20 bar eps 40, O/F 2.6/3.4/4.0 | +0.08 to +0.12% | +0.03 to +0.04% | +0.01 to +0.03% / +0.04 to +0.06% |
| LOX/LH2 32.75 bar eps 61, O/F 4.5/5.5/6.5 | +0.08 to +0.13% | +0.01 to +0.04% | +0.002 to +0.018% / +0.04 to +0.06% |

Reading:
- Tc is high by a systematic 3-4 K in every case.
- **Inferred cause:** the two codes use different species thermo data and species sets. Cantera's mechanism files use 7-coefficient NASA fits (h2o2 is valid up to 3500 K); CEA uses the McBride 2002 9-coefficient data and a larger product set.
- This was not chased further, because it is 4x inside tolerance and does not grow with O/F.

Cost (**measured**): about 2.3 s per point, single thread. Nearly all of it is about 100 SP
equilibrium solves in the throat and exit searches. This is fine for verification and for
boundary/reference evaluations, but it is not a per-cell path. The reacting solver will need
finite-rate kinetics and tabulated or direct property calls, not equilibrium calls (item 4).

Build:
- `cmake -S . -B build -DCRUCIBLE_WITH_CANTERA=ON`, which is the default.
- `CANTERA_ROOT` points at the scons-built tree.
- The adapter links `libcantera.a` and Accelerate.

Regenerating the reference:
`~/venvs/crucible/bin/python tools/cea_reference.py`
