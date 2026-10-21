# Table A for any propellant pair: design and criteria

Status: design and criteria, 13 October 2026, before any code for them. Nothing here is verified yet. The off-line measurement below is a Cantera flame calculation, not a table result. Results will be appended below a line; they do not change the criteria.

## Why the H/O table does not carry over

Table A (`TABLE_A.md`) reads the equilibrium at (f_H, T, ln rho_r). Its single composition axis is complete only because H and O are the only reacting elements: once the reacting mass is factored out, one degree of freedom is left. With C, H and O there are two, and with N reacting as well (NTO/MMH) three.

Two propellant streams alone fill a line in that space: z(f) = z_O + f (z_F − z_O). Here z is the reacting element mass fractions over the reacting mass, z_F and z_O the streams' and f the fuel stream's share. The engine's mixture-averaged diffusion moves cells off that line.

**Measured, 13 October 2026.** `tablepairs/flame_offline.py` (log beside it) ran Cantera counterflow diffusion flames on `gri30_highT.yaml`, mixture-averaged, with streams at 300 K and 3 MPa. The off-line distance is δ = z − z(f), with f from Bilger's coupling function (below).

| Streams | mdot per stream | T max | Largest abs. component of δ | Where |
|---|---|---|---|---|
| CH4 / O2 | 0.5 kg/m^2/s | 3,445 K | 0.056 (C +0.056, H −0.023, O −0.033) | f 0.359, 1,985 K |
| CH4 / O2 | 5.0 kg/m^2/s | 3,439 K | 0.051 | f 0.355, 2,115 K |
| H2 / O2 | both | 3,565 K | 2.5e-16 (two elements: on the line by construction) | |

The largest offset is on the rich side: 21% more carbon, relative, than the line holds at that f. The likely cause (not separated in the run) is that the H2 and H made there diffuse away faster than CO and CO2. Raising the strain tenfold barely changes it. Turbulent diffusion is the same for every species, so a turbulent chamber should show less. Laminar layers at walls, at the face and in the flames themselves keep it. A table must therefore give the equilibrium at elements up to about 0.06 off the line. A first-order correction in δ would be stretched there, and an extra table axis spanning that range costs too much memory (below).

## Design

1. **Builder inputs.**
   - Any Cantera YAML mechanism with NASA-7 thermo; the core's `Medium` reads NASA-7.
   - The fuel and the oxidizer as element formulas: CEA's RP-1 is C 1, H 1.9423; MMH is C 1, H 6, N 2; NTO is N 2, O 4.
   - The axes.
2. **Mechanism for the C/H/O/N pairs: `gri30_highT.yaml`,** as shipped with Cantera 3.2.0. It holds GRI-Mech 3.0's 53 species and transport. Joseph Shepherd (Caltech) added high-temperature NASA-7 polynomials from McBride et al.'s NASA TM-4513 data (the data behind CEA), so every species spans 200–5,000 K.
   - Its thermo range covers the chamber temperatures. CEA (RocketCEA 1.2.3) gives Tc 3,607 K for LOX/CH4 (O/F 3.5, 10.0 MPa), 3,677 K for LOX/RP-1 (O/F 2.6, 6.9 MPa) and 3,059 K for NTO/MMH (O/F 1.65, 1.0 MPa). `h2o2.yaml` and `gri30.yaml` stop at 3,500 K and 3,000 K.
   - The table needs only element formulas for the propellants. A fuel species absent from the mechanism (RP-1, MMH) matters only where the flow injects it, which is the liquid lane's concern.
   - Equilibrium is gas only, with no condensed carbon, as in the engine; CEA allows C(gr). This matters only for very rich mixtures (gas-generator O/F), not at chamber O/F.
3. **Reacting and inert species.**
   - A species is inert if it carries no element shared with any other species: Ar in `gri30_highT.yaml`, or He. Every other element reacts.
   - So N reacts in `gri30_highT.yaml` (NO, HCN and NH3 share it with C, H and O) even for LOX/CH4. Ambient N2 is then just an element vector off the line, which the exact solve below handles.
4. **The exact solve (the change from today's lookup).**
   - At fixed T, the ideal gas's constant-(T, V) equilibrium concentrations are c_k = C_k(T) exp(sum_m a_km pi_m).
     - c_k is in kmol/m^3.
     - C_k(T) = (p° / (R T)) exp(−g°_k(T) / (R T)).
     - a_km is the atoms of element m in species k.
     - pi_m is element m's potential over RT.
   - Element conservation, sum_k a_km c_k = b_m (b_m = rho Z_m / W_m), is the gradient of the strictly convex Phi(pi) = sum_k c_k − sum_m b_m pi_m. So Newton on pi with a line search on Phi converges from any start, with an SPD Hessian sum_k a_km a_kn c_k.
   - T comes from the energy: sum_k c_k W_k e_k(T) plus the inert species' energy equals rho e. This is an outer scalar Newton on T with a bracket guard; the equilibrium e(T) is increasing.
   - The result is the constant-(u, v) equilibrium at the query's own elements, energy and density, to the Newton tolerance. It is the same model as `LocalEquilibrium` and Cantera's `equilibrateUV` on the same species and thermo, with no interpolation in composition.
   - The final Y is projected onto the query's exact elements, as the builder already projects nodes (criterion 1).
5. **What the table holds.**
   - The thermo, on a fine uniform T grid (spacing an input, default 2 K) over the mechanism's range: ln C_k(T) and e_k(T) for every reacting species, computed by Cantera. The solve interpolates both linearly in T.
     - The error in ln C_k is about (dT)^2 / 8 times its second derivative, of order 2e-6 per K^2 for H2O at 3,000 K. At 2 K spacing that is about 1e-6 (derived).
     - The core's T still comes from `Medium`'s inversion with the returned Y, so the flow's energy is exact whatever this error.
   - Seeds: pi_m and the reacting energy e at nodes (f, T, ln rho_r) along the two-stream line, with f clustered at f_st and the ends as today.
     - f is Bilger's mixture fraction: (beta(z) − beta_O) / (beta_F − beta_O), with beta = sum_m v_m z_m / W_m and valences C +4, H +1, O −2, N 0. A stoichiometric cell therefore reads a stoichiometric seed off the line too.
     - For H/O it equals f_H exactly (derived).
     - The lookup finds the seed's T from the seed energies, as today's bracket search does. It then interpolates pi and starts the solve there.
     - Seeds hold ne_r + 1 values per node, against today's ns + 1. On today's grid (131 x 129 x 97, 1,639,203 nodes), a C/H/O/N seed table holds 5 values per node: about 66 MB. Today's form with 53 species holds 54 values per node. One off-line axis of 17 nodes would make that about 12 GB for C/H/O, and the two axes C/H/O/N needs about 200 GB (derived).
6. **Absent elements.** An element whose b_m is zero (C in pure oxidizer) is removed from the solve with every species carrying it.
7. **Fallback.**
   - The direct Cantera call, counted, when the solve fails to converge in a set number of iterations, or when T leaves the thermo range. Below the range the state is clamped, as today.
8. **File format 4.**
   - It adds the element names, a_km, W_k, the streams, the Bilger valences, the thermo grid and the seeds.
   - Format 3 files (today's H/O tables) load and run through today's interpolation unchanged.
   - A format-4 H/O table uses the exact solve.

## Criteria

0. **H/O unchanged.** C1 `eqtt` 64 x 12 to 1 ms on 4 threads with today's format 3 table prints the same history to the last digit as before. ctest results are unchanged.
1. **Conservation.** Every lookup in criterion 2 returns Y with the query's element mass fractions to 1e-12 and |sum Y − 1| ≤ 1e-12, on and off the line.
2. **Accuracy against direct Cantera.** The reference is `equilibrateUV` at the query's own elements, e and rho, on the same mechanism. |dT| is `Medium`'s T from the returned Y against Cantera's.
   - (a) **LOX/CH4** on `gri30_highT.yaml`, 1e5 random states:
     - f uniform on [0, 1];
     - an off-line δ uniform in the C/H/O off-line plane, each component up to 0.08 (1.4 times the measured 0.056), redrawn while any element fraction would be negative;
     - an N2 diluent fraction uniform on [0, 0.5), plus 1% of states at [0.5, 1);
     - eta and ln rho_r over the table's range, as in today's check.
     - Pass: |dT| at most 3 K for the 99th percentile and for the maximum, among states with T_eq in 300–4,800 K. Fallbacks at most 1e-4 of the states, each listed.
     - Report: the median, |dY_k| maximum per species, the Newton iteration counts (median, 99th, maximum) and the thermo-grid interpolation's share of |dT| (a run at 0.5 K spacing).
   - (b) **NTO/MMH and LOX/RP-1** on `gri30_highT.yaml`: the same states and pass marks. For NTO/MMH, δ also spans N.
   - (c) **H2/O2 on `h2o2.yaml` as a format 4 table:** the same pass marks, with today's format 3 table's |dT| on the same states printed beside.
3. **Against NASA CEA (reported, not judged).** For each pair at the operating point in Design item 2, the solve's T and molar mass at CEA's chamber density and internal energy, against CEA's Tc and M. The differences come from the species sets and the thermo fits (`gri30_highT.yaml`'s 53 species against CEA's thermo library). They are recorded, not tuned.
4. **Cost.** Mean time per lookup on one thread on the Mac, over criterion 2(a)'s states, seeded only from the table (no warm start from a previous step).
   - Pass: at most 5 us.
   - Reasons: today's H/O lookup costs 0.23 us per cell with 10 species, and a 53-species flow step should cost several times the H/O step's 1.37 us (not measured). Chemistry at 5 us stays at or below the flow's share. A per-cell warm start in the flow is a later saving.
   - Report: on- and off-line separately.
5. **Builder conservation.** At every seed node, the equilibrium's elements equal the node's to 1e-12 after projection.

The resolutions (thermo grid spacing; seed nodes in f, T and ln rho_r) are settings that criteria 2 and 4 select. They are reported with the results. Table files stay out of git and are reproducible; their provenance is printed by the builder.

## Order of work

- **A.** Builder: general streams and mechanism, reacting and inert elements, Bilger f with clustering, ne_r-element node projection, the thermo grid, format 4 with format 3 still loading. Criterion 5.
- **B.** The exact solve in `EquilibriumTable::lookup` for format 4; the check mode with off-line sampling. Criteria 1, 2 and 4 on LOX/CH4. Commit.
- **C.** NTO/MMH, LOX/RP-1 and H2/O2 format 4: criterion 2(b), 2(c) and 3. Criterion 0. Commit.
- **Later.** A per-cell warm start, and the first C/H/O chamber run (a pre-vaporized LOX/CH4 development chamber). Table B for C/H/O needs a kinetic mechanism validated at chamber pressures. That is a separate choice, recorded in RESEARCH.md when made.
