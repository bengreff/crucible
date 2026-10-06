# Table A: tabulated shifting equilibrium

Status: criteria stated 5 October 2026, about 22:10, before any table code. Results are appended below the line; they do not change the criteria.

## What is tested

Table A replaces the per-cell Cantera constant-(u, v) equilibrium of `Chemistry::LocalEquilibrium` with a table built offline by the same Cantera call (TECHNICAL_PLAN, *Lightweight engine*). The model is unchanged, so the reference is the direct Cantera call at the same state. The table is built for the C1 mechanism, `h2o2.yaml`.

- **Axes.**
  - f_H = Z_H / (Z_H + Z_O), with a node at the stoichiometric f of H2O. The complete-combustion products bend there.
  - Z_N, the N2 mass fraction.
  - eta = (e − e_lo) / (e_hi − e_lo), normalised per (f_H, Z_N) node. e_lo is the energy of the complete-combustion products at 200 K; e_hi is the energy of the same products, frozen, at 6000 K.
  - ln rho.
- **Values.** Each node stores the equilibrium Y_k and T, computed by Cantera's equilibrateUV started from the complete-combustion products at the node's energy.
- **Lookup.** Multilinear interpolation, with each (f_H, Z_N) corner using its own eta.
  - Z_H, Z_O and Z_N are bilinear in (f_H, Z_N), and multilinear weights reproduce bilinear functions exactly. So the returned Y carries the query's elements to rounding, and sum Y = 1.
  - The flow keeps its total energy. The cell's T comes from the core's e-to-T inversion with the new Y, so energy is exact. The tabulated T is only a consistency check.
- **Outside the table.**
  - e below e_lo: eta is clamped to 0, which returns the complete-combustion products, the low-temperature limit of equilibrium. Counted.
  - e above e_hi, rho outside the axis, or a species not made of H, O and N2 above 1e-12 (AR): the direct Cantera call is used. Counted.

## Criteria

1. **Conservation.**
   - At every node, the equilibrium's element mass fractions equal the node's to 1e-12, and |sum Y − 1| ≤ 1e-12.
   - For every lookup in criterion 2, the returned Y has the query's element mass fractions to 1e-12, and |sum Y − 1| ≤ 1e-12.
   - In the C1 runs of criterion 3, the mass and energy budgets are reported against CHAMBER_C1's 1e-11 limits, with the same normalisation caveat.
2. **Accuracy against direct Cantera** (|dT| is the core's T from the table's Y, against Cantera's equilibrium T at the same elements, e and rho).
   - (a) Random states: 1e5 states uniform in (f_H, Z_N, eta, ln rho) over the table's range. Report the median, 99th percentile and maximum of |dT|, and the maximum |dY_k| per species. Pass: 99th percentile of |dT| ≤ 3 K for states whose equilibrium T is between 200 and 4000 K.
   - (b) The run's own states: in the C1 runs of criterion 3, an audit compares every cell on every Nth reaction call with the direct Cantera call from the same pre-reaction state (N declared per run). Pass: maximum |dT| ≤ 3 K over audited cells inside the table. Report the maximum |dY_k|, and the clamp and fallback counts.
   - Why 3 K: c* goes as sqrt(T / M). A uniform 3 K error at 3350 K shifts c* by 0.045%, inside criterion 3's 0.05% (derived).
3. **The C1 cold start.** The table (`crucible_chamber_study eqt`) is compared with direct Cantera (`eq`) on 32x6 and 64x12 to 8 ms, same tree, same machine and same threads (5).
   - Pass: settled c*, vacuum Isp and vacuum thrust agree within 0.05%. 0.05% is one eighth of the 64x12 grid error in c* (+0.40%, CHAMBER_C1), so the table's error is negligible against the discretisation's.
   - Pass: injector pressure and outlet mass flow agree within 0.1% at every printed time from 0.5 ms.
   - Pass: CHAMBER_C1's settling criterion (drift below 1e-3 over the last 1 ms).
4. **Cost.** Wall time to 8 ms on the Mac for `eq` and `eqt`, same tree and threads. `eqt` is also run on 1 thread.
   - Pass: speedup of 10 or more on 64x12.
   - Pass: the 64x12 cold start to full thrust in 5 minutes or less. Full thrust means 4 ms: the recorded 64x12 vacuum thrust is within 0.1% of its settled value from 3.5 ms (`c1/eq64_log.txt`).
   - Reported: the table's own cost per cell and its share of the step.

The table's resolution (nodes per axis) is the setting that criteria 2 and 3 select; it is reported with the result. The table file is built by `crucible_build_equilibrium_table`, is not kept in git (it is reproducible), and its provenance (mechanism, axes, Cantera version, build time, node conservation error) is printed by the builder and copied into the results below.

---

## Results
