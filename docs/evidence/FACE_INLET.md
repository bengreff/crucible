# The face inlet (a) and the plate's wall temperature

Registered 22 October 2026, about 02:00 CDT, before any code. It implements the interface spec in `docs/SESSION_HANDOFF.md` (*For inquiry2: the face inlet (a) and the face's wall*, 21 October about 02:40) for `tier1_ensemble_plan.md`'s decisions 6 and 8. The spec's properties and gates are kept; what follows fixes inquiry2's choices and states the gates as they will be judged.

## 1. The open-area supply

**The switch.** `Definition::openAreaSupplies{false}`. False runs today's ring-centre code unchanged. True applies this section to every supply of a Chamber.

**Geometry, per level.** Face ring j spans node radii [ρ_j, ρ_{j+1}] = R_0 [f_j, f_{j+1}] (R_0 the face radius, f the mesh's ring fractions).
- Supply s's share of ring j: a_sj = π (min(r_o, ρ_{j+1})² − max(r_i, ρ_j)²), or 0 if not positive. A ring the supply covers whole takes a_sj = `axialArea(0, j)` exactly, so an aligned layout keeps today's areas bit for bit.
- φ_sj = a_sj / `axialArea(0, j)`. The supply's area is `supplyArea_[s]` = Σ_j a_sj, the true area to rounding, so the mass through the face sums to ṁ_s opening(t) to rounding at any mesh.
- The wall share is w_j = 1 − Σ_s φ_sj, set to 0 when below 1e-12 (rounding of a covered ring).

**The ring's kind (inquiry2's choice, as the spec proposed).** Ring j is *plate* if w_j ≥ 1/2, else it takes the supply with the largest φ_sj (the lower index on a tie). The kind is `faceSupply_[j]`, as today, and it alone governs:
- the gradient stencil's face neighbour (a wall value at the face midpoint, or zero normal gradient);
- the wall distance's plate;
- the wall functions' prescribed first cell, and the resolved SST's wall ω cells.

**φ and w govern the fluxes.**
- Inviscid: ring j's face flux is Σ_s φ_sj `supplyFace(s, g_s, ring j's interior)` + w_j × the plate's pressure-only flux. Species: Σ_s φ_sj g_s Y_s; ρk and ρω: Σ_s φ_sj × each supply's turbulence flux. g_s = ṁ_s opening(t) / `supplyArea_[s]`. Choking and the characteristic stay per supply against ring j's interior.
- A zero-flow supply (ṁ_s = 0) is a closed port: `supplyFace` at g = 0 returns pressure only. It counts in Σ φ, so its area is not wall: no shear and no heat (the LOX slots).
- Viscous, the wall share:
  - a *plate* ring passes w_j times its whole wall flux (resolved or the law's τ_w and q_w; k and ω fluxes and the corrections' fluxes too);
  - a *supply* ring with w_j > 0 passes w_j times the wall's momentum and heat flux only. Under the wall functions it is a wall face solved from the cell's state at the face distance, its cell not prescribed. Resolved, it is the wall gradient with μ_t = 0 at the wall. Its k and ω are solved, under the supply share's zero-gradient face, and no k or ω flux to the partial wall is passed.
  - *Declared limit:* on a supply ring the partial wall's sink of k and ω (resolved) is not represented. It lies under half the face (inferred small; tier 1 runs the wall functions, where it does not arise).
- Rates: the plate share's pressure force goes to `wallAxial`, the supplies' share (closed ports included) to the inlet momentum. Every viscous end flux at i = 0 is the wall's, so it goes to `wallAxial` and `wallHeat`.

**Refusals under the switch.** "Holds no ring centre" is dropped. Refused: two supplies' true annuli overlapping; a supply's outer radius beyond R_0 by more than 1e-12 relative.

**AMR.** Each level builds its own a, φ, w and kinds from its own face rings in its constructor; `core/amr.cpp`'s copy of the base parent's supply is skipped under the switch. A coarse face area is the sum of its children's, so a coarse a_sj is too.

## 2. The plate's wall temperature

`Definition::plateTemperature`, a `std::optional<double>`. The plate's temperature T_p is its value, or `wallTemperature` when empty; 0 is adiabatic. Plate faces read T_p where they now read `wallTemperature`:
- the wall value of T in the resolved wall gradient, and the adiabatic wall's zero heat;
- the resolved SST's wall ν, ρ_w and μ_w for the plate's cells (the corner cell (0, nr − 1) stays with the side wall, as it is listed today);
- the wall law's T_w.

**The adiabatic law** (member 19). With T_w = 0 on a wall face under the wall functions: β = 0, q_w = 0, and T_w = T_1 + r u_1² / (2 c_p) inside the solve (WALL_FUNCTIONS.md, eq. 10.8 with β = 0; r = Pr_w^(1/3) at T_w, c_p the layer mean). The refusal at `core/flow.cpp:107` then admits an adiabatic wall.

## 3. Gates

- **F1, switch off:** the fast ctest gives 7 pass and 4 fail with the sorted FAIL lines identical to `docs/evidence/amr_c/ctest_mass_error_2026-10-21.txt`. rl10dev's short run (`fr`, 76 × 8, to 30 µs) is `cmp`-identical to the parent commit's (mesh, history, the three fields, the log apart from wall times).
- **F2, switch on, aligned** (every φ 0 or 1): a chamber with supplies on ring nodes, laminar and SST with the wall functions, 200 steps. The general path is run (no branch for φ = 1). *Expected:* not `cmp`-identical, since the species flux is formed as area × (g Y) in place of (area × g) Y; the state's largest relative difference under 1e-10, and the face mass and momentum rates equal to 1e-12 relative.
- **F3, switch on, tier 1's true layout** (with the LOX slots as closed ports) at 149 × 16, at 1192 × 128 and on AMR (149 × 16 with 3 levels): Flow accepts it; each supply's face mass flow equals ṁ_s opening(t) to 1e-12; and on every ring it touches, its face ρu times its true area equals ṁ_s opening(t) to 1e-12, so its velocity and J are the element's own.
- **F4, AMR:** over each coarse face ring, the sum of its children's area × flux equals its own to 1e-12 relative, for mass, momentum and energy, at a uniform state (all levels' interiors equal). In a developed state the mass and energy flux still do (they do not depend on the interior); the momentum flux, which carries the interior pressure, is reported.
- **F5, the plate temperature:** (a) empty and set equal to `wallTemperature` give `cmp`-identical runs; (b) at 0 with the wall functions, the plate's q_w is 0 to rounding on every plate face, and its T_w lies at the recovery temperature (between T_1 and T_1 + u_1²/(2 c_p)).

## 4. Results

Measured on backhouse (Release, 4 threads at most) at 566b12c, against the registration commit d9e2684 built beside it for F1. Raw outputs in `face_inlet/`: `f245.txt` (`crucible_face_inlet_tests`, ctest `face_inlet`), `f3.txt` (`CRUCIBLE_T1_OPEN_AREA=1 crucible_rl10_tier1 149x16 1192x128`), `f1_rl10dev.txt` (the `cmp` and log diff, then the exit code of the diff of tier 1's default dry run against d9e2684's build), `ctest_2026-10-22_backhouse.txt` (the fast ctest and its FAIL lines against d9e2684's and the Mac record's).

- **F1: passes against the parent on the same machine. As written, against the Mac record, it fails, and the failure is the machine's.**
  - rl10dev's short run: mesh, history and the three fields are `cmp`-identical. The log differs only in its two wall-time lines (46.1 against 48.3 µs per cell update). Tier 1's dry run without the switch prints the same text as d9e2684's (diff exit 0).
  - The fast ctest (`face_inlet/ctest_2026-10-22_backhouse.txt`) gives 8 pass and 4 fail of 12. The 4 failures are the record's (amr, amr_subcycling, transport, spray), and the 12th test, `face_inlet`, passes. The sorted FAIL lines (16) are identical, values included, to d9e2684's built and run beside it.
  - Against `amr_c/ctest_mass_error_2026-10-21.txt` (the Mac) they differ in three places, as written: a second spray failure, KH-RT liquid mass through every child at 1.044e-12 against 1e-13; and amr's level 2 and 3 closures at 3.943e-14 and 8.024e-14, against 3.938e-14 and 8.132e-14. d9e2684 gives the same three differences on backhouse, so they belong to the machine and not to this change (measured). The cause, compiler and libm, is inferred.
  - The Mac's ctest was not rerun tonight (the Director: nothing heavier than a quick unit test on the Mac while it is short of memory).
- **F2, pass.** Not `cmp`-identical, as expected. The state's largest relative difference is 4.36e-14 laminar and 3.30e-14 under SST with the wall functions (limit 1e-10); the inlet mass flow differs by 0, the inlet momentum flux by 1.6e-16 and the wall's axial force by 7.5e-16 and 3.7e-16 (limit 1e-12). The open run's budgets are 3e-16 to 4e-16.
- **F3, pass,** for every grid named: per ring, the face mass flow equals Σ a_sj ṁ_s/A_s to 5.8e-16 (149 × 16), 2.5e-14 (1192 × 128) and 7.3e-15 at most over the AMR levels; per supply, the face mass flow and ρu A_true on each ring it touches to 1.0e-15, 1.0e-15 and 3.0e-15 (38, 129 and up to 76 supply-rings); the supply areas to 8.6e-16 and 1.0e-15 of the true ones. Under the switch Flow accepts the true slots (33 supplies), which the ring-centre code refuses at both grids (G6, unchanged). One step on each grid gives a mass budget of 1e-15; its energy budgets (1.9e-7 at 149 × 16, 2.9e-8 at 1192 × 128) are reported, not gated: the run is cold and unignited and one step says nothing about the walk.
- **F4, pass.** At the uniform state the children's area × flux matches each coarse ring's to 2.9e-16 (mass), 1.7e-16 (axial momentum) and 4.6e-16 (energy). After 100 base steps (t 2.5e-5 s) mass and energy still match to 2.9e-16 and 4.6e-16; the momentum flux differs by 3.1e-2 relative, as expected: it carries each level's own interior pressure (reported).
- **F5, pass.** (a) Empty and 400 K (= `wallTemperature`) give bit-identical state, species, ρk, ρω and budgets, with the wall functions and with resolved SST. (b) With `plateTemperature` 0 and the wall functions, the largest plate |q_w| is 1.5e-12 W/m² (ring-centre, 2 plate faces) and 2.3e-12 W/m² (open-area, 3 plate faces), against the side wall's 1.8e4 W/m² (ratios 8e-17 and 1.2e-16, limit 1e-9). T_w lies between T_1 and T_1 + u_1²/(2 c_p) on every plate face (both bounds met to 0 K), at 0.892 of the dynamic temperature rise, which the law sets to Pr_w^(1/3) (Pr_w about 0.71, inferred, not printed).

## 5. The gas energy leak (23 October 2026)

Registered about 00:50 CDT, before the fix's code. The Director (00:35) relayed inquiry's finding: one plain gas step on tier 1's run Definition (149 × 16, no spray, no AMR) leaves an energy budget of −1.0e-8, and the 2-level probe reaches 2.45e-7 at 0.02 ms against R3(a)'s 1e-11, while mass and momentum close at 1e-16. Both reproduce on this tree (measured, `face_inlet/leak_2026-10-23/before/`): −1.011e-08 and 2.45e-07.

**Isolation (measured, before the fix).** F6, a new mode of `crucible_face_inlet_tests`, runs F4's layout on one level. The closed port (2.1 to 2.9 mm) leaves ring 2 a supply ring with a wall share, w = 0.470. The gas is N2 at 300 K and 1 atm, the side wall is at 400 K, and each variant prints the budgets after one step and after 200 steps (`face_inlet/leak_2026-10-23/f6_before.log`):

| variant | energy, 1 step | energy, 200 steps | momentum, 200 steps |
|---|---|---|---|
| open-area, wall functions and corrections, plate 400 K | 2.147e-08 | 4.288e-06 | −3.7e-17 |
| open-area, wall functions, plate 400 K | 2.147e-08 | 4.288e-06 | 4.9e-17 |
| open-area, wall functions, adiabatic plate | −7.5e-17 | 3.4e-16 | 3.7e-17 |
| open-area, resolved SST, plate 400 K | 1.728e-08 | 3.266e-06 | 8.038e-08 |
| open-area, laminar, plate 400 K | 1.755e-08 | 3.314e-06 | 8.002e-08 |
| open-area, inviscid (the face inlet alone) | 5.2e-17 | −5.6e-17 | −8.5e-17 |
| ring-centre on ring nodes, wall functions and corrections, plate 400 K (the walls alone) | −7.7e-18 | 2.0e-16 | 8.0e-18 |

The face inlet's convective accounting closes, and so do the walls and the wall functions on ring-centre supplies. The corrections play no part. Only open-area supplies with a plate that exchanges heat leak, under every viscous model.

**The mechanism.** `Flow::transportFluxes` forms ring j's viscous end flux at i = 0 under open-area supplies as w_j times the wall's, on every ring with w_j > 0 (section 1, *Viscous, the wall share*), and the boundary rates book it on the same rings. The cell gather (`endFace`) skips any face whose kind is Supply. So on a supply ring with 0 < w_j < 1/2, the wall share's heat and momentum are booked but never reach the cell. That is the code from 4a7178c: I changed the face flux and the rates, but not the gather.

A temporary print, reverted and not committed, gave ring 2's booked share as +0.0179 W and +0.0180 W in the first step's two stages, heat into the gas. Times dt = 2.5e-7 s that is 4.6e-9 J (derived). Over |E0|, about 0.21 J for a 2.1e-6 m³ chamber of N2 at 300 K and 1 atm (derived, e about −87 kJ/kg on the formation datum), that is 2.1e-8, the measured 2.147e-8.

E0 is negative here (inferred from the signs), so heat that is booked but not received shows as a positive error. On tier 1's run, with the plate at 199.67 K below the 300 K gas, the sign flips to the −1.0e-8 seen.

Under the wall functions the law's traction on the plate has no axial part, so the axial momentum budget still closed. Resolved, the partial wall's normal stress is lost from it too (8e-8).

**A missed signal, recorded.** Section 4's F3 printed one-step energy budgets of 1.9e-7 and 2.9e-8. I set them aside as "reported, not gated: the run is cold and unignited". That was wrong: a conservation budget closes to rounding in any state, so they were this leak, visible on 22 October. No gate had a supply ring with a wall share and a plate that exchanges heat. F2's supplies sit on ring nodes, F5(b)'s plate is adiabatic, and F4 compares face fluxes, not budgets. F6 covers it from now on.

**The fix.** `endFace` skips an end face on the same condition the boundary rates use. At i = 0 under open-area supplies it skips only a ring with no wall share; elsewhere it skips a supply face, as now. A supply ring with w_j > 0 then passes w_j times the wall's momentum and heat to its cell, as section 1 registered. Its k and ω flux stays out: it is zero on such a ring, and the gather's turbulence term keeps its Wall-kind test. Without open-area supplies, and on rings with w_j = 0 or w_j ≥ 1/2, nothing changes.

**The checks.**
- **L1:** F6 passes. Each variant's mass, energy and axial momentum budgets are below 1e-13 after one step and after 200 steps. Before the fix, 4 of 7 fail.
- **L2:**
  - F2 and F5(a) print the same text as before the fix (diff). Neither has a supply ring with a wall share.
  - F4 and F5(b) still pass. Both have such rings: F4's closed port, and F5(b)'s ring 3 at w = 0.434. Their states change.
  - Ring-centre runs never reach the change, so the C1 eqtt cases, criterion 4's runs and step D's drivers are unchanged. That is by construction (the branch needs `openAreaSupplies`) and is not rerun.
- **L3:** on tier 1, every F3 one-step energy budget is below 1e-13:
  - the dry run under `CRUCIBLE_T1_OPEN_AREA=1` at 76 × 8, 149 × 16 and 1192 × 128 (before: 1.944e-07, 1.906e-07, 2.941e-08);
  - the run Definition's, from `run development` (before: −1.011e-08).
- **L4:** inquiry's 2-level probe, with the same settings (`run development`, `CRUCIBLE_T1_LEVELS=2`, to 0.02 ms, 4 threads, Table B `rl10_s.bin`): budgets (a) gas and (c) whole system below 1e-11 at 0.02 ms (R3(a)'s limit; before: 2.45e-07).
  - If L4 fails while L1 to L3 pass, a second leak remains on the AMR side. It will be looked for next, and is not fixed under this registration.

**Predictions** (inferred, before the code):
- L1 to L3 come out at rounding, 1e-16 to 1e-15, as F6's passing variants already do.
- L4 passes. The probe's leak per base step, 2.45e-7 / 31 = 7.9e-9, is the size of the one-step leak, so nothing else is needed to explain it. Still, AMR adds plate faces at coarse-fine boundaries that one level does not test.

**Results** (measured on the Mac, 23 October, about 00:49 to 00:54 CDT; logs in `face_inlet/leak_2026-10-23/`). The fix is the one line registered. L1 to L3 pass. L4 could not be judged as registered: the fixed probe stops at 14 µs, in the spray, before 0.02 ms.
- **L1, pass.** F6: all 7 variants have their largest |budget| between 1.9e-16 and 7.7e-16 (limit 1e-13). Each energy budget after 200 steps is 7.5e-16 or less, against 4.3e-6 before (`f6_after.log`).
- **L2, pass.**
  - F2's output is identical to before. In F5's output only the two open-area lines of F5(b) differ, as expected (budgets mass 2.1e-16, energy −3.9e-16).
  - F4 and F5 pass. The whole `crucible_face_inlet_tests` passes with 0 failures (`all_after.log`).
- **L3, pass.** F3's one-step energy budgets:
  - the dry run: −6.3e-16 at 76 × 8, 5.4e-16 at 149 × 16 and −1.2e-15 at 1192 × 128;
  - the run's Definition: 5.0e-16 (before: −1.011e-08).
- **L4, not judged.** The fixed probe stops in the first base step past 14.0 µs with `Flow::exchange left a cell non-admissible`.
  - Temporary prints, reverted and not committed, located it. The gas there is finite: finest level, cell (0, 21), ρ 0.0241 kg/m³, T 227.0 K. The spray's deposit is NaN.
  - It comes from one LOX drop: T 65.085 K against T_sat 65.137 K at the cell's 2404 Pa, not boiling, d 0.49 mm. Its Heun stage gave m1, h1 and T1 NaN.
  - This is the drop just below saturation at the start-up's low pressure. The mass-transfer number's 1 − Y_s goes to zero there (inferred, not traced to the line). It is a defect in the spray's start-up path, not in the gas budget. The fix only changed the face cells' states, which moved this drop onto it. The 0.02 ms run on the unfixed build did not meet it.
- **In place of L4** (added after the registration, so not a registered check): the same probe to 0.01 ms, before the drop, on both builds.
  - The gas budget (a) is 2.01e-07 before the fix and 8.50e-16 after it. The whole system (c) is 2.01e-07 and 8.51e-16 (limit 1e-11; `before/probe2_10us.log`, `after/probe2_10us.log`).
  - So on 2-level AMR the gas budget closes to rounding once the fix is in. To 0.01 ms there is no second leak (measured).
- **L4 after the rebase onto main (8fe9ce3), pass** (measured on the Mac, 23 October, 00:57 CDT; `rebased_before/` and `rebased_after/`). Main's start amendment opens the H2 supplies at 0 and the LOX at t_L = 20 ms, so the probe to 0.02 ms now runs without drops and does not meet the drop above. It is the registered check as written: inquiry's 2-level probe, run development, to 0.02 ms.
  - Fixed: (a) gas 1.58e-15 and (c) whole system 1.58e-15 (limit 1e-11), exit 0, 32 base steps, 18 regrids.
  - Unfixed (the line reversed and rebuilt on the same rebased tree, then restored): (a) 2.39e-07 and (c) 2.39e-07, with the run's F3 one-step energy budget −1.011e-08 as before.
  - The drop's NaN is not fixed by this. It is a spray defect that waits for the first LOX at low pressure and is handed to the spray's owner (SESSION_HANDOFF.md).
