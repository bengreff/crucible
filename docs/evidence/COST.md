# Cost of a run: what sets the step and where the wall time goes (14 October 2026)

Asked by the Director at about 00:45 CDT, after the S7 Mascotte runs took hours on the Mac for about 15 ms of physical time. Nothing here is implemented; it is a note to decide from. Labels: measured, derived (stated arithmetic), inferred, guessed.

## The case measured

- S7 Mascotte A-10 (SPRAY.md S7), axisymmetric, 239 x 70 = 16,730 cells: dz 2 mm, dr 0.4 mm, uniform; R 28 mm, L 478 mm. CFL 0.8, 6 threads. HLLC, SST, Table A local equilibrium, Lagrangian LOX parcels.
- Two runs on the Mac, profiled while running: **blob** (t 3.3 to 3.6 ms) and **low** drop sizes (t 2.1 to 2.4 ms).
- The Mac was loaded (both runs, a pipe test, other users' jobs), so absolute wall times are upper bounds for a quiet machine. The shares are what matter.

## What sets the step

**The acoustic CFL across the 0.4 mm radial cells** (measured by grid scaling, first 0.1 ms of the same case):

| Grid | dz, dr | dt at 0.1 ms | Ratio |
|---|---|---|---|
| 239 x 70 | 2 mm, 0.4 mm | 7.2e-8 s | 1 |
| 239 x 35 | 2 mm, 0.8 mm | 1.23e-7 s | 1.7 |
| 478 x 70 | 1 mm, 0.4 mm | 6.17e-8 s | 0.86 |

- `Flow::stableDt` takes dt = CFL V / (sum over faces A(|u_n| + a) + 2 nu_max sum A^2 / V). With a much larger than u, the acoustic part is proportional to 1/(1/dr + 1/dz): doubling dr gives 3.0/1.75 = 1.71 (derived), against 1.7 measured. A diffusive limit would give 4. So the step is acoustic, set by dr.
- The diffusive limit is far away: dr^2/(4 nu) is about 1e-4 s for nu of a few 1e-4 m^2/s in the hot gas (derived, nu guessed).
- dt falls as the chamber heats and the sound speed rises: 7.2e-8 s at 0.1 ms, 5.7e-8 s at 1 to 2 ms, 4.55e-8 s at 4 ms (measured, blob).
- **Steps per ms:** about 13,900 at first, 17,500 at 1 to 2 ms, 22,000 once hot (measured). A 15 to 25 ms run is about 300,000 to 500,000 steps.
- **How far below the flow's own limit:** the H2 jet enters at 317 m/s. With the acoustic term removed, the step would be set by |u|: about 5e-7 s if the full 317 m/s crossed the 0.4 mm cells, about 1e-6 to 2e-6 s for that speed axially across 2 mm plus 30 to 100 m/s radially across 0.4 mm (derived; the radial speeds guessed). So the acoustic step is **7 to 30 times** shorter than the convective one. The chamber is low Mach (about 0.1 in the bulk, guessed; at most about 0.25 in the jet).

## Where the wall time goes

**Wall per step** (measured from the run logs): blob 33.7 ms; low 106.8 ms on average over 2.1 to 2.4 ms. Blob runs at about 600 to 700 s per simulated ms, so 15 ms is about 2.5 to 3 h.

**Profile** (measured: macOS `sample`, 20 s at 1 ms per thread, innermost-phase attribution by `cost/parse_sample.py`; outputs in `cost/blob_profile.txt` and `cost/low_profile.txt`):

| Phase | Blob, % of busy | Low, % of busy |
|---|---|---|
| Spray (all of `Spray::advance`, loads included) | 29.3 | 16.6 |
| Viscous and diffusive fluxes | 16.2 | 5.9 |
| Transport properties (mu, lambda, D_k) | 14.2 | 6.0 |
| Inviscid fluxes (HLLC, reconstruction) | 11.6 | 4.3 |
| Gradients | 10.2 | 4.1 |
| State refresh (T from e, sound speed) | 7.5 | 3.1 |
| SST | 5.2 | 2.6 |
| Table A lookup | 3.4 | 1.1 |
| Step limit (`stableDt`) | 1.7 | 0.5 |
| Cantera equilibrium (Table A fallback) | 0 | 55.4 |
| Busy share of all thread samples | 92% | 39% (61% waiting) |

- **Blob:** the gas step is 71% and the spray 29%. Inside the gas step, transport (properties, gradients, viscous fluxes) is 40% of busy, more than the inviscid step. The spray share rises with the parcel count (27% by the run's own timer at 5,700 parcels).
- **Low:** the fallback dominates. Cells above Table A's T axis (300 to 3,500 K, the mechanism's thermo range) call Cantera directly. One call costs about 75 ms (derived: the extra wall between the 2.4 to 2.5 ms and 2.5 to 2.6 ms intervals over the extra calls, (375 - 99) s / (5,081 - 1,412)). At 0.8 to 3.6 calls per step, the other threads wait, which is the 61%.
- **Low's fallbacks point at a defect, not only a cost.** Low's T_max reached 3,634 K at 2.4 ms. The adiabatic flame temperature at 1 MPa is about 3,400 K at most (measured now with Cantera, h2o2.yaml, H2/O2 from 300 K: 3,392 K at phi 1, 3,374 K at phi 1.2; 3,367 K at phi 1 from 150 K). The chamber pressure was 1.00 MPa then (measured, injector probe), and the run's global mass and energy budgets close to about 1e-14 (measured). So some cells are locally about 240 K too hot. The cause is not found. Blob's T_max is 3,396 K, at the flame temperature.

## Levers

The ones that pay, largest first. Speedups are on blob's wall per simulated ms.

1. **Adaptive resolution with level subcycling** (the Lightweight engine's main line; core/amr.* is inquiry2's B1). The fine radial cells are needed in the shear layers near the injector, not across the whole chamber. Each coarser level also takes a step twice as long. As an example, level 2 (0.4 x 2 mm) on r < 8 mm and z < 150 mm, level 1 on r < 16 mm and z < 300 mm, and level 0 (1.6 x 8 mm) elsewhere: 3,300 cells instead of 16,730, and 3.9e7 cell updates per ms instead of 2.9e8. That is **7.5 times on the gas step**, about 4 to 6 times after regridding and ghost-cell overhead (derived; the region sizes and the overhead are guessed). On its own it gives about **2 to 2.5 times overall**, because the spray stays at 29% (derived).
2. **Spray subcycling.** The parcels are advanced every gas step of 5e-8 s. Their own time scales are much longer: drag relaxation about 1e-5 s for a 5 um drop, a 5 um drop's lifetime about 2.5e-5 s, crossing a 0.4 mm cell about 1e-5 s (derived; the drop size and the evaporation constant guessed). Advancing them every N = 10 to 20 gas steps, with their sources applied as rates in every gas step, cuts the spray to 1.5 to 3%: **1.35 to 1.4 times alone** (derived from the measured share), and with lever 1 about **5 to 8 times overall**. The sources stay inside the shared step as rates; N is chosen by measurement, from S1 to S5 rerun at N = 1 and N.
3. **Low's over-temperature first, then the fallback.** Find why cells exceed the flame temperature; a correct solution at 1 MPa should not leave Table A. If hot cells remain legitimate, extend the T axis above 3,500 K, which needs thermo fits beyond the mechanism's. Expected: low back to about 38 ms per step, **1.5 to 6 times on low while fallbacks occur** (derived from the measured 58 to 220 ms per step at 0.8 to 3.6 calls). This is a correctness item before it is a cost item.

Smaller ones:
- **Tabulated transport properties:** at most **1.17 times** (measured share 14.2%). TECHNICAL_PLAN's Lightweight engine lists this as the first RL10 candidate pending a profile; the profile now bounds it. Gradients and viscous fluxes (26%) cannot be tabulated.
- **A coarser startup grid:** 239 x 35 gives 1.7 times the step (measured) on half the cells, about 3.4 times per ms (derived, assuming the cost per cell update is unchanged). Run there until the start is forgotten, then remap. About 1.5 to 2 times overall if two thirds of the physical time runs coarse (derived). It is lever 1's base level used alone, and needs a remap.

Not levers for this case:
- **Dual time stepping, or any implicit physical step** (TECHNICAL_PLAN step 16, parked). The physical step could grow by the acoustic-to-convective ratio, 7 to 30 times. The plan's own inner cost is 10 to 30 iterations of 2 to 5 explicit-step costs, so 20 to 150 explicit steps per physical step: **about 0.05 to 1.5 times** (derived). Dual time pays where the ratio is about 1,000, the wall strip. The decision rule names that case, and it does not hold here.
- **Low-Mach preconditioning on its own** changes the pseudo-time wave speeds, so it speeds up convergence inside dual time or towards a steady state. It does not lengthen an explicit physical step; applied to the explicit march it shortens it, like M^2 (Birken and Meister 2005, cited in `LOW_MACH.md`, which rejected it on 30 September).
- **Local time stepping towards the steady state** is the largest in principle (10 to 100 times, guessed: a steady RANS state in thousands of iterations instead of 300,000 steps). It replaces the physical march with pseudo-time. Step 1 forbids that, and Ben rejected it on 5 October. It is listed only so the choice is explicit.

## Which fits CRUCIBLE's design

**Levers 1 and 2 together** (inferred from the design documents). Adaptive resolution is the Lightweight engine's main line. It keeps the physical march of step 1 and attacks the binding constraint here: the fine radial cells set both the cell count and the step everywhere. Spray subcycling is needed with it, or the spray caps the gain at about 2 to 2.5 times. Lever 3 is a fix and goes first in any case.

**What it leaves** (derived from the rows above): blob from about 650 s per ms to about 80 to 130 s per ms, so 15 ms in about 20 to 35 minutes on the Mac. That is still 4 to 7 times over the 5-minute target. The rest is the explicit acoustic step itself, 7 to 30 times below the flow's own. Closing that needs a pressure-based (all-speed) treatment of the acoustics, pseudo-time (Ben's ruling), or the GPU. That is a later decision, once levers 1 to 3 are measured.
