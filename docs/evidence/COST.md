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
- **Low's fallbacks may point at a defect, not only a cost.** Low's T_max reached 3,634 K at 2.4 ms. The adiabatic flame temperature at 1 MPa is about 3,400 K at most (measured now with Cantera, h2o2.yaml, H2/O2 from 300 K: 3,392 K at phi 1, 3,374 K at phi 1.2; 3,367 K at phi 1 from 150 K). The chamber pressure was 1.00 MPa then (measured, injector probe), and the run's global mass and energy budgets close to about 1e-14 (measured). So some cells sit about 240 K above the flame temperature. Differential diffusion of H2 can lift a diffusion flame above it (inferred), but blob, with the same transport, peaks at 3,396 K. The cause is not found.

## Levers

The ones that pay, largest first. Speedups are on blob's wall per simulated ms.

1. **Adaptive resolution with level subcycling** (the Lightweight engine's main line; core/amr.* is inquiry2's B1). The fine radial cells are needed in the shear layers near the injector, not across the whole chamber. Each coarser level also takes a step twice as long. As an example, level 2 (0.4 x 2 mm) on r < 8 mm and z < 150 mm, level 1 on r < 16 mm and z < 300 mm, and level 0 (1.6 x 8 mm) elsewhere: 3,300 cells instead of 16,730, and 3.9e7 cell updates per ms instead of 2.9e8. That is **7.5 times on the gas step**, about 4 to 6 times after regridding and ghost-cell overhead (derived; the region sizes and the overhead are guessed). On its own it gives about **2 to 2.5 times overall**, because the spray stays at 29% (derived).
2. **Spray subcycling.** The parcels are advanced every gas step of 5e-8 s. Their own time scales are much longer: drag relaxation about 1e-5 s for a 5 um drop, a 5 um drop's lifetime about 2.5e-5 s, crossing a 0.4 mm cell about 1e-5 s (derived; the drop size and the evaporation constant guessed). Advancing them every N = 10 to 20 gas steps, with their sources applied as rates in every gas step, cuts the spray to 1.5 to 3%: **1.35 to 1.4 times alone** (derived from the measured share), and with lever 1 about **5 to 8 times overall**. The sources stay inside the shared step as rates; N is chosen by measurement, from S1 to S5 rerun at N = 1 and N.
3. **Low's over-temperature first, then the fallback.** Find why cells exceed the flame temperature. If it is a defect, its fix removes the fallbacks. If the hot cells are physical, extend the T axis above 3,500 K, which needs thermo fits beyond the mechanism's. Expected: low back to about 38 ms per step, **1.5 to 6 times on low while fallbacks occur** (derived from the measured 58 to 220 ms per step at 0.8 to 3.6 calls). This is a correctness item before it is a cost item.

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

## Lever 4, for later: implicit acoustics in physical time (a note, no code)

Added at the Director's request (14 October, about 01:25 CDT) as the candidate for the last 4 to 7 times, with the GPU as the alternative.

- **What it is.** An IMEX step in physical time: convection, diffusion, the turbulence model and the sources stay explicit; only the pressure terms that carry the acoustics are implicit. Each stage solves one elliptic (Helmholtz) equation for the pressure, with coefficients from the cell's rho a^2, then fluxes the pressure work and momentum with the new pressure. The conserved variables stay conservative. This is the all-speed, semi-implicit pressure family, for example Kwatra, Su, Gretarsson and Fedkiw (2009, J. Comput. Phys. 228), which keeps a density-based conservative scheme and a general equation of state. The march stays physical (step 1): there is no pseudo-time, and nothing is iterated to convergence except one linear solve per stage.
- **Gain on the step.** The step becomes the convective one: 7 to 30 times longer here (derived above, the radial speeds guessed). Everything else per step (transport, Table A, the spray) then runs 7 to 30 times less often.
- **Cost per step.** One Helmholtz solve per stage, two stages with SSPRK2, on 16,730 cells with a 5-point stencil. With conjugate gradients and a cheap preconditioner, or a geometric multigrid on the block grid, each solve is guessed at 0.15 to 0.5 of today's explicit step, which carries mixture transport, gradients and the spray per cell. So a step costs about 1.3 to 2 explicit steps (guessed).
- **Net:** about **3.5 to 23 times** (derived from the two ranges above), and about 5 to 10 times as a central guess, because the larger steps also cut the Strang chemistry and spray calls. On top of levers 1 and 2 that brings 15 ms to about 2 to 7 minutes on the Mac (derived, all inputs guessed).
- **What it costs in physics.** Acoustic waves at steps above the acoustic CFL are damped and dispersed, so the run would not resolve chamber acoustics or combustion instability. That fits the performance-prediction framing (honest bands, not instability CFD). The fill and pressure rise (about 1 ms time constants) are still resolved. In a throat, u is about a, so the convective and acoustic steps coincide and the gain there is about 1; with level subcycling the throat level keeps its own step.
- **How it sits with the plan.** It is an implicit physical step, so the decision rule in TECHNICAL_PLAN's *Lightweight engine* applies. That rule names the wall strip's acoustic limit; here the binding limit is the bulk's, so the rule would need restating before this is built. That is the Director's and Ben's call.
- **Against the GPU** (20 to 200 times over the Mac, guessed in TECHNICAL_PLAN's *Compute budget*): the GPU keeps every acoustic step and buys throughput; lever 4 removes steps. They multiply rather than compete, but lever 4 is the smaller port.

## Decision (the Director, 14 October, about 01:25 CDT)

1. Lever 3 first, as a correctness item: find why low's cells sit about 240 K above the flame temperature while blob's do not (spray energy coupling, parcel source clipping or differential diffusion), with a one-cell budget on a hot cell.
2. Then lever 2, spray subcycling, with N chosen by measurement on S1 to S5.
3. Lever 1 is inquiry2's B1; the Director points inquiry2 at this note.
4. Lever 4 above: a note only, no code.

## Lever 3: where low's hot cells are (the trace, 14 October, about 02:10 CDT)

A rerun of low with `CRUCIBLE_HOT_TRACE=3450` (e18c710; it reproduces the original run step for step: 38,577 steps and T_max 3,634 K at 2.4 ms) crossed 3,450 K at step 37,049, t = 2.311 ms. Evidence in `cost/`: `low_hot_trace.txt` (the hottest cell every 200 steps, and every step around the crossing, by `hot_read.py`), `low_hot_cells.csv` and `low_hot_cells.txt` (the 5 by 7 cells around it at the crossing, by `hot_cells.py`). `excess` is the specific total enthalpy above the propellants' mixing line at the cell's hydrogen element fraction f_H (H2 at 287 K, the liquid at its injected enthalpy). Under equal diffusivities, with adiabatic walls and a spray that only takes heat from the gas, it has no source. All values below are measured.

- **Where.** Every cell above 3,450 K in the trace (the crossing and the 2,000 steps after it, to 2.43 ms) is in the first axial column against the injector face (z 0 to 2 mm), at r 0 to 1.2 mm. That is inside the LOX post's radius (2.5 mm) and inside the liquid core (13.7 mm long, §3 of the prediction). The real jet has liquid there; the RR mode releases its drops from the core's surface (up to 1,174 parcels in the cell at r 2.4 to 2.8 mm) and leaves a gas volume inside the core.
- **What the gas there is.** Hydrogen-rich: f_H 0.14 to 0.29 in the face column against 0.112 stoichiometric, Y_H2 up to 0.195 at the axis corner. Hydrogen reaches it along the face across the post's lip (the face cell at r 2.4 to 2.8 mm: f_H 0.226, Y_H2 0.127, 2,611 K), and the flow there runs toward the axis (u_r −6 to −22 m/s, u_z 27 to 41 m/s).
- **How far above the line.** The hot face cells sit 3.7 to 5.6 MJ/kg above the mixing line, against 0.5 to 0.8 MJ/kg in the flame's hottest cells elsewhere before the crossing (3,397 to 3,413 K, f_H 0.10 to 0.125, 7 to 21 mm downstream). Cells by the spray sit below it (down to −1.2 MJ/kg), as the heat held in the liquid requires. Mixing alone cannot raise a cell above the highest excess around it, so the face zone has gained enthalpy: by transport with unequal diffusivities (H2 and H diffuse faster than heat), or by a numerical source at the wall and axis corner.
- **Spray coupling excluded here as well as by reading the code.** The hot cells hold no parcels and receive no deposit (`Spray::exchanged()` zero in every traced step).
- **Blob** injects through the post's area, so the gas inside its core is O2 vapour and its peak stays at about 3,400 K.

**Not yet done: the one-cell budget.** It splits the energy into each face's inviscid, conductive and species-diffusion enthalpy fluxes for a face-column cell, such as the one at r 0.4 to 0.8 mm, over the steps before the crossing. A large species-diffusion term means differential diffusion; a term that does not close means a numerical source. The fluxes live in `core/flow.*` (inquiry2's files), so the budget is either a read-only flux hook agreed with inquiry2, or a driver-side recomputation from the state through `Medium::transport`. **The remedy waits on it** (inferred, not decided): if the zone's enthalpy is physical transport, the defect is the open gas volume inside the core, and the core's volume would be blocked or fed from the core itself; if it is numerical, the corner is fixed.

**The one-cell budget's readout, registered before its code (inquiry2, 20 October, 17:15 CDT, commit 8652024; the stamp first read 17:45 in error).** For the face-column cells (i 0, j 0 to 2) over 2.20 ms to the 3,450 K crossing, with b the trace's excess and its rate split per face and mechanism as F_E − h_LOX F_m − (h_H2 − h_LOX) F_H plus V dp/dt − d(ρkV)/dt and the sources: (a) ρEV's change must match the probe's face and source terms to 1e-9 of the largest term, or the cause is numerical; (b) of the positive supply of b (the face-and-mechanism integrals that add b to the cell), 70% or more from molecular diffusion through interior faces means differential diffusion, physical within the model (say so and stop); 70% or more from turbulent diffusion means the Sc_t 0.7 against Pr_t 0.9 mismatch, a model inconsistency (report; the fix is a model choice); 30% or more through the plate or axis faces, from a source that should be zero, or from inviscid faces whose b per unit mass flux lies outside the range of the two cells beside them, means numerical (fix it); anything else is reported as mixed, with the shares.
