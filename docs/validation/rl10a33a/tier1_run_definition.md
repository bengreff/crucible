# RL10A-3-3A tier 1: the run definition (draft)

Status: drafted 21 October 2026, about 02:40 CDT, before any tier-1 run. It holds no predicted number and reads no measured value. It is drafted as far as it goes without the face inlet (a), the plan's decision 6. Items marked *proposed* are engineering, open to inquiry2 and the Director until the prediction commit; a change after any tier-1 run is dated and reruns every member already run (the plan's *Pre-registration*). It is governed by `tier1_skeleton.md`, `tier1_ensemble_plan.md` (decisions 1 to 8), `tier1_setup.md` (the dry run) and the protocol's §2 to §4.

## What the run is

`crucible_rl10_tier1` gains a run mode: the dry run's Definition, the LOX injectors and the chemistry, stepped from "valves open, ignite" until settled (protocol §3). A deviation member is the same run with one input changed, started warm from the central case's settled state (the plan, *How a member runs*). The dry run's checks G1 to G8 run first in every run, and a failed check stops it.

## Fixed by the dry run and the plan's decisions (no open choice)

- **Geometry:** the PCHIP contour, the face at the station nearest −13.1 in, the 8 rings of choice 1 and the porous face of choice 3 (`tier1_setup.md`, *What the dry run builds*, and the two member rules).
- **Flows:** oxidizer 13.9525 kg/s; fuel 2.76014 kg/s, 10% transpired at a uniform mass flux over the porous face, and the rest through the 16 H2 slots in equal shares.
- **The face inlet:** (a), every supply at its true annulus (decision 6). Each LOX slot is a closed port (zero gas flow), since the liquid enters there as parcels (the handoff's spec).
- **Temperatures:** fuel 199.67 K, oxidizer 101.17 K, chamber wall 567 K, and the face's wall (the post tips) at 199.67 K (decision 8).
- **Supply turbulence:** intensity 5%; length scale 0.07 D_h at the slots and a pitch on the porous face (decision 7). Each supply's viscosity ratio is converted at the nominal face state (G4).
  - The ratio does not depend on that nominal state's pressure. It is ρk / (ωμ) = √(3/2) I C_μ^¼ L ρu / μ, and ρu is the supply's fixed mass flux (*derived*). The dry run gives the same ratios at 2 and 4 MPa (*measured*), so the nominal 2.8 MPa sets no tier-1 input. `Flow` rebuilds ω from the ratio at the face's own state, so L holds there.
- **Turbulence and wall:** SST-2003 with the Hasan-Pecnik corrections and the wall functions (`Definition::turbulence`, as the dry run sets it).
- **The LOX:** one injector per ring across its LOX slot, at 1/8 of the oxidizer, at 101.17 K, axial (no swirl, choice 2). Its speed is the mean (`velocity` 0: ṁ / (ρ_l A_L), 6.05 m/s at the nominal state). Drops are Rosin-Rammler with SMD 447 μm and q 3.207, released over Porcheron's 6.84 mm core (G8). Wall restitution is 0.5, Mascotte's (`MASCOTTE_A10_PREDICTION.md`; *guessed* there).
- **Chemistry:** Table B (decision 1), at the RL10 range's Δη 0.1 table with the corner drop and the induction cut; CVODES takes every refused cell, counted by cause.
  - The table is rebuilt in the repo by `docs/evidence/table_b/build_rl10_table.sh`, not taken from `/tmp`. That script holds the 01:00 build's arguments and environment verbatim (they were only in `/tmp/rl10s/build.sh`).
  - **Checked at about 02:43 (measured):** the rebuild's data after its two-line provenance header is byte-identical to the table the RL10 checks used (3,591,631 bytes; SHA-256 `14f7af124e4509f22eade1a6a03904080d3fbedd3a42d9406960b2f824a1a404` for both; headers of 732 and 756 bytes). The header differs by 24 bytes, exactly the mechanism's absolute path in place of `../`, plus the build time and seconds. The node report holds the same rows (sorted, identical), written in thread order.
  - So the prediction commit records the hash of the data after the header, not of the file.
  - Table A (shifting equilibrium) is member 15.
- **Settled:** F_vac, the injector-face pressure and the outlet mass flow each drift by under 0.1% over the last 1 ms (protocol §3). The drift is (max − min) / |last| over that ms, as `chamber_study` measures it (`tests/chamber_study.cpp:805-815`). That driver also requires |outlet / supply − 1| < 1e-3; tier 1 counts the LOX injected in the supply, and the parcels' inventory is steady once settled.
  - **Replaced by the Director's decision (22 October 2026, about 22:40 CDT, recorded then, before any tier-1 output exists).** Every run goes to the 60 ms ceiling and is judged over its last 20 ms, 40 to 60 ms.
    - **Settled** if two conditions hold, each for the chamber pressure and for c*:
      - the batch-means standard error of the 40 to 60 ms mean is at most 0.2% of that mean;
      - the means of the two 10 ms halves (40 to 50 and 50 to 60 ms) differ by at most 0.5% of the 20 ms mean.
    - **The chamber pressure is the injector-face pressure,** which is protocol §3's settled pressure (*read so*). c* is p0 A_t / ṁ, from the nozzle stagnation pressure, as item 3 states.
    - **The batches (engineering, fixed now).**
      - τ is the integral time scale of each quantity's 0.5 ms-printed series over 40 to 60 ms: the autocorrelation integrated to its first zero. The run also writes the series at every step it already samples, so τ is not limited by the print interval.
      - The batch length b is the shortest of 1, 2, 4 and 5 ms that is at least 10 τ. That gives 20, 10, 5 or 4 batches.
      - The judged SE is the larger of the SE with b and with the next longer length in that list. At b = 5 ms there is no longer length, and the 5 ms SE stands.
      - If 10 τ exceeds 5 ms, the 5 ms SE is used and flagged as "batches shorter than 10 τ".
      - This mirrors S7's guard (MASCOTTE_A10_PREDICTION.md §4, *The Director's decision on low and high's settle test*).
    - **The reported values** are each item-3 quantity's 40 to 60 ms mean, with its judged SE. F_vac and the outlet flow are judged by no bar, but their SE is printed.
    - **A run that does not meet the rule at 60 ms** is reported with these values and flagged unsettled. Its SE enters the quantity's range as a component, alongside protocol §4's bracket deviations. It is not extended without a dated note.
    - Protocol §3's 0.1% test over the last 1 ms is printed beside it, not judged.
    - **Item 5's parcel-rate check, re-expressed (engineering).**
      - After 60 ms, the central case continues warm at twice the parcel rate to 85 ms, and is judged over 65 to 85 ms by the same rule.
      - The rate doubles for every run if any item-3 quantity's 20 ms mean moves by more than 0.5% (the halves bar). It also doubles if the move exceeds twice the two runs' combined SE.
      - This adds 25 ms at about twice the spray's cost to the central case alone (*derived*). It is in the restated budget.

## Proposed (engineering, stated before the code)

1. **The start.** Supplies and injectors open at t = 0 on the dry run's 0.5 ms ramp.
   - The igniter is the C1-like check's scaled torch: 7 kJ over 1 ms (`TABLE_B.md`, *The C1-like check at RL10 conditions*), which lit rl10dev under Table B. The dry run's 0.5 J placeholder (rl10dev's) is replaced.
   - The igniter's real dimensions are missing (`tier1_skeleton.md`, *Open*). It is a start device: the settled state should not depend on it if that state is unique. The plan's one cold-start check (member 4's 834 K end) tests uniqueness, not the igniter. *Proposed:* a second cold start of the central case at half the igniter energy, judged by the same 0.1% agreement. One run; it is added to the budget only if the Director wants it.
2. **The end.** Settledness is checked every 1 ms from 5 ms on. The ceiling is 60 ms of flow time (*guessed*: rl10dev settled within 20 to 25 ms pre-vaporized, and the liquid's residence adds to that). A run that has not settled by the ceiling is "not measured" (the plan).
3. **The quantities,** printed every 0.5 ms and averaged over the last settled ms:
   - F_vac;
   - Isp_vac on the chamber flow and on the engine-inlet flow;
   - the nozzle stagnation pressure and the injector-face pressure;
   - c* = p0 A_t / ṁ;
   - the hot-gas wall's heat flow, and the face's apart.
   - **The face's heat flow is not part of the hot-gas wall's** (*inferred*). Protocol §3 names the hot-gas wall's heat flow as the expander cycle's driver, which is the jacket's coolant. The transpiration-cooled face is not on the jacket. Under decision 8 the face is a cold wall, so including it would bias the wall heat flow.
   - Also printed: Table B's causes by region, the parcels in flight, the liquid mass, and the spray's share of the step's cost.
4. **The mesh.** The base is 149 × 16 (G1 passes: the face at −13.0862 in) with 3 levels (1192 × 128 at the finest).
   - **The tags are AMR step C's** (`Amr::setTags`, on `main` since af91015; its checks 2, 4 and 5 are inquiry2's), with step D's thresholds fixed on C1:
     - the throat tag (0.5 r_t, to level 2);
     - the wall strip to the y+ target (300 by default; step D and the wall functions' y+ question decide 300 or 956, which is the budget's factor of 4);
     - the dynamic tags on f_H, T, the pressure jump and vorticity, with regridding every 2 base steps and a buffer of 2.
   - **Added for tier 1:** a static face tag. Every face-column cell from the axis to R_c, out to the core length (6.84 mm, 5.4 finest cells of 1.26 mm), is held at the finest level. This is because decision 6's declared limit (streams mixed within one cell) falls with the face cell's height, and because the drops are released there. Its cost is one strip of 128 × 6 finest cells, small against the wall strip (*derived*).
     - Step C's `Tags::custom` holds it with no engine change: true on level l's cell (i, j) when its centroid lies within 6.84 mm of the face, from `level(l).mesh()` (`core/amr.cpp:364`; read 21 October, about 02:55).
5. **The parcel rate.** 5e5 parcels per second per injector, Mascotte's (*guessed* there, for cost), so each parcel carries 3.5 mg (*derived*: 1.744 kg/s / 5e5), about 69 drops at the SMD.
   - A numerical check, not a member: once settled, the central case runs 1 ms more at twice the rate. If any quantity moves by more than the settle bar (0.1%), the rate doubles for every run.
   - The parcels per occupied cell near the face are printed at the end of the first ms. A bar for them is not stated (no basis yet).
6. **Threads and machines.** 4 threads per run (the budget's basis), through the slot limiter on the Mac, or on backhouse in tmux with the thread caps of the worker rules. Each run is detached and writes a done-file.
7. **The warm start.** The central case writes its settled state (`CRUCIBLE_STATE_DUMP`'s format, or AMR's own if that differs) for the deviations. The re-meshed members (2, 3, 12, 13) start cold unless a mapping exists then (the plan).

## Frozen inputs (22 October 2026, about 22:35 CDT, before any tier-1 run)

The Director's brief of 22 October gated the freeze on S7 passing or its failure being understood. S7 did not pass, and the failure is understood (`MASCOTTE_A10_PREDICTION.md` §4, 22 October). No tier-1 run has been made, and no measured RL10 value has been read.

**Frozen.** Every physical and model input listed above under *Fixed by the dry run and the plan's decisions*, plus these items from *Proposed*, adopted as written (engineering):
- **The start (item 1):** the 0.5 ms ramp, and the 7 kJ over 1 ms torch. The second cold start at half the igniter energy is not run unless the Director asks for it.
- **The quantities (item 3)**, the face's heat flow kept apart.
- **The parcel rate (item 5):** 5e5 per second per injector, with its doubling check.
- **The threads (item 6)** and **the warm start (item 7)**.

The files that carry these inputs, pinned by SHA-256:

| File | SHA-256 |
|---|---|
| `tests/rl10_tier1.cpp` (the Definition, the LOX injectors, the env defaults) | `b8de323ebd2aaee4c96e7ccb009eb631f56b5aa4b4ca8df8c2ca74d47fa88789` |
| `tests/rl10_geometry.hpp` | `f4370ff953a13d54d8f60bf4e97212f96cad21e29cdaa09a9c52808a98e03682` |
| `tests/rl10_wall.hpp` | `44efe872e792e77277e04d278ab34414e97a82c1f5357694b6fbfff380639941` |
| `tools/mechanisms/h2o2_highT.yaml` | `b65874b8b5fb2038d85c5b751c58aabeefdf8d71879f0d0536b410dc85e9e6b5` |
| `tools/o2_saturation_nist.tsv` | `ab63f25c9690a1ebb44e41ed7884fd77589f7999076956377ebe71cc1c512e08` |
| `docs/evidence/table_b/build_rl10_table.sh` | `11f90f1f445e1b6aacbb0c4b262ba6f1f2369dded9e70e9ad90e685e2656308e` |
| Table B's data after its header | `14f7af124e4509f22eade1a6a03904080d3fbedd3a42d9406960b2f824a1a404` (checked 21 October, about 02:43) |

The run mode is still to be added to `tests/rl10_tier1.cpp`. It may change the driver file, but not one input value: the dry run's printed Definition and G2 to G4 and G8 lines must stay identical (`diff` against `tier1_dry_run_length_2026-10-21.log`), apart from the switch to face inlet (a), the LOX ports and the plate's temperature, which are decisions 6 and 8. From here, a change to any frozen input is dated and reruns every member already run (the plan's *Pre-registration*).

**Not frozen: the numerical settings, which wait on the Director.**
- **The mesh (item 4):** the y+ target (300 or 956) and the tag thresholds, which wait on step D's verdicts (D1's igniter, D2's measure; `docs/evidence/AMR.md`).
- **The settle rule's reading under a spray (below).**

**What S7 shows for tier 1.**
1. **The spray's cost (measured, S7's Rosin-Rammler runs, one injector at 5e5 parcels per second):** 33.4% and 41.3% of the wall time, so the whole step costs 1.50 and 1.70 times the gas step's. The blob legs cost 2.81 and 1.83 times.
   - The budget's guessed row (1.5 to 4) keeps its range: its low end now matches S7's measurement, but the RL10's parcels per cell are not known.
   - The central case's first ms measures the RL10's own ratio.
2. **The settle rule may not be reachable with a spray (S7 measured; the carry-over to the RL10 inferred).**
   - Over 36 one-ms windows each in S7's stationary runs (low and high), protocol §3's measure (max − min) / |last| never fell below 0.1% for the outlet flow: median 1.76% (low) and 4.85% (high).
   - It fell below 0.1% for the injector-face pressure once in 72 windows, and that with the pressure imposed at S7's outlet.
   - S7's 1 ms mean evaporation rate varies by 2.7% to 5.0%.
   - The RL10 is choked, and its face inlet mixes the streams within one cell, so its own noise is not known. S7's is mostly the flow's own, not the parcels' (inferred, §4).
   - *Proposed (protocol §3 is registered, so this is the Director's call):* the central run keeps the 0.1% test, and runs on to the 60 ms ceiling if the test is not met, so the history is there either way. The Director decides, before any tier-1 output is read, how a run that is stationary but unsteady is judged. One option is the S7 amendment's form: batch means over the stationary span, with their standard error against 0.1%.
   - The parcel-rate check's 0.1% bar needs the same reading.
3. **The face corner.** S7's low run holds a pocket of about 3,000 K in front of the LOX post (§4). Tier 1's LOX slots are closed ports with the plate beside them, and the face column is held at the finest level, so the case is not the same. The central run's first ms is checked for T_max at the face before it goes on (inferred to be worth a look, not a gate).
4. **A missing piece found while freezing: the spray does not run on AMR.**
   - `Spray` takes one `Flow` (`core/spray.hpp:116`). It locates parcels on a single uniform-dz mesh (SPRAY.md, *Mesh*).
   - `docs/evidence/AMR.md:18` leaves parcels on AMR's levels to the spray owner, which is this line, and nothing has been built.
   - Tier 1 needs 3 levels and the liquid LOX together, so **the central case cannot run as defined until the spray couples to `Amr`.** Coupling means: parcels located on the finest active level covering them, the gas read there, and the deposits booked to owned cells, with the composite budgets closed.
   - That design is registered before its code, with the spray checks S1 to S5 rerun on levels.

## The run mode: spec and acceptance (23 October 2026, about 00:11 CDT, before its code)

The Director's brief (23 October, 00:00): the driver's run mode is built under development mode, spec and acceptance first. The central run launches only when this acceptance passes, S7 low and high are scored, and the Director clears it. The blind protocol is unchanged: a development run's outputs are the model's, labelled as below, and no measured RL10 value is read.

**Invocation.** `crucible_rl10_tier1 run development` or `crucible_rl10_tier1 run central`. With no `run` argument the driver is the dry run, unchanged.
- **central:** the frozen run. The base is 149 × 16 with 3 levels on 4 threads. It runs to 60 ms, then the continuation at twice the parcel rate to 85 ms.
  - It refuses to start while the numerical settings are not frozen. A compiled constant, false until the Director freezes the y+ target and the tag thresholds, makes it exit non-zero before the Definition is built. Setting it is a dated commit.
- **development (protocol §6):** every output (each printed line's header, the file names, the series' first line) is labelled `DEVELOPMENT (protocol §6): not the reported case`.
  - The environment sets the run: `CRUCIBLE_T1_BASE` (default 149x16), `CRUCIBLE_T1_LEVELS` (3), `CRUCIBLE_T1_END` (s, 1e-3), `CRUCIBLE_T1_THREADS` (4), `CRUCIBLE_T1_OUT` (the files' prefix) and `CRUCIBLE_T1_RATE_FACTOR` (1).
  - The numerical settings are provisional: step C's `Tags` defaults (y+ 300).
  - The member inputs are read from the environment as in the dry run.
- **Table B** in both modes is `CRUCIBLE_MANIFOLD_TABLE`, whose provenance must name `h2o2_highT.yaml` (as `chamber_study` checks). The repo has no SHA-256 code, so the launch script checks the data hash against the frozen one with `shasum` and records it.

**The Definition** is the dry run's F3 Definition:
- face inlet (a) with the 8 LOX slots as closed ports (decision 6), and the plate at 199.67 K (decision 8);
- the supplies on their 0.5 ms ramp;
- the igniter's region as in the dry run, at 7 kJ from 0.2 ms over 1 ms (item 1).

**Gates, before any step,** on the base mesh:
- G1, G2, G3, G4 and G8 as in the dry run, and F3 on the run's own Definition.
- G5 and G6 test the old ring-centre rule, which decision 6 superseded; they are printed, not gated. G7 is printed as declared (decision 5).
- Any failing gate makes the run exit non-zero with no step.

**The coupling.**
- `Amr(definition, levels)` with step C's tags, plus item 4's static face tag. That tag is `Tags::custom`: true on a level below the top where the cell's centroid lies within the core length (6.84 mm) of the face.
- One `ReactingFlow` per level on Table B.
- One `Spray` on the `Amr` with the 8 injectors at 5e5 parcels per second each (times the rate factor), restitution 0.5.
- After each level's step: the spray's `afterStep(l, dt)`, then that level's `react(dt)`. On a redo: the spray's `redo()`.
- Each base step is bounded by the spray's next event and by the next sample time, so samples fall on a 10 μs grid.

**The quantities,** composite over the levels (each from its owned cells and faces), sampled every 10 μs:
- ṁ = Σ outlet flow.
- F_vac = Σ (inlet momentum flux + wall axial force + body axial force), plus the spray ledger's injection thrust and wall thrust.
- Isp_vac = F_vac / (g0 ṁ_s), on the chamber flow (ṁ_s = 16.7126 kg/s) and on the engine-inlet flow (2.7946 + 13.9525 kg/s).
- p0, by `columnStagnation` on level 0's last column before the cylinder's end (rl10dev's rule).
- The face pressure: the area average over every level's owned face-column cells. The owned area must equal the face's area to 1e-12, checked at every sample.
- c* = p0 A_t / ṁ.
- The hot-gas wall's heat flow, −(Σ wall heat − Σ plate heat), and the face's, −Σ plate heat.
  - This needs one additive engine change: `BoundaryRates` books the i = 0 share of the wall heat apart as well, exposed as `Measurements::plateHeatFlow`. Nothing else in the step changes.
- Printed every 0.5 ms, not sampled:
  - the liquid mass and the parcels in flight;
  - Table B's causes by level and by region (chamber to the cylinder's end, convergent to the throat, divergent);
  - T_max over the owned cells within the core length of the face;
  - the parcels per occupied cell within the core length (at 1 ms);
  - the spray's share of the wall time;
  - the budgets (a) to (c) as in A2.

**The settle judge** is the Director's rule as registered above. It runs at 60 ms over 40 to 60 ms, and at 85 ms over 65 to 85 ms with the continuation's move test. τ is taken from the 10 μs series.
- The continuation needs the parcel rate changed mid-run: an additive `Spray::setParcelRate(factor)`.
- **The outputs:**
  - the series (`<prefix>_series.csv`);
  - the state at 60 and 85 ms, or at the end in development (`writeAmrDump`, plus the parcels' bytes);
  - `<prefix>_done` with the exit status.
- **Not built tonight: a restart reader.** A crash in a multi-day run loses it. It is recommended before the central launch.

**Acceptance:**
- **R0:** the engine changes are additive. A0's dump md5 stays `e95251ab…`.
- **R1:** the dry run is unchanged. With no `run` argument, its output on the default meshes is byte-identical to the output of 543d50b's tree. The run mode's printed Definition, G2 to G4 and G8 lines match `tier1_dry_run_length_2026-10-21.log`'s 149 × 16 lines, apart from decisions 6 and 8 and item 1's igniter.
- **R2:** the gates stop the run. Each of these exits non-zero with no step:
  - (a) `run central`, while the numerical settings are not frozen;
  - (b) `run development` with `CRUCIBLE_T1_BASE=76x8`, where G1 fails (the face at −13.2542 in).
- **R3:** a short development run on the real configuration: 149 × 16 with 3 levels, 4 threads, through the slot limiter.
  - First a cost probe to 0.02 ms.
  - Then to 0.3 ms, through the igniter's start at 0.2 ms, if the probe puts that under 3 hours. Otherwise it goes to the longest whole 0.05 ms that fits, and the shortfall is recorded.
  - Met if:
    - (a) budgets (a) to (c) are each within 1e-11 at every print;
    - (b) every quantity is finite at every print, and the series' rows are finite on a uniform 10 μs grid (to 1e-12 s);
    - (c) the face-area check holds at every sample;
    - (d) the run reaches its end (its rejected steps printed).
- **R4:** the settle judge's selftest, on synthetic 10 μs series over 40 to 60 ms:
  - (a) a constant is settled, with SE and halves difference 0;
  - (b) a +1% step at 50 ms is unsettled, its halves differing by the expected 1% to 1e-9, and its batches flagged shorter than 10 τ;
  - (c) a sine of period 2 ms gives τ within 2% of P / (2π), so b = 4 ms, and is settled.
- **Two points found while writing the code, recorded before any run (23 October 2026, about 00:18 CDT):**
  - **R3(b) amended.** p0 is a mass-flux average over its column, and c* divides by the outlet flow. At the start the gas is at rest, so both flows are zero and p0 and c* are undefined there, not wrong.
    - They are written as NaN while their column's net forward flux, or the outlet flow, is not positive. Those samples are counted and printed.
    - R3(b) holds for every other quantity at every sample, and for p0 and c* wherever they are defined.
  - **`Spray::setParcelRate` takes the rate per second,** not a factor. The driver passes twice its rate at 60 ms.
- **R3's cost probe stopped at the first step; the cause and its fix, registered before the fix's code (23 October 2026, about 00:21 CDT).**
  - **The symptom:** with 3 levels the first step threw `Invalid wall-function layer` at a level-3 side-wall cell, with p = 0 and zero mass fractions (NaN wall density, cp −0). With 1 or 2 levels the same probe runs.
  - **The cause, in AMR's wall functions (inquiry2's lane, step D):**
    - Under B2 a level computes its primitives and properties only within distance 2 of its active cells (`Flow::transportProperties(2)`). `Flow`'s wall-solve loop (`core/transport.cpp:388-393`) still visits every wall face of the level.
    - Level 3's only active cells here are the face column (the static face tag), so most of its side-wall cells lie outside that region, and their arrays were never filled.
    - With 2 levels, level 2 holds the whole wall strip (the strip refines to level 2), so no wall face falls outside. That is why step D's runs never met it.
  - **The fix:** the loop skips a wall face whose cell is outside that region (`inRegion(q, 2)`).
    - No reader takes those faces' solutions. The fluxes use faces beside active cells, and the wall-strip tag reads active cells only (`core/amr.cpp:362-366`).
    - So no run that ran before changes. Without AMR, `inRegion` is always true.
  - **Criteria:**
    - (i) A0's md5 is unchanged;
    - (ii) the dry run is byte-identical;
    - (iii) the 2-level probe to 0.02 ms (149 × 16, 4 threads) writes a series byte-identical to its series before the fix;
    - (iv) the 3-level probe passes its first step.
  - Flagged to inquiry2 in the handoff.
- **The wall-face fix, scored (23 October 2026, about 00:23 to 00:28 CDT; `run_mode_2026-10-23/`).** All four criteria are met.
  - (i) A0's md5 is unchanged, `e95251ab…`.
  - (ii) The dry run is byte-identical.
  - (iii) The 2-level probe's series and AMR dump are byte-identical to those from before the fix.
    - The spray dump differed only in the 4 trailing padding bytes of each 152-byte `Parcel`, 63 bytes over 72 parcels. The dump wrote raw structs, so those bytes were never set.
    - The driver now zeroes them. The padding-zeroed dump differs from the old one in no other byte.
  - (iv) The 3-level probe passes the wall functions. It then stops before the first 10 μs sample, on the first of the two findings below.
  - The driver's t = 0 print now says "none before the first step" for the budgets, and leaves t = 0 out of the worst budgets. Before, it printed 1.00 there, because the budget has no reference until a step is taken.
- **Two findings block R3 (23 October 2026, about 00:28 CDT). Neither is fixed. R3 is not met, and nothing past the probes was run.**
  1. **The start-up's explicit evaporation overdraws a fine cell (the spray's coupling, this line).**
     - The 3-level probe stops before the first 10 μs sample: `Flow::exchange left a cell non-admissible`. On AMR the spray releases parcels at the end of each base step, so the first base step deposits nothing; the stop came in a later one (corrected 23 October, about 00:33 CDT; it first said "the first step").
     - The cell is level 3's (4, 28): active, its full ring 2.33e-7 m³ (1.26 mm by 1.02 mm at r of about 29 mm).
     - Before the deposit its gas density was 0.096 kg/m³. The chamber starts at 1,000 Pa and is still filling with H2.
     - One step's deposit is 7.16e-8 kg of O2 vapour, about three times the cell's gas, and −0.0323 J of energy. That leaves ρE at −1.57e5 J/m³, below the floor: the gas gave the drops more heat than it held.
     - The spray has no bound on a step's exchange against its cell's content.
     - At 2 levels each cell holds 4 times the volume, and the probe ran to 0.02 ms with 72 parcels. Whether that holds to 0.3 ms is not known.
     - Options, none chosen:
       - (a) A step bound from the exchange: the deposit in a step stays a set fraction of the owning cell's gas mass and of its thermal energy above the floor. This is numerics. It is registered like the other spray changes, and A0's md5 holds if the bound never binds in S3 and S8(c).
       - (b) The start: the chamber's initial fill, or the LOX lag behind the H2. This is a choice in the Definition, so it is the Director's.
       - Recommended: (a). It keeps the Definition, and the physics it enforces is real: in one step a drop cannot draw more heat than its cell holds. Its cost during the start-up on level 3 has to be measured.
  2. **The gas energy budget (the face inlet and the wall functions, inquiry2's lane).**
     - R3(a) needs every budget within 1e-11. At 2 levels, budget (a) is 2.45e-7 at 0.02 ms, all of it energy: 1.08e-4 J against a scale of 441 J. Mass is −3.5e-16 and momentum −1.9e-16.
     - With no spray and no AMR, the single-level gas step on this Definition (F3's run-mode line) already ends at −1.011e-8 of the initial energy (`Flow::measurements`, `energyBalanceError`). About 30 such steps give the order of 2.45e-7 (inferred).
     - So the source is the gas path's energy accounting on this Definition, not the spray or AMR. The candidates are the open-area supplies, the plate's wall-function share (`lawFlux` scaled by `plateShare_`) and the wall-function heat. Inferred, not isolated.
     - Step D's and the face inlet's checks did not score the energy budget on this Definition.
- **The Director's decision on the spray overload (23 October 2026, 00:35 by the Director's message; recorded about 00:33 CDT by the Mac's clock), and two registrations made before their code and before any RL10 output.** The only RL10 runs so far are the development probes to 0.02 ms (protocol §6), with no settled or performance value.
  - **The decision:**
    - (1) The LOX must not enter the near-empty chamber. We predict steady performance, and the 1,000 Pa start without propellant is not the engine's state when its LOX arrives. Use a gas-only converged field (or the H2-lead pre-fill), then introduce the LOX.
    - (2) Add a numerical bound that protects against the same failure in transient pockets. It may only be a time-step restriction or a local subcycling of the exchange, each sub-step conserving exactly, never a clamp on the exchanged amount.
  - **Found when applying (1):** the run definition names neither a gas-only field nor a pre-fill for the central case. Item 7's warm start is for the deviations.
    - The plan's step 1 (Ben's rule) starts the reported case "valves open, ignite", from ambient gas at rest, with a declared valve schedule. A partly developed initial state is "never the reported RL10 case".
    - So a gas-only field imported as the initial state would break that rule. An H2 lead inside the same march does not: the H2 valves open first, and the LOX opens once the H2-only flow has settled. That is the real start's order, and it is a declared valve schedule (valve travel is an input).
    - (1) is therefore done as the H2 lead within the march. If the Director meant an imported field, this is the point to say so.

- **The start amendment (item 1), registered before its code (23 October 2026, about 00:33 CDT).** Item 1 was frozen on 22 October. This changes it by the Director's decision, before any RL10 output.
  - **The ambient fill is unchanged:** 1,000 Pa, 300 K, at rest.
  - **The H2 lead:** every fuel supply (the H2 slots and the porous face) opens at t = 0 on its 0.5 ms ramp, as now. The LOX injectors stay closed.
  - **The LOX opens at t_L = 20 ms** on its 0.5 ms ramp. The igniter fires at t_L + 0.2 ms, 7 kJ over 1 ms as frozen.
    - 20 ms is *derived* from estimates: H2 at about 2.8 kg/s through the throat gives about 3 bar and a residence time of about 2 ms. An approach with a time constant near 3 ms drifts under 0.1% per ms after about 17 ms.
  - **The lead's gate:**
    - At t_L, the means of the outlet flow, p0 and the face pressure over [t_L − 2, t_L − 1] ms and [t_L − 1, t_L] ms must each differ by under 0.1%.
    - The run prints them. `run central` stops at t_L if any fails. `run development` prints the result and goes on.
    - If the gate fails, t_L is amended with a dated note before the central run, and the gate's result is not used to tune anything else.
  - **Every later time counts from t_L:** the settle window is t_L + 40 to t_L + 60 ms, and the continuation at twice the parcel rate runs to t_L + 85 ms with the move test. The end is t_L + 85 ms.
  - **In development,** `CRUCIBLE_T1_LEAD` sets t_L (default 20 ms). R3's probes may shorten it, and each states its lead.
  - **Cost:** 20 ms more simulated time, about a quarter of the central run's (*derived*). Pure-H2 cells lie outside Table B's f_H range (0.133 to 0.182), so their updates take the fallback. In the 0.02 ms probe every chemistry update took it (measured). The lead's cost is measured in development before the central launch.
  - **The deviations** start warm from the central settled state (item 7), so the lead does not touch them. Member 4's cold-start check uses the same lead.

- **Registration B: the spray-gas exchange bound, a time-step restriction (23 October 2026, about 00:33 CDT, before its code).**
  - **The defect:**
    - The spray moves its parcels against the gas as it stood at the start of the step, then deposits the step's exchange (`Flow::exchange`).
    - When the drops can draw more heat in one step than the cell's gas holds, or add more vapour than the cell holds gas, the explicit deposit overshoots. The 3-level probe's cell took three times its gas mass and went below the energy floor.
  - **The bound:** for each owning cell c, from the parcels in it at the step's start and their `Spray::rates` against the cell's gas:
    - the gas's thermal relaxation time toward the drops, τ_T = m_c c_v,c / Σ N_p G_p, with G_p the per-drop conductance `Rates::conductance` (heat flow per kelvin of T_gas − T_surface) and N_p the drops per parcel;
    - the vapour time, τ_m = m_c / Σ N_p ṁ_p.
    - On AMR the cell's level l steps at the base step over 2^l, so the base step's bound is the minimum over cells of 2^l f min(τ_T, τ_m), with f = 0.25 (*declared*, engineering). Without AMR, l = 0.
    - The physics it holds to: in a step the gas gives the drops at most a quarter of its heat above the drop surface's temperature, and gains at most a quarter of its mass as vapour, so the frozen-gas move stays close to the coupled one.
  - **Applied as a step limit only:** `Spray::exchangeDt()` (const, no side effect) returns the bound, and the driver passes it with its other limits to `Amr::step`'s or `Flow::step`'s `maxDt`.
    - The deposit, the parcels' sub-steps and the ledgers are untouched. Each step conserves exactly as now; only the step can be shorter.
    - Parcels released at the base step's end, and parcels that cross into another cell during the step, are not in that step's bound. f's margin covers them (*inferred*). The count of base steps the bound shortened is printed.
  - **Criteria:**
    - (B0) A0's md5 is unchanged (the engine change is additive).
    - (B1) The bound never binds on S3, S8(c), A1 and A2: in each, every step's dt is below the bound. The smallest ratio of bound to dt is printed.
    - (B2) Applied in A2's driver (4 threads, to 1 ms), the gas state, partial densities and parcels are bit-identical to the same run without it, in the same process.
    - (B3) The tier-1 3-level probe to 0.02 ms runs through, with t_L = 0 so the LOX meets the thin chamber. Budget (b) stays within 1e-11, and the count of steps the bound shortened is printed. Budget (a) is inquiry2's energy finding, so (a) and (c) are reported, not judged.
    - (B4) The lead's own probe (the start amendment) with t_L = 20 ms: its gas-only cost is measured, and the bound is printed at the LOX's opening.
  - **Wording corrected before any B run (about 00:38 CDT):** the amendment was first named "A1", which clashes with the spray's A1, so it is now "the start amendment". τ_T's heat capacity is c_v (the cell's gas holds its volume in the step), not c_p. B2 first read "bit-identical to 543d50b's". 543d50b's results are not kept as files, so B2 compares in-process against the same run without the bound; B0's md5 covers the engine's other paths. A3 never takes the bound, so it is not part of B2.
  - Local subcycling of the exchange is the alternative, if the restriction proves too costly in the start. It is not built now.

- **Registration B scored (23 October 2026, about 00:47 CDT; the Mac; logs in `run_mode_2026-10-23/`). B0 to B2 met, B3 not met, B4 running.**
  - **(B0) met:** A0's md5 is e95251ab17e43a961a392dcd4bd0912f, unchanged.
  - **(B1) met.** The smallest ratio of bound to dt:
    - S3 143.1 and S8(c) 143.4 (`crucible_spray_tests`);
    - A1 143.1 and 143.4 (its two chambers, AMR and uniform runs together);
    - A2 59.75, on 1 thread and on 4.
  - **(B2) met:** with the bound applied, A2's state, partial densities and parcels are bit-identical to the same run without it.
  - **R1(a) still holds:** the dry run is byte-identical. `spray_verification` fails S0(b) at 1.245e3 J/kg, the same value as on 22 October (`docs/evidence/spray_amr/a0_spray_ctest_after_2026-10-22.log`). It is not this change.
  - **(B3) not met.** With t_L = 0 the 3-level probe stops as before, in `Flow::exchange` before the first 10 μs sample (`b3_probe_3level_lead0_0.02ms.log`).
- **B3's cause, measured (about 00:45 CDT, a temporary instrumented build, reverted, not committed): the LOX flashes. A step limit cannot bound a flash.**
  - The LOX enters at 101.17 K. O2's saturation pressure there is about 2.8 bar (*derived* from NIST's 2.54 bar at 100 K); the chamber starts at 1,000 Pa.
  - On its first move each new parcel flashes the enthalpy above saturation at once (`Spray::move`, *A drop above saturation*).
  - In each of the 5 level-3 sub-steps with a flash, the vapour deposited equals the flashed mass to 4 or 5 digits (1.1e-7 kg in the last). The largest rate-driven vapour in any of the other 120 sub-steps is 1.4e-13 kg (`b3_flash_diagnostic.txt`).
  - **Not a rate.** A parcel's flashed mass is a fixed share of its mass, the injector's flow over the parcel rate. A shorter step or subcycling the exchange cannot reduce it.
  - **B is kept as built.** It never binds on the A-series, and it guards the rate-driven exchange.
- **What this means for the start amendment (*derived*):**
  - At t_L the LOX meets the lead's H2. Where that gas is above about 2.8 bar, nothing flashes. The lead's pressure was estimated at about 3 bar, a thin margin.
  - Below it the flash is small per parcel: about 1% at 2.5 bar. But a parcel at full flow carries 3.5e-6 kg, so 1% is still near the gas a level-3 face cell holds (about 5e-8 kg).
  - B4 (running) prints the gate's face pressure at t_L and the count of flashes.
- **For the Director, a design question: the flash is an instantaneous point source.**
  - (a) Rely on the start amendment, if B4 shows no flash after t_L.
  - (b) A finite flash rate: a physics model, registered like the drop model.
  - (c) A flash that stops when the cell's own vapour brings it to saturation at the drop's temperature (thermodynamic, but it acts like a limit on the exchange).
  - (d) Lighter parcels while the LOX ramps (the parcel rate is numerics).
  - Recommended: (a) first, since B4 measures whether the reported case flashes at all; then (b) if it does.

- **The Director's flash decision, registered before any of B4's output past its header is read (00:55 by the message; about 00:49 CDT by the clock).** The start amendment and registration B are accepted as built.
  - The flash is real physics, but only in the start transient, which is discarded. The prediction is steady running near 30 bar, where LOX at 101.17 K does not flash.
  - **The branch is chosen by B4's face pressure at the LOX's opening.** That is fixed here as the gate's second mean: the face pressure over [t_L − 1, t_L] ms.
    - **At least 3.1 bar** (10% above the 2.8 bar saturation, for the thin margin and local dips): **(a)**, rely on the lead, with the flash check below.
    - **Below 3.1 bar: (b),** a sourced, superheat-dependent finite flash-rate model, its source and validity range cited, built under development mode. Not (c) or (d), which are numerical stand-ins.
  - **The flash check (built now; both branches report it):**
    - The spray counts the mass flashed over the whole run (`Spray::Ledger::flashedMass`, the drops' mass times the flashed share), and the driver prints it with the flash count.
    - From the sample at t_L + 1 ms on, any further flash stops `run central` (exit 5, with a dump). `run development` prints it once and goes on.
  - **The flash check's criteria:**
    - (F1) A0's md5 is unchanged.
    - (F2) In `crucible_spray_tests`: a parcel of LOX at 95 K in O2 at 1 atm (saturation 90.2 K) flashes on its first step. The ledger's flashed mass equals the liquid mass before the step times (h(95 K) − h(T_sat)) / L(T_sat), from the liquid's own tables, within 1e-12 relative. The flash count is 1.
    - (F3) The dry run is byte-identical.

- **The flash check scored (23 October 2026, about 00:51 CDT; the Mac): F1 to F3 met.**
  - (F1) A0's md5 is e95251ab17e43a961a392dcd4bd0912f, unchanged.
  - (F2) T_sat is 90.1879 K at 101,325 Pa and the flashed share is 0.038699. The flashed mass is 2.829016e-11 kg, within 1.1e-16 relative of the expected value, from 1 flash.
  - (F3) The dry run is byte-identical.
  - `spray_verification` still fails only S0(b), at 1.245e3 J/kg, as on 22 October.
  - B4 restarted at 00:51 CDT on this build (the first launch, at 00:47, ran the build before the flash check). It has not been read past its header.
  - **Seen in the first launch's log, its last line before the stop:** a Cantera warning that `equilibrate` met 188.7 K, outside h2o2_highT's 200 to 6,000 K. Inferred: p0's equilibrium stagnation in the cold H2 (supplied at 199.67 K). The lead's gate uses p0, so if p0 is undefined or wrong there, the gate fails for that reason and not for an unsteady lead.

- **R1(b): met (23 October 2026, about 00:24 CDT).** Every line of `tier1_dry_run_length_2026-10-21.log`'s header and of its 149 × 16 section appears unchanged in the run mode's output. The run mode only adds its F3 lines, G7, the gates line and the run's own lines.

## What it waits on

- **Budget (the Director's decision, 22 October 2026, 22:40 CDT):** go for tonight's work, within Ben's raised cap. The spray-on-AMR coupling is now tier 1's critical path, and tonight's main item (below and `docs/evidence/SPRAY.md`, *Spray on AMR*). The driver's run mode follows it.

- **The face inlet (a) and the face's wall temperature,** inquiry2's lane (the handoff's spec, 21 October, about 02:40), after S7 and AMR step C. The driver then sets the switch, adds the LOX ports and the face temperature, and reruns G5 and G6 as the spec's gate F3.
- **AMR step C's remaining checks** (2, 4 and 5; the code is on `main` since af91015, and checks 1 and 3 pass by AMR.md and abdda6d's message) and **step D** (the thresholds on C1, and the y+ target).
- **AMR carrying the wall functions and the SST property corrections across levels.** `Amr` refuses both today (`core/amr.cpp:25-26`), and tier 1 needs both. This is inquiry2's AMR line, before any tier-1 run.
- **Member 19's adiabatic face** needs the wall functions' adiabatic law (`core/flow.cpp:107-108`). Only member 19 needs it.
- **The two protocol §4 amendments** (choice 1's face and decision 7's length scale), inquiry2's line.
- **The budget's two measurements:** the spray's cost on Mascotte (S7) and the gas-side cost on AMR (criterion 6). Then the central cold start, then the restated budget with the Director's go, before the deviations.
