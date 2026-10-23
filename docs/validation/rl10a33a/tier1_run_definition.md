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
- **Two points found while writing the code, recorded before any run (23 October 2026, about 00:47 CDT):**
  - **R3(b) amended.** p0 is a mass-flux average over its column, and c* divides by the outlet flow. At the start the gas is at rest, so both flows are zero and p0 and c* are undefined there, not wrong.
    - They are written as NaN while their column's net forward flux, or the outlet flow, is not positive. Those samples are counted and printed.
    - R3(b) holds for every other quantity at every sample, and for p0 and c* wherever they are defined.
  - **`Spray::setParcelRate` takes the rate per second,** not a factor. The driver passes twice its rate at 60 ms.

## What it waits on

- **Budget (the Director's decision, 22 October 2026, 22:40 CDT):** go for tonight's work, within Ben's raised cap. The spray-on-AMR coupling is now tier 1's critical path, and tonight's main item (below and `docs/evidence/SPRAY.md`, *Spray on AMR*). The driver's run mode follows it.

- **The face inlet (a) and the face's wall temperature,** inquiry2's lane (the handoff's spec, 21 October, about 02:40), after S7 and AMR step C. The driver then sets the switch, adds the LOX ports and the face temperature, and reruns G5 and G6 as the spec's gate F3.
- **AMR step C's remaining checks** (2, 4 and 5; the code is on `main` since af91015, and checks 1 and 3 pass by AMR.md and abdda6d's message) and **step D** (the thresholds on C1, and the y+ target).
- **AMR carrying the wall functions and the SST property corrections across levels.** `Amr` refuses both today (`core/amr.cpp:25-26`), and tier 1 needs both. This is inquiry2's AMR line, before any tier-1 run.
- **Member 19's adiabatic face** needs the wall functions' adiabatic law (`core/flow.cpp:107-108`). Only member 19 needs it.
- **The two protocol §4 amendments** (choice 1's face and decision 7's length scale), inquiry2's line.
- **The budget's two measurements:** the spray's cost on Mascotte (S7) and the gas-side cost on AMR (criterion 6). Then the central cold start, then the restated budget with the Director's go, before the deviations.
