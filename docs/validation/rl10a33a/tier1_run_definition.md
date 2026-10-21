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

## What it waits on

- **The face inlet (a) and the face's wall temperature,** inquiry2's lane (the handoff's spec, 21 October, about 02:40), after S7 and AMR step C. The driver then sets the switch, adds the LOX ports and the face temperature, and reruns G5 and G6 as the spec's gate F3.
- **AMR step C's remaining checks** (2, 4 and 5; the code is on `main` since af91015, and checks 1 and 3 pass by AMR.md and abdda6d's message) and **step D** (the thresholds on C1, and the y+ target).
- **AMR carrying the wall functions and the SST property corrections across levels.** `Amr` refuses both today (`core/amr.cpp:25-26`), and tier 1 needs both. This is inquiry2's AMR line, before any tier-1 run.
- **Member 19's adiabatic face** needs the wall functions' adiabatic law (`core/flow.cpp:107-108`). Only member 19 needs it.
- **The two protocol §4 amendments** (choice 1's face and decision 7's length scale), inquiry2's line.
- **The budget's two measurements:** the spray's cost on Mascotte (S7) and the gas-side cost on AMR (criterion 6). Then the central cold start, then the restated budget with the Director's go, before the deviations.
