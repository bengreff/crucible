# RL10A-3-3A tier 1: the ensemble, its pre-registration and its run budget (plan)

Status: plan written 20 October 2026, about 23:45 CDT, before any tier-1 run and before any measured performance of this engine was read. It holds no predicted number. It answers Ben's note: "engine predictions are stated as a central value plus a range, the range from rerunning across the declared input ranges (the injector dimensions etc.) and the model options. Plan the RL10 pre-registration and the run budget for that ensemble." Items marked **Decide** are for the Director; the rest is engineering and stands unless overruled. The five were decided on 21 October (*Decided*, below), recorded before any tier-1 run.

## The rule (unchanged)

The protocol (`../RL10A33A_PROTOCOL.md` §4) already fixes how a range is formed. Each member is varied alone about the central case. A quantity's judged band is the central value plus and minus the sum of the absolute one-at-a-time deviations on each side (an envelope, not a standard deviation), and the quadrature sum is printed beside it, not judged. This plan treats a model option exactly like an input bracket: one more member whose deviation joins the envelope. Nothing here changes the rule.

The predicted quantities are the protocol's §3 list: F_vac, Isp_vac (on the chamber flow and on the engine-inlet flow), the nozzle stagnation and injector-face pressures, c*, and the hot-gas wall heat flow.

## The members

| # | Member | Central | Other values (runs) | Source |
|---|---|---|---|---|
| 0 | The central case | every central value below | (1, a cold start) | protocol §3, skeleton |
| 1 | Throat shape | PCHIP | sharp vertex (1) | protocol §4 |
| 2 | Injector face station | −13.1 in | −13.4, −12.8 in (2) | protocol §4 |
| 3 | Exit station | 46.0 in | 45.29, 46.40 in (2) | protocol §4 |
| 4 | Hot-gas wall temperature | 567 K | 300, 834 K (2); the chamber wall only, the face is member 19 (21 Oct, about 02:35) | protocol §4 |
| 5 | LOX SMD | SPRAY.md's central: Burick at the ring's group 5.03, SMD 447 μm (21 Oct, about 02:28; `tier1_setup.md`) | its bracket's ends (2): 263 and 581 μm, the scatter's stated width | protocol §4, SPRAY.md |
| 6 | Injection temperatures | plenum values | ±0.1 R (0: carried from tier 0's sensitivity; decided 21 Oct) | protocol §4 |
| 7 | Supply turbulence intensity | 5% | 2%, 10% (2) | protocol §4 |
| 8 | LOX post inner diameter d_L | 3.54 mm | 2.62, 4.78 mm (2) | skeleton |
| 9 | Post tip wall t | 0.272 d_L | 0.20, 0.37 d_L (2) | skeleton |
| 10 | H2 annulus gap h | 0.268 d_L | 0.09, 0.80 d_L (2) | skeleton |
| 11 | Post recess | 0.5 d_L, as geometry: the LOX slot's exit set back inside a cup (decided 21 Oct) | 0, 1.0 d_L, as geometry (2); not run: the cup is not resolved and Burick's curve has no recess term (21 Oct, about 02:08; `tier1_setup.md`) | skeleton |
| 11b | How the recess enters | the cup in the geometry | through the SMD correlation only, no cup (1; decided 21 Oct); now the central case's treatment, not run (21 Oct, about 02:08) | skeleton |
| 12 | Outer ring radius r_8 | R_c − p/2 | R_c − p (1) | skeleton |
| 13 | Inner ring radius r_1 | r_8 − 7p | the inner ring's clearance (1) | skeleton |
| 14 | Transpired H2 fraction | 10% | 5%, 15% (2; the bracket is inferred, accepted 21 Oct) | skeleton choice 3 |
| 15 | Chemistry | Table B (decided 21 Oct; it passed its gates under the split fallback bar, 21 Oct about 01:48) | Table A to 6,000 K, shifting equilibrium (1) | skeleton, TABLE_B.md |
| 16 | SST property corrections | Hasan-Pecnik | none (1; decided 21 Oct) | WALL_FUNCTIONS.md |
| 17 | Grid | the reported AMR grid | its finest level halved (1); the Richardson estimate is the member | protocol §4 |
| 18 | Supply turbulence length scale at the slots | 0.07 D_h (inferred; one pitch on the porous face, unchanged) | 0.035, 0.14 D_h (2; guessed, ×/÷ 2; decided 21 Oct, about 02:35) | the Director's decision 7, below |
| 19 | The injector face's wall | isothermal at the transpired H2's inlet temperature, 199.67 K (inferred) | adiabatic (1; decided 21 Oct, about 02:35) | the Director's decision 8, below |

Count (derived): the central cold start and 27 deviations, one of them the grid run (with member 6 carried from tier 0), plus one cold-start check (below). Member 11b was added on 21 October with the Director's decision 5; it was 26 before. *Dated note, 21 October 2026, about 02:35 CDT:* members 11 and 11b are not run (decision 5 applied, about 02:08), and members 18 and 19 add three runs (decisions 7 and 8), so the count is again 27 deviations: 26 warm and the grid run. The wall functions' measured errors join the heat-flow band without a run (protocol §4). The skeleton's three modelling choices (annular rings, no swirl, a distributed transpiration) stay declared limits with stated directions, not members, since tier 1 cannot run their alternatives.

## How a member runs (engineering, proposed)

- **The central case** starts "valves open, ignite" and runs until settled: F_vac, injector-face pressure and outlet mass flow each drift by under 0.1% over the last 1 ms (protocol §3).
- **A deviation starts warm** from the central case's settled state, with the member's change applied, and runs until settled by the same rule. The settled state does not depend on how it was reached if it is unique; that is checked once. Member 4's 834 K end also runs from a cold start, and its two settled values must agree within 0.1%. If they do not, every member runs cold and the budget below doubles. Warm starts with this one cold-start check were accepted on 21 October.
- **Members that change the geometry (2, 3, 11, 11b, 12, 13)** are meshed again. They start warm from the central case's state mapped onto the new mesh, if that mapping exists when they run; otherwise they start cold, and the budget counts them as cold runs.
- **Exit-station members (3)** change only the divergent end. They still run the whole domain, since the wall layer is subsonic there.
- **A member that fails** (does not settle, crashes, or leaves the tables' domain) is printed as "not measured" with its cause. The band is printed without it and flagged incomplete. Fixes before the prediction commit are allowed, with dated notes; any fix is applied to every member already run.
- **P2** stays at tier 0 unless the ensemble's measured cost leaves room.

## Pre-registration, in order

1. This plan committed: the members and the rules. It may be amended with dated notes until step 3.
2. The tier-1 settings and every central input fixed in `tier1_skeleton.md` (its table), and the open items above decided.
3. The central case and every member run. `RL10A33A_PREDICTION.md` gains tier 1: per quantity, the central value, each member's signed deviation, the band, the quadrature sum, and every run's log path. That commit's hash is recorded.
4. Only then does anyone read a measured value (protocol §5).

The central case's own output is not measured data. Changing the model after seeing it is allowed before step 3, but each change is dated and reruns every member already run.

## Budget

Not yet measured for tier 1's setup (liquid LOX, AMR, Table A, SST with the corrections, wall functions). The basis:

| Row | Value | Label |
|---|---|---|
| Base mesh 150 x 16, pre-vaporized, Table A, wall functions, 4 threads: cold start to settled thrust | 3.3 min (uniform rings), 5.9 min (stretched) | measured, RL10_DEV.md |
| The same on AMR to y+ 956 at the throat, with the throat tag | 28 to 35 min | derived, RL10_DEV.md |
| The same on AMR to y+ 300 | 1.7 to 2.2 h | derived, RL10_DEV.md |
| The spray's cost on the whole step | 1.5 to 4 times the gas step's cost | guessed from S6's 3.7 times on S3 (SPRAY.md); parcels per cell on RL10 not known |
| A warm-started member against a cold start | 0.5 | guessed (no start transient, a shorter settle) |
| The grid member against the central case | 6 to 8 times | derived: the finest level's cells 4 times, its steps 2 times; it dominates the cost on the wall strip |

Per run (derived from the rows):
- **The central case, and the cold-start check:** 0.7 h (28 min at y+ 956, times 1.5 for the spray) to 8.8 h (2.2 h at y+ 300, times 4) each, on 4 threads of the Mac.
- **A warm member:** 0.35 to 4.4 h.
- **The grid member:** 4.2 to 70 h.

Total for the central case, the check, 26 warm members and the grid run (derived; 25 before member 11b): **15 to 200 run-hours on 4 threads, or up to 240 if the 9 re-meshed runs (members 2, 3, 11, 11b, 12, 13) start cold.** Concurrency: 2 runs at a time on the Mac (10 cores) and 5 on backhouse (28 threads, when it is back online), so about 7 at once. The grid member is the longest single run, so it sets the wall time: **about 5 hours at the low end** (one night), and **about 3 days at the high end**, with the other members done in about 18 to 28 hours beside it. Giving the grid member more threads on backhouse would shorten it; its scaling there is not measured.

*Dated note, 21 October 2026, about 02:35 CDT (derived from the rows above):* the count is still 26 warm members and the grid run, so the total stays **15 to 200 run-hours**. The re-meshed runs are now 6 (members 2, 3, 12, 13), so the cold-start ceiling falls from about 240 to **about 230 run-hours**. Members 18 and 19 cost 0.35 to 4.4 h each per run, 1.1 to 13 run-hours together.

What sets the range is the wall-function y+ target (300 against 956: a factor of 4) and the spray's cost (a factor of about 3). Both are measured before any member runs:
- AMR criterion 6 and the wall functions' y+ question (AMR step D) give the gas-side cost.
- The spray's ratio on the Mascotte setup (SPRAY.md's restated S6 trigger) gives the spray side.
- The central case's own cold start then replaces this table. The budget is restated, with the Director's go, before the 27 deviations start.

Token cost is small against machine time: each member is one slot job with a short log read. The overnight session's share (about 20% of the week's usage) is not the binding constraint; the machines are.

## Decided (21 October 2026, the Director's brief on fa51719)

Recorded before any tier-1 run and before any measured performance is read.
1. **Chemistry.** If Table B passes its gates before the prediction commit, the central case uses Table B: it carries finite-rate chemistry into the nozzle, the highest fidelity of the two. Table A (shifting equilibrium) is then a member. Table B's gates are those of `docs/evidence/TABLE_B.md`, with the Director's bar of 21 October: c*, Isp_vac and F_vac within 0.2% of the CVODES run on C1, and a steady fallback under 1% of updates in each region. If Table B has not passed by then, Table A is central, and Table B is not a member.
   - SST with the property corrections against without them is a member (16).
   - *Dated note, 21 October, about 00:45 CDT:* Table B passed these gates on C1 with the Δη 0.1 table, the corner drop and the induction cut (`TABLE_B.md`): c* −0.06%, Isp_vac −0.13%, F_vac −0.14%, 0% fallback in each region. It has not yet been built or checked on the RL10's own range, so this is not yet a pass for the central case.
   - *Dated note, 21 October, about 01:50 CDT:* on the RL10's range (`TABLE_B.md`, *The C1-like check at RL10 conditions*), the 0.2% bars pass on the development chamber at O/F 5.055 against CVODES, both settled: c* −0.081%, Isp_vac −0.138%, F_vac −0.139%. The chamber's 1% fallback bar fails: 26.8% of its updates, all from the face column's unburnt cells, which are off the manifold by cause. The convergent and divergent are at 0%. The bar's reading is with the Director. Until it is decided, Table A stays central.
   - *Dated note, 21 October, about 01:48 CDT (the Director's decision (b) on 30483ac):* the fallback bar is split. Defect fallbacks (refused, flagged or frozen corners, integrator failures, states inside the axes the table should hold) stay barred at 1% per region. Off-manifold fallbacks by cause are reported per region and not barred, since CVODES is exact for them and the 0.2% bars include them. Each check also reports its wall time against an all-CVODES run. Re-judged from the logs (`TABLE_B.md`, *The fallback bar split*): C1 passes (0 defect and 0 by cause in every region; 76 s against 937 s). The RL10 check passes (0 defect in every region; by cause, chamber 26.8%, the face column and two flame-column cells; 290 s against 1,973 s). **So Table B is the central case's chemistry, and Table A (shifting equilibrium) is member 15's other value.** One named check follows: the queued replay of the CVODES run's samples tests the curve at the two flame-column cells (1.8% of the chamber's updates). If it shows a defect there, this is reopened.
2. **The transpired fraction's bracket,** 5% to 15% about 10%: accepted.
3. **Warm starts** for the deviations, with the one cold-start check on member 4's 834 K end: accepted.
4. **Member 6** is carried from tier 0's sensitivity, with no tier-1 runs: accepted.
5. **The recess** enters the central case geometrically, as the LOX slot's exit set back inside a cup, if the mesh resolves it. The treatment through the SMD correlation alone (no cup) is a deviation member (11b).
   - The resolution criterion for the cup (cells across its depth and width) is stated in `tier1_skeleton.md` before any run.
   - **Applied 21 October 2026, about 02:08 CDT, before any tier-1 run** (`tier1_setup.md`, *Decision 5 applied*): the dry run's G7 fails (1.41 cells across the depth at the finest planned level, and `Mesh` has no cells behind the face plane), so the branch below holds. Burick's curve has no recess term, so member 11 would repeat the central case and is not run either. The recess's whole effect is a declared limit: c* and Isp biased low.
   - If the tier-1 mesh does not resolve the cup, the central case takes the SMD correlation alone, member 11b is not run, and the missing cup is a declared limit with its stated direction. That is decided and dated before any run.

The budget's basis rows are on Table A. Table B's cost per cell update on RL10 is measured with the central case's cold start, when the budget is restated.

## Decided (21 October 2026, about 02:35 CDT, the Director's brief on f29c38c)

Recorded before any tier-1 run and before any measured performance is read. They answer `tier1_setup.md`'s three findings, and they close the face-inlet decision that was open.

6. **The face inlet is (a): an open-area sub-cell supply in `core/flow`.**
   - Each supply declares its true annulus. Each face ring it overlaps takes the supply's mass, momentum, energy and turbulence flux in proportion to the overlap, at the stream's own velocity on its true area.
   - So each slot keeps its own mass flow, momentum and J at any mesh, and no cells are added.
   - The rest of each ring's face is the post tips (wall, decision 8) and the LOX slot, where the liquid enters as parcels.
   - It is inquiry2's lane (`core/flow.*`, `core/transport.cpp`, `core/amr.cpp`). The interface spec is in `docs/SESSION_HANDOFF.md` (21 October, about 02:35), for inquiry2 to code after S7 and AMR step C. The central case does not run before it lands.
   - **Declared limit:** a ring's cell mixes its slots' streams within one cell, so the shear layers between the H2 slots, the tips and the LOX are not resolved at any affordable mesh (the setup's G5). Its direction is not inferred here. Mixing within the cell is immediate, which speeds it. But the H2's momentum is spread over the cell's height, so the gas the drops meet is slower, which slows their acceleration. The two oppose. The size falls with the face cell's height.
7. **The supply turbulence's length scale is 0.07 D_h at the slots.**
   - D_h is the annulus's hydraulic diameter, twice the slot's radial width. The value is the fully developed duct's standard estimate: Nikuradse's pipe mixing length at the centre, 0.14 R = 0.07 D. *Inferred* to hold at the slots' exit; there are no element-level data.
   - The porous face keeps one pitch (protocol §4). Its inflow is 2.0 m/s, so its k is small either way.
   - **Measured from the driver** at the nominal face state (`tier1_dry_run_length_2026-10-21.log`), the viscosity ratios are 137 to 2,520 at the slots and 524 on the porous face, against 46,000 at the slots with a pitch.
   - **Wired into `crucible_rl10_tier1`** (`CRUCIBLE_T1_SLOT_LENGTH`, central 0.07). The dry run's output is unchanged apart from that line (measured: diff against `tier1_dry_run_g8_2026-10-21.log`).
   - **The bracket is member 18: 0.035 and 0.14 D_h,** ×/÷ 2 (*guessed*). The area mean of Nikuradse's profile, 0.04 D (*derived*), lies inside it.
     - It is cheap: two warm runs that change one Definition value and no geometry, 0.7 to 8.8 run-hours (*derived* from the budget's rows), 4 to 5% of the total.
     - Both ends are dry-run (measured): slot ratios 68.5 to 1,260, and 274 to 5,039.
   - Protocol §4's row (one pitch) needs the matching amendment, beside choice 1's (inquiry2's line, `tier1_skeleton.md`).
8. **The injector face's wall is isothermal at the transpired H2's inlet temperature.**
   - The porous Rigimesh face is transpiration-cooled, so its wall is set to the fuel injector plenum's 199.67 K (359.4 R; `tm107318_inputs_addendum.md`, protocol §2). This is *inferred*: the transpired H2 holds the face near its own temperature.
   - Under decision 6 the face's wall share is the post tips. They lie between the LOX at 101 K and the H2 at 200 K, and are taken at the same value (inferred; each tip is 0.96 mm wide at the central t).
   - The chamber wall stays at 567 K. Member 4 moves the chamber wall only.
   - **The adiabatic face is member 19** (one warm run).
   - Both need inquiry2's code (spec in the handoff):
     - `Definition` has one `wallTemperature`;
     - the wall functions refuse an adiabatic wall (`core/flow.cpp:107`: the adiabatic law is not implemented).
   - Until then the dry run's Definition carries 567 K on every wall.

## Open, for the Director (as asked on 20 October; answered above on 21 October)

1. Which model options are members: Table A against Table B (proposed: a member only if Table B passes its gates), and SST with the corrections against without (proposed: a member, since the hot-wall question was a real disagreement, 6 October).
2. The transpired fraction's bracket, 5% to 15% (inferred; the source says "about 10%").
3. Warm starts for the deviations, with the one cold-start check.
4. Member 6 carried from tier 0's sensitivity rather than two tier-1 runs (a 0.1 R change is 0.056 K).
5. How the recess enters tier 1 with continuous annular inlets: as a set-back of the LOX slot's exit inside a cup, or only through the SMD correlation.
