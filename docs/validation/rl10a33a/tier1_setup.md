# RL10A-3-3A tier 1: the central case's setup and its dry-run build

Status: stated 21 October 2026, about 02:00 CDT, before its code and before any tier-1 run. It covers inputs and geometry only. It holds no predicted number and reads no measured value. It is governed by the skeleton (`tier1_skeleton.md`: choices 1 to 3, the element ranges, the cup's criterion), the plan's decisions of 21 October (`tier1_ensemble_plan.md`), and the protocol's §2 to §4. The run itself waits (the Director, 21 October): S7/Mascotte is still held behind inquiry2's corner fix.

## What the dry run builds

`crucible_rl10_tier1` (`tests/rl10_tier1.cpp`) builds the central case and exits before any flow step. Each element dimension can be replaced by a member's value through the environment, so every member's setup is checked the same way before it runs. It reads `core/*` and changes none of it. rl10dev's contour moves to `tests/rl10_geometry.hpp`, which both drivers share; rl10dev's output is unchanged, checked by diffing a short run.

1. **The contour.**
   - rl10dev's PCHIP wall through TM-107318 Table E1's 16 nodes, the throat (0, r_t), the face node at R_c = r_t 5.13 / 2.47 and the exit node (46.0 in, r_t √61). r_t is 0.062738 m.
   - The axial spacing is uniform, so (as in rl10dev) the throat sits on a station, and the face is the station nearest −13.1 in.
2. **The face (choice 1).**
   - **The rings:** 8 at r_i = r_1 + (i − 1) p, with p = √(π R_c² / 216), r_8 = R_c − p/2 and r_1 = r_8 − 7p (central).
   - **Each ring:** the LOX slot is centred on r_i, of width A_L / (2π r_i), where A_L = 27 π d_L²/4. A post-tip wall of thickness t lies on each side of it. Beyond each tip is an H2 slot of area A_H/2, where A_H = 27 π (D_a² − D_o²)/4; the slot's far edge follows from its area.
   - **The cup:** from the inner H2 slot's inner edge to the outer H2 slot's outer edge, set back by the recess.
   - **The porous face (choice 3):** every face radius outside the cups. That includes the centre inside ring 1 (the igniter's place) and the band outside ring 8.
3. **The flows** (protocol §2, choice 3; the temperatures from protocol §2 and the addendum):
   - Oxidizer 13.9525 kg/s, and fuel 13.9525 / 5.055 = 2.76014 kg/s.
   - 10% of the fuel is transpired, at a uniform mass flux over the porous area. The rest goes through the 16 H2 slots in equal shares.
   - The LOX goes 1/8 per ring, as parcels across its slot.
   - Fuel at 199.67 K, oxidizer at 101.17 K, the wall at 567 K.
4. **The supply turbulence** (protocol §4): intensity 5%, length scale one pitch. `Supply` takes a viscosity ratio, so the length scale is converted at each stream's own state: ω = √k / (C_μ^¼ L) and μ_t/μ = ρk / (ω μ), with k = 3/2 (I u)².
5. **The mesh.** `Mesh` is built from the Definition at each candidate (nz, nr, radial stretching). Each supply's ring cells follow `Flow`'s rule: a face ring belongs to the supply that holds its centre.
6. **The build's own checks.** `Flow` is constructed from the Definition, which runs its own checks on the supplies (no overlap; each supply holds a ring centre). Then `Spray` is constructed with the 8 LOX injectors, and the program exits.

## The checks, with their bars (stated before the code)

- **G1, the contour:**
  - the throat station's radius equals r_t to 1e-12 relative;
  - the face lies within 0.03 in of −13.1 in (*guessed*: a tenth of member 2's half-bracket);
  - the face radius equals R_c and the exit radius r_t √61, each to 1e-12 relative.
- **G2, the layout:**
  - no two slots overlap, ring 1's inner H2 slot clears the axis, and ring 8's outer H2 slot clears the wall;
  - each ring's LOX area is A_L and each H2 slot's area A_H/2, to 1e-12 relative;
  - the cups' and the porous face's areas sum to π R_c², to 1e-12 relative.
- **G3, the flows:**
  - the injectors' and supplies' oxidizer sums to 13.9525 kg/s and their fuel to 2.76014 kg/s, to 1e-12 relative;
  - so the global O/F is 5.055, and each ring's is 5.055 / 0.9 = 5.6167;
  - the transpired share is 10%, and the porous mass flux is the same on every porous supply's true area, to 1e-12.
- **G4, the derived injection numbers (reported, not barred):**
  - The LOX velocity per ring, at the liquid's density at 101.17 K. The skeleton derived 5.8 m/s at d_L = 3.54 mm.
  - The H2 slot velocities and the transpired velocity, the momentum-flux ratio J, and Burick's correlating group (V_g − V_L) / (V_L MR), with MR the ring O/F.
  - These are evaluated at a nominal face pressure of 2.8 MPa: our own development runs' value, never a measured one. Tier 1 replaces it with its own face pressure.
- **G5, the inlets on the mesh:**
  - For each supply: its ring cells, and their face area against its true area. A supply's mesh velocity is its true velocity times A_true / A_mesh, and its J moves by the square of that.
  - A slot is resolved as an inlet if the two areas agree within 5% (so J within about 10%; *guessed*).
  - Reported at three meshes:
    - rl10dev's 76 × 8;
    - the AMR plan's base, 150 × 16 (`AMR.md`), with nz moved to the nearest value that passes G1;
    - that base's finest planned level, 8 times finer (3 levels), built here as a uniform mesh.
- **G6:** `Flow` accepts the layout at each mesh, or its error is reported. If it refuses the true layout, the same check is repeated on a merged layout, with each ring's two H2 slots taken as one H2 band across its cup. That band's area against the slots' areas is J's error.
- **G7, the cup:** the skeleton's criterion at the finest level, namely 4 cells across the depth and 2 across each slot and tip, and the mesh carrying the cup's walls.

## Expected (*inferred*, before the code)

- **G7 fails:** the finest planned level is about 1.25 mm axially and 1 mm radially, against the 1.77 mm depth and slots as narrow as 0.33 mm. `Mesh` is one body-fitted map of r(z) with no blocked cells, so it cannot carry a groove behind the face.
- **G6 refuses the true layout at every candidate mesh:** the outer rings' slots are narrower than a ring cell, so they hold no ring centre.
- **On the merged layout,** the outer rings' H2 velocity falls by about the cup's width over the two slots' widths (about 4 on ring 8), so J falls by about its square. That breaks choice 1's statement that the injection velocities and J are the elements' own.
- **So the expected finding:** the face needs a sub-cell inlet, which takes each slot's own mass, momentum and energy flux into the cell whose face holds it, with the rest of that face a wall. Resolving the slots would instead need about 0.1 mm rings, and so time steps about 10 times shorter.
- The inlet sits in `supplyFace` (`core/flow.cpp`, inquiry2's lane). It is proposed, not coded here.

## Results (21 October 2026, about 02:08 CDT; `tier1_dry_run_2026-10-21.log`, `tier1_dry_run_probe_2026-10-21.log`)

`crucible_rl10_tier1` builds, runs in about 1 s, and reads no measured value. rl10dev's output after the contour moved to `tests/rl10_geometry.hpp` is byte-identical to before (`fr`, 76 × 8, to 30 µs: mesh, history, the three fields, and the log apart from wall times; measured: `cmp`).

| Check | 76 × 8 | 149 × 16 | 1192 × 128 (finest planned) |
|---|---|---|---|
| G1, contour | fail: face at −13.2542 in (−0.154) | pass: −13.0862 in (+0.0138) | pass: −13.0862 in |
| G2, layout | pass | pass | pass |
| G3, flows | pass | pass | pass |
| G5, supplies within 5% of true area | 0 of 25 (19 hold no ring centre) | 0 of 25 (15 hold none) | 8 of 25 (7 hold none) |
| G6, true slots | refused | refused | refused |
| G6, merged bands | refused | refused | accepted; J at 0.048 (ring 8) to 0.34 (ring 1) of the slots' |
| G7, the cup | | | fail |

- **G1:** the throat station's radius, the face radius and the exit radius are exact to rounding (relative error 0) at every mesh. rl10dev's 76 × 8 puts the face 0.154 in from −13.1, so it is a development mesh, not tier 1's.
- **G2:**
  - the rings run from ring 1's inner H2 slot at 3.677 mm to ring 8's outer slot at 123.912 mm, against the wall at 130.302 mm;
  - the slot areas are exact to 3.7e-14;
  - the cups (12,784.04 mm²) and the porous face (40,555.84 mm²) sum to π R_c² exactly.
- **G3:** oxidizer 13.952500 kg/s, fuel 2.760138 kg/s, O/F 5.055000. Transpired: 10.0000% at 6.805773 kg/(m² s) on every porous band. Ring O/F 5.61667.
- **G4** (nominal 2.8 MPa; the same on every ring by construction):
  - LOX density 1,084.65 kg/m³ at 101.17 K (`tools/o2_saturation_nist.tsv`). The skeleton's "about 1,140" was rounder, so V_L is 6.05 m/s, not 5.8.
  - H2 density 3.400 kg/m³, so V_H2 is 176.9 m/s, J is 2.68, and Burick's group (V_g − V_L)/(V_L MR) is 5.03, inside Fig. 38's 0.4 to 12.
  - The transpired H2 enters at 2.00 m/s.
  - **The supply turbulence at a pitch length scale gives viscosity ratios of 46,000 at the slots and 524 on the porous face.** Protocol §4's length scale of one pitch (15.7 mm) was written for the uniform face. At choice 1's slots, which are 0.33 to 6.1 mm wide, it is up to 47 times the slot's width.
- **G5 and G6, the probe at 149 × 1280, 2560 and 5120** (equal rings of 0.102, 0.051 and 0.025 mm):
  - `Flow` accepts the true slots from 1,280 rings.
  - Within 5%: 17, 21 and 24 of 25 supplies. The worst is ring 8's inner H2 slot, at 1.20, 1.10 and 1.05 times its true area.
  - So the 5% bar needs more than 5,120 equal rings (rings of 0.025 mm or less). That is about 40 times the finest planned level's radial count. The radial cell is then about 50 times thinner than its axial length, which sets a time step roughly 40 times shorter.
- **G7 at 1192 × 128:**
  - 1.41 cells across the cup's 1.77 mm depth (4 needed);
  - 0.33 across the narrowest slot and 0.95 across a tip (2 needed);
  - `Mesh` has no cells behind the face plane, so it cannot carry the cup's walls at any resolution.

**Verdict:** the setup's inputs (G2, G3) and its contour (G1, at tier-1 meshes) are as stated. The face of choice 1 cannot be represented by today's `Flow`:
- the true slots are refused below 1,280 rings, and they need more than 5,120 rings to meet 5%;
- the merged bands are accepted at the finest level but cut J to 5 to 34% of the elements'. Choice 1 states that J is the elements' own.
- The cup fails G7.
These are the outcomes expected before the code, except for the numbers.

## Decision 5 applied (21 October 2026, about 02:08 CDT, before any tier-1 run)

The cup is not resolved (G7), so the plan's decision 5 branch applies:
- the central case takes the recess through the SMD correlation alone;
- member 11b is not run;
- the missing cup is a declared limit.

**Its direction (inferred):**
- A recess up to one post diameter raises mixing (SP-8089 §3.1.1.1.4, on Burick's data). Without the cup, the near-face mixing is under-represented.
- So c* and Isp are biased **low**, and the heat flux to the face and the first part of the wall is biased **low**.

**A consequence, recorded as applied and reversible:** Burick's Fig. 38 is one curve for flush and recessed posts, with no recess term (`docs/evidence/spray/burick_fig38.md`).
- So "through the SMD correlation alone" gives the recess no input in tier 1. The central case's SMD is the same at recess 0, 0.5 and 1.0 d_L.
- So member 11 (recess at its bracket's ends) would repeat the central case, and is not run.
- The recess's whole effect is the declared limit above.
- If a recess-dependent drop-size source is found before the prediction commit, member 11 returns.

## Findings for the Director and inquiry2 (not decided here)

1. **The face inlet: a DECISION NEEDED for the Director.** Choice 1's slots cannot be carried at an affordable mesh (G5, G6). The options:
   - **(a)** A sub-cell ("open-area") supply in `core/flow` (inquiry2's lane).
     - A supply declares its true annulus.
     - Each face ring it overlaps takes the supply's mass, momentum and energy flux in proportion to the overlap, at the stream's own velocity (its true area), so the slot's mass flow and momentum and J are kept at any resolution.
     - The rest of each ring's face is wall.
     - A ring's cell mixes its slots' streams within one cell height, so the near-field shear layer is unresolved at any affordable mesh. That is a declared limit, as the cup is.
   - **(b)** The merged bands that today's `Flow` accepts, with J at 5 to 34% of the elements' as a declared limit. This weakens choice 1's purpose, since J and the injection velocities set the mixing.
   - **(c)** Resolving the slots: more than 5,120 rings and time steps roughly 40 times shorter. Not affordable.
   - **Recommendation: (a).** It keeps the momentum and J that choice 1 exists to keep, and it adds no cells.
2. **The supply turbulence's length scale (engineering, to be dated before any run).** At the slots, protocol §4's pitch gives viscosity ratios of 46,000. Proposed:
   - L = 0.07 times each slot's hydraulic diameter (twice its width), the fully developed duct's estimate, at the slots;
   - one pitch on the porous face.
   - This needs the protocol §4 amendment that the skeleton already names for inquiry2's line.
3. **One wall temperature (engineering, inquiry2's lane).**
   - `Definition::wallTemperature` is one value, so the face's post tips and the chamber wall are both at 567 K.
   - The face of a transpiration-cooled injector runs near the H2's temperature.
   - Either a per-segment wall temperature, or the face's wall segments treated as adiabatic, is needed. With (a), the wall share of the face is only the post tips and the slots' complement within the rings.

## The LOX drop size and its release (stated 21 October 2026, about 02:25 CDT, before its code)

The skeleton puts the LOX in as Rosin-Rammler parcels with SPRAY.md's Burick bracket, evaluated at the element's dimensions (member 5: the central SMD, with the bracket's two ends as the member's runs). This applies the rule that Mascotte's prediction states (`docs/validation/MASCOTTE_A10_PREDICTION.md` §3) to tier 1. Burick's data are wax/GN2 cold flow, not connected to the RL10. Nothing here reads an RL10 measurement.

- **The group.** Each ring's (V_g − V_L)/(V_L MR), where MR is the ring O/F and both velocities are at G4's nominal face pressure. For the central case it is 5.03 (G4), inside the range of Burick's data (about 0.3 to 11). So the case needs no extrapolation rule, unlike Mascotte at 68.4.
- **The central diameter.**
  - D̄/D_L is read from the digitised curve (`docs/evidence/spray/burick_fig38.md` §1.4), interpolated linearly in log-log between the table's neighbouring points. D_L is the post's inner diameter d_L.
  - The ordinate is used as printed: the report names a property correction (Ingebo) but does not give it.
  - *Inferred* direction: LOX's viscosity is far below the wax's, which tends toward smaller real drops. That is a declared limit of the bracket, as at Mascotte.
- **The exponent.** q = 3.207, solved from Burick's Eq. (5), D30 = 0.455 D̄ (Mascotte §3). Then SMD = X / Γ(1 − 1/q) with X = D̄ / (ln 2)^(1/q), so SMD = 0.851 D̄.
- **The bracket (member 5's two runs).**
  - SPRAY.md's rule: "the scatter band's low and high diameters at the case's correlating group".
  - Near group 5 the digitisation holds three markers (x 3.9 to 5.4), all within 4% of the curve. Three points cannot bound a scatter band.
  - So the band is the data's scatter about the curve as the digitisation states it for the populated range: high = curve × 1.3, low = curve / 1.7 (its "+30%/−40%").
  - Readings not taken, so the choice is visible:
    - the local markers (a band of about ±4%);
    - the flagged low point's factor of 2.4 below the curve (Mascotte also declined to anchor on it).
- **The release.**
  - The drops are released from each ring's LOX slot over Porcheron et al.'s core length, as at Mascotte: L/D = 2.85 (ρ_g/ρ_l)^−0.38 Oh^0.34 M^−0.13, with D = d_L, ρ_g the H2 at G4's state, M = J, and Oh = μ_l / √(ρ_l σ d_L) at 101.17 K from the NIST table that the spray loads.
  - Lasheras and Hopfinger's 6 M^−0.5 is reported beside it.
  - The Spray releases a Rosin-Rammler injector's drops from the slot's outer radius. A ring's slot stands for 27 round jets, so the release surface is the slot's edge, not each jet's. This is a declared approximation of choice 1; the slot is 0.35 to 3.4 mm wide.
- **Members.** A member that changes d_L, t or h changes the velocities, the group and D_L, so the driver evaluates its SMD at its own group. A member's central reading is the curve. Member 5's ends are the band at the central element.

**G8 (reported, with one bar):**
- The bar: every ring's group lies inside the digitised curve's table (0.23 to 14.5), so the reading is an interpolation.
- Reported: D̄/D_L and the SMD at the centre and the band's two ends; q, with D30/D̄ from it; and both core lengths.
- `Spray` is constructed with the central SMD, q and the Porcheron length (G6 at the merged layout).
