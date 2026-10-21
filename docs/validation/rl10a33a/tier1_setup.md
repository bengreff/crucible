# RL10A-3-3A tier 1: the central case's setup and its dry-run build

Status: stated 21 October 2026, about 02:15 CDT, before its code and before any tier-1 run. It covers inputs and geometry only. It holds no predicted number and reads no measured value. It is governed by the skeleton (`tier1_skeleton.md`: choices 1 to 3, the element ranges, the cup's criterion), the plan's decisions of 21 October (`tier1_ensemble_plan.md`), and the protocol's §2 to §4. The run itself waits (the Director, 21 October): S7/Mascotte is still held behind inquiry2's corner fix.

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
