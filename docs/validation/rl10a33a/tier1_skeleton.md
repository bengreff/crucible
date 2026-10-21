# RL10A-3-3A tier 1: pre-registration skeleton and declared injector ranges

Status: skeleton written 20 October 2026, about 18:00 CDT, before any tier-1 run and before any measured performance of this engine was read. It holds no predicted number. It extends the protocol (`RL10A33A_PROTOCOL.md`, written 13 October in inquiry2's line and not yet on `main`). Its §4 brackets and §5 comparison rules govern; nothing here changes them, except the face (the Director's decision, *Declarations*). What it adds: tier 1's settings to be fixed, and declared ranges for the injector-element dimensions that the inputs lack (`injector_element_inputs.md`). Every range below is **inferred** from analogous coaxial LOX/H2 injectors (`analog_coax_elements.md`, a Sonnet research agent, 20 October), not from any RL10 document.

## Declarations

- **The analog agent's incidental exposures (agent only; nothing recorded).**
  - A generic "coaxial" search of NASA SP-8089's OCR text showed one survey-table row for "Centaur, RL10A-3" with performance figures and an element count. The agent recorded no value from it and excluded the table.
  - NASA CR-187037 sizes its test element against "RL-10 injector elements"; the agent discarded that source unread beyond that sentence.
  - The author of this file read only the agent's file, after a count-only keyword grep and a digit-masked view of every hit. No RL10 value is in it.
- **The face (the Director's decision, 20 October; replaces protocol §4's uniform face for tier 1).** Each injector ring is a continuous annular inlet carrying its ring's mass flow, and the transpired H2 is a distributed inflow over the porous face (*Tier 1's modelling choices* below). The element dimensions below enter tier 1 through the inlets' widths, the injection velocities, the spray's drop size and the supply turbulence. The protocol's §4 text (inquiry2's line) needs the matching amendment when it reaches `main`.

## Tier 1's settings (to be fixed in the prediction commit)

| Setting | To be fixed as | Source |
|---|---|---|
| Mechanism and thermo | `tools/mechanisms/h2o2_highT.yaml` (TM-4513 to 6,000 K, OH from Ruscic 2002, 1 bar standard state) | TABLE_B.md; CEA check worst 1.79 K |
| Chemistry | Table A to 6,000 K (equilibrium), unless Table B passes its gates before the prediction commit | TABLE_A.md *To 6,000 K* (passes); TABLE_B.md (gates fail, 20 October) |
| Turbulence and wall | SST with the Hasan-Pecnik corrections; wall functions; uniform wall temperature per §4 | WALL_FUNCTIONS.md |
| LOX | liquid parcels; SMD central and bracket per SPRAY.md (Burick), evaluated at the element dimensions below | SPRAY.md; protocol §4 |
| Grid | adaptive resolution; settled quantities at two resolutions (finest level halved), Richardson estimate joins the band | protocol §4 |
| Start and settling | "valves open, ignite"; settled when F_vac, injector-face pressure and outlet mass flow each drift < 0.1% over 1 ms | protocol §3 |
| Operating point | P1 (oxidizer 13.9525 kg/s, O/F 5.055) | protocol §2 |

## Tier 1's modelling choices (the Director, 20 October), each with its expected effect stated before any run

Every expected effect below is **inferred** (reasoning from the geometry and the physics, no run). Each is a declared limit of tier 1, not a band member, unless a sensitivity run is added before the prediction commit.

1. **Each ring is a continuous annular inlet carrying its ring's mass flow.** The per-ring element count sets the ring's flow (central: 216 / 8 = 27 per ring, so each ring carries 1/8 of the element flow) and not its geometry, so the 27-per-ring fit conflict (*Open*, below) no longer arises. Each ring is a LOX slot at the ring radius between two H2 slots, one on each side, with the post's tip wall t between them. The slots keep the elements' flow areas: LOX area N π d_L²/4, and the H2 annulus area N π (D_a² − D_o²)/4 split equally between the two H2 slots. So the injection velocities, the momentum-flux ratio and the inputs to Burick's SMD are the discrete elements' own. The LOX slot's radial width is w = N d_L² / (8 R): 0.85 mm at the central d_L and R = 50 mm (derived).
   - **Expected effect.** A thin annular sheet has more LOX/H2 interface per unit of LOX flow area than N round jets. The ratio of (2/w) to (4/d_L) is 4R/(N d_L), or R / 23.9 mm at the central values (derived): about 0.5 on the innermost ring and about 5 on the outermost, and about 2.8 weighted by flow if the rings carry equal flow at radii 12.5 to 122 mm. The annuli also remove the gaps between elements, where hot products recirculate against the face. Net direction (inferred): faster mixing and heat release near the face on the outer rings, so the flame sits closer to the face than the real face's. c* and Isp would be biased **up**, by as much as the real chamber's combustion is incomplete at the throat, and the heat flux near the face **up**. The size is unknown; tier 2, or a ring-resolved sensitivity run, measures it.
2. **Swirl is neglected in tier 1 (a known omission).** The 162 swirlers spin the LOX inside the post, which leaves as a thinner, wider conical sheet: faster breakup, smaller drops, a wider spray angle. Tier 1 injects axial LOX with Burick's SMD, which is for posts without swirl.
   - **Expected effect (direction inferred):** larger drops and slower vaporization and mixing than the real injector, so a longer flame. c* and Isp are biased **low**, by as much as vaporization is still incomplete well upstream of the throat (nothing if it completes early). The heat flux near the face is biased **low**, since less LOX is carried outward toward the wall. This opposes choice 1's bias near the face; the two are not expected to cancel and are not netted.
   - **Tier-2 candidate:** a 2-D option with swirl, if Mascotte (S7) shows the spray is sensitive. Note that a swirl velocity in the axisymmetric equations is swirl about the chamber's axis. It represents each element's swirl only through the sheet it makes (angle, thickness, SMD), except for an element on the axis.
3. **The transpired H2 is a distributed inflow over the porous face, at the sourced fraction.** About 10% of the H2 (the RL10A-3-3, 1966; `injector_element_inputs.md`) enters at a uniform mass flux over the face area between the ring inlets, at the H2 supply temperature of protocol §4. The ring inlets carry the other 90%, so the elements run at O/F 5.055 / 0.9 = 5.62 (derived). Stoichiometric is 7.94, so the elements burn a little hotter than the engine's mean.
   - **Expected effect (inferred):** a cool, hydrogen-rich layer on the face, which lowers the heat flux into the face and the first part of the wall. The transpired H2 enters slowly and must be entrained into the burning core: a small mixing loss, so c* and Isp are biased slightly **low** against injecting all the H2 through the elements. The elements' higher O/F raises their local flame temperature by a few tens of kelvin (not computed here).

## Known element inputs (from `injector_element_inputs.md`)

216 shear-coaxial elements in 8 rings; 162 swirlers (27 per row is inferred there); swirler ribbon D 0.103 in (2.62 mm) at pitch 0.167 in (Fig. E9); a Rigimesh face with about 10% of the H2 transpired through it (the RL10A-3-3, 1966); a central augmented spark igniter. LOX per element at P1: 13.9525 / 216 = 0.0646 kg/s.

## Declared ranges for the missing element dimensions (all inferred)

Each is varied alone about the central case, as protocol §4 requires; the deviations join the judged band by its envelope rule. The central value is the geometric mean of the bracket's ends, so that no analog is preferred. Ratios are to the LOX post's inner diameter d_L.

| Missing dimension | Central | Bracket | Rule and sources |
|---|---|---|---|
| LOX post inner diameter d_L | 3.54 mm | 2.62 to 4.78 mm | Lower: the swirler ribbon's diameter, 2.62 mm, since the swirler sits in the post's bore (inferred from Fig. E9's description; the bore cannot be smaller). Upper: the SSME main injector's post ID, 0.188 in (TM-82567, "present" design), the largest multi-element flight analog found. Single-element rigs (Mascotte A-10 5.0 mm, DLR BKN 5.5 mm) run at per-element flows far above 0.065 kg/s and are not used for the bound. Falk and Burick's (1973) range table, 1.78 to 3.45 mm, lies within or below the bracket. |
| LOX post tip wall thickness t | 0.272 d_L (0.96 mm) | 0.20 to 0.37 d_L | Falk and Burick's worked example (0.20) and the SSME "present" design (0.37, derived). J-2X subscale tip thicknesses, 0.15 to 0.44 mm, overlap. Post OD = d_L + 2t: central 5.46 mm. |
| H2 annulus gap h (annulus OD = post OD + 2h) | 0.268 d_L (0.95 mm; annulus OD 7.36 mm) | 0.09 to 0.80 d_L | DLR BKN (0.09) to the SSME "present" design (0.80, derived on the assumption that the annulus OD equals the pitch, uncertain). Falk and Burick's range spans 0.04 to 0.59, their example 0.21, Mascotte A-10 0.64. |
| LOX post recess | 0.5 d_L | 0 to 1.0 d_L | SP-8089 §3.1.1.1.4: recess "no greater than one post diameter", which Burick's data (CR-120936) put at maximum mixing. DLR BKC is flush (its post stands 3 mm proud), DLR BKN is at 0.73 d_L. The midpoint is central because 0 makes the geometric mean undefined. |
| Ring radii | 8 rings equally spaced in radius between r_1 and r_8; r_8 = R_c − p/2 with p the element pitch from face area, √(π R_c² / 216) = 15.7 mm (R_c = 5.13 in) | r_8 from R_c − p/2 to R_c − p; r_1 from the inner ring's clearance (the elements of that ring must fit its circumference at the annulus OD) to r_8 − 7p | Pitch from face area per element (derived). It agrees with the SSME's 16.0 mm, while the analogs span 16 to 51 mm. Enters tier 1 as the ring inlets' radii (choice 1) and the supply turbulence's length scale (one pitch, §4). |

**Derived checks, not inputs.** At P1, LOX at about 1,140 kg/m³ (101 K) leaves the post at 10.5 m/s for d_L = 2.62 mm, 5.8 m/s at the central value and 3.2 m/s at 4.78 mm. The H2 annulus velocity follows from the annulus area and the H2 density at the injector-face pressure, which tier 1 computes; it is not set here.

## Open before the prediction commit

- **The ring layout (resolved in part by choice 1).** With continuous annuli, 27 elements need not fit a ring's circumference, so the fit conflict no longer arises: 216/8 = 27 per ring with an annulus OD of 7.36 mm needs a ring radius of at least 32 mm (27 × 7.36 / 2π), while rings one pitch apart from r_8 = R_c − p/2 put r_1 at 12.5 mm (derived). The ring radii remain the declared rule in the table above (8 rings equally spaced, r_1 = r_8 − 7p, central). Still open: the per-ring counts, if `injector_element_inputs.md`'s fuel-only outer row and the inner row without a swirler are not 27-element coaxial rings. The central case stays 27 per ring until that is settled.
- **Swirl:** declared, choice 2.
- **The transpired H2:** declared, choice 3.
- **The igniter's dimensions** are missing; the porous face (choice 3) covers the centre inside r_1.

## The cup's resolution criterion (decision 5), stated before any tier-1 run (21 October 2026, about 02:00 CDT)

The recess enters the central case as geometry (a cup) only if the tier-1 mesh resolves the cup (the plan's decision 5). The cup is the LOX slot's exit set back by the recess, 0.5 d_L = 1.77 mm at the central values. Its walls are the outer edges of the ring's two H2 slots. Its width runs from the inner H2 slot's inner edge to the outer H2 slot's outer edge: about 2.9 mm on ring 8 to 13.9 mm on ring 1 (derived, `tier1_setup.md`).

**Criterion (engineering; the counts are *guessed* minima):**
- The cup is resolved if, at the mesh's finest level, every ring's cup has at least **4 cells across its depth**, and at least **2 cells across each of its slots and post-tip walls**, so at least 10 across its width.
- The mesh must also carry the cup's walls: cells set back behind the face plane, bounded by walls.
- The reasoning:
  - With 4 cells across the depth, the flow inside the cup has an axial development, not a one-cell buffer.
  - With 2 cells across each slot and tip, the H2 streams and the post-tip wake between them stay distinct streams with edges.
  - Fewer cells than this, and the cup is a sub-cell feature whose geometry the solution cannot carry.
- **Verdict rule:** if any ring's cup fails either count, or the mesh cannot carry the walls, the cup is not resolved. Then decision 5's branch applies, dated before any run:
  - the central case takes the recess through the SMD correlation alone;
  - member 11b is not run;
  - the missing cup is a declared limit, with its stated direction.

## The prediction (empty until the prediction commit)

F_vac, Isp_vac (chamber and engine-inlet flow), nozzle stagnation and injector-face pressures, c*, and the hot-gas wall heat flow, each with its central value, judged band and quadrature sum (protocol §3 and §4), and the commit hash recorded before any comparison reads a measured value (§5).
