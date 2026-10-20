# RL10A-3-3A tier 1: pre-registration skeleton and declared injector ranges

Status: skeleton written 20 October 2026, about 18:00 CDT, before any tier-1 run and before any measured performance of this engine was read. It holds no predicted number. It extends the protocol (`RL10A33A_PROTOCOL.md`, written 13 October in inquiry2's line and not yet on `main`). Its §4 brackets and §5 comparison rules govern; nothing here changes them. What it adds: tier 1's settings to be fixed, and declared ranges for the injector-element dimensions that the inputs lack (`injector_element_inputs.md`). Every range below is **inferred** from analogous coaxial LOX/H2 injectors (`analog_coax_elements.md`, a Sonnet research agent, 20 October), not from any RL10 document.

## Declarations

- **The analog agent's incidental exposures (agent only; nothing recorded).**
  - A generic "coaxial" search of NASA SP-8089's OCR text showed one survey-table row for "Centaur, RL10A-3" with performance figures and an element count. The agent recorded no value from it and excluded the table.
  - NASA CR-187037 sizes its test element against "RL-10 injector elements"; the agent discarded that source unread beyond that sentence.
  - The author of this file read only the agent's file, after a count-only keyword grep and a digit-masked view of every hit. No RL10 value is in it.
- **The protocol's face distribution stands.** §4 sets a uniform face (O/F and mass flux), with no bracket. The element dimensions below enter tier 1 through the injection velocities, the spray's drop size and the supply turbulence, not through a resolved face.

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
| Ring radii | 8 rings equally spaced in radius between r_1 and r_8; r_8 = R_c − p/2 with p the element pitch from face area, √(π R_c² / 216) = 15.7 mm (R_c = 5.13 in) | r_8 from R_c − p/2 to R_c − p; r_1 from the inner ring's clearance (the elements of that ring must fit its circumference at the annulus OD) to r_8 − 7p | Pitch from face area per element (derived). It agrees with the SSME's 16.0 mm, while the analogs span 16 to 51 mm. Enters tier 1 only through the supply turbulence's length scale (one pitch, §4) and a ring-resolved check, if one is run. |

**Derived checks, not inputs.** At P1, LOX at about 1,140 kg/m³ (101 K) leaves the post at 10.5 m/s for d_L = 2.62 mm, 5.8 m/s at the central value and 3.2 m/s at 4.78 mm. The H2 annulus velocity follows from the annulus area and the H2 density at the injector-face pressure, which tier 1 computes; it is not set here.

## Open before the prediction commit

- **The ring layout is under-determined.** 216/8 = 27 per ring matches the 27-per-row reading of the 162 swirlers. But 27 elements of annulus OD 7.36 mm need a ring radius of at least 32 mm (27 × 7.36 / 2π), while rings one pitch apart inward from r_8 = R_c − p/2 put r_1 at 12.5 mm. So a radial pitch of p and 27 per ring cannot both hold (derived). The ring rule above becomes a declared choice, or the ring-resolved check is dropped and the uniform face stands alone.
- **Swirl.** The 162 swirlers impart LOX swirl. The 2-D axisymmetric engine carries no swirl velocity, so the swirlers enter only through Burick's SMD, which is for non-swirled posts. Declared as a modelling limit: no bracket can be set from the inputs.
- **The transpired H2 (about 10% through the Rigimesh face):** whether tier 1 injects it through the face between elements (at the face's uniform mass flux, so already inside a uniform face) or through the annuli. To be stated.
- **The igniter's dimensions** are missing; the uniform face covers the centre.

## The prediction (empty until the prediction commit)

F_vac, Isp_vac (chamber and engine-inlet flow), nozzle stagnation and injector-face pressures, c*, and the hot-gas wall heat flow, each with its central value, judged band and quadrature sum (protocol §3 and §4), and the commit hash recorded before any comparison reads a measured value (§5).
