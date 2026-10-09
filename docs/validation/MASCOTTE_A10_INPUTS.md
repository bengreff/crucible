# Mascotte A-10: input extraction for blind validation

Date: 8 October 2026, about 22:30 CDT. Extracted by a Sonnet research agent under the blind rule (inputs, boundary conditions and the names and stations of measured quantities only). Its file was screened for outcome values before it was read in full and committed: none found. The extraction is kept as the agent wrote it, from its first section heading on.

Modelling notes for the case (mine, not the agent's): the chamber is square, 50 by 50 mm, and the engine is axisymmetric, so the case will need a declared equivalent round chamber; the windows' helium film is not in the engine's medium. Both are decisions for the case's pre-registration (docs/evidence/SPRAY.md), not settled here.

## Sources used

- [N] A. Nicole, G. Ordonneau, M. Theron (ONERA/CNES), "3D Simulation of LOX/GH2 MASCOTTE Test
  Case at 10 bar," EUCASS proceedings (undated reprint, work describes a 2008-era study).
  https://www.eucass.eu/component/docindexer/?task=download&id=6352
  — cites, as its own ref. [4]: L. Vingert, "Sélection de cas tests MASCOTTE – cas subcritique
  LOX/GH2 (10 bar) et cas supercritique LOX/GH2 (60 bar)," ONERA RT n°1/11785, November 2006
  (the RCM test-case description document; not independently retrieved — ONERA-internal report,
  not found openly online).
  — also cites: M. Habiballah, L. Vingert, J.C. Traineau, P. Vuillermoz, "MASCOTTE: A test bench
  for cryogenic combustion research," 47th IAC, Beijing, Oct 1996 (not independently retrieved).
- [R] H. Riedmann, D. Banuti, B. Ivancic, O. Knab, K. Hannemann, "Modeling of H2/O2 Single-Element
  Rocket Thrust Chamber Combustion at Sub- and Supercritical Pressures with Different CFD Tools,"
  6th EUCASS, 2015. https://www.eucass.eu/component/docindexer/?task=download&id=3819
  (its §2.1 "Mascotte A-10 Test Case" explicitly states it is a short summary of Vingert &
  Habiballah's test-case description [their ref. 4], i.e. the RCM-1/RCM-2 workshop document).
- [F] N. Fdida, L. Vingert, A. Ristori, Y. Le Sant (ONERA), "Droplet size and velocity measurements
  in a cryogenic jet flame of a rocket-type combustor using high speed imaging," draft for
  Atomization & Sprays journal, arXiv:2408.17132.
- [M] J. Martin, W. Armbruster, R. Stützer, S. General, B. Knapp, D. Suslov, J. Hardi (DLR),
  "Flame characteristics of a high-pressure LOX/H2 rocket combustor with large optical access,"
  Case Studies in Thermal Engineering 28 (2021) 101546. Used only for §7 (alternative case, DLR
  P8/BKN), not for Mascotte geometry.
- [WS] Values that only came back as web-search-engine summaries of secondary sources, which I
  could not independently verify by reading the primary PDF/table myself. Flagged explicitly;
  treat as lower confidence until checked against [N]'s ref. [4] or Habiballah et al. 2006
  (Combust. Sci. Technol. 178, 101-128) directly.

Chamber-pressure note: ~10 bar / 1.0 MPa is the nominal/target operating point by which the
"A-10" case is defined (stated identically in [N] and [R]). [F] Table 1 additionally reports an
as-recorded steady-state chamber-pressure sensor average of 1.00 MPa (fire test) / 0.96 MPa
(cold He flow test) from one specific test campaign — this is a recorded operating-condition
sensor reading (facility pressure transducer), not a flame-diagnostic outcome, so it is reported
here as an input/operating record, distinguished from the nominal set point.

---

## 1. Injector geometry

| Parameter | Value | Label | Source |
|---|---|---|---|
| LOX post (oxygen tube) diameter at exit, "d" | 5 mm | stated | [R] p.2: "The diameter of the oxygen tube is 5 mm"; [F] p.4: "central tube of diameter d" (d used as the injector's reference length) |
| LOX post inlet (upstream) diameter | 3.6 mm (tapering to 5 mm at exit) | WS, unverified | Returned by search-engine summaries of a related EUCASS LOX/CH4-vs-LOX/H2 comparison paper; not read directly from a primary table. Flag for verification. |
| LOX post outer diameter | ≈5.6 mm | derived + WS | Derived: [F] states the GH2 annulus inner diameter = 1.12·d = 1.12×5 mm = 5.6 mm; if the annulus inner wall is flush with the LOX-post outer wall (standard shear-coaxial design), post outer diameter ≈ 5.6 mm. A search-summary of a secondary source stated the same 5.6 mm figure directly; not independently verified against a primary table. |
| LOX post lip thickness | ≈0.3 mm | derived by us | (5.6 − 5.0) / 2 = 0.3 mm, from the above outer/inner diameters. Arithmetic only; not a stated value. |
| GH2 annulus inner diameter | 5.6 mm | stated | [F] p.4: outer tube fed with hydrogen, "inner diameter 1.12d" = 5.6 mm (d = 5 mm); [R] p.2: "the diameters of the hydrogen annulus are 5.6 mm and 12 mm" |
| GH2 annulus outer diameter | 12 mm | stated | [R] p.2 (A-10 case, as distinct from the A-60 case where outer diameter is reduced to 10 mm); [F] p.4: "outer diameter 2.4d" = 12 mm |
| GH2 annulus radial gap width | 3.2 mm | derived | (12 − 5.6)/2 = 3.2 mm, from the above |
| Recess | not confirmed for A-10 specifically | unresolved | Literature on later Mascotte campaigns (LOX/CH4 work) explicitly distinguishes "with recess" vs "without recess (flush)" injector variants, implying the original A-10/LOX-GH2 injector is the flush (zero-recess) baseline, but no source read directly states a recess value of 0 mm for A-10. Do not assume; verify against Vingert ONERA RT 1/11785 or Habiballah et al. 2006. |
| LOX post length | ≈50 mm (claimed) | WS, unverified — low confidence | Only found via a search-engine summary ("the injector length is equal to 50 mm"); not confirmed in any primary PDF read directly. Flag for verification. |
| Swirl | None — shear coaxial injector | stated | [N], [R], [F] all describe the injector as a "shear coaxial" element (LOX in the center, GH2 in a surrounding annulus); no swirl device mentioned anywhere. |

## 2. Propellants

| Parameter | Value | Label | Source |
|---|---|---|---|
| LOX mass flow rate | 50 g/s | stated | [N] p.2 (nominal A-10 point, citing ref. [4]); [R] p.2 ("The oxygen mass flow rate is 50 g/s") |
| LOX injection temperature | 85 K (nominal, [N],[R]); 89–90 K (as-recorded, [F]) | stated | [R] p.2: "a temperature of 85 K"; [N] p.2 implies LOX injected cold via LN2 pre-cooling; [F] Table 1: Tinj(LOX) = 89 K (fire test), 90 K (cold-flow test), as-recorded averages |
| GH2 mass flow rate | 23.7 g/s | stated | [N] p.2; [R] p.2 ("the hydrogen mass flow rate is 23.7 g/s") |
| GH2 injection temperature | 287 K (nominal, [R]); 296 K (as-recorded fire test, [F]); "room temperature" ([F] text) | stated | [R] p.2; [F] p.4 and Table 1 |
| He film-cooling flow rate (window cooling, not a propellant) | ≈10.12 g/s ([N]) or ≈6 g/s ([F]) | stated (two different test campaigns, discrepant) | [N] p.2: "10.12 g/s... for helium"; [F] p.4: "a film cooling of helium was injected with a mass flow rate of ≈6 g/s." Both are stated in their respective primary texts; the discrepancy is noted, not resolved, since these may be different test campaigns/hardware configurations both labeled "A-10." |
| Mixture ratio O/F | 2.11 (nominal, [N],[R]); 2.4 (as-recorded fire-test average, [F] Table 1); ROF≈2 ("nominal" per [F] text) | stated | [N] p.2; [R] p.2; [F] p.4 and Table 1 |
| Momentum flux ratio J | 12.63 ([N]) vs 15.4 ([R]) vs 12.99 (as-recorded, [F] Table 1) | stated (discrepant between sources — not reconciled) | [N] p.2: "the momentum flux ratio J is equal to 12.63"; [R] p.2: "This results in a momentum ratio of J = 15.4"; [F] Table 1: J = 12.99 (fire), 12.38 (cold flow) |
| Velocity ratio (GH2/LOX) | 146.3 | stated | [R] p.2: "a velocity ratio of 146.3... between fuel and oxidizer" |
| LOX post-based liquid Reynolds number Re_d | ≈62,300 (fire) / 64,300 (cold flow) | stated | [F] Table 1 |
| Gaseous Weber number We_G | ≈26,400 (fire) / 25,300 (cold flow) | stated | [F] Table 1 |
| Injection velocities (individual LOX, GH2) | not directly stated as a single pair of numbers in any source read | — | Could be back-derived from mass flow, density and post/annulus areas, but not attempted here to avoid introducing an unlabeled derived "result"; flag as an exercise for whoever builds the case, not reported as fact. |

## 3. Chamber

| Parameter | Value | Label | Source |
|---|---|---|---|
| Cross-section shape | Square | stated | [N] p.1: "50 mm x 50 mm square internal section"; [R] p.2: "square cross section with an edge length of 50 mm"; [F] p.4: "square cross section of 50 mm side" |
| Cross-section size | 50 mm × 50 mm | stated | as above |
| Chamber length, injector exit plane to throat | 478 mm ([N],[R]) vs "about 500 mm" ([F]) vs 472 mm (a web-search paraphrase of the same RCM-2 description, unverified) | stated, source-discrepant | [N] p.1: "The total length of the chamber from the injector exit section to the throat is 478 mm"; [R] p.2: "The distance between the injection plane and the throat is 478 mm"; [F] p.4: "was about 500 mm in length." The 472 mm figure came only from a search-engine paraphrase and is not independently confirmed — treat 478 mm ([N],[R], same underlying citation to Vingert & Habiballah) as the better-supported figure. |
| Injector overtake into the chamber | 10 mm (injector protrudes 10 mm into the chamber so its exit is visible through the portholes) | stated | [N] p.1 |
| Nozzle throat diameter | 15 mm | stated | [N] p.1; [R] p.2 |
| Nozzle/chamber transition | A transition segment connects the square chamber to the axisymmetric (round) nozzle | stated | [F] p.4 |
| Area contraction ratio chamber:throat, εc | 14.1 | derived (by [R]'s authors, not us) | [R] p.2 |
| Characteristic chamber length l* = Vc/At | 6.6 m | derived (by [R]'s authors) | [R] p.2 |
| Windows | Side windows, helium-film cooled, for optical diagnostics; two optical ports/visualization module at the injector exit in the imaging study [F]; three selectable axial positions of the visualization module referred to as "Upstream," "Medium," "Downstream" in [N] Fig. 8 (exact mm offsets of these three positions not given in [N]'s visible text) | stated | [N] p.1, p.6 (Fig. 8); [F] p.4 |
| Window dimensions | 100 mm long × 50 mm wide × 10 mm thick, made of fused/"melt" silica | stated | [N] p.6 |
| Chamber pressure (operating set-point) | ≈10 bar / ≈1.0 MPa (nominal, defines the "A-10" designation) | stated — set-point/nominal | [N] p.1: "the A-10 point which corresponds to... a chamber pressure of approximately 10 bar"; [R] p.2: "a subcritical chamber pressure of 10 bar" |
| Chamber pressure, as-recorded sensor average | 1.00 MPa (fire test), 0.96 MPa (cold He-flow test) | stated — recorded operating condition (facility pressure transducer reading, not a flame-diagnostic outcome) | [F] Table 1 |

## 4. Ignition

| Parameter | Value | Label | Source |
|---|---|---|---|
| Ignition method | Gaseous-oxygen/gaseous-hydrogen torch igniter, placed in the first ferrule (near the injector face) | stated | [F] p.4: "Ignition is provided by a gaseous oxygen/ gaseous hydrogen torch igniter placed in the first ferrule." |
| Timing/energy | Not stated in any source read | — | — |

## 5. Wall conditions

| Parameter | Value | Label | Source |
|---|---|---|---|
| Main chamber wall cooling | Not explicitly stated in the sources read (windows are helium-film cooled; main chamber liner cooling method not found) | unresolved | — |
| Window material / thermal spec | Fused ("melt") silica; the paper separately discusses a softening-temperature figure for that material as a design/material property, not a flame measurement | stated (material spec) | [N] p.6 |
| Simulation wall boundary condition used by [R] (NOT an experimental wall condition — modeling choice only, flagged as such) | Adiabatic no-slip walls on faceplate and liner | stated, but this is a CFD modeling assumption, not a measured experimental wall condition | [R] p.5 |

## 6. Measurements (names and axial positions only — no measured values)

| Quantity | Instrument/technique | Axial station(s) | Source |
|---|---|---|---|
| Gas-phase temperature | CARS (Coherent Anti-Stokes Raman Spectroscopy) thermometry | x = 10 mm, x = 250 mm, x = 410 mm (downstream of injector exit) | [N] p.2 (CARS), p.7–9 (station positions, figure axis labels only) |
| Flame emission imaging | OH* chemiluminescence emission imaging (average image + Abel-transformed image) | Full chamber length viewed through the optical windows ("observation plane") | [N] p.4, p.9; [R] p.3 (citing Candel et al. and Grisch et al. as the source publications for this diagnostic) |
| Species concentration (named only) | H2 concentration measurements (Candel et al.); H2O concentration measurements (Grisch et al.) — both cited as existing published diagnostics for A-10, used elsewhere for comparison | Not stated in the text read | [R] p.3 |
| Droplet size and velocity | PDPA (Phase Doppler Particle Analyzer) | Atomization-zone region, reported relative to injector/post diameter d: 6 < x/d < 12 (axial), 2 < r/d < 4 (radial, off-axis) — this is the measurement coverage region, not a measured result | [F] p.2–3 (PDPA work credited to Gicquel & Vingert, 2000); [F] p.2 (this paper's own imaging-based drop-size/velocity study covers the same region using high-speed shadowgraphy, PTV and PIV/FOLKI-SPIV instead of PDPA) |
| Pressure/temperature/flow-rate sensors (facility instrumentation) | Thermocouples, pressure transducers, flow meters on propellant lines and chamber | Not detailed beyond "on the test bench" | [F] p.5 |

## 7. Other subcritical (<5 MPa) LOX/GH2 single-element spray-flame cases — candidates

| Case | Pressure | Completeness of published inputs (no outcomes) |
|---|---|---|
| Mascotte, other operating points (e.g. a "C10" point mentioned alongside A10 in one secondary source, both at 10 bar) | ~10 bar | Only seen as a passing mention ("simulated cases correspond to the LOx/GH2 Mascotte A10 and C10, both operating at 10 bar with two-phase flow conditions"); full input set not retrieved — would need the specific paper defining C10. |
| DLR P8 test facility, single-injector research combustor "BKN" (Martin et al. 2021, [M]) | Sub- to supercritical, chamber pressures stated up to 70 bar, O/F 3.4–6 in that paper's own tests (so the exact sub-critical, Mascotte-A10-like O/F≈2 regime is not this paper's operating point, but the facility/injector is explicitly built to cover both regimes) | Very complete and clearly tabulated: full injector nomenclature (LOX post orifice/length/taper angle, fuel-annulus gap, recess length all defined as named geometric parameters — see [M] Fig. 2 and its Table 1, not reproduced here), chamber geometry (50 mm bore, 362 mm length, 16.9 mm throat), extensive instrumentation (thermocouples, static/unsteady pressure at many axial stations), ignition (torch or laser), and a 253.5 mm-long full-length optical window. One of the best-documented modern alternatives. |
| Penn State preburner single-element combustor (often itself called "RCM-1" in the RCM workshop series) | Chamber pressures 300–750 psia (~2.1–5.2 MPa) in various sub-cases | Well published (multiple CFD-validation papers reproduce the geometry and flow table), but NOTE: this is a GASEOUS O2/GASEOUS H2 (gas-gas, "preburner") case, not a liquid-oxygen spray flame — it does not meet the "LOX" requirement of this task and is listed only because the task prompt explicitly asked about it; it is not a substitute LOX spray case. |
| DLR M3 test infrastructure | Not a single fixed operating point — M3 is a general multi-purpose high-pressure test facility (per "Overview on the DLR M3 Test Infrastructure") | Facility-level overview only; would need a specific published single-injector LOX/GH2 sub-critical test case run on M3 to be usable — none identified in this search. |

---

### Blindness note

Outcome values were visible on two sources while extracting inputs, and were not recorded:
- Source [N] (A. Nicole, G. Ordonneau, M. Theron, EUCASS "3D Simulation of LOX/GH2 MASCOTTE Test
  Case at 10 bar"), pages 7–9: Figures 11–13 show CARS-measured and simulated temperature-vs-height
  profiles at x = 10, 250, 410 mm (with "validation rate" percentages annotated on the plots), and
  Figure 14 shows an OH mass-fraction/OH-emission comparison image. Not recorded.
- Source [R] (Riedmann et al., EUCASS 2015), pages 2–3: Figure 2 ("Mascotte A-10: average OH*
  emission image... and Abel-transformed emission image") and Figure 3 (equivalent A-60 images)
  display actual emission-image outcome data. Not recorded/described beyond noting that such
  images exist as a diagnostic type (already a legitimate "what was measured" entry in §6).
