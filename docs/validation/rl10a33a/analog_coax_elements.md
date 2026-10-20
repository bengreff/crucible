# Analog shear-coaxial injector element dimensions (RL10A-3-3A gaps)

Written by a Sonnet research agent, 20 October 2026. Analogs only: no RL10 document of any
variant was opened, searched for, or read. Geometry, element counts and dimensionless geometric
ratios only; no thrust, Isp, c*, efficiency, chamber-pressure or heat-flux values are recorded for
any engine. Propellant temperatures and injection velocities of the analog engines are kept in a
separate, labelled column only where printed next to geometry.

**Blindness note (read first).** A generic keyword search ("coaxial") inside the OCR'd text of
NASA SP-8089 returned a summary table (state-of-the-art survey, p.9 area) that includes one row
for "Centaur, RL10A-3" with performance figures (Isp efficiency, thrust, chamber pressure, mixture
ratio) and an element count. This was an unavoidable side effect of searching a 130-page scanned
document for a generic term, not a targeted RL10 search. No value from that row was recorded,
copied, or used anywhere below, and the surrounding table itself is not used as a source for any
other engine's numbers in this file (the J-2 and SSME element counts below come from independent,
non-RL10 sources). Flagged per the blind-validation protocol, the same way
`docs/validation/MASCOTTE_A10_INPUTS.md`'s own "Blindness note" handles an analogous incident.
One other SP-8089 line ("concentric-tube elements are used in the RL-10, M-1, and J-2 engines")
was seen in a context grep; it states only that RL-10 uses a concentric-tube (shear-coaxial)
element, which the RL10A-3-3A case record already states as known context, and carries no new
geometry value, so it is not treated as a disclosure.

Labels follow the task's scheme: *drawing/table dimension*, *text statement*, *estimate/derived by
the source*, or *derived here* (this agent's arithmetic).

---

## 1. Mascotte A-10 / C-10 (ONERA, single-element LOX/GH2 rig, ~1 MPa)

Repo source: `docs/validation/MASCOTTE_A10_PREDICTION.md` and `docs/validation/MASCOTTE_A10_INPUTS.md`
(geometry sections only; no measured/outcome value from those files is reproduced here). Original
references as cited there: [R] Riedmann et al., 6th EUCASS 2015; [F] Fdida et al., arXiv:2408.17132;
[N] Nicole, Ordonneau, Theron, EUCASS (undated reprint).

| Parameter | Value | Label | Source |
|---|---|---|---|
| LOX post inner diameter, d | 5.0 mm | text statement | [R] p.2, [F] p.4 |
| LOX post outer diameter | ≈5.6 mm | derived + WS (repo's own flag) | Repo-derived from [F]'s annulus-ID statement, assuming flush post/annulus wall; a search-summary of a secondary source gave the same figure independently but unverified against a primary table |
| LOX post lip (wall) thickness | ≈0.3 mm | derived by repo | (5.6 − 5.0)/2 |
| GH2 annulus inner diameter | 5.6 mm | text statement | [F] p.4, [R] p.2 |
| GH2 annulus outer diameter | 12.0 mm | text statement | [R] p.2 ("A-10" case), [F] p.4 |
| GH2 annulus radial gap width | 3.2 mm | derived by repo | (12.0 − 5.6)/2 |
| LOX post recess (relative to its own faceplate) | not confirmed for A-10 | unresolved in repo's own extraction | repo flags this explicitly; do not assume flush |
| Injector body protrusion into the chamber (a different parameter — not post recess) | 10 mm | text statement | [N] p.1 |
| Chamber cross-section at the injector face | 50 × 50 mm square (round-equivalent radius 28.0 mm used in the repo's model) | text statement / repo modelling choice | [N],[R],[F] |
| Swirl | none (shear-coaxial) | text statement | [N],[R],[F] |
| Mascotte C-10 (same injector family) post ID / OD | 5.0 mm / 5.6 mm | WS, secondary source, unverified against a primary table | search-engine summary of a secondary figure caption; consistent with A-10's own numbers but not independently confirmed here |

**Derived ratios here:** gap/d_LOX = 3.2/5.0 = 0.64; lip-thickness/d_LOX = 0.3/5.0 = 0.06;
annulus-OD/d_LOX = 12.0/5.0 = 2.4.

## 2. DLR BKN (single-element, large-optical-access LOX/H2 research chamber, up to 70 bar)

Source: Martin, Armbruster, Stützer, General, Knapp, Suslov, Hardi, "Flame characteristics of a
high-pressure LOX/H2 rocket combustor with large optical access," *Case Studies in Thermal
Engineering* 28 (2021) 101546 (already cited in the repo as [M] in `MASCOTTE_A10_INPUTS.md` §7,
there "used only for §7," not for Mascotte geometry — read directly here for its own Table 1).
PDF fetched from elib.dlr.de and OCR'd/text-extracted locally; Table 1, p.5.

| Symbol | Meaning | Value | Label |
|---|---|---|---|
| d₍o,or₎ | LOX post metering-orifice diameter | 2.2 mm | drawing/table dimension (Table 1) |
| l₍o,or₎ | LOX post orifice length | 4.0 mm | drawing/table dimension |
| d₍o,po₎ | LOX post diameter (main body) | 5.5 mm | drawing/table dimension |
| l₍o,po₎ | LOX post length | 127.0 mm | drawing/table dimension |
| d₍o,ex₎ | LOX post exit (tip) diameter, after an outward taper | 6.3 mm | drawing/table dimension |
| α | LOX post taper angle | 8° | drawing/table dimension |
| l₍f,an₎ | fuel annulus length | 34.5 mm | drawing/table dimension |
| d₍f,ex₎ | fuel annulus exit diameter | 8.2 mm | drawing/table dimension |
| h₍f,gap₎ | fuel annulus gap height | 0.5 mm | drawing/table dimension |
| d₍c,ex₎ | window-cooling annulus exit diameter (≈ chamber bore) | 50.0 mm | drawing/table dimension |
| h₍c,gap₎ | window-cooling annulus gap height | 1.0 mm | drawing/table dimension |
| l₍re₎ | LOX post recess length | 4.0 mm | drawing/table dimension |
| — | chamber length, faceplate to throat | 362 mm | text statement (p.5) |
| d₍throat₎ | nozzle throat diameter | 16.9 mm | text statement (p.5) |

Note on scope: the paper states this element "was designed as a hybrid element optimized for
methane but also useable for hydrogen as fuel," and the tests reported are at ROF 3.4–6 (a
staged-combustion-like, fuel-lean-biased regime, unlike RL10-class gas-generator ROF); flagged as
an operating-class caveat, not a geometry caveat.

**Derived ratios here:** recess/d₍o,po₎ = 4.0/5.5 = 0.73; fuel-gap/d₍o,po₎ = 0.5/5.5 = 0.09;
(d₍o,ex₎ − d₍o,po₎)/d₍o,po₎ = 0.145 (taper growth fraction); window-cooling-annulus-OD / d₍o,po₎ = 50.0/5.5 = 9.1.

## 3. DLR BKC (single-element LOX/H2 and LOX/CH4 research chamber, sub- to supercritical)

Source: Suslov, Hardi, Knapp, Oschwald, "High-speed imaging of the disintegration of an LOX jet..."
*Aerospace Science and Technology* 86 (2019) 444–454 (text-extracted locally from elib.dlr.de PDF).

| Parameter | Value | Label |
|---|---|---|
| Chamber inner diameter | 50 mm | text statement (p.446) |
| Contraction ratio | 8.35 | text statement |
| LOX post recess or taper | **none** (flush, cylindrical) | text statement ("the LOX post has no recess or taper") |
| Injector outlet vs. faceplate | stands proud (protrudes) by **3 mm** | text statement — a protrusion, i.e. a *negative* recess |
| LOX post diameter, D | not found in extracted text (shown only in a figure, not machine-readable here) | not found |

This is a useful negative data point: at least one modern DLR single-element LOX/H2 research
injector uses zero recess / a small protrusion rather than a recessed post, for direct comparison
against the BKN example above (recess 0.73 d) and against Falk & Burick's optimum (below, 1 d).

## 4. SSME main injector, baseline ("present") design and a proposed low-loss redesign

Source: NASA TM 82567, "Low Loss Injector for Space Shuttle Main Engine" (1984; downloaded from
NTRS, text-extracted locally with `pdftotext`). Figures/tables referenced are p.11–17 in the
source numbering as printed.

| Parameter | Present (baseline/flown) design | Proposed low-loss design | Label |
|---|---|---|---|
| LOX post outer diameter | 0.329 in = 8.36 mm | 1.25 in = 31.75 mm | text statement |
| LOX post inner diameter | 0.188 in = 4.78 mm | 0.80 in = 20.32 mm | text statement |
| LOX post wall thickness (derived) | (8.36−4.78)/2 = 1.79 mm | (31.75−20.32)/2 = 5.72 mm | derived here |
| Closest element-to-element spacing diameter, D (≈ hexagonal pitch) | 0.63 in = 16.0 mm | 2.0 in = 50.8 mm | text statement |
| Circle-row count, R | 13 | 4 | text statement |
| Chamber inner diameter, C | 17.74 in = 450.6 mm (both designs) | same | text statement |
| LOX post length (periphery / center) | 5.935 in / 3.334 in = 150.7 / 84.7 mm (both designs) | same | text statement |
| Lox-plug metering sieve (proposed design only) | — | 36 holes, 0.084 in = 2.13 mm diameter | text statement |
| Hot-gas shroud sieve (proposed design only) | — | 48 holes, 0.20 in = 5.08 mm diameter | text statement |

**Derived ratios here (present design):** wall-thickness/post-ID = 1.79/4.78 = 0.37;
(pitch − post-OD)/2 as an approximate fuel-annulus-gap estimate = (16.0−8.36)/2 = 3.82 mm, so
gap/post-ID ≈ 3.82/4.78 = 0.80 — **this gap figure assumes the annulus OD equals the hexagonal
pitch (elements just touching), which is not separately confirmed in the source; flagged as a
rough estimate, not a drawing dimension.**

SSME/RS-25 element count (600) and J-2 element count (614) are well-attested by general
secondary/public sources (Wikipedia-level, NASA historical summaries) independent of the above;
no LOX post diameter for J-2 was found in this search (see §8).

## 5. Falk and Burick (1973), NASA CR-120968 / R-8973-3, "Injector Design Guidelines for Gas/Liquid Propellant Systems"

Downloaded in full from NTRS (fulltext OCR layer) and read directly — this is also the primary
source behind NASA SP-8089's own recess guidance (SP-8089 ref. 14 is Burick's companion 1972 report,
CR-120936; this 1973 report by the same two authors is the general gas/liquid design-guideline
document and gives an explicit cold-flow range table plus a full worked coaxial-element design
example). Not an RL10 source; written for a FLOX/light-hydrocarbon design study, but the circular
coaxial element data and the worked example are propellant-agnostic geometry.

**Table II, "Range of Cold-Flow Variables for Circular Coaxial Element" (p.17 as printed):**

| Parameter | Symbol | Range | Label |
|---|---|---|---|
| Liquid jet (oxidizer) diameter | D_L | 0.070, 0.108, 0.136 in (1.78, 2.74, 3.45 mm) — three discrete sizes tested | drawing/table dimension |
| Gas gap height (fuel) | h | 0.005 to 0.041 in (0.13 to 1.04 mm) | drawing/table dimension |
| Oxidizer post recess | R | 0 to 1 × D_L | drawing/table dimension |
| Gas velocity in annulus (fuel) | V_g | 100 to 630 ft/s (30.5 to 192 m/s) — velocity, kept separate per task instructions | drawing/table dimension |
| Liquid velocity at post tip (oxidizer) | V_L | 5 to 100 ft/s (1.5 to 30.5 m/s) — velocity, kept separate | drawing/table dimension |
| Injected mixture ratio | MR | 3 to 7.5 | drawing/table dimension (not a geometry value; kept for context only) |

**Worked design example (p.62–64, §6.4.2, a 148-element, 5000-lbf-class coaxial injector in a
7.90 in = 200.7 mm injector face):**

| Parameter | Value | Label |
|---|---|---|
| Liquid jet (LOX post) inner diameter, D_L | 0.100 in = 2.54 mm | text statement (design choice) |
| LOX post wall thickness | 0.020 in = 0.508 mm | text statement |
| LOX post outer diameter = fuel annulus inner diameter | 0.140 in = 3.56 mm | derived by the source (OD = ID + 2×wall) |
| Fuel annulus outer diameter | 0.182 in = 4.62 mm | derived by the source (from required annulus area) |
| Fuel annulus gap height | 0.021 in = 0.533 mm | derived by the source |
| Recess (optimum, per Table II's correlation, applied to this example) | ≈1 × D_L ≈ 0.100 in = 2.54 mm | estimate/derived by the source |
| Element count / face diameter | 148 elements / 7.90 in = 200.7 mm | text statement |

**Key design-practice finding (text statement, p.62): "optimum mixing was obtained with an element
recess of approximately 1 × D_L"** — directly matching, and the likely origin of, NASA SP-8089's
own recess guidance (below).

**Derived ratios here (from the worked example):** wall-thickness/D_L = 0.020/0.100 = 0.20;
annulus-gap/D_L = 0.021/0.100 = 0.21; post-OD/D_L = 0.140/0.100 = 1.40; recess/D_L ≈ 1.0.

## 6. J-2X / RS-68B subscale shear- and swirl-coaxial injector testing

Source: Protz, C.; Elam, S.: "Subscale Injector Testing to Support Ares Engines Development," NASA
MSFC report MSFC-2184 (NTRS 20090014736) — **abstract only**; the full report was not retrievable
(NTRS serves only an abstract for this submission; confirmed via the NTRS API, `onlyAbstract` not
set but no downloadable PDF/fulltext link present). Numbers below are quoted from the abstract
text itself, not a secondary paraphrase.

| Parameter | Value | Label |
|---|---|---|
| Shear-coaxial injector element densities tested | 1.1, 1.6, 2.3 elements/in² | text statement (abstract) |
| Swirl-coaxial injector element density tested | 1.6 elements/in² | text statement (abstract) |
| Shear-coaxial LOX post tip thickness range | 0.006 to 0.0175 in (0.15 to 0.44 mm) | text statement (abstract) |

**Derived here**, converting element density to an equivalent hexagonal pitch (pitch =
√(1/(density × sin 60°))): 1.1/in² → pitch ≈ 1.02 in = 26.0 mm; 1.6/in² → pitch ≈ 0.85 in = 21.6 mm;
2.3/in² → pitch ≈ 0.71 in = 18.0 mm. These are this agent's own arithmetic, not the source's.

## 7. Vulcain (Ariane 5 main stage, 566-element coaxial injector)

Source: a search-summary of a secondary figure caption (ResearchGate figure "Injection plate of the
Vulcain 2 LRE, composed of 566 coaxial injectors"), not independently verified against a primary
drawing or report.

| Parameter | Value | Label |
|---|---|---|
| Full-scale injector face diameter | 415 mm | WS, secondary source, unverified |
| Element count, Vulcain 2 | 566, in 13 rows | WS, secondary source, unverified |
| Element count, original Vulcain 1 | 516 | WS, secondary source, unverified |
| LOX post diameter | not found | not found |

## 8. Not found (searched, no usable number recovered)

- **J-2 (614 elements):** oxidizer-post diameter, annulus gap, recess — not found in any source
  reached in this search window. J-2's 614-element count and its "hollow oxidizer posts... fuel
  nozzles threaded through... in concentric rings" construction are well attested, but no primary
  dimension for the post itself was located.
- **LE-5 / LE-5A / LE-5B, HM7B, Vinci:** no post diameter, annulus, or recess figures found in
  English-language open sources reached in this search; likely exist only in JAXA/MHI or
  ArianeGroup/Snecma internal documentation.
- **RD-0120:** one search summary claimed "15 individual injection elements... 10 at the outer
  circle and 5 at the inner circle," which is implausibly low for an engine of this thrust class
  and could not be corroborated or traced to a primary source in this window — **treated as
  unreliable and not used.**
- **Mayer and Tamura (1996), J. Propulsion and Power 12(6):1137–1147:** a search summary returned
  "LOX post ID 5.6 mm, wall thickness 0.3 mm, recess 0.67×d_i" — these numbers are suspiciously
  identical to the Mascotte A-10 figures already in the repo (§1 above) rather than an independent
  measurement, and could not be confirmed against the primary paper (paywalled; not reached).
  **Not used — flagged as likely a search-engine conflation with Mascotte, not a genuine Mayer &
  Tamura figure.**
- **Penn State GO2/GH2 single-element combustor (Santoro/Pal/Woodward and successors):** confirmed
  to be a real, well-published single-element shear-coaxial GO2/GH2 rig (1.5 in / 38 mm chamber
  bore), but no LOX/GH2-post diameter, annulus gap or recess number was recovered from the search
  results reached; would need direct retrieval of e.g. Marshall & Pal's or Smith's papers.
  **Not used.**
- **LE-5/Vinci/HM7B velocity ratio or momentum-flux-ratio design points:** none found.
- **A NASA CR (CR-187037, Zaller, "LOX/Hydrogen Coaxial Injector Atomization Test Program," NASA
  Lewis) was found and its full text retrieved**, and it does contain an explicit "LOX Post
  Diameter" / "Gas Annulus Diameter" table for its own single-element test hardware — **but its
  running text states the test element "is smaller than SSME main chamber injector as RL-10
  injector elements," i.e. it uses RL10 as an explicit size comparator. Per the blind rule, this
  entire source was discarded unread beyond that one sentence and no number from its geometry
  table was extracted or used, out of caution that the comparison could be embedded in the same
  table.** This is a second, separate instance from the SP-8089 table noted above; both are
  disclosed here rather than silently worked around.

## 9. SP-8089 (Gill and Nurick, 1976) design-practice ranges — read directly from OCR'd primary text

Source: NASA SP-8089, "Liquid Rocket Engine Injectors" (Gill, G. S.; Nurick, W. H., Rockwell
International, March 1976). The copy found online is a 130-page scanned PDF with no text layer;
pages were rendered and OCR'd locally (`tesseract`) for this search. Section numbers and page
numbers below are as printed in the source, read from the OCR output. (The one incidental RL10
sighting while grepping this same OCR'd text is disclosed in the Blindness note above and excluded
from everything below.)

- **§3.1.1.1.4, "Impingement Distance" (p.80), normative design criterion:** "For the
  concentric-tube injector, the center-post recess (equivalent to impingement distance) should be
  no greater than one post diameter." — *text statement.*
- **State-of-the-art discussion (p.28) behind that criterion:** "data have shown that recess depth
  equal to one center-post diameter results in a maximum value of mixing and minimum drop sizes
  (refs. 5 and 14)." Ref. 14 in SP-8089's own bibliography is Burick, R. J., NASA CR-120936 (1972);
  ref. 5 is McHale, R. and Nurick, W., R-9271, Rocketdyne (1974). "Decreasing or increasing the
  recess from this value reduces mixing uniformity and increases drop size." — *text statement.*
- **§3.1.1.1.3, "Impingement Angle" (p.79–80), normative design criterion for the concentric-tube
  injector's inner-post chamfer angle:** "should not be greater than 5°." — *text statement.*
- **State-of-the-art discussion (p.26):** "Chamfer angles of about 5° generally have been used;
  however, angles up to 90° have been evaluated." — *text statement.*
- **§3.1.1.1.2, "Orifice Diameter and Diameter Ratio" (p.77), on the fuel-annulus-gap-to-post-
  diameter ratio:** no fixed numeric range is given; the normative text states only that "the
  annulus-to-center-jet ratio should be determined from manufacturing limitations, performance, and
  pressure-drop requirements," with cold-flow characterization recommended per element (refs. 5 and
  14 again). The state-of-the-art text (p.25) adds, qualitatively only: "mixing decreases when the
  ratio of annulus width to center-jet diameter increases." — *text statement*, no number.
- **§3.2.1.5, "Posts" (concentricity/cracking):** a qualitative requirement only ("shall remain
  concentric within acceptable limits and shall resist cracking"), no numeric tolerance given.

So SP-8089's own quantitative design-practice ranges, read directly, are: **recess ≤ 1 × post
diameter** (optimum ≈ 1×) and **inner-post chamfer angle ≤ 5°** (typical practice ≈5°, up to 90°
evaluated in the literature it surveys). It gives no fixed numeric range for the annulus-gap-to-
post-diameter ratio, deferring instead to Falk and Burick's cold-flow data (§5 above), which is the
same pair of authors and the same underlying experimental programme.

## 10. Summary table — what bears on the missing RL10A-3-3A dimensions

| RL10A-3-3A missing dimension | Best analog evidence | Range found |
|---|---|---|
| LOX post recess (or protrusion) | Falk & Burick Table II + worked example; SP-8089 §3.1.1.1.4; DLR BKN (table); DLR BKC (text) | 0 (BKC, flush/protruding 3 mm) to 1 × d_LOX (Falk & Burick optimum, SP-8089 ceiling); BKN sits at 0.73 × d_LOX |
| Fuel-annulus gap / d_LOX | Falk & Burick Table II and example; DLR BKN; Mascotte A-10; SSME (rough, derived) | ≈0.09 (BKN) to ≈0.80 (SSME, derived/uncertain); Falk & Burick's own range table spans h = 0.005–0.041 in at D_L = 0.070–0.136 in, i.e. gap/D_L ≈ 0.04 to 0.59; worked example 0.21 |
| LOX post tip wall thickness / d_LOX | Falk & Burick example (0.20); SSME present design (0.37, derived); J-2X/Ares abstract (absolute thickness 0.15–0.44 mm, no ID given so no ratio) | ≈0.20 to ≈0.37 |
| Element pitch / annulus OD | SSME present/proposed (pitch only, annulus OD not separately known); J-2X/Ares (density→pitch, derived here); Mascotte, DLR (single-element rigs, no pitch) | not resolvable to a clean ratio from what was found; absolute pitches found range 16.0 mm (SSME present) to 50.8 mm (SSME proposed) to 18–26 mm (Ares/J-2X, derived from element density) |
| Chamfer/taper angle at the post tip | SP-8089 (≤5°, up to 90° in the literature it cites); DLR BKN (8°, outward taper) | about 5° to 8°, with one outlier case evaluated up to 90° |
