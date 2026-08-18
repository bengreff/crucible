# VAL-2 — Anchor Specifications

| Field | Value |
|---|---|
| **ID** | VAL-2 |
| **Family** | VAL (Validation & test) |
| **Status** | Reviewed (2026-08-14) — grows as anchors are added per wave |
| **Depends on** | VAL-1, META-3 (anchor data), FND-1 |
| **Version** | 0.2 (2026-08-14 review fixes: N16, N17, N18/D-G, N19, D-C) |

---

## 0. Purpose

VAL-2 is the **per-anchor specification**: for each validation anchor, the sourced **data-of-record**, the
**target/criterion**, and the **retrieval**. It instantiates VAL-1's ladder against concrete cases. W2 fixes
the **RL10A-3-3A system anchor** and the **unit-physics anchors** (EOS Hugoniot, compressible-Euler, radiation
analytics, stopping) that validate the constituent laws underneath it; the **nuclear anchors** (KRUSTY,
NERVA — §3.4) are specified now, ahead of the nuclear coding wave (N16).

Read after VAL-1 (the framework) and META-3 §5/§6.6 (the sourced data).

## 1. Scope & razor ruling
**Owns:** the anchor specs (data + target + retrieval) and the **pass criterion**. **Defers:** the framework
(GCI/area-metric/p-box) → **VAL-1**; the test harness that runs them → **VAL-3**; the source data ledger →
**META-3**.

**Razor ruling:** infrastructure. Its obligation is that every anchor's data is **sourced and cached**
(META-1 Principle 7) and its criterion is **honest** (overlap-band, not a false cliff — §3.1).

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Anchor spec** | VAL-3 (runs it), COUP-6 (pedigree) | data-of-record + target + retrieval, per anchor |
| **Pass criterion** (§3.1) | VAL-3 | overlap of the predicted p-box with the reference p-box + reported area metric |

**Invariant:** every anchor cites a META-3 `key` with cached source; no anchor uses an un-sourced number.

## 3. Method

### 3.1 The pass criterion — overlap band, computable *(Ben, 2026-07-21; N17)*
The reference is a **p-box**, not a point, and the criterion is a computation, not a judgment call:

- **Reference p-box construction:** the published reference interval per QoI (Isp 440–446 s, thrust
  16.41–16.50 klbf) **convolved with TM-107318's own published prediction-vs-measurement error
  distribution** — the report's chamber sub-model predicts Isp 440.3 s where the cycle model gives 445.6 s,
  a documented model-scatter that *is* the anchor's error distribution. Each interval endpoint is convolved
  with that distribution, giving lower/upper reference CDFs (F̲_ref, F̄_ref). The anchor is itself a
  validated model with published error — CRUCIBLE inherits the numbers *and* their scatter.
- **Overlap metric (Ferson p-box generalization of the area metric):** with (F̲, F̄) the lower/upper CDF
  bounds of each p-box, `d = ∫ max(0, F̲_pred(q) − F̄_ref(q), F̲_ref(q) − F̄_pred(q)) dq` over the QoI axis
  — the area by which the two boxes *fail* to overlap; `d = 0` iff they intersect at every quantile.
- **Pass = overlap** (`d = 0`); **`d` is reported as the score either way** — a miss is a quantified,
  honestly-reported model-form number, never a silent fail. The **≤2% is the reported target, not a binary
  gate** (VISION_SCOPE S1, reworded per the 2026-07-21 amendment).

Report **blind** and **calibrated** per the mode definitions in §3.2. [META-3: `area-metric`, `pbox`]

### 3.2 RL10A-3-3A — system anchor
Data-of-record (NASA TM-107318 / NTRS 19970010379; cached): expander cycle, LOX/LH₂; **Pc 475 psia, O/F 5.0,
thrust 16,500 lbf (73.4 kN), Isp 440.3 s (chamber sub-model) / 445.6 s (cycle), c\* 7824 ft/s (η_c\* 0.989),
ε 61, throat Ø 2.47 in**, 216-element coaxial injector; fuel pump 6.05 lb/s→~1100 psia @31,537 rpm η≈0.58,
LOX pump 30.8 lb/s→600 psia @12,615 rpm η≈0.64; jacket heat 7994 Btu/s, ΔP 242 psid; full 16-station cycle
table (Table 6.1.1). **Target:** predicted p-box (thrust, Isp, c\*, C_F, emergent Pc) overlaps the reference
band per §3.1; `d` reported; c\* and C_F validated **separately** (the `v_e = c*·C_F` split).
[META-3: `rl10-cycle-data`, `cstar-cf-defs`]

**Geometry-of-record (N19):** the blind config is authored from the **RL10 chamber/nozzle contour tables**
(1966 design report contour stations; TM-107318 station geometry — throat Ø 2.47 in, ε 61, contraction
ratio, injector-end plane), retrieved and **cached per the archival rule** (META-1 Principle 7). Without
this the blind config cannot be authored at all. [META-3: `rl10-geometry`]

**Blind vs calibrated — the input split (N18/D-G, per VAL-1 §3.6):**
- **Blind (open-mode feed):** design spec + geometry-of-record + universal physics/closures + **technology-
  class data from *other* hardware** — specifically the **coaxial-injector-family η_c\* prior band**
  (±1–3%, cited), pump-class efficiency envelopes. **No quantity measured on RL10 itself.**
- **Calibrated (closed-mode expander):** may additionally bind **RL10's own fitted η_c\* = 0.9892 and its
  measured pump maps** (TM-107318 component data). These anchor-derived values are **calibrated-mode only**
  and every reported score is **labeled** blind or calibrated — the blind headline never silently consumes
  the anchor's own measurements (the prior circularity this ruling retires).

**Bartz's role (D-C):** Bartz is **not** a runtime closure anywhere — the runtime wall flux is SOLV-1 §3.5's
one local wall-function law. Here Bartz serves as the **nozzle-envelope validation oracle**: the predicted
integrated nozzle heat load is cross-checked against the Bartz correlation ± its published band as an anchor
diagnostic. [META-3: `bartz` (oracle)]

### 3.3 Unit-physics anchors (constituent laws)
| Anchor | Validates | Data / target | Key |
|---|---|---|---|
| **Sod shock tube** | compressible-Euler solver (SOLV-1) | exact Riemann; shock within 1 cell, L1→0 at formal order | `sod-shock` |
| **Su-Olson / Marshak** | radiation transport (SOLV-2) | analytic T_rad,T_mat; within ~1–3% | `su-olson` |
| **Principal Hugoniot** | EOS (FND-7/OFFL-5) | Al 0.3–12 Mbar laser-shock; within experimental bars | `hugoniot-anchor` |
| **NIST PSTAR/ASTAR** | stopping (FND-7/SOLV-3) | proton/alpha dE/dx + range; within ~1–5% (Bethe regime) | `stopping-astar` |
| **RP-1311 / CEA cases** | equilibrium chemistry (OFFL-3) | CEA example cases within manual tolerances | `nasa-cea` |

### 3.4 Nuclear anchors *(N16 — written now; milestone-1 scope per R1)*

**Criterion doctrine first:** S2's **≤300 pcm and ≤10% are reported targets over a p-box, not binary
gates** — mirroring S1's overlap-band doctrine (§3.1). Each nuclear QoI is compared as predicted-p-box vs
reference band (measurement uncertainty + report scatter), `d` reported, the S2 numbers quoted alongside as
targets.

**KRUSTY — criticality + coupled-transient anchor.**
- **Configuration of record:** the ICSBEP **HEU-MET-FAST-101** benchmark configuration (KRUSTY cold-critical
  geometry/materials as evaluated); model retrieval from the **MOOSE-VTB KRUSTY model** (open, incl.
  geometry + materials), cached per the archival rule.
- **Per-quantity criteria (reported targets):** k-eff within **≤300 pcm** of the benchmark value; **warm
  criticals** (temperature-defect reactivity points) reproduced within the benchmark's stated bands; the
  **coupled full-power transient** reproduces the published power/temperature **shape** (feedback
  self-regulation), compared as time-series band overlap, not point equality. [META-3: `krusty`]

**NERVA — hot-fire state-point anchors** (state-point definitions from the cited reports; data retrieved
and cached): [META-3: `nerva-pewee`, `nrx-a6`, `xe-prime`]
- **Pewee** (NASA-CR-184270, NTRS 19920005899; LA-4217-MS): >500 MW; fuel-exit 2556 K; chamber 1833 K at
  18.6 kg/s, 4275 kPa; peak ideal vac Isp 901 s. Per-quantity target **≤10%** on chamber T, flow state,
  Isp.
- **NRX-A6** (ibid.; WANL-TNR-223/224): 1120 MW × 60 min endurance state point; Pc 4089 kPa at 32.7 kg/s.
  Target **≤10%** on the held state point.
- **XE-Prime** (ibid.; Aerojet RN-S-510): 1140 MW; 244.75 kN thrust-stand thrust; Isp ~710 s; 2272 K;
  35.8 kg/s; plus the **operating map** — target **≤10%** on state points **and qualitative map-shape
  agreement** (throttle-line trends), the E-2 control-state axis being what makes the map traversable.
- Each state point is a **declared (power, ṁ, Pc, T) tuple with report provenance** — the anchor spec is
  the tuple + citation + cached retrieval, so a run's comparison set is data, not prose.

*(NSTAR, VASIMR, antimatter, BEAVRS anchors are added in the plasma/pulsed waves — already seeded in
META-3 §5.)*

## 4. Coupling relationships
- **VAL-1** provides the p-box/area-metric machinery the criterion uses; **VAL-3** runs each anchor as a test/
  CI gate; **COUP-6** folds anchor results into pedigree; the consuming modules (SOLV-1/2, FND-7, OFFL-3/5,
  SOLV-7) are validated here.
- **META-3** holds the cached source data (RL10 primary source cached locally per the archival rule).

## 5. Uncertainty & validity
The anchor *reference* itself carries uncertainty (RL10's model/measurement scatter; experimental Hugoniot
bars), which is **part of the criterion** (overlap of p-boxes) — not ignored. Validity: each anchor validates
its module only within the anchor's regime; extrapolation beyond is governed by VAL-1's widening model-form
band. Ladder: unit (Sod/Su-Olson/Hugoniot/stopping) + system (RL10).

## 6. Validation plan
This doc *is* the validation plan for W2. VAL-3 executes it: unit anchors on every commit (cheap/analytic),
the RL10 system anchor at milestones (blind then calibrated), each reported as a p-box overlap + area metric.

## 7. References
META-3 keys: `rl10-cycle-data`, `rl10-tm107318`, `rl10-geometry` *(new)*, `cstar-cf-defs`, `sod-shock`,
`su-olson`, `hugoniot-anchor`, `stopping-astar`, `nasa-cea`, `area-metric`, `pbox`, `krusty`, `nerva-pewee`,
`nrx-a6`, `xe-prime`, `bartz` (oracle). Depends on VAL-1 (framework + §3.6 blind rule), FND-1 (p-box).

*(No open questions — pass criterion (overlap band, no hard cutoff; blind + calibrated) resolved by Ben,
2026-07-21; S1 reworded in VISION_SCOPE per the 2026-07-21 amendment; blind-mode definition ruled D-G,
2026-08-14.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-17 | 0.2.1 | Transcription erratum: §3.2 data-of-record c\* unit corrected `in/s` → `ft/s` (TM-107318 reports 7824 ft/s ≈ 2385 m/s; the in/s value would be 12× low). No criterion or method change. Same fix in META-3 `rl10-cycle-data`. |
| 2026-08-14 | 0.2 | Review fixes (N16, N17, N18/D-G, N19, D-C). **N17:** S1 criterion made computable — reference p-box = reference interval convolved with TM-107318's published prediction-vs-measurement error distribution (Isp 440.3 chamber sub-model vs 445.6 cycle); overlap metric `d = ∫max(0, F̲_pred−F̄_ref, F̲_ref−F̄_pred)dq`; pass = overlap (d=0), d reported as the score (§3.1). **N18/D-G:** blind = coax-injector-class prior band + design spec + geometry-of-record, nothing measured on RL10; fitted η_c\* = 0.9892 + measured pump maps = calibrated mode only, labeled (§3.2). **N19:** RL10 geometry-of-record retrieval added (1966 design-report contours / TM-107318 stations, cached). **D-C:** Bartz recorded as nozzle-envelope validation oracle only (runtime = SOLV-1 §3.5 wall-function law). **N16:** nuclear anchor specs written — KRUSTY (HEU-MET-FAST-101 config, MOOSE-VTB retrieval, ≤300 pcm + warm criticals + coupled-transient shape) and NERVA Pewee/NRX-A6/XE-Prime (state-point tuples from NASA-CR-184270, LA-4217-MS, WANL-TNR-223/224, RN-S-510; per-quantity ≤10% + map shape); S2 thresholds stated as reported targets over a p-box, mirroring S1 (§3.4). |
| 2026-07-21 | 0.1 | Initial draft (W2 anchors). RL10A-3-3A system anchor with full sourced cycle data (TM-107318); overlap-band pass criterion + Ferson area metric, blind (open-mode) and calibrated (closed-mode), no hard cutoff (Ben 2026-07-21; S1 "≤2%" flagged as a reported target, not a gate); c\*/C_F separate validation; unit-physics anchors (Sod, Su-Olson, Hugoniot, PSTAR/ASTAR, RP-1311). Nuclear/plasma/pulsed anchors deferred to later waves. |
