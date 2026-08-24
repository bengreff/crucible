# CRUCIBLE Station 5 Certificate — the Blind RL10

Goal-B station 5 (VAL-2 §3.1/§3.2; SOLV-7; COUP-3 §3.5; COUP-7 §3.2): the RL10A-3-3A exists in this repository ONLY as data — the geometry-of-record contour CSV, the preset TOMLs, and the pinned `lox_lh2_v0.4.0` equilibrium surface — assembled by the one config-driven engine path and marched to a settled state at the coarse tier (`cells_across_throat = 5`, 41×119, 2667 gas cells, 12 flow-throughs, ~155 s laptop wall clock per member on the S2 SDC-IMEX spine — implicit class-D liner conduction + Robin-Robin wall exchange inside the step, the COUP-2 conservation audit armed on every one of the ~15,800 steps of every member, zero violations). Every reported number is a plane integral of the conserved field (SOLV-7); chamber pressure and thrust are **emergent, never imposed**. **Spine of record for the scores below: S2/S3** — see KNOWN LIMITS for the measured S4 delta and why the scores are not re-derived here.

**Labels (VAL-2 §3.2, VISION_SCOPE §9 v1.4.1):** every score below is labeled blind or calibrated, and every score carries `development-observed: yes` — the RL10 campaign is declared OPEN development (anchor deltas were watched while the instrument was built); blind remains a mechanical INPUT property: the blind members consume the design spec, the geometry of record, universal closures, and the coax-injector-family η_c\* prior band measured on other hardware — no quantity measured on the RL10 itself. The calibrated members additionally bind the anchor's own fitted η_c\* = 0.9892 and TM-107318 component data.

## The reference p-box (VAL-2 §3.1, N17)

Published interval per QoI convolved with the anchor's own documented model scatter (chamber sub-model Isp 440.3 s vs cycle 445.6 s ⇒ relative half-width 0.598% applied per QoI): Isp [440.3, 445.6] s; F [16.41, 16.50] klbf = [73.00, 73.40] kN; p_c [475 (cycle station), 482 (Table 2.5.1 injector-face)] psia = [3.2750, 3.3233] MPa; c\* 7824 ft/s = 2384.8 m/s; C_F derived F/(p_c·A_t) ∈ [1.7763, 1.8124] with A_t = π·(2.47 in)² = 0.012365 m² (the radius erratum, VAL-2 0.2.3). Pass = overlap (d = 0); d is reported either way — a miss is a quantified model-form number, never a silent fail.

## BLIND score — open mode (label: blind, development-observed: yes)

Members: the band-mid nominal (η_c\* target 0.98, realized 0.9786 = 2208.8/2257.0 against the full-equilibrium baseline) and the four declared-band corners — coax-family η_c\* edge (realized 0.966–0.969 and 0.988–0.991; S18 source-level knockdown, coarse-calibrated slope −0.847% c\* per −3×10⁵ J/kg) × wall-law band edge (the SOLV-1 §3.5 ±20–30% Colburn band, realized as h × 0.75 / × 1.25). Boxes widened by the declared ±0.5% numeric band (limit-cycle + flow-closure + integration class).

| QoI | predicted box | reference (pre-convolution) | Ferson d | d / ref-mid |
|---|---|---|---|---|
| thrust F | [70.5887, 73.7505] kN | [72.9953, 73.3957] kN | **0 (overlap)** | — |
| Isp | [424.7456, 443.7276] s | [440.3000, 445.6000] s | **0 (overlap)** | — |
| c* | [2169.7965, 2247.6825] m/s | [2384.7552, 2384.7552] m/s | 137.0735 m/s | 5.75% |
| C_F | [1.9101, 1.9457]  | [1.7763, 1.8124]  | 0.0977  | 5.45% |
| p_c | [2.9737, 3.0807] MPa | [3.2750, 3.3233] MPa | 0.1943 MPa | 5.89% |

F and Isp OVERLAP the reference at the band's upper edge; p_c, c\*, and C_F miss coherently (p_c and c\* ~6% low, C_F correspondingly high — one discretization signature, not three physics errors: at the coarse tier the under-resolved throat/wall region under-produces chamber pressure while the exit momentum integral is nearly converged, and c\* = p_c·A_t/ṁ inherits p_c's deficit while C_F = F/(p_c·A_t) inherits its inverse). The indicative dial-8 refinement moved exactly these: p_c +4.7%, c\* +4.7%, C_F −4.6%, F +0.06% (§limits for why that tier is not certified tonight).

## CALIBRATED score — closed expander (label: calibrated, development-observed: yes)

The COUP-3 §3.5 fixed point closes the cycle each step (jacket pickup → turbine power → pump map → delivered ṁ; Aitken-relaxed, residual ≤ 1e-8, engaged after the establishment window): the engine finds its OWN operating point. Nominal: delivered ṁ = 17.71 kg/s (+4.5% over design — the wall law's +16% jacket pickup driven through the declared ṁ³ impedance line), turbine 671.7 kW, T_turbine_in 235.2 K. Members: nominal + the wall-band corners (h × 0.75 ⇒ ṁ 16.74; h × 1.25 ⇒ ṁ 18.55 kg/s — the wall-function band IS the dominant p_c spread, exactly as COUP-7 §3.4 predicted).

| QoI | predicted box | reference (pre-convolution) | Ferson d | d / ref-mid |
|---|---|---|---|---|
| thrust F | [71.9764, 80.2534] kN | [72.9953, 73.3957] kN | **0 (overlap)** | — |
| Isp | [436.8746, 443.1849] s | [440.3000, 445.6000] s | **0 (overlap)** | — |
| c* | [2219.1485, 2247.1800] m/s | [2384.7552, 2384.7552] m/s | 137.5779 m/s | 5.77% |
| C_F | [1.9209, 1.9438]  | [1.7763, 1.8124]  | 0.1086  | 6.05% |
| p_c | [3.0096, 3.3617] MPa | [3.2750, 3.3233] MPa | **0 (overlap)** | — |

**F, Isp, and the emergent p_c all OVERLAP the record** (p_c box [3.01, 3.36] MPa spans the published 3.27–3.32; the nominal alone reads 463.8 psia vs the 475–482 record). Isp is nearly flat across the whole wall band (439.1–441.0 s vs record 440.3–445.6): the closed cycle trades ṁ against p_c at almost constant specific impulse — real expander-cycle self-regulation, reproduced by the coupled instrument, not imposed. c\*/C_F carry the same coarse-tier discretization signature as the blind score.

## Consistency checks

- `v_e = c\*·C_F` closes to 4.22e-4 relative on the nominals (SOLV-7 §6-2; the residual is the ṁ_exit vs ṁ_inj closure, ≤ 0.2%).
- `p_c·A_t ≈ ṁ·c\*` closes by construction of the c\* readout (SOLV-7 N11: the identity is the harness check, never the definition).
- Delivered-flow closure: the measured inflow-plane ṁ is within 0.2% of the declared/solved value on every member (sonic startup cap inactive at readout — the session-12 honesty signal).
- Steadiness: every member's residual is a stationary limit cycle (max |Δρ|/ρ 0.6–4.3×10⁻² over the probe window — a plume-fringe cell-wise max; the integral readouts' F oscillation is ≤ ±0.15% on every member, measured ±0.08% on the noisiest), inside the declared numeric band.

## Declared bands, model form, and declared devices

- **Wall function (SOLV-1 §3.5):** the one Colburn-class law, ±20–30% declared band — realized as the bracket corners; in closed mode it dominates the p_c spread (as designed). Jacket pickup at the calibrated nominal: 9.76 MW vs the record 8.43 MW (+16%, inside the band). Bartz nozzle-envelope oracle scoring: deferred (recorded), rides the next wave.
- **η_c\* prior (COUP-7 §3.2.1):** blind = the coax-family band applied as the S18 source-level knockdown (never output-side); realized η re-measured per member and reported above. Calibrated = the anchor's fitted 0.9892.
- **Equilibrium surface v0.4.0:** gas-only METASTABLE products (declared plume model — real plumes supersaturate; deck-stamped), Z narrowed to the premixed operating class, h ∈ [−1.23×10⁷, +1.2×10⁷] J/kg with rule-space (relative) interp bounds incl. envelope-edge holdout (session-12 review fixes); the projection acceptance uses the density column's own log-space bound. v0.4.0 is the plan-S6 ignition-headroom STRICT extension of v0.3.2 (OFFL-3 0.6.2): the recorded marches ran on v0.3.2, whose every in-envelope value and bound is bit-identical in v0.4.0 (verified column-by-column at regeneration).
- **Geometry:** FND-3 analytic partial fractions + apertures (zero sampling error) with State Redistribution (κ < 0.5, Berger–Giuliani) — the session-11 stair-transpiration/starvation class is retired; wall heat runs on the closure-vector (smooth) interface area, not the stair overcount.
- **Integrator (S2):** the ONE deterministic SDC-IMEX step (COUP-3 §3.1) — explicit hyperbolic class + implicit class-D liner conduction (fixed-cycle CG) with the Robin-Robin wall exchange inside each sweep (COUP-2 §3.5) and the COUP-2 conservation audit armed every step; Δt is the gas CFL alone. The liner ρc_p remains the declared steady-state continuation device until the plan's S4 gives it the physical value (now legal under the implicit class).
- **UQ:** these boxes are declared-band corner brackets (epistemic intervals), NOT the full COUP-5 ensemble p-box — that machinery is a later wave; the boxes are never collapsed to points.

## KNOWN LIMITS (recorded, with owners)

- **Certified tier = dial 5 (coarse).** Dials ≥ 8 do not survive the ESTABLISHMENT march under the honest rule-space projection acceptance: the startup transient (drain/shear at the bell wall; injector piston; or overexpanded-bell backflow — five schedules probed, each halting loudly at a different envelope edge) manufactures mixture states genuinely outside the equilibrium surface's representable set. The pre-review dial-8 'success' rode on the vacuous absolute acceptance bound the session-12 review retired. The designated cure (VISION_SCOPE v1.5, session 13 — accelerated convergence is DELETED; establishment stays a physical march) is the **PLAN_CHEMICAL_SANDBOX Phase 1–2 physics**: gas diffusion (plan S3) + the implicit integrator (**LANDED — this S2 rerun**) + the cold/unburnt chemistry branch (plan S5–S6), which together make these transient states representable and conductively coupled instead of refused; grid-sequenced restart from settled physical states (FND-6) is the legal warm-start. The indicative dial-8 deltas quoted above are direction-of-refinement evidence only.
- **Discretization band is therefore declared, not swept:** the scored boxes carry closure bands + the numeric band; the coarse-tier truncation error is visibly ~5–6% on p_c/c\*/C_F (the coherent signature above) and is NOT hidden inside the boxes.
- `interp_error_bound_log` rides outside digest v3 (recorded digest-v4 deferral). The station-4-fixture rewire deferral is DISCHARGED (S2): the fixture and the engine now march the same one integrator.
- **The recorded readouts below predate plan S4 (session 16).** They were earned on the S2/S3 spine with the wall law's chamber-fitted constant transport and the ρc_p continuation device. S4 replaced both — the FND-7 §3.3 tabulated spine now supplies per-cell μ/k/c_p/c_v and the liner carries physical areal thermal capacitance — and the preset TOMLs cited under Reproduction are the S4 ones. The **measured** S4 shift on the η = 1 baseline member (`configs/rl10_coarse.toml`, 15 819 steps, audit clean): thrust +0.30%, Isp +0.25%, c\* −0.07%, C_F +0.32%, p_c −0.02%, and **jacket heat −19.8%** (liner T_max 410.1 K) — every scored quantity moves an order of magnitude less than the declared bands the boxes carry, while the wall term (which no scored quantity reports, but which drives the closed cycle) moves by two orders of magnitude more. That last number is why the CALIBRATED members cannot simply be relabeled: the closed cycle converts jacket pickup into delivered ṁ, so a −19.8% wall term is a first-order change to those scores and **only** to those scores. The scores are therefore NOT re-derived here: re-scoring all eight members on the S4 spine is ~5 h of laptop march and it rides the **COUP-5 ensemble wave (plan S18/S19)**, which replaces these hand-run corner brackets with sampled p-boxes anyway. Two consequences are recorded rather than papered over: the wall-law band corners are now the declared `band_factor` (h × 0.75/1.25 directly) instead of the `cp_j_per_kg_k` proxy that reached the same h; and the physical liner clock (~37 ms) exceeds the 12-flow-through settle budget (~11 ms), so an S4 re-scoring must also grow that budget — run length, which plan ruling #4 accepted, not schedule tuning.

## Reproduction

`cargo run --release -p crucible -- run configs/rl10_coarse.toml` (η = 1 baseline), `configs/rl10_calibrated.toml` (calibrated nominal). Bracket members = the coarse preset + the overrides stated in the recorded-readout comments of `station5_rl10_certificate.rs` (note: those comments name the S2/S3-spine overrides; S4 reaches the wall-law band corners through `[mechanisms.wall] band_factor` instead). Table regen: `offline/scripts/make_station5_tables.py` (~2 min). Every run writes its resolved-config manifest + fields CSV under `runs/<name>/`; halts write `crash_fields.csv` (raw conserved state — the session-12 instrument that diagnosed every establishment pathology above).

## Sensitivity (uncertified indication): the blind score under the dial-8 refinement deltas

Applying the indicative deltas (p_c +4.7%, c* +4.7%, C_F -4.6%, F +0.1%, Isp +0.1%) to the blind box edges: 
- p_c: d = 0.0495 MPa (1.50% of ref-mid)
- c*: d = 31.4299 m/s (1.32% of ref-mid)
- C_F: d = 0.0099  (0.55% of ref-mid)
- thrust F: d = 0 (overlap)
- Isp: d = 0 (overlap)

The residual blind miss under refinement is the wall-law + η-prior model form the bands already carry — the certified sweep that would promote this from indication to score awaits the PLAN_CHEMICAL_SANDBOX Phase 1–2 wave (physical-march establishment with the full diffusion + cold-branch physics).

