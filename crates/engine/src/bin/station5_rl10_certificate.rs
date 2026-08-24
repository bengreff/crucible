//! Regenerates `certificates/station5_rl10_certificate.md` — GOAL-B
//! station 5: the RL10A-3-3A assembly scored against the TM-107318
//! reference p-box (VAL-2 §3.1/§3.2).
//!
//! Unlike stations 1–4, the field marches behind this certificate cost
//! minutes each, so gate 5 cannot rerun them. The house rule extends the
//! same way META-3 handles anchor data: the run READOUTS are recorded
//! constants below, each traceable to a config (the committed presets, or
//! the documented bracket variants — base preset + the stated overrides)
//! and reproducible with `crucible run <config>` — the resolved-config
//! manifest + table pins make every march a pure function of its inputs.
//! This bin recomputes the SCORING (reference-box construction, predicted
//! boxes, Ferson overlap metric) deterministically from those recorded
//! inputs; gate 5 diffs the artifact.

use std::fmt::Write as _;

// --- Recorded run readouts (SOLV-7 performance objects) --------------------
// Session 14 (2026-08-19, plan S2 — re-run on the SDC-IMEX spine: implicit
// class-D liner + Robin-Robin exchange + the COUP-2 audit armed every step;
// readouts moved < 0.1% from the session-12 records), table lox_lh2_v0.3.2
// (pins sidecar), dial 5 (cells_across_throat), 12 flow-throughs, the
// committed engine tree. Bracket variants = the named base preset + the
// stated overrides only. Plan-S6 note: the presets now pin lox_lh2_v0.4.0,
// a STRICT h-ceiling extension of v0.3.2 (old 43 h-nodes bit-exact, nodes
// appended above; verified all-columns bit-identical in the old region and
// every interp bound unchanged — OFFL-3 0.6.2), so these recorded v0.3.2
// readouts are exact for v0.4.0 marches too.

#[derive(Clone, Copy)]
#[allow(dead_code)] // full recorded readout of record; the scoring consumes a subset
struct Readout {
    f: f64,        // N
    isp: f64,      // s
    c_star: f64,   // m/s
    c_f: f64,      // 1
    p_c: f64,      // Pa
    mdot_inj: f64, // kg/s (delivered/declared)
    jacket: f64,   // W
    resid: f64,    // steadiness
}

/// configs/rl10_coarse.toml as committed (η = 1 full equilibrium — the
/// sweep/calibration baseline, NOT a scored member).
const R_COARSE_ETA1: Readout = Readout {
    f: 74162.9,
    isp: 446.25,
    c_star: 2257.0,
    c_f: 1.9389,
    p_c: 3.0933e6,
    mdot_inj: 16.947,
    jacket: 9.737e6,
    resid: 6.16e-3,
};

/// BLIND nominal: coarse preset + eta_cstar 0.98 / h_offset −7.16e5
/// (the coax-family prior band mid; realized η_c* = 2208.8/2257.0 = 0.9786).
const R_BLIND_NOM: Readout = Readout {
    f: 72179.0,
    isp: 434.21,
    c_star: 2208.8,
    c_f: 1.9277,
    p_c: 3.0280e6,
    mdot_inj: 16.947,
    jacket: 9.246e6,
    resid: 9.50e-3,
};

/// BLIND corners: η edge (offset −1.079e6 → realized 0.966..0.969, or
/// −3.56e5 → realized 0.988..0.991) × wall-law band edge.
///
/// The band edge was reached, on the S2/S3 spine these readouts come from, by
/// setting the wall law's then-private `cp_j_per_kg_k` to 3750/6250 at fixed
/// Pr and μ — which scales k = μc_p/Pr with c_p and so scales h by
/// 0.75/1.25 exactly. S4 retired that block (transport is a spine query) and
/// named the band directly: `[mechanisms.wall] band_factor = 0.75 / 1.25`.
/// The two are not identical — the proxy also moved the recovery temperature
/// and, after S3, the resolved viscous fluxes — which is one reason an S4
/// re-scoring is a re-run, not a relabel (see KNOWN LIMITS).
const R_B_E97_WLO: Readout = Readout {
    f: 71317.0,
    isp: 429.12,
    c_star: 2186.9,
    c_f: 1.9243,
    p_c: 2.9972e6,
    mdot_inj: 16.947,
    jacket: 7.741e6,
    resid: 2.31e-2,
};
const R_B_E97_WHI: Readout = Readout {
    f: 70943.4,
    isp: 426.88,
    c_star: 2180.7,
    c_f: 1.9197,
    p_c: 2.9886e6,
    mdot_inj: 16.947,
    jacket: 10.158e6,
    resid: 1.27e-2,
};
const R_B_E99_WLO: Readout = Readout {
    f: 73383.6,
    isp: 441.52,
    c_star: 2236.5,
    c_f: 1.9360,
    p_c: 3.0654e6,
    mdot_inj: 16.947,
    jacket: 8.162e6,
    resid: 4.30e-2,
};
const R_B_E99_WHI: Readout = Readout {
    f: 73020.9,
    isp: 439.26,
    c_star: 2230.4,
    c_f: 1.9314,
    p_c: 3.0575e6,
    mdot_inj: 16.947,
    jacket: 10.726e6,
    resid: 6.81e-3,
};

/// CALIBRATED nominal: configs/rl10_calibrated.toml as committed (closed
/// expander, TM component data, fitted η_c* = 0.9892 via offset −3.85e5).
const R_CAL_NOM: Readout = Readout {
    f: 76382.3,
    isp: 439.91,
    c_star: 2233.4,
    c_f: 1.9316,
    p_c: 3.1978e6,
    mdot_inj: 17.713,
    jacket: 9.759e6,
    resid: 2.49e-2,
};

/// CALIBRATED wall-band corners (S2/S3 spine, cp 3750/6250 = the h × 0.75/1.25
/// proxy — S4's direct form is `band_factor`, see the blind corners above):
/// the wall law drives the cycle, so the ±25% h band sweeps delivered
/// ṁ 16.74 → 18.56 kg/s.
const R_C_WLO: Readout = Readout {
    f: 72338.1,
    isp: 440.98,
    c_star: 2236.0,
    c_f: 1.9341,
    p_c: 3.0247e6,
    mdot_inj: 16.737,
    jacket: 8.077e6,
    resid: 3.99e-2,
};
const R_C_WHI: Readout = Readout {
    f: 79854.1,
    isp: 439.07,
    c_star: 2230.3,
    c_f: 1.9306,
    p_c: 3.3450e6,
    mdot_inj: 18.546,
    jacket: 11.357e6,
    resid: 6.12e-3,
};

/// Declared numeric half-band applied to every predicted box edge:
/// steadiness limit-cycle (probe F oscillation ≤ ±0.15%), mass-flow
/// closure (inflow-plane vs declared ≤ 0.2%), plane-integration class.
const NUMERIC_REL: f64 = 0.005;

/// INDICATIVE (not certified) refinement deltas, dial 5 → 8, measured
/// 2026-08-18 UNDER THE PRE-REVIEW acceptance gate (the vacuous absolute
/// bound) and zero-Roe wavespeeds — the establishment at dial 8 does not
/// survive the honest gate (see §limits), so these are direction-of-
/// refinement indications only, never folded into the scored boxes:
/// p_c +4.7%, c* +4.7%, C_F −4.6%, F +0.06%, Isp +0.14%.
const INDICATIVE_REFINEMENT: [(&str, f64); 5] = [
    ("p_c", 0.047),
    ("c*", 0.047),
    ("C_F", -0.046),
    ("F", 0.0006),
    ("Isp", 0.0014),
];

// --- Reference data of record (VAL-2 §3.2; META-3 rl10-cycle-data) ---------

/// Throat area πr_t², r_t = 2.47 in (the VAL-2 0.2.3 radius erratum).
const A_T: f64 = std::f64::consts::PI * (2.47 * 0.0254) * (2.47 * 0.0254);

/// Published reference intervals per QoI. Isp: chamber sub-model 440.3 s
/// vs cycle 445.6 s (TM's documented model scatter — N17 makes it the
/// anchor's error distribution). Thrust: 16.41–16.50 klbf (VAL-2 §3.1).
/// p_c: the cycle-table station 475 psia and Table 2.5.1's injector-face
/// 482 psia — the anchor's own two published chamber pressures. c*:
/// 7824 ft/s (the 0.2.1 unit erratum). C_F is DERIVED (no direct published
/// value): F/(p_c·A_t) over the published extremes.
const REF_ISP: (f64, f64) = (440.3, 445.6);
const REF_F: (f64, f64) = (16410.0 * 4.448_221_6, 16500.0 * 4.448_221_6);
const REF_PC: (f64, f64) = (475.0 * 6894.757, 482.0 * 6894.757);
const REF_CSTAR: (f64, f64) = (7824.0 * 0.3048, 7824.0 * 0.3048);

/// TM's model-scatter half-width, relative: half the chamber-vs-cycle Isp
/// spread over its midpoint (N17 — the report's documented prediction
/// scatter, applied per QoI as the convolution kernel half-width).
const SCATTER_REL: f64 = (445.6 - 440.3) / 2.0 / ((445.6 + 440.3) / 2.0);

// --- The Ferson overlap metric (VAL-2 §3.1) --------------------------------

/// Reference p-box: interval [a, b] convolved with uniform ±w·q scatter →
/// upper CDF ramps over [a(1−w), a(1+w)], lower CDF over [b(1−w), b(1+w)].
/// Predicted p-box: sharp interval [pa, pb] (an epistemic interval — the
/// declared-band corners; no distribution is claimed inside it).
/// d = ∫ max(0, F̲_pred − F̄_ref, F̲_ref − F̄_pred) dq on a fixed
/// 20001-point grid (deterministic quadrature of a piecewise-linear
/// integrand — error far below reporting precision).
fn ferson_d(pred: (f64, f64), refi: (f64, f64), w: f64) -> f64 {
    let (pa, pb) = pred;
    let (a, b) = refi;
    let (ra0, ra1) = (a * (1.0 - w), a * (1.0 + w));
    let (rb0, rb1) = (b * (1.0 - w), b * (1.0 + w));
    let lo = pa.min(ra0) - 1e-9 * pa.abs();
    let hi = pb.max(rb1) + 1e-9 * pb.abs();
    let n = 20_001usize;
    let dq = (hi - lo) / (n as f64 - 1.0);
    let ramp = |q: f64, x0: f64, x1: f64| -> f64 {
        if q <= x0 {
            0.0
        } else if q >= x1 {
            1.0
        } else {
            (q - x0) / (x1 - x0)
        }
    };
    let mut acc = 0.0f64;
    for k in 0..n {
        let q = lo + k as f64 * dq;
        let f_pred_lo = if q >= pb { 1.0 } else { 0.0 };
        let f_pred_hi = if q >= pa { 1.0 } else { 0.0 };
        let f_ref_hi = ramp(q, ra0, ra1);
        let f_ref_lo = ramp(q, rb0, rb1);
        let miss = (f_pred_lo - f_ref_hi).max(f_ref_lo - f_pred_hi).max(0.0);
        let wgt = if k == 0 || k == n - 1 { 0.5 } else { 1.0 };
        acc += miss * wgt;
    }
    acc * dq
}

/// Predicted box: [min, max] over the member set, widened by the declared
/// numeric half-band.
fn pred_box(members: &[f64]) -> (f64, f64) {
    let lo = members.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = members.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (lo * (1.0 - NUMERIC_REL), hi * (1.0 + NUMERIC_REL))
}

struct Score {
    qoi: &'static str,
    unit: &'static str,
    scale: f64,
    pred: (f64, f64),
    refi: (f64, f64),
    d: f64,
    d_rel: f64,
}

fn score(
    qoi: &'static str,
    unit: &'static str,
    scale: f64,
    members: &[f64],
    refi: (f64, f64),
) -> Score {
    let pred = pred_box(members);
    let d = ferson_d(pred, refi, SCATTER_REL);
    let mid = 0.5 * (refi.0 + refi.1);
    Score {
        qoi,
        unit,
        scale,
        pred,
        refi,
        d,
        d_rel: d / mid,
    }
}

fn score_table(w: &mut String, scores: &[Score]) {
    writeln!(
        w,
        "| QoI | predicted box | reference (pre-convolution) | Ferson d | d / ref-mid |"
    )
    .unwrap();
    writeln!(w, "|---|---|---|---|---|").unwrap();
    for s in scores {
        writeln!(
            w,
            "| {} | [{:.4}, {:.4}] {u} | [{:.4}, {:.4}] {u} | {} | {} |",
            s.qoi,
            s.pred.0 * s.scale,
            s.pred.1 * s.scale,
            s.refi.0 * s.scale,
            s.refi.1 * s.scale,
            if s.d == 0.0 {
                "**0 (overlap)**".to_string()
            } else {
                format!("{:.4} {u}", s.d * s.scale, u = s.unit)
            },
            if s.d == 0.0 {
                "—".to_string()
            } else {
                format!("{:.2}%", s.d_rel * 100.0)
            },
            u = s.unit,
        )
        .unwrap();
    }
    writeln!(w).unwrap();
}

fn render(w: &mut String) {
    let ref_cf_lo = REF_F.0 / (REF_PC.1 * A_T);
    let ref_cf_hi = REF_F.1 / (REF_PC.0 * A_T);

    let blind = [
        R_BLIND_NOM,
        R_B_E97_WLO,
        R_B_E97_WHI,
        R_B_E99_WLO,
        R_B_E99_WHI,
    ];
    let cal = [R_CAL_NOM, R_C_WLO, R_C_WHI];
    let get = |set: &[Readout], f: fn(&Readout) -> f64| -> Vec<f64> { set.iter().map(f).collect() };

    let blind_scores = [
        score("thrust F", "kN", 1e-3, &get(&blind, |r| r.f), REF_F),
        score("Isp", "s", 1.0, &get(&blind, |r| r.isp), REF_ISP),
        score("c*", "m/s", 1.0, &get(&blind, |r| r.c_star), REF_CSTAR),
        score(
            "C_F",
            "",
            1.0,
            &get(&blind, |r| r.c_f),
            (ref_cf_lo, ref_cf_hi),
        ),
        score("p_c", "MPa", 1e-6, &get(&blind, |r| r.p_c), REF_PC),
    ];
    let cal_scores = [
        score("thrust F", "kN", 1e-3, &get(&cal, |r| r.f), REF_F),
        score("Isp", "s", 1.0, &get(&cal, |r| r.isp), REF_ISP),
        score("c*", "m/s", 1.0, &get(&cal, |r| r.c_star), REF_CSTAR),
        score(
            "C_F",
            "",
            1.0,
            &get(&cal, |r| r.c_f),
            (ref_cf_lo, ref_cf_hi),
        ),
        score("p_c", "MPa", 1e-6, &get(&cal, |r| r.p_c), REF_PC),
    ];

    writeln!(w, "# CRUCIBLE Station 5 Certificate — the Blind RL10\n").unwrap();
    writeln!(
        w,
        "Goal-B station 5 (VAL-2 §3.1/§3.2; SOLV-7; COUP-3 §3.5; COUP-7 §3.2): the \
         RL10A-3-3A exists in this repository ONLY as data — the geometry-of-record \
         contour CSV, the preset TOMLs, and the pinned `lox_lh2_v0.4.0` equilibrium \
         surface — assembled by the one config-driven engine path and marched to a \
         settled state at the coarse tier (`cells_across_throat = 5`, 41×119, 2667 gas \
         cells, 12 flow-throughs, ~155 s laptop wall clock per member on the S2 \
         SDC-IMEX spine — implicit class-D liner conduction + Robin-Robin wall \
         exchange inside the step, the COUP-2 conservation audit armed on every one \
         of the ~15,800 steps of every member, zero violations). Every reported \
         number is a plane integral of the conserved field (SOLV-7); chamber pressure \
         and thrust are **emergent, never imposed**. **Spine of record for the scores \
         below: S2/S3** — see KNOWN LIMITS for the measured S4 delta and why the \
         scores are not re-derived here.\n"
    )
    .unwrap();
    writeln!(
        w,
        "**Labels (VAL-2 §3.2, VISION_SCOPE §9 v1.4.1):** every score below is labeled \
         blind or calibrated, and every score carries `development-observed: yes` — the \
         RL10 campaign is declared OPEN development (anchor deltas were watched while \
         the instrument was built); blind remains a mechanical INPUT property: the blind \
         members consume the design spec, the geometry of record, universal closures, \
         and the coax-injector-family η_c\\* prior band measured on other hardware — \
         no quantity measured on the RL10 itself. The calibrated members additionally \
         bind the anchor's own fitted η_c\\* = 0.9892 and TM-107318 component data.\n"
    )
    .unwrap();

    writeln!(w, "## The reference p-box (VAL-2 §3.1, N17)\n").unwrap();
    writeln!(
        w,
        "Published interval per QoI convolved with the anchor's own documented \
         model scatter (chamber sub-model Isp 440.3 s vs cycle 445.6 s ⇒ relative \
         half-width {:.3}% applied per QoI): Isp [440.3, 445.6] s; F [16.41, 16.50] \
         klbf = [{:.2}, {:.2}] kN; p_c [475 (cycle station), 482 (Table 2.5.1 \
         injector-face)] psia = [{:.4}, {:.4}] MPa; c\\* 7824 ft/s = {:.1} m/s; C_F \
         derived F/(p_c·A_t) ∈ [{:.4}, {:.4}] with A_t = π·(2.47 in)² = {:.6} m² \
         (the radius erratum, VAL-2 0.2.3). Pass = overlap (d = 0); d is reported \
         either way — a miss is a quantified model-form number, never a silent fail.\n",
        SCATTER_REL * 100.0,
        REF_F.0 * 1e-3,
        REF_F.1 * 1e-3,
        REF_PC.0 * 1e-6,
        REF_PC.1 * 1e-6,
        REF_CSTAR.0,
        ref_cf_lo,
        ref_cf_hi,
        A_T,
    )
    .unwrap();

    writeln!(
        w,
        "## BLIND score — open mode (label: blind, development-observed: yes)\n"
    )
    .unwrap();
    writeln!(
        w,
        "Members: the band-mid nominal (η_c\\* target 0.98, realized 0.9786 = \
         2208.8/2257.0 against the full-equilibrium baseline) and the four declared-band \
         corners — coax-family η_c\\* edge (realized 0.966–0.969 and 0.988–0.991; \
         S18 source-level knockdown, coarse-calibrated slope −0.847% c\\* per \
         −3×10⁵ J/kg) × wall-law band edge (the SOLV-1 §3.5 ±20–30% Colburn band, \
         realized as h × 0.75 / × 1.25). Boxes widened by the declared ±{:.1}% \
         numeric band (limit-cycle + flow-closure + integration class).\n",
        NUMERIC_REL * 100.0,
    )
    .unwrap();
    score_table(w, &blind_scores);
    writeln!(
        w,
        "F and Isp OVERLAP the reference at the band's upper edge; p_c, c\\*, and C_F \
         miss coherently (p_c and c\\* ~6% low, C_F correspondingly high — one \
         discretization signature, not three physics errors: at the coarse tier the \
         under-resolved throat/wall region under-produces chamber pressure while the \
         exit momentum integral is nearly converged, and c\\* = p_c·A_t/ṁ inherits \
         p_c's deficit while C_F = F/(p_c·A_t) inherits its inverse). The indicative \
         dial-8 refinement moved exactly these: p_c +4.7%, c\\* +4.7%, C_F −4.6%, \
         F +0.06% (§limits for why that tier is not certified tonight).\n"
    )
    .unwrap();

    writeln!(
        w,
        "## CALIBRATED score — closed expander (label: calibrated, development-observed: yes)\n"
    )
    .unwrap();
    writeln!(
        w,
        "The COUP-3 §3.5 fixed point closes the cycle each step (jacket pickup → \
         turbine power → pump map → delivered ṁ; Aitken-relaxed, residual ≤ 1e-8, \
         engaged after the establishment window): the engine finds its OWN operating \
         point. Nominal: delivered ṁ = 17.71 kg/s (+4.5% over design — the wall law's \
         +16% jacket pickup driven through the declared ṁ³ impedance line), turbine \
         671.7 kW, T_turbine_in 235.2 K. Members: nominal + the wall-band corners \
         (h × 0.75 ⇒ ṁ 16.74; h × 1.25 ⇒ ṁ 18.55 kg/s — the wall-function band IS \
         the dominant p_c spread, exactly as COUP-7 §3.4 predicted).\n"
    )
    .unwrap();
    score_table(w, &cal_scores);
    writeln!(
        w,
        "**F, Isp, and the emergent p_c all OVERLAP the record** (p_c box \
         [3.01, 3.36] MPa spans the published 3.27–3.32; the nominal alone reads \
         463.8 psia vs the 475–482 record). Isp is nearly flat across the whole \
         wall band (439.1–441.0 s vs record 440.3–445.6): the closed cycle trades \
         ṁ against p_c at almost constant specific impulse — real expander-cycle \
         self-regulation, reproduced by the coupled instrument, not imposed. c\\*/C_F \
         carry the same coarse-tier discretization signature as the blind score.\n"
    )
    .unwrap();

    writeln!(w, "## Consistency checks\n").unwrap();
    let id_worst = [&R_BLIND_NOM, &R_CAL_NOM]
        .iter()
        .map(|r| {
            let ve = r.c_star * r.c_f;
            ((r.f / (r.mdot_inj * 9.80665) * 9.80665 - ve) / ve).abs()
        })
        .fold(0.0f64, f64::max);
    writeln!(
        w,
        "- `v_e = c\\*·C_F` closes to {:.2e} relative on the nominals (SOLV-7 §6-2; \
         the residual is the ṁ_exit vs ṁ_inj closure, ≤ 0.2%).\n\
         - `p_c·A_t ≈ ṁ·c\\*` closes by construction of the c\\* readout (SOLV-7 N11: \
         the identity is the harness check, never the definition).\n\
         - Delivered-flow closure: the measured inflow-plane ṁ is within 0.2% of the \
         declared/solved value on every member (sonic startup cap inactive at \
         readout — the session-12 honesty signal).\n\
         - Steadiness: every member's residual is a stationary limit cycle \
         (max |Δρ|/ρ 0.6–4.3×10⁻² over the probe window — a plume-fringe cell-wise \
         max; the integral readouts' F oscillation is ≤ ±0.15% on every member, \
         measured ±0.08% on the noisiest), inside the declared numeric band.\n",
        id_worst
    )
    .unwrap();

    writeln!(w, "## Declared bands, model form, and declared devices\n").unwrap();
    writeln!(
        w,
        "- **Wall function (SOLV-1 §3.5):** the one Colburn-class law, ±20–30% declared \
         band — realized as the bracket corners; in closed mode it dominates the p_c \
         spread (as designed). Jacket pickup at the calibrated nominal: 9.76 MW vs the \
         record 8.43 MW (+16%, inside the band). Bartz nozzle-envelope oracle scoring: \
         deferred (recorded), rides the next wave.\n\
         - **η_c\\* prior (COUP-7 §3.2.1):** blind = the coax-family band applied as the \
         S18 source-level knockdown (never output-side); realized η re-measured per \
         member and reported above. Calibrated = the anchor's fitted 0.9892.\n\
         - **Equilibrium surface v0.4.0:** gas-only METASTABLE products (declared plume \
         model — real plumes supersaturate; deck-stamped), Z narrowed to the premixed \
         operating class, h ∈ [−1.23×10⁷, +1.2×10⁷] J/kg with rule-space (relative) \
         interp bounds incl. envelope-edge holdout (session-12 review fixes); the \
         projection acceptance uses the density column's own log-space bound. v0.4.0 \
         is the plan-S6 ignition-headroom STRICT extension of v0.3.2 (OFFL-3 0.6.2): \
         the recorded marches ran on v0.3.2, whose every in-envelope value and bound \
         is bit-identical in v0.4.0 (verified column-by-column at regeneration).\n\
         - **Geometry:** FND-3 analytic partial fractions + apertures (zero sampling \
         error) with State Redistribution (κ < 0.5, Berger–Giuliani) — the session-11 \
         stair-transpiration/starvation class is retired; wall heat runs on the \
         closure-vector (smooth) interface area, not the stair overcount.\n\
         - **Integrator (S2):** the ONE deterministic SDC-IMEX step (COUP-3 §3.1) — \
         explicit hyperbolic class + implicit class-D liner conduction (fixed-cycle \
         CG) with the Robin-Robin wall exchange inside each sweep (COUP-2 §3.5) and \
         the COUP-2 conservation audit armed every step; Δt is the gas CFL alone. \
         The liner ρc_p remains the declared steady-state continuation device until \
         the plan's S4 gives it the physical value (now legal under the implicit \
         class).\n\
         - **UQ:** these boxes are declared-band corner brackets (epistemic intervals), \
         NOT the full COUP-5 ensemble p-box — that machinery is a later wave; the \
         boxes are never collapsed to points.\n"
    )
    .unwrap();

    writeln!(w, "## KNOWN LIMITS (recorded, with owners)\n").unwrap();
    writeln!(
        w,
        "- **Certified tier = dial 5 (coarse).** Dials ≥ 8 do not survive the \
         ESTABLISHMENT march under the honest rule-space projection acceptance: the \
         startup transient (drain/shear at the bell wall; injector piston; or \
         overexpanded-bell backflow — five schedules probed, each halting loudly at a \
         different envelope edge) manufactures mixture states genuinely outside the \
         equilibrium surface's representable set. The pre-review dial-8 'success' rode \
         on the vacuous absolute acceptance bound the session-12 review retired. The \
         designated cure (VISION_SCOPE v1.5, session 13 — accelerated convergence is \
         DELETED; establishment stays a physical march) is the \
         **PLAN_CHEMICAL_SANDBOX Phase 1–2 physics**: gas diffusion (plan S3) + the \
         implicit integrator (**LANDED — this S2 rerun**) + the cold/unburnt \
         chemistry branch (plan S5–S6), which together make these transient states \
         representable and conductively coupled instead of refused; grid-sequenced \
         restart from settled physical states (FND-6) is the legal warm-start. The \
         indicative dial-8 deltas quoted above are direction-of-refinement evidence \
         only.\n\
         - **Discretization band is therefore declared, not swept:** the scored boxes \
         carry closure bands + the numeric band; the coarse-tier truncation error is \
         visibly ~5–6% on p_c/c\\*/C_F (the coherent signature above) and is NOT \
         hidden inside the boxes.\n\
         - `interp_error_bound_log` rides outside digest v3 (recorded digest-v4 \
         deferral). The station-4-fixture rewire deferral is DISCHARGED (S2): the \
         fixture and the engine now march the same one integrator.\n\
         - **The recorded readouts below predate plan S4 (session 16).** They were \
         earned on the S2/S3 spine with the wall law's chamber-fitted constant \
         transport and the ρc_p continuation device. S4 replaced both — the FND-7 \
         §3.3 tabulated spine now supplies per-cell μ/k/c_p/c_v and the liner carries \
         physical areal thermal capacitance — and the preset TOMLs cited under \
         Reproduction are the S4 ones. The **measured** S4 shift on the η = 1 baseline \
         member (`configs/rl10_coarse.toml`, 15 819 steps, audit clean): thrust \
         +0.30%, Isp +0.25%, c\\* −0.07%, C_F +0.32%, p_c −0.02%, and **jacket heat \
         −19.8%** (liner T_max 410.1 K) — every scored quantity moves an order of \
         magnitude less than the declared bands the boxes carry, while the wall term \
         (which no scored quantity reports, but which drives the closed cycle) moves by \
         two orders of magnitude more. That last number is why the CALIBRATED members \
         cannot simply be relabeled: the closed cycle converts jacket pickup into \
         delivered ṁ, so a −19.8% wall term is a first-order change to those scores and \
         **only** to those scores. The \
         scores are therefore NOT re-derived here: re-scoring all eight members on the \
         S4 spine is ~5 h of laptop march and it rides the **COUP-5 ensemble wave \
         (plan S18/S19)**, which replaces these hand-run corner brackets with sampled \
         p-boxes anyway. Two consequences are recorded rather than papered over: the \
         wall-law band corners are now the declared `band_factor` (h × 0.75/1.25 \
         directly) instead of the `cp_j_per_kg_k` proxy that reached the same h; and \
         the physical liner clock (~37 ms) exceeds the 12-flow-through settle budget \
         (~11 ms), so an S4 re-scoring must also grow that budget — run length, which \
         plan ruling #4 accepted, not schedule tuning.\n"
    )
    .unwrap();

    writeln!(w, "## Reproduction\n").unwrap();
    writeln!(
        w,
        "`cargo run --release -p crucible -- run configs/rl10_coarse.toml` (η = 1 \
         baseline), `configs/rl10_calibrated.toml` (calibrated nominal). Bracket \
         members = the coarse preset + the overrides stated in the recorded-readout \
         comments of `station5_rl10_certificate.rs` (note: those comments name the \
         S2/S3-spine overrides; S4 reaches the wall-law band corners through \
         `[mechanisms.wall] band_factor` instead). Table regen: \
         `offline/scripts/make_station5_tables.py` (~2 min). Every run writes its \
         resolved-config manifest + fields CSV under `runs/<name>/`; halts write \
         `crash_fields.csv` (raw conserved state — the session-12 instrument that \
         diagnosed every establishment pathology above).\n"
    )
    .unwrap();

    // Uncertified-indication sensitivity: d under the indicative refinement.
    writeln!(
        w,
        "## Sensitivity (uncertified indication): the blind score under the dial-8 refinement deltas\n"
    )
    .unwrap();
    let shifted: Vec<(usize, f64)> = vec![
        (4, 0.047),
        (2, 0.047),
        (3, -0.046),
        (0, 0.0006),
        (1, 0.0014),
    ];
    writeln!(
        w,
        "Applying the indicative deltas ({}) to the blind box edges: ",
        INDICATIVE_REFINEMENT
            .iter()
            .map(|(q, d)| format!("{q} {:+.1}%", d * 100.0))
            .collect::<Vec<_>>()
            .join(", ")
    )
    .unwrap();
    for (idx, delta) in shifted {
        let s = &blind_scores[idx];
        let pred = (s.pred.0 * (1.0 + delta), s.pred.1 * (1.0 + delta));
        let d = ferson_d(pred, s.refi, SCATTER_REL);
        writeln!(
            w,
            "- {}: d = {}",
            s.qoi,
            if d == 0.0 {
                "0 (overlap)".to_string()
            } else {
                format!(
                    "{:.4} {} ({:.2}% of ref-mid)",
                    d * s.scale,
                    s.unit,
                    d / (0.5 * (s.refi.0 + s.refi.1)) * 100.0
                )
            }
        )
        .unwrap();
    }
    writeln!(
        w,
        "\nThe residual blind miss under refinement is the wall-law + η-prior model \
         form the bands already carry — the certified sweep that would promote this \
         from indication to score awaits the PLAN_CHEMICAL_SANDBOX Phase 1–2 wave \
         (physical-march establishment with the full diffusion + cold-branch physics).\n"
    )
    .unwrap();

    let _ = (
        R_COARSE_ETA1,
        ferson_d as fn((f64, f64), (f64, f64), f64) -> f64,
    );
}

fn main() {
    let mut md = String::new();
    let w = &mut md;
    render(w);
    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../certificates/station5_rl10_certificate.md"
    );
    std::fs::write(path, md).expect("write certificate");
    println!("wrote {path}");
}
