//! Renders the Goal-B Station-1 certificate to
//! `certificates/station1_sod_certificate.md` — the committed, regenerable
//! record of the bursting-diaphragm anchor. Criteria come from the SAME
//! named constants the test suite asserts
//! (`station1_sod::SHOCK_POS_TOL_CELLS` etc.), and gate 4 of
//! `scripts/check.sh` regenerates this file and diffs it against the
//! committed copy. Regenerate: `cargo run --release --bin
//! station1_sod_certificate` (deterministic — no RNG, no wall-clock).

use crucible_solvers::euler::{EULER_FIELDS, NCOMP};
use crucible_solvers::euler_mms::{
    MMS_LEVELS, MMS_ORDER_MAX, MMS_ORDER_MIN, mms_axisym_swirl, mms_theta_mode,
};
use crucible_solvers::station1_sod::{
    ADV_C_ORDER_MAX, ADV_C_ORDER_MIN, ADV_LEVELS, ADV_MEAN_ORDER_MIN, ADV_RAMP_WIDTH,
    ADV_RHO_STEP_ORDER_MIN, ADV_T_FINAL, C_BOUNDS_TOL, CONSERVATION_T_FINAL, CONSERVATION_TOL_REL,
    CONTACT_POS_TOL_CELLS, FAN_ORDER_MIN, FAN_WINDOW, SHOCK_POS_TOL_CELLS, SOD_FINEST_L1_MAX,
    SOD_L1_DECREASE_MIN, SOD_T_FINAL, STAR_L_WINDOW, STAR_PLATEAU_TOL, STAR_WINDOW_ORDER_MIN,
    UNIFORM_STEPS, advection_order, closed_tube_conservation, sod_convergence, sod_exact,
    sod_waves, uniform_state_is_bitwise_fixed_point,
};
use std::fmt::Write as _;

fn main() {
    let mut md = String::new();
    let w = &mut md;

    writeln!(
        w,
        "# CRUCIBLE Station 1 Certificate — the Bursting Diaphragm (Sod)"
    )
    .unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "Goal-B station 1 (SOLV-1 §6-1; VAL-2 §3.3 `sod-shock`): compressible Euler \
         on the unified cylindrical grid at N_θ = 1 — PPM reconstruction, HLLC flux \
         with Batten wavespeeds, well-balanced geometric sources, one passive \
         advected composition riding the contact. The tube is a full cylinder \
         (r = 0 axis inside the domain), radially uniform, marched by the ONE \
         production integrator — COUP-3's SDC-IMEX step (S2; explicit hyperbolic \
         class, fixed sweeps, the COUP-2 conservation audit armed every step; the \
         session-7 SSP-RK2 scaffolding is retired). Oracle: the \
         exact Riemann solution (Toro exact solver; star state verified against \
         Toro Table 4.2 in the unit tests). Criteria are the named constants in \
         `crates/solvers/src/station1_sod.rs`, asserted by \
         `crates/solvers/tests/solv1_station1_sod.rs`, and enforced against this \
         committed file by gate 4 of `scripts/check.sh`."
    )
    .unwrap();

    // --- Sod ladder ----------------------------------------------------------
    let study = sod_convergence();
    let po_g = study.observed_orders(|l| l.l1_rho);
    let po_f = study.observed_orders(|l| l.l1_fan);
    let po_s = study.observed_orders(|l| l.l1_star_l);
    writeln!(w).unwrap();
    writeln!(
        w,
        "## Sod tube vs exact Riemann solution (t = {SOD_T_FINAL})"
    )
    .unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "| n_z | steps | L1(ρ) global | order | L1(ρ) fan {FAN_WINDOW:?} | order | \
         L1(ρ) star-L {STAR_L_WINDOW:?} | order |"
    )
    .unwrap();
    writeln!(w, "|---|---|---|---|---|---|---|---|").unwrap();
    for (i, l) in study.levels.iter().enumerate() {
        let ord = |v: &[f64]| {
            if i == 0 {
                "—".to_string()
            } else {
                format!("{:.3}", v[i - 1])
            }
        };
        writeln!(
            w,
            "| {} | {} | {:.4e} | {} | {:.4e} | {} | {:.4e} | {} |",
            l.n_z,
            l.steps,
            l.l1_rho,
            ord(&po_g),
            l.l1_fan,
            ord(&po_f),
            l.l1_star_l,
            ord(&po_s),
        )
        .unwrap();
    }
    writeln!(w).unwrap();
    writeln!(
        w,
        "**Criteria:** global L1 falls ≥{SOD_L1_DECREASE_MIN}× per refinement and is \
         < {SOD_FINEST_L1_MAX:e} at the finest level; star-plateau window order ≥ \
         {STAR_WINDOW_ORDER_MIN} (formal order in smooth regions — measured 2.2–2.6); \
         fan-interior order ≥ {FAN_ORDER_MIN} (centered-rarefaction startup \
         singularity gives the known ~1st-order fan interior; gated against further \
         degradation, reported honestly)."
    )
    .unwrap();

    // --- Wave structure ------------------------------------------------------
    let wav = sod_waves();
    let ex = sod_exact();
    writeln!(w).unwrap();
    writeln!(w, "## Wave structure at n_z = {}", wav.n_z).unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "Exact star state: p* = {:.5}, u* = {:.5}, ρ*L = {:.5}, ρ*R = {:.5} \
         (Toro Table 4.2).",
        ex.p_star, ex.u_star, ex.rho_star_l, ex.rho_star_r
    )
    .unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "- **Shock position:** off by {:.3} cells — criterion ≤ {SHOCK_POS_TOL_CELLS} \
         (VAL-2: shock within 1 cell).",
        wav.shock_pos_err_cells
    )
    .unwrap();
    writeln!(
        w,
        "- **Contact position** (composition midpoint): off by {:.3} cells — \
         criterion ≤ {CONTACT_POS_TOL_CELLS} (the Batten-restored wave).",
        wav.contact_pos_err_cells
    )
    .unwrap();
    writeln!(
        w,
        "- **Star plateaus:** max errors ρ*L {:.2e}, ρ*R {:.2e}, u* {:.2e}, p* {:.2e} \
         — criterion ≤ {STAR_PLATEAU_TOL:e} each.",
        wav.rho_star_l_err, wav.rho_star_r_err, wav.u_star_err, wav.p_star_err
    )
    .unwrap();
    writeln!(
        w,
        "- **Composition bounds:** C ∈ [{:.2e}, 1 + {:.2e}] over the whole tube — \
         criterion within [−{C_BOUNDS_TOL:e}, 1 + {C_BOUNDS_TOL:e}] (species ride \
         the contact, no new extrema).",
        wav.c_min,
        wav.c_max - 1.0
    )
    .unwrap();
    writeln!(
        w,
        "- **Radial uniformity:** a radially-uniform tube stays radially uniform \
         **bitwise** through the whole shock evolution: {} (one law on the one \
         grid; the axis metric and geometric sources cancel exactly).",
        wav.radially_uniform_bitwise
    )
    .unwrap();

    // --- Smooth advection ----------------------------------------------------
    let adv = advection_order();
    let ao_r = adv.observed_orders(|l| l.l1_rho);
    let ao_c = adv.observed_orders(|l| l.l1_c);
    writeln!(w).unwrap();
    writeln!(
        w,
        "## Smooth formal-order study (advected tanh ramp, width {ADV_RAMP_WIDTH}, \
         t = {ADV_T_FINAL})"
    )
    .unwrap();
    writeln!(w).unwrap();
    writeln!(w, "| n_z | L1(ρ) | order | L1(C) | order |").unwrap();
    writeln!(w, "|---|---|---|---|---|").unwrap();
    for (i, l) in adv.levels.iter().enumerate() {
        let ord = |v: &[f64]| {
            if i == 0 {
                "—".to_string()
            } else {
                format!("{:.3}", v[i - 1])
            }
        };
        writeln!(
            w,
            "| {} | {:.4e} | {} | {:.4e} | {} |",
            l.n_z,
            l.l1_rho,
            ord(&ao_r),
            l.l1_c,
            ord(&ao_c),
        )
        .unwrap();
    }
    writeln!(w).unwrap();
    writeln!(
        w,
        "**Criteria:** mean ρ order across the ladder ≥ {ADV_MEAN_ORDER_MIN} \
         (measured {:.2}; per-step ≥ {ADV_RHO_STEP_ORDER_MIN} — ρ couples to the \
         acoustic families and oscillates per step); composition order per step in \
         [{ADV_C_ORDER_MIN}, {ADV_C_ORDER_MAX}] (pure contact family — textbook \
         2nd order). Levels {ADV_LEVELS:?}.",
        adv.mean_order(|l| l.l1_rho)
    )
    .unwrap();

    // --- Whole-operator MMS (SOLV-1 §6-2) ------------------------------------
    writeln!(w).unwrap();
    writeln!(w, "## Whole-operator MMS (SOLV-1 §6-2)").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "A manufactured smooth field with radial flow, swirl, and axial flow all \
         active — every flux direction and all three geometric source terms \
         (pressure, centrifugal ρu_θ², swirl advection ρu_ru_θ) carry nonzero \
         operands, which the Sod run structurally cannot exercise (its u_r ≡ 0). \
         The analytic residual enters through the operator's source intake; the \
         computed field must recover the manufactured one at formal order in \
         **L1 per conserved component** (the shock-capturing verification norm: \
         the limiter's clipping at the θ-mode's smooth extrema is locally \
         1st-order over an O(h) measure, which L2 amplifies to a ~O(h^1.6) tail \
         while L1 retains the formal order — measured and recorded here as the \
         honest caveat). Levels {MMS_LEVELS:?}."
    )
    .unwrap();
    for study in [mms_axisym_swirl(), mms_theta_mode()] {
        writeln!(w).unwrap();
        writeln!(w, "### {}", study.label).unwrap();
        writeln!(w).unwrap();
        write!(w, "| n |").unwrap();
        for name in EULER_FIELDS {
            write!(w, " L1({name}) | order |").unwrap();
        }
        writeln!(w).unwrap();
        write!(w, "|---|").unwrap();
        for _ in 0..NCOMP {
            write!(w, "---|---|").unwrap();
        }
        writeln!(w).unwrap();
        for (i, lvl) in study.levels.iter().enumerate() {
            write!(w, "| {} |", lvl.n).unwrap();
            for k in 0..NCOMP {
                let p = if i == 0 {
                    "—".to_string()
                } else {
                    format!("{:.2}", study.observed_orders(k)[i - 1])
                };
                write!(w, " {:.2e} | {p} |", lvl.l1[k]).unwrap();
            }
            writeln!(w).unwrap();
        }
        writeln!(w).unwrap();
        writeln!(
            w,
            "**Criterion:** every component's observed order in \
             [{MMS_ORDER_MIN}, {MMS_ORDER_MAX}] (formal = 2) at every refinement."
        )
        .unwrap();
    }

    // --- Well-balance, conservation, determinism -----------------------------
    let (dm, de) = closed_tube_conservation().expect("conservation run");
    writeln!(w).unwrap();
    writeln!(w, "## Well-balance, conservation, determinism").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "- **Uniform gas at rest is a bitwise fixed point** over {UNIFORM_STEPS} \
         steps — axis config (N_θ = 1, r = 0 inside): {}; annulus config \
         (N_θ = 8): {}. The geometric sources use the same A·p products as the \
         face fluxes, so the balance is exact, not approximate.",
        uniform_state_is_bitwise_fixed_point(true),
        uniform_state_is_bitwise_fixed_point(false)
    )
    .unwrap();
    writeln!(
        w,
        "- **Closed reflecting tube through wall reflections** (t = \
         {CONSERVATION_T_FINAL}): mass drift {dm:.3e}, energy drift {de:.3e} \
         relative — criterion < {CONSERVATION_TOL_REL:e} each (flux-form \
         telescoping; measured at round-off). Momentum is exchanged with the \
         walls — that ledger is COUP-2 §3.1.2's mount-reaction term, a later \
         wave."
    )
    .unwrap();
    writeln!(
        w,
        "- **Determinism:** reruns are byte-identical (asserted in the test suite \
         on the full conserved state)."
    )
    .unwrap();

    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let dir = format!("{root}/certificates");
    std::fs::create_dir_all(&dir).expect("create certificates/");
    let path = format!("{dir}/station1_sod_certificate.md");
    std::fs::write(&path, &md).expect("write certificate");
    println!("{md}");
    println!("written: {path}");
}
