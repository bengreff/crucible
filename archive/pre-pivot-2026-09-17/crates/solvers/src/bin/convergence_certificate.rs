//! Renders the Goal-A convergence certificate to
//! `certificates/convergence_certificate.md` — the committed, regenerable
//! artifact. Criteria come from the SAME named constants the test suite
//! asserts (`certificate::MMS_ORDER_*` etc.), and gate 4 of
//! `scripts/check.sh` regenerates this file and diffs it against the
//! committed copy, so artifact and code cannot silently diverge.
//! Regenerate with: `cargo run --bin convergence_certificate`.

use crucible_solvers::certificate::{
    ANNULUS_CELLS, ANNULUS_TOL_REL, BESSEL_TOL_K, CONSERVATION_STEPS, CONSERVATION_TOL_REL,
    MMS_ORDER_MAX, MMS_ORDER_MIN, annulus_anchor, bessel_cylinder_anchor, conservation_drift,
    mms_axisymmetric, mms_theta_mode,
};
use std::fmt::Write as _;

fn main() {
    let mut md = String::new();
    let w = &mut md;

    writeln!(w, "# CRUCIBLE Convergence Certificate (Goal A)").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "Heat conduction on the unified cylindrical world-state grid, driven \
         config → loader → registry → grid → operator → this data. Criteria are \
         the named constants in `crates/solvers/src/certificate.rs`, asserted by \
         `crates/solvers/tests/goal_a_certificate.rs` and enforced against this \
         committed file by gate 4 of `scripts/check.sh` (regenerate-and-diff). \
         Regenerate: `cargo run --bin convergence_certificate` (deterministic — \
         no RNG, no wall-clock; the git commit records provenance)."
    )
    .unwrap();

    for study in [mms_axisymmetric(), mms_theta_mode()] {
        writeln!(w).unwrap();
        writeln!(w, "## {}", study.label).unwrap();
        writeln!(w).unwrap();
        writeln!(w, "| n | h | L2 error | observed order |").unwrap();
        writeln!(w, "|---|---|---|---|").unwrap();
        let orders = study.observed_orders();
        for (i, lvl) in study.levels.iter().enumerate() {
            let p = if i == 0 {
                "—".to_string()
            } else {
                format!("{:.3}", orders[i - 1])
            };
            writeln!(
                w,
                "| {} | {:.5} | {:.6e} | {} |",
                lvl.n, lvl.h, lvl.l2_error, p
            )
            .unwrap();
        }
        writeln!(w).unwrap();
        writeln!(
            w,
            "**Criterion:** observed order in [{MMS_ORDER_MIN}, {MMS_ORDER_MAX}] (formal = 2)."
        )
        .unwrap();
    }

    let (annulus_rel, _) = annulus_anchor();
    writeln!(w).unwrap();
    writeln!(w, "## Analytic anchors").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "- **Steady annulus (log profile), {ANNULUS_CELLS} radial cells:** max error \
         {annulus_rel:.3e} relative to ΔT = 200 K — criterion < {ANNULUS_TOL_REL:e}."
    )
    .unwrap();
    let bessel_abs = bessel_cylinder_anchor();
    writeln!(
        w,
        "- **Transient cylinder vs 5-term Bessel series at t̃ = 0.1 (crosses the \
         r = 0 axis):** max error {bessel_abs:.3e} K on T₀ = 100 K — criterion < \
         {BESSEL_TOL_K} K."
    )
    .unwrap();
    let drift = conservation_drift(CONSERVATION_STEPS);
    writeln!(
        w,
        "- **Closed insulated sweep, {CONSERVATION_STEPS} steps:** total-energy drift \
         {drift:.3e} relative — criterion < {CONSERVATION_TOL_REL:e} (flux-form \
         telescoping)."
    )
    .unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "Determinism: reruns are byte-identical (asserted in the test suite: field \
         bits and config content hash)."
    )
    .unwrap();

    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let dir = format!("{root}/certificates");
    std::fs::create_dir_all(&dir).expect("create certificates/");
    let path = format!("{dir}/convergence_certificate.md");
    std::fs::write(&path, &md).expect("write certificate");
    println!("{md}");
    println!("written: {path}");
}
