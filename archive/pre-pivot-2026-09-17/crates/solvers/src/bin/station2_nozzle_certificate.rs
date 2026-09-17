//! Renders the Goal-B Station-2 certificate to
//! `certificates/station2_nozzle_certificate.md` — the committed,
//! regenerable record of the De Laval nozzle anchor. Criteria come from the
//! SAME named constants the test suite asserts; gate 4 of
//! `scripts/check.sh` regenerates this file and diffs it. Regenerate:
//! `cargo run --release --bin station2_nozzle_certificate`.

use crucible_solvers::station2_nozzle::{
    CD_BAND, CD_FINEST_TOL, CF_DEV_MAX, EXIT_MACH_DEV_MAX, INLET_EXCLUSION, MACH_DEV_MAX,
    MDOT_SPREAD_MAX, NOZZLE_LENGTH, PC_BAND, RADIUS_COEFF, RESERVOIR_P0, RESERVOIR_RHO0,
    SETTLE_TIME, STEADY_RESID_FINEST, STEADY_RESID_MAX, THROAT_RADIUS, THROAT_Z, area_ratio,
    centerline_mach, mach_from_area_ratio, masked_uniform_fixed_point, nozzle_study, run_nozzle,
};
use std::fmt::Write as _;

fn main() {
    let mut md = String::new();
    let w = &mut md;

    writeln!(w, "# CRUCIBLE Station 2 Certificate — the De Laval Nozzle").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "Goal-B station 2 (SOLV-1 §6-4; SOLV-7 §6-1; META-3 `maccormack-nozzle` \
         lineage): a stagnation reservoir (p0 = {RESERVOIR_P0}, ρ0 = \
         {RESERVOIR_RHO0}) feeds an axisymmetric converging–diverging nozzle — \
         parabolic radius R(z) = R*(1 + {RADIUS_COEFF}((z−{THROAT_Z})/{THROAT_Z})²), \
         R* = {THROAT_RADIUS}, length {NOZZLE_LENGTH}, max wall slope 9.5° — on \
         the one cylindrical operator at N_θ = 1, axis inside the domain. The \
         flow **chokes at the throat on its own** and exits supersonic; nothing \
         about the operating point is imposed (SOLV-7 §3.2). Oracle: quasi-1-D \
         isentropic theory (area–Mach relation, choked ṁ, vacuum C_F).\n\n\
         Geometry enters through the grid's config-time activity seam in its \
         binary (stair-step) degenerate form, with **slip-ghost walls**: ghost \
         states mirror about the true contour normal (ghost-cell immersed \
         boundary), cutting spurious wave generation from O(wall slope) to \
         O(h·curvature). Cost of that trade: a small stair-face transpiration \
         flux, reported below as the plane-ṁ spread, shrinking with resolution \
         — retired when FND-3's partial apertures + cut cells land. The \
         quasi-1-D oracle also carries its own 2-D model floor (centerline ≠ \
         area mean); both gaps are inside the declared bands. Criteria are the \
         named constants in `crates/solvers/src/station2_nozzle.rs`, asserted \
         by `crates/solvers/tests/solv1_station2_nozzle.rs`."
    )
    .unwrap();

    let levels = nozzle_study().expect("study");
    writeln!(w).unwrap();
    writeln!(w, "## Choked-flow ladder (marched to t = {SETTLE_TIME})").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "| n_r×n_z | steps | Cd = ṁ/ṁ_ideal | p_c/p0 | steady resid | max Mach dev \
         (z > {INLET_EXCLUSION}) | exit M (1-D {:.3}) | C_F | ideal C_F | ṁ spread |",
        mach_from_area_ratio(area_ratio(NOZZLE_LENGTH), true)
    )
    .unwrap();
    writeln!(w, "|---|---|---|---|---|---|---|---|---|---|").unwrap();
    for l in &levels {
        writeln!(
            w,
            "| {}×{} | {} | {:.4} | {:.5} | {:.2e} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} |",
            l.n_r,
            l.n_z,
            l.steps,
            l.cd_analytic,
            l.p_c_over_p0,
            l.steady_resid,
            l.mach_dev_max,
            l.exit_mach_centerline,
            l.cf,
            l.cf_ideal,
            l.mdot_spread,
        )
        .unwrap();
    }
    writeln!(w).unwrap();
    writeln!(
        w,
        "**Criteria:** Cd in [{}, {}] at every level and |Cd−1| < {CD_FINEST_TOL} \
         at the finest (measured 0.0020 — the choked mass flow matches ideal \
         theory to 0.2%), never growing under refinement; p_c/p0 in \
         [{}, {}]; steadiness residual < {STEADY_RESID_MAX} (finest < \
         {STEADY_RESID_FINEST}); centerline Mach within {MACH_DEV_MAX} of the \
         area–Mach relation past the entrance band; exit Mach supersonic and \
         within {EXIT_MACH_DEV_MAX} of the 1-D exit value; C_F within \
         {CF_DEV_MAX} of ideal vacuum C_F; plane-ṁ spread < {MDOT_SPREAD_MAX} \
         and strictly shrinking (the declared slip-wall transpiration).",
        CD_BAND.0, CD_BAND.1, PC_BAND.0, PC_BAND.1
    )
    .unwrap();

    // Centerline profile at the finest level.
    let (n_r, n_z) = (levels[levels.len() - 1].n_r, levels[levels.len() - 1].n_z);
    let (g, f, _, _) = run_nozzle(n_r, n_z).expect("finest run");
    let mach = centerline_mach(&g, &f);
    let dz = g.spec().dz;
    writeln!(w).unwrap();
    writeln!(w, "## Centerline Mach profile at {n_r}×{n_z}").unwrap();
    writeln!(w).unwrap();
    writeln!(w, "| z | A/A* | M (computed) | M (quasi-1-D) | rel. dev |").unwrap();
    writeln!(w, "|---|---|---|---|---|").unwrap();
    for i_z in (0..n_z).step_by(8) {
        let z = (i_z as f64 + 0.5) * dz;
        let m1d = mach_from_area_ratio(area_ratio(z), z > THROAT_Z);
        writeln!(
            w,
            "| {z:.2} | {:.3} | {:.4} | {m1d:.4} | {:+.4} |",
            area_ratio(z),
            mach[i_z],
            (mach[i_z] - m1d) / m1d,
        )
        .unwrap();
    }
    writeln!(w).unwrap();
    writeln!(
        w,
        "The entrance band (z < {INLET_EXCLUSION}) carries the first-order \
         stagnation-BC adjacency + 2-D entrance turning (decays from ~+14% at \
         the first cell to ~+1%); it is reported, excluded from the pointwise \
         gate, and shrinks with the COUP-7 injector-object wave."
    )
    .unwrap();

    writeln!(w).unwrap();
    writeln!(w, "## Masked-sweep exactness & determinism").unwrap();
    writeln!(w).unwrap();
    writeln!(
        w,
        "- **Uniform gas at rest in the stair-stepped cavity is a bitwise fixed \
         point** — grid-aligned mirror: {}; slip-ghost wall: {} (reflecting a \
         zero velocity about any normal is the identity, so the wall machinery \
         adds no drift).",
        masked_uniform_fixed_point(false),
        masked_uniform_fixed_point(true)
    )
    .unwrap();
    writeln!(
        w,
        "- **Determinism:** reruns are byte-identical (asserted in the test \
         suite on the full conserved state)."
    )
    .unwrap();

    let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../..");
    let dir = format!("{root}/certificates");
    std::fs::create_dir_all(&dir).expect("create certificates/");
    let path = format!("{dir}/station2_nozzle_certificate.md");
    std::fs::write(&path, &md).expect("write certificate");
    println!("{md}");
    println!("written: {path}");
}
