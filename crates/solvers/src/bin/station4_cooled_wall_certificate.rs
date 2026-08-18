//! Regenerates `certificates/station4_cooled_wall_certificate.md` (GOAL-B
//! station 4). Criteria are CI-enforced in
//! `crates/solvers/tests/solv1_station4_cooled_wall.rs`; check.sh diffs the
//! committed artifact against this regeneration.

use crucible_solvers::station4_cooled_wall::*;
use std::fmt::Write as _;

fn main() {
    let (duct, rec, steps, resid) = run_duct().expect("coupled march");
    let r = duct_report(&duct, &rec, steps, resid);

    let (mut cav, cav_op) = build_stepped_cavity();
    let dirs: std::collections::BTreeSet<String> =
        cav.faces.iter().map(|f| format!("{:?}", f.dir)).collect();
    let (e_gas0, e_solid0) = cavity_energies(&cav);
    let mut cav_rec = ExchangeRecord::default();
    let mut t = 0.0;
    for _ in 0..500 {
        t += stepped_cavity_step(&mut cav, &cav_op, t, &mut cav_rec).expect("cavity step");
    }
    let (e_gas1, e_solid1) = cavity_energies(&cav);
    let lost = e_gas0 - e_gas1;
    let gained = e_solid1 - e_solid0;
    let mismatch = ((lost - gained) / gained).abs();

    let mut md = String::new();
    let w = &mut md;
    writeln!(w, "# CRUCIBLE Station 4 Certificate — the Cooled Wall\n").unwrap();
    writeln!(
        w,
        "Goal-B station 4 (SOLV-1 §3.5; COUP-2 §3.5 subset; META-3 `wall-function-heat`, \
         `robin-robin-cht` lineage): supersonic hot gas (M = {MACH_IN}, {T_IN:.0} K static — \
         free-stream recovery ≈ 1374 K) flows through a straight cylindrical duct \
         (R = {R_GAS} m, L = {DUCT_LEN} m) whose {LINER} m steel-class liner \
         (κ = {KAPPA_S} W/m·K) is regeneratively cooled outside \
         (h_cool = {H_COOL:.0} W/m²·K film at {T_COOL:.0} K). The gas-side flux comes from \
         the **one local wall-function law** — Colburn-class, near-wall state only, \
         **±20–30% declared band** — evaluated once per face per step and applied with \
         opposite signs to both sides (interface conservation by construction). The liner \
         conducts on the `Solid` region of the same world grid through the Goal-A-certified \
         operator; regions (gas/solid/exterior) are config-time data through the widened \
         FND-2 §3.6 ingest seam. Criteria are the named constants in \
         `station4_cooled_wall.rs`, asserted by `solv1_station4_cooled_wall.rs`.\n"
    )
    .unwrap();

    writeln!(
        w,
        "## The coupled duct at steady state ({} steps)\n",
        r.steps
    )
    .unwrap();
    writeln!(w, "| quantity | value |").unwrap();
    writeln!(w, "|---|---|").unwrap();
    writeln!(
        w,
        "| steadiness residual (ρ and T_solid, check window) | {:.2e} |",
        r.resid
    )
    .unwrap();
    writeln!(w, "| ṁ | {:.4} kg/s |", r.mdot).unwrap();
    writeln!(w, "| T0 in → out | {:.2} → {:.2} K |", r.t0_in, r.t0_out).unwrap();
    writeln!(
        w,
        "| gas enthalpy deficit ṁ·c_p·ΔT0 | {:.1} W |",
        r.gas_watts
    )
    .unwrap();
    writeln!(w, "| wall-face exchange Σq·A | {:.1} W |", r.wall_watts).unwrap();
    writeln!(w, "| coolant extraction | {:.1} W |", r.coolant_watts).unwrap();
    writeln!(
        w,
        "| ledger closure (wall vs gas, wall vs coolant) | {:.2e}, {:.2e} |",
        ((r.wall_watts - r.gas_watts) / r.wall_watts).abs(),
        ((r.wall_watts - r.coolant_watts) / r.wall_watts).abs()
    )
    .unwrap();
    writeln!(w, "| mid-duct film h | {:.1} W/m²·K |", r.h_mid).unwrap();
    writeln!(
        w,
        "| mid-duct wall flux q | {:.3e} W/m² (rocket-scale) |",
        r.q_mid
    )
    .unwrap();
    writeln!(
        w,
        "| liner ΔT at mid-duct (inner → outer) | {:.1} → {:.1} K |",
        r.t_liner_inner_mid, r.t_liner_outer_mid
    )
    .unwrap();
    writeln!(
        w,
        "| **series-resistance oracle, worst dev past entrance band** | **{:.2e}** |",
        r.oracle_worst
    )
    .unwrap();
    writeln!(
        w,
        "| T_aw(wall cell) / free-stream recovery at mid-duct | {:.3} |\n",
        r.t_aw_ratio_mid
    )
    .unwrap();
    writeln!(
        w,
        "**The oracle** (VAL-1 rung i): at steady state the coupled system must reproduce \
         the cylindrical film + ln-annulus + coolant-film series-resistance solution, built \
         from the simulated *gas* state and declared coolant data only — the simulated solid \
         field never enters. Pointwise agreement to {:.2e} (gate {ORACLE_REL_TOL:.0e}) past the {}-cell \
         entrance band says the wall law, the Robin faces, the region-masked conduction, and \
         the explicit flux-matched exchange compose into exactly the textbook conjugate \
         solution.\n",
        r.oracle_worst, ENTRANCE_BAND
    )
    .unwrap();

    writeln!(w, "## Stair-interface conservation (the stepped cavity)\n").unwrap();
    writeln!(
        w,
        "Hot gas at rest in a closed stepped cavity against a cold stair liner (face \
         directions exercised: {}), insulated exterior: after 500 coupled steps the gas \
         lost {:.6e} J, the liner gained {:.6e} J — relative mismatch **{:.2e}** \
         (round-off; the same q·A is applied to both sides and the shared face areas are \
         bitwise identical by the session-7 `face_radius` single-owner guarantee).\n",
        dirs.iter().cloned().collect::<Vec<_>>().join(", "),
        lost,
        gained,
        mismatch
    )
    .unwrap();

    writeln!(w, "## Exactness & determinism\n").unwrap();
    writeln!(
        w,
        "- **Uniform rest at coolant temperature is a bitwise fixed point** of the full \
         coupled step (gas fields and solid field byte-identical after 200 steps; the \
         wall machinery adds no drift — the exchange at equal temperatures is exactly \
         zero, unit-tested in `wall_heat`).\n\
         - **Reruns are bit-identical** on every field (asserted on a real march).\n\
         - **Robin-face analytic anchor:** the annulus with Dirichlet inner / Robin outer \
         reproduces the exact steady `T(r) = T1 + (T∞−T1)·ln(r/r1)/(ln(r2/r1)+κ/(h·r2))` \
         to < {ANNULUS_ROBIN_TOL:.0e} relative (32 radial cells).\n"
    )
    .unwrap();

    writeln!(w, "## Declared bands & honest scaffolding\n").unwrap();
    writeln!(
        w,
        "- The wall law carries its **±20–30% declared closure band** (SOLV-1 §3.5) — the \
         certificate's oracle validates the *coupling*, not the closure; the closure's truth \
         enters station 5 as a declared band in the p-box.\n\
         - **Near-wall sampling:** the wall-adjacent cell is itself cooled, so its recovery \
         temperature reads {:.1}% below free-stream at Δr = {DR} m. This cell-size \
         dependence of the un-resolved-boundary-layer closure is inside the declared band \
         and shrinks when FND-3 geometry + finer wall cells arrive.\n\
         - **Explicit flux-matched splitting** at the gas CFL dt is honest scaffolding \
         (guarded each step against both thermal stability limits, fail-loud): COUP-3's \
         class-`D` implicit solve with Robin-Robin Picard/Aitken sweeps supersedes it. \
         The liner ρc_p = {RHO_CP_S} J/m³·K is a **declared steady-state continuation \
         device** (steady solution independent of ρc_p; physical value ~3.6e6 only slows \
         settling).\n\
         - **Coolant side** is a configured Robin film — COUP-7's cooling-jacket boundary \
         object (channel correlation, coolant return state) supersedes it.\n\
         - **Bartz nozzle-envelope cross-check** (VAL-2, `bartz` oracle) rides with \
         station 5, where the nozzle + liner assembly exists; the stair-interface machinery \
         it needs is certified here by the stepped cavity.\n",
        (1.0 - r.t_aw_ratio_mid) * 100.0
    )
    .unwrap();

    let path = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../certificates/station4_cooled_wall_certificate.md"
    );
    std::fs::write(path, md).expect("write certificate");
    println!("wrote {path}");
}
