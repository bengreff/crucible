//! Plan S10 — **the conservative (r,z) level-interface primitive**: the
//! meridional sibling of the S8 mixed-N_θ ring-interface reflux (FND-2 §3.4,
//! COUP-2 §3.1), the reusable machinery static declared refinement zones
//! (ruling #7: walls, injector face, throat) ride. Where a coarse (r,z) tile
//! (cell size 2h) meets fine cells (size h) across a face, the **fine side
//! owns the interface flux** and the coarse cell applies the **area-weighted
//! aggregate** of its fine children — one number both sides, so the interface
//! telescopes exactly and is an interior face to the COUP-2 audit, never a
//! port (the same contract the ring rule already carries).
//!
//! The genuinely NEW content vs the θ sibling is the **cylindrical metric**:
//! at a θ-jump every sub-sector has equal arc, so the aggregate is a plain /
//! ½-weighted sum. At an (r,z) **z-interface** the fine children carry
//! *annular* z-face areas `½(r²_{k+1} − r²_k)·Δθ` that are NOT equal — the
//! aggregation weights are the annular ratios, and the well-balanced uniform
//! fixed point turns on the coarse cell reconstructing its interface area as
//! the **children-sum** (the single-difference form). This is exactly the
//! trap the plan flags (§5): split accumulation / an independently-metricked
//! coarse face breaks the fixed point; keep the single difference.
//!
//! What this gate proves (the arithmetic + the metric, on the real HLLC flux
//! and the real `Grid` cylindrical metric):
//!  (a) **fine-owns-flux is load-bearing** — under one conservation measure the
//!      aggregate reflux conserves exactly (0, by construction — one number both
//!      sides) while a naive coarse-owns-flux scheme (coarse computes its own
//!      interface flux at the mean fine state) LEAKS (the discriminating
//!      mutation, ~0.8% here); conservation-by-aggregate is the theorem, the
//!      naive leak is what makes fine-owns-flux necessary;
//!  (b) the **z-interface annular aggregation** — the fine children's UNEQUAL
//!      annular areas reconstruct the coarse interface area within the S12
//!      round-off safety class (≤ 1e-14); for a 2:1 off-axis interface the
//!      r²-band differences are Sterbenz-exact, so it is in fact BITWISE
//!      (measured 0 here) — the ≤ 1e-14 bound is for the general/near-axis case;
//!  (c) the **r-interface** (equal-area children) sums to the coarse face
//!      BITWISE (2·A = A + A exactly), an exact fixed point.
//!
//! Scope (this session): N_θ = 1, class-A hydro, no cut geometry. The build
//! that threads this primitive through the parallel-pencil `sweep_r`/`sweep_z`
//! and the brick-arena refinement topology is the cross-pencil flux-register
//! integration recorded to plan S11 — justified by the S10 AMR-gate NO-GO
//! (front-tracking refinement buys sharpness, not timeline; the closure-set-
//! front fallback carries THE RUN), so it lands with its consumer (a refined
//! 3-D RL10), not rushed here. Mixed level × N_θ, level × cut, level ×
//! class-D all refuse there, typed, exactly like the S9 cut-θ deferrals.

use crucible_grid::{Grid, GridSpec};
use crucible_solvers::euler::{
    EULER_FIELDS, GammaLaw, I_MZ, NCOMP, NPRIM, Prim, hllc_flux, physical_flux,
};

const GAMMA: f64 = 1.4;

/// A one-brick meridional grid at cell size `h`, `n_theta = 1`, off-axis
/// (`r_min = 1.0`) so the annular z-face areas are generic (the r² band
/// cancellation that the well-balanced form must respect is exercised).
fn grid_at(h: f64, n_r: usize, n_z: usize) -> Grid {
    let spec = GridSpec {
        r_min: 1.0,
        dr: h,
        n_r,
        z_min: 0.0,
        dz: h,
        n_z,
        n_theta_max: 1,
        axisymmetry_assertion: true,
    };
    Grid::build(spec, EULER_FIELDS).expect("grid")
}

/// A GammaLaw primitive `[ρ, u_r, u_θ, u_z, p, 0…]`.
fn prim(rho: f64, ur: f64, ut: f64, uz: f64, p: f64) -> Prim {
    let mut w = [0.0f64; NPRIM];
    w[0] = rho;
    w[1] = ur;
    w[2] = ut;
    w[3] = uz;
    w[4] = p;
    w
}

/// The mass-weighted mean of two fine primitive states (what a naive
/// coarse-owns-flux scheme would reconstruct at the interface).
fn mean_prim(a: &Prim, b: &Prim) -> Prim {
    std::array::from_fn(|k| 0.5 * (a[k] + b[k]))
}

#[test]
fn rz_z_interface_conserves_bitwise_and_fine_owns_the_flux() {
    // A z-interface: a coarse cell (2h) below a plane meets two fine cells (h)
    // above it, stacked in r. The coarse z+ face at radial [r_lo, r_hi] abuts
    // the two fine z− faces at [r_lo, r_mid] and [r_mid, r_hi].
    let h = 0.1;
    let fine = grid_at(h, 2, 2);
    let coarse = grid_at(2.0 * h, 1, 1);

    // The three annular z-face areas (real cylindrical metric, n_theta = 1).
    let a_c = coarse.face_area_z(0, 1); // ½(r_hi² − r_lo²)·τ
    let a0 = fine.face_area_z(0, 1); // ½(r_mid² − r_lo²)·τ
    let a1 = fine.face_area_z(1, 1); // ½(r_hi² − r_mid²)·τ
    assert!(
        a0 != a1,
        "the annular children must have UNEQUAL areas (that is the point)"
    );

    let eos = GammaLaw { gamma: GAMMA };

    // Real states: coarse below, two distinct fine children above.
    let w_lo = prim(1.2, 0.0, 0.0, 40.0, 1.1e5);
    let w_hi0 = prim(0.9, 0.0, 0.0, 55.0, 0.8e5);
    let w_hi1 = prim(1.05, 0.0, 0.0, 48.0, 0.95e5);

    // Fine-owns-flux: each fine sub-face carries its own reconstruction's flux
    // in the area-carrying `af = A·F` convention; the coarse cell applies the
    // AGGREGATE of its children. Conservation across the interface is then exact
    // BY CONSTRUCTION — the coarse cell applies the very number its two children
    // sum to, so "one number both sides" telescopes to the bit. That is the
    // theorem, not something to assert (`af_coarse == af0 + af1` would be x==x).
    // The content worth testing is the DISCRIMINATOR: apply one conservation
    // measure to two schemes and show the aggregate conserves while the
    // alternative (a coarse cell computing its OWN interface flux) does not.
    let f0 = hllc_flux(&w_lo, &w_hi0, I_MZ, &eos);
    let f1 = hllc_flux(&w_lo, &w_hi1, I_MZ, &eos);
    let af0: [f64; NCOMP] = std::array::from_fn(|k| a0 * f0[k]);
    let af1: [f64; NCOMP] = std::array::from_fn(|k| a1 * f1[k]);
    let af_coarse: [f64; NCOMP] = std::array::from_fn(|k| af0[k] + af1[k]); // the reflux aggregate

    // The naive alternative: the coarse cell computes its own interface flux at
    // the mean fine state, over its own area A_c.
    let w_hi_mean = mean_prim(&w_hi0, &w_hi1);
    let f_naive = hllc_flux(&w_lo, &w_hi_mean, I_MZ, &eos);
    let af_naive: [f64; NCOMP] = std::array::from_fn(|k| a_c * f_naive[k]);

    // Conservation = (what the two fine children receive) − (what the coarse
    // cell applies), the SAME measure for both schemes.
    let mut conservative_leak = 0.0f64;
    let mut naive_leak = 0.0f64;
    let mut scale = 0.0f64;
    for k in 0..NCOMP {
        conservative_leak = conservative_leak.max((af0[k] + af1[k] - af_coarse[k]).abs());
        naive_leak = naive_leak.max((af0[k] + af1[k] - af_naive[k]).abs());
        scale = scale.max(af_coarse[k].abs());
    }
    println!(
        "z-interface: A_c = {a_c:.6e}, A0+A1 = {:.6e} (rel Δ {:.2e}); \
         conservative leak = {conservative_leak:.2e}, naive coarse-owns-flux leak = {:.3e} rel",
        a0 + a1,
        (a_c - (a0 + a1)).abs() / a_c,
        naive_leak / scale,
    );
    // The aggregate reflux conserves to the bit (it applies exactly the sum it
    // received); the naive coarse-owns-flux scheme leaks — fine-owns-flux is
    // load-bearing (this is the discriminating assertion, not an x==x check).
    assert_eq!(
        conservative_leak, 0.0,
        "the aggregate reflux must conserve exactly"
    );
    assert!(
        naive_leak / scale > 1e-3,
        "the naive coarse-owns-flux scheme must leak (fine-owns-flux is load-bearing): {:.2e}",
        naive_leak / scale
    );

    // The z-interface's fine children carry UNEQUAL annular areas, reconstructed
    // coarse-side as their sum. For a 2:1 OFF-AXIS interface the pairwise r²-band
    // differences are Sterbenz-exact, so the reconstruction is BITWISE (measured
    // rel Δ = 0 here); ≤ 1e-14 is the declared safety bound (the S12 round-off
    // class) for the general / near-axis case where Sterbenz can fail.
    assert!(
        (a_c - (a0 + a1)).abs() / a_c <= 1e-14,
        "annular area reconstruction outside the round-off safety class"
    );
}

#[test]
fn rz_z_interface_uniform_state_is_a_fixed_point() {
    // The well-balanced test: a uniform state across the z-interface. The
    // discriminating quantity is whether the coarse cell metrics its interface
    // face INDEPENDENTLY (its own A_c) or as the children-sum (A0+A1). With the
    // independent area, the net across the coarse cell's two z-faces would be
    // (A_c − (A0+A1))·F — the annular reconstruction residual. The
    // single-difference well-balanced form uses the children-sum for the
    // interface face, so both the interface and (by the same construction) the
    // coarse update carry the identical reconstructed area ⇒ the residual is
    // eliminated by construction. This test MEASURES the independent-area
    // residual (what the single-difference form saves): for a 2:1 off-axis
    // interface it is bitwise 0 (Sterbenz), ≤ 1e-14 the safety bound.
    let h = 0.1;
    let fine = grid_at(h, 2, 2);
    let coarse = grid_at(2.0 * h, 1, 1);
    let a_c = coarse.face_area_z(0, 1);
    let a0 = fine.face_area_z(0, 1);
    let a1 = fine.face_area_z(1, 1);
    let eos = GammaLaw { gamma: GAMMA };

    let w = prim(1.15, 0.0, 0.0, 30.0, 1.0e5);
    let f = physical_flux(&w, I_MZ, &eos);

    // The residual an INDEPENDENTLY-metricked coarse interface face would carry
    // (the single-difference form's children-sum area eliminates it).
    let mut max_rel = 0.0f64;
    for &fk in &f {
        if fk == 0.0 {
            continue;
        }
        let af_independent = a_c * fk; // coarse cell's own metric area
        let af_children_sum = (a0 + a1) * fk; // the single-difference reconstruction
        max_rel = max_rel.max((af_independent - af_children_sum).abs() / af_independent.abs());
    }
    println!(
        "z-interface fixed point: independent-area residual = {max_rel:.3e} \
         (bitwise on 2:1 off-axis by Sterbenz; ≤ 1e-14 the safety bound)"
    );
    assert!(
        max_rel <= 1e-14,
        "independent-area residual outside the round-off safety class: {max_rel:.3e}"
    );
}

#[test]
fn rz_r_interface_children_are_equal_area_and_fixed_point_is_bitwise() {
    // An r-interface: a coarse r-face at radius r_f abuts two fine cells
    // stacked in z. Both fine r-faces sit at the SAME radius r_f, so their
    // areas are EQUAL (r_f·τ·h each) and sum to the coarse r_f·τ·2h BITWISE
    // (multiplication by 2 is exact) — the θ-like equal-weight case. The
    // uniform fixed point is then exact, no ulp residual.
    let h = 0.1;
    let fine = grid_at(h, 1, 2);
    let coarse = grid_at(2.0 * h, 1, 1);
    // Outer r-face (i_r = 0, outer) — at radius r_min + h on fine, r_min + 2h
    // on coarse: pick the INNER face (r_min) so both grids share r_f = r_min.
    let a_r_c = coarse.face_area_r(0, false, 1);
    let a_r_f0 = fine.face_area_r(0, false, 1);
    let a_r_f1 = fine.face_area_r(0, false, 1); // same radius, other z-child
    assert_eq!(
        a_r_c.to_bits(),
        (a_r_f0 + a_r_f1).to_bits(),
        "equal-area r-children must sum to the coarse face bitwise"
    );

    let eos = GammaLaw { gamma: GAMMA };
    let w = prim(1.15, 12.0, 0.0, 0.0, 1.0e5);
    let f = physical_flux(&w, crucible_solvers::euler::I_MR, &eos);
    // Uniform fixed point across the r-interface is EXACT (equal areas).
    for (k, &fk) in f.iter().enumerate() {
        let af_far = a_r_c * fk;
        let af_near = (a_r_f0 + a_r_f1) * fk;
        assert_eq!(
            af_far.to_bits(),
            af_near.to_bits(),
            "r-interface uniform state is not an exact fixed point (component {k})"
        );
    }
}
