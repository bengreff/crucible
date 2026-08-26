//! FND-3 §3.3 — PLIC interface planes: volume-exact plane placement in the
//! cell's local Cartesian box via the closed-form Scardovelli–Zaleski
//! plane/box volume relation (Scardovelli & Zaleski, J. Comput. Phys. 164
//! (2000) 228 — the standard clamped-cubic corner-sum form), solved here by
//! DETERMINISTIC bracketed bisection with a fixed iteration count, and the
//! interface area from the exact identity `area = dV/dd` (the derivative of
//! the same corner sum — never polygon clipping, never Marching-Cubes,
//! whose areas are ~8% biased per §3.3).

use super::Geom3dError;

/// Fixed bisection count for the volume-matching offset. 80 halvings of
/// the support bracket reach 2⁻⁸⁰ of the box diagonal — far below f64
/// resolution — and a FIXED count (no residual-dependent early exit) keeps
/// the solve a pure function of its inputs (§3.5 determinism).
pub const N_PLIC_BISECT: usize = 80;

/// Relative floor below which a normal component is treated as exactly
/// zero (its dimension "flat"): `|n_i|·L_i ≤ PLIC_DIM_FLOOR · Σ|n_j|·L_j`.
/// Rationale: the corner sum divides by Π n_i, so a near-zero component
/// amplifies round-off in the difference of near-equal cubes by
/// `~ Σ|n_j|L_j / (|n_i|L_i)`, while flattening the dimension perturbs the
/// matched volume by at most `~ ½|n_i|L_i/Σ|n_j|L_j` of the box — the two
/// errors cross near `√ε_f64 ≈ 1.5e-8`, so 1e-7 keeps BOTH below ~1e-7 of
/// the box (far under the ε_α sampling error the volume match consumes;
/// measured in the §6.6 battery: a 7.9e-13 un-snapped component moved d by
/// 2e-7 — this floor is what removes that). A fixed threshold is a fixed
/// tie-break (§3.5), not a physics clamp.
const PLIC_DIM_FLOOR: f64 = 1e-7;

/// One PLIC interface plane of a cut cell (FND-3 §3.3).
///
/// Frame convention (documented once): `normal` is the unit outward normal
/// OF THE SOLID INTERFACE pointing INTO the gas (the direction of
/// increasing classifier field), expressed in the cell's LOCAL Cartesian
/// frame components (r̂, θ̂, ẑ at the cell's metric midpoint). The plane is
/// `{ξ : normal·ξ = d}` with ξ the local offset from that midpoint (m);
/// the solid occupies `normal·ξ ≤ d`. `centroid` is the S9 simple form —
/// the cell midpoint projected onto the plane, in GLOBAL Cartesian
/// coordinates (the exact polygon centroid is a diagnostic upgrade that
/// rides later). `area` is the exact plane-in-box polygon area in physical
/// m² (of the local-box model).
#[derive(Debug, Clone, PartialEq)]
pub struct PlicPlane {
    pub normal: [f64; 3],
    pub d: f64,
    pub centroid: [f64; 3],
    pub area: f64,
}

/// Place the plane: given the (unnormalized) local classifier-field
/// gradient, the local box dimensions `l`, and the target SOLID fraction,
/// return `(unit local normal, d, area)`. The gradient is normalized, its
/// sub-floor components snapped to exact zero (see [`PLIC_DIM_FLOOR`]) and
/// renormalized so the reported normal and the volume solve agree bit-for-
/// bit; a zero/non-finite gradient refuses loudly.
pub(crate) fn fit_plane(
    grad_local: [f64; 3],
    l: [f64; 3],
    solid_frac: f64,
    cell: usize,
) -> Result<([f64; 3], f64, f64), Geom3dError> {
    let g2 = grad_local[0] * grad_local[0]
        + grad_local[1] * grad_local[1]
        + grad_local[2] * grad_local[2];
    if !g2.is_finite() || g2 == 0.0 {
        return Err(Geom3dError::DegeneratePlicNormal { cell });
    }
    let gn = g2.sqrt();
    let mut n = [grad_local[0] / gn, grad_local[1] / gn, grad_local[2] / gn];
    // Snap flat dimensions, then renormalize so (normal, d) self-agree.
    let s0: f64 = (0..3).map(|i| n[i].abs() * l[i]).sum();
    for i in 0..3 {
        if n[i].abs() * l[i] <= PLIC_DIM_FLOOR * s0 {
            n[i] = 0.0;
        }
    }
    let n2 = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
    if n2 == 0.0 {
        return Err(Geom3dError::DegeneratePlicNormal { cell });
    }
    let nn = n2.sqrt();
    n = [n[0] / nn, n[1] / nn, n[2] / nn];

    // Mirror to non-negative components: with x_i measured from the box
    // corner (mirrored where n_i < 0), n·ξ ≤ d ⇔ |n|·x ≤ d + ½Σ|n_i|L_i.
    let na = [n[0].abs(), n[1].abs(), n[2].abs()];
    let span: f64 = na[0] * l[0] + na[1] * l[1] + na[2] * l[2];
    let v_box = l[0] * l[1] * l[2];
    let target = solid_frac * v_box;

    // Bracket = plane offsets touching the box support corners: [0, span].
    let (mut lo, mut hi) = (0.0f64, span);
    for _ in 0..N_PLIC_BISECT {
        let mid = 0.5 * (lo + hi);
        if box_halfspace_volume(na, l, mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let d_pos = 0.5 * (lo + hi);
    let area = box_halfspace_area(na, l, d_pos);
    Ok((n, d_pos - 0.5 * span, area))
}

/// Split into active dimensions (n_i > 0) and flat ones; returns the
/// active (n, l) pairs, their count k, and the product of flat lengths.
fn reduce(n: [f64; 3], l: [f64; 3]) -> ([f64; 3], [f64; 3], usize, f64) {
    let mut an = [0.0f64; 3];
    let mut al = [0.0f64; 3];
    let mut k = 0usize;
    let mut flat = 1.0f64;
    for i in 0..3 {
        if n[i] > 0.0 {
            an[k] = n[i];
            al[k] = l[i];
            k += 1;
        } else {
            flat *= l[i];
        }
    }
    (an, al, k, flat)
}

/// Volume of `{x ∈ [0, L] : n·x ≤ d}` for n_i ≥ 0 — the Scardovelli–
/// Zaleski corner sum over the 2^k active-corner subsets:
/// `V(d) = flat · (1/(k!·Π n_j)) · Σ_s (−1)^{|s|} max(0, d − Σ_{j∈s} n_j L_j)^k`,
/// with flat dimensions (n_i = 0) contributing their lengths as a factor.
/// (Each relu-power term is one repeated integration of the half-space
/// indicator; the alternating subset sum is the inclusion–exclusion over
/// the box corners.) Verified: unit cube, n = (1,1,1): V(1) = 1/6,
/// V(1.5) = 1/2, V(2) = 5/6, V(3) = 1.
pub(crate) fn box_halfspace_volume(n: [f64; 3], l: [f64; 3], d: f64) -> f64 {
    let (an, al, k, flat) = reduce(n, l);
    if k == 0 {
        return if d >= 0.0 { flat } else { 0.0 };
    }
    let mut denom = 1.0f64;
    let mut kfact = 1.0f64;
    for (j, &nj) in an.iter().enumerate().take(k) {
        denom *= nj;
        kfact *= (j + 1) as f64;
    }
    let mut acc = 0.0f64;
    for s in 0..(1u32 << k) {
        let mut off = 0.0f64;
        let mut sign = 1.0f64;
        for j in 0..k {
            if (s >> j) & 1 == 1 {
                off += an[j] * al[j];
                sign = -sign;
            }
        }
        let t = d - off;
        if t > 0.0 {
            acc += sign * t.powi(k as i32);
        }
    }
    flat * acc / (kfact * denom)
}

/// The exact interface area identity `area = dV/dd` (plane cross-section
/// inside the box; exact when `n` is the UNIT normal, so `d` is physical
/// offset): the analytic derivative of the same corner sum — exponents
/// drop by one, the factorial follows. Verified: unit cube, unit normal
/// (1,1,1)/√3 at the mid-diagonal plane gives the regular-hexagon section
/// area (3√3/2)·(√2/2)² ≈ 1.29904.
pub(crate) fn box_halfspace_area(n: [f64; 3], l: [f64; 3], d: f64) -> f64 {
    let (an, al, k, flat) = reduce(n, l);
    if k == 0 {
        return 0.0;
    }
    let mut denom = 1.0f64;
    let mut kfact = 1.0f64;
    for (j, &nj) in an.iter().enumerate().take(k) {
        denom *= nj;
        if j + 1 < k {
            kfact *= (j + 1) as f64; // (k−1)!
        }
    }
    let mut acc = 0.0f64;
    for s in 0..(1u32 << k) {
        let mut off = 0.0f64;
        let mut sign = 1.0f64;
        for j in 0..k {
            if (s >> j) & 1 == 1 {
                off += an[j] * al[j];
                sign = -sign;
            }
        }
        let t = d - off;
        if t > 0.0 {
            acc += sign * t.powi(k as i32 - 1);
        }
    }
    flat * acc / (kfact * denom)
}
