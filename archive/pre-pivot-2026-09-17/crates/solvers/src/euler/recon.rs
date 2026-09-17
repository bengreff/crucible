//! PPM reconstruction (SOLV-1 §3.2; META-3 `castro-source` — Colella &
//! Woodward 1984 interface values + monotonization) on a 1-D pencil of
//! primitive states, method-of-lines form: the limited parabola's edge
//! values become the face input states for the Riemann solve. No
//! characteristic tracing — that is the CTU predictor of the split scheme;
//! on the MOL path (which is also Castro's SDC path, the COUP-3 target)
//! reconstruction feeds the flux directly.
//!
//! Pencil layout: `w[0..3]` and `w[n+3..n+6]` are ghost cells (3 each side),
//! interior cells at `w[3..n+3]`. Face `f ∈ 0..=n` sits between pencil
//! indices `f+2` and `f+3`; `faces_l[f]`/`faces_r[f]` are the states
//! immediately left/right of that face.

use super::{NPRIM, Prim};

/// Ghost cells required on each end of a pencil (PPM edge values of the
/// first ghost cell need two more neighbors beyond it).
pub const NGHOST: usize = 3;

/// Monotonized-central (van Leer) slope of one scalar (CW84 eq. 1.8).
#[inline]
fn mc_slope(wm: f64, w0: f64, wp: f64) -> f64 {
    let dl = w0 - wm;
    let dr = wp - w0;
    if dl * dr <= 0.0 {
        return 0.0;
    }
    let dc = 0.5 * (dl + dr);
    dc.abs().min(2.0 * dl.abs()).min(2.0 * dr.abs()) * dc.signum()
}

/// PPM face states for every component of a pencil. `w` has `n + 2·NGHOST`
/// entries; writes `n + 1` entries into `faces_l`/`faces_r`.
pub fn ppm_faces(w: &[Prim], n: usize, faces_l: &mut [Prim], faces_r: &mut [Prim]) {
    debug_assert_eq!(w.len(), n + 2 * NGHOST);
    debug_assert!(faces_l.len() > n && faces_r.len() > n);

    // Slopes for pencil indices 1..n+5, interface values at i+1/2 for
    // i in 1..n+4, edge pairs for cells 2..n+4 — everything face 0..=n needs.
    let m = w.len();
    let mut slope = vec![[0.0f64; NPRIM]; m];
    for i in 1..m - 1 {
        for k in 0..NPRIM {
            slope[i][k] = mc_slope(w[i - 1][k], w[i][k], w[i + 1][k]);
        }
    }
    // iface[i] = value at interface i+1/2 (CW84 eq. 1.6).
    let mut iface = vec![[0.0f64; NPRIM]; m - 1];
    for i in 1..m - 2 {
        for k in 0..NPRIM {
            iface[i][k] = 0.5 * (w[i][k] + w[i + 1][k]) - (slope[i + 1][k] - slope[i][k]) / 6.0;
        }
    }
    // Monotonized left/right edge values per cell (CW84 eq. 1.10), then the
    // MOL face states: face f takes cell f+2's right edge and cell f+3's
    // left edge.
    let edge = |i: usize| -> ([f64; NPRIM], [f64; NPRIM]) {
        let mut lo = iface[i - 1];
        let mut hi = iface[i];
        for k in 0..NPRIM {
            let c = w[i][k];
            if (hi[k] - c) * (c - lo[k]) <= 0.0 {
                lo[k] = c;
                hi[k] = c;
            } else {
                let d = hi[k] - lo[k];
                let six = 6.0 * (c - 0.5 * (lo[k] + hi[k]));
                if d * six > d * d {
                    lo[k] = 3.0 * c - 2.0 * hi[k];
                } else if d * six < -(d * d) {
                    hi[k] = 3.0 * c - 2.0 * lo[k];
                }
            }
        }
        (lo, hi)
    };
    for f in 0..=n {
        let (_, hi) = edge(f + 2);
        let (lo, _) = edge(f + 3);
        faces_l[f] = hi;
        faces_r[f] = lo;
    }
}

// Covers SOLV-1 §3.2 (PPM properties) — named per META-2 §4.
#[cfg(test)]
mod tests {
    use super::*;

    fn pencil(vals: &[f64]) -> Vec<Prim> {
        vals.iter().map(|&v| [v; NPRIM]).collect()
    }

    #[test]
    fn solv1_s32_ppm_preserves_a_constant_exactly() {
        let n = 6;
        let w = pencil(&[7.5; 12]);
        let mut l = vec![[0.0; NPRIM]; n + 1];
        let mut r = vec![[0.0; NPRIM]; n + 1];
        ppm_faces(&w, n, &mut l, &mut r);
        for f in 0..=n {
            assert_eq!(l[f][0], 7.5);
            assert_eq!(r[f][0], 7.5);
        }
    }

    #[test]
    fn solv1_s32_ppm_reproduces_a_linear_profile() {
        // On linear data every limiter is inactive and the parabola is the
        // line: face values equal the exact interface values, from both
        // sides (no jump ⇒ the Riemann solve sees a smooth field).
        let n = 6;
        let vals: Vec<f64> = (0..n + 6).map(|i| 2.0 + 0.5 * i as f64).collect();
        let w = pencil(&vals);
        let mut l = vec![[0.0; NPRIM]; n + 1];
        let mut r = vec![[0.0; NPRIM]; n + 1];
        ppm_faces(&w, n, &mut l, &mut r);
        for f in 0..=n {
            let exact = 2.0 + 0.5 * (f as f64 + 2.5);
            assert!((l[f][0] - exact).abs() < 1e-13, "face {f}: {}", l[f][0]);
            assert!((r[f][0] - exact).abs() < 1e-13);
        }
    }

    #[test]
    fn solv1_s32_ppm_face_states_stay_within_neighbor_range() {
        // Monotonicity: around a jump, no face state may overshoot the
        // adjacent cell values (this is what keeps reconstructed ρ, p
        // positive when the cells are).
        let n = 6;
        let w = pencil(&[
            1.0, 1.0, 1.0, 1.0, 1.0, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125,
        ]);
        let mut l = vec![[0.0; NPRIM]; n + 1];
        let mut r = vec![[0.0; NPRIM]; n + 1];
        ppm_faces(&w, n, &mut l, &mut r);
        for f in 0..=n {
            for v in [l[f][0], r[f][0]] {
                assert!((0.125..=1.0).contains(&v), "face {f} overshoots: {v}");
            }
        }
    }
}
