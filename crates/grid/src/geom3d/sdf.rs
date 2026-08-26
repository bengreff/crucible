//! FND-3 §3.1 — the CSG analytic signed-distance tree: primitives with
//! exact signed distances, booleans (union = min, intersection = max,
//! subtraction A−B = max(a, −b)), and the revolved-profile leaf the engine
//! contour corresponds to. Pure data + one `eval`.
//!
//! Sign: negative INSIDE the solid (the module-level convention). Every
//! primitive below is an exact distance except where noted; the booleans
//! keep the **sign** correct everywhere while the **magnitude** is only a
//! bound near concave joins (§3.1 caveat) — which is exactly why occupancy
//! is sampled rather than trusted to the magnitude. Every node is
//! 1-Lipschitz (exact distances are; a 1-Lipschitz coordinate map composed
//! with a 2-D distance is; min/max/negation preserve it), so
//! `|eval(c)| > R_circumscribed` is a sound pure-cell test even where the
//! magnitude under-estimates the true surface distance.

use super::Geom3dError;

/// The analytic SDF tree (FND-3 §3.1). Pure data: build it, [`Sdf::eval`]
/// it — no state, no caches, bit-reproducible.
#[derive(Debug, Clone, PartialEq)]
pub enum Sdf {
    /// Exact: `|p − center| − radius`.
    Sphere { center: [f64; 3], radius: f64 },
    /// Axis-aligned box; exact outside distance (Euclidean to the nearest
    /// face/edge/corner) AND exact inside distance (largest per-axis
    /// signed depth).
    BoxAxis { min: [f64; 3], max: [f64; 3] },
    /// z-aligned capped cylinder about the axis through `center`; exact
    /// both sides via the 2-D (radial, axial) box distance.
    CylinderZ {
        center: [f64; 3],
        radius: f64,
        half_height: f64,
    },
    /// Finite closed solid cone, z-aligned: apex at `apex`, opening
    /// DOWNWARD (−z) over `height`, base radius `height·tan(half_angle)`
    /// (the documented parametrization). Evaluated as the revolved
    /// (s, ζ)-profile triangle `(0,0)–(R,−h)–(0,−h)` about the apex axis:
    /// sign exact; magnitude is the exact surface distance except near the
    /// interior axis segment, where the profile's s = 0 edge makes it a
    /// lower bound — 1-Lipschitz either way, so occupancy and the
    /// pure-cell test stay sound.
    ConeZ {
        apex: [f64; 3],
        half_angle: f64,
        height: f64,
    },
    /// z-axis-normal torus about `center`; exact:
    /// `|(√(x²+y²) − major_r, z)| − minor_r`. Requires `minor_r < major_r`
    /// (no spindle self-intersection).
    Torus {
        center: [f64; 3],
        major_r: f64,
        minor_r: f64,
    },
    /// A CLOSED polygon in the (s, z) half-plane (s = √(x²+y²) ≥ 0, last
    /// vertex implicitly joined to the first) revolved about the GRID
    /// z-axis — the CSG revolved-profile leaf the engine contour
    /// corresponds to. Distance = signed distance to the polygon in
    /// (s, z): minimum distance over edges, sign by the even-odd crossing
    /// rule. Exact for the revolved solid up to the (s, z) metric, which
    /// is exact for revolution about the axis (the nearest surface point
    /// lies in the query's own meridian plane since s, s₀ ≥ 0); profile
    /// edges lying ON the axis (s = 0) revolve to interior degenerate
    /// surfaces and make the magnitude a lower bound there — sign and
    /// Lipschitz bound unaffected.
    Revolved { profile: Vec<(f64, f64)> },
    /// Union = min(a, b) — sign always correct (§3.1).
    Union(Box<Sdf>, Box<Sdf>),
    /// Intersection = max(a, b).
    Intersection(Box<Sdf>, Box<Sdf>),
    /// Subtraction A − B = max(a, −b).
    Subtraction(Box<Sdf>, Box<Sdf>),
}

impl Sdf {
    /// Signed distance (bound) at `p` — negative inside the solid.
    pub fn eval(&self, p: [f64; 3]) -> f64 {
        match self {
            Sdf::Sphere { center, radius } => {
                let d = sub(p, *center);
                norm(d) - radius
            }
            Sdf::BoxAxis { min, max } => {
                // Per-axis signed depth q_i = max(min_i − p_i, p_i − max_i):
                // positive components measure the outside offset (exact
                // Euclidean via their norm); if all are ≤ 0 the point is
                // inside and the distance is the largest (least deep) one.
                let mut out2 = 0.0f64;
                let mut inner = f64::NEG_INFINITY;
                for i in 0..3 {
                    let q = (min[i] - p[i]).max(p[i] - max[i]);
                    if q > 0.0 {
                        out2 += q * q;
                    } else {
                        inner = inner.max(q);
                    }
                }
                if out2 > 0.0 { out2.sqrt() } else { inner }
            }
            Sdf::CylinderZ {
                center,
                radius,
                half_height,
            } => {
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let qr = (dx * dx + dy * dy).sqrt() - radius;
                let qz = (p[2] - center[2]).abs() - half_height;
                let (or_, oz) = (qr.max(0.0), qz.max(0.0));
                let outside = (or_ * or_ + oz * oz).sqrt();
                let inside = qr.max(qz).min(0.0);
                outside + inside
            }
            Sdf::ConeZ {
                apex,
                half_angle,
                height,
            } => {
                let dx = p[0] - apex[0];
                let dy = p[1] - apex[1];
                let s = (dx * dx + dy * dy).sqrt();
                let z = p[2] - apex[2];
                let r_base = height * half_angle.tan();
                let profile = [(0.0, 0.0), (r_base, -height), (0.0, -height)];
                polygon_signed_distance(&profile, s, z)
            }
            Sdf::Torus {
                center,
                major_r,
                minor_r,
            } => {
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let s = (dx * dx + dy * dy).sqrt();
                let a = s - major_r;
                let b = p[2] - center[2];
                (a * a + b * b).sqrt() - minor_r
            }
            Sdf::Revolved { profile } => {
                let s = (p[0] * p[0] + p[1] * p[1]).sqrt();
                polygon_signed_distance(profile, s, p[2])
            }
            Sdf::Union(a, b) => a.eval(p).min(b.eval(p)),
            Sdf::Intersection(a, b) => a.eval(p).max(b.eval(p)),
            Sdf::Subtraction(a, b) => a.eval(p).max(-b.eval(p)),
        }
    }

    /// Refuse non-finite or degenerate parameters (META-1 P6) — recursive
    /// over the tree, called once before any sampling.
    pub fn validate(&self) -> Result<(), Geom3dError> {
        let bad = |m: String| Err(Geom3dError::BadSolid(m));
        match self {
            Sdf::Sphere { center, radius } => {
                if !finite3(*center) || !radius.is_finite() || *radius <= 0.0 {
                    return bad(format!("sphere center {center:?} radius {radius}"));
                }
            }
            Sdf::BoxAxis { min, max } => {
                if !finite3(*min) || !finite3(*max) {
                    return bad(format!("box min {min:?} max {max:?} non-finite"));
                }
                for i in 0..3 {
                    if min[i] >= max[i] {
                        return bad(format!(
                            "box inverted on axis {i}: {} >= {}",
                            min[i], max[i]
                        ));
                    }
                }
            }
            Sdf::CylinderZ {
                center,
                radius,
                half_height,
            } => {
                if !finite3(*center)
                    || !radius.is_finite()
                    || *radius <= 0.0
                    || !half_height.is_finite()
                    || *half_height <= 0.0
                {
                    return bad(format!("cylinder r {radius} hh {half_height}"));
                }
            }
            Sdf::ConeZ {
                apex,
                half_angle,
                height,
            } => {
                if !finite3(*apex)
                    || !half_angle.is_finite()
                    || *half_angle <= 0.0
                    || *half_angle >= std::f64::consts::FRAC_PI_2
                    || !height.is_finite()
                    || *height <= 0.0
                {
                    return bad(format!("cone half_angle {half_angle} height {height}"));
                }
            }
            Sdf::Torus {
                center,
                major_r,
                minor_r,
            } => {
                if !finite3(*center)
                    || !major_r.is_finite()
                    || !minor_r.is_finite()
                    || *minor_r <= 0.0
                    || *minor_r >= *major_r
                {
                    return bad(format!("torus major {major_r} minor {minor_r}"));
                }
            }
            Sdf::Revolved { profile } => {
                if profile.len() < 3 {
                    return bad(format!(
                        "revolved profile has {} < 3 vertices",
                        profile.len()
                    ));
                }
                for &(s, z) in profile {
                    if !s.is_finite() || !z.is_finite() || s < 0.0 {
                        return bad(format!("revolved profile vertex ({s}, {z})"));
                    }
                }
                for i in 0..profile.len() {
                    let a = profile[i];
                    let b = profile[(i + 1) % profile.len()];
                    if a == b {
                        return bad(format!("revolved profile zero-length edge at vertex {i}"));
                    }
                }
            }
            Sdf::Union(a, b) | Sdf::Intersection(a, b) | Sdf::Subtraction(a, b) => {
                a.validate()?;
                b.validate()?;
            }
        }
        Ok(())
    }
}

/// Signed distance to a closed polygon in a 2-D plane: min distance over
/// edges, sign by the even-odd crossing rule (a +s ray; strict comparisons
/// = the fixed tie-break of §3.5). Negative inside.
fn polygon_signed_distance(poly: &[(f64, f64)], s: f64, z: f64) -> f64 {
    let n = poly.len();
    let mut d2 = f64::INFINITY;
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (sa, za) = poly[j];
        let (sb, zb) = poly[i];
        // Distance to the segment (sa,za)–(sb,zb).
        let (ex, ey) = (sb - sa, zb - za);
        let (px, py) = (s - sa, z - za);
        let t = ((px * ex + py * ey) / (ex * ex + ey * ey)).clamp(0.0, 1.0);
        let (dx, dy) = (px - t * ex, py - t * ey);
        d2 = d2.min(dx * dx + dy * dy);
        // Even-odd crossing of the horizontal ray toward +s.
        if (za > z) != (zb > z) {
            let s_cross = sa + (sb - sa) * (z - za) / (zb - za);
            if s < s_cross {
                inside = !inside;
            }
        }
        j = i;
    }
    let d = d2.sqrt();
    if inside { -d } else { d }
}

#[inline]
fn finite3(v: [f64; 3]) -> bool {
    v[0].is_finite() && v[1].is_finite() && v[2].is_finite()
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
