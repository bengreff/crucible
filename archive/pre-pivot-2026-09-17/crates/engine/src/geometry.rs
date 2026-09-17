//! Contour-of-revolution geometry (FND-3): the cited r(z) station table →
//! the wall function, the ternary region field (gas / liner / exterior)
//! through FND-2 §3.6's ingest seam, the true wall normal for the
//! slip-ghost walls (station-2 machinery), and — session 12 — the **FND-3
//! §3.3 partial fractions and face apertures**, computed by ANALYTIC
//! clipping of the revolved piecewise-linear profile against each ring
//! cell in the cylindrical volume measure. A revolved profile is the CSG
//! revolved-profile leaf on FND-3's analytic path: the clipped integrals
//! are closed-form (piecewise polynomial), so the emitted fractions carry
//! **zero sampling error** (the §3.1 jittered-sampling bound degenerates
//! to round-off; no `N_FRAC_SAMPLES` is consumed) and are exactly
//! deterministic. Ties (wall touching a face) resolve by strict `>` —
//! fixed tie-breaking per §3.5.

use crucible_grid::Region;

/// A parsed, SI, strictly-z-increasing contour with its liner declaration.
#[derive(Debug, Clone)]
pub struct Contour {
    /// (z, r) stations in metres, z strictly increasing.
    stations: Vec<(f64, f64)>,
    pub liner_thickness_m: f64,
    pub r_throat_m: f64,
    /// Axial location of the minimum radius (the throat plane).
    pub z_throat_m: f64,
}

impl Contour {
    /// Build from stations already parsed by the config loader
    /// (`crucible_config::parse_contour_csv`) — one parser, two consumers.
    pub fn new(stations: Vec<(f64, f64)>, liner_thickness_m: f64) -> Result<Self, String> {
        if stations.len() < 2 {
            return Err("a contour needs at least 2 stations".to_string());
        }
        let (mut r_throat, mut z_throat) = (f64::INFINITY, 0.0);
        for &(z, r) in &stations {
            if r < r_throat {
                r_throat = r;
                z_throat = z;
            }
        }
        Ok(Self {
            stations,
            liner_thickness_m,
            r_throat_m: r_throat,
            z_throat_m: z_throat,
        })
    }

    pub fn z_min(&self) -> f64 {
        self.stations[0].0
    }

    pub fn z_max(&self) -> f64 {
        self.stations[self.stations.len() - 1].0
    }

    /// Wall radius at `z` — piecewise-linear between stations (the declared
    /// interpretation rule of the geometry-of-record CSV; a smoother
    /// reconstruction must not overshoot station radii). Outside the axial
    /// span the end radii extend (only ghost centroids ever ask).
    pub fn r_wall(&self, z: f64) -> f64 {
        let s = &self.stations;
        if z <= s[0].0 {
            return s[0].1;
        }
        if z >= s[s.len() - 1].0 {
            return s[s.len() - 1].1;
        }
        // partition_point: first station with z_station > z.
        let hi = s.partition_point(|&(zs, _)| zs <= z).min(s.len() - 1);
        let (z0, r0) = s[hi - 1];
        let (z1, r1) = s[hi];
        r0 + (r1 - r0) * (z - z0) / (z1 - z0)
    }

    /// Local wall slope dr/dz (piecewise-constant between stations).
    fn slope(&self, z: f64) -> f64 {
        let s = &self.stations;
        let hi = s.partition_point(|&(zs, _)| zs <= z).clamp(1, s.len() - 1);
        let (z0, r0) = s[hi - 1];
        let (z1, r1) = s[hi];
        (r1 - r0) / (z1 - z0)
    }

    /// Unit wall normal `(n_r, n_z)` of the true (smooth) wall at `(r, z)` —
    /// the slip-ghost operand (`Euler::wall_normal`). The wall tangent is
    /// `(dr, dz) = (r'(z), 1)`; any ± of the normal serves (reflection is
    /// symmetric in n̂).
    pub fn wall_normal(&self, _r: f64, z: f64) -> (f64, f64) {
        let rp = self.slope(z);
        let inv = 1.0 / (1.0 + rp * rp).sqrt();
        (inv, -rp * inv)
    }

    /// FND-2 §3.6 region classification at a cell center: gas inside the
    /// contour, liner in the declared ring outside it, exterior beyond.
    /// `liner_thickness_m = 0` ⇒ no solid anywhere (cold-flow geometry).
    pub fn classify(&self, r: f64, z: f64) -> Region {
        let rw = self.r_wall(z);
        if r < rw {
            Region::Gas
        } else if self.liner_thickness_m > 0.0 && r < rw + self.liner_thickness_m {
            Region::Solid
        } else {
            Region::Exterior
        }
    }

    /// Subdivide `[z0, z1]` at interior stations: returns the ordered
    /// breakpoints (z0, …, z1) between which `r_wall` is exactly linear.
    fn linear_spans(&self, z0: f64, z1: f64) -> Vec<f64> {
        let mut cuts = vec![z0];
        for &(zs, _) in &self.stations {
            if zs > z0 && zs < z1 {
                cuts.push(zs);
            }
        }
        cuts.push(z1);
        cuts
    }

    /// FND-3 §3.3 — exact gas volume fraction κ of the ring cell
    /// `[r0, r1] × [z0, z1]` in the cylindrical measure:
    /// `κ = ∫ ½(clamp(r_w(z), r0, r1)² − r0²) dz / (½(r1²−r0²)(z1−z0))`.
    /// Piecewise-linear `r_w` ⇒ the integrand is piecewise quadratic and
    /// the integral closed-form: on a span where `w` is linear and inside
    /// `[r0, r1]`, `∫w²dz = len·(wa² + wa·wb + wb²)/3` (exact, no division
    /// by the slope). Clamp crossings are found per linear span.
    pub fn gas_volume_fraction(&self, r0: f64, r1: f64, z0: f64, z1: f64) -> f64 {
        debug_assert!(r1 > r0 && z1 > z0);
        let mut acc = 0.0f64;
        let cuts = self.linear_spans(z0, z1);
        for pair in cuts.windows(2) {
            let (za, zb) = (pair[0], pair[1]);
            let (wa, wb) = (self.r_wall(za), self.r_wall(zb));
            // Split [za, zb] at the (at most two) points where the linear
            // w crosses r0 or r1, then integrate each piece by regime.
            let mut pts = vec![za];
            for rc in [r0, r1] {
                if (wa - rc) * (wb - rc) < 0.0 {
                    pts.push(za + (zb - za) * (rc - wa) / (wb - wa));
                }
            }
            pts.push(zb);
            pts.sort_by(|a, b| a.partial_cmp(b).expect("finite z"));
            for piece in pts.windows(2) {
                let (pa, pb) = (piece[0], piece[1]);
                let len = pb - pa;
                if len <= 0.0 {
                    continue;
                }
                let wm = self.r_wall(0.5 * (pa + pb));
                if wm <= r0 {
                    // fully covered piece
                } else if wm >= r1 {
                    acc += 0.5 * (r1 * r1 - r0 * r0) * len;
                } else {
                    let (va, vb) = (self.r_wall(pa), self.r_wall(pb));
                    acc += 0.5 * (len * (va * va + va * vb + vb * vb) / 3.0 - r0 * r0 * len);
                }
            }
        }
        (acc / (0.5 * (r1 * r1 - r0 * r0) * (z1 - z0))).clamp(0.0, 1.0)
    }

    /// FND-3 §3.3 — open-area fraction of the radial face at radius `r_f`
    /// spanning `[z0, z1]`: the length fraction where `r_wall(z) > r_f`.
    pub fn r_face_aperture(&self, r_f: f64, z0: f64, z1: f64) -> f64 {
        let mut open = 0.0f64;
        let cuts = self.linear_spans(z0, z1);
        for pair in cuts.windows(2) {
            let (za, zb) = (pair[0], pair[1]);
            let (wa, wb) = (self.r_wall(za), self.r_wall(zb));
            let len = zb - za;
            open += if wa > r_f && wb > r_f {
                len
            } else if wa <= r_f && wb <= r_f {
                0.0
            } else {
                let t = (r_f - wa) / (wb - wa);
                if wa > r_f { t * len } else { (1.0 - t) * len }
            };
        }
        (open / (z1 - z0)).clamp(0.0, 1.0)
    }

    /// FND-3 §3.3 — open-area fraction of the z-face annulus `[r0, r1]` at
    /// `z_f` (annular measure): `(clamp(r_w, r0, r1)² − r0²)/(r1² − r0²)`.
    pub fn z_face_aperture(&self, r0: f64, r1: f64, z_f: f64) -> f64 {
        let x = self.r_wall(z_f).clamp(r0, r1);
        ((x * x - r0 * r0) / (r1 * r1 - r0 * r0)).clamp(0.0, 1.0)
    }

    /// Geometric throat area πr_t² (the c\*/C_F reference, SOLV-7 §3.3).
    pub fn throat_area(&self) -> f64 {
        std::f64::consts::PI * self.r_throat_m * self.r_throat_m
    }

    /// Injector-plane (z_min) flow area πr_inj² — the ṁ/A of the declared
    /// injector mass flow.
    pub fn injector_area(&self) -> f64 {
        let r = self.stations[0].1;
        std::f64::consts::PI * r * r
    }
}
