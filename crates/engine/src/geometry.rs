//! Contour-of-revolution geometry (FND-3, stair degenerate form): the cited
//! r(z) station table → the wall function, the ternary region field
//! (gas / liner / exterior) through FND-2 §3.6's ingest seam, and the true
//! wall normal for the slip-ghost walls (station-2 machinery). Partial
//! apertures / cut cells (FND-3 proper) supersede the stair form.

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
