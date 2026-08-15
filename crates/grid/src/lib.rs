//! Implements the structural core of **FND-2 v0.5** — the world-state grid:
//! the natively cylindrical index space `(i_r, i_θ, i_z)` with exact ring
//! metrics (§3.2), the static-(r,z)-topology brick arena in Morton order
//! with per-field SoA storage and dynamic per-brick azimuthal resolution
//! N_θ (§3.2/§3.9), conservative θ-coarsening/refinement with thermalized-
//! ΔKE accounting (§3.4), the azimuthal symmetry indicator with guard,
//! hysteresis and dwell (§3.4), and fixed-shape deterministic reductions
//! (§3.7).
//!
//! Deferred to their own waves (each is additive to this layout):
//! - the full §3.3 multi-material cell model (needs FND-1 `M` + FND-7 spine;
//!   this session stores registered named `f64` fields — the conserved-state
//!   components register as fields when SOLV-1 lands);
//! - §3.5 uniform-grouping tiles (an optimization: everything materialized
//!   is the correct-but-larger representation);
//! - §3.6 FND-3 voxelization ingest + Löhner refinement (construction here
//!   takes explicit extents/spacings — enough for the conduction anchors);
//! - §3.3(7) dormant PLIC sharp-interface fields (reserved by doc, built
//!   with the liquid-interface wave).
//!
//! Determinism (§3.7): brick traversal is a Morton-sorted `Vec`, never a
//! hash map; reductions use a fixed brick-chunk decomposition with a
//! fixed-shape pairwise tree combine, so totals are bit-identical no matter
//! which thread computed which partial.

mod theta;

pub use theta::{
    CollapseLog, MomentumFields, N_DWELL, N_SYM_CADENCE, TAU_COLLAPSE, TAU_EXPAND, ThetaAction,
    ThetaController,
};

/// Leaf-brick edge length in the (r,z) plane (§3.2: small dense bricks,
/// e.g. 8×8). 8×8 = 64 cells ⇒ the per-brick active mask is one `u64`.
pub const BRICK: usize = 8;
pub const BRICK_CELLS: usize = BRICK * BRICK;

/// Coarsest adaptive azimuthal resolution (§3.4, S4): 4 is the coarsest
/// ring carrying both phases of the first symmetry-breaking mode m = 1.
/// `N_θ = 1` is reachable only via the recorded axisymmetry assertion.
pub const N_THETA_GUARD: u32 = 4;

#[derive(Debug, Clone, PartialEq)]
pub struct GridSpec {
    /// Inner radius of the innermost ring [m]; 0 puts the axis in the grid
    /// and activates the §3.2 reflecting-axis treatment.
    pub r_min: f64,
    pub dr: f64,
    pub n_r: usize,
    pub z_min: f64,
    pub dz: f64,
    pub n_z: usize,
    /// Config-declared finest azimuthal resolution; must be `4·2^n`
    /// (FND-4 §3.4-5b θ-ladder) or `1` under the axisymmetry assertion.
    pub n_theta_max: u32,
    /// §3.4: full N_θ = 1 is never an adaptive decision — only a recorded,
    /// pedigree-visible config assertion (carried into the results bundle).
    pub axisymmetry_assertion: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GridError {
    BadSpec(String),
    UnknownField(String),
    BadThetaResolution {
        requested: u32,
        reason: String,
    },
    /// ΔKE accounting needs ρ > 0 in every merged cell (META-1 P6: halt,
    /// never guess past non-physical state).
    NonPositiveDensity {
        brick: usize,
    },
}

impl std::fmt::Display for GridError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BadSpec(m) => write!(f, "bad grid spec: {m}"),
            Self::UnknownField(n) => write!(f, "unknown field {n:?}"),
            Self::BadThetaResolution { requested, reason } => {
                write!(f, "N_θ = {requested} refused: {reason}")
            }
            Self::NonPositiveDensity { brick } => {
                write!(
                    f,
                    "brick {brick}: non-positive density in ΔKE accounting — halt"
                )
            }
        }
    }
}

impl std::error::Error for GridError {}

/// Handle to a registered SoA field (index into every brick's data arrays).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldId(pub(crate) usize);

/// One 8×8 (r,z) leaf brick carrying its azimuthal ring as `n_theta`
/// contiguous θ-planes (§3.2). Storage per field: `n_theta × 64` values,
/// θ-plane-major (`idx = i_theta·64 + local_rz`), SoA across fields (§3.9).
#[derive(Debug, Clone, PartialEq)]
pub struct Brick {
    pub br: u32,
    pub bz: u32,
    pub morton: u64,
    /// Active-cell mask over the 8×8 (r,z) plane (bit `lr·8 + lz`).
    pub mask: u64,
    pub n_theta: u32,
    /// §3.4: geometry-driven floor from FND-3 (defaults to the guard);
    /// flow-adaptive collapse never goes below it.
    pub n_theta_geom_floor: u32,
    pub(crate) data: Vec<Vec<f64>>,
}

impl Brick {
    #[inline]
    pub fn field(&self, f: FieldId) -> &[f64] {
        &self.data[f.0]
    }

    #[inline]
    pub fn field_mut(&mut self, f: FieldId) -> &mut [f64] {
        &mut self.data[f.0]
    }

    #[inline]
    pub fn cell_index(&self, i_theta: u32, local_rz: usize) -> usize {
        i_theta as usize * BRICK_CELLS + local_rz
    }
}

/// Morton interleave of two 32-bit coordinates (§3.2 deterministic layout).
pub fn morton2(x: u32, y: u32) -> u64 {
    fn part1by1(v: u32) -> u64 {
        let mut v = u64::from(v);
        v = (v | (v << 16)) & 0x0000_ffff_0000_ffff;
        v = (v | (v << 8)) & 0x00ff_00ff_00ff_00ff;
        v = (v | (v << 4)) & 0x0f0f_0f0f_0f0f_0f0f;
        v = (v | (v << 2)) & 0x3333_3333_3333_3333;
        v = (v | (v << 1)) & 0x5555_5555_5555_5555;
        v
    }
    part1by1(x) | (part1by1(y) << 1)
}

/// Is `n` on the θ-ladder `{4·2^k}` (FND-2 §3.4 / FND-4 §3.4-5b)?
pub fn theta_ladder_aligned(n: u32) -> bool {
    n >= N_THETA_GUARD && n.is_multiple_of(N_THETA_GUARD) && (n / N_THETA_GUARD).is_power_of_two()
}

#[derive(Debug, Clone, PartialEq)]
pub struct Grid {
    pub spec: GridSpec,
    field_names: Vec<String>,
    /// Morton-sorted active-brick array — THE canonical traversal (§3.7).
    pub bricks: Vec<Brick>,
}

impl Grid {
    /// Build a fully-active grid over the spec'd extents with every brick at
    /// the finest resolution `n_theta_max`. (§3.6 ingest of FND-3 fractions
    /// and refinement land with the voxelization wave; edge bricks get
    /// partial masks.) Field set is fixed at construction — config-time.
    pub fn build(spec: GridSpec, field_names: &[&str]) -> Result<Grid, GridError> {
        // NaN checked explicitly: `< 0.0` alone would let NaN through.
        if spec.r_min.is_nan() || spec.r_min < 0.0 || spec.dr <= 0.0 || spec.dz <= 0.0 {
            return Err(GridError::BadSpec("need r_min ≥ 0, dr > 0, dz > 0".into()));
        }
        if spec.n_r == 0 || spec.n_z == 0 {
            return Err(GridError::BadSpec("need n_r ≥ 1, n_z ≥ 1".into()));
        }
        if spec.n_theta_max == 1 {
            if !spec.axisymmetry_assertion {
                return Err(GridError::BadThetaResolution {
                    requested: 1,
                    reason: "N_θ = 1 requires the recorded axisymmetry assertion (§3.4, S4)".into(),
                });
            }
        } else if !theta_ladder_aligned(spec.n_theta_max) {
            return Err(GridError::BadThetaResolution {
                requested: spec.n_theta_max,
                reason: format!(
                    "must be 4·2^n so coarsening lands on N_θ^guard = {N_THETA_GUARD} \
                     (FND-2 §3.4), or 1 under the axisymmetry assertion"
                ),
            });
        }
        let n_fields = field_names.len();
        let (nbr, nbz) = (spec.n_r.div_ceil(BRICK), spec.n_z.div_ceil(BRICK));
        let mut bricks = Vec::with_capacity(nbr * nbz);
        for br in 0..nbr as u32 {
            for bz in 0..nbz as u32 {
                let mut mask = 0u64;
                for lr in 0..BRICK {
                    for lz in 0..BRICK {
                        let (ir, iz) = (br as usize * BRICK + lr, bz as usize * BRICK + lz);
                        if ir < spec.n_r && iz < spec.n_z {
                            mask |= 1u64 << (lr * BRICK + lz);
                        }
                    }
                }
                if mask == 0 {
                    continue;
                }
                let len = spec.n_theta_max as usize * BRICK_CELLS;
                bricks.push(Brick {
                    br,
                    bz,
                    morton: morton2(br, bz),
                    mask,
                    n_theta: spec.n_theta_max,
                    n_theta_geom_floor: N_THETA_GUARD.min(spec.n_theta_max),
                    data: vec![vec![0.0; len]; n_fields],
                });
            }
        }
        bricks.sort_by_key(|b| b.morton);
        Ok(Grid {
            spec,
            field_names: field_names.iter().map(|s| s.to_string()).collect(),
            bricks,
        })
    }

    pub fn field_id(&self, name: &str) -> Result<FieldId, GridError> {
        self.field_names
            .iter()
            .position(|n| n == name)
            .map(FieldId)
            .ok_or_else(|| GridError::UnknownField(name.to_string()))
    }

    // --- Cylindrical metric (§3.2): exact ring factors, never Cartesian ---

    /// Inner/outer radii of ring `i_r`.
    #[inline]
    pub fn ring_radii(&self, i_r: usize) -> (f64, f64) {
        let r_i = self.spec.r_min + i_r as f64 * self.spec.dr;
        (r_i, r_i + self.spec.dr)
    }

    /// Exact cell volume `½(r_o²−r_i²)·Δθ·Δz` for a ring at resolution
    /// `n_theta` (∝ r̄).
    #[inline]
    pub fn cell_volume(&self, i_r: usize, n_theta: u32) -> f64 {
        let (r_i, r_o) = self.ring_radii(i_r);
        let dtheta = std::f64::consts::TAU / f64::from(n_theta);
        0.5 * (r_o * r_o - r_i * r_i) * dtheta * self.spec.dz
    }

    /// Radial-face area `r·Δθ·Δz` at the inner (`false`) or outer (`true`)
    /// face. At `r = 0` this is exactly zero: the axis face drops out of any
    /// flux stencil geometrically (§3.2 reflecting-axis treatment).
    #[inline]
    pub fn face_area_r(&self, i_r: usize, outer: bool, n_theta: u32) -> f64 {
        let (r_i, r_o) = self.ring_radii(i_r);
        let r = if outer { r_o } else { r_i };
        r * (std::f64::consts::TAU / f64::from(n_theta)) * self.spec.dz
    }

    /// θ-face area `Δr·Δz` (a constant-θ plane section).
    #[inline]
    pub fn face_area_theta(&self) -> f64 {
        self.spec.dr * self.spec.dz
    }

    /// z-face area `½(r_o²−r_i²)·Δθ` (annular sector).
    #[inline]
    pub fn face_area_z(&self, i_r: usize, n_theta: u32) -> f64 {
        let (r_i, r_o) = self.ring_radii(i_r);
        0.5 * (r_o * r_o - r_i * r_i) * (std::f64::consts::TAU / f64::from(n_theta))
    }

    /// §3.2 axis parity pairing: the cross-axis neighbor of θ-index `i` in
    /// the innermost ring is θ+π.
    #[inline]
    pub fn axis_pair(i_theta: u32, n_theta: u32) -> u32 {
        (i_theta + n_theta / 2) % n_theta
    }

    // --- Deterministic access & reduction (§3.7) --------------------------

    /// Brick containing (i_r, i_z), by binary search on the Morton key.
    pub fn brick_index(&self, i_r: usize, i_z: usize) -> Option<usize> {
        let key = morton2((i_r / BRICK) as u32, (i_z / BRICK) as u32);
        self.bricks.binary_search_by_key(&key, |b| b.morton).ok()
    }

    #[inline]
    pub fn local_rz(i_r: usize, i_z: usize) -> usize {
        (i_r % BRICK) * BRICK + (i_z % BRICK)
    }

    /// Volume-weighted global sum `Σ V·q` of a field — the conservation-
    /// ledger reduction. Fixed brick-chunk decomposition (one partial per
    /// brick, Morton order) + fixed-shape pairwise tree combine: the result
    /// is a pure function of the data and the brick array, independent of
    /// which thread computed which partial (§3.7).
    pub fn reduce_volume_weighted(&self, f: FieldId) -> f64 {
        let partials: Vec<f64> = self
            .bricks
            .iter()
            .map(|b| self.brick_partial(b, f))
            .collect();
        tree_combine(&partials)
    }

    fn brick_partial(&self, b: &Brick, f: FieldId) -> f64 {
        let data = b.field(f);
        let mut acc = 0.0f64;
        for i_theta in 0..b.n_theta {
            for local in 0..BRICK_CELLS {
                if b.mask & (1u64 << local) != 0 {
                    let i_r = b.br as usize * BRICK + local / BRICK;
                    let v = self.cell_volume(i_r, b.n_theta);
                    acc += v * data[b.cell_index(i_theta, local)];
                }
            }
        }
        acc
    }
}

/// Fixed-shape pairwise tree combine: reduction shape depends only on
/// `partials.len()`, never on scheduling (§3.7; META-3 `repro-sum`).
pub fn tree_combine(partials: &[f64]) -> f64 {
    match partials.len() {
        0 => 0.0,
        1 => partials[0],
        n => {
            let mid = n / 2;
            tree_combine(&partials[..mid]) + tree_combine(&partials[mid..])
        }
    }
}
