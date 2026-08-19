//! Implements the structural core of **FND-2 v0.5** — the world-state grid:
//! the natively cylindrical index space `(i_r, i_θ, i_z)` with exact ring
//! metrics (§3.2), the static-(r,z)-topology brick arena in Morton order
//! with per-field SoA storage and dynamic per-brick azimuthal resolution
//! N_θ (§3.2/§3.9), conservative θ-coarsening/refinement with thermalized-
//! ΔKE accounting (§3.4), the azimuthal symmetry indicator with guard,
//! hysteresis and dwell (§3.4), and fixed-shape deterministic reductions
//! (§3.7).
//!
//! **Sealed surface (post-review):** brick storage, Morton order, masks and
//! per-brick N_θ are private; reads go through `bricks()`/`brick()` and the
//! `for_each_active_cell` visitor, writes through `brick_field_mut`/
//! `fill_field`. The §3.7 determinism invariants (sorted Morton array,
//! mask ⇔ in-bounds, `n_theta` ⇔ data length) hold by construction because
//! no external code can mutate the structures that carry them.
//!
//! Deferred to their own waves (each is additive to this layout):
//! - the full §3.3 multi-material cell model (needs FND-1 `M` + FND-7 spine;
//!   this session stores registered named `f64` fields);
//! - §3.5 uniform-grouping tiles (an optimization);
//! - §3.6 FND-3 voxelization ingest + Löhner refinement;
//! - §3.3(7) dormant PLIC sharp-interface fields.

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
    /// Non-finite or non-positive density in ΔKE accounting (META-1 P6:
    /// halt, never guess past non-physical state).
    NonPositiveDensity {
        brick: usize,
    },
    /// A computation encountered NaN/inf where the doctrine demands a halt
    /// with diagnosis (e.g. the symmetry indicator over corrupted data).
    NonFinite {
        context: String,
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
                    "brick {brick}: non-finite or non-positive density in ΔKE accounting — halt"
                )
            }
            Self::NonFinite { context } => {
                write!(
                    f,
                    "non-finite value in {context} — halt with diagnosis (META-1 P6)"
                )
            }
        }
    }
}

impl std::error::Error for GridError {}

/// Handle to a registered SoA field (index into every brick's data arrays).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldId(pub(crate) usize);

/// Per-cell medium region — the §3.6 ingest surface in ternary degenerate
/// form (FND-3 fractions/apertures refine this same seam later). `Gas` is
/// the flow-active domain (the mask every existing consumer reads);
/// `Solid` cells carry conductive media (the station-4 liner); `Exterior`
/// cells are outside the world. The classification is config-time data —
/// an operator's domain comes from here, never from an `if(material)`
/// (META-1 Rule 12).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Region {
    Exterior,
    Gas,
    Solid,
}

/// Direction from a cell to a face-adjacent neighbor, in the fixed sweep
/// order (r−, r+, z−, z+) every deterministic enumeration uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceDir {
    RMinus,
    RPlus,
    ZMinus,
    ZPlus,
}

impl FaceDir {
    /// Index into per-cell face arrays, in the fixed enumeration order.
    #[inline]
    pub const fn index(self) -> usize {
        match self {
            FaceDir::RMinus => 0,
            FaceDir::RPlus => 1,
            FaceDir::ZMinus => 2,
            FaceDir::ZPlus => 3,
        }
    }
}

/// FND-3 §3.3 cut geometry of one cell, as ingested through the §3.6 seam:
/// the gas volume fraction κ (cylindrical measure) and the open-area
/// fraction of each face in [`FaceDir`] order. The supplier must compute a
/// shared face's aperture from the face coordinates alone so both cells
/// carry the identical bits — validated at build.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellGeom {
    pub region: Region,
    pub kappa: f64,
    pub aperture: [f64; 4],
}

/// Per-brick cut geometry ((r,z)-shaped, θ-uniform, like the masks):
/// κ and the four face apertures per local cell. Present on every brick of
/// a grid built through [`Grid::build_with_geometry`]; absent (`None`) on
/// stair/box worlds, whose accessors degenerate to κ = 1 / aperture = 1 —
/// the arithmetic identity that keeps every certified full-box operator
/// bit-identical.
#[derive(Debug, Clone, PartialEq)]
pub struct BrickGeom {
    pub kappa: Vec<f64>,
    /// Indexed `[FaceDir::index()][local_rz]`.
    pub aperture: [Vec<f64>; 4],
}

/// One config-time-identified gas↔solid interface face (SOLV-1 §3.5's
/// "wall faces identified geometrically from config data"). `dir` points
/// from the gas cell toward the solid cell; faces are θ-uniform (regions
/// are (r,z)-shaped, like activity).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InterfaceFace {
    pub gas: (usize, usize),
    pub solid: (usize, usize),
    pub dir: FaceDir,
}

/// One 8×8 (r,z) leaf brick carrying its azimuthal ring as `n_theta`
/// contiguous θ-planes (§3.2). Storage per field: `n_theta × 64` values,
/// θ-plane-major (`idx = i_theta·64 + local_rz`), SoA across fields (§3.9).
/// Structure fields are read-only outside the crate — the Morton order,
/// mask, and `n_theta` ⇔ storage-length invariants are load-bearing (§3.7).
#[derive(Debug, Clone, PartialEq)]
pub struct Brick {
    pub(crate) br: u32,
    pub(crate) bz: u32,
    pub(crate) morton: u64,
    pub(crate) mask: u64,
    pub(crate) solid_mask: u64,
    pub(crate) n_theta: u32,
    pub(crate) n_theta_geom_floor: u32,
    pub(crate) geom: Option<BrickGeom>,
    pub(crate) data: Vec<Vec<f64>>,
}

impl Brick {
    #[inline]
    pub fn br(&self) -> u32 {
        self.br
    }

    #[inline]
    pub fn bz(&self) -> u32 {
        self.bz
    }

    #[inline]
    pub fn morton(&self) -> u64 {
        self.morton
    }

    /// Flow-active (gas) cell mask over the 8×8 (r,z) plane (bit `lr·8 + lz`).
    #[inline]
    pub fn mask(&self) -> u64 {
        self.mask
    }

    /// Solid-region mask over the same plane; disjoint from `mask` by
    /// construction (a cell has exactly one [`Region`]).
    #[inline]
    pub fn solid_mask(&self) -> u64 {
        self.solid_mask
    }

    #[inline]
    pub fn n_theta(&self) -> u32 {
        self.n_theta
    }

    /// §3.4: geometry-driven floor from FND-3 (defaults to the guard).
    #[inline]
    pub fn n_theta_geom_floor(&self) -> u32 {
        self.n_theta_geom_floor
    }

    /// Whether FND-3 cut geometry is attached to this brick.
    #[inline]
    pub fn has_geom(&self) -> bool {
        self.geom.is_some()
    }

    /// Gas volume fraction κ of a local cell (FND-3 §3.3). Full-box worlds
    /// (no geometry attached) answer exactly 1.0 for active cells — the
    /// arithmetic-identity default (×1.0 and /1.0 are exact), so consumers
    /// need no mode branch.
    #[inline]
    pub fn kappa_rz(&self, local: usize) -> f64 {
        self.geom.as_ref().map_or(1.0, |g| g.kappa[local])
    }

    /// Open-area fraction of a local cell's face (FND-3 §3.3); 1.0 when no
    /// geometry is attached (see [`Brick::kappa_rz`]).
    #[inline]
    pub fn aperture_rz(&self, dir: FaceDir, local: usize) -> f64 {
        self.geom
            .as_ref()
            .map_or(1.0, |g| g.aperture[dir.index()][local])
    }

    #[inline]
    pub fn field(&self, f: FieldId) -> &[f64] {
        &self.data[f.0]
    }

    #[inline]
    pub(crate) fn field_mut(&mut self, f: FieldId) -> &mut [f64] {
        &mut self.data[f.0]
    }

    #[inline]
    pub fn cell_index(&self, i_theta: u32, local_rz: usize) -> usize {
        i_theta as usize * BRICK_CELLS + local_rz
    }

    /// Global (i_r, i_z) of a local cell — the single inverse of
    /// `Grid::local_rz` (review finding: this decode was hand-copied at 7
    /// sites; a `/` vs `%` flip silently transposes r and z).
    #[inline]
    pub fn global_rz(&self, local: usize) -> (usize, usize) {
        (
            self.br as usize * BRICK + local / BRICK,
            self.bz as usize * BRICK + local % BRICK,
        )
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

/// One visited cell: indices, flat storage index, and centroid coordinates.
#[derive(Debug, Clone, Copy)]
pub struct CellRef {
    pub bi: usize,
    pub i_r: usize,
    pub i_z: usize,
    pub i_theta: u32,
    pub local: usize,
    /// Flat index into any of this brick's field slices.
    pub idx: usize,
    pub r: f64,
    pub theta: f64,
    pub z: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Grid {
    spec: GridSpec,
    field_names: Vec<String>,
    /// Morton-sorted active-brick array — THE canonical traversal (§3.7).
    bricks: Vec<Brick>,
}

impl Grid {
    /// Build a fully-active grid over the spec'd extents with every brick at
    /// the finest resolution `n_theta_max`. (§3.6 ingest of FND-3 fractions
    /// and refinement land with the voxelization wave; edge bricks get
    /// partial masks.) Field set is fixed at construction — config-time.
    pub fn build(spec: GridSpec, field_names: &[&str]) -> Result<Grid, GridError> {
        Self::build_with_activity(spec, field_names, |_, _| true)
    }

    /// Build with a per-(i_r, i_z) activity predicate — the §3.6 ingest
    /// surface in its binary (occupancy) degenerate form: FND-3's
    /// voxelization will supply partial fractions and apertures through
    /// this same config-time seam; until then a certificate fixture may
    /// activate cells from an analytic contour. Activity is (r,z)-shaped
    /// (uniform in θ — a non-axisymmetric mask is FND-3 §3.4's `N_θ^geom`
    /// floor territory, not this seam). Inactive cells are excluded from
    /// every mask; fully-inactive bricks are not allocated at all.
    pub fn build_with_activity(
        spec: GridSpec,
        field_names: &[&str],
        active: impl Fn(usize, usize) -> bool,
    ) -> Result<Grid, GridError> {
        Self::build_with_regions(spec, field_names, |i_r, i_z| {
            if active(i_r, i_z) {
                Region::Gas
            } else {
                Region::Exterior
            }
        })
    }

    /// Build with a per-(i_r, i_z) [`Region`] classifier — the §3.6 ingest
    /// surface in ternary degenerate form (gas / solid / exterior; FND-3
    /// fractions and apertures arrive through this same config-time seam).
    /// Regions are (r,z)-shaped, uniform in θ. Bricks with no gas and no
    /// solid cells are not allocated.
    pub fn build_with_regions(
        spec: GridSpec,
        field_names: &[&str],
        classify: impl Fn(usize, usize) -> Region,
    ) -> Result<Grid, GridError> {
        // `is_finite` everywhere: `<= 0.0` is false for NaN AND +inf, so
        // comparisons alone admit non-finite worlds (review finding).
        if !spec.r_min.is_finite()
            || spec.r_min < 0.0
            || !spec.dr.is_finite()
            || spec.dr <= 0.0
            || !spec.dz.is_finite()
            || spec.dz <= 0.0
            || !spec.z_min.is_finite()
        {
            return Err(GridError::BadSpec(
                "need finite z_min, finite r_min ≥ 0, finite dr > 0, finite dz > 0".into(),
            ));
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
                let mut solid_mask = 0u64;
                for lr in 0..BRICK {
                    for lz in 0..BRICK {
                        let (ir, iz) = (br as usize * BRICK + lr, bz as usize * BRICK + lz);
                        if ir < spec.n_r && iz < spec.n_z {
                            match classify(ir, iz) {
                                Region::Gas => mask |= 1u64 << (lr * BRICK + lz),
                                Region::Solid => solid_mask |= 1u64 << (lr * BRICK + lz),
                                Region::Exterior => {}
                            }
                        }
                    }
                }
                if mask == 0 && solid_mask == 0 {
                    continue;
                }
                let len = spec.n_theta_max as usize * BRICK_CELLS;
                bricks.push(Brick {
                    br,
                    bz,
                    morton: morton2(br, bz),
                    mask,
                    solid_mask,
                    n_theta: spec.n_theta_max,
                    n_theta_geom_floor: N_THETA_GUARD.min(spec.n_theta_max),
                    geom: None,
                    data: vec![vec![0.0; len]; n_fields],
                });
            }
        }
        // An all-Exterior classifier is a config/geometry mistake (a
        // contour that misses every cell center), not a world — refuse at
        // the build seam rather than letting the first operator panic on
        // `brick(0)` (META-1 P6; review finding).
        if bricks.is_empty() {
            return Err(GridError::BadSpec(
                "region classifier marked every cell Exterior — an empty world cannot build \
                 (check contour/extents units)"
                    .into(),
            ));
        }
        bricks.sort_by_key(|b| b.morton);
        Ok(Grid {
            spec,
            field_names: field_names.iter().map(|s| s.to_string()).collect(),
            bricks,
        })
    }

    /// Build through the full FND-3 §3.3/§3.6 ingest seam: per-cell partial
    /// fractions + face apertures. The invariants a legal voxelization must
    /// satisfy are validated loudly (META-1 P6, FND-3's mask-disjointness
    /// validation class):
    /// - `Gas ⇔ κ > 0`; `Solid`/`Exterior` ⇒ κ = 0;
    /// - κ and apertures finite in [0, 1];
    /// - a shared face carries bitwise-identical apertures on both sides
    ///   (the supplier computes it from face coordinates alone);
    /// - a face adjoining a κ = 0 cell is fully covered (aperture 0).
    ///
    /// Restricted to `n_theta_max = 1`: a cut θ-face needs the azimuthal
    /// fraction field of the sampled voxelization path — that arrives with
    /// the 3-D wave; refuse rather than guess.
    pub fn build_with_geometry(
        spec: GridSpec,
        field_names: &[&str],
        geom_of: impl Fn(usize, usize) -> CellGeom,
    ) -> Result<Grid, GridError> {
        if spec.n_theta_max != 1 {
            return Err(GridError::BadSpec(
                "cut geometry at N_θ > 1 needs the sampled azimuthal-fraction path \
                 (FND-3 3-D wave); refusing rather than guessing"
                    .into(),
            ));
        }
        let mut g =
            Self::build_with_regions(spec, field_names, |i_r, i_z| geom_of(i_r, i_z).region)?;
        for b in &mut g.bricks {
            let mut bg = BrickGeom {
                kappa: vec![0.0; BRICK_CELLS],
                aperture: std::array::from_fn(|_| vec![0.0; BRICK_CELLS]),
            };
            for local in 0..BRICK_CELLS {
                let (i_r, i_z) = b.global_rz(local);
                if i_r >= g.spec.n_r || i_z >= g.spec.n_z {
                    continue;
                }
                let c = geom_of(i_r, i_z);
                let legal = c.kappa.is_finite()
                    && (0.0..=1.0).contains(&c.kappa)
                    && match c.region {
                        Region::Gas => c.kappa > 0.0,
                        Region::Solid | Region::Exterior => c.kappa == 0.0,
                    };
                if !legal {
                    return Err(GridError::BadSpec(format!(
                        "cell ({i_r}, {i_z}): κ = {} inconsistent with region {:?} \
                         (Gas ⇔ κ > 0; κ ∈ [0, 1])",
                        c.kappa, c.region
                    )));
                }
                for (d, &a) in c.aperture.iter().enumerate() {
                    if !a.is_finite() || !(0.0..=1.0).contains(&a) {
                        return Err(GridError::BadSpec(format!(
                            "cell ({i_r}, {i_z}) face {d}: aperture {a} outside [0, 1]"
                        )));
                    }
                    bg.aperture[d][local] = a;
                }
                bg.kappa[local] = c.kappa;
            }
            b.geom = Some(bg);
        }
        // Shared-face coherence: identical bits from both sides; covered
        // against any κ = 0 cell (incl. unallocated bricks / domain-interior
        // exterior cells).
        for i_r in 0..g.spec.n_r {
            for i_z in 0..g.spec.n_z {
                let here = [
                    (
                        FaceDir::RMinus,
                        i_r.checked_sub(1).map(|r| (r, i_z)),
                        FaceDir::RPlus,
                    ),
                    (
                        FaceDir::ZMinus,
                        i_z.checked_sub(1).map(|z| (i_r, z)),
                        FaceDir::ZPlus,
                    ),
                ];
                let a_here = |dir: FaceDir| {
                    g.brick_index(i_r, i_z).map_or(0.0, |bi| {
                        g.bricks[bi].aperture_rz(dir, Self::local_rz(i_r, i_z))
                    })
                };
                let kappa_here = g.kappa(i_r, i_z);
                for (dir, nbr, opp) in here {
                    let a = if g.brick_index(i_r, i_z).is_some() {
                        a_here(dir)
                    } else {
                        0.0
                    };
                    let (a_nbr, kappa_nbr) = match nbr {
                        Some((nr, nz)) => (
                            g.brick_index(nr, nz).map_or(0.0, |bi| {
                                g.bricks[bi].aperture_rz(opp, Self::local_rz(nr, nz))
                            }),
                            g.kappa(nr, nz),
                        ),
                        None => continue, // domain-edge face: no partner
                    };
                    if a.to_bits() != a_nbr.to_bits() {
                        return Err(GridError::BadSpec(format!(
                            "shared face ({i_r}, {i_z})↔{nbr:?}: apertures differ \
                             ({a} vs {a_nbr}) — the supplier must compute a face from its \
                             coordinates alone"
                        )));
                    }
                    if (kappa_here == 0.0 || kappa_nbr == 0.0) && a != 0.0 {
                        return Err(GridError::BadSpec(format!(
                            "face ({i_r}, {i_z})↔{nbr:?}: open aperture {a} against a \
                             κ = 0 cell — geometry incoherent"
                        )));
                    }
                }
            }
        }
        Ok(g)
    }

    /// Whether this world carries FND-3 cut geometry (fractions/apertures).
    pub fn has_cut_geometry(&self) -> bool {
        self.bricks.iter().any(|b| b.geom.is_some())
    }

    /// Gas volume fraction κ of cell (i_r, i_z): 0 outside the gas mask,
    /// 1 for active cells of full-box worlds.
    pub fn kappa(&self, i_r: usize, i_z: usize) -> f64 {
        if !self.is_active(i_r, i_z) {
            return 0.0;
        }
        self.brick_index(i_r, i_z)
            .map_or(0.0, |bi| self.bricks[bi].kappa_rz(Self::local_rz(i_r, i_z)))
    }

    /// Face aperture of cell (i_r, i_z) toward `dir` (1.0 on full-box
    /// worlds' active cells; 0.0 outside the mask).
    pub fn aperture(&self, i_r: usize, i_z: usize, dir: FaceDir) -> f64 {
        if !self.is_active(i_r, i_z) {
            return 0.0;
        }
        self.brick_index(i_r, i_z).map_or(0.0, |bi| {
            self.bricks[bi].aperture_rz(dir, Self::local_rz(i_r, i_z))
        })
    }

    /// The discrete wall-closure vector `(W_r, W_z)` of a cut cell — THE
    /// single owner of the identity that defines the embedded interface
    /// from apertures on the cylindrical metric (SOLV-1 §3.3 well-balance
    /// extended to cut cells):
    ///
    /// `W_r = a_r₊·A_r₊ − a_r₋·A_r₋ − κ·(A_r₊ − A_r₋)`,
    /// `W_z = (a_z₊ − a_z₋)·A_z`.
    ///
    /// The interface pressure force on the gas is `+p·(W_r, W_z)` (so the
    /// momentum rate gains `p·W/(κV)`), chosen so a uniform state at rest
    /// is preserved to round-off — the discrete divergence theorem. `|W|`
    /// is the interface area (exact for the true geometry under uniform
    /// p; the smooth-wall area the stair form overcounts), and `−W/|W|`
    /// the outward (gas→wall) interface normal. Identically (0, 0) for
    /// uncut interior cells.
    pub fn wall_closure(&self, i_r: usize, i_z: usize, n_theta: u32) -> (f64, f64) {
        let Some(bi) = self.brick_index(i_r, i_z) else {
            return (0.0, 0.0);
        };
        let b = &self.bricks[bi];
        let local = Self::local_rz(i_r, i_z);
        let a_in = self.face_area_r(i_r, false, n_theta);
        let a_out = self.face_area_r(i_r, true, n_theta);
        let a_z = self.face_area_z(i_r, n_theta);
        let kappa = b.kappa_rz(local);
        let w_r = b.aperture_rz(FaceDir::RPlus, local) * a_out
            - b.aperture_rz(FaceDir::RMinus, local) * a_in
            - kappa * (a_out - a_in);
        let w_z =
            (b.aperture_rz(FaceDir::ZPlus, local) - b.aperture_rz(FaceDir::ZMinus, local)) * a_z;
        (w_r, w_z)
    }

    #[inline]
    pub fn spec(&self) -> &GridSpec {
        &self.spec
    }

    pub fn field_id(&self, name: &str) -> Result<FieldId, GridError> {
        self.field_names
            .iter()
            .position(|n| n == name)
            .map(FieldId)
            .ok_or_else(|| GridError::UnknownField(name.to_string()))
    }

    // --- Brick access (read-only structure; writes are Grid-mediated) -----

    #[inline]
    pub fn bricks(&self) -> &[Brick] {
        &self.bricks
    }

    #[inline]
    pub fn n_bricks(&self) -> usize {
        self.bricks.len()
    }

    #[inline]
    pub fn brick(&self, bi: usize) -> &Brick {
        &self.bricks[bi]
    }

    /// Mutable access to one field's storage in one brick — the only
    /// external write path (structure stays sealed).
    #[inline]
    pub fn brick_field_mut(&mut self, bi: usize, f: FieldId) -> &mut [f64] {
        self.bricks[bi].field_mut(f)
    }

    /// Two distinct fields of one brick, mutably (e.g. `T += dt·rate`).
    pub fn brick_fields_mut2(
        &mut self,
        bi: usize,
        a: FieldId,
        b: FieldId,
    ) -> (&mut [f64], &mut [f64]) {
        assert_ne!(a.0, b.0, "brick_fields_mut2 requires distinct fields");
        let [x, y] = self.bricks[bi]
            .data
            .get_disjoint_mut([a.0, b.0])
            .expect("distinct in-range field ids");
        (x, y)
    }

    // --- Cylindrical metric (§3.2): exact ring factors, never Cartesian ---

    /// Radius of radial face `f` (the inner face of ring `f`) — the single
    /// owner of the face-radius rounding. Ring `i`'s outer face and ring
    /// `i+1`'s inner face MUST be the same number bitwise, or flux-form
    /// telescoping (COUP-2) and the well-balanced geometric sources
    /// (SOLV-1 §3.3) pick up one-ulp seams at particular radii (found by
    /// the Station-1 uniform-fixed-point certificate: the former
    /// `r_i + dr` rounds differently from `r_min + (i+1)·dr`).
    #[inline]
    pub fn face_radius(&self, f: usize) -> f64 {
        self.spec.r_min + f as f64 * self.spec.dr
    }

    /// Inner/outer radii of ring `i_r`.
    #[inline]
    pub fn ring_radii(&self, i_r: usize) -> (f64, f64) {
        (self.face_radius(i_r), self.face_radius(i_r + 1))
    }

    /// Cell-center coordinates — the single owner of the centroid
    /// convention (review finding: this arithmetic was re-derived at 8
    /// sites; solver and verifier must sample identical points).
    #[inline]
    pub fn r_center(&self, i_r: usize) -> f64 {
        self.spec.r_min + (i_r as f64 + 0.5) * self.spec.dr
    }

    #[inline]
    pub fn z_center(&self, i_z: usize) -> f64 {
        self.spec.z_min + (i_z as f64 + 0.5) * self.spec.dz
    }

    #[inline]
    pub fn theta_center(i_theta: u32, n_theta: u32) -> f64 {
        (f64::from(i_theta) + 0.5) * std::f64::consts::TAU / f64::from(n_theta)
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

    // --- Deterministic access, traversal & reduction (§3.7) ---------------

    /// Brick containing cell (i_r, i_z), by binary search on the Morton key.
    pub fn brick_index(&self, i_r: usize, i_z: usize) -> Option<usize> {
        self.brick_index_by_coords((i_r / BRICK) as u32, (i_z / BRICK) as u32)
    }

    /// Brick at brick-coordinates (br, bz), if active.
    pub fn brick_index_by_coords(&self, br: u32, bz: u32) -> Option<usize> {
        let key = morton2(br, bz);
        self.bricks.binary_search_by_key(&key, |b| b.morton).ok()
    }

    #[inline]
    pub fn local_rz(i_r: usize, i_z: usize) -> usize {
        (i_r % BRICK) * BRICK + (i_z % BRICK)
    }

    /// Is cell (i_r, i_z) active — in bounds, in an allocated brick, mask
    /// bit set? (Out-of-range indices are simply inactive.)
    pub fn is_active(&self, i_r: usize, i_z: usize) -> bool {
        i_r < self.spec.n_r
            && i_z < self.spec.n_z
            && self
                .brick_index(i_r, i_z)
                .is_some_and(|bi| self.bricks[bi].mask & (1u64 << Self::local_rz(i_r, i_z)) != 0)
    }

    /// Is cell (i_r, i_z) solid? (Out-of-range indices are exterior.)
    pub fn is_solid(&self, i_r: usize, i_z: usize) -> bool {
        i_r < self.spec.n_r
            && i_z < self.spec.n_z
            && self.brick_index(i_r, i_z).is_some_and(|bi| {
                self.bricks[bi].solid_mask & (1u64 << Self::local_rz(i_r, i_z)) != 0
            })
    }

    /// The cell's [`Region`] (out-of-range ⇒ `Exterior`).
    pub fn region(&self, i_r: usize, i_z: usize) -> Region {
        if self.is_active(i_r, i_z) {
            Region::Gas
        } else if self.is_solid(i_r, i_z) {
            Region::Solid
        } else {
            Region::Exterior
        }
    }

    /// Config-time enumeration of every gas↔solid interface face (SOLV-1
    /// §3.5 wall faces), in the fixed deterministic order: gas cell by
    /// (i_r, i_z) lexicographic, then face direction (r−, r+, z−, z+).
    pub fn gas_solid_faces(&self) -> Vec<InterfaceFace> {
        let mut out = Vec::new();
        for i_r in 0..self.spec.n_r {
            for i_z in 0..self.spec.n_z {
                if !self.is_active(i_r, i_z) {
                    continue;
                }
                let mut probe = |nbr: Option<(usize, usize)>, dir: FaceDir| {
                    if let Some((nr, nz)) = nbr
                        && self.is_solid(nr, nz)
                    {
                        out.push(InterfaceFace {
                            gas: (i_r, i_z),
                            solid: (nr, nz),
                            dir,
                        });
                    }
                };
                probe(i_r.checked_sub(1).map(|r| (r, i_z)), FaceDir::RMinus);
                probe(Some((i_r + 1, i_z)), FaceDir::RPlus);
                probe(i_z.checked_sub(1).map(|z| (i_r, z)), FaceDir::ZMinus);
                probe(Some((i_r, i_z + 1)), FaceDir::ZPlus);
            }
        }
        out
    }

    /// Interface face area per θ-plane slice (the θ-uniform region shape
    /// makes one area serve all `n_theta` planes): the r-face ring area or
    /// the z-face annular sector, evaluated on the GAS side's ring index.
    /// Both sides of the exchange integrate the identical area because the
    /// shared face is one face: for r-faces the gas ring's outer radius IS
    /// the solid ring's inner radius (`face_radius` single owner, session
    /// 7), and z-face partners share a ring. When FND-3 partial apertures
    /// make per-side areas genuinely differ, this function is the one
    /// place the shared-face definition changes.
    pub fn interface_area_per_theta(&self, face: &InterfaceFace, n_theta: u32) -> f64 {
        match face.dir {
            // Gas at (i_r, i_z), solid at (i_r∓1, i_z): the shared face is
            // the gas cell's r∓ face.
            FaceDir::RMinus => self.face_area_r(face.gas.0, false, n_theta),
            FaceDir::RPlus => self.face_area_r(face.gas.0, true, n_theta),
            FaceDir::ZMinus | FaceDir::ZPlus => self.face_area_z(face.gas.0, n_theta),
        }
    }

    /// THE canonical active-cell traversal (§3.7): Morton brick order, then
    /// θ-planes ascending, then local cells ascending. Every sweep and
    /// reduction that isn't hand-optimized should go through here so the
    /// mask test and index decode exist in exactly one place.
    pub fn for_each_active_cell(&self, mut f: impl FnMut(CellRef)) {
        for (bi, b) in self.bricks.iter().enumerate() {
            for i_theta in 0..b.n_theta {
                let theta = Self::theta_center(i_theta, b.n_theta);
                for local in 0..BRICK_CELLS {
                    if b.mask & (1u64 << local) == 0 {
                        continue;
                    }
                    let (i_r, i_z) = b.global_rz(local);
                    f(CellRef {
                        bi,
                        i_r,
                        i_z,
                        i_theta,
                        local,
                        idx: b.cell_index(i_theta, local),
                        r: self.r_center(i_r),
                        theta,
                        z: self.z_center(i_z),
                    });
                }
            }
        }
    }

    /// Set one field from a function of the cell centroid (canonical order).
    pub fn fill_field(&mut self, f: FieldId, func: impl Fn(f64, f64, f64) -> f64) {
        for b in &mut self.bricks {
            let (br, bz, mask, nt) = (b.br, b.bz, b.mask, b.n_theta);
            let data = b.field_mut(f);
            for i_theta in 0..nt {
                let theta = Self::theta_center(i_theta, nt);
                for local in 0..BRICK_CELLS {
                    if mask & (1u64 << local) == 0 {
                        continue;
                    }
                    let i_r = br as usize * BRICK + local / BRICK;
                    let i_z = bz as usize * BRICK + local % BRICK;
                    let r = self.spec.r_min + (i_r as f64 + 0.5) * self.spec.dr;
                    let z = self.spec.z_min + (i_z as f64 + 0.5) * self.spec.dz;
                    data[i_theta as usize * BRICK_CELLS + local] = func(r, theta, z);
                }
            }
        }
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

    /// Gas-volume-weighted global sum `Σ κV·q` — the stored-total reduction
    /// of the COUP-2 conservation ledger on cut worlds. Same fixed-shape
    /// structure as [`Grid::reduce_volume_weighted`]; on full-box worlds
    /// κ = 1.0 exactly, so the two reductions agree bitwise.
    pub fn reduce_kappa_volume_weighted(&self, f: FieldId) -> f64 {
        let partials: Vec<f64> = self
            .bricks
            .iter()
            .map(|b| {
                let data = b.field(f);
                let mut acc = 0.0f64;
                for i_theta in 0..b.n_theta {
                    for local in 0..BRICK_CELLS {
                        if b.mask & (1u64 << local) != 0 {
                            let i_r = b.br as usize * BRICK + local / BRICK;
                            let v = b.kappa_rz(local) * self.cell_volume(i_r, b.n_theta);
                            acc += v * data[b.cell_index(i_theta, local)];
                        }
                    }
                }
                acc
            })
            .collect();
        tree_combine(&partials)
    }

    /// Volume-weighted global sum `Σ V·q` over the SOLID region — the
    /// solid-side stored total of the COUP-2 ledger (solid cells are uncut
    /// in the current geometry class; the full cell volume is theirs).
    pub fn reduce_solid_volume_weighted(&self, f: FieldId) -> f64 {
        let partials: Vec<f64> = self
            .bricks
            .iter()
            .map(|b| {
                let data = b.field(f);
                let mut acc = 0.0f64;
                for i_theta in 0..b.n_theta {
                    for local in 0..BRICK_CELLS {
                        if b.solid_mask & (1u64 << local) != 0 {
                            let i_r = b.br as usize * BRICK + local / BRICK;
                            let v = self.cell_volume(i_r, b.n_theta);
                            acc += v * data[b.cell_index(i_theta, local)];
                        }
                    }
                }
                acc
            })
            .collect();
        tree_combine(&partials)
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
