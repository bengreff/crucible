//! FND-3 §3.3/§3.4/§3.5 — Voxelize(): cylindrical ring cells about the
//! grid z-axis → per-cell gas fraction κ, six face apertures, PLIC
//! interface planes, and the per-(r,z) azimuthal geometry floor.
//!
//! As-built sampling scheme (the §3.3/§3.5 record):
//! - **Volume**: `N_FRAC_SAMPLES` jittered-stratified strata per cut cell
//!   in the CYLINDRICAL volume measure — uniform in (r², θ, z), i.e.
//!   `r = √(r_i² + u·(r_o² − r_i²))` — evaluated in Cartesian, where the
//!   SDF / winding number is authored (the mapping is exact).
//! - **Apertures**: `N_APERTURE_SAMPLES` (8×8) per face in each face's OWN
//!   measure: r-faces uniform (θ, z) at r_f; z-faces uniform (r², θ);
//!   θ-faces uniform (r, z) — the planar rectangle.
//! - **Jitter**: keyed SplitMix64 — a pure function of (seed, entity, tag,
//!   stratum), no RNG state (§3.5). The entity key is the (r,z)-plane
//!   linear index of the cell (volume) or of the canonical face
//!   (apertures), deliberately EXCLUDING the θ index: sample patterns in
//!   different θ-sectors are exact PARAMETER rotations of one another, so
//!   an axisymmetric solid yields per-θ-identical fractions/apertures up
//!   to the classifier's own rotational round-off (the field is evaluated
//!   at rotated Cartesian points; only a sample within an ulp of the
//!   surface could classify differently — S9 review note: the §3.4
//!   zero-variance behavior is that strong, not fully structural). Every
//!   shared face IS computed once (canonical face arrays), so the two
//!   adjacent cells carry identical bits by construction (the FND-2 §3.6
//!   shared-face validation contract). The seed stays live through every
//!   key.
//! - **Pure cells** are decided by bound tests, not sampled (§3.1):
//!   [`Solid::pure_cell_hint`] on the sector's conservative Cartesian
//!   AABB (circumradius = half its diagonal).
//!
//! Cell order (documented, fixed): `i_r`-major, then `i_theta`, then
//! `i_z` — `index = (i_r·n_theta + i_theta)·n_z + i_z`. Face order per
//! cell: `[r−, r+, z−, z+, θ−, θ+]` — the grid's **`FaceDir::index` order**,
//! the single face-order owner (S10 unification: the S9 seam wart, where the
//! voxelizer emitted §3.3(1)'s {r,θ,z} order and the ingest permuted, is
//! retired — the `CellCut → CellGeomTheta` copy is now identity).
//!
//! κ is the GAS fraction (1 − solid fraction), matching the grid's
//! `CellGeom` convention.

use super::plic::{self, PlicPlane};
use super::{Geom3dError, Solid};
use std::f64::consts::{FRAC_PI_2, TAU};

/// FND-3 §3.1: fixed jittered-stratified pattern, 8³ = 512 strata per cut
/// cell (named constant; the declared per-cell fraction-error bound
/// ε_α = C_jitter·N^(−2/3) is derived from it).
pub const N_FRAC_SAMPLES: usize = 512;

/// Per-face aperture pattern, 8×8 = 64 strata in the face's own measure
/// (the §3.3 as-built record; faces are 2-D, so the same jittered-
/// stratified argument gives an N^(−3/4) rate there — apertures are
/// gate-checked directly rather than carrying their own declared bound at
/// S9).
pub const N_APERTURE_SAMPLES: usize = 64;
const N_AP_AXIS: usize = 8;

/// The §3.1 jitter constant of `ε_α = C_JITTER · N_FRAC_SAMPLES^(−2/3)`,
/// CALIBRATED against the analytic axis-centered sphere of the §6.1
/// battery (`sphere_fractions_converge_at_the_jittered_rate`): measured
/// per-cell MAX |κ − exact| at N = 512 was 1.418e-2 ⇒ C = 0.907; declared
/// with the ×1.5 margin, rounded up: 1.4 (the measurement is reprinted by
/// the gate on every run).
pub const C_JITTER: f64 = 1.4;

/// PLIC gradient probe step, as a fraction of min(Δr, Δz): the classifier
/// field is 1-Lipschitz-scaled O(cell), so h = 1e-3·min(Δr, Δz) keeps the
/// central-difference truncation O(h²) small while staying far above the
/// f64 cancellation floor (ε/h ~ 1e-13·cell); NOT smaller — the field is
/// only C⁰ at CSG joins and fp cancellation dominates below ~1e-8·cell.
pub const H_GRAD_STEP_FRAC: f64 = 1e-3;

// Sampling-domain tags of the jitter key (fixed, documented — §3.5).
const TAG_VOLUME: u64 = 1;
const TAG_FACE_R: u64 = 2;
const TAG_FACE_THETA: u64 = 3;
const TAG_FACE_Z: u64 = 4;

/// Input spec of the cylindrical voxelization target — a plain struct,
/// deliberately NOT the grid's types (this kernel is pre-grid; FND-2's
/// θ-ladder legality is the grid's own gate, not repeated here).
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelSpec {
    pub r_min: f64,
    pub dr: f64,
    pub n_r: usize,
    pub z_min: f64,
    pub dz: f64,
    pub n_z: usize,
    pub n_theta: u32,
    pub seed: u64,
}

impl VoxelSpec {
    /// Refuse non-finite or degenerate worlds (META-1 P6).
    pub fn validate(&self) -> Result<(), Geom3dError> {
        let ok = self.r_min.is_finite()
            && self.dr.is_finite()
            && self.z_min.is_finite()
            && self.dz.is_finite()
            && self.r_min >= 0.0
            && self.dr > 0.0
            && self.dz > 0.0
            && self.n_r >= 1
            && self.n_z >= 1
            && self.n_theta >= 1;
        if ok {
            Ok(())
        } else {
            Err(Geom3dError::BadSpec(format!("{self:?}")))
        }
    }
}

/// One cell's cut geometry: gas fraction κ, the six face apertures in
/// **`FaceDir::index` order** `[r−, r+, z−, z+, θ−, θ+]` (the single
/// face-order owner, S10), and the PLIC plane for sampled-cut cells
/// (0 < κ < 1).
#[derive(Debug, Clone, PartialEq)]
pub struct CellCut {
    pub kappa: f64,
    pub aperture: [f64; 6],
    pub plic: Option<PlicPlane>,
}

/// The voxelization product: cells in the documented fixed order
/// (`i_r`-major, then `i_theta`, then `i_z`) plus the declared per-cell
/// fraction-error bound ε_α of §3.1, destined for the run manifest.
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelWorld {
    pub spec: VoxelSpec,
    pub cells: Vec<CellCut>,
    pub eps_alpha: f64,
}

impl VoxelWorld {
    /// The documented cell order, in one place.
    #[inline]
    pub fn cell_index(&self, i_r: usize, i_theta: usize, i_z: usize) -> usize {
        (i_r * self.spec.n_theta as usize + i_theta) * self.spec.n_z + i_z
    }
}

/// SplitMix64 finalizer (Steele/Lea/Flood 2014) — the one mixing
/// primitive; hand-rolled, no dependency.
#[inline]
fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[inline]
fn to_unit(x: u64) -> f64 {
    // Top 53 bits → [0, 1); exact dyadic, deterministic.
    (x >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
}

#[inline]
fn jitter_state(seed: u64, entity: u64, tag: u64, stratum: u64) -> u64 {
    let mut s = splitmix64(seed);
    s = splitmix64(s ^ entity);
    s = splitmix64(s ^ tag);
    splitmix64(s ^ stratum)
}

/// The §3.5 as-built jitter: a PURE function of (seed, entity, tag,
/// stratum) → three jitters in [0, 1). No RNG state exists anywhere.
pub fn jitter01(seed: u64, entity: u64, tag: u64, stratum: u64) -> (f64, f64, f64) {
    let a = splitmix64(jitter_state(seed, entity, tag, stratum));
    let b = splitmix64(a);
    let c = splitmix64(b);
    (to_unit(a), to_unit(b), to_unit(c))
}

/// 2-D variant for face sampling — same key chain, two outputs.
pub fn jitter01_2d(seed: u64, entity: u64, tag: u64, stratum: u64) -> (f64, f64) {
    let a = splitmix64(jitter_state(seed, entity, tag, stratum));
    let b = splitmix64(a);
    (to_unit(a), to_unit(b))
}

/// Classify one sample through the ONE classifier scalar, refusing NaN
/// (META-1 P6: never classify a non-finite field value).
#[inline]
fn sample_is_solid(solid: &Solid, p: [f64; 3]) -> Result<bool, Geom3dError> {
    let f = solid.field(p);
    if f.is_finite() {
        Ok(f < 0.0)
    } else {
        Err(Geom3dError::NonFinite {
            context: format!("classifier field at {p:?}"),
        })
    }
}

/// Conservative Cartesian AABB of the ring sector
/// `[r0, r1] × [θ0, θ1] × [z0, z1]` (r0 = r1 or z0 = z1 give face
/// patches): extrema of r·cos θ / r·sin θ occur at the θ endpoints and at
/// the cardinal angles inside the arc, at one of the two radii.
fn sector_aabb(r0: f64, r1: f64, th0: f64, th1: f64, z0: f64, z1: f64) -> ([f64; 3], [f64; 3]) {
    let mut xs = (f64::INFINITY, f64::NEG_INFINITY);
    let mut ys = (f64::INFINITY, f64::NEG_INFINITY);
    let mut consider = |th: f64| {
        let (c, s) = (th.cos(), th.sin());
        for r in [r0, r1] {
            let (x, y) = (r * c, r * s);
            xs = (xs.0.min(x), xs.1.max(x));
            ys = (ys.0.min(y), ys.1.max(y));
        }
    };
    consider(th0);
    consider(th1);
    let k_lo = (th0 / FRAC_PI_2).ceil() as i64;
    let k_hi = (th1 / FRAC_PI_2).floor() as i64;
    for k in k_lo..=k_hi {
        consider(k as f64 * FRAC_PI_2);
    }
    ([xs.0, ys.0, z0], [xs.1, ys.1, z1])
}

/// Pure-region short-circuit over an AABB: `Some(true)` = all solid,
/// `Some(false)` = all gas, `None` = sample.
fn aabb_hint(solid: &Solid, lo: [f64; 3], hi: [f64; 3]) -> Option<bool> {
    let center = [
        0.5 * (lo[0] + hi[0]),
        0.5 * (lo[1] + hi[1]),
        0.5 * (lo[2] + hi[2]),
    ];
    let dx = hi[0] - lo[0];
    let dy = hi[1] - lo[1];
    let dz = hi[2] - lo[2];
    let circumradius = 0.5 * (dx * dx + dy * dy + dz * dz).sqrt();
    solid.pure_cell_hint(center, circumradius)
}

/// Geometry of one (r,z) cell of the spec, in one place.
struct RingCell {
    r0: f64,
    r1: f64,
    z0: f64,
    z1: f64,
}

fn ring_cell(spec: &VoxelSpec, i_r: usize, i_z: usize) -> RingCell {
    RingCell {
        r0: spec.r_min + i_r as f64 * spec.dr,
        r1: spec.r_min + (i_r + 1) as f64 * spec.dr,
        z0: spec.z_min + i_z as f64 * spec.dz,
        z1: spec.z_min + (i_z + 1) as f64 * spec.dz,
    }
}

/// Volume fraction of one cell sector: pure hint, else jittered-stratified
/// sampling in the cylindrical volume measure. Returns the GAS fraction κ.
fn cell_kappa(
    solid: &Solid,
    spec: &VoxelSpec,
    na: usize,
    i_r: usize,
    i_t: usize,
    i_z: usize,
) -> Result<f64, Geom3dError> {
    let c = ring_cell(spec, i_r, i_z);
    let dth = TAU / f64::from(spec.n_theta);
    let th0 = i_t as f64 * dth;
    let (lo, hi) = sector_aabb(c.r0, c.r1, th0, th0 + dth, c.z0, c.z1);
    if let Some(all_solid) = aabb_hint(solid, lo, hi) {
        return Ok(if all_solid { 0.0 } else { 1.0 });
    }
    // θ-independent entity key (see the module doc: exact θ-congruence).
    let entity = (i_r * spec.n_z + i_z) as u64;
    let r0sq = c.r0 * c.r0;
    let r1sq = c.r1 * c.r1;
    let naf = na as f64;
    let mut solid_cnt = 0usize;
    for a in 0..na {
        for b in 0..na {
            for cc in 0..na {
                let stratum = ((a * na + b) * na + cc) as u64;
                let (ju, jt, jz) = jitter01(spec.seed, entity, TAG_VOLUME, stratum);
                let u = (a as f64 + ju) / naf;
                let r = (r0sq + u * (r1sq - r0sq)).sqrt();
                let th = th0 + (b as f64 + jt) / naf * dth;
                let z = c.z0 + (cc as f64 + jz) / naf * spec.dz;
                if sample_is_solid(solid, [r * th.cos(), r * th.sin(), z])? {
                    solid_cnt += 1;
                }
            }
        }
    }
    Ok(1.0 - solid_cnt as f64 / (naf * naf * naf))
}

/// r-face aperture at radius `r_f` over sector `i_t`, span `[z0, z1]` —
/// uniform in the face's own (θ, z) measure. At r_f = 0 (the axis) the
/// face has zero area; the sampled value (every sample is the axis point)
/// is still emitted deterministically — consumers weight by the zero area.
fn r_face_aperture(
    solid: &Solid,
    spec: &VoxelSpec,
    i_fr: usize,
    i_t: usize,
    i_z: usize,
) -> Result<f64, Geom3dError> {
    let r_f = spec.r_min + i_fr as f64 * spec.dr;
    let c = ring_cell(spec, i_fr.min(spec.n_r - 1), i_z);
    let dth = TAU / f64::from(spec.n_theta);
    let th0 = i_t as f64 * dth;
    let (lo, hi) = sector_aabb(r_f, r_f, th0, th0 + dth, c.z0, c.z1);
    if let Some(all_solid) = aabb_hint(solid, lo, hi) {
        return Ok(if all_solid { 0.0 } else { 1.0 });
    }
    let entity = (i_fr * spec.n_z + i_z) as u64;
    let naf = N_AP_AXIS as f64;
    let mut solid_cnt = 0usize;
    for a in 0..N_AP_AXIS {
        for b in 0..N_AP_AXIS {
            let stratum = (a * N_AP_AXIS + b) as u64;
            let (j1, j2) = jitter01_2d(spec.seed, entity, TAG_FACE_R, stratum);
            let th = th0 + (a as f64 + j1) / naf * dth;
            let z = c.z0 + (b as f64 + j2) / naf * spec.dz;
            if sample_is_solid(solid, [r_f * th.cos(), r_f * th.sin(), z])? {
                solid_cnt += 1;
            }
        }
    }
    Ok(1.0 - solid_cnt as f64 / (naf * naf))
}

/// θ-face aperture at angle `θ_f = i_ft·Δθ` over cell (i_r, i_z) — the
/// planar rectangle `[r0, r1] × [z0, z1]`, uniform in (r, z) (the planar
/// measure — a θ-face is flat). The jitter entity excludes the face angle,
/// so every θ-face of a ring carries the same (r, z) pattern (exact
/// θ-congruence); at n_theta = 1 the single seam face at θ = 0 serves as
/// both θ− and θ+ of the full ring.
fn theta_face_aperture(
    solid: &Solid,
    spec: &VoxelSpec,
    i_r: usize,
    i_ft: usize,
    i_z: usize,
) -> Result<f64, Geom3dError> {
    let c = ring_cell(spec, i_r, i_z);
    let dth = TAU / f64::from(spec.n_theta);
    let th_f = i_ft as f64 * dth;
    let (ct, st) = (th_f.cos(), th_f.sin());
    let lo = [(c.r0 * ct).min(c.r1 * ct), (c.r0 * st).min(c.r1 * st), c.z0];
    let hi = [(c.r0 * ct).max(c.r1 * ct), (c.r0 * st).max(c.r1 * st), c.z1];
    if let Some(all_solid) = aabb_hint(solid, lo, hi) {
        return Ok(if all_solid { 0.0 } else { 1.0 });
    }
    let entity = (i_r * spec.n_z + i_z) as u64;
    let naf = N_AP_AXIS as f64;
    let mut solid_cnt = 0usize;
    for a in 0..N_AP_AXIS {
        for b in 0..N_AP_AXIS {
            let stratum = (a * N_AP_AXIS + b) as u64;
            let (j1, j2) = jitter01_2d(spec.seed, entity, TAG_FACE_THETA, stratum);
            let r = c.r0 + (a as f64 + j1) / naf * spec.dr;
            let z = c.z0 + (b as f64 + j2) / naf * spec.dz;
            if sample_is_solid(solid, [r * ct, r * st, z])? {
                solid_cnt += 1;
            }
        }
    }
    Ok(1.0 - solid_cnt as f64 / (naf * naf))
}

/// z-face aperture at `z_f` over the annular sector `[r0, r1] × sector` —
/// uniform in the face's own (r², θ) annular measure.
fn z_face_aperture(
    solid: &Solid,
    spec: &VoxelSpec,
    i_r: usize,
    i_t: usize,
    i_fz: usize,
) -> Result<f64, Geom3dError> {
    let c = ring_cell(spec, i_r, i_fz.min(spec.n_z - 1));
    let z_f = spec.z_min + i_fz as f64 * spec.dz;
    let dth = TAU / f64::from(spec.n_theta);
    let th0 = i_t as f64 * dth;
    let (lo, hi) = sector_aabb(c.r0, c.r1, th0, th0 + dth, z_f, z_f);
    if let Some(all_solid) = aabb_hint(solid, lo, hi) {
        return Ok(if all_solid { 0.0 } else { 1.0 });
    }
    let entity = (i_r * (spec.n_z + 1) + i_fz) as u64;
    let r0sq = c.r0 * c.r0;
    let r1sq = c.r1 * c.r1;
    let naf = N_AP_AXIS as f64;
    let mut solid_cnt = 0usize;
    for a in 0..N_AP_AXIS {
        for b in 0..N_AP_AXIS {
            let stratum = (a * N_AP_AXIS + b) as u64;
            let (j1, j2) = jitter01_2d(spec.seed, entity, TAG_FACE_Z, stratum);
            let u = (a as f64 + j1) / naf;
            let r = (r0sq + u * (r1sq - r0sq)).sqrt();
            let th = th0 + (b as f64 + j2) / naf * dth;
            if sample_is_solid(solid, [r * th.cos(), r * th.sin(), z_f])? {
                solid_cnt += 1;
            }
        }
    }
    Ok(1.0 - solid_cnt as f64 / (naf * naf))
}

/// PLIC for one sampled-cut cell (FND-3 §3.3: planes are fitted in the
/// cell's LOCAL Cartesian frame). Frame origin = the cell's metric
/// midpoint (r̄ = (r_i + r_o)/2, θ mid, z mid; the r²-measure volume
/// centroid differs by O(Δr²/r̄) — within the box model's own first-order
/// honesty); axes = (r̂, θ̂, ẑ) there; local box dims (Δr, r̄·Δθ, Δz), whose
/// product equals the ring-sector volume exactly. At n_theta = 1 the θ
/// dimension is the FULL ring (Δθ = 2π): a local Cartesian box is a poor
/// model of the annulus, but the same formula is emitted and DOCUMENTED as
/// the declared local-frame approximation rather than special-cased.
///
/// The normal operand splits by authoring path: CSG uses the central-
/// difference gradient of the SDF at the midpoint (step [`H_GRAD_STEP_FRAC`]);
/// a MESH cannot — the winding number of a watertight mesh is piecewise
/// constant, so the field's central difference is analytically zero away
/// from the surface — and takes the area-weighted outward facet normal
/// over the cell's AABB instead ([`super::TriMesh::interface_normal_in_box`]),
/// exact for flat facets and first-order in curvature, the same honesty
/// class as PLIC itself.
fn fit_cell_plic(
    solid: &Solid,
    spec: &VoxelSpec,
    i_r: usize,
    i_t: usize,
    i_z: usize,
    kappa: f64,
    cell: usize,
) -> Result<PlicPlane, Geom3dError> {
    let c = ring_cell(spec, i_r, i_z);
    let dth = TAU / f64::from(spec.n_theta);
    let rbar = 0.5 * (c.r0 + c.r1);
    let thc = (i_t as f64 + 0.5) * dth;
    let zc = 0.5 * (c.z0 + c.z1);
    let (ct, st) = (thc.cos(), thc.sin());
    let p0 = [rbar * ct, rbar * st, zc];
    let er = [ct, st, 0.0];
    let et = [-st, ct, 0.0];
    let ez = [0.0, 0.0, 1.0];
    let h = H_GRAD_STEP_FRAC * spec.dr.min(spec.dz);

    let g = match solid {
        Solid::Csg(_) => {
            let probe = |dir: [f64; 3], sgn: f64| -> Result<f64, Geom3dError> {
                let p = [
                    p0[0] + sgn * h * dir[0],
                    p0[1] + sgn * h * dir[1],
                    p0[2] + sgn * h * dir[2],
                ];
                let v = solid.field(p);
                if v.is_finite() {
                    Ok(v)
                } else {
                    Err(Geom3dError::NonFinite {
                        context: format!("PLIC gradient probe at {p:?}"),
                    })
                }
            };
            let mut g = [0.0f64; 3];
            for (k, dir) in [er, et, ez].into_iter().enumerate() {
                g[k] = (probe(dir, 1.0)? - probe(dir, -1.0)?) / (2.0 * h);
            }
            g
        }
        Solid::Mesh(mesh) => {
            let th0 = i_t as f64 * dth;
            let (lo, hi) = sector_aabb(c.r0, c.r1, th0, th0 + dth, c.z0, c.z1);
            let ng = mesh
                .interface_normal_in_box(lo, hi)
                .ok_or(Geom3dError::DegeneratePlicNormal { cell })?;
            // Global Cartesian → local (r̂, θ̂, ẑ) components.
            [
                ng[0] * er[0] + ng[1] * er[1],
                ng[0] * et[0] + ng[1] * et[1],
                ng[2],
            ]
        }
    };

    let l = [spec.dr, rbar * dth, spec.dz];
    let (n, d, area) = plic::fit_plane(g, l, 1.0 - kappa, cell)?;
    let n_global = [
        n[0] * er[0] + n[1] * et[0],
        n[0] * er[1] + n[1] * et[1],
        n[2],
    ];
    Ok(PlicPlane {
        normal: n,
        d,
        centroid: [
            p0[0] + d * n_global[0],
            p0[1] + d * n_global[1],
            p0[2] + d * n_global[2],
        ],
        area,
    })
}

/// The FND-3 §3.3 voxelizer at the declared `N_FRAC_SAMPLES` pattern.
pub fn voxelize(solid: &Solid, spec: &VoxelSpec) -> Result<VoxelWorld, Geom3dError> {
    voxelize_with_samples(solid, spec, N_FRAC_SAMPLES)
}

/// Sample-count-parametrized variant — the §6.1 convergence battery's
/// knob (the C_JITTER calibration sweeps N ∈ {64, 512, 4096}); production
/// callers use [`voxelize`]. `n_frac` must be a perfect cube.
#[doc(hidden)]
pub fn voxelize_with_samples(
    solid: &Solid,
    spec: &VoxelSpec,
    n_frac: usize,
) -> Result<VoxelWorld, Geom3dError> {
    spec.validate()?;
    solid.validate()?;
    let na = (n_frac as f64).cbrt().round() as usize;
    if na < 1 || na * na * na != n_frac {
        return Err(Geom3dError::BadSpec(format!(
            "N_FRAC_SAMPLES = {n_frac} is not a perfect cube"
        )));
    }
    let nt = spec.n_theta as usize;
    let (n_r, n_z) = (spec.n_r, spec.n_z);

    // Canonical face arrays — each physical face computed exactly once, so
    // shared faces are bitwise coherent by construction (FND-2 §3.6).
    let mut r_ap = vec![0.0f64; (n_r + 1) * nt * n_z];
    for i_fr in 0..=n_r {
        for i_t in 0..nt {
            for i_z in 0..n_z {
                r_ap[(i_fr * nt + i_t) * n_z + i_z] = r_face_aperture(solid, spec, i_fr, i_t, i_z)?;
            }
        }
    }
    let mut th_ap = vec![0.0f64; n_r * nt * n_z];
    for i_r in 0..n_r {
        for i_ft in 0..nt {
            for i_z in 0..n_z {
                th_ap[(i_r * nt + i_ft) * n_z + i_z] =
                    theta_face_aperture(solid, spec, i_r, i_ft, i_z)?;
            }
        }
    }
    let mut z_ap = vec![0.0f64; n_r * nt * (n_z + 1)];
    for i_r in 0..n_r {
        for i_t in 0..nt {
            for i_fz in 0..=n_z {
                z_ap[(i_r * nt + i_t) * (n_z + 1) + i_fz] =
                    z_face_aperture(solid, spec, i_r, i_t, i_fz)?;
            }
        }
    }

    let mut cells = Vec::with_capacity(n_r * nt * n_z);
    for i_r in 0..n_r {
        for i_t in 0..nt {
            for i_z in 0..n_z {
                let kappa = cell_kappa(solid, spec, na, i_r, i_t, i_z)?;
                let idx = cells.len();
                let plic = if kappa > 0.0 && kappa < 1.0 {
                    Some(fit_cell_plic(solid, spec, i_r, i_t, i_z, kappa, idx)?)
                } else {
                    None
                };
                // Emitted in **FaceDir::index order** {r−, r+, z−, z+, θ−, θ+}
                // — the single face-order owner (FND-2 §3.3(1)/lib.rs). The
                // S9 seam wart (the voxelizer once emitted {r,θ,z}) is retired
                // here (S10): the grid ingest is now a plain field-by-field
                // copy, no permutation to get wrong.
                let aperture = [
                    r_ap[(i_r * nt + i_t) * n_z + i_z],
                    r_ap[((i_r + 1) * nt + i_t) * n_z + i_z],
                    z_ap[(i_r * nt + i_t) * (n_z + 1) + i_z],
                    z_ap[(i_r * nt + i_t) * (n_z + 1) + i_z + 1],
                    th_ap[(i_r * nt + i_t) * n_z + i_z],
                    th_ap[(i_r * nt + (i_t + 1) % nt) * n_z + i_z],
                ];
                cells.push(CellCut {
                    kappa,
                    aperture,
                    plic,
                });
            }
        }
    }

    Ok(VoxelWorld {
        spec: spec.clone(),
        cells,
        // The declared §3.1 bound, destined for the run manifest (FND-6).
        eps_alpha: C_JITTER * (n_frac as f64).powf(-2.0 / 3.0),
    })
}

/// The coarsest azimuthal resolution `N_θ' ∈ {n_theta, n_theta/2, …}` (the
/// power-of-two ladder down from the sampled `n_theta`) whose piecewise-
/// constant θ-coarsening **reproduces** the per-sector `values` within
/// `eps`: for the tested `N_θ'`, every fine sector's value differs from its
/// coarse group's mean by ≤ `eps`. Searches coarse→fine and returns the
/// FIRST (coarsest) level that passes — 1 when the quantity is θ-uniform,
/// `nt` when even the second-coarsest level cannot represent the variation.
/// Only halving factors that divide `nt` are considered (a non-power-of-two
/// sampled `n_theta`, never produced by the FND-4 θ-ladder, stops at its odd
/// base — conservative, never coarser than it can halve).
fn coarsest_reproducing(values: &[f64], nt: usize, eps: f64) -> u32 {
    // Candidate levels, coarsest first: 1, 2, 4, … up to nt, keeping only
    // those that divide nt evenly (so each group is a contiguous equal split).
    let mut level = 1usize;
    let mut candidates = Vec::new();
    while level <= nt {
        if nt.is_multiple_of(level) {
            candidates.push(level);
        }
        level *= 2;
    }
    for &nth in &candidates {
        let group = nt / nth; // fine sectors per coarse sector
        let mut ok = true;
        'grp: for g in 0..nth {
            let base = g * group;
            let mut mean = 0.0f64;
            for k in 0..group {
                mean += values[base + k];
            }
            mean /= group as f64;
            for k in 0..group {
                if (values[base + k] - mean).abs() > eps {
                    ok = false;
                    break 'grp;
                }
            }
        }
        if ok {
            return nth as u32;
        }
    }
    nt as u32
}

/// FND-3 §3.4, **S10 coarsest-reproducing form**: the per-(r,z) azimuthal
/// geometry floor. Returns one `u32` per (i_r, i_z) — indexed `i_r·n_z + i_z`
/// — equal to the **coarsest N_θ that reproduces the cell's κ and all six
/// face apertures within the declared ε_α** (the max over the seven
/// quantities of [`coarsest_reproducing`]): 1 for a θ-uniform (axisymmetric)
/// cell, and only as fine as the geometry's own azimuthal structure demands
/// — never the blanket `n_theta` the S9 conservative form pinned on ANY
/// variation. The θ-coarsening projection this needs is the plan-S10
/// refinement machinery (FND-2 §3.5); the result holds independent of the
/// flow state and is a FLOOR (the FND-2 consumer takes `max(N_θ^geom,
/// N_θ^guard)` for adaptive decisions). Reducing per-region floors to
/// per-BRICK floors (max over the brick's cells) remains the FND-2
/// consumer's job.
pub fn theta_geom_floor(world: &VoxelWorld) -> Vec<u32> {
    let s = &world.spec;
    let nt = s.n_theta as usize;
    let mut out = vec![1u32; s.n_r * s.n_z];
    let mut scratch = vec![0.0f64; nt];
    for i_r in 0..s.n_r {
        for i_z in 0..s.n_z {
            let mut floor = 1u32;
            for q in 0..7usize {
                for (i_t, v) in scratch.iter_mut().enumerate() {
                    let cell = &world.cells[(i_r * nt + i_t) * s.n_z + i_z];
                    *v = if q == 0 {
                        cell.kappa
                    } else {
                        cell.aperture[q - 1]
                    };
                }
                floor = floor.max(coarsest_reproducing(&scratch, nt, world.eps_alpha));
            }
            out[i_r * s.n_z + i_z] = floor;
        }
    }
    out
}

#[cfg(test)]
mod floor_tests {
    use super::coarsest_reproducing;

    // The coarsest-reproducing search (S10) over the power-of-two θ-ladder.
    #[test]
    fn coarsest_reproducing_finds_the_ladder_rung() {
        let eps = 1e-3;
        // θ-uniform ⇒ no floor.
        assert_eq!(coarsest_reproducing(&[0.5; 8], 8, eps), 1);
        // A pure m=1 half/half split reproduces at N_θ = 2 (each half is
        // internally uniform), never coarser (the two halves differ).
        assert_eq!(
            coarsest_reproducing(&[0.2, 0.2, 0.2, 0.2, 0.9, 0.9, 0.9, 0.9], 8, eps),
            2
        );
        // Four distinct quadrant values ⇒ N_θ = 4.
        assert_eq!(
            coarsest_reproducing(&[0.1, 0.1, 0.4, 0.4, 0.7, 0.7, 0.95, 0.95], 8, eps),
            4
        );
        // A single localized sector (the sharp-feature case) ⇒ the full 8.
        assert_eq!(
            coarsest_reproducing(&[0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0, 0.0], 8, eps),
            8
        );
        // Within-tolerance ripple below ε_α is reproduced at the coarse rung
        // (the fixed-count projection absorbs it — never a spurious floor).
        assert_eq!(
            coarsest_reproducing(&[0.5, 0.5 + 5e-4, 0.5 - 5e-4, 0.5], 4, eps),
            1
        );
    }
}
