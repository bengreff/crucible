//! `geom3d` — the FND-3 §3 S9 geometry kernel: CSG analytic SDF trees
//! (§3.1), STL import + the exact generalized winding number (§3.2), and
//! `voxelize()` onto cylindrical ring cells — jittered-stratified partial
//! fractions, six face apertures, and PLIC interface planes (§3.3) — plus
//! the per-(r,z) azimuthal geometry floor (§3.4) and §3.5 determinism
//! (a pure function of {geometry, resolution, seed}; fixed tie-breaks;
//! keyed SplitMix64 jitter — no RNG state anywhere).
//!
//! Self-contained by design: nothing here touches `CellGeom`/`BrickGeom`/
//! `Grid::build_with_geometry` or any other existing (certified) path.
//! Wiring the emitted [`VoxelWorld`] through the FND-2 §3.6 ingest seam is
//! the consumer's job in a later session.
//!
//! **The one sign convention (stated once):** [`Solid::field`] is negative
//! inside the solid, positive in the gas, and its gradient points OUT OF
//! the solid — INTO the gas. Every classification and every PLIC normal in
//! this module derives from that single scalar; the tie `field == 0.0`
//! counts as NOT inside (fixed tie-breaking, §3.5).

mod plic;
mod sdf;
mod stl;
mod voxel;

pub use plic::{N_PLIC_BISECT, PlicPlane};
pub use sdf::Sdf;
pub use stl::{TriMesh, tri_touches_box};
pub use voxel::{
    C_JITTER, CellCut, H_GRAD_STEP_FRAC, N_APERTURE_SAMPLES, N_FRAC_SAMPLES, VoxelSpec, VoxelWorld,
    jitter01, jitter01_2d, theta_geom_floor, voxelize, voxelize_with_samples,
};

/// Typed refusal surface of the geometry kernel (META-1 P6: fail loud,
/// halt clean, never guess — out-of-envelope inputs refuse, nothing is
/// clamped or silently repaired).
#[derive(Debug, Clone, PartialEq)]
pub enum Geom3dError {
    /// A [`VoxelSpec`] that cannot describe a cylindrical world (non-finite
    /// or non-positive extents, zero cells, negative `r_min`, a sample
    /// count that is not a perfect cube).
    BadSpec(String),
    /// A [`Solid`] whose parameters are non-finite or geometrically
    /// degenerate (non-positive radius, inverted box, empty mesh, a
    /// revolved profile reaching s < 0 or with a zero-length edge, …).
    BadSolid(String),
    /// STL bytes that satisfy neither the binary record-count size
    /// identity nor the ASCII `solid` grammar, or that carry non-finite
    /// vertex data.
    StlParse(String),
    /// A NaN/inf surfaced mid-computation (e.g. the classifier field at a
    /// sample point) — halt with diagnosis, never classify a NaN.
    NonFinite { context: String },
    /// A cut cell whose classifier-field gradient vanishes or is
    /// non-finite at the centroid: no PLIC normal is derivable — refused,
    /// not defaulted (§3.5: no silent fallback orientation).
    DegeneratePlicNormal { cell: usize },
}

impl std::fmt::Display for Geom3dError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BadSpec(m) => write!(f, "bad voxel spec: {m}"),
            Self::BadSolid(m) => write!(f, "bad solid: {m}"),
            Self::StlParse(m) => write!(f, "STL parse refused: {m}"),
            Self::NonFinite { context } => {
                write!(
                    f,
                    "non-finite value in {context} — halt with diagnosis (META-1 P6)"
                )
            }
            Self::DegeneratePlicNormal { cell } => {
                write!(
                    f,
                    "cell {cell}: degenerate PLIC normal (zero/non-finite field gradient)"
                )
            }
        }
    }
}

impl std::error::Error for Geom3dError {}

/// The classifier seam over both authoring paths (FND-3 §2): a CSG SDF
/// tree (§3.1) or a triangle mesh classified by the generalized winding
/// number (§3.2). Both produce identical downstream cell data — the
/// authoring path sets investment emphasis, not capability.
#[derive(Debug, Clone, PartialEq)]
pub enum Solid {
    Csg(Sdf),
    Mesh(TriMesh),
}

impl Solid {
    /// Refuse degenerate/non-finite geometry before any sampling begins
    /// (META-1 P6). Called by [`voxelize`]; callers constructing solids by
    /// hand may call it directly.
    pub fn validate(&self) -> Result<(), Geom3dError> {
        match self {
            Solid::Csg(sdf) => sdf.validate(),
            Solid::Mesh(mesh) => mesh.validate(),
        }
    }

    /// The one classifier scalar (see the module-level sign convention):
    /// negative inside the solid, positive in the gas, gradient pointing
    /// out of the solid — the PLIC normal direction. For CSG it is the SDF
    /// value itself; for a mesh it is `0.5 − w` (the winding number `w`
    /// decreases outward through the surface, so `0.5 − w` increases).
    pub fn field(&self, p: [f64; 3]) -> f64 {
        match self {
            Solid::Csg(sdf) => sdf.eval(p),
            Solid::Mesh(mesh) => 0.5 - mesh.winding_number(p),
        }
    }

    /// Point classification with the fixed tie-break: `field == 0.0`
    /// (SDF surface / winding exactly ½) counts as NOT inside (§3.5).
    pub fn inside(&self, p: [f64; 3]) -> bool {
        self.field(p) < 0.0
    }

    /// Provable pure-cell test (FND-3 §3.1: interior/exterior cells are
    /// decided by bound tests, not sampled). `Some(true)` = the ball of
    /// radius `circumradius` about `center` is entirely INSIDE the solid;
    /// `Some(false)` = entirely gas; `None` = cannot prove either — sample.
    ///
    /// CSG: every [`Sdf`] node is 1-Lipschitz (exact primitive distances;
    /// min/max preserve the property), so `|eval(center)| > circumradius`
    /// proves the whole ball is one-signed — sound even where booleans
    /// make the magnitude only a lower bound (§3.1 caveat). Mesh: if no
    /// triangle touches the axis-aligned box `center ± circumradius`
    /// (which contains the circumscribed ball, hence the cell), the
    /// winding number is threshold-stable across the cell **for a
    /// WATERTIGHT mesh** (w is piecewise constant away from the surface) —
    /// classify by the center. For an imperfect mesh (the graceful-
    /// degradation class §3.2 exists for) w is merely harmonic away from
    /// triangles and can cross 0.5 inside a triangle-free cell near a
    /// hole's throat, so the hint may disagree with what sampling would
    /// have measured there — an S9 review-recorded limit; a declared
    /// watertightness gate on the hint rides the production-STL wave with
    /// Barnes-Hut.
    pub fn pure_cell_hint(&self, center: [f64; 3], circumradius: f64) -> Option<bool> {
        match self {
            Solid::Csg(sdf) => {
                let v = sdf.eval(center);
                if v.is_finite() && v.abs() > circumradius {
                    Some(v < 0.0)
                } else {
                    None
                }
            }
            Solid::Mesh(mesh) => {
                let bmin = [
                    center[0] - circumradius,
                    center[1] - circumradius,
                    center[2] - circumradius,
                ];
                let bmax = [
                    center[0] + circumradius,
                    center[1] + circumradius,
                    center[2] + circumradius,
                ];
                // Fixed iteration order (mesh order) — §3.5.
                if mesh.tris.iter().any(|t| tri_touches_box(t, bmin, bmax)) {
                    None
                } else {
                    Some(mesh.winding_number(center) > 0.5)
                }
            }
        }
    }
}
