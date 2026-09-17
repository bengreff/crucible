//! SOLV-1 §3.1 — **`F_visc`, the missing forces (plan S3; the full θ-stress
//! tensor at uniform N_θ landed at plan S9 — COUP-3 0.4.5's fixed design,
//! built)**: compressible viscous stress + Fourier heat conduction +
//! Fickian species diffusion on the exact cylindrical metric. At N_θ > 1
//! (uniform per-brick N_θ only) every component gains its **θ-θ implicit
//! core** on the periodic within-brick ring stencil (μ for u_r/u_z,
//! (4/3)μ·r̄² for the ω angular-momentum form, k for T, ρD for C) and the
//! **θ curvature/cross couplings** ride the same fixed Picard lag as the
//! S3 cross-stress terms (the τ_rθ/τ_θz θ-limbs, the ∂ω/∂θ dilatation limb
//! in e_θθ, and the θ work fluxes through the same total-energy
//! bookkeeping). The COMPLETE tensor is the point: the θ-θ core alone
//! would spuriously damp m = 1 translation (u_r = U cosθ, u_θ = −U sinθ
//! has zero true stress; its curvature partners cancel the core — the
//! transverse-flow gate measures that cancellation at O(Δθ²)). **Mixed
//! per-brick N_θ refuses → plan S11** (the ring-interface coupling of an
//! implicit operator is ◆C3-wave content), and **cut geometry at N_θ > 1
//! refuses → plan S11** with it (θ-face apertures + per-θ wall patches).
//! One flux-form operator over the gas state —
//! no material or regime branch; transport (μ, k, c_v, ρD, ∂h/∂Z) is a
//! **per-cell query on the FND-7 §3.3 spine** ([`crate::transport`], S4),
//! the same one provider the wall law next door reads. This operator holds
//! no transport constant: the caller refreshes a [`GasTransportField`] from
//! the spine once per Picard iterate, from that iterate's lag state
//! (COUP-3 §3.1, 0.4.3 — never inside the CG, which must stay the solve of
//! one fixed linear operator).
//!
//! ## The operator (continuous form, ∂/∂θ = 0)
//!
//! Strain: `e_rr = ∂u_r/∂r`, `e_θθ = u_r/r`, `e_zz = ∂u_z/∂z`,
//! `e_rz = ½(∂u_r/∂z + ∂u_z/∂r)`, `e_rθ = ½ r ∂(u_θ/r)/∂r`,
//! `e_θz = ½ ∂u_θ/∂z`; Stokes hypothesis `λ = −⅔μ`:
//! `τ_ij = 2μ e_ij − ⅔μ (∇·u) δ_ij`.
//!
//! - r-momentum: `(1/r)∂_r(r τ_rr) + ∂_z τ_rz − τ_θθ/r`
//! - θ-momentum: `(1/r²)∂_r(r² τ_rθ) + ∂_z τ_θz` — the **angular-momentum
//!   form**: solved in ω = u_θ/r, where the whole operator collapses to a
//!   symmetric two-point diffusion with r³-class weights (rigid rotation
//!   ω = const is **discretely stress-free**, and angular momentum
//!   telescopes exactly; the linear θ-momentum it induces is ledgered as
//!   an applied source, exactly like the flow operator's swirl source)
//! - z-momentum: `(1/r)∂_r(r τ_rz) + ∂_z τ_zz`
//! - energy (total-energy flux form): `∇·(τ·u + k∇T + Σ_k h_k j_k)` —
//!   dissipation is not a separate term; it emerges from the
//!   KE/internal-energy bookkeeping, which is what conserves total energy
//!   by construction. The **species-enthalpy flux** `Σ_k h_k j_k` (S4)
//!   reduces exactly to `(∂h/∂Z)|_{p,T}·j_Z` under one composition
//!   coordinate; its ∇T limb is NOT here — that limb is already inside the
//!   spine's effective conductivity, and adding it again would double-count
//!   one flux (`crate::transport` module doc has the decomposition)
//! - species: `∇·(ρD ∇C)` with ρD = μ/Sc from the spine (the Fickian
//!   gradient-diffusion occupant; the LES subgrid flux of SOLV-1 §3.4 is a
//!   later occupant of this same F_visc-class seam)
//!
//! ## Class-D treatment (COUP-3 §3.1) — implicit cores + fixed Picard
//!
//! Each solved variable (u_r, ω, u_z, T, C) gets a **symmetric two-point
//! implicit core** (its normal-gradient stress/flux terms + the negative-
//! definite geometric diagonal), solved matrix-free by the same
//! deterministic fixed-structure Jacobi-CG as the solid class-D solve
//! (`EPS_CG_RESID`/`N_CG_ITERS_MAX`, δ-form warm start). The cross-stress
//! couplings (the τ_rz cross-derivatives and the −⅔μ∇·u compressible
//! corrections between u_r and u_z) are converged by the **fixed Picard
//! sweeps of the SDC orchestrator** (the COUP-2 §3.5 discipline: fixed
//! count + named residual acceptance → `COUPLING_RESIDUAL`); their fixed
//! point is the fully-implicit solution of the complete F_visc operator,
//! so the unconditional class-D stability claim holds *at acceptance*, not
//! by hope. The contraction is strong at any Δt: the lagged terms' row
//! sums are dominated by the implicit diagonals (AM-GM bounds the stiff-
//! limit gain near 1/12). ω has **no** lagged remainder — the swirl solve
//! is exactly implicit. T is solved after the velocities (its work fluxes
//! use the same-iterate accepted velocities); C is independent.
//!
//! ## Ownership at walls (SOLV-1 §3.5 / COUP-2 §3.5 — no double count)
//!
//! Resolved diffusion flows **only through gas↔gas faces** (aperture-
//! weighted, FND-3): a face against solid or exterior contributes nothing
//! here — the one wall-function law owns config-time wall faces, and its
//! heat enters as the exchange debit in the SDC step. Domain edges take
//! the **declared viscous BCs** below (data, not code): no-slip walls with
//! declared wall velocity (a moving wall does ledgered work — the Couette
//! drive), free-slip, isothermal or adiabatic; species is zero-flux at
//! every domain edge (non-catalytic — the only occupant this session).
//!
//! ## The temperature solve's slope (S4)
//!
//! `fill_mass` gives the `T` component the mass `ρ·c_v·κV`, and `c_v` is now
//! the spine's **equilibrium** `∂e/∂T|_ρ` per cell — the derivative of the
//! very surface the Picard re-derives `T` from each iterate
//! (`sdc::derive_gas_operands` reads the composed state through
//! `EosLaw::prim_checked` + the FND-7 temperature query). That is what
//! "TableEos-consistent T refresh" means concretely: the linearization
//! slope and the surface that closes the iterate are the same object.
//! `c_v` never enters the fixed point — at convergence the right-hand side
//! is zero and the mass divides nothing — so this is about *reaching* the
//! fixed point inside the fixed sweep count. It matters: across the RL10's
//! state range the equilibrium `c_v` spans more than an order of magnitude
//! (recombination stores energy the frozen value cannot see), so the
//! constant slope S3 used would leave the truncated solve far off wherever
//! dissociation runs.
//!
//! ## Recorded deferrals (owners named)
//! - **Wall-function skin-friction momentum debit**: the wall law is a
//!   heat law (SOLV-1 §3.5); its tangential-force leg rides the COUP-2
//!   §3.1.2 mount-reaction ledger (the verdict wave). Until then wall-law
//!   faces are momentum-slip, declared.
//! - θ-diffusion fluxes and per-θ operand keying (plan S9; the S8 split) —
//!   refused, not guessed.
//! - Near-wall/boundary lagged-cross stencils are one-sided (first-order
//!   locally — the MMS battery verifies the composed order).
//!
//! Determinism: serial fixed-order assembly (brick/local/face order — the
//! same recorded perf deferral as the solid assembly; the GPU wave
//! parallelizes both), `tree_combine` reductions in the CG, fixed sweep
//! structure. Bit-identical at any thread count.

use crate::euler::{Cons, FlowLedger, I_EN, I_MR, I_MT, I_MZ, I_RC, NCOMP};
use crate::sdc::{EPS_CG_RESID, N_CG_ITERS_MAX};
use crate::transport::TransportProps;
use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, Grid, tree_combine};

pub(crate) type BufF = Vec<Vec<f64>>;

#[derive(Debug, Clone, PartialEq)]
pub enum GasDiffError {
    /// Mixed per-brick N_θ in the gas class-D solve rides plan S11 (the
    /// ring-interface coupling of an implicit operator is the ◆C3 wave's
    /// content — the explicit flow reflux does not transfer to the CG
    /// stencil unmodified, COUP-3 0.4.5). Refuse.
    AzimuthalResolution,
    /// Cut geometry at N_θ > 1 in the gas class-D solve rides plan S11
    /// (θ-face apertures + per-θ wall patches — FND-2 0.5.3 makes the
    /// flow legal there first). Refuse.
    CutThetaGeometry,
    /// Non-finite operand in the assembly — halt with diagnosis.
    NonFinite {
        i_r: usize,
        i_z: usize,
        what: &'static str,
    },
    /// The fixed-structure CG missed its acceptance (COUP-3 §3.7 rule) —
    /// the caller raises `COUPLING_RESIDUAL`.
    CgUnconverged { comp: &'static str, resid: f64 },
}

impl std::fmt::Display for GasDiffError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AzimuthalResolution => write!(
                f,
                "gas diffusion at MIXED per-brick N_θ rides plan S11 (the \
                 ring-interface coupling of an implicit operator, COUP-3 0.4.5); \
                 refusing rather than guessing"
            ),
            Self::CutThetaGeometry => write!(
                f,
                "gas diffusion with cut geometry at N_θ > 1 rides plan S11 \
                 (θ-face apertures + per-θ wall patches); refusing rather than \
                 guessing"
            ),
            Self::NonFinite { i_r, i_z, what } => write!(
                f,
                "non-finite {what} at (i_r={i_r}, i_z={i_z}) in the gas-diffusion \
                 assembly — halt with diagnosis"
            ),
            Self::CgUnconverged { comp, resid } => write!(
                f,
                "gas class-D CG ({comp}) residual {resid:.3e} exceeds {EPS_CG_RESID:.1e} \
                 after the fixed iteration structure"
            ),
        }
    }
}

impl std::error::Error for GasDiffError {}

/// Viscous velocity condition on one domain-edge face (data, not code —
/// COUP-7's boundary objects replace these closures, the same tracked
/// deferral as the conduction/flow BCs). Closures receive `(r, θ, z, t)`
/// at the face centroid and return the wall velocity `(u_r, u_θ, u_z)`.
pub enum VelocityBc<'a> {
    /// No-slip against a wall moving at the declared velocity (zero for a
    /// stationary wall). The wall does `τ·u_wall` work — ledgered as an
    /// energy port (the Couette drive).
    NoSlip(&'a dyn Fn(f64, f64, f64, f64) -> (f64, f64, f64)),
    /// Zero viscous flux (slip, no penetration handled by the hyperbolic
    /// class's reflecting BC). A true free-slip wall: τ·n̂ = 0 AT the face.
    FreeSlip,
    /// Zero-normal-gradient continuation (symmetry/open plane): the
    /// two-point normal-gradient terms vanish, but the tangential-stress
    /// and dilatation pieces are carried one-sided at the face — an
    /// infinite-channel/outflow-plane end, NOT a wall (a FreeSlip end
    /// would truncate the real τ shear and drive edge vortices).
    Continuative,
}

/// Thermal condition on one domain-edge face.
pub enum ThermalBc<'a> {
    /// Prescribed wall temperature (half-cell two-point flux).
    Isothermal(&'a dyn Fn(f64, f64, f64, f64) -> f64),
    /// Zero conductive flux.
    Adiabatic,
}

/// Species condition on one domain-edge face. The physical occupant is
/// `ZeroFlux` (non-catalytic wall — module doc); `Prescribed` exists for
/// MMS verification, exactly like `FlowBc::Prescribed`.
pub enum SpeciesBc<'a> {
    ZeroFlux,
    Prescribed(&'a dyn Fn(f64, f64, f64, f64) -> f64),
}

pub struct FaceGasBc<'a> {
    pub velocity: VelocityBc<'a>,
    pub thermal: ThermalBc<'a>,
    pub species: SpeciesBc<'a>,
}

impl FaceGasBc<'_> {
    /// The all-fluxes-zero face (adiabatic free-slip, non-catalytic).
    pub const fn free() -> Self {
        FaceGasBc {
            velocity: VelocityBc::FreeSlip,
            thermal: ThermalBc::Adiabatic,
            species: SpeciesBc::ZeroFlux,
        }
    }
}

pub struct GasDiffBcs<'a> {
    /// Ignored where the axis face has zero area (r_min = 0).
    pub r_inner: FaceGasBc<'a>,
    pub r_outer: FaceGasBc<'a>,
    pub z_lo: FaceGasBc<'a>,
    pub z_hi: FaceGasBc<'a>,
}

/// The gas-phase diffusion operator (module doc). It holds **only** its
/// declared boundary conditions: transport arrives per cell in a
/// [`GasTransportField`] the caller refreshes from the FND-7 spine once per
/// Picard iterate (COUP-3 §3.1, 0.4.3). Before S4 this struct cached four
/// scalars copied out of the wall law; the spine is now the one owner and
/// nothing is copied.
pub struct GasDiffusion<'a> {
    pub bcs: GasDiffBcs<'a>,
}

impl<'a> GasDiffusion<'a> {
    pub fn new(bcs: GasDiffBcs<'a>) -> Self {
        GasDiffusion { bcs }
    }
}

/// Per-cell spine transport for one gas class-`D` solve, per brick at
/// N_θ = 1 — the coefficient fields of COUP-3 §3.1's variable-coefficient
/// rule. `cp` is not carried: the operator needs only the coefficients that
/// multiply a gradient (μ, k, ρD, ∂h/∂Z) and the T-solve's slope c_v.
pub struct GasTransportField {
    pub(crate) mu: BufF,
    pub(crate) k: BufF,
    pub(crate) cv: BufF,
    pub(crate) rho_d: BufF,
    pub(crate) dh_dz: BufF,
}

/// Per-brick buffer sized to the brick's own θ-plane count (`n_theta ×
/// BRICK_CELLS`, θ-plane-major — the grid's `cell_index` layout). At
/// N_θ = 1 this is the pre-S9 `BRICK_CELLS` layout bit-for-bit.
fn alloc_plane(g: &Grid) -> BufF {
    (0..g.n_bricks())
        .map(|bi| vec![0.0f64; g.brick(bi).n_theta() as usize * BRICK_CELLS])
        .collect()
}

impl GasTransportField {
    pub fn alloc(g: &Grid) -> Self {
        let mk = || alloc_plane(g);
        GasTransportField {
            mu: mk(),
            k: mk(),
            cv: mk(),
            rho_d: mk(),
            dh_dz: mk(),
        }
    }

    /// Write one gas cell's spine reading, validating it at the point of
    /// production (META-1 P6 — a coefficient that is not finite and
    /// positive must never reach an assembly, where its provenance is
    /// gone). `dh_dz` is signed and only checked finite.
    pub fn set(
        &mut self,
        bi: usize,
        local: usize,
        i_r: usize,
        i_z: usize,
        tr: &TransportProps,
    ) -> Result<(), GasDiffError> {
        for (what, v) in [
            ("spine mu", tr.mu),
            ("spine k", tr.k),
            ("spine c_v", tr.cv),
            ("spine rho*D", tr.rho_d),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(GasDiffError::NonFinite { i_r, i_z, what });
            }
        }
        if !tr.dh_dz.is_finite() {
            return Err(GasDiffError::NonFinite {
                i_r,
                i_z,
                what: "spine dh/dZ",
            });
        }
        self.mu[bi][local] = tr.mu;
        self.k[bi][local] = tr.k;
        self.cv[bi][local] = tr.cv;
        self.rho_d[bi][local] = tr.rho_d;
        self.dh_dz[bi][local] = tr.dh_dz;
        Ok(())
    }

    /// Masked (non-gas) slots: values that keep the arithmetic finite and
    /// are never read (the assembly skips masked cells and suppresses every
    /// face against one). `1.0`, not `0.0`, so a hypothetical read produces
    /// an obviously-wrong number rather than a silent zero coefficient.
    pub fn set_masked(&mut self, bi: usize, local: usize) {
        for f in [&mut self.mu, &mut self.k, &mut self.cv, &mut self.rho_d] {
            f[bi][local] = 1.0;
        }
        self.dh_dz[bi][local] = 0.0;
    }
}

/// One face's transport, from the visiting cell's side. Interior faces take
/// the **arithmetic mean** of the two cells' spine readings (SOLV-1 §3.1,
/// 0.4.1): exact in the constant-coefficient limit — which is what keeps
/// every constant-occupant fixture and certificate bit-identical — and
/// second order for the smooth transport fields a gas has. There is no
/// material discontinuity to cross, because a wall face carries no resolved
/// diffusion at all. Boundary faces take the cell's own value, one-sided,
/// exactly as their lagged cross terms already do.
#[derive(Clone, Copy)]
struct FaceTr {
    mu: f64,
    k: f64,
    rho_d: f64,
    dh_dz: f64,
}

impl FaceTr {
    fn at(tr: &GasTransportField, bi: usize, local: usize) -> Self {
        FaceTr {
            mu: tr.mu[bi][local],
            k: tr.k[bi][local],
            rho_d: tr.rho_d[bi][local],
            dh_dz: tr.dh_dz[bi][local],
        }
    }

    /// The face average. Commutative in the two cells, so the coefficient
    /// seen from either side is bit-identical — the symmetry the CG's SPD
    /// contract and the flux telescoping both rest on.
    fn between(tr: &GasTransportField, a: (usize, usize), b: (usize, usize)) -> Self {
        let avg = |f: &BufF| 0.5 * (f[a.0][a.1] + f[b.0][b.1]);
        FaceTr {
            mu: avg(&tr.mu),
            k: avg(&tr.k),
            rho_d: avg(&tr.rho_d),
            dh_dz: avg(&tr.dh_dz),
        }
    }
}

/// Per-cell working-variable fields of one gas class-D solve, per brick at
/// N_θ = 1 (plane = `BRICK_CELLS`): frozen density, the three solved
/// velocity variables (ω = u_θ/r̄), temperature, and composition.
pub struct GasOperands {
    pub rho: BufF,
    pub ur: BufF,
    pub om: BufF,
    pub uz: BufF,
    pub tt: BufF,
    pub cc: BufF,
}

impl GasOperands {
    pub fn alloc(g: &Grid) -> Self {
        let mk = || alloc_plane(g);
        GasOperands {
            rho: mk(),
            ur: mk(),
            om: mk(),
            uz: mk(),
            tt: mk(),
            cc: mk(),
        }
    }

    pub fn clone_from(&mut self, other: &GasOperands) {
        for (dst, src) in [
            (&mut self.rho, &other.rho),
            (&mut self.ur, &other.ur),
            (&mut self.om, &other.om),
            (&mut self.uz, &other.uz),
            (&mut self.tt, &other.tt),
            (&mut self.cc, &other.cc),
        ] {
            for (d, s) in dst.iter_mut().zip(src) {
                d.copy_from_slice(s);
            }
        }
    }
}

/// Scratch of the gas class-D solve: lag-state cell-centered velocity
/// gradients (the Picard-lagged cross-term operands) and the CG vectors.
/// The θ-derivative buffers (`d*_dth`, plain ∂/∂θ per cell — consumers
/// divide by the cell's own radius, the `e_thth_f` pattern) and the ω
/// meridional gradients (`dom_dr`/`dom_dz`, the τ_rθ/τ_θz θ-face limbs)
/// are written only at N_θ > 1 and never read at N_θ = 1 (the structural
/// bit-identity guard — S9).
pub struct GasWork {
    pub(crate) dur_dr: BufF,
    pub(crate) dur_dz: BufF,
    pub(crate) duz_dr: BufF,
    pub(crate) duz_dz: BufF,
    pub(crate) dom_dr: BufF,
    pub(crate) dom_dz: BufF,
    pub(crate) dur_dth: BufF,
    pub(crate) duz_dth: BufF,
    pub(crate) dom_dth: BufF,
    pub(crate) b: BufF,
    pub(crate) r: BufF,
    pub(crate) z: BufF,
    pub(crate) q: BufF,
    pub(crate) p: BufF,
    pub(crate) diag: BufF,
    pub(crate) delta: BufF,
    pub(crate) mass: BufF,
    pub(crate) ke_base: BufF,
}

impl GasWork {
    pub fn alloc(g: &Grid) -> Self {
        let mk = || alloc_plane(g);
        GasWork {
            dur_dr: mk(),
            dur_dz: mk(),
            duz_dr: mk(),
            duz_dz: mk(),
            dom_dr: mk(),
            dom_dz: mk(),
            dur_dth: mk(),
            duz_dth: mk(),
            dom_dth: mk(),
            b: mk(),
            r: mk(),
            z: mk(),
            q: mk(),
            p: mk(),
            diag: mk(),
            delta: mk(),
            mass: mk(),
            ke_base: mk(),
        }
    }
}

/// The five solved components, in the fixed solve order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum GasComp {
    Ur,
    Uz,
    Om,
    T,
    C,
}

impl GasComp {
    fn name(self) -> &'static str {
        match self {
            Self::Ur => "u_r",
            Self::Uz => "u_z",
            Self::Om => "omega (u_theta/r)",
            Self::T => "T",
            Self::C => "C",
        }
    }
}

/// Where a face's neighbor lives (the conduction assembly's pattern).
#[derive(Clone, Copy)]
enum Nbr {
    InBrick(isize),
    Cross { bi: usize, local: usize },
    Missing,
    Edge,
}

/// One face's classification for this operator.
enum FaceKind<'b> {
    /// gas↔gas with open aperture: resolved two-point + cross fluxes.
    Interior { bi: usize, local: usize, ap: f64 },
    /// Domain-edge face: the declared viscous BC (aperture-weighted).
    Boundary { bc: &'b FaceGasBc<'b>, ap: f64 },
    /// Wall-law face (solid/exterior neighbor) or sealed aperture:
    /// F_visc suppressed — the wall function owns it (SOLV-1 §3.5).
    Suppressed,
}

/// Everything constant about one face visit, from the visiting cell's side.
struct FaceGeom {
    /// Full face area (per full ring, N_θ = 1), before aperture weighting.
    area: f64,
    /// Center-to-neighbor distance (dr or dz); BC faces use half of it.
    dist: f64,
    /// Face radius (r-faces: the face's own radius; z-faces: the ring r̄).
    r_face: f64,
    /// +1 if the neighbor/wall sits on the high side of the cell.
    high: f64,
    /// Face centroid (r, z) for BC closures.
    pos: (f64, f64),
    /// True for r-direction faces (selects the 4/3-vs-1 coefficient split).
    radial: bool,
}

impl GasDiffusion<'_> {
    /// Uniform per-brick N_θ (any value) is legal since plan S9; MIXED
    /// N_θ and cut geometry at N_θ > 1 refuse (→ plan S11, module doc).
    fn validate(&self, g: &Grid) -> Result<(), GasDiffError> {
        let nt0 = g.brick(0).n_theta();
        if g.bricks().iter().any(|b| b.n_theta() != nt0) {
            return Err(GasDiffError::AzimuthalResolution);
        }
        if nt0 > 1 && g.has_cut_geometry() {
            return Err(GasDiffError::CutThetaGeometry);
        }
        Ok(())
    }

    /// The validated uniform N_θ (callers run `validate` first).
    fn n_theta(g: &Grid) -> usize {
        g.brick(0).n_theta() as usize
    }

    /// Implicit two-point face coefficient of a component [per unit
    /// gradient·area]: the symmetric core the CG matrix carries, at this
    /// face's transport.
    fn face_coef(comp: GasComp, geom: &FaceGeom, tr: &FaceTr) -> f64 {
        match comp {
            GasComp::Ur => {
                if geom.radial {
                    (4.0 / 3.0) * tr.mu
                } else {
                    tr.mu
                }
            }
            GasComp::Uz => {
                if geom.radial {
                    tr.mu
                } else {
                    (4.0 / 3.0) * tr.mu
                }
            }
            // The angular-momentum form: flux = μ·A·r_f²·Δω/d (r-faces) or
            // μ·A·r̄²·Δω/d (z-faces, both cells share r̄).
            GasComp::Om => tr.mu * geom.r_face * geom.r_face,
            GasComp::T => tr.k,
            GasComp::C => tr.rho_d,
        }
    }

    /// Whether a BC face carries a two-point term for this component
    /// (the Linear-mode matrix question — no closure is evaluated).
    fn bc_active(comp: GasComp, bc: &FaceGasBc<'_>) -> bool {
        match comp {
            GasComp::Ur | GasComp::Uz | GasComp::Om => {
                matches!(bc.velocity, VelocityBc::NoSlip(_))
            }
            GasComp::T => matches!(bc.thermal, ThermalBc::Isothermal(_)),
            GasComp::C => matches!(bc.species, SpeciesBc::Prescribed(_)),
        }
    }

    /// The wall value of a component at a BC face (`None` ⇒ zero-flux for
    /// this component). ω walls convert the declared u_θ at the face
    /// radius.
    fn bc_value(
        &self,
        comp: GasComp,
        bc: &FaceGasBc<'_>,
        geom: &FaceGeom,
        theta: f64,
        time: f64,
    ) -> Option<f64> {
        let (r, z) = geom.pos;
        match comp {
            GasComp::Ur | GasComp::Uz | GasComp::Om => match &bc.velocity {
                VelocityBc::NoSlip(f) => {
                    let (w_r, w_t, w_z) = f(r, theta, z, time);
                    Some(match comp {
                        GasComp::Ur => w_r,
                        GasComp::Uz => w_z,
                        GasComp::Om => w_t / geom.r_face,
                        _ => unreachable!(),
                    })
                }
                VelocityBc::FreeSlip | VelocityBc::Continuative => None,
            },
            GasComp::T => match &bc.thermal {
                ThermalBc::Isothermal(f) => Some(f(r, theta, z, time)),
                ThermalBc::Adiabatic => None,
            },
            GasComp::C => match &bc.species {
                SpeciesBc::Prescribed(f) => Some(f(r, theta, z, time)),
                SpeciesBc::ZeroFlux => None,
            },
        }
    }

    /// Resolve the four (r,z) faces of `(i_r, i_z)` for the visiting gas
    /// cell. Fixed order r−, r+, z−, z+ (the deterministic contract).
    /// Areas are per θ-sector (`nt` = the validated uniform N_θ — exactly
    /// the pre-S9 values at N_θ = 1); the classification itself is
    /// θ-independent (masks and apertures are per (r,z); cut geometry at
    /// N_θ > 1 refuses upstream).
    #[allow(clippy::too_many_lines, clippy::too_many_arguments)]
    fn classify_faces<'b>(
        &'b self,
        g: &Grid,
        bi: usize,
        local: usize,
        i_r: usize,
        i_z: usize,
        nt: u32,
        nbrs: &[Nbr; 4],
    ) -> [(FaceKind<'b>, FaceGeom); 4] {
        let spec = g.spec();
        let (dr, dz) = (spec.dr, spec.dz);
        let (r0, z0) = (spec.r_min, spec.z_min);
        let rbar = g.r_center(i_r);
        let zbar = g.z_center(i_z);
        let b = g.brick(bi);
        let dirs = [
            FaceDir::RMinus,
            FaceDir::RPlus,
            FaceDir::ZMinus,
            FaceDir::ZPlus,
        ];
        let geoms = [
            FaceGeom {
                area: g.face_area_r(i_r, false, nt),
                dist: dr,
                r_face: r0 + i_r as f64 * dr,
                high: -1.0,
                pos: (r0 + i_r as f64 * dr, zbar),
                radial: true,
            },
            FaceGeom {
                area: g.face_area_r(i_r, true, nt),
                dist: dr,
                r_face: r0 + (i_r + 1) as f64 * dr,
                high: 1.0,
                pos: (r0 + (i_r + 1) as f64 * dr, zbar),
                radial: true,
            },
            FaceGeom {
                area: g.face_area_z(i_r, nt),
                dist: dz,
                r_face: rbar,
                high: -1.0,
                pos: (rbar, z0 + i_z as f64 * dz),
                radial: false,
            },
            FaceGeom {
                area: g.face_area_z(i_r, nt),
                dist: dz,
                r_face: rbar,
                high: 1.0,
                pos: (rbar, z0 + (i_z + 1) as f64 * dz),
                radial: false,
            },
        ];
        let bcs = [
            &self.bcs.r_inner,
            &self.bcs.r_outer,
            &self.bcs.z_lo,
            &self.bcs.z_hi,
        ];
        let mut out: [Option<(FaceKind<'b>, FaceGeom)>; 4] = [None, None, None, None];
        for (k, geom) in geoms.into_iter().enumerate() {
            let ap = b.aperture_rz(dirs[k], local);
            let kind = match nbrs[k] {
                Nbr::Edge => FaceKind::Boundary { bc: bcs[k], ap },
                Nbr::InBrick(off) => {
                    let nl = (local as isize + off) as usize;
                    if b.mask() & (1u64 << nl) != 0 && ap > 0.0 {
                        FaceKind::Interior { bi, local: nl, ap }
                    } else {
                        FaceKind::Suppressed
                    }
                }
                Nbr::Cross { bi: nbi, local: nl } => {
                    if g.brick(nbi).mask() & (1u64 << nl) != 0 && ap > 0.0 {
                        FaceKind::Interior {
                            bi: nbi,
                            local: nl,
                            ap,
                        }
                    } else {
                        FaceKind::Suppressed
                    }
                }
                Nbr::Missing => FaceKind::Suppressed,
            };
            out[k] = Some((kind, geom));
        }
        out.map(|o| o.expect("filled"))
    }

    /// Neighbor addressing for one cell (the conduction pattern: adjacent
    /// bricks resolved by the caller once per brick).
    #[allow(clippy::similar_names)]
    fn neighbors(
        g: &Grid,
        i_r: usize,
        i_z: usize,
        lr: usize,
        lz: usize,
        nb: &[Option<usize>; 4],
    ) -> [Nbr; 4] {
        let (n_r, n_z) = (g.spec().n_r, g.spec().n_z);
        let n_rm = if i_r == 0 {
            Nbr::Edge
        } else if lr > 0 {
            Nbr::InBrick(-(BRICK as isize))
        } else {
            match nb[0] {
                Some(bi) => Nbr::Cross {
                    bi,
                    local: (BRICK - 1) * BRICK + lz,
                },
                None => Nbr::Missing,
            }
        };
        let n_rp = if i_r + 1 >= n_r {
            Nbr::Edge
        } else if lr + 1 < BRICK {
            Nbr::InBrick(BRICK as isize)
        } else {
            match nb[1] {
                Some(bi) => Nbr::Cross { bi, local: lz },
                None => Nbr::Missing,
            }
        };
        let n_zm = if i_z == 0 {
            Nbr::Edge
        } else if lz > 0 {
            Nbr::InBrick(-1)
        } else {
            match nb[2] {
                Some(bi) => Nbr::Cross {
                    bi,
                    local: lr * BRICK + (BRICK - 1),
                },
                None => Nbr::Missing,
            }
        };
        let n_zp = if i_z + 1 >= n_z {
            Nbr::Edge
        } else if lz + 1 < BRICK {
            Nbr::InBrick(1)
        } else {
            match nb[3] {
                Some(bi) => Nbr::Cross {
                    bi,
                    local: lr * BRICK,
                },
                None => Nbr::Missing,
            }
        };
        [n_rm, n_rp, n_zm, n_zp]
    }

    /// Cell-centered lag-state velocity gradients (the Picard-lagged
    /// cross-term operands): central over gas-connected neighbors,
    /// one-sided at walls/edges (first-order locally — module doc). At
    /// N_θ > 1 the meridional gradients are computed per θ-plane, the ω
    /// meridional gradients (`dom_dr`/`dom_dz` — the τ_rθ/τ_θz θ-face
    /// limbs) join them, and the plain-∂/∂θ ring gradients are central
    /// over the periodic within-brick ring (`(f[j+1] − f[j−1])/(2Δθ)`).
    /// The θ buffers are untouched at N_θ = 1 (structural bit-identity —
    /// they are also never read there).
    pub(crate) fn fill_lag_gradients(
        &self,
        g: &Grid,
        lag: &GasOperands,
        work: &mut GasWork,
    ) -> Result<(), GasDiffError> {
        self.validate(g)?;
        let nt = Self::n_theta(g);
        let (dr, dz) = (g.spec().dr, g.spec().dz);
        let dtheta = std::f64::consts::TAU / nt as f64;
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            let (br, bz) = (b.br(), b.bz());
            let nb = [
                (br > 0)
                    .then(|| g.brick_index_by_coords(br - 1, bz))
                    .flatten(),
                g.brick_index_by_coords(br + 1, bz),
                (bz > 0)
                    .then(|| g.brick_index_by_coords(br, bz - 1))
                    .flatten(),
                g.brick_index_by_coords(br, bz + 1),
            ];
            for jt in 0..nt {
                for local in 0..BRICK_CELLS {
                    let idx = jt * BRICK_CELLS + local;
                    work.dur_dr[bi][idx] = 0.0;
                    work.dur_dz[bi][idx] = 0.0;
                    work.duz_dr[bi][idx] = 0.0;
                    work.duz_dz[bi][idx] = 0.0;
                    work.dom_dr[bi][idx] = 0.0;
                    work.dom_dz[bi][idx] = 0.0;
                    if b.mask() & (1u64 << local) == 0 {
                        continue;
                    }
                    let (lr, lz) = (local / BRICK, local % BRICK);
                    let (i_r, i_z) = b.global_rz(local);
                    let nbrs = Self::neighbors(g, i_r, i_z, lr, lz, &nb);
                    let dirs = [
                        FaceDir::RMinus,
                        FaceDir::RPlus,
                        FaceDir::ZMinus,
                        FaceDir::ZPlus,
                    ];
                    // Gas-connected neighbor value of a field (same
                    // θ-plane — masks/apertures are per (r,z)), else None.
                    let val = |n: Nbr, f: &BufF, dir: FaceDir| -> Option<f64> {
                        if b.aperture_rz(dir, local) <= 0.0 {
                            return None;
                        }
                        match n {
                            Nbr::InBrick(off) => {
                                let nl = (local as isize + off) as usize;
                                (b.mask() & (1u64 << nl) != 0).then(|| f[bi][jt * BRICK_CELLS + nl])
                            }
                            Nbr::Cross { bi: nbi, local: nl } => {
                                (g.brick(nbi).mask() & (1u64 << nl) != 0)
                                    .then(|| f[nbi][jt * BRICK_CELLS + nl])
                            }
                            _ => None,
                        }
                    };
                    let deriv = |f: &BufF, k_lo: usize, k_hi: usize, d: f64| -> f64 {
                        let c = f[bi][idx];
                        let lo = val(nbrs[k_lo], f, dirs[k_lo]);
                        let hi = val(nbrs[k_hi], f, dirs[k_hi]);
                        match (lo, hi) {
                            (Some(a), Some(bv)) => (bv - a) / (2.0 * d),
                            (None, Some(bv)) => (bv - c) / d,
                            (Some(a), None) => (c - a) / d,
                            // No gas neighbor either way in this direction.
                            // Correct (not a guess) for the symmetric cases
                            // this reaches today — a quasi-1-D fixture or a
                            // free-slip single-cell span, where the gradient
                            // IS zero. A genuinely under-resolved gas island
                            // (one-cell gap between walls) would also land
                            // here and silently lose its dilatation term; the
                            // FND-3 PLIC/refinement wave owns that geometry
                            // class and should refuse it at build time
                            // (S3 review finding — recorded, not cured here).
                            (None, None) => 0.0,
                        }
                    };
                    work.dur_dr[bi][idx] = deriv(&lag.ur, 0, 1, dr);
                    work.dur_dz[bi][idx] = deriv(&lag.ur, 2, 3, dz);
                    work.duz_dr[bi][idx] = deriv(&lag.uz, 0, 1, dr);
                    work.duz_dz[bi][idx] = deriv(&lag.uz, 2, 3, dz);
                    if nt > 1 {
                        work.dom_dr[bi][idx] = deriv(&lag.om, 0, 1, dr);
                        work.dom_dz[bi][idx] = deriv(&lag.om, 2, 3, dz);
                        // Periodic ring central differences, plain ∂/∂θ
                        // (consumers divide by their own radius — the
                        // `e_thth_f` per-cell-radius pattern).
                        let jm = (jt + nt - 1) % nt;
                        let jp = (jt + 1) % nt;
                        let ring = |f: &BufF| -> f64 {
                            (f[bi][jp * BRICK_CELLS + local] - f[bi][jm * BRICK_CELLS + local])
                                / (2.0 * dtheta)
                        };
                        work.dur_dth[bi][idx] = ring(&lag.ur);
                        work.duz_dth[bi][idx] = ring(&lag.uz);
                        work.dom_dth[bi][idx] = ring(&lag.om);
                    }
                }
            }
        }
        Ok(())
    }

    /// Implicit θ-θ two-point face coefficient of a component (COUP-3
    /// 0.4.5's fixed design): μ for u_r/u_z (the τ_rθ/τ_θz implicit
    /// limbs), (4/3)μ·r̄² for ω (the τ_θθ limb in the angular-momentum
    /// form — both ring cells share r̄ exactly, so the coefficient is
    /// symmetric), k for T, ρD for C. Same units contract as `face_coef`
    /// with the θ-face distance r̄·Δθ.
    fn face_coef_theta(comp: GasComp, rbar: f64, tr: &FaceTr) -> f64 {
        match comp {
            GasComp::Ur | GasComp::Uz => tr.mu,
            GasComp::Om => (4.0 / 3.0) * tr.mu * rbar * rbar,
            GasComp::T => tr.k,
            GasComp::C => tr.rho_d,
        }
    }

    /// The matrix apply of one component's implicit core: `out = L·x` in
    /// solve units (force for u_r/u_z, angular-momentum torque for ω, W
    /// for T/C-flux), constants dropped exactly (the Linear discipline of
    /// the solid assembly). `diag` receives ∂out_i/∂x_i (≤ 0). At
    /// N_θ > 1 each cell additionally carries its two periodic θ-θ core
    /// faces (`face_coef_theta`); the whole θ block is structurally
    /// skipped at N_θ = 1 (bit identity — module doc).
    pub(crate) fn apply_linear(
        &self,
        g: &Grid,
        comp: GasComp,
        tr: &GasTransportField,
        x: &BufF,
        out: &mut BufF,
        mut diag: Option<&mut BufF>,
    ) {
        let nt = Self::n_theta(g);
        let ntu = nt as u32;
        let a_th = g.face_area_theta();
        let dtheta = std::f64::consts::TAU / nt as f64;
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            let (br, bz) = (b.br(), b.bz());
            let nb = [
                (br > 0)
                    .then(|| g.brick_index_by_coords(br - 1, bz))
                    .flatten(),
                g.brick_index_by_coords(br + 1, bz),
                (bz > 0)
                    .then(|| g.brick_index_by_coords(br, bz - 1))
                    .flatten(),
                g.brick_index_by_coords(br, bz + 1),
            ];
            out[bi].fill(0.0);
            if let Some(d) = diag.as_deref_mut() {
                d[bi].fill(0.0);
            }
            for jt in 0..nt {
                for local in 0..BRICK_CELLS {
                    if b.mask() & (1u64 << local) == 0 {
                        continue;
                    }
                    let idx = jt * BRICK_CELLS + local;
                    let (lr, lz) = (local / BRICK, local % BRICK);
                    let (i_r, i_z) = b.global_rz(local);
                    let nbrs = Self::neighbors(g, i_r, i_z, lr, lz, &nb);
                    let faces = self.classify_faces(g, bi, local, i_r, i_z, ntu, &nbrs);
                    let x_c = x[bi][idx];
                    let mut acc = 0.0f64;
                    let mut dg = 0.0f64;
                    for (kind, geom) in &faces {
                        // The face's transport must be formed EXACTLY as
                        // `assemble_rates` forms it, or the CG solves a
                        // different operator than the composition applies (the
                        // S3 review wave verified that identity by
                        // finite-differencing the true Jacobian; variable
                        // coefficients must not break it).
                        let ftr = match kind {
                            FaceKind::Interior {
                                bi: nbi, local: nl, ..
                            } => FaceTr::between(tr, (bi, idx), (*nbi, jt * BRICK_CELLS + nl)),
                            _ => FaceTr::at(tr, bi, idx),
                        };
                        let coef = Self::face_coef(comp, geom, &ftr);
                        match kind {
                            FaceKind::Interior {
                                bi: nbi,
                                local: nl,
                                ap,
                            } => {
                                let a = geom.area * ap;
                                let x_n = x[*nbi][jt * BRICK_CELLS + nl];
                                acc += coef * a * (x_n - x_c) / geom.dist;
                                dg -= coef * a / geom.dist;
                            }
                            FaceKind::Boundary { bc, ap } => {
                                // Linear mode: the wall VALUE is a constant —
                                // only the −x_c part of the two-point form
                                // survives; zero-flux BCs contribute nothing.
                                if geom.area > 0.0 && Self::bc_active(comp, bc) {
                                    let a = geom.area * ap;
                                    let c = coef * a / (0.5 * geom.dist);
                                    acc -= c * x_c;
                                    dg -= c;
                                }
                            }
                            FaceKind::Suppressed => {}
                        }
                    }
                    if nt > 1 {
                        // The two periodic θ-θ core faces (θ−, θ+ — fixed
                        // order). Always interior (no θ domain edge; cut
                        // geometry at N_θ > 1 refuses upstream, so the
                        // aperture is identically 1).
                        let rbar = g.r_center(i_r);
                        let dist = rbar * dtheta;
                        let jm = (jt + nt - 1) % nt;
                        let jp = (jt + 1) % nt;
                        for jn in [jm, jp] {
                            let nidx = jn * BRICK_CELLS + local;
                            let ftr = FaceTr::between(tr, (bi, idx), (bi, nidx));
                            let coef = Self::face_coef_theta(comp, rbar, &ftr);
                            acc += coef * a_th * (x[bi][nidx] - x_c) / dist;
                            dg -= coef * a_th / dist;
                        }
                    }
                    if comp == GasComp::Ur {
                        // The negative-definite geometric diagonal of −τ_θθ/r:
                        // −(4/3)μ·u_r/r̄ · geo·κV with the metric-consistent
                        // geo = (A_out − A_in)/V = 1/r̄ (SOLV-1 §3.3 pattern).
                        let vol = g.cell_volume(i_r, ntu);
                        let geo =
                            (g.face_area_r(i_r, true, ntu) - g.face_area_r(i_r, false, ntu)) / vol;
                        let kv = b.kappa_rz(local) * vol;
                        let c = (4.0 / 3.0) * tr.mu[bi][idx] * geo / g.r_center(i_r) * kv;
                        acc -= c * x_c;
                        dg -= c;
                    }
                    out[bi][idx] = acc;
                    if let Some(d) = diag.as_deref_mut() {
                        d[bi][idx] = dg;
                    }
                }
            }
        }
    }

    /// The full Affine physics assembly: all five conserved rates per cell
    /// (mass slot exactly 0), in **conserved-density-rate units** (the
    /// flow-ledger convention: `rate = (Σ face terms + sources)/(κV)`;
    /// θ-momentum from the angular-momentum form divided by r̄κV), plus
    /// the COUP-2 ledger lines. Implicit-core parts read `sol`; the
    /// Picard-lagged cross terms read the gradients in `work` (filled from
    /// the lag state) and `lag`'s fields; work fluxes read `sol`
    /// velocities. `time` feeds the BC closures.
    #[allow(clippy::too_many_lines, clippy::too_many_arguments)]
    pub(crate) fn assemble_rates(
        &self,
        g: &Grid,
        sol: &GasOperands,
        lag: &GasOperands,
        work: &GasWork,
        tr: &GasTransportField,
        time: f64,
        rates: &mut [Vec<Cons>],
        mut ledger: Option<&mut FlowLedger>,
    ) -> Result<(), GasDiffError> {
        self.validate(g)?;
        let nt = Self::n_theta(g);
        let ntu = nt as u32;
        let a_th = g.face_area_theta();
        let dtheta = std::f64::consts::TAU / nt as f64;
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            let (br, bz) = (b.br(), b.bz());
            let nb = [
                (br > 0)
                    .then(|| g.brick_index_by_coords(br - 1, bz))
                    .flatten(),
                g.brick_index_by_coords(br + 1, bz),
                (bz > 0)
                    .then(|| g.brick_index_by_coords(br, bz - 1))
                    .flatten(),
                g.brick_index_by_coords(br, bz + 1),
            ];
            for cell in rates[bi].iter_mut() {
                *cell = [0.0; NCOMP];
            }
            for jt in 0..nt {
                for local in 0..BRICK_CELLS {
                    if b.mask() & (1u64 << local) == 0 {
                        continue;
                    }
                    let idx = jt * BRICK_CELLS + local;
                    let (lr, lz) = (local / BRICK, local % BRICK);
                    let (i_r, i_z) = b.global_rz(local);
                    let rbar = g.r_center(i_r);
                    let vol = g.cell_volume(i_r, ntu);
                    let kv = b.kappa_rz(local) * vol;
                    let nbrs = Self::neighbors(g, i_r, i_z, lr, lz, &nb);
                    let faces = self.classify_faces(g, bi, local, i_r, i_z, ntu, &nbrs);

                    // The visiting cell's operands.
                    let at = |f: &BufF| f[bi][idx];
                    let (ur_c, om_c, uz_c, tt_c, cc_c) = (
                        at(&sol.ur),
                        at(&sol.om),
                        at(&sol.uz),
                        at(&sol.tt),
                        at(&sol.cc),
                    );
                    for (what, v, positive) in [
                        ("velocity", ur_c.abs() + om_c.abs() + uz_c.abs(), false),
                        // Temperature must be POSITIVE, not merely finite: a
                        // T ≤ 0 out of the implicit solve would otherwise flow
                        // into the k∇T fluxes and downstream seams unremarked
                        // (META-1 P6 — refuse, never carry; S3 review finding).
                        ("temperature", tt_c, true),
                        ("composition operand", cc_c.abs(), false),
                    ] {
                        if !v.is_finite() || (positive && v <= 0.0) {
                            return Err(GasDiffError::NonFinite { i_r, i_z, what });
                        }
                    }
                    // Lag-state cell value for the cross terms. At N_θ > 1 the
                    // full e_θθ carries the ∂ω/∂θ dilatation limb (COUP-3
                    // 0.4.5): e_θθ = u_r/r + (1/r)∂u_θ/∂θ = u_r/r + ∂ω/∂θ —
                    // structurally skipped at N_θ = 1 (bit identity).
                    let mut e_thth_c = lag.ur[bi][idx] / rbar;
                    if nt > 1 {
                        e_thth_c += work.dom_dth[bi][idx];
                    }

                    // Accumulators in total units; θ-momentum in λ (angular-
                    // momentum) units. Fixed face order = the deterministic
                    // accumulation contract.
                    let mut tot = [0.0f64; NCOMP];
                    let mut tot_lam = 0.0f64;
                    let mut lam_abs = 0.0f64;
                    // The domain-edge (port) part of `tot`, kept separately so
                    // the ledger's gross scale does not double-count it.
                    let mut bc_port = [0.0f64; NCOMP];

                    for (kind, geom) in &faces {
                        match kind {
                            FaceKind::Interior {
                                bi: nbi,
                                local: nl,
                                ap,
                            } => {
                                let nidx = jt * BRICK_CELLS + nl;
                                let ftr = FaceTr::between(tr, (bi, idx), (*nbi, nidx));
                                let two_thirds_mu = (2.0 / 3.0) * ftr.mu;
                                let a = geom.area * ap;
                                let d = geom.dist;
                                let (ur_n, om_n, uz_n, tt_n, cc_n) = (
                                    sol.ur[*nbi][nidx],
                                    sol.om[*nbi][nidx],
                                    sol.uz[*nbi][nidx],
                                    sol.tt[*nbi][nidx],
                                    sol.cc[*nbi][nidx],
                                );
                                // Signed two-point gradients along the face
                                // normal, oriented +r/+z (exact negation seen
                                // from the other side ⇒ exact telescoping).
                                let s = geom.high;
                                let g_ur = s * (ur_n - ur_c) / d;
                                let g_om = s * (om_n - om_c) / d;
                                let g_uz = s * (uz_n - uz_c) / d;
                                let g_tt = s * (tt_n - tt_c) / d;
                                let g_cc = s * (cc_n - cc_c) / d;
                                // Face-averaged lag quantities for the cross
                                // pieces (identical from both sides). For
                                // radial faces the neighbor's e_θθ uses its
                                // own ring radius; z-face neighbors share r̄.
                                let avg_w = |w: &BufF| 0.5 * (w[bi][idx] + w[*nbi][nidx]);
                                let e_thth_f = {
                                    let (n_ir, _) = g.brick(*nbi).global_rz(*nl);
                                    let r_n = if geom.radial { g.r_center(n_ir) } else { rbar };
                                    let mut e_n = lag.ur[*nbi][nidx] / r_n;
                                    if nt > 1 {
                                        // The neighbor's own ∂ω/∂θ dilatation
                                        // limb (its e_θθ, like e_thth_c above).
                                        e_n += work.dom_dth[*nbi][nidx];
                                    }
                                    0.5 * (e_thth_c + e_n)
                                };
                                let e_rr_f = avg_w(&work.dur_dr);
                                let e_zz_f = avg_w(&work.duz_dz);
                                let dur_dz_f = avg_w(&work.dur_dz);
                                let duz_dr_f = avg_w(&work.duz_dr);
                                // Face-averaged velocities for the work flux.
                                let ur_f = 0.5 * (ur_c + ur_n);
                                let uz_f = 0.5 * (uz_c + uz_n);
                                let ut_f = if geom.radial {
                                    let (n_ir, _) = g.brick(*nbi).global_rz(*nl);
                                    0.5 * (om_c * rbar + sol.om[*nbi][nidx] * g.r_center(n_ir))
                                } else {
                                    0.5 * (om_c + om_n) * rbar
                                };
                                // The stress components at the face (implicit
                                // normal-gradient parts at sol; cross at lag).
                                // At N_θ > 1 the meridional-face θ-limbs join
                                // as lag terms (COUP-3 0.4.5): τ_rθ gains
                                // (1/r)∂u_r/∂θ, τ_θz gains (1/r)∂u_z/∂θ — the
                                // per-cell-radius face average, e_θθ pattern.
                                let (f_mr, f_mz, tau_rr_or_zz, tau_rz, tau_th);
                                if geom.radial {
                                    let tau_rr = (4.0 / 3.0) * ftr.mu * g_ur
                                        - two_thirds_mu * (e_thth_f + e_zz_f);
                                    let t_rz = ftr.mu * (g_uz + dur_dz_f);
                                    let mut t_rth = ftr.mu * geom.r_face * g_om;
                                    if nt > 1 {
                                        let (n_ir, _) = g.brick(*nbi).global_rz(*nl);
                                        let dur_dth_r_f = 0.5
                                            * (work.dur_dth[bi][idx] / rbar
                                                + work.dur_dth[*nbi][nidx] / g.r_center(n_ir));
                                        t_rth += ftr.mu * dur_dth_r_f;
                                    }
                                    f_mr = a * tau_rr;
                                    f_mz = a * t_rz;
                                    tau_rr_or_zz = tau_rr;
                                    tau_rz = t_rz;
                                    tau_th = t_rth;
                                } else {
                                    let tau_zz = (4.0 / 3.0) * ftr.mu * g_uz
                                        - two_thirds_mu * (e_rr_f + e_thth_f);
                                    let t_rz = ftr.mu * (g_ur + duz_dr_f);
                                    let mut t_thz = ftr.mu * rbar * g_om;
                                    if nt > 1 {
                                        t_thz += ftr.mu * avg_w(&work.duz_dth) / rbar;
                                    }
                                    f_mr = a * t_rz;
                                    f_mz = a * tau_zz;
                                    tau_rr_or_zz = tau_zz;
                                    tau_rz = t_rz;
                                    tau_th = t_thz;
                                }
                                // Signed accumulation: G is the flux vector
                                // component along +r/+z; a high face adds
                                // +G·A, a low face −G·A (divergence).
                                tot[I_MR] += s * f_mr;
                                tot[I_MZ] += s * f_mz;
                                // θ: the angular-momentum (λ = ρu_θr) flux is
                                // A·r·τ_(rθ|θz) with the face's own radius —
                                // `geom.r_face` is the face radius on r-faces
                                // (exact) and r̄ on z-faces (the consistent
                                // second-order lumped form of ∫r τ dA — the
                                // exact z-face moment is (r̄² + Δr²/12)ΔrΔθ;
                                // the r̄² form matches the ω-solve's ρr̄²κV
                                // inertia and telescopes exactly — S9 review
                                // note).
                                let f_lam = a * geom.r_face * tau_th;
                                tot_lam += s * f_lam;
                                lam_abs += (f_lam / rbar).abs();
                                // Energy: work + conduction.
                                let g_e = if geom.radial {
                                ur_f * tau_rr_or_zz + ut_f * tau_th + uz_f * tau_rz
                            } else {
                                ur_f * tau_rz + ut_f * tau_th + uz_f * tau_rr_or_zz
                            } + ftr.k * g_tt
                                // The species-enthalpy diffusion flux
                                // `Σ h_k j_k` (SOLV-1 §3.1, 0.4.1): with one
                                // composition coordinate it is exactly
                                // (∂h/∂Z)|_{p,T}·j_Z, and it is identically
                                // zero on a single-composition gas — which is
                                // what kept the S3 constant occupant honest.
                                // Its ∇T limb is NOT here: that limb is already
                                // inside the spine's effective conductivity
                                // (crate::transport module doc), so adding it
                                // again would double-count one flux.
                                + ftr.rho_d * ftr.dh_dz * g_cc;
                                tot[I_EN] += s * a * g_e;
                                // Species.
                                tot[I_RC] += s * a * ftr.rho_d * g_cc;
                            }
                            FaceKind::Boundary { bc, ap } => {
                                if geom.area == 0.0 {
                                    continue; // the axis face drops out
                                }
                                let ftr = FaceTr::at(tr, bi, idx);
                                let two_thirds_mu = (2.0 / 3.0) * ftr.mu;
                                let a = geom.area * ap;
                                let half = 0.5 * geom.dist;
                                let theta = Grid::theta_center(jt as u32, ntu);
                                let s = geom.high;
                                // One-sided lag pieces at the boundary (module
                                // doc: first-order locally).
                                let e_zz_f = work.duz_dz[bi][idx];
                                let e_rr_f = work.dur_dr[bi][idx];
                                let e_thth_f = e_thth_c;
                                let dur_dz_f = work.dur_dz[bi][idx];
                                let duz_dr_f = work.duz_dr[bi][idx];
                                let mut port = [0.0f64; NCOMP];
                                let mut port_lam = 0.0f64;
                                // Viscous terms per the declared velocity BC:
                                // NoSlip = two-point normal gradients against
                                // the wall values + one-sided cross pieces,
                                // work at the WALL velocity; Continuative =
                                // zero normal gradient, cross pieces only,
                                // work at the CELL velocity; FreeSlip = none.
                                let visc = match &bc.velocity {
                                    VelocityBc::FreeSlip => None,
                                    VelocityBc::NoSlip(f) => {
                                        let (r, z) = geom.pos;
                                        let (wr, wt, wz) = f(r, theta, z, time);
                                        Some((
                                            s * (wr - ur_c) / half,
                                            s * (wz - uz_c) / half,
                                            s * (wt / geom.r_face - om_c) / half,
                                            wr,
                                            wt,
                                            wz,
                                        ))
                                    }
                                    VelocityBc::Continuative => {
                                        Some((0.0, 0.0, 0.0, ur_c, om_c * rbar, uz_c))
                                    }
                                };
                                if let Some((g_ur, g_uz, g_om, wr, wt, wz)) = visc {
                                    let (tau_nn, tau_rz, mut tau_th);
                                    if geom.radial {
                                        tau_nn = (4.0 / 3.0) * ftr.mu * g_ur
                                            - two_thirds_mu * (e_thth_f + e_zz_f);
                                        tau_rz = ftr.mu * (g_uz + dur_dz_f);
                                        tau_th = ftr.mu * geom.r_face * g_om;
                                        if nt > 1 {
                                            // One-sided τ_rθ θ-limb at the
                                            // boundary (cell radius — module
                                            // doc: first-order locally).
                                            tau_th += ftr.mu * work.dur_dth[bi][idx] / rbar;
                                        }
                                        port[I_MR] += s * a * tau_nn;
                                        port[I_MZ] += s * a * tau_rz;
                                    } else {
                                        tau_nn = (4.0 / 3.0) * ftr.mu * g_uz
                                            - two_thirds_mu * (e_rr_f + e_thth_f);
                                        tau_rz = ftr.mu * (g_ur + duz_dr_f);
                                        tau_th = ftr.mu * rbar * g_om;
                                        if nt > 1 {
                                            // One-sided τ_θz θ-limb.
                                            tau_th += ftr.mu * work.duz_dth[bi][idx] / rbar;
                                        }
                                        port[I_MZ] += s * a * tau_nn;
                                        port[I_MR] += s * a * tau_rz;
                                    }
                                    port_lam += s * a * geom.r_face * tau_th;
                                    // The face's work at its declared/continued
                                    // velocity (the Couette drive when a NoSlip
                                    // wall moves) — an energy port.
                                    let g_e_work = if geom.radial {
                                        wr * tau_nn + wt * tau_th + wz * tau_rz
                                    } else {
                                        wr * tau_rz + wt * tau_th + wz * tau_nn
                                    };
                                    port[I_EN] += s * a * g_e_work;
                                }
                                // Thermal + species two-point terms (declared
                                // independently of the velocity condition).
                                if let Some(w_tt) = self.bc_value(GasComp::T, bc, geom, theta, time)
                                {
                                    let g_tt = s * (w_tt - tt_c) / half;
                                    port[I_EN] += s * a * ftr.k * g_tt;
                                }
                                if let Some(wc) = self.bc_value(GasComp::C, bc, geom, theta, time) {
                                    let g_cc = s * (wc - cc_c) / half;
                                    port[I_RC] += s * a * ftr.rho_d * g_cc;
                                    // ...and the enthalpy that flux carries.
                                    port[I_EN] += s * a * ftr.rho_d * ftr.dh_dz * g_cc;
                                }
                                for ((t_k, b_k), p_k) in
                                    tot.iter_mut().zip(bc_port.iter_mut()).zip(&port)
                                {
                                    *t_k += p_k;
                                    *b_k += p_k;
                                }
                                tot_lam += port_lam;
                                lam_abs += (port_lam / rbar).abs();
                                if let Some(l) = ledger.as_deref_mut() {
                                    #[allow(clippy::needless_range_loop)]
                                    // kk indexes two ledger arrays
                                    for kk in 0..NCOMP {
                                        if kk == I_MT {
                                            continue; // θ goes through the src lines
                                        }
                                        l.port_net[kk] += port[kk];
                                        l.port_abs[kk] += port[kk].abs();
                                    }
                                }
                            }
                            FaceKind::Suppressed => {}
                        }
                    }

                    if nt > 1 {
                        // The two periodic θ-faces (θ−, θ+ — fixed order after
                        // the r/z quartet; COUP-3 0.4.5's full θ-stress
                        // tensor). Always interior (periodic ring; no cut
                        // geometry at N_θ > 1), aperture identically 1. The
                        // per-unit-arc two-point gradients read `sol` (the
                        // implicit limbs — arithmetic mirrors `apply_linear`'s
                        // θ core); the curvature/cross limbs read the lag
                        // gradients, face-averaged (both cells share r̄, so
                        // every face quantity is identical from either side —
                        // exact telescoping).
                        let dist = rbar * dtheta;
                        let jm = (jt + nt - 1) % nt;
                        let jp = (jt + 1) % nt;
                        for (jn, s) in [(jm, -1.0f64), (jp, 1.0f64)] {
                            let nidx = jn * BRICK_CELLS + local;
                            let ftr = FaceTr::between(tr, (bi, idx), (bi, nidx));
                            let two_thirds_mu = (2.0 / 3.0) * ftr.mu;
                            let a = a_th;
                            // Signed per-unit-arc gradients along +θ:
                            // g_X = (1/r)∂X/∂θ at the face.
                            let g_ur = s * (sol.ur[bi][nidx] - ur_c) / dist;
                            let g_om = s * (sol.om[bi][nidx] - om_c) / dist;
                            let g_uz = s * (sol.uz[bi][nidx] - uz_c) / dist;
                            let g_tt = s * (sol.tt[bi][nidx] - tt_c) / dist;
                            let g_cc = s * (sol.cc[bi][nidx] - cc_c) / dist;
                            // Face-averaged lag pieces (same (r,z) cell pair).
                            let avg_w = |w: &BufF| 0.5 * (w[bi][idx] + w[bi][nidx]);
                            let e_rr_f = avg_w(&work.dur_dr);
                            let e_zz_f = avg_w(&work.duz_dz);
                            let dom_dr_f = avg_w(&work.dom_dr);
                            let dom_dz_f = avg_w(&work.dom_dz);
                            let ur_f_lag = 0.5 * (lag.ur[bi][idx] + lag.ur[bi][nidx]);
                            // The θ-column stresses at the face:
                            // τ_rθ = μ[(1/r)∂u_r/∂θ + r ∂ω/∂r]
                            // τ_θθ = (4/3)μ[(1/r)∂u_θ/∂θ + u_r/r] − ⅔μ(e_rr+e_zz)
                            //        with (1/r)∂u_θ/∂θ = ∂ω/∂θ = r̄·g_om
                            // τ_θz = μ[(1/r)∂u_z/∂θ + r ∂ω/∂z]
                            let tau_rth = ftr.mu * (g_ur + rbar * dom_dr_f);
                            let tau_thth = (4.0 / 3.0) * ftr.mu * (rbar * g_om + ur_f_lag / rbar)
                                - two_thirds_mu * (e_rr_f + e_zz_f);
                            let tau_thz = ftr.mu * (g_uz + rbar * dom_dz_f);
                            tot[I_MR] += s * a * tau_rth;
                            tot[I_MZ] += s * a * tau_thz;
                            // Angular momentum: the θ-face torque flux is
                            // A_θ·r̄·τ_θθ (force τ_θθ × arm r̄).
                            let f_lam = a * rbar * tau_thth;
                            tot_lam += s * f_lam;
                            lam_abs += (f_lam / rbar).abs();
                            // Energy: work + conduction + species enthalpy —
                            // the same total-energy bookkeeping as the r/z
                            // faces (dissipation from the KE ledger).
                            let ur_f = 0.5 * (ur_c + sol.ur[bi][nidx]);
                            let ut_f = 0.5 * (om_c + sol.om[bi][nidx]) * rbar;
                            let uz_f = 0.5 * (uz_c + sol.uz[bi][nidx]);
                            let g_e = ur_f * tau_rth
                                + ut_f * tau_thth
                                + uz_f * tau_thz
                                + ftr.k * g_tt
                                + ftr.rho_d * ftr.dh_dz * g_cc;
                            tot[I_EN] += s * a * g_e;
                            tot[I_RC] += s * a * ftr.rho_d * g_cc;
                        }
                    }

                    // −τ_θθ/r volume source of the r-momentum (implicit
                    // diagonal + lagged compressible correction), with the
                    // metric-consistent 1/r̄ (SOLV-1 §3.3 pattern). At N_θ > 1
                    // e_θθ's lagged ∂ω/∂θ dilatation limb joins (the
                    // "(2μ/r²)∂u_θ/∂θ-class basis term" — the m = 1 curvature
                    // partner of the θ-θ core).
                    let geo =
                        (g.face_area_r(i_r, true, ntu) - g.face_area_r(i_r, false, ntu)) / vol;
                    let e_rr_c = work.dur_dr[bi][idx];
                    let e_zz_c = work.duz_dz[bi][idx];
                    let mu_c = tr.mu[bi][idx];
                    let mut tau_thth =
                        (4.0 / 3.0) * mu_c * (ur_c / rbar) - (2.0 / 3.0) * mu_c * (e_rr_c + e_zz_c);
                    if nt > 1 {
                        tau_thth += (4.0 / 3.0) * mu_c * work.dom_dth[bi][idx];
                    }
                    let src_mr = -tau_thth * geo * kv;
                    tot[I_MR] += src_mr;

                    // Commit: conserved-density rates (θ from λ/(r̄κV)).
                    let inv_kv = 1.0 / kv;
                    let rate = &mut rates[bi][idx];
                    rate[I_MR] = tot[I_MR] * inv_kv;
                    rate[I_MT] = tot_lam / (rbar * kv);
                    rate[I_MZ] = tot[I_MZ] * inv_kv;
                    rate[I_EN] = tot[I_EN] * inv_kv;
                    rate[I_RC] = tot[I_RC] * inv_kv;

                    if let Some(l) = ledger.as_deref_mut() {
                        // Interior fluxes telescope for r/z-momentum, energy,
                        // species — only BC ports (accumulated in the Boundary
                        // arm above) and the volume/θ sources are ledgered.
                        l.src_net[I_MR] += src_mr;
                        l.src_abs[I_MR] += src_mr.abs();
                        l.src_net[I_MT] += tot_lam / rbar;
                        l.src_abs[I_MT] += lam_abs;
                        // Gross magnitude for S[q]: the cell's NET applied
                        // increment. The boundary ports are already in
                        // `port_abs`; adding `tot` (which contains them) would
                        // double-count and silently LOOSEN `TOL_AUDIT` on every
                        // edge-touching cell (S3 review finding) — charge only
                        // the interior part here.
                        for kk in [I_MR, I_MZ, I_EN, I_RC] {
                            l.port_abs[kk] += (tot[kk] - bc_port[kk]).abs();
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// One component's fixed-structure Jacobi-CG solve in δ-form (the S2
    /// solid-solve pattern, on module-owned buffers): solves
    /// `(mass − wqnew·L)·δ = b`, then `x += δ` over gas cells. `mass` must
    /// hold the component's per-cell mass (module doc); returns
    /// (iterations, terminal relative residual). NaN-safe acceptance.
    pub(crate) fn cg_solve(
        &self,
        g: &Grid,
        comp: GasComp,
        tr: &GasTransportField,
        wqnew: f64,
        x: &mut BufF,
        work: &mut GasWork,
    ) -> Result<(usize, f64), GasDiffError> {
        let nb = g.n_bricks();
        let masked = |bi: usize| g.brick(bi).mask();
        // Plane-wide iteration (θ-plane-major; the (r,z) gas mask repeats
        // per plane — `idx % BRICK_CELLS`). At N_θ = 1 this is the pre-S9
        // loop bit-for-bit (same length, same order). The dot's partials
        // are per (brick, θ-plane) into the fixed-shape tree: at N_θ = 1
        // that is exactly the per-brick list of before, and on an
        // axisymmetric N_θ = 2^k world the equal per-plane partials
        // combine EXACTLY (pairwise doubling is exact), so the CG's
        // accept/iterate decisions — and therefore the whole solve — stay
        // bit-identical per plane to the N_θ = 1 march (the S9 symmetry
        // gate's arithmetic basis).
        let dot = |a: &BufF, c: &BufF| -> f64 {
            let partials: Vec<f64> = (0..nb)
                .map(|bi| {
                    let m = masked(bi);
                    let nt_b = a[bi].len() / BRICK_CELLS;
                    // Per-plane sums combined by the same fixed tree
                    // WITHIN the brick (N_θ is a power of two on the
                    // ladder, so equal planes double pairwise exactly),
                    // then the per-brick partials enter the outer tree —
                    // whose shape is therefore N_θ-independent.
                    let planes: Vec<f64> = (0..nt_b)
                        .map(|jt| {
                            let mut acc = 0.0f64;
                            for local in 0..BRICK_CELLS {
                                if m & (1u64 << local) != 0 {
                                    let idx = jt * BRICK_CELLS + local;
                                    acc += a[bi][idx] * c[bi][idx];
                                }
                            }
                            acc
                        })
                        .collect();
                    tree_combine(&planes)
                })
                .collect();
            tree_combine(&partials)
        };

        // A's diagonal: mass − wqnew·diag(L)  (diag(L) ≤ 0 ⇒ positive).
        self.apply_linear(g, comp, tr, x, &mut work.q, Some(&mut work.diag));
        for bi in 0..nb {
            let m = masked(bi);
            for idx in 0..work.diag[bi].len() {
                if m & (1u64 << (idx % BRICK_CELLS)) == 0 {
                    work.diag[bi][idx] = 1.0;
                } else {
                    work.diag[bi][idx] = work.mass[bi][idx] - wqnew * work.diag[bi][idx];
                }
            }
        }
        for bi in 0..nb {
            work.delta[bi].fill(0.0);
            work.r[bi].copy_from_slice(&work.b[bi]);
        }
        let b_norm2 = dot(&work.b, &work.b);
        let eps2 = EPS_CG_RESID * EPS_CG_RESID * b_norm2;
        for bi in 0..nb {
            for idx in 0..work.z[bi].len() {
                work.z[bi][idx] = work.r[bi][idx] / work.diag[bi][idx];
            }
            work.p[bi].copy_from_slice(&work.z[bi]);
        }
        let mut rz = dot(&work.r, &work.z);
        let mut r_norm2 = dot(&work.r, &work.r);
        let mut iters = 0usize;
        while iters < N_CG_ITERS_MAX && r_norm2 > eps2 && rz > 0.0 {
            // q = (mass − wqnew·L)·p.
            self.apply_linear(g, comp, tr, &work.p, &mut work.q, None);
            for bi in 0..nb {
                let m = masked(bi);
                for idx in 0..work.q[bi].len() {
                    if m & (1u64 << (idx % BRICK_CELLS)) == 0 {
                        work.q[bi][idx] = 0.0;
                    } else {
                        work.q[bi][idx] =
                            work.mass[bi][idx] * work.p[bi][idx] - wqnew * work.q[bi][idx];
                    }
                }
            }
            let pq = dot(&work.p, &work.q);
            if pq.is_nan() || pq <= 0.0 {
                break; // SPD breakdown ⇒ the acceptance below decides
            }
            let alpha = rz / pq;
            for bi in 0..nb {
                for idx in 0..work.delta[bi].len() {
                    work.delta[bi][idx] += alpha * work.p[bi][idx];
                    work.r[bi][idx] -= alpha * work.q[bi][idx];
                }
            }
            for bi in 0..nb {
                for idx in 0..work.z[bi].len() {
                    work.z[bi][idx] = work.r[bi][idx] / work.diag[bi][idx];
                }
            }
            let rz_new = dot(&work.r, &work.z);
            let beta = rz_new / rz;
            rz = rz_new;
            for bi in 0..nb {
                for idx in 0..work.p[bi].len() {
                    work.p[bi][idx] = work.z[bi][idx] + beta * work.p[bi][idx];
                }
            }
            r_norm2 = dot(&work.r, &work.r);
            iters += 1;
        }
        let resid = (r_norm2 / b_norm2.max(f64::MIN_POSITIVE)).sqrt();
        if b_norm2.is_nan() || (b_norm2 > 0.0 && (resid.is_nan() || resid > EPS_CG_RESID)) {
            return Err(GasDiffError::CgUnconverged {
                comp: comp.name(),
                resid,
            });
        }
        for (bi, xb) in x.iter_mut().enumerate() {
            let m = masked(bi);
            for (idx, xv) in xb.iter_mut().enumerate() {
                if m & (1u64 << (idx % BRICK_CELLS)) != 0 {
                    *xv += work.delta[bi][idx];
                }
            }
        }
        Ok((iters, resid))
    }

    /// Fill `work.mass` with one component's per-cell mass:
    /// u_r/u_z: ρκV; ω: ρr̄²κV (angular-momentum form); T: ρc_vκV; C: ρκV.
    pub(crate) fn fill_mass(
        &self,
        g: &Grid,
        comp: GasComp,
        rho: &BufF,
        tr: &GasTransportField,
        work: &mut GasWork,
    ) {
        for (bi, mb) in work.mass.iter_mut().enumerate() {
            let b = g.brick(bi);
            let nt = b.n_theta();
            for (idx, mv) in mb.iter_mut().enumerate() {
                let local = idx % BRICK_CELLS;
                if b.mask() & (1u64 << local) == 0 {
                    *mv = 1.0;
                    continue;
                }
                let (i_r, _) = b.global_rz(local);
                let kv = b.kappa_rz(local) * g.cell_volume(i_r, nt);
                let rho_c = rho[bi][idx];
                *mv = match comp {
                    GasComp::Ur | GasComp::Uz | GasComp::C => rho_c * kv,
                    GasComp::Om => {
                        let r = g.r_center(i_r);
                        rho_c * r * r * kv
                    }
                    GasComp::T => rho_c * tr.cv[bi][idx] * kv,
                };
            }
        }
    }
}

/// S13b GPU-residency cross-check support (doc-hidden; **additive** — it runs
/// the real `cg_solve` and changes no production number, the same pattern as
/// S13's `EulerWorkspace::rates()`). One in-crate source of truth so the
/// `crates/gpu` harness constructs no module `BufF`: it returns the CPU oracle
/// solution *and* every dense operand the resident GPU CG
/// (`crates/gpu/cuda/residency_diffusion.cu`) needs. N_θ = 1 only (the class-D
/// residency scope; θ/mixed-N_θ is S13c). Dense layout `c = i_r*n_z + i_z`.
#[doc(hidden)]
pub struct GasCgDense {
    /// CPU oracle: `x0 + δ` over gas cells (masked cells hold `x0`).
    pub x_final: Vec<f64>,
    pub iters: usize,
    pub resid: f64,
    pub rho: Vec<f64>,
    pub cv: Vec<f64>,
    pub mu: Vec<f64>,
    pub k: Vec<f64>,
    pub rhod: Vec<f64>,
    /// 1.0 on gas cells, 0.0 elsewhere.
    pub gas: Vec<f64>,
    pub b: Vec<f64>,
    pub x0: Vec<f64>,
    pub n_r: usize,
    pub n_z: usize,
}

impl GasDiffusion<'_> {
    /// Build one component's fixed-structure CG problem from closures, solve it
    /// on the CPU (the bit-exact oracle), and return the solution + the dense
    /// operands for the device cross-check. `comp_id`: 0=Ur 1=Uz 2=Om 3=T 4=C.
    #[doc(hidden)]
    #[allow(clippy::too_many_arguments)]
    pub fn xcheck_cg_dense(
        &self,
        g: &Grid,
        comp_id: usize,
        tr: &GasTransportField,
        rho_of: impl Fn(usize, usize) -> f64,
        b_of: impl Fn(usize, usize) -> f64,
        x0_of: impl Fn(usize, usize) -> f64,
        wq: f64,
    ) -> GasCgDense {
        let comp = match comp_id {
            0 => GasComp::Ur,
            1 => GasComp::Uz,
            2 => GasComp::Om,
            3 => GasComp::T,
            _ => GasComp::C,
        };
        let mut rho = alloc_plane(g);
        let mut work = GasWork::alloc(g);
        let mut x = alloc_plane(g);
        g.for_each_active_cell(|cell| {
            rho[cell.bi][cell.idx] = rho_of(cell.i_r, cell.i_z);
            work.b[cell.bi][cell.idx] = b_of(cell.i_r, cell.i_z);
            x[cell.bi][cell.idx] = x0_of(cell.i_r, cell.i_z);
        });
        let x0_buf = x.clone();
        self.fill_mass(g, comp, &rho, tr, &mut work);
        let (iters, resid) = self
            .cg_solve(g, comp, tr, wq, &mut x, &mut work)
            .expect("xcheck cg_solve converged");
        let (n_r, n_z) = (g.spec().n_r, g.spec().n_z);
        let ncell = n_r * n_z;
        let mut out = GasCgDense {
            x_final: vec![0.0; ncell],
            iters,
            resid,
            rho: vec![0.0; ncell],
            cv: vec![0.0; ncell],
            mu: vec![0.0; ncell],
            k: vec![0.0; ncell],
            rhod: vec![0.0; ncell],
            gas: vec![0.0; ncell],
            b: vec![0.0; ncell],
            x0: vec![0.0; ncell],
            n_r,
            n_z,
        };
        g.for_each_active_cell(|cell| {
            let c = cell.i_r * n_z + cell.i_z;
            out.x_final[c] = x[cell.bi][cell.idx];
            out.rho[c] = rho[cell.bi][cell.idx];
            out.cv[c] = tr.cv[cell.bi][cell.idx];
            out.mu[c] = tr.mu[cell.bi][cell.idx];
            out.k[c] = tr.k[cell.bi][cell.idx];
            out.rhod[c] = tr.rho_d[cell.bi][cell.idx];
            out.gas[c] = 1.0;
            out.b[c] = work.b[cell.bi][cell.idx];
            out.x0[c] = x0_buf[cell.bi][cell.idx];
        });
        out
    }

    /// S13c GPU cross-check for the affine class-D FORCING (doc-hidden,
    /// **additive** — runs the real `fill_lag_gradients` + `assemble_rates`,
    /// changes no production number). Builds the operator's `sol` (current
    /// iterate) and `lag` (cross-term lag state) from closures, plus a per-cell
    /// transport field, and returns the dense 5-component diffusion rate the
    /// device `gpu_class_d_assemble` must reproduce, alongside every dense
    /// input. N_θ = 1 only (the S13c box scope). `c = i_r*n_z + i_z`.
    #[doc(hidden)]
    #[allow(clippy::too_many_arguments)]
    pub fn xcheck_assemble_dense(
        &self,
        g: &Grid,
        tr: &GasTransportField,
        sol_ur: impl Fn(usize, usize) -> f64,
        sol_om: impl Fn(usize, usize) -> f64,
        sol_uz: impl Fn(usize, usize) -> f64,
        sol_tt: impl Fn(usize, usize) -> f64,
        sol_cc: impl Fn(usize, usize) -> f64,
        lag_ur: impl Fn(usize, usize) -> f64,
        lag_uz: impl Fn(usize, usize) -> f64,
    ) -> GasAssembleDense {
        let mut sol = GasOperands::alloc(g);
        let mut lag = GasOperands::alloc(g);
        g.for_each_active_cell(|cell| {
            let (i_r, i_z) = (cell.i_r, cell.i_z);
            sol.ur[cell.bi][cell.idx] = sol_ur(i_r, i_z);
            sol.om[cell.bi][cell.idx] = sol_om(i_r, i_z);
            sol.uz[cell.bi][cell.idx] = sol_uz(i_r, i_z);
            sol.tt[cell.bi][cell.idx] = sol_tt(i_r, i_z);
            sol.cc[cell.bi][cell.idx] = sol_cc(i_r, i_z);
            lag.ur[cell.bi][cell.idx] = lag_ur(i_r, i_z);
            lag.uz[cell.bi][cell.idx] = lag_uz(i_r, i_z);
        });
        let mut work = GasWork::alloc(g);
        let mut rates = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
        self.fill_lag_gradients(g, &lag, &mut work)
            .expect("xcheck fill_lag_gradients");
        self.assemble_rates(g, &sol, &lag, &work, tr, 0.0, &mut rates, None)
            .expect("xcheck assemble_rates");
        let (n_r, n_z) = (g.spec().n_r, g.spec().n_z);
        let ncell = n_r * n_z;
        let mut out = GasAssembleDense {
            rate: vec![[0.0f64; NCOMP]; ncell],
            ur: vec![0.0; ncell],
            om: vec![0.0; ncell],
            uz: vec![0.0; ncell],
            tt: vec![0.0; ncell],
            cc: vec![0.0; ncell],
            lag_ur: vec![0.0; ncell],
            lag_uz: vec![0.0; ncell],
            mu: vec![0.0; ncell],
            k: vec![0.0; ncell],
            rhod: vec![0.0; ncell],
            dhdz: vec![0.0; ncell],
            gas: vec![0.0; ncell],
            n_r,
            n_z,
        };
        g.for_each_active_cell(|cell| {
            let c = cell.i_r * n_z + cell.i_z;
            out.rate[c] = rates[cell.bi][cell.idx];
            out.ur[c] = sol.ur[cell.bi][cell.idx];
            out.om[c] = sol.om[cell.bi][cell.idx];
            out.uz[c] = sol.uz[cell.bi][cell.idx];
            out.tt[c] = sol.tt[cell.bi][cell.idx];
            out.cc[c] = sol.cc[cell.bi][cell.idx];
            out.lag_ur[c] = lag.ur[cell.bi][cell.idx];
            out.lag_uz[c] = lag.uz[cell.bi][cell.idx];
            out.mu[c] = tr.mu[cell.bi][cell.idx];
            out.k[c] = tr.k[cell.bi][cell.idx];
            out.rhod[c] = tr.rho_d[cell.bi][cell.idx];
            out.dhdz[c] = tr.dh_dz[cell.bi][cell.idx];
            out.gas[c] = 1.0;
        });
        out
    }
}

/// Dense operands + oracle rate for the S13c affine-forcing cross-check
/// (doc-hidden). Dense layout `c = i_r*n_z + i_z`.
#[doc(hidden)]
pub struct GasAssembleDense {
    /// CPU oracle: the 5-component (+ zero mass/burn slots) diffusion rate.
    pub rate: Vec<Cons>,
    pub ur: Vec<f64>,
    pub om: Vec<f64>,
    pub uz: Vec<f64>,
    pub tt: Vec<f64>,
    pub cc: Vec<f64>,
    pub lag_ur: Vec<f64>,
    pub lag_uz: Vec<f64>,
    pub mu: Vec<f64>,
    pub k: Vec<f64>,
    pub rhod: Vec<f64>,
    pub dhdz: Vec<f64>,
    pub gas: Vec<f64>,
    pub n_r: usize,
    pub n_z: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crucible_grid::{GridSpec, Region};
    use crucible_units::{dynamic_viscosity_pa_s, specific_heat_capacity_j_per_kg_k};

    const MU: f64 = 0.01;
    const CP: f64 = 1000.0;
    const PR: f64 = 0.5;
    const GAMMA: f64 = 1.4;
    const SC: f64 = 0.8;

    fn op() -> GasDiffusion<'static> {
        GasDiffusion::new(GasDiffBcs {
            r_inner: FaceGasBc::free(),
            r_outer: FaceGasBc::free(),
            z_lo: FaceGasBc::free(),
            z_hi: FaceGasBc::free(),
        })
    }

    /// The declared-constant spine occupant these fixtures run on — the
    /// same one a `[mechanisms.transport] type = "transport_constant"`
    /// block builds, so the unit tests exercise the production path.
    fn props() -> TransportProps {
        crate::transport::ConstantTransport::new(
            specific_heat_capacity_j_per_kg_k(CP),
            dynamic_viscosity_pa_s(MU),
            PR,
            GAMMA,
            SC,
        )
        .expect("transport")
        .into_props()
    }

    /// A uniform transport field over every gas cell.
    fn tr(g: &Grid) -> GasTransportField {
        tr_with(g, &props())
    }

    fn tr_with(g: &Grid, t: &TransportProps) -> GasTransportField {
        let mut f = GasTransportField::alloc(g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for local in 0..BRICK_CELLS {
                if b.mask() & (1u64 << local) == 0 {
                    f.set_masked(bi, local);
                } else {
                    let (i_r, i_z) = b.global_rz(local);
                    f.set(bi, local, i_r, i_z, t).expect("valid transport");
                }
            }
        }
        f
    }

    fn spec() -> GridSpec {
        GridSpec {
            r_min: 0.5,
            dr: 0.1,
            n_r: 8,
            z_min: 0.0,
            dz: 0.1,
            n_z: 8,
            n_theta_max: 1,
            axisymmetry_assertion: true,
        }
    }

    fn assemble(
        g: &Grid,
        op: &GasDiffusion<'_>,
        sol: &GasOperands,
    ) -> (Vec<Vec<Cons>>, FlowLedger) {
        let mut work = GasWork::alloc(g);
        let mut rates = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
        let mut ledger = FlowLedger::default();
        let t = tr(g);
        op.fill_lag_gradients(g, sol, &mut work).expect("gradients");
        op.assemble_rates(g, sol, sol, &work, &t, 0.0, &mut rates, Some(&mut ledger))
            .expect("assembly");
        (rates, ledger)
    }

    /// One transport owner: the fixtures' coefficients come from the
    /// spine's declared-constant occupant, and the operator restates
    /// nothing (it holds no transport at all since S4).
    #[test]
    fn transport_derivation_by_hand() {
        let t = props();
        assert_eq!(t.mu, MU);
        assert_eq!(t.k, MU * CP / PR);
        assert_eq!(t.cv, CP / GAMMA);
        assert_eq!(t.rho_d, MU / SC);
        // A single-composition gas carries no species-enthalpy flux.
        assert_eq!(t.dh_dz, 0.0);
    }

    /// **The S4 blindness test.** Every S3 fixture runs a UNIFORM transport
    /// field, so a face coefficient that used only the visiting cell's
    /// value — instead of the two-cell average — would be invisible to all
    /// of them, and would also break the CG's symmetry silently (the
    /// coefficient seen from A would differ from the one seen from B). Give
    /// the spine a spatially varying reading and check every core against
    /// the face average, by hand.
    ///
    /// Mutation-proven: replacing `FaceTr::between` with `FaceTr::at` in
    /// `assemble_rates` leaves the whole rest of the battery green and
    /// fails only here.
    #[test]
    fn face_coefficients_are_the_two_cell_average_of_a_varying_spine() {
        let g = Grid::build(spec(), &["dummy"]).expect("grid");
        let o = op();
        // A transport field that varies in BOTH directions and is nowhere
        // symmetric about the probe cell — a one-sided read must show up.
        let base = props();
        let scale = |i_r: usize, i_z: usize| 1.0 + 0.31 * (i_r as f64) + 0.17 * (i_z as f64);
        let mut tf = GasTransportField::alloc(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                if b.mask() & (1u64 << l) == 0 {
                    tf.set_masked(bi, l);
                    continue;
                }
                let (i_r, i_z) = b.global_rz(l);
                let f = scale(i_r, i_z);
                tf.set(
                    bi,
                    l,
                    i_r,
                    i_z,
                    &TransportProps {
                        mu: base.mu * f,
                        k: base.k * f,
                        cv: base.cv,
                        rho_d: base.rho_d * f,
                        dh_dz: 0.0,
                        cp: base.cp,
                        cp_film: base.cp_film,
                        pr: base.pr,
                    },
                )
                .expect("valid");
            }
        }

        let (i_r, i_z) = (4usize, 4usize);
        let local = (i_r % BRICK) * BRICK + (i_z % BRICK);
        let vol = g.cell_volume(i_r, 1);
        let a_z = g.face_area_z(i_r, 1);
        let dz = g.spec().dz;

        // T(z) with a varying k: the two z-faces now carry DIFFERENT
        // coefficients, each the average of its two cells.
        let t_at = |iz: usize| 300.0 + 5.0 * (iz as f64) * (iz as f64);
        let mut sol = GasOperands::alloc(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (_, iz) = b.global_rz(l);
                sol.rho[bi][l] = 1.0;
                sol.tt[bi][l] = t_at(iz);
                sol.cc[bi][l] = 0.5;
            }
        }
        let mut work = GasWork::alloc(&g);
        let mut rates = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
        o.fill_lag_gradients(&g, &sol, &mut work)
            .expect("gradients");
        o.assemble_rates(&g, &sol, &sol, &work, &tf, 0.0, &mut rates, None)
            .expect("assembly");

        let k_hi = base.k * 0.5 * (scale(i_r, i_z) + scale(i_r, i_z + 1));
        let k_lo = base.k * 0.5 * (scale(i_r, i_z) + scale(i_r, i_z - 1));
        let expect = a_z
            * (k_hi * (t_at(i_z + 1) - t_at(i_z)) - k_lo * (t_at(i_z) - t_at(i_z - 1)))
            / (dz * vol);
        let got = rates[0][local][I_EN];
        assert!(
            (got - expect).abs() <= 1e-12 * expect.abs(),
            "varying-k T core: got {got}, expect {expect}"
        );

        // The CG operator must see the SAME face coefficients, or it solves
        // a different matrix than the composition applies (the S3 review's
        // Jacobian identity, now with variable coefficients). Check the
        // implicit apply against the same hand-built expression.
        let mut x = vec![vec![0.0f64; BRICK_CELLS]; g.n_bricks()];
        for (bi, xb) in x.iter_mut().enumerate() {
            let b = g.brick(bi);
            for (l, xv) in xb.iter_mut().enumerate() {
                let (_, iz) = b.global_rz(l);
                *xv = t_at(iz);
            }
        }
        let mut out = vec![vec![0.0f64; BRICK_CELLS]; g.n_bricks()];
        o.apply_linear(&g, GasComp::T, &tf, &x, &mut out, None);
        let expect_lin =
            a_z * (k_hi * (t_at(i_z + 1) - t_at(i_z)) - k_lo * (t_at(i_z) - t_at(i_z - 1))) / dz;
        assert!(
            (out[0][local] - expect_lin).abs() <= 1e-9 * expect_lin.abs(),
            "apply_linear must carry the same face coefficients as the assembly: \
             got {}, expect {expect_lin}",
            out[0][local]
        );

        // Symmetry: the flux across one face, seen from either side, must
        // be the exact negation — otherwise the CG matrix is not symmetric
        // and conservation stops telescoping.
        let mut probe = GasOperands::alloc(&g);
        for bi in 0..g.n_bricks() {
            for l in 0..BRICK_CELLS {
                probe.rho[bi][l] = 1.0;
                probe.tt[bi][l] = 300.0;
                probe.cc[bi][l] = 0.5;
            }
        }
        probe.tt[0][local] = 700.0;
        let mut rates2 = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
        o.fill_lag_gradients(&g, &probe, &mut work)
            .expect("gradients");
        o.assemble_rates(&g, &probe, &probe, &work, &tf, 0.0, &mut rates2, None)
            .expect("assembly");
        let nb_local = (i_r % BRICK) * BRICK + ((i_z + 1) % BRICK);
        let out_flux = rates2[0][local][I_EN] * vol; // κ = 1 in a full box
        let in_flux = rates2[0][nb_local][I_EN] * vol;
        assert!(
            out_flux < 0.0 && in_flux > 0.0,
            "the hot cell must lose what its neighbours gain"
        );
    }

    /// **S4 review finding 4.1.** Deleting the BOUNDARY arm's
    /// species-enthalpy port left the whole battery green — the interior
    /// arm has its own test, the boundary one had none, and the omission is
    /// conservation-neutral (it transports the wrong physics without
    /// breaking the ledger). A `Prescribed` species wall with a non-zero
    /// `∂h/∂Z` must carry enthalpy through that face, by hand.
    #[test]
    fn the_boundary_species_port_carries_its_enthalpy_too() {
        let g = Grid::build(spec(), &["dummy"]).expect("grid");
        const C_WALL: f64 = 0.9;
        const DH_DZ: f64 = 6.0e7;
        let wall_c = |_: f64, _: f64, _: f64, _: f64| C_WALL;
        let o = GasDiffusion::new(GasDiffBcs {
            r_inner: FaceGasBc::free(),
            r_outer: FaceGasBc::free(),
            z_lo: FaceGasBc {
                velocity: VelocityBc::FreeSlip,
                thermal: ThermalBc::Adiabatic,
                species: SpeciesBc::Prescribed(&wall_c),
            },
            z_hi: FaceGasBc::free(),
        });
        let mut sol = GasOperands::alloc(&g);
        for bi in 0..g.n_bricks() {
            for l in 0..BRICK_CELLS {
                sol.rho[bi][l] = 1.0;
                sol.tt[bi][l] = 300.0; // uniform: no Fourier flux anywhere
                sol.cc[bi][l] = 0.1;
            }
        }
        let base = props();
        let tf = tr_with(
            &g,
            &TransportProps {
                dh_dz: DH_DZ,
                ..base
            },
        );
        let mut work = GasWork::alloc(&g);
        let mut rates = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
        let mut ledger = FlowLedger::default();
        o.fill_lag_gradients(&g, &sol, &mut work)
            .expect("gradients");
        o.assemble_rates(
            &g,
            &sol,
            &sol,
            &work,
            &tf,
            0.0,
            &mut rates,
            Some(&mut ledger),
        )
        .expect("assembly");

        // The z_lo boundary cell of the first ring: one half-cell two-point
        // species flux against the wall value, and the enthalpy it carries.
        let (i_r, i_z) = (0usize, 0usize);
        let local = (i_r % BRICK) * BRICK + (i_z % BRICK);
        let vol = g.cell_volume(i_r, 1);
        let a_z = g.face_area_z(i_r, 1);
        let half = 0.5 * g.spec().dz;
        let g_cc = (C_WALL - 0.1) / half;
        let expect_c = a_z * base.rho_d * g_cc / vol;
        let expect_e = a_z * base.rho_d * DH_DZ * g_cc / vol;
        assert!(
            (rates[0][local][I_RC] - expect_c).abs() <= 1e-12 * expect_c.abs(),
            "boundary species flux: got {}, expect {expect_c}",
            rates[0][local][I_RC]
        );
        assert!(
            (rates[0][local][I_EN] - expect_e).abs() <= 1e-12 * expect_e.abs(),
            "boundary species-ENTHALPY port: got {}, expect {expect_e}",
            rates[0][local][I_EN]
        );
        // And it must be ledgered as a port, or the audit's port/source
        // books stop matching the applied increment.
        assert!(
            ledger.port_net[I_EN] > 0.0,
            "the enthalpy entering through a species wall is a PORT"
        );
    }

    /// **S4 review finding 4.1, second mutation.** A one-sided read of
    /// `dh_dz` at an interior face breaks conservation, and nothing in the
    /// battery could see it. Assert the ledger identity directly on a
    /// STRONGLY varying spine: interior fluxes must telescope, so the sum
    /// of applied increments equals ports + sources exactly.
    #[test]
    fn a_varying_spine_still_telescopes_exactly() {
        let g = Grid::build(spec(), &["dummy"]).expect("grid");
        let o = op();
        let base = props();
        let mut tf = GasTransportField::alloc(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                if b.mask() & (1u64 << l) == 0 {
                    tf.set_masked(bi, l);
                    continue;
                }
                let (i_r, i_z) = b.global_rz(l);
                // Independent, strongly-varying fields — a shared scale
                // factor would let a one-sided read of one coefficient hide
                // behind the correct averaging of another.
                let f = 1.0 + 0.7 * i_r as f64;
                let q = 1.0 + 0.5 * i_z as f64;
                tf.set(
                    bi,
                    l,
                    i_r,
                    i_z,
                    &TransportProps {
                        mu: base.mu * f,
                        k: base.k * q,
                        cv: base.cv,
                        rho_d: base.rho_d * (1.0 + 0.3 * i_z as f64),
                        dh_dz: 1.0e7 * (1.0 + 0.9 * i_r as f64),
                        cp: base.cp,
                        cp_film: base.cp_film,
                        pr: base.pr,
                    },
                )
                .expect("valid");
            }
        }
        // A state with gradients in every solved variable.
        let mut sol = GasOperands::alloc(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (i_r, i_z) = b.global_rz(l);
                let (x, y) = (i_r as f64, i_z as f64);
                sol.rho[bi][l] = 1.0 + 0.05 * x;
                sol.ur[bi][l] = 3.0 * (0.1 * x).sin() + 0.7 * y;
                sol.om[bi][l] = 0.4 + 0.03 * x * y;
                sol.uz[bi][l] = 2.0 + 0.5 * x - 0.3 * y;
                sol.tt[bi][l] = 300.0 + 11.0 * x + 7.0 * y;
                sol.cc[bi][l] = 0.2 + 0.01 * x + 0.02 * y;
            }
        }
        let mut work = GasWork::alloc(&g);
        let mut rates = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
        let mut ledger = FlowLedger::default();
        o.fill_lag_gradients(&g, &sol, &mut work)
            .expect("gradients");
        o.assemble_rates(
            &g,
            &sol,
            &sol,
            &work,
            &tf,
            0.0,
            &mut rates,
            Some(&mut ledger),
        )
        .expect("assembly");

        // Σ_cells rate[k]·κV == port_net[k] + src_net[k], to round-off.
        for k in [I_MR, I_MZ, I_EN, I_RC] {
            let mut applied = 0.0f64;
            let mut scale = 0.0f64;
            for (bi, rb) in rates.iter().enumerate() {
                let b = g.brick(bi);
                for (l, cell) in rb.iter().enumerate() {
                    if b.mask() & (1u64 << l) == 0 {
                        continue;
                    }
                    let (i_r, _) = b.global_rz(l);
                    let kv = b.kappa_rz(l) * g.cell_volume(i_r, 1);
                    applied += cell[k] * kv;
                    scale += (cell[k] * kv).abs();
                }
            }
            let books = ledger.port_net[k] + ledger.src_net[k];
            assert!(
                (applied - books).abs() <= 1e-11 * scale.max(f64::MIN_POSITIVE),
                "component {k}: applied {applied} vs ports+sources {books} \
                 (scale {scale}) — interior fluxes must telescope on a \
                 varying spine too"
            );
        }
    }

    /// The temperature solve's mass must be `ρ·c_v·κV` with the **cell's
    /// own** spine `c_v` (module doc). Every S3 fixture had a single c_v,
    /// so a stale constant here would be invisible to all of them.
    #[test]
    fn temperature_solve_mass_uses_the_per_cell_spine_slope() {
        let g = Grid::build(spec(), &["dummy"]).expect("grid");
        let o = op();
        let base = props();
        let cv_at = |i_r: usize, i_z: usize| base.cv * (1.0 + 0.4 * i_r as f64 + 0.9 * i_z as f64);
        let mut tf = GasTransportField::alloc(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                if b.mask() & (1u64 << l) == 0 {
                    tf.set_masked(bi, l);
                    continue;
                }
                let (i_r, i_z) = b.global_rz(l);
                tf.set(
                    bi,
                    l,
                    i_r,
                    i_z,
                    &TransportProps {
                        cv: cv_at(i_r, i_z),
                        ..base
                    },
                )
                .expect("valid");
            }
        }
        let mut rho = vec![vec![0.0f64; BRICK_CELLS]; g.n_bricks()];
        for (bi, rb) in rho.iter_mut().enumerate() {
            let b = g.brick(bi);
            for (l, v) in rb.iter_mut().enumerate() {
                let (i_r, _) = b.global_rz(l);
                *v = 0.7 + 0.05 * i_r as f64;
            }
        }
        let mut work = GasWork::alloc(&g);
        o.fill_mass(&g, GasComp::T, &rho, &tf, &mut work);
        for (bi, rb) in rho.iter().enumerate() {
            let b = g.brick(bi);
            for (l, rv) in rb.iter().enumerate() {
                if b.mask() & (1u64 << l) == 0 {
                    continue;
                }
                let (i_r, i_z) = b.global_rz(l);
                let kv = b.kappa_rz(l) * g.cell_volume(i_r, 1);
                let expect = rv * cv_at(i_r, i_z) * kv;
                assert!(
                    (work.mass[bi][l] - expect).abs() <= 1e-12 * expect,
                    "T mass at ({i_r}, {i_z}): got {}, expect {expect}",
                    work.mass[bi][l]
                );
            }
        }
        // The velocity/species components take ρκV — no heat capacity.
        o.fill_mass(&g, GasComp::Ur, &rho, &tf, &mut work);
        let b = g.brick(0);
        let (i_r, _) = b.global_rz(9);
        assert_eq!(
            work.mass[0][9],
            rho[0][9] * b.kappa_rz(9) * g.cell_volume(i_r, 1)
        );
    }

    /// The species-enthalpy diffusion flux `Σ h_k j_k` (S4). It is
    /// identically zero on a single-composition gas — which is why S3 could
    /// defer it honestly — so nothing in the S3 battery can see it. With a
    /// non-zero `∂h/∂Z` and a composition gradient, energy must move by
    /// exactly `ρD·(∂h/∂Z)·∇C` per face, on top of the Fourier flux.
    #[test]
    fn species_enthalpy_flux_carries_energy_down_a_composition_gradient() {
        let g = Grid::build(spec(), &["dummy"]).expect("grid");
        let o = op();
        let (i_r, i_z) = (4usize, 4usize);
        let local = (i_r % BRICK) * BRICK + (i_z % BRICK);
        let vol = g.cell_volume(i_r, 1);
        let a_z = g.face_area_z(i_r, 1);
        let dz = g.spec().dz;
        let c_at = |iz: usize| 0.1 + 0.02 * (iz as f64) * (iz as f64);

        let mut sol = GasOperands::alloc(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (_, iz) = b.global_rz(l);
                sol.rho[bi][l] = 1.0;
                sol.tt[bi][l] = 300.0; // uniform: no Fourier flux at all
                sol.cc[bi][l] = c_at(iz);
            }
        }
        let mut work = GasWork::alloc(&g);
        o.fill_lag_gradients(&g, &sol, &mut work)
            .expect("gradients");

        // With dh_dz = 0 (a single-composition gas) the energy rate is
        // exactly zero — the S3 arithmetic, unchanged.
        let quiet = tr(&g);
        let mut rates = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
        o.assemble_rates(&g, &sol, &sol, &work, &quiet, 0.0, &mut rates, None)
            .expect("assembly");
        assert_eq!(
            rates[0][local][I_EN], 0.0,
            "a single-composition gas carries no species enthalpy"
        );

        // Now switch the spine to a reacting mixture.
        const DH_DZ: f64 = 6.0e7; // J/kg, the LOX/LH2 class
        let base = props();
        let hot = tr_with(
            &g,
            &TransportProps {
                dh_dz: DH_DZ,
                ..base
            },
        );
        let mut rates = vec![vec![[0.0f64; NCOMP]; BRICK_CELLS]; g.n_bricks()];
        o.assemble_rates(&g, &sol, &sol, &work, &hot, 0.0, &mut rates, None)
            .expect("assembly");
        let expect =
            base.rho_d * DH_DZ * a_z * ((c_at(i_z + 1) - c_at(i_z)) - (c_at(i_z) - c_at(i_z - 1)))
                / (dz * vol);
        let got = rates[0][local][I_EN];
        assert!(
            (got - expect).abs() <= 1e-12 * expect.abs(),
            "species-enthalpy flux: got {got}, expect {expect}"
        );
        // The species flux itself is untouched by the enthalpy limb.
        let expect_c =
            base.rho_d * a_z * ((c_at(i_z + 1) - c_at(i_z)) - (c_at(i_z) - c_at(i_z - 1)))
                / (dz * vol);
        assert!((rates[0][local][I_RC] - expect_c).abs() <= 1e-12 * expect_c.abs());
        // Ratio check, stated as the physics: the energy the flux carries
        // is exactly ∂h/∂Z per unit of composition it moves.
        assert!(
            (got / rates[0][local][I_RC] - DH_DZ).abs() <= 1e-9 * DH_DZ,
            "the enthalpy carried per unit diffused composition must be ∂h/∂Z"
        );
    }

    /// Uniform translation + RIGID ROTATION + uniform T/C: every rate is
    /// identically zero — rigid rotation is discretely stress-free by the
    /// ω-form construction (module doc), and a uniform state produces no
    /// diffusion of anything.
    #[test]
    fn rigid_rotation_and_uniform_state_are_stress_free() {
        let g = Grid::build(spec(), &["dummy"]).expect("grid");
        let o = op();
        let mut sol = GasOperands::alloc(&g);
        for bi in 0..g.n_bricks() {
            for local in 0..BRICK_CELLS {
                sol.rho[bi][local] = 1.2;
                sol.ur[bi][local] = 0.0;
                sol.om[bi][local] = 0.4; // rigid: u_θ = 0.4·r
                sol.uz[bi][local] = 0.7; // uniform translation
                sol.tt[bi][local] = 300.0;
                sol.cc[bi][local] = 0.4;
            }
        }
        let (rates, ledger) = assemble(&g, &o, &sol);
        for b in &rates {
            for cell in b {
                for (k, v) in cell.iter().enumerate() {
                    assert_eq!(*v, 0.0, "component {k} rate not exactly zero");
                }
            }
        }
        for k in 0..NCOMP {
            assert_eq!(ledger.port_net[k], 0.0);
            assert_eq!(ledger.src_net[k], 0.0);
        }
    }

    /// The discrete coefficients, checked by hand at one interior cell for
    /// every solved component (the wall_heat `colburn_algebra_by_hand`
    /// pattern): single-variable profiles isolate each two-point core.
    #[test]
    fn face_coefficients_by_hand_at_one_cell() {
        let g = Grid::build(spec(), &["dummy"]).expect("grid");
        let o = op();
        let (i_r, i_z) = (4usize, 4usize);
        let local = (i_r % BRICK) * BRICK + (i_z % BRICK);
        let vol = g.cell_volume(i_r, 1);
        let (a_in, a_out) = (g.face_area_r(i_r, false, 1), g.face_area_r(i_r, true, 1));
        let a_z = g.face_area_z(i_r, 1);
        let (dr, dz) = (g.spec().dr, g.spec().dz);
        let base = |g: &Grid| {
            let mut s = GasOperands::alloc(g);
            for bi in 0..g.n_bricks() {
                for l in 0..BRICK_CELLS {
                    s.rho[bi][l] = 1.0;
                    s.tt[bi][l] = 300.0;
                    s.cc[bi][l] = 0.5;
                }
            }
            s
        };

        // T(z) parabola: rate_E = k·A_z·(T₊ − 2T + T₋)/(dz·κV).
        let mut sol = base(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (_, iz) = b.global_rz(l);
                sol.tt[bi][l] = 300.0 + 5.0 * (iz as f64) * (iz as f64);
            }
        }
        let (rates, _) = assemble(&g, &o, &sol);
        let t_at = |iz: usize| 300.0 + 5.0 * (iz as f64) * (iz as f64);
        let expect = props().k * a_z * ((t_at(i_z + 1) - t_at(i_z)) - (t_at(i_z) - t_at(i_z - 1)))
            / (dz * vol);
        let got = rates[0][local][I_EN];
        assert!(
            (got - expect).abs() <= 1e-12 * expect.abs(),
            "T core: got {got}, expect {expect}"
        );

        // C(z) parabola: rate_C = ρD·A_z·Δ²C/(dz·κV).
        let mut sol = base(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (_, iz) = b.global_rz(l);
                sol.cc[bi][l] = 0.1 + 0.02 * (iz as f64) * (iz as f64);
            }
        }
        let (rates, _) = assemble(&g, &o, &sol);
        let c_at = |iz: usize| 0.1 + 0.02 * (iz as f64) * (iz as f64);
        let expect =
            props().rho_d * a_z * ((c_at(i_z + 1) - c_at(i_z)) - (c_at(i_z) - c_at(i_z - 1)))
                / (dz * vol);
        let got = rates[0][local][I_RC];
        assert!(
            (got - expect).abs() <= 1e-12 * expect.abs(),
            "C core: got {got}, expect {expect}"
        );

        // u_z(r): rate_MZ = μ·(A₊Δu₊ − A₋Δu₋)/(dr·κV) — coefficient μ (the
        // τ_rz core), exact cylindrical areas.
        let mut sol = base(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (ir, _) = b.global_rz(l);
                sol.uz[bi][l] = 3.0 * (ir as f64) * (ir as f64);
            }
        }
        let (rates, _) = assemble(&g, &o, &sol);
        let u_at = |ir: usize| 3.0 * (ir as f64) * (ir as f64);
        let expect = props().mu
            * (a_out * (u_at(i_r + 1) - u_at(i_r)) - a_in * (u_at(i_r) - u_at(i_r - 1)))
            / (dr * vol);
        let got = rates[0][local][I_MZ];
        assert!(
            (got - expect).abs() <= 1e-12 * expect.abs(),
            "u_z core: got {got}, expect {expect}"
        );

        // ω(r): the angular-momentum form — rate_MT = μ·(A₊r_f₊²Δω₊ −
        // A₋r_f₋²Δω₋)/(dr·r̄·κV).
        let mut sol = base(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (ir, _) = b.global_rz(l);
                sol.om[bi][l] = 2.0 + 0.5 * (ir as f64) * (ir as f64);
            }
        }
        let (rates, _) = assemble(&g, &o, &sol);
        let om_at = |ir: usize| 2.0 + 0.5 * (ir as f64) * (ir as f64);
        let (rf_in, rf_out) = (g.face_radius(i_r), g.face_radius(i_r + 1));
        let expect = props().mu
            * (a_out * rf_out * rf_out * (om_at(i_r + 1) - om_at(i_r))
                - a_in * rf_in * rf_in * (om_at(i_r) - om_at(i_r - 1)))
            / (dr * g.r_center(i_r) * vol);
        let got = rates[0][local][I_MT];
        assert!(
            (got - expect).abs() <= 1e-12 * expect.abs(),
            "ω core: got {got}, expect {expect}"
        );

        // u_r(r) linear in r-index: the (4/3)μ radial core + the −τ_θθ/r
        // geometric source at the cell (lag = sol here, so e_rr is the
        // cell-centered central derivative and e_θθ = u_r/r̄).
        let mut sol = base(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (ir, _) = b.global_rz(l);
                sol.ur[bi][l] = 0.2 * (ir as f64);
            }
        }
        let (rates, _) = assemble(&g, &o, &sol);
        let ur_at = |ir: usize| 0.2 * (ir as f64);
        let rbar = g.r_center(i_r);
        let geo = (a_out - a_in) / vol;
        let e_rr = (ur_at(i_r + 1) - ur_at(i_r - 1)) / (2.0 * dr);
        // Radial-face compressible correction −⅔μ(e_θθ+e_zz)|face (e_zz=0):
        let e_thth_out = 0.5 * (ur_at(i_r) / rbar + ur_at(i_r + 1) / g.r_center(i_r + 1));
        let e_thth_in = 0.5 * (ur_at(i_r) / rbar + ur_at(i_r - 1) / g.r_center(i_r - 1));
        let tau_out = (4.0 / 3.0) * props().mu * (ur_at(i_r + 1) - ur_at(i_r)) / dr
            - (2.0 / 3.0) * props().mu * e_thth_out;
        let tau_in = (4.0 / 3.0) * props().mu * (ur_at(i_r) - ur_at(i_r - 1)) / dr
            - (2.0 / 3.0) * props().mu * e_thth_in;
        let tau_thth =
            (4.0 / 3.0) * props().mu * ur_at(i_r) / rbar - (2.0 / 3.0) * props().mu * e_rr;
        let expect = (a_out * tau_out - a_in * tau_in) / vol - tau_thth * geo;
        let got = rates[0][local][I_MR];
        assert!(
            (got - expect).abs() <= 1e-12 * expect.abs().max(1e-30),
            "u_r core: got {got}, expect {expect}"
        );
    }

    /// The compressible (dilatation) and cross-shear terms, on a field
    /// whose FOUR lag gradients `∂u_r/∂r`, `∂u_r/∂z`, `∂u_z/∂r`,
    /// `∂u_z/∂z` all take DISTINCT values — so any confusion among them
    /// is caught. This is a deliberate complement to the coefficient test
    /// above (which drives one velocity component at a time, leaving the
    /// other's gradients identically zero) and to the MMS study (whose
    /// shared mode gives `u_r` and `u_z` identical gradient fields):
    /// under both of those, swapping `duz_dz` for `dur_dz` in the −⅔μ∇·u
    /// corrections is INVISIBLE. Verified by planting exactly that
    /// mutation: the whole battery stayed green, this test fails loudly
    /// (S3 review finding).
    #[test]
    fn dilatation_and_cross_shear_discriminate_independent_gradients() {
        let g = Grid::build(spec(), &["dummy"]).expect("grid");
        let o = op();
        let (i_r, i_z) = (4usize, 4usize);
        let local = (i_r % BRICK) * BRICK + (i_z % BRICK);
        let (dr, dz) = (g.spec().dr, g.spec().dz);
        let vol = g.cell_volume(i_r, 1);
        let (a_in, a_out) = (g.face_area_r(i_r, false, 1), g.face_area_r(i_r, true, 1));
        let a_z = g.face_area_z(i_r, 1);
        let rbar = g.r_center(i_r);

        // Independent shapes with BILINEAR cross terms. The cross terms are
        // essential, not decoration: a dilatation error that is spatially
        // UNIFORM cancels exactly out of the r-momentum — the face term
        // carries it with weight (A_out−A_in)/V and the −τ_θθ/r source
        // with `geo`, which are the same number (a uniform isotropic
        // stress exerts no net force — real physics, and the reason a
        // simpler field is blind here). A3/B3 make every gradient vary
        // across the stencil so nothing cancels.
        const A1: f64 = 0.7;
        const A2: f64 = 0.3;
        const A3: f64 = 0.23;
        const B1: f64 = 0.11;
        const B2: f64 = 1.3;
        const B3: f64 = 0.37;
        let ur_at = |ir: usize, iz: usize| A1 * ir as f64 + A2 * iz as f64 + A3 * (ir * iz) as f64;
        let uz_at =
            |ir: usize, iz: usize| B1 * (iz * iz) as f64 + B2 * ir as f64 + B3 * (ir * iz) as f64;
        let mut sol = GasOperands::alloc(&g);
        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            for l in 0..BRICK_CELLS {
                let (ir, iz) = b.global_rz(l);
                sol.rho[bi][l] = 1.0;
                sol.tt[bi][l] = 300.0;
                sol.cc[bi][l] = 0.5;
                sol.ur[bi][l] = ur_at(ir, iz);
                sol.uz[bi][l] = uz_at(ir, iz);
            }
        }
        let (rates, _) = assemble(&g, &o, &sol);

        // Cell-centered lag gradients, from the analytic central-difference
        // formulas (exact for these polynomials) — derived here, never read
        // back from the operator.
        let f = |i: usize| i as f64;
        let dur_dr = |_ir: usize, iz: usize| (A1 + A3 * f(iz)) / dr;
        let dur_dz = |ir: usize, _iz: usize| (A2 + A3 * f(ir)) / dz;
        let duz_dr = |_ir: usize, iz: usize| (B2 + B3 * f(iz)) / dr;
        let duz_dz = |ir: usize, iz: usize| (2.0 * B1 * f(iz) + B3 * f(ir)) / dz;
        for (a, b) in [
            (dur_dr(i_r, i_z), dur_dz(i_r, i_z)),
            (dur_dr(i_r, i_z), duz_dr(i_r, i_z)),
            (dur_dr(i_r, i_z), duz_dz(i_r, i_z)),
            (dur_dz(i_r, i_z), duz_dr(i_r, i_z)),
            (dur_dz(i_r, i_z), duz_dz(i_r, i_z)),
            (duz_dr(i_r, i_z), duz_dz(i_r, i_z)),
        ] {
            assert!(
                (a - b).abs() > 1e-6 * a.abs().max(b.abs()),
                "the four lag gradients must be pairwise distinct to discriminate"
            );
        }
        let two_thirds = (2.0 / 3.0) * props().mu;
        let four_thirds = (4.0 / 3.0) * props().mu;
        let geo = (a_out - a_in) / vol;
        let e_thth = |ir: usize, iz: usize| ur_at(ir, iz) / g.r_center(ir);
        let avg = |x: f64, y: f64| 0.5 * (x + y);

        // --- r-momentum ------------------------------------------------
        // τ_rr on the two radial faces (each carries the face-averaged
        // e_θθ and e_zz — the u_z → u_r dilatation channel), the τ_rz
        // z-faces (carrying the face-averaged duz_dr — cross shear), and
        // the −τ_θθ/r volume source (carrying e_rr AND e_zz).
        let tau_rr_out = four_thirds * dur_dr(i_r, i_z)
            - two_thirds
                * (avg(e_thth(i_r, i_z), e_thth(i_r + 1, i_z))
                    + avg(duz_dz(i_r, i_z), duz_dz(i_r + 1, i_z)));
        let tau_rr_in = four_thirds * dur_dr(i_r, i_z)
            - two_thirds
                * (avg(e_thth(i_r, i_z), e_thth(i_r - 1, i_z))
                    + avg(duz_dz(i_r, i_z), duz_dz(i_r - 1, i_z)));
        let tau_rz_zp =
            props().mu * (dur_dz(i_r, i_z) + avg(duz_dr(i_r, i_z), duz_dr(i_r, i_z + 1)));
        let tau_rz_zm =
            props().mu * (dur_dz(i_r, i_z) + avg(duz_dr(i_r, i_z), duz_dr(i_r, i_z - 1)));
        let tau_thth =
            four_thirds * e_thth(i_r, i_z) - two_thirds * (dur_dr(i_r, i_z) + duz_dz(i_r, i_z));
        let expect_mr = (a_out * tau_rr_out - a_in * tau_rr_in + a_z * (tau_rz_zp - tau_rz_zm))
            / vol
            - tau_thth * geo;
        let got_mr = rates[0][local][I_MR];
        assert!(
            (got_mr - expect_mr).abs() <= 1e-11 * expect_mr.abs(),
            "r-momentum dilatation/cross-shear: got {got_mr}, expect {expect_mr}"
        );

        // --- z-momentum ------------------------------------------------
        // τ_rz on the radial faces (face-averaged dur_dz — cross shear)
        // and τ_zz on the z-faces (face-averaged e_rr and e_θθ — the
        // u_r → u_z dilatation channel).
        let tau_rz_rp =
            props().mu * (duz_dr(i_r, i_z) + avg(dur_dz(i_r, i_z), dur_dz(i_r + 1, i_z)));
        let tau_rz_rm =
            props().mu * (duz_dr(i_r, i_z) + avg(dur_dz(i_r, i_z), dur_dz(i_r - 1, i_z)));
        let g_uz_zp = (uz_at(i_r, i_z + 1) - uz_at(i_r, i_z)) / dz;
        let g_uz_zm = (uz_at(i_r, i_z) - uz_at(i_r, i_z - 1)) / dz;
        let tau_zz_p = four_thirds * g_uz_zp
            - two_thirds
                * (avg(dur_dr(i_r, i_z), dur_dr(i_r, i_z + 1))
                    + avg(ur_at(i_r, i_z), ur_at(i_r, i_z + 1)) / rbar);
        let tau_zz_m = four_thirds * g_uz_zm
            - two_thirds
                * (avg(dur_dr(i_r, i_z), dur_dr(i_r, i_z - 1))
                    + avg(ur_at(i_r, i_z), ur_at(i_r, i_z - 1)) / rbar);
        let expect_mz = (a_out * tau_rz_rp - a_in * tau_rz_rm + a_z * (tau_zz_p - tau_zz_m)) / vol;
        let got_mz = rates[0][local][I_MZ];
        assert!(
            (got_mz - expect_mz).abs() <= 1e-11 * expect_mz.abs(),
            "z-momentum cross-shear/dilatation: got {got_mz}, expect {expect_mz}"
        );
    }

    /// SOLV-1 §3.5 ownership: the resolved diffusion never reads through a
    /// gas↔solid face — the gas rates are bitwise independent of the solid
    /// cells' operand values (the wall function owns those faces), and the
    /// suppressed world still telescopes to zero net energy.
    #[test]
    fn wall_law_faces_are_suppressed() {
        let g = Grid::build_with_regions(spec(), &["dummy"], |i_r, _| {
            if i_r >= 6 { Region::Solid } else { Region::Gas }
        })
        .expect("region world");
        let o = op();
        let fill = |solid_value: f64| {
            let mut s = GasOperands::alloc(&g);
            for bi in 0..g.n_bricks() {
                let b = g.brick(bi);
                for l in 0..BRICK_CELLS {
                    let (ir, iz) = b.global_rz(l);
                    if b.mask() & (1u64 << l) != 0 {
                        s.rho[bi][l] = 1.0;
                        s.ur[bi][l] = 0.0;
                        s.om[bi][l] = 0.0;
                        s.uz[bi][l] = 0.1 * (ir as f64) + 0.05 * (iz as f64);
                        s.tt[bi][l] = 300.0 + 10.0 * (ir as f64);
                        s.cc[bi][l] = 0.5;
                    } else {
                        // Garbage the wall function would never let near
                        // the gas: MUST be unread.
                        s.rho[bi][l] = solid_value;
                        s.ur[bi][l] = solid_value;
                        s.om[bi][l] = solid_value;
                        s.uz[bi][l] = solid_value;
                        s.tt[bi][l] = solid_value;
                        s.cc[bi][l] = solid_value;
                    }
                }
            }
            s
        };
        let (ra, la) = assemble(&g, &o, &fill(1.0e6));
        let (rb, lb) = assemble(&g, &o, &fill(-5.0e8));
        for (ba, bb) in ra.iter().zip(&rb) {
            for (ca, cb) in ba.iter().zip(bb) {
                for k in 0..NCOMP {
                    assert_eq!(
                        ca[k].to_bits(),
                        cb[k].to_bits(),
                        "gas rate depends on a solid cell's operands — suppression broken"
                    );
                }
            }
        }
        assert_eq!(la.port_net[I_EN].to_bits(), lb.port_net[I_EN].to_bits());
        // All faces are either interior gas↔gas (telescope) or suppressed/
        // free-slip: the net applied energy is zero to round-off.
        let mut net = 0.0f64;
        let mut gross = 0.0f64;
        for (bi, rb) in ra.iter().enumerate() {
            let b = g.brick(bi);
            for (l, cell) in rb.iter().enumerate() {
                if b.mask() & (1u64 << l) != 0 {
                    let (ir, _) = b.global_rz(l);
                    let kv = b.kappa_rz(l) * g.cell_volume(ir, 1);
                    net += kv * cell[I_EN];
                    gross += (kv * cell[I_EN]).abs();
                }
            }
        }
        assert!(
            net.abs() <= 1e-12 * gross.max(1e-30),
            "suppressed world leaks energy: net {net:.3e} vs gross {gross:.3e}"
        );
    }

    /// **The S9 trap gate (COUP-3 0.4.5): uniform transverse flow is not
    /// spuriously damped.** The Cartesian translation u_r = U·cosθ,
    /// u_θ = −U·sinθ has ZERO true stress; the θ-θ core alone would damp
    /// it at the full μU/r² scale (the partial-tensor wrong physics the S8
    /// split refused to ship). With the curvature/cross partners present
    /// the discrete residual must sit at truncation (O(Δθ²) + O(Δr²)),
    /// orders below that scale, and shrink as N_θ refines.
    #[test]
    fn transverse_flow_is_not_spuriously_damped() {
        const U: f64 = 1.0;
        let residual = |nt: u32| -> f64 {
            let spec = GridSpec {
                r_min: 1.0,
                dr: 0.05,
                n_r: 8,
                z_min: 0.0,
                dz: 0.05,
                n_z: 8,
                n_theta_max: nt,
                axisymmetry_assertion: false,
            };
            let g = Grid::build(spec, &["dummy"]).expect("grid");
            let o = op();
            let mut sol = GasOperands::alloc(&g);
            for bi in 0..g.n_bricks() {
                let b = g.brick(bi);
                let n = b.n_theta() as usize;
                for jt in 0..n {
                    let th = Grid::theta_center(jt as u32, b.n_theta());
                    for l in 0..BRICK_CELLS {
                        let idx = jt * BRICK_CELLS + l;
                        let (ir, _) = b.global_rz(l);
                        let r = g.r_center(ir);
                        sol.rho[bi][idx] = 1.2;
                        sol.ur[bi][idx] = U * th.cos();
                        sol.om[bi][idx] = -U * th.sin() / r;
                        sol.uz[bi][idx] = 0.3;
                        sol.tt[bi][idx] = 300.0;
                        sol.cc[bi][idx] = 0.4;
                    }
                }
            }
            let mut work = GasWork::alloc(&g);
            let mut rates: Vec<Vec<Cons>> = (0..g.n_bricks())
                .map(|bi| vec![[0.0f64; NCOMP]; g.brick(bi).n_theta() as usize * BRICK_CELLS])
                .collect();
            let t = {
                let mut f = GasTransportField::alloc(&g);
                let p = props();
                for bi in 0..g.n_bricks() {
                    let b = g.brick(bi);
                    for jt in 0..b.n_theta() as usize {
                        for l in 0..BRICK_CELLS {
                            let (ir, iz) = b.global_rz(l);
                            f.set(bi, jt * BRICK_CELLS + l, ir, iz, &p).expect("tr");
                        }
                    }
                }
                f
            };
            o.fill_lag_gradients(&g, &sol, &mut work).expect("grads");
            o.assemble_rates(&g, &sol, &sol, &work, &t, 0.0, &mut rates, None)
                .expect("assembly");
            let (n_r, n_z) = (g.spec().n_r, g.spec().n_z);
            let mut worst = 0.0f64;
            for (bi, rb) in rates.iter().enumerate() {
                let b = g.brick(bi);
                for (idx, cell) in rb.iter().enumerate() {
                    let (ir, iz) = b.global_rz(idx % BRICK_CELLS);
                    // Interior cells only: at domain edges the declared
                    // zero-flux BC leaves the interior face's O(Δθ²)
                    // truncation uncancelled with a 1/Δ amplification —
                    // the module-doc first-order-local boundary effect,
                    // not the m = 1 damping this gate measures.
                    if ir == 0 || iz == 0 || ir + 1 == n_r || iz + 1 == n_z {
                        continue;
                    }
                    if b.mask() & (1u64 << (idx % BRICK_CELLS)) != 0 {
                        // C is uniform: every species gradient is an exact
                        // zero (no dilatation limb enters species).
                        assert_eq!(cell[I_RC], 0.0, "species θ flux on uniform C");
                        // I_MZ joins the residual: at z-edge cells the
                        // zero-flux boundary face cannot cancel the
                        // interior face's O(Δθ²) dilatation-limb τ_zz —
                        // truncation, not damping.
                        worst = worst
                            .max(cell[I_MR].abs())
                            .max(cell[I_MT].abs())
                            .max(cell[I_MZ].abs())
                            .max(cell[I_EN].abs());
                    }
                }
            }
            worst
        };
        // The partial-tensor damping scale the trap names: ρ·du_r/dt would
        // be ~μU/r² at r_min if the curvature partners were missing.
        let scale = MU * U / (1.0 * 1.0);
        let (r8, r32) = (residual(8), residual(32));
        println!(
            "transverse-flow residual: N_θ=8 {r8:.3e}, N_θ=32 {r32:.3e}, \
             partial-tensor scale {scale:.3e}"
        );
        assert!(
            r32 <= 0.025 * scale,
            "N_θ=32 residual {r32:.3e} not ≪ the partial-tensor damping {scale:.3e}"
        );
        assert!(
            r8 <= 0.25 * scale,
            "N_θ=8 residual {r8:.3e} not below the partial-tensor damping {scale:.3e}"
        );
        assert!(
            r32 < 0.5 * r8,
            "residual does not shrink with θ refinement: {r8:.3e} → {r32:.3e}"
        );
    }

    /// The S13b GPU cross-check accessor (`xcheck_cg_dense`) must run the real
    /// `cg_solve` and CONVERGE for every one of the five components on the same
    /// N_θ=1 box fixture the device harness uses — a CPU-side de-risk of the
    /// (Ben-gated) GPU trip: it proves the accessor + fixture are healthy and
    /// records the iteration counts the box feeds the resident CG. Uniform
    /// transport here (the varying-coefficient face-average path is the
    /// device harness's job); convergence + dense-mapping is what this guards.
    #[test]
    fn xcheck_cg_dense_converges_for_all_components() {
        let (n_r, n_z) = (48usize, 96usize);
        let spec = GridSpec {
            r_min: 0.5,
            dr: 1.0 / n_r as f64,
            n_r,
            z_min: 0.0,
            dz: 1.0 / n_z as f64,
            n_z,
            n_theta_max: 1,
            axisymmetry_assertion: true,
        };
        let g = Grid::build(spec, &["dummy"]).expect("grid");
        let tf = tr(&g);
        let o = op();
        let rho_of =
            |i_r: usize, i_z: usize| 1.0 + 0.2 * (1.7 * i_r as f64 + 0.9 * i_z as f64).sin();
        let b_of = |i_r: usize, i_z: usize| {
            0.7 * (0.30 * i_r as f64).sin() * (0.21 * i_z as f64).cos()
                + 0.15 * (0.9 * i_r as f64 - 0.5 * i_z as f64).sin()
        };
        let x0_of =
            |i_r: usize, i_z: usize| 0.4 + 0.1 * (0.5 * i_r as f64 + 0.3 * i_z as f64).cos();
        let wq = 5.0e-2;
        for comp in 0..5usize {
            let d = o.xcheck_cg_dense(&g, comp, &tf, rho_of, b_of, x0_of, wq);
            assert!(
                d.resid <= EPS_CG_RESID && d.iters >= 1 && d.iters < N_CG_ITERS_MAX,
                "component {comp}: resid {:.3e} iters {} — CG did not converge cleanly",
                d.resid,
                d.iters
            );
            // The dense map covers every cell of the box.
            assert_eq!(d.gas.iter().filter(|&&v| v == 1.0).count(), n_r * n_z);
            println!("comp {comp}: iters {}, resid {:.3e}", d.iters, d.resid);
        }
    }

    /// The S13c full-resident-iterate accessor (`sdc::xcheck_class_d_iterate_dense`)
    /// must run the real assemble → fill_rhs → CG sweep for all five components
    /// and produce a finite, actually-updated `sol` — a CPU-side de-risk of the
    /// (box-only) GPU `gpu_class_d_iterate` cross-check: it proves the CPU
    /// reference + fixture are healthy before the box run.
    #[test]
    fn xcheck_class_d_iterate_runs_and_updates() {
        let (n_r, n_z) = (48usize, 96usize);
        let spec = GridSpec {
            r_min: 0.5,
            dr: 1.0 / n_r as f64,
            n_r,
            z_min: 0.0,
            dz: 1.0 / n_z as f64,
            n_z,
            n_theta_max: 1,
            axisymmetry_assertion: true,
        };
        let g = Grid::build(spec, &["dummy"]).expect("grid");
        let tf = tr(&g);
        let o = op();
        let sol_rho = |r: usize, z: usize| 1.0 + 0.2 * (1.7 * r as f64 + 0.9 * z as f64).sin();
        let sol_ur = |r: usize, z: usize| 0.30 * (1.1 * z as f64 - 0.4 * r as f64).sin();
        let sol_om = |r: usize, z: usize| 0.25 * (0.8 * r as f64 + 0.5 * z as f64).cos();
        let sol_uz = |r: usize, z: usize| 0.35 * (0.9 * r as f64 - 1.3 * z as f64).sin();
        let sol_tt = |r: usize, z: usize| 300.0 + 20.0 * (0.3 * r as f64 + 0.2 * z as f64).sin();
        let sol_cc = |r: usize, _z: usize| 0.5 + 0.1 * (0.25 * r as f64).sin();
        let lag_ur = |r: usize, z: usize| sol_ur(r, z) + 0.05 * (0.7 * r as f64).cos();
        let lag_uz = |r: usize, z: usize| sol_uz(r, z) - 0.04 * (0.5 * r as f64).sin();
        let dlag_of = |_r: usize, _z: usize| [0.0f64; NCOMP];
        let d = crate::sdc::xcheck_class_d_iterate_dense(
            &g, &o, &tf, 5.0e-2, sol_rho, sol_ur, sol_om, sol_uz, sol_tt, sol_cc, lag_ur, lag_uz,
            dlag_of,
        );
        // Every updated component is finite and the solve moved the state.
        let mut moved = false;
        for c in 0..n_r * n_z {
            for (fin, init) in [
                (d.tt[c], d.init_tt[c]),
                (d.cc[c], d.init_cc[c]),
                (d.ur[c], d.init_ur[c]),
                (d.uz[c], d.init_uz[c]),
                (d.om[c], d.init_om[c]),
            ] {
                assert!(fin.is_finite(), "non-finite iterate output at cell {c}");
                if (fin - init).abs() > 1e-12 {
                    moved = true;
                }
            }
        }
        assert!(
            moved,
            "the resident iterate left the state entirely unchanged"
        );
    }
}
