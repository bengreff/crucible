//! SOLV-1 §3.1 — **`F_visc`, the missing forces (plan S3)**: compressible
//! viscous stress + Fourier heat conduction + Fickian species diffusion on
//! the exact cylindrical metric, axisymmetric-with-swirl (N_θ = 1; the θ
//! diffusion fluxes and per-θ re-keying arrive with the 3-D wave, plan S8
//! — refused, never guessed). One flux-form operator over the gas state —
//! no material or regime branch; transport (μ, Pr → k, Sc → ρD, c_p/c_v)
//! is **pure config data** (Rule 13), shared with the wall law (one owner:
//! [`WallLaw`]'s constant set — the degenerate FND-7 spine occupant until
//! the S4 transport tables land).
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
//! - energy (total-energy flux form): `∇·(τ·u + k∇T)` — dissipation is not
//!   a separate term; it emerges from the KE/internal-energy bookkeeping,
//!   which is what conserves total energy by construction
//! - species: `∇·(ρD ∇C)` with constant ρD = μ/Sc (the Fickian
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
//! ## Recorded deferrals (owners named)
//! - **Wall-function skin-friction momentum debit**: the wall law is a
//!   heat law (SOLV-1 §3.5); its tangential-force leg rides the COUP-2
//!   §3.1.2 mount-reaction ledger (the verdict wave). Until then wall-law
//!   faces are momentum-slip, declared.
//! - **Species-enthalpy diffusion flux** `Σ h_k j_k` and a **TableEos-
//!   consistent T refresh** in the Picard: both need the S4 spine
//!   (per-cell c_p/c_v/partial enthalpies). The constant-c_v occupant here
//!   is exact for the gamma-law class; S4 owns the general case.
//! - **COUP-8 registry row + config grammar**: lands with S4's engine
//!   wiring (a manifest without its `from_loaded` path would be half-wired).
//! - Near-wall/boundary lagged-cross stencils are one-sided (first-order
//!   locally — the MMS battery verifies the composed order).
//!
//! Determinism: serial fixed-order assembly (brick/local/face order — the
//! same recorded perf deferral as the solid assembly; the GPU wave
//! parallelizes both), `tree_combine` reductions in the CG, fixed sweep
//! structure. Bit-identical at any thread count.

use crate::euler::{Cons, FlowLedger, I_EN, I_MR, I_MT, I_MZ, I_RC, NCOMP};
use crate::sdc::{EPS_CG_RESID, N_CG_ITERS_MAX};
use crate::wall_heat::WallLaw;
use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, Grid, tree_combine};

pub(crate) type BufF = Vec<Vec<f64>>;

#[derive(Debug, Clone, PartialEq)]
pub enum GasDiffError {
    /// Gas diffusion at N_θ > 1 arrives with the 3-D wave (plan S8):
    /// θ-direction diffusion fluxes + per-θ operand keying. Refuse.
    AzimuthalResolution,
    /// A non-physical transport coefficient (META-1 P6).
    BadCoefficient(&'static str),
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
                "gas diffusion at N_θ > 1 arrives with the 3-D wave (plan S8); \
                 refusing rather than guessing"
            ),
            Self::BadCoefficient(which) => write!(
                f,
                "{which} must be finite and positive (fail loud, META-1 P6)"
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

/// The gas-phase diffusion operator (module doc). Constructed from the one
/// transport owner ([`WallLaw`]) — the constants are never restated.
pub struct GasDiffusion<'a> {
    /// Dynamic viscosity μ [Pa·s].
    pub(crate) mu: f64,
    /// Thermal conductivity k = μ·c_p/Pr [W/(m·K)].
    pub(crate) k_gas: f64,
    /// Specific heats [J/(kg·K)]: c_v = c_p/γ is the implicit T-solve's
    /// linearization slope (exact for the gamma-law class).
    pub(crate) cv: f64,
    /// Species diffusion coefficient ρD = μ/Sc [kg/(m·s)] (constant-ρD
    /// Fickian occupant).
    pub(crate) rho_d: f64,
    pub bcs: GasDiffBcs<'a>,
}

impl<'a> GasDiffusion<'a> {
    /// Build from the one transport owner: k, ρD, c_v are DERIVED from the
    /// wall law's (c_p, μ, Pr) + the declared γ and Schmidt number — no
    /// second statement of any constant (Rule 13, one owner).
    pub fn from_transport(
        law: &WallLaw,
        gamma: f64,
        schmidt: f64,
        bcs: GasDiffBcs<'a>,
    ) -> Result<Self, GasDiffError> {
        if !(gamma.is_finite() && gamma > 1.0) {
            return Err(GasDiffError::BadCoefficient("gamma (> 1)"));
        }
        if !(schmidt.is_finite() && schmidt > 0.0) {
            return Err(GasDiffError::BadCoefficient("schmidt"));
        }
        Ok(GasDiffusion {
            mu: law.mu(),
            k_gas: law.k_gas(),
            cv: law.cp() / gamma,
            rho_d: law.mu() / schmidt,
            bcs,
        })
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
        let mk = || vec![vec![0.0f64; BRICK_CELLS]; g.n_bricks()];
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
pub struct GasWork {
    pub(crate) dur_dr: BufF,
    pub(crate) dur_dz: BufF,
    pub(crate) duz_dr: BufF,
    pub(crate) duz_dz: BufF,
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
        let mk = || vec![vec![0.0f64; BRICK_CELLS]; g.n_bricks()];
        GasWork {
            dur_dr: mk(),
            dur_dz: mk(),
            duz_dr: mk(),
            duz_dz: mk(),
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
    fn validate(&self, g: &Grid) -> Result<(), GasDiffError> {
        if g.brick(0).n_theta() != 1 || g.bricks().iter().any(|b| b.n_theta() != 1) {
            return Err(GasDiffError::AzimuthalResolution);
        }
        for (which, v) in [
            ("mu", self.mu),
            ("k_gas", self.k_gas),
            ("cv", self.cv),
            ("rho_d", self.rho_d),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(GasDiffError::BadCoefficient(which));
            }
        }
        Ok(())
    }

    /// Implicit two-point face coefficient of a component [per unit
    /// gradient·area]: the symmetric core the CG matrix carries.
    fn face_coef(&self, comp: GasComp, geom: &FaceGeom) -> f64 {
        match comp {
            GasComp::Ur => {
                if geom.radial {
                    (4.0 / 3.0) * self.mu
                } else {
                    self.mu
                }
            }
            GasComp::Uz => {
                if geom.radial {
                    self.mu
                } else {
                    (4.0 / 3.0) * self.mu
                }
            }
            // The angular-momentum form: flux = μ·A·r_f²·Δω/d (r-faces) or
            // μ·A·r̄²·Δω/d (z-faces, both cells share r̄).
            GasComp::Om => self.mu * geom.r_face * geom.r_face,
            GasComp::T => self.k_gas,
            GasComp::C => self.rho_d,
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
    #[allow(clippy::too_many_lines)]
    fn classify_faces<'b>(
        &'b self,
        g: &Grid,
        bi: usize,
        local: usize,
        i_r: usize,
        i_z: usize,
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
                area: g.face_area_r(i_r, false, 1),
                dist: dr,
                r_face: r0 + i_r as f64 * dr,
                high: -1.0,
                pos: (r0 + i_r as f64 * dr, zbar),
                radial: true,
            },
            FaceGeom {
                area: g.face_area_r(i_r, true, 1),
                dist: dr,
                r_face: r0 + (i_r + 1) as f64 * dr,
                high: 1.0,
                pos: (r0 + (i_r + 1) as f64 * dr, zbar),
                radial: true,
            },
            FaceGeom {
                area: g.face_area_z(i_r, 1),
                dist: dz,
                r_face: rbar,
                high: -1.0,
                pos: (rbar, z0 + i_z as f64 * dz),
                radial: false,
            },
            FaceGeom {
                area: g.face_area_z(i_r, 1),
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
    /// one-sided at walls/edges (first-order locally — module doc).
    pub(crate) fn fill_lag_gradients(
        &self,
        g: &Grid,
        lag: &GasOperands,
        work: &mut GasWork,
    ) -> Result<(), GasDiffError> {
        self.validate(g)?;
        let (dr, dz) = (g.spec().dr, g.spec().dz);
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
            for local in 0..BRICK_CELLS {
                work.dur_dr[bi][local] = 0.0;
                work.dur_dz[bi][local] = 0.0;
                work.duz_dr[bi][local] = 0.0;
                work.duz_dz[bi][local] = 0.0;
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
                // Gas-connected neighbor value of a field, else None.
                let val = |n: Nbr, f: &BufF, dir: FaceDir| -> Option<f64> {
                    if b.aperture_rz(dir, local) <= 0.0 {
                        return None;
                    }
                    match n {
                        Nbr::InBrick(off) => {
                            let nl = (local as isize + off) as usize;
                            (b.mask() & (1u64 << nl) != 0).then(|| f[bi][nl])
                        }
                        Nbr::Cross { bi: nbi, local: nl } => {
                            (g.brick(nbi).mask() & (1u64 << nl) != 0).then(|| f[nbi][nl])
                        }
                        _ => None,
                    }
                };
                let deriv = |f: &BufF, k_lo: usize, k_hi: usize, d: f64| -> f64 {
                    let c = f[bi][local];
                    let lo = val(nbrs[k_lo], f, dirs[k_lo]);
                    let hi = val(nbrs[k_hi], f, dirs[k_hi]);
                    match (lo, hi) {
                        (Some(a), Some(bv)) => (bv - a) / (2.0 * d),
                        (None, Some(bv)) => (bv - c) / d,
                        (Some(a), None) => (c - a) / d,
                        (None, None) => 0.0,
                    }
                };
                work.dur_dr[bi][local] = deriv(&lag.ur, 0, 1, dr);
                work.dur_dz[bi][local] = deriv(&lag.ur, 2, 3, dz);
                work.duz_dr[bi][local] = deriv(&lag.uz, 0, 1, dr);
                work.duz_dz[bi][local] = deriv(&lag.uz, 2, 3, dz);
            }
        }
        Ok(())
    }

    /// The matrix apply of one component's implicit core: `out = L·x` in
    /// solve units (force for u_r/u_z, angular-momentum torque for ω, W
    /// for T/C-flux), constants dropped exactly (the Linear discipline of
    /// the solid assembly). `diag` receives ∂out_i/∂x_i (≤ 0).
    pub(crate) fn apply_linear(
        &self,
        g: &Grid,
        comp: GasComp,
        x: &BufF,
        out: &mut BufF,
        mut diag: Option<&mut BufF>,
    ) {
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
            for local in 0..BRICK_CELLS {
                if b.mask() & (1u64 << local) == 0 {
                    continue;
                }
                let (lr, lz) = (local / BRICK, local % BRICK);
                let (i_r, i_z) = b.global_rz(local);
                let nbrs = Self::neighbors(g, i_r, i_z, lr, lz, &nb);
                let faces = self.classify_faces(g, bi, local, i_r, i_z, &nbrs);
                let x_c = x[bi][local];
                let mut acc = 0.0f64;
                let mut dg = 0.0f64;
                for (kind, geom) in &faces {
                    let coef = self.face_coef(comp, geom);
                    match kind {
                        FaceKind::Interior {
                            bi: nbi,
                            local: nl,
                            ap,
                        } => {
                            let a = geom.area * ap;
                            let x_n = x[*nbi][*nl];
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
                if comp == GasComp::Ur {
                    // The negative-definite geometric diagonal of −τ_θθ/r:
                    // −(4/3)μ·u_r/r̄ · geo·κV with the metric-consistent
                    // geo = (A_out − A_in)/V = 1/r̄ (SOLV-1 §3.3 pattern).
                    let vol = g.cell_volume(i_r, 1);
                    let geo = (g.face_area_r(i_r, true, 1) - g.face_area_r(i_r, false, 1)) / vol;
                    let kv = b.kappa_rz(local) * vol;
                    let c = (4.0 / 3.0) * self.mu * geo / g.r_center(i_r) * kv;
                    acc -= c * x_c;
                    dg -= c;
                }
                out[bi][local] = acc;
                if let Some(d) = diag.as_deref_mut() {
                    d[bi][local] = dg;
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
        time: f64,
        rates: &mut [Vec<Cons>],
        mut ledger: Option<&mut FlowLedger>,
    ) -> Result<(), GasDiffError> {
        self.validate(g)?;
        let two_thirds_mu = (2.0 / 3.0) * self.mu;
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
            for local in 0..BRICK_CELLS {
                if b.mask() & (1u64 << local) == 0 {
                    continue;
                }
                let (lr, lz) = (local / BRICK, local % BRICK);
                let (i_r, i_z) = b.global_rz(local);
                let rbar = g.r_center(i_r);
                let vol = g.cell_volume(i_r, 1);
                let kv = b.kappa_rz(local) * vol;
                let nbrs = Self::neighbors(g, i_r, i_z, lr, lz, &nb);
                let faces = self.classify_faces(g, bi, local, i_r, i_z, &nbrs);

                // The visiting cell's operands.
                let at = |f: &BufF| f[bi][local];
                let (ur_c, om_c, uz_c, tt_c, cc_c) = (
                    at(&sol.ur),
                    at(&sol.om),
                    at(&sol.uz),
                    at(&sol.tt),
                    at(&sol.cc),
                );
                for (what, v) in [
                    ("velocity", ur_c.abs() + om_c.abs() + uz_c.abs()),
                    ("temperature", tt_c),
                    ("composition operand", cc_c.abs()),
                ] {
                    if !v.is_finite() {
                        return Err(GasDiffError::NonFinite { i_r, i_z, what });
                    }
                }
                // Lag-state cell value for the cross terms.
                let e_thth_c = lag.ur[bi][local] / rbar;

                // Accumulators in total units; θ-momentum in λ (angular-
                // momentum) units. Fixed face order = the deterministic
                // accumulation contract.
                let mut tot = [0.0f64; NCOMP];
                let mut tot_lam = 0.0f64;
                let mut lam_abs = 0.0f64;

                for (kind, geom) in &faces {
                    match kind {
                        FaceKind::Interior {
                            bi: nbi,
                            local: nl,
                            ap,
                        } => {
                            let a = geom.area * ap;
                            let d = geom.dist;
                            let (ur_n, om_n, uz_n, tt_n, cc_n) = (
                                sol.ur[*nbi][*nl],
                                sol.om[*nbi][*nl],
                                sol.uz[*nbi][*nl],
                                sol.tt[*nbi][*nl],
                                sol.cc[*nbi][*nl],
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
                            let avg_w = |w: &BufF| 0.5 * (w[bi][local] + w[*nbi][*nl]);
                            let e_thth_f = if geom.radial {
                                let (n_ir, _) = g.brick(*nbi).global_rz(*nl);
                                0.5 * (e_thth_c + lag.ur[*nbi][*nl] / g.r_center(n_ir))
                            } else {
                                0.5 * (e_thth_c + lag.ur[*nbi][*nl] / rbar)
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
                                0.5 * (om_c * rbar + sol.om[*nbi][*nl] * g.r_center(n_ir))
                            } else {
                                0.5 * (om_c + om_n) * rbar
                            };
                            // The stress components at the face (implicit
                            // normal-gradient parts at sol; cross at lag).
                            let (f_mr, f_mz, tau_rr_or_zz, tau_rz, tau_th);
                            if geom.radial {
                                let tau_rr = (4.0 / 3.0) * self.mu * g_ur
                                    - two_thirds_mu * (e_thth_f + e_zz_f);
                                let t_rz = self.mu * (g_uz + dur_dz_f);
                                let t_rth = self.mu * geom.r_face * g_om;
                                f_mr = a * tau_rr;
                                f_mz = a * t_rz;
                                tau_rr_or_zz = tau_rr;
                                tau_rz = t_rz;
                                tau_th = t_rth;
                            } else {
                                let tau_zz = (4.0 / 3.0) * self.mu * g_uz
                                    - two_thirds_mu * (e_rr_f + e_thth_f);
                                let t_rz = self.mu * (g_ur + duz_dr_f);
                                let t_thz = self.mu * rbar * g_om;
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
                            // θ: angular-momentum flux = A·r_f·τ_(rθ|θz)
                            // (r-faces carry r_f = the face radius;
                            // z-faces r̄ — folded into tau_th? No: tau_th
                            // here is the plain stress; the λ-flux weight
                            // is r_face for both kinds).
                            let f_lam = a * geom.r_face * tau_th;
                            tot_lam += s * f_lam;
                            lam_abs += (f_lam / rbar).abs();
                            // Energy: work + conduction.
                            let g_e = if geom.radial {
                                ur_f * tau_rr_or_zz + ut_f * tau_th + uz_f * tau_rz
                            } else {
                                ur_f * tau_rz + ut_f * tau_th + uz_f * tau_rr_or_zz
                            } + self.k_gas * g_tt;
                            tot[I_EN] += s * a * g_e;
                            // Species.
                            tot[I_RC] += s * a * self.rho_d * g_cc;
                        }
                        FaceKind::Boundary { bc, ap } => {
                            if geom.area == 0.0 {
                                continue; // the axis face drops out
                            }
                            let a = geom.area * ap;
                            let half = 0.5 * geom.dist;
                            let theta = Grid::theta_center(0, 1);
                            let s = geom.high;
                            // One-sided lag pieces at the boundary (module
                            // doc: first-order locally).
                            let e_zz_f = work.duz_dz[bi][local];
                            let e_rr_f = work.dur_dr[bi][local];
                            let e_thth_f = e_thth_c;
                            let dur_dz_f = work.dur_dz[bi][local];
                            let duz_dr_f = work.duz_dr[bi][local];
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
                                let (tau_nn, tau_rz, tau_th);
                                if geom.radial {
                                    tau_nn = (4.0 / 3.0) * self.mu * g_ur
                                        - two_thirds_mu * (e_thth_f + e_zz_f);
                                    tau_rz = self.mu * (g_uz + dur_dz_f);
                                    tau_th = self.mu * geom.r_face * g_om;
                                    port[I_MR] += s * a * tau_nn;
                                    port[I_MZ] += s * a * tau_rz;
                                } else {
                                    tau_nn = (4.0 / 3.0) * self.mu * g_uz
                                        - two_thirds_mu * (e_rr_f + e_thth_f);
                                    tau_rz = self.mu * (g_ur + duz_dr_f);
                                    tau_th = self.mu * rbar * g_om;
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
                            if let Some(w_tt) = self.bc_value(GasComp::T, bc, geom, theta, time) {
                                let g_tt = s * (w_tt - tt_c) / half;
                                port[I_EN] += s * a * self.k_gas * g_tt;
                            }
                            if let Some(wc) = self.bc_value(GasComp::C, bc, geom, theta, time) {
                                let g_cc = s * (wc - cc_c) / half;
                                port[I_RC] += s * a * self.rho_d * g_cc;
                            }
                            for (t_k, p_k) in tot.iter_mut().zip(&port) {
                                *t_k += p_k;
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

                // −τ_θθ/r volume source of the r-momentum (implicit
                // diagonal + lagged compressible correction), with the
                // metric-consistent 1/r̄ (SOLV-1 §3.3 pattern).
                let geo = (g.face_area_r(i_r, true, 1) - g.face_area_r(i_r, false, 1)) / vol;
                let e_rr_c = work.dur_dr[bi][local];
                let e_zz_c = work.duz_dz[bi][local];
                let tau_thth =
                    (4.0 / 3.0) * self.mu * (ur_c / rbar) - two_thirds_mu * (e_rr_c + e_zz_c);
                let src_mr = -tau_thth * geo * kv;
                tot[I_MR] += src_mr;

                // Commit: conserved-density rates (θ from λ/(r̄κV)).
                let inv_kv = 1.0 / kv;
                let rate = &mut rates[bi][local];
                rate[I_MR] = tot[I_MR] * inv_kv;
                rate[I_MT] = tot_lam / (rbar * kv);
                rate[I_MZ] = tot[I_MZ] * inv_kv;
                rate[I_EN] = tot[I_EN] * inv_kv;
                rate[I_RC] = tot[I_RC] * inv_kv;

                if let Some(l) = ledger.as_deref_mut() {
                    // Interior fluxes telescope for r/z-momentum, energy,
                    // species — only BC ports (above) and the volume/θ
                    // sources are ledgered. Gross magnitudes feed S[q].
                    l.src_net[I_MR] += src_mr;
                    l.src_abs[I_MR] += src_mr.abs();
                    l.src_net[I_MT] += tot_lam / rbar;
                    l.src_abs[I_MT] += lam_abs;
                    for kk in [I_MR, I_MZ, I_EN, I_RC] {
                        l.port_abs[kk] += tot[kk].abs();
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
        wqnew: f64,
        x: &mut BufF,
        work: &mut GasWork,
    ) -> Result<(usize, f64), GasDiffError> {
        let nb = g.n_bricks();
        let masked = |bi: usize| g.brick(bi).mask();
        let dot = |a: &BufF, c: &BufF| -> f64 {
            let partials: Vec<f64> = (0..nb)
                .map(|bi| {
                    let m = masked(bi);
                    let mut acc = 0.0f64;
                    for local in 0..BRICK_CELLS {
                        if m & (1u64 << local) != 0 {
                            acc += a[bi][local] * c[bi][local];
                        }
                    }
                    acc
                })
                .collect();
            tree_combine(&partials)
        };

        // A's diagonal: mass − wqnew·diag(L)  (diag(L) ≤ 0 ⇒ positive).
        self.apply_linear(g, comp, x, &mut work.q, Some(&mut work.diag));
        for bi in 0..nb {
            let m = masked(bi);
            for local in 0..BRICK_CELLS {
                if m & (1u64 << local) == 0 {
                    work.diag[bi][local] = 1.0;
                } else {
                    work.diag[bi][local] = work.mass[bi][local] - wqnew * work.diag[bi][local];
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
            for local in 0..BRICK_CELLS {
                work.z[bi][local] = work.r[bi][local] / work.diag[bi][local];
            }
            work.p[bi].copy_from_slice(&work.z[bi]);
        }
        let mut rz = dot(&work.r, &work.z);
        let mut r_norm2 = dot(&work.r, &work.r);
        let mut iters = 0usize;
        while iters < N_CG_ITERS_MAX && r_norm2 > eps2 && rz > 0.0 {
            // q = (mass − wqnew·L)·p.
            self.apply_linear(g, comp, &work.p, &mut work.q, None);
            for bi in 0..nb {
                let m = masked(bi);
                for local in 0..BRICK_CELLS {
                    if m & (1u64 << local) == 0 {
                        work.q[bi][local] = 0.0;
                    } else {
                        work.q[bi][local] =
                            work.mass[bi][local] * work.p[bi][local] - wqnew * work.q[bi][local];
                    }
                }
            }
            let pq = dot(&work.p, &work.q);
            if pq.is_nan() || pq <= 0.0 {
                break; // SPD breakdown ⇒ the acceptance below decides
            }
            let alpha = rz / pq;
            for bi in 0..nb {
                for local in 0..BRICK_CELLS {
                    work.delta[bi][local] += alpha * work.p[bi][local];
                    work.r[bi][local] -= alpha * work.q[bi][local];
                }
            }
            for bi in 0..nb {
                for local in 0..BRICK_CELLS {
                    work.z[bi][local] = work.r[bi][local] / work.diag[bi][local];
                }
            }
            let rz_new = dot(&work.r, &work.z);
            let beta = rz_new / rz;
            rz = rz_new;
            for bi in 0..nb {
                for local in 0..BRICK_CELLS {
                    work.p[bi][local] = work.z[bi][local] + beta * work.p[bi][local];
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
            for (local, xv) in xb.iter_mut().enumerate() {
                if m & (1u64 << local) != 0 {
                    *xv += work.delta[bi][local];
                }
            }
        }
        Ok((iters, resid))
    }

    /// Fill `work.mass` with one component's per-cell mass:
    /// u_r/u_z: ρκV; ω: ρr̄²κV (angular-momentum form); T: ρc_vκV; C: ρκV.
    pub(crate) fn fill_mass(&self, g: &Grid, comp: GasComp, rho: &BufF, work: &mut GasWork) {
        for (bi, mb) in work.mass.iter_mut().enumerate() {
            let b = g.brick(bi);
            for (local, mv) in mb.iter_mut().enumerate() {
                if b.mask() & (1u64 << local) == 0 {
                    *mv = 1.0;
                    continue;
                }
                let (i_r, _) = b.global_rz(local);
                let kv = b.kappa_rz(local) * g.cell_volume(i_r, 1);
                let rho_c = rho[bi][local];
                *mv = match comp {
                    GasComp::Ur | GasComp::Uz | GasComp::C => rho_c * kv,
                    GasComp::Om => {
                        let r = g.r_center(i_r);
                        rho_c * r * r * kv
                    }
                    GasComp::T => rho_c * self.cv * kv,
                };
            }
        }
    }

    /// The T-solve's linearization slope c_v [J/(kg·K)] (exact for the
    /// gamma-law class; the S4 spine owns the general case).
    pub fn cv(&self) -> f64 {
        self.cv
    }
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
        let law = WallLaw::new(
            specific_heat_capacity_j_per_kg_k(CP),
            dynamic_viscosity_pa_s(MU),
            PR,
        )
        .expect("transport");
        GasDiffusion::from_transport(
            &law,
            GAMMA,
            SC,
            GasDiffBcs {
                r_inner: FaceGasBc::free(),
                r_outer: FaceGasBc::free(),
                z_lo: FaceGasBc::free(),
                z_hi: FaceGasBc::free(),
            },
        )
        .expect("operator")
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
        op.fill_lag_gradients(g, sol, &mut work).expect("gradients");
        op.assemble_rates(g, sol, sol, &work, 0.0, &mut rates, Some(&mut ledger))
            .expect("assembly");
        (rates, ledger)
    }

    /// One transport owner: the derived constants restate NOTHING.
    #[test]
    fn transport_derivation_by_hand() {
        let o = op();
        assert_eq!(o.mu, MU);
        assert_eq!(o.k_gas, MU * CP / PR);
        assert_eq!(o.cv, CP / GAMMA);
        assert_eq!(o.rho_d, MU / SC);
        let law = WallLaw::new(
            specific_heat_capacity_j_per_kg_k(CP),
            dynamic_viscosity_pa_s(MU),
            PR,
        )
        .unwrap();
        assert!(
            GasDiffusion::from_transport(
                &law,
                1.0,
                SC,
                GasDiffBcs {
                    r_inner: FaceGasBc::free(),
                    r_outer: FaceGasBc::free(),
                    z_lo: FaceGasBc::free(),
                    z_hi: FaceGasBc::free(),
                }
            )
            .is_err()
        );
        assert!(
            GasDiffusion::from_transport(
                &law,
                GAMMA,
                0.0,
                GasDiffBcs {
                    r_inner: FaceGasBc::free(),
                    r_outer: FaceGasBc::free(),
                    z_lo: FaceGasBc::free(),
                    z_hi: FaceGasBc::free(),
                }
            )
            .is_err()
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
        let expect = o.k_gas * a_z * ((t_at(i_z + 1) - t_at(i_z)) - (t_at(i_z) - t_at(i_z - 1)))
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
        let expect = o.rho_d * a_z * ((c_at(i_z + 1) - c_at(i_z)) - (c_at(i_z) - c_at(i_z - 1)))
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
        let expect = o.mu
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
        let expect = o.mu
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
        let tau_out = (4.0 / 3.0) * o.mu * (ur_at(i_r + 1) - ur_at(i_r)) / dr
            - (2.0 / 3.0) * o.mu * e_thth_out;
        let tau_in = (4.0 / 3.0) * o.mu * (ur_at(i_r) - ur_at(i_r - 1)) / dr
            - (2.0 / 3.0) * o.mu * e_thth_in;
        let tau_thth = (4.0 / 3.0) * o.mu * ur_at(i_r) / rbar - (2.0 / 3.0) * o.mu * e_rr;
        let expect = (a_out * tau_out - a_in * tau_in) / vol - tau_thth * geo;
        let got = rates[0][local][I_MR];
        assert!(
            (got - expect).abs() <= 1e-12 * expect.abs().max(1e-30),
            "u_r core: got {got}, expect {expect}"
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
}
