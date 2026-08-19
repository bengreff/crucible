//! SOLV-1 §3.1–§3.3 — the **conserved-variable field operator**, reacting-gas
//! subset, first inhabitant: compressible Euler on the exact cylindrical
//! metric. One flux-form law over the state vector
//! `U = (ρ, ρu_r, ρu_θ, ρu_z, ρE, ρC)` — no material or regime branch; the
//! axis face has zero area and drops out geometrically; at `N_θ = 1` the two
//! θ-face fluxes of the ring cell cancel identically (same computed bits),
//! so the axisymmetric-with-swirl corner is the same operator, not a special
//! case (SOLV-1 §3.3).
//!
//! Discretization (SOLV-1 §3.2): Godunov FV — PPM reconstruction in
//! primitives (`recon`), HLLC flux with Batten wavespeeds (`hllc`),
//! well-balanced geometric sources (§3.3: the radial-momentum pressure
//! source uses the same `A·p` products as the face fluxes, so a uniform
//! state at rest is preserved to round-off, axis included). The composition
//! slot `ρC` is one passive advected mass fraction — the elemental block of
//! §3.1 at its smallest (source-free advection; the OFFL-3 equilibrium
//! projection arrives with station 3).
//!
//! Time integration note (honest scaffolding, the session-5 pattern): the
//! explicit MOL SSP-RK2 here is the fixed-order reference integrator that
//! drives the Station-1 certificate; it is superseded — not extended — by
//! COUP-3's SDC-IMEX schedule, which consumes the same flux-form spatial
//! operator (Castro's SDC path uses exactly this MOL reconstruction, hence
//! no characteristic tracing — that is the split scheme's predictor).
//!
//! EOS seam (SOLV-1 §3.4, the FND-7 spine boundary): the operator is generic
//! over [`EosLaw`] — monomorphized, no dynamic dispatch in hot loops (FND-2
//! §3.9). Two occupants: [`GammaLaw`] (the degenerate spine occupant; γ is
//! config data) and [`TableEos`] (shifting-equilibrium mode — the OFFL-3
//! (p, h, Z) surface with the fixed-count per-cell equilibrium projection).
//! The primitive vector carries two **auxiliary EOS slots** past the six
//! physical ones — specific internal energy `e` and the effective adiabatic
//! exponent `Γ₁ = ρa²/p` — reconstructed componentwise like every primitive
//! and read only by the general-EOS closures at faces (the Castro/PeleC
//! treatment of a general convex EOS in the Riemann solve: Batten wavespeeds
//! from local/Roe-averaged Γ₁). `GammaLaw` never reads the aux slots, so the
//! gamma-law path is arithmetically identical to the pre-seam operator
//! (certificates byte-identical — asserted by gate 5).
//!
//! Deferred, loud (owners named): `r_min = 0` with `N_θ > 1` refuses (the
//! cross-axis θ↔θ+π parity-pair gather, FND-2 §3.2, lands with the first
//! 3-D-across-the-axis wave); mixed per-brick N_θ refuses (refluxing =
//! COUP-2/COUP-3); apertures/cut cells (FND-3) not yet consumed — worlds
//! are full boxes.

mod exact;
mod hllc;
mod recon;
mod table_eos;

pub use exact::{RiemannSide, RiemannSolution, solve as solve_riemann};
pub use hllc::{hllc_flux, physical_flux};
use recon::{NGHOST, ppm_faces};
pub use table_eos::{EPS_P_PROJECTION, N_INFLOW_ITER, N_P_ITER_MAX, TableEos};

use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, FieldId, Grid, GridError};
use rayon::prelude::*;

/// SOLV-1 §3.6 / FND-3 §3.4 — the State Redistribution small-cell
/// threshold (Berger & Giuliani 2020, META-3 `state-redistribution`): a
/// cut cell with gas volume fraction κ below this merges its provisional
/// update into a flow-connected neighborhood whose summed κ reaches it.
/// The paper's standard half-cell target; named constant, never tuned
/// per geometry.
pub const KAPPA_SRD: f64 = 0.5;

/// A merge neighborhood: member cells as global (i_r, i_z) with their gas
/// volumes κ·V (per θ-plane) — the SRD weights.
pub type SrdNeighborhood = Vec<((usize, usize), f64)>;

/// The deterministic SRD merging neighborhood of cell (i_r, i_z): itself,
/// then flow-connected face neighbors (shared aperture > 0) added in
/// descending-κ order (ties: the fixed FaceDir enumeration order) until
/// `Σκ ≥ KAPPA_SRD`. `None` for a regular cell (κ ≥ threshold). The
/// under-resolved case — no reachable neighborhood — refuses loudly.
/// Public because the engine's wall-exchange debit must deposit into the
/// same merged control volume SRD stabilizes (single owner of the rule).
/// Returned weights are the members' gas volumes κ·V (per θ-plane).
pub fn srd_neighborhood(
    g: &Grid,
    i_r: usize,
    i_z: usize,
) -> Result<Option<SrdNeighborhood>, FlowError> {
    let nt = g.brick(0).n_theta();
    if g.kappa(i_r, i_z) >= KAPPA_SRD {
        return Ok(None);
    }
    let members = neighborhood_cells(g, i_r, i_z)?;
    Ok(Some(
        members
            .into_iter()
            .map(|(r, z)| ((r, z), g.kappa(r, z) * g.cell_volume(r, nt)))
            .collect(),
    ))
}

/// The one neighborhood-construction rule (see [`srd_neighborhood`]).
fn neighborhood_cells(g: &Grid, i_r: usize, i_z: usize) -> Result<Vec<(usize, usize)>, FlowError> {
    let mut members = vec![(i_r, i_z)];
    let mut sum = g.kappa(i_r, i_z);
    if sum >= KAPPA_SRD {
        return Ok(members);
    }
    let nbrs = [
        (FaceDir::RMinus, i_r.checked_sub(1).map(|r| (r, i_z))),
        (FaceDir::RPlus, Some((i_r + 1, i_z))),
        (FaceDir::ZMinus, i_z.checked_sub(1).map(|z| (i_r, z))),
        (FaceDir::ZPlus, Some((i_r, i_z + 1))),
    ];
    let mut cand: Vec<((usize, usize), f64)> = Vec::new();
    for (dir, nbr) in nbrs {
        if let Some((nr, nz)) = nbr
            && g.aperture(i_r, i_z, dir) > 0.0
            && g.is_active(nr, nz)
        {
            cand.push(((nr, nz), g.kappa(nr, nz)));
        }
    }
    // Stable sort: descending κ, FaceDir order breaking ties.
    cand.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("finite κ"));
    for ((nr, nz), k) in cand {
        members.push((nr, nz));
        sum += k;
        if sum >= KAPPA_SRD {
            return Ok(members);
        }
    }
    Err(FlowError::CutCellUnmergeable { i_r, i_z })
}

/// Components of `U` (SOLV-1 §3.1) and of the primitive view
/// `W = (ρ, u_r, u_θ, u_z, p, C | e, Γ₁)`. Slots 1–3 are the velocity/
/// momentum directions, so a face's normal is named by its slot index.
/// `W` carries `NPRIM − NCOMP = 2` auxiliary EOS slots (specific internal
/// energy, effective Γ₁) with no conserved counterpart — filled by
/// `EosLaw::prim_checked`, reconstructed componentwise, read only by the
/// general-EOS face closures (`GammaLaw` ignores them).
pub const NCOMP: usize = 6;
pub const NPRIM: usize = 8;
pub const I_RHO: usize = 0;
pub const I_MR: usize = 1;
pub const I_MT: usize = 2;
pub const I_MZ: usize = 3;
pub const I_EN: usize = 4;
pub const I_RC: usize = 5;
/// Aux primitive slot: specific internal energy `e` (J/kg).
pub const I_EI: usize = 6;
/// Aux primitive slot: effective adiabatic exponent `Γ₁ = ρa²/p`.
pub const I_G1: usize = 7;

pub type Cons = [f64; NCOMP];
pub type Prim = [f64; NPRIM];

/// Build a primitive state from the six physical slots, aux slots zero.
/// Correct for `GammaLaw` (which never reads aux); a `TableEos` primitive
/// must come from `TableEos::prim_checked` (which fills them).
pub const fn prim6(rho: f64, u_r: f64, u_t: f64, u_z: f64, p: f64, c: f64) -> Prim {
    [rho, u_r, u_t, u_z, p, c, 0.0, 0.0]
}

/// SOLV-1 §3.4 — the constitutive closure the operator is generic over (the
/// FND-7 spine boundary). Implementations read only the state (`U`/`W`),
/// never a material or regime label (Rule 12); which occupant runs is config
/// data (Rule 13). Monomorphized into the sweeps — no hot-loop dispatch.
pub trait EosLaw {
    /// Conserved → primitive (+aux slots), with the META-1 P6 checks: a
    /// non-physical or off-table state is a halt diagnosis (the caller
    /// attaches the cell location), never a clamp.
    fn prim_checked(&self, u: &Cons) -> Result<Prim, &'static str>;
    /// [`EosLaw::prim_checked`] with an optional warm-start hint (the
    /// cell's previous projected pressure). Occupants whose conversion has
    /// no iterative solve ignore it (the default — `GammaLaw` stays
    /// arithmetically identical); `TableEos` uses it to tighten the
    /// projection bracket. The hint is an acceleration, never a physics
    /// input: any valid state converts identically with or without it (to
    /// the projection tolerance).
    fn prim_checked_hinted(&self, u: &Cons, _hint: Option<f64>) -> Result<Prim, &'static str> {
        self.prim_checked(u)
    }
    /// Primitive → conserved.
    fn prim_to_cons(&self, w: &Prim) -> Cons;
    /// Total energy density ρE at a (possibly face-reconstructed) state.
    fn total_energy(&self, w: &Prim) -> f64;
    /// Sound speed at a (possibly face-reconstructed) state.
    fn sound_speed_w(&self, w: &Prim) -> f64;
    /// Roe-average sound speed for the Batten wavespeed bounds; receives the
    /// already-formed Roe enthalpy/velocity ingredients so every occupant
    /// shares one Roe algebra.
    #[allow(clippy::too_many_arguments)]
    fn roe_sound_speed(
        &self,
        wl: &Prim,
        wr: &Prim,
        h_roe: f64,
        q2_roe: f64,
        sql: f64,
        sqr: f64,
        inv: f64,
    ) -> f64;
    /// Reservoir-isentrope ghost for [`FlowBc::StagnationInflow`] — an
    /// EOS-specific closed form. Occupants without one refuse loudly (the
    /// COUP-7 injector object owns inflow for the table EOS).
    fn stagnation_ghost(
        &self,
        p0: f64,
        rho0: f64,
        c_frac: f64,
        u_n: f64,
        normal: usize,
    ) -> Result<Prim, FlowError>;

    /// Declared-mass-flux inflow ghost for [`FlowBc::MassFlowInflow`]:
    /// interior static pressure `p_int`, declared `(ṁ/A, h_total, c_frac)`,
    /// `sign` = +1 entering through a low face, −1 through a high face,
    /// `normal` the velocity slot. Default: refuse (occupants opt in with
    /// their own deterministic solve).
    fn mass_flow_inflow_ghost(
        &self,
        _mdot_per_area: f64,
        _h_total: f64,
        _c_frac: f64,
        _p_int: f64,
        _sign: f64,
        _normal: usize,
    ) -> Result<Prim, FlowError> {
        Err(FlowError::BcUnsupportedByEos {
            bc: "MassFlowInflow",
        })
    }
}

/// Grid field names for `U`, in component order (FND-2 §3.4).
pub const EULER_FIELDS: &[&str] = &["rho", "mom_r", "mom_theta", "mom_z", "rho_e", "rho_c"];

/// Gamma-law EOS — the first, degenerate occupant of the FND-7 spine seam.
/// γ is pure data (config parameter); nothing here reads a material label.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GammaLaw {
    pub gamma: f64,
}

impl GammaLaw {
    #[inline]
    pub fn sound_speed(&self, rho: f64, p: f64) -> f64 {
        (self.gamma * p / rho).sqrt()
    }

    /// Total energy density ρE from a primitive state.
    #[inline]
    pub fn total_energy(&self, w: &Prim) -> f64 {
        let rho = w[I_RHO];
        w[4] / (self.gamma - 1.0) + 0.5 * rho * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3])
    }

    #[inline]
    pub fn prim_to_cons(&self, w: &Prim) -> Cons {
        let rho = w[I_RHO];
        [
            rho,
            rho * w[1],
            rho * w[2],
            rho * w[3],
            self.total_energy(w),
            rho * w[I_RC],
        ]
    }

    /// Conserved → primitive with the META-1 P6 checks: non-finite or
    /// non-positive ρ/p is a halt (the caller attaches the cell diagnosis),
    /// never a clamp.
    #[inline]
    pub fn prim_checked(&self, u: &Cons) -> Result<Prim, &'static str> {
        let rho = u[I_RHO];
        if !rho.is_finite() || rho <= 0.0 {
            return Err("non-positive or non-finite density");
        }
        let inv = 1.0 / rho;
        let (ur, ut, uz) = (u[I_MR] * inv, u[I_MT] * inv, u[I_MZ] * inv);
        let ke = 0.5 * rho * (ur * ur + ut * ut + uz * uz);
        let p = (self.gamma - 1.0) * (u[I_EN] - ke);
        if !p.is_finite() || p <= 0.0 {
            return Err("non-positive or non-finite pressure");
        }
        let c = u[I_RC] * inv;
        if !c.is_finite() {
            return Err("non-finite composition");
        }
        Ok([rho, ur, ut, uz, p, c, 0.0, 0.0])
    }
}

impl EosLaw for GammaLaw {
    #[inline]
    fn prim_checked(&self, u: &Cons) -> Result<Prim, &'static str> {
        GammaLaw::prim_checked(self, u)
    }

    #[inline]
    fn prim_to_cons(&self, w: &Prim) -> Cons {
        GammaLaw::prim_to_cons(self, w)
    }

    #[inline]
    fn total_energy(&self, w: &Prim) -> f64 {
        GammaLaw::total_energy(self, w)
    }

    #[inline]
    fn sound_speed_w(&self, w: &Prim) -> f64 {
        // Same arithmetic as `sound_speed(rho, p)` — γ from config, never
        // from the (unread) aux slots.
        (self.gamma * w[4] / w[I_RHO]).sqrt()
    }

    #[inline]
    fn roe_sound_speed(
        &self,
        _wl: &Prim,
        _wr: &Prim,
        h_roe: f64,
        q2_roe: f64,
        _sql: f64,
        _sqr: f64,
        _inv: f64,
    ) -> f64 {
        ((self.gamma - 1.0) * (h_roe - 0.5 * q2_roe))
            .max(0.0)
            .sqrt()
    }

    /// Reservoir-isentrope ghost state from the interior-extrapolated
    /// normal velocity (see [`FlowBc::StagnationInflow`]).
    fn stagnation_ghost(
        &self,
        p0: f64,
        rho0: f64,
        c_frac: f64,
        u_n: f64,
        normal: usize,
    ) -> Result<Prim, FlowError> {
        let ga = self.gamma;
        let c0sq = ga * p0 / rho0;
        let csq = c0sq - 0.5 * (ga - 1.0) * u_n * u_n;
        if csq <= 0.0 {
            return Err(FlowError::InflowBeyondVacuumLimit);
        }
        let p = p0 * (csq / c0sq).powf(ga / (ga - 1.0));
        let rho = ga * p / csq;
        let mut m = prim6(rho, 0.0, 0.0, 0.0, p, c_frac);
        m[normal] = u_n;
        Ok(m)
    }
}

/// Boundary condition on one grid face (COUP-7 replaces these with
/// declarative boundary objects; same tracked deferral as conduction's BCs).
pub enum FlowBc<'a> {
    /// Solid wall: mirror ghost with the normal velocity negated.
    Reflecting,
    /// Zero-gradient outflow (first-order extrapolation).
    Transmissive,
    /// Prescribed primitive state at `(r, θ, z, t)` — MMS verification and
    /// supersonic inflow.
    Prescribed(&'a (dyn Fn(f64, f64, f64, f64) -> Prim + Sync)),
    /// Subsonic reservoir inflow: the ghost state sits on the isentrope of
    /// a stagnation reservoir `(p0, ρ0)` at the interior-extrapolated
    /// normal velocity — `c² = c0² − ½(γ−1)u²`, `p = p0·(c²/c0²)^(γ/(γ−1))`
    /// (the standard total-condition inflow; the COUP-7 injector object
    /// supersedes this with declared provenance/envelope/band).
    StagnationInflow { p0: f64, rho0: f64, c_frac: f64 },
    /// COUP-7 §3.2.1 prior-tier injector inflow: declared mass flux
    /// `ṁ/A` [kg s⁻¹ m⁻²], total enthalpy `h_total` [J/kg] (table
    /// coordinate), composition `c_frac` (= Z in shifting mode). Pressure
    /// extrapolates from the interior (standard subsonic inflow); the EOS
    /// occupant solves the face state by its own fixed-count iteration
    /// (`EosLaw::mass_flow_inflow_ghost`). The emergent-quantity rule lives
    /// here: the boundary states ṁ and inlet enthalpy; `p_c` is whatever
    /// the field produces (COUP-7 §3.2).
    MassFlowInflow {
        mdot_per_area: f64,
        h_total: f64,
        c_frac: f64,
    },
    /// Ambient-pressure outflow (the vacuum-plume seam at the table-envelope
    /// floor): while the interior normal flow at the face is subsonic the
    /// ghost carries the interior state with its pressure replaced by
    /// `p_ambient` (the standard subsonic pressure outlet — this is what
    /// lets a quiescent fill column drain and the nozzle establish); once
    /// the exit runs supersonic the branch condition makes it pure
    /// extrapolation, and the imposed value has no upstream influence at
    /// all (characteristics all leave). The branch is a fixed comparison on
    /// the data — deterministic. The ambient is a declared function of time
    /// (a deterministic startup schedule — e.g. the altitude-cell pump-down
    /// that lets a vacuum nozzle establish quasi-statically instead of
    /// through a violent drain); a constant closure is the steady form.
    PressureOutflow(&'a (dyn Fn(f64) -> f64 + Sync)),
}

pub struct FlowBcs<'a> {
    /// Ignored when `r_min = 0`: the axis is geometric (zero-area face +
    /// through-axis mirror ghosts with u_r AND u_θ negated — the N_θ = 1
    /// degenerate form of the θ↔θ+π parity pairing, FND-2 §3.2).
    pub r_inner: FlowBc<'a>,
    pub r_outer: FlowBc<'a>,
    pub z_lo: FlowBc<'a>,
    pub z_hi: FlowBc<'a>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FlowError {
    /// Same session rule as conduction: refluxing across an N_θ jump
    /// arrives with COUP-2/COUP-3; refuse rather than guess.
    MixedThetaResolution,
    /// Cross-axis parity-pair gather at N_θ > 1 is a deferred wave.
    AxisWithAzimuthalResolution,
    /// META-1 P6: halt with diagnosis (mechanism, location), never clamp.
    NonPhysicalState {
        i_r: usize,
        i_z: usize,
        i_theta: u32,
        what: &'static str,
    },
    /// A `StagnationInflow` face saw an interior velocity beyond the
    /// reservoir's vacuum limit `u² < 2c0²/(γ−1)` — the isentrope has no
    /// state there; halt, never clamp (META-1 P6).
    InflowBeyondVacuumLimit,
    /// A configured BC has no closed form under the selected EOS occupant
    /// (e.g. `StagnationInflow` under `TableEos` — the COUP-7 injector
    /// object owns inflow there). Refuse at use, never approximate.
    BcUnsupportedByEos { bc: &'static str },
    /// A cut cell below [`KAPPA_SRD`] has no flow-connected neighborhood
    /// reaching the merge target — the geometry has a feature thinner than
    /// the cell (FND-2 §5 under-resolved class). Refine or fix the
    /// contour; never smear (META-1 P6).
    CutCellUnmergeable { i_r: usize, i_z: usize },
}

impl std::fmt::Display for FlowError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MixedThetaResolution => write!(
                f,
                "mixed per-brick N_θ in one flow sweep — conservative flux aggregation \
                 across an N_θ jump arrives with COUP-2/COUP-3; refusing rather than guessing"
            ),
            Self::AxisWithAzimuthalResolution => write!(
                f,
                "r_min = 0 with N_θ > 1: the cross-axis θ↔θ+π parity-pair gather \
                 (FND-2 §3.2) is a deferred wave; refusing rather than guessing"
            ),
            Self::NonPhysicalState {
                i_r,
                i_z,
                i_theta,
                what,
            } => write!(
                f,
                "{what} at (i_r={i_r}, i_z={i_z}, i_θ={i_theta}) — halt with diagnosis"
            ),
            Self::InflowBeyondVacuumLimit => write!(
                f,
                "stagnation-inflow face saw interior velocity beyond the reservoir's \
                 vacuum limit — no state on the isentrope; halt, never clamp"
            ),
            Self::BcUnsupportedByEos { bc } => write!(
                f,
                "{bc} has no closed form under the selected EOS occupant — refusing \
                 rather than approximating (the COUP-7 boundary object owns this inflow)"
            ),
            Self::CutCellUnmergeable { i_r, i_z } => write!(
                f,
                "cut cell (i_r={i_r}, i_z={i_z}) below the SRD threshold has no \
                 flow-connected merge neighborhood — geometry under-resolved at this \
                 cell size (feature thinner than a cell); refine, never smear"
            ),
        }
    }
}

impl std::error::Error for FlowError {}

/// Resolved field handles for the six components of `U`.
#[derive(Debug, Clone, Copy)]
pub struct EulerFields {
    ids: [FieldId; NCOMP],
}

impl EulerFields {
    pub fn resolve(g: &Grid) -> Result<Self, GridError> {
        let mut ids = [g.field_id(EULER_FIELDS[0])?; NCOMP];
        for (k, name) in EULER_FIELDS.iter().enumerate() {
            ids[k] = g.field_id(name)?;
        }
        Ok(Self { ids })
    }

    #[inline]
    pub fn ids(&self) -> [FieldId; NCOMP] {
        self.ids
    }
}

/// Set every component of `U` from a primitive-state function of the cell
/// centroid (initial conditions; canonical traversal order).
pub fn fill_from_prim(
    g: &mut Grid,
    f: &EulerFields,
    eos: &impl EosLaw,
    func: impl Fn(f64, f64, f64) -> Prim,
) {
    for k in 0..NCOMP {
        let id = f.ids[k];
        g.fill_field(id, |r, th, z| eos.prim_to_cons(&func(r, th, z))[k]);
    }
}

/// The operator (SOLV-1 §1.1): reconstruction + Riemann flux + geometric
/// sources, writing `U` in place. `source` is the external per-component
/// volumetric source intake (MMS verification now; the SOLV-4 reaction
/// intake is this same slot). Generic over the [`EosLaw`] occupant
/// (monomorphized; `GammaLaw` default keeps existing fixtures unchanged).
pub struct Euler<'a, E: EosLaw = GammaLaw> {
    pub eos: E,
    pub source: &'a (dyn Fn(f64, f64, f64, f64) -> Cons + Sync),
    pub bcs: FlowBcs<'a>,
    /// Outward unit wall normal `(n_r, n_z)` of the true (smooth) wall at
    /// `(r, z)`, for masked stair-step walls: when `Some`, ghost states at
    /// **interior** run-boundary faces mirror the velocity about the true
    /// wall tangent (`v − 2(v·n̂)n̂` — ghost-cell immersed-boundary slip
    /// wall) instead of the grid-aligned stair face, cutting spurious wave
    /// generation from O(wall slope) to O(h·curvature). `None` → grid-
    /// aligned mirror (walls that lie exactly on grid faces). Superseded by
    /// FND-3's partial apertures + cut cells when that wave lands.
    pub wall_normal: Option<&'a (dyn Fn(f64, f64) -> (f64, f64) + Sync)>,
    /// Whether the slip ghost also applies at stair-step **z-faces** (the
    /// axial faces where the wall column changes). At a step, the slip
    /// ghost lets near-tangent flow glide *through* the face (the declared
    /// stair transpiration, station 2); in hypersonic flow that flux can
    /// STARVE the cell just downstream of the step (density runaway → CFL
    /// collapse — observed at the RL10 exit lip at M ≈ 4.4). `false` keeps
    /// slip on wall-parallel r-faces (the spurious-wave fix that matters)
    /// while step z-faces take the non-transpiring grid-aligned mirror.
    /// Numerics policy, pure data (Rule 13); retired with FND-3's cut cells
    /// + State Redistribution.
    pub slip_wall_z_faces: bool,
}

/// Per-brick scratch: primitives (AoS is fine for the CPU reference path;
/// the grid state itself stays SoA per FND-2 §3.9) and the RK stage data.
struct Scratch {
    /// `NPRIM`-wide primitive (+aux) states per cell.
    prim: Vec<Vec<Prim>>,
    rate: Vec<Vec<Cons>>,
    u0: Vec<Vec<Cons>>,
    /// Brick index by (br·nbz + bz) — resolved once, not per cell; `None`
    /// where a fully-inactive brick was never allocated (masked worlds).
    bmap: Vec<Option<usize>>,
    nbz: usize,
    /// Cell activity by (i_r·n_z + i_z): the pencil sweeps decompose each
    /// line into maximal active runs against this map; a run boundary in
    /// the interior is a stair-step wall face (reflecting), a run boundary
    /// at the domain edge takes the configured BC.
    act: Vec<bool>,
    /// State-Redistribution data — present iff the grid carries FND-3 cut
    /// geometry. Neighborhoods are geometry-time data (fixed for the run).
    cut: Option<CutScratch>,
    /// Whether `prim` holds a previous fill (⇒ usable warm-start hints).
    primed: bool,
}

/// Cells addressed as (brick index, local index) — the sweeps' native form.
type BrickLocalCells = Vec<(usize, usize)>;

/// Opaque persistent stepping workspace (see [`Euler::workspace`]).
pub struct EulerWorkspace(Scratch);

/// SRD bookkeeping (Berger & Giuliani 2020): the small-cell neighborhoods
/// in fixed lexicographic owner order and the overlap counts n_j.
struct CutScratch {
    /// Per small cell: members as (bi, local) with members[0] = owner, and
    /// the members' gas volumes κ·V (per θ-plane) — the merge weights.
    small: Vec<(BrickLocalCells, Vec<f64>)>,
    /// Overlap count n_j per (bi·BRICK_CELLS + local): 1 + the number of
    /// OTHER cells' neighborhoods containing j (every cell's own
    /// neighborhood — trivial for regular cells — is the 1).
    counts: Vec<u32>,
}

impl<E: EosLaw + Sync> Euler<'_, E> {
    fn validate(&self, g: &Grid) -> Result<u32, FlowError> {
        let nt = g.brick(0).n_theta();
        if g.bricks().iter().any(|b| b.n_theta() != nt) {
            return Err(FlowError::MixedThetaResolution);
        }
        if g.spec().r_min == 0.0 && nt > 1 {
            return Err(FlowError::AxisWithAzimuthalResolution);
        }
        Ok(nt)
    }

    fn scratch(&self, g: &Grid, nt: u32) -> Result<Scratch, FlowError> {
        let plane = nt as usize * BRICK_CELLS;
        let nb = g.n_bricks();
        let (n_r, n_z) = (g.spec().n_r, g.spec().n_z);
        let nbr = n_r.div_ceil(BRICK);
        let nbz = n_z.div_ceil(BRICK);
        let mut bmap = vec![None; nbr * nbz];
        for br in 0..nbr {
            for bz in 0..nbz {
                bmap[br * nbz + bz] = g.brick_index_by_coords(br as u32, bz as u32);
            }
        }
        let mut act = vec![false; n_r * n_z];
        for b in g.bricks() {
            for local in 0..BRICK_CELLS {
                if b.mask() & (1u64 << local) != 0 {
                    let (i_r, i_z) = b.global_rz(local);
                    act[i_r * n_z + i_z] = true;
                }
            }
        }
        let cut = if g.has_cut_geometry() {
            let mut small = Vec::new();
            let mut counts = vec![1u32; nb * BRICK_CELLS];
            let to_bl = |(r, z): (usize, usize)| {
                let bi = bmap[(r / BRICK) * nbz + z / BRICK].expect("active cell's brick");
                (bi, (r % BRICK) * BRICK + z % BRICK)
            };
            for i_r in 0..n_r {
                for i_z in 0..n_z {
                    if !act[i_r * n_z + i_z] || g.kappa(i_r, i_z) >= KAPPA_SRD {
                        continue;
                    }
                    let cells = neighborhood_cells(g, i_r, i_z)?;
                    let kv: Vec<f64> = cells
                        .iter()
                        .map(|&(r, z)| g.kappa(r, z) * g.cell_volume(r, nt))
                        .collect();
                    let members: Vec<(usize, usize)> = cells.into_iter().map(to_bl).collect();
                    for &(bi, local) in &members[1..] {
                        counts[bi * BRICK_CELLS + local] += 1;
                    }
                    small.push((members, kv));
                }
            }
            Some(CutScratch { small, counts })
        } else {
            None
        };
        Ok(Scratch {
            prim: vec![vec![[0.0; NPRIM]; plane]; nb],
            rate: vec![vec![[0.0; NCOMP]; plane]; nb],
            u0: vec![vec![[0.0; NCOMP]; plane]; nb],
            bmap,
            nbz,
            act,
            cut,
            primed: false,
        })
    }

    /// Build a persistent workspace for repeated stepping: reuse across
    /// steps keeps the previous stage's projected pressures as warm-start
    /// hints for the equilibrium projection (`EosLaw::prim_checked_hinted`)
    /// and retires the per-step scratch allocation churn. One-shot
    /// [`Euler::step`] builds a fresh (cold) one each call.
    pub fn workspace(&self, g: &Grid) -> Result<EulerWorkspace, FlowError> {
        let nt = self.validate(g)?;
        Ok(EulerWorkspace(self.scratch(g, nt)?))
    }

    /// One SSP-RK2 (Heun) step of size `dt` at time `t` (one-shot form —
    /// fresh workspace; see [`Euler::step_ws`]).
    pub fn step(&self, g: &mut Grid, f: &EulerFields, t: f64, dt: f64) -> Result<(), FlowError> {
        let mut ws = self.workspace(g)?;
        self.step_ws(g, f, &mut ws, t, dt)
    }

    /// One SSP-RK2 (Heun) step of size `dt` at time `t`. Fixed traversal and
    /// accumulation order (r-sweep, θ-sweep, z-sweep, geometric, external);
    /// bit-reproducible at any thread count (single-threaded reference).
    /// On cut-geometry worlds every stage ends with a State-Redistribution
    /// pass (Berger & Giuliani apply SRD stagewise inside SSP-RK): the
    /// provisional divide-by-κV update on a small cell is merged into its
    /// neighborhood — conservative by construction, and the reason Δt keeps
    /// the UNCUT CFL.
    pub fn step_ws(
        &self,
        g: &mut Grid,
        f: &EulerFields,
        ws: &mut EulerWorkspace,
        t: f64,
        dt: f64,
    ) -> Result<(), FlowError> {
        let nt = self.validate(g)?;
        let s = &mut ws.0;
        let ids = f.ids();

        // Snapshot U⁰.
        for (bi, u0) in s.u0.iter_mut().enumerate() {
            let b = g.brick(bi);
            for k in 0..NCOMP {
                let src = b.field(ids[k]);
                for (cell, u) in u0.iter_mut().enumerate() {
                    u[k] = src[cell];
                }
            }
        }

        // Stage 1: U¹ = SRD(U⁰ + dt·L(U⁰, t)).
        self.rhs(g, f, nt, s, t)?;
        for bi in 0..g.n_bricks() {
            for (k, &id) in ids.iter().enumerate() {
                let dst = g.brick_field_mut(bi, id);
                for (cell, v) in dst.iter_mut().enumerate() {
                    *v = s.u0[bi][cell][k] + dt * s.rate[bi][cell][k];
                }
            }
        }
        Self::srd(g, f, nt, s);

        // Stage 2: Uⁿ⁺¹ = SRD(½(U⁰ + U¹ + dt·L(U¹, t+dt))).
        self.rhs(g, f, nt, s, t + dt)?;
        for bi in 0..g.n_bricks() {
            for (k, &id) in ids.iter().enumerate() {
                let dst = g.brick_field_mut(bi, id);
                for (cell, v) in dst.iter_mut().enumerate() {
                    *v = 0.5 * (s.u0[bi][cell][k] + *v + dt * s.rate[bi][cell][k]);
                }
            }
        }
        Self::srd(g, f, nt, s);
        Ok(())
    }

    /// The State-Redistribution pass (META-3 `state-redistribution`),
    /// applied to the freshly-written stage state. With neighborhoods
    /// `M_i` (trivial `{i}` for regular cells) and overlap counts
    /// `n_j = #{i : j ∈ M_i}`:
    ///
    /// `Q_i = Σ_{j∈M_i} (κ_j V_j / n_j)·Û_j / Σ_{j∈M_i} (κ_j V_j / n_j)`,
    /// then `U_j = (1/n_j)·Σ_{i : j∈M_i} Q_i`.
    ///
    /// `Σ κVU` is preserved exactly (each j appears in n_j neighborhoods);
    /// only cells touched by a small neighborhood change. Fixed owner and
    /// member order ⇒ deterministic.
    fn srd(g: &mut Grid, f: &EulerFields, nt: u32, s: &Scratch) {
        let Some(cut) = &s.cut else { return };
        if cut.small.is_empty() {
            return;
        }
        let ids = f.ids();
        for j_theta in 0..nt {
            // (bi, local) → (Σ Q contributions, is-small-owner).
            let mut acc: std::collections::BTreeMap<(usize, usize), ([f64; NCOMP], bool)> =
                std::collections::BTreeMap::new();
            for (members, _) in &cut.small {
                acc.entry(members[0]).or_insert(([0.0; NCOMP], false)).1 = true;
            }
            let read = |g: &Grid, bi: usize, local: usize| -> Cons {
                let b = g.brick(bi);
                let idx = b.cell_index(j_theta, local);
                std::array::from_fn(|k| b.field(ids[k])[idx])
            };
            // Neighborhood averages against the CURRENT (provisional) state.
            let mut qs: Vec<[f64; NCOMP]> = Vec::with_capacity(cut.small.len());
            for (members, kv) in &cut.small {
                let mut num = [0.0f64; NCOMP];
                let mut den = 0.0f64;
                for (&(bi, local), &kvj) in members.iter().zip(kv) {
                    let w = kvj / f64::from(cut.counts[bi * BRICK_CELLS + local]);
                    let u = read(g, bi, local);
                    for k in 0..NCOMP {
                        num[k] += w * u[k];
                    }
                    den += w;
                }
                qs.push(std::array::from_fn(|k| num[k] / den));
            }
            // Base term: every affected non-owner keeps its own trivial
            // neighborhood's Q = Û; owners' own Q is the merged one.
            for (members, _) in &cut.small {
                for &(bi, local) in members {
                    acc.entry((bi, local)).or_insert(([0.0; NCOMP], false));
                }
            }
            for ((bi, local), (a, owner)) in &mut acc {
                if !*owner {
                    let u = read(g, *bi, *local);
                    for k in 0..NCOMP {
                        a[k] += u[k];
                    }
                }
            }
            for ((members, _), q) in cut.small.iter().zip(&qs) {
                for &(bi, local) in members {
                    let (a, _) = acc.get_mut(&(bi, local)).expect("member entry");
                    for k in 0..NCOMP {
                        a[k] += q[k];
                    }
                }
            }
            // Write back U_j = acc_j / n_j in fixed (bi, local) order.
            for ((bi, local), (a, _)) in &acc {
                let n = f64::from(cut.counts[*bi * BRICK_CELLS + *local]);
                let idx = g.brick(*bi).cell_index(j_theta, *local);
                for (k, &id) in ids.iter().enumerate() {
                    g.brick_field_mut(*bi, id)[idx] = a[k] / n;
                }
            }
        }
    }

    /// March `n_steps` of size `dt` from `t0`; returns the final time.
    pub fn advance(
        &self,
        g: &mut Grid,
        f: &EulerFields,
        t0: f64,
        dt: f64,
        n_steps: usize,
    ) -> Result<f64, FlowError> {
        let mut t = t0;
        for _ in 0..n_steps {
            self.step(g, f, t, dt)?;
            t += dt;
        }
        Ok(t)
    }

    /// CFL timestep `cfl / max Σ_d (|u_d|+c)/Δ_d` over active cells — a
    /// deterministic fixed rule over the wave speeds (COUP-3 §2 contract;
    /// COUP-3 owns the production Δt schedule). Cold projections (no
    /// workspace); [`Euler::stable_dt_ws`] is the warm-started form.
    pub fn stable_dt(&self, g: &Grid, f: &EulerFields, cfl: f64) -> Result<f64, FlowError> {
        self.stable_dt_inner(g, f, cfl, None)
    }

    /// [`Euler::stable_dt`] with warm-start hints from a stepping
    /// workspace's last primitive fill (an acceleration, never physics).
    pub fn stable_dt_ws(
        &self,
        g: &Grid,
        f: &EulerFields,
        ws: &EulerWorkspace,
        cfl: f64,
    ) -> Result<f64, FlowError> {
        self.stable_dt_inner(g, f, cfl, ws.0.primed.then_some(&ws.0.prim))
    }

    /// Parallel per-brick partial maxima (exact for f64 max — no rounding,
    /// so order-free) with the first bad cell in brick order raised — the
    /// serial sweep's diagnosis.
    fn stable_dt_inner(
        &self,
        g: &Grid,
        f: &EulerFields,
        cfl: f64,
        hints: Option<&Vec<Vec<Prim>>>,
    ) -> Result<f64, FlowError> {
        let nt = self.validate(g)?;
        let ids = f.ids();
        let dtheta = std::f64::consts::TAU / f64::from(nt);
        let (dr, dz) = (g.spec().dr, g.spec().dz);
        let partials: Vec<Result<f64, FlowError>> = (0..g.n_bricks())
            .into_par_iter()
            .map(|bi| {
                let b = g.brick(bi);
                let fields: [&[f64]; NCOMP] = std::array::from_fn(|k| b.field(ids[k]));
                let mut max_sig = 0.0f64;
                for j in 0..nt {
                    for local in 0..BRICK_CELLS {
                        if b.mask() & (1u64 << local) == 0 {
                            continue;
                        }
                        let idx = j as usize * BRICK_CELLS + local;
                        let u: Cons = std::array::from_fn(|k| fields[k][idx]);
                        let (i_r, i_z) = b.global_rz(local);
                        let hint = hints.map(|h| h[bi][idx][4]);
                        match self.eos.prim_checked_hinted(&u, hint) {
                            Ok(w) => {
                                let c = self.eos.sound_speed_w(&w);
                                let mut sig = (w[1].abs() + c) / dr + (w[3].abs() + c) / dz;
                                if nt > 1 {
                                    sig += (w[2].abs() + c) / (g.r_center(i_r) * dtheta);
                                }
                                max_sig = max_sig.max(sig);
                            }
                            Err(what) => {
                                return Err(FlowError::NonPhysicalState {
                                    i_r,
                                    i_z,
                                    i_theta: j,
                                    what,
                                });
                            }
                        }
                    }
                }
                Ok(max_sig)
            })
            .collect();
        let mut max_sig = 0.0f64;
        for p in partials {
            max_sig = max_sig.max(p?);
        }
        Ok(cfl / max_sig)
    }

    /// L(U): the flux-divergence + source contribution (SOLV-1 §2 contract),
    /// into `s.rate`.
    fn rhs(
        &self,
        g: &Grid,
        f: &EulerFields,
        nt: u32,
        s: &mut Scratch,
        t: f64,
    ) -> Result<(), FlowError> {
        self.fill_prims(g, f, nt, s)?;
        for rate in &mut s.rate {
            for cell in rate.iter_mut() {
                *cell = [0.0; NCOMP];
            }
        }
        self.sweep_r(g, nt, s, t)?;
        self.sweep_theta(g, nt, s);
        self.sweep_z(g, nt, s, t)?;
        self.add_sources(g, nt, s, t);
        Ok(())
    }

    /// Parallel over bricks (each brick's prim block is its own task —
    /// FND-2 §3.7: ownership partition, so results are thread-count-
    /// independent). Errors are gathered per brick and the FIRST in brick
    /// order is raised — the same diagnosis the serial sweep chose.
    fn fill_prims(
        &self,
        g: &Grid,
        f: &EulerFields,
        nt: u32,
        s: &mut Scratch,
    ) -> Result<(), FlowError> {
        let ids = f.ids();
        let primed = s.primed;
        let results: Vec<Result<(), FlowError>> = s
            .prim
            .par_iter_mut()
            .enumerate()
            .map(|(bi, prim)| {
                let b = g.brick(bi);
                let fields: [&[f64]; NCOMP] = std::array::from_fn(|k| b.field(ids[k]));
                for j in 0..nt {
                    for local in 0..BRICK_CELLS {
                        if b.mask() & (1u64 << local) == 0 {
                            continue;
                        }
                        let idx = j as usize * BRICK_CELLS + local;
                        let u: Cons = std::array::from_fn(|k| fields[k][idx]);
                        // Warm-start hint: the previous stage's projected p
                        // for this cell (an acceleration, never physics).
                        let hint = if primed { Some(prim[idx][4]) } else { None };
                        match self.eos.prim_checked_hinted(&u, hint) {
                            Ok(w) => prim[idx] = w,
                            Err(what) => {
                                let (i_r, i_z) = b.global_rz(local);
                                return Err(FlowError::NonPhysicalState {
                                    i_r,
                                    i_z,
                                    i_theta: j,
                                    what,
                                });
                            }
                        }
                    }
                }
                Ok(())
            })
            .collect();
        for r in results {
            r?;
        }
        s.primed = true;
        Ok(())
    }

    /// Solid-wall mirror ghosts on the low side of a run (`normal` is the
    /// velocity slot to negate) — used for reflecting domain faces AND for
    /// interior stair-step wall faces when no true wall normal is declared.
    fn mirror_low(w: &mut [Prim], n: usize, normal: usize) {
        for k in 1..=NGHOST {
            let mut m = w[NGHOST + (k - 1).min(n - 1)];
            m[normal] = -m[normal];
            w[NGHOST - k] = m;
        }
    }

    fn mirror_high(w: &mut [Prim], n: usize, normal: usize) {
        for k in 1..=NGHOST {
            let mut m = w[NGHOST + n.saturating_sub(k).min(n - 1)];
            m[normal] = -m[normal];
            w[NGHOST + n - 1 + k] = m;
        }
    }

    /// Reflect a primitive state's meridional velocity about the true wall
    /// tangent: `v ← v − 2(v·n̂)n̂` on (u_r, u_z); u_θ is tangential to any
    /// axisymmetric wall and untouched.
    #[inline]
    fn slip_reflect(mut m: Prim, n_hat: (f64, f64)) -> Prim {
        let vn = m[I_MR] * n_hat.0 + m[I_MZ] * n_hat.1;
        m[I_MR] -= 2.0 * vn * n_hat.0;
        m[I_MZ] -= 2.0 * vn * n_hat.1;
        m
    }

    /// Interior-wall ghosts on the low side of a run: slip-mirror about the
    /// true wall normal at `(r, z)` when declared, else the grid-aligned
    /// mirror in slot `normal`.
    fn wall_ghosts_low(&self, w: &mut [Prim], n: usize, normal: usize, r: f64, z: f64) {
        let slip = normal != I_MZ || self.slip_wall_z_faces;
        match self.wall_normal {
            Some(nf) if slip => {
                let n_hat = nf(r, z);
                for k in 1..=NGHOST {
                    w[NGHOST - k] = Self::slip_reflect(w[NGHOST + (k - 1).min(n - 1)], n_hat);
                }
            }
            _ => Self::mirror_low(w, n, normal),
        }
    }

    fn wall_ghosts_high(&self, w: &mut [Prim], n: usize, normal: usize, r: f64, z: f64) {
        let slip = normal != I_MZ || self.slip_wall_z_faces;
        match self.wall_normal {
            Some(nf) if slip => {
                let n_hat = nf(r, z);
                for k in 1..=NGHOST {
                    w[NGHOST + n - 1 + k] =
                        Self::slip_reflect(w[NGHOST + n.saturating_sub(k).min(n - 1)], n_hat);
                }
            }
            _ => Self::mirror_high(w, n, normal),
        }
    }

    /// Fill the low-side ghosts of a run from a domain BC. `normal` is the
    /// velocity slot to mirror; `pos(k)` gives the k-th ghost centroid.
    fn fill_ghosts_low(
        &self,
        w: &mut [Prim],
        n: usize,
        bc: &FlowBc<'_>,
        normal: usize,
        pos: impl Fn(usize) -> (f64, f64, f64),
        t: f64,
    ) -> Result<(), FlowError> {
        for k in 1..=NGHOST {
            w[NGHOST - k] = match bc {
                FlowBc::Reflecting => {
                    let mut m = w[NGHOST + (k - 1).min(n - 1)];
                    m[normal] = -m[normal];
                    m
                }
                FlowBc::Transmissive => w[NGHOST],
                FlowBc::Prescribed(f) => {
                    let (r, th, z) = pos(k);
                    f(r, th, z, t)
                }
                FlowBc::StagnationInflow { p0, rho0, c_frac } => {
                    self.eos
                        .stagnation_ghost(*p0, *rho0, *c_frac, w[NGHOST][normal], normal)?
                }
                FlowBc::MassFlowInflow {
                    mdot_per_area,
                    h_total,
                    c_frac,
                } => self.eos.mass_flow_inflow_ghost(
                    *mdot_per_area,
                    *h_total,
                    *c_frac,
                    w[NGHOST][4],
                    1.0,
                    normal,
                )?,
                FlowBc::PressureOutflow(p_amb) => {
                    let m = w[NGHOST];
                    // Outflow through a LOW face is velocity toward −normal.
                    let out_mach = (-m[normal]) / self.eos.sound_speed_w(&m);
                    if out_mach >= 1.0 {
                        m
                    } else {
                        let mut g = m;
                        g[4] = p_amb(t);
                        g
                    }
                }
            };
        }
        Ok(())
    }

    fn fill_ghosts_high(
        &self,
        w: &mut [Prim],
        n: usize,
        bc: &FlowBc<'_>,
        normal: usize,
        pos: impl Fn(usize) -> (f64, f64, f64),
        t: f64,
    ) -> Result<(), FlowError> {
        for k in 1..=NGHOST {
            w[NGHOST + n - 1 + k] = match bc {
                FlowBc::Reflecting => {
                    let mut m = w[NGHOST + n.saturating_sub(k).min(n - 1)];
                    m[normal] = -m[normal];
                    m
                }
                FlowBc::Transmissive => w[NGHOST + n - 1],
                FlowBc::Prescribed(f) => {
                    let (r, th, z) = pos(k);
                    f(r, th, z, t)
                }
                FlowBc::StagnationInflow { p0, rho0, c_frac } => self.eos.stagnation_ghost(
                    *p0,
                    *rho0,
                    *c_frac,
                    w[NGHOST + n - 1][normal],
                    normal,
                )?,
                FlowBc::MassFlowInflow {
                    mdot_per_area,
                    h_total,
                    c_frac,
                } => self.eos.mass_flow_inflow_ghost(
                    *mdot_per_area,
                    *h_total,
                    *c_frac,
                    w[NGHOST + n - 1][4],
                    -1.0,
                    normal,
                )?,
                FlowBc::PressureOutflow(p_amb) => {
                    let m = w[NGHOST + n - 1];
                    let out_mach = m[normal] / self.eos.sound_speed_w(&m);
                    if out_mach >= 1.0 {
                        m
                    } else {
                        let mut g = m;
                        g[4] = p_amb(t);
                        g
                    }
                }
            };
        }
        Ok(())
    }

    /// Radial sweep. Face f (0..=n_r) sits at radius `r_min + f·dr`; the
    /// axis face (f = 0 when r_min = 0) has zero area and drops out
    /// geometrically (FND-2 §3.2). Each pencil line is decomposed into
    /// maximal active runs (`Scratch::act`): an interior run boundary is a
    /// stair-step wall face (reflecting mirror — the binary-aperture
    /// degenerate wall; FND-3 cut cells supersede), a domain-edge boundary
    /// takes the configured BC.
    fn sweep_r(&self, g: &Grid, nt: u32, s: &mut Scratch, t: f64) -> Result<(), FlowError> {
        let spec = g.spec();
        let (n, n_z) = (spec.n_r, spec.n_z);
        let (r0, dr) = (spec.r_min, spec.dr);
        let on_axis = r0 == 0.0;
        let area: Vec<f64> = (0..=n)
            .map(|fi| {
                if fi < n {
                    g.face_area_r(fi, false, nt)
                } else {
                    g.face_area_r(n - 1, true, nt)
                }
            })
            .collect();
        let vol: Vec<f64> = (0..n).map(|i| g.cell_volume(i, nt)).collect();

        // Parallel by brick z-row: a pencil at i_z only touches bricks with
        // bz = i_z/BRICK, so rows are an ownership partition of the rate
        // buffers (FND-2 §3.7: per-cell accumulation order identical at any
        // thread count). Rows report errors in order; the first is raised —
        // the same cell the serial sweep would have named.
        let nbr = n.div_ceil(BRICK);
        let nbz = s.nbz;
        let mut rows: Vec<Vec<(usize, &mut Vec<Cons>)>> = (0..nbz).map(|_| Vec::new()).collect();
        for (bi, rv) in s.rate.iter_mut().enumerate() {
            rows[g.brick(bi).bz() as usize].push((bi, rv));
        }
        let (prim, act, bmap) = (&s.prim, &s.act, &s.bmap);
        let results: Vec<Result<(), FlowError>> = rows
            .into_par_iter()
            .enumerate()
            .map(|(bz, mut row)| {
                let mut slot: Vec<Option<usize>> = vec![None; nbr];
                for (k, (bi, _)) in row.iter().enumerate() {
                    slot[g.brick(*bi).br() as usize] = Some(k);
                }
                let mut w = vec![[0.0f64; NPRIM]; n + 2 * NGHOST];
                let mut fl = vec![[0.0f64; NPRIM]; n + 1];
                let mut fr = vec![[0.0f64; NPRIM]; n + 1];
                let mut af = vec![[0.0f64; NCOMP]; n + 1];
                let mut ap = vec![1.0f64; n + 1];
                let mut kap = vec![1.0f64; n];
                for i_z in bz * BRICK..((bz + 1) * BRICK).min(n_z) {
                    let lz = i_z % BRICK;
                    let z = g.z_center(i_z);
                    for j in 0..nt {
                        let theta = Grid::theta_center(j, nt);
                        let mut i = 0usize;
                        while i < n {
                            if !act[i * n_z + i_z] {
                                i += 1;
                                continue;
                            }
                            let start = i;
                            while i < n && act[i * n_z + i_z] {
                                i += 1;
                            }
                            let len = i - start;
                            for (q, ii) in (start..start + len).enumerate() {
                                let bi =
                                    bmap[(ii / BRICK) * nbz + bz].expect("active cell's brick");
                                let local = (ii % BRICK) * BRICK + lz;
                                w[NGHOST + q] = prim[bi][j as usize * BRICK_CELLS + local];
                                let b = g.brick(bi);
                                kap[q] = b.kappa_rz(local);
                                ap[q] = b.aperture_rz(FaceDir::RMinus, local);
                                if q + 1 == len {
                                    ap[len] = b.aperture_rz(FaceDir::RPlus, local);
                                }
                            }
                            if start == 0 && on_axis {
                                // Through-axis mirror: ê_r and ê_θ both flip
                                // (the N_θ=1 degenerate parity pairing).
                                for k in 1..=NGHOST {
                                    let mut m = w[NGHOST + (k - 1).min(len - 1)];
                                    m[I_MR] = -m[I_MR];
                                    m[I_MT] = -m[I_MT];
                                    w[NGHOST - k] = m;
                                }
                            } else if start == 0 {
                                self.fill_ghosts_low(
                                    &mut w,
                                    len,
                                    &self.bcs.r_inner,
                                    I_MR,
                                    |k| (r0 - (k as f64 - 0.5) * dr, theta, z),
                                    t,
                                )?;
                            } else {
                                self.wall_ghosts_low(&mut w, len, I_MR, r0 + start as f64 * dr, z);
                            }
                            if start + len == n {
                                self.fill_ghosts_high(
                                    &mut w,
                                    len,
                                    &self.bcs.r_outer,
                                    I_MR,
                                    |k| (r0 + (n as f64 + k as f64 - 0.5) * dr, theta, z),
                                    t,
                                )?;
                            } else {
                                self.wall_ghosts_high(
                                    &mut w,
                                    len,
                                    I_MR,
                                    r0 + (start + len) as f64 * dr,
                                    z,
                                );
                            }
                            ppm_faces(&w[..len + 2 * NGHOST], len, &mut fl, &mut fr);
                            for fi in 0..=len {
                                let flux = hllc_flux(&fl[fi], &fr[fi], I_MR, &self.eos);
                                // Aperture-weighted open area (FND-3 §3.3);
                                // ap = 1.0 exactly on full-box worlds, so
                                // `(A·1.0)·F ≡ A·F` bitwise — no mode branch.
                                let aa = area[start + fi] * ap[fi];
                                for k in 0..NCOMP {
                                    af[fi][k] = aa * flux[k];
                                }
                            }
                            for q in 0..len {
                                let ii = start + q;
                                let k_slot = slot[ii / BRICK].expect("active brick in row");
                                let idx = j as usize * BRICK_CELLS + (ii % BRICK) * BRICK + lz;
                                let rate = &mut row[k_slot].1[idx];
                                for k in 0..NCOMP {
                                    rate[k] += (af[q][k] - af[q + 1][k]) / (kap[q] * vol[ii]);
                                }
                            }
                        }
                    }
                }
                Ok(())
            })
            .collect();
        for r in results {
            r?;
        }
        Ok(())
    }

    /// Azimuthal sweep: periodic ring pencils, each inside one brick. At
    /// N_θ = 1 both faces carry the same computed flux and cancel exactly —
    /// geometry, not a code branch.
    fn sweep_theta(&self, g: &Grid, nt: u32, s: &mut Scratch) {
        let n = nt as usize;
        let a_th = g.face_area_theta();
        let mut w = vec![[0.0f64; NPRIM]; n + 2 * NGHOST];
        let mut fl = vec![[0.0f64; NPRIM]; n + 1];
        let mut fr = vec![[0.0f64; NPRIM]; n + 1];
        let mut af = vec![[0.0f64; NCOMP]; n + 1];

        for bi in 0..g.n_bricks() {
            let mask = g.brick(bi).mask();
            for local in 0..BRICK_CELLS {
                if mask & (1u64 << local) == 0 {
                    continue;
                }
                let (i_r, _) = g.brick(bi).global_rz(local);
                // One per-ring metric ratio `A_θ/V` (= 1/(r̄·Δθ) exactly),
                // applied to the flux difference — the same conditioning
                // rule as the z sweep.
                let inv = a_th / g.cell_volume(i_r, nt);
                // Periodic gather with wrapped ghosts: pencil index pi maps
                // to ring cell (pi − NGHOST) mod n, offset by NGHOST·n so
                // the subtraction cannot underflow at any n ≥ 1.
                for (pi, cell) in w.iter_mut().enumerate() {
                    let jj = (pi + NGHOST * n - NGHOST) % n;
                    *cell = s.prim[bi][jj * BRICK_CELLS + local];
                }
                ppm_faces(&w, n, &mut fl, &mut fr);
                for fi in 0..=n {
                    af[fi] = hllc_flux(&fl[fi], &fr[fi], I_MT, &self.eos);
                }
                for j in 0..n {
                    let rate = &mut s.rate[bi][j * BRICK_CELLS + local];
                    for k in 0..NCOMP {
                        rate[k] += (af[j][k] - af[j + 1][k]) * inv;
                    }
                }
            }
        }
    }

    /// Axial sweep (the Sod direction). The metric identity `A_z/V = 1/dz`
    /// holds exactly for every ring, so the divergence is accumulated as
    /// `(F_i − F_{i+1})/dz` — the well-conditioned form (META-1 §3). This is
    /// also what keeps a radially-uniform state radially uniform **bitwise**:
    /// with per-ring `A_z·F` products, identical physics on different rings
    /// rounds differently by an ulp (found by the Station-1 certificate).
    fn sweep_z(&self, g: &Grid, nt: u32, s: &mut Scratch, t: f64) -> Result<(), FlowError> {
        let spec = g.spec();
        let (n, n_r) = (spec.n_z, spec.n_r);
        let (z0, dz) = (spec.z_min, spec.dz);
        let inv_dz = 1.0 / dz;
        // Parallel by brick column (see sweep_r's row rationale, transposed:
        // a pencil at i_r touches only bricks with br = i_r/BRICK).
        let nbz = s.nbz;
        let nbr = n_r.div_ceil(BRICK);
        let mut cols: Vec<Vec<(usize, &mut Vec<Cons>)>> = (0..nbr).map(|_| Vec::new()).collect();
        for (bi, rv) in s.rate.iter_mut().enumerate() {
            cols[g.brick(bi).br() as usize].push((bi, rv));
        }
        let (prim, act, bmap) = (&s.prim, &s.act, &s.bmap);
        let results: Vec<Result<(), FlowError>> = cols
            .into_par_iter()
            .enumerate()
            .map(|(br, mut col)| {
                let mut slot: Vec<Option<usize>> = vec![None; nbz];
                for (k, (bi, _)) in col.iter().enumerate() {
                    slot[g.brick(*bi).bz() as usize] = Some(k);
                }
                let mut w = vec![[0.0f64; NPRIM]; n + 2 * NGHOST];
                let mut fl = vec![[0.0f64; NPRIM]; n + 1];
                let mut fr = vec![[0.0f64; NPRIM]; n + 1];
                let mut af = vec![[0.0f64; NCOMP]; n + 1];
                let mut ap = vec![1.0f64; n + 1];
                let mut kap = vec![1.0f64; n];
                for i_r in br * BRICK..((br + 1) * BRICK).min(n_r) {
                    let lr = i_r % BRICK;
                    let r = g.r_center(i_r);
                    for j in 0..nt {
                        let theta = Grid::theta_center(j, nt);
                        let mut i = 0usize;
                        while i < n {
                            if !act[i_r * n + i] {
                                i += 1;
                                continue;
                            }
                            let start = i;
                            while i < n && act[i_r * n + i] {
                                i += 1;
                            }
                            let len = i - start;
                            for (q, ii) in (start..start + len).enumerate() {
                                let bi = bmap[br * nbz + ii / BRICK].expect("active cell's brick");
                                let local = lr * BRICK + (ii % BRICK);
                                w[NGHOST + q] = prim[bi][j as usize * BRICK_CELLS + local];
                                let b = g.brick(bi);
                                kap[q] = b.kappa_rz(local);
                                ap[q] = b.aperture_rz(FaceDir::ZMinus, local);
                                if q + 1 == len {
                                    ap[len] = b.aperture_rz(FaceDir::ZPlus, local);
                                }
                            }
                            if start == 0 {
                                self.fill_ghosts_low(
                                    &mut w,
                                    len,
                                    &self.bcs.z_lo,
                                    I_MZ,
                                    |k| (r, theta, z0 - (k as f64 - 0.5) * dz),
                                    t,
                                )?;
                            } else {
                                self.wall_ghosts_low(&mut w, len, I_MZ, r, z0 + start as f64 * dz);
                            }
                            if start + len == n {
                                self.fill_ghosts_high(
                                    &mut w,
                                    len,
                                    &self.bcs.z_hi,
                                    I_MZ,
                                    |k| (r, theta, z0 + (n as f64 + k as f64 - 0.5) * dz),
                                    t,
                                )?;
                            } else {
                                self.wall_ghosts_high(
                                    &mut w,
                                    len,
                                    I_MZ,
                                    r,
                                    z0 + (start + len) as f64 * dz,
                                );
                            }
                            ppm_faces(&w[..len + 2 * NGHOST], len, &mut fl, &mut fr);
                            for fi in 0..=len {
                                let flux = hllc_flux(&fl[fi], &fr[fi], I_MZ, &self.eos);
                                // ap = 1.0 exactly on full-box worlds
                                // (`F·1.0 ≡ F`).
                                for k in 0..NCOMP {
                                    af[fi][k] = flux[k] * ap[fi];
                                }
                            }
                            for q in 0..len {
                                let ii = start + q;
                                let k_slot = slot[ii / BRICK].expect("active brick in column");
                                let idx = j as usize * BRICK_CELLS + lr * BRICK + (ii % BRICK);
                                let rate = &mut col[k_slot].1[idx];
                                for k in 0..NCOMP {
                                    rate[k] += (af[q][k] - af[q + 1][k]) * inv_dz / kap[q];
                                }
                            }
                        }
                    }
                }
                Ok(())
            })
            .collect();
        for r in results {
            r?;
        }
        Ok(())
    }

    /// Geometric sources (SOLV-1 §3.3) + the external source intake.
    /// Well-balance: the radial-momentum pressure source is written with the
    /// same `A·p` products the face fluxes use, so for a uniform state at
    /// rest it is the exact negation of the pressure-flux difference —
    /// preserved to round-off, no equilibrium drift. The centrifugal
    /// `ρu_θ²` term and the swirl-advection term `−ρu_ru_θ` use the
    /// metric-consistent `1/r̄ = (A_out−A_in)/V`.
    fn add_sources(&self, g: &Grid, nt: u32, s: &mut Scratch, t: f64) {
        let prim = &s.prim;
        s.rate.par_iter_mut().enumerate().for_each(|(bi, rate_v)| {
            let mask = g.brick(bi).mask();
            for local in 0..BRICK_CELLS {
                if mask & (1u64 << local) == 0 {
                    continue;
                }
                let (i_r, i_z) = g.brick(bi).global_rz(local);
                let a_in = g.face_area_r(i_r, false, nt);
                let a_out = g.face_area_r(i_r, true, nt);
                let vol = g.cell_volume(i_r, nt);
                let geo = (a_out - a_in) / vol;
                let (r, z) = (g.r_center(i_r), g.z_center(i_z));
                // Embedded-interface pressure closure (cut cells only —
                // gated on geometry presence so full-box worlds stay
                // bit-identical; `+0.0` could flip a −0.0 rate bit).
                let wall = if g.brick(bi).has_geom() {
                    let (w_r, w_z) = g.wall_closure(i_r, i_z, nt);
                    let inv_kv = 1.0 / (g.brick(bi).kappa_rz(local) * vol);
                    Some((w_r * inv_kv, w_z * inv_kv))
                } else {
                    None
                };
                for j in 0..nt {
                    let idx = j as usize * BRICK_CELLS + local;
                    let wc = prim[bi][idx];
                    let (rho, ur, ut, p) = (wc[I_RHO], wc[1], wc[2], wc[4]);
                    let rate = &mut rate_v[idx];
                    rate[I_MR] += (a_out * p - a_in * p) / vol + rho * ut * ut * geo;
                    rate[I_MT] -= rho * ur * ut * geo;
                    if let Some((wr_kv, wz_kv)) = wall {
                        rate[I_MR] += p * wr_kv;
                        rate[I_MZ] += p * wz_kv;
                    }
                    let src = (self.source)(r, Grid::theta_center(j, nt), z, t);
                    for k in 0..NCOMP {
                        rate[k] += src[k];
                    }
                }
            }
        });
    }
}
