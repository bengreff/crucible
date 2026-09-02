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
//! Time integration (S2): this operator is COUP-3 §3.1's **explicit
//! hyperbolic class `A`** — it evaluates `L(U)` ([`Euler::eval_rhs`]) with
//! a per-evaluation COUP-2 port/source ledger; the one production advance
//! is the SDC-IMEX step (`crate::sdc`), which owns node weights, state
//! composition, the stagewise SRD passes, and the Δt rule. The session-7
//! MOL SSP-RK2 scaffolding is retired (superseded, not extended — Castro's
//! SDC path uses exactly this MOL reconstruction, hence no characteristic
//! tracing).
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
//! Azimuthal capability (S8): the sweeps run at **each brick's own N_θ**.
//! `r_min = 0` with `N_θ > 1` uses the FND-2 §3.2 **θ↔θ+π parity-pair
//! gather** for the cross-axis ghosts (`u_r` and `u_θ` negated — the basis
//! flip; at N_θ = 1 the partner is the cell itself, so the certified
//! axisymmetric corner is arithmetically identical). Mixed per-brick N_θ
//! marches by the FND-2 §3.4 **ring-interface rule**: r/z pencils decompose
//! into maximal uniform-N_θ segments; at an N_θ jump the **fine side owns
//! the interface flux** (its own reconstruction against piecewise-constant-
//! prolonged coarse ghosts) and the coarse cell applies the area-weighted
//! aggregate of the same numbers — bit-exact telescoping, interior to the
//! COUP-2 ledger. Constraints (validated loudly): 2:1 ladder adjacency
//! between face-adjacent bricks; ≥ NGHOST cells of uniform N_θ on each side
//! of a jump (proper nesting); combustion at mixed N_θ refuses → plan S11.
//! Cut geometry (S9, FND-3 §3.3 3-D apertures): legal at **uniform** N_θ —
//! the sweeps and SRD consume per-θ-sector κ and six face apertures
//! (`kappa_cell`/`aperture_cell`; at N_θ = 1 plane 0 keeps every certified
//! contour world bit-identical); mixed-N_θ cut worlds and combustion on a
//! cut world at N_θ > 1 refuse → plan S11. The Δt rule carries
//! the per-brick θ-CFL member and, when combustion is scheduled, the
//! SOLV-4 §3.6 front-carrier signal `σ_front` with its `S_T_MACH_LIMIT`
//! model-form scale-separation refusal (COUP-3 §3.4, S8).

mod blend_eos;
mod combustion;
mod exact;
mod hllc;
mod recon;
mod table_eos;

pub use blend_eos::{BurnBlendEos, EPS_B_PURE_BURNT, EPS_B_PURE_UNBURNT};
pub use combustion::{
    BURN_COMPLETE, C_NAGUMO_SLOPE, CombMarshal, Combustion, EPS_BURN_BOUND, EPS_IGNITED,
    IgnitionColumns, S_T_MACH_LIMIT, THETA_CELLS, consumption_rate, reacting_measure,
};
pub use exact::{RiemannSide, RiemannSolution, solve as solve_riemann};
pub use hllc::{hllc_flux, physical_flux};
use recon::{NGHOST, ppm_faces};
pub use table_eos::{EPS_P_PROJECTION, N_INFLOW_ITER, N_P_ITER_MAX, TableEos, TableEosMarshal};

use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, FieldId, Grid, GridError};
use rayon::prelude::*;

/// SOLV-1 §3.6 / FND-3 §3.4 — the State Redistribution small-cell
/// threshold (Berger & Giuliani 2020, META-3 `state-redistribution`): a
/// cut cell with gas volume fraction κ below this merges its provisional
/// update into a flow-connected neighborhood whose summed κ reaches it.
/// The paper's standard half-cell target; named constant, never tuned
/// per geometry.
pub const KAPPA_SRD: f64 = 0.5;

/// A merge neighborhood in its (r,z)-keyed view: member cells as global
/// (i_r, i_z) with their gas volumes κ·V (per θ-plane) — the SRD weights.
pub type SrdNeighborhood = Vec<((usize, usize), f64)>;

/// The deterministic SRD merging neighborhood of cell (i_r, i_z) in its
/// **(r,z)-keyed sector-0 view**: itself, then flow-connected face
/// neighbors (shared aperture > 0) added in descending-κ order (ties: the
/// fixed candidate order, see [`neighborhood_cells`]) until
/// `Σκ ≥ KAPPA_SRD`. `None` for a regular cell (κ ≥ threshold). The
/// under-resolved case — no reachable neighborhood — refuses loudly.
/// Public because the engine's wall-exchange debit must deposit into the
/// same merged control volume SRD stabilizes (single owner of the rule).
/// Returned weights are the members' gas volumes κ·V (per θ-plane).
/// Its consumers (the wall patches) march N_θ = 1 worlds — a neighborhood
/// that recruits a θ-neighbor (possible only at N_θ > 1, S9) cannot be
/// expressed in this view and refuses rather than mis-keying.
pub fn srd_neighborhood(
    g: &Grid,
    i_r: usize,
    i_z: usize,
) -> Result<Option<SrdNeighborhood>, FlowError> {
    let nt = g.brick(0).n_theta();
    if g.kappa(i_r, i_z) >= KAPPA_SRD {
        return Ok(None);
    }
    let members = neighborhood_cells(g, i_r, i_z, 0)?;
    if members.iter().any(|&(_, _, j)| j != 0) {
        return Err(FlowError::ThetaCutUnsupported {
            what: "the (r,z)-keyed SRD debit neighborhood (a θ-recruited member; \
                   per-θ wall patches)",
        });
    }
    Ok(Some(
        members
            .into_iter()
            .map(|(r, z, _)| ((r, z), g.kappa(r, z) * g.cell_volume(r, nt)))
            .collect(),
    ))
}

/// θ-mapping of a neighbor segment's cell state for an N_θ-interface ghost
/// (FND-2 §3.4, S8): piecewise-constant prolongation from a coarser
/// neighbor (`j >> 1`), equal-volume pair-mean restriction from a finer
/// one. `nt_nbr ∈ {nts/2, 2·nts}` by the 2:1 adjacency rule (equal
/// resolutions never form an interface — segments are maximal).
#[inline]
fn theta_mapped(prim: &[Prim], local: usize, j: u32, nts: u32, nt_nbr: u32) -> Prim {
    if nt_nbr < nts {
        prim[(j >> 1) as usize * BRICK_CELLS + local]
    } else {
        let a = prim[2 * j as usize * BRICK_CELLS + local];
        let b = prim[(2 * j as usize + 1) * BRICK_CELLS + local];
        std::array::from_fn(|k| 0.5 * (a[k] + b[k]))
    }
}

/// One uniform-N_θ span of a pencil run (S8): its stored per-θ face values
/// (`af[j·(len+1) + fi]`) and per-sector cell κ (`kap[j·len + q]`, S9 —
/// cut geometry is per-θ-sector), for the two-pass sweeps.
struct SweepSeg {
    start: usize,
    len: usize,
    nts: u32,
    af: Vec<Cons>,
    kap: Vec<f64>,
}

impl SweepSeg {
    fn new(start: usize, len: usize, nts: u32) -> Self {
        SweepSeg {
            start,
            len,
            nts,
            af: vec![[0.0; NCOMP]; nts as usize * (len + 1)],
            kap: vec![1.0; nts as usize * len],
        }
    }
}

/// FND-2 §3.4 interface-flux fix-up: at each N_θ jump, replace the COARSE
/// side's boundary face entry with `agg(child_a, child_b)` of the FINE
/// side's two children — one computed number on both sides, so the
/// interface telescopes exactly. `agg` is the plain sum where `af` carries
/// area (the r sweep) and the ½-weighted sum in metric-ratio form (the z
/// sweep, where `A_zf/A_zc = 1/2` exactly). Adjacent segments differ by
/// exactly one ladder factor (`Euler::validate`).
fn reflux_fixup(segs: &mut [SweepSeg], agg: impl Fn(&Cons, &Cons) -> Cons) {
    for si in 1..segs.len() {
        let (lhs, rhs) = segs.split_at_mut(si);
        let left = &mut lhs[si - 1];
        let right = &mut rhs[0];
        let (ll, rl) = (left.len, right.len);
        if left.nts == 2 * right.nts {
            // Left finer: its high-face fluxes aggregate into the right
            // (coarse) segment's low entries.
            for jc in 0..right.nts {
                let a = left.af[(2 * jc) as usize * (ll + 1) + ll];
                let b = left.af[(2 * jc + 1) as usize * (ll + 1) + ll];
                right.af[jc as usize * (rl + 1)] = agg(&a, &b);
            }
        } else {
            // Right finer: its low-face fluxes aggregate into the left
            // (coarse) segment's high entries.
            for jc in 0..left.nts {
                let a = right.af[(2 * jc) as usize * (rl + 1)];
                let b = right.af[(2 * jc + 1) as usize * (rl + 1)];
                left.af[jc as usize * (ll + 1) + ll] = agg(&a, &b);
            }
        }
    }
}

/// The one neighborhood-construction rule (see [`srd_neighborhood`]),
/// sector-resolved (S9): the neighborhood of θ-sector cell
/// (i_r, i_theta, i_z) is itself, then flow-connected face neighbors
/// (shared aperture > 0) added in descending-κ order until
/// `Σκ ≥ KAPPA_SRD`. Ties break on the fixed candidate order
/// **r−, r+, θ−, θ+, z−, z+** (the S9 extension of the S2 FaceDir tie
/// rule; documented order, stable sort). θ-candidates are the two ring
/// neighbors (j ∓ 1 mod N_θ) of the SAME (r,z) cell, offered only at
/// N_θ > 1 — at N_θ = 1 the ring neighbor is the cell itself, so it is
/// excluded STRUCTURALLY and the N_θ = 1 neighborhoods are bit-identical
/// to the pre-S9 rule.
fn neighborhood_cells(
    g: &Grid,
    i_r: usize,
    i_z: usize,
    i_theta: u32,
) -> Result<Vec<(usize, usize, u32)>, FlowError> {
    let nt = g.brick(0).n_theta();
    let mut members = vec![(i_r, i_z, i_theta)];
    let mut sum = g.kappa_at(i_r, i_theta, i_z);
    if sum >= KAPPA_SRD {
        return Ok(members);
    }
    let mut cand: Vec<((usize, usize, u32), f64)> = Vec::new();
    let push_rz = |dir: FaceDir, nbr: Option<(usize, usize)>, cand: &mut Vec<_>| {
        if let Some((nr, nz)) = nbr
            && g.aperture_at(i_r, i_theta, i_z, dir) > 0.0
            && g.is_active(nr, nz)
        {
            cand.push(((nr, nz, i_theta), g.kappa_at(nr, i_theta, nz)));
        }
    };
    push_rz(
        FaceDir::RMinus,
        i_r.checked_sub(1).map(|r| (r, i_z)),
        &mut cand,
    );
    push_rz(FaceDir::RPlus, Some((i_r + 1, i_z)), &mut cand);
    if nt > 1 {
        let ring = [
            (FaceDir::ThetaMinus, (i_theta + nt - 1) % nt),
            (FaceDir::ThetaPlus, (i_theta + 1) % nt),
        ];
        for (dir, jn) in ring {
            if g.aperture_at(i_r, i_theta, i_z, dir) > 0.0 {
                cand.push(((i_r, i_z, jn), g.kappa_at(i_r, jn, i_z)));
            }
        }
    }
    push_rz(
        FaceDir::ZMinus,
        i_z.checked_sub(1).map(|z| (i_r, z)),
        &mut cand,
    );
    push_rz(FaceDir::ZPlus, Some((i_r, i_z + 1)), &mut cand);
    // Stable sort: descending κ, the fixed candidate order breaking ties.
    cand.sort_by(|a, b| b.1.partial_cmp(&a.1).expect("finite κ"));
    for ((nr, nz, nj), k) in cand {
        members.push((nr, nz, nj));
        sum += k;
        if sum >= KAPPA_SRD {
            return Ok(members);
        }
    }
    Err(FlowError::CutCellUnmergeable { i_r, i_z, i_theta })
}

/// Components of `U` (SOLV-1 §3.1) and of the primitive view
/// `W = (ρ, u_r, u_θ, u_z, p, C, b | e, Γ₁)`. Slots 1–3 are the velocity/
/// momentum directions, so a face's normal is named by its slot index.
/// `W` carries `NPRIM − NCOMP = 2` auxiliary EOS slots (specific internal
/// energy, effective Γ₁) with no conserved counterpart — filled by
/// `EosLaw::prim_checked`, reconstructed componentwise, read only by the
/// general-EOS face closures (`GammaLaw` ignores them).
///
/// **Slot 6 (`I_RB`) is the S6 burn-progress fixed `+1` widening** (SOLV-1
/// §3.4, SOLV-4 §3.6): conserved `ρb`, `b ∈ [0,1]` the burnt mass fraction —
/// **the doc's burn-progress `c`** (the code's `c`/`I_RC` slot is the
/// composition element `Z`, a pre-existing name). It is advected as a passive
/// scalar exactly like `Z` and is **inert (source-free) unless a combustion
/// occupant is scheduled**, so a plain shifting run never reads it and the
/// stations stay byte-for-byte (it defaults to 0 and 0 advects to 0). This is
/// a single compile-time-fixed component — NOT the S5b config-variable
/// `{ρX_k}` widening (SOLV-1 §3.4 change log).
pub const NCOMP: usize = 7;
pub const NPRIM: usize = 9;
pub const I_RHO: usize = 0;
pub const I_MR: usize = 1;
pub const I_MT: usize = 2;
pub const I_MZ: usize = 3;
pub const I_EN: usize = 4;
pub const I_RC: usize = 5;
/// Conserved burn-progress slot `ρb` (SOLV-4 §3.6); see the module const doc.
pub const I_RB: usize = 6;
/// Aux primitive slot: specific internal energy `e` (J/kg).
pub const I_EI: usize = 7;
/// Aux primitive slot: effective adiabatic exponent `Γ₁ = ρa²/p`.
pub const I_G1: usize = 8;

pub type Cons = [f64; NCOMP];
pub type Prim = [f64; NPRIM];

/// Build a primitive state from the physical slots, burn + aux slots zero
/// (the inert default — `b = 0`). Correct for `GammaLaw` (which never reads
/// burn or aux); a `TableEos`/blended primitive must come from its own
/// `prim_checked` (which fills them). `c` here is the composition `Z`.
pub const fn prim6(rho: f64, u_r: f64, u_t: f64, u_z: f64, p: f64, c: f64) -> Prim {
    [rho, u_r, u_t, u_z, p, c, 0.0, 0.0, 0.0]
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

/// Grid field names for `U`, in component order (FND-2 §3.4). `rho_b` is the
/// S6 burn-progress `ρb` (SOLV-4 §3.6; the doc's `c`).
pub const EULER_FIELDS: &[&str] = &[
    "rho",
    "mom_r",
    "mom_theta",
    "mom_z",
    "rho_e",
    "rho_c",
    "rho_b",
];

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
            rho * w[I_RB],
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
        // Burn progress advects as a passive scalar (inert under GammaLaw —
        // never a combustion occupant; carried so the state width is uniform).
        let b = u[I_RB] * inv;
        if !b.is_finite() {
            return Err("non-finite burn progress");
        }
        Ok([rho, ur, ut, uz, p, c, b, 0.0, 0.0])
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

#[derive(Debug, Clone, PartialEq)]
pub enum FlowError {
    /// A brick's N_θ is off the ladder (`4·2^k`, or 1 under the recorded
    /// axisymmetry assertion) — the ring-interface nesting premise fails.
    ThetaOffLadder { brick: usize, n_theta: u32 },
    /// Face-adjacent bricks differ by more than one ladder factor — the
    /// FND-2 §3.4 2:1 adjacency rule; re-tile the N_θ profile.
    ThetaAdjacency {
        a: (u32, u32),
        b: (u32, u32),
        nt_a: u32,
        nt_b: u32,
    },
    /// An N_θ jump sits within NGHOST cells of a run boundary or another
    /// jump along a pencil — the interface reconstruction has no uniform
    /// stencil (proper nesting, FND-2 §3.4); re-tile the N_θ profile.
    ThetaNesting {
        dir: &'static str,
        line: usize,
        at: usize,
    },
    /// A capability that requires uniform N_θ was scheduled on a mixed-N_θ
    /// world (combustion, cut geometry/SRD → plan S11).
    MixedThetaUnsupported { what: &'static str },
    /// A capability not yet carried by cut geometry at N_θ > 1 (S9 opened
    /// the Euler sweeps + SRD only; the rest rides plan S11).
    ThetaCutUnsupported { what: &'static str },
    /// SOLV-4 §3.6 (S8, 0.4.7): `S_T` crossed `S_T_MACH_LIMIT ×` the
    /// cell's own sound speed — the quasi-isobaric flamelet premises are
    /// broken (fast-deflagration/DDT class, out of the declared model
    /// form). A model-form limit on velocities, never on mesh rates.
    FrontCarrierScaleSeparation {
        i_r: usize,
        i_z: usize,
        i_theta: u32,
        s_t: f64,
        sound_speed: f64,
    },
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
    /// contour; never smear (META-1 P6). Sector-resolved at N_θ > 1 (S9;
    /// `i_theta` = 0 on N_θ = 1 worlds).
    CutCellUnmergeable {
        i_r: usize,
        i_z: usize,
        i_theta: u32,
    },
}

impl std::fmt::Display for FlowError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ThetaOffLadder { brick, n_theta } => write!(
                f,
                "brick {brick} carries N_θ = {n_theta}, off the ladder (4·2^k, or 1 under \
                 the recorded axisymmetry assertion) — ring-interface nesting undefined"
            ),
            Self::ThetaAdjacency { a, b, nt_a, nt_b } => write!(
                f,
                "face-adjacent bricks {a:?} (N_θ = {nt_a}) and {b:?} (N_θ = {nt_b}) differ \
                 by more than one ladder factor — FND-2 §3.4 2:1 adjacency; re-tile the \
                 N_θ profile"
            ),
            Self::ThetaNesting { dir, line, at } => write!(
                f,
                "N_θ jump at {dir}-index {at} on pencil line {line} sits within NGHOST \
                 cells of a run boundary or another jump — proper nesting (FND-2 §3.4); \
                 re-tile the N_θ profile"
            ),
            Self::MixedThetaUnsupported { what } => write!(
                f,
                "{what} requires uniform N_θ — refusing on this mixed-N_θ world rather \
                 than guessing"
            ),
            Self::ThetaCutUnsupported { what } => write!(
                f,
                "{what} is not yet carried by cut geometry at N_θ > 1 — refusing rather \
                 than guessing (S9 opened the Euler sweeps + SRD; the rest rides plan S11)"
            ),
            Self::FrontCarrierScaleSeparation {
                i_r,
                i_z,
                i_theta,
                s_t,
                sound_speed,
            } => write!(
                f,
                "front-carrier scale separation lost at (i_r={i_r}, i_z={i_z}, \
                 i_θ={i_theta}): S_T = {s_t:.3e} m/s exceeds S_T_MACH_LIMIT × the local \
                 sound speed {sound_speed:.3e} m/s — the quasi-isobaric flamelet model \
                 form does not cover this fast-deflagration/DDT-class state \
                 (SOLV-4 §3.6, S8)"
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
            Self::CutCellUnmergeable { i_r, i_z, i_theta } => write!(
                f,
                "cut cell (i_r={i_r}, i_z={i_z}, i_θ={i_theta}) below the SRD threshold \
                 has no flow-connected merge neighborhood — geometry under-resolved at \
                 this cell size (feature thinner than a cell); refine, never smear"
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
    /// The SOLV-4 §3.6 burn-progress source (S6). `Some` ⇒ `eval_rhs`
    /// accumulates the SOLV-4.4 rate law into the `ρb` slot at every node;
    /// `None` ⇒ `ρb` is an inert passive scalar (the shifting stations). The
    /// occupant carries the blended EOS + ignition closures (an occupant of
    /// this operator, selected as config data — Rule 13).
    pub combustion: Option<&'a Combustion<'a>>,
}

/// COUP-2 §3.1 per-evaluation ledger of one rhs evaluation, in conserved
/// units per second (κV-weighted): net and gross boundary-**port** fluxes
/// (run-boundary faces — domain BCs, wall faces incl. the declared stair
/// transpiration; interior faces telescope and are never ledgered) and
/// net/gross applied volumetric **sources** (geometric, wall-closure,
/// external intake). The SDC step composes these with its node weights;
/// gross magnitudes feed the §3.1.1 `S[q]` tolerance scale.
#[derive(Debug, Clone, Copy, Default)]
pub struct FlowLedger {
    pub port_net: [f64; NCOMP],
    pub port_abs: [f64; NCOMP],
    pub src_net: [f64; NCOMP],
    pub src_abs: [f64; NCOMP],
}

/// Per-brick scratch: primitives (AoS is fine for the CPU reference path;
/// the grid state itself stays SoA per FND-2 §3.9) and the rate-evaluation
/// data the SDC step composes.
pub(crate) struct Scratch {
    /// `NPRIM`-wide primitive (+aux) states per cell — the SDC step's gas
    /// class-D operand source after a `refresh_prims`/`eval_rhs` (S3).
    pub(crate) prim: Vec<Vec<Prim>>,
    pub(crate) rate: Vec<Vec<Cons>>,
    pub(crate) u0: Vec<Vec<Cons>>,
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
    /// COUP-2 ledger of the most recent [`Euler::eval_rhs`].
    pub(crate) ledger: FlowLedger,
}

/// Cells addressed as (brick index, θ-plane-major storage index
/// `i_theta·BRICK_CELLS + local`) — the field storage's native form (at
/// N_θ = 1 the index IS the local (r,z) index).
type BrickLocalCells = Vec<(usize, usize)>;

/// One sweep's (net, gross) port-flux partial (COUP-2 ledger).
type SweepPorts = ([f64; NCOMP], [f64; NCOMP]);

/// Opaque persistent stepping workspace (see [`Euler::workspace`]).
pub struct EulerWorkspace(pub(crate) Scratch);

impl EulerWorkspace {
    /// The COUP-2 ledger of the most recent [`Euler::eval_rhs`].
    pub fn ledger(&self) -> &FlowLedger {
        &self.0.ledger
    }

    /// The class-A rate buffer of the most recent [`Euler::eval_rhs`],
    /// per brick (cell index `idx = i_theta·BRICK_CELLS + local`). Exposed
    /// read-only for the S13 GPU residency cross-check (`crates/gpu`): the
    /// CPU stays the bit-exact reference the device kernels are scored
    /// against. No physics path reads it.
    pub fn rates(&self) -> &[Vec<Cons>] {
        &self.0.rate
    }
}

/// SRD bookkeeping (Berger & Giuliani 2020): the small-cell neighborhoods
/// in fixed owner order and the overlap counts n_j. Sector-resolved (S9):
/// each θ-plane's neighborhoods come from that sector's own κ and
/// apertures (`neighborhood_cells`), and members address cells by their
/// θ-plane-major storage index — at N_θ = 1 `idx = local`, so the pre-S9
/// bookkeeping is reproduced bit-identically.
struct CutScratch {
    /// Per small sector-cell, in the fixed build order (θ-plane-major:
    /// sector 0's (r,z) sweep first, then sector 1, …): members as
    /// (bi, idx) with members[0] = owner, and the members' gas volumes
    /// κ_j·V — the merge weights.
    small: Vec<(BrickLocalCells, Vec<f64>)>,
    /// Overlap count n_j per (bi·N_θ·BRICK_CELLS + idx): 1 + the number of
    /// OTHER cells' neighborhoods containing j (every cell's own
    /// neighborhood — trivial for regular cells — is the 1).
    counts: Vec<u32>,
}

impl<E: EosLaw + Sync> Euler<'_, E> {
    /// S8/S9 N_θ legality (module doc): per-brick ladder membership, 2:1
    /// face-adjacency, uniform-N_θ-only capabilities (cut geometry/SRD and
    /// combustion at mixed N_θ → S11; combustion on cut worlds at
    /// N_θ > 1 → S11). Pure structure checks — a pure function of the
    /// grid, run before every evaluation (cheap: O(bricks)).
    fn validate(&self, g: &Grid) -> Result<(), FlowError> {
        let nt0 = g.brick(0).n_theta();
        let mixed = g.bricks().iter().any(|b| b.n_theta() != nt0);
        for (bi, b) in g.bricks().iter().enumerate() {
            let nt = b.n_theta();
            if nt != 1 && !crucible_grid::theta_ladder_aligned(nt) {
                return Err(FlowError::ThetaOffLadder {
                    brick: bi,
                    n_theta: nt,
                });
            }
            // 2:1 adjacency toward the high-side neighbors (each pair
            // checked once).
            for (dbr, dbz) in [(1u32, 0u32), (0, 1)] {
                if let Some(ni) = g.brick_index_by_coords(b.br() + dbr, b.bz() + dbz) {
                    let ntn = g.brick(ni).n_theta();
                    let (lo, hi) = (nt.min(ntn), nt.max(ntn));
                    if hi > 2 * lo {
                        return Err(FlowError::ThetaAdjacency {
                            a: (b.br(), b.bz()),
                            b: (b.br() + dbr, b.bz() + dbz),
                            nt_a: nt,
                            nt_b: ntn,
                        });
                    }
                }
            }
        }
        if mixed {
            if g.has_cut_geometry() {
                return Err(FlowError::MixedThetaUnsupported {
                    what: "cut geometry / State Redistribution (mixed-N_θ cut worlds \
                           ride plan S11; uniform N_θ is legal — S9)",
                });
            }
            if self.combustion.is_some() {
                return Err(FlowError::MixedThetaUnsupported {
                    what: "the SOLV-4 §3.6 combustion operator (ring-interface \
                           c-diffusion rides plan S11)",
                });
            }
        } else if nt0 > 1
            && g.has_cut_geometry()
            && self.combustion.is_some()
            && !g.geometry_is_theta_uniform()
        {
            // S9 opened the sweeps + SRD to per-sector cut geometry; S11
            // extended the combustion operator's D_c face stencil to read
            // per-sector apertures too (kv, meridional faces, and the
            // θ-face aperture), so a REVOLVED (θ-uniform) cut wall — every
            // ◆C3 engine world — is legal. A genuinely θ-varying wall
            // (CSG/STL) still rides the θ-varying-geometry wave.
            return Err(FlowError::ThetaCutUnsupported {
                what: "the SOLV-4 §3.6 combustion operator on θ-VARYING cut geometry \
                       (revolved cut walls are legal since S11; CSG/STL θ-variation \
                       rides the θ-varying-geometry wave)",
            });
        }
        Ok(())
    }

    fn scratch(&self, g: &Grid) -> Result<Scratch, FlowError> {
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
            // Uniform N_θ guaranteed by `validate` (mixed-N_θ cut worlds
            // refuse → plan S11); the SRD weights key on it. Sector-resolved
            // (S9): each θ-plane's small cells and neighborhoods come from
            // that sector's own κ/apertures, in the fixed θ-plane-major
            // build order — at N_θ = 1 this is exactly the pre-S9 (r,z)
            // sweep, bit-identical.
            let nt = g.brick(0).n_theta();
            let plane = nt as usize * BRICK_CELLS;
            let mut small = Vec::new();
            let mut counts = vec![1u32; nb * plane];
            let to_bl = |(r, z, j): (usize, usize, u32)| {
                let bi = bmap[(r / BRICK) * nbz + z / BRICK].expect("active cell's brick");
                (
                    bi,
                    j as usize * BRICK_CELLS + (r % BRICK) * BRICK + z % BRICK,
                )
            };
            for j in 0..nt {
                for i_r in 0..n_r {
                    for i_z in 0..n_z {
                        if !act[i_r * n_z + i_z] || g.kappa_at(i_r, j, i_z) >= KAPPA_SRD {
                            continue;
                        }
                        let cells = neighborhood_cells(g, i_r, i_z, j)?;
                        let kv: Vec<f64> = cells
                            .iter()
                            .map(|&(r, z, jj)| g.kappa_at(r, jj, z) * g.cell_volume(r, nt))
                            .collect();
                        let members: Vec<(usize, usize)> = cells.into_iter().map(to_bl).collect();
                        for &(bi, idx) in &members[1..] {
                            counts[bi * plane + idx] += 1;
                        }
                        small.push((members, kv));
                    }
                }
            }
            Some(CutScratch { small, counts })
        } else {
            None
        };
        // Per-brick plane sizes (S8): each brick's buffers span its own
        // N_θ × 64 cells.
        let plane = |bi: usize| g.brick(bi).n_theta() as usize * BRICK_CELLS;
        Ok(Scratch {
            prim: (0..nb).map(|bi| vec![[0.0; NPRIM]; plane(bi)]).collect(),
            rate: (0..nb).map(|bi| vec![[0.0; NCOMP]; plane(bi)]).collect(),
            u0: (0..nb).map(|bi| vec![[0.0; NCOMP]; plane(bi)]).collect(),
            bmap,
            nbz,
            act,
            cut,
            primed: false,
            ledger: FlowLedger::default(),
        })
    }

    /// Build a persistent workspace for repeated stepping: reuse across
    /// steps keeps the previous stage's projected pressures as warm-start
    /// hints for the equilibrium projection (`EosLaw::prim_checked_hinted`)
    /// and retires the per-step scratch allocation churn. One-shot
    /// [`Euler::step`] builds a fresh (cold) one each call.
    pub fn workspace(&self, g: &Grid) -> Result<EulerWorkspace, FlowError> {
        self.validate(g)?;
        Ok(EulerWorkspace(self.scratch(g)?))
    }

    /// Refresh the workspace's primitive cache from the CURRENT grid
    /// state — no rate evaluation, no ledger touch. The SDC step's gas
    /// class-D solve derives its operands from this (S3); warm-start
    /// hints ride along exactly as in `eval_rhs`.
    pub fn refresh_prims(
        &self,
        g: &Grid,
        f: &EulerFields,
        ws: &mut EulerWorkspace,
    ) -> Result<(), FlowError> {
        self.validate(g)?;
        self.fill_prims(g, f, &mut ws.0)
    }

    /// Evaluate `L(U)` — the flux-divergence + source contribution (SOLV-1
    /// §2 contract) — into the workspace's rate buffer at time `t`, with
    /// the COUP-2 per-evaluation port/source ledger. The caller (COUP-3's
    /// SDC step, `crate::sdc`) owns node weights, state composition, the
    /// stagewise SRD passes, and the Δt rule.
    pub fn eval_rhs(
        &self,
        g: &Grid,
        f: &EulerFields,
        ws: &mut EulerWorkspace,
        t: f64,
    ) -> Result<(), FlowError> {
        self.validate(g)?;
        self.rhs(g, f, &mut ws.0, t)?;
        // The burn-progress source (S6) rides the class-A rate: after `rhs`
        // has filled the primitive cache, rate, and ledger, SOLV-4.4 adds its
        // ρb source + source-ledger, so the SDC step composes and audits it
        // exactly like the flow's own geometric source (no SDC-step change).
        if let Some(comb) = self.combustion {
            comb.accumulate(g, &f.ids(), &mut ws.0)?;
        }
        Ok(())
    }

    /// Apply the State-Redistribution pass to the CURRENT field state —
    /// the SDC step calls this after each node-state composition (Berger &
    /// Giuliani apply SRD stagewise): the provisional divide-by-κV update
    /// on a small cell is merged into its neighborhood — conservative by
    /// construction, and the reason Δt keeps the UNCUT CFL. A no-op on
    /// worlds without cut geometry (bit-identity preserved).
    pub fn apply_srd(&self, g: &mut Grid, f: &EulerFields, ws: &EulerWorkspace) {
        Self::srd(g, f, &ws.0);
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
    /// member order ⇒ deterministic. Sector-resolved (S9): the
    /// neighborhoods carry θ-plane-major member indices and per-sector
    /// κ_j·V weights, so one pass covers every θ-plane; at N_θ = 1 the
    /// index IS the (r,z) local index and the pass is the pre-S9
    /// arithmetic, bit-identical.
    fn srd(g: &mut Grid, f: &EulerFields, s: &Scratch) {
        let Some(cut) = &s.cut else { return };
        if cut.small.is_empty() {
            return;
        }
        // Uniform N_θ on cut worlds (`validate`); the counts stride.
        let plane = g.brick(0).n_theta() as usize * BRICK_CELLS;
        let ids = f.ids();
        // (bi, idx) → (Σ Q contributions, is-small-owner).
        let mut acc: std::collections::BTreeMap<(usize, usize), ([f64; NCOMP], bool)> =
            std::collections::BTreeMap::new();
        for (members, _) in &cut.small {
            acc.entry(members[0]).or_insert(([0.0; NCOMP], false)).1 = true;
        }
        let read = |g: &Grid, bi: usize, idx: usize| -> Cons {
            let b = g.brick(bi);
            std::array::from_fn(|k| b.field(ids[k])[idx])
        };
        // Neighborhood averages against the CURRENT (provisional) state.
        let mut qs: Vec<[f64; NCOMP]> = Vec::with_capacity(cut.small.len());
        for (members, kv) in &cut.small {
            let mut num = [0.0f64; NCOMP];
            let mut den = 0.0f64;
            for (&(bi, idx), &kvj) in members.iter().zip(kv) {
                let w = kvj / f64::from(cut.counts[bi * plane + idx]);
                let u = read(g, bi, idx);
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
            for &(bi, idx) in members {
                acc.entry((bi, idx)).or_insert(([0.0; NCOMP], false));
            }
        }
        for ((bi, idx), (a, owner)) in &mut acc {
            if !*owner {
                let u = read(g, *bi, *idx);
                for k in 0..NCOMP {
                    a[k] += u[k];
                }
            }
        }
        for ((members, _), q) in cut.small.iter().zip(&qs) {
            for &(bi, idx) in members {
                let (a, _) = acc.get_mut(&(bi, idx)).expect("member entry");
                for k in 0..NCOMP {
                    a[k] += q[k];
                }
            }
        }
        // Write back U_j = acc_j / n_j in fixed (bi, idx) order.
        for ((bi, idx), (a, _)) in &acc {
            let n = f64::from(cut.counts[*bi * plane + *idx]);
            for (k, &id) in ids.iter().enumerate() {
                g.brick_field_mut(*bi, id)[*idx] = a[k] / n;
            }
        }
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
        self.validate(g)?;
        let ids = f.ids();
        let (dr, dz) = (g.spec().dr, g.spec().dz);
        // The front-carrier member's operands (SOLV-4 §3.6, S8): the same Δ
        // the rate law's matched coefficients use, and the meridional
        // 1/Δ² sum of its explicit diffusion (the θ term is per cell).
        let delta = (dr * dz).sqrt();
        let inv_sq_rz = 1.0 / (dr * dr) + 1.0 / (dz * dz);
        let partials: Vec<Result<f64, FlowError>> = (0..g.n_bricks())
            .into_par_iter()
            .map(|bi| {
                let b = g.brick(bi);
                let nt = b.n_theta();
                let dtheta = std::f64::consts::TAU / f64::from(nt);
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
                                // COUP-3 §3.4 (S8): the front-carrier signal
                                // joins the reduction; SOLV-4's separation
                                // guard refuses at the sonic end.
                                if let Some(comb) = self.combustion {
                                    let inv_sq = if nt > 1 {
                                        let arc = g.r_center(i_r) * dtheta;
                                        inv_sq_rz + 1.0 / (arc * arc)
                                    } else {
                                        inv_sq_rz
                                    };
                                    let (sf, s_t) = comb
                                        .front_carrier_signal(&w, delta, inv_sq)
                                        .map_err(|what| FlowError::NonPhysicalState {
                                            i_r,
                                            i_z,
                                            i_theta: j,
                                            what,
                                        })?;
                                    // The model-form guard compares
                                    // VELOCITIES (S_T vs the local c) —
                                    // mesh rates belong only to the Δt
                                    // member (SOLV-4 0.4.7).
                                    if s_t > S_T_MACH_LIMIT * c {
                                        return Err(FlowError::FrontCarrierScaleSeparation {
                                            i_r,
                                            i_z,
                                            i_theta: j,
                                            s_t,
                                            sound_speed: c,
                                        });
                                    }
                                    sig += sf;
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
    /// into `s.rate`, with the per-evaluation COUP-2 ledger (ports from the
    /// r/z sweeps' run boundaries; the θ sweep is periodic — its ring
    /// fluxes telescope exactly; sources from the geometric/closure/
    /// external pass).
    fn rhs(&self, g: &Grid, f: &EulerFields, s: &mut Scratch, t: f64) -> Result<(), FlowError> {
        self.fill_prims(g, f, s)?;
        for rate in &mut s.rate {
            for cell in rate.iter_mut() {
                *cell = [0.0; NCOMP];
            }
        }
        s.ledger = FlowLedger::default();
        self.sweep_r(g, s, t)?;
        self.sweep_theta(g, s);
        self.sweep_z(g, s, t)?;
        self.add_sources(g, s, t);
        Ok(())
    }

    /// Parallel over bricks (each brick's prim block is its own task —
    /// FND-2 §3.7: ownership partition, so results are thread-count-
    /// independent). Errors are gathered per brick and the FIRST in brick
    /// order is raised — the same diagnosis the serial sweep chose.
    fn fill_prims(&self, g: &Grid, f: &EulerFields, s: &mut Scratch) -> Result<(), FlowError> {
        let ids = f.ids();
        let primed = s.primed;
        let results: Vec<Result<(), FlowError>> = s
            .prim
            .par_iter_mut()
            .enumerate()
            .map(|(bi, prim)| {
                let b = g.brick(bi);
                let nt = b.n_theta();
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
    /// takes the configured BC. Each run further decomposes into maximal
    /// **uniform-N_θ segments** (S8, module doc), processed in TWO passes:
    /// pass 1 reconstructs and stores every segment's fluxes; the fix-up
    /// replaces each coarse side's N_θ-jump face entry with the exact
    /// aggregate of the fine side's children (the fine side owns the flux);
    /// pass 2 accumulates rates as the SINGLE difference
    /// `(af[q] − af[q+1])/(κV)` — the same one-rounding well-balanced form
    /// as a uniform run, which is what keeps a uniform state a bitwise
    /// fixed point across the interface (a split accumulation breaks the
    /// exact cancellation against the geometric pressure source — measured,
    /// the first cut of this sweep did exactly that).
    fn sweep_r(&self, g: &Grid, s: &mut Scratch, t: f64) -> Result<(), FlowError> {
        let spec = g.spec();
        let (n, n_z) = (spec.n_r, spec.n_z);
        let (r0, dr) = (spec.r_min, spec.dr);
        let on_axis = r0 == 0.0;

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
        let results: Vec<Result<SweepPorts, FlowError>> = rows
            .into_par_iter()
            .enumerate()
            .map(|(bz, mut row)| {
                let mut slot: Vec<Option<usize>> = vec![None; nbr];
                for (k, (bi, _)) in row.iter().enumerate() {
                    slot[g.brick(*bi).br() as usize] = Some(k);
                }
                // Metric caches per distinct N_θ in this row — the same
                // function calls (same bits) as the former per-call vectors.
                let mut metrics: std::collections::BTreeMap<u32, (Vec<f64>, Vec<f64>)> =
                    std::collections::BTreeMap::new();
                for (bi, _) in &row {
                    let nt = g.brick(*bi).n_theta();
                    metrics.entry(nt).or_insert_with(|| {
                        let area = (0..=n)
                            .map(|fi| {
                                if fi < n {
                                    g.face_area_r(fi, false, nt)
                                } else {
                                    g.face_area_r(n - 1, true, nt)
                                }
                            })
                            .collect();
                        let vol = (0..n).map(|i| g.cell_volume(i, nt)).collect();
                        (area, vol)
                    });
                }
                // Brick N_θ of an active cell in this row.
                let nt_of = |ii: usize| -> u32 {
                    let bi = bmap[(ii / BRICK) * nbz + bz].expect("active cell's brick");
                    g.brick(bi).n_theta()
                };
                let mut w = vec![[0.0f64; NPRIM]; n + 2 * NGHOST];
                let mut fl = vec![[0.0f64; NPRIM]; n + 1];
                let mut fr = vec![[0.0f64; NPRIM]; n + 1];
                let mut ap = vec![1.0f64; n + 1];
                let mut port_net = [0.0f64; NCOMP];
                let mut port_abs = [0.0f64; NCOMP];
                for i_z in bz * BRICK..((bz + 1) * BRICK).min(n_z) {
                    let lz = i_z % BRICK;
                    let z = g.z_center(i_z);
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
                        // Maximal uniform-N_θ segments of this run; proper
                        // nesting: a segment owning an interface spans at
                        // least NGHOST uniform cells (module doc).
                        let mut segs: Vec<SweepSeg> = Vec::new();
                        let mut seg_start = start;
                        while seg_start < start + len {
                            let nts = nt_of(seg_start);
                            let mut seg_end = seg_start + 1;
                            while seg_end < start + len && nt_of(seg_end) == nts {
                                seg_end += 1;
                            }
                            let seg_len = seg_end - seg_start;
                            if (seg_start > start || seg_end < start + len) && seg_len < NGHOST {
                                return Err(FlowError::ThetaNesting {
                                    dir: "r",
                                    line: i_z,
                                    at: seg_start,
                                });
                            }
                            segs.push(SweepSeg::new(seg_start, seg_len, nts));
                            seg_start = seg_end;
                        }
                        // Pass 1: reconstruct + flux every segment, stored.
                        for si in 0..segs.len() {
                            let (seg_start, seg_len, nts) =
                                (segs[si].start, segs[si].len, segs[si].nts);
                            let seg_end = seg_start + seg_len;
                            let left_nt = (si > 0).then(|| segs[si - 1].nts);
                            let right_nt = (si + 1 < segs.len()).then(|| segs[si + 1].nts);
                            let area = &metrics[&nts].0;
                            for j in 0..nts {
                                let theta = Grid::theta_center(j, nts);
                                // Per-sector cut geometry (FND-3 §3.3, S9):
                                // κ and the face apertures of THIS θ-plane.
                                // At N_θ = 1 plane 0 carries the identical
                                // bits the pre-S9 (r,z) reads produced.
                                for (q, ii) in (seg_start..seg_end).enumerate() {
                                    let bi =
                                        bmap[(ii / BRICK) * nbz + bz].expect("active cell's brick");
                                    let b = g.brick(bi);
                                    let local = (ii % BRICK) * BRICK + lz;
                                    segs[si].kap[j as usize * seg_len + q] = b.kappa_cell(j, local);
                                    ap[q] = b.aperture_cell(FaceDir::RMinus, j, local);
                                    if q + 1 == seg_len {
                                        ap[seg_len] = b.aperture_cell(FaceDir::RPlus, j, local);
                                    }
                                    w[NGHOST + q] = prim[bi][j as usize * BRICK_CELLS + local];
                                }
                                // Low ghosts.
                                if let Some(lnt) = left_nt {
                                    // N_θ-interface ghosts from the left
                                    // segment: prolong (coarser) / restrict
                                    // (finer) — FND-2 §3.4.
                                    for k in 1..=NGHOST {
                                        let ii = seg_start - k;
                                        let bi = bmap[(ii / BRICK) * nbz + bz]
                                            .expect("active cell's brick");
                                        let local = (ii % BRICK) * BRICK + lz;
                                        w[NGHOST - k] = theta_mapped(&prim[bi], local, j, nts, lnt);
                                    }
                                } else if start == 0 && on_axis {
                                    // FND-2 §3.2 θ↔θ+π parity-pair gather:
                                    // ê_r and ê_θ both flip. At N_θ = 1 the
                                    // partner is this cell — the identical
                                    // mirror arithmetic.
                                    let jp = Grid::axis_pair(j, nts) as usize;
                                    for k in 1..=NGHOST {
                                        let ii = start + (k - 1).min(seg_len - 1);
                                        let bi = bmap[(ii / BRICK) * nbz + bz]
                                            .expect("active cell's brick");
                                        let local = (ii % BRICK) * BRICK + lz;
                                        let mut m = prim[bi][jp * BRICK_CELLS + local];
                                        m[I_MR] = -m[I_MR];
                                        m[I_MT] = -m[I_MT];
                                        w[NGHOST - k] = m;
                                    }
                                } else if start == 0 {
                                    self.fill_ghosts_low(
                                        &mut w,
                                        seg_len,
                                        &self.bcs.r_inner,
                                        I_MR,
                                        |k| (r0 - (k as f64 - 0.5) * dr, theta, z),
                                        t,
                                    )?;
                                } else {
                                    self.wall_ghosts_low(
                                        &mut w,
                                        seg_len,
                                        I_MR,
                                        r0 + start as f64 * dr,
                                        z,
                                    );
                                }
                                // High ghosts.
                                if let Some(rnt) = right_nt {
                                    for k in 1..=NGHOST {
                                        let ii = seg_start + seg_len + k - 1;
                                        let bi = bmap[(ii / BRICK) * nbz + bz]
                                            .expect("active cell's brick");
                                        let local = (ii % BRICK) * BRICK + lz;
                                        w[NGHOST + seg_len - 1 + k] =
                                            theta_mapped(&prim[bi], local, j, nts, rnt);
                                    }
                                } else if start + len == n {
                                    self.fill_ghosts_high(
                                        &mut w,
                                        seg_len,
                                        &self.bcs.r_outer,
                                        I_MR,
                                        |k| (r0 + (n as f64 + k as f64 - 0.5) * dr, theta, z),
                                        t,
                                    )?;
                                } else {
                                    self.wall_ghosts_high(
                                        &mut w,
                                        seg_len,
                                        I_MR,
                                        r0 + (start + len) as f64 * dr,
                                        z,
                                    );
                                }
                                ppm_faces(&w[..seg_len + 2 * NGHOST], seg_len, &mut fl, &mut fr);
                                for fi in 0..=seg_len {
                                    let flux = hllc_flux(&fl[fi], &fr[fi], I_MR, &self.eos);
                                    // Aperture-weighted open area (FND-3
                                    // §3.3); ap = 1.0 exactly on full-box
                                    // worlds (`(A·1.0)·F ≡ A·F` bitwise).
                                    let aa = area[seg_start + fi] * ap[fi];
                                    let dst = &mut segs[si].af[j as usize * (seg_len + 1) + fi];
                                    for k in 0..NCOMP {
                                        dst[k] = aa * flux[k];
                                    }
                                }
                            }
                        }
                        // Fix-up (FND-2 §3.4): the FINE side owns each
                        // N_θ-jump face's flux; the coarse side's boundary
                        // entry becomes the aggregate of the fine children
                        // (each fine af already carries its own area, so
                        // the aggregate is the plain sum — exact where the
                        // children agree, by the power-of-two metric split).
                        reflux_fixup(&mut segs, |a, b| std::array::from_fn(|k| a[k] + b[k]));
                        // Pass 2: accumulate rates as the SINGLE difference
                        // (af[q] − af[q+1])/(κV) — the one-rounding
                        // well-balanced form — plus the COUP-2 ledger
                        // (run-boundary faces only; interfaces telescope).
                        for (si, seg) in segs.iter().enumerate() {
                            let vol = &metrics[&seg.nts].1;
                            let first = si == 0;
                            let last = si + 1 == segs.len();
                            let stride = seg.len + 1;
                            for j in 0..seg.nts as usize {
                                let af = &seg.af[j * stride..(j + 1) * stride];
                                for q in 0..seg.len {
                                    let ii = seg.start + q;
                                    let k_slot = slot[ii / BRICK].expect("active brick in row");
                                    let idx = j * BRICK_CELLS + (ii % BRICK) * BRICK + lz;
                                    let rate = &mut row[k_slot].1[idx];
                                    let kap = seg.kap[j * seg.len + q];
                                    for k in 0..NCOMP {
                                        rate[k] += (af[q][k] - af[q + 1][k]) / (kap * vol[ii]);
                                    }
                                }
                                if first && last {
                                    for k in 0..NCOMP {
                                        port_net[k] += af[0][k] - af[seg.len][k];
                                        port_abs[k] += af[0][k].abs() + af[seg.len][k].abs();
                                    }
                                } else if first {
                                    for k in 0..NCOMP {
                                        port_net[k] += af[0][k];
                                        port_abs[k] += af[0][k].abs();
                                    }
                                } else if last {
                                    for k in 0..NCOMP {
                                        port_net[k] -= af[seg.len][k];
                                        port_abs[k] += af[seg.len][k].abs();
                                    }
                                }
                            }
                        }
                    }
                }
                Ok((port_net, port_abs))
            })
            .collect();
        for r in results {
            let (net, abs) = r?;
            for k in 0..NCOMP {
                s.ledger.port_net[k] += net[k];
                s.ledger.port_abs[k] += abs[k];
            }
        }
        Ok(())
    }

    /// Azimuthal sweep: periodic ring pencils, each inside one brick at its
    /// own N_θ (S8). At N_θ = 1 both faces carry the same computed flux and
    /// cancel exactly — geometry, not a code branch. Cut geometry (S9,
    /// FND-3 §3.3): each θ-face flux is aperture-weighted and the rate
    /// divide uses the sector's own κ. **One canonical aperture side per
    /// face**: face `fi` (between ring cells fi−1 and fi, mod N_θ) reads
    /// the θ+ aperture of the cell BEFORE it, so both accumulation
    /// directions use the identical bits and the periodic face 0 ≡ face
    /// N_θ (their fluxes already agree by the periodic gather; the shared
    /// bits keep the ring telescoping exact). Full-box worlds: ap = 1.0
    /// and κ = 1.0 exactly, so `F·1.0 ≡ F` and `x/1.0 ≡ x` — the S12
    /// arithmetic-identity idiom, every pre-S9 world bit-identical.
    fn sweep_theta(&self, g: &Grid, s: &mut Scratch) {
        let a_th = g.face_area_theta();
        let n_max = g.bricks().iter().map(|b| b.n_theta()).max().unwrap_or(1) as usize;
        let mut w = vec![[0.0f64; NPRIM]; n_max + 2 * NGHOST];
        let mut fl = vec![[0.0f64; NPRIM]; n_max + 1];
        let mut fr = vec![[0.0f64; NPRIM]; n_max + 1];
        let mut af = vec![[0.0f64; NCOMP]; n_max + 1];
        let mut ap = vec![1.0f64; n_max + 1];

        for bi in 0..g.n_bricks() {
            let b = g.brick(bi);
            let mask = b.mask();
            let n = b.n_theta() as usize;
            for local in 0..BRICK_CELLS {
                if mask & (1u64 << local) == 0 {
                    continue;
                }
                let (i_r, _) = b.global_rz(local);
                // One per-ring metric ratio `A_θ/V` (= 1/(r̄·Δθ) exactly),
                // applied to the flux difference — the same conditioning
                // rule as the z sweep.
                let inv = a_th / g.cell_volume(i_r, b.n_theta());
                // Periodic gather with wrapped ghosts: pencil index pi maps
                // to ring cell (pi − NGHOST) mod n, offset by NGHOST·n so
                // the subtraction cannot underflow at any n ≥ 1.
                for (pi, cell) in w[..n + 2 * NGHOST].iter_mut().enumerate() {
                    let jj = (pi + NGHOST * n - NGHOST) % n;
                    *cell = s.prim[bi][jj * BRICK_CELLS + local];
                }
                // Canonical θ-face apertures (module note above): face fi
                // reads the θ+ aperture of ring cell (fi − 1) mod n.
                for (fi, a) in ap[..=n].iter_mut().enumerate() {
                    let jp = ((fi + n - 1) % n) as u32;
                    *a = b.aperture_cell(FaceDir::ThetaPlus, jp, local);
                }
                ppm_faces(&w[..n + 2 * NGHOST], n, &mut fl, &mut fr);
                for fi in 0..=n {
                    let flux = hllc_flux(&fl[fi], &fr[fi], I_MT, &self.eos);
                    for k in 0..NCOMP {
                        af[fi][k] = flux[k] * ap[fi];
                    }
                }
                for j in 0..n {
                    let kap = b.kappa_cell(j as u32, local);
                    let rate = &mut s.rate[bi][j * BRICK_CELLS + local];
                    for k in 0..NCOMP {
                        rate[k] += (af[j][k] - af[j + 1][k]) * inv / kap;
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
    /// Mixed N_θ (S8): pencils decompose into uniform-N_θ segments exactly
    /// as in `sweep_r`; in this sweep's metric-ratio form the fine→coarse
    /// aggregate carries the exact area ratio `A_zf/A_zc = 1/2`.
    fn sweep_z(&self, g: &Grid, s: &mut Scratch, t: f64) -> Result<(), FlowError> {
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
        let results: Vec<Result<SweepPorts, FlowError>> = cols
            .into_par_iter()
            .enumerate()
            .map(|(br, mut col)| {
                let mut slot: Vec<Option<usize>> = vec![None; nbz];
                for (k, (bi, _)) in col.iter().enumerate() {
                    slot[g.brick(*bi).bz() as usize] = Some(k);
                }
                let nt_of = |ii: usize| -> u32 {
                    let bi = bmap[br * nbz + ii / BRICK].expect("active cell's brick");
                    g.brick(bi).n_theta()
                };
                let mut w = vec![[0.0f64; NPRIM]; n + 2 * NGHOST];
                let mut fl = vec![[0.0f64; NPRIM]; n + 1];
                let mut fr = vec![[0.0f64; NPRIM]; n + 1];
                let mut ap = vec![1.0f64; n + 1];
                let mut port_net = [0.0f64; NCOMP];
                let mut port_abs = [0.0f64; NCOMP];
                for i_r in br * BRICK..((br + 1) * BRICK).min(n_r) {
                    let lr = i_r % BRICK;
                    let r = g.r_center(i_r);
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
                        // Maximal uniform-N_θ segments (see sweep_r).
                        let mut segs: Vec<SweepSeg> = Vec::new();
                        let mut seg_start = start;
                        while seg_start < start + len {
                            let nts = nt_of(seg_start);
                            let mut seg_end = seg_start + 1;
                            while seg_end < start + len && nt_of(seg_end) == nts {
                                seg_end += 1;
                            }
                            let seg_len = seg_end - seg_start;
                            if (seg_start > start || seg_end < start + len) && seg_len < NGHOST {
                                return Err(FlowError::ThetaNesting {
                                    dir: "z",
                                    line: i_r,
                                    at: seg_start,
                                });
                            }
                            segs.push(SweepSeg::new(seg_start, seg_len, nts));
                            seg_start = seg_end;
                        }
                        // Pass 1: reconstruct + flux every segment, stored.
                        for si in 0..segs.len() {
                            let (seg_start, seg_len, nts) =
                                (segs[si].start, segs[si].len, segs[si].nts);
                            let seg_end = seg_start + seg_len;
                            let left_nt = (si > 0).then(|| segs[si - 1].nts);
                            let right_nt = (si + 1 < segs.len()).then(|| segs[si + 1].nts);
                            for j in 0..nts {
                                let theta = Grid::theta_center(j, nts);
                                // Per-sector cut geometry (FND-3 §3.3, S9)
                                // — see sweep_r.
                                for (q, ii) in (seg_start..seg_end).enumerate() {
                                    let bi =
                                        bmap[br * nbz + ii / BRICK].expect("active cell's brick");
                                    let b = g.brick(bi);
                                    let local = lr * BRICK + (ii % BRICK);
                                    segs[si].kap[j as usize * seg_len + q] = b.kappa_cell(j, local);
                                    ap[q] = b.aperture_cell(FaceDir::ZMinus, j, local);
                                    if q + 1 == seg_len {
                                        ap[seg_len] = b.aperture_cell(FaceDir::ZPlus, j, local);
                                    }
                                    w[NGHOST + q] = prim[bi][j as usize * BRICK_CELLS + local];
                                }
                                if let Some(lnt) = left_nt {
                                    for k in 1..=NGHOST {
                                        let ii = seg_start - k;
                                        let bi = bmap[br * nbz + ii / BRICK]
                                            .expect("active cell's brick");
                                        let local = lr * BRICK + (ii % BRICK);
                                        w[NGHOST - k] = theta_mapped(&prim[bi], local, j, nts, lnt);
                                    }
                                } else if start == 0 {
                                    self.fill_ghosts_low(
                                        &mut w,
                                        seg_len,
                                        &self.bcs.z_lo,
                                        I_MZ,
                                        |k| (r, theta, z0 - (k as f64 - 0.5) * dz),
                                        t,
                                    )?;
                                } else {
                                    self.wall_ghosts_low(
                                        &mut w,
                                        seg_len,
                                        I_MZ,
                                        r,
                                        z0 + start as f64 * dz,
                                    );
                                }
                                if let Some(rnt) = right_nt {
                                    for k in 1..=NGHOST {
                                        let ii = seg_start + seg_len + k - 1;
                                        let bi = bmap[br * nbz + ii / BRICK]
                                            .expect("active cell's brick");
                                        let local = lr * BRICK + (ii % BRICK);
                                        w[NGHOST + seg_len - 1 + k] =
                                            theta_mapped(&prim[bi], local, j, nts, rnt);
                                    }
                                } else if start + len == n {
                                    self.fill_ghosts_high(
                                        &mut w,
                                        seg_len,
                                        &self.bcs.z_hi,
                                        I_MZ,
                                        |k| (r, theta, z0 + (n as f64 + k as f64 - 0.5) * dz),
                                        t,
                                    )?;
                                } else {
                                    self.wall_ghosts_high(
                                        &mut w,
                                        seg_len,
                                        I_MZ,
                                        r,
                                        z0 + (start + len) as f64 * dz,
                                    );
                                }
                                ppm_faces(&w[..seg_len + 2 * NGHOST], seg_len, &mut fl, &mut fr);
                                for fi in 0..=seg_len {
                                    let flux = hllc_flux(&fl[fi], &fr[fi], I_MZ, &self.eos);
                                    // ap = 1.0 exactly on full-box worlds
                                    // (`F·1.0 ≡ F`).
                                    let dst = &mut segs[si].af[j as usize * (seg_len + 1) + fi];
                                    for k in 0..NCOMP {
                                        dst[k] = flux[k] * ap[fi];
                                    }
                                }
                            }
                        }
                        // Fix-up (FND-2 §3.4): in this sweep's metric-ratio
                        // form the fine→coarse aggregate carries the exact
                        // area ratio A_zf/A_zc = 1/2.
                        reflux_fixup(&mut segs, |a, b| {
                            std::array::from_fn(|k| 0.5 * (a[k] + b[k]))
                        });
                        // Pass 2: single-difference accumulation + ledger
                        // (see sweep_r); the segment's ring z-face area
                        // restores conserved units on the ledger side.
                        for (si, seg) in segs.iter().enumerate() {
                            let a_z = g.face_area_z(i_r, seg.nts);
                            let first = si == 0;
                            let last = si + 1 == segs.len();
                            let stride = seg.len + 1;
                            for j in 0..seg.nts as usize {
                                let af = &seg.af[j * stride..(j + 1) * stride];
                                for q in 0..seg.len {
                                    let ii = seg.start + q;
                                    let k_slot = slot[ii / BRICK].expect("active brick in column");
                                    let idx = j * BRICK_CELLS + lr * BRICK + (ii % BRICK);
                                    let rate = &mut col[k_slot].1[idx];
                                    let kap = seg.kap[j * seg.len + q];
                                    for k in 0..NCOMP {
                                        rate[k] += (af[q][k] - af[q + 1][k]) * inv_dz / kap;
                                    }
                                }
                                if first && last {
                                    for k in 0..NCOMP {
                                        port_net[k] += a_z * (af[0][k] - af[seg.len][k]);
                                        port_abs[k] +=
                                            a_z * (af[0][k].abs() + af[seg.len][k].abs());
                                    }
                                } else if first {
                                    for k in 0..NCOMP {
                                        port_net[k] += a_z * af[0][k];
                                        port_abs[k] += a_z * af[0][k].abs();
                                    }
                                } else if last {
                                    for k in 0..NCOMP {
                                        port_net[k] -= a_z * af[seg.len][k];
                                        port_abs[k] += a_z * af[seg.len][k].abs();
                                    }
                                }
                            }
                        }
                    }
                }
                Ok((port_net, port_abs))
            })
            .collect();
        for r in results {
            let (net, abs) = r?;
            for k in 0..NCOMP {
                s.ledger.port_net[k] += net[k];
                s.ledger.port_abs[k] += abs[k];
            }
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
    fn add_sources(&self, g: &Grid, s: &mut Scratch, t: f64) {
        let prim = &s.prim;
        let partials: Vec<([f64; NCOMP], [f64; NCOMP])> = s
            .rate
            .par_iter_mut()
            .enumerate()
            .map(|(bi, rate_v)| {
                let mask = g.brick(bi).mask();
                let nt = g.brick(bi).n_theta();
                // COUP-2 ledger partials: the exact κV-weighted increments
                // applied here (geometric + wall-closure + external), net
                // and gross — per brick, combined in fixed brick order.
                let mut net = [0.0f64; NCOMP];
                let mut abs = [0.0f64; NCOMP];
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
                    let has_geom = g.brick(bi).has_geom();
                    for j in 0..nt {
                        let idx = j as usize * BRICK_CELLS + local;
                        let wc = prim[bi][idx];
                        let (rho, ur, ut, p) = (wc[I_RHO], wc[1], wc[2], wc[4]);
                        let rate = &mut rate_v[idx];
                        // Per-sector κV ledger weight (S9): each θ-plane's
                        // applied increment weighs by its own sector's gas
                        // volume. At N_θ = 1 plane 0 = the pre-S9 bits.
                        let kv = g.brick(bi).kappa_cell(j, local) * vol;
                        let s_mr = (a_out * p - a_in * p) / vol + rho * ut * ut * geo;
                        let s_mt = rho * ur * ut * geo;
                        rate[I_MR] += s_mr;
                        rate[I_MT] -= s_mt;
                        net[I_MR] += kv * s_mr;
                        abs[I_MR] += (kv * s_mr).abs();
                        net[I_MT] -= kv * s_mt;
                        abs[I_MT] += (kv * s_mt).abs();
                        // Embedded-interface pressure closure (cut cells
                        // only — gated on geometry presence so full-box
                        // worlds stay bit-identical; `+0.0` could flip a
                        // −0.0 rate bit). Per-sector (S9): the sector's own
                        // closure vector, SOLV-1 §3.3 / FND-2 §3.4.
                        if has_geom {
                            let (w_r, w_theta, w_z) = g.wall_closure_cell(i_r, j, i_z, nt);
                            let inv_kv = 1.0 / (g.brick(bi).kappa_cell(j, local) * vol);
                            let (wr_kv, wz_kv) = (w_r * inv_kv, w_z * inv_kv);
                            rate[I_MR] += p * wr_kv;
                            rate[I_MZ] += p * wz_kv;
                            net[I_MR] += kv * (p * wr_kv);
                            abs[I_MR] += (kv * (p * wr_kv)).abs();
                            net[I_MZ] += kv * (p * wz_kv);
                            abs[I_MZ] += (kv * (p * wz_kv)).abs();
                            // θ-limb (S9): applied only when nonzero — on
                            // the certified N_θ = 1 revolved path W_θ is
                            // exactly 0 (both θ-apertures the same bits),
                            // and adding its +0.0 could flip a −0.0 rate
                            // bit (the same reasoning as the has_geom
                            // gate). A pure-data comparison, not a mode.
                            if w_theta != 0.0 {
                                let wt_kv = w_theta * inv_kv;
                                rate[I_MT] += p * wt_kv;
                                net[I_MT] += kv * (p * wt_kv);
                                abs[I_MT] += (kv * (p * wt_kv)).abs();
                            }
                        }
                        let src = (self.source)(r, Grid::theta_center(j, nt), z, t);
                        for k in 0..NCOMP {
                            rate[k] += src[k];
                            net[k] += kv * src[k];
                            abs[k] += (kv * src[k]).abs();
                        }
                    }
                }
                (net, abs)
            })
            .collect();
        for (net, abs) in partials {
            for k in 0..NCOMP {
                s.ledger.src_net[k] += net[k];
                s.ledger.src_abs[k] += abs[k];
            }
        }
    }
}
