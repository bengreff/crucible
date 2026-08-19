//! The conduction operator: flux-form ∇·(k∇T) + S on the exact cylindrical
//! metric (FND-2 §3.2), one law for every cell (Rule 12 — no material or
//! regime branch anywhere; k and ρc_p arrive as data).
//!
//! Discretization: cell-centered finite volume, two-point flux per face —
//! heat into a cell through a face is `k·A_face·(T_neighbor − T_cell)/d`
//! with the true face area and center-to-center distance (θ distance is the
//! arc r̄·Δθ). Boundary Dirichlet faces use the half-cell distance `d/2`;
//! Neumann faces impose the outward heat flux directly. At `r = 0` the axis
//! face has zero area and drops out geometrically; at `N_θ = 1` the θ-flux
//! vanishes identically (the neighbor is the cell itself) — the same
//! operator is the axisymmetric operator, no special case.
//!
//! **Time integration (S2, COUP-3 §3.1): this operator is class `D` —
//! spatially-coupled implicit diffusion.** The explicit `step`/`advance`
//! scaffolding of sessions 5–12 is retired; the one production advance is
//! the SDC-IMEX step (`crate::sdc`), whose fixed-cycle CG solve applies
//! this module's [`Conduction::assemble_heat`] — the single owner of the
//! spatial discretization (affine in T at frozen operands). Gas↔solid wall
//! exchange enters as Robin interface data ([`GasFaceRobin`], linear in the
//! solid cell's T — COUP-2 §3.5's Robin-Robin placement); the coolant side
//! stays a [`FaceBc::Robin`] on exterior faces (COUP-7 owns that closure).
//! [`Conduction::stable_dt`] survives as the *explicit stability bound* —
//! a reference quantity the stiffness tests measure against, no longer a
//! step controller.
//!
//! Sweep structure (post-review): neighbor bricks are resolved ONCE per
//! brick (≤4 Morton lookups), and every in-brick access — the center, both
//! θ-neighbors, and all (r,z) neighbors of the 36/64 interior cells — is
//! direct index arithmetic.
//!
//! Determinism (§3.7): fixed Morton-brick / θ-plane / cell sweep order,
//! plain f64, and the per-cell accumulation order (r−, r+, θ−, θ+, z−, z+,
//! source) is part of the certified bit-identical behavior.

use std::collections::BTreeMap;

use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, FieldId, Grid};

#[derive(Debug, Clone, PartialEq)]
pub enum SolverError {
    /// This sweep requires one uniform N_θ across bricks; the N_θ-jump
    /// refluxing lands with the plan's 3-D wave (fail loud, no guess).
    MixedThetaResolution,
    NonFiniteState {
        i_r: usize,
        i_z: usize,
    },
    /// A face against a cell outside this operator's domain, with no
    /// interior-face treatment configured — refuse, never guess (the
    /// Goal-A full-box fixtures declare `InteriorFaces::refuse()`).
    UnhandledInteriorFace {
        i_r: usize,
        i_z: usize,
    },
    /// A non-physical coefficient (κ, ρc_p, or a Robin h must be finite
    /// and positive).
    BadCoefficient(&'static str),
}

impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MixedThetaResolution => write!(
                f,
                "mixed per-brick N_θ in one conduction sweep — flux aggregation across an \
                 N_θ jump arrives with the 3-D refluxing wave; refusing rather than guessing"
            ),
            Self::NonFiniteState { i_r, i_z } => {
                write!(
                    f,
                    "non-finite temperature at (i_r={i_r}, i_z={i_z}) — halt with diagnosis"
                )
            }
            Self::UnhandledInteriorFace { i_r, i_z } => {
                write!(
                    f,
                    "cell (i_r={i_r}, i_z={i_z}) has a face against a cell outside the \
                     operator's domain but no interior-face treatment is configured — \
                     refusing rather than guessing (COUP-2 §3.5)"
                )
            }
            Self::BadCoefficient(which) => {
                write!(
                    f,
                    "{which} must be finite and positive (fail loud, META-1 P6)"
                )
            }
        }
    }
}

impl std::error::Error for SolverError {}

/// Boundary condition on one grid face. Functions receive `(r, θ, z, t)` at
/// the face centroid — pure functions, so determinism is preserved.
/// (COUP-7 will replace these closures with declarative boundary-object
/// data; recorded as a tracked deferral in CLAUDE.md.)
pub enum FaceBc<'a> {
    /// Prescribed temperature at the face.
    Dirichlet(&'a dyn Fn(f64, f64, f64, f64) -> f64),
    /// Prescribed outward heat flux [W/m²]; 0.0 = insulated.
    HeatFlux(f64),
    /// Convective exchange against an ambient at `t_inf` through film
    /// coefficient `h` [W/(m²·K)], discretized as the film + half-cell
    /// conduction series resistance (the consistent 2nd-order face form;
    /// COUP-2 §3.5's coolant side until COUP-7's jacket object lands).
    Robin { h: f64, t_inf: f64 },
}

/// COUP-2 §3.5 — one gas-facing solid face's Robin interface data, as the
/// class-`D` implicit solve consumes it: heat INTO the solid cell is
/// `area_scale · A_face · (t_aw − T_cell) / (1/h_film + half_d/κ_solid)` —
/// **linear in the solid cell's T** (the wall-function `h` is the Robin
/// coefficient), with the gas-side operands (h_film, T_aw) frozen per
/// Picard sweep by the SDC orchestrator. `area_scale` maps stair-face area
/// onto the embedded-interface area |W| on cut worlds (1.0 on box worlds).
#[derive(Debug, Clone, Copy)]
pub struct GasFaceRobin {
    pub h_film: f64,
    pub t_aw: f64,
    pub area_scale: f64,
}

/// Key of a gas-facing solid face in the exchange map: the SOLID cell's
/// (i_r, i_z) and the face direction *as seen from the solid cell*
/// (`FaceDir::index()`).
pub type ExchangeKey = (usize, usize, u8);

/// What a face against a cell *outside the operator's domain* does — the
/// interior counterpart of `Bcs` (domain edges). Gas-exchange faces are no
/// longer configured here: they arrive per-solve as [`GasFaceRobin`] data
/// through [`Conduction::assemble_heat`] (COUP-2 §3.5 sweep placement,
/// owned by the SDC orchestrator).
pub struct InteriorFaces<'a> {
    /// Face against an exterior cell: `None` ⇒ hard error; `Some` ⇒ the BC
    /// (insulated `HeatFlux(0.0)`, coolant `Robin`, …). Position-dependent
    /// `Dirichlet` closures see the face centroid as usual.
    pub exterior: Option<FaceBc<'a>>,
}

impl InteriorFaces<'_> {
    /// The full-box declaration: any interior face is a loud error (the
    /// Goal-A fixtures — their domain has no interior boundary at all).
    pub fn refuse() -> Self {
        InteriorFaces { exterior: None }
    }
}

/// Which region mask is this operator's domain (data from the config-time
/// region classification — never an `if(material)` in the sweep).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Domain {
    FlowActive,
    Solid,
}

pub struct Bcs<'a> {
    /// Ignored when `r_min = 0` (the axis face has zero area).
    pub r_inner: FaceBc<'a>,
    pub r_outer: FaceBc<'a>,
    pub z_lo: FaceBc<'a>,
    pub z_hi: FaceBc<'a>,
}

/// The operator: uniform conductivity k [W/(m·K)] and volumetric heat
/// capacity ρc_p [J/(m³·K)] as pure data (per-cell spine-derived k arrives
/// with FND-7), a volumetric source S(r, θ, z, t) [W/m³], and the face BCs.
pub struct Conduction<'a> {
    pub kappa: f64,
    pub rho_cp: f64,
    pub source: &'a dyn Fn(f64, f64, f64, f64) -> f64,
    pub bcs: Bcs<'a>,
    /// The region mask this operator sweeps (config-time data).
    pub domain: Domain,
    /// Treatment of faces against cells outside that domain.
    pub interior: InteriorFaces<'a>,
}

/// Assembly mode of [`Conduction::assemble_heat`]. The operator is affine
/// in T at frozen operands: `heat(T) = Ã·T + b̃`. `Affine` evaluates the
/// full form (physics); `Linear` evaluates `Ã·T` exactly — every constant
/// term (Dirichlet values, imposed fluxes, Robin ambients, exchange T_aw,
/// the source) is dropped, so the CG matrix apply carries no cancellation
/// error against b̃.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssembleMode {
    Affine,
    Linear,
}

/// COUP-2 ledger lines of one assembly: the port subset (everything that
/// does NOT telescope — domain-edge BC heats, exterior-face heats, the
/// gas-exchange heats) plus the volumetric source, and the gross magnitude
/// scale (Σ|every term summed|, the §3.1.1 S[q] throughput ingredient).
/// All in W at the assembled state.
#[derive(Debug, Clone, Copy, Default)]
pub struct HeatLedger {
    pub bc_w: f64,
    pub exterior_w: f64,
    pub exchange_w: f64,
    pub source_w: f64,
    pub gross_w: f64,
}

impl HeatLedger {
    /// Total non-telescoping (port + source) heat rate [W] — what the
    /// COUP-2 identity charges against the solid's stored-energy change.
    pub fn applied_w(&self) -> f64 {
        self.bc_w + self.exterior_w + self.exchange_w + self.source_w
    }
}

/// Where a face's neighbor value comes from: same brick (index offset), a
/// specific adjacent brick (resolved once per brick), an in-bounds position
/// whose brick is not allocated (exterior), or the grid edge.
#[derive(Clone, Copy)]
enum Nbr {
    InBrick(isize),
    Cross { bi: usize, local: usize },
    Missing,
    Edge,
}

impl Conduction<'_> {
    /// The one spatial discretization, as heat into each domain cell [W]
    /// (module doc): `out[bi][idx] = Σ_faces heat + S·V`, affine in the
    /// `t_field` values at frozen operands. Consumers: the class-`D`
    /// fixed-cycle CG solve (`crate::sdc`) applies it in `Linear` mode per
    /// iteration and `Affine` mode for right-hand sides and ledgers.
    ///
    /// - `exchange`: gas-facing solid faces' Robin data (Solid domain only;
    ///   a gas face with no entry / no map refuses loudly).
    /// - `diag`: when present, receives `∂heat_i/∂T_i` (≤ 0) per cell — the
    ///   Jacobi preconditioner ingredient (operand-frozen, T-independent).
    /// - `ledger`: when present, accumulates the COUP-2 lines (fixed
    ///   accumulation order — part of the deterministic contract).
    /// - `exchange_heats`: when present, records each gas-face heat term
    ///   [W] in visit order keyed by [`ExchangeKey`] — the numbers the gas
    ///   side must debit exactly (conservation by construction).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn assemble_heat(
        &self,
        g: &Grid,
        t_field: FieldId,
        time: f64,
        mode: AssembleMode,
        exchange: Option<&BTreeMap<ExchangeKey, GasFaceRobin>>,
        out: &mut [Vec<f64>],
        mut diag: Option<&mut [Vec<f64>]>,
        mut ledger: Option<&mut HeatLedger>,
        mut exchange_heats: Option<&mut Vec<(ExchangeKey, f64)>>,
    ) -> Result<(), SolverError> {
        self.validate()?;
        let nt = g.brick(0).n_theta();
        if g.bricks().iter().any(|b| b.n_theta() != nt) {
            return Err(SolverError::MixedThetaResolution);
        }
        let dtheta = std::f64::consts::TAU / f64::from(nt);
        let (dr, dz) = (g.spec().dr, g.spec().dz);
        let (n_r, n_z) = (g.spec().n_r, g.spec().n_z);
        let z0 = g.spec().z_min;
        let r0 = g.spec().r_min;
        let linear = mode == AssembleMode::Linear;

        for bi in 0..g.n_bricks() {
            let (br, bz, dmask, gmask) = {
                let b = g.brick(bi);
                let dm = match self.domain {
                    Domain::FlowActive => b.mask(),
                    Domain::Solid => b.solid_mask(),
                };
                (b.br(), b.bz(), dm, b.mask())
            };
            out[bi].fill(0.0);
            if let Some(d) = diag.as_deref_mut() {
                d[bi].fill(0.0);
            }
            if dmask == 0 {
                continue;
            }
            // Adjacent bricks, resolved once per brick (≤4 Morton lookups).
            let nb_rm = (br > 0)
                .then(|| g.brick_index_by_coords(br - 1, bz))
                .flatten();
            let nb_rp = g.brick_index_by_coords(br + 1, bz);
            let nb_zm = (bz > 0)
                .then(|| g.brick_index_by_coords(br, bz - 1))
                .flatten();
            let nb_zp = g.brick_index_by_coords(br, bz + 1);

            for j in 0..nt {
                let theta = (f64::from(j) + 0.5) * dtheta;
                for local in 0..BRICK_CELLS {
                    if dmask & (1u64 << local) == 0 {
                        continue;
                    }
                    let (lr, lz) = (local / BRICK, local % BRICK);
                    let (i_r, i_z) = (br as usize * BRICK + lr, bz as usize * BRICK + lz);
                    let rbar = r0 + (i_r as f64 + 0.5) * dr;
                    let zbar = z0 + (i_z as f64 + 0.5) * dz;
                    let vol = g.cell_volume(i_r, nt);

                    let here = g.brick(bi);
                    let t_here = here.field(t_field);
                    let idx = here.cell_index(j, local);
                    let t_c = t_here[idx];
                    if !t_c.is_finite() {
                        return Err(SolverError::NonFiniteState { i_r, i_z });
                    }

                    // Face neighbor resolution — pure index arithmetic. A
                    // `None` neighbor brick is an unallocated (fully
                    // exterior) region: the face is an interior boundary.
                    let n_rm = if i_r == 0 {
                        Nbr::Edge
                    } else if lr > 0 {
                        Nbr::InBrick(-(BRICK as isize))
                    } else {
                        match nb_rm {
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
                        match nb_rp {
                            Some(bi) => Nbr::Cross { bi, local: lz },
                            None => Nbr::Missing,
                        }
                    };
                    let n_zm = if i_z == 0 {
                        Nbr::Edge
                    } else if lz > 0 {
                        Nbr::InBrick(-1)
                    } else {
                        match nb_zm {
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
                        match nb_zp {
                            Some(bi) => Nbr::Cross {
                                bi,
                                local: lr * BRICK,
                            },
                            None => Nbr::Missing,
                        }
                    };

                    let mut heat_in = 0.0f64; // W
                    let mut diag_c = 0.0f64; // ∂heat_in/∂t_c [W/K]
                    let mut gross = 0.0f64;

                    // Per non-edge face: in-domain neighbor ⇒ two-point
                    // conductive flux; out-of-domain ⇒ the interior-face
                    // treatment (gas-exchange Robin or exterior BC),
                    // classified from the neighbor's masks — data, not
                    // `if(material)`.
                    macro_rules! interior_face {
                        ($n:expr, $area:expr, $dist:expr, $dir:expr, $pos:expr) => {{
                            let (nbr_t, in_domain, is_gas) = match $n {
                                Nbr::InBrick(off) => {
                                    let nl = (local as isize + off) as usize;
                                    (
                                        t_here[(idx as isize + off) as usize],
                                        dmask & (1u64 << nl) != 0,
                                        gmask & (1u64 << nl) != 0,
                                    )
                                }
                                Nbr::Cross { bi: nbi, local: nl } => {
                                    let nb = g.brick(nbi);
                                    let ndm = match self.domain {
                                        Domain::FlowActive => nb.mask(),
                                        Domain::Solid => nb.solid_mask(),
                                    };
                                    (
                                        nb.field(t_field)[nb.cell_index(j, nl)],
                                        ndm & (1u64 << nl) != 0,
                                        nb.mask() & (1u64 << nl) != 0,
                                    )
                                }
                                Nbr::Missing => (f64::NAN, false, false),
                                Nbr::Edge => unreachable!("edge faces take the BC path"),
                            };
                            if in_domain {
                                let term = self.kappa * $area * (nbr_t - t_c) / $dist;
                                heat_in += term;
                                diag_c -= self.kappa * $area / $dist;
                                gross += term.abs();
                            } else if is_gas && self.domain == Domain::Solid {
                                // Gas-exchange Robin face (COUP-2 §3.5):
                                // linear in t_c at frozen gas operands.
                                let Some(map) = exchange else {
                                    return Err(SolverError::UnhandledInteriorFace { i_r, i_z });
                                };
                                let key: ExchangeKey = (i_r, i_z, $dir.index() as u8);
                                let Some(gr) = map.get(&key) else {
                                    return Err(SolverError::UnhandledInteriorFace { i_r, i_z });
                                };
                                let resist = 1.0 / gr.h_film + 0.5 * $dist / self.kappa;
                                let t_drive = if linear { 0.0 } else { gr.t_aw };
                                let term = gr.area_scale * $area * (t_drive - t_c) / resist;
                                heat_in += term;
                                diag_c -= gr.area_scale * $area / resist;
                                gross += term.abs();
                                if let Some(rec) = exchange_heats.as_deref_mut() {
                                    rec.push((key, term));
                                }
                                if let Some(l) = ledger.as_deref_mut() {
                                    l.exchange_w += term;
                                }
                            } else {
                                match &self.interior.exterior {
                                    Some(bc) => {
                                        let (term, d) = self.face_bc_heat(
                                            bc,
                                            $area,
                                            t_c,
                                            0.5 * $dist,
                                            $pos,
                                            linear,
                                        );
                                        heat_in += term;
                                        diag_c += d;
                                        gross += term.abs();
                                        if let Some(l) = ledger.as_deref_mut() {
                                            l.exterior_w += term;
                                        }
                                    }
                                    None => {
                                        return Err(SolverError::UnhandledInteriorFace {
                                            i_r,
                                            i_z,
                                        });
                                    }
                                }
                            }
                        }};
                    }
                    macro_rules! edge_face {
                        ($bc:expr, $area:expr, $half_d:expr, $pos:expr) => {{
                            let (term, d) =
                                self.face_bc_heat($bc, $area, t_c, $half_d, $pos, linear);
                            heat_in += term;
                            diag_c += d;
                            gross += term.abs();
                            if let Some(l) = ledger.as_deref_mut() {
                                l.bc_w += term;
                            }
                        }};
                    }

                    // r− face. a_in == 0.0 at the r = 0 axis: drops out.
                    let a_in = g.face_area_r(i_r, false, nt);
                    match n_rm {
                        Nbr::Edge => {
                            if a_in > 0.0 {
                                edge_face!(
                                    &self.bcs.r_inner,
                                    a_in,
                                    0.5 * dr,
                                    (r0, theta, zbar, time)
                                );
                            }
                        }
                        n => interior_face!(
                            n,
                            a_in,
                            dr,
                            FaceDir::RMinus,
                            (r0 + i_r as f64 * dr, theta, zbar, time)
                        ),
                    }

                    // r+ face.
                    let a_out = g.face_area_r(i_r, true, nt);
                    match n_rp {
                        Nbr::Edge => {
                            edge_face!(
                                &self.bcs.r_outer,
                                a_out,
                                0.5 * dr,
                                (r0 + n_r as f64 * dr, theta, zbar, time)
                            );
                        }
                        n => interior_face!(
                            n,
                            a_out,
                            dr,
                            FaceDir::RPlus,
                            (r0 + (i_r + 1) as f64 * dr, theta, zbar, time)
                        ),
                    }

                    // θ faces: periodic, arc distance r̄·Δθ, always in-brick.
                    // At N_θ = 1 the neighbor is the cell itself ⇒ flux ≡ 0
                    // (and ∂/∂t_c ≡ 0: both terms carry the cell's own T).
                    let a_th = g.face_area_theta();
                    let arc = rbar * dtheta;
                    let jm = (j + nt - 1) % nt;
                    let jp = (j + 1) % nt;
                    let th_m = self.kappa * a_th * (t_here[here.cell_index(jm, local)] - t_c) / arc;
                    let th_p = self.kappa * a_th * (t_here[here.cell_index(jp, local)] - t_c) / arc;
                    heat_in += th_m;
                    heat_in += th_p;
                    gross += th_m.abs() + th_p.abs();
                    if jm != j {
                        diag_c -= 2.0 * self.kappa * a_th / arc;
                    }

                    // z faces.
                    let a_z = g.face_area_z(i_r, nt);
                    match n_zm {
                        Nbr::Edge => {
                            edge_face!(&self.bcs.z_lo, a_z, 0.5 * dz, (rbar, theta, z0, time));
                        }
                        n => interior_face!(
                            n,
                            a_z,
                            dz,
                            FaceDir::ZMinus,
                            (rbar, theta, z0 + i_z as f64 * dz, time)
                        ),
                    }
                    match n_zp {
                        Nbr::Edge => {
                            edge_face!(
                                &self.bcs.z_hi,
                                a_z,
                                0.5 * dz,
                                (rbar, theta, z0 + n_z as f64 * dz, time)
                            );
                        }
                        n => interior_face!(
                            n,
                            a_z,
                            dz,
                            FaceDir::ZPlus,
                            (rbar, theta, z0 + (i_z + 1) as f64 * dz, time)
                        ),
                    }

                    let mut cell_heat = heat_in;
                    if !linear {
                        let sv = (self.source)(rbar, theta, zbar, time) * vol;
                        cell_heat += sv;
                        gross += sv.abs();
                        if let Some(l) = ledger.as_deref_mut() {
                            l.source_w += sv;
                        }
                    }
                    out[bi][idx] = cell_heat;
                    if let Some(d) = diag.as_deref_mut() {
                        d[bi][idx] = diag_c;
                    }
                    if let Some(l) = ledger.as_deref_mut() {
                        l.gross_w += gross;
                    }
                }
            }
        }
        Ok(())
    }

    /// Explicit-stability step bound `dt ≤ C·ρc_p/(k·Σ 2/d_i²)` with the
    /// smallest distances on the grid (θ arc at the innermost ring). No
    /// longer a step controller (class `D` is implicit — COUP-3 §3.1): this
    /// is the reference bound the stiffness tests measure the implicit
    /// solve against, and a diagnostic scale.
    pub fn stable_dt(&self, g: &Grid, safety: f64) -> f64 {
        let nt = g.brick(0).n_theta();
        let dtheta = std::f64::consts::TAU / f64::from(nt);
        let arc_min = (g.spec().r_min + 0.5 * g.spec().dr) * dtheta;
        let mut inv = 2.0 / (g.spec().dr * g.spec().dr) + 2.0 / (g.spec().dz * g.spec().dz);
        if nt > 1 {
            inv += 2.0 / (arc_min * arc_min);
        }
        safety * self.rho_cp / (self.kappa * inv)
    }

    /// One face-BC heat term and its `∂/∂T_cell` [W/K]. In `linear` mode
    /// the constant part (prescribed value/flux/ambient) is dropped exactly.
    #[inline]
    fn face_bc_heat(
        &self,
        bc: &FaceBc<'_>,
        area: f64,
        t_c: f64,
        half_d: f64,
        pos: (f64, f64, f64, f64),
        linear: bool,
    ) -> (f64, f64) {
        match bc {
            FaceBc::Dirichlet(f) => {
                let t_face = if linear {
                    0.0
                } else {
                    f(pos.0, pos.1, pos.2, pos.3)
                };
                let coeff = self.kappa * area / half_d;
                (coeff * (t_face - t_c), -coeff)
            }
            FaceBc::HeatFlux(q_out) => {
                let term = if linear { 0.0 } else { -q_out * area };
                (term, 0.0)
            }
            // Film + half-cell conduction in series — the consistent face
            // form; bounded above by the Dirichlet coefficient.
            FaceBc::Robin { h, t_inf } => {
                let drive = if linear { 0.0 } else { *t_inf };
                let coeff = area / (1.0 / h + half_d / self.kappa);
                (coeff * (drive - t_c), -coeff)
            }
        }
    }

    /// Fail-loud coefficient checks, once per assembly (META-1 P6).
    fn validate(&self) -> Result<(), SolverError> {
        if !(self.kappa.is_finite() && self.kappa > 0.0) {
            return Err(SolverError::BadCoefficient("kappa"));
        }
        if !(self.rho_cp.is_finite() && self.rho_cp > 0.0) {
            return Err(SolverError::BadCoefficient("rho_cp"));
        }
        let robin_ok = |bc: &FaceBc<'_>| match bc {
            FaceBc::Robin { h, t_inf } => h.is_finite() && *h > 0.0 && t_inf.is_finite(),
            _ => true,
        };
        for bc in [
            &self.bcs.r_inner,
            &self.bcs.r_outer,
            &self.bcs.z_lo,
            &self.bcs.z_hi,
        ] {
            if !robin_ok(bc) {
                return Err(SolverError::BadCoefficient(
                    "Robin (h, T∞) on a domain edge",
                ));
            }
        }
        if let Some(bc) = &self.interior.exterior
            && !robin_ok(bc)
        {
            return Err(SolverError::BadCoefficient(
                "Robin (h, T∞) on an exterior face",
            ));
        }
        Ok(())
    }
}
