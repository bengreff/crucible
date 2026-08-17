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
//! Sweep structure (post-review): neighbor bricks are resolved ONCE per
//! brick (≤4 Morton lookups), and every in-brick access — the center, both
//! θ-neighbors, and all (r,z) neighbors of the 36/64 interior cells — is
//! direct index arithmetic. The v1 sweep did up to 8 binary searches per
//! cell per step, exactly the pointer-chasing FND-2 §3.9 forbids; this
//! sweep is the template the flow solver copies, so the pattern matters
//! more than this operator's own cost.
//!
//! Determinism (§3.7): fixed Morton-brick / θ-plane / cell sweep order,
//! two-pass rate-then-update (Jacobi form), plain f64, and the per-cell
//! accumulation order (r−, r+, θ−, θ+, z−, z+, source) is part of the
//! certified bit-identical behavior.

use crucible_grid::{BRICK, BRICK_CELLS, FaceDir, FieldId, Grid};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolverError {
    /// This session's sweep requires one uniform N_θ across bricks; the
    /// N_θ-jump refluxing lands with COUP-2/COUP-3 (fail loud, no guess).
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
                 N_θ jump arrives with COUP-2/COUP-3; refusing rather than guessing"
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

/// What a face against a cell *outside the operator's domain* does — the
/// interior counterpart of `Bcs` (domain edges), COUP-2 §3.5's interface
/// delineation. Both fields are explicit at every construction site (no
/// hidden defaults, FND-4 §3.5).
pub struct InteriorFaces<'a> {
    /// Face against a flow-active (gas) cell: heat INTO this domain cell
    /// [W/m²], from the coupler's single per-face evaluation (conservation
    /// by construction — the gas side applies the same number negated).
    /// `None` ⇒ such a face is a hard error.
    pub gas: Option<&'a dyn Fn(usize, usize, u32, FaceDir) -> f64>,
    /// Face against an exterior cell: `None` ⇒ hard error; `Some` ⇒ the BC
    /// (insulated `HeatFlux(0.0)`, coolant `Robin`, …). Position-dependent
    /// `Dirichlet` closures see the face centroid as usual.
    pub exterior: Option<FaceBc<'a>>,
}

impl InteriorFaces<'_> {
    /// The full-box declaration: any interior face is a loud error (the
    /// Goal-A fixtures — their domain has no interior boundary at all).
    pub fn refuse() -> Self {
        InteriorFaces {
            gas: None,
            exterior: None,
        }
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
    /// One explicit step: fills `rate` with dT/dt at time `t`, then applies
    /// `T += dt·rate`. Fixed traversal order; bit-reproducible.
    pub fn step(
        &self,
        g: &mut Grid,
        t_field: FieldId,
        rate: FieldId,
        t: f64,
        dt: f64,
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

        // Reused rate buffer: rates are staged per brick, then written once
        // (keeps reads of neighbor bricks and the write disjoint).
        let mut buf = vec![0.0f64; nt as usize * BRICK_CELLS];

        // Pass 1: rates, brick by brick in Morton order.
        for bi in 0..g.n_bricks() {
            let (br, bz, dmask, gmask) = {
                let b = g.brick(bi);
                let dm = match self.domain {
                    Domain::FlowActive => b.mask(),
                    Domain::Solid => b.solid_mask(),
                };
                (b.br(), b.bz(), dm, b.mask())
            };
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

                    // Per non-edge face: in-domain neighbor ⇒ two-point
                    // conductive flux; out-of-domain ⇒ the interior-face
                    // treatment (gas exchange or exterior BC), classified
                    // from the neighbor's masks — data, not `if(material)`.
                    let face = |n: Nbr,
                                area: f64,
                                dist: f64,
                                dir: FaceDir,
                                pos: (f64, f64, f64, f64)|
                     -> Result<f64, SolverError> {
                        let (nbr_t, in_domain, is_gas) = match n {
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
                            return Ok(self.kappa * area * (nbr_t - t_c) / dist);
                        }
                        // Interior boundary. Gas side ⇒ coupler flux; the
                        // domain being FlowActive makes a solid/exterior
                        // neighbor exterior-like by the same rule.
                        if is_gas && self.domain == Domain::Solid {
                            match self.interior.gas {
                                Some(q_in) => Ok(q_in(i_r, i_z, j, dir) * area),
                                None => Err(SolverError::UnhandledInteriorFace { i_r, i_z }),
                            }
                        } else {
                            match &self.interior.exterior {
                                Some(bc) => Ok(self.face_bc_heat(bc, area, t_c, 0.5 * dist, pos)),
                                None => Err(SolverError::UnhandledInteriorFace { i_r, i_z }),
                            }
                        }
                    };

                    let mut heat_in = 0.0f64; // W

                    // r− face. a_in == 0.0 at the r = 0 axis: drops out.
                    let a_in = g.face_area_r(i_r, false, nt);
                    match n_rm {
                        Nbr::Edge => {
                            if a_in > 0.0 {
                                heat_in += self.face_bc_heat(
                                    &self.bcs.r_inner,
                                    a_in,
                                    t_c,
                                    0.5 * dr,
                                    (r0, theta, zbar, t),
                                );
                            }
                        }
                        n => {
                            heat_in += face(
                                n,
                                a_in,
                                dr,
                                FaceDir::RMinus,
                                (r0 + i_r as f64 * dr, theta, zbar, t),
                            )?;
                        }
                    }

                    // r+ face.
                    let a_out = g.face_area_r(i_r, true, nt);
                    match n_rp {
                        Nbr::Edge => {
                            heat_in += self.face_bc_heat(
                                &self.bcs.r_outer,
                                a_out,
                                t_c,
                                0.5 * dr,
                                (r0 + n_r as f64 * dr, theta, zbar, t),
                            );
                        }
                        n => {
                            heat_in += face(
                                n,
                                a_out,
                                dr,
                                FaceDir::RPlus,
                                (r0 + (i_r + 1) as f64 * dr, theta, zbar, t),
                            )?;
                        }
                    }

                    // θ faces: periodic, arc distance r̄·Δθ, always in-brick.
                    // At N_θ = 1 the neighbor is the cell itself ⇒ flux 0.
                    let a_th = g.face_area_theta();
                    let arc = rbar * dtheta;
                    let jm = (j + nt - 1) % nt;
                    let jp = (j + 1) % nt;
                    heat_in += self.kappa * a_th * (t_here[here.cell_index(jm, local)] - t_c) / arc;
                    heat_in += self.kappa * a_th * (t_here[here.cell_index(jp, local)] - t_c) / arc;

                    // z faces.
                    let a_z = g.face_area_z(i_r, nt);
                    match n_zm {
                        Nbr::Edge => {
                            heat_in += self.face_bc_heat(
                                &self.bcs.z_lo,
                                a_z,
                                t_c,
                                0.5 * dz,
                                (rbar, theta, z0, t),
                            );
                        }
                        n => {
                            heat_in += face(
                                n,
                                a_z,
                                dz,
                                FaceDir::ZMinus,
                                (rbar, theta, z0 + i_z as f64 * dz, t),
                            )?;
                        }
                    }
                    match n_zp {
                        Nbr::Edge => {
                            heat_in += self.face_bc_heat(
                                &self.bcs.z_hi,
                                a_z,
                                t_c,
                                0.5 * dz,
                                (rbar, theta, z0 + n_z as f64 * dz, t),
                            );
                        }
                        n => {
                            heat_in += face(
                                n,
                                a_z,
                                dz,
                                FaceDir::ZPlus,
                                (rbar, theta, z0 + (i_z + 1) as f64 * dz, t),
                            )?;
                        }
                    }

                    let s = (self.source)(rbar, theta, zbar, t);
                    buf[idx] = (heat_in + s * vol) / (self.rho_cp * vol);
                }
            }
            let len = nt as usize * BRICK_CELLS;
            g.brick_field_mut(bi, rate)[..len].copy_from_slice(&buf[..len]);
        }

        // Pass 2: apply.
        for bi in 0..g.n_bricks() {
            let (t_slice, r_slice) = g.brick_fields_mut2(bi, t_field, rate);
            for (tv, rv) in t_slice.iter_mut().zip(r_slice.iter()) {
                *tv += dt * rv;
            }
        }
        Ok(())
    }

    /// March `n_steps` of size `dt` from `t0`; returns the final time.
    pub fn advance(
        &self,
        g: &mut Grid,
        t_field: FieldId,
        rate: FieldId,
        t0: f64,
        dt: f64,
        n_steps: usize,
    ) -> Result<f64, SolverError> {
        let mut t = t0;
        for _ in 0..n_steps {
            self.step(g, t_field, rate, t, dt)?;
            t += dt;
        }
        Ok(t)
    }

    /// Explicit-stability step bound `dt ≤ C·ρc_p/(k·Σ 2/d_i²)` with the
    /// smallest distances on the grid (θ arc at the innermost ring). The
    /// bound is spectrally sharp including Dirichlet boundaries (a boundary
    /// face adds to the diagonal but has no off-diagonal partner, leaving
    /// the Gershgorin radius unchanged — review-verified, with dt at
    /// 0.999× stable and 1.02× divergent).
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

    #[inline]
    fn face_bc_heat(
        &self,
        bc: &FaceBc<'_>,
        area: f64,
        t_c: f64,
        half_d: f64,
        pos: (f64, f64, f64, f64),
    ) -> f64 {
        match bc {
            FaceBc::Dirichlet(f) => {
                self.kappa * area * (f(pos.0, pos.1, pos.2, pos.3) - t_c) / half_d
            }
            FaceBc::HeatFlux(q_out) => -q_out * area,
            // Film + half-cell conduction in series — the consistent face
            // form; bounded above by the Dirichlet coefficient, so the
            // `stable_dt` Gershgorin bound continues to cover it.
            FaceBc::Robin { h, t_inf } => area * (t_inf - t_c) / (1.0 / h + half_d / self.kappa),
        }
    }

    /// Fail-loud coefficient checks, once per step (META-1 P6).
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
