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
//! Determinism (§3.7): fixed Morton-brick / θ-plane / cell sweep order,
//! two-pass rate-then-update (Jacobi form), plain f64.

use crucible_grid::{BRICK, BRICK_CELLS, FieldId, Grid};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SolverError {
    /// This session's sweep requires one uniform N_θ across bricks; the
    /// N_θ-jump refluxing lands with COUP-2/COUP-3 (fail loud, no guess).
    MixedThetaResolution,
    NonFiniteState {
        i_r: usize,
        i_z: usize,
    },
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
        }
    }
}

impl std::error::Error for SolverError {}

/// Boundary condition on one grid face. Functions receive `(r, θ, z, t)` at
/// the face centroid — pure functions, so determinism is preserved.
pub enum FaceBc<'a> {
    /// Prescribed temperature at the face.
    Dirichlet(&'a dyn Fn(f64, f64, f64, f64) -> f64),
    /// Prescribed outward heat flux [W/m²]; 0.0 = insulated.
    HeatFlux(f64),
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
        let nt = g.bricks[0].n_theta;
        if g.bricks.iter().any(|b| b.n_theta != nt) {
            return Err(SolverError::MixedThetaResolution);
        }
        let dtheta = std::f64::consts::TAU / f64::from(nt);
        let (dr, dz) = (g.spec.dr, g.spec.dz);
        let (n_r, n_z) = (g.spec.n_r, g.spec.n_z);
        let z0 = g.spec.z_min;
        let r0 = g.spec.r_min;

        // Pass 1: rates, brick by brick in Morton order.
        for bi in 0..g.bricks.len() {
            let (br, bz, mask) = (g.bricks[bi].br, g.bricks[bi].bz, g.bricks[bi].mask);
            for j in 0..nt {
                let theta = (f64::from(j) + 0.5) * dtheta;
                for local in 0..BRICK_CELLS {
                    if mask & (1u64 << local) == 0 {
                        continue;
                    }
                    let i_r = br as usize * BRICK + local / BRICK;
                    let i_z = bz as usize * BRICK + local % BRICK;
                    let rbar = r0 + (i_r as f64 + 0.5) * dr;
                    let zbar = z0 + (i_z as f64 + 0.5) * dz;
                    let vol = g.cell_volume(i_r, nt);
                    let t_c = self.read(g, t_field, i_r, j, i_z);
                    if !t_c.is_finite() {
                        return Err(SolverError::NonFiniteState { i_r, i_z });
                    }

                    let mut heat_in = 0.0f64; // W

                    // r− face.
                    let a_in = g.face_area_r(i_r, false, nt);
                    if i_r > 0 {
                        heat_in +=
                            self.kappa * a_in * (self.read(g, t_field, i_r - 1, j, i_z) - t_c) / dr;
                    } else if a_in > 0.0 {
                        heat_in += self.face_bc_heat(
                            &self.bcs.r_inner,
                            a_in,
                            t_c,
                            0.5 * dr,
                            (r0, theta, zbar, t),
                        );
                    } // a_in == 0.0: the r = 0 axis face — drops out geometrically.

                    // r+ face.
                    let a_out = g.face_area_r(i_r, true, nt);
                    if i_r + 1 < n_r {
                        heat_in +=
                            self.kappa * a_out * (self.read(g, t_field, i_r + 1, j, i_z) - t_c)
                                / dr;
                    } else {
                        heat_in += self.face_bc_heat(
                            &self.bcs.r_outer,
                            a_out,
                            t_c,
                            0.5 * dr,
                            (r0 + n_r as f64 * dr, theta, zbar, t),
                        );
                    }

                    // θ faces: periodic, arc distance r̄·Δθ. At N_θ = 1 the
                    // neighbor is the cell itself ⇒ flux exactly zero.
                    let a_th = g.face_area_theta();
                    let arc = rbar * dtheta;
                    let jm = (j + nt - 1) % nt;
                    let jp = (j + 1) % nt;
                    heat_in +=
                        self.kappa * a_th * (self.read(g, t_field, i_r, jm, i_z) - t_c) / arc;
                    heat_in +=
                        self.kappa * a_th * (self.read(g, t_field, i_r, jp, i_z) - t_c) / arc;

                    // z faces.
                    let a_z = g.face_area_z(i_r, nt);
                    if i_z > 0 {
                        heat_in +=
                            self.kappa * a_z * (self.read(g, t_field, i_r, j, i_z - 1) - t_c) / dz;
                    } else {
                        heat_in += self.face_bc_heat(
                            &self.bcs.z_lo,
                            a_z,
                            t_c,
                            0.5 * dz,
                            (rbar, theta, z0, t),
                        );
                    }
                    if i_z + 1 < n_z {
                        heat_in +=
                            self.kappa * a_z * (self.read(g, t_field, i_r, j, i_z + 1) - t_c) / dz;
                    } else {
                        heat_in += self.face_bc_heat(
                            &self.bcs.z_hi,
                            a_z,
                            t_c,
                            0.5 * dz,
                            (rbar, theta, z0 + n_z as f64 * dz, t),
                        );
                    }

                    let s = (self.source)(rbar, theta, zbar, t);
                    let rate_val = (heat_in + s * vol) / (self.rho_cp * vol);
                    let bi_w = g
                        .brick_index(i_r, i_z)
                        .expect("cell exists by construction");
                    let idx = g.bricks[bi_w].cell_index(j, crucible_grid::Grid::local_rz(i_r, i_z));
                    g.bricks[bi_w].field_mut(rate)[idx] = rate_val;
                }
            }
        }

        // Pass 2: apply.
        for b in &mut g.bricks {
            let n = b.n_theta as usize * BRICK_CELLS;
            for i in 0..n {
                let r = b.field(rate)[i];
                b.field_mut(t_field)[i] += dt * r;
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
    /// smallest distances on the grid (θ arc at the innermost ring).
    pub fn stable_dt(&self, g: &Grid, safety: f64) -> f64 {
        let nt = g.bricks[0].n_theta;
        let dtheta = std::f64::consts::TAU / f64::from(nt);
        let arc_min = (g.spec.r_min + 0.5 * g.spec.dr) * dtheta;
        let mut inv = 2.0 / (g.spec.dr * g.spec.dr) + 2.0 / (g.spec.dz * g.spec.dz);
        if nt > 1 {
            inv += 2.0 / (arc_min * arc_min);
        }
        safety * self.rho_cp / (self.kappa * inv)
    }

    #[inline]
    fn read(&self, g: &Grid, f: FieldId, i_r: usize, j: u32, i_z: usize) -> f64 {
        let bi = g
            .brick_index(i_r, i_z)
            .expect("neighbor within active grid");
        let b = &g.bricks[bi];
        b.field(f)[b.cell_index(j, crucible_grid::Grid::local_rz(i_r, i_z))]
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
        }
    }
}
