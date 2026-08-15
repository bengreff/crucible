//! FND-2 §3.4 — adaptive azimuthal resolution: conservative θ-coarsening /
//! refinement (the same operator runs unchanged at every N_θ; these are pure
//! conservative projections of stored fields), the thermalized-ΔKE ledger
//! entry per collapse event, the normalized azimuthal-variance symmetry
//! indicator, and the collapse/expand decision machinery (guard, hysteresis,
//! dwell). All loops are fixed-order (§3.7).

use crate::{BRICK, BRICK_CELLS, FieldId, Grid, GridError, N_THETA_GUARD};

/// §3.4 named defaults (config-documented; rationale in the doc):
/// τ_collapse sits orders below any declared physics band; τ_expand gives
/// 10² hysteresis against thrashing.
pub const TAU_COLLAPSE: f64 = 1e-6;
pub const TAU_EXPAND: f64 = 1e-4;
/// Indicator evaluation cadence in steps, and consecutive over/under
/// evaluations required before any change.
pub const N_SYM_CADENCE: u32 = 32;
pub const N_DWELL: u32 = 4;

/// Which registered fields are the density/momentum components, for the
/// collapse ΔKE accounting (§3.4: conserving ρ𝐮 and ρE through a merge
/// turns sub-ring kinetic energy into internal energy — logged, not hidden).
#[derive(Debug, Clone, Copy)]
pub struct MomentumFields {
    pub rho: FieldId,
    pub mom_r: FieldId,
    pub mom_theta: FieldId,
    pub mom_z: FieldId,
}

/// One collapse event's ledger entry (COUP-2 diagnostics + the declared
/// truncation bound, §5).
#[derive(Debug, Clone, PartialEq)]
pub struct CollapseLog {
    pub brick: usize,
    pub n_theta_from: u32,
    pub n_theta_to: u32,
    /// Resolved kinetic energy [J] thermalized by the merge (≥ 0 by
    /// convexity); 0.0 when no momentum fields were declared.
    pub thermalized_ke: f64,
    /// True when the collapse was the recorded axisymmetry assertion, not
    /// an adaptive decision.
    pub by_assertion: bool,
}

impl Grid {
    /// One ladder step of conservative θ-coarsening on brick `bi`
    /// (N_θ → N_θ/2): merged θ-pairs get the volume-weighted mean of every
    /// field (cells in a merged pair share a ring ⇒ equal volumes ⇒ plain
    /// pair mean), preserving every ring integral to round-off. Adaptive
    /// collapse never drops below `max(N_θ^guard, N_θ^geom)` (§3.4).
    pub fn coarsen_theta(
        &mut self,
        bi: usize,
        momentum: Option<MomentumFields>,
    ) -> Result<CollapseLog, GridError> {
        let floor = self.bricks[bi].n_theta_geom_floor.max(N_THETA_GUARD);
        let from = self.bricks[bi].n_theta;
        let to = from / 2;
        if to < floor {
            return Err(GridError::BadThetaResolution {
                requested: to,
                reason: format!(
                    "adaptive collapse floor is max(N_θ^guard = {N_THETA_GUARD}, N_θ^geom = {}) \
                     — full N_θ = 1 is only a recorded config assertion (§3.4, S4)",
                    self.bricks[bi].n_theta_geom_floor
                ),
            });
        }
        let dke = self.merge_theta_pairs(bi, momentum)?;
        Ok(CollapseLog {
            brick: bi,
            n_theta_from: from,
            n_theta_to: to,
            thermalized_ke: dke,
            by_assertion: false,
        })
    }

    /// The §3.4 recorded axisymmetry assertion: collapse brick `bi` all the
    /// way to N_θ = 1. Permitted only when the spec carries the assertion —
    /// this is a config act, pedigree-visible, never an adaptive decision.
    pub fn assert_axisymmetric(
        &mut self,
        bi: usize,
        momentum: Option<MomentumFields>,
    ) -> Result<CollapseLog, GridError> {
        if !self.spec.axisymmetry_assertion {
            return Err(GridError::BadThetaResolution {
                requested: 1,
                reason: "spec carries no axisymmetry assertion (§3.4: N_θ = 1 must be a \
                         recorded, pedigree-visible config assertion)"
                    .into(),
            });
        }
        let from = self.bricks[bi].n_theta;
        let mut dke = 0.0f64;
        while self.bricks[bi].n_theta > 1 {
            dke += self.merge_theta_pairs(bi, momentum)?;
        }
        Ok(CollapseLog {
            brick: bi,
            n_theta_from: from,
            n_theta_to: 1,
            thermalized_ke: dke,
            by_assertion: true,
        })
    }

    /// Merge θ-pairs (2j, 2j+1) → j for every field; returns the resolved
    /// kinetic energy thermalized by this halving.
    fn merge_theta_pairs(
        &mut self,
        bi: usize,
        momentum: Option<MomentumFields>,
    ) -> Result<f64, GridError> {
        let from = self.bricks[bi].n_theta;
        let to = from / 2;
        debug_assert!(to >= 1);

        // ΔKE per merged pair, before overwriting (fixed iteration order).
        let mut dke = 0.0f64;
        if let Some(m) = momentum {
            let b = &self.bricks[bi];
            for local in 0..BRICK_CELLS {
                if b.mask & (1u64 << local) == 0 {
                    continue;
                }
                let i_r = b.br as usize * BRICK + local / BRICK;
                let vf = self.cell_volume(i_r, from);
                for j in 0..to {
                    let (ia, ib) = (b.cell_index(2 * j, local), b.cell_index(2 * j + 1, local));
                    let (ra, rb) = (b.field(m.rho)[ia], b.field(m.rho)[ib]);
                    if ra <= 0.0 || rb <= 0.0 {
                        return Err(GridError::NonPositiveDensity { brick: bi });
                    }
                    let ke = |i: usize| {
                        let (mr, mt, mz) = (
                            b.field(m.mom_r)[i],
                            b.field(m.mom_theta)[i],
                            b.field(m.mom_z)[i],
                        );
                        (mr * mr + mt * mt + mz * mz) / (2.0 * b.field(m.rho)[i])
                    };
                    let fine = vf * (ke(ia) + ke(ib));
                    let (mr, mt, mz) = (
                        0.5 * (b.field(m.mom_r)[ia] + b.field(m.mom_r)[ib]),
                        0.5 * (b.field(m.mom_theta)[ia] + b.field(m.mom_theta)[ib]),
                        0.5 * (b.field(m.mom_z)[ia] + b.field(m.mom_z)[ib]),
                    );
                    let rbar = 0.5 * (ra + rb);
                    let coarse = 2.0 * vf * (mr * mr + mt * mt + mz * mz) / (2.0 * rbar);
                    dke += fine - coarse;
                }
            }
        }

        // Conservative projection: pair mean, every field.
        let b = &mut self.bricks[bi];
        for data in &mut b.data {
            let mut out = vec![0.0f64; to as usize * BRICK_CELLS];
            for j in 0..to as usize {
                for local in 0..BRICK_CELLS {
                    let a = data[(2 * j) * BRICK_CELLS + local];
                    let c = data[(2 * j + 1) * BRICK_CELLS + local];
                    out[j * BRICK_CELLS + local] = 0.5 * (a + c);
                }
            }
            *data = out;
        }
        b.n_theta = to;
        Ok(dke)
    }

    /// One ladder step of conservative θ-refinement on brick `bi`
    /// (N_θ → 2·N_θ): limited piecewise-linear prolongation in θ (minmod —
    /// introduces no new extrema), children `q_j ∓ σ_j/4`, whose mean is
    /// exactly `q_j` ⇒ ring integrals preserved to round-off (§3.4).
    pub fn refine_theta(&mut self, bi: usize) -> Result<(), GridError> {
        let from = self.bricks[bi].n_theta;
        if from == 1 {
            return Err(GridError::BadThetaResolution {
                requested: 2,
                reason: "an asserted-axisymmetric region stays at N_θ = 1 (the assertion is \
                         recorded config intent; the indicator is structurally zero there)"
                    .into(),
            });
        }
        let to = from * 2;
        if to > self.spec.n_theta_max {
            return Err(GridError::BadThetaResolution {
                requested: to,
                reason: format!(
                    "above the config-declared finest N_θ^max = {}",
                    self.spec.n_theta_max
                ),
            });
        }
        let b = &mut self.bricks[bi];
        let n = from as usize;
        for data in &mut b.data {
            let mut out = vec![0.0f64; to as usize * BRICK_CELLS];
            for j in 0..n {
                let (jm, jp) = ((j + n - 1) % n, (j + 1) % n);
                for local in 0..BRICK_CELLS {
                    let q = data[j * BRICK_CELLS + local];
                    let (qm, qp) = (
                        data[jm * BRICK_CELLS + local],
                        data[jp * BRICK_CELLS + local],
                    );
                    let sigma = minmod(qp - q, q - qm);
                    out[(2 * j) * BRICK_CELLS + local] = q - 0.25 * sigma;
                    out[(2 * j + 1) * BRICK_CELLS + local] = q + 0.25 * sigma;
                }
            }
            *data = out;
        }
        b.n_theta = to;
        Ok(())
    }

    /// §3.4 (S7) symmetry indicator for brick `bi`:
    /// `A_q = Σ V·(q − ⟨q⟩_ring)² / (Σ V·⟨q⟩_ring² + V_b·floor_q²)`, maximum
    /// over the given conserved fields (each with its absolute floor —
    /// config-documented default 10⁻⁶ × the field's global reference
    /// magnitude). Fixed-order reductions; a pure function of the data.
    pub fn symmetry_indicator(&self, bi: usize, fields: &[(FieldId, f64)]) -> f64 {
        let b = &self.bricks[bi];
        let nt = b.n_theta as usize;
        let mut a_max = 0.0f64;
        for &(f, floor) in fields {
            let data = b.field(f);
            let mut num = 0.0f64;
            let mut den = 0.0f64;
            let mut v_brick = 0.0f64;
            for local in 0..BRICK_CELLS {
                if b.mask & (1u64 << local) == 0 {
                    continue;
                }
                let i_r = b.br as usize * BRICK + local / BRICK;
                let v = self.cell_volume(i_r, b.n_theta);
                let mut mean = 0.0f64;
                for j in 0..nt {
                    mean += data[j * BRICK_CELLS + local];
                }
                mean /= nt as f64;
                for j in 0..nt {
                    let d = data[j * BRICK_CELLS + local] - mean;
                    num += v * d * d;
                }
                den += v * (nt as f64) * mean * mean;
                v_brick += v * nt as f64;
            }
            let a = num / (den + v_brick * floor * floor);
            a_max = a_max.max(a);
        }
        a_max
    }
}

#[inline]
fn minmod(a: f64, b: f64) -> f64 {
    if a * b <= 0.0 {
        0.0
    } else if a.abs() < b.abs() {
        a
    } else {
        b
    }
}

/// §3.4 collapse/expand decision state for one brick: hysteresis
/// (τ_collapse ≪ τ_expand) plus a dwell of `N_DWELL` consecutive
/// over/under-threshold evaluations before any change. Driven by COUP-3
/// every `N_SYM_CADENCE` steps.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ThetaController {
    under: u32,
    over: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThetaAction {
    Hold,
    Collapse,
    Expand,
}

impl ThetaController {
    pub fn observe(&mut self, indicator: f64) -> ThetaAction {
        if indicator < TAU_COLLAPSE {
            self.under += 1;
            self.over = 0;
            if self.under >= N_DWELL {
                self.under = 0;
                return ThetaAction::Collapse;
            }
        } else if indicator > TAU_EXPAND {
            self.over += 1;
            self.under = 0;
            if self.over >= N_DWELL {
                self.over = 0;
                return ThetaAction::Expand;
            }
        } else {
            self.under = 0;
            self.over = 0;
        }
        ThetaAction::Hold
    }
}
