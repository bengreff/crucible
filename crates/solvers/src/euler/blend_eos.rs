//! SOLV-4 §3.6 — the **burn-progress blended EOS occupant**. The chemical
//! regime's thermochemistry is the `b`-blend (`b` = the doc's burn-progress
//! `c`, the code's `I_RB` slot) of two OFFL-3 branches, both bound as
//! [`TableEos`] on the same `(p, h, Z)` schema: **unburnt** (the S5 frozen-
//! reactant surface) and **burnt** (the shifting-equilibrium surface).
//!
//! **The energy-conserving flamelet closure (SOLV-4 §3.6).** Both sub-states
//! share the cell pressure `p` **and** the cell's actual static specific
//! enthalpy `h = e + p/ρ` — the adiabatic constant-`p` identity
//! `h_b = h_u = h` on S5's shared CEA formation reference (burning conserves
//! enthalpy), so the two branches are interrogated at the **same** `(p, h, Z)`
//! and differ only in composition. The cell's one mean density is partitioned
//! by **mass-weighted specific volume**
//! `1/ρ = (1−b)/ρ_u(p, h, Z) + b/ρ_b(p, h + h_offset, Z)`, and the pressure is
//! the root of that identity — the [`TableEos`] equilibrium projection
//! generalized to interrogate both branches (one deterministic Illinois
//! root find). It recovers the limits exactly: `b = 1` is the pure burnt
//! projection **bit-for-bit** (shifting mode), `b = 0` the pure unburnt one.
//!
//! `T_u = T_unburnt(p, h, Z)` (the unburnt branch's temperature at the cell's
//! actual `h`) is the reaction coordinate the SOLV-4 §3.6 rate law keys on:
//! cold at an energy-conserving flame front (`h` ≈ uniform), and **hot where
//! an igniter deposit or compression raised `h`** — which is exactly why the
//! igniter fires `τ_ign` (`T_u` sees the deposited enthalpy).

use super::{
    Cons, EosLaw, FlowError, I_EI, I_EN, I_G1, I_MR, I_MT, I_MZ, I_RB, I_RC, I_RHO, Prim, TableEos,
};

/// **Pure-BURNT threshold** on the `[0,1]` progress variable `b`: at
/// `b ≥ 1 − EPS_B_PURE_BURNT` the cell is read as pure burnt (the unburnt
/// branch is skipped). **Derived from the one owner,
/// [`super::BURN_COMPLETE`] (SOLV-4 0.4.3):** the combustion source zeroes
/// its rate (and its closure queries: post-flame gas is not a reactant, the
/// §3.6 domain guard) at `b ≥ 1 − BURN_COMPLETE`; the S6 explicit tier
/// **asymptoted toward that fixed point from below**, and the S7 class-`R`
/// implicit solve **parks exactly AT it** (SOLV-4 0.4.4) (S6 measured:
/// cells park ~2e-7 under it). The blend's pure-burnt region must therefore
/// *strictly contain* that pinning attractor — an equal threshold is a
/// knife edge on whose wrong side every burnt cell keeps interrogating the
/// unburnt branch at ~10⁻³ weight forever, and hot burnt gas (whose `h`
/// legitimately exceeds the unburnt surface's ceiling) refuses on a branch
/// describing 0.1 % of its mass. Factor 2 gives the attractor ~100× margin
/// over the largest per-step creep (`Δb ~ (1−b)·Δt/τ`). The crossing step
/// this skip declares is `EPS·(v_u/v_b)` ≈ **2.6e-4** relative at flame
/// states (the dropped unburnt branch is the *dense* one, `v_u ≪ v_b`) —
/// far inside the density column's own interp bound. **Not** a physics
/// threshold; the physical `b ∈ [0,1]` bound is owned by the combustion
/// source's loud positivity guard (SOLV-4 §3.6), never clamped here.
pub const EPS_B_PURE_BURNT: f64 = 2.0 * super::BURN_COMPLETE;

/// **Pure-UNBURNT threshold**: at `b ≤ EPS_B_PURE_UNBURNT` the cell is read
/// as pure unburnt (the burnt branch is skipped — the envelope-robustness
/// motive: a cold `b = 0`⁺ cell projects without a spurious burnt-branch
/// query). **Deliberately asymmetric to [`EPS_B_PURE_BURNT`] (review
/// finding, S6 close):** the reaction pins an attractor only at the BURNT
/// end — on the cold side the metastable fringe (`b < a`) decays smoothly
/// through every value to 0, parking nowhere — and the cold-side crossing
/// step scales as `EPS·(v_b/v_u)` ≈ **7.8·EPS** at flame states (here the
/// dropped burnt branch is the *light* one), so a 2e-3 threshold would step
/// the specific volume ~1.5 % — *outside* the density column's own ~0.8 %
/// bound. At 1e-9 the step is ~8e-9: a true regularization.
pub const EPS_B_PURE_UNBURNT: f64 = 1.0e-9;

/// Inward relative shrink on the closed-form h-derived pressure bracket so
/// endpoint queries cannot fall an ulp outside either surface's envelope.
/// Wider than `TableEos`'s 1e-12 (S7): the blend's partition arithmetic
/// (the mass-weighted floor bound and the ÷b balance form) accumulates a
/// few ulps between the bracket bound and the branch query coordinate, so
/// an endpoint can land a rounding error outside the envelope; 1e-9
/// relative on h is ~10⁻² J/kg — far below every declared bound.
const H_BRACKET_MARGIN: f64 = 1.0e-9;

/// Fixed log-scan resolution of the non-monotone fallback (mirrors
/// [`TableEos`]'s `N_P_SCAN`); sized to match the tables' own p-axis density.
const N_P_SCAN: usize = 64;

/// Outcome of the fixed mid-`b` log-scan ([`scan_first_crossing`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum PScan {
    /// `g` hit exactly zero at a scan node before any sign change — the node
    /// IS the accepted root.
    ExactRoot(f64),
    /// The FIRST sign change over the sweep: the crossing's bracket `(a, b)`
    /// with its endpoint residuals `(ga, gb)`.
    Bracket { a: f64, b: f64, ga: f64, gb: f64 },
    /// No sign change anywhere: the `|g|`-best sample over the whole sweep
    /// (the tangency-acceptance seed).
    NoCrossing { p_best: f64 },
}

/// The mid-transition scan (SOLV-4 §3.6 v0.4.9): the fixed `N_P_SCAN`
/// log-scan of `g` over `[lo, hi]` from the cold end, returning the FIRST
/// sign crossing's bracket — **multi-root resolves by CONTINUITY** (Ben
/// ruling 2026-08-31, PLAN §8 v1.16 (3)): a folded `g` with a further
/// admissible root past the first crossing is a degenerate-lookup corner of
/// the tabulated surface, resolved by keeping the branch reached from the
/// bracket's cold end — the same rule `TableEos`'s warm path and the resident
/// GPU projection carry. The v0.4.8 `>1 crossing ⇒ typed refusal` (and the
/// sweep completion that was only its counter) is retired: the loop breaks at
/// the first crossing, the pre-v0.4.8 structure.
///
/// **Bit-identity contract:** the scan grid (`lo·(hi/lo)^(k/N)`), the
/// `g`-evaluation order, and the `best`-sample tracking are byte-for-byte the
/// pre-v0.4.8 loop's, so every single-root bracket and every no-crossing
/// tangency seed are the identical floats — every previously-accepted root is
/// bit-identical. An `Err` from `g` propagates from any visited node.
///
/// `ga0` must be the caller's already-computed `g(lo)` (nonzero — the caller
/// returns on exact-zero endpoints before scanning).
pub(crate) fn scan_first_crossing(
    lo: f64,
    hi: f64,
    ga0: f64,
    g: &impl Fn(f64) -> Result<f64, &'static str>,
) -> Result<PScan, &'static str> {
    let ratio = hi / lo;
    let mut prev_p = lo;
    let mut prev_g = ga0;
    let mut best = (lo, ga0.abs());
    for k in 1..=N_P_SCAN {
        let pk = lo * ratio.powf(k as f64 / N_P_SCAN as f64);
        let gk = g(pk)?;
        if gk == 0.0 {
            return Ok(PScan::ExactRoot(pk));
        }
        if gk.abs() < best.1 {
            best = (pk, gk.abs());
        }
        if prev_g * gk < 0.0 {
            return Ok(PScan::Bracket {
                a: prev_p,
                b: pk,
                ga: prev_g,
                gb: gk,
            });
        }
        prev_p = pk;
        prev_g = gk;
    }
    Ok(PScan::NoCrossing { p_best: best.0 })
}

/// The SOLV-4 §3.6 blended occupant: two `(p, h, Z)` branches + the S18
/// `h_offset` knockdown (applied to the **burnt** interrogation only — the
/// unburnt reactants carry no combustion-completeness deficit).
#[derive(Debug, Clone)]
pub struct BurnBlendEos<'t> {
    unburnt: TableEos<'t>,
    burnt: TableEos<'t>,
    /// S18 η_c\* knockdown on the burnt branch (J/kg); 0 = full equilibrium.
    pub h_offset: f64,
}

impl<'t> BurnBlendEos<'t> {
    /// Bind the blend to its two already-bound branches. The unburnt branch
    /// must be the S5 frozen-reactant surface, the burnt branch the
    /// equilibrium surface — both `(p, h, Z)` (checked by `TableEos::bind`).
    pub fn new(unburnt: TableEos<'t>, burnt: TableEos<'t>) -> Self {
        Self {
            unburnt,
            burnt,
            h_offset: 0.0,
        }
    }

    /// The burnt branch (S18 knockdown carrier + shifting-mode fast paths).
    pub fn burnt(&self) -> &TableEos<'t> {
        &self.burnt
    }

    /// Conserved state from `(p, h, Z, b)` and a velocity — the forward blend
    /// (density = the mass-weighted specific volume at the common `(p, h, Z)`),
    /// for constructor initial conditions. `b = 0` is the pure unburnt state,
    /// `b = 1` the pure burnt — both at the *same* `h` (the adiabatic-flamelet
    /// identity), so a flame front is initialized by varying only `b`.
    pub fn cons_from_phzb(
        &self,
        p: f64,
        h: f64,
        z: f64,
        b: f64,
        vel: [f64; 3],
    ) -> Result<Cons, &'static str> {
        let b_c = b.clamp(0.0, 1.0);
        let rho = 1.0 / self.blend_inv_rho(p, h, z, b_c)?;
        let e = h - p / rho;
        let ke = 0.5 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
        Ok([
            rho,
            rho * vel[0],
            rho * vel[1],
            rho * vel[2],
            rho * (e + ke),
            rho * z,
            rho * b,
        ])
    }

    /// The sub-state enthalpy partition `(h_u, h_b)` at a trial cell
    /// enthalpy (SOLV-4 §3.6, extended S7). On the shared-`h` flamelet
    /// domain (`h ≥` the unburnt branch's h-envelope floor) both branches
    /// read the cell's own `h` — the adiabatic constant-`p` identity.
    /// **Below the floor** that identity's premise (adiabatic burning) has
    /// been broken by real heat loss — the quench trajectory: wall cooling
    /// drives a mid-transition cell's mixture enthalpy under any
    /// representable reactant state. The **declared cold-side extension**:
    /// the reactant sub-state pins AT the floor (the coldest representable
    /// reactant gas — physically its condensation edge) and the products
    /// carry the balance, `h_b = (h − (1−b)·h_floor)/b` — continuous at the
    /// floor (`h_b = h` there), mass-consistent by construction, and
    /// self-limiting: as cooling continues `h_b` walks down the burnt
    /// branch until ITS envelope refuses — the blend's true cold edge —
    /// while a low-`b` sub-floor cell (mostly unrepresentably-cold
    /// reactants) refuses through the same gate immediately (the huge
    /// `1/b` deficit lands far off the burnt surface). No clamp: every
    /// state either projects on declared physics or refuses loudly.
    fn partition_h(&self, h: f64, b: f64) -> (f64, f64) {
        let (hu_floor, hu_ceil) = self.unburnt.envelopes()[1];
        if (h >= hu_floor && h <= hu_ceil) || b <= EPS_B_PURE_UNBURNT {
            return (h, h);
        }
        // TWO-SIDED (S7 ◆C2 shake-out): the reactant sub-state pins at the
        // NEAR edge of its envelope — the floor under heat loss (quench),
        // the CEILING under spark/blast superheat (a mid-transition kernel
        // cell whose h passed the metastable-reactant ceiling is sub-µs
        // from burnt; pinning keeps the closures live at the ceiling values
        // so the class-R reaction finishes it — the model self-heals
        // through the reaction rather than refusing).
        let hu_pin = h.clamp(hu_floor, hu_ceil);
        // The WINDOW form (S8 review wave, SOLV-4 0.4.7 — supersedes the S7
        // balance↔trace two-form fallback): the products carry the balance,
        // taken into their own h-window by a declared clamp. Continuous in
        // BOTH h and b — the S7 switch was measured discontinuous inside
        // the projection's own scan bracket (a pseudo-root site), and its
        // B_PARTITION_MIN crossing stepped the mixture volume by
        // (1−b)·δ·∂v_b/∂h: the 1/b amplification cancels the b weight
        // exactly, so the step was NOT trace-bounded (the ◆C2 purge cells
        // sat at the measured knife edge b = 0.010). Where the clamp
        // engages, the residual mixture deficit is declared UNATTRIBUTED —
        // it is exactly the state's distance beyond every representable
        // mixture at this b — and the projection's root-residual acceptance
        // refuses when it exceeds the declared interpolation error (the
        // blend's true cold edge, enforced at the acceptance).
        let (hb_lo, hb_hi) = {
            let e = self.burnt.envelopes()[1];
            (e.0 - self.h_offset, e.1 - self.h_offset)
        };
        let balance = (h - (1.0 - b) * hu_pin) / b;
        (hu_pin, balance.clamp(hb_lo, hb_hi))
    }

    /// `1/ρ` of the blend at a trial `(p, h, Z)` and burn fraction `b`: the
    /// mass-weighted specific volume on the partitioned sub-state
    /// enthalpies ([`Self::partition_h`]), skipping a branch of ~zero
    /// weight (so the pure limits never query the other surface's
    /// envelope).
    fn blend_inv_rho(&self, p: f64, h: f64, z: f64, b: f64) -> Result<f64, &'static str> {
        let (h_u, h_b) = self.partition_h(h, b);
        let mut v = 0.0;
        if b < 1.0 - EPS_B_PURE_BURNT {
            let ru = self
                .unburnt
                .rho_at(p, h_u, z)
                .map_err(|_| "unburnt branch query failed inside the blend projection bracket")?;
            v += (1.0 - b) / ru;
        }
        if b > EPS_B_PURE_UNBURNT {
            let rb = self
                .burnt
                .rho_at(p, h_b + self.h_offset, z)
                .map_err(|_| "burnt branch query failed inside the blend projection bracket")?;
            v += b / rb;
        }
        Ok(v)
    }

    /// Project the pressure: root of `blend_inv_rho(p) − 1/ρ = 0` over the
    /// admissible bracket (the intersection of the *required* branches'
    /// envelopes). Deterministic Illinois; a mild non-monotone corner falls
    /// back to a fixed log-scan for a sign change, then (S7) a fixed-count
    /// golden-section **tangency acceptance** in the branches' own
    /// mass-weighted rule-space bound — the RL10-plume feature.
    ///
    /// **Pure limits delegate to the branch's own `TableEos` projection
    /// (S7):** `b = 1` runs the burnt surface's projection (at the knocked
    /// coordinate `e + h_offset`) and `b = 0` the unburnt's — so the pure
    /// corners carry exactly the single-surface semantics, near-vacuum
    /// tangency acceptance included, with one owner for that code path.
    fn project_pressure(
        &self,
        rho: f64,
        e: f64,
        z: f64,
        b: f64,
        hint: Option<f64>,
    ) -> Result<f64, &'static str> {
        let inv = 1.0 / rho;
        let need_u = b < 1.0 - EPS_B_PURE_BURNT;
        let need_b = b > EPS_B_PURE_UNBURNT;
        if !need_b {
            return match hint {
                Some(ph) if ph.is_finite() && ph > 0.0 => {
                    self.unburnt.project_pressure_hinted(rho, e, z, ph)
                }
                _ => self.unburnt.project_pressure(rho, e, z),
            };
        }
        if !need_u {
            let e_q = e + self.h_offset;
            return match hint {
                Some(ph) if ph.is_finite() && ph > 0.0 => {
                    self.burnt.project_pressure_hinted(rho, e_q, z, ph)
                }
                _ => self.burnt.project_pressure(rho, e_q, z),
            };
        }
        let [(pu_lo, pu_hi), (hu_lo, hu_hi), _] = self.unburnt.envelopes();
        let [(pb_lo, pb_hi), (hb_lo, hb_hi), _] = self.burnt.envelopes();

        // h-window (of h = e + p/ρ) and p-window over the two branches (the
        // mid-b path always needs both — the pure limits delegated above).
        // The burnt window is shifted by −h_offset (that branch is queried
        // at h_b + h_offset), and the LOWER bound carries the cold-side
        // partition extension ([`Self::partition_h`]): the cell h may run
        // down to the mass-weighted floor `(1−b)·hu_lo + b·(hb_lo − off)` —
        // reactant sub-state pinned at its floor, products at theirs —
        // which reduces to `max(hu_lo, hb_lo − off)`-style shared-h logic
        // continuously at `b → 1⁻` and to `hu_lo` at `b → 0⁺`.
        let (mut h_lo, mut h_hi) = (f64::NEG_INFINITY, f64::INFINITY);
        let (mut p_lo_env, mut p_hi_env) = (0.0_f64, f64::INFINITY);
        if need_u && need_b {
            // The mid-b admissible h-window is the PRODUCTS branch's own
            // (S7 ◆C2, the supersonic-purge finding): the partition's
            // pin/balance/trace machinery makes every cell-h queryable
            // wherever h_b lands inside the products window — the reactant
            // sub-state pins at its near envelope edge outside its own
            // band, the balance form attributes the difference where
            // representable, the trace form where not. The products window
            // strictly contains the reactant band, so it IS the window; a
            // mass-weighted bracket floor excluded real roots (measured: a
            // 1060 m/s nozzle purge cell at c ≈ 0.01 with e_int below the
            // mass-weighted floor and its root at ~1e4–5e4 Pa).
            h_lo = hb_lo - self.h_offset;
            h_hi = hb_hi - self.h_offset;
            p_lo_env = pu_lo.max(pb_lo);
            p_hi_env = pu_hi.min(pb_hi);
        } else if need_u {
            h_lo = h_lo.max(hu_lo);
            h_hi = h_hi.min(hu_hi);
            p_lo_env = p_lo_env.max(pu_lo);
            p_hi_env = p_hi_env.min(pu_hi);
        } else if need_b {
            h_lo = h_lo.max(hb_lo - self.h_offset);
            h_hi = h_hi.min(hb_hi - self.h_offset);
            p_lo_env = p_lo_env.max(pb_lo);
            p_hi_env = p_hi_env.min(pb_hi);
        }
        let p_from_h_lo = rho * (h_lo - e);
        let p_from_h_hi = rho * (h_hi - e);
        let mut lo = p_lo_env
            .max(p_from_h_lo * (1.0 + H_BRACKET_MARGIN))
            .max(0.0);
        let hi = p_hi_env.min(p_from_h_hi * (1.0 - H_BRACKET_MARGIN));
        if !(lo.is_finite() && hi.is_finite()) || lo >= hi || hi <= 0.0 {
            return Err("no admissible pressure bracket: blended state outside the table envelope");
        }
        lo = lo.max(0.0);

        let g = |p: f64| -> Result<f64, &'static str> {
            Ok(self.blend_inv_rho(p, e + p * inv, z, b)? - inv)
        };
        let ga0 = g(lo)?;
        let gb0 = g(hi)?;
        if ga0 == 0.0 {
            return Ok(lo);
        }
        if gb0 == 0.0 {
            return Ok(hi);
        }
        // The mid-b bracket is the PRODUCTS branch's whole h-window (above)
        // — up to five decades of p with wildly asymmetric endpoint
        // magnitudes (a near-vacuum endpoint's specific volume dwarfs the
        // dense end's), where regula-falsi creeps from the flat end and can
        // exhaust its fixed budget (measured on the ◆C2 purge state). The
        // fixed log-scan below therefore ALWAYS runs first: it brackets the
        // first sign change to a single scan cell (~a fifth of a decade),
        // where Illinois converges in a handful of iterations; no sign
        // change anywhere falls through to the tangency acceptance. A
        // further root past the first crossing resolves by CONTINUITY
        // ([`scan_first_crossing`], SOLV-4 §3.6 v0.4.9) — the first
        // crossing from the cold end is the accepted branch, never a halt.
        let ratio = hi / lo;
        let p_best = match scan_first_crossing(lo, hi, ga0, &g)? {
            PScan::ExactRoot(pk) => return Ok(pk),
            PScan::Bracket {
                a: pa,
                b: pb,
                ga,
                gb,
            } => {
                // Root-residual acceptance (S8 review wave, SOLV-4 0.4.7):
                // Illinois/bisection accepts on BRACKET WIDTH, so a fold or
                // residual partition corner could hand back a pseudo-root
                // with a finite density mismatch. Every accepted root must
                // meet the blended surface within the branches' own
                // volume-weighted rule-space bound — the same acceptance
                // the tangency path applies. (An artifact without stamped
                // bounds keeps the bracket-width acceptance — the pre-S8
                // semantics; every production surface stamps them.)
                let p_root = TableEos::illinois_root(pa, pb, ga, gb, &g)?;
                let g_root = g(p_root)?;
                return match self.root_within_bound(p_root, g_root, rho, e, z, b)? {
                    None | Some(true) => Ok(p_root),
                    Some(false) => Err(
                        "blend root fails the volume-weighted rule-space acceptance — \
                         pseudo-root or a state beyond the declared interpolation error \
                         (the blend's cold edge refuses here)",
                    ),
                };
            }
            PScan::NoCrossing { p_best } => p_best,
        };
        // Near-vacuum tangency acceptance (S7, mirroring `TableEos`): the
        // constraint line h = e + p/ρ runs almost parallel to the blended
        // ρ-contour in the deep plume, so the crossing can be a tangency the
        // scan never sign-changes on. Golden-section minimize |g| around the
        // best sample, then accept iff the blended surface is met within the
        // branches' own MASS-WEIGHTED rule-space density bounds — the blend's
        // declared interpolation error, data not a tunable (META-1 P6). Both
        // production branches stamp log bounds; an artifact without one keeps
        // the refusal (pre-0.3.2 behavior, preserved bit-for-bit).
        let (Some(lb_u), Some(lb_b)) = (self.unburnt.rho_bound_log(), self.burnt.rho_bound_log())
        else {
            return Err(
                "blended state off both surfaces over the admissible bracket — no sign \
                 change and a branch carries no rule-space bound for tangency acceptance",
            );
        };
        let phi = 0.618_033_988_749_894_9_f64;
        let (mut x0, mut x3) = (
            (p_best / ratio.powf(1.0 / N_P_SCAN as f64)).max(lo),
            (p_best * ratio.powf(1.0 / N_P_SCAN as f64)).min(hi),
        );
        let mut x1 = x3 - phi * (x3 - x0);
        let mut x2 = x0 + phi * (x3 - x0);
        let mut g1 = g(x1)?; // signed: v_blend − 1/ρ
        let mut g2 = g(x2)?;
        for _ in 0..super::table_eos::N_P_ITER_MAX {
            if g1.abs() < g2.abs() {
                x3 = x2;
                x2 = x1;
                g2 = g1;
                x1 = x3 - phi * (x3 - x0);
                g1 = g(x1)?;
            } else {
                x0 = x1;
                x1 = x2;
                g1 = g2;
                x2 = x0 + phi * (x3 - x0);
                g2 = g(x2)?;
            }
        }
        let (p_best, g_best) = if g1.abs() < g2.abs() {
            (x1, g1)
        } else {
            (x2, g2)
        };
        // g is a specific-volume mismatch; ratio = v_blend/v_cell = 1 + g·ρ,
        // so |ln ratio| is the same rule-space measure the branches' density
        // bounds are stamped in. The bound composes VOLUME-weighted (review
        // finding): v = (1−b)v_u + b·v_b is volume-additive, so the
        // first-order sensitivity of ln v to each branch's own ln-ρ bound
        // carries that branch's volume share — at flame states the light
        // burnt branch dominates v long before b → 1.
        let vratio = 1.0 + g_best * rho;
        let h_best = e + p_best * inv;
        let (h_u, h_b) = self.partition_h(h_best, b);
        let vu = (1.0 - b)
            / self
                .unburnt
                .rho_at(p_best, h_u, z)
                .map_err(|_| "unburnt branch query failed at the tangency acceptance")?;
        let vb = b / self
            .burnt
            .rho_at(p_best, h_b + self.h_offset, z)
            .map_err(|_| "burnt branch query failed at the tangency acceptance")?;
        let bound = (vu * lb_u + vb * lb_b) / (vu + vb);
        if vratio > 0.0 && vratio.ln().abs() <= bound {
            return Ok(p_best);
        }
        Err(
            "blended state off both surfaces beyond their volume-weighted declared \
             interpolation-error bound — no admissible projection",
        )
    }

    /// The S8 root-residual test (SOLV-4 0.4.7): is the accepted root's
    /// specific-volume mismatch within the branches' volume-weighted
    /// rule-space density bounds? `Ok(None)` when a branch carries no
    /// stamped bound (non-production artifact — the caller keeps the
    /// bracket-width acceptance, the pre-S8 semantics).
    fn root_within_bound(
        &self,
        p: f64,
        g_val: f64,
        rho: f64,
        e: f64,
        z: f64,
        b: f64,
    ) -> Result<Option<bool>, &'static str> {
        if g_val == 0.0 {
            return Ok(Some(true));
        }
        let (Some(lb_u), Some(lb_b)) = (self.unburnt.rho_bound_log(), self.burnt.rho_bound_log())
        else {
            return Ok(None);
        };
        let vratio = 1.0 + g_val * rho;
        if vratio <= 0.0 {
            return Ok(Some(false));
        }
        let h_p = e + p / rho;
        let (h_u, h_b) = self.partition_h(h_p, b);
        let vu = (1.0 - b)
            / self
                .unburnt
                .rho_at(p, h_u, z)
                .map_err(|_| "unburnt branch query failed at the root acceptance")?;
        let vb = b / self
            .burnt
            .rho_at(p, h_b + self.h_offset, z)
            .map_err(|_| "burnt branch query failed at the root acceptance")?;
        let bound = (vu * lb_u + vb * lb_b) / (vu + vb);
        Ok(Some(vratio.ln().abs() <= bound))
    }

    /// `true` when the cell's enthalpy sits below the unburnt branch's own
    /// h-envelope floor — **colder than any representable reactant** (the
    /// gas-phase metastable model's condensation edge, ~100 K class). The
    /// second face of the SOLV-4 §3.6 cold-side non-reactive floor (S7):
    /// such a cell cannot be interrogated for `T_u` at all, and the honest
    /// physical statement is that near-condensing reactants do not burn —
    /// the rate law reads zero there by declaration, never a refusal.
    /// (Reached in practice by strong wall cooling of nearly-burnt cells,
    /// whose `b` parks just under the reaction's `1 − BURN_COMPLETE` fixed
    /// point and therefore keeps its closure queries live — the quench_box
    /// trajectory.)
    pub fn below_unburnt_floor(&self, w: &Prim) -> bool {
        let h = w[I_EI] + w[4] / w[I_RHO];
        h < self.unburnt.envelopes()[1].0
    }

    /// The unburnt-branch temperature `T_u(p, h, Z)` at a projected primitive
    /// — the SOLV-4 §3.6 rate-law coordinate (`S_L`, `τ_ign`). `h` is the
    /// cell's actual static enthalpy `e + p/ρ`, so `T_u` rises with an igniter
    /// deposit.
    pub fn unburnt_temperature(&self, w: &Prim) -> Result<f64, &'static str> {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        let h = w[I_EI] + p / rho;
        // The rate-law coordinate reads the PINNED reactant sub-state (the
        // two-sided partition): a superheated mid-transition cell's T_u is
        // the ceiling temperature — its closures stay live (τ_ign sub-µs
        // there) so the reaction can finish it.
        let (h_u, _) = self.partition_h(h, w[I_RB].clamp(0.0, 1.0));
        self.unburnt
            .temp_at(p, h_u, z)
            .map_err(|_| "unburnt temperature query failed (T_u off the ignition-surface envelope)")
    }

    /// The unburnt-branch density `ρ_u(p, h, Z)` — the `ρ_u` factor of the
    /// propagation source `ρ_u·S_T·|∇c|` (SOLV-4.4).
    pub fn unburnt_density(&self, w: &Prim) -> Result<f64, &'static str> {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        let h = w[I_EI] + p / rho;
        let (h_u, _) = self.partition_h(h, w[I_RB].clamp(0.0, 1.0));
        self.unburnt
            .rho_at(p, h_u, z)
            .map_err(|_| "unburnt density query failed inside the blend")
    }

    /// The **blended diagnostic temperature** `T = (1−b)·T_u + b·T_b` at a
    /// projected primitive — the SOLV-4 §3.6 declared mass-weighted model
    /// form, recovering each pure limit. The blend↔class-D seam (S7): this is
    /// the temperature the gas-diffusion operator and the wall law read on a
    /// blended run (the same closure-injection seam as `TableEos`), so a
    /// mid-transition cell conducts on its mixture temperature, not one
    /// branch's. Skips a ~zero-weight branch exactly as the projection does.
    pub fn temperature_w(&self, w: &Prim) -> Result<f64, &'static str> {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        let b = w[I_RB].clamp(0.0, 1.0);
        let h = w[I_EI] + p / rho;
        let (h_u, h_b) = self.partition_h(h, b);
        let mut t = 0.0;
        if b < 1.0 - EPS_B_PURE_BURNT {
            let tu = self
                .unburnt
                .temp_at(p, h_u, z)
                .map_err(|_| "unburnt temperature query failed at the projected state")?;
            t += (1.0 - b) * tu;
        }
        if b > EPS_B_PURE_UNBURNT {
            let tb = self
                .burnt
                .temp_at(p, h_b + self.h_offset, z)
                .map_err(|_| "burnt temperature query failed at the projected state")?;
            t += b * tb;
        }
        Ok(t)
    }

    /// The FND-7 §3.3 spine interrogation coordinate at a blended primitive
    /// (the `TableEos::interrogation_php` mirror, S7): `(p, h + b·h_offset,
    /// Z)`. The knockdown offsets only the **burnt** branch's coordinate, so
    /// the transport interrogation carries it mass-weighted — reducing to the
    /// `TableEos` coordinate at `b = 1` and the true enthalpy at `b = 0`. The
    /// spine's chemical-regime surface is generated on the EQUILIBRIUM
    /// composition; reading it for unburnt/mid-transition cells is a
    /// **declared model-form gap** (FND-7 0.5.3 — the unburnt-composition
    /// transport rides S5b's per-species work), inside the spine's own
    /// declared 10–20 % band for the states a startup march visits.
    pub fn interrogation_php(&self, w: &Prim) -> crate::transport::MediumState {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        let b = w[I_RB].clamp(0.0, 1.0);
        crate::transport::MediumState {
            p,
            h: w[I_EI] + b * self.h_offset + p / rho,
            z,
        }
    }

    /// The declared (p, h, Z) envelope box of the BLEND: the intersection of
    /// the two branches' envelopes (the burnt h-window shifted by
    /// −`h_offset`, since that branch is interrogated at `h + h_offset`) —
    /// the states a mid-transition cell can legally occupy. Assembly-time
    /// sizing/refusals (the `TableEos::envelopes` mirror, S7).
    pub fn envelopes(&self) -> [(f64, f64); 3] {
        let u = self.unburnt.envelopes();
        let b = self.burnt.envelopes();
        [
            (u[0].0.max(b[0].0), u[0].1.min(b[0].1)),
            (
                u[1].0.max(b[1].0 - self.h_offset),
                u[1].1.min(b[1].1 - self.h_offset),
            ),
            (u[2].0.max(b[2].0), u[2].1.min(b[2].1)),
        ]
    }

    /// S13c GPU marshaling (doc-hidden, additive): the UNBURNT branch's
    /// temperature + density columns (the combustion source's `T_u`, `ρ_u`
    /// direct interps) and that surface's own (p, h, Z) envelopes. Changes no
    /// production number.
    #[doc(hidden)]
    pub fn xcheck_unburnt_marshal(
        &self,
    ) -> (
        crucible_tables::ColumnMarshal,
        crucible_tables::ColumnMarshal,
        [(f64, f64); 3],
    ) {
        let m = self.unburnt.xcheck_marshal();
        (m.temperature, m.rho, self.unburnt.envelopes())
    }

    /// S13c GPU marshaling (doc-hidden, additive): both branches' DENSITY
    /// columns + their (p, h, Z) envelopes + `h_offset` — everything the device
    /// blend projection `gpu_blend_project` needs to reproduce
    /// `project_pressure`/`blend_inv_rho` bit-for-formula. Changes no production
    /// number.
    #[doc(hidden)]
    #[allow(clippy::type_complexity)]
    pub fn xcheck_blend_marshal(
        &self,
    ) -> (
        crucible_tables::ColumnMarshal,
        crucible_tables::ColumnMarshal,
        [(f64, f64); 3],
        [(f64, f64); 3],
        f64,
    ) {
        (
            self.unburnt.xcheck_marshal().rho,
            self.burnt.xcheck_marshal().rho,
            self.unburnt.envelopes(),
            self.burnt.envelopes(),
            self.h_offset,
        )
    }

    /// The UNION of the two branches' (p, h, Z) envelopes (the burnt
    /// h-window shifted by −`h_offset`) — the box this occupant may
    /// interrogate a co-keyed surface (the FND-7 spine) over. Distinct from
    /// [`Self::envelopes`] (the intersection — the mid-transition
    /// PROJECTION's admissible box): a pure-burnt cell legally roams the
    /// whole burnt envelope (S7 review finding).
    pub fn interrogation_envelopes(&self) -> [(f64, f64); 3] {
        let u = self.unburnt.envelopes();
        let b = self.burnt.envelopes();
        [
            (u[0].0.min(b[0].0), u[0].1.max(b[0].1)),
            (
                u[1].0.min(b[1].0 - self.h_offset),
                u[1].1.max(b[1].1 - self.h_offset),
            ),
            (u[2].0.min(b[2].0), u[2].1.max(b[2].1)),
        ]
    }

    fn prim_checked_impl(&self, u: &Cons, hint: Option<f64>) -> Result<Prim, &'static str> {
        let rho = u[I_RHO];
        if !rho.is_finite() || rho <= 0.0 {
            return Err("non-positive or non-finite density");
        }
        let inv = 1.0 / rho;
        let (ur, ut, uz) = (u[I_MR] * inv, u[I_MT] * inv, u[I_MZ] * inv);
        let e = u[I_EN] * inv - 0.5 * (ur * ur + ut * ut + uz * uz);
        if !e.is_finite() {
            return Err("non-finite internal energy");
        }
        let z = u[I_RC] * inv;
        if !z.is_finite() {
            return Err("non-finite elemental mixture fraction");
        }
        let b_raw = u[I_RB] * inv;
        if !b_raw.is_finite() {
            return Err("non-finite burn progress");
        }
        // The blend math needs b ∈ [0,1]; a tiny advective over/undershoot is
        // a numerical fact about a [0,1] fraction, so it is bounded HERE for
        // the projection only — the STORED b (and the audited ρb) stay the
        // true advected value, and the physical [0,1] violation guard belongs
        // to the combustion source (SOLV-4 §3.6), not this projection.
        let b = b_raw.clamp(0.0, 1.0);
        let p = self.project_pressure(rho, e, z, b, hint)?;
        let h = e + p * inv;
        let (h_u, h_b) = self.partition_h(h, b);
        // Mass-weighted blended sound speed a² = (1−b)a_u² + b·a_b² (declared
        // model-form; skips a ~zero-weight branch, matching the projection).
        let mut a2 = 0.0;
        if b < 1.0 - EPS_B_PURE_BURNT {
            let au = self
                .unburnt
                .sound_at(p, h_u, z)
                .map_err(|_| "unburnt sound-speed query failed at the projected state")?;
            a2 += (1.0 - b) * au * au;
        }
        if b > EPS_B_PURE_UNBURNT {
            let ab = self
                .burnt
                .sound_at(p, h_b + self.h_offset, z)
                .map_err(|_| "burnt sound-speed query failed at the projected state")?;
            a2 += b * ab * ab;
        }
        let g1 = rho * a2 / p;
        Ok([rho, ur, ut, uz, p, z, b_raw, e, g1])
    }
}

impl EosLaw for BurnBlendEos<'_> {
    fn prim_checked(&self, u: &Cons) -> Result<Prim, &'static str> {
        self.prim_checked_impl(u, None)
    }

    /// Warm start (S7): the hint reaches the PURE-LIMIT delegated
    /// projections (where `TableEos`'s root-uniqueness guard makes it a pure
    /// acceleration); the mid-transition generalized projection ignores it —
    /// front cells are a thin minority, so the pure limits are where the
    /// march's cost lives.
    fn prim_checked_hinted(&self, u: &Cons, hint: Option<f64>) -> Result<Prim, &'static str> {
        self.prim_checked_impl(u, hint)
    }

    fn prim_to_cons(&self, w: &Prim) -> Cons {
        let rho = w[I_RHO];
        let ke = 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]);
        [
            rho,
            rho * w[1],
            rho * w[2],
            rho * w[3],
            rho * (w[I_EI] + ke),
            rho * w[I_RC],
            rho * w[I_RB],
        ]
    }

    fn total_energy(&self, w: &Prim) -> f64 {
        let rho = w[I_RHO];
        rho * (w[I_EI] + 0.5 * (w[1] * w[1] + w[2] * w[2] + w[3] * w[3]))
    }

    fn sound_speed_w(&self, w: &Prim) -> f64 {
        (w[I_G1] * w[4] / w[I_RHO]).sqrt()
    }

    fn roe_sound_speed(
        &self,
        wl: &Prim,
        wr: &Prim,
        _h_roe: f64,
        _q2_roe: f64,
        sql: f64,
        sqr: f64,
        inv: f64,
    ) -> f64 {
        // The datum-free general-EOS Roe sound speed (as `TableEos`): average
        // the squared sound speeds `c² = Γ₁·p/ρ` from each side's own
        // projected (blended) state — always positive, exact for identical
        // states, and formation-reference-safe.
        let cl2 = wl[I_G1] * wl[4] / wl[I_RHO];
        let cr2 = wr[I_G1] * wr[4] / wr[I_RHO];
        ((sql * cl2 + sqr * cr2) * inv).sqrt()
    }

    fn stagnation_ghost(
        &self,
        _p0: f64,
        _rho0: f64,
        _c_frac: f64,
        _u_n: f64,
        _normal: usize,
    ) -> Result<Prim, FlowError> {
        // No closed-form isentrope on tabulated surfaces; the COUP-7 injector
        // object owns inflow. Refuse, never approximate (META-1 P6).
        Err(FlowError::BcUnsupportedByEos {
            bc: "StagnationInflow",
        })
    }

    /// Injector inflow: the injected propellant is **unburnt** (`b = 0`), so
    /// the face state is the unburnt branch's mass-flow ghost (which already
    /// sets the burn slot to 0). The reactants burn in the chamber.
    fn mass_flow_inflow_ghost(
        &self,
        mdot_per_area: f64,
        h_total: f64,
        c_frac: f64,
        p_int: f64,
        sign: f64,
        normal: usize,
    ) -> Result<Prim, FlowError> {
        self.unburnt
            .mass_flow_inflow_ghost(mdot_per_area, h_total, c_frac, p_int, sign, normal)
    }
}

/// SOLV-4 §3.6 v0.4.9 — the mid-transition scan's first-crossing CONTINUITY
/// resolution on synthetic closures (the production test surfaces are
/// monotone in the projection variable, so no genuine two-crossing blend
/// state is constructible from them; the resolution is pinned here instead,
/// and `project_pressure` consumes the same helper).
#[cfg(test)]
mod tests {
    use super::{N_P_SCAN, PScan, scan_first_crossing};

    /// The exact scan grid the helper walks.
    fn nodes(lo: f64, hi: f64) -> Vec<f64> {
        let ratio = hi / lo;
        (1..=N_P_SCAN)
            .map(|k| lo * ratio.powf(k as f64 / N_P_SCAN as f64))
            .collect()
    }

    /// The bracket the first-hit loop stores for a root at `root`.
    fn first_bracket(lo: f64, hi: f64, root: f64, g: impl Fn(f64) -> f64) -> PScan {
        let ns = nodes(lo, hi);
        let k = ns.iter().position(|&p| p > root).expect("crossing node");
        let a = if k == 0 { lo } else { ns[k - 1] };
        PScan::Bracket {
            a,
            b: ns[k],
            ga: g(a),
            gb: g(ns[k]),
        }
    }

    #[test]
    fn single_crossing_returns_the_first_hit_bracket_bit_exactly() {
        // Monotone g with one sign change: the returned bracket must be the
        // SAME floats the first-hit loop stores (bit-identity of every
        // accepted root rides on this).
        let (lo, hi) = (1.0e2, 1.0e6);
        let g = |p: f64| -> Result<f64, &'static str> { Ok(p - 3.7e4) };
        let out = scan_first_crossing(lo, hi, g(lo).unwrap(), &g).expect("single root scans");
        let want = first_bracket(lo, hi, 3.7e4, |p| p - 3.7e4);
        assert_eq!(out, want, "first-crossing bracket must be bit-identical");
    }

    #[test]
    fn two_crossings_resolve_to_the_first_by_continuity() {
        // A folded g (two admissible roots): v0.4.9 resolves by CONTINUITY —
        // the first crossing from the cold end is the accepted branch, never
        // a refusal (Ben ruling 2026-08-31; PLAN §8 v1.16 (3)).
        let (lo, hi) = (1.0, 1.0e4);
        let gf = |p: f64| (p - 30.0) * (p - 3.0e3);
        let g = |p: f64| -> Result<f64, &'static str> { Ok(gf(p)) };
        let out = scan_first_crossing(lo, hi, g(lo).unwrap(), &g)
            .expect("a folded g resolves by continuity (SOLV-4 §3.6 v0.4.9)");
        assert_eq!(out, first_bracket(lo, hi, 30.0, gf));
    }

    #[test]
    fn three_crossings_also_resolve_to_the_first() {
        // Multiplicity beyond two: the same first-crossing resolution.
        let (lo, hi) = (1.0, 1.0e4);
        let gf = |p: f64| (p - 12.0) * (p - 350.0) * (p - 4.7e3);
        let g = |p: f64| -> Result<f64, &'static str> { Ok(gf(p)) };
        let out = scan_first_crossing(lo, hi, g(lo).unwrap(), &g).expect("three roots resolve");
        assert_eq!(out, first_bracket(lo, hi, 12.0, gf));
    }

    #[test]
    fn no_crossing_returns_the_best_sample() {
        // One-signed g: the tangency seed is the |g|-minimal node over the
        // WHOLE sweep.
        let (lo, hi) = (1.0, 1.0e4);
        let g =
            |p: f64| -> Result<f64, &'static str> { Ok((p.ln() - 100.0f64.ln()).powi(2) + 0.5) };
        let out = scan_first_crossing(lo, hi, g(lo).unwrap(), &g).expect("no-crossing scans");
        let ns = nodes(lo, hi);
        let p_want = ns
            .iter()
            .copied()
            .fold((lo, g(lo).unwrap().abs()), |acc, p| {
                let gp = g(p).unwrap().abs();
                if gp < acc.1 { (p, gp) } else { acc }
            })
            .0;
        assert_eq!(out, PScan::NoCrossing { p_best: p_want });
    }

    #[test]
    fn exact_grid_zero_before_any_crossing_is_the_root() {
        // g exactly zero AT a scan node with no prior sign change: returned
        // as the root.
        let (lo, hi) = (1.0, 1.0e4);
        let p0 = nodes(lo, hi)[17];
        let g = move |p: f64| -> Result<f64, &'static str> { Ok(p - p0) };
        let out = scan_first_crossing(lo, hi, g(lo).unwrap(), &g).expect("exact zero scans");
        assert_eq!(out, PScan::ExactRoot(p0));
    }

    #[test]
    fn exact_grid_zero_after_a_crossing_is_never_visited() {
        // A sign change, then an exact grid zero further on: the scan stops
        // at the first crossing (continuity) — the later root is a further
        // admissible branch the resolution does not take.
        let (lo, hi) = (1.0, 1.0e4);
        let p1 = nodes(lo, hi)[40];
        let gf = move |p: f64| (p - 5.0) * (p - p1);
        let g = move |p: f64| -> Result<f64, &'static str> { Ok(gf(p)) };
        let out = scan_first_crossing(lo, hi, g(lo).unwrap(), &g).expect("first crossing accepted");
        assert_eq!(out, first_bracket(lo, hi, 5.0, gf));
    }

    #[test]
    fn g_error_propagates_from_any_visited_node() {
        // An Err from g refuses the scan when it lands BEFORE the first
        // crossing (a visited node); a node past the accepted bracket is
        // never visited — the sweep stops at the crossing.
        let (lo, hi) = (1.0, 1.0e4);
        let ns = nodes(lo, hi);
        let (early, late) = (ns[3], ns[50]);
        let g_pre = move |p: f64| -> Result<f64, &'static str> {
            if p == early {
                Err("bad point")
            } else {
                Ok(p - 5.0e3)
            }
        };
        assert_eq!(
            scan_first_crossing(lo, hi, g_pre(lo).unwrap(), &g_pre),
            Err("bad point")
        );
        let g_post = move |p: f64| -> Result<f64, &'static str> {
            if p == late {
                Err("bad point")
            } else {
                Ok(p - 5.0)
            }
        };
        assert_eq!(
            scan_first_crossing(lo, hi, g_post(lo).unwrap(), &g_post),
            Ok(first_bracket(lo, hi, 5.0, |p| p - 5.0))
        );
    }
}
