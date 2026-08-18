//! SOLV-1 §3.4 — **shifting-equilibrium mode**: the [`EosLaw`] occupant
//! backed by the OFFL-3 `(p, h, Z)` equilibrium surface. The composition
//! slot of `U` carries the **elemental mixture fraction Z** (advected,
//! source-free — element conservation exact by construction); every
//! conserved→primitive conversion performs the **per-cell equilibrium
//! projection**: solve for the pressure `p` at which the tabulated
//! equilibrium density at `(p, h = e + p/ρ, Z)` matches the cell's conserved
//! density. Species, temperature, and heat release are whatever the surface
//! says at that state — combustion happens *in the EOS*, never as a separate
//! `if(burning)` branch (Rule 12).
//!
//! **The p-iteration** (SOLV-1 §3.4): `h = e + p/ρ` is linear in `p` at
//! fixed conserved state, so the admissible pressure bracket (the p-range
//! keeping both `p` and `h` inside the declared table envelope) is closed
//! form; the root of `ρ_tab(p, h(p), Z) − ρ` is then located by a
//! **deterministic Illinois regula-falsi** with a fixed relative tolerance
//! [`EPS_P_PROJECTION`] and fixed max iteration count [`N_P_ITER_MAX`] —
//! the META-1 §2.2 sanctioned convergence form (absolute + deterministic;
//! never wall-clock). Failure to bracket or converge is a halt diagnosis,
//! never a clamp: a state off the tabulated surface refuses (FND-5 §3.5).
//!
//! **η_c\* knockdown hook (S18, COUP-7 §3.2.1):** [`TableEos::h_offset`] is
//! the source-term-level combustion-completeness knockdown — an enthalpy
//! offset added to the `h`-coordinate of every equilibrium interrogation
//! (δh < 0 ⇒ the surface is read as if that much reaction enthalpy had not
//! been released, so T, a, Γ₁ and the emergent p_c all carry the deficit
//! self-consistently). The conserved energy is untouched (the deficit is
//! unreleased chemical energy riding in `U`); output-side multiplication of
//! c\*/Isp is forbidden (SOLV-1 §3.4). The COUP-7 injector prior tier sets
//! δh (calibrated so delivered c\* = η_c\*·c\*_ideal at the anchor state);
//! 0.0 = full equilibrium.

use super::{Cons, EosLaw, FlowError, I_EI, I_EN, I_G1, I_MR, I_MT, I_MZ, I_RC, I_RHO, Prim};
use crucible_tables::{BoundColumn, Table, TableError};

/// Fixed maximum iteration count of the equilibrium pressure projection.
/// Illinois regula-falsi typically converges in < 15; plain bisection would
/// need ~38 to reach [`EPS_P_PROJECTION`] over the full table p-envelope, so
/// 48 is headroom, not a tunable.
pub const N_P_ITER_MAX: usize = 48;

/// Relative bracket-width tolerance of the pressure projection — orders
/// below the table's own interpolation-error bound, so the projection never
/// contributes to the physics error budget.
pub const EPS_P_PROJECTION: f64 = 1e-11;

/// Inward relative shrink applied to the closed-form h-derived pressure
/// bracket so endpoint queries cannot fall an ulp outside the envelope.
const H_BRACKET_MARGIN: f64 = 1e-12;

/// The shifting-equilibrium EOS occupant: bound columns of one OFFL-3
/// `(p, h, Z)` equilibrium surface (columns resolved once — `BoundColumn`,
/// FND-5 §3.3 — for the ~10⁸ projections of an anchor run).
#[derive(Debug, Clone)]
pub struct TableEos<'t> {
    rho: BoundColumn<'t>,
    sound: BoundColumn<'t>,
    temperature: BoundColumn<'t>,
    /// Declared envelopes of the (p, h, Z) axes, cached at bind.
    p_env: (f64, f64),
    h_env: (f64, f64),
    z_env: (f64, f64),
    /// S18 knockdown: enthalpy offset (J/kg) added to the h-coordinate of
    /// every equilibrium interrogation. See the module header.
    pub h_offset: f64,
}

impl<'t> TableEos<'t> {
    /// Bind the occupant to a loaded equilibrium surface. Refuses (never
    /// guesses) if the table is not the `(p, h, Z)` schema or a column is
    /// missing/relabeled (META-2 §4 ★ units gate runs per column).
    pub fn bind(table: &'t Table) -> Result<Self, String> {
        let rho = table
            .bind("density", "kg/m^3")
            .map_err(|e| format!("equilibrium surface: {e}"))?;
        let sound = table
            .bind("sound_speed", "m/s")
            .map_err(|e| format!("equilibrium surface: {e}"))?;
        let temperature = table
            .bind("temperature", "K")
            .map_err(|e| format!("equilibrium surface: {e}"))?;
        let expected = ["p", "h", "Z"];
        for (i, want) in expected.iter().enumerate() {
            let got = rho.axis_name(i);
            if got != *want {
                return Err(format!(
                    "equilibrium surface axis {i} is {got:?}, expected {want:?} — \
                     not a (p, h, Z) local-state surface (OFFL-3 §3.3); refusing the bind"
                ));
            }
        }
        Ok(Self {
            p_env: rho.axis_envelope(0),
            h_env: rho.axis_envelope(1),
            z_env: rho.axis_envelope(2),
            rho,
            sound,
            temperature,
            h_offset: 0.0,
        })
    }

    /// The equilibrium pressure projection: root of
    /// `ρ_tab(p, e_q + p/ρ, Z) − ρ` over the admissible bracket, where
    /// `e_q = e + h_offset`. Deterministic Illinois regula-falsi (module
    /// header). Returns the located pressure.
    fn project_pressure(&self, rho: f64, e_q: f64, z: f64) -> Result<f64, &'static str> {
        let inv = 1.0 / rho;
        // Closed-form admissible bracket: h(p) = e_q + p/ρ must lie in the
        // h envelope AND p in the p envelope; shrink the h-derived bounds
        // inward so endpoint rounding cannot leave the envelope.
        let p_from_h_lo = rho * (self.h_env.0 - e_q);
        let p_from_h_hi = rho * (self.h_env.1 - e_q);
        let mut lo = self.p_env.0.max(p_from_h_lo * (1.0 + H_BRACKET_MARGIN));
        let hi = self.p_env.1.min(p_from_h_hi * (1.0 - H_BRACKET_MARGIN));
        if !(lo.is_finite() && hi.is_finite()) || lo >= hi || hi <= 0.0 {
            return Err("no admissible pressure bracket: state outside the table envelope");
        }
        lo = lo.max(0.0);

        let g = |p: f64| -> Result<f64, &'static str> {
            self.rho
                .interpolate(&[p, e_q + p * inv, z])
                .map(|r| r - rho)
                .map_err(|_| "equilibrium surface query failed inside the projection bracket")
        };

        let mut ga = g(lo)?;
        let mut gb = g(hi)?;
        if ga == 0.0 {
            return Ok(lo);
        }
        if gb == 0.0 {
            return Ok(hi);
        }
        if ga * gb > 0.0 {
            return Err("no equilibrium state in the table envelope for this (rho, e, Z)");
        }
        let (mut a, mut b) = (lo, hi);
        for _ in 0..N_P_ITER_MAX {
            if (b - a).abs() <= EPS_P_PROJECTION * a.abs().max(b.abs()) {
                return Ok(0.5 * (a + b));
            }
            // Regula-falsi iterate, bisection fallback if it degenerates.
            let denom = gb - ga;
            let mut p = if denom != 0.0 {
                b - gb * (b - a) / denom
            } else {
                0.5 * (a + b)
            };
            let (loe, hie) = (a.min(b), a.max(b));
            if !(p > loe && p < hie) {
                p = 0.5 * (a + b);
            }
            let gp = g(p)?;
            if gp == 0.0 {
                return Ok(p);
            }
            if gp * gb < 0.0 {
                a = b;
                ga = gb;
            } else {
                // Illinois: halve the stale endpoint's residual so the
                // bracket cannot stagnate on one side.
                ga *= 0.5;
            }
            b = p;
            gb = gp;
        }
        if (b - a).abs() <= EPS_P_PROJECTION * a.abs().max(b.abs()) * 10.0 {
            return Ok(0.5 * (a + b));
        }
        Err("equilibrium pressure projection did not converge in the fixed iteration budget")
    }

    /// Equilibrium temperature at a primitive state produced by this
    /// occupant (diagnostics + the wall-law operand feed).
    pub fn temperature_w(&self, w: &Prim) -> Result<f64, TableError> {
        let (rho, p, z) = (w[I_RHO], w[4], w[I_RC]);
        let h = w[I_EI] + self.h_offset + p / rho;
        self.temperature.interpolate(&[p, h, z])
    }

    /// Conserved state on the equilibrium surface from `(p, h, Z)` and a
    /// velocity — the constructor initial conditions and the COUP-7 injector
    /// inflow use (the exact inverse of the projection at `h_offset = 0`;
    /// with a knockdown in force, `h` is the *query* coordinate and the
    /// stored energy carries the deficit consistently).
    pub fn cons_from_phz(&self, p: f64, h: f64, z: f64, vel: [f64; 3]) -> Result<Cons, TableError> {
        let rho = self.rho.interpolate(&[p, h, z])?;
        let e_true = h - p / rho - self.h_offset;
        let ke = 0.5 * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
        Ok([
            rho,
            rho * vel[0],
            rho * vel[1],
            rho * vel[2],
            rho * (e_true + ke),
            rho * z,
        ])
    }

    /// The declared (p, h, Z) envelopes (assembly-time sizing/refusals).
    pub fn envelopes(&self) -> [(f64, f64); 3] {
        [self.p_env, self.h_env, self.z_env]
    }
}

impl EosLaw for TableEos<'_> {
    fn prim_checked(&self, u: &Cons) -> Result<Prim, &'static str> {
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
        if z < self.z_env.0 || z > self.z_env.1 {
            return Err("elemental mixture fraction outside the table envelope");
        }
        let e_q = e + self.h_offset;
        let p = self.project_pressure(rho, e_q, z)?;
        let h = e_q + p * inv;
        let a = self
            .sound
            .interpolate(&[p, h, z])
            .map_err(|_| "sound-speed query failed at the projected state")?;
        let g1 = rho * a * a / p;
        Ok([rho, ur, ut, uz, p, z, e, g1])
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
        h_roe: f64,
        q2_roe: f64,
        sql: f64,
        sqr: f64,
        inv: f64,
    ) -> f64 {
        // Roe-averaged Γ₁ — the standard general-convex-EOS extension of the
        // Batten bounds (the aux slot is reconstructed like every primitive).
        let g1_roe = (sql * wl[I_G1] + sqr * wr[I_G1]) * inv;
        ((g1_roe - 1.0) * (h_roe - 0.5 * q2_roe)).max(0.0).sqrt()
    }

    fn stagnation_ghost(
        &self,
        _p0: f64,
        _rho0: f64,
        _c_frac: f64,
        _u_n: f64,
        _normal: usize,
    ) -> Result<Prim, FlowError> {
        // No closed-form isentrope on a tabulated surface; the COUP-7
        // injector object owns inflow in shifting mode. Refuse, never
        // approximate (META-1 P6).
        Err(FlowError::BcUnsupportedByEos {
            bc: "StagnationInflow",
        })
    }
}
