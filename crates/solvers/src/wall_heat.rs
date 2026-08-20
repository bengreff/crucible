//! SOLV-1 §3.5 — the **one local wall-function heat law**, every wall,
//! every engine (D-C). Gas-side convective flux from the local near-wall
//! state only: a Reynolds-analogy/Colburn-class wall function
//! [META-3: `wall-function-heat`, Kays–Crawford–Weigand lineage]:
//!
//! - `Re_y = ρ·|u_t|·y/μ` on the wall distance `y` (half the wall-adjacent
//!   cell's normal size);
//! - `C_f/2 = CF2_COEFF · Re_y^{−1/4}` (the 1/7-power-law wall function);
//! - `St = (C_f/2) · Pr^{−2/3}` (the Colburn analogy);
//! - `h = k_gas/y + St·ρ·|u_t|·c_p` — the first term is the molecular
//!   conduction floor (`Nu_y = 1`), which makes the law smooth and defined
//!   down to `u_t = 0` with no regime branch (Rule 12); the convective term
//!   vanishes as `u_t^{3/4}`;
//! - `T_aw = T + Pr^{1/3}·u_t²/(2 c_p)` (turbulent recovery factor);
//! - the exchanged flux couples through the film + solid half-cell series
//!   resistance: `q = (T_aw − T_s)/(1/h + d_s/κ_s)` with `T_s` the solid
//!   surface cell's center value and `d_s` its center-to-face distance —
//!   the standard consistent CHT flux form (COUP-2 §3.5), evaluated ONCE
//!   per face so both sides exchange the identical number (conservation by
//!   construction).
//!
//! **Declared band: ±20–30%** on `h` (SOLV-1 §3.5; the closure family
//! spread — recorded in the station-4 certificate, consumed by COUP-5 as a
//! closure band, never folded into the returned value). Transport
//! (μ, c_p, Pr) is constant config data this session — the degenerate
//! FND-7 spine occupant, exactly the gamma-law pattern (session 7); the
//! spine supplies per-cell transport when OFFL-5 lands.

use crucible_units as units;

/// The classic power-law wall-function skin-friction coefficient
/// (Kays–Crawford; the 1/7-velocity-profile constant).
pub const CF2_COEFF: f64 = 0.0225;
/// `C_f/2 ∝ Re_y^{CF2_RE_EXP}`.
pub const CF2_RE_EXP: f64 = -0.25;
/// Colburn: `St·Pr^{2/3} = C_f/2`.
pub const COLBURN_PR_EXP: f64 = -2.0 / 3.0;
/// Turbulent recovery factor `r = Pr^{1/3}`.
pub const RECOVERY_PR_EXP: f64 = 1.0 / 3.0;
/// The declared closure band on `h` (SOLV-1 §3.5): ±20–30%.
pub const DECLARED_BAND: (f64, f64) = (0.20, 0.30);

#[derive(Debug, Clone, PartialEq)]
pub enum WallHeatError {
    /// An operand outside the law's domain (non-finite, or ≤ 0 where
    /// positivity is physical) — refuse, never clamp (META-1 P6).
    BadOperand { which: &'static str, value: f64 },
}

impl std::fmt::Display for WallHeatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BadOperand { which, value } => write!(
                f,
                "wall-function operand {which} = {value} outside the law's domain — \
                 refusing, never clamping (SOLV-1 §3.5, META-1 P6)"
            ),
        }
    }
}

impl std::error::Error for WallHeatError {}

/// The law's constant transport data (degenerate spine occupant).
#[derive(Debug, Clone, Copy)]
pub struct WallLaw {
    cp: f64, // J/(kg·K)
    mu: f64, // Pa·s
    pr: f64, // 1
}

/// Local near-wall gas state, SI (assembled by the coupler from `U`).
#[derive(Debug, Clone, Copy)]
pub struct NearWallGas {
    pub rho: f64,
    /// Wall-tangential speed |u_t| ≥ 0.
    pub u_t: f64,
    pub temperature: f64,
    /// Wall distance y (half the wall-adjacent cell's normal size).
    pub y: f64,
}

impl WallLaw {
    /// Unit-typed at the boundary (META-2 §4 ★); Pr is dimensionless.
    pub fn new(
        cp: units::SpecificHeatCapacity,
        mu: units::DynamicViscosity,
        pr: f64,
    ) -> Result<WallLaw, WallHeatError> {
        let law = WallLaw {
            cp: units::si(cp),
            mu: units::si(mu),
            pr,
        };
        for (which, v) in [("cp", law.cp), ("mu", law.mu), ("pr", law.pr)] {
            if !(v.is_finite() && v > 0.0) {
                return Err(WallHeatError::BadOperand { which, value: v });
            }
        }
        Ok(law)
    }

    fn check_gas(gas: &NearWallGas) -> Result<(), WallHeatError> {
        for (which, v, positive) in [
            ("rho", gas.rho, true),
            ("u_t", gas.u_t, false),
            ("temperature", gas.temperature, true),
            ("y", gas.y, true),
        ] {
            if !(v.is_finite() && (v > 0.0 || (!positive && v >= 0.0))) {
                return Err(WallHeatError::BadOperand { which, value: v });
            }
        }
        Ok(())
    }

    /// Gas conductivity from the constant transport set: k = μ·c_p/Pr.
    pub fn k_gas(&self) -> f64 {
        self.mu * self.cp / self.pr
    }

    /// The constant transport set, read back (SI). The wall law is the ONE
    /// owner of these constants (Rule 13); the S3 gas-diffusion operator
    /// derives its μ/k/ρD/c_v from here — never a second statement.
    pub fn cp(&self) -> f64 {
        self.cp
    }

    pub fn mu(&self) -> f64 {
        self.mu
    }

    pub fn pr(&self) -> f64 {
        self.pr
    }

    /// Film coefficient h [W/(m²·K)] from the local state — the law itself.
    pub fn film_coefficient(&self, gas: &NearWallGas) -> Result<f64, WallHeatError> {
        Self::check_gas(gas)?;
        let floor = self.k_gas() / gas.y;
        if gas.u_t == 0.0 {
            return Ok(floor);
        }
        let re_y = gas.rho * gas.u_t * gas.y / self.mu;
        let cf2 = CF2_COEFF * re_y.powf(CF2_RE_EXP);
        let st = cf2 * self.pr.powf(COLBURN_PR_EXP);
        Ok(floor + st * gas.rho * gas.u_t * self.cp)
    }

    /// Adiabatic-wall (recovery) temperature [K].
    pub fn adiabatic_wall_temperature(&self, gas: &NearWallGas) -> Result<f64, WallHeatError> {
        Self::check_gas(gas)?;
        Ok(gas.temperature + self.pr.powf(RECOVERY_PR_EXP) * gas.u_t * gas.u_t / (2.0 * self.cp))
    }

    /// The exchanged wall flux q [W/m²] (positive gas→solid) through the
    /// film + solid half-cell series resistance, with the diagnostics the
    /// coupler and certificate consume. `t_solid` is the surface solid
    /// cell's center temperature, `d_solid` its center-to-face distance,
    /// `kappa_solid` its conductivity.
    pub fn wall_exchange(
        &self,
        gas: &NearWallGas,
        t_solid: f64,
        d_solid: f64,
        kappa_solid: f64,
    ) -> Result<WallExchange, WallHeatError> {
        for (which, v) in [
            ("t_solid", t_solid),
            ("d_solid", d_solid),
            ("kappa_solid", kappa_solid),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(WallHeatError::BadOperand { which, value: v });
            }
        }
        let h = self.film_coefficient(gas)?;
        let t_aw = self.adiabatic_wall_temperature(gas)?;
        let u_series = 1.0 / (1.0 / h + d_solid / kappa_solid);
        Ok(WallExchange {
            q: u_series * (t_aw - t_solid),
            h,
            t_aw,
            u_series,
        })
    }

    /// The flux alone (see [`WallLaw::wall_exchange`]).
    pub fn wall_flux(
        &self,
        gas: &NearWallGas,
        t_solid: f64,
        d_solid: f64,
        kappa_solid: f64,
    ) -> Result<f64, WallHeatError> {
        Ok(self.wall_exchange(gas, t_solid, d_solid, kappa_solid)?.q)
    }
}

/// One evaluated wall exchange: the flux, and the pieces the coupler's
/// stability limit and the certificate's diagnostics need.
#[derive(Debug, Clone, Copy)]
pub struct WallExchange {
    /// Exchanged flux [W/m²], positive gas→solid.
    pub q: f64,
    /// Gas-side film coefficient [W/(m²·K)].
    pub h: f64,
    /// Adiabatic-wall (recovery) temperature [K].
    pub t_aw: f64,
    /// The series coefficient q/(T_aw − T_s) [W/(m²·K)].
    pub u_series: f64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crucible_units::{dynamic_viscosity_pa_s, specific_heat_capacity_j_per_kg_k};

    fn law() -> WallLaw {
        WallLaw::new(
            specific_heat_capacity_j_per_kg_k(1004.0),
            dynamic_viscosity_pa_s(4.0e-5),
            0.72,
        )
        .expect("valid transport set")
    }

    #[test]
    fn colburn_algebra_by_hand() {
        // One case computed independently: ρ=2, u=300, y=1e-3, μ=4e-5,
        // cp=1004, Pr=0.72.
        let l = law();
        let gas = NearWallGas {
            rho: 2.0,
            u_t: 300.0,
            temperature: 800.0,
            y: 1.0e-3,
        };
        let re_y: f64 = 2.0 * 300.0 * 1.0e-3 / 4.0e-5; // 1.5e4
        let cf2 = 0.0225 * re_y.powf(-0.25);
        let st = cf2 * 0.72f64.powf(-2.0 / 3.0);
        let floor = 4.0e-5 * 1004.0 / 0.72 / 1.0e-3;
        let expect = floor + st * 2.0 * 300.0 * 1004.0;
        assert_eq!(l.film_coefficient(&gas).unwrap(), expect);
        let t_aw = 800.0 + 0.72f64.powf(1.0 / 3.0) * 300.0 * 300.0 / (2.0 * 1004.0);
        assert_eq!(l.adiabatic_wall_temperature(&gas).unwrap(), t_aw);
    }

    #[test]
    fn static_gas_reduces_to_molecular_conduction() {
        let l = law();
        let gas = NearWallGas {
            rho: 1.0,
            u_t: 0.0,
            temperature: 500.0,
            y: 2.0e-3,
        };
        assert_eq!(l.film_coefficient(&gas).unwrap(), l.k_gas() / 2.0e-3);
        // No kinetic recovery at rest: T_aw = T exactly.
        assert_eq!(l.adiabatic_wall_temperature(&gas).unwrap(), 500.0);
    }

    #[test]
    fn equal_temperatures_at_rest_exchange_exactly_zero() {
        // The coupled fixed-point requirement: uniform T ⇒ q = 0 exactly.
        let l = law();
        let gas = NearWallGas {
            rho: 1.0,
            u_t: 0.0,
            temperature: 400.0,
            y: 1.0e-3,
        };
        assert_eq!(l.wall_flux(&gas, 400.0, 5.0e-4, 20.0).unwrap(), 0.0);
    }

    #[test]
    fn hot_gas_drives_flux_into_the_wall_and_scales_with_speed() {
        let l = law();
        let slow = NearWallGas {
            rho: 2.0,
            u_t: 50.0,
            temperature: 900.0,
            y: 1.0e-3,
        };
        let fast = NearWallGas { u_t: 400.0, ..slow };
        let q_slow = l.wall_flux(&slow, 400.0, 5.0e-4, 20.0).unwrap();
        let q_fast = l.wall_flux(&fast, 400.0, 5.0e-4, 20.0).unwrap();
        assert!(q_slow > 0.0);
        assert!(
            q_fast > q_slow,
            "convection must grow with tangential speed"
        );
    }

    #[test]
    fn refusals_fire_on_bad_operands() {
        let l = law();
        let good = NearWallGas {
            rho: 1.0,
            u_t: 10.0,
            temperature: 500.0,
            y: 1.0e-3,
        };
        for bad in [
            NearWallGas { rho: -1.0, ..good },
            NearWallGas { u_t: -3.0, ..good },
            NearWallGas {
                temperature: f64::NAN,
                ..good
            },
            NearWallGas { y: 0.0, ..good },
        ] {
            assert!(l.film_coefficient(&bad).is_err());
        }
        assert!(l.wall_flux(&good, -5.0, 5.0e-4, 20.0).is_err());
        assert!(
            WallLaw::new(
                specific_heat_capacity_j_per_kg_k(1004.0),
                dynamic_viscosity_pa_s(0.0),
                0.72
            )
            .is_err()
        );
    }
}
