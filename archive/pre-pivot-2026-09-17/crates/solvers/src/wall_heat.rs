//! SOLV-1 §3.5 — the **one local wall-function heat law**, every wall,
//! every engine (D-C). Gas-side convective flux from the local near-wall
//! state only: a Reynolds-analogy/Colburn-class wall function
//! [META-3: `wall-function-heat`, Kays–Crawford–Weigand lineage]:
//!
//! - `Re_y = ρ·|u_t|·y/μ` on the wall distance `y` (half the wall-adjacent
//!   cell's normal size);
//! - `C_f/2 = CF2_COEFF · Re_y^{−1/4}` (the 1/7-power-law wall function);
//! - `St = (C_f/2) · Pr^{−2/3}` (the Colburn analogy);
//! - `h = k/y + St·ρ·|u_t|·c_p,film` — the first term is the conduction
//!   floor (`Nu_y = 1`), which makes the law smooth and defined down to
//!   `u_t = 0` with no regime branch (Rule 12); the convective term
//!   vanishes as `u_t^{3/4}`. **The convective limb drives on the
//!   film-mean slope `c_p,film`, not the local `c_p`** — the analogy
//!   transports enthalpy (`q_w = St·ρu·(h_aw − h_w)`), and `c_p·ΔT` is
//!   `Δh` only for the mean slope across the film. The local equilibrium
//!   `c_p` is the peak of a strongly-peaked curve and using it here
//!   overpredicts `q_w` by ~1.9× at RL10 chamber conditions (S4 review
//!   finding E1; the spine owns the distinction — `crate::transport`);
//! - `T_aw = T + Pr^{1/3}·u_t²/(2 c_p)` — recovery uses the **local**
//!   slope, because it converts a recovery enthalpy into a temperature at
//!   the boundary-layer EDGE state, where the gas is at its own
//!   equilibrium composition. Two different c_p's, two different
//!   questions, both from the one provider;
//! - the exchanged flux couples through the film + solid half-cell series
//!   resistance: `q = (T_aw − T_s)/(1/h + d_s/κ_s)` with `T_s` the solid
//!   surface cell's center value and `d_s` its center-to-face distance —
//!   the standard consistent CHT flux form (COUP-2 §3.5), evaluated ONCE
//!   per face so both sides exchange the identical number (conservation by
//!   construction).
//!
//! **Declared band: ±20–30%** on `h` (SOLV-1 §3.5; the closure family
//! spread — recorded in the station-4 certificate, consumed by COUP-5 as a
//! closure band, never folded into the returned value).
//!
//! ## Transport (S4 — the law states no constant of its own)
//!
//! μ, k, c_p and Pr arrive as a [`TransportProps`] **operand**, queried from
//! the FND-7 §3.3 spine at the near-wall gas cell's own state
//! ([`crate::transport`]) — the same one provider the resolved `F_visc`
//! next door reads. Before S4 this law privately held `(c_p, μ, Pr)` as
//! config constants; that was the degenerate spine occupant, and it is now
//! [`crate::transport::ConstantTransport`], stated once.
//!
//! Evaluating the operands at the **gas** state rather than at a film or
//! Eckert reference temperature is a declared choice inside the ±20–30%
//! band (SOLV-1 §3.5, 0.4.1) — the band is the closure-family spread and
//! reference-temperature variants sit within it. The refinement is a
//! recorded deferral, not a silent approximation.
//!
//! `WallLaw` is therefore **stateless**: pure algebra over declared
//! constants and its operands. That is the point — a wall law that carried
//! properties would be a second owner of them.

use crate::transport::TransportProps;

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

/// The one wall-function heat law. Since S4 it holds no transport — that
/// arrives as a [`TransportProps`] operand from the FND-7 spine (module
/// doc) — only its **declared band coordinate**: `band` multiplies `h`,
/// realizing the ±20–30% Colburn-family band of [`DECLARED_BAND`] so a
/// p-box corner member is a labeled config value rather than a proxy knob.
/// `band = 1` is the law itself.
#[derive(Debug, Clone, Copy)]
pub struct WallLaw {
    band: f64,
}

impl Default for WallLaw {
    fn default() -> Self {
        WallLaw::new()
    }
}

/// Local near-wall gas state, SI (assembled by the coupler from `U`), with
/// the spine's transport at that same state. One struct, because the law's
/// contract is "everything local, nothing remembered" — a transport set
/// that arrived separately could silently belong to a different cell.
#[derive(Debug, Clone, Copy)]
pub struct NearWallGas {
    pub rho: f64,
    /// Wall-tangential speed |u_t| ≥ 0.
    pub u_t: f64,
    pub temperature: f64,
    /// Wall distance y (half the wall-adjacent cell's normal size).
    pub y: f64,
    /// Spine transport at this cell's state (FND-7 §3.3).
    pub tr: TransportProps,
}

impl WallLaw {
    /// The law at its nominal band coordinate.
    pub fn new() -> WallLaw {
        WallLaw { band: 1.0 }
    }

    /// The law at a declared band corner (module doc). Refuses a factor
    /// that is not a positive finite number — a band coordinate that is
    /// not a number would silently zero or NaN every wall flux.
    pub fn with_band(band: f64) -> Result<WallLaw, WallHeatError> {
        if !(band.is_finite() && band > 0.0) {
            return Err(WallHeatError::BadOperand {
                which: "band_factor",
                value: band,
            });
        }
        Ok(WallLaw { band })
    }

    /// The declared band coordinate this law is evaluated at.
    pub fn band(&self) -> f64 {
        self.band
    }

    fn check_gas(gas: &NearWallGas) -> Result<(), WallHeatError> {
        for (which, v, positive) in [
            ("rho", gas.rho, true),
            ("u_t", gas.u_t, false),
            ("temperature", gas.temperature, true),
            ("y", gas.y, true),
            // The spine validates its own outputs (it refuses rather than
            // returning a bad one), but this law is also called directly by
            // fixtures that hand-build a `TransportProps` — so the operands
            // it divides and raises to fractional powers are checked here
            // too. Defense in depth at a coupler-rate call is free.
            ("cp", gas.tr.cp, true),
            ("cp_film", gas.tr.cp_film, true),
            ("mu", gas.tr.mu, true),
            ("k", gas.tr.k, true),
            ("pr", gas.tr.pr, true),
        ] {
            if !(v.is_finite() && (v > 0.0 || (!positive && v >= 0.0))) {
                return Err(WallHeatError::BadOperand { which, value: v });
            }
        }
        Ok(())
    }

    /// Film coefficient h [W/(m²·K)] from the local state — the law itself.
    pub fn film_coefficient(&self, gas: &NearWallGas) -> Result<f64, WallHeatError> {
        Self::check_gas(gas)?;
        let floor = gas.tr.k / gas.y;
        if gas.u_t == 0.0 {
            return Ok(self.band * floor);
        }
        let re_y = gas.rho * gas.u_t * gas.y / gas.tr.mu;
        let cf2 = CF2_COEFF * re_y.powf(CF2_RE_EXP);
        let st = cf2 * gas.tr.pr.powf(COLBURN_PR_EXP);
        // The band multiplies the WHOLE film coefficient — conduction floor
        // included — because it is a band on h, not on one of its limbs.
        Ok(self.band * (floor + st * gas.rho * gas.u_t * gas.tr.cp_film))
    }

    /// Adiabatic-wall (recovery) temperature [K].
    pub fn adiabatic_wall_temperature(&self, gas: &NearWallGas) -> Result<f64, WallHeatError> {
        Self::check_gas(gas)?;
        Ok(gas.temperature
            + gas.tr.pr.powf(RECOVERY_PR_EXP) * gas.u_t * gas.u_t / (2.0 * gas.tr.cp))
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
        WallLaw::new()
    }

    /// The transport operand the fixtures use — the declared-constant
    /// spine occupant, exactly as a config would build it.
    fn tr() -> TransportProps {
        crate::transport::ConstantTransport::new(
            specific_heat_capacity_j_per_kg_k(1004.0),
            dynamic_viscosity_pa_s(4.0e-5),
            0.72,
            1.4,
            0.7,
        )
        .expect("valid transport set")
        .into_props()
    }

    fn gas(rho: f64, u_t: f64, temperature: f64, y: f64) -> NearWallGas {
        NearWallGas {
            rho,
            u_t,
            temperature,
            y,
            tr: tr(),
        }
    }

    #[test]
    fn colburn_algebra_by_hand() {
        // One case computed independently: ρ=2, u=300, y=1e-3, μ=4e-5,
        // cp=1004, Pr=0.72.
        let l = law();
        let g = gas(2.0, 300.0, 800.0, 1.0e-3);
        let re_y: f64 = 2.0 * 300.0 * 1.0e-3 / 4.0e-5; // 1.5e4
        let cf2 = 0.0225 * re_y.powf(-0.25);
        let st = cf2 * 0.72f64.powf(-2.0 / 3.0);
        let floor = 4.0e-5 * 1004.0 / 0.72 / 1.0e-3;
        let expect = floor + st * 2.0 * 300.0 * 1004.0;
        assert_eq!(l.film_coefficient(&g).unwrap(), expect);
        let t_aw = 800.0 + 0.72f64.powf(1.0 / 3.0) * 300.0 * 300.0 / (2.0 * 1004.0);
        assert_eq!(l.adiabatic_wall_temperature(&g).unwrap(), t_aw);
    }

    #[test]
    fn static_gas_reduces_to_molecular_conduction() {
        let l = law();
        let g = gas(1.0, 0.0, 500.0, 2.0e-3);
        assert_eq!(l.film_coefficient(&g).unwrap(), tr().k / 2.0e-3);
        // No kinetic recovery at rest: T_aw = T exactly.
        assert_eq!(l.adiabatic_wall_temperature(&g).unwrap(), 500.0);
    }

    #[test]
    fn equal_temperatures_at_rest_exchange_exactly_zero() {
        // The coupled fixed-point requirement: uniform T ⇒ q = 0 exactly.
        let l = law();
        let g = gas(1.0, 0.0, 400.0, 1.0e-3);
        assert_eq!(l.wall_flux(&g, 400.0, 5.0e-4, 20.0).unwrap(), 0.0);
    }

    #[test]
    fn hot_gas_drives_flux_into_the_wall_and_scales_with_speed() {
        let l = law();
        let slow = gas(2.0, 50.0, 900.0, 1.0e-3);
        let fast = NearWallGas { u_t: 400.0, ..slow };
        let q_slow = l.wall_flux(&slow, 400.0, 5.0e-4, 20.0).unwrap();
        let q_fast = l.wall_flux(&fast, 400.0, 5.0e-4, 20.0).unwrap();
        assert!(q_slow > 0.0);
        assert!(
            q_fast > q_slow,
            "convection must grow with tangential speed"
        );
    }

    /// A hotter, more-conductive spine reading must move `h` — the whole
    /// point of S4 is that these are operands, not constants. Guards the
    /// refactor against a stale-transport regression (a law that ignored
    /// `gas.tr` would still pass every test above, which all share one
    /// transport set).
    #[test]
    fn film_coefficient_tracks_the_spine_operands() {
        let l = law();
        let base = gas(2.0, 300.0, 800.0, 1.0e-3);
        let h0 = l.film_coefficient(&base).unwrap();
        let hotter = NearWallGas {
            tr: TransportProps {
                k: base.tr.k * 4.0,
                ..base.tr
            },
            ..base
        };
        assert!(l.film_coefficient(&hotter).unwrap() > h0);
        let fatter = NearWallGas {
            tr: TransportProps {
                cp_film: base.tr.cp_film * 3.0,
                ..base.tr
            },
            ..base
        };
        assert!(l.film_coefficient(&fatter).unwrap() > h0);
        let recovering = NearWallGas {
            tr: TransportProps {
                cp: base.tr.cp * 3.0,
                ..base.tr
            },
            ..base
        };
        assert!(
            l.adiabatic_wall_temperature(&recovering).unwrap()
                < l.adiabatic_wall_temperature(&base).unwrap(),
            "a larger LOCAL c_p absorbs the same kinetic energy into a smaller \
             temperature rise"
        );
    }

    /// **S4 review finding E1, as a test.** The convective limb must drive
    /// on the film-mean slope and the recovery on the local one. A law that
    /// used a single c_p for both would pass every other test here, because
    /// they all run a constant occupant where the two coincide.
    #[test]
    fn the_two_heat_capacities_drive_different_limbs() {
        let l = law();
        let base = gas(2.0, 300.0, 800.0, 1.0e-3);
        // Raising the LOCAL c_p alone must not move h at all...
        let local_only = NearWallGas {
            tr: TransportProps {
                cp: base.tr.cp * 2.0,
                ..base.tr
            },
            ..base
        };
        assert_eq!(
            l.film_coefficient(&local_only).unwrap(),
            l.film_coefficient(&base).unwrap(),
            "the Colburn limb transports enthalpy across the film; the local \
             equilibrium slope is not the film mean"
        );
        // ...and raising the FILM-mean slope alone must not move T_aw.
        let film_only = NearWallGas {
            tr: TransportProps {
                cp_film: base.tr.cp_film * 2.0,
                ..base.tr
            },
            ..base
        };
        assert_eq!(
            l.adiabatic_wall_temperature(&film_only).unwrap(),
            l.adiabatic_wall_temperature(&base).unwrap(),
            "recovery is an edge-state property"
        );
    }

    #[test]
    fn refusals_fire_on_bad_operands() {
        let l = law();
        let good = gas(1.0, 10.0, 500.0, 1.0e-3);
        for bad in [
            NearWallGas { rho: -1.0, ..good },
            NearWallGas { u_t: -3.0, ..good },
            NearWallGas {
                temperature: f64::NAN,
                ..good
            },
            NearWallGas { y: 0.0, ..good },
            // A hand-built transport operand is checked too (module doc).
            NearWallGas {
                tr: TransportProps { mu: 0.0, ..good.tr },
                ..good
            },
            NearWallGas {
                tr: TransportProps {
                    k: f64::NAN,
                    ..good.tr
                },
                ..good
            },
        ] {
            assert!(l.film_coefficient(&bad).is_err());
        }
        assert!(l.wall_flux(&good, -5.0, 5.0e-4, 20.0).is_err());
    }
}
