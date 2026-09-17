//! META-2 §4 ★ — units are type-checked at interface boundaries (config
//! parse, mechanism-setup structs, public signatures, results); plain `f64`
//! in documented SI base units inside hot kernels (FND-1 §3.2).
//!
//! This crate is the **one owner** of the `uom` dependency (pinned in
//! Cargo.toml): boundary code names quantities through these re-exports and
//! constructors, so swapping the underlying library is a one-crate event.
//! Two rules of use:
//!
//! 1. **Construct** from a documented-SI `f64` with the `*_si` constructors
//!    below (or uom's `Quantity::new::<unit>` for non-base units).
//! 2. **Extract** at the kernel boundary with [`si`]: the returned `f64` is
//!    the SI base-unit value (uom's `si::f64` storage is base-normalized,
//!    zero-cost). Kernels never see a quantity type.
//!
//! Dimensionless parameters (γ, Pr, mixture fraction Z) stay plain `f64` —
//! a unitless quantity type adds noise, not safety (META-2 §4).
//!
//! The table seam cannot carry compile-time types (units cross as data,
//! FND-5 §3.1) — its gate is `Table::expect_units` (session 9); this crate
//! is the in-process half of the same review-finding-9 closure.

pub use uom::si::f64::{
    DynamicViscosity, HeatTransfer, Length, MassDensity, Pressure, SpecificHeatCapacity,
    TemperatureInterval, ThermalConductivity, ThermodynamicTemperature, Time, Velocity,
    VolumetricHeatCapacity,
};

use uom::si::Quantity;
use uom::si::dynamic_viscosity::pascal_second;
use uom::si::heat_transfer::watt_per_square_meter_kelvin;
use uom::si::length::meter;
use uom::si::mass_density::kilogram_per_cubic_meter;
use uom::si::pressure::pascal;
use uom::si::specific_heat_capacity::joule_per_kilogram_kelvin;
use uom::si::thermal_conductivity::watt_per_meter_kelvin;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::time::second;
use uom::si::velocity::meter_per_second;
use uom::si::volumetric_heat_capacity::joule_per_cubic_meter_kelvin;

/// Extract the SI base-unit `f64` at the kernel boundary (zero-cost).
#[inline]
pub fn si<D, U>(q: Quantity<D, U, f64>) -> f64
where
    D: uom::si::Dimension + ?Sized,
    U: uom::si::Units<f64> + ?Sized,
{
    q.value
}

macro_rules! ctor {
    ($(#[$doc:meta] $name:ident -> $ty:ident, $unit:ident;)*) => {
        $(
            #[$doc]
            #[inline]
            pub fn $name(v: f64) -> $ty {
                $ty::new::<$unit>(v)
            }
        )*
    };
}

ctor! {
    /// Length from meters.
    length_m -> Length, meter;
    /// Temperature from kelvin.
    temperature_k -> ThermodynamicTemperature, kelvin;
    /// Pressure from pascals.
    pressure_pa -> Pressure, pascal;
    /// Velocity from m/s.
    velocity_m_per_s -> Velocity, meter_per_second;
    /// Time from seconds.
    time_s -> Time, second;
    /// Mass density from kg/m³.
    mass_density_kg_per_m3 -> MassDensity, kilogram_per_cubic_meter;
    /// Thermal conductivity from W/(m·K).
    thermal_conductivity_w_per_m_k -> ThermalConductivity, watt_per_meter_kelvin;
    /// Volumetric heat capacity from J/(m³·K).
    volumetric_heat_capacity_j_per_m3_k -> VolumetricHeatCapacity, joule_per_cubic_meter_kelvin;
    /// Convective film coefficient from W/(m²·K).
    heat_transfer_w_per_m2_k -> HeatTransfer, watt_per_square_meter_kelvin;
    /// Dynamic viscosity from Pa·s.
    dynamic_viscosity_pa_s -> DynamicViscosity, pascal_second;
    /// Specific heat capacity from J/(kg·K).
    specific_heat_capacity_j_per_kg_k -> SpecificHeatCapacity, joule_per_kilogram_kelvin;
}

/// Temperature interval (ΔT) from kelvin. Distinct from
/// [`ThermodynamicTemperature`] by uom's `TemperatureKind` — the guard that
/// keeps absolute temperatures out of difference slots (and Celsius bugs
/// out of the codebase).
#[inline]
pub fn temperature_interval_k(v: f64) -> TemperatureInterval {
    TemperatureInterval::new::<uom::si::temperature_interval::kelvin>(v)
}

/// ΔT = hot − cold as the interval type the flux laws multiply by.
#[inline]
pub fn delta_t(
    hot: ThermodynamicTemperature,
    cold: ThermodynamicTemperature,
) -> TemperatureInterval {
    temperature_interval_k(hot.value - cold.value)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn constructors_are_si_base_normalized() {
        assert_eq!(si(length_m(2.5)), 2.5);
        assert_eq!(si(temperature_k(300.0)), 300.0);
        assert_eq!(si(pressure_pa(32.75e5)), 32.75e5);
        assert_eq!(si(thermal_conductivity_w_per_m_k(20.0)), 20.0);
        assert_eq!(si(volumetric_heat_capacity_j_per_m3_k(3.6e6)), 3.6e6);
        assert_eq!(si(heat_transfer_w_per_m2_k(1.5e4)), 1.5e4);
    }

    #[test]
    fn dimensioned_arithmetic_carries_units() {
        // q = h · ΔT must come out as a heat flux [W/m²]: the type system
        // does the bookkeeping this crate exists for. ΔT is the interval
        // kind — `h * temperature_k(...)` does not compile (the guard).
        let h = heat_transfer_w_per_m2_k(1.0e4);
        let dt = delta_t(temperature_k(350.0), temperature_k(300.0));
        let q = h * dt;
        assert_eq!(q.value, 5.0e5);
    }
}
