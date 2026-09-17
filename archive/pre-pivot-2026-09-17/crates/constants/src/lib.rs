//! The FND-1 (v0.x §3.2) typed constants module, implementing the META-3 §4
//! physical-constants table. Single pinned source: **CODATA 2022 recommended
//! values** (NIST CUU, physics.nist.gov), retrieved 2026-08-14.
//!
//! Rules (META-2 §2.2, FND-1 §3.2):
//! - SI base units, `f64`, values at full source precision — never re-typed
//!   per solver; cite by META-3 key.
//! - Every value carries its provenance: exact-by-definition constants are
//!   flagged as such; measured values carry their CODATA 1σ standard
//!   uncertainty. The uncertainty is *metadata for the ledger and UQ
//!   declarations* (FND-1) — runtime arithmetic uses the scalar value.
//! - `REGISTRY` iterates in the fixed META-3 §4 table order (deterministic
//!   output ordering, META-1 §2 / META-2 §4 ★).
//!
//! Units typing (review finding 9) is closed at the boundaries by
//! `crucible-units` (session 10); this crate deliberately remains the
//! documented-SI `f64` kernel layer — a constant gains a typed wrapper at
//! the boundary that consumes it, when one appears.

/// Where a constant's value comes from, and what uncertainty it carries.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Provenance {
    /// Exact by definition (2019 SI redefinition, or a defined convention
    /// such as standard gravity). Zero uncertainty.
    Exact,
    /// CODATA 2022 measured value; `sigma` is the 1σ standard uncertainty in
    /// the same SI units as the value.
    Measured { sigma: f64 },
}

/// One ledger entry: a constant plus the metadata that makes it citable
/// (META-1 Principle 7 — an input a paper cannot cite is an input we do not use).
#[derive(Debug, Clone, Copy)]
pub struct Constant {
    /// META-3 §4 table key (`g0`, `c`, `kB`, …) — the citation handle.
    pub key: &'static str,
    pub symbol: &'static str,
    /// Value in SI base/coherent units, full CODATA source precision.
    pub value: f64,
    pub provenance: Provenance,
    pub units: &'static str,
    pub note: &'static str,
}

/// The pinned source, stamped: cite this string in provenance output.
pub const SOURCE: &str =
    "CODATA 2022 recommended values (NIST CUU), retrieved 2026-08-14; META-3 \u{a7}4";

/// Standard gravity g₀ [m/s²] — exact by definition (Isp convention). META-3:`g0`.
pub const G0: f64 = 9.806_65;
/// Speed of light in vacuum c [m/s] — exact by definition. META-3:`c`.
pub const C: f64 = 299_792_458.0;
/// Boltzmann constant k_B [J/K] — exact (2019 SI). META-3:`kB`.
pub const K_B: f64 = 1.380_649e-23;
/// Avogadro constant N_A [1/mol] — exact (2019 SI). META-3:`NA`.
pub const N_A: f64 = 6.022_140_76e23;
/// Elementary charge e [C] — exact (2019 SI). META-3:`e`.
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Electron mass m_e [kg] — CODATA 2022, 1σ = 2.8e-40 kg. META-3:`me`.
pub const M_E: f64 = 9.109_383_713_9e-31;
/// Proton mass m_p [kg] — CODATA 2022, 1σ = 5.2e-37 kg. META-3:`mp`.
pub const M_P: f64 = 1.672_621_925_95e-27;
/// Neutron mass m_n [kg] — CODATA 2022, 1σ = 8.5e-37 kg. META-3:`mn`.
pub const M_N: f64 = 1.674_927_500_56e-27;
/// Atomic mass constant m_u [kg] — CODATA 2022, 1σ = 5.2e-37 kg. META-3:`mu`.
pub const M_U: f64 = 1.660_539_068_92e-27;
/// Electron-volt → joule conversion [J] — exact (numerically = e). META-3:`eV`.
pub const EV: f64 = 1.602_176_634e-19;
/// Barn → m² conversion [m²] — exact by definition. META-3:`barn`.
pub const BARN: f64 = 1e-28;

/// The full ledger, in META-3 §4 table order (fixed — deterministic iteration).
pub const REGISTRY: &[Constant] = &[
    Constant {
        key: "g0",
        symbol: "g\u{2080}",
        value: G0,
        provenance: Provenance::Exact,
        units: "m/s^2",
        note: "standard gravity, exact by definition (Isp convention)",
    },
    Constant {
        key: "c",
        symbol: "c",
        value: C,
        provenance: Provenance::Exact,
        units: "m/s",
        note: "speed of light in vacuum, exact by definition",
    },
    Constant {
        key: "kB",
        symbol: "k_B",
        value: K_B,
        provenance: Provenance::Exact,
        units: "J/K",
        note: "Boltzmann constant, exact (2019 SI redefinition)",
    },
    Constant {
        key: "NA",
        symbol: "N_A",
        value: N_A,
        provenance: Provenance::Exact,
        units: "1/mol",
        note: "Avogadro constant, exact (2019 SI redefinition)",
    },
    Constant {
        key: "e",
        symbol: "e",
        value: E_CHARGE,
        provenance: Provenance::Exact,
        units: "C",
        note: "elementary charge, exact (2019 SI redefinition)",
    },
    Constant {
        key: "me",
        symbol: "m_e",
        value: M_E,
        provenance: Provenance::Measured { sigma: 2.8e-40 },
        units: "kg",
        note: "electron mass, CODATA 2022",
    },
    Constant {
        key: "mp",
        symbol: "m_p",
        value: M_P,
        provenance: Provenance::Measured { sigma: 5.2e-37 },
        units: "kg",
        note: "proton mass, CODATA 2022",
    },
    Constant {
        key: "mn",
        symbol: "m_n",
        value: M_N,
        provenance: Provenance::Measured { sigma: 8.5e-37 },
        units: "kg",
        note: "neutron mass, CODATA 2022",
    },
    Constant {
        key: "mu",
        symbol: "m_u",
        value: M_U,
        provenance: Provenance::Measured { sigma: 5.2e-37 },
        units: "kg",
        note: "atomic mass constant, CODATA 2022",
    },
    Constant {
        key: "eV",
        symbol: "eV",
        value: EV,
        provenance: Provenance::Exact,
        units: "J",
        note: "electron-volt to joule, exact (numerically equal to e)",
    },
    Constant {
        key: "barn",
        symbol: "b",
        value: BARN,
        provenance: Provenance::Exact,
        units: "m^2",
        note: "barn to square metre, exact by definition",
    },
];

/// Look up a constant by its META-3 §4 key (linear scan over the fixed-order
/// slice — no hash iteration in output paths, META-2 §4 ★).
pub fn by_key(key: &str) -> Option<&'static Constant> {
    REGISTRY.iter().find(|c| c.key == key)
}

// Tests cover META-3 §4 (constants table, CODATA 2022 stamp) and FND-1 §3.2
// (typed constants module) — named to the doc IDs per META-2 §4 / VAL-3 §3.1.
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn meta3_s4_exact_constants_match_their_definitions() {
        // 2019 SI defining constants + defined conventions, bit-exact.
        assert_eq!(G0, 9.806_65);
        assert_eq!(C, 299_792_458.0);
        assert_eq!(K_B, 1.380_649e-23);
        assert_eq!(N_A, 6.022_140_76e23);
        assert_eq!(E_CHARGE, 1.602_176_634e-19);
        assert_eq!(EV, 1.602_176_634e-19);
        assert_eq!(BARN, 1e-28);
    }

    #[test]
    fn meta3_s4_ev_is_numerically_the_elementary_charge() {
        assert_eq!(EV, E_CHARGE);
    }

    #[test]
    fn meta3_s4_registry_is_complete_ordered_and_unique() {
        let expected_order = [
            "g0", "c", "kB", "NA", "e", "me", "mp", "mn", "mu", "eV", "barn",
        ];
        let keys: Vec<&str> = REGISTRY.iter().map(|c| c.key).collect();
        assert_eq!(
            keys, expected_order,
            "META-3 \u{a7}4 table order is the contract"
        );
        for c in REGISTRY {
            assert!(!c.units.is_empty(), "{}: bare value with no unit", c.key);
            assert!(!c.note.is_empty(), "{}: uncited value", c.key);
            assert!(
                c.value.is_finite() && c.value > 0.0,
                "{}: non-physical",
                c.key
            );
        }
    }

    #[test]
    fn meta3_s4_registry_agrees_with_named_constants() {
        assert_eq!(by_key("g0").unwrap().value, G0);
        assert_eq!(by_key("c").unwrap().value, C);
        assert_eq!(by_key("kB").unwrap().value, K_B);
        assert_eq!(by_key("NA").unwrap().value, N_A);
        assert_eq!(by_key("e").unwrap().value, E_CHARGE);
        assert_eq!(by_key("me").unwrap().value, M_E);
        assert_eq!(by_key("mp").unwrap().value, M_P);
        assert_eq!(by_key("mn").unwrap().value, M_N);
        assert_eq!(by_key("mu").unwrap().value, M_U);
        assert_eq!(by_key("eV").unwrap().value, EV);
        assert_eq!(by_key("barn").unwrap().value, BARN);
        assert!(by_key("nonexistent").is_none());
    }

    #[test]
    fn fnd1_s32_measured_masses_are_sane() {
        // Compile-time facts: physical ordering m_u (C-12/12) < m_p < m_n, and
        // the n–p mass difference ≈ 2.3e-30 kg (≈ 1.29 MeV/c²), loose bracket.
        const {
            assert!(M_U < M_P && M_P < M_N);
            let dm = M_N - M_P;
            assert!(dm > 2.2e-30 && dm < 2.4e-30);
        }
        // CODATA 2022 relative standard uncertainties are all ≤ ~5.1e-10.
        for c in REGISTRY {
            if let Provenance::Measured { sigma } = c.provenance {
                let rel = sigma / c.value;
                assert!(rel > 0.0 && rel < 1e-9, "{}: implausible \u{3c3}", c.key);
            }
        }
    }
}
