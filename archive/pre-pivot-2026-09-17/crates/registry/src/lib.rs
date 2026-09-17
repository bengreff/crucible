//! Implements the loader-facing subset of **COUP-8 v0.4 §3.1–§3.4**: the
//! static, const-evaluable mechanism `Manifest` and the explicit registry
//! table (`id → Manifest`) that the FND-4 loader dispatches config blocks
//! against — validation before construction, the MOOSE `validParams()`
//! pattern as data.
//!
//! Not yet carried (arrives with its consumer, per the skeleton split):
//! `grid_fields` (needs FND-2's conserved-state component set), `halts[]`
//! (COUP-4), table validity envelopes on `TableReq` (FND-5 session), and the
//! mechanism *constructor* column (first real operator, the conduction
//! session). Closed enums below extend **by doc amendment only** (their
//! vocabularies are owned by VISION_SCOPE §7.4/§7.5/§7.6 and COUP-7/COUP-8).
//!
//! Determinism (COUP-8 §3.2/§3.5): the registry is one explicit table sorted
//! by `id` — construction *asserts* sorted-unique order, lookups are binary
//! search, iteration is slice order. No distributed registration, no hash
//! maps, no link-order dependence.

/// COUP-8 §3.4 — the interface contract SemVer a mechanism was written against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct InterfaceVersion {
    pub major: u16,
    pub minor: u16,
}

/// The core's current interface major, and the minor range it supports.
/// Mechanisms migrate independently within the minor range (COUP-8 §3.4).
pub const CORE_INTERFACE_MAJOR: u16 = 1;
pub const CORE_SUPPORTED_MINOR_MIN: u16 = 0;
pub const CORE_SUPPORTED_MINOR_MAX: u16 = 0;

impl InterfaceVersion {
    /// The §3.3 check-4 predicate: does the core support this version?
    pub fn supported_by_core(self) -> bool {
        self.major == CORE_INTERFACE_MAJOR
            && (CORE_SUPPORTED_MINOR_MIN..=CORE_SUPPORTED_MINOR_MAX).contains(&self.minor)
    }
}

/// Typed parameter schema entry (COUP-8 §3.1 `params`).
///
/// `required` is encoded structurally: `default: None` ⇔ required. An
/// optional parameter *must* carry a default — FND-4 §3.5 puts every default
/// in the registry descriptor, so "optional without default" would be a
/// hidden default and is unrepresentable here.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParamSpec {
    pub name: &'static str,
    pub ty: ParamType,
    pub default: Option<ParamValue>,
    /// Inclusive validity range for numeric parameters (SI units documented
    /// in the parameter name per META-2 §2.1).
    pub range: Option<(f64, f64)>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParamType {
    Bool,
    Int,
    Float,
    Str,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ParamValue {
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(&'static str),
}

/// Boundary-object port declaration (COUP-8 §3.1 `ports`, COUP-7 objects).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PortSpec {
    pub name: &'static str,
    pub role: PortRole,
    pub kind: PortKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortRole {
    Require,
    Provide,
}

/// Closed port-kind vocabulary — the COUP-7 v1 boundary-object set
/// (VISION_SCOPE §7.5). Extends by doc amendment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PortKind {
    Pump,
    PowerSupply,
    Beam,
    CoilSet,
    ReactivitySchedule,
    AntiprotonDelivery,
    Injector,
    Radiator,
}

/// Coupler participation (COUP-8 §3.1 `couplers`); the 7 kinds are the
/// exhaustive VISION_SCOPE §7.4 list — a new kind requires a doc amendment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CouplerSpec {
    pub kind: CouplerKind,
    pub quantity: Quantity,
    pub direction: Direction,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CouplerKind {
    Deposition,
    WallExchange,
    Recession,
    FieldSampling,
    KineticsFeedback,
    EventInjection,
    PortAccounting,
}

/// Closed conserved/exchanged-quantity vocabulary for coupler edges.
/// Seed set; extends by doc amendment as couplers land (COUP-2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Quantity {
    Mass,
    Momentum,
    Energy,
    Power,
    Heat,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Source,
    Sink,
    Bidirectional,
}

/// Operating-regime vocabulary for the per-regime chaotic classification
/// (COUP-8 §3.1, O21). Derived from VISION_SCOPE §7.6's execution modes;
/// refined when the operating-profile grammar lands (FND-4 deferred block).
/// Until then the FND-4 loader is conservative: ANY `Chaotic` entry on a
/// selected mechanism blocks `relaxed` determinism mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Regime {
    Steady,
    Transient,
    Pulsed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChaoticClass {
    Chaotic,
    NonChaotic,
}

/// Required-table declaration (COUP-8 §3.1 `tables`). The consumer envelope
/// field arrives with the FND-5 loader session (envelope-coverage check §3.3-2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TableReq {
    pub semantic_name: &'static str,
    pub required: bool,
}

/// COUP-8 §3.1 — the static capability declaration, validated before any
/// construction. One per mechanism (and per FND-7 material type).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Manifest {
    /// Stable registry key; config selects by this string.
    pub id: &'static str,
    pub interface_version: InterfaceVersion,
    pub tables: &'static [TableReq],
    pub couplers: &'static [CouplerSpec],
    pub ports: &'static [PortSpec],
    /// O21: explicit per-regime classification — never by omission.
    pub chaotic_class: &'static [(Regime, ChaoticClass)],
    pub params: &'static [ParamSpec],
}

impl Manifest {
    /// Conservative O21 predicate used by FND-4 until the operating-profile
    /// grammar can name the configured regime.
    pub fn any_chaotic(&self) -> bool {
        self.chaotic_class
            .iter()
            .any(|(_, c)| *c == ChaoticClass::Chaotic)
    }
}

/// The explicit registry table (COUP-8 §3.2): mechanisms and material types,
/// each a slice sorted by `id`. Construction asserts the ordering invariants
/// so lookup order is canonical, not incidental.
#[derive(Debug, Clone, Copy)]
pub struct Registry {
    mechanisms: &'static [&'static Manifest],
    materials: &'static [&'static Manifest],
}

impl Registry {
    /// Build a registry, asserting COUP-8 §3.2/§3.5 invariants: ids sorted
    /// and unique; chaotic classification explicitly declared (never by
    /// omission); no param named `type` (reserved by the FND-4 block
    /// grammar); no `.` in ids or port names (reserved by the O20
    /// `instance.port` binding grammar).
    ///
    /// # Panics
    /// On any violated invariant — a malformed registry is a build defect,
    /// not a runtime condition.
    pub fn new(
        mechanisms: &'static [&'static Manifest],
        materials: &'static [&'static Manifest],
    ) -> Registry {
        for set in [mechanisms, materials] {
            for pair in set.windows(2) {
                assert!(
                    pair[0].id < pair[1].id,
                    "registry ids must be sorted and unique: {:?} !< {:?}",
                    pair[0].id,
                    pair[1].id
                );
            }
            for m in set {
                assert!(
                    !m.chaotic_class.is_empty(),
                    "{}: chaotic_class must be declared explicitly, never by omission (O21)",
                    m.id
                );
                assert!(
                    !m.id.contains('.'),
                    "{}: `.` reserved by O20 binding grammar",
                    m.id
                );
                for p in m.ports {
                    assert!(
                        !p.name.contains('.'),
                        "{}.{}: `.` reserved by O20 binding grammar",
                        m.id,
                        p.name
                    );
                }
                for p in m.params {
                    assert!(
                        p.name != "type",
                        "{}: param name `type` is reserved by the FND-4 block grammar",
                        m.id
                    );
                }
            }
        }
        Registry {
            mechanisms,
            materials,
        }
    }

    pub fn mechanism(&self, id: &str) -> Option<&'static Manifest> {
        self.mechanisms
            .binary_search_by(|m| m.id.cmp(id))
            .ok()
            .map(|i| self.mechanisms[i])
    }

    pub fn material(&self, id: &str) -> Option<&'static Manifest> {
        self.materials
            .binary_search_by(|m| m.id.cmp(id))
            .ok()
            .map(|i| self.materials[i])
    }

    /// Valid mechanism ids in canonical (sorted) order — for fail-loud
    /// "unknown mechanism, valid ids are …" diagnostics.
    pub fn mechanism_ids(&self) -> impl Iterator<Item = &'static str> + '_ {
        self.mechanisms.iter().map(|m| m.id)
    }

    pub fn material_ids(&self) -> impl Iterator<Item = &'static str> + '_ {
        self.materials.iter().map(|m| m.id)
    }
}

// Covers COUP-8 §3.2 (registry invariants) — named to doc ID per META-2 §4.
#[cfg(test)]
mod tests {
    use super::*;

    static M_A: Manifest = Manifest {
        id: "a-mech",
        interface_version: InterfaceVersion { major: 1, minor: 0 },
        tables: &[],
        couplers: &[],
        ports: &[],
        chaotic_class: &[(Regime::Steady, ChaoticClass::NonChaotic)],
        params: &[],
    };
    static M_B: Manifest = Manifest {
        id: "b-mech",
        interface_version: InterfaceVersion { major: 1, minor: 0 },
        tables: &[],
        couplers: &[],
        ports: &[],
        chaotic_class: &[(Regime::Transient, ChaoticClass::Chaotic)],
        params: &[],
    };

    static SORTED: [&Manifest; 2] = [&M_A, &M_B];
    static UNSORTED: [&Manifest; 2] = [&M_B, &M_A];

    #[test]
    fn coup8_s32_lookup_and_canonical_order() {
        let r = Registry::new(&SORTED, &[]);
        assert_eq!(r.mechanism("a-mech").unwrap().id, "a-mech");
        assert!(r.mechanism("zz").is_none());
        let ids: Vec<_> = r.mechanism_ids().collect();
        assert_eq!(ids, ["a-mech", "b-mech"]);
    }

    #[test]
    #[should_panic(expected = "sorted")]
    fn coup8_s32_unsorted_registry_is_a_build_defect() {
        let _ = Registry::new(&UNSORTED, &[]);
    }

    #[test]
    fn coup8_s31_chaotic_predicate() {
        assert!(!M_A.any_chaotic());
        assert!(M_B.any_chaotic());
    }

    #[test]
    fn coup8_s34_interface_version_gate() {
        assert!(InterfaceVersion { major: 1, minor: 0 }.supported_by_core());
        assert!(!InterfaceVersion { major: 2, minor: 0 }.supported_by_core());
        assert!(!InterfaceVersion { major: 1, minor: 9 }.supported_by_core());
    }
}
