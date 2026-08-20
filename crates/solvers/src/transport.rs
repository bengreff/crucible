//! FND-7 §3.3 — **the constitutive spine's transport seam**: the ONE
//! provider of the diffusive-flux closure over the medium state, shared by
//! the resolved `F_visc` operator ([`crate::gas_diffusion`]) and the one
//! wall-function heat law ([`crate::wall_heat`]). Neither states a
//! transport constant of its own; both ask here (Rule 13 — one owner).
//!
//! ## Two occupants, one law
//!
//! Which occupant runs is **config data, never a code branch** — the
//! [`crate::euler::EosLaw`] pattern:
//!
//! - [`ConstantTransport`] — the declared-constant occupant (c_p, μ, Pr, γ,
//!   Sc). The degenerate spine of sessions 7–15: a single-composition gas
//!   whose properties do not vary. Still the right occupant for every
//!   analytic fixture, where a constant-coefficient exact solution is the
//!   point.
//! - [`TabulatedTransport`] — the OFFL-5 §3.1a chemical-regime surface:
//!   mixture-averaged Chapman-Enskog transport and CEA caloric derivatives
//!   over the **local state (p, h, Z)**, the same coordinate (and the same
//!   envelope) as the equilibrium EOS surface, so one cell coordinate
//!   serves both and the two refuse together.
//!
//! ## What is tabulated versus what is closed here
//!
//! The surface ships **molecular** transport (μ, `k_frozen`) and the
//! **caloric** derivatives (`cp_frozen`, `cp_equilibrium`,
//! `cv_equilibrium`, `dh_dz`). It deliberately does *not* ship an
//! equilibrium conductivity, because that would be a second model of a flux
//! the operator already carries. On the equilibrium manifold
//! `Y_k(T, p, Z)` the species-enthalpy flux decomposes exactly:
//!
//! ```text
//! Σ_k h_k j_k = −ρD [ (c_p,eq − c_p,fr)·∇T + (∂h/∂Z)|_{p,T}·∇Z ]
//! ```
//!
//! so its temperature limb folds into the Fourier flux as the **effective
//! (equilibrium) conductivity** `k = k_fr + ρD·(c_p,eq − c_p,fr)` — which
//! is the classical equilibrium conductivity — and its composition limb is
//! the resolved species-enthalpy flux the operator adds. **One flux,
//! decomposed by the chain rule, never counted twice.** That composition is
//! done *here*, once, so no consumer re-derives it.
//!
//! **The ∇p limb is neglected — and NOT because it is small.** Measured,
//! `(∂h/∂p)|_{T,Z}·∇p` runs **22–45% of the retained ∇T limb** along an
//! isentrope (it reduces to `−(RT/2ΔH_diss)·dlnp/dlnT`, a structural ~⅓,
//! not a small parameter). It is neglected because of *where* it lives: in
//! a boundary layer `∂p/∂n ≈ 0` to leading order, so the limb vanishes
//! exactly where diffusive flux matters, and streamwise — where ∇p is
//! large — Pe ≫ 1 and the whole diffusive flux is negligible against
//! advection. Stated this way the neglect is falsifiable: it fails in a
//! low-Péclet corner with a streamwise pressure gradient. (This is ordinary
//! Fickian diffusion of the equilibrium composition's pressure response —
//! *not* Stefan-Maxwell baro-diffusion, which is a different term and rides
//! plan S5's species-vector state along with Soret/Dufour.)
//!
//! ## Which c_p — the S4 review's first finding
//!
//! Two different questions wear the same symbol, and conflating them is a
//! **factor 1.7–2.4** on wall heat flux wherever dissociation runs:
//!
//! - `cp` is the **local** slope `∂h/∂T|_p` at this state. It is what
//!   converts a recovery *enthalpy* `r·u²/2` into a recovery *temperature*
//!   at the boundary-layer edge, and it is what pairs with `k` to form the
//!   Prandtl group.
//! - `cp_film` is the **film-mean** slope `(h_aw − h_w)/(T_aw − T_w)`. It
//!   is what a Reynolds/Colburn analogy must drive on, because that
//!   analogy transports *enthalpy*: `q_w = St·ρu·(h_aw − h_w)`, and
//!   `c_p·ΔT` equals `Δh` only for the mean slope over the film.
//!
//! The equilibrium `c_p` is the **peak** of a strongly-peaked curve — at
//! the RL10 chamber it is ~7970 J/(kg·K) against a film mean of ~4130 — so
//! using it as the driving slope overpredicts `q_w` by ~1.9×, one-signed,
//! against a declared ±20–30% band. (The pre-S4 constant `c_p = 5000` was
//! accidentally *inside* that band at 1.21×: S4 made the property more
//! accurate and, until this was fixed, the flux less so.)
//!
//! `cp_film` ships the **frozen** c_p as the declared proxy: measured
//! against the true film mean over the RL10 operating set it is 0.94–1.06
//! at chamber and throat states, degrading to 0.76–0.86 below ~0.2 MPa.
//! That is inside the closure's own band where the wall load lives.
//! **Recorded deferral:** the exact form drives on `h_aw − h_w` directly —
//! the spine is already keyed on `h` and `TableEos` can invert `T(p,h,Z)`
//! at the wall temperature, so this is a root-solve per wall patch, not a
//! new model. It rides the mount-reaction/verdict wave with the
//! skin-friction debit, which touches the same faces.
//!
//! `Pr = μ c_p / k` is likewise **derived, never a second datum**: the
//! constant occupant declares Pr and derives k; the tabulated occupant
//! tabulates k and derives Pr. Both are exact, and a test asserts the two
//! directions agree.
//!
//! ## The Schmidt number
//!
//! `ρD = μ/Sc` is a **declared closure constant**, carried here so it has
//! one owner, not a tabulated property: one effective diffusion coefficient
//! over one composition coordinate is a modeling choice, and it carries its
//! own band (META-3 `schmidt-combustion-gas`, ~0.5 ± 0.3 for H/O combustion
//! gases). Per-species diffusion coefficients become spine outputs when the
//! species-vector state lands (plan S5).
//!
//! ## Declared band
//!
//! **10–20%** on the tabulated occupant's μ and k (FND-7 §3.3/§5): the
//! mixture-averaged-vs-multicomponent spread plus the pure-species
//! Lennard-Jones fits' own. It is an `UncertainInput` declaration consumed
//! by COUP-5, **never folded into a returned value** (FND-5 §2) — as with
//! the wall law's ±20–30%, this module returns scalars.

use crucible_tables::{BoundColumn, Table, TableError};

/// The declared band on the tabulated occupant's transport (FND-7 §3.3):
/// mixture-averaged Chapman-Enskog vs multicomponent, plus the LJ fits'
/// spread at combustion temperatures. Recorded for COUP-5, never applied.
pub const DECLARED_BAND: (f64, f64) = (0.10, 0.20);

/// The HDF5 group's column names and their SI units — the cross-language
/// contract with `offline/crucible_offl/transport.py`'s `_TR_COLUMNS`. The
/// units strings are gated once at bind (`Table::bind`), so a producer-side
/// relabel is a load refusal, never a silent misread.
pub const COLUMNS: [(&str, &str); 6] = [
    ("viscosity", "Pa*s"),
    ("conductivity_frozen", "W/(m*K)"),
    ("cp_frozen", "J/(kg*K)"),
    ("cp_equilibrium", "J/(kg*K)"),
    ("cv_equilibrium", "J/(kg*K)"),
    ("dh_dz", "J/kg"),
];

/// The axis names the surface must carry, in order — the (p, h, Z) local
/// state of FND-7 §3.7 / OFFL-3 §3.3 (S22).
pub const AXES: [&str; 3] = ["p", "h", "Z"];

/// Declared operand rails, shared by BOTH occupants (S4 review): the
/// constant occupant enforces them through its COUP-8 manifest ranges, the
/// tabulated one through [`TabulatedTransport::props`]. Sourced in
/// `crucible_solvers::mechanism`'s `TRANSPORT_CONSTANT_MANIFEST` doc —
/// stated there once, referenced here.
pub const MU_RANGE: (f64, f64) = (1e-7, 1e-2);
/// Heat-capacity rail, applied to c_p and c_v alike (both are specific
/// heats of the same gas; the equilibrium ones run far above the frozen
/// value where dissociation stores energy, which is why the ceiling is the
/// generous 1e5 rather than a combustion-products figure).
pub const CP_RANGE: (f64, f64) = (50.0, 1e5);
/// How far below `k_frozen` the composed effective conductivity may sit
/// before it is refused. Sized to admit generator round-off in the
/// undissociated corners (measured worst case on the shipped surface:
/// 2e-13 relative) and nothing physical.
pub const K_ORDERING_SLACK: f64 = 1e-9;

/// The medium state a spine query is keyed on: the runtime local state.
/// Formed by the EOS occupant that owns the enthalpy coordinate (see
/// [`crate::euler::TableEos::interrogation_php`]) so the transport surface
/// and the EOS surface are always interrogated at the identical triple.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MediumState {
    /// Pressure [Pa].
    pub p: f64,
    /// Specific enthalpy [J/kg] in the EOS surface's own coordinate.
    pub h: f64,
    /// Elemental mixture fraction Z.
    pub z: f64,
}

/// The closed diffusive-flux coefficients at one state — everything both
/// consumers need, with the derivations of the module doc already done.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransportProps {
    /// Dynamic viscosity μ [Pa·s].
    pub mu: f64,
    /// **Effective (equilibrium) thermal conductivity** [W/(m·K)]:
    /// `k_fr + ρD·(c_p,eq − c_p,fr)` (module doc).
    pub k: f64,
    /// Equilibrium specific heat at constant pressure [J/(kg·K)] — the
    /// **local** enthalpy slope `∂h/∂T|_p` at this state.
    pub cp: f64,
    /// **Frozen-composition** specific heat at constant pressure
    /// [J/(kg·K)] — the declared proxy for the **film-mean** slope
    /// `(h_aw − h_w)/(T_aw − T_w)` that a Reynolds-analogy wall law must
    /// drive on. It is NOT a second statement of `cp`: the two answer
    /// different questions, and the difference is a factor ~2 wherever
    /// dissociation runs. See the module doc's "Which c_p" note.
    pub cp_film: f64,
    /// Equilibrium specific heat at constant volume [J/(kg·K)] — the
    /// class-`D` temperature solve's linearization slope `∂e/∂T|_ρ`.
    pub cv: f64,
    /// Species diffusion coefficient ρD = μ/Sc [kg/(m·s)].
    pub rho_d: f64,
    /// `∂h/∂Z|_{p,T}` [J/kg] — the species-enthalpy diffusion-flux
    /// coefficient. Identically zero on a single-composition gas.
    pub dh_dz: f64,
    /// Prandtl number μc_p/k on the EFFECTIVE conductivity — derived,
    /// never a second datum. The transport-consistent group: as the
    /// reactive limb of `k` grows, `pr → Sc` (Le → 1), which is the
    /// physically right asymptote for an equilibrium boundary layer.
    pub pr: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum TransportError {
    /// A declared constant outside the law's domain (META-1 P6).
    BadCoefficient { which: &'static str, value: f64 },
    /// The surface refused the query — off-envelope, off-domain, or a
    /// schema mismatch. Carries the table's own diagnosis.
    OffSurface(String),
    /// A shipped column produced a non-physical value: the table is wrong,
    /// and clamping it would launder that into a plausible run.
    NonPhysical {
        column: &'static str,
        value: f64,
        state: MediumState,
    },
}

impl std::fmt::Display for TransportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::BadCoefficient { which, value } => write!(
                f,
                "transport constant {which} = {value} outside the spine's domain — \
                 refusing, never clamping (FND-7 §3.3, META-1 P6)"
            ),
            Self::OffSurface(e) => write!(f, "spine transport surface refused: {e}"),
            Self::NonPhysical {
                column,
                value,
                state,
            } => write!(
                f,
                "spine transport column {column} = {value} at (p={:.6e}, h={:.6e}, \
                 Z={:.6}) is not physical — the surface is wrong; halting rather \
                 than clamping (META-1 P6)",
                state.p, state.h, state.z
            ),
        }
    }
}

impl std::error::Error for TransportError {}

impl From<TableError> for TransportError {
    fn from(e: TableError) -> Self {
        Self::OffSurface(e.to_string())
    }
}

/// The declared-constant occupant: a single-composition gas whose
/// properties do not vary with state. Its fundamental data is the
/// `(c_p, μ, Pr, γ, Sc)` set the configs have declared since session 7 —
/// `k` follows from Pr, `c_v` from γ, `ρD` from Sc, and `∂h/∂Z` is
/// identically zero because there is nothing for composition to change.
#[derive(Debug, Clone, Copy)]
pub struct ConstantTransport {
    cp: f64,
    mu: f64,
    pr: f64,
    gamma: f64,
    schmidt: f64,
}

impl ConstantTransport {
    /// Unit-typed at the boundary (META-2 §4 ★); Pr, γ, Sc are
    /// dimensionless. Every operand is validated here — the one place
    /// these constants enter the runtime.
    pub fn new(
        cp: crucible_units::SpecificHeatCapacity,
        mu: crucible_units::DynamicViscosity,
        pr: f64,
        gamma: f64,
        schmidt: f64,
    ) -> Result<Self, TransportError> {
        let law = ConstantTransport {
            cp: crucible_units::si(cp),
            mu: crucible_units::si(mu),
            pr,
            gamma,
            schmidt,
        };
        for (which, v) in [
            ("cp", law.cp),
            ("mu", law.mu),
            ("pr", law.pr),
            ("schmidt", law.schmidt),
        ] {
            if !(v.is_finite() && v > 0.0) {
                return Err(TransportError::BadCoefficient { which, value: v });
            }
        }
        // γ ≤ 1 degenerates the gamma-law energy relation (the `flow`
        // manifest's range carries the same bound and the same reason).
        if !(law.gamma.is_finite() && law.gamma > 1.0) {
            return Err(TransportError::BadCoefficient {
                which: "gamma (> 1)",
                value: law.gamma,
            });
        }
        Ok(law)
    }

    /// The closed bundle. Public so a fixture can build the operand a
    /// consumer would otherwise get through [`TransportSpine::eval`] —
    /// same expression, one owner.
    pub fn into_props(&self) -> TransportProps {
        TransportProps {
            mu: self.mu,
            // The session-7..15 `WallLaw::k_gas()` expression, unchanged:
            // a single-composition gas has c_p,eq = c_p,fr, so the module
            // doc's effective-conductivity correction is exactly zero.
            k: self.mu * self.cp / self.pr,
            cp: self.cp,
            // A single-composition gas has one c_p: local and film-mean
            // coincide because there is no composition shift to peak.
            cp_film: self.cp,
            cv: self.cp / self.gamma,
            rho_d: self.mu / self.schmidt,
            dh_dz: 0.0,
            pr: self.pr,
        }
    }
}

/// The tabulated occupant: the OFFL-5 §3.1a surface, bound column-by-column
/// with the units gate, plus the one declared Schmidt number.
#[derive(Debug)]
pub struct TabulatedTransport<'t> {
    mu: BoundColumn<'t>,
    k_frozen: BoundColumn<'t>,
    cp_frozen: BoundColumn<'t>,
    cp_equilibrium: BoundColumn<'t>,
    cv_equilibrium: BoundColumn<'t>,
    dh_dz: BoundColumn<'t>,
    schmidt: f64,
}

impl<'t> TabulatedTransport<'t> {
    /// Bind the six columns and check the axis schema. Fails loud on a
    /// units drift, a missing column, or wrong axes — at bind, once.
    pub fn bind(table: &'t Table, schmidt: f64) -> Result<Self, TransportError> {
        if !(schmidt.is_finite() && schmidt > 0.0) {
            return Err(TransportError::BadCoefficient {
                which: "schmidt",
                value: schmidt,
            });
        }
        let col = |i: usize| -> Result<BoundColumn<'t>, TransportError> {
            let (name, units) = COLUMNS[i];
            Ok(table.bind(name, units)?)
        };
        let bound = Self {
            mu: col(0)?,
            k_frozen: col(1)?,
            cp_frozen: col(2)?,
            cp_equilibrium: col(3)?,
            cv_equilibrium: col(4)?,
            dh_dz: col(5)?,
            schmidt,
        };
        for (i, want) in AXES.iter().enumerate() {
            let got = bound.mu.axis_name(i);
            if got != *want {
                return Err(TransportError::OffSurface(format!(
                    "spine transport surface axis {i} is {got:?}, expected {want:?} — \
                     the runtime local-state coordinate is (p, h, Z) (FND-7 §3.7)"
                )));
            }
        }
        Ok(bound)
    }

    /// The declared validity envelope per axis, in [`AXES`] order — the
    /// engine cross-checks it against the EOS surface's at assembly.
    pub fn envelopes(&self) -> [(f64, f64); 3] {
        [
            self.mu.axis_envelope(0),
            self.mu.axis_envelope(1),
            self.mu.axis_envelope(2),
        ]
    }

    /// The measured interpolation bounds (FND-5 §3.4) of each column, in
    /// [`COLUMNS`] order — declared metadata for COUP-5's epistemic
    /// interval, never folded into a value.
    pub fn interp_error_bounds(&self) -> [f64; 6] {
        [
            self.mu.interp_error_bound(),
            self.k_frozen.interp_error_bound(),
            self.cp_frozen.interp_error_bound(),
            self.cp_equilibrium.interp_error_bound(),
            self.cv_equilibrium.interp_error_bound(),
            self.dh_dz.interp_error_bound(),
        ]
    }

    fn props(&self, st: MediumState) -> Result<TransportProps, TransportError> {
        let q = [st.p, st.h, st.z];
        let mu = self.mu.interpolate(&q)?;
        let k_fr = self.k_frozen.interpolate(&q)?;
        let cp_fr = self.cp_frozen.interpolate(&q)?;
        let cp_eq = self.cp_equilibrium.interpolate(&q)?;
        let cv_eq = self.cv_equilibrium.interpolate(&q)?;
        let dh_dz = self.dh_dz.interpolate(&q)?;
        for (column, v, positive) in [
            ("viscosity", mu, true),
            ("conductivity_frozen", k_fr, true),
            ("cp_frozen", cp_fr, true),
            ("cp_equilibrium", cp_eq, true),
            ("cv_equilibrium", cv_eq, true),
            ("dh_dz", dh_dz, false),
        ] {
            if !v.is_finite() || (positive && v <= 0.0) {
                return Err(TransportError::NonPhysical {
                    column,
                    value: v,
                    state: st,
                });
            }
        }
        let rho_d = mu / self.schmidt;
        // The effective (equilibrium) conductivity — the module doc's
        // chain-rule limb. Composition shift can only ADD enthalpy
        // response, so physically `cp_eq >= cp_fr` and this can only raise
        // k. That ordering is NOT guaranteed by the positivity loop above,
        // which tests each column alone (S4 review): 2.3% of the shipped
        // surface's nodes carry `cp_eq - cp_fr` slightly negative — CEA
        // round-off, magnitude ≤ 7e-10 J/(kg·K), in the cold corners where
        // nothing is dissociated and there is no reactive conduction to
        // find. Those are the data as generated and are left alone
        // (doctoring a table to satisfy an inequality is worse than the
        // 1e-13 W/(m·K) it would "fix"). What guards the real failure is
        // the explicit floor below: a surface whose ordering were wrong by
        // any physically meaningful amount would drive k below k_fr, and
        // that refuses.
        let k = k_fr + rho_d * (cp_eq - cp_fr);
        if !(k.is_finite() && k > 0.0) {
            return Err(TransportError::NonPhysical {
                column: "conductivity_effective",
                value: k,
                state: st,
            });
        }
        // The reactive limb may not SUBTRACT conduction. `K_ORDERING_SLACK`
        // admits the round-off band above and nothing wider.
        if k < k_fr * (1.0 - K_ORDERING_SLACK) {
            return Err(TransportError::NonPhysical {
                column: "conductivity_effective (below the frozen conductivity)",
                value: k,
                state: st,
            });
        }
        // The tabulated occupant gets the SAME declared rails the constant
        // occupant's manifest carries (S4 review): a table is not more
        // trustworthy than a config value, and the reasons those bounds
        // exist — a metal-like μ gives the gas a metal-like k, a Pr outside
        // [0.05, 5] is not a gas — apply to whatever produced the number.
        for (column, v, lo, hi) in [
            ("viscosity", mu, MU_RANGE.0, MU_RANGE.1),
            ("cp_equilibrium", cp_eq, CP_RANGE.0, CP_RANGE.1),
            ("cp_frozen", cp_fr, CP_RANGE.0, CP_RANGE.1),
            ("cv_equilibrium", cv_eq, CP_RANGE.0, CP_RANGE.1),
        ] {
            if v < lo || v > hi {
                return Err(TransportError::NonPhysical {
                    column,
                    value: v,
                    state: st,
                });
            }
        }
        Ok(TransportProps {
            mu,
            k,
            cp: cp_eq,
            cp_film: cp_fr,
            cv: cv_eq,
            rho_d,
            dh_dz,
            pr: mu * cp_eq / k,
        })
    }
}

/// The spine's transport slot. One `eval`, two occupants, no branch at any
/// consumer (module doc).
///
/// The occupants differ in size (six bound columns vs five scalars) and
/// that is fine: exactly one of these is constructed per run, at assembly,
/// and it is then borrowed — never moved, never collected. Boxing the large
/// variant would buy nothing but a pointer chase on the hot query.
#[allow(clippy::large_enum_variant)]
#[derive(Debug)]
pub enum TransportSpine<'t> {
    Constant(ConstantTransport),
    Tabulated(TabulatedTransport<'t>),
}

impl TransportSpine<'_> {
    /// The one query: closed diffusive-flux coefficients at a medium state.
    pub fn eval(&self, st: MediumState) -> Result<TransportProps, TransportError> {
        match self {
            // The constant occupant is state-independent by construction —
            // that IS its physics, not a shortcut (a single-composition gas
            // whose properties do not vary).
            Self::Constant(c) => Ok(c.into_props()),
            Self::Tabulated(t) => t.props(st),
        }
    }

    /// Whether this occupant's properties vary with state — the only thing
    /// a consumer may ask *about* the occupant, and it asks so it can size
    /// a diagnostic, never to pick a code path.
    pub fn is_tabulated(&self) -> bool {
        matches!(self, Self::Tabulated(_))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crucible_units::{dynamic_viscosity_pa_s, specific_heat_capacity_j_per_kg_k};

    const CP: f64 = 5000.0;
    const MU: f64 = 1.0e-4;
    const PR: f64 = 0.6;
    const GAMMA: f64 = 1.2;
    const SC: f64 = 0.5;

    fn constant() -> TransportSpine<'static> {
        TransportSpine::Constant(
            ConstantTransport::new(
                specific_heat_capacity_j_per_kg_k(CP),
                dynamic_viscosity_pa_s(MU),
                PR,
                GAMMA,
                SC,
            )
            .expect("valid constant transport set"),
        )
    }

    /// The bit-identity contract with sessions 7–15: `k = μc_p/Pr` is the
    /// retired `WallLaw::k_gas()` expression, `c_v = c_p/γ` and
    /// `ρD = μ/Sc` are the retired `GasDiffusion::from_transport` ones.
    /// Certificates regenerate byte-identically only if these hold to the
    /// last bit.
    #[test]
    fn constant_occupant_reproduces_the_retired_expressions_bitwise() {
        let s = constant();
        let st = MediumState {
            p: 3.0e6,
            h: -1.0e6,
            z: 1.0 / 6.0,
        };
        let t = s.eval(st).unwrap();
        assert_eq!(t.mu, MU);
        assert_eq!(t.k, MU * CP / PR);
        assert_eq!(t.cp, CP);
        assert_eq!(t.cv, CP / GAMMA);
        assert_eq!(t.rho_d, MU / SC);
        assert_eq!(t.pr, PR);
        assert_eq!(t.dh_dz, 0.0);
    }

    #[test]
    fn constant_occupant_is_state_independent() {
        let s = constant();
        let a = s
            .eval(MediumState {
                p: 1.0e2,
                h: -1.2e7,
                z: 0.155,
            })
            .unwrap();
        let b = s
            .eval(MediumState {
                p: 7.0e6,
                h: 3.0e6,
                z: 0.185,
            })
            .unwrap();
        assert_eq!(a, b);
    }

    /// Pr is a derived group either way round: the constant occupant
    /// declares Pr and derives k, so μc_p/k must return the declared Pr.
    #[test]
    fn prandtl_is_consistent_with_the_derived_conductivity() {
        let s = constant();
        let t = s
            .eval(MediumState {
                p: 1.0,
                h: 0.0,
                z: 0.0,
            })
            .unwrap();
        let round_trip = t.mu * t.cp / t.k;
        assert!(
            (round_trip - t.pr).abs() <= 8.0 * f64::EPSILON * t.pr,
            "Pr round-trip {round_trip} vs declared {}",
            t.pr
        );
    }

    #[test]
    fn bad_constants_refuse_rather_than_clamp() {
        let ok = |cp, mu, pr, gamma, sc| {
            ConstantTransport::new(
                specific_heat_capacity_j_per_kg_k(cp),
                dynamic_viscosity_pa_s(mu),
                pr,
                gamma,
                sc,
            )
        };
        assert!(ok(CP, MU, PR, GAMMA, SC).is_ok());
        assert!(ok(0.0, MU, PR, GAMMA, SC).is_err());
        assert!(ok(CP, -MU, PR, GAMMA, SC).is_err());
        assert!(ok(CP, MU, f64::NAN, GAMMA, SC).is_err());
        // γ = 1 is the singular isothermal limit, not a valid gas.
        assert!(ok(CP, MU, PR, 1.0, SC).is_err());
        assert!(ok(CP, MU, PR, GAMMA, 0.0).is_err());
    }
}
