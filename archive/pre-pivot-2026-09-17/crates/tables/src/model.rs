//! FND-5 §3.1 — the in-memory table model after a validated load, and the
//! closed error set. Everything here is owned `f64` data: the HDF5 handle is
//! dropped at the end of `Table::open`, so queries touch no C library.

use std::fmt;

/// Table-format schema major version (§3.1 `schema_version`); the loader
/// refuses any other major (§3.2).
pub const TABLE_SCHEMA_MAJOR: i64 = 1;

/// §3.1 mandatory provenance attributes — the set that flows into the run
/// manifest (FND-4/FND-6) and walks back to META-3 keys.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Provenance {
    pub producer: String,
    pub producer_version: String,
    pub input_deck_hash: String,
    pub source_library: String,
    pub generator_commit: String,
    /// Present iff the generator was stochastic.
    pub rng_seed: Option<i64>,
}

/// One coordinate axis: explicit grid points (never an implied linspace) and
/// the declared validity envelope, which may be tighter than the grid.
#[derive(Debug, Clone, PartialEq)]
pub struct Axis {
    pub name: String,
    /// Strictly increasing grid coordinates.
    pub points: Vec<f64>,
    pub envelope_min: f64,
    pub envelope_max: f64,
}

/// One value dataset (quantity) with its per-column semantics (§3.1).
#[derive(Debug, Clone, PartialEq)]
pub struct TableValue {
    pub name: String,
    /// Row-major over the axes in table order.
    pub data: Vec<f64>,
    pub units: String,
    /// `lin`/`log` tokens, one per axis then one for the value, `-`-separated
    /// (e.g. `"log-log"` for a 1-axis log-log table) — the concrete grammar
    /// realizing §3.1's `interp_rule`.
    pub interp_rule: String,
    /// §3.4: producer-measured bound, table metadata for COUP-5 — never
    /// folded into a runtime return value.
    pub interp_error_bound: f64,
    /// Session-12 review: the same producer-measured bound in the VALUE'S
    /// RULE SPACE (|Δ ln value|, ≈ relative) for log-valued columns — the
    /// scale-honest form an acceptance gate must use where an absolute
    /// bound is attained at the dense end of a log-valued range. `None`
    /// for linear-valued columns and pre-0.3.2 artifacts (consumers fall
    /// back to the absolute bound — old-artifact behavior preserved).
    /// RECORDED DEFERRAL: rides OUTSIDE digest v3 (the pinned
    /// cross-language field set); digest v4 folds it in.
    pub interp_error_bound_log: Option<f64>,
    /// Sibling `sigma_<name>` dataset (per-point 1σ), if the producer
    /// supplied one; a scalar band may come as `sigma_<name>` attribute.
    pub sigma: Option<Vec<f64>>,
    pub sigma_scalar: Option<f64>,
}

/// §3.5 out-of-envelope policy: `Refuse` (default for headline runs) errors;
/// `Flag` computes the value **iff still inside the grid domain** and
/// reports the excursion to the caller (COUP-5 records per-member hits).
/// Outside the grid domain is always a hard error — the runtime never
/// extrapolates, under either policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum EnvelopePolicy {
    #[default]
    Refuse,
    Flag,
}

/// A fully-validated, immutable table.
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    pub group_path: String,
    pub kind: String,
    pub schema_version: String,
    pub data_version: String,
    pub interp_method: String,
    pub provenance: Provenance,
    pub axes: Vec<Axis>,
    pub values: Vec<TableValue>,
    /// `sha256:<hex>` canonical content digest (see `digest` module) —
    /// verified against the pin at load, recorded into the manifest.
    pub content_digest: String,
}

/// One envelope excursion under `EnvelopePolicy::Flag`.
#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopeHit {
    pub axis: String,
    pub query: f64,
    pub envelope_min: f64,
    pub envelope_max: f64,
}

/// The closed load/query error set — every variant is a halt-with-diagnosis,
/// never a guess (META-1 P6).
#[derive(Debug, Clone, PartialEq)]
pub enum TableError {
    Hdf5(String),
    MissingAttribute {
        path: String,
        attribute: String,
    },
    MissingDataset {
        path: String,
    },
    SchemaMajorMismatch {
        found: String,
        supported: i64,
    },
    PinVersionMismatch {
        pinned: String,
        found: String,
    },
    PinDigestMismatch {
        pinned: String,
        found: String,
    },
    NonMonotonicAxis {
        axis: String,
    },
    EnvelopeOutsideGrid {
        axis: String,
    },
    ShapeMismatch {
        value: String,
        expected: usize,
        found: usize,
    },
    BadInterpRule {
        value: String,
        rule: String,
        reason: String,
    },
    LogScaleRequiresPositive {
        value: String,
        axis_or_value: String,
    },
    UnknownValue {
        name: String,
    },
    /// META-2 §4 ★ boundary units check: the consumer's expected units
    /// string does not match the column's declared `units` metadata.
    UnitsMismatch {
        value: String,
        expected: String,
        found: String,
    },
    QueryArity {
        expected: usize,
        found: usize,
    },
    OutOfEnvelope {
        axis: String,
        query: f64,
        envelope_min: f64,
        envelope_max: f64,
    },
    OutOfDomain {
        axis: String,
        query: f64,
        grid_min: f64,
        grid_max: f64,
    },
    NotYetImplemented {
        what: String,
        arrives_with: String,
    },
    UnknownKind {
        kind: String,
    },
    UnknownInterpMethod {
        method: String,
    },
    /// A dataset/group the schema does not account for: unlisted axis
    /// dataset, unknown group member, or a `sigma_<x>` with no value `<x>`.
    /// Refused so the content digest covers every byte the file carries
    /// (same-label/different-bytes must fail the pin, §3.2).
    UnexpectedMember {
        path: String,
        reason: String,
    },
}

impl fmt::Display for TableError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Hdf5(m) => write!(f, "HDF5 error: {m}"),
            Self::MissingAttribute { path, attribute } => {
                write!(
                    f,
                    "{path}: mandatory attribute {attribute:?} missing — table cannot load (FND-5 §3.1)"
                )
            }
            Self::MissingDataset { path } => write!(f, "{path}: dataset missing"),
            Self::SchemaMajorMismatch { found, supported } => {
                write!(
                    f,
                    "table schema_version {found:?} major-incompatible with supported major {supported} (§3.2)"
                )
            }
            Self::PinVersionMismatch { pinned, found } => {
                write!(
                    f,
                    "data_version {found:?} does not match the config pin {pinned:?} (§3.2, S6)"
                )
            }
            Self::PinDigestMismatch { pinned, found } => {
                write!(
                    f,
                    "content digest {found} does not match the pin {pinned} — same-label/different-bytes refused (§3.2)"
                )
            }
            Self::NonMonotonicAxis { axis } => {
                write!(f, "axis {axis:?} is not strictly increasing")
            }
            Self::EnvelopeOutsideGrid { axis } => {
                write!(
                    f,
                    "axis {axis:?}: declared envelope exceeds the grid domain — nothing may be valid where nothing is tabulated"
                )
            }
            Self::ShapeMismatch {
                value,
                expected,
                found,
            } => {
                write!(
                    f,
                    "value {value:?}: {found} elements, axes imply {expected}"
                )
            }
            Self::BadInterpRule {
                value,
                rule,
                reason,
            } => {
                write!(f, "value {value:?}: interp_rule {rule:?} invalid: {reason}")
            }
            Self::LogScaleRequiresPositive {
                value,
                axis_or_value,
            } => {
                write!(
                    f,
                    "value {value:?}: log scale on {axis_or_value:?} requires strictly positive data"
                )
            }
            Self::UnknownValue { name } => write!(f, "no value dataset named {name:?}"),
            Self::UnitsMismatch {
                value,
                expected,
                found,
            } => {
                write!(
                    f,
                    "value {value:?}: consumer expects units {expected:?} but the table declares \
                     {found:?} — refusing the bind (META-2 §4 ★); a relabeled column must never \
                     be silently misread"
                )
            }
            Self::QueryArity { expected, found } => {
                write!(
                    f,
                    "query has {found} coordinates, table has {expected} axes"
                )
            }
            Self::OutOfEnvelope {
                axis,
                query,
                envelope_min,
                envelope_max,
            } => {
                write!(
                    f,
                    "OutOfEnvelope: axis {axis:?} query {query} outside declared validity [{envelope_min}, {envelope_max}] — refusing, never clamping (§3.5)"
                )
            }
            Self::OutOfDomain {
                axis,
                query,
                grid_min,
                grid_max,
            } => {
                write!(
                    f,
                    "axis {axis:?} query {query} outside tabulated domain [{grid_min}, {grid_max}] — the runtime never extrapolates (§3.5)"
                )
            }
            Self::NotYetImplemented { what, arrives_with } => {
                write!(
                    f,
                    "{what} is not yet implemented (arrives with {arrives_with}) — refusing to load rather than guessing"
                )
            }
            Self::UnknownKind { kind } => {
                write!(
                    f,
                    "unknown table kind {kind:?}; known: regular, thermo_potential, sample_set"
                )
            }
            Self::UnknownInterpMethod { method } => {
                write!(
                    f,
                    "unknown interp_method {method:?}; known: multilinear, pchip"
                )
            }
            Self::UnexpectedMember { path, reason } => {
                write!(
                    f,
                    "{path}: unexpected member — {reason}; the digest must cover every byte \
                     (§3.2), so unaccounted content refuses to load"
                )
            }
        }
    }
}

impl std::error::Error for TableError {}
