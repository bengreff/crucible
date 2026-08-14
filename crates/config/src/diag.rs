//! FND-4 §3.4 error accumulation: the semantic pass collects **every**
//! problem and reports them together, in deterministic order (sorted by
//! dotted path, then message — the COUP-8 §3.3 `(mechanism-id, field)`
//! ordering generalized to the whole document). Never first-error-only,
//! never hash order.

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Diagnostic {
    /// Dotted TOML path of the offending item, e.g. `mechanisms.hydro.type`.
    /// Empty for document-level faults (parse failure, bad schema_version).
    pub path: String,
    pub message: String,
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Diagnostics(Vec<Diagnostic>);

impl Diagnostics {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, path: impl Into<String>, message: impl Into<String>) {
        self.0.push(Diagnostic {
            path: path.into(),
            message: message.into(),
        });
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Canonical order: (path, message). Called once at emission so the
    /// report order is a pure function of the config, not of check order.
    pub(crate) fn into_sorted(mut self) -> Self {
        self.0.sort_by(|a, b| {
            (a.path.as_str(), a.message.as_str()).cmp(&(b.path.as_str(), b.message.as_str()))
        });
        self
    }

    pub fn iter(&self) -> impl Iterator<Item = &Diagnostic> {
        self.0.iter()
    }
}

impl fmt::Display for Diagnostics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "config refused to load ({} fault(s)):", self.0.len())?;
        for d in &self.0 {
            if d.path.is_empty() {
                writeln!(f, "  - {}", d.message)?;
            } else {
                writeln!(f, "  - [{}] {}", d.path, d.message)?;
            }
        }
        Ok(())
    }
}

impl std::error::Error for Diagnostics {}
