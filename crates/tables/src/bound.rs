//! FND-5 §3.3 — `BoundColumn`: a column handle resolved **once, at bind
//! time**, for kernel-rate interrogation (SOLV-1 §3.4 makes ~10⁸ equilibrium
//! projections per anchor run; the general [`Table::interpolate`] re-parses
//! the `interp_rule` grammar and allocates per query, which is fine for
//! diagnostics and forbidden in a per-cell loop).
//!
//! Binding performs, once: the META-2 §4 ★ units gate (`expect_units`), the
//! value lookup, the `interp_rule` parse, stride precomputation, and
//! pre-linearization of `log`-scaled axes/values (their `ln` is taken at
//! bind). Queries then run the **identical fixed-order multilinear kernel**
//! as the general path — same partition points, same corner order, same
//! accumulation order — so a bound query is **bit-identical** to
//! `Table::interpolate` (asserted by test), just without the per-query parse
//! and heap traffic. Envelope policy is `Refuse` (the headline-run default);
//! `Flag`-policy ensemble bookkeeping stays on the general path (COUP-5).

use crate::interp::{Scale, parse_rule};
use crate::model::{Table, TableError};

/// Hard cap on bound-query dimensionality: query scratch lives on the stack
/// (no per-query allocation). Every FND-5 surface in the project is ≤ 4-D
/// ((p, h, Z) + one spare); a wider table binds with a load-style refusal.
pub const MAX_BOUND_AXES: usize = 8;

/// **S13c GPU marshaling** (doc-hidden): the flat arrays a device kernel needs
/// to reproduce [`BoundColumn::interpolate`] bit-for-formula. Owned (copied out
/// of the borrowed table) so it can outlive the bind for upload.
#[doc(hidden)]
#[derive(Debug, Clone)]
pub struct ColumnMarshal {
    /// Grid points per axis, in table order.
    pub axis_points: Vec<Vec<f64>>,
    /// Per-axis: is the interp rule `log` on that axis?
    pub axis_is_log: Vec<bool>,
    /// Row-major strides over the axes (table order).
    pub strides: Vec<usize>,
    /// The value data block, row-major.
    pub data: Vec<f64>,
    /// Is the value column `log`-scaled?
    pub value_is_log: bool,
}

/// A pre-resolved, allocation-free view of one value column of a [`Table`].
/// Borrows the table (a `Table` is immutable after load, FND-5 §3.6).
#[derive(Debug, Clone)]
pub struct BoundColumn<'t> {
    table: &'t Table,
    value_index: usize,
    axis_scales: Vec<Scale>,
    value_scale: Scale,
    /// Row-major strides over the axes in table order.
    strides: Vec<usize>,
    /// Per log axis: `ln` of its grid points (`None` for lin axes).
    ln_points: Vec<Option<Vec<f64>>>,
    /// For a log-scaled value: `ln` of the whole data block, taken at bind —
    /// the same `ln(f)` the general kernel takes per corner, just hoisted.
    ln_data: Option<Vec<f64>>,
}

impl Table {
    /// Bind a column by `(name, expected SI units)` — the META-2 §4 ★ units
    /// gate runs here, once, so a relabeled column refuses at bind, never
    /// misreads at kernel rate.
    pub fn bind(&self, value_name: &str, units: &str) -> Result<BoundColumn<'_>, TableError> {
        self.expect_units(value_name, units)?;
        let value_index = self
            .values
            .iter()
            .position(|v| v.name == value_name)
            .expect("expect_units verified existence");
        let value = &self.values[value_index];
        let n = self.axes.len();
        if n > MAX_BOUND_AXES {
            return Err(TableError::BadInterpRule {
                value: value_name.to_string(),
                rule: value.interp_rule.clone(),
                reason: format!("{n} axes exceed the {MAX_BOUND_AXES}-axis bound-query cap"),
            });
        }
        let (axis_scales, value_scale) = parse_rule(&value.name, &value.interp_rule, n)?;

        let mut strides = vec![1usize; n];
        for i in (0..n.saturating_sub(1)).rev() {
            strides[i] = strides[i + 1] * self.axes[i + 1].points.len();
        }

        let ln_points = self
            .axes
            .iter()
            .zip(&axis_scales)
            .map(|(a, s)| match s {
                Scale::Lin => None,
                Scale::Log => Some(a.points.iter().map(|p| p.ln()).collect()),
            })
            .collect();
        let ln_data = match value_scale {
            Scale::Lin => None,
            Scale::Log => Some(value.data.iter().map(|f| f.ln()).collect()),
        };

        Ok(BoundColumn {
            table: self,
            value_index,
            axis_scales,
            value_scale,
            strides,
            ln_points,
            ln_data,
        })
    }
}

impl BoundColumn<'_> {
    /// The bound column's name (diagnostics).
    pub fn name(&self) -> &str {
        &self.table.values[self.value_index].name
    }

    /// §3.4: the producer-measured interpolation-error bound of this column
    /// (table metadata for the UQ budget — never folded into a return value).
    pub fn interp_error_bound(&self) -> f64 {
        self.table.values[self.value_index].interp_error_bound
    }

    /// The rule-space (|Δ ln|, ≈ relative) bound for log-valued columns —
    /// `None` on linear-valued columns and pre-0.3.2 artifacts (see the
    /// model's field doc; the recorded digest-v4 deferral rides there).
    pub fn interp_error_bound_log(&self) -> Option<f64> {
        self.table.values[self.value_index].interp_error_bound_log
    }

    /// The grid-domain span of one axis (for consumers that must size
    /// iteration brackets inside the tabulated domain, e.g. the SOLV-1 §3.4
    /// p-iteration).
    pub fn axis_domain(&self, axis: usize) -> (f64, f64) {
        let pts = &self.table.axes[axis].points;
        (pts[0], pts[pts.len() - 1])
    }

    /// The declared validity envelope of one axis (§3.5; may be tighter than
    /// the grid domain — queries outside it refuse).
    pub fn axis_envelope(&self, axis: usize) -> (f64, f64) {
        let a = &self.table.axes[axis];
        (a.envelope_min, a.envelope_max)
    }

    /// The axis name (bind-time schema checks by consumers).
    pub fn axis_name(&self, axis: usize) -> &str {
        &self.table.axes[axis].name
    }

    /// **S13c GPU marshaling** (doc-hidden, additive). Flatten the bound
    /// column into the plain arrays a device kernel needs to reproduce
    /// `interpolate` bit-for-formula: per-axis grid points + a log-scale flag,
    /// the row-major value data + a log-scale flag, and the strides. The
    /// device interp re-does the same fixed 8-corner reduction in the same
    /// rule-space, so it is a per-cell gather (CPU↔GPU = FMA-order ECT). The
    /// CPU path here is untouched — this only reads it out.
    #[doc(hidden)]
    pub fn marshal(&self) -> ColumnMarshal {
        ColumnMarshal {
            axis_points: self.table.axes.iter().map(|a| a.points.clone()).collect(),
            axis_is_log: self.axis_scales.iter().map(|s| matches!(s, Scale::Log)).collect(),
            strides: self.strides.clone(),
            data: self.table.values[self.value_index].data.clone(),
            value_is_log: matches!(self.value_scale, Scale::Log),
        }
    }

    /// §3.3 `Interpolate(query)` under `Refuse` policy — deterministic
    /// scalar or hard refusal, allocation-free on the success path, and
    /// bit-identical to [`Table::interpolate`] on the same query.
    pub fn interpolate(&self, query: &[f64]) -> Result<f64, TableError> {
        let axes = &self.table.axes;
        let n = axes.len();
        if query.len() != n {
            return Err(TableError::QueryArity {
                expected: n,
                found: query.len(),
            });
        }

        // Envelope + domain checks, axes in fixed table order (§3.5),
        // `Refuse` policy: any excursion is a hard error.
        for (a, &q) in axes.iter().zip(query) {
            let (grid_min, grid_max) = (a.points[0], a.points[a.points.len() - 1]);
            if !q.is_finite() || q < grid_min || q > grid_max {
                return Err(TableError::OutOfDomain {
                    axis: a.name.clone(),
                    query: q,
                    grid_min,
                    grid_max,
                });
            }
            if q < a.envelope_min || q > a.envelope_max {
                return Err(TableError::OutOfEnvelope {
                    axis: a.name.clone(),
                    query: q,
                    envelope_min: a.envelope_min,
                    envelope_max: a.envelope_max,
                });
            }
        }

        // Identical kernel to `Table::multilinear`, scratch on the stack.
        let mut cell = [0usize; MAX_BOUND_AXES];
        let mut t = [0.0f64; MAX_BOUND_AXES];
        for (d, ((a, &q), scale)) in axes.iter().zip(query).zip(&self.axis_scales).enumerate() {
            let hi = a
                .points
                .partition_point(|&p| p <= q)
                .min(a.points.len() - 1);
            let i = hi.max(1) - 1;
            let frac = match scale {
                Scale::Lin => {
                    let (x0, x1) = (a.points[i], a.points[i + 1]);
                    (q - x0) / (x1 - x0)
                }
                Scale::Log => {
                    let lnp = self.ln_points[d].as_ref().expect("log axis has ln_points");
                    (q.ln() - lnp[i]) / (lnp[i + 1] - lnp[i])
                }
            };
            cell[d] = i;
            t[d] = frac;
        }

        let value = &self.table.values[self.value_index];
        let data: &[f64] = match self.value_scale {
            Scale::Lin => &value.data,
            Scale::Log => self.ln_data.as_ref().expect("log value has ln_data"),
        };

        // Corners in ascending index order — the fixed reduction order.
        let mut acc = 0.0f64;
        for corner in 0..(1usize << n) {
            let mut w = 1.0f64;
            let mut idx = 0usize;
            for d in 0..n {
                let up = (corner >> d) & 1 == 1;
                w *= if up { t[d] } else { 1.0 - t[d] };
                idx += (cell[d] + usize::from(up)) * self.strides[d];
            }
            acc += w * data[idx];
        }
        Ok(match self.value_scale {
            Scale::Lin => acc,
            Scale::Log => acc.exp(),
        })
    }
}
