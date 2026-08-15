//! FND-5 §3.3/§3.5 — deterministic multilinear interpolation on regular
//! (rectilinear) grids, in the linearizing space named by `interp_rule`.
//!
//! Multilinear is the default for every conservation- or positivity-critical
//! quantity: within each cell the result is a convex combination of the
//! 2^N corner values, so it can never overshoot the tabulated bounds.
//! Evaluation order is fixed (axes in table order, corners in ascending
//! index), all `f64`, no explicit FMA — bit-reproducible on a target (§3.6).

use crate::model::{EnvelopeHit, EnvelopePolicy, Table, TableError, TableValue};

/// Per-axis/value scale parsed from the `interp_rule` grammar
/// (`lin`/`log` tokens, axes in order then the value).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Scale {
    Lin,
    Log,
}

pub(crate) fn parse_rule(
    value_name: &str,
    rule: &str,
    n_axes: usize,
) -> Result<(Vec<Scale>, Scale), TableError> {
    let toks: Vec<&str> = rule.split('-').collect();
    if toks.len() != n_axes + 1 {
        return Err(TableError::BadInterpRule {
            value: value_name.to_string(),
            rule: rule.to_string(),
            reason: format!(
                "{} tokens, need {} (one per axis, then the value)",
                toks.len(),
                n_axes + 1
            ),
        });
    }
    let mut scales = Vec::with_capacity(toks.len());
    for t in &toks {
        scales.push(match *t {
            "lin" => Scale::Lin,
            "log" => Scale::Log,
            other => {
                return Err(TableError::BadInterpRule {
                    value: value_name.to_string(),
                    rule: rule.to_string(),
                    reason: format!("unknown token {other:?}; expected \"lin\" or \"log\""),
                });
            }
        });
    }
    let value_scale = scales.pop().expect("len >= 1 by construction");
    Ok((scales, value_scale))
}

impl Table {
    /// §3.3 `Interpolate(query)`: a deterministic **scalar**, or a hard
    /// refusal — never a distribution, never a clamped/extrapolated value.
    /// `Refuse` policy (the headline-run default) errors on any envelope
    /// excursion; use [`Table::interpolate_flagged`] for `Flag`.
    pub fn interpolate(&self, value_name: &str, query: &[f64]) -> Result<f64, TableError> {
        let (v, hits) = self.query_inner(value_name, query, EnvelopePolicy::Refuse)?;
        debug_assert!(hits.is_empty(), "Refuse policy cannot produce hits");
        Ok(v)
    }

    /// §3.5 `Flag` policy: computes the value when the query is inside the
    /// tabulated domain, reporting each envelope excursion for COUP-5 to
    /// record per-member. Outside the grid domain still hard-errors.
    pub fn interpolate_flagged(
        &self,
        value_name: &str,
        query: &[f64],
    ) -> Result<(f64, Vec<EnvelopeHit>), TableError> {
        self.query_inner(value_name, query, EnvelopePolicy::Flag)
    }

    fn query_inner(
        &self,
        value_name: &str,
        query: &[f64],
        policy: EnvelopePolicy,
    ) -> Result<(f64, Vec<EnvelopeHit>), TableError> {
        let value = self
            .values
            .iter()
            .find(|v| v.name == value_name)
            .ok_or_else(|| TableError::UnknownValue {
                name: value_name.to_string(),
            })?;
        if query.len() != self.axes.len() {
            return Err(TableError::QueryArity {
                expected: self.axes.len(),
                found: query.len(),
            });
        }

        // Envelope + domain checks, axes in fixed table order (§3.5).
        let mut hits = Vec::new();
        for (a, &q) in self.axes.iter().zip(query) {
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
                match policy {
                    EnvelopePolicy::Refuse => {
                        return Err(TableError::OutOfEnvelope {
                            axis: a.name.clone(),
                            query: q,
                            envelope_min: a.envelope_min,
                            envelope_max: a.envelope_max,
                        });
                    }
                    EnvelopePolicy::Flag => hits.push(EnvelopeHit {
                        axis: a.name.clone(),
                        query: q,
                        envelope_min: a.envelope_min,
                        envelope_max: a.envelope_max,
                    }),
                }
            }
        }

        Ok((self.multilinear(value, query)?, hits))
    }

    /// Fixed-order multilinear kernel in the `interp_rule` space.
    fn multilinear(&self, value: &TableValue, query: &[f64]) -> Result<f64, TableError> {
        let n = self.axes.len();
        let (axis_scales, value_scale) = parse_rule(&value.name, &value.interp_rule, n)?;

        // Per axis: containing interval index and the local coordinate
        // t ∈ [0,1], computed in the axis's declared space.
        let mut cell = Vec::with_capacity(n);
        let mut t = Vec::with_capacity(n);
        for ((a, &q), scale) in self.axes.iter().zip(query).zip(&axis_scales) {
            // partition_point: first index with point > q; interval is [i-1, i].
            let hi = a
                .points
                .partition_point(|&p| p <= q)
                .min(a.points.len() - 1);
            let i = hi.max(1) - 1;
            let (x0, x1) = (a.points[i], a.points[i + 1]);
            let frac = match scale {
                Scale::Lin => (q - x0) / (x1 - x0),
                Scale::Log => (q.ln() - x0.ln()) / (x1.ln() - x0.ln()),
            };
            cell.push(i);
            t.push(frac);
        }

        // Row-major strides, axes in table order.
        let mut strides = vec![1usize; n];
        for i in (0..n.saturating_sub(1)).rev() {
            strides[i] = strides[i + 1] * self.axes[i + 1].points.len();
        }

        // Corners in ascending index order — the fixed reduction order.
        let mut acc = 0.0f64;
        for corner in 0..(1usize << n) {
            let mut w = 1.0f64;
            let mut idx = 0usize;
            for d in 0..n {
                let up = (corner >> d) & 1 == 1;
                w *= if up { t[d] } else { 1.0 - t[d] };
                idx += (cell[d] + usize::from(up)) * strides[d];
            }
            let f = value.data[idx];
            let fv = match value_scale {
                Scale::Lin => f,
                Scale::Log => f.ln(),
            };
            acc += w * fv;
        }
        Ok(match value_scale {
            Scale::Lin => acc,
            Scale::Log => acc.exp(),
        })
    }
}
