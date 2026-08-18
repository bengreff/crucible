//! FND-5 §3.2 — the load-time contract. `Table::open` verifies the pin
//! (data_version + content digest), the schema major version, the full
//! mandatory-metadata set, and the envelope declarations — any failure
//! **halts with a diagnosis** (META-1 §5: reproducibility by construction).
//!
//! Concrete on-disk schema (the §3.1 realization; the Python OFFL side
//! mirrors this exactly):
//!
//! - group attributes: `kind`, `schema_version` (`"major.minor"`),
//!   `data_version`, `interp_method`, `axis_order` (axis names joined by
//!   `\n`), and the provenance set `producer`, `producer_version`,
//!   `input_deck_hash`, `source_library`, `generator_commit`
//!   (+ optional i64 `rng_seed` iff the generator was stochastic);
//! - `axes/<name>`: 1-D f64 grid points, attrs `envelope_min`/`envelope_max`;
//! - `values/<name>`: f64 data (any HDF5 shape; read row-major flat, element
//!   count must equal the axis-length product), attrs `units`,
//!   `interp_rule`, `interp_error_bound`, optional scalar attr `sigma`;
//!   optional per-point dataset `values/sigma_<name>`.

use crate::digest::content_digest;
use crate::interp::parse_rule;
use crate::model::{Axis, Provenance, TABLE_SCHEMA_MAJOR, Table, TableError, TableValue};
use hdf5::types::VarLenUnicode;

/// The config-side pin (FND-4 `[tables]`): the loader refuses any mismatch.
/// `content_digest: None` means "first load" — the caller pins the digest
/// the load reports (it is always computed and returned on the `Table`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Pin {
    pub data_version: String,
    pub content_digest: Option<String>,
}

fn h5err(e: hdf5::Error) -> TableError {
    TableError::Hdf5(e.to_string())
}

fn read_str_attr(loc: &hdf5::Group, path: &str, name: &str) -> Result<String, TableError> {
    let names = loc.attr_names().map_err(h5err)?;
    if !names.iter().any(|n| n == name) {
        return Err(TableError::MissingAttribute {
            path: path.to_string(),
            attribute: name.to_string(),
        });
    }
    let attr = loc.attr(name).map_err(h5err)?;
    let v: VarLenUnicode = attr.read_scalar().map_err(h5err)?;
    Ok(v.as_str().to_string())
}

fn read_f64_attr(ds: &hdf5::Dataset, path: &str, name: &str) -> Result<f64, TableError> {
    let names = ds.attr_names().map_err(h5err)?;
    if !names.iter().any(|n| n == name) {
        return Err(TableError::MissingAttribute {
            path: path.to_string(),
            attribute: name.to_string(),
        });
    }
    ds.attr(name).map_err(h5err)?.read_scalar().map_err(h5err)
}

fn read_str_attr_ds(ds: &hdf5::Dataset, path: &str, name: &str) -> Result<String, TableError> {
    let names = ds.attr_names().map_err(h5err)?;
    if !names.iter().any(|n| n == name) {
        return Err(TableError::MissingAttribute {
            path: path.to_string(),
            attribute: name.to_string(),
        });
    }
    let v: VarLenUnicode = ds.attr(name).map_err(h5err)?.read_scalar().map_err(h5err)?;
    Ok(v.as_str().to_string())
}

impl Table {
    /// Open and fully validate one table group against its pin (§3.2).
    pub fn open(file_path: &str, group_path: &str, pin: &Pin) -> Result<Table, TableError> {
        let file = hdf5::File::open(file_path).map_err(h5err)?;
        let group = file.group(group_path).map_err(h5err)?;

        // Mandatory identity + provenance attributes (§3.1).
        let kind = read_str_attr(&group, group_path, "kind")?;
        let schema_version = read_str_attr(&group, group_path, "schema_version")?;
        let data_version = read_str_attr(&group, group_path, "data_version")?;
        let interp_method = read_str_attr(&group, group_path, "interp_method")?;
        let provenance = Provenance {
            producer: read_str_attr(&group, group_path, "producer")?,
            producer_version: read_str_attr(&group, group_path, "producer_version")?,
            input_deck_hash: read_str_attr(&group, group_path, "input_deck_hash")?,
            source_library: read_str_attr(&group, group_path, "source_library")?,
            generator_commit: read_str_attr(&group, group_path, "generator_commit")?,
            rng_seed: {
                let names = group.attr_names().map_err(h5err)?;
                if names.iter().any(|n| n == "rng_seed") {
                    Some(
                        group
                            .attr("rng_seed")
                            .map_err(h5err)?
                            .read_scalar()
                            .map_err(h5err)?,
                    )
                } else {
                    None
                }
            },
        };

        // Schema major gate (§3.2).
        let major: i64 = schema_version
            .split('.')
            .next()
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| TableError::SchemaMajorMismatch {
                found: schema_version.clone(),
                supported: TABLE_SCHEMA_MAJOR,
            })?;
        if major != TABLE_SCHEMA_MAJOR {
            return Err(TableError::SchemaMajorMismatch {
                found: schema_version,
                supported: TABLE_SCHEMA_MAJOR,
            });
        }

        // Kind / method gates — deferred kinds refuse loudly (lib.rs).
        match kind.as_str() {
            "regular" => {}
            "thermo_potential" => {
                return Err(TableError::NotYetImplemented {
                    what: "kind = \"thermo_potential\" (Helmholtz potential tables, S11)".into(),
                    arrives_with: "the FND-7/OFFL-5 constitutive-spine session".into(),
                });
            }
            "sample_set" => {
                return Err(TableError::NotYetImplemented {
                    what: "kind = \"sample_set\" (empirical sample-sets, FND-1 §3.9)".into(),
                    arrives_with: "the OFFL-2 nuclear-data session".into(),
                });
            }
            other => {
                return Err(TableError::UnknownKind {
                    kind: other.to_string(),
                });
            }
        }
        match interp_method.as_str() {
            "multilinear" => {}
            "pchip" => {
                return Err(TableError::NotYetImplemented {
                    what: "interp_method = \"pchip\" (monotone-cubic shipped coefficients)".into(),
                    arrives_with: "the first C¹ consumer (EOS/reactivity differencing)".into(),
                });
            }
            other => {
                return Err(TableError::UnknownInterpMethod {
                    method: other.to_string(),
                });
            }
        }

        // Strict membership (review finding): every member of the table
        // group must be accounted for, or the digest would not cover the
        // file's bytes and same-label/different-bytes could pass the pin.
        let mut top = group.member_names().map_err(h5err)?;
        top.sort();
        for m in &top {
            if m != "axes" && m != "values" {
                return Err(TableError::UnexpectedMember {
                    path: format!("{group_path}/{m}"),
                    reason: "a table group holds exactly the `axes` and `values` subgroups".into(),
                });
            }
        }

        // Axes, in the declared order.
        let axis_order = read_str_attr(&group, group_path, "axis_order")?;
        let axes_group = group.group("axes").map_err(h5err)?;
        let declared: Vec<&str> = axis_order.split('\n').filter(|s| !s.is_empty()).collect();
        let mut axis_members = axes_group.member_names().map_err(h5err)?;
        axis_members.sort();
        for m in &axis_members {
            if !declared.contains(&m.as_str()) {
                return Err(TableError::UnexpectedMember {
                    path: format!("{group_path}/axes/{m}"),
                    reason: "dataset not listed in `axis_order`".into(),
                });
            }
        }
        let mut axes = Vec::new();
        for name in axis_order.split('\n').filter(|s| !s.is_empty()) {
            let path = format!("{group_path}/axes/{name}");
            let ds = axes_group
                .dataset(name)
                .map_err(|_| TableError::MissingDataset { path: path.clone() })?;
            let points: Vec<f64> = ds.read_raw().map_err(h5err)?;
            if points.len() < 2 || points.windows(2).any(|w| w[0] >= w[1]) {
                return Err(TableError::NonMonotonicAxis {
                    axis: name.to_string(),
                });
            }
            let envelope_min = read_f64_attr(&ds, &path, "envelope_min")?;
            let envelope_max = read_f64_attr(&ds, &path, "envelope_max")?;
            if envelope_min > envelope_max
                || envelope_min < points[0]
                || envelope_max > points[points.len() - 1]
            {
                return Err(TableError::EnvelopeOutsideGrid {
                    axis: name.to_string(),
                });
            }
            axes.push(Axis {
                name: name.to_string(),
                points,
                envelope_min,
                envelope_max,
            });
        }
        let expected_len: usize = axes.iter().map(|a| a.points.len()).product();

        // Value datasets, sorted by name (deterministic §3.6; digest order).
        // The `sigma_` prefix is RESERVED for uncertainty companions
        // (§3.1: "a sibling sigma_<name> uncertainty dataset"); a
        // `sigma_<x>` with no value `<x>` is refused loudly — the v1 reader
        // silently dropped it, so a physics value named e.g. `sigma_t`
        // vanished and the writer/reader digests disagreed (review finding).
        // Producers must name values outside the reserved prefix.
        let values_group = group.group("values").map_err(h5err)?;
        let mut member_names = values_group.member_names().map_err(h5err)?;
        member_names.sort();
        for m in member_names.iter().filter(|n| n.starts_with("sigma_")) {
            let base = &m["sigma_".len()..];
            if !member_names.iter().any(|n| n == base) {
                return Err(TableError::UnexpectedMember {
                    path: format!("{group_path}/values/{m}"),
                    reason: format!(
                        "`sigma_` is reserved for uncertainty companions and no value {base:?} \
                         exists; rename the dataset if it is a physical value"
                    ),
                });
            }
        }
        let mut values = Vec::new();
        for name in member_names.iter().filter(|n| !n.starts_with("sigma_")) {
            let path = format!("{group_path}/values/{name}");
            let ds = values_group.dataset(name).map_err(h5err)?;
            let data: Vec<f64> = ds.read_raw().map_err(h5err)?;
            if data.len() != expected_len {
                return Err(TableError::ShapeMismatch {
                    value: name.clone(),
                    expected: expected_len,
                    found: data.len(),
                });
            }
            let units = read_str_attr_ds(&ds, &path, "units")?;
            let interp_rule = read_str_attr_ds(&ds, &path, "interp_rule")?;
            let interp_error_bound = read_f64_attr(&ds, &path, "interp_error_bound")?;
            let sigma = {
                let sname = format!("sigma_{name}");
                if member_names.contains(&sname) {
                    let s: Vec<f64> = values_group
                        .dataset(&sname)
                        .map_err(h5err)?
                        .read_raw()
                        .map_err(h5err)?;
                    if s.len() != expected_len {
                        return Err(TableError::ShapeMismatch {
                            value: sname,
                            expected: expected_len,
                            found: s.len(),
                        });
                    }
                    Some(s)
                } else {
                    None
                }
            };
            let sigma_scalar = {
                let names = ds.attr_names().map_err(h5err)?;
                if names.iter().any(|n| n == "sigma") {
                    Some(
                        ds.attr("sigma")
                            .map_err(h5err)?
                            .read_scalar()
                            .map_err(h5err)?,
                    )
                } else {
                    None
                }
            };
            // Exactly one sigma form: the digest hashes one marker, so a
            // file carrying both holds bytes the pin does not cover —
            // refuse it (same-label/different-bytes, §3.2; review finding).
            if sigma.is_some() && sigma_scalar.is_some() {
                return Err(TableError::UnexpectedMember {
                    path: path.clone(),
                    reason: "both per-point sigma_ dataset and scalar sigma attribute \
                             present — the digest covers exactly one sigma form, so the \
                             other is unpinned content"
                        .into(),
                });
            }

            // Interp-rule + log-positivity validation (§3.3) at load, not
            // first query — a bad table never loads.
            let (axis_scales, value_scale) = parse_rule(name, &interp_rule, axes.len())?;
            for (a, s) in axes.iter().zip(&axis_scales) {
                if *s == crate::interp::Scale::Log && a.points[0] <= 0.0 {
                    return Err(TableError::LogScaleRequiresPositive {
                        value: name.clone(),
                        axis_or_value: a.name.clone(),
                    });
                }
            }
            if value_scale == crate::interp::Scale::Log && data.iter().any(|&x| x <= 0.0) {
                return Err(TableError::LogScaleRequiresPositive {
                    value: name.clone(),
                    axis_or_value: "value".into(),
                });
            }

            values.push(TableValue {
                name: name.clone(),
                data,
                units,
                interp_rule,
                interp_error_bound,
                sigma,
                sigma_scalar,
            });
        }

        // Content digest + pin verification (§3.2, S6).
        let digest = content_digest(
            &kind,
            &schema_version,
            &data_version,
            &interp_method,
            &provenance,
            &axes,
            &values,
        );
        if data_version != pin.data_version {
            return Err(TableError::PinVersionMismatch {
                pinned: pin.data_version.clone(),
                found: data_version,
            });
        }
        if let Some(pinned) = &pin.content_digest
            && *pinned != digest
        {
            return Err(TableError::PinDigestMismatch {
                pinned: pinned.clone(),
                found: digest,
            });
        }

        Ok(Table {
            group_path: group_path.to_string(),
            kind,
            schema_version,
            data_version,
            interp_method,
            provenance,
            axes,
            values,
            content_digest: digest,
        })
    }
}
