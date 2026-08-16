//! Test/dev table writer. **Production tables come from the Python OFFL
//! pipelines** (the h5py side of the seam); this writer exists so the §6-1
//! round-trip test and grid/solver fixtures are hermetic in Rust CI, and so
//! the canonical digest has two independent call sites on day one.
//!
//! Deliberately permissive about `kind`/`interp_method` values: producer
//! honesty is enforced by the *reader* (a deferred or unknown kind must
//! refuse at load, and tests need to write such tables).

use crate::digest::content_digest;
use crate::model::{Axis, Provenance, TableError, TableValue};
use hdf5::types::VarLenUnicode;
use std::str::FromStr;

fn h5err(e: hdf5::Error) -> TableError {
    TableError::Hdf5(e.to_string())
}

fn vlu(s: &str) -> Result<VarLenUnicode, TableError> {
    VarLenUnicode::from_str(s).map_err(|e| TableError::Hdf5(format!("bad string: {e}")))
}

fn write_str_attr(loc: &hdf5::Group, name: &str, value: &str) -> Result<(), TableError> {
    loc.new_attr::<VarLenUnicode>()
        .create(name)
        .map_err(h5err)?
        .write_scalar(&vlu(value)?)
        .map_err(h5err)
}

fn write_str_attr_ds(ds: &hdf5::Dataset, name: &str, value: &str) -> Result<(), TableError> {
    ds.new_attr::<VarLenUnicode>()
        .create(name)
        .map_err(h5err)?
        .write_scalar(&vlu(value)?)
        .map_err(h5err)
}

/// Everything a producer must supply — exactly the §3.1 mandatory set.
#[derive(Debug, Clone)]
pub struct WriteSpec {
    pub kind: String,
    pub schema_version: String,
    pub data_version: String,
    pub interp_method: String,
    pub provenance: Provenance,
    pub axes: Vec<Axis>,
    pub values: Vec<TableValue>,
}

/// Write one table group; returns the canonical `sha256:` content digest for
/// pinning. Values are stored sorted by name; value data is written flat
/// row-major (the reader accepts any shape with the right element count).
pub fn write_table(
    file_path: &str,
    group_path: &str,
    spec: &WriteSpec,
) -> Result<String, TableError> {
    // Mirror the reader's schema contract so writer and reader can never
    // disagree about a file's contents (review finding): the `sigma_`
    // prefix is reserved for uncertainty companions, and axis names must
    // survive the `\n`-joined `axis_order` encoding.
    for v in &spec.values {
        if v.name.starts_with("sigma_") {
            return Err(TableError::UnexpectedMember {
                path: format!("{group_path}/values/{}", v.name),
                reason: "`sigma_` is reserved for uncertainty companions (§3.1); name the \
                         value outside the reserved prefix"
                    .into(),
            });
        }
    }
    for a in &spec.axes {
        if a.name.contains('\n') || a.name.is_empty() {
            return Err(TableError::UnexpectedMember {
                path: format!("{group_path}/axes/{}", a.name),
                reason: "axis names must be non-empty and newline-free (the `axis_order` \
                         attribute is `\\n`-joined)"
                    .into(),
            });
        }
    }
    let mut values = spec.values.clone();
    values.sort_by(|a, b| a.name.cmp(&b.name));

    let file = hdf5::File::append(file_path).map_err(h5err)?;
    let group = file.create_group(group_path).map_err(h5err)?;

    write_str_attr(&group, "kind", &spec.kind)?;
    write_str_attr(&group, "schema_version", &spec.schema_version)?;
    write_str_attr(&group, "data_version", &spec.data_version)?;
    write_str_attr(&group, "interp_method", &spec.interp_method)?;
    write_str_attr(&group, "producer", &spec.provenance.producer)?;
    write_str_attr(
        &group,
        "producer_version",
        &spec.provenance.producer_version,
    )?;
    write_str_attr(&group, "input_deck_hash", &spec.provenance.input_deck_hash)?;
    write_str_attr(&group, "source_library", &spec.provenance.source_library)?;
    write_str_attr(
        &group,
        "generator_commit",
        &spec.provenance.generator_commit,
    )?;
    if let Some(seed) = spec.provenance.rng_seed {
        group
            .new_attr::<i64>()
            .create("rng_seed")
            .map_err(h5err)?
            .write_scalar(&seed)
            .map_err(h5err)?;
    }
    let order: Vec<&str> = spec.axes.iter().map(|a| a.name.as_str()).collect();
    write_str_attr(&group, "axis_order", &order.join("\n"))?;

    let axes_group = group.create_group("axes").map_err(h5err)?;
    for a in &spec.axes {
        let ds = axes_group
            .new_dataset_builder()
            .with_data(&a.points)
            .create(a.name.as_str())
            .map_err(h5err)?;
        for (attr, v) in [
            ("envelope_min", a.envelope_min),
            ("envelope_max", a.envelope_max),
        ] {
            ds.new_attr::<f64>()
                .create(attr)
                .map_err(h5err)?
                .write_scalar(&v)
                .map_err(h5err)?;
        }
    }

    let values_group = group.create_group("values").map_err(h5err)?;
    for v in &values {
        let ds = values_group
            .new_dataset_builder()
            .with_data(&v.data)
            .create(v.name.as_str())
            .map_err(h5err)?;
        write_str_attr_ds(&ds, "units", &v.units)?;
        write_str_attr_ds(&ds, "interp_rule", &v.interp_rule)?;
        ds.new_attr::<f64>()
            .create("interp_error_bound")
            .map_err(h5err)?
            .write_scalar(&v.interp_error_bound)
            .map_err(h5err)?;
        if let Some(s) = &v.sigma {
            values_group
                .new_dataset_builder()
                .with_data(s)
                .create(format!("sigma_{}", v.name).as_str())
                .map_err(h5err)?;
        }
        if let Some(s) = v.sigma_scalar {
            ds.new_attr::<f64>()
                .create("sigma")
                .map_err(h5err)?
                .write_scalar(&s)
                .map_err(h5err)?;
        }
    }

    Ok(content_digest(
        &spec.kind,
        &spec.data_version,
        &spec.interp_method,
        &spec.provenance,
        &spec.axes,
        &values,
    ))
}
