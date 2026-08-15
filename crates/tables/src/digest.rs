//! FND-5 §3.2 content addressing — the canonical digest a pin verifies.
//!
//! **This algorithm is the cross-language contract.** The Python OFFL
//! pipelines must reproduce it byte-for-byte to stamp `content_digest`s;
//! keep this doc comment and the OFFL implementation in lockstep.
//!
//! Digest = `"sha256:" + hex(SHA-256(stream))` where `stream` is:
//!
//! 1. the ASCII tag `crucible-table-digest-v1\n`;
//! 2. for each **axis in table order**: `axis\n`, name, `\n`, point count as
//!    u64 LE, points as f64 LE, envelope min then max as f64 LE;
//! 3. for each **value dataset sorted by name**: `value\n`, name, `\n`,
//!    units, `\n`, interp_rule, `\n`, interp_error_bound as f64 LE, element
//!    count as u64 LE, data as f64 LE, then the sigma marker: `0u8` (none),
//!    or `1u8` + per-point sigma as f64 LE, or `2u8` + scalar sigma f64 LE;
//! 4. the trailer `meta\n`, kind, `\n`, data_version, `\n`, interp_method,
//!    `\n`, then each provenance field (producer, producer_version,
//!    input_deck_hash, source_library, generator_commit) each followed by
//!    `\n`, and the rng-seed marker: `0u8`, or `1u8` + seed as i64 LE.
//!
//! Floats hash by their IEEE-754 bit pattern (`to_le_bytes`), so the digest
//! is exact — no formatting, no rounding. Provenance is included
//! deliberately: *any* difference from the pinned bytes must fail the pin
//! (the DVC same-label/different-bytes lesson, §3.6 of FND-4).

use crate::model::{Axis, Provenance, TableValue};
use sha2::{Digest, Sha256};

pub(crate) fn content_digest(
    kind: &str,
    data_version: &str,
    interp_method: &str,
    provenance: &Provenance,
    axes: &[Axis],
    values_sorted_by_name: &[TableValue],
) -> String {
    let mut h = Sha256::new();
    h.update(b"crucible-table-digest-v1\n");
    for a in axes {
        h.update(b"axis\n");
        h.update(a.name.as_bytes());
        h.update(b"\n");
        h.update((a.points.len() as u64).to_le_bytes());
        for p in &a.points {
            h.update(p.to_le_bytes());
        }
        h.update(a.envelope_min.to_le_bytes());
        h.update(a.envelope_max.to_le_bytes());
    }
    debug_assert!(
        values_sorted_by_name
            .windows(2)
            .all(|w| w[0].name < w[1].name),
        "digest requires values sorted by name"
    );
    for v in values_sorted_by_name {
        h.update(b"value\n");
        h.update(v.name.as_bytes());
        h.update(b"\n");
        h.update(v.units.as_bytes());
        h.update(b"\n");
        h.update(v.interp_rule.as_bytes());
        h.update(b"\n");
        h.update(v.interp_error_bound.to_le_bytes());
        h.update((v.data.len() as u64).to_le_bytes());
        for x in &v.data {
            h.update(x.to_le_bytes());
        }
        match (&v.sigma, v.sigma_scalar) {
            (Some(s), _) => {
                h.update([1u8]);
                for x in s {
                    h.update(x.to_le_bytes());
                }
            }
            (None, Some(s)) => {
                h.update([2u8]);
                h.update(s.to_le_bytes());
            }
            (None, None) => h.update([0u8]),
        }
    }
    h.update(b"meta\n");
    for s in [kind, data_version, interp_method] {
        h.update(s.as_bytes());
        h.update(b"\n");
    }
    for s in [
        provenance.producer.as_str(),
        provenance.producer_version.as_str(),
        provenance.input_deck_hash.as_str(),
        provenance.source_library.as_str(),
        provenance.generator_commit.as_str(),
    ] {
        h.update(s.as_bytes());
        h.update(b"\n");
    }
    match provenance.rng_seed {
        Some(seed) => {
            h.update([1u8]);
            h.update(seed.to_le_bytes());
        }
        None => h.update([0u8]),
    }
    let mut hex = String::with_capacity(71);
    hex.push_str("sha256:");
    for b in h.finalize() {
        use std::fmt::Write;
        write!(hex, "{b:02x}").expect("writing to String cannot fail");
    }
    hex
}
