//! FND-5 §3.2 content addressing — the canonical digest a pin verifies.
//!
//! **This algorithm is the cross-language contract.** The Python OFFL
//! pipelines must reproduce it byte-for-byte to stamp `content_digest`s;
//! keep this doc comment and the OFFL implementation in lockstep.
//!
//! Digest = `"sha256:" + hex(SHA-256(stream))` where `stream` is:
//!
//! 1. the ASCII tag `crucible-table-digest-v3\n`;
//! 2. for each **axis in table order**: the tag `axis\n`, then STR(name),
//!    point count as u64 LE, points as f64 LE, envelope min then max as
//!    f64 LE;
//! 3. for each **value dataset sorted by name**: the tag `value\n`, then
//!    STR(name), STR(units), STR(interp_rule), interp_error_bound as
//!    f64 LE, element count as u64 LE, data as f64 LE, then the sigma
//!    marker: `0u8` (none), or `1u8` + per-point sigma as f64 LE, or
//!    `2u8` + scalar sigma as f64 LE — **exactly one form per value**
//!    (writer and reader both refuse a value carrying both: the stream
//!    hashes one marker, so a second form would be unpinned bytes), and
//!    the per-point sigma count is defined equal to the (prefixed) data
//!    count — enforced at every stamping site, keeping the stream
//!    prefix-unambiguous without a second length field;
//! 4. the trailer `meta\n`, then STR(kind), STR(schema_version) — v3:
//!    the reader gates on it, so it is reader-consequential bytes the pin
//!    must cover (review finding; v2 left it unpinned) — STR(data_version),
//!    STR(interp_method), STR(producer), STR(producer_version),
//!    STR(input_deck_hash), STR(source_library), STR(generator_commit),
//!    and the rng-seed marker: `0u8`, or `1u8` + seed as i64 LE.
//!
//! `STR(s)` = UTF-8 byte length as u64 LE, then the bytes — **length-
//! prefixed, never delimiter-terminated**: the v1 encoding used `\n`
//! terminators, so a newline *inside* one field shifted bytes into the next
//! and two different tables could hash identically (review finding). v2
//! makes every field boundary explicit. Floats hash by their IEEE-754 bit
//! pattern (`to_le_bytes`) — exact, no formatting. Provenance is included
//! deliberately: *any* difference from the pinned bytes must fail the pin
//! (the DVC same-label/different-bytes lesson).

use crate::model::{Axis, Provenance, TableValue};
use sha2::{Digest, Sha256};

fn put_str(h: &mut Sha256, s: &str) {
    h.update((s.len() as u64).to_le_bytes());
    h.update(s.as_bytes());
}

pub(crate) fn content_digest(
    kind: &str,
    schema_version: &str,
    data_version: &str,
    interp_method: &str,
    provenance: &Provenance,
    axes: &[Axis],
    values_sorted_by_name: &[TableValue],
) -> String {
    let mut h = Sha256::new();
    h.update(b"crucible-table-digest-v3\n");
    for a in axes {
        h.update(b"axis\n");
        put_str(&mut h, &a.name);
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
        put_str(&mut h, &v.name);
        put_str(&mut h, &v.units);
        put_str(&mut h, &v.interp_rule);
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
    for s in [
        kind,
        schema_version,
        data_version,
        interp_method,
        provenance.producer.as_str(),
        provenance.producer_version.as_str(),
        provenance.input_deck_hash.as_str(),
        provenance.source_library.as_str(),
        provenance.generator_commit.as_str(),
    ] {
        put_str(&mut h, s);
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
