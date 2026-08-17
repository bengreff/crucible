"""FND-5 table seam, Python side: model, canonical digest v2, h5py writer.

This module is the OFFL half of the cross-language contract. The Rust half
is ``crates/tables/src/digest.rs`` (the digest) and ``reader.rs`` (the
on-disk schema); **keep the three in lockstep**. A table this writer stamps
must open under ``Table::open`` with the digest it returned as the pin —
that round-trip is CI-enforced by the committed seam fixture
(``tests/test_seam_fixture.py`` ↔ ``crates/tables/tests/fnd5_python_seam.rs``).

Digest = ``"sha256:" + hex(SHA-256(stream))`` where ``stream`` is:

1. the ASCII tag ``crucible-table-digest-v2\\n``;
2. for each **axis in table order**: the tag ``axis\\n``, then STR(name),
   point count as u64 LE, points as f64 LE, envelope min then max as f64 LE;
3. for each **value dataset sorted by name** (byte-wise on UTF-8): the tag
   ``value\\n``, then STR(name), STR(units), STR(interp_rule),
   interp_error_bound as f64 LE, element count as u64 LE, data as f64 LE,
   then the sigma marker: ``0u8`` (none), or ``1u8`` + per-point sigma as
   f64 LE, or ``2u8`` + scalar sigma as f64 LE (per-point wins if both);
4. the trailer ``meta\\n``, then STR(kind), STR(data_version),
   STR(interp_method), STR(producer), STR(producer_version),
   STR(input_deck_hash), STR(source_library), STR(generator_commit), and
   the rng-seed marker: ``0u8``, or ``1u8`` + seed as i64 LE.

``STR(s)`` = UTF-8 byte length as u64 LE, then the bytes — length-prefixed,
never delimiter-terminated (the v1 newline lesson). Floats hash by their
IEEE-754 bit pattern (``struct.pack('<d', x)``) — exact, no formatting.
"""

from __future__ import annotations

import hashlib
import struct
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Provenance:
    """§3.1 mandatory provenance — identical fields to the Rust struct."""

    producer: str
    producer_version: str
    input_deck_hash: str
    source_library: str
    generator_commit: str
    rng_seed: int | None = None


@dataclass(frozen=True)
class Axis:
    """One coordinate axis: explicit strictly-increasing grid + envelope."""

    name: str
    points: tuple[float, ...]
    envelope_min: float
    envelope_max: float


@dataclass(frozen=True)
class TableValue:
    """One value dataset with its per-column semantics (§3.1)."""

    name: str
    data: tuple[float, ...]
    units: str
    interp_rule: str
    interp_error_bound: float
    sigma: tuple[float, ...] | None = None
    sigma_scalar: float | None = None


@dataclass(frozen=True)
class WriteSpec:
    """Everything a producer must supply — exactly the §3.1 mandatory set."""

    kind: str
    schema_version: str
    data_version: str
    interp_method: str
    provenance: Provenance
    axes: tuple[Axis, ...] = field(default_factory=tuple)
    values: tuple[TableValue, ...] = field(default_factory=tuple)


def _put_str(h, s: str) -> None:
    b = s.encode("utf-8")
    h.update(struct.pack("<Q", len(b)))
    h.update(b)


def _put_f64(h, x: float) -> None:
    h.update(struct.pack("<d", x))


def content_digest(
    kind: str,
    data_version: str,
    interp_method: str,
    provenance: Provenance,
    axes: tuple[Axis, ...],
    values_sorted_by_name: tuple[TableValue, ...],
) -> str:
    """The canonical content digest a pin verifies (digest.rs mirror)."""
    names = [v.name.encode("utf-8") for v in values_sorted_by_name]
    assert names == sorted(names), "digest requires values sorted by name"

    h = hashlib.sha256()
    h.update(b"crucible-table-digest-v2\n")
    for a in axes:
        h.update(b"axis\n")
        _put_str(h, a.name)
        h.update(struct.pack("<Q", len(a.points)))
        for p in a.points:
            _put_f64(h, p)
        _put_f64(h, a.envelope_min)
        _put_f64(h, a.envelope_max)
    for v in values_sorted_by_name:
        h.update(b"value\n")
        _put_str(h, v.name)
        _put_str(h, v.units)
        _put_str(h, v.interp_rule)
        _put_f64(h, v.interp_error_bound)
        h.update(struct.pack("<Q", len(v.data)))
        for x in v.data:
            _put_f64(h, x)
        if v.sigma is not None:
            h.update(b"\x01")
            for x in v.sigma:
                _put_f64(h, x)
        elif v.sigma_scalar is not None:
            h.update(b"\x02")
            _put_f64(h, v.sigma_scalar)
        else:
            h.update(b"\x00")
    h.update(b"meta\n")
    for s in (
        kind,
        data_version,
        interp_method,
        provenance.producer,
        provenance.producer_version,
        provenance.input_deck_hash,
        provenance.source_library,
        provenance.generator_commit,
    ):
        _put_str(h, s)
    if provenance.rng_seed is not None:
        h.update(b"\x01")
        h.update(struct.pack("<q", provenance.rng_seed))
    else:
        h.update(b"\x00")
    return "sha256:" + h.hexdigest()


def _sorted_values(spec: WriteSpec) -> tuple[TableValue, ...]:
    return tuple(sorted(spec.values, key=lambda v: v.name.encode("utf-8")))


def spec_digest(spec: WriteSpec) -> str:
    """Digest of a spec without writing it (values sorted here)."""
    return content_digest(
        spec.kind,
        spec.data_version,
        spec.interp_method,
        spec.provenance,
        spec.axes,
        _sorted_values(spec),
    )


def write_table(file_path: str, group_path: str, spec: WriteSpec) -> str:
    """Write one table group to the reader.rs on-disk schema; return the
    canonical ``sha256:`` digest for pinning.

    Imports h5py lazily so the digest half of the module stays stdlib-pure
    (usable anywhere, e.g. to re-derive a pin without the HDF5 stack).
    """
    import h5py
    import numpy as np

    for v in spec.values:
        if v.name.startswith("sigma_"):
            raise ValueError(
                f"{group_path}/values/{v.name}: `sigma_` is reserved for "
                "uncertainty companions (§3.1); name the value outside the "
                "reserved prefix"
            )
    for a in spec.axes:
        if not a.name or "\n" in a.name:
            raise ValueError(
                f"{group_path}/axes/{a.name!r}: axis names must be non-empty "
                "and newline-free (the `axis_order` attribute is \\n-joined)"
            )
    values = _sorted_values(spec)

    utf8 = h5py.string_dtype(encoding="utf-8")

    def put_attr(obj, name: str, s: str) -> None:
        obj.attrs.create(name, s, dtype=utf8)

    with h5py.File(file_path, "a") as f:
        g = f.create_group(group_path)
        put_attr(g, "kind", spec.kind)
        put_attr(g, "schema_version", spec.schema_version)
        put_attr(g, "data_version", spec.data_version)
        put_attr(g, "interp_method", spec.interp_method)
        put_attr(g, "producer", spec.provenance.producer)
        put_attr(g, "producer_version", spec.provenance.producer_version)
        put_attr(g, "input_deck_hash", spec.provenance.input_deck_hash)
        put_attr(g, "source_library", spec.provenance.source_library)
        put_attr(g, "generator_commit", spec.provenance.generator_commit)
        if spec.provenance.rng_seed is not None:
            g.attrs.create("rng_seed", np.int64(spec.provenance.rng_seed))
        put_attr(g, "axis_order", "\n".join(a.name for a in spec.axes))

        axes_g = g.create_group("axes")
        for a in spec.axes:
            ds = axes_g.create_dataset(a.name, data=np.asarray(a.points, dtype=np.float64))
            ds.attrs.create("envelope_min", np.float64(a.envelope_min))
            ds.attrs.create("envelope_max", np.float64(a.envelope_max))

        values_g = g.create_group("values")
        for v in values:
            ds = values_g.create_dataset(v.name, data=np.asarray(v.data, dtype=np.float64))
            put_attr(ds, "units", v.units)
            put_attr(ds, "interp_rule", v.interp_rule)
            ds.attrs.create("interp_error_bound", np.float64(v.interp_error_bound))
            if v.sigma is not None:
                values_g.create_dataset(
                    f"sigma_{v.name}", data=np.asarray(v.sigma, dtype=np.float64)
                )
            if v.sigma_scalar is not None:
                ds.attrs.create("sigma", np.float64(v.sigma_scalar))

    return content_digest(
        spec.kind,
        spec.data_version,
        spec.interp_method,
        spec.provenance,
        spec.axes,
        values,
    )
