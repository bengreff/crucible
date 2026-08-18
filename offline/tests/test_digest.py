"""Digest v3 mirror tests. The golden vector is the cross-language contract:
it must equal the one asserted in
``crates/tables/tests/fnd5_s6.rs::fnd5_s32_digest_golden_vector_guards_the_python_contract``.
If either side changes, both must change together (and the tag string bumps).
"""

from crucible_offl.tables import (
    Axis,
    Provenance,
    TableValue,
    WriteSpec,
    spec_digest,
)

GOLDEN = "sha256:a3b0bc5987ad11ba166d006ab7103c16ff1d3a8f17f64a248d4b11d8091023f2"


def provenance(**overrides) -> Provenance:
    base = dict(
        producer="crucible-test",
        producer_version="0.1.0",
        input_deck_hash="sha256:0000",
        source_library="analytic",
        generator_commit="deadbeef",
        rng_seed=None,
    )
    base.update(overrides)
    return Provenance(**base)


def spec(axes, values, **overrides) -> WriteSpec:
    base = dict(
        kind="regular",
        schema_version="1.0",
        data_version="1.0.0",
        interp_method="multilinear",
        provenance=provenance(),
        axes=tuple(axes),
        values=tuple(values),
    )
    base.update(overrides)
    return WriteSpec(**base)


def axis(name, points) -> Axis:
    return Axis(name, tuple(points), points[0], points[-1])


def value(name, data, rule="lin-lin", bound=0.5, **kw) -> TableValue:
    return TableValue(name, tuple(data), "K", rule, bound, **kw)


def test_golden_vector_matches_rust():
    s = spec([axis("x", [0.0, 1.0])], [value("t", [10.0, 20.0])])
    assert spec_digest(s) == GOLDEN


def test_newline_in_provenance_changes_the_digest():
    # The v1 lesson: \n-terminated strings let two different tables hash
    # identically; v2 length-prefixes make these distinct.
    base = spec([axis("x", [0.0, 1.0])], [value("t", [1.0, 2.0], bound=1e-9)])
    a = WriteSpec(
        **{
            **base.__dict__,
            "provenance": provenance(producer="cea\nv2", producer_version="1.0"),
        }
    )
    b = WriteSpec(
        **{
            **base.__dict__,
            "provenance": provenance(producer="cea", producer_version="v2\n1.0"),
        }
    )
    assert spec_digest(a) != spec_digest(b)


def test_sigma_markers_are_distinct():
    ax = [axis("x", [0.0, 1.0])]
    none = spec(ax, [value("t", [1.0, 2.0])])
    per_point = spec(ax, [value("t", [1.0, 2.0], sigma=(0.1, 0.2))])
    scalar = spec(ax, [value("t", [1.0, 2.0], sigma_scalar=0.1)])
    digests = {spec_digest(none), spec_digest(per_point), spec_digest(scalar)}
    assert len(digests) == 3


def test_rng_seed_marker_changes_the_digest():
    ax = [axis("x", [0.0, 1.0])]
    vals = [value("t", [1.0, 2.0])]
    unseeded = spec(ax, vals)
    seeded = spec(ax, vals, provenance=provenance(rng_seed=42))
    zero_seeded = spec(ax, vals, provenance=provenance(rng_seed=0))
    assert len({spec_digest(unseeded), spec_digest(seeded), spec_digest(zero_seeded)}) == 3


def test_value_order_is_by_name_bytes():
    ax = [axis("x", [0.0, 1.0])]
    ab = spec(ax, [value("a", [1.0, 2.0]), value("b", [3.0, 4.0])])
    ba = spec(ax, [value("b", [3.0, 4.0]), value("a", [1.0, 2.0])])
    assert spec_digest(ab) == spec_digest(ba)


def test_float_bit_pattern_negative_zero_is_distinct():
    # Floats hash by IEEE-754 bits, not by value: -0.0 != +0.0 in the stream.
    a = spec([axis("x", [0.0, 1.0])], [value("t", [0.0, 1.0])])
    b = spec([axis("x", [0.0, 1.0])], [value("t", [-0.0, 1.0])])
    assert spec_digest(a) != spec_digest(b)


def test_both_sigma_forms_are_refused():
    # The digest hashes exactly one sigma marker; both set would leave
    # unpinned bytes in the file (review finding).
    import pytest

    s = spec(
        [axis("x", [0.0, 1.0])],
        [value("t", [1.0, 2.0], sigma=(0.1, 0.2), sigma_scalar=0.1)],
    )
    with pytest.raises(ValueError, match="both per-point sigma"):
        spec_digest(s)


def test_shape_mismatches_are_refused_at_digest_time():
    import pytest

    wrong_len = spec([axis("x", [0.0, 1.0])], [value("t", [1.0, 2.0, 3.0])])
    with pytest.raises(ValueError, match="axes imply"):
        spec_digest(wrong_len)
    short_sigma = spec([axis("x", [0.0, 1.0])], [value("t", [1.0, 2.0], sigma=(0.1,))])
    with pytest.raises(ValueError, match="sigma has 1 elements"):
        spec_digest(short_sigma)
