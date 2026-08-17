"""The Python half of the seam-fixture pin: regenerating the fixture spec
must reproduce PIN_DIGEST exactly. The Rust half
(crates/tables/tests/fnd5_python_seam.rs) opens the committed file with the
same pin — together they prove h5py-written bytes and the Rust reader agree
about every field the digest covers."""

from crucible_offl.seam_fixture import GROUP_PATH, PIN_DIGEST, fixture_spec
from crucible_offl.tables import spec_digest, write_table


def test_fixture_spec_digest_matches_pin():
    assert spec_digest(fixture_spec()) == PIN_DIGEST


def test_written_fixture_stamps_the_pin(tmp_path):
    assert write_table(str(tmp_path / "seam.h5"), GROUP_PATH, fixture_spec()) == PIN_DIGEST
