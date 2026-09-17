"""h5py writer round-trip: the file carries exactly the reader.rs schema and
the returned digest is stable across writes (Tier-2 determinism)."""

import h5py
import numpy as np
import pytest

from crucible_offl.seam_fixture import GROUP_PATH, fixture_spec
from crucible_offl.tables import Axis, TableValue, write_table


def test_writer_layout_matches_reader_schema(tmp_path):
    path = str(tmp_path / "seam.h5")
    write_table(path, GROUP_PATH, fixture_spec())
    with h5py.File(path, "r") as f:
        g = f[GROUP_PATH]
        assert set(g.keys()) == {"axes", "values"}
        assert g.attrs["kind"] == "regular"
        assert g.attrs["schema_version"] == "1.0"
        assert g.attrs["axis_order"] == "p\nh"
        assert g.attrs["source_library"] == "Glenn-α (UTF-8 probe)"
        assert g.attrs["rng_seed"] == 20260817
        assert g.attrs["rng_seed"].dtype == np.int64
        assert list(g["axes"].keys()) == ["h", "p"]  # h5py sorts; both present
        p = g["axes/p"]
        assert p.dtype == np.float64
        assert p.attrs["envelope_min"] == 2.0e4
        vals = g["values"]
        assert set(vals.keys()) == {
            "temperature",
            "sigma_temperature",
            "gamma_eff",
            "sound_speed",
        }
        t = vals["temperature"]
        assert t.attrs["units"] == "K"
        assert t.attrs["interp_rule"] == "log-lin-lin"
        assert t.attrs["interp_error_bound"] == 0.5
        assert "sigma" in vals["gamma_eff"].attrs
        assert "sigma" not in t.attrs
        assert vals["sound_speed"].shape == (12,)


def test_write_is_deterministic(tmp_path):
    d1 = write_table(str(tmp_path / "a.h5"), GROUP_PATH, fixture_spec())
    d2 = write_table(str(tmp_path / "b.h5"), GROUP_PATH, fixture_spec())
    assert d1 == d2


def test_writer_refuses_reserved_sigma_prefix(tmp_path):
    spec = fixture_spec()
    bad = spec.__class__(
        **{
            **spec.__dict__,
            "values": (
                TableValue("sigma_t", (1.0,) * 12, "K", "log-lin-lin", 0.1),
            ),
        }
    )
    with pytest.raises(ValueError, match="reserved"):
        write_table(str(tmp_path / "bad.h5"), GROUP_PATH, bad)


def test_writer_refuses_newline_axis_name(tmp_path):
    spec = fixture_spec()
    bad = spec.__class__(
        **{
            **spec.__dict__,
            "axes": (Axis("p\nq", (1.0, 2.0), 1.0, 2.0),),
        }
    )
    with pytest.raises(ValueError, match="newline"):
        write_table(str(tmp_path / "bad.h5"), GROUP_PATH, bad)
