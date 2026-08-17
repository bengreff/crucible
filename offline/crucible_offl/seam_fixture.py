"""The cross-language seam fixture: a table written by *this* package that
the Rust loader must open and digest-verify. One spec, three assertions:

- ``tests/test_seam_fixture.py`` regenerates it and asserts ``PIN_DIGEST``;
- ``scripts/make_seam_fixture.py`` writes the committed copy at
  ``crates/tables/tests/fixtures/python_seam.h5``;
- ``crates/tables/tests/fnd5_python_seam.rs`` opens that copy with the same
  pin — if h5py and the Rust reader ever disagree about bytes, attributes,
  or the digest stream, that test halts the battery.

The spec deliberately exercises every schema feature at once: a log-spaced
axis, an envelope strictly inside the grid, per-point sigma, scalar sigma,
no sigma, an rng_seed, and a non-ASCII provenance string (UTF-8 byte-length
prefixing, not char count).
"""

from __future__ import annotations

from .tables import Axis, Provenance, TableValue, WriteSpec

GROUP_PATH = "/chem/seam_probe"
DATA_VERSION = "0.1.0"
# Pinned on first generation (2026-08-17); any change to the writer or the
# digest that moves this is a seam break and must bump both sides together.
PIN_DIGEST = "sha256:221f4c4014b78664d844f74859dcd38980dbc624d89c94f4428ad7340b3bbb3c"


def fixture_spec() -> WriteSpec:
    p = (1.0e4, 1.0e5, 1.0e6, 1.0e7)  # Pa, log-spaced
    hgrid = (-2.0e6, 1.0e6, 4.0e6)  # J/kg, linear
    n = len(p) * len(hgrid)
    # Row-major over (p, h): index = i_p * len(h) + i_h.
    temperature = tuple(300.0 + 10.0 * i_p + 1.5 * i_h for i_p in range(4) for i_h in range(3))
    gamma_eff = tuple(1.4 - 0.01 * i for i in range(n))
    sound_speed = tuple(340.0 * (1.0 + 0.05 * i) for i in range(n))
    return WriteSpec(
        kind="regular",
        schema_version="1.0",
        data_version=DATA_VERSION,
        interp_method="multilinear",
        provenance=Provenance(
            producer="crucible-offl/seam-fixture",
            producer_version="0.1.0",
            input_deck_hash="sha256:seam-fixture-has-no-deck",
            source_library="Glenn-α (UTF-8 probe)",
            generator_commit="worktree",
            rng_seed=20260817,
        ),
        axes=(
            Axis("p", p, envelope_min=2.0e4, envelope_max=1.0e7),
            Axis("h", hgrid, envelope_min=-2.0e6, envelope_max=3.5e6),
        ),
        values=(
            TableValue(
                "temperature",
                temperature,
                units="K",
                interp_rule="log-lin-lin",
                interp_error_bound=0.5,
                sigma=tuple(0.01 * t for t in temperature),
            ),
            TableValue(
                "gamma_eff",
                gamma_eff,
                units="1",
                interp_rule="log-lin-lin",
                interp_error_bound=1e-3,
                sigma_scalar=0.005,
            ),
            TableValue(
                "sound_speed",
                sound_speed,
                units="m/s",
                interp_rule="log-lin-log",
                interp_error_bound=0.1,
            ),
        ),
    )
