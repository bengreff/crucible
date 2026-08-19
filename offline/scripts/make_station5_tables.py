"""Generate the station-5 production tables (OFFL-3, LOX/LH2 wide-envelope).

Usage (from repo root):
    offline/.venv/bin/python offline/scripts/make_station5_tables.py

Same pipeline as make_station3_tables.py with the station-5 envelope grid
(pressure axis widened down to the vacuum-plume fringe — an envelope
setting, OFFL-3 §3.3): writes tables/chem/lox_lh2_v<DATA_VERSION>.h5 + the
pins sidecar. The station-3 v0.1.0 artifact and its certificate are
untouched (a new data_version is a new artifact, never a mutation).
"""

import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from crucible_offl.surfaces import station5_envelope_grid, write_station3_tables

DATA_VERSION = "0.3.2"


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True, check=True
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "status", "--porcelain"], cwd=repo, capture_output=True, text=True, check=True
    ).stdout.strip()
    if dirty:
        commit += "-dirty"  # git-describe convention: data hygiene over vanity
    out = repo / "tables" / "chem" / f"lox_lh2_v{DATA_VERSION}.h5"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.unlink(missing_ok=True)
    digests = write_station3_tables(
        str(out), DATA_VERSION, commit, eq_grid=station5_envelope_grid(), gas_only=True
    )
    pins = out.with_suffix(".pins.toml")
    with open(pins, "w") as f:
        f.write(
            "# Table pins (FND-5 §3.2), consumed by the FND-4 §6-4 [tables] grammar\n"
            "# (the machine-written single pin owner).\n"
        )
        for group, digest in sorted(digests.items()):
            f.write(f'\n["{group}"]\ndata_version = "{DATA_VERSION}"\ncontent_digest = "{digest}"\n')
    print(f"wrote {out}")
    for group, digest in sorted(digests.items()):
        print(f"  {group}: {digest}")
    print(f"pins recorded in {pins}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
