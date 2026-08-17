"""Generate the station-3 production tables (OFFL-3, LOX/LH2 design window).

Usage (from repo root):
    offline/.venv/bin/python offline/scripts/make_station3_tables.py

Writes tables/chem/lox_lh2_v<DATA_VERSION>.h5 (two groups: the (p, h, Z)
equilibrium surface and the (p_c, MR) performance reference) plus a pins
sidecar TOML, and prints the digests. Regeneration on the same machine and
tool pins must reproduce the digests bit-for-bit (Tier-2 determinism,
META-1 §2.1) — CI-asserted by test_station3_surfaces.py on a subgrid.
"""

import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from crucible_offl.surfaces import write_station3_tables

DATA_VERSION = "0.1.0"


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
    digests = write_station3_tables(str(out), DATA_VERSION, commit)
    pins = out.with_suffix(".pins.toml")
    with open(pins, "w") as f:
        f.write(
            "# Table pins (FND-5 §3.2): the config-side [tables] wiring (FND-4 §6-4)\n"
            "# consumes these when it lands; until then this sidecar is the record.\n"
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
