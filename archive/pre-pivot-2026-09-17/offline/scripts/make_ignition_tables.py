"""Generate the ignition/flame closure surface — S_L(p, T_u, Z) and
tau_ign(p, T_u, Z) for the SOLV-4 §3.6 burn-progress rate law (OFFL-3 §3.3,
plan S6).

Usage (from repo root):
    offline/.venv/bin/python offline/scripts/make_ignition_tables.py

Writes tables/chem/lox_lh2_ignition_v<DATA_VERSION>.h5 + the pins sidecar.
Takes a few minutes (each grid node is one Cantera FreeFlame solve + one 0-D
reactor integration on h2o2.yaml). The two columns share the (p, T_u, Z)
axes:

  * laminar_flame_speed [m/s], LINEAR-valued — carries 0 at the flammability
    limit (a non-lighting FreeFlame, or a runaway auto-ignitive solve above
    the physical ceiling, ⇒ 0). This is the propagation term's S_L.
  * ignition_delay [s], LOG-valued and capped at TAU_MAX (10 s) — a reactor
    that does not ignite by the horizon returns the cap, so tau -> "very large"
    smoothly, never inf. This is the auto-ignition term's tau_ign.

Extinction is thus the columns' own values (Rule 12: no runtime branch). The
Cantera pin is part of the record (data_version event on a pin bump); the
provenance deck stamps cantera {version} + the mechanism.

Determinism: FreeFlame + reactor are bit-reproducible cross-process given the
declared fixed solver settings (see crucible_offl/ignition.py); the regen
probe (offline/tests/test_regen_determinism.py) pins it.
"""

import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from crucible_offl.ignition import IgnitionEngine, ignition_grid, write_ignition_table

DATA_VERSION = "0.3.0"


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
    out = repo / "tables" / "chem" / f"lox_lh2_ignition_v{DATA_VERSION}.h5"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.unlink(missing_ok=True)
    eng = IgnitionEngine()
    digests = write_ignition_table(
        str(out), DATA_VERSION, commit, engine=eng, grid=ignition_grid()
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
