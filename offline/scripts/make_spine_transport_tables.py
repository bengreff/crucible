"""Generate the spine's chemical-regime transport surface (OFFL-5 §3.1a).

Usage (from repo root):
    offline/.venv/bin/python offline/scripts/make_spine_transport_tables.py

Writes tables/spine/lox_lh2_transport_v<DATA_VERSION>.h5 + the pins
sidecar. Takes ~10 s and prints CEA's "Singular update matrix" chatter from
the constant-temperature Z-difference solves — those solves are checked for
convergence and the derivative is step-independent to ~1e-4 relative
(OFFL-5 §6), so the chatter is noise, not a fault. Suppressing a solver's
warnings to make output tidy is exactly the habit META-1 P6 forbids.

**Grid density is a measured choice.** 57 x 37 x 7 puts every column's
rule-space interpolation bound at or below ~3.3% against the declared
10-20% physical band (FND-7 §3.3) — a factor >= 3 of margin, which is what
FND-5 §3.4 asks for. The next step up (85 x 55 x 9) halves the bound for
~3x the nodes and generation time; that is the refinement to reach for if a
consumer ever needs it, not a default.
"""

import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from crucible_offl.chemistry import EquilibriumEngine
from crucible_offl.surfaces import station5_envelope_grid
from crucible_offl.transport import matched_transport_grid, write_transport_table

DATA_VERSION = "0.1.0"
#: Node counts (module doc). The grid itself is DERIVED from the
#: equilibrium surface's declared envelope so the two cannot drift.
N_P, N_H, N_Z = 57, 37, 7


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
    out = repo / "tables" / "spine" / f"lox_lh2_transport_v{DATA_VERSION}.h5"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.unlink(missing_ok=True)
    grid = matched_transport_grid(
        station5_envelope_grid(), n_p=N_P, n_h=N_H, n_z=N_Z
    )
    digests = write_transport_table(
        str(out),
        DATA_VERSION,
        commit,
        grid,
        engine=EquilibriumEngine(gas_only=True),
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
