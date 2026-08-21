"""Generate the unburnt-reactant surface — the burn-progress c = 0 branch
(OFFL-3 §3.3, plan S5).

Usage (from repo root):
    offline/.venv/bin/python offline/scripts/make_unburnt_tables.py

Writes tables/chem/lox_lh2_unburnt_v<DATA_VERSION>.h5 + the pins sidecar.
Takes ~4 s. The surface is the **gas-phase ideal-gas frozen reactant
mixture** (gaseous H2 + O2 at the elemental proportions Z), on the SAME CEA
enthalpy reference as the equilibrium surface so the SOLV-4 §3.6 blend
`h = (1-c)h_u + c*h_b` adds on one reference. It uses the equilibrium
surface's FND-5 schema (density/sound_speed/temperature/gamma_eff/mbar), so
`TableEos` binds it with no new occupant (SOLV-1 §3.4).

**Grid density is a measured choice.** 9 x 121 x 15 puts every column's
holdout bound well inside any reasonable model band: temperature 1.66 K abs
(~0.08% hot, ~1.7% at the 100 K cold floor), density 0.81% rule-space,
sound speed 0.90 m/s, gamma_eff 1.2e-3. The p-axis is coarse on purpose —
the ideal-gas density is exactly log-linear in p and the other columns are
p-invariant, so p resolution costs nothing and buys nothing beyond bracketing
the envelope. The h-axis carries the T(h) curvature (worst at the H2-rich
hot corner, where H2's vibrational modes activate), the Z-axis the M-bar and
mixture-cp curvature — both measured, refine there if a consumer needs it.
"""

import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from crucible_offl.chemistry import FrozenReactantEngine
from crucible_offl.surfaces import unburnt_reactant_grid, write_unburnt_table

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
    out = repo / "tables" / "chem" / f"lox_lh2_unburnt_v{DATA_VERSION}.h5"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.unlink(missing_ok=True)
    eng = FrozenReactantEngine()
    digests = write_unburnt_table(
        str(out), DATA_VERSION, commit, engine=eng, grid=unburnt_reactant_grid(eng)
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
