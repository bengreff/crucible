"""Regenerate the committed cross-language seam fixture.

Usage (from repo root):
    offline/.venv/bin/python offline/scripts/make_seam_fixture.py

Writes crates/tables/tests/fixtures/python_seam.h5 and prints the digest,
which must equal crucible_offl.seam_fixture.PIN_DIGEST (and the pin in
crates/tables/tests/fnd5_python_seam.rs). If it doesn't, the seam moved:
update both pins in the same commit and say why.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from crucible_offl.seam_fixture import DATA_VERSION, GROUP_PATH, PIN_DIGEST, fixture_spec
from crucible_offl.tables import write_table


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    out = repo / "crates" / "tables" / "tests" / "fixtures" / "python_seam.h5"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.unlink(missing_ok=True)
    digest = write_table(str(out), GROUP_PATH, fixture_spec())
    print(f"wrote {out}")
    print(f"data_version   {DATA_VERSION}")
    print(f"content_digest {digest}")
    if digest != PIN_DIGEST:
        print(f"MISMATCH: pinned {PIN_DIGEST}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
