#!/usr/bin/env bash
# VAL-3 §3.2 per-commit gate battery (fast tier), run in fixed order.
# Any red gate breaks the build; every session ends with this green (VISION_SCOPE §12).
#
# What's armed lives in the test files themselves (names cite their doc §);
# certificates/ carries the certified claims. Not yet armed (arrive with
# their subject): full COUP-2 port-accounting audit; 1-vs-N-thread byte
# identity (first parallel sweep); Su-Olson anchor (SOLV-2).
# Battery budget: the station-2 ladder and station-4 coupled duct are the
# heavy items (~1 min each in gate 3; ~1 min + ~15 s in gate 5); next
# growth splits to the VAL-3 §3.2 milestone tier.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== gate 1/5: format (cargo fmt --check)"
cargo fmt --all --check

echo "== gate 2/5: lints (cargo clippy -D warnings)"
cargo clippy --workspace --all-targets -- -D warnings

echo "== gate 3/5: tests (cargo test)"
cargo test --workspace --quiet

echo "== gate 4/5: offline seam + OFFL-3 battery (pytest)"
# The Python half of the two-language rule (VISION_SCOPE §6): digest v2
# golden vector, h5py writer/fixture, RP-1311 reproduction, CEA↔Cantera
# cross-check, production-table holdout bounds, coordinate consistency.
if [ ! -x offline/.venv/bin/python ]; then
  echo "offline/.venv missing — bootstrap: /usr/local/bin/python3.13 -m venv offline/.venv \\"
  echo "  && offline/.venv/bin/pip install -e 'offline[test]'" >&2
  exit 1
fi
(cd offline && .venv/bin/python -m pytest -q)

echo "== gate 5/5: certificate artifacts match regeneration (no silent drift)"
# Release profile: the Sod ladder is a real solver march (~6 s optimized,
# ~20× that unoptimized). Optimization does not change f64 results in Rust
# (no fast-math — META-1 §2.2), so dev/test/release agree bitwise.
cargo run -q --release --bin convergence_certificate >/dev/null
cargo run -q --release --bin station1_sod_certificate >/dev/null
cargo run -q --release --bin station2_nozzle_certificate >/dev/null
offline/.venv/bin/python offline/scripts/station3_flame_certificate.py >/dev/null
cargo run -q --release --bin station4_cooled_wall_certificate >/dev/null
git diff --exit-code certificates/
# A brand-new (untracked) certificate is invisible to git diff — refuse it
# too, or a first regeneration ships uncommitted (review finding).
if [ -n "$(git ls-files --others --exclude-standard certificates/)" ]; then
  echo "untracked certificate artifact(s) present — commit them:" >&2
  git ls-files --others --exclude-standard certificates/ >&2
  exit 1
fi

echo "== all gates green"
