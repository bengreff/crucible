#!/usr/bin/env bash
# VAL-3 §3.2 per-commit gate battery (fast tier), run in fixed order.
# Any red gate breaks the build; every session ends with this green (VISION_SCOPE §12).
#
# Gates not yet armed (placeholders — armed when their first subject lands):
#   - MMS order-of-accuracy        (first PDE operator, FND-2/SOLV-1)
#   - conservation audit           (COUP-2 harness)
#   - determinism 1-vs-N threads   (first parallel numeric path)
#   - envelope / fail-loud         (FND-5 loader)
#   - analytic anchors             (Sod, Su-Olson, … — VAL-2)
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== gate 1/3: format (cargo fmt --check)"
cargo fmt --all --check

echo "== gate 2/3: lints (cargo clippy -D warnings)"
cargo clippy --workspace --all-targets -- -D warnings

echo "== gate 3/3: tests (cargo test)"
cargo test --workspace --quiet

echo "== all gates green"
