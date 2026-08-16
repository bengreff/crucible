#!/usr/bin/env bash
# VAL-3 §3.2 per-commit gate battery (fast tier), run in fixed order.
# Any red gate breaks the build; every session ends with this green (VISION_SCOPE §12).
#
# Armed inside gate 3 (test names cite their doc §): MMS order-of-accuracy
# (cert_mms_*), envelope/fail-loud (fnd5_s6_4, review_* regressions),
# analytic anchors (cert_annulus_*, cert_bessel_*), closed-sweep
# conservation (cert_closed_sweep_*), determinism/rerun byte-identity
# (cert_config_drives_*, fnd2_s6_3, fnd5_s6_6).
# Not yet armed (arrive with their subject): full COUP-2 port-accounting
# audit; 1-vs-N-thread byte identity (first parallel sweep); Sod/Su-Olson
# anchors (SOLV-1/SOLV-2).
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== gate 1/4: format (cargo fmt --check)"
cargo fmt --all --check

echo "== gate 2/4: lints (cargo clippy -D warnings)"
cargo clippy --workspace --all-targets -- -D warnings

echo "== gate 3/4: tests (cargo test)"
cargo test --workspace --quiet

echo "== gate 4/4: certificate artifact matches regeneration (no silent drift)"
cargo run -q --bin convergence_certificate >/dev/null
git diff --exit-code certificates/

echo "== all gates green"
