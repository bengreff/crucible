#!/usr/bin/env bash
# VAL-3 §3.2 per-commit gate battery (fast tier), run in fixed order.
# Any red gate breaks the build; every session ends with this green (VISION_SCOPE §12).
#
# Armed inside gate 3 (test names cite their doc §): MMS order-of-accuracy
# (cert_mms_*, solv1_s62_mms_*), envelope/fail-loud (fnd5_s6_4,
# station1_gamma_*, review_* regressions), analytic anchors (cert_annulus_*,
# cert_bessel_*, the station1_* Sod battery vs the exact Riemann oracle, the
# station2_* choked-nozzle battery vs isentropic theory), closed-sweep
# conservation (cert_closed_sweep_*, station1_closed_tube_*),
# determinism/rerun byte-identity (cert_config_drives_*, station1_sod_rerun_*,
# station2_rerun_*, fnd2_s6_3, fnd5_s6_6).
# Not yet armed (arrive with their subject): full COUP-2 port-accounting
# audit; 1-vs-N-thread byte identity (first parallel sweep); Su-Olson
# anchor (SOLV-2).
# Battery budget note: the station-2 ladder is the heavy item (~1 min in
# gate 3 + ~1 min in gate 4); if it grows further, split it to the VAL-3
# §3.2 milestone tier.
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== gate 1/4: format (cargo fmt --check)"
cargo fmt --all --check

echo "== gate 2/4: lints (cargo clippy -D warnings)"
cargo clippy --workspace --all-targets -- -D warnings

echo "== gate 3/4: tests (cargo test)"
cargo test --workspace --quiet

echo "== gate 4/4: certificate artifacts match regeneration (no silent drift)"
# Release profile: the Sod ladder is a real solver march (~6 s optimized,
# ~20× that unoptimized). Optimization does not change f64 results in Rust
# (no fast-math — META-1 §2.2), so dev/test/release agree bitwise.
cargo run -q --release --bin convergence_certificate >/dev/null
cargo run -q --release --bin station1_sod_certificate >/dev/null
cargo run -q --release --bin station2_nozzle_certificate >/dev/null
git diff --exit-code certificates/

echo "== all gates green"
