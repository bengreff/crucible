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
# station2_rerun_*, station4_rerun_*, fnd2_s6_3, fnd5_s6_6), cross-language
# seam (fnd5_python_seam, station3_tables vs the Python-stamped pins),
# conjugate-wall battery (station4_*: series-resistance oracle, stair-
# interface round-off ledger, Robin annulus anchor, coupled fixed point).
# Armed inside gate 4 (offline pytest): digest v2 golden vector + fixture
# pin, RP-1311 example-8 reproduction, CEA↔Cantera cross-check,
# frozen/shifting bracket, production-table fresh-holdout bounds,
# regeneration digest-identity, (p,h,Z)↔(p_c,MR) coordinate consistency.
# Not yet armed (arrive with their subject): full COUP-2 port-accounting
# audit; 1-vs-N-thread byte identity (first parallel sweep); Su-Olson
# anchor (SOLV-2).
# Battery budget note: the station-2 ladder and the station-4 coupled duct
# are the heavy items (~1 min each in gate 3; station-2 ~1 min + station-4
# ~15 s in gate 5); next growth splits to the VAL-3 §3.2 milestone tier.
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

echo "== all gates green"
