# CRUCIBLE Station 3 Certificate — the Flame Seam

Goal-B station 3 (OFFL-3 §6; FND-5 §3.2; VISION_SCOPE §6 two-language rule): the LOX/LH2 flame is computed **offline** by NASA CEA (equilibrium free-energy minimization over the full Glenn species set), tabulated as local-state surfaces, and crosses the project's one cross-language seam as versioned, digest-pinned HDF5 that the Rust runtime loads, verifies, and interpolates. The physical system: the RL10 combustion chamber at p_c = 32.75 bar burning liquid hydrogen (20.27 K) with liquid oxygen (90.17 K) at MR = 5 — the flame the station-5 blind run will feed on. Python never runs at simulation time; the tables are the only crossing.

Criteria are CI-enforced in `offline/tests/` (pytest: digest golden vector, RP-1311 reproduction, CEA↔Cantera cross-check, frozen/shifting bracket, holdout error bounds, regeneration determinism, coordinate consistency) and `crates/tables/tests/` (`fnd5_python_seam.rs`: the h5py-written fixture opens under its Python-stamped pin; `station3_tables.rs`: the production tables load, interpolate identically to the Python reference evaluator, and refuse out-of-envelope). Regenerate: `offline/scripts/station3_flame_certificate.py`.

## The seam itself (digest v2, byte-for-byte)

- Golden vector (both languages assert the same constant): `sha256:a3b0bc5987ad11ba166d006ab7103c16ff1d3a8f17f64a248d4b11d8091023f2`
- Cross-language fixture pin (Python stamps, Rust verifies): `sha256:5f91b7ab6cbd34c836fa214175553de5f18cee3c6842837443e0f4ace49a78b4`
- Production table pins (`tables/chem/lox_lh2_v0.1.0.pins.toml`):
  - `/chem/lox_lh2/equilibrium`: `sha256:8cc3f8b8d810745d9f6ef39671b0dd4acccac974c06a4a99d873e5221093d6b3`
  - `/chem/lox_lh2/performance`: `sha256:aafc17e702383b8c4a96189358d0f525cbb209fbec80340de9fc8b3afe2a1d7d`
- Toolchain pins: Python 3.13.7, cea 3.3.2 (libcea 3.3.2), cantera 3.2.0, h5py 3.16.0, numpy 2.5.2, crucible-offl 0.1.0; generator commit `0847979f4d31`.

## OFFL-3 §6-1 — RP-1311 example 8 reproduced through the pipeline

LOX/LH2 IAC rocket, p_c = 53.3172 bar, o/f = 5.55157 (McBride & Gordon, RP-1311, NTRS 19960044559). Reference = the published example output; tolerance = its printed precision (the pytest gates).

| quantity | published | pipeline | rel. dev |
|---|---|---|---|
| T_c [K] | 3383.845 | 3383.84 | 1.1e-07 |
| rho_c [kg/m^3] | 2.41 | 2.40968 | 1.3e-04 |
| h_c [kJ/kg] | -1026.05 | -1026.05 | 4.8e-06 |
| s_c [kJ/kg-K] | 18.659 | 18.6586 | 2.3e-05 |
| Mbar [kg/kmol] | 12.716 | 12.7157 | 2.7e-05 |
| gamma_s | 1.145 | 1.1447 | 2.7e-04 |
| a_c [m/s] | 1591.47 | 1591.47 | 1.6e-06 |
| X_H2O | 0.63456 | 0.634556 | 6.0e-06 |
| X_H2 | 0.29479 | 0.294794 | 1.5e-05 |
| X_H | 0.033498 | 0.0334979 | 3.4e-06 |
| X_OH | 0.033341 | 0.0333414 | 1.2e-05 |
| c* [m/s] | 2332.34 | 2332.34 | 1.8e-06 |
| T_throat [K] | 3185.673 | 3185.67 | 4.8e-08 |
| Isp_vac(eps=25) [m/s] | 4348.51 | 4348.51 | 8.8e-07 |

## OFFL-3 §6-3 — frozen/shifting bracket (the model-form error bar)

Vacuum Isp at area ratio 25: shifting 4348.5 m/s, frozen-from-chamber 4178.7 m/s — a strict bracket of width 3.91% that contains the JANNAF kinetic-efficiency knockdown (~0.8–1% of shifting for LOX/LH2 at MR≈5, META-3 `jannaf-eff`). Real delivered performance sits between the two tables; COUP-5 samples this as a discrete epistemic dimension (S19) — the disagreement **is** the error bar, never averaged away.

## OFFL-3 §6-2 — two independent solvers on the chamber state

CEA (Glenn NASA9 fits) vs Cantera 3.2 `h2o2.yaml` (NASA7 fits), no shared implementation, same (p, h, Z) coordinate:

| p_c [bar] | MR | T: CEA [K] | T: Cantera [K] | dT | dM̄ | dX_H2O |
|---|---|---|---|---|---|---|
| 32.75 | 5.0 | 3225.41 | 3228.91 | +0.108% | +0.036% | +0.00130 |
| 53.32 | 5.55157 | 3383.84 | 3388.17 | +0.128% | +0.045% | +0.00167 |
| 10.00 | 4.0 | 2865.57 | 2867.22 | +0.057% | +0.018% | +0.00060 |

Gates (2× the measured band): |dT| ≤ 0.3%, |dM̄| ≤ 0.15%, |dX| ≤ 5e-3, |dcp_eq| ≤ 1%. The residual disagreement is thermodynamic-data uncertainty (§3.4), carried, not hidden.

## The tables (FND-5 schema, committed at `tables/chem/`)

**Equilibrium surface** `(p, h, Z) → T, ρ, γ_eff, a, M̄, condensed_fraction, X_k` (41×31×13, log-p): the S22 local-state coordinate of SOLV-1's shifting mode. **Performance reference** `(p_c, MR) → c*_ideal, T_c, γ, M̄` (11×11): SOLV-7's anchor and the SOLV-1 §3.4 knockdown reference. Declared envelopes (refusal-enforced, grid strictly wider):

- `p`: grid [1500, 8e+06], envelope [2000, 7e+06]
- `h`: grid [-1.18e+07, -100000], envelope [-1.15e+07, -200000]
- `Z`: grid [0.1, 0.26], envelope [0.111111, 0.25]

Stored `interp_error_bound` per column = **measured** holdout max (midpoints + ¼-offsets, direct CEA solves) × 1.5 declared sampling margin (FND-5 §3.4); verified in CI by a fresh disjoint ⅜-offset sweep against the committed artifact — a violation means *refine the grid*, never *relax the gate*. Full-envelope bounds are dominated by the H₂O condensation kink in the deep-cold corners (reported honestly by the `condensed_fraction` column; the certificate's gas-region numbers below show the working regime is ~50× tighter):

| column | stored bound | | column | stored bound |
|---|---|---|---|---|
| X_H | 0.00149 | condensed_fraction | 0.0306 |
| X_H2 | 0.00357 | density | 0.36 |
| X_H2O | 0.03 | gamma_eff | 0.119 |
| X_O | 0.000316 | mbar | 0.277 |
| X_O2 | 0.00155 | sound_speed | 49.3 |
| X_OH | 0.0017 | temperature | 36.5 |

Performance reference bounds: T_c 10.4, c_star_ideal 2.69, gamma 0.000956, mbar 0.0135.

## OFFL-3 §6-6 — one flame, two parameterizations (S22)

Along the design line h = h_inj(Z), the (p, h, Z) surface and the (p_c, MR) performance reference must describe the same chamber. Direct CEA solution alongside:

| p_c [bar] | MR | T surface [K] | T perf-ref [K] | T direct [K] | surf−direct |
|---|---|---|---|---|---|
| 32.75 | 5.0 | 3225.00 | 3225.15 | 3225.41 | -0.42 |
| 20.00 | 4.5 | 3059.47 | 3060.19 | 3060.19 | -0.72 |
| 45.00 | 6.0 | 3433.21 | 3435.56 | 3435.56 | -2.36 |
| 15.00 | 3.8 | 2811.66 | 2806.18 | 2811.81 | -0.15 |
| 50.00 | 7.0 | 3529.42 | 3531.93 | 3531.93 | -2.51 |

## The RL10 chamber through the whole seam (the station-3 headline)

At (p_c = 32.75 bar, h = h_inj, Z = 1/6): the Rust loader verifies the pin, interpolates T = 3225.00 K against CEA's direct 3225.41 K (-0.42 K in the gas region, vs the kink-dominated stored bound ±36.5 K), and c*_ideal = 2362.48 m/s against the direct rocket solve 2362.53 m/s (-0.05 m/s). The same numbers are asserted in `station3_tables.rs` — Python and Rust provably interpolate the same surface the same way (≤1e-6 abs).

## Deferred, loudly (owners)

- **Frozen-path surface vs (p, h, {X_k})** — its axis set belongs to SOLV-1's frozen-advection consumer; arrives with that wave (the frozen *bracket* is validated above; `chemistry.py` header).
- **Quasi-1-D expansion oracles** (`oracle`-labeled, never field-interpolated) — SOLV-7/VAL-2 C_F cross-check wave, station 5.
- **B′ ablation tables** — SOLV-8 wave. The **transport feed is DISCHARGED** (plan S4): `chemistry.py` now returns the caloric companions (c_p,fr, c_v,eq) from the same CEA solve, and `crucible_offl/transport.py` assembles them with Cantera's mixture-averaged transport into the spine's chemical-regime surface (OFFL-5 §3.1a, `tables/spine/lox_lh2_transport_v0.1.0.h5`) — written in the runtime local-state (p, h, Z) coordinate, not the (T, p, Z) feed one. OFFL-3 still ships no runtime transport table; the spine (FND-7) is the sole runtime provider.
- **Config→tables pin wiring** (FND-4 §6-4) — the `[tables]` grammar consumes `lox_lh2_v0.1.0.pins.toml` when it lands; until then the sidecar + the Rust test pins are the record.
- **Per-point sigma columns** — the thermodynamic-data band enters with the COUP-5 UQ wave; the seam already carries sigma companions (proven by the fixture).

Honest scaffolding note: the surfaces are *the equilibrium*, including condensed H₂O where the rectangular grid's deep-cold corners demand it (`condensed_fraction` column); no engine trajectory enters that region, and the gas-only working regime is where all gas-region numbers above live.
