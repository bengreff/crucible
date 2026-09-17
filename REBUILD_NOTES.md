# Lessons and references for the rebuild

17 September 2026. This is a selective extraction from the old project, not a renewed implementation plan. `VISION_SCOPE.md` is authoritative. Archive-relative references below point into `archive/pre-pivot-2026-09-17/`; their contents were preserved unchanged. Line numbers refer to that frozen snapshot.

The extraction inspected repository structure, selected numerical routines, certificate machinery, and relevant session records. It did not independently rerun or validate the legacy solver. Statements about historical results below are attributed to those records.

**Lessons worth carrying forward**

1. **Build the human experiment loop early.** The previous plan explicitly deferred visualization until a large certification run. That conflicts with the new product. The new instrument must let a person see the arrangement, observe the calculated state, ask where energy goes, and compare a change. An inaccessible solver is an incomplete sandbox. This is a product decision from the pivot, not a measured failure of an interface that did not yet exist.

2. **Keep the physical scope distinct from resolution.** The reacting chamber belongs inside the simulator. Turbulence and chemistry can be represented by responsive models. This does not justify assuming away their influence, nor does inclusion require resolving every eddy or collision. Avoid inheriting the old escalation from a physical question to mandatory three-dimensional startup simulation.

3. **A good headline number can hide the wrong internal behavior.** The [session record](archive/pre-pivot-2026-09-17/SESSION_LOG.md#L703) reports that adding/fixing transport changed scored performance by at most about 0.32% in one coarse study while wall heat changed substantially. An apparent 3.27% heat change became 19.8% after correcting a driving potential; the earlier small difference involved cancellation. Inspect fields, local budgets, heat loads, and intermediate quantities, not only thrust/Isp. This is a reason visualization belongs in the scientific workflow.

4. **Identify exactly what a certificate proves.** The [RL10 certificate generator](archive/pre-pivot-2026-09-17/crates/engine/src/bin/station5_rl10_certificate.rs#L1) explicitly recomputes scoring from recorded readouts rather than rerunning the underlying flow. Its constants include different model assumptions and calibrated versus prior-band cases. A green certificate diff establishes that score generation is unchanged; it does not establish that the current executable reproduces the run or that combustion startup is validated.

5. **A numerically completed run is not a physical feasibility verdict.** The [later GPU startup record](archive/pre-pivot-2026-09-17/SESSION_LOG.md#L2240) reports a conservation-audited trajectory that ignited a small sustained kernel but reached only about 22.5 kN against a 75.6 kN target. The log attributes the shortfall to the flame-spreading/modeling gap at that setup. Preserve this as a limitation of that run and model, not evidence that the real engine cannot work. The new interface must distinguish numerical failure, an unmet target, and missing physical capability.

6. **Different processes settle on different clocks.** The [wall-response discussion](archive/pre-pivot-2026-09-17/SESSION_LOG.md#L719) compares an approximately 37 ms wall thermal timescale with an 11 ms run. Settled gas did not mean a settled wall. Show what has equilibrated, what is still changing, and why a run ended.

7. **Geometry and conservation need shared definitions.** [Grid metrics](archive/pre-pivot-2026-09-17/crates/grid/src/lib.rs#L1026) assign each shared face one radius calculation. The comments record that algebraically equivalent radius formulas rounded differently and broke exact cancellation. Preserve this ownership principle in axisymmetric volumes, face areas, source terms, and diagnostics; do not inherit the entire adaptive-azimuthal grid just to keep it.

8. **Table validity and thermodynamic consistency are part of the physics.** Keep explicit units, data provenance, interpolation coordinates, and applicability checks. The [old blend model](archive/pre-pivot-2026-09-17/crates/solvers/src/euler/blend_eos.rs#L1) contains branch thresholds, pressure-root selection, and warm-start behavior that deserve a fresh physical review before reuse. Continuous interpolation or selecting the first root is not, by itself, proof that a closure describes the new medium correctly.

9. **Measure the expensive operation.** The [GPU profiling record](archive/pre-pivot-2026-09-17/SESSION_LOG.md#L2263) reports that reaction-node projections consumed 79% of kernel time in that case; flux sweeps were about 5%. Reusing pressure hints reduced the measured composed step from roughly 0.137 to 0.028 seconds. These are historical case-specific timings. The reusable lesson is to profile a representative coupled experiment rather than optimize an assumed bottleneck.

10. **Save enough state to explain and reproduce a run.** Configuration and table identities, model version, controls, field state, numerical caches needed for continuation, and diagnostic histories all matter. Preserve the restart concept. Replace opaque failure labels with useful explanations and retain the last inspectable state. Reproducibility does not require universal cross-device bit identity.

**Specific code to consult**

These are references for extraction or redesign, not approved drop-in components. Read their tests and assumptions before carrying code into the new implementation.

| Need | Specific archived chunk | Evidence and caution |
|---|---|---|
| Axisymmetric geometric metrics | [grid/src/lib.rs](archive/pre-pivot-2026-09-17/crates/grid/src/lib.rs#L1036): `face_radius`, `ring_radii`, `cell_volume`, `face_area_r`, `face_area_z` | Shared-face ownership and cylindrical measures. Extract from the larger 3D/brick architecture. |
| Compressible-gas numerical flux | [euler/hllc.rs](archive/pre-pivot-2026-09-17/crates/solvers/src/euler/hllc.rs#L18): `physical_flux`, `hllc_flux`; [recon.rs](archive/pre-pivot-2026-09-17/crates/solvers/src/euler/recon.rs#L34): `ppm_faces` | Review EOS and fixed state-component assumptions. This is gas-dynamics code, not an electromagnetic/plasma solver. |
| Primitive-state recovery from thermochemistry | [euler/table_eos.rs](archive/pre-pivot-2026-09-17/crates/solvers/src/euler/table_eos.rs#L119): `TableEos::bind`, `cons_from_phz` | Useful lookup/projection pattern; medium coordinates and admissible states must fit the new model. |
| Fast checked table access | [tables/src/bound.rs](archive/pre-pivot-2026-09-17/crates/tables/src/bound.rs#L177): `BoundColumn::interpolate`; [interp.rs](archive/pre-pivot-2026-09-17/crates/tables/src/interp.rs#L60): reference API | Checks table domain and narrower validity envelope, then interpolates in declared coordinates. Do not assume table interpolation error equals physical-model uncertainty. |
| Content-addressed data | [tables/src/digest.rs](archive/pre-pivot-2026-09-17/crates/tables/src/digest.rs#L1): v3 encoding; [Python tables.py](archive/pre-pivot-2026-09-17/offline/crucible_offl/tables.py#L139): `content_digest` | Length-prefixed metadata, units, data, and provenance; check [cross-language fixture tests](archive/pre-pivot-2026-09-17/crates/tables/tests/fnd5_python_seam.rs#L26). Old logs mention earlier digest versions; consult current archived code. |
| Chemistry-library interface | [chemistry.py](archive/pre-pivot-2026-09-17/offline/crucible_offl/chemistry.py#L135): `EquilibriumEngine`, `state_php`, `chamber`, `performance`; `FrozenReactantEngine` at line 404 | Explicit SI conversions and reactant/equilibrium distinction. Equilibrium performance is a reference, not validation of mixing or finite-rate combustion. |
| Audit of coupled transfers | [sdc.rs](archive/pre-pivot-2026-09-17/crates/solvers/src/sdc.rs#L1928): `audit_stored`, `audit_check`, `audit_check_n` | Stored changes compared with transfers using the integration weights actually applied. Extend accounting for fields, reaction products, and external work; do not inherit all SDC interfaces. |
| Save/resume | [gpu_engine_run.rs](archive/pre-pivot-2026-09-17/crates/gpu/src/gpu_engine_run.rs#L94): `write_checkpoint`, `read_checkpoint` | Includes physical state, hints, clock, trackers, and input digests. Temporary-write/rename is useful; inspect crash-durability and format compatibility requirements afresh. |
| Measurements and spatial export | [engine/run.rs](archive/pre-pivot-2026-09-17/crates/engine/src/run.rs#L1515): `plane_area`, `plane_mdot`, `plane_thrust`; `fields_csv` at line 1716 | Useful reference for measurements and first field views. Gas-plane thrust alone does not supply complete magnetic, radiation, or transient force accounting. Avoid coupling the new UI to the old run controller. |
| Small numerical reference problems | [Sod tests](archive/pre-pivot-2026-09-17/crates/solvers/tests/solv1_station1_sod.rs#L19); [manufactured-solution tests](archive/pre-pivot-2026-09-17/crates/solvers/tests/solv1_euler_mms.rs#L11); [table holdouts](archive/pre-pivot-2026-09-17/offline/tests/test_station3_surfaces.py#L59) | Preserve the physical questions and independently expected answers. Reconsider exact-output and architecture-specific assertions when the model changes. |
| GPU implementation references | [gpu/src/engine_host.rs](archive/pre-pivot-2026-09-17/crates/gpu/src/engine_host.rs); [residency_diffusion.cu](archive/pre-pivot-2026-09-17/crates/gpu/cuda/residency_diffusion.cu); [residency_engine_xcheck.rs](archive/pre-pivot-2026-09-17/crates/gpu/src/residency_engine_xcheck.rs) | Device-resident state and CPU/device comparisons are worth studying. These depend on the old state layout and coupling; they are not the foundation of the new architecture by default. |

**Three small patterns, extracted verbatim**

Shared radial-face arithmetic, from `crates/grid/src/lib.rs:1036`:

```rust
pub fn face_radius(&self, f: usize) -> f64 {
    self.spec.r_min + f as f64 * self.spec.dr
}
```

The value is that neighboring cells refer to this same definition, not that the formula is elaborate. Preserve that invariant when redesigning the grid.

Unambiguous string hashing, from `crates/tables/src/digest.rs`:

```rust
fn put_str(h: &mut Sha256, s: &str) {
    h.update((s.len() as u64).to_le_bytes());
    h.update(s.as_bytes());
}
```

This avoids using a separator that can also appear inside metadata. The full schema and ordering must remain explicit; copying this helper alone does not define a reproducible file format.

The energy flux in the old gas solver, from `crates/solvers/src/euler/hllc.rs:28`:

```rust
f[n] += p;
f[I_EN] = u_n * (e_tot + p);
```

This compactly exposes pressure's contribution to momentum and energy transport. Preserve the governing accounting when redesigning interfaces; do not reuse this gas-only expression as the total energy flux for electromagnetic fields and energetic particles.

**Design decisions to retain, reconsider, and retire**

Retain explicit units and material/data provenance; one owner for each exchange; conservative spatial measures; model-domain checks; independently specified reference problems; reproducible experiments; useful saved state; and performance work guided by measurements. Expose these through understandable displays and diagnostics rather than requiring users to read internal logs.

Reconsider the state representation, EOS/reaction closure, grid and time integrator, input format, data storage, and CPU/GPU division against the supported physics and the human workflow. The archived Rust/Python/CUDA work is available, but the language split and HDF5-only interface are not governing constraints. Prefer a shared experiment and result representation for the interface and batch work so they cannot silently run different physics.

Retire mandatory full 3D, automatic axisymmetric-to-3D ambition, universal single-law blending, mandatory physical cold startup of an entire engine, numerical-verdict-as-real-feasibility language, unconditional certificate byte equality, per-document crate architecture, and the requirement to finish physics before visualization. Detailed structures, service lifetime, pumps, whole-vehicle simulation, and arbitrary liquid injection are not automatic obligations of the new medium/reaction scope.

**What remains genuinely unresolved**

The plasma/fluid/particle treatment, appropriate transport closures, first antimatter species and configuration, defensible post-ICF data/handoff, and specific validation anchors need selection. Chemical reacting-chamber models must be assessed in their own right; the old burn-progress closure is not accepted just because code exists. None of these is settled by writing the new vision.

The new instrument's usefulness also needs a literature comparison and human testing. The first implementation should let a person define, run, inspect, and compare a small supported experiment, while magnetic and particle cases test the physical design early. A chosen investigation can then test the tool's research value without becoming the only question the instrument can ask.

**Archive and reset verification**

The archive manifest records the pre-pivot revision, branch, original paths, hashes, and modes for 236 tracked files and four untracked planning documents. The archive operation preserved them; no legacy source was rewritten. Ignored local artifacts were moved but are not included in Git. The prior active CI and agent-specific settings are archived, so they cannot impose the previous roadmap on the root project.

This reset was checked for archive integrity and active-document consistency. Legacy solver tests, physical benchmarks, and GPU runs were not repeated, because this change neither alters nor adopts that implementation. The old environment may need recreation after relocation. No claim of newly validated physics follows from the reset.
