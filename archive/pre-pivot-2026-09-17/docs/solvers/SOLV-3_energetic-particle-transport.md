# SOLV-3 — Energetic-Particle-Transport Operator

| Field | Value |
|---|---|
| **ID** | SOLV-3 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, FND-7 (collision coefficients), COUP-3, COUP-2; SOLV-2, SOLV-4 |
| **Version** | 0.1 |

---

## 0. Purpose

SOLV-3 is the **one transport operator for all fast charged particles** — fission fragments, fusion α/p,
annihilation-produced charged fragments/protons, sputtered ions, fast electrons — evolving their phase-space
distribution, self-heat deposition, and orbits over the medium-state vector `M`. It carries **no constitutive
law**: the drag, velocity-space diffusion, and effective charge are **queries into the FND-7 constitutive
spine** (the same stopping physics for a fragment in a wall and an α in a plasma — one law, only `M` differs).
This is the single biggest seam risk in the project, closed by construction.

Read after FND-7 (the collision-coefficient bundle it queries), SOLV-4 (the birth sources it consumes),
SOLV-2 (the photon field it emits into), and META-1 Rule 12.

## 1. Scope & razor ruling
**Owns:** the *transport* of fast charged particles — the phase-space evolution (drag + diffusion), the
range-test deposition split, and the orbit integration. **Defers:** the drag/diffusion/`Z_eff` **coefficients**
→ **FND-7** (`spine(M).stopping` = the collision bundle `{drag_e, drag_ion, D_∥, D_⊥/ν_pitch, Z_eff}`, FND-7
§3.5); reaction/annihilation *sources* → **SOLV-4/OFFL-4**; neutrons/photons → **SOLV-2**; time integration →
**COUP-3**.

**Razor ruling:** slowing-down *is* self-heat, so SOLV-3 is full reaction physics — but the medium properties
are the spine's. A private plasma-only friction law, or a separate "fragment" transport, would be a Rule-12
seam and is forbidden.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Collision-coefficient query** | FND-7 | `{drag_e, drag_ion, D_∥, D_⊥/ν_pitch, Z_eff(v;M)}` from the one AA dielectric response — no material label |
| **Birth-source intake** | SOLV-4/OFFL-4 | fast-charged birth spectra |
| **Self-heat write** | FND-2 (`U`) | fast-species number/momentum/energy as **conserved fields** (COUP-2 telescopes them) + self-heat to **T_e and T_i** channels (split at v_c) |
| **Radiative emission** | SOLV-2 | brems/synchrotron photons as birth sources (emit-once/transport-once) |
| **Distribution export** `f(E)` | SOLV-4 | slowing-down `f` for in-flight secondary-reaction rates |

**Invariant:** every coefficient is a spine query (no private law); the same operator handles every fast
species (source spectrum differs, physics does not); terms self-zero **by physics** (ion synchrotron ∝1/m⁴),
never by an `if(species)` branch.

## 3. Method & governing structure

### 3.1 The adaptive-reduction ladder (one operator, nested conservative reductions)
Analogous to FND-2's adaptive azimuthal resolution (N_θ), SOLV-3 is **one kinetic operator with a nested ladder** sharing
the *same* spine coefficients — switched by a **continuous physical-validity indicator**, never a regime label:
1. **Gaffey analytic slowing-down** `f(v) ∝ Θ(v_b−v)/(v_c³+v³)` — the steady-collisional self-heat path;
   electron/ion split at the critical velocity v_c (`stix-critical-energy`). [META-3: `gaffey-slowing`]
2. **Continuum reduced-Fokker-Planck moments** on the grid — the **runtime default** (RABBIT-class:
   source + slowing-down + pitch-angle scatter; ~1000× faster than Monte-Carlo at benchmarked fidelity).
   [META-3: `rabbit-fast-ion`]
3. **δf markers** — a bounded runtime *fallback* for orbit-width/anisotropy-critical cases (first-orbit loss
   in a compact magnetic nozzle, collimated fragment beams); counter-based RNG, deterministic.
4. **Monte-Carlo guiding-center** (NUBEAM/ASCOT-class) — **offline oracle only** (VISION_SCOPE §5.3).
   [META-3: `nubeam-oracle`]

Rungs 1↔2 switch on τ_s/τ_loss; 2↔3 on an orbit-width/anisotropy indicator with hysteresis (mirroring FND-2's
`symmetry-indicator`). Short-range products (fragments: R ≪ ℓ) never enter transport at all — the **range
test** deposits them locally (a smooth split, not a toggle).

### 3.2 Orbits & radiative losses
Coupled **guiding-center↔Boris** integration switched on the adiabaticity/magnetization indicator ε_gc
(`gc-boris`) in marker mode; a GC-drift fluid limit in moment mode. Radiative cooling is a drag on `f(v)` that
**emits into SOLV-2's photon field**: bremsstrahlung ∝1/m² (`brems-emc`), synchrotron ∝1/m⁴
(`synchrotron-larmor`) — the ion terms vanish by mass scaling, by physics.

### 3.3 Self-heat and two-temperature deposition
The slowing-down moment integral (§3.1 rung 1) gives the self-heat, split into **electron and ion** channels
at v_c (`stix-critical-energy`) and written to the T_e/T_i state (the split is load-bearing for the p-¹¹B
ignition margin). Momentum deposition is written as a conserved field so COUP-2 audits it as an internal
transfer summing to zero.

### 3.4 Determinism
Marker draws use counter-based RNG keyed on stable indices; moment updates and depositions use fixed-order
reductions (FND-2 §3.7). The rung switch is a deterministic function of the field (indicator + hysteresis).

## 4. Coupling relationships
- **FND-7** supplies the collision bundle; **SOLV-4** supplies birth spectra and consumes `f(E)`;
  **SOLV-2** receives brems/synchrotron photons; **COUP-3** advances drag/self-heat as stiff sources;
  **COUP-2** audits the internal energy/momentum transfer; **FND-2** stores the fast-species fields;
  **SOLV-1** (W4 two-fluid) receives the T_e/T_i self-heat split.

## 5. Uncertainty & validity
Consumes the spine's stopping band (wide in WDM, FND-7 §5) + the fast-ion charge-state band (`zeff-betz`); the
moment-closure truncation vs marker/MC is a **declared model-form band** (measured by rung-to-rung
comparison). Validity indicators flag when a rung is out of its regime. Ladder: unit (Gaffey/Spitzer
analytics, Stix split) → benchmark (NIST stopping via the spine, plasma stopping data) → reference-code
(ASCOT/NUBEAM oracle); plasma modules permanently rung-(ii) flagged.

## 6. Validation plan
1. **Slowing-down:** recover Gaffey/Spitzer analytics + the Stix electron/ion heating split. [META-3:
   `gaffey-slowing`, `spitzer-slowing`, `stix-critical-energy`]
2. **Stopping:** NIST PSTAR/ASTAR cold ranges (via the spine) + plasma stopping data — continuous across.
3. **Orbits:** guiding-center↔Boris switch reproduces a known drift orbit; **ASCOT/NUBEAM** differential
   oracle on a magnetized first-orbit-loss case. [META-3: `gc-boris`, `nubeam-oracle`]
4. **Determinism:** identical results at 1 vs N threads.

## 7. References
META-3 keys: `rabbit-fast-ion`, `gaffey-slowing`, `nubeam-oracle`, `gc-boris`, `spitzer-slowing`,
`stix-critical-energy`, `brems-emc`, `synchrotron-larmor`, `zeff-betz`, `stopping-rpa-lda`,
`stopping-li-petrasso`, `stopping-bps`, `radiation-partition`. Depends on FND-7 (collision bundle), SOLV-4,
SOLV-2, COUP-3, COUP-2, FND-2.

*(No open questions — friction coefficients are FND-7 spine queries, not a private law (FND-7 §3.5 broadened
to the collision bundle 2026-07-21); moments-vs-markers is a hysteretic validity switch on one operator, not
an if(regime) branch.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-21 | 0.1 | Initial draft. One transport operator for all fast charged particles; **all coefficients are FND-7 spine queries** (collision bundle {drag, D_∥, D_⊥/ν_pitch, Z_eff}), no private friction law. Nested adaptive-reduction ladder (Gaffey ⊂ continuum reduced-FP [runtime default] ⊂ δf markers [fallback] ⊂ offline MC-GC oracle) switched by continuous validity indicators; range test deposits short-range fragments locally. Guiding-center↔Boris orbits; brems/synchrotron emit into SOLV-2; two-temperature self-heat split at v_c. |
