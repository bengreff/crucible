# COUP-2 — Conservation Audit & Operator Coupling

| Field | Value |
|---|---|
| **ID** | COUP-2 |
| **Family** | COUP (Coupling & orchestration) |
| **Status** | Reviewed (2026-08-14) |
| **Depends on** | FND-2, COUP-3; FND-1, COUP-8 |
| **Version** | 0.3.1 (2026-08-25 plan S8: mixed-N_θ ring-interface faces are interior to the one flux ledger — computed once, bit-exact telescoping, no new row). 0.3 (2026-08-24 plan S7: Robin-Robin acceptance `EPS_ROBIN_RESID` relaxed 1e-6 → 1e-4 for the near-vacuum cold-fill startup regime + new absolute floor `EPS_ROBIN_Q_FLOOR_W` = 1e-2 W on the residual's relative scale; the audit's `burn_progress` row carries the class-`R` applied increments). 0.2 (2026-08-14 review fix wave) |

---

## 0. Purpose

COUP-2 owns two things on the one grid: (1) the **every-step global conservation/consistency audit** — mass,
species, momentum, energy (and charge where relevant) — computed as `Δ(stored) = Σ(port fluxes) + sources`,
halting on any violation beyond tolerance; and (2) the **operator-coupling contract** — which operator reads
and writes which field of the conserved state `U`, and the **radiation-partition invariant** (emit-once/
transport-once) that keeps SOLV-2 honest. The audit is *the accountant made structural*: because the field
update is flux-form, the interior fluxes telescope and the global balance collapses to **port accounting**.

Read after FND-2 (the ledger + ports it audits) and COUP-3 (the step it audits after).

## 1. Scope & razor ruling
**Owns:** the conservation-audit algorithm + tolerance; the operator read/write contract on `U`; the
radiation-partition invariant; the **wall-exchange (Robin-Robin) coupling** contract and its deterministic
fixed-sweep solve. **Defers:** the time-integration schedule → **COUP-3**; the operators → **SOLV-\***; the
boundary-object port set → **COUP-7**; the `Manifest.couplers[]` declaration → **COUP-8**; the reduced-mesh
transfer of the retired COUP-1 (there is one discretization now, FND-2 §3.4).

**Razor ruling:** infrastructure, but load-bearing — a conservation violation is a **bug that halts**
(META-1 Principle 6), never a tolerance to live with. "Roughly conserved" is forbidden anywhere.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Conservation ledger** | COUP-3 (post-step) | per-quantity stored totals + port-flux accumulation; audited every step |
| **Operator-coupling map** | SOLV-\*, COUP-3 | the declared reads/writes of each operator on `U` (closed enums, COUP-8); a write outside the declaration is a defect |
| **Radiation-partition invariant** | SOLV-2/4 | energy emitted as radiation is subtracted from local deposit; M1↔Sₙ band split conserves total radiative energy |
| **Wall-exchange coupling** | SOLV-1 (gas-side wall function) ↔ solid conduction; COUP-7 (coolant side only) | Robin-Robin interface with a deterministic fixed-sweep (Aitken/IQN-accelerated) solve, placed per §3.5 |

**Invariant:** `Δ(stored) = Σ(port fluxes) + Σ(sources)` holds to the **named per-quantity tolerance
`TOL_AUDIT[q]`** (§3.1.1, derived + manifest-recorded) every step, with `Σ(sources)` the integrator-applied
increments (§3.1); a violation halts with a diagnosis (quantity, location, magnitude). Ports obey a **Modelica across/through**
discipline (through-flows sum to zero at a node).

## 3. Method

### 3.1 Flux-telescoping audit — a bookkeeping identity over integrator-applied increments
Because SOLV-1's update is **flux-form**, interior face fluxes cancel pairwise and the global sum of any
conserved quantity collapses to its **boundary/port fluxes** (LeVeque; Lax-Wendroff conservative form). The
audit therefore accumulates port fluxes and compares to the change in stored totals — an `O(ports)` check, not
an `O(cells)` re-integration. [META-3: `fv-telescoping`]

**`Σ(sources)` is defined as what the integrator applied (O8).** Each operator accumulates, at commit time,
the **exact increment it added to `U`** — the SDC-node-weighted explicit and implicit contributions as
committed (COUP-3 §3.1) — into the per-quantity ledger. The audit is a **bookkeeping identity** over these
integrator-applied increments; it is **never** an independent `rate × Δt` recomputation, whose `O(Δt²)`
disagreement with the SDC splitting structure would swamp any usable tolerance. What the audit certifies is
that *no increment escaped the ledger*, not that a re-derived rate matches.

**Mixed-N_θ ring-interface faces are interior faces (S8, 0.3.1).** Where bricks of differing azimuthal
resolution meet (FND-2 §3.4), the interface flux is computed **once** at the fine resolution and the coarse
cell applies the area-weighted aggregate of the same numbers — so the exchange telescopes **exactly** (to
the bit, not to a tolerance) and contributes nothing to the port ledger. The audit needs no new row and no
widened tolerance for a mixed-N_θ world; a scheme that computed the two sides independently would be a
defect this section forbids, not a tolerance to absorb.

**Precomputed-kernel deposition (the E-1 kernel mode):** the fraction of a kernel's source that deposits
outside the domain is an **explicit boundary-ledger (escape) term** per step, so kernel deposition closes the
identity exactly. Kernel **normalization** (deposit + escape = source, per kernel) is verified **offline** in
the kernel pipeline's verification (OFFL-1/VAL-3), never re-integrated per step.

### 3.1.1 Audit tolerance — `TOL_AUDIT` *(named, derived, manifest-recorded — O7)*
Per conserved quantity `q`, per step:
`TOL_AUDIT[q] = max( K_AUDIT · ε_machine · √N_cells · S[q] , TOL_AUDIT_FLOOR[q] )`, where
- `ε_machine` = f64 unit roundoff (~2.2e-16); `√N_cells` is the stochastic accumulation scaling of the
  fixed-shape tree reduction over the active cells (Higham) [META-3: `higham-rounding`, `repro-sum`];
- `S[q]` = the per-quantity magnitude scale: the **gross stored magnitude `Σ κV·|q|`, summed at both step endpoints** (the after-state's gross content carries the never-ledgered interior telescoping fluxes' composed magnitudes — the θ sweep's, above all), plus cumulative gross
  port/source throughput of `q` this step (so cancellation-heavy balances are scaled by what was actually
  summed, not the small net — on BOTH sides of the identity: the stored term was the |net total| until S8,
  when the first mirror-symmetric-swirl world measured the net θ-momentum cancelling to ~0 while the stored
  reduction's own rounding did not, collapsing the tolerance below round-off);
- `K_AUDIT` = named safety factor (default 100 — headroom for flux-aggregation and ledger arithmetic beyond
  the bare reduction);
- `TOL_AUDIT_FLOOR[q]` = per-quantity absolute floor (default: `K_AUDIT · ε_machine · √N_cells ×` a declared
  reference scale for `q` from the config's operating profile), so an identically-zero stored total cannot
  make the tolerance vanish.
All four are **named constants recorded in the run manifest** (FND-4). At RL10-class grids (~10⁵ cells) this
lands near `1e-10` **relative to `S[q]`** — the §6 test's earlier "~1e-10 relative" phrasing was this same
number stated loosely; the binding definition is the absolute, derived `TOL_AUDIT[q]` above.

### 3.1.2 Momentum audit — wall-force / mount-reaction ledger *(O9)*
Momentum cannot telescope to ports alone: flow pushes on **anchored material** (walls, structure). Rule:
**anchored-solid momentum is not integrated** — anchored material is held by an implicit mount reaction (the
ship, outside the sim, VISION_SCOPE §4.2). The per-step **impulse applied to anchored material** — pressure +
wall-function tangential force integrated over the config-time-identified wall faces (§3.5) — accumulates
into an explicit **mount-reaction ledger term**, closing the momentum identity:
`Δ(fluid momentum) = Σ(port momentum fluxes) − impulse_on_anchored + Σ(applied increments)`.
The accumulated mount reaction is the thrust-frame force; at steady state it must agree with **SOLV-7's**
exit-plane momentum-flux thrust integral (SOLV-7 owns thrust; the agreement is a §6 cross-check, not a second
thrust definition).

### 3.2 Ports as Modelica connectors
Each port carries an **across** variable (potential — p, T; equal at a node) and a **through**/`flow` variable
(rate — mass/energy/momentum; summing to zero at a node). This makes the audit a **sum-to-zero** check per
node and matches the boundary-object contract (COUP-7): a boundary object is a port term, never an interior
source that escapes telescoping. [META-3: `modelica-connector`]

### 3.3 Operator-coupling contract
Every operator declares (COUP-8 `Manifest.grid_fields`, `couplers[]`) exactly which components of `U` it reads
and writes, over **closed enums**. **Enforcement mechanism (O11), three layers, never conflated:** (1)
**test-time property checks** (VAL-3) — fixtures that drive each operator and assert no undeclared field
changed; (2) **debug-build write asserts** — every field write path asserts the writer's declaration in debug
builds; (3) **release runtime carries no per-write checks** — it relies on **Manifest-granted typed field
access**: an operator receives accessor views only for its declared fields, so an undeclared write has no code
path. A violation surfacing anywhere is a defect, not a tolerated event. This is what lets the audit be
*sound* — there is no hidden coupling to slip past it. The seven coupler kinds of VISION_SCOPE §7.4 are the closed
vocabulary (unchanged — no cycle-power coupler, COUP-7 §3.2).

### 3.4 Radiation-partition invariant
Energy leaving a cell as radiation (SOLV-2) is **subtracted from that cell's local deposit** (OpenMC KERMA
discipline), so energy is neither double-counted nor lost across the emit→transport→absorb chain; the M1↔Sₙ
band split conserves total radiative energy at the aligned edge. COUP-2 audits this as part of the energy
balance. [META-3: `radiation-partition`, `sn-deposition`]

### 3.5 Wall-exchange (Robin-Robin) coupling *(delineated per D-C/O10 — one owner per quantity)*
Gas↔solid heat exchange is a **conjugate-heat-transfer** interface; Dirichlet-Neumann coupling is unstable
for Biot>1 (thin-fluid/solid-conductive — the rocket-wall regime), so COUP-2 uses a **Robin-Robin** interface
with an optimized coupling coefficient (unconditionally stable), solved by a **fixed number of Picard sweeps**
with **Aitken/IQN-ILS acceleration**. The fixed count (not a variable tolerance stop) preserves Tier-1
determinism (META-1 §2.2). [META-3: `robin-robin-cht`, `aitken-iqnils`]

**Delineation:**
- **SOLV-1** evaluates the **one local wall-function law** `h` (SOLV-1 §3.5 — the closure is stated there
  once, not here) at **config-time geometrically-identified wall faces** — face selection is data from
  apertures/volume fractions, never an `if(material)`; the resolved `F_visc` is **suppressed at those faces**
  (the wall function replaces, never adds).
- **COUP-2** owns the **sweep placement**: the wall faces enter each SDC node's **spatially-coupled implicit
  diffusion solve** (COUP-3 §3.1 class `D`) as Robin interface conditions with the wall-function `h` as the
  Robin coefficient; the fixed Picard/Aitken sweeps run **inside that class-`D` solve, once per SDC sweep**.
  The post-sweep converged wall-heat integral is the jacket enthalpy rise **COUP-3 §3.5**'s expander solve
  consumes after the sweeps.
- **COUP-7**'s cooling-jacket boundary object owns **only the coolant side** (channel correlation, coolant
  return state).

### 3.6 Determinism
The audit reduction uses fixed brick-chunk order + fixed-shape tree combine (bit-identical at any thread
count, FND-2 §3.7). The wall-exchange solve uses a fixed sweep count. No part of the audit or coupling depends
on schedule or hash order.

## 4. Coupling relationships
- **COUP-3** calls the audit after each advance; a violation halts the step (COUP-4 verdict).
- **SOLV-1** provides face fluxes (the telescoping quantity); **SOLV-2** the radiation exchange; **COUP-7**
  the boundary-object port terms; **COUP-8** the declared coupler/grid-field enums the audit enforces.
- **FND-2** provides the stored ledger + port hooks; **COUP-6** may note a persistent near-tolerance imbalance
  in pedigree.

## 5. Uncertainty & validity
COUP-2 originates **no physics uncertainty**; it is a correctness gate. Its tolerance is the derived,
manifest-recorded named constant of §3.1.1 (not a magic number). A run that cannot hold conservation is
halted, not reported with a caveat. Validation-ladder status: infrastructure (rung i), CI-gated.

## 6. Validation plan
1. **Closed-box conservation:** a no-port field-physics sweep conserves mass/energy to `TOL_AUDIT` (§3.1.1).
2. **Port closure:** `Δ(stored) = Σ fluxes` closes to `TOL_AUDIT[q]` with boundary objects present (the
   derived absolute tolerance of §3.1.1 — supersedes the earlier "~1e-10 relative" phrasing).
3. **Radiation partition:** transported energy = subtracted deposit; band-split conserves radiative energy.
4. **CHT stability:** the Robin-Robin solve is stable for Biot>1 and converges within the fixed sweep count.
5. **Declared-or-forbidden:** an undeclared `U` write is caught by the test-time property check and the
   debug-build assert (§3.3).
6. **Determinism:** identical audit results at 1 vs N threads.
7. **Momentum ledger:** a thrusting fixture closes the momentum identity with the mount-reaction term
   (§3.1.2); the steady-state mount reaction matches SOLV-7's thrust integral within `TOL_AUDIT`.
8. **Kernel escape:** a kernel-deposition fixture with partial escape closes the energy identity via the
   boundary-ledger escape term (§3.1); a planted mis-normalized kernel is caught offline, not per step.

## 7. References
META-3 keys: `fv-telescoping`, `modelica-connector`, `radiation-partition`, `sn-deposition`,
`robin-robin-cht`, `aitken-iqnils`, `repro-sum`, `higham-rounding`. Depends on FND-2, COUP-3 (§3.1 operator
classes, §3.5 expander solve), COUP-8, COUP-7, SOLV-1 (§3.5 wall function), SOLV-2, SOLV-7 (thrust
cross-check).

*(No open questions.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-08-25 | 0.3.1 | **Plan S8 (landed with the code).** (1) **Mixed-N_θ interface faces stated interior to the one flux ledger** (§3.1 note): the FND-2 §3.4 ring-interface exchange computes each N_θ-jump face's flux once (fine side) and applies the aggregate coarse-side — bit-exact telescoping, no port term, no new audit row. (2) **§3.1.1 `S[q]`'s stored term is the gross magnitude `Σ κV·|q|`**, not the |net total|: the first mirror-symmetric-swirl world (the S8 `axis_pulse_3d` gate) measured net θ-momentum cancelling to ~0 while the stored reduction's cancellation rounding did not — the net-scaled tolerance collapsed below round-off and the audit halted on its own arithmetic. The gross form is what §3.1.1's own rationale ("scaled by what was actually summed") always meant; a tolerance change only (halt-gate), no accepted number moves. |
| 2026-08-24 | 0.3 | **Plan S7 (landed with the code): the Robin-Robin acceptance meets the cold start.** (1) **`EPS_ROBIN_RESID` relaxed 1e-6 → 1e-4** (measured: the startup march's near-vacuum cold fill (~10² Pa) weakens the exchange map's contraction — the wall-adjacent gas cell's thermal mass is ~10⁴× smaller than the stations' dense fills, and the fixed five sweeps land at ~3e-5 relative on a ~10 W exchange = 0.3 mW of staleness). The acceptance is a **halt-gate, not a solution modifier**, so no accepted number moves anywhere and the stations' certificates stay **byte-identical**; 1e-4 still sits 3+ orders below the wall law's ±20–30 % band — the COUP-3 §3.5 sizing principle. (2) New **absolute floor `EPS_ROBIN_Q_FLOOR_W` = 1e-2 W** on the residual's relative scale: a ~zero-heat cold start cannot fail on noise over nothing; a micro-thruster-class config whose *real* exchange lives at that floor must revisit the constant, loudly. (3) The audit's **`burn_progress` row now carries the class-`R` applied increments** (net + gross) exactly as applied — the §3.1 integrator-applied-increment rule extended to the new implicit occupant (COUP-3 0.4.4). |
| 2026-08-14 | 0.2 | **Review fix wave (O7, O8, O9, O10/D-C, O11).** §3.1: `Σ(sources)` defined as the **integrator-applied increments** (SDC commit-time ledger; never a rate×Δt recomputation); kernel-mode **escape fraction** as an explicit boundary-ledger term, normalization verified offline (O8). New §3.1.1: **`TOL_AUDIT[q]`** — derived named tolerance (`K_AUDIT·ε_machine·√N_cells·S[q]`, absolute floor), manifest-recorded; §6.2's "~1e-10 relative" reconciled as the same number at RL10 scale, binding form absolute (O7). New §3.1.2: **momentum audit** — anchored-solid rule + per-step mount-reaction impulse ledger term + SOLV-7 thrust cross-check (O9). §3.3: declared-or-forbidden **enforcement mechanism** split into test-time property checks / debug-build write asserts / release Manifest-granted typed field access (O11). §3.5: wall-exchange **delineation** — SOLV-1 evaluates the one wall-function `h` at config-time geometric wall faces with `F_visc` suppressed (closure stated in SOLV-1 §3.5 only); COUP-2 owns Robin-Robin placement inside the SDC class-`D` solve (COUP-3 §3.1), feeding COUP-3 §3.5; COUP-7 owns coolant side only (O10/D-C). §2/§5/§6 aligned; validation items 7–8 added. |
| 2026-07-21 | 0.1 | Initial draft. Flux-telescoping every-step conservation audit (Δstored = Σ port fluxes) with Modelica across/through ports; operator-coupling declared-or-forbidden contract over closed enums (seven couplers, unchanged); radiation-partition invariant (emit-once/transport-once); Robin-Robin wall-exchange CHT coupling with deterministic fixed-sweep Aitken/IQN acceleration (produces the expander jacket ΔH); fixed-order-reduction determinism. |
