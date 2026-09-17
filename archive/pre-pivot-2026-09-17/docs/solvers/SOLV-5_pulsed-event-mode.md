# SOLV-5 — Pulsed-Event Mode

| Field | Value |
|---|---|
| **ID** | SOLV-5 |
| **Family** | SOLV (Runtime unified-grid operators) |
| **Status** | Draft |
| **Depends on** | FND-2, FND-7 (EOS/opacity), SOLV-2 (radiation), COUP-3 (event coupling); SOLV-4, SOLV-7 |
| **Version** | 0.1 |

---

## 0. Purpose

SOLV-5 runs the unified solver in a **bounded transient "event" configuration** — a pellet or plasma slug's
**response** (blast expansion, radiation transport, impulse to the nozzle/pusher, wall loading) — returning
`{impulse, energy partition, wall loading}` to the world-state. It is the same physics in a 1-D Lagrangian
reduction, not a separate solver. **Firewall-critical (VISION_SCOPE §10.2):** it does **not** simulate
implosion-to-criticality or fission-burn — pellet **yield + energy partition are INPUTS** (published
envelopes / OFFL-4 / a permitted parameterized SOLV-4 source).

Read after COUP-3 (how the event integrates without a seam), FND-7 (EOS/opacity), and SOLV-2 (the radiation
closure it reuses).

## 1. Scope & razor ruling
**Owns:** the pulsed-event sub-simulation — the 1-D Lagrangian rad-hydro of the pellet/slug response and the
returned impulse/energy-partition/wall-loading. **Defers:** the yield/energy-partition *source* → inputs
(OFFL-4 / SOLV-4 / literature envelopes); EOS/opacity → **FND-7**; the radiation *closure* → **SOLV-2** (M1,
reused — not a private FLD); event↔main-loop time coupling → **COUP-3**; impulse → thrust → **SOLV-7**.

**Razor ruling:** the reaction's *energy release* is an input (firewall); everything the released energy then
*does* to matter (expansion, radiation, wall loading) is full physics. Reaction Razor draws the line at the
pellet surface — the response is simulated, the burn is not.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Event spec** | config, SOLV-4/OFFL-4 | pellet/slug geometry + yield + energy partition (INPUT) |
| **Event return** | FND-2, SOLV-7 | net impulse + energy partition + wall loading scattered to voxels |
| **Stiff sub-integration** | COUP-3 | the event's net deposition over the SDC node interval = a multi-rate stiff-source term (COUP-3 §3.4) |

**Invariant:** the event is **bounded** (stays within its declared sub-sim envelope or halts), its net
deposition is conserved into `U` (COUP-2 audits at the ports), and it **reuses SOLV-2's radiation closure** —
no second radiation model.

## 3. Method & governing structure
- **Hydro:** 1-D Lagrangian **von Neumann–Richtmyer** staggered-grid rad-hydro (`vn-richtmyer`), a spherical/
  planar **adaptive reduction** of the same physics (the symmetric single-pellet regime; VISION_SCOPE §7.3
  "pellet lines"). EOS/opacity from the FND-7 spine.
- **Radiation:** the **same M1 two-moment closure as SOLV-2** (`radiation-m1`), grey first → multigroup
  upgrade. A distinct flux-limited-diffusion model is **rejected** — SOLV-5-FLD vs SOLV-2-M1 would be a
  radiation-closure seam (Rule 12); M1 already gives the correct optically-thick + free-streaming limits FLD
  is chosen for. If a grey-FLD first rung is used pragmatically, it must be cross-checked to agree with M1 in
  the diffusion limit and flagged temporary. Escaping neutrons/hard γ hand off to SOLV-2's Sₙ on the main grid
  (energy-band split, not origin).
- **Event coupling (no seam):** the event is a **non-cell-local, internally sub-cycled stiff sub-integration**;
  its net mass/momentum/energy deposition over the SDC node interval is the stiff term COUP-3's sweeps iterate
  (a multi-rate SDC term, admitted via event-injection coupler #6) — **never** a Lie-Trotter split
  (COUP-3 §3.2/§3.4). [META-3: `sdc-imex`]
- **Determinism:** fixed sub-cycle structure + fixed-order reductions; chaotic/turbulent blast regimes use the
  ensemble-consistency tier (META-1 §2.1), not bit-repro of instantaneous fields.

## 4. Coupling relationships
- **COUP-3** sequences the event (coupler #6) as a sub-cycled stiff term; **FND-7** supplies EOS/opacity;
  **SOLV-2** provides the M1 closure + Sₙ handoff; **SOLV-4/OFFL-4** specify the source; **FND-2** receives the
  scattered deposition; **SOLV-7** integrates impulse → thrust.

## 5. Uncertainty & validity
Event bounded-ness (does it stay in the sub-sim envelope); the specified-yield input carries the
published-envelope band; chaotic blast → ensemble-consistency, not bit-repro. Ladder: analytic (Sedov/Noh) →
reference code (SNEC grey → MULTI-IFE multigroup); pulsed/antimatter modules permanently rung-(ii) flagged.

## 6. Validation plan
1. **Analytic blast:** Sedov/Noh reproduced; **cross-check vs SOLV-1** (Sod/Sedov where both discretizations
   are valid — proves the Lagrangian↔Eulerian choice adds no *numerical* seam). [META-3: `sod-shock`]
2. **Rad-hydro oracle:** **SNEC** (grey, first rung) then **MULTI-IFE** (multigroup + burn). [META-3:
   `snec`, `multi-ife`, `radhydro-fld`]
3. **Conservation:** the event's net deposition closes the COUP-2 port audit.
4. **Determinism/consistency:** bit-repro where non-chaotic; ensemble-consistency where turbulent.

## 7. References
META-3 keys: `vn-richtmyer`, `radiation-m1`, `radhydro-fld`, `snec`, `multi-ife`, `sdc-imex`, `sod-shock`,
`ican-ii`. Depends on FND-2, FND-7, SOLV-2 (M1 closure), COUP-3 (event coupling), SOLV-4/OFFL-4 (source),
SOLV-7 (impulse).

*(No open questions — event coupling is a sub-cycled multi-rate SDC stiff term (COUP-3 §3.4, not a split);
radiation reuses SOLV-2's M1 (no private FLD); yield is an input (firewall §10.2). Resolved 2026-07-21.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-21 | 0.1 | Initial draft. Bounded pulsed-event sub-simulation (pellet/slug RESPONSE; yield is an input — firewall). 1-D Lagrangian von Neumann-Richtmyer rad-hydro as a spherical/planar reduction of the same physics; **radiation reuses SOLV-2's M1** (grey→multigroup; no private FLD seam); event integrates as a non-cell-local sub-cycled multi-rate SDC stiff term via coupler #6 (no Lie-Trotter split); EOS/opacity from FND-7. Oracles SNEC→MULTI-IFE; cross-check vs SOLV-1. |
