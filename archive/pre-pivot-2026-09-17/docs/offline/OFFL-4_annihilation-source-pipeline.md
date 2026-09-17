# OFFL-4 — Annihilation Source Pipeline

| Field | Value |
|---|---|
| **ID** | OFFL-4 |
| **Family** | OFFL (Offline pipeline) |
| **Status** | Draft |
| **Depends on** | FND-5 (schema), FND-1 (model-form band); feeds SOLV-4 |
| **Version** | 0.2 |

---

## 0. Purpose

OFFL-4 is the **offline Geant4 pipeline** that produces antiproton-annihilation **source terms** — product
spectra, energy partition, and per-nucleus fission probabilities — for SOLV-4 to emit. It moves **no
particles**: annihilation products are transported by SOLV-2/SOLV-3 at runtime. Read after FND-5 (the table
contract) and SOLV-4 (the consumer).

## 1. Scope & razor ruling
**Owns:** the Geant4 annihilation runs and the emitted source-term tables + the declared model-form band.
**Defers:** all transport of the products → SOLV-2/SOLV-3; the emission at runtime → SOLV-4; the schema →
FND-5.

**Razor ruling & firewall (§10.2/§10.3):** published spectra / energy partitions / fission probabilities
**only** — never implosion-to-criticality or device-physics optimization. Geant4 (permissive); **FLUKA only
if licensable** for an independent student, else Geant4-only with the wider band.

## 2. Interfaces & contracts
| Interface | Consumed by | Contract |
|---|---|---|
| **Annihilation source tables** | SOLV-4 (via FND-5) | product spectra (pions/charged fragments, neutrons) + energy partition + per-nucleus fission probabilities; **source terms only** |
| **Per-quantity model-form ensemble** | FND-1, COUP-6 | the inter-list spread per emitted quantity (FTFP/INCL++/CHIPS ensemble members, §3) — not a single blanket band |

**Invariant:** OFFL-4 emits sources, never transports; every table carries FND-5 provenance; the pedigree is
permanently rung-(ii) flagged (no hardware anchor beyond thin-target benchmarks).

## 3. Method
- **Physics list:** Geant4 **FTFP_INCLXX** — **INCL++** is now the default/recommended model for antiproton
  annihilation at-rest/in-flight — ~1e6 events per configuration (minutes–hours; effectively free per §8).
  [META-3: `incl++`]
- **Model-form disputes become UQ (META-1 §4):** run **FTFP and CHIPS as model-form ensemble members**, not
  just a central ± band — the inter-list spread *is* the uncertainty. Anchor to data: `pbar-annihilation`
  (⟨n_π⟩≈5, ~16 n/annihilation on U, ~52% via fission, fragments ~160 MeV), `pbar-fission-prob`
  (Bocquet/PS177), `pbar-modern-data` (ASACUSA/AEgIS 2024–25).
- **Cross-check:** **SHIELD-HIT12A** (independent, non-Geant4) for an energy-partition second opinion.
  [META-3: `shieldhit`]

## 4. Coupling relationships
- **SOLV-4** emits the tabulated sources into SOLV-2/SOLV-3; **FND-5** stores them; **FND-1** carries the band;
  **COUP-6** flags the pedigree; **VAL-2** holds the `antimatter-ps177` anchor.

## 5. Uncertainty & validity
**Model-form uncertainty = the per-quantity inter-list spread**, carried as ensemble members (INCL++ central;
FTFP/CHIPS as edges), **never collapsed to a single blanket number**: FTFP underestimates heavy fragments
~12×, prong-multiplicity discrepancies reach ~4×, neutron-dose spread ~50%, while prong multiplicity itself is
pinned to a few % by PS177/AEgIS. A single headline band — whether the legacy "15–30%" *or* a blanket
"factor-of-several" — both **over-states** the well-measured quantities and **under-states** the divergent
ones; only the per-quantity spread is honest. **Conversion efficiency and waste-heat load are not emitted here
at all** — they are computed at runtime by transporting these products (neutrinos escape; γ/charged deposit
per geometry, SOLV-2/3). C3's claim is a *relative* dose-reduction factor, robust to correlated model error.
[META-3: `pbar-list-band`]

> **Resolved (R3, 2026-08-13):** the VISION_SCOPE "15–30%" blanket was retired via the §15 ceremony (§5.2/§9
> amended). Per Ben's ruling, neither "15–30%" nor a blanket "factor-of-several" is the right frame — model-form
> is the **per-quantity physics-list spread** (above), and efficiency/waste-heat are **computed by product
> transport** (SOLV-2/3), not assigned. No doctrine change; the source-only firewall and rung-(ii) pedigree stand.

## 6. Validation plan
1. **Spectra/multiplicities** vs PS177 / arXiv:2407.06721 within the declared band.
2. **Inter-list spread** captured as ensemble members (FTFP/INCL++/CHIPS).
3. **SHIELD-HIT12A** energy-partition cross-check.
4. **Provenance:** physics-list + Geant4 version + seed pinned; regenerable.

## 7. References
META-3 keys: `incl++`, `pbar-list-band`, `shieldhit`, `pbar-annihilation`, `pbar-fission-prob`,
`pbar-modern-data`, `antimatter-ps177`. Depends on FND-5, FND-1; feeds SOLV-4, COUP-6, VAL-2.

*(No open questions — physics list (FTFP_INCLXX central, FTFP/CHIPS as ensemble members) resolved 2026-07-21;
model-form carried as the **per-quantity inter-list spread**, not a blanket band; VISION_SCOPE §5.2/§9 "15–30%"
retired via §15 on 2026-08-13 per Ben's R3 ruling; efficiency/waste-heat computed downstream by SOLV-2/3.)*

## 8. Change log
| Date | Version | Change |
|---|---|---|
| 2026-07-21 | 0.1 | Initial draft. Geant4 FTFP_INCLXX (INCL++ central) annihilation source-term pipeline; FTFP/CHIPS as model-form ensemble members; SHIELD-HIT12A cross-check; **factor-of-several** declared band (supersedes 15–30%; VISION_SCOPE reconciliation flagged for Ben). Emits source terms only; products transported by SOLV-2/3; firewall-clean. |
| 2026-08-13 | 0.2 | R3 reframe: blanket "factor-of-several" headline demoted to the **per-quantity inter-list spread** it already computes; efficiency/waste-heat clarified as **computed downstream by SOLV-2/3 product transport**, not emitted here. VISION_SCOPE §5.2/§9 "15–30%" retired via §15 (2026-08-13, Ben's R3 ruling). |
