# The Space Physics Computational Laboratory
## The Final Vision — Post-Audit Synthesis (2026-07-04)

This document is the complete vision in one coherent read, synthesized from all 18 design
documents after the Phase 0 audit and repairs. It exists to answer one question: *does the
whole thing hold together?* Findings from writing it are in the last section.

---

## 1. THE ONE-PARAGRAPH VISION

A terminal-based physics simulator in which you fly ships through a physically honest
universe. One gravity rule carries you from a launch pad on a rotating Earth, through
max-Q, into orbit, past the planets, out of the Sun's wind bubble, across light-years of
interstellar medium at relativistic speed, and down to the last stable orbit around a black
hole — with no seams, no patched conics, no scripted content. Every one of ~130 billion stars
has a position derived from galactic dynamics; every planet an atmosphere derived from
formation physics; every biosphere exists because a deterministic seed rolled favorably
under literature-calibrated probabilities. The ship is not a marker but a machine: tanks,
reactors, shields, and plumbing that heat up, wear out, break, and cascade. Physics in,
telemetry out. Bit-exact reproducible, always.

## 2. WHAT MAKES IT DEFENSIBLE (the five invariants)

Everything else is negotiable; these five are the architecture:

1. **One force rule everywhere.** F = −∇Φ_mean_field + Σ(point sources above adaptive
   threshold) − double-count corrections, evaluated through the velocity-complete weak-field
   law dp/dt = −γm(1+β²)∇Φ with v = (pc²/E)(1+2Φ/c²). Newtonian at rest, correct photon
   deflection at β→1, Kerr geodesics inside 100 r_s of compact objects (integrator-restart
   event, hysteresis exit at 110 r_s). No hierarchy, no patched regimes, no special cases.

2. **Bit-exact determinism as a testable property.** Position of anything = pure function of
   (seed, T). No integration history in the universe. Physics path and display path fully
   separated. All summations ID-ordered, all transcendentals via libm, background Monte
   Carlo bound to sim-time apply points (deferred-apply protocol). Hardware speed can change
   how long you wait — never where the ship goes. CI enforces this from the first commit.

3. **Expensive physics offline, tables at runtime.** Cantera, OpenMC, Geant4, FEniCS run in
   the characterization tool; the simulator interpolates. The only runtime ODEs are three
   tiny fixed-size systems (thermal network, cabin atmosphere, reactor point kinetics).
   Per-tick cost ~30-115 μs on one core.

4. **The mean field matches the observed Milky Way by construction.** Analytic
   potential-density pairs (McMillan 2017 + rigid bar at 27°/38 km s⁻¹ kpc⁻¹ + Cox-Gomez
   spiral + warp) pinned to measured parameters, plus one precomputed basis expansion for
   the LMC's halo wake. The density that sources gravity is the density that makes stars —
   self-consistency is structural, not aspirational.

5. **The paper is the physics.** Where a subsystem rests on original research
   (black_hole_paper for the drive; the ISM co-evolution paper for the forward shield), the
   paper and its calculation scripts are authoritative and the design docs are summaries.

## 3. THE UNIVERSE

**Window:** 10 million years from J2000.0, galactocentric double-double coordinates
(sub-nm at 30 kpc), double-double time (sub-ns at 10 Myr).

**Galaxy:** the hybrid mean field above; ~170 globular + ~50 open clusters as orbiting
analytic profiles; Sgr A* (4.297×10⁶ M_☉) + NSC + NSD + CMZ at the center; LMC/SMC/Sgr dSph.
Rotation curve 233±3 km/s at R₀ = 8.275 kpc — a validation gate, not a tuning target.

**Stars:** ~130 billion — DERIVED from the mass model + Kroupa IMF (physics.md §3.2.9), not
chosen — generated on demand (~55 ns each) from time-invariant guiding-center
cells, discarded after use. Kroupa IMF, population-correct ages/metallicities/velocity
dispersions, remnants (10¹⁰ white dwarfs, 10⁹ neutron stars with kicks, 10⁸ black holes),
binaries with environment-clipped separations, deterministic flare/supernova schedules via
hierarchical bucket seeding (O(1) at any T). Spiral arms are material overdensities frozen
into guiding-phase space at epoch — drift vs. the potential's rigid pattern is bounded and
documented (≤3° at R₀ over the window). Queries unwind guiding centers (φ_g = φ − Ω(R)T)
before cell lookup. ~1M real Gaia stars overlay the procedural field: complete to 100 pc,
anchors and recorded-fraction samples beyond, real exoplanet systems augmented in stable gaps.

**Planetary systems:** occurrence from Kepler/TESS/RV statistics (giants mass-first from the
mass function; small planets radius-first via Chen-Kipping), peas-in-a-pod architectures,
moons with tidal heating, belts with individually addressable >1 km members, rings, comets,
stellar winds with astropauses. All motion analytic: Kepler + first-order Laplace-Lagrange
secular rates; the solar system uses SPICE inside ±14 kyr, C¹-blended to the same secular
form beyond — a documented self-consistent fiction, because chaos (Lyapunov ~5 Myr) makes
ground truth impossible for ANY method.

**Life:** sequential Poisson stages (abiogenesis λ=1.5 Gyr⁻¹ → photosynthesis → oxygenation
→ eukaryogenesis (hard step) → multicellularity → complex life; no intelligence by scope).
Expected: 2-10 billion microbial worlds, 20-200 million complex biospheres, pigments set by
host-star spectra, biosignatures (O₂, CH₄ disequilibrium, vegetation edges) queryable.
Substellar census: 25-100 billion brown dwarfs; ~2-3 TRILLION rogue planets (~20 per star,
MOA 2023 microlensing — most sub-Earth-mass).

**Environment:** 5-phase ISM, HII regions around every O/B star, SNRs on Sedov-Taylor
clocks, molecular clouds in arms, GCR modulated inside astropauses, trapped belts around
magnetized planets — evaluated per tick at ~0.5 μs.

## 4. THE SHIP

A component graph (primitives + materials + connections) is the single source of truth,
from which derive: a voxel octree (radiation, "what's at X?"), a thermal network (implicit),
a reduced structural model with stress recovery, resource/electrical graphs, and an analytic
surface mesh (aero). Damage is per-component parameterized (ablation depth maps, dpa
profiles, activation inventories) — the octree re-rasterizes exactly from it, so saves stay
1-10 MB and round-trip bit-exactly.

Failures propagate physically: stress exceeds degraded yield → connection breaks → thermal
paths and pressure boundaries change → leaks, brownouts, explosions (blast model), ship
splitting with momentum conservation. Multi-object tracking, docking, collisions. Attitude
control is quaternion-PD with wheel/RCS allocation; maneuvers execute align→burn→cutoff;
open-loop scripted profiles (pitch programs, entry schedules) make validation flights
repeatable. Ground contact: pad constraint + leg spring-dampers + tip-over check.

The standalone Monte Carlo radiation transport crate is both a simulator subsystem and a
research instrument: analytic EM physics + Geant4-derived element tables, mutable octree,
and the MC → 1D-thermal-relief → ablation macro-loop with an explicit fluence↔time mapping.
It runs on the deferred-apply schedule; timewarp is honest (requested vs. granted, with the
binding constraint shown).

## 5. THE THREE SHIPS

| | Falcon 9 / Dragon | NTP Mars Ship | BH-Drive Starship |
|---|---|---|---|
| Regime | Pad → LEO → entry | Interplanetary | Interstellar → compact objects |
| Physics exercised | Atmosphere, aero, staging, ground | Reactor kinetics, shadow shield, aerocapture | Everything, incl. GR + relativistic ISM |
| Key validation | Max-Q ~80 s, MECO ~160 s, ~16 t LEO | ΔV 3.6 km/s, dose <50 mSv/yr | 55.3 MN at 10⁹ kg, τ_unfed 15.8 yr, v_eff 0.276c |

The BH drive is the black_hole_paper architecture: CFL strange-quark-matter shell absorbing
83% of the Hawking spectrum and re-emitting keV pairs (Usov equilibrium, β = 0.332),
anti-Helmholtz confinement (the 2.8 s lateral instability is the ship's defining emergency),
no DRIVE radiators (ship radiators serve reactor/life-support heat only), no off switch —
thrust throttles only through mass drift over months. The
forward ISM shield comes from research paper #1 and erodes over the mission.

## 6. THE RESEARCH PROGRAM (year 1)

1. **Shield-flux co-evolution at 0.3-0.9c** — the transport crate's dynamic-geometry loop
   applied to sustained relativistic ISM bombardment. The cascade physics is well-known;
   the time-resolved shield/damage co-evolution is not. Submits ~month 10-11.
2. **Second paper from the architecture itself** — either the velocity-complete weak-field
   formulation or deterministic procedural galactic dynamics (material-arm scheme + its
   error budget). Drafts in Phase 5.

Year 2+: laser sail (positioned against the Starshot literature), ACMF ignition (positioned
against ICAN-II/AIMStar), antimatter production, curved-spacetime ship characterization.

## 7. TIMELINE (canonical in SCOPE.md)

Gate-based, ~15 months, band 13-20, at ~20 h/wk + AI-assisted implementation, laptop-only
compute (Geant4 tables and the LMC run start as background jobs in month 1):

G1 (M2) determinism + gravity gates → G2 (M5) Falcon 9 validated → G3 (M9) transport engine
benchmarked + NTP → G4 (M12) galaxy + relativistic + BH ship → G5 (M15) every validation row
green and the full mission flown end-to-end. Gate slips trigger re-planning, never silent
schedule erosion. Scope is fixed; time flexes.

---

## 8. SYNTHESIS FINDINGS — what writing this revealed

Writing the vision as one story surfaced six residual issues. Two are real inconsistencies
(fixed in the docs alongside this document), one is a genuine vision-level decision that
needs Ben's explicit call, three are minor.

**F1 — THE MORTALITY PROBLEM (RESOLVED 2026-07-04, Ben's decision: operator framing).**
The flagship narrative — "returns home, finding that millions of years have passed" — is
physically incompatible with a living crew: galaxy crossing costs 10 Myr of proper time at
0.01c and ~200,000 yr even at the drive's realistic 0.45-0.5c ceiling; ISCO time dilation
(1.4×) doesn't rescue it. **Adopted:** the player is an operator, not a passenger — timewarp
is the operator's time machine. Crewed missions are a bounded mode with real life-support/
dose constraints (nearby stars, years-to-decades of proper time — exactly the paper's
mission table); galaxy-scale voyages fly uncrewed, and it is the SHIP that returns to a
changed Earth. Written into SCOPE.md's vision statement. Physics untouched.

**F2 — interface.md contradicted the drive physics (FIXED).** `bh_drive <on/off/...>`
offered an off switch that Hawking radiation does not permit (hawking.md rewrite). Command
set corrected to feed-rate control only.

**F3 — "Time acceleration is unlimited" in SCOPE.md overstated (FIXED).** ship.md §10 grants
warp subject to per-subsystem dt validity and compute/transport ceilings. SCOPE's universe
section now says "unbounded request, physics-limited grant" and points to ship.md §10.

**F4 — Ship 3 has no fuel store (FIXED, minor).** The paper's missions carry 0.7×-33× dry
mass in bulk matter (the sprint's fuel is a 300 m asteroid), but the component list had no
bulk-mass store; the railgun had nothing to feed from. Added a bulk propellant store
component (grappled mass hopper) to Ship 3, sized per mission.

**F5 — Ship 3's design point (RESOLVED 2026-07-04, Ben's decision: mass is a free dial).**
BH mass is fully commandable in operation — one ship spans 55 MN to 21,800 MN. The honest
dynamics are specced in hawking.md "Mass maneuvering": thrust-up (evaporate down) takes
~15 yr from the 10⁹ kg reference but only days at the low-mass end, so the cruise band is a
mission-planning choice; thrust-down is railgun-feed-limited; and the FIXED shell must be
sized for the minimum planned operating mass (kT ∝ 1/M against the stopping budget). The
simulator enforces the stopping-margin gate rather than a hard interlock.

**F7 — THE STAR COUNT WAS INTERNALLY INCONSISTENT (RESOLVED 2026-07-04).** "200 billion"
sat inside the popular 100-400 billion range but contradicted the project's own McMillan
2017 mass model by ~50%: (4.9×10¹⁰ M_☉ stellar − ~8×10⁹ M_☉ remnants) / 0.35 M_☉ Kroupa
mean = ~1.2×10¹¹ living stars, ~1.3×10¹¹ with remnants. Since the architecture's principle
is "the density that sources Φ is the density that makes stars," the count MUST be an
output, not a target. All docs now derive it (~130 billion; physics.md §3.2.9 census table
added with sourced populations). Same pass updated the rogue-planet abundance from the stale
"1-2 per star" to the 2023 MOA microlensing measurement (~20 per star, ~2-3 trillion, most
sub-Earth) and scaled the biosphere expectations (2-10 B microbial, 20-200 M complex).

**F6 — Cosmetics (RESOLVED 2026-07-05 in the consistency re-read).** HUD date fixed to the
J2000-era epoch (campaign start defined per scenario); SCOPE's document table refreshed with
current line counts + AUDIT/VISION rows; physics.md §8 storage total corrected to ~630 MB.

**Verdict (updated 2026-07-04):** F1, F5, and F7 are resolved; F2-F4 fixed; F6 cosmetic.
The vision is internally consistent end-to-end — every subsystem's inputs are produced by
another subsystem, every claim traces to a spec section, every population number derives
from the same mass model that sources gravity, every spec section has a validation row, and
the timeline funds it all.
