# Ship Radiation Module

The ship-side integration of the radiation transport tool (radiation_transport.md). Manages the cached dose-rate map, triggers re-computation, accumulates dose per-tick, and reports crew/electronics exposure.

This module is the CONSUMER of the radiation transport tool's output. It does NOT run Monte Carlo transport -- it reads cached results.

---

## ROLE IN THE PER-TICK LOOP

Phase 3f of ship.md:
```
crew_dose += cached_dose_rate_at_crew × dt
electronics_dose[i] += cached_dose_rate_at_component[i] × dt
component.damage_dpa += cached_dpa_rate[component] × dt
```

Cost: **~0.1 μs per tick.** Pure multiply-accumulate from the cached map.

## THE CACHED DOSE-RATE MAP

Produced by the radiation transport tool (background thread). Contains:

```rust
struct DoseRateMap {
    // Per-voxel (for thermal heating from radiation)
    voxel_dose_rate: Vec<f64>,        // Gy/s per octree leaf
    voxel_dpa_rate: Vec<f64>,         // dpa/s per octree leaf
    
    // Per sub-component (for crew and electronics)
    crew_dose_rate: Vec<f64>,         // Sv/s per crew station
    electronics_dose_rate: Vec<f64>,  // Gy/s per electronics component
    electronics_see_rate: Vec<f64>,   // SEU/bit/s per electronics component
    
    // Metadata (scheduling fields — see the deferred-apply protocol below)
    job_index: u64,                   // which TransportJob produced this map
    computed_from: DoubleDouble,      // T_trig of the producing job
    applied_at: DoubleDouble,         // T_apply — when this map became physics-visible
}
```

## RE-COMPUTATION TRIGGERS

Triggers are evaluated in the per-tick event-detection phase (ship.md Phase 4) as
deterministic conditions on sim state. Each trigger belongs to a class with a fixed
**apply delay Δt_apply (sim time)**, a fixed primary count, and a minimum re-trigger
interval (sim time):

| Class | Triggers | Δt_apply | N_primaries | Min re-trigger |
|---|---|---|---|---|
| ACUTE | SPE onset (seed-determined flare at current star) | 10 s | 2×10⁵ | 60 s |
| GEOMETRY | Shield damage: structural failure or ablation modifies shielding (octree version change) | 30 s | 10⁶ | 60 s |
| DRIFT | Velocity change > 10% since last map; enter/exit magnetosphere or atmosphere; reactor power change > 10%; BH mass crosses a table band | 300 s | 10⁶ | 300 s |
| PERIODIC | Absolute sim-time grid: every 86,400 s of coordinate time, aligned to multiples of 86,400 s from J2000.0 (catches gradual change during timewarp) | at grid tick | 10⁶ | — |

Between applications the cached rates are used (rate × dt extrapolation, deterministic).
Staleness is bounded by Δt_apply + trigger-detection granularity — physically acceptable for
every class (radiation environments change on longer timescales than their class delays).

## DETERMINISTIC SCHEDULING PROTOCOL (deferred-apply)

The wall clock never influences which tick a dose map becomes physics-visible on. The
protocol makes both the CONTENT and the APPLICATION TIME of every map a pure function of
sim state:

```
On trigger at coordinate time T_trig (detected in Phase 4, tick N):
  job = TransportJob {
      job_index:    monotonic counter (deterministic: increments per job created),
      class:        trigger class (table above),
      inputs:       snapshot at end of tick N — octree version, source set,
                    environment (velocity, position-derived fields, reactor/BH state),
      seed:         hash(universe_seed, "radiation", job_index),
      n_primaries:  fixed per class,
      T_apply:      T_trig + Δt_apply(class),
  }
  submit to background thread pool.

At every tick with start time ≥ job.T_apply (first such tick, before Phase 1):
  if job result ready:  swap dose-rate map atomically; job retired.
  else:                 BLOCK the tick loop until the result is ready, then swap.
                        (Blocking costs wall-clock time only — zero physics effect.)
```

Rules that keep the schedule deterministic:
- **One job in flight per class.** While a class has a pending job, new triggers of that
  class are suppressed; the condition is re-evaluated in the first Phase 4 after T_apply.
- **Same-tick collisions:** if multiple classes trigger on one tick, jobs are created in
  fixed class order (ACUTE, GEOMETRY, DRIFT, PERIODIC) so job_index assignment is deterministic.
- **Inputs are snapshots**, taken at the trigger tick — never read live by the worker.
- **Content determinism:** counter-based PRNG per primary + fixed-order tally merge
  (radiation_transport.md); same job → bit-identical map on any hardware and thread count.
- **Application determinism:** T_apply is sim time; tick boundaries are deterministic;
  therefore the map swap happens on the same tick in every run. A slow machine blocks
  longer at the swap point; it does not get a different trajectory.

Timewarp interaction: at high warp, PERIODIC jobs fall due faster in wall-clock terms than
the ~10-60 s Monte Carlo runtime, and the tick loop will block on them — the transport
throughput sets the effective WALL-CLOCK rate in radiation-active regimes. CRITICAL RULE:
granted dt NEVER depends on wall-clock. If compute or transport cannot keep up, the sim
simply runs slower than the requested rate — the tick/dt sequence is identical on any
hardware; only the wall-clock duration differs. (Auto-reducing warp on wall-time would make
dt hardware-dependent and break bit-determinism; any automatic warp adjustment must be
triggered by deterministic sim state — e.g., pending-job count at a tick boundary — never
by elapsed wall time.)

```rust
struct DoseRateMap {
    // ... per-voxel and per-sub-component rates as above ...
    job_index: u64,            // which job produced this map
    computed_from: DoubleDouble, // T_trig of the producing job (staleness bookkeeping)
    applied_at: DoubleDouble,    // T_apply — the sim time this map became physics-visible
}
```

(The previous `valid: bool` + wall-clock polling loop is deleted: map validity is now a
scheduling property, not a runtime flag.)

## DOSE TRACKING

### Crew dose

Accumulated in Sieverts. Quality factors convert absorbed dose (Gy) to dose-equivalent (Sv):
- Photons, electrons: Q = 1
- Protons: Q = 2-5
- Neutrons: Q = 5-20 (energy-dependent)
- Heavy ions: Q = 20-40

Dose limits (NASA-STD-3001):
- 30-day limit: 250 mSv
- Annual limit: 500 mSv (recently updated)
- Career limit: 600 mSv

### Electronics dose

Total ionizing dose (TID) accumulated in Gray. SEE rate accumulated in upsets/bit.

Failure thresholds depend on the device:
- Radiation-hardened (rad-hard): TID > 1-10 kGy (100 krad - 1 Mrad; note 100 krad = 1 kGy, NOT 100 kGy)
- Radiation-tolerant: TID > 0.3-1 kGy
- Commercial: TID > 0.05-0.3 kGy
- SEU rate: > ~10^-7 upsets/bit/day triggers error correction concern

### Material degradation

Displacement damage dose (dpa) accumulated per component. Feeds into material property evaluation:
- yield_strength(T, dpa) degrades with increasing dpa
- thermal_conductivity(T, dpa) degrades
- swelling(dpa, T) may change geometry

## TOOLS

| Tool | When | Purpose |
|---|---|---|
| radiation_transport crate | Background thread | Monte Carlo transport to produce dose-rate map |
| None at runtime per-tick | -- | Pure multiply-accumulate from cached map |
