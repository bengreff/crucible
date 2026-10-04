# Determinism Specification

Every computation in the simulator must be bit-exact reproducible. Same initial state → same trajectory → same outcome. Every time, every platform. This document specifies how determinism is guaranteed across all systems.

---

## THE GUARANTEE

Given the same (ship_initial_state, coordinate_time T, universe_seed), the simulation produces BIT-IDENTICAL results regardless of:
- What the player queried, in what order, at what time
- Which platform the code runs on (x86, ARM, etc.)
- Which run this is (first or thousandth)
- Compiler version (given the same Rust edition and flags)

---

## COORDINATE SYSTEM

**Galactocentric frame.** All positions stored as double-double (two f64s, ~31 significant digits) in meters from Sgr A*. Sub-nanometer precision at 30 kpc.

**Ship-local frame.** Octree voxel positions in plain f64 meters from ship center of mass. 22 femtometer precision at 100 m.

**Time.** Double-double seconds from J2000.0 (2000-01-01 12:00:00 TDB). Sub-nanosecond precision at T = 10 million years. Proper time accumulated separately: dτ = dt / γ per tick.

No frame switching. The galactocentric frame is used for all trajectory integration. The ship-local frame is used only for the octree. Conversion between them is a single offset (ship's galactocentric position), applied once per tick.

---

## FLOATING-POINT RULES

### 1. Deterministic math library

```rust
// NEVER: system libm (platform-dependent results for transcendentals)
let x = f64::sin(theta);

// ALWAYS: portable libm crate (bit-identical everywhere)
let x = libm::sin(theta);
```

The `libm` crate provides pure-Rust implementations of sin, cos, exp, log, sqrt, atan2, etc. Bit-reproducible across all platforms. Required for ALL transcendental function calls in the simulator.

IEEE 754 basic operations (+, -, ×, ÷) are already bit-exact across platforms with the same rounding mode (round-to-nearest-even, the Rust default).

### 2. No implicit FMA

Rust does not auto-insert FMA. Explicit `f64::mul_add(a, b, c)` is permitted but must be used CONSISTENTLY -- either always or never for a given computation. Mixed use (FMA on some platforms, separate mul+add on others) breaks determinism.

Policy: use `mul_add` explicitly where it improves accuracy (double-double algorithms). Never rely on the compiler to insert it.

### 3. Fixed summation order

Floating-point addition is non-associative: (a + b) + c ≠ a + (b + c).

Every sum over a variable-length collection must use a DETERMINISTIC ORDER independent of discovery order, query history, or hash table state.

```rust
// NEVER: sum over unsorted collection
let force: Vec3 = sources.iter().map(|s| gravity(s)).sum();

// ALWAYS: sort by permanent ID first
sources.sort_by_key(|s| s.id);
let force: Vec3 = sources.iter().map(|s| gravity(s)).sum();
```

This applies to:
- Active point source gravity summation (sort by source ID)
- Mean field basis function summation (fixed index order, 0..N)
- Radiation transport tallies (accumulate by voxel index order)
- Any reduction over generated stars (sort by cell_index + star_index)

### 4. No non-deterministic parallelism in physics

The per-tick physics loop is SINGLE-THREADED. At ~30-55 μs per tick (physics.md Section 7), parallelism is unnecessary and would introduce non-deterministic reduction ordering.

Parallel computation is permitted ONLY for:
- Radiation transport Monte Carlo (counter-based PRNG, see below). Its results DO feed
  physics (dose → dpa → material strength), so both its content AND its application time
  must be deterministic: jobs are triggered at deterministic sim times, computed from
  snapshotted inputs with fixed primary counts, and become physics-visible at a
  predetermined sim time T_apply — the tick loop blocks if the result is late. Wall-clock
  speed affects only how long the block lasts, never which tick the map lands on. Full
  protocol: radiation_shielding.md "Deterministic scheduling protocol."
- Background data loading (no physics impact)
- Player display queries (read-only, no physics impact)

For parallel Monte Carlo: use a counter-based PRNG (Philox). Each primary particle's RNG stream is seed + particle_index. Same particle always gets same random numbers regardless of thread assignment. Tallies accumulated in fixed voxel-index order after all particles complete.

### 5. Double-double arithmetic

The TwoSum and TwoProd algorithms are fixed sequences of f64 operations. Deterministic by construction. The `twofloat` crate or a custom implementation must use explicit operation ordering (no compiler reordering).

```rust
// TwoSum: exact sum of two f64s as (hi, lo)
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}
```

Each line is one f64 operation in a fixed sequence. The compiler must not reorder. In Rust, this is guaranteed for non-SIMD scalar operations at standard optimization levels.

---

## PROCEDURAL GENERATION DETERMINISM

### No shared mutable state

Every generated object's properties are derived from a hash-based seed, NOT from a sequential RNG that advances with each generation call.

```rust
fn generate_star(cell: CellIndex, index: u32, T: f64) -> Star {
    let seed = deterministic_hash(cell, index);  // pure function of coordinates
    let mut rng = Rng::from_seed(seed);           // independent per star
    
    let mass = sample_imf(&mut rng);
    let age = sample_age(&mut rng, population);
    // ... all properties from this independent RNG
}
```

Generating star A does not affect star B. Query region 1 first or region 2 first -- identical results. The RNG for each star is seeded from its spatial coordinates, not from a global counter.

### Positions are pure functions

Every star's position at coordinate time T is a deterministic function of (orbital_params, T):

```rust
fn star_position(params: &EpicyclicParams, T: f64) -> Position {
    let R = params.R_guide + params.A_R * libm::cos(params.kappa * T + params.phi_R);
    let phi = params.phi_guide + params.Omega * T   // phi_guide: azimuth of the guiding
              + (2.0 * params.Omega / params.kappa) // center at epoch (was missing)
              * params.A_R * libm::sin(params.kappa * T + params.phi_R) / params.R_guide;
    let z = params.A_z * libm::cos(params.nu * T + params.phi_z);
    galactocentric_from_cylindrical(R, phi, z)
}
```

No integration history. No persistent state. Same (params, T) always returns the same position. Uses `libm` for all trig.

### Hash function

The hash function mapping cell coordinates to seeds must be:
- Deterministic (same input → same output, always)
- Well-distributed (uniform coverage of the u64 seed space)
- Fast (~1-5 ns)
- Platform-independent (no pointer hashing, no random seeding)

Use a portable integer hash (e.g., SipHash with a fixed key, or xxHash with a fixed seed). Rust's `DefaultHasher` is SipHash but its seed varies by process. Use a FIXED-SEED hasher:

```rust
use std::hash::{BuildHasherDefault};
use siphasher::SipHasher13;

const UNIVERSE_SEED: u64 = 0x5A7E_F1C2_D3B4_A690;  // fixed, part of the universe definition

fn cell_seed(cell: CellIndex) -> u64 {
    let mut hasher = SipHasher13::new_with_keys(UNIVERSE_SEED, 0);
    cell.x.hash(&mut hasher);
    cell.y.hash(&mut hasher);
    cell.z.hash(&mut hasher);
    hasher.finish()
}
```

---

## PHYSICS PATH vs. DISPLAY PATH

Two completely separate data paths. The display path CANNOT affect the physics path.

```
PHYSICS PATH (deterministic, runs every tick):
  ┌─────────────────────────────────────────────┐
  │ 1. Evaluate mean field at ship position     │
  │ 2. Sum gravity from active sources (sorted) │
  │ 3. Compute thrust, drag, other forces       │
  │ 4. Integrate trajectory (DOP853)            │
  │ 5. Update ship internal state               │
  │ 6. Accumulate proper time                   │
  │ 7. Periodic refresh: update active sources  │
  └─────────────────────────────────────────────┘
  Input: ship state at tick N
  Output: ship state at tick N+1
  No player interaction. No query results. No display state.

DISPLAY PATH (player-driven, read-only):
  ┌─────────────────────────────────────────────┐
  │ Player queries a region                     │
  │ → Generate stars (stateless, independent)   │
  │ → Compute positions at current T            │
  │ → Apply filters                             │
  │ → Return results for display                │
  └─────────────────────────────────────────────┘
  Input: query parameters + current coordinate time
  Output: list of objects for display
  CANNOT write to: active source list, ship state, trajectory, any physics state.
```

The display path reads the universe (star seeds, positions at T) but never writes to physics state. A player looking at the star map, inspecting a planet, or querying a region has ZERO effect on the ship's trajectory.

---

## ACTIVE SOURCE LIST MANAGEMENT

The active source list determines which bodies contribute point-source gravity to the ship. It must be a deterministic function of (ship_position, coordinate_time T).

### Refresh procedure

```rust
fn refresh_active_sources(&mut self, ship_pos: &Position, T: f64) {
    // 1. Threshold from PREVIOUS tick's total acceleration (breaks circularity)
    let threshold = EPSILON * self.a_total_previous_tick;
    
    // 2. Scan spatial index for candidates
    //    The scan is deterministic: fixed octree structure, fixed traversal order
    let search_radius = max_influence_radius(threshold);
    let candidates = self.spatial_index.query_sphere(ship_pos, search_radius, T);
    
    // 3. Generate stars for all candidate cells (stateless, independent seeds)
    let mut sources: Vec<Source> = Vec::new();
    for cell in candidates.cells_in_order() {  // deterministic cell enumeration
        for star in cell.generate_stars(T) {     // deterministic per-cell generation
            let accel = star.gravitational_acceleration(ship_pos);
            if accel > threshold {
                sources.push(star.as_source());
            }
        }
    }
    
    // 4. Add solar system bodies (always active within solar system)
    sources.extend(self.solar_system_bodies(T));
    
    // 5. Sort by permanent ID for deterministic summation
    sources.sort_by_key(|s| s.permanent_id);
    
    // 6. Apply hard cap
    if sources.len() > MAX_ACTIVE_SOURCES {
        sources.sort_by(|a, b| b.acceleration.partial_cmp(&a.acceleration).unwrap());
        sources.truncate(MAX_ACTIVE_SOURCES);
        sources.sort_by_key(|s| s.permanent_id);  // re-sort by ID for summation
    }
    
    self.active_sources = sources;
}
```

The procedure is called every N ticks or when the ship moves more than X meters from the last refresh position. Both triggers are deterministic (depend only on the trajectory, which is deterministic).

### System entry

When a star crosses the threshold where it becomes the dominant non-mean-field source, Layer 2 generation triggers (planets, moons, etc.). The trigger is deterministic: same trajectory → same crossing time → same system generated at the same tick.

Generated system bodies (planets, moons) enter the active source list at the NEXT refresh. They don't retroactively affect the trajectory.

---

## REGIME TRANSITIONS

### 1PN to full GR (near compact objects)

Transition uses hysteresis to prevent oscillation:

```rust
const GR_ENTER_THRESHOLD: f64 = 100.0;   // Enter GR at 100 r_s
const GR_EXIT_THRESHOLD: f64 = 110.0;    // Exit GR at 110 r_s

fn update_gr_mode(&mut self, r: f64, r_s: f64) {
    if !self.gr_mode && r < GR_ENTER_THRESHOLD * r_s {
        self.gr_mode = true;
        // Switch integrator to Kerr geodesic equations
        // Convert state: (x, p) galactocentric → (t, r, θ, φ, p_t, p_r, p_θ, p_φ) Boyer-Lindquist
    }
    if self.gr_mode && r > GR_EXIT_THRESHOLD * r_s {
        self.gr_mode = false;
        // Switch back to the galactocentric weak-field geodesic law
        // (velocity-complete, physics.md Section 2.5) + Φ²/EIH corrections
        // Convert state: Boyer-Lindquist → (x, p) galactocentric
    }
}
```

The transition is deterministic: same trajectory → same r at each tick → same transition decision.

At the transition boundary (100 r_s), the 1PN and Kerr formulations agree to ~3×10⁻⁵ fractionally (residual terms are O((GM/rc²)²) = (1/200)² and O((v/c)⁴), both ~2.5×10⁻⁵ there — NOT percent-level). The residual step is handled by treating the transition as an integration-step boundary: the integrator stops exactly at the transition event and restarts with the new formulation, so no discontinuity ever appears inside a step. The accumulated trajectory error from the formulation change is bounded by the 1PN truncation error itself.

### Atmospheric / vacuum transition

Atmospheric drag activates when altitude < atmosphere_top. The atmosphere model (NRLMSISE or equivalent) returns density = 0 above the exobase. No discontinuity -- drag force goes to zero smoothly.

---

## TRAJECTORY INTEGRATOR DETERMINISM

The DOP853 adaptive integrator computes an error estimate at each step and adjusts dt accordingly. This is deterministic IF:
- Force evaluation is deterministic (guaranteed by sorted summation, libm, no parallelism)
- Error tolerance is fixed (atol, rtol set at initialization, never modified)
- Error norm uses a deterministic reduction (max over fixed-order iteration)

The integrator state is:
```rust
struct IntegratorState {
    position: Position,        // double-double
    momentum: Vec3,            // f64 (relativistic 3-momentum)
    proper_time: (f64, f64),   // double-double
    rest_mass: f64,
    dt: f64,                   // current adaptive timestep
    // DOP853 internal state (Nordsieck vector, error history)
}
```

All fields are deterministic functions of the initial state and the force history. No randomness, no external input.

---

## RADIATION TRANSPORT DETERMINISM

The Monte Carlo transport is inherently random, but deterministically so:

```rust
fn run_transport(octree: &Octree, sources: &[&dyn RadiationSource], n_primaries: u64, 
                 base_seed: u64) -> Tallies {
    let mut tallies = Tallies::new(octree);
    
    for i in 0..n_primaries {
        // Each primary gets a unique, deterministic RNG stream
        let particle_seed = deterministic_hash(base_seed, i);
        let mut rng = Rng::from_seed(particle_seed);
        
        let primary = sources.sample_particle(&mut rng);
        transport_particle(octree, primary, &mut rng, &mut tallies);
    }
    
    tallies  // accumulated in voxel-index order, deterministic
}
```

Each primary particle's entire transport history (interactions, secondaries, energy deposition) is determined by its particle_seed. Particle 0 always gets the same trajectory. Particle 999,999 always gets the same trajectory. The total tallies (sum over all particles in voxel-index order) are deterministic.

For parallel execution: partition particles into chunks, each chunk processed on one thread. Each particle still uses its particle_seed-derived RNG. Tallies from each thread are merged in fixed chunk order (not completion order). Result is identical to single-threaded execution.

---

## WHAT CANNOT AFFECT THE PHYSICS

| Action | Physics impact | Why |
|---|---|---|
| Player queries a region | None | Display path only, read-only |
| Player inspects a body | None | Generates Layer 3 detail for display, not for physics |
| Player queries in different order | None | Stateless generation, independent seeds |
| Player queries at different wall-clock time | None | Sim-time is independent of wall-clock time |
| Running on faster/slower hardware | None | dt is sim-time, not wall-clock time |
| Background radiation-transport completion timing | None | Deferred-apply protocol: results bind at predetermined sim time T_apply; tick loop blocks if late (radiation_shielding.md) |
| Running on different platform (x86 vs ARM) | None | libm + IEEE 754 + no FMA = bit-exact |
| Player issues a command (throttle, orientation) | Yes (intentional) | This is input that changes the ship state |
| Player advances time (time jump) | Yes (intentional) | This changes coordinate time T |

Only explicit player commands (throttle, orientation, time) affect the physics. Everything else is observation-only.

---

## TESTING DETERMINISM

### Regression test

Run a reference trajectory (Falcon 9 launch, 10 minutes, ~600,000 ticks). Record the final ship state (position, velocity, fuel mass, proper time) to full f64 precision. On every commit, re-run and compare. Any bit-level difference is a determinism regression.

### Cross-platform test

Run the same reference trajectory on x86_64 and aarch64. Compare final states. Must be bit-identical.

### Query-independence test

Run a trajectory twice:
- Run 1: player queries nothing
- Run 2: player queries 1000 random regions during the flight

Compare final ship states. Must be bit-identical.

### Order-independence test

Generate all stars within 100 pc. Then generate them again, starting from a different region. Compare all star properties. Must be bit-identical.
