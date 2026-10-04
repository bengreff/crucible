# Procedural Generation Architecture

Complete specification for generating the galaxy's ~130 billion stars (a count DERIVED from the mass model — physics.md Section 3.2.9), their planetary systems, biospheres, minor bodies, and environmental fields. Everything is deterministic from seed. Generation is layered -- cheap summary data for galaxy-scale queries, full detail on demand.

**Simulation timescale:** 10 million years from J2000.0 epoch. Sufficient to cross the galaxy at 0.01c. All generation must produce consistent results across this window.

**Time epoch:** J2000.0 (2000-01-01 12:00:00 TDB). Gaia DR3 catalog (epoch 2016.0) is propagated back 16 years to J2000.0 using measured proper motions and radial velocities.

**Determinism:** All generation uses independent per-object hash-derived seeds -- no shared mutable RNG state. Generation order cannot affect results. All transcendentals via `libm` for cross-platform bit-exactness. Player queries are read-only and cannot affect physics. See [determinism.md](determinism.md) for the complete specification.

---

## 1. DESIGN PRINCIPLES

**Determinism:** Every property of every object is a pure function of its seed and the coordinate time T. No runtime state, no integration history, no randomness outside the seed. Visiting a system at T=0, leaving, and returning at T=5,000,000 years produces the same system in a state consistent with 5 million years of deterministic evolution.

**Layered detail:** The galaxy map draws on ~130 billion stars at ~55 ns each. Planetary systems are generated only when inspected (~1-10 ms). Geology, atmosphere, and biosphere only on body inspection (~0.1-1 ms). Surface physics only on approach. Each layer adds detail without regenerating previous layers.

**Physical realism within scope:** Every generated property uses the best available observational distributions and physical models. Star masses from the IMF, planet occurrence from Kepler, biospheres from abiogenesis probability research. Where models are uncertain, we use the median estimate and acknowledge the uncertainty.

**Omniscient observer:** The player sees all objects regardless of detectability. No magnitude limits, no survey bias, no detection thresholds. Every brown dwarf, rogue planet, and cold neutron star is visible when queried. This is a physics simulator, not a detection simulator.

**Current state only:** The player sees the state of objects at the current coordinate time T. Future events (supernovae, flares) are not predicted -- they happen when T reaches their seed-determined time. The player observes and infers, like a real observer with perfect instruments.

---

## 2. GENERATION LAYERS

| Layer | Triggered by | Per-object cost | Data generated |
|---|---|---|---|
| **0: Galaxy structure** | Precomputed | N/A | Mean field potential, density model, cluster orbits |
| **1: Star catalog** | Spatial query or nearby-star refresh | ~55 ns | Mass, age, spectral type, luminosity, position at T, velocity |
| **2: Star system** | System inspection or ship entering system | ~1-10 ms | Planets, moons, asteroid belts, rings, stellar wind, system architecture |
| **3: Body detail** | Body inspection | ~0.1-1 ms | Geology, atmosphere, magnetic field, rotation, weather, biosphere |
| **4: Surface/environment** | Ship approach (atmospheric flight range) | ~1 ms | Atmospheric profile, radiation environment, magnetic field topology |

Each layer is generated from the same seed chain. Layer N includes all data from layers 0 through N-1.

---

## 3. STAR CATALOG GENERATION (Layer 1)

### 3.1 Spatial Indexing

The galaxy is divided into generation cells using a log-radius cylindrical grid in (log R, φ, z) coordinates, where R is galactocentric radius. Log-radius naturally gives smaller cells near the galactic center (high density) and larger cells in the outer disk (low density).

Base grid:
- log₁₀ R: 100 bins from -3 to 2 (0.001 pc to 100 kpc)
- φ: 3600 bins (0.1° per bin)
- z: 200 bins from -10 kpc to +10 kpc (100 pc per bin, finer near midplane)

Cells that exceed ~1000 stars after density evaluation get octree-subdivided until each sub-cell contains ≤ 1000 stars.

Within Sgr A*'s sphere of influence (~2 pc): cells are defined in (a, e) orbital element space instead of (R, φ, z). Semi-major axis in 50 logarithmic bins from 0.001 to 2 pc. Eccentricity in 20 linear bins from 0 to 1.

Each cell's seed is a hash of its integer grid indices. Same cell always produces the same stars.

**Query mapping (guiding-center unwinding) + drift buffer:** Disk stars drift Ω(R_g)·T in
azimuth — ~2.4 kpc of arc at R₀ over 10 Myr. A query at position (R, φ) at time T must
therefore be UNWOUND into guiding-center space before cell lookup:

    φ_g,query(R_bin) = φ − Ω(R_bin) · T     (evaluated per radius bin, since Ω varies with R)

then expand the unwound search region by the residual excursions:
- Disk stars: epicyclic amplitude A_R (typically 0.1-1 kpc thin disk) + cell shear
  ΔΩ(across the cell's R extent)·T·R + dispersion drift (σ_φ·T ~ 225 pc over 10 Myr)
- Galactic center stars: apocenter - pericenter range (can be large for eccentric orbits)

Without the unwinding term the query misses the bulk of the drift entirely (the 225 pc
dispersion buffer is only the SPREAD around the ~2.4 kpc systematic shift).

### 3.2 Star Count Per Cell

The number of stars generated in a cell is:

N_cell = ∫_cell n_stars(R, φ, z) dV

where n_stars is the total stellar number density:

n_stars(R, φ, z) = n_thin(R, z) + n_thick(R, z) + n_bulge(r) + n_halo(r)

Each component uses the density profiles from physics.md Section 3.2:

**Thin disk:** n_thin(R, z) = n_0_thin × exp(-R/R_d) × exp(-|z|/h_z), with R_d = 2.50 kpc, h_z = 0.30 kpc.

**Thick disk:** n_thick(R, z) = n_0_thick × exp(-R/R_d) × exp(-|z|/h_z), with R_d = 3.02 kpc, h_z = 0.90 kpc.

**Bulge:** n_bulge(r) from the power-law + exponential cutoff profile (r_cut = 2.1 kpc).

**Halo:** n_halo(r) ∝ r^(-3.5) (power-law, from the stellar halo density profile).

The normalization constants n_0 are NOT free parameters: n_0 = ρ_component(from the mass model, physics.md Section 3.2.1) / ⟨m⟩_population. The total count is then an output — ~1.3×10¹¹ stars (physics.md Section 3.2.9) — guaranteeing that the density making stars is exactly the density sourcing the potential.

**Spiral arm modulation (material arms, frozen at epoch):** Cell counts must be
TIME-INVARIANT — "same cell always produces the same stars" (Section 3.1) forbids any
T-dependence in N_cell. The arm enhancement is therefore written into the cells in
guiding-center coordinates, evaluated once at the J2000.0 epoch:

N_cell = ∫_cell n_disk(R_g, z) × [1 + (f_arm − 1) × S(R_g, φ_g − φ_arm,0(R_g))] dV

where φ_arm,0(R_g) is the logarithmic-spiral locus (physics.md Section 3.2.3) at T = 0, and
f_arm ~ 2-3 for old stars, ~5-10 for young/massive stars (which are also tightly concentrated:
narrower S). No T appears — determinism and order-independence hold exactly.

Stars inherit the enhancement through their guiding centers, so at time T the material
pattern sits at φ_arm,0(R_g) + Ω(R_g)·T: it shears with differential rotation rather than
rotating rigidly at Ω_spiral like the potential's Cox-Gomez pattern. **Bounded mismatch:**
the drift is (Ω(R_g) − Ω_spiral)·T ≈ 3° at R₀ and up to ~14° (~1.2 kpc) in the inner disk at
the 10 Myr window edge — below the arm FWHM (0.5-1 kpc) at the solar circle, comparable to it
inside ~5 kpc. Acknowledged limitation (Section 14.3); no fixed-parameter kinematic scheme
produces an exactly rigid pattern at all radii, so this error is bounded, not eliminable.
For young stars the material picture is also the physically correct one (age ≪ shear time:
they are born in the arm and have not left it within the window).

### 3.3 Stellar Property Assignment

For each star, the seed generates a chain of properties:

**Step 1: Population assignment.** Based on the cell's galactic position, assign the star to a population:
- Thin disk (probability proportional to n_thin / n_total at cell center)
- Thick disk (n_thick / n_total)
- Bulge (n_bulge / n_total)
- Halo (n_halo / n_total)

This determines age distribution, metallicity, and velocity dispersion.

**Step 2: Mass.** Draw from the Kroupa/Chabrier IMF:

dN/dm ∝ m^(-α), where:
- α = 0.3 for 0.01 ≤ m/M_sun < 0.08 (brown dwarfs, sub-filter by default)
- α = 1.3 for 0.08 ≤ m/M_sun < 0.5
- α = 2.3 for 0.5 ≤ m/M_sun

The IMF is sampled via inverse CDF. ~75% of stars are M-dwarfs (0.08-0.5 M_sun). ~10% are K-dwarfs. ~3% are G-dwarfs. Stars > 8 M_sun are ~0.3% of the population.

Brown dwarfs (0.013-0.08 M_sun) are generated from the sub-stellar tail of the IMF. Rogue planets are NOT an IMF tail draw — the measured abundance (~20 per star, Sumi+ 2023 MOA microlensing, rising steeply toward sub-Earth masses) far exceeds any IMF extrapolation because most are ejected from planetary systems; they get their own occurrence model (Section 8.2). Both are hidden in default queries and shown with an explicit "substellar objects" filter.

**Step 3: Age.** From the star formation history of the host population:
- Thin disk: roughly uniform 0-10 Gyr, with slight increase toward older ages
- Thick disk: 8-12 Gyr (formed early, little recent star formation)
- Bulge: 8-13 Gyr (mostly old, some younger populations in the nuclear disk)
- Halo: 10-13 Gyr (oldest population)

Massive stars (> 2 M_sun) in the thin disk are preferentially young (concentrated in spiral arms, star-forming regions). The seed determines whether a massive star is young enough to still be on the main sequence.

**Step 4: Metallicity.** From the age-metallicity relation and galactic gradient:
- [Fe/H] = [Fe/H]_0 - 0.06 × (R - R_0) / kpc + 0.04 × age / Gyr
- Scatter: σ_[Fe/H] ~ 0.15 dex (drawn from seed)
- [Fe/H]_0 ~ 0.0 at the solar position, ~+0.3 in the inner disk, ~-0.5 in the outer disk
- Halo: [Fe/H] ~ -1.5 ± 0.5 (metal-poor)

**Step 5: Evolutionary state.** Compare age to main-sequence lifetime:

t_MS ≈ 10 Gyr × (M/M_sun)^(-2.5)

| Condition | State | Properties |
|---|---|---|
| age < t_MS | Main sequence | L, T_eff, R from mass via power-law relations |
| age ≈ t_MS (within 10%) | Subgiant / RGB | Expanding, cooling. L increases, T_eff decreases. |
| age > t_MS, M < 8 M_sun | White dwarf | M_WD = 0.08 × M_initial + 0.48 M_sun (Cummings+ 2018; the 0.39 intercept underpredicted — 1 M_sun → 0.47 vs. observed ~0.55). T_eff from cooling age: T ∝ t^(-0.35) (Mestel's t^(-1.4) is the LUMINOSITY exponent; L ∝ T⁴ ⇒ T exponent ≈ -0.35). |
| age > t_MS, 8 < M < 25 M_sun | Neutron star | M_NS ~ 1.4 M_sun. Kick velocity 200-500 km/s applied to orbit. |
| age > t_MS, M > 25 M_sun | Black hole | M_BH from fallback prescriptions. Kick ~50-100 km/s. |

**Main-sequence relations (used for current L, T_eff, R):**

| Mass range | L/L_sun | T_eff (K) | R/R_sun |
|---|---|---|---|
| 0.08-0.45 | 0.23 × (M/M_sun)^2.3 | 3000 + 2100×(M/M_sun) | 0.91 × (M/M_sun)^0.8 |
| 0.45-2 | (M/M_sun)^4.0 | 5800 × (M/M_sun)^0.65 | (M/M_sun)^0.8 |
| 2-55 | 1.4 × (M/M_sun)^3.5 | 5800 × (M/M_sun)^0.6 | 1.3 × (M/M_sun)^0.55 |
| > 55 | 32000 × (M/M_sun)^1.0 (near-Eddington) | 40000-50000 | ~15-25 |

(Branch boundary at ~55 M_sun, where the two power laws cross — the previous boundary at
20 M_sun made 20-55 M_sun stars discontinuously ~13× overluminous.)

These are approximate fits. More precise: use pre-tabulated MIST/MESA isochrones indexed by (mass, age, [Fe/H]).

**Spectral type:** From T_eff:

| T_eff range (K) | Spectral type |
|---|---|
| > 30,000 | O |
| 10,000-30,000 | B |
| 7,500-10,000 | A |
| 6,000-7,500 | F |
| 5,200-6,000 | G |
| 3,900-5,200 | K |
| 2,400-3,900 | M |
| 1,300-2,400 | L (brown dwarf) |
| 600-1,300 | T (brown dwarf) |
| < 600 | Y (brown dwarf) |

Subtype (0-9) interpolated within each class.

**Step 6: Velocity and orbital parameters.** Assign velocity from the rotation curve + population-dependent velocity dispersion (physics.md Section 3.2.8):

v_R = draw from Gaussian(0, σ_R)
v_φ = V_circ(R) - v_asymmetric_drift + draw from Gaussian(0, σ_φ)
v_z = draw from Gaussian(0, σ_z)

Convert to epicyclic orbital parameters:
- R_guide from angular momentum: R_guide × V_circ(R_guide) = R × v_φ
- A_R from radial energy: A_R ≈ v_R / κ(R_guide)
- A_z from vertical energy: A_z ≈ v_z / ν(R_guide)
- Phases from position at epoch

For neutron stars: add natal kick (200-500 km/s in random direction) to the progenitor's velocity before computing orbital parameters. This can put NS on highly eccentric or unbound orbits.

For black holes: smaller kick (~50-100 km/s).

### 3.4 Stellar Activity and Events

**Rotation period:** From gyrochronology (Skumanich relation): P_rot ∝ age^0.5. Young stars rotate fast (P ~ 1-5 days). Sun-like at solar age: ~25 days. M-dwarfs: slow spindown, many retain fast rotation.

**Magnetic activity:** From rotation period and spectral type. Fast rotators are more active. Activity index R'_HK correlates with rotation. Determines flare rate, spot coverage, chromospheric emission.

**Flare schedule (seed-determined):**

Flare rate depends on spectral type and activity:
- Active M-dwarf: ~1 major flare per day
- Quiet M-dwarf: ~1 per month
- Active G-star: ~1 per week
- Quiet G-star (Sun-like): ~1 per month
- F/A/B stars: minimal flaring

Flare times use HIERARCHICAL INTERVAL SEEDING, not a cumulative sum. (A cumulative sum of
inter-flare intervals is O(N) to evaluate at arbitrary T — an active M-dwarf accumulates
~3.6×10⁹ draws over 10 Myr, making "is it flaring at T?" unanswerable. A Poisson process is
independent across disjoint intervals, so it can be seeded per-interval with no loss of
statistical correctness.)

    1. Time is divided into fixed buckets aligned to the J2000.0 epoch. Bucket length is
       chosen per star from its flare rate so the expected count per bucket is ~1-10
       (active M-dwarf: 1-day buckets; quiet G star: 30-day buckets; length is a pure
       function of the star's seed-derived rate).
    2. Per-bucket seed: hash(star_seed, "flare", bucket_index).
    3. Count in bucket: Poisson(rate × bucket_length) drawn from the bucket seed.
    4. Flare times uniform within the bucket, energies from the power-law, durations from
       the energy — all from the same bucket seed.

Evaluation at time T touches only the O(1) buckets overlapping the query window, at any T in
the 10 Myr range, in any order — deterministic and order-independent by construction. The
same scheme applies to every recurring stochastic event stream: SPEs at other stars,
microflaring, discrete mass-loss events.

Flare energies from the observed power-law: dN/dE ∝ E^(-1.8). Range: 10^24 J (microflare) to 10^33 J (superflare). Seed determines each flare's energy.

**Supernova timing:** For stars with M > 8 M_sun, the supernova time is:

t_SN = t_MS + δ_SN

where δ_SN is a seed-determined offset (~0-10,000 years, representing the duration of the final evolutionary stages). If the current coordinate time T falls between t_SN and t_SN + t_remnant_visible (~100,000 years for SNR visibility), the star is replaced by a neutron star or black hole surrounded by an expanding supernova remnant (see Section 9.2).

**Variable stars:** Stars in the instability strip (Cepheids, RR Lyrae, Mira variables) are identified by mass + evolutionary state. Pulsation period from the period-luminosity relation. Phase at time T from the seed. Luminosity varies sinusoidally (simplified) with the pulsation period.

### 3.5 Binary and Multiple Systems

When a star is generated, the seed determines its multiplicity.

**Binary fraction** depends on primary mass and environment:

| Primary mass | Solar neighborhood binary fraction | Source |
|---|---|---|
| O/B (> 8 M_sun) | ~70% | Sana+ 2012 |
| A (2-8 M_sun) | ~50% | De Rosa+ 2014 |
| FGK (0.5-2 M_sun) | ~46% | Raghavan+ 2010 |
| M (0.08-0.5 M_sun) | ~26% | Winters+ 2019 |

In denser environments, the binary fraction decreases according to the encounter survival limit (physics.md Section 3.3). The maximum binary separation is clipped to the local limit. If the drawn separation exceeds the limit, the system is generated as a single star.

**Companion mass ratio:** For solar-type primaries: roughly flat distribution in q = M2/M1 from 0.1 to 1.0 (Raghavan+ 2010). For massive stars: strong preference for q → 1 (twin excess). For M-dwarfs: roughly flat.

**Orbital period:** Log-normal distribution with:
- Solar-type: peak at log P(days) ≈ 5.0 (P ~ 300 years), σ ≈ 2.3
- M-dwarfs: peak at log P ≈ 4.5 (P ~ 90 years)
- Massive stars: bimodal (tight + wide)

**Eccentricity:** Thermal distribution f(e) = 2e for wide binaries. Circularized (e ~ 0) for close binaries with P < 10 days (tidal circularization).

**Triple and higher-order systems:** ~10% of binaries are triples. Generated as hierarchical: close inner pair + distant outer companion. Inner period << outer period for stability (P_outer / P_inner > 4.7 for stability, Mardling & Aarseth 2001).

The binary system is one atomic generation unit. Both stars share the same galactic orbital parameters (center of mass). Individual positions at time T = CM galactic position + Kepler offset from the mutual orbit.

---

## 4. STAR SYSTEM GENERATION (Layer 2)

Triggered when the player inspects a star or the ship's gravitational threshold activation reaches the system.

### 4.1 Architecture Overview

A star system contains:
- Host star(s) (from Layer 1, plus any binary companions)
- Planets (rocky, sub-Neptune, gas giant, ice giant)
- Moons (per planet)
- Asteroid belts (density fields with individually generated members > 1 km)
- Ring systems (density fields around gas/ice giants)
- Debris disk / Kuiper belt analog
- Oort cloud (extent from environment limit)
- Stellar wind structure
- Any active comets

Everything is deterministic from the star's seed.

### 4.2 Planet Occurrence Rates

From Kepler, K2, TESS, and radial velocity surveys:

| Planet type | Radius range | Period range | Occurrence per FGK star | Source |
|---|---|---|---|---|
| Hot Jupiters | > 6 R_earth | < 10 days | 0.5-1% | Wright+ 2012 |
| Warm Jupiters | > 6 R_earth | 10-200 days | 1-3% | Petigura+ 2018 |
| Cold Jupiters | > 6 R_earth | 200d - 20 yr | 10-15% | Wittenmyer+ 2020, RV surveys |
| Sub-Neptunes | 1.7-3.5 R_earth | 1-400 days | 30-50% | Fulton+ 2017, Hsu+ 2019 |
| Super-Earths | 1.0-1.7 R_earth | 1-400 days | 20-30% | Hsu+ 2019 |
| Earth-sized | 0.7-1.0 R_earth | 50-400 days | 10-20% | Bryson+ 2021 |
| Long-period giants | > 0.3 M_jup | 1-20 yr | 15-20% | Cumming+ 2008, Fernandes+ 2019 |

For M-dwarfs: higher occurrence of small planets (~2-3× for super-Earths), lower occurrence of gas giants (~3-5× lower). For A/F stars: higher occurrence of gas giants.

**Metallicity correlation:** Giant planet occurrence increases strongly with host star metallicity: P(giant) ∝ 10^(2×[Fe/H]) (Fischer & Valenti 2005). Metal-poor stars rarely host gas giants.

**Total planets per system:** Average ~2-3 detectable Kepler-size planets per star. Including undetected long-period and small planets: estimated 5-10 per system. Systems with gas giants tend to have fewer small inner planets (giant migration disruption).

### 4.3 System Architecture Generation

The generation pipeline:

**Step 1: Determine giant planet presence.** Roll for cold Jupiter, warm Jupiter, hot Jupiter based on occurrence rates modified by metallicity and stellar mass. Giant planets are placed first because they dominate system architecture.

**Step 2: Place giants.** If present, draw semi-major axis from the observed distribution (peaked at ~1-3 AU for cold Jupiters). Draw mass from the planetary mass function (dN/dm ∝ m^(-1.1) for 0.3-10 M_jup, Cumming+ 2008). Draw eccentricity from the observed distribution (Beta distribution with mean e ~ 0.25 for single giants, higher for multi-giant systems).

**Step 3: Place inner system.** Generate small planets (super-Earths, sub-Neptunes, Earth-sized) in the region interior to any giant. Follow the "peas in a pod" pattern (Weiss+ 2018):
- Adjacent planets have similar radii (within ~30%)
- Spacing: 10-30 mutual Hill radii between adjacent pairs
- R_Hill = a × (M_planet / 3M_star)^(1/3)
- Draw number of inner planets from Poisson distribution with mean ~3-5
- Place from inside out, ensuring stability

**Step 4: Place outer system.** Beyond the giant(s), generate ice giants and distant rocky/icy bodies:
- Ice giant occurrence: ~10-20% per system
- Kuiper belt analog: ~50% of systems, extent from ~30-100 AU scaled by stellar mass

**Step 5: Stability check.** Verify that all adjacent planet pairs are separated by > 10 mutual Hill radii. Remove any planet that violates stability (it would have been ejected during formation).

**Step 6: System extent.** Clip the outer boundary of the system to the local environment maximum (encounter survival limit from physics.md Section 3.3).

### 4.4 Known Exoplanet Systems (Augmentation, Not Override)

For stars with confirmed exoplanets from the NASA Exoplanet Archive (~5,500 systems as of 2025):

1. Load the known planets with their measured orbital elements, masses, and radii.
2. Identify dynamically stable gaps between known planets where additional planets could exist.
3. Procedurally generate planets in those gaps using the same architecture rules (Step 3-5 above), but constrained to not conflict with the known planets.
4. Generate moons, rings, and asteroid belts for both known and procedural planets.

This means visiting TRAPPIST-1 gives you the 7 known planets with their real properties, PLUS any procedurally generated bodies in the gaps (moons, asteroids, outer ice giants beyond the known system).

### 4.5 Planet Properties (Summary Level)

For each planet, the seed determines:

**Radius:** Drawn from the observed radius distribution within the planet's type category.

**Mass:** Direction of inference depends on regime — the Chen & Kipping (2017) relation is
DEGENERATE in the Jovian regime (R ≈ const from 0.3 to 13 M_jup), so mass-from-radius is
ill-defined there and would assign every giant ~0.3 M_jup:

| Regime | Rule | Scatter |
|---|---|---|
| Giants (> 6 R_earth) | MASS FIRST from the planetary mass function dN/dm ∝ m^(-1.1), 0.3-13 M_jup (Cumming+ 2008; already drawn in Step 2 of Section 4.3) → radius from the Jovian M-R relation with inflation for hot Jupiters | ~0.1 in log R |
| Neptunian (1.23-14.3 R_earth, non-giant) | Radius first → M = 2.0 × (R/R_earth)^1.7 M_earth | ~0.25 dex |
| Terran (< 1.23 R_earth) | Radius first → M = (R/R_earth)^3.7 M_earth | ~0.15 dex |

**Bulk density:** From mass and radius: ρ = 3M / (4π R³).

**Composition class:**
- ρ > 4 g/cm³: iron-rich rocky (Mercury-like)
- 3-5 g/cm³: silicate rocky (Earth/Venus-like)
- 1.5-3 g/cm³: icy/volatile-rich (ocean world or sub-Neptune with envelope)
- 0.5-1.5 g/cm³: gas/ice giant (H/He dominated)
- < 0.5 g/cm³: inflated hot Jupiter

**Orbital elements:** Semi-major axis, eccentricity, inclination (from Step 3-5). Argument of pericenter and mean anomaly from seed.

**Equilibrium temperature:** T_eq = T_eff_star × sqrt(R_star / (2 × a)) × (1 - albedo)^0.25. Default albedo: 0.3 for rocky, 0.5 for gas giants.

### 4.6 Moon Systems

**Gas giants:** Each gas giant generates a moon system:
- Number of major moons (M > 10^20 kg): 1-8, scaling with planet mass. Jupiter analog: ~4 Galilean-class moons. Saturn analog: ~1-2 Titan-class.
- Number of minor moons: 10-80, scaling with planet mass.
- Major moon masses: drawn from the satellite-to-planet mass ratio distribution (~10^(-4) total moon mass / planet mass for regular satellites, Canup & Ward 2006).
- Orbital radii: placed from ~5 to ~60 planet radii. Inner moons closer. Regular satellites follow Laplace resonance patterns (like Io-Europa-Ganymede).
- Compositions: inner moons rocky/volcanic, outer moons icy. Transition at the planet's ice line analog.

**Rocky planets:** 0-2 moons.
- Earth-mass in HZ: ~50% chance of a large moon (from giant impact, like Earth-Moon). Moon mass ~1-5% of planet mass.
- Mars-sized: ~30% chance of 1-2 small moons (captured asteroids).

**Tidal heating:** Moons with non-zero eccentricity orbiting within ~15 planet radii experience tidal heating:

Q_tidal = (21/2) × (k_2 / Q) × (n^5 × R_moon^5 × e^2) / G

where k_2 is the Love number (~0.3 for icy, ~0.1 for rocky), Q is the tidal quality factor (~100 for ice, ~10-50 for rock with liquid interior), n is the orbital mean motion, R_moon is the moon radius, and e is orbital eccentricity.

Compare Q_tidal to the moon's radiogenic heat production (~5 × 10^(-12) W/kg for silicates). If Q_tidal dominates: the moon is tidally heated (Io analog if rocky, Europa analog if icy with subsurface ocean).

Tidally heated icy moons with subsurface oceans are candidate habitable environments (assessed in Section 6).

### 4.7 Asteroid Belts

Each system gets 0-3 asteroid belts from the architecture:
- **Inner belt** (between rocky and giant planets, like the Main Belt): present if there's a giant that prevented planet formation in that zone. Total mass ~10^(-4) to 10^(-3) M_earth.
- **Outer belt / Kuiper analog** (beyond the outermost planet): present in ~50% of systems. Mass ~0.01-0.1 M_earth.
- **Trojan populations** at L4/L5 of gas giants: ~10% of giant planet mass in Trojans (scaled from Jupiter's).

Belt properties:
- Inner/outer radius, surface density profile (typically ∝ r^(-1.5))
- Size-frequency distribution: cumulative N(>D) ∝ D^(-2.5) (collisional equilibrium)
- Total number > 1 km: estimated from total mass and SFD. Solar Main Belt: ~10^6 asteroids > 1 km. Scales with belt mass.
- Composition: inner belts rocky (S/C-type), outer belts icy (C/D-type)

Individual asteroids > 1 km are generated on query from belt_seed + size_rank_index:
- Size from the SFD (inverse CDF sampling)
- Density: 1.5-5 g/cm³ (from composition type)
- Shape: triaxial ellipsoid with axis ratios a:b:c drawn from observed distribution (mean 1:0.8:0.6)
- Rotation period: from the observed spin distribution (2-20 hours for D > 1 km; spin barrier at ~2.2 hours for rubble piles)
- Orbital elements within the belt boundaries
- Designation: `{star_id}-Belt{letter}-{size_rank}` (deterministic from seed)

### 4.8 Ring Systems

Gas and ice giants: ~50% have visible ring systems (drawn from seed).

Ring properties:
- Inner radius: from the Roche limit: r_Roche = 2.46 × R_planet × (ρ_planet / ρ_ring_particle)^(1/3)
- Outer radius: 1.5-3 × R_planet for substantial rings (Saturn-like), up to 5 × R_planet for tenuous
- Optical depth profile τ(r): stepped/structured with gaps at resonances with major moons
- Gap locations: r_gap where n_particle / n_moon = p/q (integer ratio, strongest at 2:1, 3:2)
- Composition: ice (Saturn-like), rock/dust (Uranus-like)
- Particle size distribution: 1 cm to 10 m (power law, dN/dD ∝ D^(-3))

Rings are continuous density fields, not individual particles. The ship passing through experiences:
- Particle flux: F = n_particles × v_rel × σ_cross_section
- Micrometeorite impact rate from the particle density
- Optical depth affects illumination (shadow of the planet, ring shadow on planet)

### 4.9 Stellar Wind

Within each star system, the stellar wind creates a bubble (astropause) where the wind dominates over the ISM.

**Wind parameters at distance r from the star:**

| Property | Formula | Solar value at 1 AU |
|---|---|---|
| Proton density | n(r) = M_dot / (4π r² v_wind m_p) | ~5 cm^(-3) |
| Velocity | v_wind (constant beyond ~10 R_star) | ~400 km/s |
| Temperature | ~10^5-10^6 K (coronal) | ~10^5 K |
| Magnetic field | B_r ∝ r^(-2), B_φ ∝ (Ω_star r sin θ) / v_wind × r^(-1) | ~5 nT |

**Mass loss rate M_dot** scales with stellar type:
- Sun-like: ~2 × 10^(-14) M_sun/yr
- Active young star: 10-100× solar
- M-dwarf: ~0.1-1× solar (but concentrated, stronger per unit area)
- Red giant: 10^(-8) to 10^(-6) M_sun/yr (massive winds)
- O/B star: 10^(-6) to 10^(-4) M_sun/yr

**Astropause radius** (where wind ram pressure = ISM pressure):

R_astropause = sqrt(M_dot × v_wind / (4π × n_ISM × m_p × v_ISM²))

Solar value: ~120 AU. Scales with sqrt(M_dot × v_wind / n_ISM).

Inside the astropause: the ship is in the stellar wind. Outside: the ISM.

### 4.10 Orbital Evaluation Over the 10 Myr Window

ALL orbital motion is analytic: every body's position at coordinate time T is a pure function
of (elements, T). No runtime N-body, no integration history — required by determinism
(query-independence: a system's state at T must not depend on when it was first generated).

**Planets:** Kepler propagation of osculating elements + first-order Laplace-Lagrange secular
theory. At generation time (Layer 2), the secular apsidal and nodal rates ϖ̇_i, Ω̇_i are
computed deterministically from the system's generated masses and semi-major axes (standard
L-L coefficients; the stability step in Section 4.3 already enforces the non-resonant spacing
that makes first-order theory applicable). Position at T = Kepler(a, e, i, ϖ₀ + ϖ̇T,
Ω₀ + Ω̇T, M₀ + nT). Near-resonant configurations get the same formulas — an acknowledged
first-order approximation (Section 14.3).

**Moons:** Kepler around the planet + nodal/apsidal precession from the planet's J2 (derived
from its rotation and oblateness at Layer 3) and the stellar tide.

**Host star reflex:** r_star = −Σ m_i r_i / M_star, added analytically. This is what makes
the "solar system barycenter wobbles ~1 solar radius" validation target pass — the wobble is
bookkeeping, not runtime N-body.

**The solar system:** SPICE (DE441) is authoritative within its validity span (roughly
±14,000 years around J2000). Outside the span, the SAME Kepler + secular representation is
used, with elements and secular rates fitted to the SPICE arc and blended C¹ over a ±100 yr
window at each boundary. Beyond the span, positions are a self-consistent fiction and
documented as such: the real solar system is chaotic (Lyapunov time ~5 Myr), so no method —
including a full N-body integration — yields ground truth at T = +1 Myr. What the player gets
is deterministic, physically plausible, and exactly reproducible.

**Asteroids, comets, belt members:** Kepler + the same planetary secular precession applied
to their elements. Belt density fields precess coherently.

**Validation:** short-arc offline N-body cross-checks (REBOUND) of the secular representation
against SPICE inside the span, and against L-L theory for generated systems over 10⁵-10⁶ yr
arcs. N-body is a validation tool ONLY — it never runs in the simulator.

---

## 5. PLANET DETAIL GENERATION (Layer 3)

Triggered when the player inspects a specific planet.

### 5.1 Interior Structure

**Rocky planets (< 1.7 R_earth):**
- Iron core: mass fraction 0.25-0.40 (drawn from seed, mean 0.33 like Earth)
- Silicate mantle: bulk of remaining mass
- Volatile layer: water/ice if equilibrium temperature allows and formation location was beyond ice line
- Core size: R_core = R_planet × (M_core / M_planet)^(1/3) × (ρ_mantle / ρ_core)^(1/3)

**Sub-Neptunes (1.7-3.5 R_earth):**
- Rocky/icy core: 50-95% of total mass
- H/He envelope: 0.1-10% by mass (dramatically affects radius)
- Envelope mass fraction from mass-radius relation inversion (Lopez & Fortney 2014)

**Gas giants (> 3.5 R_earth):**
- H/He dominated atmosphere
- Heavy element core: 5-50 M_earth (correlates with planet mass and metallicity)
- Internal structure: molecular H₂ → metallic H → core (for Jupiter-class)

### 5.2 Atmosphere

**Atmospheric presence check:**
- No atmosphere if: T_eq > 2000K (stripped), mass < 0.1 M_earth (can't retain), no magnetic field AND around active M-dwarf (stripped by stellar wind)
- Thin atmosphere if: mass 0.1-0.3 M_earth or heavily irradiated
- Substantial atmosphere otherwise

**Composition (for rocky planets):**

Primary atmosphere (accreted H/He) is lost for planets < 10 M_earth in the HZ. Secondary atmosphere from outgassing:
- Default: N₂ + CO₂ (like early Earth/Venus)
- CO₂ fraction: from the carbon-silicate cycle equilibrium. Hot planets: CO₂ dominated (Venus). Temperate with plate tectonics: moderate CO₂ (~300-1000 ppm if life absent)
- If life exists and has reached the photosynthetic stage: O₂ accumulates. Post-GOE equivalent: 1-21% O₂
- Water vapor: from surface water inventory and temperature

**Surface pressure:** 0.01-100 bar (drawn from distribution). Correlation with mass: more massive planets tend to have thicker atmospheres.

**Temperature profile:**
- Surface temperature from T_eq + greenhouse warming: T_surface = T_eq × (1 + τ_IR)^0.25 where τ_IR is the infrared optical depth (from CO₂ + H₂O content)
- Lapse rate: dT/dz = -g / c_p (adiabatic) or modified by condensation
- Scale height: H = k_B T / (μ g), where μ is mean molecular weight

**Gas giant atmospheres:** H₂/He with trace CH₄, NH₃, H₂O. Cloud decks at condensation levels. Hot Jupiters: TiO, VO, alkali metals in upper atmosphere.

### 5.3 Magnetic Field

Dipole moment from the Christensen+ 2009 scaling law:

B_dipole ∝ (ρ_core × F_convective × R_core)^(1/3)

where F_convective is the convective heat flux in the core.

Simplified: scale from Earth's field using mass, rotation, and core state:

| Planet type | Typical B_surface (Gauss) | Condition |
|---|---|---|
| Earth-like, fast rotator, liquid core | 0.25-1.0 | Needs liquid iron core (age < ~10 Gyr for Earth-mass) |
| Super-Earth, fast rotator | 0.5-5.0 | Stronger due to larger core |
| Tidally locked (M-dwarf HZ) | 0-0.1 | Weak dynamo from slow/no rotation |
| Gas giant | 4-14 (Jupiter-like) | Strong, from metallic H convection |
| Ice giant | 0.2-0.5 | Moderate, from ionic water convection |

No magnetic field if: planet is tidally locked with no differential rotation, or core is fully solidified (old, small planet).

### 5.4 Rotation

**Rotation period:**
- Draw from the observed distribution (for solar system analogs)
- Rocky planets: 10-100 hours (if not tidally locked)
- Gas giants: 8-16 hours (spin-up from accretion)
- Tidally locked if: tidal synchronization timescale < system age
  - t_sync ∝ a^6 / (M_star × R_planet^3 × k_2/Q) (Goldreich & Soter 1966)
  - Most HZ planets around M-dwarfs are tidally locked
  - Mercury-like 3:2 resonances possible for moderate eccentricity

**Obliquity:** Drawn from uniform distribution 0-90° (isotropic). Earth's 23.5° is not special. Tidally locked planets: obliquity → 0 (or Cassini state equilibrium).

### 5.5 Surface and Geology

**Tectonic state:** From mass, age, internal heat budget:
- Mantle convective vigor: Ra = ρgαΔT d³ / (κη). Higher Ra → more vigorous convection
- Small planets (< 0.5 M_earth) or old planets (> 8 Gyr for Earth-mass): stagnant lid. Low volcanism, cratered surface.
- ~Earth mass, moderate age: plate tectonics possible (~30% chance, drawn from seed). Active volcanism, recycled surface.
- Large super-Earths (> 3 M_earth): stagnant lid (too viscous) or episodic resurfacing.

**Surface features (abstract description):**
- Cratered: old surface, no resurfacing (Mercury, Moon analog)
- Volcanic: active resurfacing (Io, young Venus analog)
- Tectonic: plate boundaries, mountain chains, ocean trenches (Earth analog)
- Eroded: atmosphere + liquid → weathering, sedimentary features (Earth, possibly early Mars)
- Oceanic: substantial water coverage (> 50%). Draw ocean fraction from seed.

**Water inventory:** From formation location relative to ice line, volatile delivery during late heavy bombardment analog, and atmospheric retention. Ocean depth: ~1-100 km for water worlds.

### 5.6 Weather (Static Climate)

Climate map as a function of latitude, from:
- Stellar flux: S(latitude) = S_0 × cos(latitude) (modified by obliquity and eccentricity)
- Greenhouse effect: ΔT_greenhouse from atmospheric composition (Section 5.2)
- Albedo: ice-albedo feedback at high latitudes
- Rotation rate: fast → more uniform temperature, strong Coriolis → jet streams. Slow/locked → large day-night contrast.

**Output per latitude band:**
- Surface temperature T(lat)
- Cloud coverage fraction
- Wind speed estimate (from pressure gradient and Coriolis)
- Precipitation type (liquid water, CO₂ snow, methane rain depending on T and composition)

For tidally locked planets: subsolar point temperature, antistellar point temperature, terminator ring climate. Atmospheric heat redistribution depends on surface pressure.

---

## 6. BIOSPHERE GENERATION

### 6.1 Habitability Assessment

**Habitable zone check** using Kopparapu+ 2013 (full formula with T_eff-dependent coefficients, see plan document for exact coefficients):

Conservative HZ: between runaway greenhouse (inner) and maximum greenhouse (outer).
Optimistic HZ: between recent Venus (inner) and early Mars (outer).

**Planet suitability:**

| Criterion | Requirement | Weight |
|---|---|---|
| In conservative HZ | Yes | Required |
| Mass | 0.3-10 M_earth | Required |
| Surface pressure | > 100 mbar | Required |
| Liquid water possible | T_surface 273-373 K at surface pressure | Required |
| Magnetic field present | From Section 5.3 | ×0.7 if absent |
| Plate tectonics active | From Section 5.5 | ×0.8 if absent |
| Not tidally locked | From Section 5.4 | ×0.5 if locked (habitable but reduced) |
| Host star is M-dwarf | Flare activity reduces habitability | ×0.6 |

Habitability score H = product of all weights, ∈ [0, 1] for qualifying planets, 0 for non-qualifying.

**Moon habitability:** Tidally heated moons with subsurface oceans (Section 4.6) get a separate habitability pathway:
- H_moon = 0.3 × (Q_tidal sufficient for liquid water) × (ice shell thickness allows chemical exchange)
- Lower H than surface habitable planets but non-zero. Microbial life possible, complex multicellular unlikely.

### 6.2 Life Stage Transitions (Sequential Poisson Model)

Each transition is a Poisson process with rate λ × H per Gyr. The seed generates a uniform random number U per transition. The transition occurred if U < P(transition by available time).

P(transition) = 1 - exp(-λ × H × max(0, t_available))

where t_available = system_age - t_onset (minimum time before transition can begin).

| Stage | Transition | λ (Gyr^-1) | t_onset (Gyr after formation) | Source |
|---|---|---|---|---|
| 0 → 1 | Prebiotic → microbial (abiogenesis) | 1.5 | 0.5 | Kipping 2020, Lineweaver 2002 |
| 1 → 2 | Microbial → photosynthetic | 0.5 | 1.0 after abiogenesis | Earth: ~1.3 Gyr gap |
| 2 → 3 | Photosynthetic → oxygenated (GOE) | Threshold | 0.5 Gyr after photosynthesis | Environmental, not evolutionary |
| 1 → 4 | Prokaryotic → eukaryotic | 0.3 | 1.5 after abiogenesis | Carter 2008 (hard step) |
| 4 → 5 | Eukaryotic → multicellular | 0.5 | 1.0 after eukaryotes | Evolved ~25× on Earth |
| 5 → 6 | Multicellular → complex (Cambrian equivalent) | 0.7 | Requires O₂ > 10% PAL | Oxygen-gated |

Stages are sequential: a planet must pass through each in order. Later stages require earlier ones to have occurred.

**No intelligent life.** Stage 6 (complex multicellular) is the maximum. This is a scope decision, not a probability judgment.

### 6.3 Biosphere Properties

For each inhabited planet, generate:

**Biomass:** Scale from Earth's ~2 × 10^15 kg:
- Biomass ∝ planet_surface_area × surface_coverage × H × (available_energy / Earth_energy)
- Available energy: from stellar flux in the photosynthetically active range (400-700 nm for Earth, shifted for other star types)

**Surface coverage:** Fraction of land/ocean with visible biomass:
- Microbial stage: ~5-20% (ocean surfaces, hydrothermal vents, soil)
- Photosynthetic stage: ~20-50% (ocean surface blooms, coastal mats)
- Multicellular stage: ~50-90% (land plants, ocean ecosystems)

**Dominant biome types (from seed):**

| Biome | Conditions | Characteristics |
|---|---|---|
| Deep-ocean chemosynthetic | Any stage, subsurface ocean or deep sea | Dark, hot vent communities. No photosynthesis. |
| Ocean surface photosynthetic | Photosynthetic+, ocean > 10% | Cyanobacteria/algae analogs. Surface blooms. |
| Coastal/intertidal | Multicellular+, ocean + land | Diverse, high productivity. |
| Terrestrial photosynthetic | Multicellular+, land areas | Plant analogs. Cover depends on climate. |
| Arid/extremophile | Any stage, marginal conditions | Sparse, stress-tolerant organisms. |

**Photosynthetic pigment:** Determined by host star spectrum:
- G-star (Sun-like): peak photosynthesis at ~680 nm → green/dark green pigment (chlorophyll analog)
- M-dwarf: peak flux in near-IR → pigment absorbs broadly, appears black/very dark red
- K-star: intermediate → dark green to brown
- F-star: more blue flux → possibly blue-green pigment

**Maximum organism size:** Scales with time since multicellularity onset:
- At onset: millimeters
- After 100 Myr: centimeters
- After 500 Myr: meters (Cambrian-level diversity)
- After 1 Gyr: tens of meters (large animals/trees)

### 6.4 Biosignature Outputs

| Biosignature | Stage required | Signal strength |
|---|---|---|
| Biogenic CH₄ | Microbial (1) | ~1000 ppm (pre-oxygenation), ~1.8 ppm (post-oxygenation) |
| Atmospheric O₂ | Oxygenated (3) | 0.1-21% PAL |
| O₃ (ozone layer) | Oxygenated (3) | UV absorption, protects surface |
| O₂ + CH₄ disequilibrium | Oxygenated + microbial | Simultaneous presence requires continuous biological source |
| Vegetation spectral edge | Multicellular (5) | Sharp reflectance change at pigment-dependent wavelength |
| N₂O | Microbial (1) denitrification | ~0.3 ppm (Earth-like) |

### 6.5 Expected Populations

With ~130 billion stars (physics.md Section 3.2.9), ~50% with planets, ~20 billion rocky HZ candidates:

| Life stage | Expected count in MW | Notes |
|---|---|---|
| Microbial | 2-10 billion planets | P(abiogenesis) ~ 10-50% for suitable planets |
| Photosynthetic | 0.7-3 billion | ~50% of microbial worlds develop photosynthesis given time |
| Oxygenated | 0.3-2 billion | Requires ~2 Gyr after photosynthesis |
| Eukaryotic | 70-700 million | Hard step, λ = 0.3 Gyr^-1 |
| Complex multicellular | 20-200 million | Requires O₂ + eukaryotes + time |
| Subsurface (moons) | 0.7-7 billion | Lower complexity, mostly microbial |

---

## 7. MINOR BODY DETAIL (Layer 3)

### 7.1 Individual Asteroids

When queried (player selects region within a belt):

- **Designation:** `{star_id}-Belt{letter}-{size_rank}` from seed
- **Size:** From the belt's size-frequency distribution (inverse CDF)
- **Mass:** From size and density (composition-dependent)
- **Shape:** Triaxial ellipsoid, axis ratios from observed distribution
- **Composition:** C-type (carbonaceous, ρ ~ 1.5-2.5), S-type (silicate, ρ ~ 2.5-3.5), M-type (metallic, ρ ~ 4-5). Fraction depends on belt location.
- **Rotation:** Period 2-20 hours. Axis orientation from seed. Spin barrier at ~2.2 hours for rubble piles.
- **Orbit:** Within belt boundaries. Keplerian elements from seed. Proper elements perturbed by resonances with planets.
- **Surface features:** Abstract (cratered, rubble pile, monolithic). Binary asteroids: ~15% of NEA analogs.
- **Visitable:** The ship can fly to any asteroid > 1 km and orbit it. Gravity computed from mass. Surface gravity: ~0.001-0.01 m/s² for km-class.

### 7.2 Comets

Generated from the Oort cloud / Kuiper belt analog:
- **Long-period comets:** Highly eccentric orbits (e → 1) from the Oort cloud. New comets scattered by passing star perturbation or galactic tides. ~1-10 active comets per system at any time.
- **Short-period comets:** From the Kuiper belt, scattered inward by giant planets. Jupiter-family comets with periods 5-20 years.

Comet properties:
- Nucleus: 1-50 km, icy (H₂O, CO₂, CO, organics), low density (~0.5 g/cm³)
- Activity: sublimation rate from distance to star. Active within ~3-5 AU for water ice. CO/CO₂ sublimation extends activity further.
- Coma: gas and dust cloud, radius ~10^4-10^5 km when active
- Tail: ion tail (anti-sunward) and dust tail (curved). Not simulated in detail but noted in properties.

### 7.3 Rings (Detail)

On close approach to a ringed planet:
- Particle density profile n(r, z): exponential in z (scale height ~10-100 m for dense rings), power-law or structured in r
- Gap structure from moon resonances (computed from moon orbital periods)
- Particle composition and albedo
- Ship interaction: micrometeorite flux from particle density and relative velocity

---

## 8. FREE-FLOATING OBJECTS

### 8.1 Brown Dwarfs

Objects with mass 13-80 M_jupiter (deuterium fusion threshold to hydrogen fusion threshold). Generated from the sub-stellar tail of the IMF in the same cells as stars.

Estimated population: ~25-100 billion in the galaxy.

Properties:
- Mass, age from seed
- T_eff: cools with age. Birth: ~2500 K. At 10 Gyr: ~300-500 K (Y-dwarf)
- Spectral type: L (1300-2400 K), T (600-1300 K), Y (< 600 K)
- Radius: ~1 R_jupiter (roughly constant with mass due to degeneracy)
- Luminosity: ~10^(-4) to 10^(-6) L_sun (mostly thermal IR)
- Can host tight planetary systems (0.01-1 AU)

### 8.2 Rogue Planets

Sub-stellar objects (< 13 M_jupiter) not bound to any star. Ejected during planetary system formation.

Estimated population: ~20 per star → ~2-3 TRILLION in the galaxy (Sumi+ 2023, MOA 9-year
microlensing — the first measurement sensitive below Earth mass; supersedes the older
"1-2 per star" estimates by an order of magnitude). The mass function rises steeply toward
sub-Earth masses: Earth-mass and smaller rogues dominate the count; Jupiter-mass rogues are
rare (~1 per ~20 stars).

Properties:
- Mass: drawn from the measured (steep) mass function, ~0.01-13 M_earth dominant range
- No stellar heating. Surface temperature from internal heat only: T_surface ~ 30-50 K for Earth-mass after 5 Gyr
- Retained H₂ atmosphere for massive ones (> 3 M_earth): could provide modest greenhouse warming
- Subsurface liquid water theoretically possible from radiogenic heating for massive rogue planets with thick H₂ atmospheres (Stevenson 1999). Extremely speculative.

### 8.3 Visibility

Both brown dwarfs and rogue planets are hidden by default in star queries. Shown with an explicit "substellar objects" filter. When the ship approaches within gravitational threshold, they become active point sources like any star.

---

## 9. ENVIRONMENTAL FIELDS (Per-Tick)

These fields are evaluated every tick and feed into the ship's physics modules.

### 9.1 Interstellar Medium

The ISM density, temperature, and composition at the ship's position, from the galactic gas disk model:

**Ambient density:**
- n_HI(R, z) from the HI gas disk (R_d = 7.0 kpc, h_z = 0.085 kpc + flare)
- n_H2(R, z) from the H2 gas disk (R_d = 1.5 kpc, h_z = 0.045 kpc)
- Total n_ISM = n_HI + 2×n_H2 (factor 2 for molecular)
- Solar neighborhood: n_ISM ~ 0.5-1.0 cm^(-3)
- Include helium: multiply by 1.1 (10% He by number)

**ISM phases (from position-dependent model):**

| Phase | n (cm^-3) | T (K) | Filling fraction | Where |
|---|---|---|---|---|
| Cold neutral (CNM) | 20-100 | 50-100 | ~2% | Disk midplane, spiral arms |
| Warm neutral (WNM) | 0.2-0.5 | 6,000-10,000 | ~40% | Disk, widespread |
| Warm ionized (WIM) | 0.1-0.3 | 8,000 | ~25% | Disk, near OB stars |
| Hot ionized (HIM) | 0.003-0.01 | 10^6 | ~30% | Superbubbles, halo |
| Molecular | 100-10^6 | 10-50 | ~1% | Molecular clouds, spiral arms |

The local phase is determined by proximity to hot stars (→ WIM), supernova remnants (→ HIM), and molecular clouds (→ molecular). Far from any of these: default to WNM.

**ISM affects the ship:**
- At v < 0.01c: negligible drag, minor radiation background
- At 0.01-0.3c: Bethe-Bloch energy deposition, heating, sputtering
- At > 0.3c: hadronic cascades, pion production (see physics.md Section 4.1)
- Radiation dose from cosmic rays: modulated by stellar wind within a system

### 9.2 Nebulae (Localized ISM Structures)

Nebulae override the ambient ISM within their boundaries.

**Generation:**

| Nebula type | Generated from | Radius | n (cm^-3) | T (K) |
|---|---|---|---|---|
| HII region | Every O/B star on MS | R_Strömgren = (3 N_Lyc / (4π α n²))^(1/3), typically 1-100 pc | 10-10^4 | 8,000-12,000 |
| Planetary nebula | Stars in PN phase (0.5-8 M_sun, age within ~20,000 yr of PN formation) | 0.1-1 pc (expanding) | 10^2-10^4 | 10,000-20,000 |
| Supernova remnant | Stars that went SN (age within ~100,000 yr of t_SN) | R_SNR = 0.3 × (E_SN / n_0)^(1/5) × t^(2/5) pc (Sedov-Taylor) | 0.1-100 (shell) | 10^6-10^8 (shell), cooler interior |
| Molecular cloud | Dense regions in spiral arms (from H2 density model) | 10-100 pc | 10^2-10^5 | 10-50 |

**Known nebulae:** Pre-loaded catalog of ~3000 from SIMBAD/VizieR (positions, sizes, types).

**Procedural nebulae:** Generated automatically from stellar properties. Every O/B star creates an HII region. Every recent SN creates an SNR. Molecular clouds are placed along spiral arms in the densest H2 regions.

**Ship interaction:** When inside a nebula, the ship experiences the nebula's density and temperature instead of ambient ISM. HII regions and SNR can produce significant radiation. Molecular clouds increase drag and absorb starlight (visual extinction).

### 9.3 Stellar Wind (System-Level)

Within a star system's astropause (see Section 4.9):
- Ship is immersed in the stellar wind, not the ISM
- Wind density, velocity, temperature, and magnetic field from Section 4.9
- Solar particle events (flares) propagate through the wind as density/velocity enhancements
- Modulates cosmic ray flux (heliosphere analog shields against GCR)

Transition from wind to ISM at the astropause is a discontinuity in density and temperature.

### 9.4 Radiation Environment

Position-dependent radiation field:

**Galactic cosmic rays (GCR):**
- Base flux: ~4 particles/cm²/s (at solar minimum equivalent)
- Energy spectrum: power-law, peaking at ~300 MeV/nucleon
- Modulated by stellar wind: reduced by ~50% at solar maximum equivalent (inside astropause)
- Higher in the halo (less shielding), lower in dense molecular clouds

**Stellar radiation:**
- Luminosity from host star(s), inverse square law
- UV flux: critical for atmospheric photochemistry and biological damage
- X-ray/EUV from active stars: higher for young/active stars

**Trapped radiation belts:**
- Near planets with magnetic fields: Van Allen belt analogs
- Proton and electron flux from dipole trapping geometry
- Dose rate can be extreme inside belts (100-1000× ambient)
- Belt inner edge: ~1.2 R_planet. Outer edge: ~6-10 R_planet (scales with magnetic field strength)

### 9.5 Per-Tick Evaluation Cost

| Field | Evaluation | Cost |
|---|---|---|
| Ambient ISM density | Gas disk model lookup | ~0.1 μs |
| ISM phase | Proximity check to active nebulae | ~0.1 μs |
| Stellar wind (if in system) | Analytic formula at distance r | ~0.05 μs |
| Nebula override | Check active nebula cache | ~0.1 μs |
| GCR flux | Base rate × modulation factors | ~0.05 μs |
| Stellar radiation | L / (4π r²) | ~0.02 μs |
| **Total** | | **~0.5 μs** |

Added to the per-tick budget. Total per-tick: ~30-115 μs → ~31-116 μs. Negligible impact.

---

## 10. QUERY SYSTEM

### 10.1 Query Types

**Region query:** "Show all objects within X parsecs of position P" (or of the ship).
**Filter query:** "Show all G-type stars within X parsecs" or "planets with oxygen atmospheres."
**System query:** "Show contents of this star system."
**Body query:** "Show detail for this body."

### 10.2 Performance Target

10,000 bodies listed in < 1 second.

### 10.3 Query Pipeline

```
1. Player specifies: center, radius, filters                    [input]
2. Estimate object count from density model                      [~0.1 μs]
3. If count > 100,000: shrink radius or require tighter filter   [~0.1 μs]
4. Determine overlapping spatial cells                           [~1 μs]
5. For each cell: generate Layer 1 data for all stars            [~55 ns/star]
6. Check Gaia index: substitute real stars, suppress duplicates  [~0.01 ms]
7. Apply Layer 1 filters (spectral type, mass, age, remnant)     [~5 ns/star]
8. If Layer 2+ filters: generate systems for L1 survivors        [~1 ms/system]
9. Sort by distance                                              [~0.1 ms for 10K]
10. Return top N with minimal data                               [~0.1 ms]
```

For 10,000 stars with Layer 1 filters only: steps 5-7 take ~0.6 ms. Total < 2 ms. Well within target.

For Layer 2 filters (e.g., "has habitable zone planet"): pre-filter at Layer 1 (only check stars with T_eff 3000-7000K, age > 1 Gyr) to reduce Layer 2 generation count. If 1000 stars pass L1 filter: 1000 × 1 ms = 1 second. At the limit. For expensive filters, the query system restricts the search radius to keep total computation under 1 second.

### 10.4 Filter Hierarchy

| Filter | Layer | Cost per candidate | Pre-filter strategy |
|---|---|---|---|
| Spectral type, mass, luminosity | 1 | ~55 ns | None needed |
| Binary system, remnant type | 1 | ~55 ns | None needed |
| Distance, galactic position | 1 | ~55 ns | None needed |
| Has planets / planet count | 2 | ~1 ms | L1: any main-sequence star |
| Has HZ planet | 2 | ~1 ms | L1: T_eff 3000-7000K, age > 1 Gyr |
| Has rings, has moons | 2 | ~1 ms | L1: any star with giant planet candidates |
| Has life | 3 | ~1 ms (L2) + ~0.1 ms (L3) | L1: as above + L2: has HZ rocky planet |
| Atmosphere composition | 3 | ~1.1 ms | L1 + L2 pre-filter |
| Surface conditions | 3 | ~1.1 ms | L1 + L2 pre-filter |

### 10.5 Query Radius Scaling

Filters that eliminate most candidates allow larger search radii for the same computation budget:

Corrected star-count arithmetic (at 0.1 stars/pc³, a 50 pc sphere holds ~52,000 stars — the
previous ladder's "~500 at 50 pc" was off ~100×). Layer-1-only filters are cheap at any
radius (52,000 × 55 ns ≈ 3 ms); the binding cost is Layer-2 generation (~1 ms/candidate):

- No filter: ~50 pc (~52,000 stars, ~3 ms L1; result list display-capped, sorted by distance)
- "G-type stars": ~150 pc (~4% of stars; L1-only, still milliseconds)
- "Has HZ planet": ~15-25 pc per second of compute budget — the L1 pre-filter (T_eff, age)
  passes most stars, so ~1,000 L2 generations ≈ 1 s sets the radius
- "Has life": same L2 bound (L3 adds ~0.1 ms only for L2 survivors); wide-radius rare-object
  searches are NOT feasible synchronously (2,000 pc = 3×10⁹ stars) — they run as PROGRESSIVE
  background queries: expanding-shell search in the display path, streaming results as found,
  fully interruptible, zero physics impact

Exact radius computed dynamically: estimate N_L2 = n_density × V_sphere × L1_pass_rate. Set radius so N_L2 × 1 ms < budget; beyond that, progressive mode.

---

## 11. THRESHOLD ACTIVATION AND INTEGRATION

### 11.1 The Uniform Rule

One formula everywhere: a body becomes a gravitational point source on the ship when:

a_body = GM_body / r² > ε × |a_total_ship|

where ε = 10^-6 and |a_total_ship| is the ship's current total gravitational acceleration.

This same formula activates Earth in LEO, Jupiter during interplanetary flight, nearby stars in interstellar space, and Sgr A* across the inner galaxy.

### 11.2 Activation Triggers Generation

When a body crosses the threshold, the engine needs its mass and position. This triggers the appropriate generation layer:

| Event | What happens |
|---|---|
| Star crosses threshold from interstellar distance | Already generated at Layer 1 (from periodic refresh). Mass and position known. No new generation. |
| Ship approaches a star system (star becomes dominant non-mean-field source) | Trigger Layer 2: generate full system architecture. Planets, moons, belts, rings. |
| Planet crosses threshold within system | Already generated at Layer 2. Orbital elements known. |
| Moon crosses threshold near planet | Already generated at Layer 2. |
| Large asteroid on very close approach | Generate individual asteroid from belt seed + index. |

The ship never "discovers" a body by threshold -- it knows about everything (omniscient). But the engine only computes gravitational effects when the threshold is crossed, which drives when detailed generation is needed.

### 11.3 Tidal Forces (Threshold-Activated)

Tidal tensor (gradient of gravity) is computed when the ship is within ~10× a body's Roche radius:

r_tidal_activation = 10 × R_body × (2 × M_body / M_ship)^(1/3)

For practical purposes: tidal forces are computed from the nearest 0-3 massive bodies (the ones whose tidal influence is non-negligible). Applied to the ship's structural model and to any nearby moons/rings for physics accuracy.

### 11.4 System Entry Detection

"Entering a star system" is defined as: the star's individual gravitational acceleration on the ship exceeds the mean field gradient AND is the dominant individual source. This naturally defines the transition from interstellar to system-level physics.

At system entry:
1. Layer 2 generation triggers (planets, moons, belts, rings, wind)
2. The astropause position is computed
3. If the ship crosses the astropause: ISM transitions to stellar wind
4. Planets begin appearing as individual point sources as the ship penetrates deeper

### 11.5 Mean Field Integration

The mean field provides:
1. **Gravitational potential** at the ship's position (for gravity computation)
2. **Stellar density** at any position (for determining how many stars to generate per cell)
3. **Gas density** at the ship's position (for ISM evaluation)
4. **Cluster positions** at current time T (from precomputed orbits)

Potential and density come from the same analytic potential-density pairs (McMillan
axisymmetric + bar + Cox-Gomez spiral + warp; physics.md Section 3.1). The density ρ that
sources Φ via Poisson's equation is the same density that determines star counts —
self-consistent by construction, with no SCF/analytic mismatch. The LMC halo-response
perturbation term contributes potential only (its density correction to local star counts
is negligible inside 30 kpc).

### 11.6 Gaia Integration

1. Gaia stars are stored in a k-d tree indexed by guiding center position (propagated from catalog epoch 2016.0 to J2000.0)
2. On each nearby-star refresh and on each query: the Gaia index is checked alongside the procedural generator
3. Where Gaia stars exist: the procedural star count for that cell is reduced by the number of Gaia entries
4. Gaia stars have priority: their real properties (parallax, proper motion, T_eff, [Fe/H]) override procedural generation
5. For Gaia stars with known exoplanets: the real planetary system is loaded and AUGMENTED with procedural planets in stable gaps (Section 4.4)
6. Within ~500 pc: most bright stars are Gaia entries. Beyond: mostly procedural with scattered Gaia landmarks.

---

## 12. OBJECT DESIGNATION SYSTEM

All objects get systematic IDs from their generation coordinates:

| Object type | Designation format | Example |
|---|---|---|
| Star | `GCS-R{R_kpc}-P{phi_deg}-Z{z_pc}-{index}` | `GCS-R8.247-P000.00-Z+021-0042` |
| Gaia star | Gaia DR3 source_id | `Gaia DR3 4111834567890` |
| Named star | Common name | `Proxima Centauri` |
| Planet | `{star_id} {letter}` (by orbital distance) | `GCS-...-0042 b` |
| Moon | `{planet_id}.{number}` | `GCS-...-0042 b.1` |
| Asteroid | `{star_id}-Belt{letter}-{rank}` | `GCS-...-0042-BeltA-00001` |
| Brown dwarf | `GCS-BD-{coords}-{index}` | `GCS-BD-R8.1-P359-Z-002-0001` |
| Rogue planet | `GCS-RP-{coords}-{index}` | `GCS-RP-R8.3-P001-Z+005-0003` |
| Nebula | Catalog name or `NEB-{coords}` | `Orion Nebula` or `NEB-R7.2-P210` |

Designations are deterministic from seed. Re-querying produces the same designation. Players can bookmark any designation and return to it.

---

## 13. THE COMPLETE ENGINE LOOP

### 13.1 Per-Tick (~31-116 μs)

```
Gravity:
  1. Mean field gradient at ship position                        [1-3 μs]
  2. Point source gravity from all active sources                [1-10 μs]
  3. Tidal tensor from nearest massive bodies (if activated)     [0-5 μs]

Environment:
  4. ISM/stellar wind density and composition                    [0.5 μs]
  5. Radiation field (GCR + stellar)                             [0.1 μs]

Ship physics:
  6. Assemble total force (gravity + thrust + drag + ...)        [1 μs]
  7. Integrate trajectory (DOP853 step)                          [10 μs]
  8. Update thermal, resources, atmosphere, structural           [10 μs]
  9. Update proper time and coordinate time                      [0.1 μs]
```

### 13.2 Periodic Refresh (~0.1-1 ms, every few seconds sim-time)

```
  1. Query spatial index for candidate stars                     [0.1 ms]
  2. Evaluate positions at current T                             [~50 ns × N]
  3. Compute acceleration from each on ship                      [~10 ns × N]
  4. Apply threshold: update active point source list            [~1 μs]
  5. Check for system entry/exit                                 [~1 μs]
  6. If new system entered: trigger Layer 2 generation           [1-10 ms, once]
  7. Update active nebula cache                                  [~0.1 ms]
```

### 13.3 On Player Query (< 1 second)

```
  1. Parse query (center, radius, filters)
  2. Generate Layer 1 data for region
  3. Apply filters
  4. Generate higher layers if needed
  5. Return sorted results with minimal data
```

### 13.4 On Body Inspection (~0.1-1 ms)

```
  1. Generate Layer 3 for the body
  2. Return: geology, atmosphere, magnetic field, rotation,
     weather, biosphere assessment, full orbital elements
```

### 13.5 Data Flow

```
PRECOMPUTED DATA (18 MB + Gaia 500 MB)
  │
  ├── Hybrid analytic mean field (+ LMC-response coefficients) ──→ Gravity on ship (per-tick)
  │                  └──→ Same potential-density pairs ──→ Star counts per cell
  │
  ├── Cluster orbits ──→ Cluster positions at T
  │
  ├── Epicyclic tables ──→ Star position evaluation
  │
  └── Gaia catalog ──→ Real star overlay + known exoplanets

PROCEDURAL GENERATION (on-demand, deterministic from seed)
  │
  ├── Layer 1: Star catalog ──→ Active point sources ──→ Gravity
  │
  ├── Layer 2: Systems ──→ Planets, moons, belts, rings, wind
  │
  ├── Layer 3: Body detail ──→ Geology, atmosphere, biosphere
  │
  └── Layer 4: Surface ──→ Atmospheric profile ──→ Ship drag/heating

ENVIRONMENTAL FIELDS (continuous per-tick)
  │
  ├── ISM density from gas disk model
  ├── Nebula overrides from procedural generation
  ├── Stellar wind from current system
  └── Radiation from stellar + cosmic sources
```

---

## 14. CONSISTENCY AND VALIDATION

### 14.1 Determinism

Every path through the generation pipeline is a pure function of (seed, T). No runtime state, no mutable caches (caches are purely performance optimizations that could be recomputed from scratch), no randomness outside the seed.

Test: generate a system at T=0, advance to T=5,000,000, regenerate the same system. All seed-derived properties (mass, composition, biosphere stage) identical. Time-dependent properties (stellar position, planet orbital position, flare occurrence, SNR expansion) consistent with 5 million years of deterministic evolution.

### 14.2 Population Checks

- Total star count from density model integration: should yield ~1.3×10¹¹ (the mass-model-derived value, physics.md Section 3.2.9 — NOT an independent target)
- IMF integration: should yield ~75% M-dwarfs, ~10% K, ~3% G, ~1% F, ~0.3% A/B/O
- Binary fraction: ~46% for solar-type, decreasing toward galactic center
- Planet occurrence: ~50% of FGK stars with at least one planet
- Life occurrence: ~2-10 billion microbial worlds, ~20-200 million complex biospheres (scaled to the derived ~1.3×10¹¹ star count)

### 14.3 Scope Boundaries

Within scope (100% realistic, no approximations):
- Stellar orbits in the mean field
- Planetary system architecture and dynamics
- Planet properties from physical models
- Biosphere from Poisson life model with literature-calibrated rates
- Environmental fields from galactic gas/radiation models
- Gravitational physics on the ship

Out of scope (acknowledged limitations):
- Star-star gravitational interactions (encounters)
- Galactic center orbital element drift (~0.3% over 10 Myr)
- Rigid spiral-pattern rotation in star counts: arms are material, frozen at epoch in guiding-center space; drift vs. the potential's rigid pattern is ≤3° at R₀, ≤14° inner disk over the window (Section 3.2)
- Resonant/chaotic planetary dynamics: secular theory is first-order Laplace-Lagrange; the solar system beyond the SPICE span is a self-consistent fiction (Section 4.10 — the real system is chaotic with ~5 Myr Lyapunov time, so no method yields ground truth)
- Intelligent life
- Stellar interior dynamics (convection, magnetic cycles -- simplified to observable correlations)
- Detailed atmospheric chemistry (simplified to equilibrium + biosignature overlays)
- Plate tectonics dynamics (simplified to active/inactive state)
