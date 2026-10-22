# The face inlet (a) and the plate's wall temperature

Registered 22 October 2026, about 02:00 CDT, before any code. It implements the interface spec in `docs/SESSION_HANDOFF.md` (*For inquiry2: the face inlet (a) and the face's wall*, 21 October about 02:40) for `tier1_ensemble_plan.md`'s decisions 6 and 8. The spec's properties and gates are kept; what follows fixes inquiry2's choices and states the gates as they will be judged.

## 1. The open-area supply

**The switch.** `Definition::openAreaSupplies{false}`. False runs today's ring-centre code unchanged. True applies this section to every supply of a Chamber.

**Geometry, per level.** Face ring j spans node radii [ρ_j, ρ_{j+1}] = R_0 [f_j, f_{j+1}] (R_0 the face radius, f the mesh's ring fractions).
- Supply s's share of ring j: a_sj = π (min(r_o, ρ_{j+1})² − max(r_i, ρ_j)²), or 0 if not positive. A ring the supply covers whole takes a_sj = `axialArea(0, j)` exactly, so an aligned layout keeps today's areas bit for bit.
- φ_sj = a_sj / `axialArea(0, j)`. The supply's area is `supplyArea_[s]` = Σ_j a_sj, the true area to rounding, so the mass through the face sums to ṁ_s opening(t) to rounding at any mesh.
- The wall share is w_j = 1 − Σ_s φ_sj, set to 0 when below 1e-12 (rounding of a covered ring).

**The ring's kind (inquiry2's choice, as the spec proposed).** Ring j is *plate* if w_j ≥ 1/2, else it takes the supply with the largest φ_sj (the lower index on a tie). The kind is `faceSupply_[j]`, as today, and it alone governs:
- the gradient stencil's face neighbour (a wall value at the face midpoint, or zero normal gradient);
- the wall distance's plate;
- the wall functions' prescribed first cell, and the resolved SST's wall ω cells.

**φ and w govern the fluxes.**
- Inviscid: ring j's face flux is Σ_s φ_sj `supplyFace(s, g_s, ring j's interior)` + w_j × the plate's pressure-only flux. Species: Σ_s φ_sj g_s Y_s; ρk and ρω: Σ_s φ_sj × each supply's turbulence flux. g_s = ṁ_s opening(t) / `supplyArea_[s]`. Choking and the characteristic stay per supply against ring j's interior.
- A zero-flow supply (ṁ_s = 0) is a closed port: `supplyFace` at g = 0 returns pressure only. It counts in Σ φ, so its area is not wall: no shear and no heat (the LOX slots).
- Viscous, the wall share:
  - a *plate* ring passes w_j times its whole wall flux (resolved or the law's τ_w and q_w; k and ω fluxes and the corrections' fluxes too);
  - a *supply* ring with w_j > 0 passes w_j times the wall's momentum and heat flux only. Under the wall functions it is a wall face solved from the cell's state at the face distance, its cell not prescribed. Resolved, it is the wall gradient with μ_t = 0 at the wall. Its k and ω are solved, under the supply share's zero-gradient face, and no k or ω flux to the partial wall is passed.
  - *Declared limit:* on a supply ring the partial wall's sink of k and ω (resolved) is not represented. It lies under half the face (inferred small; tier 1 runs the wall functions, where it does not arise).
- Rates: the plate share's pressure force goes to `wallAxial`, the supplies' share (closed ports included) to the inlet momentum. Every viscous end flux at i = 0 is the wall's, so it goes to `wallAxial` and `wallHeat`.

**Refusals under the switch.** "Holds no ring centre" is dropped. Refused: two supplies' true annuli overlapping; a supply's outer radius beyond R_0 by more than 1e-12 relative.

**AMR.** Each level builds its own a, φ, w and kinds from its own face rings in its constructor; `core/amr.cpp`'s copy of the base parent's supply is skipped under the switch. A coarse face area is the sum of its children's, so a coarse a_sj is too.

## 2. The plate's wall temperature

`Definition::plateTemperature`, a `std::optional<double>`. The plate's temperature T_p is its value, or `wallTemperature` when empty; 0 is adiabatic. Plate faces read T_p where they now read `wallTemperature`:
- the wall value of T in the resolved wall gradient, and the adiabatic wall's zero heat;
- the resolved SST's wall ν, ρ_w and μ_w for the plate's cells (the corner cell (0, nr − 1) stays with the side wall, as it is listed today);
- the wall law's T_w.

**The adiabatic law** (member 19). With T_w = 0 on a wall face under the wall functions: β = 0, q_w = 0, and T_w = T_1 + r u_1² / (2 c_p) inside the solve (WALL_FUNCTIONS.md, eq. 10.8 with β = 0; r = Pr_w^(1/3) at T_w, c_p the layer mean). The refusal at `core/flow.cpp:107` then admits an adiabatic wall.

## 3. Gates

- **F1, switch off:** the fast ctest gives 7 pass and 4 fail with the sorted FAIL lines identical to `docs/evidence/amr_c/ctest_mass_error_2026-10-21.txt`. rl10dev's short run (`fr`, 76 × 8, to 30 µs) is `cmp`-identical to the parent commit's (mesh, history, the three fields, the log apart from wall times).
- **F2, switch on, aligned** (every φ 0 or 1): a chamber with supplies on ring nodes, laminar and SST with the wall functions, 200 steps. The general path is run (no branch for φ = 1). *Expected:* not `cmp`-identical, since the species flux is formed as area × (g Y) in place of (area × g) Y; the state's largest relative difference under 1e-10, and the face mass and momentum rates equal to 1e-12 relative.
- **F3, switch on, tier 1's true layout** (with the LOX slots as closed ports) at 149 × 16, at 1192 × 128 and on AMR (149 × 16 with 3 levels): Flow accepts it; each supply's face mass flow equals ṁ_s opening(t) to 1e-12; and on every ring it touches, its face ρu times its true area equals ṁ_s opening(t) to 1e-12, so its velocity and J are the element's own.
- **F4, AMR:** over each coarse face ring, the sum of its children's area × flux equals its own to 1e-12 relative, for mass, momentum and energy, at a uniform state (all levels' interiors equal). In a developed state the mass and energy flux still do (they do not depend on the interior); the momentum flux, which carries the interior pressure, is reported.
- **F5, the plate temperature:** (a) empty and set equal to `wallTemperature` give `cmp`-identical runs; (b) at 0 with the wall functions, the plate's q_w is 0 to rounding on every plate face, and its T_w lies at the recovery temperature (between T_1 and T_1 + u_1²/(2 c_p)).

## 4. Results

Measured on backhouse (Release, 4 threads at most) at 566b12c, against the registration commit d9e2684 built beside it for F1. Raw outputs in `face_inlet/`: `f245.txt` (`crucible_face_inlet_tests`, ctest `face_inlet`), `f3.txt` (`CRUCIBLE_T1_OPEN_AREA=1 crucible_rl10_tier1 149x16 1192x128`), `f1_rl10dev.txt` (the `cmp` and log diff, then the exit code of the diff of tier 1's default dry run against d9e2684's build), `ctest_2026-10-22_backhouse.txt` (the fast ctest and its FAIL lines against d9e2684's and the Mac record's).

- **F1: passes against the parent on the same machine. As written, against the Mac record, it fails, and the failure is the machine's.**
  - rl10dev's short run: mesh, history and the three fields are `cmp`-identical. The log differs only in its two wall-time lines (46.1 against 48.3 µs per cell update). Tier 1's dry run without the switch prints the same text as d9e2684's (diff exit 0).
  - The fast ctest (`face_inlet/ctest_2026-10-22_backhouse.txt`) gives 8 pass and 4 fail of 12. The 4 failures are the record's (amr, amr_subcycling, transport, spray), and the 12th test, `face_inlet`, passes. The sorted FAIL lines (16) are identical, values included, to d9e2684's built and run beside it.
  - Against `amr_c/ctest_mass_error_2026-10-21.txt` (the Mac) they differ in three places, as written: a second spray failure, KH-RT liquid mass through every child at 1.044e-12 against 1e-13; and amr's level 2 and 3 closures at 3.943e-14 and 8.024e-14, against 3.938e-14 and 8.132e-14. d9e2684 gives the same three differences on backhouse, so they belong to the machine and not to this change (measured). The cause, compiler and libm, is inferred.
  - The Mac's ctest was not rerun tonight (the Director: nothing heavier than a quick unit test on the Mac while it is short of memory).
- **F2, pass.** Not `cmp`-identical, as expected. The state's largest relative difference is 4.36e-14 laminar and 3.30e-14 under SST with the wall functions (limit 1e-10); the inlet mass flow differs by 0, the inlet momentum flux by 1.6e-16 and the wall's axial force by 7.5e-16 and 3.7e-16 (limit 1e-12). The open run's budgets are 3e-16 to 4e-16.
- **F3, pass,** for every grid named: per ring, the face mass flow equals Σ a_sj ṁ_s/A_s to 5.8e-16 (149 × 16), 2.5e-14 (1192 × 128) and 7.3e-15 at most over the AMR levels; per supply, the face mass flow and ρu A_true on each ring it touches to 1.0e-15, 1.0e-15 and 3.0e-15 (38, 129 and up to 76 supply-rings); the supply areas to 8.6e-16 and 1.0e-15 of the true ones. Under the switch Flow accepts the true slots (33 supplies), which the ring-centre code refuses at both grids (G6, unchanged). One step on each grid gives a mass budget of 1e-15; its energy budgets (1.9e-7 at 149 × 16, 2.9e-8 at 1192 × 128) are reported, not gated: the run is cold and unignited and one step says nothing about the walk.
- **F4, pass.** At the uniform state the children's area × flux matches each coarse ring's to 2.9e-16 (mass), 1.7e-16 (axial momentum) and 4.6e-16 (energy). After 100 base steps (t 2.5e-5 s) mass and energy still match to 2.9e-16 and 4.6e-16; the momentum flux differs by 3.1e-2 relative, as expected: it carries each level's own interior pressure (reported).
- **F5, pass.** (a) Empty and 400 K (= `wallTemperature`) give bit-identical state, species, ρk, ρω and budgets, with the wall functions and with resolved SST. (b) With `plateTemperature` 0 and the wall functions, the largest plate |q_w| is 1.5e-12 W/m² (ring-centre, 2 plate faces) and 2.3e-12 W/m² (open-area, 3 plate faces), against the side wall's 1.8e4 W/m² (ratios 8e-17 and 1.2e-16, limit 1e-9). T_w lies between T_1 and T_1 + u_1²/(2 c_p) on every plate face (both bounds met to 0 K), at 0.892 of the dynamic temperature rise, which the law sets to Pr_w^(1/3) (Pr_w about 0.71, inferred, not printed).
