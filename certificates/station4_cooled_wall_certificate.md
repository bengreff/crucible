# CRUCIBLE Station 4 Certificate — the Cooled Wall

Goal-B station 4 (SOLV-1 §3.5; COUP-2 §3.5 subset; META-3 `wall-function-heat`, `robin-robin-cht` lineage): supersonic hot gas (M = 2, 800 K static — free-stream recovery ≈ 1374 K) flows through a straight cylindrical duct (R = 0.05 m, L = 0.4 m) whose 0.01 m steel-class liner (κ = 20 W/m·K) is regeneratively cooled outside (h_cool = 20000 W/m²·K film at 300 K). The gas-side flux comes from the **one local wall-function law** — Colburn-class, near-wall state only, **±20–30% declared band** — entering each SDC sweep's class-D implicit solve as Robin interface data (COUP-2 §3.5) with the gas debited exactly what the solid received, i.e. applied with opposite signs to both sides (interface conservation by construction). The liner conducts on the `Solid` region of the same world grid through the Goal-A-certified operator; regions (gas/solid/exterior) are config-time data through the widened FND-2 §3.6 ingest seam. Criteria are the named constants in `station4_cooled_wall.rs`, asserted by `solv1_station4_cooled_wall.rs`.

## The coupled duct at steady state (12492 steps)

| quantity | value |
|---|---|
| steadiness residual (ρ and T_solid, check window) | 6.16e-4 |
| ṁ | 7.7596 kg/s |
| T0 in → out | 1439.94 → 1422.83 K |
| gas enthalpy deficit ṁ·c_p·ΔT0 | 133295.0 W |
| wall-face exchange Σq·A | 133744.1 W |
| coolant extraction | 133744.0 W |
| ledger closure (wall vs gas, wall vs coolant) | 3.36e-3, 1.05e-6 |
| mid-duct film h | 2298.9 W/m²·K |
| mid-duct wall flux q | 1.063e6 W/m² (rocket-scale) |
| liner ΔT at mid-duct (inner → outer) | 762.4 → 399.6 K |
| **series-resistance oracle, worst dev past entrance band** | **9.44e-4** |
| T_aw(wall cell) / free-stream recovery at mid-duct | 0.940 |

**The oracle** (VAL-1 rung i): at steady state the coupled system must reproduce the cylindrical film + ln-annulus + coolant-film series-resistance solution, built from the simulated *gas* state and declared coolant data only — the simulated solid field never enters. Pointwise agreement to 9.44e-4 (gate 5e-3) past the 8-cell entrance band says the wall law, the Robin faces, the region-masked conduction, and the Robin-Robin exchange inside the class-D solve compose into exactly the textbook conjugate solution.

## Stair-interface conservation (the stepped cavity)

Hot gas at rest in a closed stepped cavity against a cold stair liner (face directions exercised: RPlus, ZPlus), insulated exterior: after 500 coupled steps the gas lost 2.306330e-1 J, the liner gained 2.306330e-1 J — relative mismatch **3.39e-12** (round-off; the same q·A is applied to both sides and the shared face areas are bitwise identical by the session-7 `face_radius` single-owner guarantee).

## Exactness & determinism

- **Uniform rest at coolant temperature is a bitwise fixed point** of the full coupled step (gas fields and solid field byte-identical after 200 steps; the wall machinery adds no drift — the exchange at equal temperatures is exactly zero, unit-tested in `wall_heat`).
- **Reruns are bit-identical** on every field (asserted on a real march).
- **Robin-face analytic anchor:** the annulus with Dirichlet inner / Robin outer reproduces the exact steady `T(r) = T1 + (T∞−T1)·ln(r/r1)/(ln(r2/r1)+κ/(h·r2))` to < 2e-3 relative (32 radial cells).

## Declared bands & deferrals

- The wall law carries its **±20–30% declared closure band** (SOLV-1 §3.5) — the certificate's oracle validates the *coupling*, not the closure; the closure's truth enters station 5 as a declared band in the p-box.
- **Near-wall sampling:** the wall-adjacent cell is itself cooled, so its recovery temperature reads 6.0% below free-stream at Δr = 0.0025 m. This cell-size dependence of the un-resolved-boundary-layer closure is inside the declared band and shrinks when FND-3 geometry + finer wall cells arrive.
- **Integrator (S2):** the ONE SDC-IMEX step — Δt is the gas CFL alone; the solid + exchange are inside the class-`D` implicit solve (fixed-cycle CG, Robin-Robin Picard/Aitken fixed sweeps, COUP-2 audit every step). The old explicit-splitting thermal stability guards are gone by construction. The liner ρc_p = 200 J/m³·K remains a **declared steady-state continuation device** (steady solution independent of ρc_p; the physical value ~3.6e6 — now legal under the implicit class — lands with the plan's S4 wave).
- **Coolant side** is a configured Robin film — COUP-7's cooling-jacket boundary object (channel correlation, coolant return state) supersedes it.
- **Bartz nozzle-envelope cross-check** (VAL-2, `bartz` oracle) rides with station 5, where the nozzle + liner assembly exists; the stair-interface machinery it needs is certified here by the stepped cavity.

