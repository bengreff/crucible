# Table B: a tabulated finite-rate progress manifold

Status (6 October 2026, about 02:34 CDT): design and criteria, stated before any code. Nothing here is built yet. A research probe of the design's two derived claims ran at 02:38 to 02:47 and amended the snap and criterion 3 (*Research probe*, below). Probe (4) (03:01) found that burnt gas mixed into fresh gas below about 70% burnt is on no curve the design builds. The design was amended at 03:10 with a counted CVODES fallback and gate 2(c) restated (*Probe (4)*, below).

Table B replaces the per-cell CVODES integration of `Chemistry::FiniteRate`, with and without the PaSR closure, by a table built offline with Cantera and CVODES on the full mechanism (TECHNICAL_PLAN, *Lightweight engine*). Unlike Table A, it changes the model: the composition is constrained to a one-parameter family of states for each element mix, energy and density. The reference is therefore the direct finite-rate integration, and the criteria test what that constraint costs.

The literature behind the choices is a Sonnet agent's search on 6 October. Its summaries were read; none of the papers was read in full tonight. Where a step rests on my own reasoning, it is marked *derived*, and where on a number not yet measured, *guessed*.

## What is tested

### The manifold

- **Trajectories.** Each node of (f_H, Z_N, e, ln ρ) is a constant-(e, ρ) homogeneous reactor integrated by CVODES on the full mechanism, from a start state to equilibrium. A homogeneous reactor at the cell's state is what PaSR's fine structure is, so the table adds no second turbulence-chemistry closure (TECHNICAL_PLAN). The precedent for PaSR with tabulated homogeneous-reactor trajectories is Golovitchev, Nordin, Jarnicki and Chomiak (SAE 2000-01-1891) and Kärrholm's thesis (Chalmers, 2008), from the search.
- **Two branches per node, which meet at equilibrium.**
  - The lower branch approaches equilibrium from below in c.
    - Its start is the reactants (the node's elements as H2, O2 and N2) where they can carry the node's e at 200 K or above.
    - Otherwise, the start is the equilibrium composition at the node's elements and ρ at the highest temperature whose frozen composition carries e at 200 K or above. That covers gas that has lost energy since it burned: the nozzle's expanding flow and the cooled wall layer.
  - The upper branch approaches equilibrium from above in c. Its start is the complete-combustion products at the node's e and ρ, so it covers compressed or reheated gas that must dissociate.
- **Why one variable can carry the late phase** (*derived*). The bimolecular chain reactions (H + O2 = OH + O, O + H2 = OH + H, OH + H2 = H2O + H) conserve the number of moles and are fast. Only the three-body recombinations change the number of moles, and they are slow. Once the bimolecular steps are near equilibrium, the composition at given elements, e and ρ has one degree of freedom left. This is the partial-equilibrium picture of H2/O2 recombination; my source for it is textbook memory, not re-read tonight. So nozzle and wall-layer states lie on the same one-parameter curve whatever path brought them there, and the start of the lower branch only matters until the trajectory reaches that curve.
- **Heat loss needs no extra axis.** e is an independent axis, not tied to an adiabatic flamelet as in a flamelet table. Gas cooled at the wall has a lower e and is found on that node's lower branch. The literature's non-adiabatic flamelet tables add an enthalpy-defect axis (Ma, Wu, Ihme and Hickey, AIAA J. 2018; Perakis, Haidn and Ihme 2021) because a flamelet has none. For H2/O2, Betti et al. (AIAA J. 54(5), 2016) put near-wall recombination at 7 to 14% of the wall heat flux against frozen flow, so it matters for the RL10's expander heat.

### The progress variable

- **c = Y_H2O.** The cell carries its own Y_H2O, so c needs no new transported variable.
- **Why it resolves the induction period** (*derived*). During induction, H2O forms through OH + H2 = H2O + H, one of the three chain steps. The net branching cycle is H + O2 + 3 H2 → 2 H2O + 3 H, so Y_H2O grows in proportion to the radical pool, by the same exponential. Then it rises through the heat release and through recombination (H + OH + M → H2O + M).
- **The c axis.** Nodes are geometric from 1e-9 of the node's equilibrium Y_H2O up to a tenth of it, then uniform to equilibrium on the lower branch. The upper branch is uniform. c = 0 is a node.
- **Monotonicity is checked, not assumed.** The literature warns that a progress variable can turn back (Ihme, Shunn and Zhang, J. Comput. Phys. 231, 2012; Prüfert et al., Flow Turbul. Combust. 2015). The builder checks every trajectory: c must move monotonically from its start to within 1e-4 of its equilibrium value. If Y_H2O fails anywhere, the next candidate is used, in this order:
  1. Y_H2O + Y_HO2 + Y_H2O2 (the low-temperature path binds O2 as HO2 and H2O2 before H2O);
  2. the O2 consumed.
  If none passes on every trajectory, the builder stops and the design comes back here. It is not tuned.
- **Preferential diffusion** makes H2 progress variables non-monotonic in flamelets (Najafi-Yazdi, Cuenot and Mongeau 2012; Sci. Rep. 2024). That is a diffusion effect. A homogeneous reactor has none, so the risk here is chemistry alone. The free flame (criterion 4) is where diffusion is tested.

### Stored values and lookup

- **Per point:** the H/O species' mass fractions (N2 and AR are fixed by their elements), and the relaxation rate λ = (dc/dt) / (c_eq − c). λ is positive on both branches and finite at equilibrium, where the reaction becomes a linear relaxation.
- **Lookup.**
  - Y is interpolated linearly in c at each corner of the other four axes. So every corner returns Y_H2O = c, and the corners' own elements are kept.
  - The other four axes are interpolated multilinearly. As in Table A, the element mass fractions are bilinear in (f_H, Z_N), so the query's elements are reproduced to rounding.
  - λ is interpolated linearly in c at each corner, then ln λ multilinearly across the other axes. Linear in c is exact for a rate proportional to c (pure branching) or affine in c (initiation plus branching), and the logarithm follows the Arrhenius dependence on e.
  - c_eq comes from the same corners.
  - e uses Table A's first-design normalisation, η between the complete-combustion products at 200 K and the same products frozen at 6,000 K, per (f_H, Z_N) corner. Table A moved to a T axis because equilibrium is sharp in T. Here e must be the axis, because each trajectory is at constant e.
- **The reaction step.**
  - It advances c at constant (e, ρ) by sub-steps of at most one c node. Each sub-step is the exact exponential update c_eq − (c_eq − c) exp(−κ_eff λ Δt), with λ at the sub-step's midpoint by predictor-corrector, and κ_eff = 1 without PaSR.
  - It then sets the cell's composition to the manifold's at the new c.
  - Stiff relaxation therefore lands on c_eq exactly at any Δt. The table's fast-chemistry limit is local equilibrium, consistent with Table A.
  - T comes from the core's exact e-to-T inversion, so energy is exact.
- **States off the manifold: snap and measure.** Setting the composition to the manifold's is the standard flamelet and FGM practice: the fast directions are assumed to relax at once. Mixing in a cell (PaSR, numerical diffusion) produces compositions that are not on the curve. The snap moves them onto it at the same c, e and ρ, and T changes with it.
  - The change |ΔT_snap| is the measured departure, reported as a distribution in criterion 5.
  - The remedy if criterion 5 fails is the ISAT precedent (Pope 1997; Lu, Lu and Pope 2009): direct CVODES above a stated |ΔT_snap|. It is not built until it is needed.
- **Outside the table** (e or ρ off the axes, AR above 1e-12, or a node whose trajectory did not reach equilibrium): direct CVODES. Counted.

### Size and build (*guessed*, to be set by criteria 2 to 5)

- **Axes.**
  - f_H: 33 nodes, clustered at stoichiometric and both ends as in Table A.
  - Z_N: 5 nodes (0, 0.1, 0.3, 0.6, 1). N2 matters only for C1's fill, and the RL10 uses the Z_N = 0 slice.
  - η: 48 nodes.
  - ln ρ: 21 nodes, 1e-3 to 1e2 kg/m³, 4 per decade.
  - c: 40 nodes.
- **Size.** About 6.7 million points, at 9 doubles each about 480 MB. Doubles are needed for element conservation to 1e-12. Limit: 1 GB, so the Mac (4.5 GB free tonight) can hold it.
- **Build.** About 330,000 CVODES integrations: 166,000 nodes with two branches each, for the 10-species `h2o2.yaml`. At 5 to 20 ms each that is 30 min to 2 h on one core, or a few minutes across backhouse's threads.
  - Cold nodes (reactants near 300 K) have astronomically long induction. They are integrated with no time limit and an absolute tolerance of 1e-30, because only the compositions along c are kept, not the time. A trajectory that does not reach equilibrium is reported, and its node falls back to direct CVODES.
- **Mechanism.** The first table is `h2o2.yaml`, for C1. The RL10's production mechanism (the Burke et al. 2012 or Li et al. 2004 class, TECHNICAL_PLAN *Chemistry*) gets its own table, built and checked the same way.

## Research probe: the two derived claims, tested (6 October 2026, 02:38 to 02:47 CDT)

Exploratory, not a criterion. `crucible_table_b_probe` (tests/table_b_probe.cpp) runs the engine's own CVODES reactor (rtol 1e-10) on h2o2.yaml at O/F 6 by mass and a 3 MPa chamber (adiabatic equilibrium 3,485 K, c = 0.8544). It builds each node's lower branch as stated above. Output: `table_b/probe_2026-10-06.txt`.

- **(1) Y_H2O is monotone in ignition.** From the reactants at 900, 1,000, 1,200 and 1,500 K and 3 MPa, Y_H2O never turns back (largest fall 0). This covers four starts at one O/F and one pressure; the builder's check stays.
- **(2) The collapse onto one curve holds, but late.** The chamber's composition, frozen, is expanded isentropically to ρ/ρ_c 0.3, 0.1 and 0.03, or cooled at constant ρ to 2,500, 2,000 and 1,500 K. It then reacts at constant (e, ρ) and is compared at equal c with its node's lower branch, at φ = 0, 0.1, 0.5 and 0.9 of the way from its start to c_eq.
  - At φ = 0, which is the snap at equal c, T differs by −71 to −477 K and the largest |ΔY| is 0.017 to 0.071.
  - At φ = 0.5 the states above 2,100 K agree within 0.5 K. The others are still 8 to 187 K apart.
  - At φ = 0.9 all agree within 1.6 K, and the departure from chain-reaction equilibrium, |ln(forward/reverse)|, is the same on both trajectories (0.02 to 0.76). So both trajectories end on the same one-dimensional slow curve. Strict partial equilibrium (departure 0) does not hold on it below about 2,500 K.
  - These starts relax half-way to c_eq in 12 ns to 1.3 µs.
- **(3) A parcel expanding gradually, as in a nozzle.** The chamber's state is expanded as ρ = ρ_c exp(−t/τ), split as the engine splits it: a frozen adiabatic expansion, then reaction at fixed (e, ρ), 200 splits per e-folding. At ρ/ρ_c 0.3, 0.1, 0.03 and 0.01 the parcel is compared with its node's lower branch. The C1 throat sets τ of order 10 µs (*derived*: r_t / u ≈ 0.01 m / 1.5 km/s).

  | τ | snap at equal c: ΔT at ρ/ρ_c 0.3, 0.1, 0.03, 0.01 | snap at equal moles: ΔT at the same | moles' rate, branch / parcel − 1, at equal moles: largest |
  |---|---|---|---|
  | 3 µs | +0.32, +0.02, −11.1, −245 K | +0.01, +0.03, +0.49, +3.2 K | 2.0% |
  | 30 µs | +0.02, +0.05, −0.46, −11.6 K | 0.00, 0.00, +0.01, +0.23 K | 0.26% |
  | 300 µs | 0.00, 0.00, −0.01, −0.53 K | 0.00, 0.00, 0.00, +0.01 K | 0.03% |

  - *What the snap at equal c gets wrong* (*derived*). As T falls, the fast chain reactions shift H2O (OH + H2 = H2O + H moves toward H2O). This changes c, but it conserves the moles per kilogram, N = Σ Y_k / W_k. A snap at fixed c undoes that shift. A snap at fixed N keeps it, since only the slow recombinations change N. In the freezing nozzle (low ρ, fast expansion) the shift is large, and the snap at equal c is 245 K off where the snap at equal N is 3.2 K off.
  - *What λ measures.* Near partial equilibrium, dY_H2O/dt is dominated by the small net rate of the fast chain reaction OH + H2 = H2O + H. So the parcel's instantaneous λ differs from the branch's by up to 98%, even where T agrees within 0.05 K. The branch's λ is the rate along the slow curve, and the update is consistent with it by construction. The well-conditioned comparison is the moles' rate dN/dt, which agrees within 2.0% (above).
  - Moles along the ignition branches (from (1)) never rise. They fall by 1% of their total by the time c reaches 0.5 to 0.8% of c_eq. Before that (the induction period) N is flat, and only c resolves the state.
- **What this changes in the design.**
  - The claim in *Why one variable can carry the late phase* holds for the slow curve, not for strict partial equilibrium. The start of the lower branch matters until φ is about 0.5 above 2,100 K and about 0.9 below.
  - **The snap, amended before any code (6 October 2026, 02:47 CDT).** The snap keeps the cell's N, e, ρ and elements, not its c, wherever the node's branch has lost at least 1% of its fall in N. Before that point (induction) it keeps c, as stated above. The handover between the two is a design item settled in the builder: a single rule continuous in the cell's state, written here before it is coded. N is linear in the transported Y, so it needs no new transported variable. The table still indexes and advances by c, which is monotone and resolves induction.
  - **The handover, stated before code (6 October 2026, 02:59 CDT).**
    - *One curve per node.* The lower branch runs from its start c_0 up to c_eq, and the upper branch from c_eq up to complete combustion, c_cp. Joined at c_eq they make one curve, parametrised by c. Along it, N is flat in induction and then falls to N_eq at c_eq. It keeps falling on the upper branch, since more H2O means fewer moles. N needs no storage: it is linear in the stored Y.
    - *The weight.* φ_N(c) = (N_0 − N(c)) / (N_0 − N_eq) is the share of the lower branch's fall in N reached at c (it is above 1 on the upper branch). The cell's weight is w = clamp((φ_N(c_cell) − 0.005) / 0.005, 0, 1). So w is 0 until the branch has lost 0.5% of its fall, and 1 from 1%, as the amendment states.
    - *The target.* c_N is the point on the curve where N equals the cell's N. It is sought only beyond c_a, the point of 0.5% fall, and clamped to [c_a, c_cp]. The snap goes to c* = (1 − w) c_cell + w c_N, and takes the curve's composition there.
    - Every step is continuous in the cell's c and N, so the snap is continuous in the cell's state. The cell's e, ρ and elements are kept, and T follows from e.
    - *A builder gate, added to 2(a).* N must fall strictly along the curve from c_a to c_cp at every node. If it does not, the builder stops and the design comes back here.
    - *Reported in criterion 5:* the share of snaps with 0 < w < 1, beside the |ΔT_snap| distribution.
    - *Derived, then tested:* this rule also covers burnt gas mixed into fresh gas (PaSR, numerical diffusion), which leaves the curve at an intermediate c and N. In hot gas the chain reactions that carry the mixture back to the curve conserve N, so the equal-N point is where it rejoins. Probe (4), below, found this holds only for mixtures above about 70% burnt.
  - **Criterion 3 gains (d), stated now (02:47 CDT).** Gradual-expansion parcels: 100 chamber states (equilibrium at 2,500 to 3,600 K, 0.5 to 10 MPa, O/F 2 to 16), each expanded as above at τ = 3, 30 and 300 µs to ρ/ρ_c 0.01. Table B's update with the engine's splitting is compared with CVODES with the same splitting. Pass: at ρ/ρ_c 0.1, 0.03 and 0.01, T within 3 K and Y_H2O within 1e-3, at the 99th percentile. 1e-3 in Y_H2O is about 0.04% in Isp_vac (*derived*, linearly: from frozen to shifting, Y_H2O at ρ/ρ_c 0.01 moves by about 0.11 in the probe (0.854 to 0.965) and C1's 1-D Isp_vac by 3.9% (419.4 to 436.0 s)).
  - 3(b)'s sudden expansion stays as stated. The probe's (2) suggests that 3(b) will fail at φ near 0 under any one-curve table. If it does, that is reported as the method's limit for an instantaneous jump, which a flow step of about 70 ns does not make.

## Probe (4): burnt gas mixed into fresh gas (6 October 2026, 03:01 CDT)

Exploratory, like (1) to (3). The chamber's equilibrium gas (3,485 K) is mixed with fresh reactants (300 K, 3 MPa) at a burnt mass share m. Y, e and 1/ρ are mass-weighted. The truth is CVODES from the mixture at constant (e, ρ). It is compared with the lower branch of the mixture's own node (the same elements, e and ρ), built as stated. Output: `table_b/probe4_mixing_2026-10-06.txt`.

| m | T_mix | c_mix | the node's lower branch | truth: time to φ 0.9 | snap at equal c: ΔT | snap at equal N: ΔT; largest \|ΔY\|; T − T_truth at the same c | equal N: time to φ 0.9 / truth's | branch − truth in T at φ 0.9 |
|---|---|---|---|---|---|---|---|---|
| 0.03 | 421 K | 0.026 | reactants at 271 K; not ignited by 1e15 s | 6e8 s | c above the branch | c above the branch | | |
| 0.10 | 699 K | 0.085 | reactants at 201 K; not ignited by 1e15 s | 17 ms | c above the branch | c above the branch | | |
| 0.20 | 1,080 K | 0.171 | frozen equilibrium at 200 K, c_0 0.490 | 276 ns | c below c_0 | −66 K; 0.61; −1,723 K | 0.017 | −184 K |
| 0.30 | 1,437 K | 0.256 | the same, c_0 0.498 | 30 ns | c below c_0 | −35 K; 0.51; −1,351 K | 0.27 | −125 K |
| 0.50 | 2,088 K | 0.427 | the same, c_0 0.513 | 28 ns | c below c_0 | +69 K; 0.33; −749 K | 0.67 | −35 K |
| 0.70 | 2,676 K | 0.598 | the same, c_0 0.527 | 41 ns | −1,469 K | +27 K; 0.20; −331 K | 0.86 (equal c: 1.20) | −5.4 K |
| 0.90 | 3,223 K | 0.769 | the same, c_0 0.541 | 65 ns | −574 K | +1.9 K; 0.068; −79 K | 0.89 (equal c: 1.52) | −0.6 K |

At m 0.7 and 0.9 the stated rule's weight is 1, so its snap is the equal-N snap. Below 0.7 the rule is undefined, because the cell's c is off the curve.

- **Mixtures below about 70% burnt are on no curve the design builds.** There are two ways this happens.
  - *Cold nodes (m 0.03 and 0.1).* The node's e can be carried by reactants at 200 K or above, so its lower branch starts from bare reactants, at 201 and 271 K. That branch never ignites: c stays at 0 through 1e15 s. The mixture carries products and radicals and does ignite, though slowly (17 ms at m 0.1; 6e8 s at m 0.03). Its c is above everything the branch reaches.
  - *Low-e nodes (m 0.2 to 0.5).* The node's e is too low for bare reactants at 200 K, so its lower branch starts from a frozen equilibrium at 200 K. That start is radical-laden, with c_0 of 0.49 to 0.51, above the mixture's c. Both states have the same elements, e and ρ. The mixture holds its chemical energy in unburnt H2 and O2; the branch's start holds it in radicals. The equal-N snap moves mass fractions by up to 0.61. Compared at the same c, it is 750 to 1,720 K from the truth, and the branch is still 35 to 184 K off at φ 0.9.
- **Hot mixtures (m 0.7 and 0.9): the stated rule works where the equal-c snap does not.**
  - T at the snap is within 27 and 1.9 K (equal c: −1,469 and −574 K).
  - The time to φ 0.9 is within 14% and 11% (equal c: +20% and +52%).
  - At φ 0.9 the branch and the truth agree within 5.4 and 0.6 K.
  - The snap moves c forward, from 0.60 to 0.76 and from 0.77 to 0.82. So the state is labelled ahead of the truth: compared with the truth at the same c, it is 331 and 79 K cooler. The derived claim holds in hot gas only in this weaker sense.
- **The largest |ΔY| separates the two kinds of snap.** The nozzle's equal-N snaps (probe 3) move mass fractions by 3e-10 to 3.6e-3. The mixtures' snaps move them by 0.068 to 0.61. The smallest mixture jump is 19 times the largest nozzle jump.
- **The cold-node build assumption is wrong.** *Size and build* expects cold nodes to reach equilibrium if integrated with no time limit. The probe's two reactant branches, at 201 and 271 K, did not ignite by 1e15 s. The injected reactants at 300 K are expected to behave the same (*derived*: the 3% mixture at 421 K, which already carries radicals, takes 6e8 s). Gate 2(c), which requires zero such trajectories in C1's range, would therefore fail on C1's own injected gas.

**The design, amended before any code (6 October 2026, 03:10 CDT).**
- **When a cell is represented.** The table handles a cell only when (i) its c is within the range its node's curve reaches, from c_0 to c_cp, and (ii) the snap, by the stated rule, moves no mass fraction by more than ΔY_max = 0.01. Otherwise the cell goes to direct CVODES, counted by cause.
  - ΔY_max is *guessed*. It sits between the nozzle's largest jump (3.6e-3) and the hottest mixture's (0.068). It is not tuned after a run: criterion 5's audit reports the |ΔT| of the snaps just under it, and a fail there brings the design back here.
  - This replaces the remedy in *States off the manifold*, which was to build the fallback only if needed. The fallback is now in the design from the start.
- **Frozen branches.** A lower branch that has not reached equilibrium by 1e15 s is kept as far as it went and flagged. Within its reach, the table handles the cell: the gas is frozen, and the branch's λ says so. Above its reach, rule (i) sends the cell to CVODES.
- **Gate 2(c), restated.** The builder counts and reports the frozen branches, by node. They are no longer required to be zero in C1's range, because C1's injected gas sits on them. 2(c) now requires that every branch either reaches equilibrium or is flagged frozen, with no integrator failure.
- **Criterion 3 gains (e), reported and not judged.** Burnt-fresh mixtures: 200 of 3(d)'s chamber states, each mixed with fresh reactants at 300 K at m drawn uniformly from 0 to 1. The table's path, with its CVODES fallback, is compared with CVODES. Reported: the fallback share against m, the time to φ 0.9, and T at φ 0.9.
- **Criterion 5 also reports** the CVODES share by cause: c off the curve, ΔY above ΔY_max, a frozen branch exceeded, or outside the table. It also reports the share of the step's cost each cause takes. Criterion 6 then decides whether the fallback is affordable.
- **The alternative, not built unless criterion 6 fails on the fallback's cost.** Add N, the moles per kilogram, as a second table coordinate. Each node would carry trajectories started along the burnt-fresh mixing line, giving a two-parameter manifold. It is about ten times the size at the present axes (*derived*: about ten starts per node), which is over the 1 GB limit, so other axes would have to be coarsened.

## Criteria

Stated before code, 6 October 2026, about 02:34 CDT.

1. **Conservation.**
   - At every node, the stored composition has the node's element mass fractions to 1e-12, and |Σ Y − 1| ≤ 1e-12.
   - Every lookup in criteria 3 and 5 returns the query's element mass fractions to 1e-12, with |Σ Y − 1| ≤ 1e-12.
   - In the C1 runs of criterion 5, the mass and energy budgets are reported against CHAMBER_C1's limits.
2. **The build** (gates, checked by the builder before the table is written).
   - (a) The progress variable passes the monotonicity check on every trajectory. The variable used is reported.
   - (b) Every node's equilibrium end state agrees with direct Cantera equilibrium at the same elements, e and ρ to |ΔT| ≤ 3 K. This is Table A's criterion and reason: c* to 0.05%.
   - (c) The count of trajectories that did not reach equilibrium is reported, and must be zero for the C1 range (f_H at O/F 4 to 6, all of Z_N, the e and ρ that C1 visits).
3. **Homogeneous reactors against CVODES.** The table's own update (the runtime path) and CVODES on the full mechanism are run from the same start state, at constant (e, ρ). States are drawn off the nodes, uniformly in the stated ranges.
   - (a) **Ignition.** 1,000 reactant states with O/F 2 to 16, initial T 900 to 2,000 K, p 0.1 to 10 MPa, and Z_N 0 and 0.3. The ignition delay is the time of the largest dT/dt. Pass: 99th percentile of |error| ≤ 5%. The largest error is reported.
     - 5% is small against the mechanisms' own disagreement at rocket pressures, which is tens of percent (*guessed*; the mechanism papers will confirm or correct it).
   - (b) **Recombination.** 1,000 states taken from equilibrium at 2,500 to 3,600 K and 0.5 to 10 MPa, frozen, then set to an expanded (e, ρ) with ρ falling 2 to 20 times isentropically, as in a nozzle. Pass: the time for c to cover half the distance to c_eq within 5% (99th percentile), and T at 10 of those half-times within 3 K.
   - (c) Reported: the largest |ΔT| along each trajectory after ignition, per range.
4. **The free flame** (FLAME_C2's case: h2o2.yaml, mixture-averaged transport, dz = 10 µm). The table's S_c is compared with the direct-chemistry run's S_c.
   - Pass: within 5%.
   - A fail is reported as the method's limit for resolved flames: preferential diffusion takes the flame off the homogeneous-reactor curve. Criterion 5 then still decides use in the chamber, where cells are a thousand times the flame thickness.
5. **The C1 finite-rate cold start**, 64x12, from the same tree on the same machine and threads. The table runs against direct CVODES, both with the PaSR closure (chamber_study `frt`, PASR_C2 criterion 4) and both without it.
   - Pass: settled c*, vacuum Isp and vacuum thrust within 0.1%.
   - Pass: light-off time (as defined in PASR_C2) within 5%.
   - Pass: injector pressure within 1% at every printed time from 1 ms.
   - Reported: an audit of every Nth reaction call against CVODES from the same pre-reaction state (|ΔT| distribution, N stated per run), the |ΔT_snap| distribution, and the fallback counts.
6. **Cost.** Table B's 64x12 finite-rate cold start to 4 ms is timed against direct CVODES.
   - Pass: 5 minutes or less on the Mac (the *Lightweight engine* target).
   - Reported: the speedup, the table's share of the step, the memory, and the build time on backhouse.

The resolution per axis is the setting that criteria 2 to 5 select, reported with the result. The table file is reproducible and not kept in git. The builder prints its provenance: mechanism, axes, progress variable, Cantera version, build time and the node conservation error.

## Not covered

- Liquid oxygen (the RL10's LOX enters below 200 K, as a liquid). Table B is gas-phase, like the engine.
- Mechanisms with nitrogen chemistry. As in Table A, N2 is a diluent and a third body only. The builder refuses a mechanism in which N2 reacts.
- Turbulent fluctuations of c inside a cell beyond what PaSR's κ represents. There is no presumed-PDF table over c.
