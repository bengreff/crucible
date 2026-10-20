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

## The builder, stated before its code (14 October 2026, about 01:00 CDT)

What the design above leaves to the builder, settled here before the builder is written. Nothing here changes a criterion.

- **Integration.** The engine's own CVODES reactor (`ReactionStep`, the same rates as the fallback) gains a one-step mode that shows the builder every internal step of one uninterrupted integration. Every branch, not only the cold ones, uses rtol 1e-10 and atol 1e-30, up to 1e15 s. A branch stops when |c − c_eq| ≤ 1e-6 of its span |c_eq − c_start|; one still short at 1e15 s is frozen (as amended at 03:10).
- **The branches' starts** are as stated: the lower branch from the reactants (H2, O2 and N2 at the node's elements) at the node's (e, ρ) when they carry e at 200 K or above, otherwise from the equilibrium at the node's elements and ρ at the highest T* whose frozen composition carries e at 200 K (bisection, as in the probe); the upper branch from the complete-combustion products. c_eq comes from `equilibrateUV` of the start.
- **The c nodes, relative to the branch's own start.** The lower branch's nodes sit at s = (c − c_0)/(c_eq − c_0): s = 0, then 12 nodes geometric from 1e-9 to 0.1, then 18 uniform to 1 (c_eq). With the reactants' start (c_0 = 0) this is the stated axis. With the other start, c_0 > 0, so the stated axis would put nodes below the curve's reach. The upper branch has 9 nodes uniform from c_eq to c_cp. 40 points per node, as guessed. A frozen lower branch places the same fractions on [c_0, c_reach].
- **Values at a node** are interpolated linearly in c between the two internal steps that bracket it (c is monotone there, by gate 2(a)). λ = (dc/dt)/(c_eq − c) uses the exact rate at each step. At c_eq, where λ is a 0/0, the node takes the λ of the node before it (the relaxation rate is finite and smooth there, derived).
- **Exact elements.** Each stored composition is projected onto its node's elements as in Table A (smallest change weighted by Y_k). The conservation error is printed before and after.
- **The progress-variable candidates** are linear in Y: Y_H2O, then Y_H2O + Y_HO2 + Y_H2O2, then the O2 consumed (−Y_O2 plus a constant). Gate 2(a) measures each trajectory's largest turn-back against the approach, relative to its span, until c is within 1e-4 of the span from c_eq. Anything above 1e-9 (CVODES's noise at rtol 1e-10) is a turn-back. The builder uses the first candidate with no turn-back anywhere, or stops.
- **The N gate** (03:10, 2(a)) is checked on the stored points of each node's joined curve: N strictly falls from c_a to c_cp.
- **Nodes c cannot describe.** At f_H = 0 (oxygen alone) and f_H = 1 (hydrogen alone) no H2O can form, so c ≡ 0 while the gas still dissociates. Such a node stores the equilibrium at its (e, ρ) at every point, with c_0 = c_eq = c_cp = 0, and is flagged degenerate. A branch of zero span (|c_eq − c_start| < 1e-12) is stored the same way. How the runtime treats these flags is stated with the runtime, before its code.
- **The thermo range.** A node whose branch passes above the mechanism's highest thermo temperature is built anyway, since the NASA fits extrapolate smoothly. It is flagged outside, so the runtime sends it to CVODES (as *Outside the table* states).
- **Output.** Per node: c_0, c_eq, c_cp, the flags (frozen, degenerate, zero span, outside), the gates' numbers and the end-state |ΔT| against `equilibrateUV`. Per point: the H/O species' Y and λ. The file layout is `core/manifold_table.hpp`. A small table (f_H 5, Z_N 2, η 6, ln ρ 5) exercises every gate first. The full axes are then built on backhouse.

## The 3,500 K ceiling: thermo to 6,000 K (20 October 2026, about 18:00 CDT)

The first small table (300 nodes, 14 October) failed every gate, and the failures cluster where the energy axis puts states above h2o2.yaml's 3,500 K thermo range (the builder agent's diagnosis, *Stopped at* in the handoff of that date). η = 1 is the complete-combustion products frozen at 6,000 K, so the axis reaches far past 3,500 K by design. **Decision (the Director, 20 October): extend the thermo to 6,000 K with sourced data**, not shrink the η range. Table A's T axis gets the same extension.

- **The data.** `tools/mechanisms/h2o2_highT.yaml`, written by `tools/mechanisms/make_h2o2_highT.py`: Cantera's h2o2.yaml (reactions, transport) with every species' thermo replaced by the NASA 7-coefficient fits, 200 to 1,000 to 6,000 K, of McBride, Gordon and Reno, NASA TM-4513 (1993), as shipped in Cantera's `nasa_gas.yaml`. Each species carries its source in its `note`. The engine reads NASA-7 only (`adapters/reaction.cpp`), so nothing in the core changes. Ar's single range is written as two identical ranges (the same function).
- **Species check against CEA's own data** (measured; `tools/mechanisms/check_h2o2_highT.py` against the NASA-9 `thermo.inp` shipped with RocketCEA 1.2.3, i.e. McBride, Zehe and Gordon, NASA TP-2002-211556; output `table_b/thermo_highT_check_2026-10-20.txt`):
  - h2o2_highT, 3,500 to 5,000 K: cp within 0.15% for every species; h within 0.05 kJ/mol for H2, H, O, O2, H2O, N2, Ar. 5,000 to 6,000 K: cp within 0.3%.
  - h2o2.yaml (GRI) extrapolated past its range: cp off by up to 15% at 3,500 to 5,000 K and 62% at 5,000 to 6,000 K (H2, O2, H2O, OH). This is what the first table's outside nodes were built on.
  - **Shared by both files: OH's enthalpy is 2.07 kJ/mol from CEA's at every temperature** (g/RT 0.07 at 4,000 K, about 7% in OH's formation constant), and HO2 and H2O2 by 0.4 to 0.56 kJ/mol. CEA's OH cites Ruscic (2002) for D0(H-OH), the revised heat of formation; GRI and TM-4513 predate it. So the engine has carried this offset since Cantera was first matched to CEA. It is not corrected here: the equilibrium check below measures what it costs, and the correction (adopting the revised heat of formation) would be a stated decision of its own.
- **The equilibrium check, stated before it runs.** CEA (RocketCEA, HP chamber problem) against Cantera with h2o2_highT, same reactants: gaseous H2 and O2 at 300, 1,500 and 3,000 K, O/F 1 to 12, 0.1 to 20 MPa, keeping the cases whose CEA chamber temperature is 3,000 to 5,000 K. Pass: |ΔT| ≤ 3 K (Table A's criterion and reason, c* to 0.05%). The major mole fractions are reported. The same cases with h2o2.yaml are printed beside, not judged.
- **Then**, in order: the builders take the mechanism as an argument (default h2o2.yaml, printed in the provenance); Table A's T edges extend to 6,000 K and Table A is rebuilt with h2o2_highT and re-checked by its criterion 2(a); the small Table B is rebuilt and its gate table reported against the gates as written. A run that uses either table must run its CVODES fallback and core thermo on the same mechanism, so a driver that loads a table built on h2o2_highT must load h2o2_highT too.

## OH's heat of formation, the reference pressure, and the equilibrium check (20 October 2026, about 17:30 CDT)

- **Decision (the Director, 20 October): adopt the revised OH heat of formation**, Ruscic et al. (2002), J. Phys. Chem. A 106, 2727, which CEA uses: the better-measured value, and what the engine is validated against; TM-4513's is an outdated number, not a modelling choice. `make_h2o2_highT.py` shifts a6 of both OH ranges so that h(298.15 K) = 37,278.206 J/mol, the value in CEA's `thermo.inp` (TM-4513's is 39,346.882 J/mol); cp and s are unchanged. The species check rerun (`table_b/thermo_highT_check_ohfix_2026-10-20.txt`): OH's h is now within 0.034 kJ/mol of CEA's below 3,000 K and 0.011 above. HO2 and H2O2 keep their 0.4 to 0.56 kJ/mol offsets (not in the decision; minor species, see the check below).
- **Found while running the check: the fits were read at the wrong reference pressure** (measured; `table_b/cea_check_reference_pressure_probe.txt`). With OH corrected, the first run failed: Cantera 1.5 to 3.6 K below CEA, worst −3.62 K, 17 of 153 cases above 3 K, all from 3,000 K reactants at CEA Tc above 4,000 K. Replacing every species' thermo with CEA's own NASA-9 (`diag_nasa9_swap.py ALL`) left the same gap (−1.70 K, −3.32 K on two probes), so it was not the data. At CEA's own state (1,500 K reactants, O/F 6, 1 MPa) the enthalpy agreed to 0.17 kJ/kg, but Cantera's TP equilibrium at CEA's T was more dissociated (X_H 0.10922 against 0.10850). NASA Glenn's data (TM-4513 and CEA's TP-2002-211556) are for a 1 bar standard state, and Cantera's `nasa_gas.yaml` declares no reference pressure, so Cantera read them at its 1 atm default: each Kp off by 1.01325^Δn. Declared at 1 bar, CEA's own data reproduce CEA's T to 0.04 K on all four probes. `make_h2o2_highT.py` now writes `reference-pressure: 1.0e5` on every species. The engine's core never evaluates entropy (only `adapters/thermo.cpp` reads it, through Cantera), so the declaration reaches equilibrium and the reverse rates through Cantera and nothing else changes.
  - **Closed 20 October (TABLE_A.md, *What h2o2.yaml's 1-atm reading moved*):** its reacting species are 1 bar data read at 1 atm. Correcting that moves C1 by about 1e-4 relative, and moves the CEA check's worst from +6.13 to +7.60 K. The note below is kept as written. **Not checked here, an open item:** Cantera's GRI `h2o2.yaml` (the engine's mechanism for every run so far and Table A) also declares no reference pressure. The species check's g/RT differences hint that its fits mix the two standard states (AR and N2 differ from CEA's by about ln 1.01325, H2, O2 and H2O by 3e-4 or less; inferred, not measured). It is beside the check below at worst +6.13 K.
- **The equilibrium check as stated (passes)**, `tools/cea_equilibrium_check.py` driving `crucible_equilibrium_hp` (new, `tools/equilibrium_hp.cpp`), output `table_b/cea_equilibrium_check_2026-10-20.txt`. 198 cases (reactants at 300, 1,500 and 3,000 K; O/F 1 to 12; 0.1 to 20 MPa), 153 with CEA Tc 3,000 to 5,000 K, judged on those.
  - h2o2_highT (TM-4513, OH from Ruscic 2002, 1 bar): **worst |ΔT| 1.79 K** (1,500 K reactants, O/F 2, 20 MPa), bar 3 K. By reactant temperature: worst 0.76 K at 300 K, 1.79 K at 1,500 K, 0.33 K at 3,000 K; 4 cases above 1 K. Worst mole-fraction difference 2.4e-4 (H).
  - **The OH decision's effect**, beside, not judged: the same file with TM-4513's OH, worst +7.42 K (300 K reactants, O/F 12, 20 MPa), mean +5.4 K at 300 K reactants, and X_OH off by up to 1.03e-2. The file as committed in acab338 (TM-4513's OH, read at 1 atm): worst +5.95 K. GRI's h2o2.yaml: worst +6.13 K (the two errors partly cancel in both).

## The small table rebuilt on h2o2_highT: the gate table (20 October 2026, about 17:45 CDT)

The same small axes as on 14 October, now with the command recorded: `crucible_build_manifold_table <out> 2 2 0 0,0.3 6 5 1e-3 30 2 ../tools/mechanisms/h2o2_highT.yaml` (f_H 0, 0.056, 0.112 = f_st, 0.556, 1; Z_N 0 and 0.3; η 0 to 1 in 6; ln ρ over 1e-3 to 30 kg/m³ in 5: 300 nodes, 7 s per candidate on 2 threads). The 14 October command was not recorded, so its ρ range may differ. Output: `table_b/small_highT_2026-10-20.txt`, per node `table_b/small_highT_nodes_2026-10-20.csv`.

**Two builder changes, neither touching a gate.**
- It now prints gates 2(b), 2(c) and criterion 1 for every candidate. Before, it printed them only for the last one tried (the O2 consumed), the worst on the rich side. `CRUCIBLE_TABLE_B_NODE_REPORT=<path>` writes one row per node and candidate.
- Criterion 1's sum is checked with the inert species at their exact values (N2 = Z_N, Ar 0). The table stores only the reacting species and the runtime rebuilds the inert ones from Z_N, so the integration's own N2 is not what is stored. It drifts by about rtol: every one of the 16 nodes with |Σ Y − 1| above 1e-12 (up to 7.7e-11) was at Z_N 0.3 (`table_b/small_highT_before_inert_fix_2026-10-20.txt`).

| Gate (as written) | Y_H2O | Y_H2O + Y_HO2 + Y_H2O2 | The O2 consumed | 14 October (h2o2.yaml, as reported) |
|---|---|---|---|---|
| 2(a) worst turn-back (bar 1e-9) | 5.78e-2 | 5.35e-2 | 2.15e6 | 2.05e5, 2.05e5, 2.67e6 |
| 2(a) N-gate failures (bar 0) | 25 | 25 | 27 | 26 to 27 |
| 2(b) worst \|ΔT\| (bar 3 K) | 64.9 K | 64.9 K | 593 K | 503 K |
| 2(c) integrator failures (bar 0) | 40 | 40 | 42 | 39 |
| 2(c) frozen branches, lower and upper (reported) | 8, 13 | 8, 13 | 6, 6 | 2 of 11 |
| Criterion 1, \|Σ Y − 1\| and elements after projection (bar 1e-12) | 2.2e-16, 2.5e-16 | 4.4e-16, 2.5e-16 | 2.2e-16, 2.5e-16 | 3.6e-11, 2.8e-16 |
| Nodes outside the thermo range | 0 | 0 | 1 | 146 |

The table is not written: gates 2(a), 2(b) and 2(c) fail for every candidate. The thermo ceiling was one cause, not the main one. Nodes outside the range fell from 146 to 0, and the worst Y_H2O turn-back from 2e5 to 0.058, but the failures that remain are where Y_H2O's node report puts them (120 nodes at f_H 0 and 1 are degenerate by design):

- **2(c), integrator failures (40): CVODES flags −3 (error test) and −4 (corrector convergence) at t = 0, at the dense end.** By ρ: 18 at 30 kg/m³, 17 at 2.28 and 5 at 0.17. By f_H: 20 at f_st, 15 at 0.056 and 5 at 0.556. By η: 17 at η = 0, the rest at 0.2 to 0.8. A rerun at atol 1e-20 instead of the stated 1e-30 still has 38, so the tolerance is not the cause (measured; the edit was reverted).
  - Inferred, not yet measured: these are the lower branches that start from the frozen high-T* composition at about 200 K, where a dense radical pool recombines in picoseconds. It is a defect of the builder's start or integration, not of the manifold's physics. The next step is to tag the failing branch and print its start state and the first right-hand side.
- **2(a), turn-back (25 nodes above 1e-9; worst 5.8% of the span).** 17 are at f_H 0.556 (O/F 0.8; η 0.2 to 1) and 8 at f_st (η 0.4 to 1). The N gate fails at f_H 0.556, η 0.6 to 1 (23) and f_H 0.056, η 1 (2). These are hot or very rich states, where water and the radical pool exchange, so a single linear-in-Y variable is not monotone. That is physics, not the ceiling (inferred, as the handoff expected).
- **2(b) (Y_H2O): 2 nodes above 3 K, both at f_H 0.556, η 0.6, ρ 1e-3** (32 and 65 K). There the lower branch reaches c_eq while hydrogen recombination has not finished, so c has stopped describing the state. On the O2 consumed it is 593 K at the same node, because at O/F 0.8 the O2 is gone long before the end.
- **No node of this small table lies in the C1 range** (O/F 4 to 6, f_H 0.143 to 0.2), which criterion 2(c) names. Its f_H nodes are at 0.112 and 0.556. Whether the gates pass where C1 and the RL10 run is still unmeasured.

**For the Director (a design question, not settled here):** the hot and very rich failures of 2(a) are a property of linear-in-Y progress variables. The options are to state the table's domain (f_H and η ranges) with everything outside sent to CVODES, as *Outside the table* already does for the thermo range, or to choose a different progress variable. The 2(c) defect is engineering and is fixed first, then a small table with f_H nodes inside O/F 4 to 6 measures the C1 range.

## The integrator failures: two causes, measured (20 October 2026, about 18:15 CDT)

The Director's order: fix gate 2(c)'s 40 integrator failures first, measuring the cause. The builder now tags a failing branch with its name, how it started (reactants, or the frozen start and its T*), the start's T, first right-hand side and Y, and the furthest step reached (the node report's `failure` column). Same small table and command as above. Tagged run: `table_b/small_highT_tagged_nodes_2026-10-20.csv`. **Every one of the 40 is a lower branch from the frozen start, failing before its first step. They have two causes:**

- **(A) 17 at η = 0: the frozen start's bisection had no sign change (a builder defect, fixed).** At η = 0 the node's energy is the complete-combustion products' at 200 K. That is the frozen family's low end, so g(T*) = e(200 K, frozen at T*) − e is ≥ 0 over the whole bracket [300, 6,000] K. The bisection then drifts to its upper bound: 15 nodes started at T* = 6,000 K, i.e. from the 6,000 K equilibrium at T = 6,000 K, with c_0 = c_eq to 7 digits. The other 2, at f_H 0.556, found round-off roots at 599 and 669 K and started from equilibrium at 200 K. CVODES fails (flags −3 and −4) on both kinds of start.
  - **The fix (an engineering change to the stated start; the design's text is unchanged in intent):** g is evaluated at the bracket's ends before bisecting. If g(300 K) ≥ −1e-9 of the node's energy span, the start is the complete-combustion products at the node's energy (at η = 0, at 200 K exactly; the equilibrium at 300 K differs from them only by trace dissociation). If g(6,000 K) < 0, the builder throws a named error ("no frozen start up to 6000 K carries the node's energy") instead of drifting; it never fired here.
  - **Result** (`table_b/small_highT_bracketfix_2026-10-20.txt`, nodes `..._bracketfix_nodes_2026-10-20.csv`): 0 failures among the 30 non-degenerate η = 0 nodes (10 are zero-span points, 20 frozen at 200 K, as they should be). Frozen lower branches rise from 8 to 25 for Y_H2O for that reason.
- **(B) 23 at η 0.2 to 0.8: the cold dense-radical start (the Director's first suspect; open).** All at f_H 0.056 and f_st, ρ 2.28 and 30 kg/m³, T* 3,400 to 5,460 K. They start at 200 K with the T* radical pool (start Y in the tagged report), where dT/dt is 9e11 to 4e15 K/s, so the recombination time T/|dT/dt| is 5e-14 to 2e-10 s. CVODES fails with flag −4 (corrector convergence) before taking a step.
  - **Probe (measured, reverted): CVODES's own first step is the trigger.** Forcing the initial step to 1e-18 s (`table_b/probe_initial_step_1e-18_2026-10-20.txt`; 1e-22 s beside) lets all but 1 of the 23 integrate. But that step makes an upper branch fail that passed before (f_H 0.056, η 1, ρ 30: complete-combustion products at 6,000 K, dT/dt −1e14 K/s, radicals exactly 0), and 1e-22 s makes two slow reactant starts at ρ 1e-3 fail (−3 after 54 to 60 steps). So a fixed initial step is not the fix.
  - **Inferred, not measured: the difference-quotient Jacobian.** CVODES's dense Jacobian perturbs each species by about √(unit round-off) / its error weight, and at atol 1e-30 a species near 0 is perturbed by about 1e-38 against right-hand sides of 1e7 to 1e15 s⁻¹. Its columns are then round-off divided by the perturbation, and Newton cannot converge where the chemistry is fastest. That fits both the cold starts and the new upper-branch failure, and the 14 October result that atol 1e-20 changed little (38 failures).
  - **Next (engineering, not a design question):** give the builder's `follow()` an analytic Jacobian (Cantera 3.2's analytic derivatives of the net production rates, plus T and the constant-volume energy), or scale the perturbation; then rerun the small table. If that is not enough, an initial step tied to the start's own time scale (a stated rule, for example 1e-6 of T/|dT/dt|). The runtime's CVODES fallback uses the same Jacobian, so the change is checked with `ctest -R reaction` and C1 before it is kept.
- **Gates otherwise unchanged by fix A:** Y_H2O's worst turn-back 5.78e-2 and 25 N-gate failures, 2(b) 64.9 K; the O2 consumed's worst turn-back falls from 2.15e6 to 7.04 (its worst sat on an η = 0 node).

## The runtime, stated before its code (14 October 2026, about 01:20 CDT)

- **Why the corners cannot be read at the cell's c** (*derived*). Lean of stoichiometric, c_cp = 9 Y_H (all the hydrogen as water), so it rises with f_H. In cold burnt gas, c_eq is close to c_cp. A cell at its own c_eq, between f_H nodes f_1 < f < f_2, then has c above c_cp(f_1): the leaner corner's whole curve ends below the cell's c. The rich side is the mirror image. So *Y is interpolated linearly in c at each corner* would extrapolate at the nozzle exit and in the cooled wall layer, which are the places the table is for.
- **The rule instead.** Every node stores its points at the same fractions: s = (c − c_0)/(c_eq − c_0) on the lower branch (s from 0 to 1), and s = 1 + (c − c_eq)/(c_cp − c_eq) on the upper branch (s from 1 to 2). The 16 corners are blended point by point at these common s with the multilinear weights of (f_H, Z_N, η, ln ρ). That gives the cell one curve of 40 points. Y, λ and c_eq on it are linear in s between points. The cell's s is where the blended curve's c equals the cell's c. The blended c is monotone in s, since each corner's is and the weights are positive, so the inversion is a segment search. The design's properties hold: the curve returns Y_H2O = c (or the chosen c) exactly, and the elements bilinearly. The query's η uses its own complete-combustion products, as in Table A's check.
- **When the table handles a cell** (the 03:10 rule, made concrete).
  - The table handles a cell when: every corner exists and none is flagged outside; e, ρ and the elements are inside the axes; AR ≤ 1e-12; the cell's c lies on the blended curve; and the snap moves no Y by more than ΔY_max.
  - A frozen corner places its lower-branch fractions on [c_0, c_reach]. If any corner is frozen and the cell's s, or a sub-step, would pass 1 on the lower branch, the cell goes to CVODES ("frozen reach").
  - Degenerate corners (f_H = 0 or 1) enter the blend with their equilibrium at every s. A cell whose every weighted corner is degenerate (pure O2 or pure H2) is set to the blend, i.e. equilibrium. This is counted, and its share of cells and the largest T reached are reported in criterion 5.
- **The step.** Sub-steps never cross a point of the blended curve. Each is the stated exponential update toward the blended c_eq, with λ at the midpoint by predictor-corrector.
  - Under PaSR, λ is multiplied by κ_eff = 1 − s_mix (1 − κ), as the engine's closure multiplies dz/dt. κ comes from the same chemical time τ_c = max_i Y_i/|dY_i/dt| over the same species. Here dY_i/dt = (dY_i/dc) λ (c_eq − c), with the slope taken along the blended curve's segment.
  - The composition is then the curve's at the new c, with e and ρ kept, and T comes from the core's exact inversion.
- **The snap** is the 02:59 handover, on the blended curve. φ_N, c_a and c_N come from N = Σ Y_k/W_k of the blended points.
- **Counted per call:** table, and CVODES by cause (outside the axes, outside flag, AR, off the curve, ΔY above ΔY_max, frozen reach), plus the cells set by degenerate corners alone. This is criterion 5's breakdown.

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
