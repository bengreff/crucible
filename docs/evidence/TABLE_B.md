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

## The analytic Jacobian, stated before its code (20 October 2026, about 18:15 CDT)

The Director's order: give the integration an exact Jacobian in place of CVODES's difference quotients, and check it before keeping it. `ReactionStep` (`advance` and `follow`) is the builder's integrator and the runtime's CVODES fallback, so the change reaches both.

- **The Jacobian of z = [T, Y_1..Y_K] at fixed ρ** (*derived*). With C_j = ρ Y_j / W_j, f_k = W_k ω_k / ρ and f_T = −Σ u_k ω_k / (ρ c_v):
  - ∂f_k/∂Y_j = (W_k / W_j)(D_kj + R T ∂ω_k/∂P), with D = Cantera's `netProductionRates_ddCi()` (each C_j alone, T and P held), and the P term because at fixed volume P = R T Σ C.
  - ∂f_k/∂T = (W_k / ρ)(∂ω_k/∂T + (P/T) ∂ω_k/∂P), from `getNetProductionRates_ddT` (concentrations held) and `_ddP`.
  - ∂f_T/∂Y_j = −(1 / (W_j c_v)) [Σ_k u_k (D_kj + R T ∂ω_k/∂P) + f_T c_v,j], with c_v,j molar; c_v = Σ Y_k c_v,k / W_k is linear in each Y (Cantera's unnormalised mass fractions).
  - ∂f_T/∂T = −(Σ_k c_v,k ω_k + Σ_k u_k ∂ω_k/∂T|_V) / (ρ c_v) − f_T (∂c_v/∂T) / c_v, with ∂c_v,k/∂T from the NASA-7 coefficients.
  - Cantera's settings: `skip-third-bodies` and `skip-falloff` false, so third-body and falloff terms are in. Cantera itself takes the falloff [M] derivative and the T derivatives of falloff and of the equilibrium constants by a relative 1e-8 perturbation (well scaled: T and [M] are never near zero). Everything else is analytic.
- **Under the PaSR closure** (segregation > 0) the Newton matrix is κ_eff(z) times the Jacobian above: κ_eff is held over the Newton iteration and its own derivative is left out. The right-hand side is unchanged, so this changes convergence, not the solution (CVODES needs only an approximate Newton matrix). With segregation 0 it is the Jacobian above.
- **The difference quotient stays available** (a constructor option) for the comparisons below. The analytic one becomes the default.

**Checks, stated before the code:**

1. **It matches finite differences at well-scaled states.** States with every species at Y ≥ 1e-3 (renormalised random mixes, 6 per point), at T 800, 1,500, 2,500 and 3,500 K and ρ 0.05, 1 and 20 kg/m³, on h2o2.yaml and h2o2_highT. Reference: central differences of `rates` with steps 1e-6 of |z_j|. Pass: for every column j, max_i |J_ij − J_ij^FD| ≤ 1e-6 max_i |J_ij^FD|. The worst column and the worst entry above 1e-3 of its column's largest are reported. In `crucible_reaction_tests`.
2. **The cold dense-radical starts integrate.** The small table as before (`crucible_build_manifold_table <out> 2 2 0 0,0.3 6 5 1e-3 30 2 ../tools/mechanisms/h2o2_highT.yaml`): the 23 lower branches that failed (f_H 0.056 and f_st, η 0.2 to 0.8, ρ 2.28 and 30 kg/m³; `small_highT_bracketfix_nodes_2026-10-20.csv`) and the 6,000 K upper branch the forced first step broke (f_H 0.056, η 1, ρ 30) all integrate, and the table has 0 integrator failures in all. Its gate table is reported in full.
3. **Nothing that passed regresses.** `ctest -R "reaction|pasr"` passes as before.
4. **C1 is unchanged within the integrator's tolerance.** C1 finite-rate, inviscid, 32x6 (`crucible_chamber_study fr 32 6 <end> 2 <prefix>`, h2o2.yaml, rtol 1e-6, atol 1e-12), the analytic build against the difference-quotient build of the parent commit, same machine and threads. Pass: c*, vacuum Isp and vacuum thrust at the end within 1e-4 relative, the light-off time within 1%, and the injector pressure within 1e-3 relative at every printed time. Reported: the reaction's wall time and steps in both. The end time is chosen from a timed short run and stated before the comparison.
   - **End time, stated before the comparison (18:12 CDT):** 4 ms, with C1's 300 J igniter (`crucible_chamber_study fr 32 6 4e-3 2 <prefix> 300`). The timed run (difference quotient, to 1 ms) took 221 s per ms once lit, so about 15 min per run.

### Results (20 October, 18:25 CDT)

**Check 1, as stated: FAIL** (measured). The worst column was 1.402e-6 (h2o2.yaml, 1 column of 792) and 2.911e-6 (h2o2_highT, 3 columns), against 1e-6, all in the T row. Two findings on the way:

- **A Cantera 3.2 bug (measured, reproduced in plain Cantera).** After a composition change at unchanged T and ρ, `getNetProductionRates_ddT` returns stale values. The falloff species' rows (H, O2, OH, HO2, H2O2) were wrong by up to 5e4 times, and the first run of the check failed at 1.3e7 to 2.4e7. `invalidateCache()` before the derivative calls gives the fresh values exactly, and the code now does that.
- **What remained was round-off, in two places (measured by scans).**
  - The reference's own: its worst column falls with the step as round-off over the step does. At steps of 1e-7, 1e-6 and 1e-5 it was 1.9e-5, 1.4e-6 and 0 on h2o2, and 1.3e-5, 2.9e-6 and 1.0e-6 on highT.
  - Cantera's: the reference was then made converged (Richardson-extrapolated central differences at 1e-5, whose values at h, 2h, 4h and 8h change by a factor of 4 each time). The one column left was ∂f_T/∂T at 2,500 K and ρ 20 (1.02e-6), and it moves with Cantera's one-sided step `rtol-delta`:

| `rtol-delta` | 1e-10 | 1e-9 | 1e-8 (default) | 3e-8 | 1e-7 | 3e-7 | 1e-6 |
|---|---|---|---|---|---|---|---|
| highT worst column | 1.4e-4 | 1.6e-5 | 1.02e-6 | 6.1e-7 | 4.4e-7 | 1.3e-6 | 4.4e-6 |
| h2o2 worst column | | | 4.9e-7 | 2.3e-7 | 3.1e-7 | 9.4e-7 | |

- **Restated after the result (flagged): Cantera's step is 1e-7, not 1e-8, and the reference is the Richardson one at 1e-5.** The bar stays at 1e-6. The worst column is then 3.09e-7 (h2o2) and 4.38e-7 (highT): pass.

**Check 2, as stated: FAIL** (measured; `table_b/small_highT_analyticjac_2026-10-20.txt`, nodes `..._analyticjac_nodes_2026-10-20.csv`). The 6,000 K upper branch integrates, but the same 23 lower branches fail (the same set as the bracket-fix run). Each fails with flag −4 at t = 0, and the last step sizes (1.25e-12 to 9.34e-12 s) match the difference-quotient run's to every printed digit. The gate table:

| Gate | Y_H2O | Y_H2O + Y_HO2 + Y_H2O2 | The O2 consumed |
|---|---|---|---|
| 2(a) worst turn-back (bar 1e-9) | 5.78e-2 | 5.78e-2 | 7.04 |
| 2(a) N-gate failures (bar 0) | 25 | 25 | 25 |
| 2(b) worst \|ΔT\| (bar 3 K) | 65.4 K | 65.4 K | 132.4 K |
| 2(c) integrator failures (bar 0) | 23 | 23 | 23 |
| 2(c) frozen branches, lower and upper | 25, 20 | 25, 20 | 15, 10 |
| Criterion 1 after projection (bar 1e-12) | 2.2e-16, 2.4e-16 | 4.4e-16, 2.6e-16 | 4.4e-16, 3.4e-16 |

- The O2 consumed's 2(b) moved from 592.7 K (f_H 0.556, Z_N 0.3, η 0.6, ρ 1e-3) to 132.4 K (f_H 0.556, Z_N 0, η 0.2, ρ 0.013). Inferred: at O/F 0.8 the O2 runs out long before the end, so that variable's end point is ill-conditioned and moves with the integrator's path.

**The mechanism, written down before more code (rule 4: two fixes have missed it, the forced first step and the Jacobian).** *Inferred, from CVODES 5.3's step logic:*
- The Jacobian was not the cause. The last step sizes are the same with either Jacobian because the failure path does not depend on it.
- CVODES estimates a first step. At every corrector failure it cuts h by 4, and it gives up at the 10th failure (`maxncf`, default 10). So the last h is the first h × 4⁻⁹, and the first step was about 3e-7 to 2.4e-6 s (*derived*). The recombination time is 5e-14 to 2e-10 s, so that is 3 to 7 orders too long.
- At the last h, h/τ is still 10 to 1e2 at the fastest starts. One implicit step must then carry T from 200 K to thousands of kelvin, and 3 Newton iterations from the predicted state cannot converge on Arrhenius terms, however exact the Jacobian.
- This fits the probe: a forced first step of 1e-18 s let 22 of the 23 through.

**Next (an engineering change; checked, not assumed):** allow more corrector failures (`CVodeSetMaxConvFails`), so that h keeps falling by 4 until Newton converges (from 1e-6 s, 30 failures reach 1e-18 s). This changes nothing on any path that never reached 10 failures in one step: CVODES's path is identical up to that point. So every branch and run that passed before is reproduced bit for bit, which the reaction tests and C1 will show. Before that change, the first step and the failure count are measured on one failing start.

**The mechanism, measured (18:35 CDT; a probe build that reports CVODES's counters on failure, not kept).** All 23 starts have a first step of 3.3e-7 to 2.4e-6 s, against T/|dT/dt| of 4.9e-14 to 2.2e-10 s. Each starts the same way: exactly 10 corrector failures with 10 Jacobian evaluations (a fresh exact Jacobian at every retry), 12 to 22 Newton iterations, no step taken. This is the inferred mechanism, confirmed.

**The change: 40 corrector failures allowed per step, not 10** (`CVodeSetMaxConvFails` in `ReactionStep`; with the analytic Jacobian, the code that is kept).
- **Check 2 rerun after the change (flagged as after the result):** 0 integrator failures in all three candidates, so gate 2(c) passes (`table_b/small_highT_maxconvfails_2026-10-20.txt`, nodes `..._maxconvfails_nodes_2026-10-20.csv`). The 6,000 K upper branch still integrates.
- **The gate table after the change:**

| Gate | Y_H2O | Y_H2O + Y_HO2 + Y_H2O2 | The O2 consumed |
|---|---|---|---|
| 2(a) worst turn-back (bar 1e-9) | 5.78e-2 | 1.21e-1 | 7.04 |
| 2(a) N-gate failures (bar 0) | 27 | 25 | 25 |
| 2(b) worst \|ΔT\| (bar 3 K) | 65.4 K | 65.4 K | 132.4 K |
| 2(c) integrator failures (bar 0) | **0** | **0** | **0** |
| 2(c) frozen branches, lower and upper | 25, 20 | 25, 20 | 15, 10 |
| Criterion 1 after projection (bar 1e-12) | 2.2e-16, 2.4e-16 | 4.4e-16, 2.7e-16 | 4.4e-16, 3.4e-16 |

- **The 23 branches now integrated add to 2(a), and nothing else changed.** Y_H2O's N-gate failures rise from 25 to 27. The second candidate's worst turn-back rises to 0.121, at f_H 0.056, Z_N 0.3, η 0.2 and ρ 30, which is one of those branches. 2(b) and the frozen counts are unchanged.
- **Branch time:** mean 9.6 ms and worst 89 ms, against 4.7 and 20 ms, over 316 branches against 270. The new ones are the cold dense-radical starts.
- **Gates 2(a) and 2(b) still fail.** Their failures are at f_H 0.556 and f_st and at f_H 0.056, all outside O/F 4 to 6. That is item 2's question.

**Check 3: pass** (measured). `ctest -R "reaction|pasr"` passes with the analytic Jacobian and 40 corrector failures (`reaction_verification` 2 s, `pasr_verification` 391 s).
- **The 40 failures change nothing there:** a `crucible_pasr_tests` built with the analytic Jacobian and the default 10 prints the same output to every digit, timing lines aside.
- **The Jacobian does move the PaSR test's integration errors, within its limits:** against a build that predates it, for example 1.441e-8 to 2.318e-8 in T at rtol 1e-10, and 2.255e-6 to 1.461e-6 at rtol 1e-7. A different Newton matrix takes a different path inside the tolerance. That is expected, not a regression, and every case still passes its stated limit.
- `reaction_verification`'s ignition and gri30 outputs are unchanged apart from wall times.

**Check 4: pass** (measured, 18:47 CDT; histories and logs in `/tmp/c1jac`, not kept). Both builds ran to 4 ms with the 300 J igniter. Against the difference-quotient build of the parent commit, the analytic build with 40 corrector failures gives:
- c* 2,477.28 m/s, vacuum thrust 1,665.07 N and vacuum Isp 424.53 s in both, identical to every printed digit, so within 1e-4 relative.
- At every printed time, the injector pressure and the outlet mass flow differ by 0, and the vacuum thrust by at most 1.34e-7 relative. The injector-pressure bar was 1e-3.
- 34,701 flow steps in both.
- The light-off time is not printed in finite-rate mode, so that part of the check was not judged. The identical pressure history through ignition is the evidence that it did not move.
- **Reaction wall time:** 923.4 s against 778.2 s (+19%, or 138.6 against 116.8 µs per cell update; measured, on a shared machine at load 17 to 23, so the comparison is loose). The study does not print CVODES's own step count. On these states the exact Jacobian buys no speed. Its use is the 23 cold starts above, through the corrector failures it shares with the quotient.

## The small table in O/F 4 to 6 and a chamber pressure band (20 October 2026, 18:43 CDT)

**The table:** f_H nodes at O/F 6, 5 and 4 (1/7, 1/6, 1/5; the builder's new `CRUCIBLE_TABLE_B_FH_NODES`); Z_N {0, 0.3}; η 6 nodes; ln ρ 5 nodes over [0.3, 80] kg/m³; on h2o2_highT. Command: `CRUCIBLE_TABLE_B_FH_NODES=0.142857,0.166667,0.2 crucible_build_manifold_table <out> 0 0 0 0,0.3 6 5 0.3 80 2 ../tools/mechanisms/h2o2_highT.yaml`.

**The density band (declared; no RL10 measurement read).** The RL10's chamber pressure is one of its blind predictions, so the band is a wide declared one: 1 to 10 MPa, which contains C1 (injector pressure 1.55 MPa at 0.5 ms and 3.3 MPa at 2.5 ms, measured this session) and any engine of the RL10's class. In density (*derived*, ideal gas at O/F 5, W about 12 kg/kmol): hot products at 1 MPa and 3,600 K are 0.37 kg/m³, and cold products at 10 MPa and 200 K are 73. The table is rectangular in (η, ln ρ), so its dense hot corner lies far above the band.

**The gate table** (measured; `table_b/of4to6_highT_2026-10-20.txt`, nodes `..._nodes_2026-10-20.csv`). The first candidate, Y_H2O, passes 2(a), so the builder did not try the others.

| Gate | Y_H2O |
|---|---|
| 2(a) worst turn-back (bar 1e-9) | **0** (pass) |
| 2(a) N-gate failures (bar 0) | **0** (pass) |
| 2(b) worst \|ΔT\| (bar 3 K) | **0.010 K** (pass; f_H 0.2, Z_N 0, η 0.6, ρ 0.3) |
| 2(c) integrator failures (bar 0) | **2** (FAIL) |
| 2(c) frozen branches, lower and upper | 15, 20 |
| Criterion 1 after projection (bar 1e-12) | 4.4e-16, 3.6e-16 (pass) |
| Nodes: degenerate, zero-span, outside the thermo range | 0, 20, 0 |

- **The 2 failures are not integrator failures but the builder's named refusal** ("no frozen start up to 6000 K carries the node's energy"). Both are at η 0.6 and ρ 80 kg/m³ (f_H 1/7 and 1/6, Z_N 0).
  - At that density the 6,000 K equilibrium is so little dissociated that no composition frozen from it holds the node's energy at 200 K.
  - The node's equilibrium pressure is of order 200 MPa (*inferred*: about 3,500 K at W 12), 20 times the band's top. So the failures lie outside the band, in the table's rectangular corner.
- **Inside O/F 4 to 6, where C1 and the RL10 run, the progress variable Y_H2O is monotone and 2(b) holds to 0.01 K.** The 2(a) and 2(b) failures of the earlier small table were at f_H 0.056, f_st and 0.556 (O/F 16.9, 7.9 and 0.8), outside this range.

## The fallback rate in a C1-like run, stated before it is measured (20 October 2026, 18:50 CDT)

The Director's rule: a declared domain with a CVODES fallback is acceptable if, in a C1-like run, the fallback takes under about 1% of cell updates; otherwise the progress variable is reopened.

- **What is counted.** The runtime is not yet written, so the count is the cell updates that the *domain* sends to CVODES. A cell update is counted if its state is outside the table's axes (f_H, Z_N, η, ρ), or if any of its 16 corners has failed a gate: a 2(a) turn-back above 1e-9 or an N-gate failure, a 2(b) |ΔT| above 3 K, a 2(c) failure, or the outside flag. The runtime's other causes need the runtime and are not counted: off the curve, ΔY above ΔY_max, frozen reach and AR. So the measured rate is a **lower bound** on the fallback rate.
- **The table:** as the O/F 4 to 6 table above, but with ln ρ over [1e-3, 80] kg/m³ in 9 nodes, so that the axes reach C1's nozzle (its exit density is about 0.03 kg/m³, *inferred*).
- **The run:** C1 finite-rate 32x6 to 4 ms with the 300 J igniter, on h2o2_highT (the table's mechanism, through `CRUCIBLE_MECHANISM`), with the analytic Jacobian. Every 2 µs, each cell's (f_H, Z_N, η, ρ) is written through the chamber study's new `CRUCIBLE_STATE_SAMPLES`, with η taken from the cell's own complete-combustion products as in the builder. Each sample is weighted by the flow steps since the previous one, so the fraction is of cell updates.
- **Reported:** the fraction over the whole run (the stated figure, judged against 1%), and over the last millisecond (the burning chamber); each broken down by cause. The causes are Z_N above the axis (the ambient N2 the start purges), ρ, η or f_H outside the axes, and a failed corner.
- **The decision follows the rule.** Under 1%: a declared domain, named by its axes. Over 1%: the cause breakdown says whether a wider axis would fix it or the progress variable must be reopened.

### The table it is measured against (measured, 18:46 CDT; `table_b/c1range_highT_2026-10-20.txt`, nodes `..._nodes_2026-10-20.csv`)

`CRUCIBLE_TABLE_B_FH_NODES=0.142857,0.166667,0.2 crucible_build_manifold_table <out> 0 0 0 0,0.3 6 9 1e-3 80 2 ../tools/mechanisms/h2o2_highT.yaml`: 324 nodes, 536 branches (mean 8.7 ms, worst 82 ms).

| Gate | Y_H2O | Y_H2O + Y_HO2 + Y_H2O2 | The O2 consumed |
|---|---|---|---|
| 2(a) worst turn-back (bar 1e-9) | 0.253 (f_H 1/6, Z_N 0, η 0.6, ρ 1e-3) | 0.198 (same node) | 4.76 |
| 2(a) N-gate failures (bar 0) | 8 | 8 | 8 |
| 2(b) worst \|ΔT\| (bar 3 K) | 0.126 K | 0.126 K | 1.05 K |
| 2(c) integrator failures (bar 0) | 2 | 2 | 2 |
| 2(c) frozen branches, lower and upper | 27, 28 | 27, 28 | 9, 0 |
| Criterion 1 after projection (bar 1e-12) | 4.4e-16, 3.2e-16 | 4.4e-16, 3.3e-16 | 4.4e-16, 3.1e-16 |

- **2(a) fails at low density, in every candidate.** For Y_H2O, 85 of 324 nodes fail. 83 of them are at ρ ≤ 0.069 kg/m³ (all four of the lowest density nodes: 25, 25, 23 and 10 failures), and they are spread over η 0.2 to 1. The other 2 are the ρ 80 refusals.
  - *Inferred:* at 1e-3 to 0.07 kg/m³ and these energies, the equilibrium is so dissociated that the path from the frozen start makes H2O and then gives some back as the radical pool grows. A product mass fraction then turns back.
  - From 0.28 to 19.5 kg/m³ (the band's interior) there is no 2(a) failure, as in the 0.3 to 80 table.
- **2(b) passes in all candidates** (worst 1.05 K). **2(c)'s 2 failures** are the same two ρ 80 refusals as in the O/F 4 to 6 table.
- **So no progress variable passes on the C1-range axes,** and the fallback rate is measured against Y_H2O, the declared one (candidate 0). A cell is counted as fallback if any of its 16 corners fails a gate.
- **Deviation from the statement above (flagged).** The C1 run on h2o2_highT failed at 0.206 ms, as the igniter turned on: "Temperature iteration did not converge", in the flow's energy-to-temperature inversion (`/tmp/c1fb/fb_highT_failed.log`). C1 on h2o2.yaml passes the same point. That failure is an open finding, not fixed this session. **Its mechanism (measured, 19:00 CDT, in a probe build that prints the failing inversion, not kept):**
  - Six µs into the igniter, one cell holds e = 79.5 MJ/kg at 12,176 K the step before, almost fully dissociated (Y_O 0.81, Y_H 0.16, Y_OH 0.024).
  - h2o2_highT's 1,000 to 6,000 K fits, extrapolated past their range, turn cv negative: from 10,630 K for OH, 11,340 K for O and 12,140 K for H2O. O's u then peaks at 23.6 MJ/kg. The mixture's e(T) never reaches the cell's energy, and Newton walks off to 1e33 K.
  - GRI h2o2.yaml's fits (to 3,500 K) extrapolate monotonically to 40,000 K, so the same run passes on them. Its temperatures there are equally outside its fits.
  - So it is the igniter driving cells to about 12,000 K, beyond any fit in the tree.
  - **The question for the fix (engineering, open):** the thermo above 6,000 K (for example TM-4513's 6,000 to 20,000 K range, or a constant cp beyond the top of the fit), or a gentler igniter deposit.
  - **So the sampled run is on h2o2.yaml.**
  - f_H, Z_N and ρ do not depend on the mechanism file, but η does. Its top end is u(6,000 K) of the products, and h2o2.yaml's polynomials end at 3,500 K. Extrapolated, its u(6,000 K) of the O/F 4 to 6 products is 8.0 to 8.6% higher than h2o2_highT's, so the run's η reads up to 0.079 low near η 1. At 200 K the two energies agree. (Measured, in a probe.)
  - So each sample's η is re-expressed in h2o2_highT's energies before it is classified. The cell's energy e is recovered from η on h2o2.yaml's ends, and then η = (e − u_b(200 K)) / (u_b(6,000 K) − u_b(200 K)), using the products from (f_H, Z_N) (new `tools/state_samples_eta.cpp`, built as `crucible_state_samples_eta`).
  - What is left: the C1 states are h2o2.yaml's, whose thermo below 3,500 K differs from h2o2_highT's by the OH and 1-bar corrections (about 1e-4 in C1's c*, TABLE_A.md).

### The fallback rate: measured (20 October 2026, 19:08 CDT; `table_b/c1_fallback_2026-10-20.txt`)

The run: C1 finite-rate 32x6 to 4 ms, 300 J igniter, h2o2.yaml, analytic Jacobian. Its end values match the analytic C1 of check 4 to every printed digit (c* 2,477.28 m/s), so the sampling does not disturb the run. 2,001 samples of 192 cells, 6,662,784 cell updates. η re-expressed in h2o2_highT's energies, classified against the C1-range table with Y_H2O.

| | Whole run (the stated figure) | Last ms |
|---|---|---|
| **Fallback, share of cell updates** | **47.97%** | **41.21%** |
| A corner failed a gate | 41.86% | 41.21% |
| η above the axis (the igniter's cells, η up to about 15) | 4.57% | 0 |
| Z_N above the axis (the ambient N2 the start purges) | 1.45% | 0 |
| No H or O (pure N2) | 0.09% | 0 |

With η left on h2o2.yaml's scale the whole run is 47.87%, so the conversion moves it by 0.1 points.

**Where it is (last ms, measured):**
- Every cell of the divergent nozzle from column 19 to the exit (throat at 16), and 19% of column 18. These are the cells with ρ below 0.28 kg/m³: 41.2% of the run's cell updates.
- Every one of their failing corners fails 2(a)'s turn-back (Y_H2O turning back on the low-density branches, the 83 nodes at ρ ≤ 0.069). None fails 2(b), 2(c) or the outside flag.
- The chamber and the convergent section (ρ 0.28 to 4.8 kg/m³, 58.8% of updates) are on the table at every update.

**The decision, by the Director's rule: over 1%, so the progress variable is reopened.**
- A wider axis would not help: the fallback is failed corners inside the axes, not cells outside them.
- A declared domain that leaves out the nozzle would hand CVODES 41% of the steady run, about 40 times the rule's bar.
- The measured rate is a lower bound (the runtime's own causes are not counted), so the conclusion does not depend on what is left out.

**A candidate for the reopened variable (a design choice, for the Director; *derived*, not tested):** the specific entropy s.
- Every branch the builder integrates is an adiabatic reactor at fixed volume, so u and ρ are held, and dS/dt = Σ_r (r_f − r_r) R ln(r_f / r_r) / ρ ≥ 0.
- That holds term by term when every reaction is reversible with its reverse rate from the equilibrium constant. All 29 reactions in h2o2_highT and h2o2.yaml are written `<=>` with no explicit reverse rate (checked).
- So s is non-decreasing on every branch by the second law, and 2(a)'s turn-back is zero to round-off by construction, at every density.
- What it would still have to pass, by measurement:
  - the N gate and the resolution near equilibrium, where ds/dt → 0;
  - the runtime cost of s per cell (it is not linear in Y: one log per species);
  - its span on branches that barely react.

## The specific entropy as progress variable: its checks, stated before its code (20 October 2026, 19:25 CDT)

The Director approved s, on the second-law argument above. These checks are stated before any code, which comes next session. Nothing below is measured yet.

- **The variable.** c = s, the mixture's specific entropy, s = Σ Y_k (s_k°(T) − R_k ln(X_k p / p°)).
  - At build it comes from Cantera along each branch.
  - At run time it comes from the cell's own T, ρ and Y on the same NASA-7 data. The 7th coefficient is already in `Species`. Above 6,000 K it uses the constant-cp continuation of the section below, s° continued as cp ln(T / 6,000 K).
- **A consequence of the argument, found while stating it (derived): both branches now approach equilibrium from below.**
  - At fixed (u, ρ, elements), equilibrium is the entropy maximum. So the upper branch, from the complete-combustion products, also rises in s toward s_eq. It does not fall toward it, as Y_H2O does.
  - The runtime's coordinate (lower branch 0 to 1, upper 1 to 2, "The runtime" below) therefore cannot tell the branches apart by c alone.
  - **The rule, stated now (engineering):** a cell takes the branch whose blended curve, at the cell's s, is nearer the cell's composition, by the largest |ΔY_k| (the snap's own measure). Ties and both-above-ΔY_max cases go to CVODES, counted.
  - **Check E0:** on the C1 samples, the share of cell updates that the rule sends to CVODES. It is reported and counted in the fallback rate.
- **E1. 2(a)'s turn-back: zero by construction, measured anyway.** The bar stays as stated: no decrease above 1e-9 of the branch's span, until within 1e-4 of the span from s_eq, on every branch of three tables: the small one, the O/F 4 to 6 one, and the C1-range one (whose 83 low-density nodes fail for every linear candidate). The N gate is measured on the same branches.
  - A turn-back above the bar would mean an irreversible reaction or a thermo inconsistency between the builder's s and the integrator's. It would be reported, not tuned.
- **E2. Resolution near equilibrium: the c nodes are clustered toward s_eq.** The entropy production is quadratic in the distance from equilibrium, so s_eq − s ∝ |Y − Y_eq|² there (derived). Nodes uniform in s would leave the last approach with almost no points.
  - **The node rule:** nodes are uniform in the distance, at s_eq − s = (s_eq − s_0)(1 − ξ)², with ξ uniform over the existing 30 lower-branch and 9 upper-branch slots. The early geometric nodes in Y_H2O for the induction are kept, mapped through s.
  - λ = (ds/dt) / (s_eq − s) stays finite at equilibrium. With a linear relaxation of rate k it tends to 2k, so the exponential update of the stated step still lands on s_eq exactly.
  - **The check:** at the midpoint between every pair of c nodes on every branch, the table's T and Y (interpolated linearly in s, as the lookup does) are compared with the trajectory's own state at the same s. For that, the builder integrates to the midpoints as extra output points.
  - **Bars:** |ΔT| ≤ 3 K (2(b)'s bar and reason) and max_k |ΔY_k| ≤ ΔY_max = 0.01 (the snap's bar). Both are reported separately for the last quarter of the span, where s flattens. A fail there adds nodes toward s_eq. It does not relax a bar.
- **E3. Span on weakly reacting branches.** On frozen or barely reacting branches, s_eq − s_0 can sit near round-off. Every branch's span is reported against 1e-9 |s| (the turn-back bar's scale). A branch under it is treated as frozen, as now, and counted.
- **E4. Cost per cell.** s needs one logarithm per species plus ln T at the cell's state, once per cell per reaction call, since the lookup inverts the blended curve at the cell's s.
  - It is timed in isolation (ns per cell, `h2o2_highT`, 10 species), and as its share of the C1 32x6 step at the same threads.
  - **Bar (guessed): under 5% of the flow step.** Over it is reported against criterion 6's 5 minutes, not tuned.
- **E5. The fallback rate, rerun.** The same samples (`/tmp/c1fb/states_highT.csv`, h2o2.yaml's run re-expressed), then the samples of the h2o2_highT run in the next section, both classified against the C1-range table with s, by `tools/fallback_rate.py`. The Director's 1% bar applies, with E0's rejects counted.

## The entropy candidate in the builder and the E0 replay, stated before their code (20 October 2026, about 22:05 CDT)

What E0 to E5 leave to the code, settled before it is written. No check or bar above changes.

- **c and its rate.** c = s from Cantera's `entropy_mass` at the step's (T, ρ, Y), and ds/dt = −Σ_k μ_k (dY_k/dt) / (W_k T), with μ_k from `getChemPotentials`: the Gibbs relation at fixed u and ρ. This is the integrator's own thermo (the step integrates T on Cantera's), so a turn-back cannot come from two thermo sets. It is candidate 3; `CRUCIBLE_TABLE_B_CANDIDATES` (a comma list, default 0,1,2,3) selects which candidates are tried, so s can be built where Y_H2O already passes.
- **The stop rule changes for s (*derived*).** "|c − c_eq| ≤ 1e-6 of the span" would stop an s-branch at about 1e-3 of the composition distance, since s_eq − s is quadratic in it, and 2(b) would then read the stop and not the table. Finer in s is below what s resolves (s is about 1e4 J/(kg K), the integration's rtol is 1e-10). So an s-branch stops when max_k |Y_k − Y_k,eq| ≤ 1e-6 of max_k |Y_k,start − Y_k,eq|. E1's turn-back window stays in s, as stated.
- **Zero span for s:** |s_eq − s_start| < 1e-9 |s_eq| (E3's scale). Such a branch is stored as a zero-span branch is now, and counted under E3.
- **The stored points.** Lower point p: s = s_eq − (s_eq − s_0)(1 − ξ_p)², with ξ_p the existing fractions (0, 12 geometric from 1e-9 to 0.1, 18 uniform to 1). So the induction stays geometric: s − s_0 ≈ 2ξ (s_eq − s_0), from 2e-9 to 0.19 of the span. This is my reading of E2's "kept, mapped through s" (engineering). Upper point p, from 0 at equilibrium to 8 at complete combustion: s = s_eq − (s_eq − s_cp)(p/8)². A frozen branch uses its reach in place of s_eq.
- **The file.** Each point now also stores its c (layout version 2: Y, λ, c), and the table stores the mechanism's reference pressure, so that the runtime's s uses the builder's p°. Version 1 files are refused (none is in use).
- **E2's reference.** For each pair of neighbouring stored points on a branch, the midpoint s_m = (s_p + s_p+1)/2. The trajectory's own state at s_m is found by integrating from the last internal step below it (`ReactionStep::advance`, the builder's rtol and atol), with an Illinois iteration on the time, until |s − s_m| ≤ 1e-6 of the pair's s gap. The table's state there is the stored Y interpolated linearly in s, with T from the node's e on that Y. Per node: the worst |ΔT| and max_k |ΔY_k|, over all midpoints and over those with s_eq − s ≤ span / 4 (the last quarter in s).
- **E0's replay** (`tools/table_b_branch_rule.cpp`, new), on samples that now carry each cell's T and Y (`CRUCIBLE_STATE_SAMPLES`). Per sample:
  - s_cell from the cell's T, ρ and Y on the mechanism's NASA-7 data, with the table's p° (the runtime's evaluation, also E4's timing);
  - the 16 corners blended point by point, as "The runtime" states, separately for each branch;
  - on each branch, the cell's s located on the blended curve by segment search, and Y there;
  - the nearer branch taken by max_k |ΔY_k|; CVODES, counted by cause, if the cell's s is on neither curve ("off the curve", including above the blended s_eq), if both ΔY exceed ΔY_max = 0.01, or on an exact tie.
- **A prediction for E0, stated before it is measured (*derived*, the size *inferred*).** The blended curve's s is a weighted mean of the corners' s, and s is curved in η and ln ρ: at fixed composition ∂²s/∂e² = −1/(T² c_v). At the C1-range table's spacing (Δη 0.2, about 4.5 MJ/kg; Δ ln ρ 1.4) the blend's s_eq can sit tens of J/(kg K) from the cell's own. Near equilibrium s_eq − s is quadratic in the distance, so that offset is a large share of the composition distance. So E0 may reject many near-equilibrium cells (s above the blended s_eq, or ΔY above 0.01), in the chamber as well as the nozzle. If it does, the cause is measured first: the blend's s_eq against each sampled cell's own s_eq from `equilibrateUV`. The fix is not chosen before that.

### E2 measured on the small table: linear in s fails, by a square root (20 October 2026, 22:15 CDT)

The small table with s (`table_b/small_s_2026-10-20.txt`): E2 fails at 142 of 180 non-degenerate nodes. Worst |ΔT| 72 K (f_H 0.112, η 1, ρ 1e-3), median 26 K. Every worst case is in the last quarter of the span.
- **The mechanism (*derived*, matches the measurement):** near equilibrium s_eq − s ∝ d², with d the composition distance. So Y is linear in d = √(s_eq − s), and has an infinite slope in s at s_eq.
  - Linear interpolation in s on the last interval (from d_1 to 0) puts the s-midpoint's Y at d_1/2, where the trajectory is at d_1/√2. That is an error of 0.207 d_1.
  - On the upper branch d_1 is 1/8 of the span. At η 1 the T span is 6,000 − 2,666 K, so the expected error is about 86 K, against 72 K measured.
  - Adding nodes cannot fix this efficiently: about 230 uniform nodes per branch would be needed to reach 3 K.
- **The alternative, stated before it is measured (engineering, bars unchanged):** the lookup interpolates the stored Y linearly in d = √(s_ref − s), with s_ref the branch's last point (s_eq, or a frozen branch's reach). Since every corner's points sit at s_eq − span (1 − ξ)², the blend's points keep the common ξ, and the cell's ξ comes from one square root.
  - Checked at the d-midpoints, where interpolation in d errs most, with the same bars (3 K, ΔY 0.01, all midpoints and the last quarter). The builder now reports both.
- **2(a) on the small table: one node over the bar, by round-off (measured):** a fall of 1.1e-11 J/(kg K) (1.7 ulp of s = 28,769 J/(kg K)) on an upper branch of span 0.007 J/(kg K), at f_H 0.556 (O/F 0.8), Z_N 0.3, η 0.2, ρ 30. That is not a second-law turn-back. The 1e-9 bar is below s's round-off for spans under about 0.1 J/(kg K). It is reported, not re-barred; the node is outside O/F 4 to 6.
- The N gate fails at 26 nodes (f_H 0.556: 23; 0.056: 3; mostly η ≥ 0.6), outside O/F 4 to 6, as with Y_H2O. 2(b) passes (0.015 K), 2(c) has 0 integrator failures, and E3 finds 30 nodes with a branch under 1e-9 |s_eq| (stored as zero span).

### E2 in the distance measured, the mechanism, and per-node refinement, stated before its code (20 October 2026, about 22:16 CDT)

**Measured.** On the small table, interpolating in d = √(s_ref − s) at the d-midpoints gives a worst |ΔT| of 25.9 K (f_H 0.056, η 1, ρ 0.013) and a worst ΔY of 0.030. 117 of 180 nodes fail, against 142 when linear in s. Median 4.1 K, 90% 10.6 K.

A per-pair dump (`CRUCIBLE_TABLE_B_E2_DUMP`) on f_H 0.056, Z_N 0, 6 η by 3 ρ:
- **The stored points are right.** Each matches the trajectory at its own s to 0.010 K at worst. The error is the curve bending between points, not the points themselves.
- **(i) The lower branch's last two or three gaps at high η.** At η 1, ρ 1e-3, T goes 2,523, 2,509, 2,573, 2,677 K over ξ 0.94 to 1. *Inferred:* the shuffle reactions equilibrate first, then three-body recombination raises T by about 170 K while s gains under 1% of its span. This is the slow approach that a nozzle samples. The worst gap is 19.4 K at the d-midpoint.
- **(ii) The upper branch's 8 gaps,** each about 400 K of T at η 1: 2 to 7 K.

**Why uniform refinement is a poor fix (*inferred*).** Linear interpolation errs as the square of the gap. A uniform fix sized to the worst node needs about 2.5 times the lower uniform points everywhere, for a bend whose place moves with (η, ρ).

**The change (engineering, bars unchanged).**
- **Per-node refinement.** After the 40 base points, each gap on a branch whose d-midpoint fails 3 K or ΔY 0.01 gets a new point: the trajectory's own state at that d-midpoint. That state comes from CVODES with the Illinois iteration, projected onto the node's elements, with λ from the rates there. This repeats until every gap passes, or the node holds 400 points; a node at the cap fails E2 and is reported.
- **The independent check.** The midpoints are now the refinement criterion, so E2 is re-measured where it was not used: at the quarter points of every final gap (d at 1/4 and 3/4), against the same bars, over every gap and over the last quarter of the span. The report also gives the points added (median and maximum per node) and the nodes at the cap.
- **Layout version 3.** Each node carries its first point, its point count and its lower branch's count. The base 40 points keep the common ξ.
  - The runtime takes the cell's ξ from one square root on the corners' blended s_eq and spans: ξ = 1 − √((s̄_eq − s_cell)/s̄pan). Every stored point satisfies s = s_eq − span (1 − ξ)² on its own corner, so this ξ is exact for the blend.
  - Each corner's segment then comes from its own search in d: 16 searches of up to 400 points, about 9 comparisons each (*inferred*: under 1 µs against the 155 µs cell update).
- **The N gate runs on the final points.**

### Midpoint refinement measured; the quarter points fail by an off-centre bend; refinement on three interior points, stated before its code (20 October 2026, 22:20 CDT)

**Measured** on the small table (`table_b/small_s_refined_2026-10-20.txt`):
- Refinement adds a median of 1 point per node and at most 11. No node reaches the cap. The final midpoints are at worst 2.996 K and ΔY 8.7e-3, under the bars by construction.
- The independent quarter points fail. 8 of 12,064 checks are over a bar: worst |ΔT| 4.32 K (f_H 0.556, Z_N 0.3, η 0.6, ρ 1e-3) and worst ΔY 0.0124 (f_H 0.112, Z_N 0.3, η 1, ρ 1e-3). All 8 are on the lower branch.
- **The mechanism (measured from the per-gap dump):** each failing gap is lopsided. One quarter point is 3.0 to 4.3 K off and the other 0.4 to 1.5 K.
  - The deviation from the chord is not the symmetric parabola a midpoint test assumes. The curvature is concentrated near one end of the gap: the turn into the recombination approach, ignition, or the last gap into s_eq.
  - So the midpoint underestimates the gap's worst error, by up to 1.5 times here.

This is the second fix of this bug to fail, so the mechanism is written down above before any further code.

**The change (engineering, bars unchanged).** A gap is refined if any of its three interior points (d at 1/4, 1/2 and 3/4) fails 3 K or ΔY 0.01. The new point goes at the worst of the three, with the trajectory's own state there. The independent check moves to the eighth points (1/8, 3/8, 5/8, 7/8) of every final gap. The cap and the report are otherwise unchanged.

### Three-point refinement measured; the eighth points exceed by sampling; a refinement target under the bar, stated before its code (20 October 2026, 22:24 CDT)

**Measured** on the small table (`table_b/small_s_refined3_2026-10-20.txt`):
- Refinement adds a median of 1 point per node and at most 10. None reaches the cap. The final 1/4, 1/2 and 3/4 points are at worst 2.994 K and ΔY 8.8e-3.
- The independent eighth points: 6 of 24,100 checks are over a bar, at worst 3.19 K. ΔY is at worst 9.6e-3, under its bar. 4 are on the upper branch and 2 on the lower.
- **The mechanism (measured from the dump):** in each of the 6 gaps, the three samples sat just under 3 K. The error profile across the gap is smooth (for example 0.67, 1.98, 3.19, 1.48 K at the eighths) and peaks between two samples, up to 6.7% above the highest sample. This is not a new bend. A refinement that stops exactly at the bar will always be exceeded by a finer check of its own gaps.

**The change (engineering; E2's bars unchanged at 3 K and 0.01).**
- The refinement's target is 2/3 of the bars: 2 K and ΔY 0.0067. The eighth-point check stays at the full bars.
- The margin covers a smooth profile's peak between samples (measured at up to 6.7% here). 2/3 leaves room for sharper profiles than these; that factor is *guessed*.

### Refinement to 2/3 measured; capped nodes from a time coordinate that stands still; the search seeded from the rate, stated before its code (20 October 2026, about 22:24 CDT)

**Measured.** On the small table, the eighth points pass: worst 2.17 K and ΔY 6.3e-3. Refinement adds a median of 4 points per node and at most 23, and no node reaches the cap.

The O/F 4 to 6 table (0.3 to 80 kg/m³) and the C1-range table (1e-3 to 80) both fail. They have 10 and 13 nodes at the cap, each with 360 points added, and the eighth points reach 8.5 and 9.9 K.

The capped nodes are f_H 1/7 at η 0.8 and f_H 0.2 at η 0.6, both Z_N 0, at every ρ. On one of them (f_H 0.2, η 0.6, ρ 80), the refinement inserted the same state 341 times, at a c below its left neighbour's.

**The mechanism, measured.** The E2 search returned a state far from its target without saying so. It now throws if it does not converge, and the throw names its numbers.
- The lower branch starts at 467 K and induces for 8.3e7 s (about 2.6 years) before it ignites. Ignition needs internal steps h below the double-precision spacing of t, which is about 1.5e-8 s there. So t + h == t.
- CVODES advances the state while the recorded time stands still: 2,105 of 10,895 raw steps share a time with the next.
- The system is autonomous at fixed (u, ρ), so the state path, and hence the stored points, is unaffected. 2(b) passes at 0.011 K.
- The search, however, integrates from a raw step over the recorded time to the next. That interval is 0, so the bracket collapses.

**The fix.**
- The search keeps integrating locally from the raw step's state, starting at t = 0.
- Where the recorded interval is not positive, the first trial step comes from the rate there: (s_target − s_a)/(ds/dt), with ds/dt = λ (s_eq − s_a). The bracket then grows by doubling, up to 100 times.
- A search that still fails throws, and the node is counted as an E2 failure; it is never stored.

### The search seeded from the rate measured; the remaining failures sit at s's integration noise floor; a counted fallback to the raw steps, stated before its code (20 October 2026, about 22:25 CDT)

**Measured.** The nodes at the cap are gone from both tables: none at the cap, at most 15 points added. The eighth points are at worst 1.98 K and ΔY 3.4e-3 (O/F 4 to 6), and 2.09 K and ΔY 7.1e-3 (C1 range).

Two E2 searches fail in each table, with CVODES flag −3, at f_H 1/6 and 0.2, Z_N 0.3, η 0.2, ρ 80.
- **The mechanism (measured on f_H 0.2):** the node's upper branch is frozen, with a span of 0.0029 J/(kg K) (2e-7 of s). Its last raw steps sit at t of about 1.6e14 s, with λ of 1.7e-16 /s, and differ in s by less than CVODES's reproducibility of s there (rtol 1e-10 of Y).
- So re-integrating from one raw step over the recorded interval falls short of the next raw step's s. The bracket then doubles until the integration fails at 6.2e14 s.
- The T change along the whole of such a branch is a small fraction of a kelvin (*inferred* from its span), so E2 there tests the noise, not the curve.

**The fix.** If a search fails (no convergence, or a CVODES failure), the trajectory's state at the target comes from the two bracketing raw steps, interpolated linearly in s. Those are CVODES internal steps, the source of the stored points themselves.
- Each such fallback is counted and reported: the number, the nodes, and the largest span among their branches. It is not an E2 pass on its own: the report shows where E2 was not independent.

### The raw-step fallback measured: E2 passes on both full tables; what still fails as written (20 October 2026, 22:26 CDT)

**Measured** (`table_b/of46_s_2026-10-20.txt`, `table_b/c1r_s_2026-10-20.txt`, nodes `..._nodes_2026-10-20.csv`):

| Check | O/F 4 to 6 (0.3 to 80 kg/m³, 180 nodes) | C1 range (1e-3 to 80 kg/m³, 324 nodes) |
|---|---|---|
| E2 at the eighth points (bars 3 K, ΔY 0.01) | 1.98 K, 3.35e-3: pass | 2.09 K, 7.08e-3: pass |
| Points added per node (median, max; at the cap) | 5, 15; 0 | 4, 15; 0 |
| Searches from the raw steps | 2, on 2 nodes (largest span 2.947e-3 J/(kg K)) | the same 2 |
| E1 turn-back (bar 1e-9 of the span) | 1.384e-9: FAIL | 1.384e-9: FAIL |
| N-gate failures (bar 0) | 0 | 8 |
| 2(b) worst \|ΔT\| (bar 3 K) | 0.011 K | 0.030 K |
| 2(c) integrator failures (bar 0) | 2 | 2 |
| Criterion 1 after projection (bar 1e-12) | 2.8e-16 | 3.6e-16 |
| E3: nodes with a branch under 1e-9 \|s_eq\| | 30 | 54 |

The gates as written fail, on three counts. Each is named here; none is tuned.
- **The turn-back, 1.384e-9 of the span.** It is at one node in both tables: f_H 1/6, Z_N 0.3, η 0.2, ρ 80. Its upper branch is frozen, with a span of 2.63e-3 J/(kg K), about 2e-7 of s.
  - So the turn-back is 3.6e-12 J/(kg K), about 2 ulp of s there (*derived*). On a span this small, 1e-9 of the span is under the resolution of s in double precision. E3's zero-span scale, 1e-9 \|s\|, does not catch it.
  - T_eq is 1,738 K at 80 kg/m³. That is far above 10 MPa (*inferred*), outside the band.
- **2(c): the two named refusals** at η 0.6, ρ 80, Z_N 0 (f_H 1/7 and 1/6): no frozen start up to 6,000 K carries the node's energy. These are the two documented under the earlier O/F 4 to 6 table, at about 200 MPa, outside the band.
- **The N gate, C1 range only: 8 nodes, all at η 1 and ρ ≤ 0.069 kg/m³.** They are f_H 1/6 and 0.2, at T_eq 2,760 to 3,451 K. η 1 is the products' energy frozen at 6,000 K, the top of the axis. The same 8 failed with Y_H2O (the table above at 18:46), since N along a branch does not depend on the progress variable. Whether C1's cells reach them is measured by E5.

## E0, E4 and E5 measured on h2o2_highT's C1 run: the domain passes everywhere, the branch rule fails at equilibrium (20 October 2026, about 22:52 CDT)

**The run.** C1 32x6 finite-rate on h2o2_highT, the torch igniter (50 kW for 1 ms), 2 threads, to 4 ms (`/tmp/c1s/`, not in git). It ended at 4 ms after 32,976 steps in 451 s wall; the reaction took 427.6 s (67.5 µs per cell update). Drift over the last ms: outlet mass flow 1.3e-3, injector pressure 1.1e-3, vacuum thrust 1.6e-3. Its state samples carry T and Y.

**The tools.** `tools/table_b_branch_rule.cpp` (`crucible_table_b_branch_rule`) replays E0's rule as stated at 22:05, with the runtime's ξ from one square root (22:16). `tools/fallback_rate.py` now splits the chamber (columns 0 to 8), the convergent (9 to 15) and the divergent (16 on), and joins the replay's rejects.
- s_cell on the core's NASA-7 data agrees with Cantera's `entropy_mass` to 1.8e-11 J/(kg K) at worst, over 3,800 samples.

**E5 without E0, the domain alone** (`table_b/e5_torch50_c1r_s_domain_only_2026-10-20.txt`; the axes and the gates of the C1-range table with s, as before):

| Share of the region's cell updates | Whole run | Last ms |
|---|---|---|
| Chamber (columns 0 to 8) | 0.746% | **0.000%** |
| Convergent (9 to 15) | 1.058% | **0.000%** |
| Divergent (16 on) | 2.063% | **0.000%** |
| All | 1.473% (1.269% is the purged ambient N2, above the Z_N axis) | **0.000%** |

With Y_H2O the divergent nozzle was 82.3% over the last ms, through 2(a) at low density. With s, no corner of a nozzle cell fails a gate. So the domain is what entropy was for, and it now passes.

**E0, the nearer-branch rule** (`table_b/e0_torch50_c1r_s_2026-10-20.txt`, `table_b/e5_torch50_c1r_s_2026-10-20.txt`): it fails, everywhere.

| Share of the region's cell updates sent to CVODES | Whole run | Last ms |
|---|---|---|
| Chamber | 98.221% | **100.000%** |
| Convergent | 97.875% | **100.000%** |
| Divergent | 85.313% | **92.708%** |
| All | 91.691% | 96.354% |

Over the last ms, every reject is the same cause: the cell's s is above the blended s_eq on both branches, so it is on neither curve.

**The cause, measured as the prediction stated** (920 samples, every 400th inside the axes, each with its own s_eq from `equilibrateUV`):

| Last ms | The cell's own deficit, s_eq,own − s_cell (median; 95%) | The blend's shortfall, s_eq,own − s̄_eq (median; 5% to 95%) | Samples where the shortfall exceeds the deficit |
|---|---|---|---|
| Chamber (61) | under 1e-4 J/(kg K); 1e-5 | 13.05; 13.05 to 13.07 J/(kg K) | 61 of 61 |
| Convergent (60) | under 1e-4; 1e-4 | 13.72; 13.04 to 17.69 | 60 of 60 |
| Divergent (120) | 13.07; 51.1 | 61.33; 17.4 to 86.5 | 120 of 120 |

- The chamber and the convergent gas is at its own equilibrium, to under 1e-4 J/(kg K). The nozzle's gas is a median 13 J/(kg K) short of its own, which is the freezing expansion.
- The blended s_eq is always below the cell's own, by 13 J/(kg K) in the chamber and 61 in the nozzle. Over all 920 samples, the shortfall is negative only down to −1.1 J/(kg K).
- *Derived:* s_eq(e, ρ, elements) is concave, since it is the maximum of a concave function under linear constraints, the thermodynamic stability condition. So a multilinear blend of it is a chord, and the chord lies below the curve. At fixed composition, ∂²s/∂e² = −1/(T² c_v), as the prediction stated.
- So a cell at or near equilibrium has its s above the end of both blended curves. s is the maximum at equilibrium, and both branches end there from below, so nothing covers the range above s̄_eq. The prediction is confirmed in sign and size ("tens of J/(kg K)").

**E4, measured on one thread** (368,183 cells, the machine shared at a load of about 10 on 10 cores; four runs):
- s_cell costs 72 to 92 ns per cell. The whole rule costs 248 to 305 ns per cell: s, the 16 corners and both branches.
- C1's step at 2 threads costs 71.2 µs of wall per cell update, of which 3.70 µs is outside the reaction (both measured from the run's log).
- So s is 1.0% (*inferred*, if it runs on both threads like the reaction) to 2.0% (on one) of the step outside the reaction, and 0.1% of the whole step. That is under E4's 5% bar: **E4 passes.** The whole rule is 3.4% to 6.7% of the step outside the reaction.

### The cell placed by its moles per kilogram, stated before its code (20 October 2026, about 22:52 CDT)

The 19:25 branch rule was stated as engineering, and it fails E0 by the predicted mechanism. Its replacement is the design's own snap at equal N (the 02:47 amendment and the 02:59 handover), with one change: φ_N is the cell's own. The bars do not change.

**Why N, physically.** N = Σ Y_k/W_k is linear in Y. Unlike s, it is not an extremum at equilibrium: along each node's joined curve it keeps falling through equilibrium onto the upper branch (N_0 to N_eq on the lower, N_eq to N_cp on the upper; the N gate, measured at every node of the C1-range table but the 8 at η 1). So the blended curve's N covers the equilibrium point from both sides, and a cell near equilibrium lands on one branch or the other near τ = 1, with no gap. The 02:47 probe already found the snap at equal N 3.2 K off in a freezing nozzle where the snap at equal c was 245 K off. N also decides the branch, so the nearer-branch rule by max|ΔY| is dropped.

**The rule.**
- **The blended curve.** It is parametrised by τ, from 0 at the start to 1 at equilibrium (or the reach) on the lower branch, and from 1 to 2 at complete combustion on the upper branch. On the lower branch τ = ξ; on the upper, τ = 1 + u, with s = s_ref − span u² on each corner. Each corner gives its Y at its own s at that τ by the search in d, and the blend is the multilinear weights' sum. N̄(τ) is N of the blend, over the reacting species; N2 is the cell's own Z_N on both sides, so it cancels.
- **φ_N, from the cell's own N:** φ_N = (N̄(0) − N_cell)/(N̄(0) − N̄(1⁻)), with 1⁻ the lower branch's end. If N̄(0) ≤ N̄(1⁻) the curve has no fall in N, and the cell goes to CVODES ("no fall in N").
  - The weight is as at 02:59: w = clamp((φ_N − 0.005)/0.005, 0, 1).
- **w = 0 (induction):** the cell's s places it on the lower branch, with ξ from one square root as now. If s is outside the lower curve, the cell is off the curve in s.
- **w = 1:** τ_N, where N̄(τ_N) = N_cell.
  - The branch comes from N: lower if N_cell ≥ N̄(1⁻), upper if N_cell ≤ N̄(1⁺).
  - Between the two is a gap, possible only where a corner is frozen, since N̄(1⁻) = N̄(1⁺) at equilibrium. A cell there goes to CVODES ("frozen reach", the 03:10 cause). Below N̄(2) it is off the curve in N.
  - τ_N comes from an Illinois iteration on τ within the branch, to |N̄ − N_cell| ≤ 1e-12 N_cell, or 60 iterations.
- **0 < w < 1:** s* = (1 − w) s_cell + w s̄(τ_N), with s̄ the lower branch's blended s, and τ* placed from s* on the lower branch, as 02:59 states.
- **The snap's bar.** If max_k |Ȳ(τ*) − Y_cell,k| > ΔY_max = 0.01, the cell goes to CVODES ("ΔY above ΔY_max").
- **A check on the search's premise.** N̄ must not rise along the branch searched, between τ = 0 and the end. It is sampled at the corners' 40 base τ values, and a rise above 1e-12 N̄ beyond the first point that has lost 0.5% of the fall is a cause of its own ("blended N not monotone"). The N gate guarantees this per corner, and so for their blend, but the replay measures it.

**What is measured with it, before any runtime code.**
- E0: the share of cell updates this rule sends to CVODES, by cause, by region, over the whole run and the last ms. The split by w (0, between, 1) and by branch is also given.
- E5: the domain plus E0's rejects, against the Director's 1% bar in the nozzle.
- E4: the rule's cost per cell, on one thread.
- The nearer-branch rule's result stays in the replay's output, for comparison.
- *Expected (inferred, not a bar):* the chamber's equilibrium cells land at τ near 1. Their ΔY is then the blend's interpolation error in Y_eq at this table's spacing (Δη 0.2, Δ ln ρ 1.4), which may be near 0.01. A fail there is the axes' coarseness, not the placement. It is reported against the table's spacing, not tuned.

### The N placement measured: from 0.5 ms on, under the bar in every region; the start fails it (20 October 2026, about 23:10 CDT)

The same samples, table and joins as above; the replay's `n_rule` column (`table_b/e0_n_torch50_c1r_s_2026-10-20.txt`, `table_b/e5_n_torch50_c1r_s_2026-10-20.txt`; the run's log is `table_b/c1_torch50_highT_run_2026-10-20.log`).

**E5 with the N placement** (the domain plus the placement's rejects):

| Share of the region's cell updates | Whole run | 0 to 0.5 ms | 0.5 to 1 ms | Each 0.5 ms window from 1 to 4 ms | Last ms |
|---|---|---|---|---|---|
| Chamber | 2.867% | 32.400% | 0.000% | 0.000% | **0.000%** |
| Convergent | 2.789% | 31.515% | 0.000% | 0.000% | **0.000%** |
| Divergent | 3.143% | 35.514% | 0.004% | 0.000% | **0.000%** |
| All | 2.988% | 33.763% | 0.002% | 0.000% | **0.000%** |

- **Against the Director's 1% bar:** from 0.5 ms to the end the nozzle is at 0.004% at worst, and 0 from 1 ms on. So the burning engine passes, in every region.
- **The whole run fails, in every region,** and the fallback is all in the first 0.5 ms: the purge, the cold propellant and the ignition.
- The divergent's whole-run 3.143% splits into: the domain 2.082% (Z_N above the axis 1.890%, pure N2 0.173%, η 0.019%), and the placement's rejects 1.061%.

**The placement, unweighted** (368,183 samples inside the axes):
- 357,366 on the table (97.06%): 166,634 on the lower branch at w 1, 190,728 on the upper at w 1, and 4 on the lower at w 0.
- Not one sample in "no fall in N", "N below the blended N_cp" or "blended N not monotone".
- Illinois converges in 1 to 3 iterations for 99.3% of the searched cells, at most 9; none ran to 60.
- **The placed point's ΔY:**
  - upper branch, w 1: median 8.4e-4, max 2.2e-3;
  - lower branch, w 1: median 3.3e-3, 99% 5.7e-3. The expectation stated above (near 0.01, from the axes' spacing) is measured lower.
- The nearer-branch rule's result is unchanged in the same file, for comparison: 29,824 on the table.

**The rejects, each measured** (unweighted count; whole-run share of the divergent's updates):
- **Off the curve in s, below the blended start, at w 0** (9,539; 0.824%). *Measured:*
  - The cells are cold propellant before ignition: t 0.022 to 0.29 ms, T 126 to 465 K (median 250 K), η 0.59 to 0.78, Z_N median 0.002.
  - φ_N is −0.56 to −0.02 (median −0.11): the cell has more moles than the blended start. Its s is 2 to 2,412 J/(kg K) below the start (median 1,423).
  - The corners at f_H 1/7 and 1/6, Z_N 0 and η 0.6 start from the frozen family, not from the reactants. Their energy is below the reactants' at 200 K, so the builder's rule applies. Their start's N is 0.936 to 0.980 of the reactants', so it is partly burnt. At η 0.8 every start is the reactants (1.000).
  - So between η 0.6 and 0.8 the blended start is partly burnt, while the cell is still pure reactants. The cell is on neither side of the curve, and that is the axis's coarseness at the reactants' 200 K floor.
- **ΔY above ΔY_max, lower branch, w 1** (713; 0.215%):
  - The ignition front: t 0.222 to 0.306 ms, T 649 to 3,136 K, φ_N 0.01 to 0.999 (median 0.34), ΔY median 0.070.
  - *Inferred:* these are burning cells whose gas is a mixture of burnt and unburnt, off any one homogeneous reactor's path (probe 4's case).
- **ΔY above ΔY_max, lower branch, w 0** (328; 0.023%):
  - Before ignition: t 0.074 to 0.29 ms, T 294 to 2,061 K, Z_N median 0.22 (the purge's edge), ΔY median 0.15.
  - Not measured further.
- **Frozen reach** (185; 0.016%):
  - t 0.27 to 0.76 ms, T 2,881 to 3,579 K, φ_N 1.00 to 1.045.
  - Each has a frozen corner of nonzero weight: f_H 1/6, Z_N 0.3, η 0.6, whose lower branch from the reactants does not ignite by 1e15 s.
  - The gap (N̄_lo(1) − N̄_up(1))/N̄_lo(1) is 3e-7 to 1.5e-2 (median 1.3e-3). That is the 03:10 cause, as designed.

**E4 for the placement** (one thread, the machine at a load of about 2; `table_b/e0_n_...`):

| Per cell | ns |
|---|---|
| s_cell | 71 |
| s_cell and the nearer-branch rule | 248 |
| s_cell and the N placement, the runtime's form (each branch end read from its stored point, Illinois on the whole branch) | **1,255** |
| The same with the replay's monotone check (N̄ at the 40 base τ, and Illinois in their bracket) | 5,051 |

- The first version evaluated the four branch ends by the search, and cost 2,021 ns. Its outcomes were identical, sample for sample, by cause count.
- **Against the step:** C1's step outside the reaction is 3.70 µs per cell update. So the placement is 17% (*inferred*, on both threads) to 34% (on one) of it. E4's bar was stated for s, which passes at 1% to 2%. Applied to the whole placement, the 5% bar fails.
- **Against what it replaces:** CVODES costs 67.5 µs per cell update, 54 times the placement.
- **Against criterion 6, as E4 states for an overrun:** with the placement in place of CVODES, the step would be about 4.3 to 5.0 µs per cell update. C1 to 4 ms would then take about 27 to 31 s of wall, against 451 s now (*inferred*). This leaves out the update after the placement (the step in s and Y read at the new s, about one more blend of Y) and the fallback cells' CVODES.

**What this leaves.**
- The N placement replaces the nearer-branch rule. The steady engine is on the table: the chamber and convergent near equilibrium, the freezing nozzle on both branches.
- The start's 0 to 0.5 ms fails the bar, by four named causes, of which the largest two are the axes:
  - Z_N above 0.3: the purged ambient N2;
  - the η spacing at the reactants' 200 K floor: the cold propellant.
- Whether the 1% bar covers the start transient, and so whether these axes are extended, is the Director's.

## The runtime's step on the placed point, stated before its code (20 October 2026, about 23:10 CDT)

**The Director's decisions on 69bb304 (night brief, 20 October):**
1. The N placement replaces the nearer-branch rule. The replay keeps both, for comparison.
2. The 1% bar is a steady-state cost bar. A startup fallback is exact (CVODES), so it costs time, not accuracy. The startup share is reported separately. The two axis extensions (Z_N above 0.3; an η node at the reactants' 200 K floor) are done only if each takes under 30 minutes.
3. E4's 5% bar applies to the entropy evaluation, which passes (1 to 2%). The placement's cost is recorded as measured; optimising it is later speed work.

**The step, on the point the N placement gives** (engineering, within the design of *Stored values and lookup* and *The runtime*):
- **The coordinate.** The placement gives the branch b and τ; u = 1 − τ on the lower branch and τ − 1 on the upper. Both branches rise in s to s_eq, so on both the reaction moves u toward 0. On the blended curve s̄_b(u) = s̄_eq,b − span̄_b u² holds exactly, since it holds on every corner at the same u.
- **The rate.** λ̄_b(u) = exp(Σ′ w_q ln λ_q(u) / Σ′ w_q), with each corner's λ linear in u in the same bracket as its Y. This is the design's ln λ blended multilinearly across the other axes.
  - New: Σ′ runs over the corners of nonzero weight whose branch has nonzero span. A zero-span branch stores λ = 0 and its equilibrium Y at every u, so it carries no rate.
- **The update.** With d(s̄_eq − s̄)/dt = −λ̄ (s̄_eq − s̄), the coordinate obeys du/dt = −λ̄ u / 2. The exact exponential update is u_new = u exp(−κ_eff λ̄ Δt / 2), with κ_eff = 1 without PaSR.
- **Sub-steps** never cross a base u of the branch: the builder's 40 slots, common to every corner (lower u = 1 − ξ_p, upper u = j/8). In each sub-step:
  - λ̄ at its start, then a predictor to u_p (clamped at the next base u), then λ̄ at (u + u_p)/2, and the update with that λ̄;
  - if the update would cross the next base u, the cell stops there, the time to it, 2 ln(u/u_b) / λ̄_m, is taken from Δt, and the next sub-step starts.
  - The last interval ends at u = 0, which is never reached in finite time. So there are at most 40 sub-steps, and a stiff cell lands on the blended equilibrium at any Δt.
- **The composition.** The reacting species' Y = Ȳ_b(u_new), the blend at the new u. This is the design's snap: the fast directions are taken as relaxed. N2 and AR are the cell's own (AR ≤ 1e-12). The elements are the blend's, which are bilinear, so they are the cell's to rounding. e and ρ are kept, and T comes from the core's inversion (`Flow::setMassFractions`).
- **CVODES, counted by cause:** the placement's causes, and two new ones.
  - **"frozen corner on the branch":** a corner of nonzero weight whose placed branch is frozen. Its end is its reach, not s_eq, and its stored λ is relative to s_eq, so it does not fit the update in u. This is the 03:10 frozen-reach rule, made concrete in u.
  - **"no rate at a corner":** a corner in Σ′ whose λ at the cell's u is not positive and finite.
- **The cell's inputs.**
  - f_H and Z_N come from the element mass fractions. The species' element fractions are given by the adapter from Cantera at setup.
  - η is from the cell's own complete-combustion products at 200 K and 6,000 K, on the core's `Medium` (as the chamber study's samples).
  - s is from the core's NASA-7 data with the table's p° (the replay's evaluation, moved into the runtime's file until a `Medium` method is agreed with inquiry2).
  - e is from the conserved state, as `reactTable` takes it.
- **Where.**
  - `core/manifold_runtime.*`, with no Cantera.
  - `ReactingFlow::useManifoldTable` (FiniteRate) in the reaction adapter, with CVODES for every rejected cell.
  - `crucible_chamber_study fr` reads the table from `CRUCIBLE_MANIFOLD_TABLE`.
  - The PaSR closure under Table B (*The runtime*, κ_eff from τ_c along the curve) is not coded tonight, and `useManifoldTable` refuses a flow with the closure. C1 here is inviscid `fr` without it.

**Checks, with bars, before any of the code:**
- **R1. The runtime's placement is the replay's.** The replay's `n_rule` column is recomputed by the runtime's code on the same samples (`/tmp/c1s/states_torch50_TY.csv`, the C1-range table with s). Bar: the same cause for every sample; τ within 1e-9 where both place.
- **R2. The step at the cell, on the C1-range table,** at sampled cells placed on each branch (chamber, convergent, divergent):
  - (a) Δt = 0 returns the placed Ȳ(u_0) bit for bit.
  - (b) Δt = 1 s returns the blended equilibrium Ȳ_b(0) to 1e-14.
  - (c) The cell's H, O and N mass fractions are kept to 1e-12 relative.
  - (d) A stepped cell at w = 1, placed again, reads back its branch and τ to 1e-9. So repeated calls do not move it by the placement alone.
    - Measured, not barred: at w < 1 this cannot hold. There s_cell, the entropy of the blend's Y at the cell's (e, ρ), is not the blend's s̄, which is the chord found in E0.
  - (e) Measured, not barred: two calls of Δt/2 against one of Δt, for Δt from 1e-8 to 1e-5 s.
- **R3. C1 with Table B live,** against the CVODES run whose samples were replayed (`table_b/c1_torch50_highT_run_2026-10-20.log`): h2o2_highT, 32x6, `fr`, the 50 J torch, to 4 ms, 2 threads.
  - **Reported:**
    - c*, vacuum Isp and vacuum thrust at the end, and the injector pressure at each printed time, each against CVODES;
    - the wall time, and the reaction's µs per cell update;
    - Table B's share of cell updates, and CVODES's by cause, over the whole run and the last ms, with the share before 0.5 ms (the start) apart.
  - **Expected (*inferred*):**
    - The chamber's cells sit on the blended equilibrium of a table at Δη 0.2 and Δ ln ρ 1.4, where the placed ΔY was a median 3.3e-3 (lower branch) and 8.4e-4 (upper). If that is a temperature error of order 10 K at 3,350 K, c* moves by about 0.15%.
    - The wall time is 40 to 90 s, against 451 s.
  - **Bar (*guessed*): c*, Isp and thrust within 0.3% of the CVODES run.** A fail is measured before anything changes: the table's Y at the cells against their own Cantera equilibrium. The remedy would be the η and ln ρ spacing, not tuning.

### The runtime measured: R2 passes, R1 passes but for its τ bar, C1 live fails the 0.3% bar, and a steady fallback at the face fails the 1% bar (20 October 2026, about 23:50 CDT)

Coded as stated: `core/manifold_runtime.*`, `ReactingFlow::useManifoldTable`, and `crucible_chamber_study fr` with `CRUCIBLE_MANIFOLD_TABLE` (it checks the table's mechanism against the run's). The replay (`tools/table_b_branch_rule.cpp`) runs the runtime beside its own placement. Evidence:
- `table_b/r1r2_torch50_c1r_s_2026-10-20.txt`: R1, R2 and E4 on the CVODES run's samples.
- `table_b/c1_torch50_highT_tableB_live_2026-10-20.log`: C1 live, the timed run.
- `table_b/c1_torch50_highT_tableB_live_cells_2026-10-20.log`: the same run again, with the counts by region and by cell. Its c*, Isp, thrust, drifts and counts match the timed run's to every printed digit. Its wall time (113 s) is not a measure, since the fast ctest ran beside it.

The diagnostic below: the live run's samples (`states_live_TY.csv`) filtered to the last ms by region, every 50th kept, then `CRUCIBLE_E0_OWN_EQ_EVERY=1` on the replay. Cantera's `equilibrateUV` takes about 0.5 s a sample here (*measured*: 541 in 4 min), so the check is not run on every sample.

Command, from `build/`: `CRUCIBLE_MECHANISM=../tools/mechanisms/h2o2_highT.yaml CRUCIBLE_MANIFOLD_TABLE=/tmp/tbs/c1r_s.bin ./crucible_chamber_study fr 32 6 4e-3 2 <prefix> 50 1e-3`. The replay: `CRUCIBLE_E0_OWN_EQ_EVERY=0 ./crucible_table_b_branch_rule /tmp/tbs/c1r_s.bin ../tools/mechanisms/h2o2_highT.yaml /tmp/c1s/states_torch50_TY.csv <out>`. The own-equilibrium check, every 20th sample by default, takes most of its time otherwise.

**R1, the runtime's placement against the replay's** (357,366 CVODES samples placed by both on the same branch; measured):
- The same cause on every sample: **pass.**
- τ within 1e-9: **fail as written,** worst 1.24e-8. The worst samples are on the upper branch at τ 1.019 to 1.034, all at w = 1.
  - *Inferred:* both codes solve N to 1e-12 × N_cell by Illinois, from different brackets, and N is flat in τ near equilibrium. So 1e-8 in τ is inside the solves' own tolerance, and the bar was set below what that tolerance gives. Y moves by less than 1e-8 for it.
- On the live run's own states: the runtime adds its own cause (frozen corner on the branch), which the replay does not have. Worst |Δτ| 3.5e-7.

**R2, the step at the cell** (7,148 placed samples, every 50th; measured):
- (a) Δt = 0 returns the placed Ȳ bit for bit: **pass.** Two fixes to the check itself, not to the runtime. The samples' Y are printed to 10 digits and sum to 1 only within 1.0e-10, so y0 is normalised first. The cell is then placed again on that y0.
- (b) Δt = 1 s returns the blended equilibrium: max |ΔY| 0, **pass.**
- (c) The elements after 1e-6 s: worst relative 8.3e-16, **pass.**
- (d) At w = 1, placed again after 1e-6 s: worst |Δτ| 3.1e-12, all on the same branch, **pass.**
  - At w < 1 (1 sample): |Δτ| 1.56e-3. As stated, a step there is not idempotent. C1 live had 51 such updates in 18.5 million.
- (e) At w = 1, one step of Δt against two of Δt/2, max |ΔY|:

  | Δt | 1e-8 s | 1e-7 s | 1e-6 s | 1e-5 s |
  |---|---|---|---|---|
  | max \|ΔY\| | 3.3e-9 | 7.5e-7 | 1.45e-5 | 3.7e-5 |

  C1's step is 1.2e-7 s, so the step's own splitting error is of order 1e-6 in Y per step (*derived* from the row).
- s_cell against Cantera's `entropy_mass`: 1.8e-11 J/(kg K) at worst.

**E4, the runtime on one thread** (measured):
- place: 1,379 ns per cell;
- advance over 1e-7 s: 2,560 ns per cell;
- each includes e from T (15 ns), which the timing needs and the runtime does not.
- The entropy alone is 74 ns per cell, so it passes E4's 5% (decision 3). The placement's cost is recorded as measured.

**R3, C1 with Table B live, against the CVODES run** (h2o2_highT, 32x6, `fr`, the 50 J torch, to 4 ms, 2 threads; measured):

| | CVODES | Table B live | Table B / CVODES − 1 |
|---|---|---|---|
| c* [m/s] | 2476.61 | 2461.82 | −0.60% |
| Isp_vac [s] | 424.43 | 423.08 | −0.32% |
| F_vac [N] | 1664.67 | 1659.60 | −0.30% |
| Injector face [MPa] | 3.31243 | 3.29957 | −0.39% |
| Chamber-end equilibrium p0 [MPa] | 3.15291 | 3.13453 | −0.58% |
| Outlet mass flow [kg/s] | 0.399948 | 0.400005 | +0.01% |
| F_vac drift over the last ms | 1.58e-3 | 5.14e-3 | |
| March wall [s] | 451 | 69 | 6.5 times faster |
| Reaction [µs per cell update] | 67.5 | 9.87 | 6.8 times faster |

- **The 0.3% bar (*guessed*): fail.** c* is −0.60%. Isp (−0.32%) and thrust (−0.30%) are just past the bar.
- The injector pressure against CVODES at each printed time: −0.54% (0.5 ms), −0.08% (1 ms), −0.34%, −0.40%, −0.39%, −0.40%, −0.40% and −0.39% (4 ms).
  - So the deficit is steady at 0.4% from 2 ms on. It is not the start transient, though neither run is settled by the protocol's 0.1% rule.
- **The stated fail diagnostic: the live cells against their own Cantera equilibrium** (`equilibrateUV` at each sample's e and ρ, every 50th sample of the last ms; `table_b/c1_live_own_eq_{chamber,convergent}_2026-10-20.txt`; measured):

  | Last ms | Samples | T_cell − T_own (min; median; max) | max \|Y_cell − Y_own\| (median; max) | s_eq,own − s_cell (median) |
  |---|---|---|---|---|
  | Chamber (i ≤ 8) | 541 | −23.1; −23.0; −12.5 K | 3.27e-3; 3.30e-3 | 0.156 J/(kg K) |
  | Convergent (8 < i < 16) | 420 | −44.1; −23.3; −13.0 K | 3.32e-3; 6.08e-3 | 0.161 J/(kg K) |

  - So the chamber on the table is 12.5 to 23 K colder than its own equilibrium (the convergent section 13 to 44 K), against an expected order of 10 K. *Derived:* c* goes about as √T at a fixed molar mass, so 23 K at 3,350 K is about −0.34%. The rest of the −0.60% is the composition, not split out.
  - The CVODES run's chamber cells hold their own equilibrium to under 1e-4 J/(kg K) (*E0, E4 and E5 measured*). So the gap is the table's blended equilibrium, not the flow.
  - *Inferred cause:* the equilibrium Y is curved in η between nodes 0.2 apart, and the blend is linear. The dissociated share grows faster than linearly with T, so the chord holds too much of it, which leaves the cell colder at its fixed e. That is the stated remedy's direction: the η spacing (and ln ρ), not tuning.
- The wall time (69 s) is inside the expected 40 to 90 s. The end-of-run analysis takes another 156 s outside the march (measured from the files' times); the chamber study's own post-processing, not Table B's.
- **Table B's share of cell updates** (CVODES takes the rest, by cause):
  - **The start, to 0.5 ms: 72.96%**, by the start's causes in *The N placement measured* (the purge above the Z_N axis, cold propellant below the 200 K floor, the ignition front, frozen reach).
  - 0.5 to 1 ms: 99.96%.
  - From 1 ms the frozen-corner cause grows and holds: 1.46% (1 to 1.5 ms), 1.94%, 2.41%, 2.48%, 2.50%, and 2.50% (3.5 to 4 ms).
  - Over the run: 95.99%.
- **The last ms, by region: chamber 93.43%, convergent 97.02%, divergent 100%.** All the fallback is one cause, frozen corner on the branch. **The steady bar (1%) fails in the chamber (6.57%) and the convergent section (2.98%).** It passes in the divergent section.

**Where the steady fallback sits** (measured, then inferred):
- By cell over the last ms, 84 cells take any fallback.
  - The face column (i = 0) takes 8,500 of its 25,000 updates in every row: 34%, and exactly one update per step (8,500 steps in the last ms).
  - Column 1 takes 1,894 to 1,971 (about 8%).
- The face cells' end-of-step states (the samples) replay as table on the upper branch, at 3,354 K, η 0.574 and ρ 1.37 kg/m³. So the fallback is on a state the samples do not hold. *Inferred:* it is the second half-step's state, after the flow step brings in fresh propellant, once per step.
- The table's frozen lower-branch nodes are all at f_H 1/6, Z_N 0.3 and η 0.6, at every ρ. Its frozen upper-branch nodes are at η 0.2.
- The face cells' Z_N over the last ms is 2e-92 to 8e-70 (the purged N2's remnant; measured on the samples, columns 0 to 11 all below 1.4e-46). So any cell inside the Z_N 0 to 0.3 cell gives its Z_N 0.3 corners a weight of that size. The rule `corners[q].w > 0` then refuses the cell if any of them is frozen on the placed branch.
- *Inferred, not yet checked on the mid-step states:* the face cells after the inflow are placed on the lower branch near η 0.6. There the frozen corner (f_H 1/6, Z_N 0.3, η 0.6) has a weight of about 1e-70 and refuses them.

## The igniter and a guard above 6,000 K (20 October 2026, 19:15 to 19:45 CDT)

The Director's ruling: the 12,176 K igniter cell is unphysical. Real torch or spark igniters do not make 12,000 K bulk gas in a finite-volume cell. So both sides are fixed: the igniter is sized to a physical one, and a guard on h2o2_highT means the inversion can never fail.

**The guard (`core/medium.*`, `adapters/reaction.cpp`).**
- `Species::tHigh`: above it, cp is held at its tHigh value and h continues linearly, so e(T) keeps increasing. This continuation is labelled outside the model.
- The loader sets tHigh to the top of the fits when the top is 6,000 K or more (h2o2_highT). GRI h2o2.yaml's fits end at 3,500 K and keep their extrapolation, which is monotone to 40,000 K, so every result on h2o2.yaml is unchanged.
- The reaction step's CVODES integrates T with Cantera's own thermo, which extrapolates the fits above 6,000 K. So in a cell above 6,000 K the step and the flow disagree; the step's |ΔT| mismatch is already measured. Those cells are outside the model either way.
- `ReactingFlow` now also checks tHigh between the flow's medium and the mechanism's.

**The igniter, sized (sources from a Sonnet search, screened; no RL10 or Mascotte performance figure was read).**
- **A spark is sub-grid.**
  - Energies: 8 mJ (compact exciter) and 50 mJ (conventional) per spark (Tinker, NTRS 20210000596). Reliable multi-spark ignition needed 55 to 75 mJ (NTRS 20120011145). Generic exciters give 50 to 120 mJ at 100 to 300 sparks per second (US patent 10590887).
  - That is at most about 36 W. Its kernel (about 1 mm, *inferred*) is far below C1's 5.2 mm cells. Spread over the igniter's 3.1 cm³ it heats the gas by kilojoules per kilogram, and cannot stand for the spark.
  - So the model's igniter stands for a **torch**: hot gas at its flame temperature, about 3,100 K for an H2/O2 torch (Repas 1994, NASA TM-106493; taken from the indexed abstract, *provisional*).
- **The torch's power: a declared range, 50 to 160 kW for C1 (*inferred*; no sourced torch flow was found).**
  - The basis is a torch flow of 1 to 3% of C1's 0.4 kg/s (*guessed*) at the stoichiometric H2/O2 heating value, about 13.4 MJ/kg of mixture (120 MJ/kg of H2, one ninth by mass).
  - The duration is kept at 1 ms from 0.2 ms: a torch runs through the start transient, and C1 lights within it.
  - The run uses the low end, **50 kW (50 J over 1 ms)**.
  - The former 300 J over 1 ms is 300 kW, about twice the range's top.
- **What power cannot fix (*derived*).** At 0.2 ms the igniter's region holds propellant at about 0.007 kg/m³, because the nozzle is open to 1 kPa while the valve ramps. A volumetric deposit into that gas heats it at P / (ρ V), so any power drives a slow cell far above a torch's 3,100 K. A real torch adds its heat with its own mass, as hot gas.
  - The physical representation would be a hot-gas inflow, or a deposit that stops at the torch temperature. It is a change in `core/flow.*` (inquiry2's), and is not made here.

**The guard's unit check (`crucible_mixture_tests` 1b, stated in its header before it ran): passes.** On h2o2_highT at the igniter cell's composition:
- The top of the fits is 6,000 K.
- cv(20,000 K) equals cv(6,000 K) exactly, and cv just below 6,000 K agrees with it to 4.1e-15.
- e(20,000 K) − e(6,000 K) = cv × 14,000 K to 2.2e-16.
- The inversion at 8,000, 12,176 and 20,000 K is within 4.5e-12 K.
- gri30 (fits to 3,500 K) gets no continuation.

Fast tests: core, thermo_cea, reaction and mixture pass.

**C1 32x6 finite-rate on h2o2_highT, two runs, both past 0.206 ms (where the run failed before).** Samples in `/tmp/c1g/` (not in git).
- **The former 300 J igniter, as the guard's check:**
  - T_max 7,840 K at 0.216 ms. Above 6,000 K (on the guard) for the 8 samples from 0.212 to 0.226 ms.
  - Stopped at 0.828 ms with the igniter still on: my watcher read a half-written history line as 8.3 s. So it shows the guard carrying the old igniter through its hot peak, not to 4 ms.
- **The torch-sized igniter (50 kW for 1 ms), to 4 ms (`table_b/c1_fallback_torch_highT_2026-10-20.txt`):**
  - T_max 4,029 K at 0.288 ms, and 3,489 K after 1 ms. The guard never engages (no cell above 6,000 K). The peak is still above a torch's 3,100 K, as the derivation above expects for a volumetric deposit.
  - At 4 ms: injector pressure 3.3124 MPa, F_vac 1,663.9 N, Isp_vac 424.4 s on the supply flow. These are not compared with the h2o2.yaml run, which differs in both mechanism and igniter.

**The fallback rate on h2o2_highT, against the C1-range table with Y_H2O, by region (throat station 16):**

| Share of that region's cell updates | Whole run | Last ms |
|---|---|---|
| **Chamber and convergent (columns 0 to 15)** | **4.146%** | **0.000%** |
| — a corner failed a gate | 3.449% | 0 |
| — Z_N above the axis (the purged ambient N2) | 0.649% | 0 |
| — η above the axis (the igniter's cells) | 0.037% | 0 |
| — no H or O | 0.011% | 0 |
| Divergent (columns 16 on) | 86.040% | 82.300% |
| All cells | 45.093% | 41.150% |

- **Where the chamber's failed corners are:** all of them are in the first 0.5 ms (39.3% of that window's chamber updates, then 0.00% in every later 0.5 ms window). All are at ρ below 0.28 kg/m³, before the chamber pressurises, on the same low-density 2(a) turn-back nodes as the nozzle.
- **Against the h2o2.yaml run with the 300 J igniter (same split, e6d4378's samples):** chamber and convergent 9.237% over the whole run, of which 4.849% was the igniter's cells above the η axis. So the torch-sized igniter removes the η overflow (4.85% to 0.04%), and what remains of the chamber's transient fallback is the low-density turn-back.
- **So:** the burning chamber and the convergent section are on the table at every update. The start transient and the divergent nozzle are not, both through 2(a) at ρ < 0.28 kg/m³, which is what entropy is for (E5). The nozzle stays failing until entropy lands, as the Director said.

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
