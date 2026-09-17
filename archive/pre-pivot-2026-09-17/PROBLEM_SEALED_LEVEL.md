# The Sealed-Level Problem

*A formal statement, version 2. No prior knowledge is assumed beyond the definitions below.
Version 1 was reviewed by an independent mathematician; the defects found and the repairs made are
listed in §8.*

**The question in one sentence.** For the fine dynamical system F defined in §1, the coarse level
Π_h defined in §3, the program class Σ defined in §4, and the budgets (C, C_V): what is the sealing
number ε\*(𝒢; C, C_V) defined in §6, and which program attains a declared target, or what
obstruction bounds it from below?

---

## 0. Conventions

ℝ^d is Euclidean space. ⟨·⟩_T is the time average over an interval of length T. Functions are as
regular as needed for the equations to hold in the weak sense. "Computable at cost c" means
obtainable by a deterministic procedure using at most c elementary operations. Every map below may
take the value REFUSE; REFUSE is a legal value and is priced in §6. Norms: ‖·‖ on the performance
space ℝ³ is Euclidean; the norm on the coarse state and its dual are fixed by the verifier (§5) and
used consistently.

## 1. The fine system F

Let Ω ⊂ ℝ³ be a bounded domain invariant under rotations about a fixed axis, with boundary
∂Ω = Γ_w ⊔ Γ_in ⊔ Γ_e ⊔ Γ_c (four declared pieces: walls, inflow ports, one exit plane, conductors).
Let S be a finite label set; to each s ∈ S attach constants m_s > 0 and q_s ∈ ℝ. Let Ξ ⊂ ℝ³ \ Ω̄ be
a declared source set. Let μ₀, ε₀, g₀ > 0 be constants.

**Unknowns.** f_s : Ω × ℝ³ × [0, T] → [0, ∞) for each s ∈ S; vector fields E, B : ℝ³ × [0, T] → ℝ³.

**(F1) Transport.** For each s ∈ S,
∂_t f_s + v·∇_x f_s + (q_s/m_s)(E + v × B)·∇_v f_s = Σ_{s'∈S} C_{ss'}[f_s, f_{s'}] + R_s[f] − L_s[f],
where C_{ss'} are given bilinear operators conserving the mass, momentum, and energy of each pair;
R_s[f] = Σ over declared conversions (a, b) → (c, d) of ±∫∫ σ_{ab→cd}(|v − v'|) |v − v'| f_a(v) f_b(v')
(…) dv' are given bilinear operators with given kernels σ, each conversion releasing a given energy
Q_{ab→cd} distributed over the products by a given kernel; L_s[f] is a given nonnegative linear
operator (a sink).

**(F2) Fields.**
∇ × B = μ₀(J + J_ext) + μ₀ε₀ ∂_t E,  ∇ × E = −∂_t B,  ∇·E = ρ_c/ε₀,  ∇·B = 0,
with ρ_c = Σ_s q_s ∫ f_s dv, J = Σ_s q_s ∫ v f_s dv, and supp J_ext(·, t) ⊂ Ξ.

**(F3) Boundary conditions.** On Γ_w, for v·n < 0, f_s(x, v) = (K_w f(x, ·))_s(v), where K_w is a
given kernel acting on the outgoing part v·n > 0, with a given energy accommodation law and a given
label-conversion law. On Γ_in, for v·n < 0, f_s(x, v, t) = g_s(x, v, t), prescribed data. On Γ_e,
f_s(x, v) = 0 for v·n < 0. On Γ_c, E × n = 0, and the potential of E takes prescribed values V_k(t)
on the k-th connected component of Γ_c.

**(F4) Initial data.** z(0) = ({f_s(0)}, E(0), B(0)) ∈ Z₀, a declared set.

**The control.** G := (Ω, Ξ, J_ext, {V_k}, {g_s}) ranges over a family 𝒢 that is the image of a
compact subset of ℝ^{d_G} under a fixed continuous parametrization (so the function-valued entries
J_ext, V_k, g_s are finitely parametrized). Write z(t; G) for the solution.

**Assumption A0 (the sole external hypothesis).** F, with its given kernels C, σ, Q, L, K_w, is the
true dynamics of the system of interest. Nothing else in this document is assumed about the world.

**Assumption A1 (steady statistical operating point, uniform).** For every G ∈ 𝒢 and every
z(0) ∈ Z₀, F has a unique global weak solution, and there exist T₁(G) and a unique invariant
probability measure μ_G such that for every bounded observable Ψ, ⟨Ψ(z(t; G))⟩_{[T₁, T₁+T]} → ∫ Ψ dμ_G
as T → ∞, with covariances decaying on a correlation time τ_c(G), and **sup_{G∈𝒢} τ_c(G) < ∞ and
sup_{G∈𝒢} T₁(G) < ∞**. Where A1 fails at G, every quantity below evaluated at G equals REFUSE.

## 2. The functional to be predicted

For a state z, on the exit plane Γ_e with unit normal n pointing out of Ω, define the momentum flux
leaving the domain (kinetic flux minus the Maxwell stress),
Φ_p[z] = ∫_{Γ_e} [ Σ_s m_s ∫ v (v·n) f_s dv − T·n ] dA,  T := ε₀(E⊗E − ½|E|²I) + (B⊗B − ½|B|²I)/μ₀,
the mass flux Φ_m[z] = ∫_{Γ_e} Σ_s m_s ∫ (v·n) f_s dv dA, and the input power
P_in[z] = ∫_{Γ_in} Σ_s (m_s|v|²/2) (v·n) g_s dv dA + Σ_k V_k I_k − ∫_Ξ E·J_ext dx,
where I_k is the current into the k-th conductor. Define
p[z] := ( Φ_p·n,  Φ_p·n / (g₀ Φ_m),  (Φ_p·n)² / (2 Φ_m P_in) ),
with p := REFUSE if Φ_m ≤ 0 or P_in ≤ 0. For a declared dwell T_d, the **true performance** is
**P\*(G) := ⟨p[z(t; G)]⟩_{[T₁, T₁+T_d]}**, which under A1 equals ∫ p dμ_G up to the dwell term of §5.
P\* : 𝒢 → ℝ³ ∪ {REFUSE}.

## 3. The coarse level and the exact reduced identity

Let 𝒯_h be a mesh of Ω into cells c of diameter ≤ h. For a fine state z write, pointwise in x,
ρ_s(x) := m_s ∫ f_s dv and u_s(x) := (∫ v f_s dv)/(∫ f_s dv) (the fine-scale density and mean
velocity of label s), and ū_c := Σ_s ∫_c ρ_s u_s dx / Σ_s ∫_c ρ_s dx (the cell mean velocity).
Define the projection Π_h : z ↦ M = (M_c)_c, M_c ∈ ℝ^n, by
- the first moments ∫_c ∫ (1, v, |v|²/2) f_s dv dx for each s;
- the cell averages of the fields, written as a flux function ψ_c with B_pol = ∇ψ × ∇θ/(2πr),
  the azimuthal component b_c = (r B_θ)_c, and the potential φ_c of the curl-free part of E;
- the **sub-cell kinetic energy** k_c := ½ Σ_s ∫_c ρ_s(x) |u_s(x) − ū_c|² dx: the energy of the
  fine-scale mean-velocity field's variation within the cell (a second moment in x of a first
  moment in v; it is not the thermal energy and is not annihilated by any velocity-space lift);
- the **sub-cell dissipation rate** ε_c := −(d/dt) k_c restricted to the contribution of the
  operators C_{ss'} (a functional of z at one instant), and ω_c := ε_c/k_c (REFUSE if k_c = 0).

**The lift.** Λ_h : M ↦ z^eq(M) is the maximum-entropy state consistent with the first-moment
entries of M (a product of local Maxwell–Boltzmann states, uniform within each cell) together with
the cell-averaged fields. It factors through the first-moment and field entries of M only.

**The ideal generator.** A_h(M) := Π_h ∘ (right-hand side of F1–F2) ∘ Λ_h (M). It is fully
determined by F and 𝒯_h and contains no free element.

**Exact reduced identity (Mori–Zwanzig).** For every solution of F there exist a kernel K_G and a
process ξ_G, depending on G only through Ω, Ξ, ∂Ω and the kernels of F, such that
dM/dt = A_h(M) + ∫₀ᵗ K_G(M(τ), t − τ) dτ + ξ_G(t),   with  𝔼[ξ_G(t) | M(0)] = 0.
Under A1, ξ_G has a stationary mean m_ξ(M) and covariance Q_ξ(M) with correlation time ≤ τ_c(G).

**Information bound.** There is h_KS(G) > 0 such that no functional of M(0) determines ξ_G(t) for
t ≫ 1/h_KS, and for F (which has a continuum of scales with faster dynamics at smaller scales) the
prediction horizon for z is finite uniformly in the precision of z(0). This bounds trajectories,
not ergodic averages, and is not used as an obstruction below.

## 4. The program class Σ

A candidate σ is a finite program in a **fixed total flux-form language ℒ**: every program in ℒ
terminates, contains no randomness primitive, and can only write the reduced generator as
"A_h(M) + Σ over faces of (flux expressions) + Σ over cells of (source expressions)" with the
projected sources of F fixed; conservation (S5 i) and ∇_h·B = 0 (S5 ii) are therefore properties of
the syntax, and determinism (S5 v) is a property of ℒ. A candidate specifies:

- **(S1) Augmented state.** M^σ = (M, a_σ(z)), where a_σ is a declared finite list of additional
  functionals of z with explicit definitions, so that their true values exist.
- **(S2) Closure.** A map Φ_σ : (M^σ restricted to a stencil of radius ρ_σ about c; θ) ↦ ℝ^{n_σ},
  continuous in all arguments and Lipschitz in θ, with θ ∈ Θ_σ ⊂ ℝ^{p_σ} compact. The reduced
  generator is **dM^σ/dt = A_h(M) + Φ_σ(M^σ; θ).**
- **(S3) Side maps.** A finite list Ψ_σ^j, each mapping a lower-dimensional restriction of M^σ to
  values used by Φ_σ, each with a declared cost.
- **(S4) Validity predicate.** V_σ : M^σ_c ↦ {0, 1}, computable from M^σ alone.
- **(S5) Admissibility.** (i) discrete conservation and (ii) ∇_h·B = 0 (syntactic in ℒ);
  (iii) sealed-corner condition: there is a declared region Z_chem of state space on which a
  pre-existing operator 𝒮₀ is certified to satisfy |P_{𝒮₀} − P\*| ≤ ε_chem, and on Z_chem,
  Φ_σ(M^σ; θ) coincides with 𝒮₀'s closure within ε_chem for every θ ∈ Θ_σ; (iv) no term of Φ_σ is
  defined by a predicate on M^σ (continuity on the whole state); (v) determinism (syntactic in ℒ);
  (vi) positivity of masses and energies is preserved; **(vii) coverage: P_{σ,θ}(G) ≠ REFUSE for all
  G in a declared required subset 𝒢_req ⊂ 𝒢 and all θ ∈ Θ_σ.**
- **Cost.** cost(σ, h) := elementary operations per unit physical time per cell, side maps included.

## 5. The verifier and the certified bound

The verifier 𝒱 has access to:
(a) the **box oracle** 𝒪_box(Bx, ∂-data, T), returning at declared cost c_box(|Bx|, T) the stationary
statistics under F of any observable on a sub-box Bx ⊂ Ω of diameter ≤ ℓ with prescribed boundary
data, averaged over T;
(b) the **analytic set** 𝒜 = {(G_k, exact solutions of F)};
(c) a finite **data set** D = {(G_i, O_i, y_i, Σ_i)}: values y_i with noise covariance Σ_i of
observables O_i (each a functional of Π_h z) of F on full domains, partitioned once and permanently
into D_cal ⊔ D_val with no G in both; D may be empty on a subset of 𝒢;
(d) **the fluctuation ensemble ℰ**: a fixed rule (part of the problem, not of σ) assigning to any
(k, ω) a probability law on sub-cell velocity fields, namely the maximum-entropy Gaussian random field
with the prescribed energy k, decorrelation rate ω, and a declared spectral shape. This makes the
augmented lift Λ_h^σ : M^σ ↦ (a law on fine states) single-valued.

**All constants below are computed by 𝒱, never declared by σ.** 𝒱 holds the box-sampling seeds; σ
never sees them.

For (σ, θ, G), 𝒱 computes:
- **(V1)** the reduced prediction P_{σ,θ}(G) by running the reduced generator to its dwell, or
  REFUSE(T_max) if no stationary state is reached by a declared T_max;
- **(V2)** the feasible set Θ\*_σ := {θ ∈ Θ_σ : σ passes 𝒜, and ‖Σ_i^{−1/2}(O_i[P_{σ,θ}] − y_i)‖ ≤ κ
  for all i ∈ D_cal}; if Θ\*_σ = ∅, σ is inadmissible;
- **(V3)** the closure residual r_σ(G) := sup_j ‖Φ_σ(M^σ_j; θ) − Φ^true_j‖_\* over a sample of boxes
  {Bx_j} drawn (i) from the reduced trajectory's own states and (ii) from the true states available
  in 𝒜 and D, by 𝒱's rule and seeds; Φ^true_j := 𝔼_ℰ[ Π_h(RHS of F at the box's stationary state) ]
  − A_h(M_j), the expectation over the fluctuation ensemble; ‖·‖_\* is the norm dual to the
  amplification of (V4);
- **(V4)** the amplification Λ_σ(G) := ‖(∂p/∂M) J^{−1}‖ with J the Jacobian of A_h + Φ_σ at the
  reduced stationary state; 𝒱 also computes σ_min(J) and returns REFUSE if σ_min(J) < j₀ (a declared
  margin), so that the reduced fixed point is certified hyperbolic with margin;
- **(V5)** the invalid fraction ν_σ(G) := Σ_{c: V_σ(M_c)=0} ‖∂p/∂M_c‖ / Σ_c ‖∂p/∂M_c‖ (the share of the
  performance sensitivity carried by invalid cells), and the Lipschitz constant L_σ of Φ^true on the
  invalid set **estimated by 𝒱 from finite differences of 𝒪_box** with its own confidence bound;
- **(V6)** the discretization term e_h(σ, G) **computed by 𝒱 by mesh refinement** (Richardson
  estimate on h, h/2);
- **(V7)** the dwell term τ(T_d) := (2 ∫₀^∞ C_p(s) ds / T_d)^{1/2}, where C_p is the autocovariance of
  p along the stationary reduced trajectory augmented by the box-measured orthogonal noise
  (the Green–Kubo integrated autocorrelation, resolved fluctuations included);
- **(V8)** the remainder constant ρ_σ **computed by 𝒱** from the second derivative of A_h + Φ_σ and the
  variation of J along the perturbation path;
- **(V9)** **validation**: for all j ∈ D_val, ‖Σ_j^{−1/2}(O_j[P_{σ,θ}] − y_j)‖ ≤ κ + ‖Σ_j^{−1/2}‖ ε_σ(G_j)
  for all θ ∈ Θ\*_σ; a σ that fails (V9) is **falsified and inadmissible**.

**Certified bound.**
**ε_σ(G) := sup_{θ∈Θ\*_σ} { e_h + Λ_σ(G)·[ r_σ(G) + L_σ ν_σ(G) ] + ρ_σ r_σ(G)² + τ(T_d) }.**
If P_{σ,θ}(G) = REFUSE for some θ ∈ Θ\*_σ, then ε_σ(G) := ε_R, the declared prior range of p over 𝒢.

**Assumption A2 (coverage, quantitative).** Φ_σ is Lipschitz in θ with a modulus 𝒱 can compute, and
Θ_σ has covering number 𝒩(Θ_σ, η). Then with probability ≥ 1 − δ(N), δ(N) := 2𝒩(Θ_σ, η) e^{−2Nη²},
the residual on the N sampled boxes bounds the residual along the reduced trajectory to within η.
If K programs are evaluated adaptively against the same seeds, δ is replaced by K δ (union bound);
𝒱 refreshes seeds for the final evaluation.

**Assumption A3 (hyperbolicity of the true operating point, uniform).** For all G ∈ 𝒢, the true
stationary state of F is hyperbolic and lies in the basin of the reduced stationary state under the
perturbation Φ^true − Φ_σ for residuals up to r_σ; the constants are uniform over 𝒢. 𝒱 can certify
the reduced side of this (V4) but not the true side; A3 is therefore a hypothesis of the theorem,
stated as such (see R3 g).

**Assumption A4 (continuity in G).** G ↦ ε_σ(G) admits a modulus of continuity ω_ε computable by 𝒱,
so that sup_{G∈𝒢} is evaluable on a finite net of 𝒢 at finite cost.

**Claim that any solution must establish.** Under A0–A4, for all θ ∈ Θ\*_σ,
|P_{σ,θ}(G) − P\*(G)| ≤ ε_σ(G)  with probability ≥ 1 − δ(N).

**Certified ranking.** (G_a, G_b) ∈ 𝓡_σ iff |P_σ(G_a) − P_σ(G_b)| > ε_σ(G_a) + ε_σ(G_b) holds for
every θ ∈ Θ\*_σ. The **geometry signal** is |P_σ(G_a) − P_σ(G_b)|; the band swallows it exactly when
this inequality fails.

## 6. The question

Given F (under A0), 𝒢 with its required subset 𝒢_req, 𝒯_h, the language ℒ and class Σ, the oracle
cost c_box, the budgets C (runtime per prediction) and C_V (verification), D, 𝒜, ℰ, and the declared
constants (κ, ε_R, T_d, T_max, ℓ, N, j₀, η):

**Definition (sealing number).**
ε\*(𝒢; C, C_V) := inf { sup_{G∈𝒢} ε_σ(G) : σ ∈ Σ admissible (S5 i–vii), not falsified (V9),
cost(σ, h) ≤ C, verification cost ≤ C_V }.
Because (S5 vii) forbids refusal on 𝒢_req, the trivial refuser is not admissible, and ε\* is defined
jointly with the ranked sets below. The level Π_h is *sealed to ε* on 𝒢 at budget (C, C_V) if ε\* ≤ ε.

**(Q1) Value.** Determine ε\*(𝒢; C, C_V), and the maximal certified-ranked set
𝓡\*(𝒢; C, C_V) := ⋃ 𝓡_σ over the same set of programs.

**(Q2) Extremizer to a target.** For a declared ε_target, exhibit a program σ, admissible and not
falsified, with sup_{G∈𝒢} ε_σ(G) ≤ ε_target, together with its certificate (V1–V9) and its ranked set
𝓡_σ; **or**

**(Q3) Obstruction.** Prove a lower bound ε\*(𝒢; C, C_V) ≥ ε₀ > 0 at the given budgets, together with
the largest 𝓡\* that remains certifiable.

## 7. Remarks

**(R1) Not vacuous, not impossible.** The program σ_DNS whose augmentation a_σ carries, per cell, a
fine-grid state of F (so that Φ_σ reproduces the exact projected rate) is in Σ at a cost that grows
without bound as the fine grid is refined, and has ε = e_h + τ in the limit. Hence ε\*(C) is finite
and nonincreasing in C, with limit e_h + τ as C → ∞. The content of the problem is the value at
finite C. (The exact closure of F is not itself in Σ: it is a functional of z with memory, not a
continuous function of a finite stencil of M^σ.)

**(R2) Source of finite-cost solutions.** On any region of state space where the linearized
collision operators have a uniform spectral gap, no marginal field modes exist, and the kernel K_G
decays on a time τ_K much shorter than the relaxation time τ_M of M (a normally hyperbolic slow
manifold), a closure with a single exponential kernel realized through (k, ω) has residual
r = O((τ_K/τ_M)²) uniformly on that region, with the constant measurable by 𝒪_box as the second
continued-fraction coefficient of the memory kernel. The sealing number is small exactly to the
extent that 𝒢 keeps the reduced trajectory in such regions, which (S4) records cell by cell.

**(R3) Obstructions anticipated, and where the statement stands relative to each.**
(a) *Self-reference (incompleteness, undefinability).* Accuracy is a statement about z, not
expressible in the language of M; the validity predicate V_σ is a syntactic guard, not a proof. The
certificate is relative to an external verifier with oracle access to F (𝒪_box, 𝒜, D), never
issued by the reduced system about itself.
(b) *Undecidability over programs.* "σ attains ε on all G" is not decidable over a Turing-complete
language; ℒ is total and flux-form so that (S5 i, ii, v) are syntactic, and Q1 asks for a bound with
a finite-cost evaluator; Q3 admits an obstruction as an answer. ε\* is an infimum over a countable
set and need not be attained; Q2 asks for a declared target, not for the infimum.
(c) *Finite horizon.* P\* is an ergodic average; the dwell enters through τ(T_d); no pointwise
prediction of z is requested; the information bound of §3 bounds trajectories only.
(d) *Non-termination.* Whether the march reaches a stationary state is not decidable in general;
REFUSE(T_max) is legal and priced, and forbidden on 𝒢_req.
(e) *The calibration loop.* θ is fixed by data generated by the system it models; this is closed once
as Θ\*_σ from D_cal; D_val only falsifies (V9) and never moves θ.
(f) *Levels.* The problem measures how sealed the level Π_h is: how much of the dynamics at that
level is expressible in that level's own variables. The exact identity of §3 says every level is
sealable at the price of a memory and a noise term; the sealing number prices what remains after
finite memory and dropped noise.
(g) *Bifurcations.* The bound is a linear-response bound at a hyperbolic operating point (A3, V4). It
certifies nothing across a bifurcation of the true system in G. Geometries at which the true system
changes stability are precisely those the certificate cannot cover; the honest output there is
REFUSE, and 𝒢_req must be chosen inside a single basin with margin j₀. This is the sharpest
limitation of the statement and is not hidden.
(h) *Resolution of the oracle.* 𝒪_box sees boxes of diameter ≤ ℓ. Two admissible closures that differ
by a conservative term statistically invisible on every ℓ-box but with an O(1) effect on a
domain-scale stationary structure are indistinguishable to 𝒱: this yields a lower bound
ε\* ≥ ε₀(ℓ), decreasing in ℓ, and is the expected shape of a Q3 result.

**(R4) What a solution implies.** A solution of Q2 with sup_G ε_σ(G) ≤ ε_target yields, under
A0–A4, a computable map P_σ : 𝒢 → ℝ³ ∪ {REFUSE} with |P_σ − P\*| ≤ ε_target on 𝒢 (probability
≥ 1 − δ) and the certified ranking 𝓡_σ. A solution of Q3 yields the certified statement that no
program in Σ at budget C does better than ε₀ on 𝒢, together with the ranking 𝓡\* that is still
certifiable. In either case the question of whether the coarse instrument's output converges on the
true performance of every geometry in 𝒢_req is settled under A0 and the stated regularity
hypotheses A1–A4, of which A3 (the true operating point is hyperbolic, uniformly) is the one that
carries physics and is not verifiable from the coarse side.

## 8. Defects found in version 1 by independent review, and the repairs

| Defect | Repair |
|---|---|
| k_c was defined as fluctuation about a velocity-space maximum-entropy lift matching energy, so k_c ≡ 0 and ω_c = 0/0 | k_c redefined as the sub-cell kinetic energy of the fine-scale mean-velocity field about the cell mean (a spatial second moment); ε_c as its collisional dissipation |
| Λ_h^σ was one-to-many (two scalars do not fix a fluctuation field), making Φ^true a set and opening a gaming channel inside the verifier | the fluctuation ensemble ℰ is declared by the problem (maximum-entropy Gaussian field with fixed spectral shape); Φ^true is an expectation over ℰ |
| L_σ, ρ_σ, e_h declared by the program, hence tunable to shrink the bound | all three computed by 𝒱 (finite differences of the oracle, Richardson refinement, second derivatives) |
| box seeds visible to σ; D_val never entered the bound | 𝒱 holds seeds and refreshes them for final evaluation; (V9) validation added; falsified programs are inadmissible |
| A2 assumed uniform convergence over a class of unbounded capacity | replaced by a covering-number bound with δ(N) = 2𝒩 e^{−2Nη²} and a union bound over K adaptive evaluations |
| (S5 i, v) undecidable in a general language (Rice) | Σ restricted to a total flux-form language ℒ where conservation, ∇·B = 0, and determinism are syntactic |
| the trivial refuser attained ε = ε_R at zero cost | (S5 vii): refusal forbidden on a declared required subset 𝒢_req; ε\* defined jointly with the ranked set |
| Q2's (1 + η)ε\* was unverifiable without a matching lower bound | Q2 restated to a declared target ε_target |
| R1's witness σ = F was not in Σ | witness replaced by σ_DNS (a per-cell fine-grid augmentation) |
| ν_σ undefined; norms undeclared; τ(T_d) used only the orthogonal noise | ν_σ defined by sensitivity share; r_σ in the dual norm of the amplification; τ from the Green–Kubo integral including resolved fluctuations |
| Maxwell-stress sign ambiguous; no REFUSE for Φ_m ≤ 0, P_in ≤ 0 | Φ_p written as kinetic flux minus T·n with T defined; REFUSE clauses added |
| 𝒢 "finitely parametrized" while containing functions; sup over 𝒢 not evaluable | 𝒢 defined as a continuous image of a compact finite-dimensional set; A4 modulus of continuity added |
| A1 gave τ_c(G) < ∞ pointwise only | uniformity over 𝒢 added to A1 |
| A3 unverifiable from the coarse side and silently carrying the physics | stated as a hypothesis of the theorem with the bifurcation caveat R3(g); the reduced side certified by V4 with margin j₀ |
