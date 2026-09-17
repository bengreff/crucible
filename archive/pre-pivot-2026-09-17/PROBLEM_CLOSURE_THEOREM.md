# The Closure Theorem Problem

*Version 5 (2026-09-10). Versions 1–3 were reviewed cold by independent mathematicians; §10 lists
the third review's findings and repairs. Version 4 was submitted to a solver, which disproved it by
an admissible but unphysical fine system: a loss operator cancelling the collisional gain, a flow
at Mach ε⁻³ so that particles transit with ~10⁻⁴ collisions while the thermal Knudsen number is
small, and two non-Maxwellian inlets with identical moments but different pair relative speeds,
hence different thrust, invisible to the coarse inlet condition. Version 5 adds the three
hypotheses that make local equilibrium reachable (H2 i–iii). The statement is a single inequality
with two exponents, one for each small parameter the coarse description carries.*

**Prompt.** Decide Statement (T) in §6 and determine the sharp exponents (β\*, γ\*) of §7. No context
beyond this document is required or assumed.

---

## 0. Conventions

ℝ^d is Euclidean space; I the identity tensor; ⊗ the outer product; ∇ the gradient in x ∈ ℝ³; ê_z
the unit vector along the symmetry axis; ‖·‖_{L²(Ω)} the L² norm on Ω. Positive constants μ₀,
ε_vac (vacuum permittivity), g₀ are fixed. A finite label set S is fixed; each s ∈ S carries a mass
m_s > 0 and a charge q_s ∈ ℝ; one distinguished label e ∈ S has q_e < 0 and m_e ≤ m_s/1000 for all
s ≠ e; S_h := S \ {e} are the heavy labels; S_h⁺ := {s ∈ S_h : q_s ≠ 0}; e_c := |q_e| (the elementary
charge; the label e and the charge e_c are distinct symbols); Z_s := q_s/e_c. "The Chapman–Enskog
tensors of an operator in a field B" means the full first-order Chapman–Enskog stress tensor (five
viscosity coefficients), heat-flux tensor (parallel, perpendicular, cross), diffusion coefficients,
resistivity η, and inter-label energy exchange rates, defined by the linearized operator about its
local equilibria. sup ∅ := 0 for all suprema of nonnegative quantities below.

## 1. The fine system F

Let Ω ⊂ ℝ³ be a bounded domain invariant under rotations about ê_z, with boundary
∂Ω = Γ_w ⊔ Γ_in ⊔ Γ_e ⊔ Γ_c (walls, inflow ports, one planar exit face normal to ê_z, conductors
Γ_c = ⊔_j Γ_{c,j}). Let Ξ ⊂ ℝ³ \ Ω̄ be a bounded source set; L_G := diam Ω; ℓ_mfp(G) the supremum
over Ω of the mean free path under μ_G (H1); U_G := sup_Ω |ū| under μ_G (§2), the reference speed.

**Unknowns.** f_s : Ω × ℝ³ × [0, ∞) → [0, ∞) for s ∈ S; E, B : ℝ³ × [0, ∞) → ℝ³; surface charge
densities ς_j on Γ_{c,j}.

**(F1)** For each s ∈ S,
∂_t f_s + v·∇_x f_s + (q_s/m_s)(E + v × B)·∇_v f_s = Σ_{s'∈S} C_{ss'}[f_s, f_{s'}] + R_s[f] − L_s[f],
where the C_{ss'} are given bilinear collision operators, each conserving the mass, momentum and
energy of the pair and satisfying an H-theorem; R_s[f] is a given bilinear reaction operator built
from given kernels σ_{ab→cd}(|v − v'|) with given energy releases Q_{ab→cd} ∈ ℝ (positive for
exothermic, negative for endothermic conversions; exchanged with the kinetic energy of the
products) and product kernels, conserving total mass, momentum, and energy including Q; a product
not in S carries its share of Q out of the system and is accounted in L; L_s[f] is a given
nonnegative operator, at most bilinear in f.

**(F2)** On ℝ³: ∇ × B = μ₀(J + J_ext) + μ₀ε_vac ∂_t E, ∇ × E = −∂_t B, ∇·E = (ρ_c + Σ_j ς_j δ_{Γ_{c,j}})/ε_vac,
∇·B = 0, with ρ_c := Σ_s q_s ∫ f_s dv and J := Σ_s q_s ∫ v f_s dv on Ω (zero outside),
supp J_ext ⊂ Ξ, and the Silver–Müller radiation condition at infinity. Coulomb gauge:
E = −∇φ_F − ∂_t A, B = ∇ × A, ∇·A = 0, φ_F → 0 at infinity.

**(F3)** On Γ_w ∪ Γ_c, for v·n < 0, f_s(x, v) = (K_w f(x, ·))_s(v), with K_w a given kernel acting on
the outgoing part, carrying a given wall temperature T_w(x) and a given label-conversion law; on
Γ_in, for v·n < 0, f_s = g_s (given); on Γ_e, f_s = 0 for v·n < 0; on Γ_{c,j}: E × n = 0 and
φ_F = V_j (given), with ς_j the unknown making this consistent.

**(F4)** z(0) ∈ Z₀ (constrained in H1).

**Control.** G := (Ω, Ξ, J_ext, {V_j}, {g_s}, T_w, Z₀), with J_ext, V_j, g_s time-independent. Write
z(t; G, z₀) for the solution from z₀ ∈ Z₀.

**Fine performance functionals.** With T := ε_vac(E⊗E − ½|E|²I) + (B⊗B − ½|B|²I)/μ₀ and
I_j := ∫_{Γ_{c,j}} J·n dA,
Φ_p[z] := ∫_{Γ_e} ê_z·[ Σ_s m_s ∫ v (v·ê_z) f_s dv − T·ê_z ] dA,  Φ_m[z] := ∫_{Γ_e} Σ_s m_s ∫ (v·ê_z) f_s dv dA,
P_in[z] := ∫_{Γ_in} ∫_{v·n<0} Σ_s (m_s|v|²/2)(−v·n) g_s dv dA + Σ_j V_j I_j − ∫_Ξ E·J_ext dx,
p[z] := ( Φ_p,  Φ_p/(g₀Φ_m),  Φ_p²/(2Φ_m P_in) ) ∈ ℝ³.

**Fine performance:** the ergodic mean **P_F(G) := ∫ p dμ_G** (H1), componentwise.

**Intrinsic relative variability.** For each component i and the dwell T_d of §5,
**Δ_{F,i}(G) := ( Var_{μ_G}[ ⟨p_i⟩_{[0,T_d]} ] )^{1/2} / |P_{F,i}(G)|**: the relative root-mean-square
spread of the system's own dwell averages under its stationary measure.

## 2. The coarse variables (at scale h)

Fix a cell size h > 0 and a partition 𝒯_h of Ω into axisymmetric (annular) cells of diameter ≤ h;
⟨·⟩_c is the cell average (including the azimuthal average). For x ∈ Ω define from z the pointwise
moments n_s := ∫ f_s dv, u_s := n_s^{-1}∫ v f_s dv, ρ_s := m_s n_s, ρ := Σ_s ρ_s, u := ρ^{-1}Σ_s ρ_s u_s,
T_h := (Σ_{s∈S_h} m_s ∫ |v − u|² f_s dv)/(3 Σ_{s∈S_h} n_s), T_e := (m_e ∫ |v − u|² f_e dv)/(3 n_e),
and ū_c := ⟨ρu⟩_c/⟨ρ⟩_c. Per cell:
- **sub-cell kinetic energy per unit mass** k := ½ ⟨ Σ_s ρ_s |u_s − ū_c|² ⟩_c / ⟨ρ⟩_c;
- **its viscous dissipation** ε_k := ⟨ Π_CE(f) : ∇(u − ū_c) ⟩_c / ⟨ρ⟩_c ≥ 0;
- ω := max( ε_k/k, ν_CE/h² ) if k > 0, and ω := ν_CE/h² if k = 0, where ν_CE := μ/ρ is the local
  Chapman–Enskog kinematic viscosity (the floor is the viscous decay rate at the cell scale; it is
  computable from the coarse state alone, so the coarse system depends on G only through M).
Write B = ∇ψ × ∇θ/(2πr) + (b/r) e_θ in cylindrical coordinates about ê_z, defining ψ and b = rB_θ.
The **coarse state** is M := (ρ_s, ρu, ρe_h, ρe_e, ψ, b, k, ω) as cell-averaged fields, with
e_h := (3/2)(Σ_{s∈S_h} n_s) T_h/ρ and e_e := (3/2) n_e T_e/ρ. Π_h : z ↦ M is the resulting map. The
potential φ is determined by (C7), not stored.

**Fluctuation intensity.** Let Ω_h^int be the union of cells not adjacent to Γ_w ∪ Γ_c, and
Ω_h^flow := {c ∈ Ω_h^int : |ū_c| ≥ δ_u U_G} with δ_u := 10^{-2}. Define
𝓀(G) := sup_{c ∈ Ω_h^flow} k / U_G²  (under μ_G).

## 3. The coarse system 𝒞_θ (at scale h)

**Free constants.** θ := (C_S, C_p, C_ε, C_{ω1}, C_{ω2}, c_μ, Pr_t, Sc_t, C_sh, σ_k, κ_vK, B_log, A_J) ∈
Θ := [10^{-3}, 10^3]^{13}, independent of G, h, and ε. (κ_vK, B_log, A_J are the wall-law constants;
A_J multiplies the Jayatilleke function P(Pr) := 9.24[(Pr/Pr_t)^{3/4} − 1][1 + 0.28 exp(−0.007 Pr/Pr_t)].)

Coarse electromagnetic quantities: J := (∇ × B)/μ₀; the coarse system is quasi-neutral. All
transport coefficients are the Chapman–Enskog tensors of the C_{ss'} in the local field at the local
state: Π_CE, 𝒦_h, 𝒦_e, D_s, η, ν_eq; the rates ṅ_s, Q̇_h, Q̇_e, L̇_h, L̇_e are the moments of R_s and L_s
at the local equilibrium with moments M. b̂ := B/|B|, P_∥ := b̂⊗b̂, P_⊥ := I − P_∥ where α ≥ ε, and
isotropic tensors where α < ε; S := ½(∇u + ∇uᵀ) − ⅓(∇·u)I; α := n_e/Σ_{s∈S_h} n_s;
c_{p,h} := (5/2)(Σ_{s∈S_h} n_s)/ρ, c_{p,e} := (5/2) n_e/ρ; ν_CE := μ/ρ the kinematic viscosity.

**Turbulent coefficients.** ν_t := c_μ k/ω; κ_t := ν_t/Pr_t; D_t := ν_t/Sc_t; tensors
ν_t P_⊥ + ν_t(1 − |b̂|²)P_∥ etc.

**(C1) Mass, per label.** ∂_t ρ_s + ∇·(ρ_s u + 𝒥_s) = m_s ṅ_s + W_s,
𝒥_s := −(D_s P_∥ + (D_s + D_t)P_⊥)·∇ρ_s + (ρ_s/ρ) Σ_{s'} (D_{s'} P_∥ + (D_{s'} + D_t)P_⊥)·∇ρ_{s'} (Σ_s 𝒥_s = 0),
W_s the wall source of (C11).

**(C2) Momentum.** ∂_t(ρu) + ∇·(ρu⊗u + (p + ⅔ρk) I) = ∇·Π + J × B,
Π := Π_CE + 2ρν_t P_⊥ S P_⊥, p := Σ_{s∈S_h} n_s T_h + n_e T_e.

**(C3) Heavy internal energy.** ∂_t(ρe_h) + ∇·(ρe_h u) + p_h ∇·u = −∇·q_h + Π_h:∇u + ν_eq(T_e − T_h)
+ Q̇_h − L̇_h + ρ C_ε k ω,  q_h := −(𝒦_h + ρ c_{p,h} κ_t P_⊥)·∇T_h.

**(C4) Light internal energy.** ∂_t(ρe_e) + ∇·(ρe_e u) + p_e ∇·u = −∇·q_e + η|J|² − ν_eq(T_e − T_h)
+ Q̇_e − L̇_e,  q_e := −(𝒦_e + ρ c_{p,e} κ_t P_⊥)·∇T_e, parallel part replaced by (C12) on open lines.

**(C5) Generalized Ohm's law.** E = −u × B + ηJ + E_H, E_H := (J × B − ∇p_e)/(e_c n_e), E_H := 0 where
n_e = 0.

**(C6) Induction.** ∂_t ψ = −2πr E_θ, ∂_t b = −r (∇ × E)_θ, with E from (C5); the external field
∇ × B_ext = μ₀J_ext enters through the boundary values of ψ and b on ∂Ω (time-independent).
∇·B = 0 identically.

**(C7) Potential.** ∇²φ = −∇·E on Ω with E from (C5); on Γ_{c,j} the current condition of (C11)
with φ_w = V_j (Robin-type); on Γ_w, J·n = 0; on Γ_in ∪ Γ_e, ∂_nφ = 0; with no conductor, ∫_Ω φ = 0.

**(C8)** ∂_t(ρk) + ∇·(ρuk) = ∇·[ρ(σ_k ν_t)P_⊥·∇k] + ρΓk − ρC_ε kω − ρC_sh ω_E k,
**(C9)** ∂_t(ρω) + ∇·(ρuω) = ∇·[ρ(σ_k ν_t)P_⊥·∇ω] + ρC_{ω1}Γω − ρC_{ω2}ω², with ω ≥ ν_CE/h² enforced
as a floor (so ν_t ≤ c_μ k h²/ν_CE is bounded by coarse quantities); ω_E := |∇v_E|, v_E := E × B/|B|² where α ≥ ε and 0 where α < ε;
**Γ := C_S |S|²/ω + C_p γ_lin(M).**

**Definition of γ_lin(M).** At a point x, linearize (C1)–(C6) about M(x) with all coefficients and
background gradients ∇M(x) frozen at their values at x, ν_t, κ_t, D_t := 0, (C12) omitted, (C7) in
its Fourier form −|κ|²φ̂ = −iκ·Ê. For a plane wave ∝ exp(iκ·x + λt) this is a polynomial pencil in
λ with root set Λ(κ). γ_lin(M) := max( 0, sup{ Re λ : λ ∈ Λ(κ), 2π/h ≤ |κ| ≤ 1/ρ_L(M) } ), where
ρ_L(M) := (Σ_{s∈S_h⁺} m_s ρ_s v_{th,s}/|q_s|)/((Σ_{s∈S_h⁺} ρ_s)|B|), v_{th,s} := (T_h/m_s)^{1/2}
(ρ_L := +∞ where S_h⁺ has no mass or B = 0).

**(C10) Boundary conditions.** On Γ_in: the moments of g_s. On Γ_e: outflow (zero normal gradient
of all M; inflow forbidden). On Γ_w ∪ Γ_c: (C11).

**(C11) Material boundary** (applied at the wall-adjacent state at y := h/2). The α-weighted
combination of two laws:
(neutral law, weight 1 − α) zero normal mass flux; tangential shear stress and heat flux from the
all-y⁺ blended wall law: with y⁺ := y u_τ/ν_CE, u_τ is the solution of |u_t|/u_τ = y⁺ for y⁺ ≤ y⁺_c and
|u_t|/u_τ = κ_vK^{-1} ln(y⁺) + B_log for y⁺ > y⁺_c, y⁺_c the crossover where the two coincide;
τ_w = ρu_τ²; q_w = ρ c_{p,h} u_τ (T_h − T_w)/T⁺ with T⁺ := Pr y⁺ for y⁺ ≤ y⁺_c and
T⁺ := Pr_t[κ_vK^{-1} ln(y⁺) + B_log] + A_J P(Pr) otherwise, Pr := μ c_{p,h}/κ_h (the law reduces to the
laminar law as y⁺ → 0); ∂_n k = 0 and ω = max(u_τ/(√c_μ κ_vK y), 6ν_CE/(0.075 y²));
(charged law, weight α) normal particle flux of each s ∈ S_h equal to Γ_s^w := ρ_s c_s^w/m_s,
c_s^w := ((T_e + T_h)/m_s)^{1/2}, attenuated by (1 + Kn_w)^{-1}, Kn_w := ℓ_mfp/y; wall energy flux
Σ_{s∈S_h} Γ_s^w [ γ_sh T_e + m_s (c_s^w)²/2 ], γ_sh := 2 + e_c(φ − φ_w)/T_e; normal current
J·n = e_c Σ_{s∈S_h} Z_s Γ_s^w − e_c Γ_e^w exp(−e_c(φ − φ_w)/T_e), Γ_e^w := n_e (T_e/2πm_e)^{1/2};
φ_w = V_j on Γ_{c,j}; on Γ_w, φ_w is the value making J·n = 0; on cells crossed by open field lines
the particle and energy fluxes of this law are replaced by the end-loss fluxes of (C12). Incoming
labels are returned as prescribed by K_w's conversion law (the source W_s) at T_w.

**(C12) Along-line kinetic closure.** For each open field line ℓ (meeting ∂Ω at both ends) with arc
length s, for the heavy labels and for e separately, let F_ℓ(s, v_∥, μ), μ := v_⊥²/(2|B|), be the
stationary solution of v_∥ ∂_s F − [ μ ∂_s|B| + (q/m) ∂_s φ ] ∂_{v_∥} F = C_ℓ[F], with C_ℓ the
one-dimensional reduction of C at the local moments, φ from (C7), and data on incoming
characteristics only equal to the local equilibrium with the moments of M at that end cell; trapped
orbits are determined by C_ℓ. q_∥^{kin} := ∫ (m(v_∥² + 2μ|B|)/2) v_∥ F dv_∥ dμ; the end-loss fluxes are
the outgoing fluxes of F at the ends. On cells crossed by open lines,
q_∥ := sign(q_∥^{kin}) |q_CE||q_∥^{kin}|/(|q_CE| + |q_∥^{kin}|) (zero if q_∥^{kin} = 0).

**Coarse initial state.** M(0) := Π_h(z₀), z₀ ∈ Z₀.

**Coarse performance.** With the coarse fluxes on Γ_e computed as Σ_s ρ_s u_s(u_s·ê_z) + (p + ⅔ρk) ê_z
(u_s := u + 𝒥_s/ρ_s) and T from the coarse field, **P_𝒞(G; θ, h) := p evaluated at the unique stationary
state M\* of §6(a)**.

**Reduction property.** Where B ≡ 0 and n_e ≡ 0, (C1)–(C12) reduce to the compressible reacting
Navier–Stokes–Fourier–Fick equations with a two-equation eddy-viscosity model and an all-y⁺ wall
law.

## 4. Timescales and the response operator (constructive, at scale h)

Fix G and μ_G (H1); H := L²(μ_G); V_h the span of the components of Π_h z (all cells) and their
same-cell pairwise products; 𝒫 the orthogonal projection onto V_h. The Mori projection gives the
exact identity dY/dt = AY + ∫₀ᵗ K(τ)Y(t − τ)dτ + ξ(t), Y := 𝒫z, with
K(τ) := ⟨ξ(τ)⊗ξ(0)⟩ ⟨Y⊗Y⟩^{-1}. Define:
- τ_c(G) := max_i ∫₀^∞ |C_{p_i}(τ)|/C_{p_i}(0) dτ (0 if all C_{p_i}(0) = 0);
- τ_M(G): the **relaxation time** of the coarse observables: the smallest τ such that, for every
  z₀ ∈ Z₀ and all t ≥ 0, Σ_j ‖⟨Π_h z(t; G, z₀)⟩_{[t, t+T_d]} − ⟨Π_h z⟩_{μ_G}‖_{L²} ≤ C₀ e^{−t/τ} ×
  (the same quantity at t = 0), with C₀ := 3 fixed (well-defined for a deterministic steady state,
  where the correlation times vanish, and for a mixing one, where it is the mixing time);
- τ_K(G) := ∫₀^∞ ‖K(τ)‖/‖K(0)‖ dτ (0 if K ≡ 0 or if ⟨Y⊗Y⟩ is singular, as for a steady state);
- the response operator R_G: for bounded axisymmetric fields (a, q) on Ω (body acceleration,
  heating rate per unit mass), perturb (F1) by −a·∇_v f_s and by q (m_s/T_h)(|v − u|²/3 − T_h/m_s) f_s/2
  on each heavy label; R_G maps (a, q) to the linear response of P_F, with output norm
  max_i |δP_{F,i}|/|P_{F,i}| and input norm sup_Ω ( |a| τ_M/U_G + |q| τ_M/(c_{p,h}T_h) );
- ℓ_mfp(G); ρ̂(G) := sup over {α ≥ ε} of ρ_L/(gradient length of M); Kn(G) := sup ℓ_mfp/(gradient
  length of M); 𝓀(G) as in §2.

## 5. The class 𝒢(ε, h; C_R, δ_t, δ_η, K₀)

Fix ε ∈ (0, 1), h > 0, and class constants C_R, δ_t, δ_η > 0, K₀ ∈ (0, 1). Set T₁ := τ_M log(1/ε),
T_d := max(τ_c/ε², τ_M). G belongs to the class iff:
- **H1.** For every z₀ ∈ Z₀, F has a unique global weak solution; there is a unique invariant
  probability measure μ_G reached from Z₀, in the sense that every component of Π_h z(t) averaged
  over [T₁, T₁+T_d] lies within 3Δ_F-scale spread of its μ_G-mean; τ_c ≤ τ_M; μ_G is linearly
  responsive; and Z₀ ⊂ { z : Σ_j ‖Π_h(z)_j − ⟨Π_h z⟩_{μ_G,j}‖_{L²}/‖⟨Π_h z⟩_{μ_G,j}‖_{L²} ≤ 1 } (initial
  states within a fixed relative neighbourhood of the operating state).
- **H2.** τ_K ≤ ε τ_M; Kn ≤ ε; ρ̂ ≤ ε; ℓ_mfp/L_G ≤ ε²; ℓ_mfp/ε ≤ h ≤ ε L_G; reaction zones (where a
  conversion rate exceeds ε times its maximum) have thickness ≥ h/ε; products of conversions with
  Q > 0 have gyroradius and slowing-down length ≤ ε times the local gradient length; every open
  field line of the mean field is crossed by the mean flow in time ≥ τ_M/ε.
  **Local equilibrium is reachable (added in version 5 after a counterexample):**
  (i) *collisions per transit:* with ν_coll(x) the collision frequency of the dominant heavy label
  under μ_G and the flight length ℓ_flight(x) := (|ū| + v_th)/ν_coll, require
  sup_Ω ℓ_flight/(gradient length of M) ≤ ε and ℓ_flight/L_G ≤ ε² (Kn is measured with the flow
  speed included, so a fast flow cannot transit without relaxing);
  (ii) *inlets near equilibrium:* for every s, the inflow datum g_s is within ε of a local
  Maxwell–Boltzmann state with the same (mass, momentum, energy) moments, in relative L¹ over
  velocity, pointwise on Γ_in;
  (iii) *losses slower than collisions:* for every s and every f in the support of μ_G,
  L_s[f](x, v) ≤ ε ν_coll(x) f_s(x, v) (radiative and escape losses per collision time are small,
  so the loss operator cannot cancel the collisional relaxation the transport coefficients
  presuppose).
- **H3.** ‖R_G‖ ≤ C_R.
- **H4.** Under μ_G: Φ_m > 0, P_in > 0, Φ_p/(g₀Φ_m) ≥ δ_t, Φ_p²/(2Φ_m P_in) ∈ [δ_η, 1].
- **H5.** 𝓀(G) ≤ K₀ and sup_Ω ν_t/ν_CE ≤ K₀/ε under μ_G with ω as defined in §2 (the modelled turbulent
  transport is bounded relative to the collisional transport).

**Witness of nonemptiness (to be checked by the solver, not assumed).** A straight cylinder of
length 20 L with a uniform axial field B_z, a single neutral heavy label with a hydrogen-like
collision kernel, inflow at one end at Mach 0.3 and Knudsen 10^{-3}, outflow at the other, walls at
the inflow temperature, no reactions, no conductors, is expected to satisfy H1–H5 with 𝓀 = O(Kn)
for all sufficiently small ε. If the solver finds no witness, §8 applies.

## 6. The statement

**(T)** For every choice of class constants (C_R, δ_t, δ_η, K₀) there exist θ\* ∈ Θ, ε_⋆ > 0,
constants c₁, c₂ > 0, and exponents β, γ > 0, all depending only on the class constants (not on
G, ε, or h), such that for every ε ∈ (0, ε_⋆), every admissible h, and every G in the class:

(a) **Convergence.** For every z₀ ∈ Z₀ the coarse system 𝒞_{θ\*} at scale h from M(0) = Π_h(z₀)
converges to one stationary state M\* independent of z₀, which is nonlinearly stable with
exponential return rate at least 1/(c₁τ_M) in the norm Σ_j ‖M_j − M\*_j‖_{L²}/‖M\*_j‖_{L²} (over
components with nonzero norm), and is reached to distance ε in that norm within c₁τ_M log(1/ε).

(b) **Fidelity within intrinsic variability.** For each i ∈ {1, 2, 3},
|P_{𝒞,i}(G; θ\*, h) − P_{F,i}(G)| ≤ c₂ [ ε^β + 𝓀(G)^γ + Δ_{F,i}(G) ] |P_{F,i}(G)|.

## 7. The exponents

Define the sharp exponents
β\* := sup{ β : (a)–(b) hold for some θ\*, c₁, c₂, γ > 0 },
γ\* := sup{ γ : (a)–(b) hold for some θ\*, c₁, c₂, β > 0 },
and the **floor**, for κ₀ ∈ (0, K₀]:
ε_floor(ε, h, κ₀) := inf_{θ∈Θ} sup{ max_i |P_{𝒞,i} − P_{F,i}|/|P_{F,i}| : G in the class with 𝓀(G) ≥ κ₀ },
where for a θ at which (a) fails on some such G the inner supremum is +∞.
The matching lower-bound exponent is
γ_\* := inf{ γ' : ε_floor(ε, h, κ₀) ≥ c κ₀^{γ'} for some c > 0 and all small ε, admissible h, κ₀ }.

**(T-exp)** Determine β\* and γ\*, and decide whether γ_\* = γ\* (the 𝓀-term of (b) is sharp: no
closure of this form removes it).

## 8. The question

Decide (T). If true, prove it and determine (β\*, γ\*) and whether γ_\* = γ\*. A proof must exhibit
θ\* explicitly or give a finite procedure that computes it from the kernels of F (C, σ, Q, L, K_w)
and the class constants; a pure existence proof answers (T) but does not deliver the coarse system
as an instrument, and is to be labelled as such. If false, exhibit a class member G and prove that
no θ ∈ Θ satisfies (a)–(b) with any constants and exponents. If you
conclude the class is empty for all small ε under the definitions above, prove that, identify the
weakest relaxation of H1–H5 under which it is nonempty, and answer for that class; emptiness alone
is not an answer.

## 9. Remarks

**(R1) The two small parameters.** ε measures the distance of the fine dynamics from a smooth,
collisional, magnetized, short-memory fluid; 𝓀 measures the kinetic energy the cells cannot see.
Both are computed from the coarse state itself. (b) says the coarse description's relative error is
controlled by a power of each, plus the system's own statistical spread. β measures how fast the
kinetic and memory closures converge; γ measures how fast the error grows with sub-cell
fluctuation; γ_\* = γ\* says the growth is a property of the closure class, not of a poor choice of
constants.

**(R2) Expected mechanisms (non-binding).** Every non-Chapman–Enskog term of 𝒞 enters at first
order in ε (memory, gyroviscous anisotropy, the sub-cell normal stress), so β\* = 1 is the natural
expectation. The true sub-cell momentum flux is a tensor ⟨ρu′⊗u′⟩ whose split and alignment the
model fixes by ⅔ρk I + 2ρν_t P_⊥SP_⊥; two class members with identical M but different sub-cell
anisotropy differ in exit momentum flux by an O(1) fraction of 𝓀 ρU_G², suggesting γ\* = γ_\* = 1.

**(R3) The single external clause.** Everything above is a statement about two explicit systems of
equations. Its bearing on a physical device rests on one clause outside mathematics: that F, with
its given kernels, is the true dynamics of that device.

**(R4) Excluded by construction.** Field nulls where α ≥ ε; reaction products with orbits comparable
to the device; sub-grid reaction fronts; time-dependent external sources.

## 10. Defects found in version 3 by the third cold review, and the repairs

| Defect | Repair |
|---|---|
| the weakly fluctuating class was empty (fluctuation intensity normalized by the local flow speed, zero at walls and stagnation points), so (T_w) was vacuously true | 𝓀 defined on interior cells with |ū_c| ≥ δ_u U_G, normalized by U_G²; a nonemptiness witness given |
| the strongly fluctuating class had no upper bound on 𝓀 (worst-case error infinite) and was empty at small h; sup ∅ := −∞ made the floor −∞ | H5 bounds 𝓀 ≤ K₀; sup ∅ := 0; h quantified with the class |
| the fine-side relaxation time τ_M was used inside the coarse definition of ω (a dependence of 𝒞 on G outside M); τ_M was a correlation time, undefined for a steady state | ω floored at the coarse viscous rate ν_CE/h²; τ_M redefined as the relaxation time from Z₀ (version 4 self-review) |
| the exponents of (T_s) were not sharp ((0, 10) was a correct answer) | β\*, γ\* defined as suprema; γ_\* as the matching infimum; the question is whether γ_\* = γ\* |
| H3's input norm used the local flow speed, vacuous where ū = 0 | normalized by U_G |
| the wall law had frozen turbulent-wall constants and no laminar limit; ν_t could blow up as ω → 0 | all-y⁺ blended wall law reducing to the laminar law; its constants moved into θ; ω floored at 1/τ_M; H5 bounds ν_t/ν_CE |
| Z₀ unconstrained (a large sub-cell reservoir defeats any c₁); P_𝒞 undefined if starts settle differently | Z₀ within a fixed relative neighbourhood of the operating state; (a) asserts a z₀-independent M\*; P_𝒞 defined there |
| C_R, δ_t, δ_η, ε_⋆, and H1's constant chosen by the solver | class constants universally quantified with the existential constants depending only on them; H1's constant set to 3 |
| "understood to be nonempty" stipulated what was provably false | stipulation removed; witness supplied; §8 governs emptiness |
| two statements with three exponents, one misnamed | one statement (T), one inequality, two exponents (β, γ) each attached to its own small parameter; the floor as the matching lower bound |
| **(version 4 → 5, from the solver's counterexample)** the hypotheses permitted a fine system that never reaches local equilibrium: a loss operator equal to the collisional gain, a hypersonic transit with almost no collisions under a thermal Knudsen convention, and non-Maxwellian inlets seen by the coarse system only through their moments | H2 (i) Knudsen number measured with the flight length (flow speed included) and collisions per transit ≥ 1/ε; (ii) inlets within ε of local Maxwellians; (iii) losses per collision time ≤ ε. Each corresponds to a physical exclusion the coarse instrument must enforce: collisionless regions, beam-like injection, loss-dominated media |
