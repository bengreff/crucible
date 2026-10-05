// Verification of the PaSR turbulence-chemistry closure (TECHNICAL_PLAN step 7, "PaSR"). Criteria
// 1 to 3 stated 4 October 2026, before the closure was implemented. Mechanism h2o2.yaml; the
// chemical-time species S = {H2, O2, H2O}; Sc_t 0.7 and beta* 0.09 as in the SST model.
// The closure, per cell and reaction substep:
//   kappa = tau_c / (tau_c + tau_mix),  kappa_eff = 1 - s (1 - kappa),
//   tau_c = the largest Y_i / |dY_i/dt| over the species of S with Y_i > 0 and dY_i/dt != 0
//           (kappa = 1 if there is none), evaluated from z at every right-hand-side call;
//   tau_mix = C_mix sqrt(nu_eff / epsilon), nu_eff = (mu + mu_t) / rho, epsilon = beta* k omega;
//   s = the largest over the species of S with 0 < X_i < 1 of
//       min(1, (nu_t / Sc_t) |grad X_i|^2 / (beta* omega X_i (1 - X_i)));
//   CVODES integrates dz/dt = kappa_eff(z) f(z), with tau_mix and s frozen over the substep.
//   1. The closure's factor and inputs.
//      (a) kappa_eff from the engine against the formula above evaluated here from
//          ReactionSource::rates, at three states (stoichiometric H2/O2 unburnt at 1200 K and
//          101325 Pa; the same at fixed rho after a laminar substep that takes it through half its
//          temperature rise; its constant-(u, v) equilibrium) and three closures (tau_mix, s) =
//          (1e-7 s, 1), (1e-6 s, 0.5), (1e-5 s, 0.2): within 1e-13 relative. Pure N2 (every species
//          of S dormant): exactly 1.
//      (b) The cell inputs from the flow. A straight duct (R 5 mm, L 40 mm, nz 16, nr 4, adiabatic
//          slip walls, so the wall distance is infinite) at rest, 300 K and 101325 Pa, with
//          X_H2 = 0.1 + 20 z/m, X_O2 = 0.85 - 20 z/m, X_N2 = 0.05, uniform omega = 1e3 1/s and
//          k = 10 m^2/s^2 (s about 0.26 to 0.58 off the ends, derived), then k = 25 (s clipped to 1
//          in the outer cells), C_mix = 1. At rest there is no strain, so mu_t = rho k / omega, and
//          the least-squares gradients are exact for a field linear in z away from the open ends.
//          In every cell off the two end columns: tau_mix within 1e-12 relative of
//          C_mix sqrt(((mu + rho k / omega) / rho) / (beta* k omega)), with mu from the medium's
//          transport at the cell state, and s within 1e-12 relative of the formula with
//          |grad X_i| = 20 /m. With k = 0 in every cell: s = 0 exactly.
//      (c) ReactingFlow rejects the closure with LocalEquilibrium or Frozen chemistry and on a flow
//          without turbulence; FiniteRate with turbulence is accepted.
//   2. One substep against an independent integration. Stoichiometric H2/O2 at 1200 K and
//      101325 Pa, at fixed rho and e, over one substep of 3e-4 s with (tau_mix, s) = (1e-5 s, 1),
//      (1e-4 s, 0.5) and the laminar s = 0. Reference: classical RK4 written here, on
//      ReactionSource::rates and this file's own kappa_eff, with the step halved until T and every
//      Y_k change by less than 1e-10 (relative for T, absolute for Y) between successive halvings.
//      At 20 equally spaced times (the observer of one uninterrupted engine integration), with
//      rtol 1e-10 and atol 1e-16: T within 1e-7 relative and every Y_k within 1e-8 absolute of the
//      reference, and the largest error at rtol 1e-10 below the largest at rtol 1e-7 (the
//      difference is the integrator's). The closure must change T by more than 1e-2 relative at
//      some sample against the laminar reference, or the case does not test it.
//   3. Limits.
//      (a) s = 0: the substep equals the substep without the closure bit for bit, at the three
//          states of 1(a).
//      (b) Fast mixing: with s = 1 and tau_mix = 1e-8, 1e-9 and 1e-10 s, the largest relative T
//          difference from the laminar substep over the samples of 2 falls with an observed order
//          in tau_mix of at least 0.8 between the two smallest.
//      (c) Inert: pure N2 at 1500 K through a substep of 1e-4 s with (1e-6 s, 1): z unchanged exactly.
//      (d) Reaction-free mixing: the duct of 1(b) with N2 and AR only (both inert in h2o2.yaml;
//          X_AR = 0.1 + 20 z/m), k = 10. One ReactingFlow step with the closure against one Flow step
//          of the same length from the same state: every conserved field within 1e-14 relative.
// Criterion 4, a turbulent reacting chamber with the closure (budgets, positivity, the kappa_eff and
// s fields reported against the same run without it), is stated before its first run.
