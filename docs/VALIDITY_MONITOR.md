# Validity monitor: design note (30 September 2026, no code yet)

Purpose: mark, cell by cell and step by step, where the plasma model in use stops being
predictive, and turn that into a statement attached to every reported number ("this impulse
has X% of its force history in cells outside the model"). This follows the outside review
(`docs/ASTRA_PLASMA_MODEL_2026-09-30.md`, opinion, not canon), which asks for exactly this
end-condition. The monitor does not change the solution. It only reads it.

The model it guards is the one planned next: single-fluid resistive MHD, ideal-gas or tabulated
EOS, optional second temperature, no Hall term, no electron pressure in Ohm's law, axisymmetric.

## 1. Inputs per cell

From the state: n_i, n_e = Z n_i, T_e, T_i (equal if one-temperature), u, B (total field,
b0 + b1 in the split form), J = curl B / mu0 from the same discrete circulation the solver uses
(`jTheta` in the spike), ion mass m_i, charge state Z (from the EOS/ionisation model, or stated
if fixed).

Gradient length, per cell: L = min over q in {n_e, T_e, p, |B|} of |q| / |grad q|, with grad q
from face differences (not limited slopes, so a sharp interface such as a diamagnetic cavity
edge gives a short L instead of being hidden by the limiter). Also kept: the cell size h, because
L < 2h means the gradient is not resolved and every ratio below is then a lower bound.

## 2. Quantities, formulas and thresholds

Formulas are from the NRL Plasma Formulary (2019 edition) in SI. Coulomb logarithm: the same
one the Spitzer conductivity uses, stated once in code and in the report: electron-ion
ln Lambda = 24 - ln(sqrt(n_e[cm^-3]) / T_e[eV]) for T_e > 10 Z^2 eV, and
23 - ln(sqrt(n_e[cm^-3]) Z T_e[eV]^-1.5) below that, floored at 2 (a floor of 2 is itself a flag:
the plasma is then strongly coupled and Spitzer does not hold).

The thresholds below are my defaults, chosen as "an order of magnitude of margin". They are
engineering parameters, set in one place, and each report prints the values used.

| # | Check | Ratio | Marginal | Outside |
|---|---|---|---|---|
| 1 | Collisionality (fluid closure) | Kn_e = lambda_e / L, Kn_i = lambda_i / L, mean free paths lambda = v_th tau | > 0.03 | > 0.3 |
| 2 | Local thermal equilibrium | tau_e / t_flow, tau_i / t_flow, with t_flow = L / max(\|u\|, c_s) | > 0.1 | > 1 |
| 3 | Two-temperature need | tau_eq,ei / t_flow (electron-ion equilibration) | > 0.1 with one T | n/a with two T |
| 4 | Ion orbit size | rho_i / L (ion gyroradius) | > 0.1 | > 1 |
| 5 | Hall scale | d_i / L, d_i = c / omega_pi | > 0.1 | > 1 |
| 6 | Omitted Ohm-law terms (below) | R_Hall, R_pe, R_me | > 0.1 | > 1 |
| 7 | Resolution | L / h | < 4 | < 2 |
| 8 | Vacuum and floors | cell at density or pressure floor, or dual-energy pressure in use | any | at floor |

Magnetisation (omega_c tau) is reported but is not a validity failure by itself: when
omega_ce tau_e > 1 the heat flux is anisotropic. If the model has anisotropic transport this is
fine. If not, the cell is flagged under row 2 as "closure missing", not "model invalid".

The magnetic Reynolds number Rm = mu0 sigma u L is reported, not flagged: a resistive model is
valid at low Rm. It tells the reader whether the answer is ideal-like or diffusion-dominated.

## 3. Omitted Ohm-law terms, estimated properly

Generalised Ohm's law: E + u x B = eta J + (J x B - grad p_e) / (n_e e) + (m_e / (n_e e^2)) dJ/dt.
The model keeps the first two terms on each side. Comparing the omitted electric fields with
u x B directly is misleading in RZ. Our state has only J_theta (no B_theta), so the Hall field
(J x B)/(n e) is purely poloidal and does nothing to the theta component. What matters is the
field it would generate: its curl drives a B_theta that the model never creates. So the
monitor compares rates of change of B, not electric fields:

- R_Hall = |curl( (J x B) / (n_e e) )| / max(|curl(u x B)|, |curl(eta J)|)
- R_pe = |curl( grad p_e / (n_e e) )| / same denominator. This is the Biermann battery,
  |grad n_e x grad T_e| / (n_e e). It is non-zero wherever density and temperature gradients
  are not parallel, which is generic in a laser-ablation plume.
- R_me = (m_e / (n_e e^2)) |J| / t_flow / |E_retained| (electron inertia; normally tiny,
  kept for completeness).

All three use the discrete curl of the cell-centred vectors with the same stencil, so a
smooth field gives the correct ratio and a one-cell spike shows up as a resolution flag
(row 7) rather than a false physics flag.

## 4. Turning flags into a statement about the answer

A tenuous, dynamically irrelevant cell outside the model should not condemn the run. A
flagged cell carrying the force should. So the monitor weights flags by what the answer
depends on.

For each reported integral output Q (impulse, coil reaction force, energy partition), the solver
already books Q as a sum of cell or face contributions (for example bodyAxialForce is a sum over
cells). The monitor accumulates, over time, the part of |dQ| that came from cells in each flag
class:

  share_outside(Q) = integral of sum over outside cells |dQ_cell| dt / integral of sum over all cells |dQ_cell| dt,

and the same for marginal. Absolute values are used so that cancelling contributions cannot hide
a flagged region.

Report rule (default):

- share_outside < 1%, share_marginal < 10%: predictive, report with the normal band.
- share_outside between 1% and 10%: report, with the flagged share stated next to the number
  and widened band (the flagged share is added to the uncertainty in full, both directions).
- share_outside > 10%: the number is not reported as a prediction. The run reports where
  (cells, z-r map), when (first time over threshold) and which check failed, plus the
  momentum and energy that passed through the flagged region.

The same accounting runs for energy: Joule heat, sync energy from the dual-energy
correction, and floor energy are separate lines in the energy budget, each with its share in
flagged cells.

## 5. What it cannot see (stated in every report)

- Axisymmetry: no azimuthal modes. Rayleigh-Taylor or flute modes at a decelerating
  plasma-field interface, and lower-hybrid drift at a thin current sheet, are invisible. The
  monitor estimates the RT e-folding count, integral of sqrt(k a_eff) dt with
  a_eff = interface deceleration and k = 1/h (the fastest mode the mesh could carry), and flags
  the interface when it exceeds 3. That says "a 3-D instability may have time to grow", not what
  it does.
- Kinetic effects that a local fluid check cannot detect: non-local electron heat flux from a
  hot region into a cold one (flag 1 catches the local part only), trapped or escaping
  electrons, beam-like ion populations in the far plume.
- Errors in the inputs themselves (EOS, ionisation, conductivity tables). These go into the
  pedigree and the input bands, not the monitor.

## 6. Order of work

1. Compute and dump the per-cell ratios (rows 1 to 8) in the spike, viewable as maps. Needs only
   the state and the conductivity model from task 2.
2. Add the per-output share accounting to the force and energy budgets.
3. Only then wire the report rule into the result format.

Nothing here needs Ben unless the thresholds above are judged too loose or too tight for the
first validation case. They are parameters, recorded with every result.
