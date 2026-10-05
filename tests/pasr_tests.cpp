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
//          Amended 5 October 2026 after run 2 (before it, run 1 stopped at a duct without ambient
//          turbulence values; the duct now takes the uniform field as ambient). Run 2 measured
//          gaps of 1e-15 in axial and 8e-15 in radial momentum, against field maxima of 1.5e-4 and
//          1.0e-14: at rest both momenta are round-off, so "relative to the field's max" is
//          ill-posed for them. ReactingFlow without the closure gave the same gaps (measured), so
//          they come from the reaction substep rewriting Y at the 1e-16 level, not from the
//          closure. Criterion now: (i) the step with the closure equals the same ReactingFlow
//          step without it bit for bit; (ii) against the Flow step, density, total energy,
//          partial densities and rho k, rho omega within 1e-14 relative to the field's max, and
//          both momenta within 1e-14 of sqrt(rho_max p0) (the momentum scale of a pressure wave
//          of amplitude p0, up to sqrt(gamma)).
//      The reference of criterion 2 halves the RK4 step up to 2^20 steps per sample (2^17 in
//      run 2, where classical RK4 was unstable up to 2^17 and stable from 2^18 on this
//      mechanism, measured; the gap between 2^18 and 2^19 was 9e-13 in T and 7e-13 in Y).
//      Criterion 2 amended 5 October 2026 after run 3: the case (1e-4 s, 0.5) passed its accuracy
//      checks (T 4.1e-10, Y 2.4e-10) but changed T by only 1.5e-3 against laminar, below the
//      1e-2 that shows a case tests the closure, so it failed as stated. Mechanism (inferred):
//      with s = 0.5, kappa_eff >= 0.5, so the closure at most halves the rate, and only in the
//      heat release (tau_c << tau_mix), which lasts microseconds between samples 15 us apart;
//      in the induction tau_c >> tau_mix and kappa_eff is near 1. It is replaced by
//      (1e-5 s, 0.9), where kappa_eff can fall to 0.1 (stated before its run: the T change is
//      predicted above 1e-2, from the 0.28 measured with (1e-5 s, 1) in run 3).
//      Run 4: the prediction failed, 3.5e-3 (accuracy passed, T 5.4e-10, Y 3.1e-10). Measured
//      mechanism (docs/evidence/pasr/heat_release_2026-10-05.txt): every case ignites near 18 us,
//      and the closure stretches the 10-90% rise from 3.6 us (laminar) to 7.8 us (s 0.9) and
//      12.9 us (s 1). The check passes only if the stretched rise still runs at the 30 us
//      sample, so it measures where the samples fall, not the closure. Stopped after two
//      attempts at the fractional case; this check stays failing until it is restated.
//      Restated 5 October 2026 (03:15), before run 5, which is run once: the closure's effect is
//      read from a fine record of the two converged references, T at 1280 equal times (every
//      1/64 of a sample, 0.23 us apart), which resolves the 3.6 to 12.9 us rise wherever it
//      falls. Check, for each case with s > 0: the largest relative T difference between the
//      reference with the closure and the laminar reference over the fine record exceeds 1e-2.
//      The cases, the 20 samples and the accuracy checks are unchanged. Predicted for
//      (1e-5 s, 0.9), derived from the rise table: when the laminar rise ends, the stretched
//      rise is about half done, so the gap is a large part of the temperature rise, order 0.1.
//   4. A turbulent reacting chamber with the closure. Stated 5 October 2026, before any run of it;
//      it runs in a chamber study, not in this file. The case is the C1 chamber of
//      tests/chamber_study.cpp: contour, premixed H2/O2 at O/F 5 and 0.4 kg/s, valve ramp, igniter,
//      ambient N2 at 1 kPa and 300 K. Chemistry is FiniteRate (rtol 1e-6, atol 1e-12). It is made
//      viscous and turbulent as in tests/step_cost.cpp: mixture-averaged transport, SST-2003,
//      no-slip walls at 600 K, supply turbulence I 0.05 with viscosity ratio 10, and a
//      Spalart-Rumsey ambient. The closure is on (C_mix 1, S = {H2, O2, H2O}); the control is
//      the same run without it. Grids 32x6 and 64x12 (uniform rings, so the wall is not
//      resolved; the first-cell y+ is reported), each to 8 ms, C1's settling time.
//      Judged, for each run:
//      (a) Budgets over the whole run. Mass below 1e-11 of the initial fill, as C1. Energy below
//          1e-11 of the magnitude of the gas's total energy at the end of the run; normalizing by
//          the initial N2 fill failed C1 for a reason that is not a leak (CHAMBER_C1.md). Both
//          normalizations are reported.
//      (b) Positivity: the run reaches 8 ms without an exception, and rho, p, T, k and omega are
//          positive in every cell at every history sample (2 us apart).
//      (c) Settled, as C1 criterion 1: relative drift of outlet mass flow, injector pressure and
//          vacuum thrust below 1e-3 over the last 1 ms, and outlet mass flow within 1e-3 of the
//          supply. A run that fails only this is extended once to 12 ms and reported.
//      (d) The closure's fields are admissible: kappa_eff in (0, 1] and s in [0, 1] in every cell
//          at the field snapshots.
//      Reported, not judged: the kappa_eff and s fields at light-off and at the end; the
//      fraction of the chamber volume with s > 0.01; the largest mass fraction clipped to zero
//      after a reaction substep; and, against the control, light-off time, settled injector
//      pressure, c* and vacuum Isp.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/flow.hpp"

using namespace crucible;
using thermo::Mixing;

namespace {

int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-66s %.3e (limit %.1e)  %s\n", what, value, limit, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

constexpr double betaStar = 0.09, schmidt = 0.7, p0 = 101325;
double sq(double x) { return x * x; }

// This file's own kappa_eff (the header's formula), from z = [T, Y] and its laminar rate.
double kappaEff(const double* z, const double* rate, const std::vector<std::size_t>& set, double tauMix, double s) {
  double tauC = 0;
  bool any = false;
  for (std::size_t i : set)
    if (z[i + 1] > 0 && rate[i + 1] != 0) { tauC = std::max(tauC, z[i + 1] / std::abs(rate[i + 1])); any = true; }
  if (!any) return 1;
  return 1 - s * (1 - tauC / (tauC + tauMix));
}

// Classical RK4 on dz/dt = kappaEff(z) f(z) at fixed rho over 20 equal samples of dt, with perSample
// steps per sample; empty if the integration failed (an unstable step size). If record is given,
// T is also kept at 64 equal times per sample (perSample is a multiple of 64).
std::vector<std::vector<double>> rk4(thermo::ReactionSource& source, double rho, std::vector<double> z, double dt,
                                     long perSample, const std::vector<std::size_t>& set, double tauMix, double s,
                                     std::vector<double>* record = nullptr) {
  const std::size_t n = z.size();
  std::vector<double> k1(n), k2(n), k3(n), k4(n), w(n);
  auto f = [&](const std::vector<double>& x, std::vector<double>& out) {
    source.rates(rho, x.data(), out.data());
    const double factor = kappaEff(x.data(), out.data(), set, tauMix, s);
    for (double& v : out) v *= factor;
  };
  const double h = dt / (20.0 * double(perSample));
  std::vector<std::vector<double>> samples;
  if (record) record->clear();
  try {
    for (int sample = 0; sample < 20; ++sample) {
      for (long step = 0; step < perSample; ++step) {
        f(z, k1);
        for (std::size_t i = 0; i < n; ++i) w[i] = z[i] + 0.5 * h * k1[i];
        f(w, k2);
        for (std::size_t i = 0; i < n; ++i) w[i] = z[i] + 0.5 * h * k2[i];
        f(w, k3);
        for (std::size_t i = 0; i < n; ++i) w[i] = z[i] + h * k3[i];
        f(w, k4);
        for (std::size_t i = 0; i < n; ++i) z[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        if (!std::isfinite(z[0]) || z[0] <= 0 || z[0] > 1e4) return {};
        if (record && (step + 1) % (perSample / 64) == 0) record->push_back(z[0]);
      }
      samples.push_back(z);
    }
  } catch (const std::exception&) {
    return {};
  }
  return samples;
}

// Largest T difference (relative) and Y difference (absolute) between two sample sets.
struct Difference { double t = 0, y = 0; };
Difference difference(const std::vector<std::vector<double>>& a, const std::vector<std::vector<double>>& b) {
  Difference d;
  for (std::size_t s = 0; s < a.size(); ++s) {
    d.t = std::max(d.t, std::abs(a[s][0] - b[s][0]) / b[s][0]);
    for (std::size_t i = 1; i < a[s].size(); ++i) d.y = std::max(d.y, std::abs(a[s][i] - b[s][i]));
  }
  return d;
}

// The reference of criterion 2: RK4 with the step halved until T and every Y change by less than
// 1e-10 between successive halvings; the finer of the last pair, with its fine T record.
std::vector<std::vector<double>> reference(thermo::ReactionSource& source, double rho, const std::vector<double>& z,
                                           double dt, const std::vector<std::size_t>& set, double tauMix, double s,
                                           long& perSample, std::vector<double>& record) {
  std::vector<std::vector<double>> coarse = rk4(source, rho, z, dt, 256, set, tauMix, s);
  for (perSample = 512; perSample <= (1L << 20); perSample *= 2) {
    auto fine = rk4(source, rho, z, dt, perSample, set, tauMix, s, &record);
    if (!fine.empty() && !coarse.empty()) {
      const auto d = difference(coarse, fine);
      if (d.t < 1e-10 && d.y < 1e-10) return fine;
    }
    coarse = std::move(fine);
  }
  return {};
}

// The engine's substep with the closure, sampled at 20 equal times.
std::vector<std::vector<double>> engine(thermo::ReactionSource& source, double rho, std::vector<double> z, double dt,
                                        const Mixing& mixing, const std::vector<std::size_t>& set, double rtol) {
  thermo::ReactionStep step(source, rtol, 1e-16);
  std::vector<std::vector<double>> samples;
  step.advance(rho, z.data(), dt, mixing, set, 20,
               [&](double, const double* x) { samples.emplace_back(x, x + z.size()); });
  return samples;
}

struct Duct {
  Definition definition;
  std::vector<Primitive> cells;
  std::vector<double> y;
};

// The straight duct of 1(b) and 3(d) at rest, 300 K and p0: R 5 mm, L 40 mm, nz 16, nr 4, adiabatic
// slip walls. moles(z, x) writes the mole fractions at axial position z.
template <class Moles>
Duct duct(const Medium& medium, const TransportFits& fits, bool turbulence, Moles moles) {
  Duct d;
  auto& def = d.definition;
  def.experiment = Case::UniformDuct;
  def.nz = 16; def.nr = 4; def.length = 0.04; def.inletRadius = def.exitRadius = def.throatRadius = 0.005;
  def.species = medium.species();
  const std::size_t n = medium.size();
  def.composition.assign(n, 0.0);
  for (std::size_t k = 0; k < n; ++k) def.composition[k] = medium.species()[k].name == "N2";
  def.totalPressure = def.backPressure = p0;
  def.transport = fits;
  def.wallSlip = true; def.wallTemperature = 0;
  def.turbulence.enabled = turbulence;
  def.turbulence.ambientK = 10; def.turbulence.ambientOmega = 1e3;  // the uniform field of 3(d)
  Mesh m(def);
  d.cells.resize(m.cells.size());
  d.y.resize(m.cells.size() * n);
  std::vector<double> x(n);
  for (std::size_t q = 0; q < m.cells.size(); ++q) {
    std::fill(x.begin(), x.end(), 0.0);
    moles(m.cells[q].z, x.data());
    double mass = 0;
    for (std::size_t k = 0; k < n; ++k) mass += x[k] * medium.species()[k].molarMass;
    double* y = d.y.data() + q * n;
    for (std::size_t k = 0; k < n; ++k) y[k] = x[k] * medium.species()[k].molarMass / mass;
    d.cells[q] = {p0 / (medium.gasConstant(y) * 300), 0, 0, p0};
  }
  return d;
}

std::vector<double> uniformTurbulence(std::size_t cells, double k, double omega) {
  std::vector<double> kOmega(2 * cells);
  for (std::size_t q = 0; q < cells; ++q) { kOmega[2 * q] = k; kOmega[2 * q + 1] = omega; }
  return kOmega;
}

}  // namespace

int main() {
  const std::string mechanism = "h2o2.yaml";
  thermo::ReactionSource source(mechanism);
  Medium medium = source.medium();
  const auto fits = thermo::transportFits(mechanism);
  medium.setTransport(fits);
  const std::size_t ns = medium.size();
  auto index = [&](const char* name) {
    for (std::size_t k = 0; k < ns; ++k) if (medium.species()[k].name == name) return k;
    throw std::runtime_error(std::string("missing species ") + name);
  };
  const std::size_t iH2 = index("H2"), iO2 = index("O2"), iH2O = index("H2O"), iN2 = index("N2"), iAR = index("AR");
  const std::vector<std::size_t> set = {iH2, iO2, iH2O};
  const std::vector<std::string> setNames = {"H2", "O2", "H2O"};

  // The three states of 1(a) and 3(a), at one density.
  std::vector<double> unburnt(ns + 1);
  {
    const auto y = source.massFractions("H2:2, O2:1");
    unburnt[0] = 1200;
    std::copy(y.begin(), y.end(), unburnt.begin() + 1);
  }
  const double rho = source.density(1200, p0, unburnt.data() + 1);
  std::vector<double> equilibrium = unburnt;
  source.equilibrateUV(rho, equilibrium.data());
  std::vector<double> half = unburnt;
  {
    thermo::ReactionStep laminar(source, 1e-10, 1e-16);
    const double middle = 0.5 * (unburnt[0] + equilibrium[0]);
    for (int i = 0; i < 100000 && half[0] < middle; ++i) laminar.advance(rho, half.data(), 1e-7);
  }
  std::printf("states: rho %.6f kg/m^3; T %.2f K (unburnt), %.2f K (half rise), %.2f K (equilibrium)\n", rho,
              unburnt[0], half[0], equilibrium[0]);
  const std::vector<std::vector<double>> states = {unburnt, half, equilibrium};

  std::printf("1(a) kappa_eff against the formula\n");
  {
    double worst = 0;
    std::vector<double> rate(ns + 1);
    for (const auto& z : states)
      for (const auto& [tauMix, s] : std::vector<std::pair<double, double>>{{1e-7, 1}, {1e-6, 0.5}, {1e-5, 0.2}}) {
        source.rates(rho, z.data(), rate.data());
        const double mine = kappaEff(z.data(), rate.data(), set, tauMix, s);
        const double theirs = thermo::reactingFraction(z.data(), rate.data(), set, {tauMix, s});
        worst = std::max(worst, std::abs(theirs - mine) / mine);
      }
    check(worst < 1e-13, "3 states x 3 closures: |engine - formula| / formula", worst, 1e-13);
    std::vector<double> n2(ns + 1, 0.0);
    n2[0] = 1500; n2[1 + iN2] = 1;
    source.rates(rho, n2.data(), rate.data());
    const double f = thermo::reactingFraction(n2.data(), rate.data(), set, {1e-7, 1});
    check(f == 1, "pure N2: kappa_eff - 1 (exactly 0)", std::abs(f - 1), 0);
  }

  std::printf("1(b) cell inputs on the duct\n");
  auto linearH2O2 = [&](double z, double* x) { x[iH2] = 0.1 + 20 * z; x[iO2] = 0.85 - 20 * z; x[iN2] = 0.05; };
  for (double k : {10.0, 25.0, 0.0}) {
    const double omega = 1e3;
    Duct d = duct(medium, fits, true, linearH2O2);
    Flow flow(d.definition);
    flow.setInitialState(d.cells, d.y);
    flow.setTurbulence(uniformTurbulence(d.cells.size(), k, omega));
    std::vector<double> time, segregation;
    flow.mixingInputs(1, set, time, segregation);
    const auto& m = flow.mesh();
    if (k == 0) {
      const bool zero = std::all_of(segregation.begin(), segregation.end(), [](double v) { return v == 0; });
      check(zero, "k = 0: s exactly 0 in every cell (1 = no)", zero ? 0 : 1, 0);
      continue;
    }
    double worstTime = 0, worstS = 0, sMin = 1, sMax = 0;
    std::vector<double> diffusion(ns), work;
    for (int i = 1; i < m.nz - 1; ++i)
      for (int j = 0; j < m.nr; ++j) {
        const auto q = m.index(i, j);
        const auto w = flow.cellPrimitive(q);
        const auto y = flow.massFractions(q);
        const double mu = medium.transport(flow.temperature(q), w.p, y.data(), diffusion.data(), work).viscosity;
        const double tau = std::sqrt(((mu + w.rho * k / omega) / w.rho) / (betaStar * k * omega));
        worstTime = std::max(worstTime, std::abs(time[q] - tau) / tau);
        double s = 0;
        for (double x : {0.1 + 20 * m.cells[q].z, 0.85 - 20 * m.cells[q].z})
          s = std::max(s, std::min(1.0, k / omega / schmidt * 400 / (betaStar * omega * x * (1 - x))));
        worstS = std::max(worstS, std::abs(segregation[q] - s) / s);
        sMin = std::min(sMin, s); sMax = std::max(sMax, s);
      }
    std::printf("  k = %g: s from %.4f to %.4f off the end columns\n", k, sMin, sMax);
    check(worstTime < 1e-12, k == 10 ? "k = 10: tau_mix against the formula (relative)" : "k = 25: tau_mix against the formula (relative)", worstTime, 1e-12);
    check(worstS < 1e-12, k == 10 ? "k = 10: s against the formula (relative)" : "k = 25: s against the formula (relative)", worstS, 1e-12);
  }

  std::printf("1(c) where the closure applies\n");
  {
    Duct d = duct(medium, fits, true, linearH2O2);
    Flow turbulent(d.definition);
    turbulent.setInitialState(d.cells, d.y);
    turbulent.setTurbulence(uniformTurbulence(d.cells.size(), 10, 1e3));
    d.definition.turbulence.enabled = false;
    Flow laminar(d.definition);
    laminar.setInitialState(d.cells, d.y);
    auto accepted = [&](Flow& flow, thermo::Chemistry chemistry) {
      thermo::ReactingFlow reacting(flow, mechanism, 1, 1e-8, 1e-14, chemistry);
      try { reacting.setMixingClosure(1, setNames); } catch (const std::invalid_argument&) { return false; }
      return true;
    };
    const int wrong = accepted(turbulent, thermo::Chemistry::LocalEquilibrium) + accepted(turbulent, thermo::Chemistry::Frozen) +
                      accepted(laminar, thermo::Chemistry::FiniteRate) + !accepted(turbulent, thermo::Chemistry::FiniteRate);
    check(wrong == 0, "4 cases (equilibrium, frozen, no turbulence, FiniteRate): wrong", wrong, 0);
  }

  std::printf("2. one substep of 3e-4 s against RK4 (stoichiometric H2/O2, 1200 K, 1 atm)\n");
  const double dt = 3e-4;
  long laminarLevels = 0, levels = 0;
  std::vector<double> laminarRecord, record;
  const auto laminarReference = reference(source, rho, unburnt, dt, set, 1, 0, laminarLevels, laminarRecord);
  for (const auto& [tauMix, s] : std::vector<std::pair<double, double>>{{1e-5, 1}, {1e-5, 0.9}, {1e-5, 0}}) {
    std::printf("  tau_mix %.0e s, s %.1f\n", tauMix, s);
    const auto ref = s == 0 ? laminarReference : reference(source, rho, unburnt, dt, set, tauMix, s, levels, record);
    if (s == 0) levels = laminarLevels;
    if (ref.empty()) { check(false, "reference converged (1 = no)", 1, 0); continue; }
    std::printf("    reference: %ld steps per sample; T at the end %.4f K\n", levels, ref.back()[0]);
    const auto fine = engine(source, rho, unburnt, dt, {tauMix, s}, set, 1e-10);
    const auto loose = engine(source, rho, unburnt, dt, {tauMix, s}, set, 1e-7);
    const auto d = difference(fine, ref), dl = difference(loose, ref);
    check(d.t < 1e-7, "rtol 1e-10: T against the reference (relative)", d.t, 1e-7);
    check(d.y < 1e-8, "rtol 1e-10: Y_k against the reference (absolute)", d.y, 1e-8);
    check(std::max(d.t, d.y) < std::max(dl.t, dl.y), "rtol 1e-7 error (the rtol 1e-10 error must be below it)",
          std::max(dl.t, dl.y), std::max(d.t, d.y));
    if (s > 0) {
      // Restated 5 October: over the fine record (1280 times), not the 20 samples.
      double effect = record.size() == laminarRecord.size() && record.size() == 1280 ? 0 : -1;
      std::size_t at = 0;
      for (std::size_t j = 0; effect >= 0 && j < record.size(); ++j) {
        const double gap = std::abs(record[j] - laminarRecord[j]) / laminarRecord[j];
        if (gap > effect) { effect = gap; at = j; }
      }
      std::printf("    fine record: largest T gap at t = %.4e s (T %.2f K with the closure, %.2f K laminar)\n",
                  (at + 1) * dt / 1280, record.empty() ? 0.0 : record[at], laminarRecord.empty() ? 0.0 : laminarRecord[at]);
      check(effect > 1e-2, "the closure's largest T change against laminar, fine record (must exceed)", effect, 1e-2);
    }
  }

  std::printf("3. limits\n");
  {
    bool same = true;
    for (const auto& z : states) {
      std::vector<double> a = z, b = z;
      thermo::ReactionStep one(source, 1e-10, 1e-16), two(source, 1e-10, 1e-16);
      one.advance(rho, a.data(), dt);
      two.advance(rho, b.data(), dt, {1e-6, 0}, set);
      same = same && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
    }
    check(same, "(a) s = 0 against the laminar substep, 3 states: bitwise (1 = differs)", same ? 0 : 1, 0);

    const auto laminar = engine(source, rho, unburnt, dt, {}, set, 1e-10);
    std::vector<double> gaps;
    for (double tauMix : {1e-8, 1e-9, 1e-10}) {
      gaps.push_back(difference(engine(source, rho, unburnt, dt, {tauMix, 1}, set, 1e-10), laminar).t);
      std::printf("    tau_mix %.0e s: largest relative T difference from laminar %.4e\n", tauMix, gaps.back());
    }
    const double order = std::log10(gaps[1] / gaps[2]);
    check(order >= 0.8, "(b) observed order in tau_mix, 1e-9 to 1e-10 (at least)", order, 0.8);

    std::vector<double> n2(ns + 1, 0.0), before;
    n2[0] = 1500; n2[1 + iN2] = 1;
    before = n2;
    const double rhoN2 = source.density(1500, p0, n2.data() + 1);
    thermo::ReactionStep inert(source, 1e-10, 1e-16);
    inert.advance(rhoN2, n2.data(), 1e-4, {1e-6, 1}, set);
    const bool unchanged = std::memcmp(n2.data(), before.data(), n2.size() * sizeof(double)) == 0;
    check(unchanged, "(c) pure N2 through 1e-4 s with (1e-6 s, 1): z changed (1 = yes)", unchanged ? 0 : 1, 0);

    auto inertMix = [&](double z, double* x) { x[iAR] = 0.1 + 20 * z; x[iN2] = 0.9 - 20 * z; };
    Duct d = duct(medium, fits, true, inertMix);
    Flow a(d.definition), b(d.definition), c(d.definition);
    for (Flow* f : {&a, &b, &c}) {
      f->setInitialState(d.cells, d.y);
      f->setTurbulence(uniformTurbulence(d.cells.size(), 10, 1e3));
    }
    thermo::ReactingFlow reacting(a, mechanism, 2, 1e-10, 1e-16, thermo::Chemistry::FiniteRate);
    reacting.setMixingClosure(1, setNames);
    const double taken = reacting.step();
    thermo::ReactingFlow without(c, mechanism, 2, 1e-10, 1e-16, thermo::Chemistry::FiniteRate);
    without.step();
    b.step(taken);
    auto flat = [](const Flow& f) {
      std::vector<double> v;
      for (const auto& u : f.state()) v.insert(v.end(), u.begin(), u.end());
      return v;
    };
    auto bitwise = [](const std::vector<double>& x, const std::vector<double>& y) {
      return x.size() == y.size() && std::memcmp(x.data(), y.data(), x.size() * sizeof(double)) == 0;
    };
    const bool identical = c.time() == a.time() && bitwise(flat(a), flat(c)) &&
                           bitwise(a.partialDensities(), c.partialDensities()) && bitwise(a.turbulence(), c.turbulence());
    check(identical, "(d)(i) with the closure against without it: bitwise (1 = differs)", identical ? 0 : 1, 0);
    double rhoMax = 0;
    for (const auto& u : b.state()) rhoMax = std::max(rhoMax, u[0]);
    const double momentumScale = std::sqrt(rhoMax * p0);
    double worst = 0;
    // fixed[f] > 0 is the scale of field f; otherwise the field's own max.
    auto compare = [&](const std::vector<double>& x, const std::vector<double>& y, std::size_t stride,
                       const std::vector<double>& fixed) {
      for (std::size_t f = 0; f < stride; ++f) {
        double scale = 0, gap = 0;
        for (std::size_t q = 0; q < x.size() / stride; ++q) {
          scale = std::max(scale, std::abs(y[q * stride + f]));
          gap = std::max(gap, std::abs(x[q * stride + f] - y[q * stride + f]));
        }
        if (f < fixed.size() && fixed[f] > 0) scale = fixed[f];
        worst = std::max(worst, scale > 0 ? gap / scale : (gap > 0 ? 1.0 : 0.0));
      }
    };
    compare(flat(a), flat(b), a.state()[0].size(), {0, momentumScale, momentumScale, 0});
    compare(a.partialDensities(), b.partialDensities(), ns, {});
    compare(a.turbulence(), b.turbulence(), 2, {});
    check(worst < 1e-14, "(d)(ii) against Flow (relative; momenta to sqrt(rho_max p0))", worst, 1e-14);
  }

  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
