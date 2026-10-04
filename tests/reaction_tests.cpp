// Item-4 verification: CRUCIBLE-owned stiff reaction integration (CVODES BDF over Cantera
// local rates) for a closed adiabatic constant-volume parcel. Tolerances stated before the run:
//   ignition delay (T = T0 + 400 K) within 0.1% of Cantera's IdealGasReactor/ReactorNet
//     integrated at rtol 1e-12, both sampled on the same time grid;
//   T(t_end) within 0.01% and every Y_k within 1e-6 of that reference;
//   element mass fractions conserved to 1e-8 relative; internal energy to 1e-6 of cv0*T0;
//   final state within 0.5 K and 1e-4 in every Y_k of Cantera's UV equilibrium.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "adapters/reaction.hpp"
#include "cantera/core.h"
#include "cantera/zeroD/ReactorFactory.h"
#include "cantera/zeroD/ReactorNet.h"

using namespace crucible::thermo;

namespace {

struct Case {
  const char* name;
  const char* mechanism;
  const char* moles;
  double t0, p0, tEnd;
};

// First crossing of tCross by linear interpolation between samples.
double crossing(const std::vector<double>& times, const std::vector<double>& temps, double tCross) {
  for (std::size_t i = 1; i < temps.size(); ++i)
    if (temps[i - 1] < tCross && temps[i] >= tCross)
      return times[i - 1] + (tCross - temps[i - 1]) / (temps[i] - temps[i - 1]) * (times[i] - times[i - 1]);
  return NAN;
}

double seconds(std::chrono::steady_clock::time_point a) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
}

}  // namespace

int main() {
  const std::vector<Case> cases = {
      {"H2/O2 1000K 1atm", "h2o2.yaml", "H2:2, O2:1", 1000.0, 101325.0, 0.02},
      {"CH4/O2 1400K 20bar", "gri30.yaml", "CH4:1, O2:2", 1400.0, 20e5, 0.02},
  };
  const int samples = 200000;
  int failures = 0;
  auto check = [&](bool ok, const char* what, double value, double tol) {
    std::printf("  %-34s %12.4e  (tol %.1e) %s\n", what, value, tol, ok ? "ok" : "FAIL");
    failures += !ok;
  };

  for (const Case& c : cases) {
    std::printf("%s (%s)\n", c.name, c.mechanism);
    ReactionSource source(c.mechanism);
    const std::size_t n = source.nSpecies();
    std::vector<double> z(n + 1);
    z[0] = c.t0;
    const auto y0 = source.massFractions(c.moles);
    std::copy(y0.begin(), y0.end(), z.begin() + 1);
    const double rho = source.density(c.t0, c.p0, y0.data());
    const double u0 = source.internalEnergy(c.t0, rho, y0.data());
    const double cv0 = source.cv(c.t0, rho, y0.data());
    const auto e0 = source.elementMassFractions(y0.data());

    // CRUCIBLE integration.
    std::vector<double> times(samples + 1), temps(samples + 1), refTemps(samples + 1);
    times[0] = 0.0;
    temps[0] = refTemps[0] = c.t0;
    int sample = 0;
    ReactionStep step(source, 1e-8, 1e-14);
    auto clock = std::chrono::steady_clock::now();
    step.advance(rho, z.data(), c.tEnd, samples, [&](double t, const double* s) {
      ++sample;
      times[sample] = t;
      temps[sample] = s[0];
    });
    const double ours = seconds(clock);

    // Reference: Cantera's own reactor network on an independent Solution.
    auto solution = Cantera::newSolution(c.mechanism, "", "none");
    solution->thermo()->setState_TPX(c.t0, c.p0, c.moles);
    auto reactor = Cantera::newReactor4("IdealGasReactor", solution);
    Cantera::ReactorNet net(reactor);
    net.setTolerances(1e-12, 1e-20);
    net.setMaxSteps(1000000);
    clock = std::chrono::steady_clock::now();
    for (int i = 1; i <= samples; ++i) {
      net.advance(times[i]);
      refTemps[i] = reactor->temperature();
    }
    const double theirs = seconds(clock);
    const double* yRef = reactor->massFractions();

    const double tau = crossing(times, temps, c.t0 + 400.0);
    const double tauRef = crossing(times, refTemps, c.t0 + 400.0);
    double dyRef = 0;
    for (std::size_t k = 0; k < n; ++k) dyRef = std::max(dyRef, std::abs(z[k + 1] - yRef[k]));
    const auto e1 = source.elementMassFractions(z.data() + 1);
    double de = 0;
    for (std::size_t m = 0; m < e0.size(); ++m)
      if (e0[m] > 0) de = std::max(de, std::abs(e1[m] - e0[m]) / e0[m]);
    const double du = std::abs(source.internalEnergy(z[0], rho, z.data() + 1) - u0) / (cv0 * c.t0);
    std::vector<double> eq = z;
    source.equilibrateUV(rho, eq.data());
    double dyEq = 0;
    for (std::size_t k = 0; k < n; ++k) dyEq = std::max(dyEq, std::abs(z[k + 1] - eq[k + 1]));

    std::printf("  tau %.6e s (reference %.6e s); T_end %.3f K (reference %.3f, UV equilibrium %.3f)\n",
                tau, tauRef, z[0], reactor->temperature(), eq[0]);
    std::printf("  CVODES steps %ld, wall %.2f s (reference wall %.2f s, %d samples)\n", step.steps(),
                ours, theirs, samples);
    check(std::abs(tau / tauRef - 1) <= 1e-3, "ignition delay rel. difference", std::abs(tau / tauRef - 1), 1e-3);
    check(std::abs(z[0] / reactor->temperature() - 1) <= 1e-4, "T_end rel. difference",
          std::abs(z[0] / reactor->temperature() - 1), 1e-4);
    check(dyRef <= 1e-6, "max |Y - Y_ref| at t_end", dyRef, 1e-6);
    check(de <= 1e-8, "element mass fraction drift", de, 1e-8);
    check(du <= 1e-6, "internal energy drift / (cv0 T0)", du, 1e-6);
    check(std::abs(z[0] - eq[0]) <= 0.5, "|T_end - T_eq| [K]", std::abs(z[0] - eq[0]), 0.5);
    check(dyEq <= 1e-4, "max |Y - Y_eq|", dyEq, 1e-4);
  }

  // Mixture-averaged transport: the core's evaluation from exported fits against Cantera's own
  // MixTransport over temperature, pressure and composition, including trace and absent species.
  for (const char* mechanism : {"h2o2.yaml", "gri30.yaml"}) {
    std::printf("mixture-averaged transport (%s)\n", mechanism);
    ReactionSource source(mechanism);
    auto medium = source.medium();
    medium.setTransport(transportFits(mechanism));
    auto solution = Cantera::newSolution(mechanism, "", "mixture-averaged");
    auto gas = solution->thermo();
    auto transport = solution->transport();
    const std::size_t n = source.nSpecies();
    const char* mixtures[] = {"H2:2, O2:1", "H2:1, O2:1, H2O:2, OH:0.1, H:0.05, O:0.02", "N2:1", "O2:1, N2:3.76, H2:1e-12",
                              "H2O:1, H2:0.3, O2:0.1, HO2:1e-6"};
    double worstMu = 0, worstLambda = 0, worstD = 0;
    std::vector<double> d(n), dRef(n), work;
    for (const char* x : mixtures)
      for (double t : {250.0, 300.0, 800.0, 1500.0, 2500.0, 3400.0})
        for (double p : {1e4, 101325.0, 3e6, 2e7}) {
          gas->setState_TPX(t, p, x);
          const double* y = gas->massFractions();
          const auto mine = medium.transport(t, p, y, d.data(), work);
          transport->getMixDiffCoeffs(dRef.data());
          worstMu = std::max(worstMu, std::abs(mine.viscosity / transport->viscosity() - 1));
          worstLambda = std::max(worstLambda, std::abs(mine.conductivity / transport->thermalConductivity() - 1));
          // D_km carries the factor 1 - Y_k, which Cantera evaluates as (W - X_k W_k)/W: its round-off
          // relative to 1 - Y_k exceeds 1e-13 once 1 - Y_k < 1e-3. There the difference is scaled by
          // the largest coefficient of the state.
          const double scale = *std::max_element(dRef.begin(), dRef.end());
          for (std::size_t k = 0; k < n; ++k) {
            double e = 1 - y[k] >= 1e-3 ? std::abs(d[k] / dRef[k] - 1) : std::abs(d[k] - dRef[k]) / scale;
            worstD = std::max(worstD, e);
          }
        }
    check(worstMu <= 1e-12, "viscosity max rel. difference", worstMu, 1e-12);
    check(worstLambda <= 1e-12, "conductivity max rel. difference", worstLambda, 1e-12);
    check(worstD <= 1e-12, "D_km max rel. difference (see note)", worstD, 1e-12);
  }
  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
