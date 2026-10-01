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
  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
