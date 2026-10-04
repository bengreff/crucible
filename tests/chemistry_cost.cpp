// Cost of one cell's chemistry call for each mode, on states a rocket chamber visits
// (h2o2.yaml, H2/O2 at O/F 6 by mass, 3 MPa). Measured on the machine that runs it; used to size
// chamber meshes (TECHNICAL_PLAN step 15). Not a verification test.
//
// Usage: crucible_chemistry_cost [calls per case, 2000]
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "adapters/reaction.hpp"
#include "cantera/core.h"

using namespace crucible::thermo;

namespace {
double seconds(auto&& f, int calls) {
  auto t0 = std::chrono::steady_clock::now();
  for (int i = 0; i < calls; ++i) f(i);
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() / calls;
}
}  // namespace

int main(int argc, char** argv) {
  const int calls = argc > 1 ? std::atoi(argv[1]) : 2000;
  ReactionSource source("h2o2.yaml");
  const std::size_t n = source.nSpecies();
  // Reactants at 300 K, O/F 6 by mass (H2 mass fraction 1/7).
  std::vector<double> fresh(n + 1, 0.0);
  auto y = source.massFractions("H2:1");
  auto yo = source.massFractions("O2:1");
  for (std::size_t k = 0; k < n; ++k) fresh[k + 1] = y[k] / 7 + yo[k] * 6 / 7;
  fresh[0] = 300;
  const double rho = source.density(300, 3e6, fresh.data() + 1);
  // The burned chamber state: constant-(h, p) equilibrium at 3 MPa, taken from Cantera directly.
  auto sol = Cantera::newSolution("h2o2.yaml", "", "none");
  auto& gas = *sol->thermo();
  gas.setMassFractions(fresh.data() + 1);
  gas.setState_TP(300, 3e6);
  gas.equilibrate("HP");
  std::vector<double> burned(n + 1);
  burned[0] = gas.temperature();
  gas.getMassFractions(burned.data() + 1);
  const double rhoBurned = gas.density();
  std::printf("burned state: T = %.1f K, p = %.3f MPa\n", burned[0], gas.pressure() / 1e6);

  std::vector<double> z(n + 1);
  // Near-equilibrium cell: the burned state with a small temperature disturbance, as after one
  // flow step (|dT| ~ 1 K).
  auto near = [&](int i) { z = burned; z[0] += 1.0 + 1e-3 * (i % 7); };
  for (const char* solver : {"auto", "element_potential", "gibbs", "vcs"}) {
    for (const char* pair : {"UV", "TV"}) {
      double t = 0;
      try {
        t = seconds([&](int i) {
          near(i);
          gas.setMassFractions(z.data() + 1);
          gas.setState_TD(z[0], rhoBurned);
          gas.equilibrate(pair, solver);
        }, calls / 20 + 1);
        std::printf("Cantera %s equilibrium, solver %s, near-equilibrium cell: %.1f us\n", pair, solver, t * 1e6);
      } catch (const std::exception& e) {
        std::printf("Cantera %s equilibrium, solver %s: failed\n", pair, solver);
      }
    }
  }
  for (double rtol : {1e-8, 1e-6}) {
    ReactionStep s(source, rtol, rtol * 1e-6);
    long steps = 0;
    double t = seconds([&](int i) { near(i); s.advance(rhoBurned, z.data(), 5e-8); steps += s.steps(); }, calls);
    std::printf("finite rate 50 ns, rtol %.0e, near-equilibrium cell: %.1f us per call, %.1f CVODES steps\n",
                rtol, t * 1e6, double(steps) / calls);
    steps = 0;
    t = seconds([&](int) { z = fresh; z[0] = 1500; s.advance(rho, z.data(), 5e-8); steps += s.steps(); }, calls / 10 + 1);
    std::printf("finite rate 50 ns, rtol %.0e, fresh at 1500 K:       %.1f us per call, %.1f CVODES steps\n",
                rtol, t * 1e6, double(steps) / (calls / 10 + 1));
  }
  double rateCall = seconds([&](int i) { near(i); std::vector<double> d(n + 1); source.rates(rhoBurned, z.data(), d.data()); }, calls);
  std::printf("one rates() evaluation: %.2f us\n", rateCall * 1e6);
  return 0;
}
