// Item-1 verification: the Cantera-backed ideal rocket against NASA CEA (via RocketCEA,
// tools/cea_reference.py -> tools/cea_reference.csv). Tolerance stated before the run:
// chamber temperature, c* and vacuum Isp each within 0.5% of CEA, equilibrium and frozen.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "adapters/thermo.hpp"

using namespace crucible::thermo;

int main(int argc, char** argv) {
  const std::string csv = argc > 1 ? argv[1] : CRUCIBLE_SOURCE_DIR "/tools/cea_reference.csv";
  std::ifstream in(csv);
  if (!in) { std::fprintf(stderr, "cannot open %s\n", csv.c_str()); return 2; }

  struct Pair { std::string mechanism; Reactant fuel, oxidizer; };
  const std::map<std::string, Pair> pairs = {
      {"H2/O2", {"h2o2.yaml", {"H2:1"}, {"O2:1"}}},
      {"CH4/O2", {"gri30.yaml", {"CH4:1"}, {"O2:1"}}},
      {"LH2/LOX", {"h2o2.yaml", {"H2:1", 20.27, -9012.0}, {"O2:1", 90.17, -12979.0}}},
  };
  std::map<std::string, std::unique_ptr<Mixture>> mixtures;

  const double tol = 0.005;
  int failures = 0, rows = 0;
  double worst = 0;
  std::string line;
  std::getline(in, line);
  std::printf("%-8s %4s %-11s %9s %9s %7s | %8s %8s %7s | %7s %7s %7s\n", "pair", "O/F", "mode",
              "Tc", "Tc_CEA", "dTc%", "c*", "c*_CEA", "dc*%", "Isp", "Isp_CEA", "dIsp%");
  while (std::getline(in, line)) {
    std::stringstream ss(line);
    std::string pair, mode, f;
    double v[7];
    std::getline(ss, pair, ',');
    for (int i = 0; i < 3; ++i) { std::getline(ss, f, ','); v[i] = std::stod(f); }
    std::getline(ss, mode, ',');
    for (int i = 3; i < 6; ++i) { std::getline(ss, f, ','); v[i] = std::stod(f); }
    const Pair& p = pairs.at(pair);
    auto& mix = mixtures[pair];
    if (!mix) mix = std::make_unique<Mixture>(p.mechanism);
    const auto r = mix->idealRocket(p.fuel, p.oxidizer, v[0], v[1] * 1e5, v[2],
                                    mode == "frozen" ? Expansion::FrozenAtChamber : Expansion::Equilibrium);
    const double dT = r.tc / v[3] - 1, dC = r.cstar / v[4] - 1, dI = r.ispVac / v[5] - 1;
    const double m = std::max({std::abs(dT), std::abs(dC), std::abs(dI)});
    worst = std::max(worst, m);
    const bool ok = m <= tol;
    failures += !ok;
    ++rows;
    std::printf("%-8s %4.1f %-11s %9.2f %9.2f %+7.3f | %8.2f %8.2f %+7.3f | %7.2f %7.2f %+7.3f %s\n",
                pair.c_str(), v[0], mode.c_str(), r.tc, v[3], 100 * dT, r.cstar, v[4], 100 * dC,
                r.ispVac, v[5], 100 * dI, ok ? "" : "FAIL");
  }
  std::printf("%d rows, worst |delta| %.3f%% (tolerance %.1f%%), %d failures\n", rows, 100 * worst,
              100 * tol, failures);
  return failures == 0 && rows > 0 ? 0 : 1;
}
