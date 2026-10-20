#include "cantera/core.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
// Rewrites a CRUCIBLE_STATE_SAMPLES file's eta from mechanism a's energies to mechanism b's, so that a run on
// one mechanism can be placed in a Table B built on another (docs/evidence/TABLE_B.md, "The fallback rate").
// Each cell's complete-combustion products follow from (f_H, Z_N) alone (zH = f_H (1 - Z_N),
// zO = (1 - f_H)(1 - Z_N)), as in tests/chamber_study.cpp and the builder's productsOf.
//
// Usage: crucible_state_samples_eta <mechanism a> <mechanism b> <samples in> <samples out>
struct Ends { double lo, hi; };
static Ends ends(Cantera::ThermoPhase& g, double fH, double zN) {
  std::vector<double> y(g.nSpecies(), 0.0);
  const std::size_t kH2O = g.speciesIndex("H2O"), kH2 = g.speciesIndex("H2"), kO2 = g.speciesIndex("O2"),
                    kN2 = g.speciesIndex("N2");
  const double hInWater = 2 * g.atomicWeight(g.elementIndex("H")) / g.molecularWeight(kH2O);
  const double zH = fH * (1 - zN), zO = (1 - fH) * (1 - zN);
  if (zH * (1 - hInWater) >= zO * hInWater) {
    y[kH2O] = zO / (1 - hInWater); y[kH2] = std::max(zH - y[kH2O] * hInWater, 0.0);
  } else {
    y[kH2O] = zH / hInWater; y[kO2] = std::max(zO - y[kH2O] * (1 - hInWater), 0.0);
  }
  y[kN2] = zN;
  g.setMassFractions_NoNorm(y.data());
  g.setState_TD(200, 1.0); const double lo = g.intEnergy_mass();
  g.setState_TD(6000, 1.0); const double hi = g.intEnergy_mass();
  return {lo, hi};
}
int main(int, char** argv) {
  auto a = Cantera::newSolution(argv[1], "", "none")->thermo();
  auto b = Cantera::newSolution(argv[2], "", "none")->thermo();
  std::ifstream in(argv[3]);
  std::FILE* out = std::fopen(argv[4], "w");
  std::string line;
  std::getline(in, line);
  std::fprintf(out, "%s\n", line.c_str());
  double worst = 0;
  while (std::getline(in, line)) {
    std::vector<std::string> c;
    std::stringstream s(line); std::string f;
    while (std::getline(s, f, ',')) c.push_back(f);
    const double fH = std::atof(c[4].c_str()), zN = std::atof(c[5].c_str()), eta = std::atof(c[6].c_str());
    if (std::isfinite(fH) && std::isfinite(eta)) {
      const Ends ea = ends(*a, fH, zN), eb = ends(*b, fH, zN);
      const double e = ea.lo + eta * (ea.hi - ea.lo), etaB = (e - eb.lo) / (eb.hi - eb.lo);
      worst = std::max(worst, std::abs(etaB - eta));
      char buf[32]; std::snprintf(buf, sizeof buf, "%.9e", etaB); c[6] = buf;
    }
    for (std::size_t k = 0; k < c.size(); ++k) std::fprintf(out, k ? ",%s" : "%s", c[k].c_str());
    std::fprintf(out, "\n");
  }
  std::fclose(out);
  std::printf("eta rewritten from %s to %s energies; worst change %.4e\n", argv[1], argv[2], worst);
}
