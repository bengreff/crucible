// C1 eqtt as tests/chamber_study.cpp builds it, and its chamber pressure, for the AMR drivers
// (tests/amr_budget.cpp, tests/amr_c4.cpp). Moved here from amr_budget.cpp, the same arithmetic.
#pragma once

#include <cmath>
#include <vector>

#include "adapters/reaction.hpp"
#include "cantera/core.h"
#include "core/flow.hpp"
#include "tests/c1_geometry.hpp"

namespace crucible::c1_eqtt {
using namespace crucible::c1_geometry;

constexpr double kMixtureRatio = 5, kMassFlow = 0.4, kAmbientP = 1000, kAmbientT = 300, kIgniter = 0.5;

// C1 eqtt as tests/chamber_study.cpp builds it, on nz x nr.
inline Definition c1eqtt(thermo::ReactionSource& source, int nz, int nr) {
  Definition d;
  d.experiment = Case::Chamber;
  d.nz = nz;
  d.nr = nr;
  d.contour = geometry().contour;
  d.species = source.medium().species();
  d.composition = source.massFractions("N2:1");
  d.backPressure = kAmbientP;
  d.ambientTemperature = kAmbientT;
  Supply s;
  s.outerRadius = kRc;
  s.massFlow = kMassFlow;
  s.totalTemperature = 300;
  auto h2 = source.massFractions("H2:1"), o2 = source.massFractions("O2:1");
  s.composition.resize(h2.size());
  for (std::size_t k = 0; k < h2.size(); ++k) s.composition[k] = (h2[k] + kMixtureRatio * o2[k]) / (1 + kMixtureRatio);
  s.opens = 0;
  s.ramp = 5e-4;
  s.turbulenceIntensity = 0.05;
  s.viscosityRatio = 10;
  d.supplies = {s};
  d.igniter = {0.005, 0.015, 0.010, kIgniter, 2e-4, 2e-4};
  d.transport = thermo::transportFits("h2o2.yaml");
  d.wallTemperature = 600;
  d.turbulence.enabled = true;
  const Medium medium = d.medium();
  std::vector<double> diffusion(d.species.size()), work;
  const double a0 = medium.soundSpeed(kAmbientT, d.composition.data());
  const double rho0 = kAmbientP / (medium.gasConstant(d.composition.data()) * kAmbientT);
  const double mu0 = medium.transport(kAmbientT, kAmbientP, d.composition.data(), diffusion.data(), work).viscosity;
  d.turbulence.ambientK = 9e-9 * a0 * a0;
  d.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
  return d;
}

// Equilibrium isentropic stagnation pressure of a cell, as tests/chamber_study.cpp: Newton on ln p0 along the
// isentrope. It stops at the accuracy equilibrate's relative tolerance (1e-9, its default) allows the step,
// 10 x 1e-9 x (1 + |s|/R) (AMR.md, *Amendment: the stagnation loop exits at the equilibrium solver's accuracy*);
// `oldExit` keeps the exit before it, 1e-13, for that amendment's check. `iterations` gets the SP solves made.
inline double stagnationPressure(Cantera::ThermoPhase& gas, double t, double p, const double* y, double kinetic,
                                 bool oldExit = false, int* iterations = nullptr) {
  constexpr double rtol = 1e-9;
  gas.setState_TPY(t, p, y);
  gas.equilibrate("TP");
  const double s = gas.entropy_mass(), h0 = gas.enthalpy_mass() + kinetic;
  double p0 = p;
  int solves = 0;
  for (int n = 0; n < 50; ++n) {
    gas.setState_SP(s, p0);
    gas.equilibrate("SP");
    ++solves;
    double step = (h0 - gas.enthalpy_mass()) / (p0 / gas.density());
    p0 *= std::exp(step);
    const double exit =
        oldExit ? 1e-13 : 10 * rtol * (1 + std::abs(s) * gas.meanMolecularWeight() / Cantera::GasConstant);
    if (std::abs(step) < exit) break;
  }
  if (iterations) *iterations = solves;
  return p0;
}
// One cell's mass flux through its downstream axial face, and that times its stagnation pressure.
inline void addStagnation(const Flow& flow, Cantera::ThermoPhase& gas, int i, int j, double& p0, double& flux,
                          bool oldExit = false, int* iterations = nullptr) {
  const auto& mesh = flow.mesh();
  auto q = mesh.index(i, j);
  auto w = flow.cellPrimitive(q);
  auto y = flow.massFractions(q);
  double g = w.rho * w.uz * mesh.axialArea(i + 1, j);
  p0 += g * stagnationPressure(gas, flow.temperature(q), w.p, y.data(), 0.5 * (w.uz * w.uz + w.ur * w.ur), oldExit,
                               iterations);
  flux += g;
}
// Mass-flux-averaged equilibrium stagnation pressure of a column (chamber_study's chamber pressure).
inline double columnStagnation(const Flow& flow, Cantera::ThermoPhase& gas, int column) {
  double p0 = 0, flux = 0;
  for (int j = 0; j < flow.mesh().nr; ++j) addStagnation(flow, gas, column, j, p0, flux);
  return p0 / flux;
}
}  // namespace crucible::c1_eqtt
