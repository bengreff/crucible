#include "adapters/thermo.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "cantera/core.h"

namespace crucible::thermo {

namespace {
constexpr double g0 = 9.80665;
}

std::string dataDirectory() { return CRUCIBLE_CANTERA_DATA; }

struct Mixture::Impl {
  std::shared_ptr<Cantera::Solution> solution;
  std::shared_ptr<Cantera::ThermoPhase> gas;
};

Mixture::Mixture(const std::string& mechanism) : impl_(std::make_unique<Impl>()) {
  Cantera::addDataDirectory(dataDirectory());
  impl_->solution = Cantera::newSolution(mechanism, "", "none");
  impl_->gas = impl_->solution->thermo();
}

Mixture::~Mixture() = default;

double Mixture::reactantEnthalpy(const Reactant& r, double pressure) {
  auto& gas = *impl_->gas;
  // Condensed feeds carry their own enthalpy; the gas phase only supplies composition.
  gas.setState_TPX(r.molarEnthalpy ? 298.15 : r.temperature, pressure, r.composition);
  if (r.molarEnthalpy) return *r.molarEnthalpy * 1000.0 / gas.meanMolecularWeight();
  return gas.enthalpy_mass();
}

RocketPoint Mixture::idealRocket(const Reactant& fuel, const Reactant& oxidizer, double mixtureRatio,
                                 double pc, double areaRatio, Expansion mode) {
  auto& gas = *impl_->gas;
  const std::size_t n = gas.nSpecies();
  const double wf = 1.0 / (1.0 + mixtureRatio), wo = mixtureRatio / (1.0 + mixtureRatio);

  // Reactant mass fractions and the mixture's specific enthalpy as fed.
  std::vector<double> y(n), yf(n), yo(n);
  const double hf = reactantEnthalpy(fuel, pc);
  gas.getMassFractions(yf.data());
  const double ho = reactantEnthalpy(oxidizer, pc);
  gas.getMassFractions(yo.data());
  for (std::size_t k = 0; k < n; ++k) y[k] = wf * yf[k] + wo * yo[k];
  const double h0 = wf * hf + wo * ho;

  // Chamber: adiabatic equilibrium at Pc (stagnation, infinite-area chamber).
  // Cryogenic feeds put h0 below any gas-phase state of the unreacted mixture, so start
  // from a hot equilibrium (same elements) and move it onto h0 before the HP solve.
  gas.setState_TPY(3000.0, pc, y.data());
  gas.equilibrate("TP");
  gas.setState_HP(h0, pc);
  gas.equilibrate("HP");
  RocketPoint out;
  out.tc = gas.temperature();
  out.pc = pc;
  out.gammaChamber = gas.cp_mass() / gas.cv_mass();
  out.mwChamber = gas.meanMolecularWeight();
  const double s0 = gas.entropy_mass();
  std::vector<double> yc(n);
  gas.getMassFractions(yc.data());

  // State on the isentrope at pressure p; returns mass flux rho*u.
  struct Station { double g, u, rho, t, p; };
  auto expand = [&](double p) {
    if (mode == Expansion::FrozenAtChamber) {
      gas.setMassFractions(yc.data());
      gas.setState_SP(s0, p);
    } else {
      gas.setState_SP(s0, p);
      gas.equilibrate("SP");
    }
    const double dh = h0 - gas.enthalpy_mass();
    const double u = std::sqrt(std::max(2.0 * dh, 0.0));
    return Station{gas.density() * u, u, gas.density(), gas.temperature(), p};
  };

  // Throat: maximum mass flux along the isentrope (golden section in ln p).
  double a = std::log(0.3 * pc), b = std::log(0.9 * pc);
  const double r = 0.5 * (std::sqrt(5.0) - 1.0);
  double c = b - r * (b - a), d = a + r * (b - a);
  double gc = expand(std::exp(c)).g, gd = expand(std::exp(d)).g;
  while (b - a > 1e-7) {
    if (gc > gd) { b = d; d = c; gd = gc; c = b - r * (b - a); gc = expand(std::exp(c)).g; }
    else { a = c; c = d; gc = gd; d = a + r * (b - a); gd = expand(std::exp(d)).g; }
  }
  const Station throat = expand(std::exp(0.5 * (a + b)));
  out.cstar = pc / throat.g;

  // Exit: supersonic branch where rho*u*eps equals the throat mass flux.
  double lo = std::log(throat.p * 1e-7), hi = std::log(throat.p);
  // Below the throat pressure the flux falls monotonically: flux*eps above the throat
  // value means the local area is still smaller than the exit, so expand further.
  for (int i = 0; i < 200 && hi - lo > 1e-10; ++i) {
    const double m = 0.5 * (lo + hi);
    if (expand(std::exp(m)).g * areaRatio > throat.g) hi = m; else lo = m;
  }
  const Station exit = expand(std::exp(0.5 * (lo + hi)));
  const double vacuumVelocity = exit.u + exit.p / (exit.rho * exit.u);
  out.ispVac = vacuumVelocity / g0;
  out.cf = vacuumVelocity / out.cstar;
  out.tExit = exit.t;
  out.pExit = exit.p;
  return out;
}

}  // namespace crucible::thermo
