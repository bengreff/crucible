#pragma once
// Real-gas thermochemistry adapter (TECHNICAL_PLAN: Cantera behind an adapter, one
// context per thread). The core never sees Cantera types; each Mixture owns its own
// Cantera Solution and must not be shared across threads.
#include <memory>
#include <optional>
#include <string>

namespace crucible::thermo {

// A propellant as fed: mole composition ("H2:1"), feed temperature, and optionally an
// explicit molar enthalpy [J/mol] for condensed feeds whose phase the mechanism lacks
// (liquid H2/O2: values from CEA thermo.lib so both codes start from the same h0).
struct Reactant {
  std::string composition;
  double temperature = 298.15;
  std::optional<double> molarEnthalpy;
};

enum class Expansion { Equilibrium, FrozenAtChamber };

// One-dimensional ideal rocket (infinite-area chamber, isentropic expansion), the
// same problem NASA CEA solves for its "rocket" case.
struct RocketPoint {
  double tc = 0;            // chamber temperature [K]
  double pc = 0;            // chamber pressure [Pa]
  double cstar = 0;         // characteristic velocity [m/s]
  double ispVac = 0;        // vacuum specific impulse [s]
  double cf = 0;            // vacuum thrust coefficient
  double tExit = 0;         // exit temperature [K]
  double pExit = 0;         // exit pressure [Pa]
  double gammaChamber = 0;  // frozen cp/cv in the chamber
  double mwChamber = 0;     // chamber mean molar mass [kg/kmol]
};

class Mixture {
 public:
  explicit Mixture(const std::string& mechanism);
  ~Mixture();
  Mixture(const Mixture&) = delete;
  Mixture& operator=(const Mixture&) = delete;

  // Specific enthalpy of a reactant as fed [J/kg].
  double reactantEnthalpy(const Reactant& r, double pressure);

  RocketPoint idealRocket(const Reactant& fuel, const Reactant& oxidizer, double mixtureRatio,
                          double pc, double areaRatio, Expansion mode);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// Directory holding the Cantera mechanism files (set at build time).
std::string dataDirectory();

}  // namespace crucible::thermo
