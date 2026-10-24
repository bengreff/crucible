#pragma once

#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "core/medium.hpp"

namespace crucible::thermo {

// Cantera-backed local chemistry for a closed, adiabatic, constant-volume parcel.
// State vector z = [T, Y_1..Y_K]; density is held fixed. Construct one per executing
// thread (TECHNICAL_PLAN: independent Cantera contexts per thread).
class ReactionSource {
 public:
  explicit ReactionSource(const std::string& mechanism);
  ~ReactionSource();
  ReactionSource(const ReactionSource&) = delete;
  ReactionSource& operator=(const ReactionSource&) = delete;

  std::size_t nSpecies() const;
  std::size_t nElements() const;
  std::string elementName(std::size_t m) const;  // in the order of elementMassFractions
  std::string speciesName(std::size_t k) const;
  // The mechanism's species as core thermally perfect species (NASA 7-coefficient data).
  Medium medium() const;

  // Mass fractions of a mole-fraction composition such as "CH4:1, O2:2".
  std::vector<double> massFractions(const std::string& moles);
  double density(double t, double p, const double* y);
  double pressure(double t, double rho, const double* y);
  double internalEnergy(double t, double rho, const double* y);  // [J/kg]
  double cv(double t, double rho, const double* y);              // [J/(kg K)]

  // Element mass fractions of y (no normalisation applied).
  std::vector<double> elementMassFractions(const double* y);

  // Each reaction's species with a nonzero reactant, product and net coefficient (a third-body collider is
  // among the first two only if the equation names it, as in H + O2 + N2 <=> HO2 + N2), and whether it runs in
  // reverse too.
  struct ReactionSpecies {
    std::vector<std::size_t> reactants, products, changed;
    bool reversible = false;
  };
  std::vector<ReactionSpecies> reactionSpecies() const;

  // dz/dt at fixed density: dY_k/dt = w_k W_k / rho, dT/dt = -sum(u_k w_k) / (rho cv).
  void rates(double rho, const double* z, double* dzdt);
  // The exact Jacobian of rates() at fixed density, column-major: jac[j * n + i] = d(dz_i/dt)/dz_j with
  // n = nSpecies() + 1 (docs/evidence/TABLE_B.md, "The analytic Jacobian"). Built from Cantera's analytic
  // net-production derivatives (third bodies and falloff included) and the NASA-7 heat capacities.
  // dzdt, if given, returns rates() at the same state.
  void jacobian(double rho, const double* z, double* jac, double* dzdt = nullptr);

  // Adiabatic constant-(u, v) equilibrium of the same elements; z holds the start
  // state and returns the equilibrium.
  void equilibrateUV(double rho, double* z);
  // Isothermal constant-(T, v) equilibrium at T = z[0] (built tables).
  void equilibrateTV(double rho, double* z);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  void equilibrate(double rho, double* z, const char* constraints);
};

// Cantera's mixture-averaged transport fits for the mechanism's species, for the core's own
// evaluation (Medium::transport, no Cantera calls in the flow step). Throws if the mechanism has
// no transport data or uses the CHEMKIN fit form. Cantera fits from the highest of the species' lower
// thermo limits (300 K for h2o2.yaml); below that the polynomials are extrapolations and fail (H2's
// conductivity turns negative below about 97 K). minTemperature > 0 fits from there instead, for flows
// with gas that cold (liquid oxygen; docs/evidence/SPRAY.md). The kinetic theory's collision integrals
// are tabulated from T* = kT/epsilon = 0.1 (60 K for H2O).
TransportFits transportFits(const std::string& mechanism, double minTemperature = 0);

// The species with lowest ranges, 20 to 200 K (extrapolated below 20 K), for H2 (normal hydrogen), O2 and H2O
// (docs/evidence/LOW_T_THERMO.md; tools/mechanisms/make_low_t.py), each joined at 200 K to TM-4513's low range
// (h2o2_highT.yaml) with cp, h and s continuous. At and above 200 K nothing changes. Other species are returned
// as they are (their low ranges continue below 200 K, as before). Throws if H2, O2 or H2O has another low range.
// For the flow's medium only: the kinetics keep the mechanism's own thermo below 200 K.
std::vector<Species> withLowTemperatureRanges(std::vector<Species> species);

// The PaSR closure's inputs for one cell (TECHNICAL_PLAN step 7, "PaSR"), from the flow
// (Flow::mixingInputs) and frozen over a reaction substep: the mixing time tau_mix [s] and the
// segregation s in [0, 1]. s = 0 is the mean-state (laminar) rate.
struct Mixing {
  double time = std::numeric_limits<double>::infinity();
  double segregation = 0;
};

// The PaSR factor kappa_eff = 1 - s (1 - kappa) on the rate dz/dt of z = [T, Y_1..Y_K], with
// kappa = tau_c / (tau_c + tau_mix) and tau_c the largest Y_i / |dY_i/dt| over the given species
// (indices into Y) with Y_i > 0 and dY_i/dt != 0; kappa = 1 if there is none. Exactly 1 when s = 0.
double reactingFraction(const double* z, const double* dzdt, const std::vector<std::size_t>& species,
                        const Mixing& mixing);

// CRUCIBLE-owned stiff integration of the reaction sources (CVODES variable-order BDF,
// dense direct linear solve). The Newton matrix uses ReactionSource::jacobian (under the PaSR closure,
// times kappa_eff held over the iteration), or CVODES's own difference quotients if asked.
class ReactionStep {
 public:
  enum class Jacobian { Analytic, DifferenceQuotient };
  ReactionStep(ReactionSource& source, double rtol, double atol, Jacobian jacobian = Jacobian::Analytic);
  ~ReactionStep();
  ReactionStep(const ReactionStep&) = delete;
  ReactionStep& operator=(const ReactionStep&) = delete;

  using Observer = std::function<void(double t, const double* z)>;
  // Advance z over dt at fixed density. With samples > 1 the observer sees the state at
  // dt*i/samples, i = 1..samples, from one uninterrupted integration.
  void advance(double rho, double* z, double dt, int samples = 1, const Observer& observe = {});
  // The same under the PaSR closure: dz/dt = reactingFraction(z, f(z), species, mixing) f(z). With
  // mixing.segregation = 0 this is the integration above, bit for bit.
  void advance(double rho, double* z, double dt, const Mixing& mixing, const std::vector<std::size_t>& species,
               int samples = 1, const Observer& observe = {});
  long steps() const;  // internal steps taken by the last advance

  // One uninterrupted CVODES integration at fixed density (CVodeReInit from z, CVodeSetStopTime(tEnd),
  // then CVode(..., CV_ONE_STEP) repeatedly), showing the observer every internal step reached. Stops
  // when the observer returns false or t reaches tEnd. Returns the time reached; z returns the state
  // there. Throws if more than 2e6 internal steps are taken (docs/evidence/TABLE_B.md, "The builder,
  // stated before its code": the offline table's cold, no-time-limit branches).
  using StepObserver = std::function<bool(double t, const double* z)>;
  double follow(double rho, double* z, double tEnd, const StepObserver& observe);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace crucible::thermo
