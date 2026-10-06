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

  // dz/dt at fixed density: dY_k/dt = w_k W_k / rho, dT/dt = -sum(u_k w_k) / (rho cv).
  void rates(double rho, const double* z, double* dzdt);

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
// no transport data or uses the CHEMKIN fit form.
TransportFits transportFits(const std::string& mechanism);

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
// dense direct linear solve with a difference-quotient Jacobian).
class ReactionStep {
 public:
  ReactionStep(ReactionSource& source, double rtol, double atol);
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

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace crucible::thermo
