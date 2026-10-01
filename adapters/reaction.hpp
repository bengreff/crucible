#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

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

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

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
  long steps() const;  // internal steps taken by the last advance

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace crucible::thermo
