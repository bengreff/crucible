#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "core/flow.hpp"

namespace crucible::thermo {

// Reacting flow by symmetric (Strang) splitting: half reaction step, one flow step, half reaction
// step. Reaction acts per cell at fixed density and total energy (closed adiabatic constant-volume
// parcel, CVODES BDF over Cantera rates). The flow's medium must be the mechanism's species in the
// mechanism's order (ReactionSource::medium()). Each worker thread owns its Cantera context.
//
// Chemistry selects the reaction model applied in every cell alike:
// FiniteRate integrates the mechanism's kinetics (the physical model).
// LocalEquilibrium relaxes each cell to constant-(u, v) equilibrium after every flow step: the
// infinitely fast limit, used to verify against CEA. It burns any cold premixed gas at once, so an
// igniter has no role in it.
// Frozen leaves composition unchanged.
enum class Chemistry { FiniteRate, LocalEquilibrium, Frozen };

class ReactingFlow {
 public:
  struct Stats {
    long steps = 0;
    long asymmetricSteps = 0;  // flow step came back shorter than planned (step rejection)
    long replans = 0;          // first half reaction redone because it lowered the CFL step
    double maxTemperatureMismatch = 0;  // |T_chemistry - T_core| after reaction substeps [K]
  };
  ReactingFlow(Flow& flow, const std::string& mechanism, int threads, double rtol = 1e-8,
               double atol = 1e-14, Chemistry chemistry = Chemistry::FiniteRate);
  ~ReactingFlow();
  ReactingFlow(const ReactingFlow&) = delete;
  ReactingFlow& operator=(const ReactingFlow&) = delete;

  // The PaSR turbulence-chemistry closure (TECHNICAL_PLAN step 7, "PaSR") in every cell alike: each
  // reaction substep integrates dz/dt = kappa_eff(z) f(z) with the cell's mixing time and
  // segregation frozen at the substep's start (Flow::mixingInputs, mixing constant cmix) and the
  // chemical time over the named species. Needs FiniteRate chemistry and a flow with turbulence.
  void setMixingClosure(double cmix, const std::vector<std::string>& species);
  // The closure's inputs of the last reaction substep, per cell (empty without the closure).
  const std::vector<double>& mixingTimes() const { return mixingTime_; }
  const std::vector<double>& segregations() const { return segregation_; }
  double step(double maxDt = std::numeric_limits<double>::infinity());
  void advanceTo(double time);
  void react(double dt);  // one reaction substep over every cell (equilibrium: dt is unused)
  Chemistry chemistry() const { return chemistry_; }
  const Stats& stats() const { return stats_; }

 private:
  struct Worker;
  Flow& flow_;
  std::vector<std::unique_ptr<Worker>> workers_;
  Stats stats_;
  Chemistry chemistry_;
  std::vector<double> saved_;
  bool closure_ = false;
  double cmix_ = 1;
  std::vector<std::size_t> closureSpecies_;
  std::vector<double> mixingTime_, segregation_;
};

}  // namespace crucible::thermo
