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
class ReactingFlow {
 public:
  struct Stats {
    long steps = 0;
    long asymmetricSteps = 0;  // flow step came back shorter than planned (step rejection)
    long replans = 0;          // first half reaction redone because it lowered the CFL step
    double maxTemperatureMismatch = 0;  // |T_chemistry - T_core| after reaction substeps [K]
  };
  ReactingFlow(Flow& flow, const std::string& mechanism, int threads, double rtol = 1e-8,
               double atol = 1e-14);
  ~ReactingFlow();
  ReactingFlow(const ReactingFlow&) = delete;
  ReactingFlow& operator=(const ReactingFlow&) = delete;

  double step(double maxDt = std::numeric_limits<double>::infinity());
  void advanceTo(double time);
  void react(double dt);  // one reaction substep over every cell
  const Stats& stats() const { return stats_; }

 private:
  struct Worker;
  Flow& flow_;
  std::vector<std::unique_ptr<Worker>> workers_;
  Stats stats_;
  std::vector<double> saved_;
};

}  // namespace crucible::thermo
