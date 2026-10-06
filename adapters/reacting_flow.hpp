#pragma once

#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "core/equilibrium_table.hpp"
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
    double maxClippedFraction = 0;      // largest -Y_k set to zero after a reaction substep
    // Equilibrium table (useEquilibriumTable): cells taken from it, those clamped to its low-energy
    // edge, and those outside it that took the direct call. The audit's cells, and its largest
    // |T_core(table Y) - T_Cantera| [K] and |Y_table - Y_Cantera| over those inside the table.
    long tableCells = 0, tableClamped = 0, tableFallbacks = 0, auditedCells = 0;
    double auditTemperature = 0, auditMassFraction = 0;
    double reactWall = 0;  // wall time inside react() [s]
  };
  // Sets the flow's thread count (Flow::setThreads), whose pool then runs the reaction loops too;
  // it must not change while this object lives.
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
  // Table A (TECHNICAL_PLAN, *Lightweight engine*; docs/evidence/TABLE_A.md): with LocalEquilibrium,
  // each cell's equilibrium comes from the table, and cells outside it take the direct Cantera call.
  // Every audit-th reaction call (0: never) also makes the direct call in every cell and records the
  // difference in the stats; the table's result is kept. The table must list the mechanism's species
  // in order and outlive this object.
  void useEquilibriumTable(const EquilibriumTable& table, long audit = 0);
  void react(double dt);  // one reaction substep over every cell (equilibrium: dt is unused)
  Chemistry chemistry() const { return chemistry_; }
  const Stats& stats() const { return stats_; }

 private:
  struct Worker;
  void checkThreads() const;
  Flow& flow_;
  std::vector<std::unique_ptr<Worker>> workers_;
  Stats stats_;
  Chemistry chemistry_;
  std::vector<double> saved_;
  bool closure_ = false;
  double cmix_ = 1;
  std::vector<std::size_t> closureSpecies_;
  std::vector<double> mixingTime_, segregation_;
  const EquilibriumTable* table_ = nullptr;
  long audit_ = 0, tableCalls_ = 0;
  void reactTable();
};

}  // namespace crucible::thermo
