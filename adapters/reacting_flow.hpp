#pragma once

#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "core/equilibrium_table.hpp"
#include "core/flow.hpp"
#include "core/manifold_runtime.hpp"

namespace crucible { class Spray; }

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
    // Table B (useManifoldTable): cell updates by cause (ManifoldRuntime::Cause; Table is the table's own,
    // every other cause took CVODES), the largest max|dY| of a placed point the table stepped, and the table's
    // updates placed at w < 1 (in part by s).
    std::array<long, static_cast<std::size_t>(ManifoldRuntime::Cause::Count)> manifoldCauses{};
    double manifoldMaxSnap = 0;
    long manifoldPartial = 0;
    // The same counts per cell: [q * Cause::Count + cause], q the flow's cell index (zeroed if the cell count changes).
    std::vector<long> manifoldCellCauses;
    // FrozenCorner updates by the weight of the heaviest refusing corner, by decade: below 1e-30, to 1e-12, to
    // 1e-6, to 1e-2, and above; and by that corner's table node (sized to the table's nodes).
    std::array<long, 5> manifoldFrozenWeight{};
    std::vector<long> manifoldFrozenNode;
    // FiniteRate cell updates the chemistry skip left unchanged (setChemistrySkip), and those CVODES integrated.
    long skippedCells = 0, integratedCells = 0;
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
  // Table B (TECHNICAL_PLAN, *Lightweight engine*; docs/evidence/TABLE_B.md, "The runtime's step on the placed
  // point"): with FiniteRate and no PaSR closure, each cell is placed on the table's blended curve and advanced
  // along it; a cell the table does not handle (counted by cause in the stats) is integrated by CVODES as
  // without the table. The table must be built with the specific entropy from the same mechanism, and outlive
  // this object.
  void useManifoldTable(const ManifoldTable& table);
  void react(double dt);  // one reaction substep over every cell (equilibrium: dt is unused)
  // The chemistry skip (docs/evidence/SPEED.md, "The chemistry skip"), on unless set off: with FiniteRate, a
  // cell whose species include no reactive one (one some reaction changes) is left unchanged, since every rate
  // is zero there; so is a cell at or below kSkipTemperature whose one reactive species is H2, or O2 (no
  // oxidizer, or no fuel). The one reaction such a cell can start is the molecule's dissociation into its atoms
  // (checked on the mechanism at construction; if it fails, that molecule is not skipped), whose equilibrium
  // bounds the change at any step length (tests/chemistry_skip_tests.cpp). A partial density that is not exactly
  // zero, negative included, counts as present. The atom may be present too, at 0 < Y <= kAtomTrace, if the mechanism
  // has no reaction direction that needs no reactive species but the molecule and the atom and changes another
  // (dissociation and recombination only): the cell then moves only toward the dissociation equilibrium, so its
  // change at any step length is at most the larger of Y_atom and the atom's equilibrium fraction. The atoms
  // come from CVODES itself, which leaves the H2 it integrates with a positive trace of H.
  static constexpr double kSkipTemperature = 500, kAtomTrace = 1e-14;
  void setChemistrySkip(bool on) { skip_ = on; }
  // A liquid spray (core/spray.hpp) inside the shared step: it moves after the flow step with the
  // same dt, before the second reaction half (docs/evidence/SPRAY.md, Coupling). Null: none.
  void setSpray(crucible::Spray* spray) { spray_ = spray; }
  Chemistry chemistry() const { return chemistry_; }
  const Stats& stats() const { return stats_; }
  // A restart checkpoint (core/checkpoint.hpp, docs/evidence/RESTART.md): the stats and the audit's call
  // count. Each reaction call starts every cell afresh, so nothing else carries from one step to the next.
  void writeCheckpoint(CheckpointWriter& out) const;
  void readCheckpoint(CheckpointReader& in);

 private:
  struct Worker;
  void checkThreads() const;
  Flow& flow_;
  crucible::Spray* spray_ = nullptr;
  std::vector<std::unique_ptr<Worker>> workers_;
  Stats stats_;
  Chemistry chemistry_;
  std::vector<double> saved_;
  bool closure_ = false;
  double cmix_ = 1;
  std::vector<std::size_t> closureSpecies_;
  std::vector<double> mixingTime_, segregation_;
  const EquilibriumTable* table_ = nullptr;
  std::unique_ptr<ManifoldRuntime> manifold_;
  long audit_ = 0, tableCalls_ = 0;
  // The chemistry skip's species classes, bits of a cell's mask: 0 a species no reaction has, then H2 and O2
  // where the mechanism's check passed, then every other reactive species, then H and O where the closure check
  // passed (their indices in atom_, else the species count).
  enum : unsigned char { kFuel = 1, kOxidizer = 2, kReactive = 4, kFuelAtom = 8, kOxidizerAtom = 16 };
  std::vector<unsigned char> skipClass_;
  std::array<std::size_t, 2> atom_{};
  bool skip_ = true;
  void reactTable();
  template <class Self, class Io> static void checkpoint(Self& self, Io& io);
};

}  // namespace crucible::thermo
