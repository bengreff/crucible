#include "adapters/reacting_flow.hpp"
#include "core/spray.hpp"
#include "core/checkpoint.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "adapters/reaction.hpp"

namespace crucible::thermo {

namespace {
constexpr std::size_t kBlock = 8;  // cells per block of the reaction loops
}

struct ReactingFlow::Worker {
  ReactionSource source;
  ReactionStep step;
  std::vector<double> z, y, out;
  double mismatch = 0, clipped = 0;
  // Table A tallies of one reaction call.
  long table = 0, clamped = 0, fallbacks = 0, audited = 0;
  double auditT = 0, auditY = 0;
  // Table B tallies of one reaction call, and its scratch.
  std::array<long, static_cast<std::size_t>(ManifoldRuntime::Cause::Count)> causes{};
  double snap = 0;
  long partial = 0;
  // The chemistry skip's tallies of one reaction call.
  long skipped = 0, integrated = 0;
  std::array<long, 5> frozenWeight{};
  std::vector<long> frozenNode;
  ManifoldRuntime::Work manifold;
  Worker(const std::string& mechanism, double rtol, double atol)
      : source(mechanism), step(source, rtol, atol), z(source.nSpecies() + 1), y(source.nSpecies()),
        out(source.nSpecies()) {}
};

ReactingFlow::ReactingFlow(Flow& flow, const std::string& mechanism, int threads, double rtol,
                           double atol, Chemistry chemistry)
    : flow_(flow), chemistry_(chemistry) {
  // One worker per thread of the flow's pool, which runs the reaction loops too.
  flow_.setThreads(threads);
  for (int i = 0; i < flow_.threads(); ++i)
    workers_.push_back(std::make_unique<Worker>(mechanism, rtol, atol));
  const auto& mine = flow.medium().species();
  const auto theirs = workers_[0]->source.medium().species();
  if (mine.size() != theirs.size())
    throw std::invalid_argument("Flow medium does not match the mechanism's species.");
  for (std::size_t k = 0; k < mine.size(); ++k)
    if (mine[k].name != theirs[k].name || mine[k].low != theirs[k].low || mine[k].high != theirs[k].high ||
        mine[k].tHigh != theirs[k].tHigh)
      throw std::invalid_argument("Flow medium species " + mine[k].name + " differs from the mechanism.");
  // The chemistry skip's classes. A species is reactive if some reaction changes it (a nonzero net coefficient).
  // Any reaction that changes something consumes some reactive species (mass is conserved), so in a cell with no
  // reactive species every rate is zero. H2 (O2) is skipped alone only if every reaction direction that needs no
  // reactive species but H2 (O2) changes no reactive species but H2 and H (O2 and O): its dissociation. H (O)
  // may then join H2 (O2) at a trace if no direction that needs no reactive species but the two changes another.
  const auto reactions = workers_[0]->source.reactionSpecies();
  skipClass_.assign(mine.size(), 0);
  for (const auto& r : reactions)
    for (auto k : r.changed) skipClass_[k] = kReactive;
  auto find = [&](const char* name) {
    std::size_t k = 0;
    while (k < mine.size() && mine[k].name != name) ++k;
    return k;
  };
  atom_.fill(mine.size());
  for (const auto& [molecule, atom, bit, atomBit, slot] :
       {std::tuple{"H2", "H", kFuel, kFuelAtom, 0}, std::tuple{"O2", "O", kOxidizer, kOxidizerAtom, 1}}) {
    const std::size_t m = find(molecule), a = find(atom);
    if (m == mine.size() || !skipClass_[m]) continue;
    auto within = [&](const std::vector<std::size_t>& set, std::size_t x, std::size_t y) {
      return std::all_of(set.begin(), set.end(), [&](std::size_t k) { return k == x || k == y || !skipClass_[k]; });
    };
    bool dissociation = true;
    for (const auto& r : reactions) {
      if (within(r.reactants, m, m) && !within(r.products, m, a)) dissociation = false;
      if (r.reversible && within(r.products, m, m) && !within(r.reactants, m, a)) dissociation = false;
    }
    if (!dissociation) continue;
    skipClass_[m] = bit;
    bool closed = a < mine.size();
    for (const auto& r : reactions) {
      if (within(r.reactants, m, a) && !within(r.products, m, a)) closed = false;
      if (r.reversible && within(r.products, m, a) && !within(r.reactants, m, a)) closed = false;
    }
    if (closed) {
      skipClass_[a] = atomBit;
      atom_[static_cast<std::size_t>(slot)] = a;
    }
  }
}

ReactingFlow::~ReactingFlow() = default;

void ReactingFlow::checkThreads() const {
  if (static_cast<std::size_t>(flow_.threads()) != workers_.size())
    throw std::logic_error("The flow's thread count changed after the ReactingFlow was made.");
}

void ReactingFlow::setMixingClosure(double cmix, const std::vector<std::string>& species) {
  if (chemistry_ != Chemistry::FiniteRate)
    throw std::invalid_argument("The PaSR closure needs FiniteRate chemistry.");
  if (!flow_.definition().turbulence.enabled) throw std::invalid_argument("The PaSR closure needs turbulence.");
  if (!(cmix > 0) || !std::isfinite(cmix)) throw std::invalid_argument("The PaSR mixing constant must be positive.");
  const auto& mine = flow_.medium().species();
  std::vector<std::size_t> indices;
  for (const auto& name : species) {
    std::size_t k = 0;
    while (k < mine.size() && mine[k].name != name) ++k;
    if (k == mine.size()) throw std::invalid_argument("PaSR species " + name + " is not in the mechanism.");
    indices.push_back(k);
  }
  if (indices.empty()) throw std::invalid_argument("The PaSR closure needs at least one species.");
  if (manifold_) throw std::invalid_argument("Table B under the PaSR closure is not coded yet.");
  closure_ = true;
  cmix_ = cmix;
  closureSpecies_ = std::move(indices);
}

void ReactingFlow::useManifoldTable(const ManifoldTable& table) {
  if (chemistry_ != Chemistry::FiniteRate) throw std::invalid_argument("Table B needs FiniteRate chemistry.");
  if (closure_) throw std::invalid_argument("Table B under the PaSR closure is not coded yet (TABLE_B.md, \"The runtime's step\").");
  ReactionSource& source = workers_[0]->source;
  const std::size_t ns = source.nSpecies(), ne = source.nElements();
  std::size_t eH = ne, eO = ne, eN = ne;
  for (std::size_t m = 0; m < ne; ++m) {
    const auto name = source.elementName(m);
    if (name == "H") eH = m;
    else if (name == "O") eO = m;
    else if (name == "N") eN = m;
  }
  if (eH == ne || eO == ne) throw std::invalid_argument("Table B needs H and O in the mechanism.");
  std::vector<ManifoldRuntime::Elements> elements(ns);
  std::vector<double> unit(ns, 0.0);
  for (std::size_t k = 0; k < ns; ++k) {
    unit[k] = 1;
    const auto z = source.elementMassFractions(unit.data());
    elements[k] = {z[eH], z[eO], eN < ne ? z[eN] : 0.0};
    unit[k] = 0;
  }
  manifold_ = std::make_unique<ManifoldRuntime>(table, flow_.medium(), std::move(elements));
  for (auto& w : workers_) {
    w->manifold = manifold_->work();
    w->frozenNode.assign(table.nodes.size(), 0);
  }
  stats_.manifoldFrozenNode.assign(table.nodes.size(), 0);
}

void ReactingFlow::useEquilibriumTable(const EquilibriumTable& table, long audit) {
  if (chemistry_ != Chemistry::LocalEquilibrium)
    throw std::invalid_argument("An equilibrium table needs LocalEquilibrium chemistry.");
  const auto& mine = flow_.medium().species();
  if (table.species.size() != mine.size()) throw std::invalid_argument("The table's species differ from the flow's.");
  for (std::size_t k = 0; k < mine.size(); ++k)
    if (table.species[k] != mine[k].name)
      throw std::invalid_argument("Table species " + table.species[k] + " differs from the flow's " + mine[k].name + ".");
  if (audit < 0) throw std::invalid_argument("The audit interval must not be negative.");
  table_ = &table;
  audit_ = audit;
}

// The same relaxation as react() with LocalEquilibrium, from the table: the cell's mass fractions
// (clipped and normalised as Flow::massFractions does) and its e and rho give the equilibrium Y,
// set at fixed density and total energy.
void ReactingFlow::reactTable() {
  const bool audit = audit_ > 0 && tableCalls_ % audit_ == 0;
  ++tableCalls_;
  const std::size_t cells = flow_.state().size();
  const std::size_t ns = flow_.medium().size(), nt = flow_.turbulence().size() / cells;
  checkThreads();
  for (const auto& w : workers_) {
    w->mismatch = w->clipped = w->auditT = w->auditY = 0;
    w->table = w->clamped = w->fallbacks = w->audited = 0;
  }
  const auto* list = flow_.stepCells();
  flow_.pool().blocks(list ? list->size() : cells, kBlock, [&](std::size_t begin, std::size_t end, int w) {
    Worker& worker = *workers_[w];
    for (std::size_t n = begin; n < end; ++n) {
      const std::size_t q = list ? (*list)[n] : n;
      const auto& u = flow_.state()[q];
      const double rho = u[0];
      const double* partial = flow_.partialDensities().data() + q * ns;
      double sum = 0;
      for (std::size_t k = 0; k < ns; ++k) sum += worker.y[k] = std::max(partial[k], 0.0);
      for (double& v : worker.y) v /= sum;
      const double e = (u[3] - (u[1] * u[1] + u[2] * u[2]) / (2 * rho)) / rho - (nt ? flow_.turbulence()[q * nt] / rho : 0.0);
      double t = 0;
      const auto result = table_->lookup(worker.y.data(), e, rho, worker.out.data(), t);
      const bool inside = result != EquilibriumTable::Result::Outside;
      if (!inside || audit) {
        worker.z[0] = flow_.temperature(q);
        std::copy(worker.y.begin(), worker.y.end(), worker.z.begin() + 1);
        worker.source.equilibrateUV(rho, worker.z.data());
        double total = 0;
        for (std::size_t k = 0; k < ns; ++k) {
          worker.clipped = std::max(worker.clipped, -worker.z[k + 1]);
          total += worker.z[k + 1] = std::max(worker.z[k + 1], 0.0);
        }
        for (std::size_t k = 0; k < ns; ++k) worker.z[k + 1] /= total;
      }
      if (!inside) {
        ++worker.fallbacks;
        flow_.setMassFractions(q, worker.z.data() + 1);
        worker.mismatch = std::max(worker.mismatch, std::abs(flow_.temperature(q) - worker.z[0]));
        continue;
      }
      ++worker.table;
      if (result == EquilibriumTable::Result::Clamped) ++worker.clamped;
      flow_.setMassFractions(q, worker.out.data());
      if (audit) {
        ++worker.audited;
        worker.auditT = std::max(worker.auditT, std::abs(flow_.temperature(q) - worker.z[0]));
        for (std::size_t k = 0; k < ns; ++k) worker.auditY = std::max(worker.auditY, std::abs(worker.out[k] - worker.z[k + 1]));
      }
    }
  });
  for (const auto& w : workers_) {
    stats_.maxTemperatureMismatch = std::max(stats_.maxTemperatureMismatch, w->mismatch);
    stats_.maxClippedFraction = std::max(stats_.maxClippedFraction, w->clipped);
    stats_.tableCells += w->table;
    stats_.tableClamped += w->clamped;
    stats_.tableFallbacks += w->fallbacks;
    stats_.auditedCells += w->audited;
    stats_.auditTemperature = std::max(stats_.auditTemperature, w->auditT);
    stats_.auditMassFraction = std::max(stats_.auditMassFraction, w->auditY);
  }
}

void ReactingFlow::react(double dt) {
  if (chemistry_ == Chemistry::Frozen) return;
  const auto clock0 = std::chrono::steady_clock::now();
  struct Timer {
    double& total;
    std::chrono::steady_clock::time_point start;
    ~Timer() { total += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); }
  } timer{stats_.reactWall, clock0};
  if (table_) {
    reactTable();
    return;
  }
  const bool equilibrium = chemistry_ == Chemistry::LocalEquilibrium;
  const std::size_t cells = flow_.state().size();
  const std::size_t ns = flow_.medium().size();
  checkThreads();
  if (closure_) flow_.mixingInputs(cmix_, closureSpecies_, mixingTime_, segregation_);
  // Workers take small blocks of cells from a shared counter, so the hot cells of a flame or a
  // light-off do not all fall to one worker. Each cell's result is independent of which worker
  // takes it (the integrator is reinitialised per cell).
  for (const auto& w : workers_) {
    w->mismatch = w->clipped = w->snap = 0;
    w->partial = w->skipped = w->integrated = 0;
    w->causes.fill(0);
    w->frozenWeight.fill(0);
    std::fill(w->frozenNode.begin(), w->frozenNode.end(), 0);
  }
  const std::size_t nt = manifold_ ? flow_.turbulence().size() / cells : 0;
  constexpr std::size_t kCauses = static_cast<std::size_t>(ManifoldRuntime::Cause::Count);
  if (manifold_ && stats_.manifoldCellCauses.size() != cells * kCauses) stats_.manifoldCellCauses.assign(cells * kCauses, 0);
  // Under AMR step B2 the cells a level steps (Flow::stepCells), else every cell.
  const auto* list = flow_.stepCells();
  flow_.pool().blocks(list ? list->size() : cells, kBlock, [&](std::size_t begin, std::size_t end, int w) {
    Worker& worker = *workers_[w];
    for (std::size_t n = begin; n < end; ++n) {
      const std::size_t q = list ? (*list)[n] : n;
      const double rho = flow_.state()[q][0];
      auto y = flow_.massFractions(q);
      if (manifold_) {
        // Table B first; CVODES below for a cell it does not handle. e as reactTable takes it.
        const auto& u = flow_.state()[q];
        const double e = (u[3] - (u[1] * u[1] + u[2] * u[2]) / (2 * rho)) / rho - (nt ? flow_.turbulence()[q * nt] / rho : 0.0);
        const auto placed = manifold_->advance(e, rho, flow_.temperature(q), y.data(), dt, worker.manifold);
        ++worker.causes[static_cast<std::size_t>(placed.cause)];
        ++stats_.manifoldCellCauses[q * kCauses + static_cast<std::size_t>(placed.cause)];  // one worker per cell
        if (placed.cause == ManifoldRuntime::Cause::FrozenCorner) {
          const double fw = placed.frozenW;
          ++worker.frozenWeight[fw < 1e-30 ? 0 : fw < 1e-12 ? 1 : fw < 1e-6 ? 2 : fw < 1e-2 ? 3 : 4];
          ++worker.frozenNode[placed.frozenNode];
        }
        if (placed.cause == ManifoldRuntime::Cause::Table) {
          worker.snap = std::max(worker.snap, placed.dY);
          if (placed.w < 1) ++worker.partial;
          double sum = 0;
          for (double v : y) sum += v;
          for (double& v : y) v /= sum;
          flow_.setMassFractions(q, y.data());
          continue;
        }
      }
      if (skip_ && !equilibrium) {
        // The chemistry skip (setChemistrySkip): a negative partial density counts as present, so the integrator's
        // path still clips it. Read from the partial densities, since massFractions has already clipped it to 0.
        const double* partial = flow_.partialDensities().data() + q * ns;
        unsigned mask = 0;
        for (std::size_t k = 0; k < ns; ++k)
          if (partial[k] != 0) mask |= skipClass_[k];
        // The molecule, with at most a positive trace of its atom.
        auto alone = [&](unsigned molecule, unsigned atom, std::size_t a) {
          return (mask & (molecule | atom)) && !(mask & ~(molecule | atom)) &&
                 (!(mask & atom) || (partial[a] > 0 && y[a] <= kAtomTrace));
        };
        if (mask == 0 || ((alone(kFuel, kFuelAtom, atom_[0]) || alone(kOxidizer, kOxidizerAtom, atom_[1])) &&
                          flow_.temperature(q) <= kSkipTemperature)) {
          ++worker.skipped;
          continue;
        }
      }
      if (!equilibrium) ++worker.integrated;
      worker.z[0] = flow_.temperature(q);
      std::copy(y.begin(), y.end(), worker.z.begin() + 1);
      if (equilibrium) worker.source.equilibrateUV(rho, worker.z.data());
      else if (closure_) worker.step.advance(rho, worker.z.data(), dt, {mixingTime_[q], segregation_[q]}, closureSpecies_);
      else worker.step.advance(rho, worker.z.data(), dt);
      double sum = 0;
      for (std::size_t k = 0; k < ns; ++k) {
        worker.clipped = std::max(worker.clipped, -worker.z[k + 1]);
        y[k] = std::max(worker.z[k + 1], 0.0);
        sum += y[k];
      }
      for (double& v : y) v /= sum;
      flow_.setMassFractions(q, y.data());
      worker.mismatch = std::max(worker.mismatch, std::abs(flow_.temperature(q) - worker.z[0]));
    }
  });
  for (const auto& w : workers_) {
    stats_.maxTemperatureMismatch = std::max(stats_.maxTemperatureMismatch, w->mismatch);
    stats_.maxClippedFraction = std::max(stats_.maxClippedFraction, w->clipped);
    for (std::size_t c = 0; c < w->causes.size(); ++c) stats_.manifoldCauses[c] += w->causes[c];
    stats_.manifoldMaxSnap = std::max(stats_.manifoldMaxSnap, w->snap);
    stats_.manifoldPartial += w->partial;
    for (std::size_t c = 0; c < w->frozenWeight.size(); ++c) stats_.manifoldFrozenWeight[c] += w->frozenWeight[c];
    for (std::size_t n = 0; n < w->frozenNode.size(); ++n) stats_.manifoldFrozenNode[n] += w->frozenNode[n];
    stats_.skippedCells += w->skipped;
    stats_.integratedCells += w->integrated;
  }
}

double ReactingFlow::step(double maxDt) {
  // The first half reaction step heats the gas, which can lower the flow's CFL limit below the
  // planned step. Re-plan from the saved composition until the flow takes the full step, so both
  // reaction halves span exactly half the flow step.
  // Steps end on the flow's schedule events (valves, igniter) so both halves stay symmetric.
  double planned = std::min({flow_.stableDt(), maxDt, flow_.nextEvent(flow_.time()) - flow_.time()});
  if (spray_) planned = std::min(planned, spray_->nextEvent(flow_.time()) - flow_.time());
  if (chemistry_ != Chemistry::FiniteRate) {
    // Equilibrium is reached after every flow step, so a split is not needed: flow, then relax.
    const double taken = flow_.step(planned);
    if (spray_) spray_->advance(taken);
    react(taken);
    ++stats_.steps;
    return taken;
  }
  saved_ = flow_.partialDensities();
  for (int attempt = 0;; ++attempt) {
    react(0.5 * planned);
    const double limit = flow_.stableDt();
    if (limit >= planned) break;
    if (attempt == 20) throw std::runtime_error("Reaction step keeps lowering the flow time step.");
    ++stats_.replans;
    flow_.setPartialDensities(saved_);
    planned = limit;
  }
  const double taken = flow_.step(planned);
  if (taken < planned) ++stats_.asymmetricSteps;  // admissibility rejection inside the flow step
  if (spray_) spray_->advance(taken);
  react(0.5 * taken);
  ++stats_.steps;
  return taken;
}

double ReactingFlow::implicitStep(double dt, const Flow::ImplicitControl& control) {
  if (chemistry_ == Chemistry::FiniteRate || spray_)
    throw std::runtime_error("The implicit prototype takes equilibrium or frozen chemistry and no spray.");
  const double taken = flow_.implicitStep(dt, control);
  react(taken);
  ++stats_.steps;
  return taken;
}

void ReactingFlow::advanceTo(double target) {
  if (!std::isfinite(target) || target < flow_.time()) throw std::invalid_argument("Invalid target time.");
  while (flow_.time() < target) step(target - flow_.time());
}

// The checkpoint's fields (docs/evidence/RESTART.md); the wall time is not kept.
template <class Self, class Io> void ReactingFlow::checkpoint(Self& self, Io& io) {
  io.mark("reacting");
  io.same(static_cast<int>(self.chemistry_), "the chemistry");
  auto& st = self.stats_;
  io.field(st.steps), io.field(st.asymmetricSteps), io.field(st.replans), io.field(st.maxTemperatureMismatch);
  io.field(st.maxClippedFraction), io.field(st.tableCells), io.field(st.tableClamped), io.field(st.tableFallbacks);
  io.field(st.auditedCells), io.field(st.auditTemperature), io.field(st.auditMassFraction);
  io.field(st.manifoldCauses), io.field(st.manifoldMaxSnap), io.field(st.manifoldPartial), io.field(st.manifoldCellCauses);
  io.field(st.manifoldFrozenWeight), io.field(st.manifoldFrozenNode);
  io.field(st.skippedCells), io.field(st.integratedCells);
  io.field(self.tableCalls_);
}
void ReactingFlow::writeCheckpoint(CheckpointWriter& out) const { checkpoint(*this, out); }
void ReactingFlow::readCheckpoint(CheckpointReader& in) { checkpoint(*this, in); }

}  // namespace crucible::thermo
