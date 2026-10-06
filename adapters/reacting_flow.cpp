#include "adapters/reacting_flow.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
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
    if (mine[k].name != theirs[k].name || mine[k].low != theirs[k].low || mine[k].high != theirs[k].high)
      throw std::invalid_argument("Flow medium species " + mine[k].name + " differs from the mechanism.");
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
  closure_ = true;
  cmix_ = cmix;
  closureSpecies_ = std::move(indices);
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
  flow_.pool().blocks(cells, kBlock, [&](std::size_t begin, std::size_t end, int w) {
    Worker& worker = *workers_[w];
    for (std::size_t q = begin; q < end; ++q) {
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
  for (const auto& w : workers_) w->mismatch = w->clipped = 0;
  flow_.pool().blocks(cells, kBlock, [&](std::size_t begin, std::size_t end, int w) {
    Worker& worker = *workers_[w];
    for (std::size_t q = begin; q < end; ++q) {
      const double rho = flow_.state()[q][0];
      auto y = flow_.massFractions(q);
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
  }
}

double ReactingFlow::step(double maxDt) {
  // The first half reaction step heats the gas, which can lower the flow's CFL limit below the
  // planned step. Re-plan from the saved composition until the flow takes the full step, so both
  // reaction halves span exactly half the flow step.
  // Steps end on the flow's schedule events (valves, igniter) so both halves stay symmetric.
  double planned = std::min({flow_.stableDt(), maxDt, flow_.nextEvent(flow_.time()) - flow_.time()});
  if (chemistry_ != Chemistry::FiniteRate) {
    // Equilibrium is reached after every flow step, so a split is not needed: flow, then relax.
    const double taken = flow_.step(planned);
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
  react(0.5 * taken);
  ++stats_.steps;
  return taken;
}

void ReactingFlow::advanceTo(double target) {
  if (!std::isfinite(target) || target < flow_.time()) throw std::invalid_argument("Invalid target time.");
  while (flow_.time() < target) step(target - flow_.time());
}

}  // namespace crucible::thermo
