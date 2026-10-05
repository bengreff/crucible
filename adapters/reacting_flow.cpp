#include "adapters/reacting_flow.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include "adapters/reaction.hpp"

namespace crucible::thermo {

struct ReactingFlow::Worker {
  ReactionSource source;
  ReactionStep step;
  std::vector<double> z;
  double mismatch = 0, clipped = 0;
  Worker(const std::string& mechanism, double rtol, double atol)
      : source(mechanism), step(source, rtol, atol), z(source.nSpecies() + 1) {}
};

ReactingFlow::ReactingFlow(Flow& flow, const std::string& mechanism, int threads, double rtol,
                           double atol, Chemistry chemistry)
    : flow_(flow), chemistry_(chemistry) {
  for (int i = 0; i < std::max(1, threads); ++i)
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

void ReactingFlow::react(double dt) {
  if (chemistry_ == Chemistry::Frozen) return;
  const bool equilibrium = chemistry_ == Chemistry::LocalEquilibrium;
  const std::size_t cells = flow_.state().size(), n = workers_.size();
  const std::size_t ns = flow_.medium().size();
  if (closure_) flow_.mixingInputs(cmix_, closureSpecies_, mixingTime_, segregation_);
  // Workers take small blocks of cells from a shared counter, so the hot cells of a flame or a
  // light-off do not all fall to one worker. Each cell's result is independent of which worker
  // takes it (the integrator is reinitialised per cell).
  constexpr std::size_t kBlock = 8;
  std::atomic<std::size_t> next{0};
  auto run = [&](std::size_t w) {
    Worker& worker = *workers_[w];
    worker.mismatch = 0;
    worker.clipped = 0;
    for (std::size_t start; (start = next.fetch_add(kBlock)) < cells;)
    for (std::size_t q = start; q < std::min(cells, start + kBlock); ++q) {
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
  };
  std::vector<std::thread> pool;
  for (std::size_t w = 1; w < n; ++w) pool.emplace_back(run, w);
  run(0);
  for (auto& t : pool) t.join();
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
