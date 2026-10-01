#include "adapters/reacting_flow.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <thread>
#include <vector>

#include "adapters/reaction.hpp"

namespace crucible::thermo {

struct ReactingFlow::Worker {
  ReactionSource source;
  ReactionStep step;
  std::vector<double> z;
  double mismatch = 0;
  Worker(const std::string& mechanism, double rtol, double atol)
      : source(mechanism), step(source, rtol, atol), z(source.nSpecies() + 1) {}
};

ReactingFlow::ReactingFlow(Flow& flow, const std::string& mechanism, int threads, double rtol,
                           double atol)
    : flow_(flow) {
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

void ReactingFlow::react(double dt) {
  const std::size_t cells = flow_.state().size(), n = workers_.size();
  const std::size_t ns = flow_.medium().size();
  // Contiguous cell blocks per worker; each cell's result is independent of the partition.
  auto run = [&](std::size_t w) {
    Worker& worker = *workers_[w];
    worker.mismatch = 0;
    for (std::size_t q = cells * w / n; q < cells * (w + 1) / n; ++q) {
      const double rho = flow_.state()[q][0];
      auto y = flow_.massFractions(q);
      worker.z[0] = flow_.temperature(q);
      std::copy(y.begin(), y.end(), worker.z.begin() + 1);
      worker.step.advance(rho, worker.z.data(), dt);
      double sum = 0;
      for (std::size_t k = 0; k < ns; ++k) { y[k] = std::max(worker.z[k + 1], 0.0); sum += y[k]; }
      for (double& v : y) v /= sum;
      flow_.setMassFractions(q, y.data());
      worker.mismatch = std::max(worker.mismatch, std::abs(flow_.temperature(q) - worker.z[0]));
    }
  };
  std::vector<std::thread> pool;
  for (std::size_t w = 1; w < n; ++w) pool.emplace_back(run, w);
  run(0);
  for (auto& t : pool) t.join();
  for (const auto& w : workers_) stats_.maxTemperatureMismatch = std::max(stats_.maxTemperatureMismatch, w->mismatch);
}

double ReactingFlow::step(double maxDt) {
  // The first half reaction step heats the gas, which can lower the flow's CFL limit below the
  // planned step. Re-plan from the saved composition until the flow takes the full step, so both
  // reaction halves span exactly half the flow step.
  double planned = std::min(flow_.stableDt(), maxDt);
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
