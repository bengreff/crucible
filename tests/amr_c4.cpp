// AMR: fixing the tag thresholds, then criterion 4 (docs/evidence/AMR.md, *Fixing the tag thresholds, then
// criterion 4*, recorded before this file was written). C1 eqtt with the wall functions (law) and the property
// corrections (hp), to 8 ms.
// Usage: crucible_amr_c4 amr scale=S [end time s, default 8e-3] [threads, default 4]: 32 x 6 plus one level, set
//          up as crucible_amr_budget's tagged case, with step C's four dynamic thresholds times S;
//        crucible_amr_c4 uniform [end] [threads]: the uniform 64 x 12 run.
// CRUCIBLE_EQ_TABLE names Table A; CRUCIBLE_C4_EXIT_CHECK=<t s> runs the stagnation loop's exit check from t and
// ends the run there. Each run prints one line starting "c4 " that tests/amr_c4_judge.sh reads:
//   - settled c*, vacuum Isp and vacuum thrust: the means of 21 samples over the last millisecond, each at the
//     first step ending at or after 7.00, 7.05, ..., 8.00 ms;
//   - c*: the mass-flux-weighted equilibrium stagnation pressure on the fine chamber-end column (on AMR, level
//     1's cells where it is active over a base row, else the base cell) times A_t over the outlet mass flow;
//   - vacuum thrust: inlet momentum flux plus wall and body axial forces, summed over the levels; Isp = F / (mdot g0);
//   - light-off: the end of the first step (base step on AMR) with the largest Y_H2O over the owned cells above 0.5;
//   - t90 (AMR.md, *Criterion 4: the Director's additions*; reported, not judged): the chamber pressure (c*'s
//     numerator) sampled at the first step ending at or after each multiple of 10 us and at the end, none where the
//     column's mass flux is not positive; after the last sample below 0.9 of the end value, the linear crossing.
// Printed, not judged: each quantity at the end, level 1's active cells, regrids, budgets, and for criterion 5(a)
// the cell updates (each level's active cells times its steps) and the wall time per update, net of the history's
// sampling, with the load.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "cantera/core.h"
#include "core/amr.hpp"
#include "core/equilibrium_table.hpp"
#include "tests/c1_eqtt.hpp"

using namespace crucible;
using namespace crucible::c1_eqtt;

namespace {
constexpr double kG0 = 9.80665, kWindow = 1e-3, kSpacing = 5e-5, kHistory = 1e-5;
constexpr int kSamples = 21;
double seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
double load() {
  double l[1] = {0};
  return getloadavg(l, 1) == 1 ? l[0] : -1;
}
struct Sample {
  double cstar, isp, fVac;
};
// The largest Y_H2O over a level's cells, the active ones only if a mask is given.
double largestWater(const Flow& f, std::size_t water, const std::vector<unsigned char>* active) {
  double y = 0;
  for (std::size_t q = 0; q < f.mesh().cells.size(); ++q)
    if (!active || (*active)[q]) y = std::max(y, f.massFractions(q)[water]);
  return y;
}
}  // namespace

int main(int argc, char** argv) {
  double end = 8e-3, scale = 1;
  int threads = 4, positional = 0;
  bool uniform = false, amrRun = false;
  for (int k = 1; k < argc; ++k) {
    const std::string a = argv[k];
    if (a == "uniform") uniform = true;
    else if (a == "amr") amrRun = true;
    else if (a.rfind("scale=", 0) == 0) scale = std::atof(a.c_str() + 6);
    else if (positional++ == 0) end = std::atof(a.c_str());
    else threads = std::max(1, std::atoi(a.c_str()));
  }
  if (uniform == amrRun || !(scale > 0)) throw std::invalid_argument("crucible_amr_c4: amr scale=S, or uniform");
  const char* path = std::getenv("CRUCIBLE_EQ_TABLE");
  if (!path) throw std::runtime_error("crucible_amr_c4 needs CRUCIBLE_EQ_TABLE (Table A's file)");
  const EquilibriumTable table = EquilibriumTable::load(path);
  thermo::ReactionSource source("h2o2.yaml");
  Definition d = c1eqtt(source, uniform ? 64 : 32, uniform ? 12 : 6);
  d.turbulence.wallFunctions = d.turbulence.propertyCorrections = true;
  std::size_t water = d.species.size();
  for (std::size_t k = 0; k < d.species.size(); ++k)
    if (d.species[k].name == "H2O") water = k;
  if (water == d.species.size()) throw std::runtime_error("crucible_amr_c4: no H2O");
  auto historySolution = Cantera::newSolution("h2o2.yaml", "", "none");
  Cantera::ThermoPhase& gas = *historySolution->thermo();
  std::printf("AMR criterion 4: C1 eqtt with law and hp, %s, to %.4g ms, %d threads; table %s\n",
              uniform ? "uniform 64 x 12" : ("32 x 6 plus one level, dynamic thresholds x " + std::to_string(scale)).c_str(),
              end * 1e3, threads, table.provenance.c_str());

  std::unique_ptr<Amr> amr;
  std::unique_ptr<Flow> single;
  std::vector<std::unique_ptr<thermo::ReactingFlow>> reacting;
  std::uint64_t updates = 0, activeSum = 0;
  std::vector<std::uint64_t> levelCells;
  if (uniform) {
    single = std::make_unique<Flow>(d);
    reacting.push_back(std::make_unique<thermo::ReactingFlow>(*single, "h2o2.yaml", threads, 1e-6, 1e-12,
                                                              thermo::Chemistry::LocalEquilibrium));
    reacting.back()->useEquilibriumTable(table, 0);
  } else {
    // As crucible_amr_budget 8e-3 4 wall tags (step D's D1): the fixed pattern first, then the tags.
    amr = std::make_unique<Amr>(d, 1);
    amr->setRestricted(true);
    for (int l = 0; l < amr->levels(); ++l) {
      reacting.push_back(std::make_unique<thermo::ReactingFlow>(amr->level(l), "h2o2.yaml", threads, 1e-6, 1e-12,
                                                                thermo::Chemistry::LocalEquilibrium));
      reacting.back()->useEquilibriumTable(table, 0);
    }
    levelCells.assign(static_cast<std::size_t>(amr->levels()), 0);
    amr->setAfterStep([&](int l, double dt) {
      if (l == 0)
        for (int k = 0; k < amr->levels(); ++k) {
          // The base's mask is empty: it is active everywhere.
          std::uint64_t n = amr->active(k).empty() ? amr->level(k).mesh().cells.size() : 0;
          for (auto a : amr->active(k)) n += a;
          levelCells[static_cast<std::size_t>(k)] = n;
        }
      updates += levelCells[static_cast<std::size_t>(l)];
      reacting[static_cast<std::size_t>(l)]->react(dt);
    });
    const Mesh& fine = amr->level(1).mesh();
    std::vector<unsigned char> active(fine.cells.size(), 0);
    for (int i = 0; i < fine.nz; ++i)
      for (int j = 0; j < fine.nr; ++j) active[fine.index(i, j)] = ((i >> 1) < 16 && (j >> 1) <= 1) || (j >> 1) == 5;
    amr->setActive(1, active);
    amr->averageDown();
    Amr::Tags tags;
    const std::size_t eH = gas.elementIndex("H");
    for (const auto& sp : d.species) {
      const std::size_t k = gas.speciesIndex(sp.name);
      tags.hydrogen.push_back(gas.nAtoms(k, eH) * gas.atomicWeight(eH) / gas.molecularWeight(k));
    }
    // S scales step C's thresholds, as registered, not Tags' defaults (now s* times them).
    tags.mixture = 0.02 * scale;
    tags.temperature = 0.05 * scale;
    tags.pressure = 0.10 * scale;
    tags.vorticity = 0.05 * scale;
    std::printf("thresholds: f_H %.6g, T %.6g, pressure %.6g, vorticity %.6g; interval %d, buffer %d, wall y+ %.6g, throat %.3g r_t to level %d\n",
                tags.mixture, tags.temperature, tags.pressure, tags.vorticity, tags.interval, tags.buffer, tags.wallYPlus,
                tags.throatHalfWidth, tags.throatLevel);
    amr->setTags(tags);
  }
  // The finest mesh: its throat and its chamber-end column.
  const Mesh& finest = uniform ? single->mesh() : amr->level(1).mesh();
  int throat = 0, fe = 0;
  for (int i = 0; i <= finest.nz; ++i)
    if (finest.radius[i] < finest.radius[throat]) throat = i;
  for (int i = 0; i < finest.nz; ++i)
    if (finest.cells[finest.index(i, 0)].z < geometry().cylinderEnd) fe = i;
  const double at = std::numbers::pi * finest.radius[throat] * finest.radius[throat];
  std::printf("fine throat station %d, fine chamber-end column %d (%s)\n", throat, fe, fe % 2 ? "odd: the base column ends with it" : "even");

  auto time = [&]() { return uniform ? single->time() : amr->time(); };
  auto flows = [&]() {
    std::vector<const Flow*> out;
    if (uniform) out.push_back(single.get());
    else
      for (int l = 0; l < amr->levels(); ++l) out.push_back(&amr->level(l));
    return out;
  };
  // The chamber-end column's flux-weighted stagnation pressure and its mass flux (columnStagnation's loop).
  // `oldExit` and `counts` (each cell's SP solves) are for the exit check below.
  auto column = [&](double& p0, double& flux, bool oldExit = false, std::vector<int>* counts = nullptr) {
    p0 = flux = 0;
    int n = 0;
    auto add = [&](const Flow& f, int i, int j) {
      addStagnation(f, gas, i, j, p0, flux, oldExit, &n);
      if (counts) counts->push_back(n);
    };
    if (uniform)
      for (int j = 0; j < finest.nr; ++j) add(*single, fe, j);
    else {
      const Flow& base = amr->level(0);
      const Flow& fineFlow = amr->level(1);
      for (int J = 0; J < base.mesh().nr; ++J)
        if (amr->active(1)[fineFlow.mesh().index(fe, 2 * J)]) {
          add(fineFlow, fe, 2 * J);
          add(fineFlow, fe, 2 * J + 1);
        } else
          add(base, fe >> 1, J);
    }
  };
  auto sample = [&]() {
    double mdot = 0, fVac = 0;
    for (const Flow* f : flows()) {
      const auto m = f->measurements();
      mdot += m.outletMassFlow;
      fVac += m.inletMomentumFlux + m.wallAxialForce + m.bodyAxialForce;
    }
    double p0, flux;
    column(p0, flux);
    return Sample{p0 / flux * at / mdot, fVac / (mdot * kG0), fVac};
  };
  // The history's chamber pressure: none (NaN) where the column carries no net outflow.
  std::vector<std::pair<double, double>> history;
  double historyWall = 0;
  // The amendment's check (AMR.md, *Amendment: the stagnation loop exits at the equilibrium solver's accuracy*):
  // from the first record at or after CRUCIBLE_C4_EXIT_CHECK seconds, the column's p0 with the new exit and the
  // old; the run ends at the first such record where some cell runs the old exit to 50 solves.
  const char* checkEnv = std::getenv("CRUCIBLE_C4_EXIT_CHECK");
  const double checkFrom = checkEnv ? std::atof(checkEnv) : -1;
  auto exitCheck = [&](double t) {
    double pNew, pOld, flux;
    std::vector<int> nNew, nOld;
    column(pNew, flux, false, &nNew);
    column(pOld, flux, true, &nOld);
    if (!(flux > 0)) return;
    pNew /= flux;
    pOld /= flux;
    std::printf("exit check at t %.6f ms, SP solves per cell (new / old):", t * 1e3);
    bool exercised = false;
    for (std::size_t k = 0; k < nNew.size(); ++k) {
      std::printf(" %d/%d", nNew[k], nOld[k]);
      exercised = exercised || nOld[k] == 50;
    }
    std::printf("\n");
    if (!exercised) return;
    const double half = 0.5 * std::pow(10.0, std::floor(std::log10(std::abs(pOld))) - 5);
    const double cstar = sample().cstar;
    const bool pass = std::abs(pNew - pOld) < half;
    std::printf("exit check: column p0 new %.12e Pa, old %.12e Pa, |difference| %.3e Pa (%.3e relative); half a unit in "
                "the 6th digit %.3e Pa: %s. Not judged: c* %.4f m/s, its printed half unit %.3e relative\n",
                pNew, pOld, std::abs(pNew - pOld), std::abs(pNew - pOld) / pOld, half, pass ? "PASS" : "FAIL", cstar,
                5e-5 / cstar);
    std::exit(pass ? 0 : 1);
  };
  auto record = [&](double t) {
    const auto start = std::chrono::steady_clock::now();
    double p0, flux;
    column(p0, flux);
    history.emplace_back(t, flux > 0 ? p0 / flux : std::nan(""));
    historyWall += seconds(start);
    if (checkFrom >= 0 && t >= checkFrom) exitCheck(t);
  };
  auto largestY = [&]() {
    if (uniform) return largestWater(*single, water, nullptr);
    // Covered base cells hold their children's mean, so every cell of the base and level 1's active cells
    // give the owned cells' largest.
    return std::max(largestWater(amr->level(0), water, nullptr), largestWater(amr->level(1), water, &amr->active(1)));
  };

  std::vector<Sample> samples;
  const double window = end - kWindow;
  double lightOff = -1, lightStep = 0;
  const double load0 = load();
  const auto clock0 = std::chrono::steady_clock::now();
  double nextPrint = 5e-4;
  while (time() < end) {
    const double t0 = time();
    if (uniform) reacting[0]->step(end - t0);
    else {
      amr->step(end - t0);
      std::uint64_t n = 0;
      for (auto a : amr->active(1)) n += a;
      activeSum += n;
    }
    const double t1 = time();
    if (lightOff < 0 && largestY() > 0.5) {
      lightOff = t1;
      lightStep = t1 - t0;
      std::printf("light-off at %.6f ms, step %.3e s\n", lightOff * 1e3, lightStep);
    }
    if (static_cast<int>(samples.size()) < kSamples - 1 && t1 >= window + kSpacing * static_cast<double>(samples.size()))
      samples.push_back(sample());
    if (t1 >= kHistory * static_cast<double>(history.size())) record(t1);
    if (t1 >= nextPrint) {
      std::printf("  t %.4f ms, %.0f s\n", t1 * 1e3, seconds(clock0));
      std::fflush(stdout);
      nextPrint += 5e-4;
    }
  }
  const double wall = seconds(clock0) - historyWall, load1 = load();
  samples.push_back(sample());
  record(time());
  // t90: after the last sample below 0.9 of the end value, the linear crossing (none if the end value is not positive).
  const double pEnd = history.back().second, p90 = 0.9 * pEnd;
  double t90 = -1;
  if (pEnd > 0) {
    std::size_t k = history.size();
    for (std::size_t n = 0; n < history.size(); ++n)
      if (!(history[n].second >= p90)) k = n;
    if (k == history.size()) t90 = history.front().first;
    else {
      const auto [ta, pa] = history[k];
      const auto [tb, pb] = history[k + 1];
      t90 = std::isnan(pa) ? tb : ta + (p90 - pa) / (pb - pa) * (tb - ta);
    }
  }
  Sample mean{0, 0, 0};
  for (const auto& s : samples) {
    mean.cstar += s.cstar / static_cast<double>(samples.size());
    mean.isp += s.isp / static_cast<double>(samples.size());
    mean.fVac += s.fVac / static_cast<double>(samples.size());
  }
  const Sample last = samples.back();
  std::uint64_t steps = 0;
  if (uniform) {
    const auto m = single->measurements();
    steps = m.steps;
    updates = m.steps * finest.cells.size();
    std::printf("uniform: %llu steps (%llu rejected); budgets (compensated): mass %.3e, energy %.3e, momentum %.3e; igniter |E - 0.5 J| %.3e J\n",
                static_cast<unsigned long long>(m.steps), static_cast<unsigned long long>(m.rejectedSteps),
                m.compensatedMassBalanceError, m.compensatedEnergyBalanceError, m.compensatedMomentumBalanceError,
                std::abs(m.igniterEnergy - kIgniter));
  } else {
    steps = amr->steps();
    const auto b = amr->budget();
    std::uint64_t n = 0;
    for (auto a : amr->active(1)) n += a;
    std::printf("AMR: %llu base steps (%llu redone), %llu regrids; level 1 active %llu of %zu at the end, %.1f on average; "
                "budgets at the end: mass %.3e of the content, energy %.3e, momentum %.3e of |initial| + gross; igniter |E - 0.5 J| %.3e J\n",
                static_cast<unsigned long long>(amr->steps()), static_cast<unsigned long long>(amr->rejectedSteps()),
                static_cast<unsigned long long>(amr->regrids()), static_cast<unsigned long long>(n), amr->active(1).size(),
                static_cast<double>(activeSum) / static_cast<double>(std::max<std::uint64_t>(amr->steps(), 1)),
                b.massResidual / b.mass, b.energyResidual / std::abs(b.energy), b.momentumResidual / (std::abs(b.momentum) + b.momentumGross),
                std::abs(b.igniterEnergy - kIgniter));
  }
  std::printf("at the end: c* %.4f m/s, Isp_vac %.5f s, F_vac %.4f N; settled over %zu samples: c* %.4f m/s, Isp_vac %.5f s, F_vac %.4f N\n",
              last.cstar, last.isp, last.fVac, samples.size(), mean.cstar, mean.isp, mean.fVac);
  std::printf("reported, not judged: chamber pressure %.6g Pa at the end; t90 %.6f ms (%zu history samples)\n", pEnd,
              t90 * 1e3, history.size());
  std::printf("printed for 5(a), not judged: %llu cell updates, %.1f s of wall net of the history's %.1f s, %.4e s per update; load %.1f at the start, %.1f at the end\n",
              static_cast<unsigned long long>(updates), wall, historyWall, wall / static_cast<double>(updates), load0, load1);
  std::printf("c4 %s scale %.9g cstar %.12e isp %.12e fvac %.12e lightoff %.12e step %.12e steps %llu samples %zu t90 %.12e p8 %.12e\n",
              uniform ? "uniform" : "amr", uniform ? 0.0 : scale, mean.cstar, mean.isp, mean.fVac, lightOff, lightStep,
              static_cast<unsigned long long>(steps), samples.size(), t90, pEnd);
  return 0;
}
