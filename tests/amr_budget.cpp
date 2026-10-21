// AMR criterion 2 on step B1 (docs/evidence/AMR.md, *Step B1: how criterion 2 is run*, recorded before this
// file was written): budgets across coarse-fine faces under subcycling, with the pattern fixed.
// C1 as `crucible_chamber_study eqtt` (Table A local equilibrium, transport, SST with resolved walls: no wall
// functions, no property corrections; one supply over the whole injector face; the igniter 0.5 J over 0.2 ms
// from 0.2 ms in z 5 to 15 mm, r up to 10 mm), base 32 x 6 plus one level (64 x 12), to 8 ms, beside the
// single-level 32 x 6 run on the same mesh. Level 1 covers base rows 0 and 1 of columns 0 to 15 and base
// row 5 in every column. Judged at 8 ms, with the gas's largest mass and largest |energy| over every base step:
//   - the composite mass and energy budget errors below 1e-12 of those, and no more than 10 times the
//     single-level run's, normalised the same way;
//   - the axial momentum balance error below 1e-12 of the gross momentum exchange, likewise;
//   - the composite igniter energy (the levels' owned deliveries summed) within 1e-14 of 0.5 J, relative.
// Printed beside, not judged: the single-level values, c* and vacuum thrust at the end on both runs, the
// largest errors over the run, and the wall time.
// Usage: crucible_amr_budget [end time s, default 8e-3] [threads, default 4]; CRUCIBLE_EQ_TABLE names Table A.
// Step B2's checks (AMR.md, Step B2: setup): b1 runs B1 (with the temperature-guess reset) instead of B2,
// ghost=N sets the ghost width, dump=PATH writes the active cells, steps and budget at the end, and amr
// stops after the AMR run (no single-level run, nothing judged).
// Step C's checks (AMR.md, Step C: setup): tags turns on the tags and regridding at Tags' defaults, with the f_H
// weights from each species' H content. It prints the regrids that changed the pattern and level 1's active cells,
// and judges check 2: the largest fill content error over every regrid below 1e-14, over the parent's content as
// stated (the first after light-off printed) and, beside it, over the children's gross content (AMR.md, Step C: results).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "cantera/core.h"
#include "core/amr.hpp"
#include "core/equilibrium_table.hpp"
#include "tests/amr_dump.hpp"
#include "tests/c1_geometry.hpp"

using namespace crucible;
using namespace crucible::c1_geometry;

namespace {
constexpr double kMixtureRatio = 5, kMassFlow = 0.4, kAmbientP = 1000, kAmbientT = 300, kIgniter = 0.5;
int failures = 0;
void check(bool ok, const std::string& what, double value, double limit) {
  std::printf("  %-74s %.3e (limit %.3g)  %s\n", what.c_str(), value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}
double seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// C1 eqtt as tests/chamber_study.cpp builds it, on nz x nr.
Definition c1eqtt(thermo::ReactionSource& source, int nz, int nr) {
  Definition d;
  d.experiment = Case::Chamber;
  d.nz = nz;
  d.nr = nr;
  d.contour = geometry().contour;
  d.species = source.medium().species();
  d.composition = source.massFractions("N2:1");
  d.backPressure = kAmbientP;
  d.ambientTemperature = kAmbientT;
  Supply s;
  s.outerRadius = kRc;
  s.massFlow = kMassFlow;
  s.totalTemperature = 300;
  auto h2 = source.massFractions("H2:1"), o2 = source.massFractions("O2:1");
  s.composition.resize(h2.size());
  for (std::size_t k = 0; k < h2.size(); ++k) s.composition[k] = (h2[k] + kMixtureRatio * o2[k]) / (1 + kMixtureRatio);
  s.opens = 0;
  s.ramp = 5e-4;
  s.turbulenceIntensity = 0.05;
  s.viscosityRatio = 10;
  d.supplies = {s};
  d.igniter = {0.005, 0.015, 0.010, kIgniter, 2e-4, 2e-4};
  d.transport = thermo::transportFits("h2o2.yaml");
  d.wallTemperature = 600;
  d.turbulence.enabled = true;
  const Medium medium = d.medium();
  std::vector<double> diffusion(d.species.size()), work;
  const double a0 = medium.soundSpeed(kAmbientT, d.composition.data());
  const double rho0 = kAmbientP / (medium.gasConstant(d.composition.data()) * kAmbientT);
  const double mu0 = medium.transport(kAmbientT, kAmbientP, d.composition.data(), diffusion.data(), work).viscosity;
  d.turbulence.ambientK = 9e-9 * a0 * a0;
  d.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
  return d;
}

// Equilibrium isentropic stagnation pressure of a cell, as tests/chamber_study.cpp.
double stagnationPressure(Cantera::ThermoPhase& gas, double t, double p, const double* y, double kinetic) {
  gas.setState_TPY(t, p, y);
  gas.equilibrate("TP");
  const double s = gas.entropy_mass(), h0 = gas.enthalpy_mass() + kinetic;
  double p0 = p;
  for (int n = 0; n < 50; ++n) {
    gas.setState_SP(s, p0);
    gas.equilibrate("SP");
    double step = (h0 - gas.enthalpy_mass()) / (p0 / gas.density());
    p0 *= std::exp(step);
    if (std::abs(step) < 1e-13) break;
  }
  return p0;
}
// Mass-flux-averaged equilibrium stagnation pressure of a column (chamber_study's chamber pressure).
double columnStagnation(const Flow& flow, Cantera::ThermoPhase& gas, int column) {
  const auto& mesh = flow.mesh();
  double p0 = 0, flux = 0;
  for (int j = 0; j < mesh.nr; ++j) {
    auto q = mesh.index(column, j);
    auto w = flow.cellPrimitive(q);
    auto y = flow.massFractions(q);
    double g = w.rho * w.uz * mesh.axialArea(column + 1, j);
    p0 += g * stagnationPressure(gas, flow.temperature(q), w.p, y.data(), 0.5 * (w.uz * w.uz + w.ur * w.ur));
    flux += g;
  }
  return p0 / flux;
}

// The largest content and the errors, normalised as the criterion reads them.
struct Tracker {
  double largestMass = 0, largestEnergy = 0, worstMass = 0, worstEnergy = 0, worstMomentum = 0;
  void content(double mass, double energy) {
    largestMass = std::max(largestMass, mass);
    largestEnergy = std::max(largestEnergy, std::abs(energy));
  }
  void errors(double mass, double energy, double momentum) {
    worstMass = std::max(worstMass, std::abs(mass));
    worstEnergy = std::max(worstEnergy, std::abs(energy));
    worstMomentum = std::max(worstMomentum, std::abs(momentum));
  }
};
}  // namespace

int main(int argc, char** argv) {
  double end = 8e-3;
  int threads = 4, ghost = 0, positional = 0;
  bool b1 = false, amrOnly = false, tagged = false;
  std::string dumpPath;
  for (int k = 1; k < argc; ++k) {
    const std::string a = argv[k];
    if (a == "b1") b1 = true;
    else if (a == "amr") amrOnly = true;
    else if (a == "tags") tagged = true;
    else if (a.rfind("ghost=", 0) == 0) ghost = std::atoi(a.c_str() + 6);
    else if (a.rfind("dump=", 0) == 0) dumpPath = a.substr(5);
    else if (positional++ == 0) end = std::atof(a.c_str());
    else threads = std::max(1, std::atoi(a.c_str()));
  }
  const char* path = std::getenv("CRUCIBLE_EQ_TABLE");
  if (!path) throw std::runtime_error("crucible_amr_budget needs CRUCIBLE_EQ_TABLE (Table A's file)");
  const EquilibriumTable table = EquilibriumTable::load(path);
  thermo::ReactionSource source("h2o2.yaml");
  const Definition d = c1eqtt(source, 32, 6);
  std::printf("AMR criterion 2 on %s: C1 eqtt (resolved walls), base 32 x 6 plus one level, to %.4g ms, %d threads%s%s\n", b1 ? "B1" : "B2",
              end * 1e3, threads, ghost > 0 ? (", ghost width " + std::to_string(ghost)).c_str() : "",
              tagged ? ", step C's tags and regridding" : "");
  std::printf("table %s: %s\n", path, table.provenance.c_str());

  Amr amr(d, 1);
  amr.setRestricted(!b1);
  if (ghost > 0) amr.setGhostWidth(ghost);
  std::vector<std::unique_ptr<thermo::ReactingFlow>> reacting;
  for (int l = 0; l < amr.levels(); ++l) {
    reacting.push_back(std::make_unique<thermo::ReactingFlow>(amr.level(l), "h2o2.yaml", threads, 1e-6, 1e-12,
                                                              thermo::Chemistry::LocalEquilibrium));
    reacting.back()->useEquilibriumTable(table, 0);
  }
  amr.setAfterStep([&](int l, double dt) { reacting[static_cast<std::size_t>(l)]->react(dt); });
  const Mesh& base = amr.level(0).mesh();
  const Mesh& fine = amr.level(1).mesh();
  std::vector<unsigned char> active(fine.cells.size(), 0);
  for (int i = 0; i < fine.nz; ++i)
    for (int j = 0; j < fine.nr; ++j) active[fine.index(i, j)] = ((i >> 1) < 16 && (j >> 1) <= 1) || (j >> 1) == 5;
  amr.setActive(1, active);
  amr.averageDown();
  int throat = 0, endColumn = 0;
  for (int i = 0; i <= base.nz; ++i)
    if (base.radius[i] < base.radius[throat]) throat = i;
  for (int i = 0; i < base.nz; ++i)
    if (base.cells[base.index(i, 0)].z < geometry().cylinderEnd) endColumn = i;
  const double at = std::numbers::pi * base.radius[throat] * base.radius[throat];
  std::size_t owned = 0, igniterFine = 0, igniterBase = 0;
  for (auto a : active) owned += a;
  for (std::size_t q = 0; q < fine.cells.size(); ++q) {
    const auto& c = fine.cells[q];
    igniterFine += c.z >= d.igniter.zMin && c.z <= d.igniter.zMax && c.r <= d.igniter.rMax;
  }
  for (std::size_t q = 0; q < base.cells.size(); ++q) {
    const auto& c = base.cells[q];
    igniterBase += c.z >= d.igniter.zMin && c.z <= d.igniter.zMax && c.r <= d.igniter.rMax;
  }
  std::printf("throat station %d, chamber-end column %d; level 1: %zu of %zu cells active; igniter centroids: %zu fine, %zu base\n",
              throat, endColumn, owned, fine.cells.size(), igniterFine, igniterBase);

  auto historySolution = Cantera::newSolution("h2o2.yaml", "", "none");
  Cantera::ThermoPhase& gas = *historySolution->thermo();
  if (tagged) {
    Amr::Tags tags;
    const std::size_t eH = gas.elementIndex("H");
    for (const auto& sp : d.species) {
      const std::size_t k = gas.speciesIndex(sp.name);
      tags.hydrogen.push_back(gas.nAtoms(k, eH) * gas.atomicWeight(eH) / gas.molecularWeight(k));
    }
    amr.setTags(tags);
  }
  auto activeCells = [&](int l) {
    std::size_t n = 0;
    for (auto a : amr.active(l)) n += a;
    return n;
  };
  std::uint64_t regrids = 0;
  double fillWorst = 0, fillGrossWorst = 0, fillFirstLit = -1;

  // The AMR run.
  // The momentum error's denominator as Flow::measurements: |initial momentum| + the gross exchange.
  Tracker composite;
  const auto start = amr.budget();
  const double momentum0 = std::abs(start.momentum);
  composite.content(start.mass, start.energy);
  auto momentumScale = [&](const Amr::Budget& b) { return std::max(momentum0 + b.momentumGross, std::numeric_limits<double>::min()); };
  const auto clock0 = std::chrono::steady_clock::now();
  double nextPrint = 5e-4;
  while (amr.time() < end) {
    amr.step(end - amr.time());
    if (amr.regrids() != regrids) {
      regrids = amr.regrids();
      fillWorst = std::max(fillWorst, amr.fillError());
      fillGrossWorst = std::max(fillGrossWorst, amr.fillGross());
      if (fillFirstLit < 0 && amr.time() > d.igniter.start + d.igniter.duration) fillFirstLit = amr.fillError();
    }
    const auto b = amr.budget();
    composite.content(b.mass, b.energy);
    composite.errors(b.massResidual / composite.largestMass, b.energyResidual / composite.largestEnergy,
                     b.momentumResidual / momentumScale(b));
    if (amr.time() >= nextPrint || amr.time() >= end) {
      std::printf("  AMR t %.4f ms, %llu base steps (%llu redone), %.0f s: mass %.2e, energy %.2e, momentum %.2e; igniter %.15f J\n",
                  amr.time() * 1e3, static_cast<unsigned long long>(amr.steps()), static_cast<unsigned long long>(amr.rejectedSteps()),
                  seconds(clock0), b.massResidual / composite.largestMass, b.energyResidual / composite.largestEnergy,
                  b.momentumResidual / momentumScale(b), b.igniterEnergy);
      if (tagged)
        std::printf("    regrids %llu, level 1 active %zu of %zu, fill error largest %.3e (gross %.3e)\n", static_cast<unsigned long long>(regrids),
                    activeCells(1), amr.active(1).size(), fillWorst, fillGrossWorst);
      std::fflush(stdout);
      nextPrint += 5e-4;
    }
  }
  const double amrWall = seconds(clock0);
  if (!dumpPath.empty()) writeAmrDump(amr, dumpPath);
  const auto final = amr.budget();
  double mdot = 0, fVac = 0;
  for (int l = 0; l < amr.levels(); ++l) {
    const auto m = amr.level(l).measurements();
    mdot += m.outletMassFlow;
    fVac += m.inletMomentumFlux + m.wallAxialForce + m.bodyAxialForce;
  }
  const double amrCstar = columnStagnation(amr.level(0), gas, endColumn) * at / mdot;
  for (int l = 0; l < amr.levels(); ++l) {
    const auto& s = reacting[static_cast<std::size_t>(l)]->stats();
    std::printf("  level %d: table cells %ld, clamped %ld, outside (direct) %ld\n", l, s.tableCells, s.tableClamped, s.tableFallbacks);
  }
  if (amrOnly) {
    std::printf("AMR only, not judged: mass %.3e, energy %.3e, momentum %.3e at the end (largest %.3e, %.3e, %.3e); igniter |E - 0.5 J| %.3e J; "
                "c* %.2f m/s, vacuum thrust %.2f N; %llu base steps, %.1f s of wall\n",
                final.massResidual / composite.largestMass, final.energyResidual / composite.largestEnergy,
                final.momentumResidual / momentumScale(final), composite.worstMass, composite.worstEnergy, composite.worstMomentum,
                std::abs(final.igniterEnergy - kIgniter), amrCstar, fVac, static_cast<unsigned long long>(amr.steps()), amrWall);
    if (tagged)
      std::printf("AMR only, not judged: %llu regrids, fill error largest %.3e (first after light-off %.3e), over the gross content %.3e\n",
                  static_cast<unsigned long long>(regrids), fillWorst, fillFirstLit, fillGrossWorst);
    return 0;
  }

  // The single-level run on the same base mesh.
  Flow single(d, Mesh(base));
  thermo::ReactingFlow singleReacting(single, "h2o2.yaml", threads, 1e-6, 1e-12, thermo::Chemistry::LocalEquilibrium);
  singleReacting.useEquilibriumTable(table, 0);
  Tracker alone;
  const auto m0 = single.measurements();
  alone.content(m0.mass, m0.energy);
  const auto clock1 = std::chrono::steady_clock::now();
  while (single.time() < end) {
    singleReacting.step(end - single.time());
    const auto m = single.measurements();
    alone.content(m.mass, m.energy);
    alone.errors(m.massBalanceError * m0.mass / alone.largestMass, m.energyBalanceError * m0.energy / alone.largestEnergy,
                 m.momentumBalanceError);
  }
  const double singleWall = seconds(clock1);
  const auto ms = single.measurements();
  std::printf("contents: AMR mass %.9e to %.9e kg, energy %.9e to %.9e J; single level mass %.9e to %.9e kg, energy %.9e to %.9e J\n",
              start.mass, final.mass, start.energy, final.energy, m0.mass, ms.mass, m0.energy, ms.energy);
  const double singleMass = ms.massBalanceError * m0.mass / alone.largestMass;
  const double singleEnergy = ms.energyBalanceError * m0.energy / alone.largestEnergy;
  const double singleCstar = columnStagnation(single, gas, endColumn) * at / ms.outletMassFlow;
  const double singleFVac = ms.inletMomentumFlux + ms.wallAxialForce + ms.bodyAxialForce;

  const double amrMass = final.massResidual / composite.largestMass, amrEnergy = final.energyResidual / composite.largestEnergy;
  const double amrMomentum = final.momentumResidual / momentumScale(final);
  std::printf("at %.4f ms (largest mass %.6e kg and |energy| %.6e J over the AMR run's base steps; %.6e and %.6e alone)\n", end * 1e3,
              composite.largestMass, composite.largestEnergy, alone.largestMass, alone.largestEnergy);
  check(std::abs(amrMass) < 1e-12, "composite mass budget error / largest mass", std::abs(amrMass), 1e-12);
  check(std::abs(amrMass) <= 10 * std::abs(singleMass), "  the same against 10 x the single-level run's",
        std::abs(amrMass), 10 * std::abs(singleMass));
  check(std::abs(amrEnergy) < 1e-12, "composite energy budget error / largest |energy|", std::abs(amrEnergy), 1e-12);
  check(std::abs(amrEnergy) <= 10 * std::abs(singleEnergy), "  the same against 10 x the single-level run's",
        std::abs(amrEnergy), 10 * std::abs(singleEnergy));
  check(std::abs(amrMomentum) < 1e-12, "composite axial momentum balance error / gross exchange", std::abs(amrMomentum), 1e-12);
  check(std::abs(amrMomentum) <= 10 * std::abs(ms.momentumBalanceError), "  the same against 10 x the single-level run's",
        std::abs(amrMomentum), 10 * std::abs(ms.momentumBalanceError));
  // Criterion 2 says "exact to 1e-14 relative"; the first run judged 1e-14 J, twice as loose (AMR.md, criterion 2 results).
  const double igniterError = std::abs(final.igniterEnergy - kIgniter) / kIgniter;
  check(igniterError <= 1e-14, "composite igniter energy, |E / 0.5 J - 1|", igniterError, 1e-14);
  if (tagged) {
    check(fillWorst <= 1e-14, "step C check 2: the largest fill content error over the regrids", fillWorst, 1e-14);
    check(fillGrossWorst <= 1e-14, "  the same over the children's gross content", fillGrossWorst, 1e-14);
    std::printf("      printed, not judged: %llu regrids changed the pattern; the first fill after light-off %.3e; level 1 ends with %zu of %zu active\n",
                static_cast<unsigned long long>(regrids), fillFirstLit, activeCells(1), amr.active(1).size());
  }
  std::printf("      printed, not judged: AMR igniter |E - 0.5 J| %.3e J; single level mass %.3e, energy %.3e, momentum %.3e, igniter |E - 0.5 J| %.3e J\n",
              std::abs(final.igniterEnergy - kIgniter), singleMass, singleEnergy, ms.momentumBalanceError,
              std::abs(ms.igniterEnergy - kIgniter));
  std::printf("      printed, not judged: largest over the run, AMR mass %.3e, energy %.3e, momentum %.3e; single level %.3e, %.3e, %.3e\n",
              composite.worstMass, composite.worstEnergy, composite.worstMomentum, alone.worstMass, alone.worstEnergy, alone.worstMomentum);
  std::printf("      printed, not judged: AMR c* %.2f m/s, vacuum thrust %.2f N, mdot out %.6f kg/s, %llu base steps, %.1f s of wall\n",
              amrCstar, fVac, mdot, static_cast<unsigned long long>(amr.steps()), amrWall);
  std::printf("      printed, not judged: single level c* %.2f m/s, vacuum thrust %.2f N, mdot out %.6f kg/s, %llu steps, %.1f s of wall\n",
              singleCstar, singleFVac, ms.outletMassFlow, static_cast<unsigned long long>(ms.steps), singleWall);
  std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
