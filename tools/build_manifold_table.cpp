// Builds Table B, the tabulated finite-rate progress manifold (core/manifold_table.hpp;
// TECHNICAL_PLAN, *Lightweight engine*; docs/evidence/TABLE_B.md), offline from the engine's own
// CVODES reactor (adapters/reaction.hpp's ReactionStep::follow) on the full mechanism. Every node
// of (f_H, Z_N, eta, ln rho) is built as two branches meeting at equilibrium: a lower branch from
// the node's reactants (or, for gas that has already lost energy, from a frozen high-T
// equilibrium) up to c_eq, and an upper branch from the complete-combustion products down to c_eq.
// Both branches are stored at the same 40 fractions s of their own span (TABLE_B.md, "The c nodes,
// relative to the branch's own start" and "The runtime": s = (c-c_0)/(c_eq-c_0) on the lower branch,
// s = 1 + (c-c_eq)/(c_cp-c_eq) on the upper), so every node uses the identical s layout and the
// runtime can blend corners point by point. A frozen branch (induction not complete by 1e15 s)
// places the same fractions on the span it actually reached.
//
// Mechanism: h2o2.yaml or a file with its species (default h2o2.yaml; tools/mechanisms/h2o2_highT.yaml
// carries its thermo to 6,000 K), printed in the provenance; the builder refuses a mechanism in which
// N2 reacts (checked directly on its reactions, not inferred).
//
// The progress variable is chosen from three linear-in-Y candidates in order (Y_H2O; Y_H2O + Y_HO2
// + Y_H2O2; the O2 consumed), each tried as a full build: gate 2(a) (turn-back) and the N gate must
// both pass on every node for the candidate to be accepted. This cannot be done cheaper by caching
// one candidate-independent integration, because the stopping time of a branch's own CVODES
// integration (|c - c_eq| <= 1e-6 of its span) is itself defined by the candidate's c; a second or
// third attempt therefore re-integrates. In practice the first candidate is expected to pass (the
// probe found Y_H2O never turns back), so this costs nothing on the runs that matter.
//
// Usage: crucible_build_manifold_table <out> <f_H nodes, lean side> <f_H nodes, rich side>
//                                       <f_H spacing at f_st> <Z_N list, comma separated>
//                                       <eta nodes> <ln rho nodes> <rho min> <rho max> <threads> [mechanism]
// f_H is built exactly as Table A's f axis (tools/build_equilibrium_table.cpp): uniform on [0, f_st]
// and [f_st, 1] with f_st the H mass fraction in H2O, plus clustering at f_st and the ends. eta is
// uniform on [0, 1]. ln rho is uniform over [log(rho min), log(rho max)].
//
// The table is written only if every gate passes. Set CRUCIBLE_TABLE_B_WRITE_ANYWAY=1 to write it
// regardless, clearly labelled in the provenance, for inspecting a failed build.
// CRUCIBLE_TABLE_B_NODE_REPORT=<path> writes one CSV row per node and candidate: its gate numbers,
// flags, c_0, c_eq, c_cp, T_eq and, for an integrator failure, the exception's text.
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "adapters/reaction.hpp"
#include "cantera/base/global.h"
#include "cantera/core.h"
#include "core/manifold_table.hpp"

using namespace crucible;

namespace {

constexpr double kTEnd = 1e15;            // a branch's integration horizon (TABLE_B.md)
constexpr double kStopFrac = 1e-6;        // a branch stops within this fraction of its span of c_eq
constexpr double kTurnBackWindow = 1e-4;  // turn-back is measured only until c is this close to c_eq
constexpr double kTurnBackGate = 1e-9;    // above this, a turn-back fails gate 2(a) (CVODES noise floor)
constexpr double kZeroSpanTol = 1e-12;    // |c_eq - c_start| below this: the branch is a point
constexpr double kPhiA = 0.005;           // c_a: the lower branch's 0.5% fall in N
constexpr double kDTGate2b = 3.0;         // K, gate 2(b)
constexpr std::size_t kNotReacting = std::numeric_limits<std::size_t>::max();

struct Indices { std::size_t h2, o2, h2o, n2; double hInWater; };

// Raw reactants (H2, O2, N2) at the given elements: no water, since nothing has burned.
std::vector<double> reactantsOf(const Indices& s, std::size_t ns, double zh, double zo, double zn) {
  std::vector<double> y(ns, 0.0);
  y[s.h2] = zh;
  y[s.o2] = zo;
  y[s.n2] = zn;
  return y;
}

// Complete-combustion products of the elements (as tools/build_equilibrium_table.cpp): all of the
// scarcer element in H2O, the rest as the excess H2 or O2.
std::vector<double> productsOf(const Indices& s, std::size_t ns, double zh, double zo, double zn) {
  std::vector<double> y(ns, 0.0);
  if (zh * (1 - s.hInWater) >= zo * s.hInWater) {
    y[s.h2o] = zo / (1 - s.hInWater);
    y[s.h2] = std::max(zh - y[s.h2o] * s.hInWater, 0.0);
  } else {
    y[s.h2o] = zh / s.hInWater;
    y[s.o2] = std::max(zo - y[s.h2o] * (1 - s.hInWater), 0.0);
  }
  y[s.n2] = zn;
  return y;
}

double bisect(const std::function<double(double)>& g, double lo, double hi) {
  double glo = g(lo);
  for (int i = 0; i < 200 && hi - lo > 1e-10 * hi; ++i) {
    double mid = 0.5 * (lo + hi), gm = g(mid);
    if ((gm < 0) == (glo < 0)) { lo = mid; glo = gm; } else { hi = mid; }
  }
  return 0.5 * (lo + hi);
}

// A node's two branches' starts, as the probe and TABLE_B.md state them.
struct Starts {
  std::vector<double> lower;  // [T, Y...]
  std::vector<double> upper;  // [T, Y...]
  bool fromReactants = false;
};

// Species indices and the element columns (Table A's style): species made of H and O alone are
// "reacting" (h+o ~ 1); every other species must be the sole carrier of its own elements.
struct SpeciesTable {
  Indices idx;
  std::size_t ns = 0;
  std::vector<std::array<double, 2>> elements;  // per species: (H, O) mass fraction, (0,0) if inert
  std::vector<std::size_t> reactingPosition;    // per species: position in reactingNames, or kNotReacting
  std::vector<std::string> reactingNames;
  std::vector<double> weight;  // molecular weights, per species
  std::size_t iH2O = 0, iHO2 = 0, iH2O2 = 0, iO2 = 0;
};

SpeciesTable classify(thermo::ReactionSource& probe, const std::string& mechanism) {
  SpeciesTable sp;
  sp.ns = probe.nSpecies();
  const std::size_t ne = probe.nElements();
  std::vector<std::vector<double>> columns;
  std::vector<double> unit(sp.ns, 0.0);
  for (std::size_t k = 0; k < sp.ns; ++k) {
    unit[k] = 1;
    columns.push_back(probe.elementMassFractions(unit.data()));
    unit[k] = 0;
  }
  auto index = [&](const char* name) {
    for (std::size_t k = 0; k < sp.ns; ++k) if (probe.speciesName(k) == name) return k;
    throw std::runtime_error(std::string("species missing: ") + name);
  };
  sp.idx.h2 = index("H2");
  sp.idx.o2 = index("O2");
  sp.idx.h2o = index("H2O");
  sp.idx.n2 = index("N2");
  sp.iH2O = sp.idx.h2o;
  sp.iHO2 = index("HO2");
  sp.iH2O2 = index("H2O2");
  sp.iO2 = sp.idx.o2;
  const auto mH = std::max_element(columns[sp.idx.h2].begin(), columns[sp.idx.h2].end()) - columns[sp.idx.h2].begin();
  const auto mO = std::max_element(columns[sp.idx.o2].begin(), columns[sp.idx.o2].end()) - columns[sp.idx.o2].begin();
  sp.elements.resize(sp.ns);
  sp.reactingPosition.assign(sp.ns, kNotReacting);
  for (std::size_t k = 0; k < sp.ns; ++k) {
    const double h = columns[k][mH], o = columns[k][mO];
    if (h + o > 1 - 1e-12) {
      sp.elements[k] = {h, o};
      sp.reactingPosition[k] = sp.reactingNames.size();
      sp.reactingNames.push_back(probe.speciesName(k));
      continue;
    }
    if (h + o > 1e-12) throw std::runtime_error("species " + probe.speciesName(k) + " mixes H/O with another element");
    sp.elements[k] = {0.0, 0.0};
    for (std::size_t m = 0; m < ne; ++m)
      for (std::size_t r = 0; r < sp.ns; ++r)
        if (r != k && columns[k][m] > 0 && columns[r][m] > 0)
          throw std::runtime_error("species " + probe.speciesName(k) + " shares an element with " + probe.speciesName(r));
  }
  sp.idx.hInWater = sp.elements[sp.idx.h2o][0];
  sp.weight.resize(sp.ns);
  auto solution = Cantera::newSolution(mechanism, "", "none");
  solution->thermo()->getMolecularWeights(sp.weight.data());
  return sp;
}

// Refuses a mechanism in which N2 appears as a net reactant or product of any reaction (a
// third-body efficiency does not count). Table A's own rule, applied directly to the kinetics.
void refuseReactingN2(const std::string& mechanism) {
  auto solution = Cantera::newSolution(mechanism, "", "none");
  auto& kin = *solution->kinetics();
  for (std::size_t i = 0; i < kin.nReactions(); ++i) {
    auto rxn = kin.reaction(i);
    double net = 0.0;
    auto r = rxn->reactants.find("N2");
    auto p = rxn->products.find("N2");
    if (r != rxn->reactants.end()) net -= r->second;
    if (p != rxn->products.end()) net += p->second;
    if (net != 0.0) throw std::runtime_error("mechanism has reacting N2: " + rxn->equation());
  }
}

// One raw CVODES step of a branch, kept only transiently while the node is built.
struct RawStep {
  double t = 0, T = 0, c = 0, lambda = 0;
  std::vector<double> y;  // full species, size ns
};

// The s fractions, identical for every node and every build (TABLE_B.md, "The runtime": "every node
// stores its points at the same fractions"). Lower: s=0, 12 geometric 1e-9..0.1, 18 uniform 0.1..1.
std::vector<double> lowerFractions(int geometric, int uniform) {
  std::vector<double> s;
  s.push_back(0.0);
  for (int g = 1; g <= geometric; ++g) {
    const double frac = geometric > 1 ? static_cast<double>(g - 1) / (geometric - 1) : 1.0;
    s.push_back(1e-9 * std::pow(0.1 / 1e-9, frac));
  }
  for (int u = 1; u <= uniform; ++u) s.push_back(0.1 + 0.9 * static_cast<double>(u) / uniform);
  return s;
}

struct NodeOutcome {
  double turnBack = 0;
  bool nGateChecked = false, nGateOk = true;
  bool hasDT2bLower = false, hasDT2bUpper = false;
  double dT2bLower = 0, dT2bUpper = 0;
  bool lowerFrozen = false, upperFrozen = false, integratorFailure = false;
  bool degenerate = false, zeroSpan = false, outside = false;
  double rawElementError = 0, projectedElementError = 0, sumError = 0;
  int branchesIntegrated = 0;
  double branchTimeSumMs = 0, branchTimeWorstMs = 0;
  double c0 = 0, cEq = 0, cCp = 0, tEq = 0;
  std::string failure;  // the exception's text when integratorFailure
};

// Builds one node into `out` at `nodeIdx`. Everything that can throw (equilibration, integration)
// is caught here: a failure leaves the node's stored data unfinished, which is fine because the
// table is only written once gate 2(c) (zero integrator failures) has passed.
NodeOutcome buildNode(thermo::ReactionSource& source, thermo::ReactionStep& step, const SpeciesTable& sp,
                       int candidate, double maxThermoT, const std::vector<double>& sLower,
                       double f, double z, double eta, double rho, double eLo, double eHi,
                       ManifoldTable& out, std::size_t nodeIdx) {
  NodeOutcome result;
  const std::size_t ns = sp.ns;
  ManifoldTable::Node& node = out.nodes[nodeIdx];
  node.progressVariable = candidate;
  const double e = eLo + eta * (eHi - eLo);
  const double zh = f * (1 - z), zo = (1 - f) * (1 - z), zn = z;

  auto energyAt = [&](double t, const std::vector<double>& y) { return source.internalEnergy(t, rho, y.data()); };
  auto temperatureAt = [&](double ee, const std::vector<double>& y) {
    return bisect([&](double t) { return energyAt(t, y) - ee; }, 100.0, 6000.0);
  };
  auto progress = [&](const double* y) -> double {
    switch (candidate) {
      case 0: return y[sp.iH2O];
      case 1: return y[sp.iH2O] + y[sp.iHO2] + y[sp.iH2O2];
      default: return zo - y[sp.iO2];
    }
  };
  auto dcdt = [&](const double* dzdt) -> double {  // dzdt size ns+1, [0]=dT/dt
    switch (candidate) {
      case 0: return dzdt[1 + sp.iH2O];
      case 1: return dzdt[1 + sp.iH2O] + dzdt[1 + sp.iHO2] + dzdt[1 + sp.iH2O2];
      default: return -dzdt[1 + sp.iO2];
    }
  };

  const int lowerCount = out.lowerPoints();
  const int upperCount = out.upperUniform;

  try {
    if (f == 0.0 || f == 1.0) {
      // Nodes c cannot describe (TABLE_B.md): store the equilibrium everywhere, c0=cEq=cCp=0.
      auto y = reactantsOf(sp.idx, ns, zh, zo, zn);
      std::vector<double> zstate(ns + 1);
      zstate[0] = temperatureAt(e, y);
      std::copy(y.begin(), y.end(), zstate.begin() + 1);
      source.equilibrateUV(rho, zstate.data());
      if (zstate[0] > maxThermoT) result.outside = true;
      node.c0 = node.cEq = node.cCp = 0.0;
      node.flags = ManifoldTable::kDegenerate | (result.outside ? ManifoldTable::kOutside : 0u);
      result.degenerate = true;
      for (int p = 0; p < out.pointsPerNode(); ++p) {
        double* py = out.pointY(nodeIdx, p);
        for (std::size_t q = 0; q < ns; ++q)
          if (sp.reactingPosition[q] != kNotReacting) py[sp.reactingPosition[q]] = zstate[q + 1];
        out.lambda(nodeIdx, p) = 0.0;
      }
      return result;
    }

    const std::vector<double> yReactants = reactantsOf(sp.idx, ns, zh, zo, zn);
    const std::vector<double> yComplete = productsOf(sp.idx, ns, zh, zo, zn);

    // The lower branch's start.
    std::vector<double> zLower(ns + 1);
    if (energyAt(200.0, yReactants) <= e) {
      zLower[0] = temperatureAt(e, yReactants);
      std::copy(yReactants.begin(), yReactants.end(), zLower.begin() + 1);
    } else {
      auto frozenAt = [&](double tStar) {
        std::vector<double> w(ns + 1);
        w[0] = tStar;
        std::copy(yReactants.begin(), yReactants.end(), w.begin() + 1);
        source.equilibrateTV(rho, w.data());
        return std::vector<double>(w.begin() + 1, w.end());
      };
      const double tStar = bisect([&](double ts) { return energyAt(200.0, frozenAt(ts)) - e; }, 300.0, 6000.0);
      auto y = frozenAt(tStar * (1 - 1e-9));
      zLower[0] = temperatureAt(e, y);
      std::copy(y.begin(), y.end(), zLower.begin() + 1);
    }

    // The shared equilibrium (both branches meet here).
    std::vector<double> zEq = zLower;
    source.equilibrateUV(rho, zEq.data());
    const double cEq = progress(zEq.data() + 1);
    result.tEq = zEq[0];
    if (zEq[0] > maxThermoT) result.outside = true;

    // The upper branch's start: complete combustion.
    std::vector<double> zUpper(ns + 1);
    zUpper[0] = temperatureAt(e, yComplete);
    std::copy(yComplete.begin(), yComplete.end(), zUpper.begin() + 1);
    const double cCp = progress(zUpper.data() + 1);
    const double c0 = progress(zLower.data() + 1);
    node.c0 = c0;
    node.cEq = cEq;
    node.cCp = cCp;
    result.c0 = c0;
    result.cEq = cEq;
    result.cCp = cCp;
    if (zLower[0] > maxThermoT || zUpper[0] > maxThermoT) result.outside = true;

    const double spanLower = std::fabs(cEq - c0), spanUpper = std::fabs(cCp - cEq);
    const bool zeroSpanLower = spanLower < kZeroSpanTol, zeroSpanUpper = spanUpper < kZeroSpanTol;
    if (zeroSpanLower || zeroSpanUpper) result.zeroSpan = true;

    // Integrates one branch from zStart toward branchCEq, recording every internal step (with its
    // exact rate, hence lambda) until within kStopFrac of the span or the time horizon.
    auto runBranch = [&](const std::vector<double>& zStart, double branchCEq, bool& frozenOut) {
      std::vector<RawStep> steps;
      const double cStart = progress(zStart.data() + 1);
      const double span = std::fabs(branchCEq - cStart);
      std::vector<double> dzdt(ns + 1);
      auto withLambda = [&](const double* z, double t) {
        source.rates(rho, z, dzdt.data());
        const double c = progress(z + 1);
        const double denom = branchCEq - c;
        RawStep s;
        s.t = t;
        s.T = z[0];
        s.c = c;
        s.lambda = denom != 0.0 ? dcdt(dzdt.data()) / denom : 0.0;
        s.y.assign(z + 1, z + 1 + ns);
        return s;
      };
      steps.push_back(withLambda(zStart.data(), 0.0));
      frozenOut = true;
      std::vector<double> z = zStart;
      const auto clock0 = std::chrono::steady_clock::now();
      step.follow(rho, z.data(), kTEnd, [&](double t, const double* zz) {
        steps.push_back(withLambda(zz, t));
        const bool close = std::fabs(branchCEq - steps.back().c) <= kStopFrac * span;
        if (close) frozenOut = false;
        return !close;
      });
      const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - clock0).count();
      result.branchTimeSumMs += ms;
      result.branchTimeWorstMs = std::max(result.branchTimeWorstMs, ms);
      ++result.branchesIntegrated;
      return steps;
    };

    std::vector<RawStep> lowerSteps, upperSteps;
    bool lowerFrozen = false, upperFrozen = false;
    if (!zeroSpanLower) lowerSteps = runBranch(zLower, cEq, lowerFrozen);
    if (!zeroSpanUpper) upperSteps = runBranch(zUpper, cEq, upperFrozen);
    result.lowerFrozen = lowerFrozen;
    result.upperFrozen = upperFrozen;
    if (lowerFrozen || upperFrozen) node.flags |= ManifoldTable::kFrozen;

    for (auto& s : lowerSteps) if (s.T > maxThermoT) result.outside = true;
    for (auto& s : upperSteps) if (s.T > maxThermoT) result.outside = true;
    node.flags |= (result.outside ? ManifoldTable::kOutside : 0u) | (result.zeroSpan ? ManifoldTable::kZeroSpan : 0u) |
                  (result.degenerate ? ManifoldTable::kDegenerate : 0u);

    // Turn-back (gate 2a): the largest fall against the approach, relative to the span, measured
    // only until c is within kTurnBackWindow of the span from c_eq.
    auto turnBackOf = [&](const std::vector<RawStep>& steps, double branchCEq, double span) {
      if (steps.empty() || span < kZeroSpanTol) return 0.0;
      const double sign = branchCEq >= steps.front().c ? 1.0 : -1.0;
      double best = steps.front().c, worst = 0.0;
      for (const auto& s : steps) {
        if (std::fabs(branchCEq - s.c) <= kTurnBackWindow * span) break;
        if (sign * (s.c - best) > 0) best = s.c;
        worst = std::max(worst, sign * (best - s.c) / span);
      }
      return worst;
    };
    result.turnBack = std::max(turnBackOf(lowerSteps, cEq, spanLower), turnBackOf(upperSteps, cEq, spanUpper));

    // Gate 2(b): the branch's own reached end state against the independent equilibrateUV call.
    if (!zeroSpanLower && !lowerFrozen && !lowerSteps.empty()) {
      result.hasDT2bLower = true;
      result.dT2bLower = std::fabs(lowerSteps.back().T - zEq[0]);
    }
    if (!zeroSpanUpper && !upperFrozen && !upperSteps.empty()) {
      result.hasDT2bUpper = true;
      result.dT2bUpper = std::fabs(upperSteps.back().T - zEq[0]);
    }

    // Interpolates a branch's recorded steps linearly in c to the target c (TABLE_B.md, "Values at
    // a node"): the two bracketing internal steps, c monotone there by gate 2(a).
    auto atC = [&](const std::vector<RawStep>& steps, double c, std::vector<double>& y, double& lambdaOut) {
      for (std::size_t i = 1; i < steps.size(); ++i) {
        const double a = steps[i - 1].c, b = steps[i].c;
        if ((c - a) * (c - b) <= 0 && a != b) {
          const double w = (c - a) / (b - a);
          y.resize(ns);
          for (std::size_t q = 0; q < ns; ++q) y[q] = steps[i - 1].y[q] + w * (steps[i].y[q] - steps[i - 1].y[q]);
          lambdaOut = steps[i - 1].lambda + w * (steps[i].lambda - steps[i - 1].lambda);
          return;
        }
      }
      const auto& nearest = (!steps.empty() && std::fabs(c - steps.front().c) <= std::fabs(c - steps.back().c))
                                 ? steps.front()
                                 : steps.back();
      y = nearest.y;
      lambdaOut = nearest.lambda;
    };

    // Table A's projection: the smallest change weighted by Y_k that restores the node's exact
    // (H, O) elements.
    auto project = [&](std::vector<double>& y, double& rawErr, double& projErr, double& sumErr) {
      // The table stores the reacting species only; the runtime rebuilds the inert ones from Z_N.
      // So they are checked at their exact values (N2 = Z_N, others 0), not at the integration's,
      // which drift by about rtol.
      for (std::size_t q = 0; q < ns; ++q)
        if (sp.reactingPosition[q] == kNotReacting) y[q] = q == sp.idx.n2 ? zn : 0.0;
      auto errors = [&](double err[2]) {
        double total = 0;
        err[0] = -zh;
        err[1] = -zo;
        for (std::size_t q = 0; q < ns; ++q) {
          total += y[q];
          err[0] += y[q] * sp.elements[q][0];
          err[1] += y[q] * sp.elements[q][1];
        }
        return total - 1.0;
      };
      double err[2];
      errors(err);
      rawErr = std::max(std::fabs(err[0]), std::fabs(err[1]));
      for (int pass = 0; pass < 2; ++pass) {
        double m[2][2] = {}, lam[2] = {-err[0], -err[1]};
        for (std::size_t q = 0; q < ns; ++q)
          for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b) m[a][b] += y[q] * sp.elements[q][a] * sp.elements[q][b];
        for (int a = 0; a < 2; ++a)
          if (m[a][a] == 0) { m[a][a] = 1; m[a][1 - a] = m[1 - a][a] = 0; lam[a] = 0; }
        const double det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        const double l0 = (lam[0] * m[1][1] - m[0][1] * lam[1]) / det;
        const double l1 = (m[0][0] * lam[1] - m[1][0] * lam[0]) / det;
        for (std::size_t q = 0; q < ns; ++q) y[q] += y[q] * (sp.elements[q][0] * l0 + sp.elements[q][1] * l1);
        errors(err);
      }
      sumErr = std::fabs(errors(err));
      projErr = std::max(std::fabs(err[0]), std::fabs(err[1]));
    };

    auto storePoint = [&](int p, std::vector<double> y, double lambdaVal) {
      double rawErr, projErr, sumErr;
      project(y, rawErr, projErr, sumErr);
      result.rawElementError = std::max(result.rawElementError, rawErr);
      result.projectedElementError = std::max(result.projectedElementError, projErr);
      result.sumError = std::max(result.sumError, sumErr);
      double* py = out.pointY(nodeIdx, p);
      for (std::size_t q = 0; q < ns; ++q)
        if (sp.reactingPosition[q] != kNotReacting) py[sp.reactingPosition[q]] = y[q];
      out.lambda(nodeIdx, p) = lambdaVal;
    };
    std::vector<double> eqY(zEq.begin() + 1, zEq.end());

    const double cReachLower = zeroSpanLower ? cEq : (lowerFrozen ? lowerSteps.back().c : cEq);
    for (int p = 0; p < lowerCount; ++p) {
      if (zeroSpanLower) { storePoint(p, eqY, 0.0); continue; }
      const double c = c0 + sLower[p] * (cReachLower - c0);
      if (p == lowerCount - 1 && !lowerFrozen) {
        storePoint(p, eqY, 0.0);  // lambda fixed below (copied from the point before)
      } else {
        std::vector<double> y;
        double lam;
        atC(lowerSteps, c, y, lam);
        storePoint(p, y, lam);
      }
    }
    if (!zeroSpanLower && !lowerFrozen) out.lambda(nodeIdx, lowerCount - 1) = out.lambda(nodeIdx, lowerCount - 2);

    const double cReachUpper = zeroSpanUpper ? cEq : (upperFrozen ? upperSteps.back().c : cEq);
    for (int p = 0; p < upperCount; ++p) {
      const int idx = lowerCount + p;
      if (zeroSpanUpper) { storePoint(idx, eqY, 0.0); continue; }
      const double c = cReachUpper + (upperCount > 1 ? static_cast<double>(p) / (upperCount - 1) : 1.0) * (cCp - cReachUpper);
      if (p == 0 && !upperFrozen) {
        storePoint(idx, eqY, 0.0);  // lambda fixed below
      } else {
        std::vector<double> y;
        double lam;
        atC(upperSteps, c, y, lam);
        storePoint(idx, y, lam);
      }
    }
    if (!zeroSpanUpper && !upperFrozen) out.lambda(nodeIdx, lowerCount) = out.lambda(nodeIdx, lowerCount - 1);

    // The N gate (03:10, 2a): N = sum Y_k/W_k strictly falls from c_a to c_cp along the stored
    // points. c_a is the first stored point past phi_N(c) = 0.5% of the lower branch's own fall.
    if (!result.degenerate) {
      auto nAt = [&](int p) {
        const double* py = out.pointY(nodeIdx, p);
        double m = 0;
        for (std::size_t q = 0; q < ns; ++q)
          if (sp.reactingPosition[q] != kNotReacting) m += py[sp.reactingPosition[q]] / sp.weight[q];
        return m;
      };
      const double n0 = nAt(0), neq = nAt(lowerCount - 1);
      int idxA = lowerCount;  // no portion to check if the fall never reaches 0.5%
      if (std::fabs(n0 - neq) > 1e-15) {
        for (int p = 0; p < lowerCount; ++p) {
          const double phi = (n0 - nAt(p)) / (n0 - neq);
          if (phi >= kPhiA) { idxA = p; break; }
        }
      }
      result.nGateChecked = true;
      bool ok = true;
      for (int p = idxA; p + 1 < out.pointsPerNode(); ++p) {
        const double na = nAt(p), nb = nAt(p + 1);
        if (!(nb < na + 1e-13 * std::max(1.0, std::fabs(na)))) { ok = false; break; }
      }
      result.nGateOk = ok;
    }
  } catch (const std::exception& ex) {
    result.integratorFailure = true;
    result.failure = ex.what();
  }
  return result;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 11 && argc != 12) {
    std::fprintf(stderr,
                 "usage: %s <out> <f_H lean> <f_H rich> <f_H spacing at f_st> <Z_N list, comma> <eta nodes> "
                 "<ln rho nodes> <rho min> <rho max> <threads> [mechanism]\n",
                 argv[0]);
    return 2;
  }
  const std::string out = argv[1];
  const int fLean = std::atoi(argv[2]), fRich = std::atoi(argv[3]);
  const double fCluster = std::atof(argv[4]);
  std::vector<double> zList;
  {
    std::string s = argv[5];
    std::size_t pos = 0;
    while (pos <= s.size()) {
      const std::size_t comma = s.find(',', pos);
      zList.push_back(std::atof(s.substr(pos, comma - pos).c_str()));
      if (comma == std::string::npos) break;
      pos = comma + 1;
    }
    std::sort(zList.begin(), zList.end());
  }
  const int nEta = std::atoi(argv[6]), nLnRho = std::atoi(argv[7]);
  const double rhoMin = std::atof(argv[8]), rhoMax = std::atof(argv[9]);
  const int threads = std::max(1, std::atoi(argv[10]));
  const std::string mechanism = argc == 12 ? argv[11] : "h2o2.yaml";
  const bool writeAnyway = std::getenv("CRUCIBLE_TABLE_B_WRITE_ANYWAY") != nullptr;
  const char* nodeReportPath = std::getenv("CRUCIBLE_TABLE_B_NODE_REPORT");
  std::FILE* nodeReport = nodeReportPath ? std::fopen(nodeReportPath, "w") : nullptr;
  if (nodeReportPath && !nodeReport) throw std::runtime_error(std::string("cannot write ") + nodeReportPath);
  if (nodeReport)
    std::fprintf(nodeReport, "candidate,f_H,Z_N,eta,rho,T_eq,c0,c_eq,c_cp,degenerate,zero_span,outside,lower_frozen,upper_frozen,"
                             "integrator_failure,turn_back,n_gate_checked,n_gate_ok,dT_2b_lower,dT_2b_upper,raw_element_error,"
                             "sum_error,failure\n");

  refuseReactingN2(mechanism);

  thermo::ReactionSource probe(mechanism);
  const SpeciesTable sp = classify(probe, mechanism);
  double minT, maxT;
  {
    auto solution = Cantera::newSolution(mechanism, "", "none");
    minT = solution->thermo()->minTemp();
    maxT = solution->thermo()->maxTemp();
  }
  (void)minT;

  ManifoldTable axes;
  for (int i = 0; i <= fLean; ++i) axes.fH.push_back(sp.idx.hInWater * i / fLean);
  for (int i = 1; i <= fRich; ++i) axes.fH.push_back(sp.idx.hInWater + (1 - sp.idx.hInWater) * i / fRich);
  std::size_t clustered = 0;
  if (fCluster > 0) {
    for (double d = fCluster; d < sp.idx.hInWater / fLean; d *= 1.5, ++clustered) axes.fH.push_back(sp.idx.hInWater - d);
    for (double d = fCluster; d < (1 - sp.idx.hInWater) / fRich; d *= 1.5, ++clustered) axes.fH.push_back(sp.idx.hInWater + d);
    for (double d = fCluster; d < sp.idx.hInWater / fLean; d *= 1.5, ++clustered) axes.fH.push_back(d);
    for (double d = fCluster; d < (1 - sp.idx.hInWater) / fRich; d *= 1.5, ++clustered) axes.fH.push_back(1 - d);
    std::sort(axes.fH.begin(), axes.fH.end());
  }
  axes.zN = zList;
  axes.nEta = nEta;
  axes.nLnRho = nLnRho;
  axes.lnRhoMin = std::log(rhoMin);
  axes.lnRhoMax = std::log(rhoMax);
  axes.reactingSpecies = sp.reactingNames;

  const std::size_t nFH = axes.fH.size(), nZN = axes.zN.size();
  axes.eLo.assign(nFH * nZN, 0.0);
  axes.eHi.assign(nFH * nZN, 0.0);
  for (std::size_t i = 0; i < nFH; ++i)
    for (std::size_t j = 0; j < nZN; ++j) {
      const double zh = axes.fH[i] * (1 - axes.zN[j]), zo = (1 - axes.fH[i]) * (1 - axes.zN[j]);
      auto y = productsOf(sp.idx, sp.ns, zh, zo, axes.zN[j]);
      axes.eLo[axes.fzIndex(i, j)] = probe.internalEnergy(200.0, 1.0, y.data());
      axes.eHi[axes.fzIndex(i, j)] = probe.internalEnergy(6000.0, 1.0, y.data());
    }

  const std::size_t totalNodes = nFH * nZN * static_cast<std::size_t>(nEta) * static_cast<std::size_t>(nLnRho);
  const std::vector<double> sLower = lowerFractions(axes.lowerGeometric, axes.lowerUniform);

  struct ThreadStats {
    double worstTurnBack = 0;
    std::array<std::size_t, 4> worstTurnBackLoc{};
    std::size_t nGateFailures = 0;
    std::array<std::size_t, 4> worstNGateLoc{};
    double worstDT2b = 0;
    std::array<std::size_t, 4> worstDT2bLoc{};
    std::size_t frozenLower = 0, frozenUpper = 0, integratorFailures = 0;
    std::size_t degenerateCount = 0, zeroSpanCount = 0, outsideCount = 0;
    double worstRaw = 0, worstProjected = 0, worstSum = 0;
    std::array<std::size_t, 4> worstProjectedLoc{};
    std::size_t branchCount = 0;
    double branchTimeSumMs = 0, branchTimeWorstMs = 0;
  };

  auto buildAttempt = [&](int candidate, ThreadStats& merged) {
    ManifoldTable table = axes;
    table.provenance.clear();
    table.nodes.assign(totalNodes, ManifoldTable::Node{});
    table.points.assign(totalNodes * static_cast<std::size_t>(table.pointsPerNode()) * table.pointWidth(), 0.0);

    std::atomic<std::size_t> next{0};
    std::mutex report;
    std::vector<std::thread> pool;
    const auto clock0 = std::chrono::steady_clock::now();
    auto run = [&]() {
      thermo::ReactionSource source(mechanism);
      thermo::ReactionStep step(source, 1e-10, 1e-30);
      ThreadStats ts;
      for (std::size_t idx; (idx = next.fetch_add(1)) < totalNodes;) {
        std::size_t rem = idx;
        const std::size_t l = rem % static_cast<std::size_t>(nLnRho);
        rem /= static_cast<std::size_t>(nLnRho);
        const std::size_t k = rem % static_cast<std::size_t>(nEta);
        rem /= static_cast<std::size_t>(nEta);
        const std::size_t j = rem % nZN;
        const std::size_t i = rem / nZN;
        const double eta = table.eta(k), rho = std::exp(table.lnRho(l));
        const auto r = buildNode(source, step, sp, candidate, maxT, sLower, table.fH[i], table.zN[j], eta, rho,
                                  table.eLo[table.fzIndex(i, j)], table.eHi[table.fzIndex(i, j)], table, idx);
        if (nodeReport) {
          std::string why = r.failure;
          for (char& ch : why) if (ch == ',' || ch == '\n' || ch == '"') ch = ' ';
          std::lock_guard<std::mutex> lock(report);
          std::fprintf(nodeReport, "%d,%.6f,%.4f,%.5f,%.6e,%.3f,%.9e,%.9e,%.9e,%d,%d,%d,%d,%d,%d,%.6e,%d,%d,%s,%s,%.3e,%.3e,%s\n",
                       candidate, table.fH[i], table.zN[j], eta, rho, r.tEq, r.c0, r.cEq, r.cCp, r.degenerate, r.zeroSpan,
                       r.outside, r.lowerFrozen, r.upperFrozen, r.integratorFailure, r.turnBack, r.nGateChecked, r.nGateOk,
                       r.hasDT2bLower ? std::to_string(r.dT2bLower).c_str() : "",
                       r.hasDT2bUpper ? std::to_string(r.dT2bUpper).c_str() : "", r.rawElementError, r.sumError,
                       why.c_str());
        }
        if (r.turnBack > ts.worstTurnBack) { ts.worstTurnBack = r.turnBack; ts.worstTurnBackLoc = {i, j, k, l}; }
        if (r.nGateChecked && !r.nGateOk) { ++ts.nGateFailures; ts.worstNGateLoc = {i, j, k, l}; }
        if (r.hasDT2bLower && r.dT2bLower > ts.worstDT2b) { ts.worstDT2b = r.dT2bLower; ts.worstDT2bLoc = {i, j, k, l}; }
        if (r.hasDT2bUpper && r.dT2bUpper > ts.worstDT2b) { ts.worstDT2b = r.dT2bUpper; ts.worstDT2bLoc = {i, j, k, l}; }
        ts.frozenLower += r.lowerFrozen;
        ts.frozenUpper += r.upperFrozen;
        ts.integratorFailures += r.integratorFailure;
        ts.degenerateCount += r.degenerate;
        ts.zeroSpanCount += r.zeroSpan;
        ts.outsideCount += r.outside;
        ts.worstRaw = std::max(ts.worstRaw, r.rawElementError);
        if (r.projectedElementError > ts.worstProjected) { ts.worstProjected = r.projectedElementError; ts.worstProjectedLoc = {i, j, k, l}; }
        ts.worstSum = std::max(ts.worstSum, r.sumError);
        ts.branchCount += static_cast<std::size_t>(r.branchesIntegrated);
        ts.branchTimeSumMs += r.branchTimeSumMs;
        ts.branchTimeWorstMs = std::max(ts.branchTimeWorstMs, r.branchTimeWorstMs);
      }
      std::lock_guard<std::mutex> lock(report);
      if (ts.worstTurnBack > merged.worstTurnBack) { merged.worstTurnBack = ts.worstTurnBack; merged.worstTurnBackLoc = ts.worstTurnBackLoc; }
      merged.nGateFailures += ts.nGateFailures;
      if (ts.nGateFailures > 0) merged.worstNGateLoc = ts.worstNGateLoc;
      if (ts.worstDT2b > merged.worstDT2b) { merged.worstDT2b = ts.worstDT2b; merged.worstDT2bLoc = ts.worstDT2bLoc; }
      merged.frozenLower += ts.frozenLower;
      merged.frozenUpper += ts.frozenUpper;
      merged.integratorFailures += ts.integratorFailures;
      merged.degenerateCount += ts.degenerateCount;
      merged.zeroSpanCount += ts.zeroSpanCount;
      merged.outsideCount += ts.outsideCount;
      merged.worstRaw = std::max(merged.worstRaw, ts.worstRaw);
      if (ts.worstProjected > merged.worstProjected) { merged.worstProjected = ts.worstProjected; merged.worstProjectedLoc = ts.worstProjectedLoc; }
      merged.worstSum = std::max(merged.worstSum, ts.worstSum);
      merged.branchCount += ts.branchCount;
      merged.branchTimeSumMs += ts.branchTimeSumMs;
      merged.branchTimeWorstMs = std::max(merged.branchTimeWorstMs, ts.branchTimeWorstMs);
    };
    for (int w = 0; w < threads; ++w) pool.emplace_back(run);
    for (auto& t : pool) t.join();
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();
    return std::make_pair(std::move(table), wall);
  };

  auto locStr = [&](const std::array<std::size_t, 4>& loc) {
    char buf[160];
    std::snprintf(buf, sizeof buf, "f_H %.5f Z_N %.4f eta %.5f rho %.4e", axes.fH[loc[0]], axes.zN[loc[1]], axes.eta(loc[2]),
                  std::exp(axes.lnRho(loc[3])));
    return std::string(buf);
  };

  const char* candidateNames[3] = {"Y_H2O", "Y_H2O + Y_HO2 + Y_H2O2", "the O2 consumed"};
  int chosen = -1;
  ManifoldTable table;
  ThreadStats stats;
  double buildWall = 0;
  for (int candidate = 0; candidate < 3; ++candidate) {
    ThreadStats attemptStats;
    auto [attemptTable, wall] = buildAttempt(candidate, attemptStats);
    const bool turnBackOk = attemptStats.worstTurnBack <= kTurnBackGate;
    const bool nGateOk = attemptStats.nGateFailures == 0;
    std::printf("candidate %d (%s): worst turn-back %.3e at %s; N-gate failures %zu%s; %.0f s\n", candidate,
                candidateNames[candidate], attemptStats.worstTurnBack, locStr(attemptStats.worstTurnBackLoc).c_str(),
                attemptStats.nGateFailures, attemptStats.nGateFailures ? (" (worst at " + locStr(attemptStats.worstNGateLoc) + ")").c_str() : "",
                wall);
    std::printf("  its 2(b) worst |dT| %.3f K at %s; 2(c) frozen lower %zu, upper %zu, integrator failures %zu; "
                "criterion 1 |sum Y - 1| %.2e, elements %.2e after projection\n",
                attemptStats.worstDT2b, locStr(attemptStats.worstDT2bLoc).c_str(), attemptStats.frozenLower,
                attemptStats.frozenUpper, attemptStats.integratorFailures, attemptStats.worstSum, attemptStats.worstProjected);
    if (turnBackOk && nGateOk) {
      chosen = candidate;
      table = std::move(attemptTable);
      stats = attemptStats;
      buildWall = wall;
      break;
    }
    if (candidate == 2) {
      table = std::move(attemptTable);
      stats = attemptStats;
      buildWall = wall;
    }
  }

  const bool gate2a = chosen >= 0;
  std::printf("gate 2(a) (progress variable): %s\n", gate2a ? ("pass, using " + std::string(candidateNames[chosen])).c_str() : "FAIL (no candidate passed on every node)");

  const bool gate2b = stats.worstDT2b <= kDTGate2b;
  std::printf("gate 2(b) (equilibrium end state vs equilibrateUV): worst |dT| %.3f K at %s (limit %.1f K): %s\n", stats.worstDT2b,
              locStr(stats.worstDT2bLoc).c_str(), kDTGate2b, gate2b ? "pass" : "FAIL");

  const bool gate2c = stats.integratorFailures == 0;
  std::printf("gate 2(c): frozen branches: lower %zu, upper %zu (reported, not required zero); integrator failures %zu (must be 0): %s\n",
              stats.frozenLower, stats.frozenUpper, stats.integratorFailures, gate2c ? "pass" : "FAIL");

  const bool criterion1 = stats.worstProjected <= 1e-12 && stats.worstSum <= 1e-12;
  std::printf("criterion 1 (node conservation): element error %.2e raw, %.2e after projection (worst at %s), |sum Y - 1| %.2e "
              "(limit 1e-12): %s\n",
              stats.worstRaw, stats.worstProjected, locStr(stats.worstProjectedLoc).c_str(), stats.worstSum,
              criterion1 ? "pass" : "FAIL");

  std::printf("nodes: %zu degenerate, %zu with a zero-span branch, %zu outside the thermo range\n", stats.degenerateCount,
              stats.zeroSpanCount, stats.outsideCount);
  if (stats.branchCount > 0)
    std::printf("branch integration: %zu branches, mean %.2f ms, worst %.2f ms\n", stats.branchCount,
                stats.branchTimeSumMs / static_cast<double>(stats.branchCount), stats.branchTimeWorstMs);

  const bool allGatesPass = gate2a && gate2b && gate2c && criterion1;
  if (nodeReport) std::fclose(nodeReport);

  char text[1536];
  const std::time_t now = std::time(nullptr);
  char date[64];
  std::strftime(date, sizeof date, "%Y-%m-%d %H:%M %Z", std::localtime(&now));
  std::string zStr;
  for (double v : axes.zN) zStr += (zStr.empty() ? "" : ",") + std::to_string(v);
  std::snprintf(text, sizeof text,
                "Table B built %s by crucible_build_manifold_table from %s, Cantera %s; f_H %zu nodes (%d lean, %d rich, "
                "%zu clustered), Z_N {%s}, eta %d over [0,1], ln rho %d over [%.3g, %.3g] kg/m3; progress variable %s; "
                "%zu nodes in %.0f s on %d threads; node element error %.2e raw, %.2e after projection, |sum Y - 1| %.2e; "
                "gates: 2a %s, 2b %s, 2c %s, criterion 1 %s%s",
                date, mechanism.c_str(), Cantera::version().c_str(), nFH, fLean, fRich, clustered, zStr.c_str(), nEta, nLnRho,
                rhoMin, rhoMax, chosen >= 0 ? candidateNames[chosen] : candidateNames[2], totalNodes, buildWall, threads,
                stats.worstRaw, stats.worstProjected, stats.worstSum, gate2a ? "pass" : "FAIL", gate2b ? "pass" : "FAIL",
                gate2c ? "pass" : "FAIL", criterion1 ? "pass" : "FAIL", allGatesPass ? "" : " -- GATE FAILURE");
  table.provenance = text;
  std::printf("%s\n", text);

  if (!allGatesPass && !writeAnyway) {
    std::printf("Not writing %s: a gate failed. Set CRUCIBLE_TABLE_B_WRITE_ANYWAY=1 to write it anyway for inspection.\n",
                out.c_str());
    return 1;
  }
  if (!allGatesPass) std::printf("WRITING ANYWAY FOR INSPECTION (a gate failed; this table is not for use).\n");
  table.save(out);
  std::printf("wrote %s (%.1f MB)\n", out.c_str(), table.points.size() * 8.0 / 1e6);
  return allGatesPass ? 0 : 1;
}
