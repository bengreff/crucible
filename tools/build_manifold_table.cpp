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
// The progress variable is chosen from four candidates in order (Y_H2O; Y_H2O + Y_HO2 + Y_H2O2; the
// O2 consumed; the specific entropy s, TABLE_B.md "The entropy candidate in the builder"), each tried
// as a full build (CRUCIBLE_TABLE_B_CANDIDATES=<comma list> selects which, in that order): gate 2(a) (turn-back) and the N gate must
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
// CRUCIBLE_TABLE_B_FH_NODES=<comma list> gives the f_H nodes directly (the lean, rich and clustering
// arguments are then ignored), for small tables that measure the gates over a chosen range.
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

constexpr int kEntropy = 3;  // the specific entropy as progress variable

namespace {

constexpr double kTEnd = 1e15;            // a branch's integration horizon (TABLE_B.md)
constexpr double kStopFrac = 1e-6;        // a branch stops within this fraction of its span of c_eq
constexpr double kTurnBackWindow = 1e-4;  // turn-back is measured only until c is this close to c_eq
constexpr double kTurnBackGate = 1e-9;    // above this, a turn-back fails gate 2(a) (CVODES noise floor)
constexpr double kTurnBackFloor = 1e-14;  // a fall within this fraction of |c_eq| is c's round-off, not counted (TABLE_B.md, 23 October)
constexpr double kZeroSpanTol = 1e-12;    // |c_eq - c_start| below this: the branch is a point
constexpr double kPhiA = 0.005;           // c_a: the lower branch's 0.5% fall in N
constexpr double kDTGate2b = 3.0;         // K, gate 2(b)
constexpr double kE2DT = 3.0, kE2DY = 0.01;  // E2's bars (TABLE_B.md, "The specific entropy's checks")
constexpr double kRefineFrac = 2.0 / 3.0;     // the refinement's target, a fraction of E2's bars (TABLE_B.md)
constexpr std::size_t kMaxPointsPerNode = 400;  // the refinement's cap (TABLE_B.md, "per-node refinement")
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
// CRUCIBLE_TABLE_B_E2_DUMP=<path>: one line per E2 pair (a diagnostic, TABLE_B.md "E2 measured").
std::FILE* gE2Dump = nullptr;
std::mutex gE2DumpLock;

struct RawStep {
  double t = 0, T = 0, c = 0, lambda = 0;
  std::vector<double> y;  // full species, size ns
};

// A stored point while a node is built: full species (projected), lambda, c.
struct CurvePoint {
  std::vector<double> y;
  double lambda = 0, c = 0;
};

// A node's final points in the table's per-point layout, and its lower branch's count.
struct NodeCurve {
  std::vector<double> points;
  std::uint32_t lowerCount = 0;
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
  bool hasE2 = false;  // E2 (s only): the table's T and Y between stored points against the trajectory
  // E2 after refinement: the final gaps' 1/4, 1/2 and 3/4 points (the refinement criterion, so under
  // the bars unless the node hit the cap) and the independent eighth points, over every gap and the
  // last quarter of the span; the points the refinement added.
  double e2MidDT = 0, e2MidDY = 0;
  double e2QDT = 0, e2QDY = 0, e2QDTLast = 0, e2QDYLast = 0;
  std::size_t addedPoints = 0;
  bool e2Capped = false;
  std::size_t e2RawFallbacks = 0;  // searches resolved from the raw steps (at s's noise floor)
  double e2FallbackSpan = 0;       // the largest span of a branch with such a fallback
  double spanLower = 0, spanUpper = 0;  // E3: |c_eq - c_start| per branch
  bool e2Failure = false;               // an integration to a midpoint threw (reported, not a 2(c) failure)
  bool nGateChecked = false, nGateOk = true;
  double fallWithinFloor = 0;  // the largest fall the round-off floor kept from counting, over |c_eq|
  bool floorDecided = false;   // such a fall was above kTurnBackGate of its span (the old gate would fail it)
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
NodeOutcome buildNode(thermo::ReactionSource& source, thermo::ReactionStep& step, Cantera::ThermoPhase& th,
                       const SpeciesTable& sp, int candidate, double maxThermoT, const std::vector<double>& sLower,
                       double f, double z, double eta, double rho, double eLo, double eHi,
                       ManifoldTable& out, std::size_t nodeIdx, NodeCurve& curve) {
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
  // c at z = [T, Y...]. The entropy is Cantera's, the integrator's own thermo.
  std::vector<double> mu(ns);
  auto setState = [&](const double* z) {
    th.setMassFractions_NoNorm(z + 1);
    th.setState_TD(z[0], rho);
  };
  auto progress = [&](const double* z) -> double {
    const double* y = z + 1;
    switch (candidate) {
      case 0: return y[sp.iH2O];
      case 1: return y[sp.iH2O] + y[sp.iHO2] + y[sp.iH2O2];
      case 2: return zo - y[sp.iO2];
      default: setState(z); return th.entropy_mass();
    }
  };
  // dc/dt at z, dzdt = rates at z (size ns+1, [0] = dT/dt). For s, the Gibbs relation at fixed u and rho:
  // ds/dt = -sum_k mu_k (dY_k/dt) / (W_k T).
  auto dcdt = [&](const double* z, const double* dzdt) -> double {
    switch (candidate) {
      case 0: return dzdt[1 + sp.iH2O];
      case 1: return dzdt[1 + sp.iH2O] + dzdt[1 + sp.iHO2] + dzdt[1 + sp.iH2O2];
      case 2: return -dzdt[1 + sp.iO2];
      default: {
        setState(z);
        th.getChemPotentials(mu.data());
        double sum = 0;
        for (std::size_t q = 0; q < ns; ++q) sum += mu[q] * dzdt[1 + q] / sp.weight[q];
        return -sum / z[0];
      }
    }
  };
  // The stored points' fractions of a branch's span: linear for the linear candidates; for s, uniform in
  // the distance to equilibrium, which is quadratic in s (E2's node rule).
  auto lowerFrac = [&](double xi) { return candidate == kEntropy ? 1 - (1 - xi) * (1 - xi) : xi; };
  auto upperFrac = [&](double u) { return candidate == kEntropy ? u * u : u; };

  const int lowerCount = out.lowerPoints();
  const int upperCount = out.upperUniform;
  // The node's points after refinement; empty until then, when the base points are the curve.
  std::vector<CurvePoint> lowerPts, upperPts;
  bool refined = false;
  auto fullY = [&](int p) {
    std::vector<double> y(ns, 0.0);
    const double* py = out.pointY(nodeIdx, p);
    for (std::size_t q = 0; q < ns; ++q)
      y[q] = sp.reactingPosition[q] != kNotReacting ? py[sp.reactingPosition[q]] : (q == sp.idx.n2 ? zn : 0.0);
    return y;
  };
  auto basePoint = [&](int p) { return CurvePoint{fullY(p), out.lambda(nodeIdx, p), out.c(nodeIdx, p)}; };
  auto writeCurve = [&]() {
    if (!refined) {
      lowerPts.clear();
      upperPts.clear();
      for (int p = 0; p < lowerCount; ++p) lowerPts.push_back(basePoint(p));
      for (int p = lowerCount; p < lowerCount + upperCount; ++p) upperPts.push_back(basePoint(p));
    }
    const std::size_t width = out.pointWidth();
    curve.lowerCount = static_cast<std::uint32_t>(lowerPts.size());
    curve.points.assign((lowerPts.size() + upperPts.size()) * width, 0.0);
    std::size_t at = 0;
    for (const auto* list : {&lowerPts, &upperPts})
      for (const auto& pt : *list) {
        double* dst = curve.points.data() + at * width;
        for (std::size_t q = 0; q < ns; ++q)
          if (sp.reactingPosition[q] != kNotReacting) dst[sp.reactingPosition[q]] = pt.y[q];
        dst[out.reactingSpecies.size()] = pt.lambda;
        dst[out.reactingSpecies.size() + 1] = pt.c;
        ++at;
      }
  };

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
      for (int p = 0; p < out.basePointsPerNode(); ++p) {
        double* py = out.pointY(nodeIdx, p);
        for (std::size_t q = 0; q < ns; ++q)
          if (sp.reactingPosition[q] != kNotReacting) py[sp.reactingPosition[q]] = zstate[q + 1];
        out.lambda(nodeIdx, p) = 0.0;
        out.c(nodeIdx, p) = 0.0;
      }
      writeCurve();
      return result;
    }

    const std::vector<double> yReactants = reactantsOf(sp.idx, ns, zh, zo, zn);
    const std::vector<double> yComplete = productsOf(sp.idx, ns, zh, zo, zn);

    // The lower branch's start.
    std::vector<double> zLower(ns + 1);
    double lowerTStar = std::numeric_limits<double>::quiet_NaN();  // NaN: the reactants' start
    double lowerAlpha = std::numeric_limits<double>::quiet_NaN();  // the reactants' share of a start beyond the frozen family
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
      // The bracket is checked before bisecting (TABLE_B.md, *The integrator failures*). At or below
      // the family's low end (eta = 0 carries the products at 200 K exactly), the equilibrium at 300 K
      // differs from the complete-combustion products only by trace dissociation, so the products
      // are the start; without this the bisection drifts to its upper bound.
      const auto g = [&](double ts) { return energyAt(200.0, frozenAt(ts)) - e; };
      if (g(300.0) >= -1e-9 * (eHi - eLo)) {
        zLower[0] = temperatureAt(e, yComplete);
        std::copy(yComplete.begin(), yComplete.end(), zLower.begin() + 1);
      } else if (g(6000.0) >= 0) {
        const double tStar = bisect(g, 300.0, 6000.0);
        lowerTStar = tStar;
        auto y = frozenAt(tStar * (1 - 1e-9));
        zLower[0] = temperatureAt(e, y);
        std::copy(y.begin(), y.end(), zLower.begin() + 1);
      } else {
        // No frozen start up to 6000 K (the thermo's end) carries e (dense nodes, where the 6000 K
        // equilibrium barely dissociates): the family continues to the reactants, mixed with the 6000 K
        // frozen composition at 200 K (TABLE_B.md, 23 October, fix 1). Both ends carry the node's elements.
        const std::vector<double> yHot = frozenAt(6000.0);
        auto mixAt = [&](double a) {
          std::vector<double> y(ns);
          for (std::size_t q = 0; q < ns; ++q) y[q] = (1 - a) * yHot[q] + a * yReactants[q];
          return y;
        };
        const double alpha = bisect([&](double a) { return energyAt(200.0, mixAt(a)) - e; }, 0.0, 1.0);
        lowerTStar = 6000.0;
        lowerAlpha = alpha;
        auto y = mixAt(alpha * (1 - 1e-9));
        zLower[0] = temperatureAt(e, y);
        std::copy(y.begin(), y.end(), zLower.begin() + 1);
      }
    }

    // The shared equilibrium (both branches meet here).
    std::vector<double> zEq = zLower;
    source.equilibrateUV(rho, zEq.data());
    const double cEq = progress(zEq.data());
    result.tEq = zEq[0];
    if (zEq[0] > maxThermoT) result.outside = true;

    // The upper branch's start: complete combustion.
    std::vector<double> zUpper(ns + 1);
    zUpper[0] = temperatureAt(e, yComplete);
    std::copy(yComplete.begin(), yComplete.end(), zUpper.begin() + 1);
    const double cCp = progress(zUpper.data());
    const double c0 = progress(zLower.data());
    node.c0 = c0;
    node.cEq = cEq;
    node.cCp = cCp;
    result.c0 = c0;
    result.cEq = cEq;
    result.cCp = cCp;
    if (zLower[0] > maxThermoT || zUpper[0] > maxThermoT) result.outside = true;

    const double spanLower = std::fabs(cEq - c0), spanUpper = std::fabs(cCp - cEq);
    result.spanLower = spanLower;
    result.spanUpper = spanUpper;
    // For s, a zero span is one below 1e-9 |s_eq| (E3's scale); for the linear candidates, 1e-12.
    const double zeroSpanTol = candidate == kEntropy ? 1e-9 * std::fabs(cEq) : kZeroSpanTol;
    const bool zeroSpanLower = spanLower < zeroSpanTol, zeroSpanUpper = spanUpper < zeroSpanTol;
    if (zeroSpanLower || zeroSpanUpper) result.zeroSpan = true;

    // Integrates one branch from zStart toward branchCEq, recording every internal step (with its
    // exact rate, hence lambda) until within kStopFrac of the span or the time horizon.
    auto runBranch = [&](const char* name, const std::vector<double>& zStart, double branchCEq, bool& frozenOut) {
      std::vector<RawStep> steps;
      const double cStart = progress(zStart.data());
      const double span = std::fabs(branchCEq - cStart);
      // An s-branch stops on composition (TABLE_B.md, "The entropy candidate in the builder"): within
      // 1e-6 of the start's largest species distance from the equilibrium.
      double yDistStart = 0;
      for (std::size_t q = 0; q < ns; ++q) yDistStart = std::max(yDistStart, std::fabs(zStart[q + 1] - zEq[q + 1]));
      auto closeEnough = [&](const double* zz, double c) {
        if (candidate != kEntropy) return std::fabs(branchCEq - c) <= kStopFrac * span;
        double d = 0;
        for (std::size_t q = 0; q < ns; ++q) d = std::max(d, std::fabs(zz[q + 1] - zEq[q + 1]));
        return d <= kStopFrac * yDistStart;
      };
      std::vector<double> dzdt(ns + 1);
      auto withLambda = [&](const double* z, double t) {
        source.rates(rho, z, dzdt.data());
        const double c = progress(z);
        const double denom = branchCEq - c;
        RawStep s;
        s.t = t;
        s.T = z[0];
        s.c = c;
        s.lambda = denom != 0.0 ? dcdt(z, dzdt.data()) / denom : 0.0;
        s.y.assign(z + 1, z + 1 + ns);
        return s;
      };
      steps.push_back(withLambda(zStart.data(), 0.0));
      frozenOut = true;
      std::vector<double> z = zStart;
      const auto clock0 = std::chrono::steady_clock::now();
      try {
        step.follow(rho, z.data(), kTEnd, [&](double t, const double* zz) {
          steps.push_back(withLambda(zz, t));
          const bool close = closeEnough(zz, steps.back().c);
          if (close) frozenOut = false;
          return !close;
        });
      } catch (const std::exception& ex) {
        // Tag the failing branch with its start (TABLE_B.md, 2(c)): which branch, how it started,
        // the start's T and first right-hand side, the furthest step reached, and the start's Y.
        source.rates(rho, zStart.data(), dzdt.data());
        std::string why = std::string(name) + " branch";
        char buf[256];
        if (name[0] == 'l')
          why += std::isnan(lowerTStar) ? " (reactants start)"
                 : !std::isnan(lowerAlpha) ? (std::snprintf(buf, sizeof buf, " (6000 K frozen start mixed with the reactants, alpha %.6g)", lowerAlpha), buf)
                                           : (std::snprintf(buf, sizeof buf, " (frozen start, T* %.6g K)", lowerTStar), buf);
        std::snprintf(buf, sizeof buf, ": start T %.6g K, dT/dt %.4e K/s, c %.6e to c_eq %.6e; %zu steps recorded, last t %.4e s T %.6g K; ",
                      zStart[0], dzdt[0], cStart, branchCEq, steps.size(), steps.back().t, steps.back().T);
        why += buf;
        why += ex.what();
        why += "; start Y";
        for (std::size_t q = 0; q < ns; ++q) { std::snprintf(buf, sizeof buf, " %.4e", zStart[q + 1]); why += buf; }
        why += "; start dY/dt";
        for (std::size_t q = 0; q < ns; ++q) { std::snprintf(buf, sizeof buf, " %.4e", dzdt[q + 1]); why += buf; }
        throw std::runtime_error(why);
      }
      const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - clock0).count();
      result.branchTimeSumMs += ms;
      result.branchTimeWorstMs = std::max(result.branchTimeWorstMs, ms);
      ++result.branchesIntegrated;
      return steps;
    };

    std::vector<RawStep> lowerSteps, upperSteps;
    bool lowerFrozen = false, upperFrozen = false;
    if (!zeroSpanLower) lowerSteps = runBranch("lower", zLower, cEq, lowerFrozen);
    if (!zeroSpanUpper) upperSteps = runBranch("upper", zUpper, cEq, upperFrozen);
    result.lowerFrozen = lowerFrozen;
    result.upperFrozen = upperFrozen;
    if (lowerFrozen || upperFrozen) node.flags |= ManifoldTable::kFrozen;

    for (auto& s : lowerSteps) if (s.T > maxThermoT) result.outside = true;
    for (auto& s : upperSteps) if (s.T > maxThermoT) result.outside = true;
    node.flags |= (result.outside ? ManifoldTable::kOutside : 0u) | (result.zeroSpan ? ManifoldTable::kZeroSpan : 0u) |
                  (result.degenerate ? ManifoldTable::kDegenerate : 0u);

    // Turn-back (gate 2a): the largest fall against the approach, relative to the span, measured
    // only until c is within kTurnBackWindow of the span from c_eq. A fall within kTurnBackFloor of
    // |c_eq| is c's own round-off and is not counted, but is reported (TABLE_B.md, 23 October, fix 2).
    auto turnBackOf = [&](const std::vector<RawStep>& steps, double branchCEq, double span) {
      if (steps.empty() || span < kZeroSpanTol) return 0.0;
      const double sign = branchCEq >= steps.front().c ? 1.0 : -1.0;
      const double floor = kTurnBackFloor * std::fabs(branchCEq);
      double best = steps.front().c, worst = 0.0;
      for (const auto& s : steps) {
        if (std::fabs(branchCEq - s.c) <= kTurnBackWindow * span) break;
        if (sign * (s.c - best) > 0) best = s.c;
        const double fall = sign * (best - s.c);
        if (fall > floor) {
          worst = std::max(worst, fall / span);
        } else if (fall > 0) {
          result.fallWithinFloor = std::max(result.fallWithinFloor, fall / std::fabs(branchCEq));
          if (fall / span > kTurnBackGate) result.floorDecided = true;
        }
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

    auto storePoint = [&](int p, std::vector<double> y, double lambdaVal, double cVal) {
      out.c(nodeIdx, p) = cVal;
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
      if (zeroSpanLower) { storePoint(p, eqY, 0.0, cEq); continue; }
      const double c = c0 + lowerFrac(sLower[p]) * (cReachLower - c0);
      if (p == lowerCount - 1 && !lowerFrozen) {
        storePoint(p, eqY, 0.0, cEq);  // lambda fixed below (copied from the point before)
      } else {
        std::vector<double> y;
        double lam;
        atC(lowerSteps, c, y, lam);
        storePoint(p, y, lam, c);
      }
    }
    if (!zeroSpanLower && !lowerFrozen) out.lambda(nodeIdx, lowerCount - 1) = out.lambda(nodeIdx, lowerCount - 2);

    const double cReachUpper = zeroSpanUpper ? cEq : (upperFrozen ? upperSteps.back().c : cEq);
    for (int p = 0; p < upperCount; ++p) {
      const int idx = lowerCount + p;
      if (zeroSpanUpper) { storePoint(idx, eqY, 0.0, cEq); continue; }
      const double c = cReachUpper + upperFrac(upperCount > 1 ? static_cast<double>(p) / (upperCount - 1) : 1.0) * (cCp - cReachUpper);
      if (p == 0 && !upperFrozen) {
        storePoint(idx, eqY, 0.0, cEq);  // lambda fixed below
      } else {
        std::vector<double> y;
        double lam;
        atC(upperSteps, c, y, lam);
        storePoint(idx, y, lam, c);
      }
    }
    if (!zeroSpanUpper && !upperFrozen) out.lambda(nodeIdx, lowerCount) = out.lambda(nodeIdx, lowerCount - 1);

    // E2 and per-node refinement (s only; TABLE_B.md, "per-node refinement" and "refinement on three
    // interior points"). The stored Y are interpolated linearly in the distance d = sqrt(s_ref - s),
    // s_ref the branch's last point (s_eq, or a frozen branch's reach). A gap is refined if any of its
    // points at d 1/4, 1/2 and 3/4 of the way differs from the trajectory's own state there by over 3 K
    // or 0.01 in Y; the worst of them, the trajectory's state, becomes a new point, until every gap
    // passes or the node holds kMaxPointsPerNode points. The refinement's target is kRefineFrac of the
    // bars; the final gaps are then checked at their eighth points, which it did not use, at the full bars.
    for (int p = 0; p < lowerCount; ++p) lowerPts.push_back(basePoint(p));
    for (int p = lowerCount; p < lowerCount + upperCount; ++p) upperPts.push_back(basePoint(p));
    refined = true;
    if (candidate == kEntropy) {
      auto searchAt = [&](const std::vector<RawStep>& steps, std::size_t i, double sm, double tol, std::vector<double>& z) {
        const RawStep& a = steps[i - 1];
        std::vector<double> za(ns + 1);
        za[0] = a.T;
        std::copy(a.y.begin(), a.y.end(), za.begin() + 1);
        auto f = [&](double dt) {
          z = za;
          try {
            if (dt > 0) step.advance(rho, z.data(), dt);
          } catch (const std::exception& ex) {
            char buf[300];
            std::snprintf(buf, sizeof buf, "E2 search for s %.9e: from raw step %zu of %zu (t %.6e, c %.9e, lambda %.3e) "
                          "over %.3e s: ", sm, i, steps.size(), a.t, a.c, a.lambda, dt);
            throw std::runtime_error(buf + std::string(ex.what()));
          }
          return progress(z.data()) - sm;
        };
        // The integration restarts at t = 0 from step a's state. After a long induction the raw steps' own
        // times can stand still (t + h == t at ignition; TABLE_B.md, "a time coordinate that stands
        // still"), so where the recorded interval is not positive the first trial comes from the rate.
        double hi = steps[i].t - a.t;
        if (!(hi > 0)) {
          const double rate = a.lambda * (cEq - a.c);
          hi = rate > 0 ? (sm - a.c) / rate : 1e-15;
          if (!(hi > 0)) hi = 1e-15;
        }
        double lo = 0, flo = a.c - sm, fhi = f(hi);
        for (int grow = 0; fhi < 0 && grow < 100; ++grow) { lo = hi; flo = fhi; hi *= 2; fhi = f(hi); }
        int side = 0;
        const double fhi0 = fhi, hi0 = hi;
        for (int it = 0; it < 60; ++it) {
          const double dt = (lo * fhi - hi * flo) / (fhi - flo);
          const double fm = f(dt);
          if (std::fabs(fm) <= tol) return;
          if (fm < 0) { lo = dt; flo = fm; if (side == -1) fhi *= 0.5; side = -1; }
          else { hi = dt; fhi = fm; if (side == 1) flo *= 0.5; side = 1; }
        }
        // Not converged: a state at another s must not become a point (TABLE_B.md, "the search").
        char buf[400];
        std::snprintf(buf, sizeof buf,
                      "E2 search for s %.9e not converged: raw step %zu of %zu at t %.6e (c %.9e to %.9e, dt %.3e; %zu raw "
                      "steps share a time with the next), first "
                      "f(hi) %.3e, bracket [%.6e, %.6e] f [%.3e, %.3e], tol %.3e",
                      sm, i, steps.size(), a.t, a.c, steps[i].c, hi0,
                      static_cast<std::size_t>(std::count_if(steps.begin() + 1, steps.end(),
                                                             [&](const RawStep& r) { return r.t == (&r - 1)->t; })),
                      fhi0, lo, hi, flo, fhi, tol);
        throw std::runtime_error(buf);
      };
      // At s's integration noise floor a search can fail (TABLE_B.md, "a counted fallback to the raw
      // steps"); the state then comes from the bracketing raw steps, linear in s, and is counted.
      double fallbackSpan = 0;
      auto trajectoryAt = [&](const std::vector<RawStep>& steps, double sm, double tol, std::vector<double>& z) {
        std::size_t i = 1;
        while (i < steps.size() && steps[i].c < sm) ++i;
        if (i >= steps.size()) i = steps.size() - 1;
        try {
          searchAt(steps, i, sm, tol, z);
        } catch (const std::exception&) {
          const RawStep& a = steps[i - 1];
          const RawStep& b = steps[i];
          const double w = b.c != a.c ? std::clamp((sm - a.c) / (b.c - a.c), 0.0, 1.0) : 0.0;
          z.assign(ns + 1, 0.0);
          z[0] = a.T + w * (b.T - a.T);
          for (std::size_t q = 0; q < ns; ++q) z[q + 1] = a.y[q] + w * (b.y[q] - a.y[q]);
          ++result.e2RawFallbacks;
          result.e2FallbackSpan = std::max(result.e2FallbackSpan, fallbackSpan);
        }
      };
      auto differ = [&](const std::vector<double>& y, const std::vector<double>& z, double& dT, double& dY) {
        dT = std::fabs(temperatureAt(e, y) - z[0]);
        dY = 0;
        for (std::size_t q = 0; q < ns; ++q)
          if (sp.reactingPosition[q] != kNotReacting) dY = std::max(dY, std::fabs(y[q] - z[q + 1]));
      };
      std::vector<double> rates(ns + 1);
      std::size_t total = lowerPts.size() + upperPts.size();
      auto refine = [&](std::vector<CurvePoint>& pts, const std::vector<RawStep>& steps, double span, double sRef,
                        int branch) {
        if (steps.size() < 2) return;
        auto dist = [&](double c) { return std::sqrt(std::max(sRef - c, 0.0)); };
        auto tolOf = [&](double ca, double cb) { return std::max(1e-6 * std::fabs(cb - ca), 1e-11 * std::fabs(ca)); };
        std::vector<double> z, ym(ns);
        // The table's Y at fraction w of a gap in d, and the trajectory's state at the same s.
        auto at = [&](std::size_t i, double w, double tol, double& sOut, std::vector<double>& zOut, double& dT,
                      double& dY) {
          const double da = dist(pts[i].c), db = dist(pts[i + 1].c), dq = da + w * (db - da);
          sOut = sRef - dq * dq;
          for (std::size_t q = 0; q < ns; ++q) ym[q] = pts[i].y[q] + w * (pts[i + 1].y[q] - pts[i].y[q]);
          trajectoryAt(steps, sOut, tol, zOut);
          differ(ym, zOut, dT, dY);
        };
        for (std::size_t i = 0; i + 1 < pts.size();) {
          const double ca = pts[i].c, cb = pts[i + 1].c;
          if (ca == cb) { ++i; continue; }
          // The worst of the three interior points, by its excess over the bars.
          double worstExcess = -1, worstDT = 0, worstDY = 0;
          std::vector<double> worstZ;
          for (const double w : {0.25, 0.5, 0.75}) {
            double sw, dT, dY;
            at(i, w, tolOf(ca, cb), sw, z, dT, dY);
            const double excess = std::max(dT / kE2DT, dY / kE2DY) / kRefineFrac;
            if (excess > worstExcess) { worstExcess = excess; worstZ = z; }
            worstDT = std::max(worstDT, dT);
            worstDY = std::max(worstDY, dY);
          }
          const bool fails = worstExcess > 1;
          if (fails && total < kMaxPointsPerNode) {
            CurvePoint n;
            n.y.assign(worstZ.begin() + 1, worstZ.end());
            double rawErr, projErr, sumErr;
            project(n.y, rawErr, projErr, sumErr);
            result.rawElementError = std::max(result.rawElementError, rawErr);
            result.projectedElementError = std::max(result.projectedElementError, projErr);
            result.sumError = std::max(result.sumError, sumErr);
            n.c = progress(worstZ.data());
            source.rates(rho, worstZ.data(), rates.data());
            n.lambda = cEq != n.c ? dcdt(worstZ.data(), rates.data()) / (cEq - n.c) : 0.0;
            pts.insert(pts.begin() + static_cast<std::ptrdiff_t>(i) + 1, std::move(n));
            ++total;
            ++result.addedPoints;
            continue;  // the left part is checked next
          }
          if (fails) result.e2Capped = true;
          result.e2MidDT = std::max(result.e2MidDT, worstDT);
          result.e2MidDY = std::max(result.e2MidDY, worstDY);
          ++i;
        }
        // The independent check: the eighth points of every final gap.
        for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
          const double ca = pts[i].c, cb = pts[i + 1].c;
          if (ca == cb) continue;
          for (const double w : {0.125, 0.375, 0.625, 0.875}) {
            double sq, dT, dY;
            at(i, w, tolOf(ca, cb), sq, z, dT, dY);
            result.e2QDT = std::max(result.e2QDT, dT);
            result.e2QDY = std::max(result.e2QDY, dY);
            if (std::fabs(cEq - sq) <= 0.25 * span) {
              result.e2QDTLast = std::max(result.e2QDTLast, dT);
              result.e2QDYLast = std::max(result.e2QDYLast, dY);
            }
            if (gE2Dump) {
              std::lock_guard<std::mutex> lock(gE2DumpLock);
              std::fprintf(gE2Dump, "%.5f,%.4f,%.5f,%.6e,%d,%zu,%.4f,%.9e,%.9e,%.9e,%.4f,%.4f,%.3e\n", f, zn, eta, rho,
                           branch, i, w, ca, cb, sRef, temperatureAt(e, pts[i].y), dT, dY);
            }
          }
        }
      };
      try {
        result.hasE2 = true;
        fallbackSpan = spanLower;
        if (!zeroSpanLower) refine(lowerPts, lowerSteps, spanLower, lowerPts.back().c, 0);
        fallbackSpan = spanUpper;
        if (!zeroSpanUpper) refine(upperPts, upperSteps, spanUpper, upperPts.front().c, 1);
      } catch (const std::exception& ex) {
        result.e2Failure = true;
        std::lock_guard<std::mutex> lock(gE2DumpLock);
        std::fprintf(stderr, "E2 failure at f_H %.5f Z_N %.4f eta %.5f rho %.4e: %s\n", f, zn, eta, rho, ex.what());
      }
    }

    // The N gate (03:10, 2a): N = sum Y_k/W_k strictly falls from c_a to c_cp along the stored
    // points. c_a is the first stored point past phi_N(c) = 0.5% of the lower branch's own fall.
    if (!result.degenerate) {
      const int finalLower = static_cast<int>(lowerPts.size());
      const int finalCount = finalLower + static_cast<int>(upperPts.size());
      auto nAt = [&](int p) {
        const auto& y = p < finalLower ? lowerPts[p].y : upperPts[p - finalLower].y;
        double m = 0;
        for (std::size_t q = 0; q < ns; ++q)
          if (sp.reactingPosition[q] != kNotReacting) m += y[q] / sp.weight[q];
        return m;
      };
      const double n0 = nAt(0), neq = nAt(finalLower - 1);
      int idxA = finalLower;  // no portion to check if the fall never reaches 0.5%
      if (std::fabs(n0 - neq) > 1e-15) {
        for (int p = 0; p < finalLower; ++p) {
          const double phi = (n0 - nAt(p)) / (n0 - neq);
          if (phi >= kPhiA) { idxA = p; break; }
        }
      }
      result.nGateChecked = true;
      bool ok = true;
      for (int p = idxA; p + 1 < finalCount; ++p) {
        const double na = nAt(p), nb = nAt(p + 1);
        if (!(nb < na + 1e-13 * std::max(1.0, std::fabs(na)))) { ok = false; break; }
      }
      result.nGateOk = ok;
    }
  } catch (const std::exception& ex) {
    result.integratorFailure = true;
    result.failure = ex.what();
    refined = false;
  }
  writeCurve();
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
  auto parseList = [](const std::string& s) {
    std::vector<double> list;
    std::size_t pos = 0;
    while (pos <= s.size()) {
      const std::size_t comma = s.find(',', pos);
      list.push_back(std::atof(s.substr(pos, comma - pos).c_str()));
      if (comma == std::string::npos) break;
      pos = comma + 1;
    }
    std::sort(list.begin(), list.end());
    return list;
  };
  const std::vector<double> zList = parseList(argv[5]);
  const char* fHNodes = std::getenv("CRUCIBLE_TABLE_B_FH_NODES");
  const int nEta = std::atoi(argv[6]), nLnRho = std::atoi(argv[7]);
  const double rhoMin = std::atof(argv[8]), rhoMax = std::atof(argv[9]);
  const int threads = std::max(1, std::atoi(argv[10]));
  const std::string mechanism = argc == 12 ? argv[11] : "h2o2.yaml";
  const bool writeAnyway = std::getenv("CRUCIBLE_TABLE_B_WRITE_ANYWAY") != nullptr;
  const char* nodeReportPath = std::getenv("CRUCIBLE_TABLE_B_NODE_REPORT");
  if (const char* dump = std::getenv("CRUCIBLE_TABLE_B_E2_DUMP")) {
    gE2Dump = std::fopen(dump, "w");
    std::fprintf(gE2Dump, "f_H,Z_N,eta,rho,branch,gap,w,c_a,c_b,s_ref,T_a,dT_check,dY_check\n");
  }
  std::FILE* nodeReport = nodeReportPath ? std::fopen(nodeReportPath, "w") : nullptr;
  if (nodeReportPath && !nodeReport) throw std::runtime_error(std::string("cannot write ") + nodeReportPath);
  if (nodeReport)
    std::fprintf(nodeReport, "candidate,f_H,Z_N,eta,rho,T_eq,c0,c_eq,c_cp,degenerate,zero_span,outside,lower_frozen,upper_frozen,"
                             "integrator_failure,turn_back,n_gate_checked,n_gate_ok,dT_2b_lower,dT_2b_upper,raw_element_error,"
                             "sum_error,span_lower,span_upper,e2q_dT,e2q_dY,e2q_dT_last_quarter,e2q_dY_last_quarter,e2_mid_dT,e2_mid_dY,added_points,e2_capped,e2_raw_fallbacks,e2_failure,failure,fall_within_floor\n");

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
  if (fHNodes) {
    axes.fH = parseList(fHNodes);
  } else if (fCluster > 0) {
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
    double worstFallWithinFloor = 0;  // fix 2's uncounted falls, over |c_eq|
    std::size_t floorDecidedNodes = 0;
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
    // E2 and E3 (s only).
    double e2QDT = 0, e2QDY = 0, e2QDTLast = 0, e2QDYLast = 0, e2MidDT = 0, e2MidDY = 0;
    std::array<std::size_t, 4> e2QDTLoc{}, e2QDYLoc{};
    std::vector<std::size_t> added;  // per node with E2
    std::size_t e2Failures = 0, e2Capped = 0, e2RawFallbacks = 0, e2FallbackNodes = 0;
    double e2FallbackSpan = 0;
    double minSpanRatio = std::numeric_limits<double>::infinity();
  };

  auto buildAttempt = [&](int candidate, ThreadStats& merged) {
    ManifoldTable table = axes;
    table.provenance.clear();
    table.nodes.assign(totalNodes, ManifoldTable::Node{});
    // The base points are built in a fixed layout; each node's final (refined) points go to curves[node]
    // and are assembled into the variable layout after the build.
    const auto base = static_cast<std::uint32_t>(table.basePointsPerNode());
    for (std::size_t n = 0; n < totalNodes; ++n) {
      table.nodes[n].firstPoint = static_cast<std::uint32_t>(n) * base;
      table.nodes[n].pointCount = base;
      table.nodes[n].lowerCount = static_cast<std::uint32_t>(table.lowerPoints());
    }
    table.points.assign(totalNodes * base * table.pointWidth(), 0.0);
    std::vector<NodeCurve> curves(totalNodes);

    std::atomic<std::size_t> next{0};
    std::mutex report;
    std::vector<std::thread> pool;
    const auto clock0 = std::chrono::steady_clock::now();
    auto run = [&]() {
      thermo::ReactionSource source(mechanism);
      thermo::ReactionStep step(source, 1e-10, 1e-30);
      auto solution = Cantera::newSolution(mechanism, "", "none");  // the entropy candidate's thermo
      Cantera::ThermoPhase* thermoPhase = solution->thermo().get();
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
        const auto r = buildNode(source, step, *thermoPhase, sp, candidate, maxT, sLower, table.fH[i], table.zN[j], eta, rho,
                                  table.eLo[table.fzIndex(i, j)], table.eHi[table.fzIndex(i, j)], table, idx, curves[idx]);
        if (nodeReport) {
          std::string why = r.failure;
          for (char& ch : why) if (ch == ',' || ch == '\n' || ch == '"') ch = ' ';
          std::lock_guard<std::mutex> lock(report);
          std::fprintf(nodeReport, "%d,%.6f,%.4f,%.5f,%.6e,%.3f,%.9e,%.9e,%.9e,%d,%d,%d,%d,%d,%d,%.6e,%d,%d,%s,%s,%.3e,%.3e,%.6e,%.6e,%s,%s,%s,%s,%s,%s,%zu,%d,%zu,%d,%s,%.3e\n",
                       candidate, table.fH[i], table.zN[j], eta, rho, r.tEq, r.c0, r.cEq, r.cCp, r.degenerate, r.zeroSpan,
                       r.outside, r.lowerFrozen, r.upperFrozen, r.integratorFailure, r.turnBack, r.nGateChecked, r.nGateOk,
                       r.hasDT2bLower ? std::to_string(r.dT2bLower).c_str() : "",
                       r.hasDT2bUpper ? std::to_string(r.dT2bUpper).c_str() : "", r.rawElementError, r.sumError,
                       r.spanLower, r.spanUpper, r.hasE2 ? std::to_string(r.e2QDT).c_str() : "",
                       r.hasE2 ? std::to_string(r.e2QDY).c_str() : "", r.hasE2 ? std::to_string(r.e2QDTLast).c_str() : "",
                       r.hasE2 ? std::to_string(r.e2QDYLast).c_str() : "", r.hasE2 ? std::to_string(r.e2MidDT).c_str() : "",
                       r.hasE2 ? std::to_string(r.e2MidDY).c_str() : "", r.addedPoints, r.e2Capped, r.e2RawFallbacks, r.e2Failure,
                       why.c_str(), r.fallWithinFloor);
        }
        if (r.turnBack > ts.worstTurnBack) { ts.worstTurnBack = r.turnBack; ts.worstTurnBackLoc = {i, j, k, l}; }
        if (r.nGateChecked && !r.nGateOk) { ++ts.nGateFailures; ts.worstNGateLoc = {i, j, k, l}; }
        ts.worstFallWithinFloor = std::max(ts.worstFallWithinFloor, r.fallWithinFloor);
        ts.floorDecidedNodes += r.floorDecided;
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
        if (r.hasE2) {
          if (r.e2QDT > ts.e2QDT) { ts.e2QDT = r.e2QDT; ts.e2QDTLoc = {i, j, k, l}; }
          if (r.e2QDY > ts.e2QDY) { ts.e2QDY = r.e2QDY; ts.e2QDYLoc = {i, j, k, l}; }
          ts.e2QDTLast = std::max(ts.e2QDTLast, r.e2QDTLast);
          ts.e2QDYLast = std::max(ts.e2QDYLast, r.e2QDYLast);
          ts.e2MidDT = std::max(ts.e2MidDT, r.e2MidDT);
          ts.e2MidDY = std::max(ts.e2MidDY, r.e2MidDY);
          ts.added.push_back(r.addedPoints);
          ts.e2Capped += r.e2Capped;
          ts.e2Failures += r.e2Failure;
          ts.e2RawFallbacks += r.e2RawFallbacks;
          ts.e2FallbackNodes += r.e2RawFallbacks > 0;
          ts.e2FallbackSpan = std::max(ts.e2FallbackSpan, r.e2FallbackSpan);
          const double sEq = std::fabs(r.cEq);
          for (double span : {r.spanLower, r.spanUpper})
            if (span >= 1e-9 * sEq) ts.minSpanRatio = std::min(ts.minSpanRatio, span / sEq);
        }
      }
      std::lock_guard<std::mutex> lock(report);
      if (ts.worstTurnBack > merged.worstTurnBack) { merged.worstTurnBack = ts.worstTurnBack; merged.worstTurnBackLoc = ts.worstTurnBackLoc; }
      merged.nGateFailures += ts.nGateFailures;
      merged.worstFallWithinFloor = std::max(merged.worstFallWithinFloor, ts.worstFallWithinFloor);
      merged.floorDecidedNodes += ts.floorDecidedNodes;
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
      if (ts.e2QDT > merged.e2QDT) { merged.e2QDT = ts.e2QDT; merged.e2QDTLoc = ts.e2QDTLoc; }
      if (ts.e2QDY > merged.e2QDY) { merged.e2QDY = ts.e2QDY; merged.e2QDYLoc = ts.e2QDYLoc; }
      merged.e2QDTLast = std::max(merged.e2QDTLast, ts.e2QDTLast);
      merged.e2QDYLast = std::max(merged.e2QDYLast, ts.e2QDYLast);
      merged.e2MidDT = std::max(merged.e2MidDT, ts.e2MidDT);
      merged.e2MidDY = std::max(merged.e2MidDY, ts.e2MidDY);
      merged.added.insert(merged.added.end(), ts.added.begin(), ts.added.end());
      merged.e2Capped += ts.e2Capped;
      merged.e2Failures += ts.e2Failures;
      merged.e2RawFallbacks += ts.e2RawFallbacks;
      merged.e2FallbackNodes += ts.e2FallbackNodes;
      merged.e2FallbackSpan = std::max(merged.e2FallbackSpan, ts.e2FallbackSpan);
      merged.minSpanRatio = std::min(merged.minSpanRatio, ts.minSpanRatio);
    };
    for (int w = 0; w < threads; ++w) pool.emplace_back(run);
    for (auto& t : pool) t.join();
    // The variable layout (version 3): each node's final points, in node order.
    std::size_t totalPoints = 0;
    for (const auto& cv : curves) totalPoints += cv.points.size() / table.pointWidth();
    table.points.clear();
    table.points.reserve(totalPoints * table.pointWidth());
    for (std::size_t n = 0; n < totalNodes; ++n) {
      table.nodes[n].firstPoint = static_cast<std::uint32_t>(table.points.size() / table.pointWidth());
      table.nodes[n].pointCount = static_cast<std::uint32_t>(curves[n].points.size() / table.pointWidth());
      table.nodes[n].lowerCount = curves[n].lowerCount;
      table.points.insert(table.points.end(), curves[n].points.begin(), curves[n].points.end());
    }
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();
    return std::make_pair(std::move(table), wall);
  };

  auto locStr = [&](const std::array<std::size_t, 4>& loc) {
    char buf[160];
    std::snprintf(buf, sizeof buf, "f_H %.5f Z_N %.4f eta %.5f rho %.4e", axes.fH[loc[0]], axes.zN[loc[1]], axes.eta(loc[2]),
                  std::exp(axes.lnRho(loc[3])));
    return std::string(buf);
  };

  const char* candidateNames[4] = {"Y_H2O", "Y_H2O + Y_HO2 + Y_H2O2", "the O2 consumed", "the specific entropy"};
  std::vector<int> candidates = {0, 1, 2, 3};
  if (const char* list = std::getenv("CRUCIBLE_TABLE_B_CANDIDATES")) {
    candidates.clear();
    for (double v : parseList(list)) candidates.push_back(static_cast<int>(v));
    for (int c : candidates)
      if (c < 0 || c > 3) throw std::runtime_error("CRUCIBLE_TABLE_B_CANDIDATES: no candidate " + std::to_string(c));
  }
  int chosen = -1;
  ManifoldTable table;
  ThreadStats stats;
  double buildWall = 0;
  for (int candidate : candidates) {
    ThreadStats attemptStats;
    auto [attemptTable, wall] = buildAttempt(candidate, attemptStats);
    const bool turnBackOk = attemptStats.worstTurnBack <= kTurnBackGate;
    const bool nGateOk = attemptStats.nGateFailures == 0;
    std::printf("candidate %d (%s): worst turn-back %.3e at %s; N-gate failures %zu%s; %.0f s\n", candidate,
                candidateNames[candidate], attemptStats.worstTurnBack, locStr(attemptStats.worstTurnBackLoc).c_str(),
                attemptStats.nGateFailures, attemptStats.nGateFailures ? (" (worst at " + locStr(attemptStats.worstNGateLoc) + ")").c_str() : "",
                wall);
    std::printf("  turn-back round-off floor (%.0e |c_eq|): %zu nodes with an uncounted fall above %.0e of the span; "
                "the largest uncounted fall %.3e of |c_eq|\n",
                kTurnBackFloor, attemptStats.floorDecidedNodes, kTurnBackGate, attemptStats.worstFallWithinFloor);
    std::printf("  its 2(b) worst |dT| %.3f K at %s; 2(c) frozen lower %zu, upper %zu, integrator failures %zu; "
                "criterion 1 |sum Y - 1| %.2e, elements %.2e after projection\n",
                attemptStats.worstDT2b, locStr(attemptStats.worstDT2bLoc).c_str(), attemptStats.frozenLower,
                attemptStats.frozenUpper, attemptStats.integratorFailures, attemptStats.worstSum, attemptStats.worstProjected);
    if (candidate == kEntropy) {
      auto added = attemptStats.added;
      std::sort(added.begin(), added.end());
      const std::size_t medianAdded = added.empty() ? 0 : added[added.size() / 2];
      const std::size_t maxAdded = added.empty() ? 0 : added.back();
      std::printf("  E2 refinement (linear in sqrt(s_ref - s), the 1/4, 1/2, 3/4 points to 2/3 of 3 K and dY 0.01, cap %zu points): points added "
                  "per node median %zu, max %zu; %zu nodes at the cap; final 1/4, 1/2, 3/4 points worst |dT| %.4f K, max|dY| %.3e; "
                  "%zu integration failures; %zu searches from the raw steps on %zu nodes (largest such branch span "
                  "%.3e J/(kg K))\n",
                  kMaxPointsPerNode, medianAdded, maxAdded, attemptStats.e2Capped, attemptStats.e2MidDT,
                  attemptStats.e2MidDY, attemptStats.e2Failures, attemptStats.e2RawFallbacks,
                  attemptStats.e2FallbackNodes, attemptStats.e2FallbackSpan);
      std::printf("  E2 at the eighth points (independent of the refinement; bars 3 K and dY 0.01): worst |dT| %.4f K at "
                  "%s, worst max|dY| %.3e at %s; last quarter of the span |dT| %.4f K, max|dY| %.3e: %s\n",
                  attemptStats.e2QDT, locStr(attemptStats.e2QDTLoc).c_str(), attemptStats.e2QDY,
                  locStr(attemptStats.e2QDYLoc).c_str(), attemptStats.e2QDTLast, attemptStats.e2QDYLast,
                  attemptStats.e2QDT <= kE2DT && attemptStats.e2QDY <= kE2DY && attemptStats.e2Capped == 0 &&
                          attemptStats.e2Failures == 0
                      ? "pass"
                      : "FAIL");
      std::printf("  E3 (span against 1e-9 |s_eq|): %zu nodes with a branch under it (stored as zero span); smallest "
                  "span above it %.3e |s_eq|\n", attemptStats.zeroSpanCount, attemptStats.minSpanRatio);
    }
    if (turnBackOk && nGateOk) {
      chosen = candidate;
      table = std::move(attemptTable);
      stats = attemptStats;
      buildWall = wall;
      break;
    }
    if (candidate == candidates.back()) {
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
  std::string zStr, fHStr;
  for (double v : axes.zN) zStr += (zStr.empty() ? "" : ",") + std::to_string(v);
  if (fHNodes) {
    for (double v : axes.fH) fHStr += (fHStr.empty() ? "{" : ",") + std::to_string(v);
    fHStr += "} as given";
  } else {
    fHStr = "(" + std::to_string(fLean) + " lean, " + std::to_string(fRich) + " rich, " + std::to_string(clustered) +
            " clustered)";
  }
  std::snprintf(text, sizeof text,
                "Table B built %s by crucible_build_manifold_table from %s, Cantera %s; f_H %zu nodes %s, Z_N {%s}, "
                "eta %d over [0,1], ln rho %d over [%.3g, %.3g] kg/m3; progress variable %s; "
                "%zu nodes in %.0f s on %d threads; node element error %.2e raw, %.2e after projection, |sum Y - 1| %.2e; "
                "gates: 2a %s, 2b %s, 2c %s, criterion 1 %s%s",
                date, mechanism.c_str(), Cantera::version().c_str(), nFH, fHStr.c_str(), zStr.c_str(), nEta, nLnRho,
                rhoMin, rhoMax, chosen >= 0 ? candidateNames[chosen] : candidateNames[candidates.back()], totalNodes, buildWall, threads,
                stats.worstRaw, stats.worstProjected, stats.worstSum, gate2a ? "pass" : "FAIL", gate2b ? "pass" : "FAIL",
                gate2c ? "pass" : "FAIL", criterion1 ? "pass" : "FAIL", allGatesPass ? "" : " -- GATE FAILURE");
  table.provenance = text;
  {
    auto solution = Cantera::newSolution(mechanism, "", "none");
    table.referencePressure = solution->thermo()->refPressure();
  }
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
