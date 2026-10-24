// RL10A-3-3A tier 1, the central case's dry-run build (docs/validation/rl10a33a/tier1_setup.md, its checks G1
// to G7 stated before this code). It builds the case's Definition: the PCHIP wall shared with rl10dev
// (tests/rl10_geometry.hpp), the face of choice 1 (8 rings, each a LOX slot between two post-tip walls and two
// H2 slots, the slots at the elements' own flow areas) and of choice 3 (10% of the H2 transpired at a uniform
// mass flux over the face outside the rings' cups), the supply turbulence at 5% with a length scale of 0.07
// hydraulic diameters at the slots and a pitch on the porous face (the plan, 21 October, about 02:40), and the
// LOX as Rosin-Rammler drops from Burick's curve at the rings' group over Porcheron's core length (G8). It
// reports the layout, the flows and the injection numbers, maps each supply onto candidate meshes by Flow's
// ring-centre rule, constructs Flow and Spray on each, and exits before any step. No measured value is read.
//
// Usage: crucible_rl10_tier1 [NZxNR ...]   (default 76x8 149x16 1192x128)
//        crucible_rl10_tier1 run development|central   (tier1_run_definition.md, *The run mode*)
//        crucible_rl10_tier1 selftest                  (the settle judge and the run length's moving window, R4)
// The run takes CRUCIBLE_MANIFOLD_TABLE (Table B), and in development CRUCIBLE_T1_BASE (149x16), _LEVELS (3),
// _LEAD [s] (t_L, 20e-3), _END [s] (the simulated time, t_L + 1e-3), _THREADS (4), _OUT (the files' prefix, rl10_t1) and
// _RATE_FACTOR (1), _CFL (Definition::cfl, 0.4), _SEED (the spray's, 1). Restart (tests/checkpoint_file.hpp, docs/evidence/RESTART.md): CRUCIBLE_CHECKPOINT=<stem> writes
// checkpoints and CRUCIBLE_RESTART=auto continues from the latest, bit for bit the uninterrupted run; _THREADS may
// change and _END be extended across a restart. On any exit the run writes <out>_done ("exit N"), <out> being
// rl10_t1_central or <_OUT>_DEVELOPMENT; a killed process writes none.
// A member's element dimensions replace the central ones through the environment (the skeleton's brackets):
//   CRUCIBLE_T1_POST_ID      d_L [mm], central 3.54
//   CRUCIBLE_T1_TIP          t / d_L, central 0.272
//   CRUCIBLE_T1_GAP          h / d_L, central 0.268
//   CRUCIBLE_T1_RECESS       recess / d_L, central 0.5
//   CRUCIBLE_T1_OUTER_RING   (R_c - r_8) / p, central 0.5 (bracket 0.5 to 1)
//   CRUCIBLE_T1_INNER_RING   pitch (r_1 = r_8 - 7p, central) or clearance (r_1 = 27 D_a / 2 pi)
//   CRUCIBLE_T1_TRANSPIRED   transpired share of the fuel, central 0.10
//   CRUCIBLE_T1_SLOT_LENGTH  the supply turbulence's length scale at the slots / their hydraulic diameter, central 0.07
//   CRUCIBLE_T1_FACE_PRESSURE  G4's nominal face pressure [Pa], 2.8e6 (our development runs', not a measurement)
//   CRUCIBLE_RADIAL_STRETCHING  Definition::radialStretching, 0
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "cantera/core.h"
#include "core/amr.hpp"
#include "core/flow.hpp"
#include "core/manifold_runtime.hpp"
#include "core/manifold_table.hpp"
#include "core/spray.hpp"
#include "tests/amr_dump.hpp"
#include "tests/checkpoint_file.hpp"
#include "tests/rl10_geometry.hpp"

using namespace crucible;

namespace {

constexpr double kPi = std::numbers::pi;
const std::string kSaturation = CRUCIBLE_SOURCE_DIR "/tools/o2_saturation_nist.tsv";
// Protocol section 2 (P1) and its addendum's plenum temperatures; the wall at the bracket's central 567 K.
constexpr double kOxidizer = 13.9525, kMixtureRatio = 5.055, kFuelT = 199.67, kOxidizerT = 101.17, kWallT = 567;
constexpr int kRings = 8, kElements = 216, kPerRing = 27;
constexpr double kAmbientP = 1000, kAmbientT = 300, kRamp = 5e-4;
constexpr double kIntensity = 0.05, kCmu = 0.09;

double env(const char* name, double fallback) {
  const char* v = std::getenv(name);
  return v ? std::atof(v) : fallback;
}

int failures = 0, gateFailures = 0;
void check(const char* name, bool ok, const std::string& detail) {
  if (!ok) ++failures;
  // The run mode's gates: G5 and G6 test the ring-centre rule decision 6 superseded, and G7 is declared.
  if (!ok && std::string(name) != "G5" && std::string(name) != "G6" && std::string(name) != "G7") ++gateFailures;
  std::printf("%-4s %s  %s\n", name, ok ? "pass" : "FAIL", detail.c_str());
}
std::string strf(const char* f, double a, double b = 0, double c = 0, double d = 0) {
  char s[256];
  std::snprintf(s, sizeof s, f, a, b, c, d);
  return s;
}
double rel(double a, double b) { return std::abs(a - b) / std::abs(b); }

// Burick (1972) Fig. 38, the digitised curve: (V_g - V_L) / (V_L MR) against D-bar / D_L, the mass median over the
// post diameter (docs/evidence/spray/burick_fig38.md, section 1.4). Wax/GN2 cold flow, used as printed.
constexpr std::array<std::array<double, 2>, 32> kBurick{{
    {0.23, 0.277}, {0.28, 0.276}, {0.35, 0.272}, {0.45, 0.265}, {0.55, 0.260}, {0.65, 0.252}, {0.75, 0.247},
    {0.85, 0.242}, {1.00, 0.241}, {1.30, 0.238}, {1.60, 0.248}, {2.00, 0.243}, {2.30, 0.225}, {2.60, 0.209},
    {2.80, 0.206}, {3.20, 0.169}, {3.60, 0.165}, {4.00, 0.159}, {4.50, 0.154}, {5.00, 0.148}, {5.50, 0.154},
    {6.00, 0.129}, {6.50, 0.133}, {7.00, 0.129}, {7.70, 0.125}, {8.50, 0.118}, {9.50, 0.113}, {10.50, 0.112},
    {11.50, 0.113}, {12.50, 0.101}, {13.50, 0.097}, {14.50, 0.094}}};
// The scatter's stated width about the curve (its section 1.5, "+30%/-40%"): member 5's two ends.
constexpr double kScatterHigh = 1.3, kScatterLow = 1 / 1.7;

// Beyond the table: a power law through its last point at the slope fitted over group 4 to 12
// (MASCOTTE_A10_PREDICTION.md section 3; tier1_setup.md, *Two member rules*).
constexpr double kBurickSlope = -0.371;

// D-bar / D_L at group x, linear in log-log between neighbours, the power law past the end; NaN below the start.
double burick(double x) {
  if (x > kBurick.back()[0]) return kBurick.back()[1] * std::pow(x / kBurick.back()[0], kBurickSlope);
  for (std::size_t i = 0; i + 1 < kBurick.size(); ++i) {
    const auto& [x0, y0] = kBurick[i];
    const auto& [x1, y1] = kBurick[i + 1];
    if (x >= x0 && x <= x1)
      return std::exp(std::log(y0) + (std::log(x / x0) / std::log(x1 / x0)) * std::log(y1 / y0));
  }
  return std::nan("");
}

// A mass-based Rosin-Rammler's D30 / D-bar; Burick's Eq. (5) sets it to 0.455.
double d30OverMedian(double q) { return std::pow(std::tgamma(1 - 3 / q), -1.0 / 3) / std::pow(std::log(2.0), 1 / q); }

// The settle judge (tier1_run_definition.md, *Settled*: the Director's decision of 22 October 2026, about 22:40
// CDT) on a series sampled every dt over a window of `span`: the batch-means standard error and the halves'
// difference, each over the window's mean. tau is the autocorrelation integrated to its first zero; the batch
// length is the shortest of 1, 2, 4 and 5 ms at least 10 tau, and the judged SE the larger of its and the next
// longer length's.
struct Settle {
  double mean{}, se{}, halves{}, tau{}, batch{};
  bool shortBatches{false}, settled{false};
};
double batchSe(const std::vector<double>& x, std::size_t batches) {
  const std::size_t n = x.size() / batches;
  std::vector<double> m(batches, 0.0);
  for (std::size_t b = 0; b < batches; ++b) {
    for (std::size_t k = 0; k < n; ++k) m[b] += x[b * n + k];
    m[b] /= static_cast<double>(n);
  }
  double mean = 0, s = 0;
  for (double v : m) mean += v;
  mean /= static_cast<double>(batches);
  for (double v : m) s += (v - mean) * (v - mean);
  return std::sqrt(s / (static_cast<double>(batches) * static_cast<double>(batches - 1)));
}
Settle settle(const std::vector<double>& x, double dt, double span) {
  Settle s;
  const std::size_t n = x.size(), half = n / 2;
  double a = 0, b = 0, var = 0;
  for (std::size_t k = 0; k < n; ++k) (k < half ? a : b) += x[k];
  s.mean = (a + b) / static_cast<double>(n);
  s.halves = std::abs(b / static_cast<double>(n - half) - a / static_cast<double>(half)) / std::abs(s.mean);
  for (double v : x) var += (v - s.mean) * (v - s.mean);
  if (var > 0) {
    double previous = 1;
    for (std::size_t k = 1; k < n; ++k) {
      double c = 0;
      for (std::size_t i = 0; i + k < n; ++i) c += (x[i] - s.mean) * (x[i + k] - s.mean);
      const double r = c / var;
      if (r <= 0) {  // to the zero between k - 1 and k, linearly
        s.tau += 0.5 * dt * previous * previous / (previous - r);
        break;
      }
      s.tau += 0.5 * dt * (previous + r);
      previous = r;
    }
  }
  constexpr std::array<double, 4> lengths{1e-3, 2e-3, 4e-3, 5e-3};
  std::size_t pick = lengths.size() - 1;
  s.shortBatches = true;
  for (std::size_t i = 0; i < lengths.size(); ++i)
    if (lengths[i] >= 10 * s.tau) { pick = i, s.shortBatches = false; break; }
  s.batch = lengths[pick];
  auto se = [&](std::size_t i) { return batchSe(x, static_cast<std::size_t>(std::lround(span / lengths[i]))); };
  s.se = se(pick);
  if (pick + 1 < lengths.size()) s.se = std::max(s.se, se(pick + 1));
  s.se /= std::abs(s.mean);
  s.settled = s.se <= 0.002 && s.halves <= 0.005;
  return s;
}

// The run length (tier1_run_definition.md, *Settled*: the Director's decision C of 24 October 2026, 01:40 CDT). In
// samples of 10 us (kSample), counted from t_L: from 30 ms on, every 1 ms, the trailing 20 ms window of the face
// pressure and of c* is judged; the first check at which both settle is t_s, and 60 ms is the ceiling (judged there,
// reported flagged if unsettled). The continuation at twice the parcel rate runs to t_s + 25 ms, judged over its last
// 20 ms.
constexpr long kWindow = 2000, kFirstCheck = 3000, kCheckEvery = 100, kCeiling = 6000, kContinuation = 2500;
enum class SettleCheck { None, NotSettled, Settled, Ceiling };
struct RunCheck {
  SettleCheck check{SettleCheck::None};
  Settle face, cstar;
};
// At sample k (at k kSample; t_L at leadK; face[n - 1] and cstar[n - 1] are sample n): None if k is not a check,
// otherwise the window of samples (k - kWindow, k] judged.
RunCheck runCheck(long k, long leadK, const std::vector<double>& face, const std::vector<double>& cstar) {
  RunCheck r;
  const long since = k - leadK;
  if (since < kFirstCheck || since > kCeiling || since % kCheckEvery != 0) return r;
  auto judge = [&](const std::vector<double>& x) {
    return settle(std::vector<double>(x.begin() + (k - kWindow), x.begin() + k), 1e-5, 20e-3);
  };
  r.face = judge(face), r.cstar = judge(cstar);
  r.check = r.face.settled && r.cstar.settled ? SettleCheck::Settled : since == kCeiling ? SettleCheck::Ceiling : SettleCheck::NotSettled;
  return r;
}

// R4 (tier1_run_definition.md, *The run mode*): the judge on synthetic 10 us series over 40 to 60 ms, (a) to (c), and
// the run length's moving window (decision C), (d) to (g).
int selftest() {
  constexpr double dt = 1e-5, span = 20e-3, period = 2e-3;
  std::vector<double> constant, step, sine;
  for (int k = 1; k <= 2000; ++k) {
    const double t = 40e-3 + k * dt;
    constant.push_back(3.0);
    step.push_back(k <= 1000 ? 3.0 : 3.03);
    sine.push_back(1 + 1e-3 * std::sin(2 * kPi * t / period));
  }
  const Settle c = settle(constant, dt, span), s = settle(step, dt, span), w = settle(sine, dt, span);
  const double expected = 0.01 / 1.005, tau = period / (2 * kPi);
  check("R4", c.settled && c.se == 0 && c.halves == 0, strf("(a) constant: settled %.0f, SE %.3g, halves %.3g", c.settled, c.se, c.halves));
  check("R4", !s.settled && std::abs(s.halves - expected) <= 1e-9 && s.shortBatches,
        strf("(b) +1%% step at 50 ms: settled %.0f, halves %.12f (expected %.12f), batches shorter than 10 tau %.0f", s.settled,
            s.halves, expected, s.shortBatches) + strf(" (tau %.4f ms, SE %.4g)", s.tau * 1e3, s.se));
  check("R4", std::abs(w.tau / tau - 1) <= 0.02 && w.batch == 4e-3 && w.settled,
        strf("(c) sine of period 2 ms: tau %.5f ms against P / 2 pi %.5f ms (%+.3f%%), batch %.0f ms", w.tau * 1e3, tau * 1e3,
            100 * (w.tau / tau - 1), w.batch * 1e3) + strf(", settled %.0f (SE %.3g, halves %.3g)", w.settled, w.se, w.halves));
  // The run length (decision C) on series of a 20 ms lead to t_L + 60 ms: the first check that settles, or the ceiling.
  constexpr long leadK = 2000;
  auto series = [&](double base, double from, double to, double slope) {  // +10% over (t_L + from, t_L + to] ms; a drift per ms
    std::vector<double> x;
    for (long n = 1; n <= leadK + kCeiling; ++n) {
      const long m = n - leadK;
      const bool in = m > std::lround(from * 100) && m <= std::lround(to * 100);
      x.push_back(base * (1 + (in ? 0.1 : 0) + slope * static_cast<double>(m) * 1e-2));
    }
    return x;
  };
  auto first = [&](const std::vector<double>& face, const std::vector<double>& cstar, long& at) {
    for (long k = 1; k <= leadK + kCeiling; ++k) {
      const RunCheck r = runCheck(k, leadK, face, cstar);
      if (r.check == SettleCheck::Settled || r.check == SettleCheck::Ceiling) {
        at = k - leadK;
        return r.check;
      }
    }
    return SettleCheck::None;
  };
  const auto flat = series(2.8e6, 0, 0, 0), cFlat = series(1700, 0, 0, 0);
  const auto cBlock = series(1700, 15, 17, 0), faceBlock = series(2.8e6, 20, 22, 0), cDrift = series(1700, 0, 0, 1e-3);
  long d = 0, e = 0, f = 0, g = 0;
  const SettleCheck dc = first(flat, cFlat, d), ec = first(flat, cBlock, e), fc = first(faceBlock, cBlock, f), gc = first(flat, cDrift, g);
  check("R4", dc == SettleCheck::Settled && d == 3000, strf("(d) constant: t_s at t_L + %.2f ms (expected 30, the first check)", d * 1e-2));
  check("R4", ec == SettleCheck::Settled && e == 3700,
        strf("(e) c* +10%% over t_L + 15 to 17 ms: t_s at t_L + %.2f ms (expected 37, the first window clear of it)", e * 1e-2));
  check("R4", fc == SettleCheck::Settled && f == 4200,
        strf("(f) and the face pressure +10%% over t_L + 20 to 22 ms: t_s at t_L + %.2f ms (expected 42: both must settle)", f * 1e-2));
  const RunCheck between = runCheck(leadK + 3050, leadK, flat, cFlat), beyond = runCheck(leadK + 6100, leadK, flat, cFlat);
  check("R4", gc == SettleCheck::Ceiling && g == kCeiling && between.check == SettleCheck::None && beyond.check == SettleCheck::None,
        strf("(g) c* drifting 0.1%% per ms: the ceiling at t_L + %.2f ms, unsettled %.0f (expected 60, 1); no check between "
             "checks or after the ceiling %.0f", g * 1e-2, gc == SettleCheck::Ceiling,
             between.check == SettleCheck::None && beyond.check == SettleCheck::None));
  std::printf("%d check(s) failed\n", failures);
  return failures ? 1 : 0;
}

// Equilibrium isentropic stagnation pressure of a cell, and its column's mass-flux average (chamber_study's, which
// rl10dev used); NaN while the column's net forward mass flux is not positive.
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
double columnStagnation(const Flow& flow, Cantera::ThermoPhase& gas, int column) {
  const auto& mesh = flow.mesh();
  double p0 = 0, flux = 0;
  for (int j = 0; j < mesh.nr; ++j) {
    auto q = mesh.index(column, j);
    auto w = flow.cellPrimitive(q);
    auto y = flow.massFractions(q);
    double g = w.rho * w.uz * mesh.axialArea(column + 1, j);
    if (g == 0) continue;
    p0 += g * stagnationPressure(gas, flow.temperature(q), w.p, y.data(), 0.5 * (w.uz * w.uz + w.ur * w.ur));
    flux += g;
  }
  return flux > 0 ? p0 / flux : std::nan("");
}

// The run (tier1_run_definition.md, *The run mode*).
constexpr bool kNumericsFrozen = false;  // the y+ target and the tag thresholds; set by a dated commit
constexpr double kG0 = 9.80665, kEngineFuel = 2.7946, kParcelRate = 5e5, kSample = 1e-5;
struct Setup {
  std::string mode, out, mechanism;
  int argc{};
  char** argv{};  // the run's signature (tests/checkpoint_file.hpp)
  int levels{3}, threads{4};
  double end{1e-3}, rateFactor{1}, core{}, cylinderEnd{}, chamberFlow{};
  double lead{20e-3};  // t_L, the H2 lead (tier1_run_definition.md, the start amendment)
  Definition definition;
  Spray::Settings spray;
  std::vector<double> o2;
};
std::string doneFile;

int runCase(const Setup& s) {
  const bool development = s.mode == "development";
  const std::string label = development ? "DEVELOPMENT (protocol §6): not the reported case" : "RL10A-3-3A tier 1, the central case";
  const char* tablePath = std::getenv("CRUCIBLE_MANIFOLD_TABLE");
  if (!tablePath) throw std::runtime_error("the run needs CRUCIBLE_MANIFOLD_TABLE (Table B)");
  const ManifoldTable manifold = ManifoldTable::load(tablePath);
  if (manifold.provenance.find(s.mechanism.substr(s.mechanism.find_last_of('/') + 1)) == std::string::npos)
    throw std::runtime_error(std::string("the table ") + tablePath + " was not built from " + s.mechanism);
  std::printf("\n%s\nTable B %s: %s\n", label.c_str(), tablePath, manifold.provenance.c_str());

  Amr amr(s.definition, s.levels);
  amr.setThreads(s.threads);
  auto solution = Cantera::newSolution(s.mechanism, "", "none");
  Cantera::ThermoPhase& gas = *solution->thermo();
  Amr::Tags tags;  // the defaults: the thresholds frozen at s* = 1/4 for tier 1 (AMR.md, Criterion 4: results, da3555f)
  const std::size_t eH = gas.elementIndex("H");
  for (const auto& sp : s.definition.species) {
    const std::size_t k = gas.speciesIndex(sp.name);
    tags.hydrogen.push_back(gas.nAtoms(k, eH) * gas.atomicWeight(eH) / gas.molecularWeight(k));
  }
  const double core = s.core;
  tags.custom = [&amr, core](int l, int i, int j) {  // item 4's static face tag
    const Mesh& m = amr.level(l).mesh();
    return m.cells[m.index(i, j)].z < core;
  };
  amr.setTags(tags);
  std::vector<std::unique_ptr<thermo::ReactingFlow>> reacting;
  for (int l = 0; l < amr.levels(); ++l) {
    reacting.push_back(std::make_unique<thermo::ReactingFlow>(amr.level(l), s.mechanism, s.threads, 1e-6, 1e-12,
                                                              thermo::Chemistry::FiniteRate));
    reacting.back()->useManifoldTable(manifold);
    // CRUCIBLE_CHEM_SKIP=0 integrates every cell, for timing the skip against the unskipped path.
    reacting.back()->setChemistrySkip(env("CRUCIBLE_CHEM_SKIP", 1) != 0);
  }
  Spray::Settings settings = s.spray;
  for (auto& j : settings.injectors) j.parcelRate = kParcelRate * s.rateFactor;
  Spray spray(amr, LiquidOxygen::load(kSaturation, amr.level(0).medium().enthalpy(298.15, s.o2.data())), settings);
  amr.setAfterStep([&](int l, double dt) {
    spray.afterStep(l, dt);
    reacting[static_cast<std::size_t>(l)]->react(dt);
  });
  amr.setRedo([&] { spray.redo(); });

  // The base's stations: the throat (A_t, the regions), p0's column, the face's area.
  const Flow& base = amr.level(0);
  const Mesh& mesh0 = base.mesh();
  int throat = 0, endColumn = 0;
  for (int i = 0; i <= mesh0.nz; ++i)
    if (mesh0.radius[i] < mesh0.radius[throat]) throat = i;
  for (int i = 0; i < mesh0.nz; ++i)
    if (mesh0.cells[mesh0.index(i, 0)].z < s.cylinderEnd) endColumn = i;
  const double at = kPi * mesh0.radius[throat] * mesh0.radius[throat], zThroat = throat * mesh0.dz;
  double faceArea = 0;
  for (int j = 0; j < mesh0.nr; ++j) faceArea += mesh0.axialArea(0, j);
  std::printf("%s: base %d x %d plus %d levels (%d x %d finest), %d threads, to %.4g ms, CFL %.3g; parcels %.3g per s per injector; "
              "throat station %d (A_t %.6e m2), p0 column %d at z %.4f mm, face area %.6e m2\n",
              s.mode.c_str(), mesh0.nz, mesh0.nr, s.levels, mesh0.nz << s.levels, mesh0.nr << s.levels, s.threads, s.end * 1e3,
              s.definition.cfl, kParcelRate * s.rateFactor, throat, at, endColumn, mesh0.cells[mesh0.index(endColumn, 0)].z * 1e3, faceArea);
  std::printf("tags (thresholds frozen at s* = 1/4, the y+ target %s): interval %d, buffer %d, throat %.2f r_t to level %d, "
              "wall y+ %.0f, f_H %.3g, T %.3g, p %.3g, vorticity %.3g; the face column to the top within %.3f mm\n",
              kNumericsFrozen ? "frozen" : "provisional (step D)",
              tags.interval, tags.buffer, tags.throatHalfWidth, tags.throatLevel, tags.wallYPlus, tags.mixture,
              tags.temperature, tags.pressure, tags.vorticity, core * 1e3);
  // The start amendment: H2 alone from t = 0, the LOX from t_L, the judge's times counted from t_L.
  const long leadK = std::lround(s.lead / kSample);
  std::printf("the start: the H2 supplies open at 0, the LOX at t_L = %.4g ms on its %.2g ms ramp, the igniter at %.4g ms; "
              "the lead's gate %s; the exchange bound (registration B) a step limit, f 0.25\n", s.lead * 1e3, kRamp * 1e3,
              s.definition.igniter.start * 1e3, leadK >= 200 ? "judged at t_L" : "not judged (a lead under 2 ms)");

  // The composite owned cells.
  auto owned = [&](int l, std::size_t q) {
    return (l == 0 || amr.active(l)[q]) && !(l + 1 < amr.levels() && amr.covered(l, q));
  };
  struct Sample {
    double t{}, mdot{}, fvac{}, ispChamber{}, ispInlet{}, p0{}, pFace{}, cstar{}, qWall{}, qFace{}, liquid{}, parcels{};
    double areaError{};
  };
  auto sample = [&]() {
    Sample x;
    x.t = amr.time();
    double wallHeat = 0, plateHeat = 0, pressureArea = 0, area = 0;
    for (int l = 0; l < amr.levels(); ++l) {
      const Flow& f = amr.level(l);
      const auto m = f.measurements();
      x.mdot += m.outletMassFlow;
      x.fvac += m.inletMomentumFlux + m.wallAxialForce + m.bodyAxialForce;
      wallHeat += m.wallHeatFlow;
      plateHeat += m.plateHeatFlow;
      for (int j = 0; j < f.mesh().nr; ++j) {
        const std::size_t q = f.mesh().index(0, j);
        if (!owned(l, q)) continue;
        area += f.mesh().axialArea(0, j);
        pressureArea += f.mesh().axialArea(0, j) * f.cellPrimitive(q).p;
      }
    }
    x.fvac += spray.ledger().injectionThrust + spray.ledger().wallThrust;
    x.ispChamber = x.fvac / (kG0 * s.chamberFlow);
    x.ispInlet = x.fvac / (kG0 * (kEngineFuel + kOxidizer));
    x.p0 = columnStagnation(base, gas, endColumn);
    x.pFace = pressureArea / area;
    x.areaError = rel(area, faceArea);
    x.cstar = x.mdot > 0 ? x.p0 * at / x.mdot : std::nan("");
    x.qWall = -(wallHeat - plateHeat);
    x.qFace = -plateHeat;
    x.liquid = spray.liquidMass();
    x.parcels = static_cast<double>(spray.parcels().size());
    return x;
  };

  // Budgets (a) to (c) as A2's (SPRAY.md, Spray on AMR), over the gas's initial content.
  double mass0 = 0, energy0 = 0;
  for (std::size_t q = 0; q < base.state().size(); ++q) {
    mass0 += base.state()[q][0] * mesh0.cells[q].volume;
    energy0 += base.state()[q][3] * mesh0.cells[q].volume;
  }
  auto budgets = [&]() {
    const auto& l = spray.ledger();
    const auto b = amr.budget();
    double exchangedMass = 0, exchangedEnergy = 0, exchangedMomentum = 0;
    for (int k = 0; k < amr.levels(); ++k) {
      const auto& m = amr.level(k).measurements();
      exchangedMass += m.exchangedMass;
      exchangedEnergy += m.exchangedEnergy;
      exchangedMomentum += m.exchangedMomentum;
    }
    const double massScale = mass0 + l.injectedMass;
    const double energyScale = std::abs(energy0) + std::abs(l.injectedEnergy) + std::abs(l.outflowEnergy) + l.wallLoss;
    const double momentumScale = b.momentumGross + std::abs(l.injectedMomentum) + std::abs(l.outflowMomentum) + std::abs(l.wallImpulse);
    const auto e = spray.balanceErrors();
    const double systemMass = b.massResidual + exchangedMass + spray.liquidMass() - (l.injectedMass - l.outflowMass);
    const double systemEnergy = b.energyResidual + exchangedEnergy + spray.liquidEnergy() - (l.injectedEnergy - l.outflowEnergy - l.wallLoss);
    const double systemMomentum =
        b.momentumResidual + exchangedMomentum + spray.liquidMomentum() - (l.injectedMomentum - l.outflowMomentum + l.wallImpulse);
    return std::array<double, 3>{
        std::max({std::abs(b.massResidual) / massScale, std::abs(b.energyResidual) / energyScale, std::abs(b.momentumResidual) / momentumScale}),
        std::max({std::abs(e[0]), std::abs(e[1]), std::abs(e[2])}),
        std::max({std::abs(systemMass) / massScale, std::abs(systemEnergy) / energyScale, std::abs(systemMomentum) / momentumScale})};
  };

  std::FILE* csv = nullptr;
  std::vector<Sample> series;
  std::array<double, 3> worstBudget{};
  double worstArea = 0, worstGrid = 0;
  long nonfinite = 0, undefined = 0;
  const auto clock0 = std::chrono::steady_clock::now();
  auto print = [&](const Sample& x) {
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();
    // The budgets have no reference before the first step.
    const bool stepped = amr.steps() > 0;
    const auto b = stepped ? budgets() : std::array<double, 3>{};
    if (stepped)
      for (int k = 0; k < 3; ++k) worstBudget[k] = std::max(worstBudget[k], b[k]);
    const double quantities[] = {x.mdot, x.fvac, x.ispChamber, x.ispInlet, x.pFace, x.qWall, x.qFace, x.liquid, b[0], b[1], b[2]};
    for (double v : quantities) nonfinite += !std::isfinite(v);
    double tMax = 0;
    for (int l = 0; l < amr.levels(); ++l) {
      const Flow& f = amr.level(l);
      for (std::size_t q = 0; q < f.mesh().cells.size(); ++q)
        if (f.mesh().cells[q].z < core && owned(l, q)) tMax = std::max(tMax, f.temperature(q));
    }
    // Parcels within the core length, by the finest level's cell they lie in.
    const Mesh& fine = amr.level(amr.levels() - 1).mesh();
    std::vector<unsigned char> occupied(fine.cells.size(), 0);
    long near = 0, cells = 0;
    for (const Parcel& p : spray.parcels()) {
      if (!(p.z < core)) continue;
      const int i = std::min(fine.nz - 1, static_cast<int>(p.z / fine.dz));
      const double eta = p.r / (fine.radius[i] + (fine.radius[i + 1] - fine.radius[i]) * (p.z / fine.dz - i));
      const int j = std::clamp(static_cast<int>(std::upper_bound(fine.fraction.begin(), fine.fraction.end(), eta) - fine.fraction.begin()) - 1, 0, fine.nr - 1);
      ++near;
      auto& o = occupied[fine.index(i, j)];
      cells += !o;
      o = 1;
    }
    std::printf("\n--- t %.3f ms  [%s]\n  %llu base steps (%llu redone), %llu regrids; wall %.0f s, the spray's share %.1f%%; active cells",
                x.t * 1e3, label.c_str(), static_cast<unsigned long long>(amr.steps()),
                static_cast<unsigned long long>(amr.rejectedSteps()), static_cast<unsigned long long>(amr.regrids()), wall,
                wall > 0 ? 100 * spray.wallTime() / wall : 0.0);
    for (int l = 1; l < amr.levels(); ++l) {
      long n = 0;
      for (const auto a : amr.active(l)) n += a;
      std::printf(" %d: %ld", l, n);
    }
    std::printf("\n  F_vac %.6g N, Isp_vac %.6g s (chamber flow) and %.6g s (engine-inlet flow); outlet %.6g kg/s\n", x.fvac,
                x.ispChamber, x.ispInlet, x.mdot);
    std::printf("  p0 %.6g Pa, face pressure %.6g Pa, c* %.6g m/s; hot-gas wall %.6g W, face %.6g W (gas to wall)\n", x.p0, x.pFace,
                x.cstar, x.qWall, x.qFace);
    std::printf("  liquid %.6g kg in %.0f parcels; within %.2f mm of the face: T_max %.1f K, %ld parcels in %ld finest cells (%.2f per cell)\n",
                x.liquid, x.parcels, core * 1e3, tMax, near, cells, cells ? static_cast<double>(near) / cells : 0.0);
    if (stepped)
      std::printf("  budgets: (a) gas %.2e, (b) liquid %.2e, (c) whole system %.2e (limit 1e-11)\n", b[0], b[1], b[2]);
    else
      std::printf("  budgets: none before the first step\n");
    std::printf("  Table B, table / all updates: ");
    for (int l = 0; l < amr.levels(); ++l) {
      const auto& st = reacting[static_cast<std::size_t>(l)]->stats();
      const Mesh& m = amr.level(l).mesh();
      constexpr auto causes = static_cast<std::size_t>(ManifoldRuntime::Cause::Count);
      std::array<long, 3> table{}, all{};
      for (std::size_t q = 0; q < m.cells.size() && (q + 1) * causes <= st.manifoldCellCauses.size(); ++q) {
        const double z = m.cells[q].z;
        const int region = z < s.cylinderEnd ? 0 : z < zThroat ? 1 : 2;
        for (std::size_t c = 0; c < causes; ++c) all[region] += st.manifoldCellCauses[q * causes + c];
        table[region] += st.manifoldCellCauses[q * causes];
      }
      std::printf("%slevel %d: chamber %ld / %ld, convergent %ld / %ld, divergent %ld / %ld", l ? "; " : "", l, table[0], all[0],
                  table[1], all[1], table[2], all[2]);
    }
    std::printf("\n  chemistry skipped / integrated (ReactingFlow::setChemistrySkip): ");
    for (int l = 0; l < amr.levels(); ++l) {
      const auto& st = reacting[static_cast<std::size_t>(l)]->stats();
      std::printf("%slevel %d: %ld / %ld", l ? "; " : "", l, st.skippedCells, st.integratedCells);
    }
    std::printf("\n");
    std::fflush(stdout);
  };
  auto dump = [&](const std::string& tag) {
    writeAmrDump(amr, s.out + "_" + tag + "_amr.bin");
    std::FILE* f = std::fopen((s.out + "_" + tag + "_spray.bin").c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write the spray dump");
    for (int l = 1; l < amr.levels(); ++l) std::fwrite(amr.active(l).data(), 1, amr.active(l).size(), f);
    const double t = amr.time();
    const std::uint64_t n = spray.parcels().size();
    std::fwrite(&t, sizeof t, 1, f);
    std::fwrite(&n, sizeof n, 1, f);
    // Each parcel with its trailing padding zeroed, so equal states give equal files.
    for (const Parcel& p : spray.parcels()) {
      unsigned char bytes[sizeof(Parcel)];
      std::memcpy(bytes, &p, sizeof(Parcel));
      constexpr std::size_t used = offsetof(Parcel, phase) + sizeof(int);
      std::memset(bytes + used, 0, sizeof(Parcel) - used);
      std::fwrite(bytes, 1, sizeof(Parcel), f);
    }
    std::fwrite(&spray.ledger(), sizeof(Spray::Ledger), 1, f);
    std::fclose(f);
    std::printf("  state written: %s_%s_amr.bin and _spray.bin (for reading out; a restart is CRUCIBLE_CHECKPOINT's file)\n",
                s.out.c_str(), tag.c_str());
  };
  auto window = [&](long last, const char* name) {  // the 20 ms ending at sample `last`
    std::printf("\nThe settle judge over %.0f to %.0f ms, t_L + %.0f to t_L + %.0f (%s) [%s]\n", (last - 2000) * kSample * 1e3,
                last * kSample * 1e3, (last - leadK - 2000) * kSample * 1e3, (last - leadK) * kSample * 1e3, name, label.c_str());
    std::array<Settle, 8> out{};
    const char* names[] = {"face pressure", "c*", "F_vac", "outlet flow", "p0", "Isp_vac (chamber)", "hot-gas wall heat", "face heat"};
    for (std::size_t k = 0; k < out.size(); ++k) {
      std::vector<double> x;
      for (long n = last - 2000; n < last; ++n) {
        const Sample& v = series[static_cast<std::size_t>(n)];
        const double values[] = {v.pFace, v.cstar, v.fvac, v.mdot, v.p0, v.ispChamber, v.qWall, v.qFace};
        x.push_back(values[k]);
      }
      out[k] = settle(x, kSample, 20e-3);
      std::printf("  %-18s mean %.7g, SE %.3f%%, halves %.3f%%, tau %.4f ms, batch %.0f ms%s%s\n", names[k], out[k].mean,
                  100 * out[k].se, 100 * out[k].halves, out[k].tau * 1e3, out[k].batch * 1e3,
                  out[k].shortBatches ? " (batches shorter than 10 tau)" : "", k < 2 ? (out[k].settled ? ": settled" : ": NOT settled") : "");
    }
    std::printf("  the run is %s\n", out[0].settled && out[1].settled ? "settled" : "NOT settled (reported with these values, flagged)");
    return out;
  };

  std::array<Settle, 8> first{};
  long k = 0, shortened = 0, flashesAt = 0, settledK = 0;  // settledK: t_s in samples (decision C), 0 before
  bool flashReported = false;
  double smallestBound = std::numeric_limits<double>::infinity();
  // Restart (tests/checkpoint_file.hpp, docs/evidence/RESTART.md): the engine's objects and this loop's state; the
  // series file is cut back to its length at the checkpoint and appended.
  CheckpointFile checkpoint(s.argc, s.argv);
  std::uint64_t csvSize = 0;
  auto fields = [&](auto& io) {
    io.object(amr);
    for (auto& r : reacting) io.object(*r);
    io.object(spray);
    io.mark("tier1");
    io.same(mass0, "the initial gas mass");
    io.same(energy0, "the initial gas energy");
    io.field(series);
    io.field(worstBudget), io.field(worstArea), io.field(worstGrid), io.field(nonfinite), io.field(undefined);
    for (auto& f : first)
      io.field(f.mean), io.field(f.se), io.field(f.halves), io.field(f.tau), io.field(f.batch), io.field(f.shortBatches),
          io.field(f.settled);
    io.field(k), io.field(shortened), io.field(flashesAt), io.field(flashReported), io.field(smallestBound), io.field(settledK);
    io.field(csvSize);
  };
  if (checkpoint.read([&](CheckpointReader& in) { fields(in); })) {
    csv = reopenAt(s.out + "_series.csv", csvSize);
    std::printf("restart: t %.6g ms, %ld samples, %llu base steps, %zu parcels\n", amr.time() * 1e3, k,
                static_cast<unsigned long long>(amr.steps()), spray.parcels().size());
  } else {
    csv = std::fopen((s.out + "_series.csv").c_str(), "w");
    if (!csv) throw std::runtime_error("cannot write " + s.out + "_series.csv");
    std::fprintf(csv, "# %s\nt_s,mdot_kg_s,fvac_N,isp_chamber_s,isp_inlet_s,p0_Pa,pface_Pa,cstar_m_s,qwall_W,qface_W,liquid_kg,parcels\n",
                 label.c_str());
    print(sample());
  }
  auto save = [&] {
    csvSize = flushedSize(csv);
    checkpoint.write([&](CheckpointWriter& out) { fields(out); });
  };
  auto stateOut = [&] {
    writeStateOut([&](CheckpointWriter& out) {
      out.object(amr);
      for (const auto& r : reacting) out.object(*r);
      out.object(spray);
    });
  };
  // The run ends at s.end, or at t_s + 25 ms once t_s is known (decision C; the central case's s.end is the ceiling's).
  auto runEnd = [&] { return settledK ? std::min(s.end, static_cast<double>(settledK + kContinuation) * kSample) : s.end; };
  while (amr.time() < runEnd() - 1e-13) {
    const double t = amr.time(), next = static_cast<double>(k + 1) * kSample;
    const double other = std::min({runEnd() - t, spray.nextEvent(t) - t, next - t}), bound = spray.exchangeDt();
    const double dt = amr.step(std::min(other, bound));
    smallestBound = std::min(smallestBound, bound);
    shortened += bound < other && dt >= bound * (1 - 1e-12);
    if (amr.time() < next - 1e-13) continue;
    ++k;
    const Sample x = sample();
    worstGrid = std::max(worstGrid, std::abs(x.t - static_cast<double>(k) * kSample));
    worstArea = std::max(worstArea, x.areaError);
    const double row[] = {x.mdot, x.fvac, x.ispChamber, x.ispInlet, x.pFace, x.qWall, x.qFace, x.liquid};
    for (double v : row) nonfinite += !std::isfinite(v);
    undefined += !std::isfinite(x.p0) || !std::isfinite(x.cstar);
    std::fprintf(csv, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.0f\n", x.t, x.mdot, x.fvac, x.ispChamber,
                 x.ispInlet, x.p0, x.pFace, x.cstar, x.qWall, x.qFace, x.liquid, x.parcels);
    series.push_back(x);
    if (k % 50 == 0 || amr.time() >= runEnd() - 1e-13) {
      std::fflush(csv);
      print(x);
      std::printf("  exchange bound: smallest %.3g s since the last line; %ld base steps shortened by it so far; %ld flashes so far, "
                  "%.4g kg flashed\n", smallestBound, shortened, spray.ledger().flashes, spray.ledger().flashedMass);
      smallestBound = std::numeric_limits<double>::infinity();
    }
    // The flash check (the Director's flash decision): no flash after t_L + 1 ms.
    if (k == leadK + 100) flashesAt = spray.ledger().flashes;
    if (k > leadK + 100 && spray.ledger().flashes > flashesAt && !flashReported) {
      flashReported = true;
      std::printf("\nA parcel flashed after t_L + 1 ms (by %.5g ms: %ld flashes since, %.4g kg flashed over the run) [%s]: %s\n",
                  x.t * 1e3, spray.ledger().flashes - flashesAt, spray.ledger().flashedMass, label.c_str(),
                  development ? "development goes on" : "run central stops");
      if (!development) {
        std::fclose(csv);
        dump("flash");
        return 5;
      }
    }
    if (k == leadK && leadK >= 200) {  // the lead's gate: the means over [t_L - 2, t_L - 1] and [t_L - 1, t_L] ms
      std::printf("\nThe lead's gate at t_L = %.4g ms [%s]\n", s.lead * 1e3, label.c_str());
      const char* names[] = {"outlet flow", "p0", "face pressure"};
      bool pass = true;
      for (int q = 0; q < 3; ++q) {
        double a = 0, b = 0;
        for (long n = leadK - 200; n < leadK; ++n) {
          const Sample& v = series[static_cast<std::size_t>(n)];
          const double value = q == 0 ? v.mdot : q == 1 ? v.p0 : v.pFace;
          (n < leadK - 100 ? a : b) += value / 100;
        }
        const double change = std::abs(b - a) / std::abs(a);
        const bool ok = change < 1e-3;  // NaN fails
        pass = pass && ok;
        std::printf("  %-14s %.7g, then %.7g: %.4f%%%s\n", names[q], a, b, 100 * change, ok ? "" : ": FAILS");
      }
      std::printf("  the gate %s\n", pass ? "passes" : development ? "FAILS (development: the run goes on)" : "FAILS: run central stops");
      if (!pass && !development) {
        std::fclose(csv);
        dump("lead");
        return 4;
      }
    }
    if (settledK == 0 && k - leadK >= kFirstCheck && k - leadK <= kCeiling && (k - leadK) % kCheckEvery == 0) {
      std::vector<double> face, cstar;
      for (const Sample& v : series) face.push_back(v.pFace), cstar.push_back(v.cstar);
      const RunCheck c = runCheck(k, leadK, face, cstar);
      std::printf("the settle check at t_L + %.0f ms (window from t_L + %.0f): face pressure SE %.3f%%, halves %.3f%%; c* SE %.3f%%, "
                  "halves %.3f%%: %s\n", static_cast<double>(k - leadK) * kSample * 1e3, static_cast<double>(k - leadK - kWindow) * kSample * 1e3,
                  100 * c.face.se, 100 * c.face.halves, 100 * c.cstar.se, 100 * c.cstar.halves,
                  c.check == SettleCheck::Settled ? "settled, t_s" : c.check == SettleCheck::Ceiling ? "NOT settled at the ceiling, t_s" : "not settled");
      if (c.check == SettleCheck::Settled || c.check == SettleCheck::Ceiling) {
        settledK = k;
        first = window(k, c.check == SettleCheck::Settled ? "the frozen parcel rate, t_s" : "the frozen parcel rate, the ceiling");
        dump("settled");
        if (s.end > static_cast<double>(k) * kSample + 1e-13) {
          spray.setParcelRate(2 * kParcelRate * s.rateFactor);
          std::printf("  from t_s = t_L + %.0f ms: the parcel rate doubled to %.3g per s per injector to t_s + 25 ms (item 5's check)\n",
                      static_cast<double>(k - leadK) * kSample * 1e3, 2 * kParcelRate * s.rateFactor);
        }
      }
    }
    if (settledK && k == settledK + kContinuation) {
      const auto second = window(k, "twice the parcel rate");
      dump("continuation");
      bool doubles = false;
      const char* names[] = {"face pressure", "c*", "F_vac", "outlet flow", "p0", "Isp_vac (chamber)", "hot-gas wall heat", "face heat"};
      for (std::size_t q = 0; q < first.size(); ++q) {
        if (q == 3) continue;  // the outlet flow is not an item-3 quantity
        const double move = std::abs(second[q].mean - first[q].mean) / std::abs(first[q].mean);
        const double combined = std::hypot(first[q].se * std::abs(first[q].mean), second[q].se * std::abs(second[q].mean)) / std::abs(first[q].mean);
        const bool moved = move > 0.005 || move > 2 * combined;
        doubles = doubles || moved;
        std::printf("  %-18s moves %.3f%% (bar 0.5%%, twice the combined SE %.3f%%)%s\n", names[q], 100 * move, 200 * combined,
                    moved ? ": MOVED" : "");
      }
      std::printf("  item 5: the parcel rate %s\n", doubles ? "DOUBLES for every run" : "holds");
    }
    checkpoint.killPoint(k);
    if (checkpoint.due(k) && amr.time() < runEnd() - 1e-13) save();
  }
  if (checkpoint.enabled()) save();
  stateOut();
  std::fclose(csv);
  if (settledK == 0) dump("end");
  std::printf("\nend [%s]: %ld samples, the grid to %.3g s; face area to %.3g at worst; %ld non-finite values; %ld samples with p0 or c* "
              "undefined (no forward flux); worst budgets (a) %.2e, (b) %.2e, (c) %.2e; %ld base steps shortened by the exchange "
              "bound; %ld flashes, %.4g kg flashed\n", label.c_str(), k, worstGrid, worstArea, nonfinite, undefined, worstBudget[0], worstBudget[1],
              worstBudget[2], shortened, spray.ledger().flashes, spray.ledger().flashedMass);
  return 0;
}

struct Ring {
  double r{}, w{};
  double h2InInner{}, h2InOuter{}, loxInner{}, loxOuter{}, h2OutInner{}, h2OutOuter{};
};
struct Stream {
  std::string name;
  double inner{}, outer{}, massFlow{}, velocity{}, ratio{};
  bool porous{false};
  [[nodiscard]] double area() const { return kPi * (outer * outer - inner * inner); }
};

}  // namespace

int main(int argc, char** argv) {
  try {
    std::vector<std::pair<int, int>> meshes;
    std::string mode;  // empty: the dry run
    if (argc > 1 && std::string(argv[1]) == "selftest") return selftest();
    if (argc > 1 && std::string(argv[1]) == "run") {
      mode = argc == 3 ? argv[2] : "";
      if (mode != "development" && mode != "central") throw std::invalid_argument("run development|central");
      if (mode == "central" && !kNumericsFrozen) {
        std::printf("run central: refused. The numerical settings (the y+ target and the tag thresholds) are not frozen "
                    "(tier1_run_definition.md, *Frozen inputs*); no Definition is built and no step is taken.\n");
        return 2;
      }
      int nz = 149, nr = 16;
      if (mode == "development" && std::getenv("CRUCIBLE_T1_BASE") &&
          std::sscanf(std::getenv("CRUCIBLE_T1_BASE"), "%dx%d", &nz, &nr) != 2)
        throw std::invalid_argument("CRUCIBLE_T1_BASE: NZxNR");
      meshes.push_back({nz, nr});
    }
    const bool run = !mode.empty();
    for (int a = run ? argc : 1; a < argc; ++a) {
      int nz = 0, nr = 0;
      if (std::sscanf(argv[a], "%dx%d", &nz, &nr) != 2) throw std::invalid_argument("mesh: NZxNR");
      meshes.push_back({nz, nr});
    }
    if (meshes.empty()) meshes = {{76, 8}, {149, 16}, {1192, 128}};

    // The element dimensions (tier1_skeleton.md, *Declared ranges*).
    const double dL = env("CRUCIBLE_T1_POST_ID", 3.54) * 1e-3;
    const double t = env("CRUCIBLE_T1_TIP", 0.272) * dL, h = env("CRUCIBLE_T1_GAP", 0.268) * dL;
    const double recess = env("CRUCIBLE_T1_RECESS", 0.5) * dL;
    const double outerOffset = env("CRUCIBLE_T1_OUTER_RING", 0.5);
    const std::string innerRule = std::getenv("CRUCIBLE_T1_INNER_RING") ? std::getenv("CRUCIBLE_T1_INNER_RING") : "pitch";
    if (innerRule != "pitch" && innerRule != "clearance") throw std::invalid_argument("CRUCIBLE_T1_INNER_RING: pitch or clearance");
    const double transpired = env("CRUCIBLE_T1_TRANSPIRED", 0.10);
    const double slotLength = env("CRUCIBLE_T1_SLOT_LENGTH", 0.07);
    const double faceP = env("CRUCIBLE_T1_FACE_PRESSURE", 2.8e6);
    const double stretching = env("CRUCIBLE_RADIAL_STRETCHING", 0);
    const double dO = dL + 2 * t, dA = dO + 2 * h;
    const double aL = kPerRing * kPi * dL * dL / 4, aH = kPerRing * kPi * (dA * dA - dO * dO) / 4;

    // The wall at the first mesh's nz gives R_c (the same at every nz).
    const double rc = rl10::wall(meshes.front().first).rc;
    const double pitch = std::sqrt(kPi * rc * rc / kElements);
    const double r8 = rc - outerOffset * pitch;
    const double r1Rule = innerRule == "pitch" ? r8 - (kRings - 1) * pitch : kPerRing * dA / (2 * kPi);
    // Ring 1 must hold its inner H2 slot: the radius at which that slot reaches the axis solves
    // r - A_L / (4 pi r) - t = sqrt(A_H / 2 pi). And every ring must clear the next: r_1 is the smallest radius at
    // or above both at which the evenly spaced rings do not overlap (tier1_setup.md, *Two member rules*).
    const double reach = t + std::sqrt(aH / (2 * kPi));
    const double r1Min = (reach + std::sqrt(reach * reach + aL / kPi)) / 2;
    auto layout = [&](double r1) {  // choice 1's rings from r_1 to r_8
      std::vector<Ring> rings(kRings);
      for (int i = 0; i < kRings; ++i) {
        Ring& g = rings[i];
        g.r = r1 + i * (r8 - r1) / (kRings - 1);
        g.w = aL / (2 * kPi * g.r);
        g.loxInner = g.r - g.w / 2;
        g.loxOuter = g.r + g.w / 2;
        g.h2InOuter = g.loxInner - t;
        const double inner2 = g.h2InOuter * g.h2InOuter - aH / (2 * kPi);
        const double roundOff = 1e-12 * g.h2InOuter * g.h2InOuter;  // at r1Min the slot reaches the axis exactly
        g.h2InInner = inner2 > roundOff ? std::sqrt(inner2) : inner2 > -roundOff ? 0 : -1;
        g.h2OutInner = g.loxOuter + t;
        g.h2OutOuter = std::sqrt(g.h2OutInner * g.h2OutInner + aH / (2 * kPi));
      }
      return rings;
    };
    auto ordered = [&](const std::vector<Ring>& rings) {
      bool ok = rings.front().h2InInner >= 0 && rings.back().h2OutOuter < rc;
      for (int i = 0; i < kRings; ++i) {
        const Ring& g = rings[i];
        ok = ok && g.h2InInner < g.h2InOuter && g.h2InOuter < g.loxInner && g.loxInner < g.loxOuter &&
             g.loxOuter < g.h2OutInner && g.h2OutInner < g.h2OutOuter;
        if (i + 1 < kRings) ok = ok && g.h2OutOuter < rings[i + 1].h2InInner;
      }
      return ok;
    };
    double r1 = std::max(r1Rule, r1Min), r1Clear = r1;
    if (!ordered(layout(r1))) {  // step out to a clear layout, then bisect back to the boundary
      double lo = r1, hi = r1;
      while (hi < r8 && !ordered(layout(hi))) lo = hi, hi += 1e-4;
      for (int it = 0; it < 60 && hi < r8; ++it) {
        const double mid = (lo + hi) / 2;
        (ordered(layout(mid)) ? hi : lo) = mid;
      }
      r1 = r1Clear = hi;
    }
    const double spacing = (r8 - r1) / (kRings - 1);
    if (run)
      std::printf("RL10A-3-3A tier 1, the run mode (%s; tier1_run_definition.md, *The run mode*): the gates on the base mesh, "
                  "then the run; no measured value is read\n", mode == "development" ? "DEVELOPMENT (protocol §6): not the reported case" : "central");
    else
      std::printf("RL10A-3-3A tier 1, central case dry run (tier1_setup.md): no step is taken, no measured value is read\n");
    std::printf("element: d_L %.4f mm, t %.4f mm (%.3f d_L), h %.4f mm (%.3f d_L), post OD %.4f mm, annulus OD %.4f mm, "
                "recess %.4f mm (%.3f d_L)\n", dL * 1e3, t * 1e3, t / dL, h * 1e3, h / dL, dO * 1e3, dA * 1e3,
                recess * 1e3, recess / dL);
    std::printf("face: R_c %.4f mm, pitch %.4f mm, r_8 %.4f mm (R_c - %.3f p), r_1 %.4f mm (%s %.4f mm, ring 1's "
                "smallest %.4f mm, clear of ring 2 from %.4f mm), ring spacing %.4f mm; per ring LOX area %.4f mm2, H2 area "
                "%.4f mm2\n", rc * 1e3, pitch * 1e3, r8 * 1e3, outerOffset, r1 * 1e3, innerRule.c_str(), r1Rule * 1e3,
                r1Min * 1e3, r1Clear * 1e3, spacing * 1e3, aL * 1e6, aH * 1e6);

    // The rings (choice 1).
    const std::vector<Ring> rings = layout(r1);
    std::printf("\nring  r [mm]    H2 inner slot [mm]   tip   LOX slot [mm] (w)          tip   H2 outer slot [mm]   cup width [mm]\n");
    for (int i = 0; i < kRings; ++i) {
      const Ring& g = rings[i];
      std::printf("%d  %8.3f   [%7.3f, %7.3f]   %.3f  [%7.3f, %7.3f] (%.4f)   %.3f  [%7.3f, %7.3f]   %.3f\n", i + 1,
                  g.r * 1e3, g.h2InInner * 1e3, g.h2InOuter * 1e3, t * 1e3, g.loxInner * 1e3, g.loxOuter * 1e3,
                  g.w * 1e3, t * 1e3, g.h2OutInner * 1e3, g.h2OutOuter * 1e3, (g.h2OutOuter - g.h2InInner) * 1e3);
    }

    // G2: the layout.
    std::printf("\n");
    check("G2", ordered(rings), strf("no overlap: ring 1's inner H2 slot from %.4f mm (axis clear), ring 8's outer to %.4f mm "
                             "(wall at %.4f mm)", rings.front().h2InInner * 1e3, rings.back().h2OutOuter * 1e3, rc * 1e3));
    double areaError = 0, cups = 0;
    for (const Ring& g : rings) {
      areaError = std::max({areaError, rel(kPi * (g.loxOuter * g.loxOuter - g.loxInner * g.loxInner), aL),
                            rel(kPi * (g.h2InOuter * g.h2InOuter - g.h2InInner * g.h2InInner), aH / 2),
                            rel(kPi * (g.h2OutOuter * g.h2OutOuter - g.h2OutInner * g.h2OutInner), aH / 2)});
      cups += kPi * (g.h2OutOuter * g.h2OutOuter - g.h2InInner * g.h2InInner);
    }
    check("G2", areaError < 1e-12, strf("slot areas: largest relative error %.3g against A_L and A_H / 2", areaError));

    // The streams: 16 H2 slots, 9 porous bands (choice 3).
    const double fuel = kOxidizer / kMixtureRatio, sweep = fuel * transpired;
    std::vector<Stream> streams;
    std::vector<std::array<double, 2>> gaps;
    gaps.push_back({0, rings.front().h2InInner});
    for (int i = 0; i + 1 < kRings; ++i) gaps.push_back({rings[i].h2OutOuter, rings[i + 1].h2InInner});
    gaps.push_back({rings.back().h2OutOuter, rc});
    double porous = 0;
    std::erase_if(gaps, [](const auto& g) { return g[1] - g[0] < 1e-9; });  // no band where a slot meets the axis or the next cup
    for (const auto& g : gaps) porous += kPi * (g[1] * g[1] - g[0] * g[0]);
    const double flux = sweep / porous;
    check("G2", rel(cups + porous, kPi * rc * rc) < 1e-12,
          strf("cups %.2f mm2 + porous %.2f mm2 = pi R_c^2 to %.3g", cups * 1e6, porous * 1e6,
              rel(cups + porous, kPi * rc * rc)));
    for (std::size_t g = 0; g < gaps.size(); ++g) {
      Stream s{"porous " + std::to_string(g), gaps[g][0], gaps[g][1]};
      s.massFlow = flux * s.area();
      s.porous = true;
      streams.push_back(s);
    }
    for (int i = 0; i < kRings; ++i) {
      const double each = (fuel - sweep) / (2 * kRings);
      streams.push_back({"H2 in " + std::to_string(i + 1), rings[i].h2InInner, rings[i].h2InOuter, each});
      streams.push_back({"H2 out " + std::to_string(i + 1), rings[i].h2OutInner, rings[i].h2OutOuter, each});
    }
    std::sort(streams.begin(), streams.end(), [](const Stream& a, const Stream& b) { return a.inner < b.inner; });

    // The medium, the liquid, and G3.
    const std::string mechanism = std::getenv("CRUCIBLE_MECHANISM") ? std::getenv("CRUCIBLE_MECHANISM")
                                                                     : CRUCIBLE_SOURCE_DIR "/tools/mechanisms/h2o2_highT.yaml";
    thermo::ReactionSource source(mechanism);
    // Fits from 60 K, as the spray's other cases (SPRAY.md, Results: S2): the drops' films and the gas near them reach
    // 61 to 67 K, below which fits from the mechanism's own lower limit turn negative (spray NaN, 23 October).
    const TransportFits fits = thermo::transportFits(mechanism, 60);
    Definition base;  // the medium with its transport data
    base.species = source.medium().species();
    base.transport = fits;
    const Medium medium = base.medium();
    const auto h2 = source.massFractions("H2:1"), o2 = source.massFractions("O2:1");
    const LiquidOxygen liquid = LiquidOxygen::load(kSaturation, medium.enthalpy(298.15, o2.data()));
    const double rhoL = liquid.at(kOxidizerT).density;
    double fuelSum = 0, sweepSum = 0, fluxSpread = 0;
    for (const Stream& s : streams) {
      fuelSum += s.massFlow;
      if (s.porous) {
        sweepSum += s.massFlow;
        fluxSpread = std::max(fluxSpread, rel(s.massFlow / s.area(), flux));
      }
    }
    const double oxSum = kRings * (kOxidizer / kRings);
    check("G3", rel(oxSum, kOxidizer) < 1e-12 && rel(fuelSum, fuel) < 1e-12,
          strf("oxidizer %.6f kg/s, fuel %.6f kg/s, O/F %.6f, chamber flow %.6f kg/s", oxSum, fuelSum, oxSum / fuelSum,
              oxSum + fuelSum));
    check("G3", rel(sweepSum / fuel, transpired) < 1e-12 && fluxSpread < 1e-12,
          strf("transpired %.6f of the fuel at %.6f kg/(m2 s) on every porous band (spread %.3g); ring O/F %.5f",
              sweepSum / fuel, flux, fluxSpread, (kOxidizer / kRings) / ((fuel - sweep) / kRings)));

    // G4: the injection numbers at the nominal face pressure, and each stream's viscosity ratio from its length scale:
    // slotLength times the hydraulic diameter 2 (r_o - r_i) at a slot or band, one pitch on the porous face.
    const double rhoG = faceP / (medium.gasConstant(h2.data()) * kFuelT);
    std::vector<double> diffusion(medium.size()), work;
    const double muG = medium.transport(kFuelT, faceP, h2.data(), diffusion.data(), work).viscosity;
    const double ringOF = (kOxidizer / kRings) / ((fuel - sweep) / kRings);
    std::printf("\nG4 (reported) at a nominal %.3g MPa face pressure, H2 at %.2f K: rho_H2 %.4f kg/m3, mu %.4g Pa s; "
                "LOX at %.2f K: rho_L %.2f kg/m3\n", faceP * 1e-6, kFuelT, rhoG, muG, kOxidizerT, rhoL);
    auto turbulence = [&](Stream& s) {
      s.velocity = s.massFlow / (rhoG * s.area());
      const double length = s.porous ? pitch : slotLength * 2 * (s.outer - s.inner);
      const double k = 1.5 * std::pow(kIntensity * s.velocity, 2);
      s.ratio = rhoG * k / (std::sqrt(k) / (std::pow(kCmu, 0.25) * length) * muG);
    };
    for (Stream& s : streams) turbulence(s);
    for (int i = 0; i < kRings; ++i) {
      const double vL = (kOxidizer / kRings) / (rhoL * aL);
      const double vG = ((fuel - sweep) / kRings) / (rhoG * aH);
      const double j = rhoG * vG * vG / (rhoL * vL * vL);
      std::printf("  ring %d: V_L %.4f m/s, V_H2 %.3f m/s, J %.3f, (V_g - V_L) / (V_L MR) %.3f\n", i + 1, vL, vG, j,
                  (vG - vL) / (vL * ringOF));
    }
    double slotLo = 1e300, slotHi = 0, wLo = 1e300, wHi = 0;  // every porous band has the same velocity and L
    for (const Stream& s : streams) {
      if (s.porous) continue;
      slotLo = std::min(slotLo, s.ratio), slotHi = std::max(slotHi, s.ratio);
      wLo = std::min(wLo, s.outer - s.inner), wHi = std::max(wHi, s.outer - s.inner);
    }
    std::printf("  transpired: %.4f m/s; viscosity ratios from I %.2f: slots %.4g to %.4g (L = %.3f D_h, %.4f to %.4f mm), "
                "porous %.4g (L = p, %.4f mm)\n", flux / rhoG, kIntensity, slotLo, slotHi, slotLength,
                slotLength * 2 * wLo * 1e3, slotLength * 2 * wHi * 1e3, streams[0].ratio, pitch * 1e3);

    // G8: the LOX drop size from Burick's curve at each ring's group, and the core length (tier1_setup.md,
    // *The LOX drop size and its release*). Every ring carries the same flows through the same areas.
    double qLo = 3.01, qHi = 10;  // d30OverMedian rises with q
    for (int it = 0; it < 200; ++it) {
      const double q = (qLo + qHi) / 2;
      (d30OverMedian(q) > 0.455 ? qHi : qLo) = q;
    }
    const double q = (qLo + qHi) / 2;
    const double smdOverMedian = 1 / (std::tgamma(1 - 1 / q) * std::pow(std::log(2.0), 1 / q));
    const double vL = (kOxidizer / kRings) / (rhoL * aL), vG = ((fuel - sweep) / kRings) / (rhoG * aH);
    const double group = (vG - vL) / (vL * ringOF), momentum = rhoG * vG * vG / (rhoL * vL * vL);
    const double medianRatio = burick(group);
    check("G8", std::isfinite(medianRatio),
          strf("every ring's group %.4f inside the digitised curve's table [%.2f, %.2f]", group, kBurick.front()[0],
              kBurick.back()[0]) +
              (group > kBurick.back()[0] ? strf(" or past its end: EXTRAPOLATED at slope %.3f, %.2f times the data's "
                                               "last group (about 11)", kBurickSlope, group / 11)
                                         : std::string()));
    const auto lox = liquid.at(kOxidizerT);
    const double oh = lox.viscosity / std::sqrt(lox.density * lox.tension * dL);
    const double core = dL * 2.85 * std::pow(rhoG / lox.density, -0.38) * std::pow(oh, 0.34) * std::pow(momentum, -0.13);
    const double coreLH = dL * 6 / std::sqrt(momentum);
    const double smd = std::isfinite(medianRatio) ? smdOverMedian * medianRatio * dL : 0;
    std::printf("  q %.4f (D30 / D-bar %.6f), SMD / D-bar %.4f; D-bar / D_L %.4f: SMD %.1f um central, %.1f um low "
                "(curve / 1.7), %.1f um high (curve x 1.3)\n", q, d30OverMedian(q), smdOverMedian, medianRatio, smd * 1e6,
                smd * kScatterLow * 1e6, smd * kScatterHigh * 1e6);
    std::printf("  core length: LOX at %.2f K: mu_l %.4g Pa s, sigma %.4g N/m; Oh %.4g, M = J %.3f; Porcheron %.2f mm "
                "(%.3f d_L, used), Lasheras-Hopfinger %.2f mm (%.3f d_L)\n", kOxidizerT, lox.viscosity, lox.tension, oh,
                momentum, core * 1e3, core / dL, coreLH * 1e3, coreLH / dL);

    // The Definition (tier1_skeleton.md, *Tier 1's settings*).
    auto define = [&](int nz, int nr, bool merged) {
      Definition d;
      d.experiment = Case::Chamber;
      d.nz = nz;
      d.nr = nr;
      d.contour = rl10::wall(nz).contour;
      d.species = medium.species();
      d.composition = source.massFractions("N2:1");
      d.backPressure = kAmbientP;
      d.ambientTemperature = kAmbientT;
      d.radialStretching = stretching;
      // Lever B's test (docs/evidence/SPEED.md): CRUCIBLE_T1_CFL sets Definition::cfl in development only; the central
      // case keeps its registered 0.4.
      if (mode == "development") d.cfl = env("CRUCIBLE_T1_CFL", d.cfl);
      d.transport = fits;
      d.wallTemperature = kWallT;  // the face's own 199.67 K (the plan, 21 October) waits on a per-segment wall temperature
      d.turbulence.enabled = true;
      d.turbulence.wallFunctions = true;
      d.turbulence.propertyCorrections = true;
      std::vector<double> dd(medium.size()), w;
      const double a0 = medium.soundSpeed(kAmbientT, d.composition.data());
      const double rho0 = kAmbientP / (medium.gasConstant(d.composition.data()) * kAmbientT);
      const double mu0 = medium.transport(kAmbientT, kAmbientP, d.composition.data(), dd.data(), w).viscosity;
      d.turbulence.ambientK = 9e-9 * a0 * a0;
      d.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
      auto supply = [&](const Stream& s) {
        Supply p;
        p.innerRadius = s.inner;
        p.outerRadius = s.outer;
        p.massFlow = s.massFlow;
        p.totalTemperature = kFuelT;
        p.composition = h2;
        p.opens = 0;
        p.ramp = kRamp;
        p.turbulenceIntensity = kIntensity;
        p.viscosityRatio = s.ratio;
        return p;
      };
      for (const Stream& s : streams)
        if (s.porous || !merged) d.supplies.push_back(supply(s));
      if (merged) {  // each ring's two H2 slots as one band across its cup, at the slots' flow
        for (int i = 0; i < kRings; ++i) {
          Stream s{"H2 band", rings[i].h2InInner, rings[i].h2OutOuter, (fuel - sweep) / kRings};
          turbulence(s);
          d.supplies.push_back(supply(s));
        }
        std::sort(d.supplies.begin(), d.supplies.end(),
                  [](const Supply& a, const Supply& b) { return a.innerRadius < b.innerRadius; });
      }
      // rl10dev's igniter, scaled from C1's (tier 1's igniter is open, tier1_skeleton.md *Open*).
      const double scale = rc / 0.025;
      d.igniter = {0.005 * scale, 0.015 * scale, 0.010 * scale, 0.5, 2e-4, 2e-4};
      return d;
    };
    // The run's Definition (tier1_run_definition.md, *The run mode*): face inlet (a) with each ring's LOX slot a
    // closed port (decision 6), the plate at the fuel's temperature (decision 8), the torch at 7 kJ over 1 ms (item 1).
    auto tier1 = [&](int nz, int nr) {
      Definition dm = define(nz, nr, false);
      dm.openAreaSupplies = true;
      for (int i = 0; i < kRings; ++i) {
        Supply p = dm.supplies.front();
        p.innerRadius = rings[i].loxInner;
        p.outerRadius = rings[i].loxOuter;
        p.massFlow = 0;
        p.composition = o2;
        dm.supplies.push_back(p);
      }
      dm.plateTemperature = kFuelT;
      dm.igniter.energy = 7e3;
      dm.igniter.duration = 1e-3;
      return dm;
    };

    // G1, G5, G6 per mesh; G7 at the finest.
    for (const auto& [nz, nr] : meshes) {
      std::printf("\n=== %d x %d (radial stretching %g)\n", nz, nr, stretching);
      std::optional<rl10::Wall> found;
      try {
        found = rl10::wall(nz);
      } catch (const std::exception& e) {
        check("G1", false, std::string("the contour: ") + e.what());
        continue;
      }
      const rl10::Wall& wall = *found;
      const Definition d = define(nz, nr, false);
      const Mesh mesh(d);
      const double throatR = mesh.radius[wall.throatStation];
      check("G1", rel(throatR, wall.rt) < 1e-12 && std::abs(wall.face - rl10::kNominalFace) <= 0.03 &&
                      rel(mesh.radius[0], rc) < 1e-12 && rel(mesh.radius[nz], wall.rt * std::sqrt(61.0)) < 1e-12,
            strf("throat station %.0f radius to %.3g; face at %.4f in (%+.4f from -13.1)", wall.throatStation,
                rel(throatR, wall.rt), wall.face, wall.face - rl10::kNominalFace) +
                strf("; face radius to %.3g; exit radius to %.3g", rel(mesh.radius[0], rc), rel(mesh.radius[nz], wall.rt * std::sqrt(61.0))) +
                strf("; dz %.4f mm, face ring heights %.4f to %.4f mm",
                    mesh.dz * 1e3, (mesh.fraction[1] - mesh.fraction[0]) * rc * 1e3,
                    (mesh.fraction[nr] - mesh.fraction[nr - 1]) * rc * 1e3));
      // Flow's rule: a face ring belongs to the supply that holds its centre.
      int resolved = 0, empty = 0, worstIndex = -1;
      double worst = 0;
      std::printf("  supply         [inner, outer] mm        true mm2   rings  mesh mm2   A_mesh/A_true  J_mesh/J\n");
      for (std::size_t s = 0; s < streams.size(); ++s) {
        const Stream& st = streams[s];
        double area = 0;
        int count = 0;
        for (int j = 0; j < nr; ++j) {
          const double centre = mesh.ringMiddle(0, j);
          if (centre >= st.inner && centre < st.outer) { area += mesh.axialArea(0, j); ++count; }
        }
        const double ratio = area / st.area();
        if (count == 0) ++empty;
        if (std::abs(ratio - 1) <= 0.05) ++resolved;
        if (count > 0 && std::abs(ratio - 1) > worst) { worst = std::abs(ratio - 1); worstIndex = static_cast<int>(s); }
        std::printf("  %-12s [%8.3f, %8.3f]  %9.2f   %4d  %9.2f   %8.4f   %s\n", st.name.c_str(), st.inner * 1e3,
                    st.outer * 1e3, st.area() * 1e6, count, area * 1e6, ratio,
                    count ? strf("%.4f", 1 / (ratio * ratio)).c_str() : "-");
      }
      check("G5", resolved == static_cast<int>(streams.size()),
            strf("%.0f of %.0f supplies within 5%% of their true area; %.0f hold no ring centre", resolved, streams.size(), empty) +
                (worstIndex >= 0 ? " (worst held: " + streams[worstIndex].name + strf(", %.3f)", 1 + worst) : std::string()));
      // The merged layout's bands: each carries its ring's two slots' flow over the cup, so its J against the slots'
      // is (A_H / A_mesh)^2.
      std::printf("  merged band    [inner, outer] mm      slots mm2   cup mm2   rings  mesh mm2   J_mesh/J\n");
      for (int i = 0; i < kRings; ++i) {
        const Ring& g = rings[i];
        double area = 0;
        int count = 0;
        for (int j = 0; j < nr; ++j) {
          const double centre = mesh.ringMiddle(0, j);
          if (centre >= g.h2InInner && centre < g.h2OutOuter) { area += mesh.axialArea(0, j); ++count; }
        }
        std::printf("  ring %d       [%8.3f, %8.3f]  %9.2f  %8.2f   %4d  %9.2f   %s\n", i + 1, g.h2InInner * 1e3,
                    g.h2OutOuter * 1e3, aH * 1e6, kPi * (g.h2OutOuter * g.h2OutOuter - g.h2InInner * g.h2InInner) * 1e6,
                    count, area * 1e6, count ? strf("%.4f", std::pow(aH / area, 2)).c_str() : "-");
      }
      for (bool merged : {false, true}) {
        const char* name = merged ? "merged (each ring's H2 slots as one band across its cup)" : "true slots";
        try {
          Definition dm = define(nz, nr, merged);
          Flow flow(dm);
          Spray::Settings settings;
          const Medium& fm = flow.medium();
          for (int i = 0; i < kRings; ++i) {
            Spray::Injector j;
            j.rInner = rings[i].loxInner;
            j.rOuter = rings[i].loxOuter;
            j.massFlow = kOxidizer / kRings;
            j.opens = 0;
            j.ramp = kRamp;
            j.temperature = kOxidizerT;
            j.size = Spray::Injector::Size::RosinRammler;
            j.sauterDiameter = smd;
            j.exponent = q;
            j.coreLength = core;
            settings.injectors.push_back(j);
          }
          settings.restitution = 0.5;
          [[maybe_unused]] Spray spray(flow, LiquidOxygen::load(kSaturation, fm.enthalpy(298.15, o2.data())), settings);
          check("G6", true, std::string(name) + strf(": Flow and Spray constructed (%.0f supplies, %.0f injectors, %.0f cells)",
                                                    dm.supplies.size(), settings.injectors.size(), flow.mesh().cells.size()));
        } catch (const std::exception& e) {
          check("G6", false, std::string(name) + ": " + e.what());
        }
      }
      if (run || std::getenv("CRUCIBLE_T1_OPEN_AREA")) {
        // F3 (docs/evidence/FACE_INLET.md): the true layout under openAreaSupplies, each ring's LOX slot a
        // closed port, the valves open (ramp 0) so the face carries every supply's flow at t = 0. Each ring's
        // face mass flow against sum_s a_sj mdot_s / A_s (a_sj and A_s the true annuli's, computed here);
        // each supply's face mass flow, sum_j a_sj (rho u)_sj, with rho u from supplyFace against ring j's
        // first cell, against mdot_s; the engine's supply areas against the true. On AMR (3 levels, every
        // cell active) the wall functions and the corrections stay on (AMR.md, step D, D4).
        auto open = [&]() {
          if (run) {  // the run's own Definition, the valves open
            Definition dm = tier1(nz, nr);
            for (Supply& p : dm.supplies) p.ramp = 0;
            return dm;
          }
          Definition dm = define(nz, nr, false);
          dm.openAreaSupplies = true;
          for (int i = 0; i < kRings; ++i) {
            Supply p = dm.supplies.front();
            p.innerRadius = rings[i].loxInner;
            p.outerRadius = rings[i].loxOuter;
            p.massFlow = 0;
            p.composition = o2;
            dm.supplies.push_back(p);
          }
          for (Supply& p : dm.supplies) p.ramp = 0;
          return dm;
        };
        auto judge = [&](Flow& flow, const std::string& label) {
          const Definition& dm = flow.definition();
          const Mesh& fm = flow.mesh();
          const double face = fm.radius[0];
          const auto flux = flow.injectorFluxes();
          auto share = [&](const Supply& p, int j) {
            const double lo = face * fm.fraction[j], hi = face * fm.fraction[j + 1];
            return std::max(0.0, kPi * (std::pow(std::min(p.outerRadius, hi), 2) - std::pow(std::max(p.innerRadius, lo), 2)));
          };
          double ringError = 0, supplyError = 0, areaError = 0;
          int touched = 0;
          for (int j = 0; j < fm.nr; ++j) {
            double expected = 0;
            for (const Supply& p : dm.supplies) expected += share(p, j) * p.massFlow * p.opening(0) / (kPi * (p.outerRadius * p.outerRadius - p.innerRadius * p.innerRadius));
            const double got = fm.axialArea(0, j) * flux[static_cast<std::size_t>(j)][0];
            ringError = std::max(ringError, expected == 0 ? (got == 0 ? 0.0 : 1.0) : rel(got, expected));
          }
          for (std::size_t s = 0; s < dm.supplies.size(); ++s) {
            const Supply& p = dm.supplies[s];
            const double area = kPi * (p.outerRadius * p.outerRadius - p.innerRadius * p.innerRadius);
            areaError = std::max(areaError, rel(flow.supplyAreas()[s], area));
            const double g = p.massFlow * p.opening(0) / flow.supplyAreas()[s];
            double sum = 0;
            for (int j = 0; j < fm.nr; ++j) {
              const double a = share(p, j);
              if (!(a > 0)) continue;
              const std::size_t q = fm.index(0, j);
              const auto y = flow.massFractions(q);
              const double rhoU = flow.supplyFace(p, g, flow.cellPrimitive(q), y.data()).flux[0];
              sum += a * rhoU;
              if (p.massFlow > 0) { supplyError = std::max(supplyError, rel(rhoU * area, p.massFlow)); ++touched; }
              else if (rhoU != 0) supplyError = 1;
            }
            if (p.massFlow > 0) supplyError = std::max(supplyError, rel(sum, p.massFlow));
          }
          check("F3", ringError < 1e-12 && supplyError < 1e-12 && areaError < 1e-12,
                label + strf(": %.0f supplies; per ring, face mass flow to %.3g; per supply, face mass flow and rho u A_true "
                            "(%.0f supply-rings)", dm.supplies.size(), ringError, touched) +
                    strf(" to %.3g; supply areas to %.3g (limit 1e-12)", supplyError, areaError));
        };
        try {
          Flow flow(open());
          judge(flow, strf("%.0f x ", nz) + strf("%.0f, wall functions and corrections", nr));
          flow.step();
          const auto m = flow.measurements();
          std::printf("     one step: dt %.3e s, budgets mass %.3e, energy %.3e\n", m.dt, m.massBalanceError, m.energyBalanceError);
        } catch (const std::exception& e) {
          check("F3", false, std::string("Flow: ") + e.what());
        }
        if (nz == 149 && nr == 16) {
          try {
            Amr amr(open(), 2);
            for (int l = 1; l < amr.levels(); ++l) amr.setActive(l, std::vector<unsigned char>(amr.level(l).mesh().cells.size(), 1));
            for (int l = 0; l < amr.levels(); ++l) judge(amr.level(l), strf("AMR level %.0f, wall functions and corrections (step D, D4)", l));
          } catch (const std::exception& e) {
            check("F3", false, std::string("AMR: ") + e.what());
          }
        }
      }
      if (nz == meshes.back().first && nr == meshes.back().second) {
        // G7 at the finest mesh given: the cup's criterion (tier1_skeleton.md).
        double minSlot = 1e9, minTip = 1e9;
        for (const Ring& g : rings) {
          auto across = [&](double a, double b) {  // the face rings' heights over [a, b]: cells across it
            double n = 0;
            for (int j = 0; j < nr; ++j) {
              const double lo = mesh.fraction[j] * rc, hi = mesh.fraction[j + 1] * rc;
              n += std::max(0.0, std::min(hi, b) - std::max(lo, a)) / (hi - lo);
            }
            return n;
          };
          minSlot = std::min({minSlot, across(g.h2InInner, g.h2InOuter), across(g.loxInner, g.loxOuter),
                              across(g.h2OutInner, g.h2OutOuter)});
          minTip = std::min({minTip, across(g.h2InOuter, g.loxInner), across(g.loxOuter, g.h2OutInner)});
        }
        check("G7", false,
              strf("the cup: %.2f cells across its depth (need 4), fewest across a slot %.2f and a tip %.2f (need 2); ",
                  recess / mesh.dz, minSlot, minTip) +
                  "Mesh is one body-fitted map of r(z) from the face plane, with no cells behind it, so it cannot carry the cup's walls");
      }
    }
    std::printf("\n%d check(s) failed\n", failures);
    if (!run) return 0;
    std::printf("gates (G1 to G4, G8, F3 on the run's Definition): %d failed; G5 and G6 test the ring-centre rule decision 6 "
                "superseded, G7 is declared (decision 5)\n", gateFailures);
    if (gateFailures) {
      std::printf("run %s: stopped at the gates, no step taken\n", mode.c_str());
      return 3;
    }
    Setup setup;
    setup.mode = mode;
    setup.mechanism = mechanism;
    const auto& [nz, nr] = meshes.front();
    setup.definition = tier1(nz, nr);
    setup.core = core;
    setup.cylinderEnd = (-10.00 - rl10::wall(nz).face) * rl10::kInch;
    setup.chamberFlow = oxSum + fuelSum;
    setup.o2 = o2;
    setup.out = "rl10_t1_central";
    if (mode == "development") {
      setup.levels = static_cast<int>(env("CRUCIBLE_T1_LEVELS", 3));
      setup.threads = static_cast<int>(env("CRUCIBLE_T1_THREADS", 4));
      setup.lead = env("CRUCIBLE_T1_LEAD", 20e-3);
      setup.end = env("CRUCIBLE_T1_END", setup.lead + 1e-3);  // the simulated time, not counted from t_L
      setup.rateFactor = env("CRUCIBLE_T1_RATE_FACTOR", 1);
      setup.out = std::string(std::getenv("CRUCIBLE_T1_OUT") ? std::getenv("CRUCIBLE_T1_OUT") : "rl10_t1") + "_DEVELOPMENT";
    } else {
      setup.end = setup.lead + 85e-3;
    }
    setup.definition.igniter.start = setup.lead + 2e-4;
    for (int i = 0; i < kRings; ++i) {
      Spray::Injector j;
      j.rInner = rings[i].loxInner;
      j.rOuter = rings[i].loxOuter;
      j.massFlow = kOxidizer / kRings;
      j.opens = setup.lead;
      j.ramp = kRamp;
      j.temperature = kOxidizerT;
      j.size = Spray::Injector::Size::RosinRammler;
      j.sauterDiameter = smd;
      j.exponent = q;
      j.coreLength = core;
      setup.spray.injectors.push_back(j);
    }
    setup.spray.restitution = 0.5;
    // Lever B's noise floor: CRUCIBLE_T1_SEED sets the spray's seed (Spray::Settings::seed, 1) in development only.
    if (setup.mode == "development") setup.spray.seed = static_cast<std::uint64_t>(env("CRUCIBLE_T1_SEED", 1));
    setup.argc = argc, setup.argv = argv;
    doneFile = setup.out + "_done";
    const int status = runCase(setup);
    if (std::FILE* f = std::fopen(doneFile.c_str(), "w")) std::fprintf(f, "exit %d\n", status), std::fclose(f);
    return status;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    if (!doneFile.empty())
      if (std::FILE* f = std::fopen(doneFile.c_str(), "w")) std::fprintf(f, "exit 1: %s\n", e.what()), std::fclose(f);
    return 1;
  }
}
