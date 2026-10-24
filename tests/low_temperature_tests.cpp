// The lowest ranges below 200 K (docs/evidence/LOW_T_THERMO.md; thermo::withLowTemperatureRanges), on
// h2o2_highT.yaml's medium. Tolerances stated before the run (registered 24 October 2026):
//   1. Nothing at or above 200 K changes, bit for bit: Medium's properties, gas constant, e, h, cv, sound speed and
//      species enthalpies, and ManifoldRuntime::entropy, for 64 compositions (the 10 pure species, then random) at
//      200 K exactly and 2,000 temperatures to 20,000 K (the 1,000 K and 6,000 K joins and both sides of them).
//   2. The temperature inversion of every e at or above e(200 K) is bit for bit the same, from guesses on both sides
//      of 200 K (20, 60, 150, 199.99 K, the target itself, 0.7 and 1.3 times it, 300, 3,000 and 20,000 K, none).
//   3. At 200 K the lowest ranges join TM-4513's: |cp jump| / cp <= 1e-14, |h jump| <= 1e-12 of |h| plus 1e-9 J/kg,
//      |s jump| <= 1e-12 of s, per species (H2, O2, H2O), taken between 200 K and the double below it.
//   4. Below 200 K the inversion round-trips (|T - T0| <= 1e-9 T0) from 20 to 199.99 K for every composition,
//      cv > 0 from 1 K, and the energy floor lies below e(20 K).
//   5. A shock tube of H2, O2 and H2O (Y 0.2, 0.3, 0.5), 10 bar 1,500 K | 1 bar 300 K, every state above 200 K:
//      the run with the lowest ranges is bit for bit the run without them (state and partial densities after the
//      same steps). The control, the right side at 120 K: the two runs differ (the check discriminates).
//   6. The lowest ranges' cp against independent tables: H2 within 1.5% of normal hydrogen's ideal-gas Cp0 (NIST
//      WebBook, Leachman et al. 2009) at 100, 150 and 199.9 K (against the 200 K entry); O2 within 0.2% of JANAF
//      (O-029) at 100 K and of the WebBook's Cp0 at 150 K; H2O within 0.5% of JANAF (H-064) at 100 K. And H2 is
//      normal, not equilibrium hydrogen: more than 15% from JANAF's H2 (H-050, equilibrium) at 100 K.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "adapters/reaction.hpp"
#include "core/flow.hpp"
#include "core/manifold_runtime.hpp"

using namespace crucible;

namespace {
int failures = 0;
void check(bool ok, const char* what, double value, double tol) {
  std::printf("  %-60s %12.4e  (tol %.1e) %s\n", what, value, tol, ok ? "ok" : "FAIL");
  failures += !ok;
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }
bool same(const std::vector<double>& a, const std::vector<double>& b) {
  return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

struct Run { std::vector<double> state, partial; int steps{}; };
Run shockTube(const std::vector<Species>& species, const std::vector<double>& y, double tRight) {
  const int n = 400;
  Definition d;
  d.nz = n; d.nr = 1; d.length = 1; d.inletRadius = d.throatRadius = d.exitRadius = 0.05;
  d.experiment = Case::ShockTube;
  d.species = species;
  d.composition = y;
  Flow f(d);
  const double r = f.medium().gasConstant(y.data());
  std::vector<Primitive> cells;
  std::vector<double> ys;
  for (int i = 0; i < n; ++i) {
    const bool left = (i + 0.5) / n < 0.5;
    const double p = left ? 1e6 : 1e5, t = left ? 1500 : tRight;
    cells.push_back({p / (r * t), 0, 0, p});
    ys.insert(ys.end(), y.begin(), y.end());
  }
  f.setInitialState(cells, ys);
  Run out;
  while (f.time() < 2e-4) { f.step(); ++out.steps; }
  for (const auto& u : f.state()) out.state.insert(out.state.end(), u.begin(), u.end());
  out.partial = f.partialDensities();
  return out;
}
}  // namespace

int main() {
  thermo::ReactionSource source(CRUCIBLE_SOURCE_DIR "/tools/mechanisms/h2o2_highT.yaml");
  const std::vector<Species> before = source.medium().species();
  const std::vector<Species> after = thermo::withLowTemperatureRanges(before);
  const Medium old(before), now(after);
  const std::size_t ns = old.size();
  std::vector<std::vector<double>> compositions;
  for (std::size_t k = 0; k < ns; ++k) {
    compositions.emplace_back(ns, 0.0);
    compositions.back()[k] = 1;
  }
  std::mt19937_64 rng(20261024);
  std::uniform_real_distribution<double> unit(0, 1);
  while (compositions.size() < 64) {
    std::vector<double> y(ns);
    double sum = 0;
    for (auto& v : y) sum += v = std::pow(unit(rng), 4);
    for (auto& v : y) v /= sum;
    compositions.push_back(y);
  }
  std::vector<double> temperatures{200.0, std::nextafter(1000.0, 0.0), 1000.0, std::nextafter(6000.0, 0.0), 6000.0,
                                   std::nextafter(6000.0, 1e9), 20000.0};
  for (int i = 0; i < 2000; ++i) temperatures.push_back(200 * std::pow(100.0, i / 1999.0));

  std::printf("1. At and above 200 K, every property bit for bit\n");
  {
    long compared = 0, differ = 0;
    std::vector<double> ha(ns), hb(ns);
    for (const auto& y : compositions)
      for (double t : temperatures) {
        const auto a = old.properties(t, y.data()), b = now.properties(t, y.data());
        old.speciesEnthalpies(t, ha.data());
        now.speciesEnthalpies(t, hb.data());
        const bool ok = same(a.r, b.r) && same(a.e, b.e) && same(a.cv, b.cv) &&
                        same(old.gasConstant(y.data()), now.gasConstant(y.data())) &&
                        same(old.internalEnergy(t, y.data()), now.internalEnergy(t, y.data())) &&
                        same(old.enthalpy(t, y.data()), now.enthalpy(t, y.data())) &&
                        same(old.cv(t, y.data()), now.cv(t, y.data())) &&
                        same(old.soundSpeed(t, y.data()), now.soundSpeed(t, y.data())) && same(ha, hb) &&
                        same(ManifoldRuntime::entropy(old, t, 1.0, y.data(), 1e5),
                             ManifoldRuntime::entropy(now, t, 1.0, y.data(), 1e5));
        ++compared;
        differ += !ok;
      }
    std::printf("  %ld composition-temperature pairs\n", compared);
    check(differ == 0, "pairs with any property not bit for bit", static_cast<double>(differ), 0);
  }

  std::printf("2. The inversion of every e at or above e(200 K), bit for bit, from guesses either side of 200 K\n");
  {
    long compared = 0, differ = 0;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto& y : compositions)
      for (std::size_t i = 0; i < temperatures.size(); i += 7) {
        const double t = temperatures[i], e = old.internalEnergy(t, y.data());
        for (double g : {20.0, 60.0, 150.0, 199.99, t, 0.7 * t, 1.3 * t, 300.0, 3000.0, 20000.0, 0.0, nan}) {
          ++compared;
          differ += !same(old.temperature(e, y.data(), g), now.temperature(e, y.data(), g));
        }
      }
    std::printf("  %ld inversions\n", compared);
    check(differ == 0, "inversions not bit for bit", static_cast<double>(differ), 0);
  }

  std::printf("3. The joins at 200 K (between 200 K and the double below it)\n");
  {
    const double below = std::nextafter(200.0, 0.0);
    for (const char* name : {"H2", "O2", "H2O"}) {
      std::size_t k = 0;
      while (after[k].name != name) ++k;
      std::vector<double> y(ns, 0.0);
      y[k] = 1;
      const double cpJump = std::abs(now.cv(below, y.data()) - now.cv(200, y.data())) / now.cv(200, y.data());
      const double h = now.enthalpy(200, y.data()), hJump = std::abs(now.enthalpy(below, y.data()) - h);
      const double s = ManifoldRuntime::entropy(now, 200, 1.0, y.data(), 1e5);
      const double sJump = std::abs(ManifoldRuntime::entropy(now, below, 1.0, y.data(), 1e5) - s) / std::abs(s);
      std::printf("  %s:\n", name);
      check(cpJump <= 1e-14, "  |cv jump| / cv", cpJump, 1e-14);
      check(hJump <= 1e-12 * std::abs(h) + 1e-9, "  |h jump| [J/kg]", hJump, 1e-12 * std::abs(h) + 1e-9);
      check(sJump <= 1e-12, "  |s jump| / s", sJump, 1e-12);
    }
  }

  std::printf("4. Below 200 K: round trips from 20 to 199.99 K, cv > 0 from 1 K, the floor below e(20 K)\n");
  {
    double worst = 0, leastCv = std::numeric_limits<double>::infinity(), floorMargin = std::numeric_limits<double>::infinity();
    for (const auto& y : compositions) {
      for (int i = 0; i <= 400; ++i) {
        const double t0 = 20 + (199.99 - 20) * i / 400.0, e = now.internalEnergy(t0, y.data());
        for (double g : {t0, 300.0, 20.0})
          worst = std::max(worst, std::abs(now.temperature(e, y.data(), g) - t0) / t0);
      }
      for (double t = 1; t < 200; t += 0.5) leastCv = std::min(leastCv, now.cv(t, y.data()));
      floorMargin = std::min(floorMargin, now.internalEnergy(20, y.data()) - now.energyFloor(y.data()));
    }
    check(worst <= 1e-9, "largest |T - T0| / T0", worst, 1e-9);
    check(leastCv > 0, "least cv from 1 to 200 K [J/(kg K)]", leastCv, 0);
    check(floorMargin > 0, "least e(20 K) - floor [J/kg]", floorMargin, 0);
  }

  std::printf("5. A shock tube above 200 K, with and without the lowest ranges; the control at 120 K\n");
  {
    std::vector<double> y(ns, 0.0);
    for (std::size_t k = 0; k < ns; ++k) {
      if (before[k].name == "H2") y[k] = 0.2;
      if (before[k].name == "O2") y[k] = 0.3;
      if (before[k].name == "H2O") y[k] = 0.5;
    }
    const Run a = shockTube(before, y, 300), b = shockTube(after, y, 300);
    std::printf("  300 K: %d and %d steps\n", a.steps, b.steps);
    check(a.steps == b.steps && same(a.state, b.state) && same(a.partial, b.partial),
          "runs above 200 K not bit for bit (1 if they differ)", a.steps == b.steps && same(a.state, b.state) &&
          same(a.partial, b.partial) ? 0 : 1, 0);
    const Run c = shockTube(before, y, 120), e = shockTube(after, y, 120);
    std::printf("  120 K: %d and %d steps\n", c.steps, e.steps);
    const bool differ = c.steps != e.steps || !same(c.state, e.state);
    check(differ, "the control differs (1 if it does)", differ ? 1 : 0, 0);
  }

  std::printf("6. The lowest ranges' cp against independent tables (J/(mol K) over R = 8.314462618)\n");
  {
    struct Reference { const char* name; double t, cp, tol; const char* source; };
    const Reference references[] = {
        {"H2", 100, 22.571, 0.015, "H2 100 K, NIST normal H2 Cp0"},
        {"H2", 150, 25.384, 0.015, "H2 150 K, NIST normal H2 Cp0"},
        {"H2", 199.9, 27.278, 0.015, "H2 199.9 K, NIST normal H2 Cp0 at 200 K"},
        {"O2", 100, 29.106, 0.002, "O2 100 K, JANAF O-029"},
        {"O2", 150, 29.112, 0.002, "O2 150 K, NIST Cp0"},
        {"H2O", 100, 33.299, 0.005, "H2O 100 K, JANAF H-064"},
    };
    auto cpOverR = [&](const char* name, double t) {
      std::vector<double> y(ns, 0.0);
      for (std::size_t k = 0; k < ns; ++k) y[k] = after[k].name == name;
      const double r = now.gasConstant(y.data());
      return (now.cv(t, y.data()) + r) / r;
    };
    for (const auto& ref : references) {
      const double error = std::abs(cpOverR(ref.name, ref.t) / (ref.cp / 8.314462618) - 1);
      check(error <= ref.tol, (std::string("|cp/table - 1|, ") + ref.source).c_str(), error, ref.tol);
    }
    const double fromEquilibrium = std::abs(cpOverR("H2", 100) / (28.154 / 8.314462618) - 1);
    check(fromEquilibrium > 0.15, "H2 100 K, |cp/JANAF H-050 - 1| (must exceed)", fromEquilibrium, 0.15);
  }

  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
