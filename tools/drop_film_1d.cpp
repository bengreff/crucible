// S2's reference (docs/evidence/SPRAY.md): the steady, spherically symmetric evaporation of a liquid
// oxygen drop at its wet bulb in stagnant gas, resolved in r with variable properties, against the
// drop model's rates (Spray::rates at Re = 0, with the 1/3 rule).
//
// The film is solved with Cantera alone (h2o2.yaml, ideal gas, mixture-averaged diffusion with a
// correction velocity, mixture conductivity); the surface condition (Raoult with NIST's vapour pressure)
// and the latent heat are the model's own LiquidOxygen table, so the comparison is of the film model.
//
// The equations, in xi = a / r (xi = 0 far away, 1 at the surface), with q = mdot / (4 pi a):
//   every species but O2 has no net flux, O2 carries mdot:  4 pi r^2 j_k = mdot (eps_k - Y_k);
//   j_k = -rho (W_k / W) D_km dX_k/dr + rho Y_k V_c, V_c = sum_m (W_m / W) D_mm dX_m/dr (sum_k j_k = 0);
//   energy: 4 pi r^2 lambda dT/dr = mdot (h_O2(T) - h_O2(T_s) + L(T_s)) (the wet bulb: no heat into the
//   liquid, and the vapour formed at T_s takes L).
// With c_k = (W_k / W) D_km and A = sum_c Y_c (1 - c_O2 / c_c) over the species c other than O2:
//   d ln X_c / dxi = -q W_c / (W rho c_c (1 - A)),   dT/dxi = -q (h_O2(T) - h_O2(T_s) + L) / lambda,
// integrated by RK4 from the far state at xi = 0. q is set by T(1) = T_s and T_s by
// sum_c X_c(1) = 1 - p_sat(T_s) / p.
//
// Usage: crucible_drop_film_1d [steps]   (default 2000; each case is also run at twice the steps)
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "adapters/reaction.hpp"
#include "cantera/core.h"
#include "core/spray.hpp"

using namespace crucible;

namespace {
constexpr double kPi = 3.14159265358979323846;
const std::string kSaturation = CRUCIBLE_SOURCE_DIR "/tools/o2_saturation_nist.tsv";

struct Film {
  std::shared_ptr<Cantera::ThermoPhase> gas;
  std::shared_ptr<Cantera::Transport> transport;
  std::size_t n{}, o2{};
  std::vector<std::size_t> carriers;  // species present far away, O2 aside
  std::vector<double> x, d, hrt;
  double p{};
  // At (T, ln X_c): the right-hand side for the state (T, ln X_c...) given q, T_s and h_O2(T_s) - L.
  void rhs(const std::vector<double>& s, double q, double hSurface, std::vector<double>& out) {
    std::fill(x.begin(), x.end(), 0.0);
    double others = 0;
    for (std::size_t i = 0; i < carriers.size(); ++i) others += x[carriers[i]] = std::exp(s[i + 1]);
    x[o2] = std::max(1 - others, 0.0);
    const double t = s[0];
    if (!(t > 20 && t < 6000) || !(others < 1 + 1e-12)) {  // a trial stage far off the solution
      std::fill(out.begin(), out.end(), std::nan(""));
      return;
    }
    gas->setState_TPX(t, p, x.data());
    const double rho = gas->density(), w = gas->meanMolecularWeight();
    const double lambda = transport->thermalConductivity();
    transport->getBinaryDiffCoeffs(n, d.data());
    gas->getEnthalpy_RT(hrt.data());
    const double hO2 = hrt[o2] * Cantera::GasConstant * t / gas->molecularWeight(o2);
    // Mixture-averaged D_km = (1 - Y_k) / sum_{j != k} X_j / D_kj, with 1 - Y_O2 summed from the others.
    auto yOf = [&](std::size_t k) { return x[k] * gas->molecularWeight(k) / w; };
    double yOthers = 0, sumO2 = 0;
    for (auto c : carriers) {
      yOthers += yOf(c);
      sumO2 += x[c] / d[n * c + o2];
    }
    const double cO2 = gas->molecularWeight(o2) / w * yOthers / sumO2;
    double a = 0;
    std::vector<double> cc(carriers.size());
    for (std::size_t i = 0; i < carriers.size(); ++i) {
      const auto c = carriers[i];
      double sum = 0;
      for (std::size_t j = 0; j < n; ++j)
        if (j != c && x[j] > 0) sum += x[j] / d[n * j + c];
      cc[i] = gas->molecularWeight(c) / w * (1 - yOf(c)) / sum;
      a += yOf(c) * (1 - cO2 / cc[i]);
    }
    out[0] = -q * (hO2 - hSurface) / lambda;
    for (std::size_t i = 0; i < carriers.size(); ++i)
      out[i + 1] = -q * gas->molecularWeight(carriers[i]) / (w * rho * cc[i] * (1 - a));
  }
  // The surface state (T(1), ln X_c(1)) for q and T_s.
  std::vector<double> shoot(const std::vector<double>& far, double q, double hSurface, int steps) {
    std::vector<double> s = far, k1(s.size()), k2(s.size()), k3(s.size()), k4(s.size()), tmp(s.size());
    const double h = 1.0 / steps;
    for (int i = 0; i < steps; ++i) {
      rhs(s, q, hSurface, k1);
      for (std::size_t m = 0; m < s.size(); ++m) tmp[m] = s[m] + 0.5 * h * k1[m];
      rhs(tmp, q, hSurface, k2);
      for (std::size_t m = 0; m < s.size(); ++m) tmp[m] = s[m] + 0.5 * h * k2[m];
      rhs(tmp, q, hSurface, k3);
      for (std::size_t m = 0; m < s.size(); ++m) tmp[m] = s[m] + h * k3[m];
      rhs(tmp, q, hSurface, k4);
      for (std::size_t m = 0; m < s.size(); ++m) s[m] += h / 6 * (k1[m] + 2 * k2[m] + 2 * k3[m] + k4[m]);
      if (!(s[0] > 20)) return {};  // overshot far below any surface temperature (or a trial stage failed)
    }
    return s;
  }
};

struct Solution {
  double ts{}, q{}, surfaceOthers{};
  bool boiling{};
};

// q for a given T_s (T(1) = T_s), by bisection then secant on ln q.
double solveQ(Film& film, const std::vector<double>& far, double ts, double hSurface, int steps) {
  auto miss = [&](double lq) {
    const auto s = film.shoot(far, std::exp(lq), hSurface, steps);
    return s.empty() ? -1e3 : s[0] - ts;
  };
  double lo = std::log(1e-12), hi = std::log(1e-2);
  while (miss(hi) > 0) hi += 2;
  while (miss(lo) < 0) lo -= 2;
  for (int it = 0; it < 200 && hi - lo > 1e-13; ++it) {
    const double mid = 0.5 * (lo + hi);
    (miss(mid) > 0 ? lo : hi) = mid;
  }
  return std::exp(0.5 * (lo + hi));
}

Solution solve(Film& film, const LiquidOxygen& liquid, const std::vector<double>& far, int steps) {
  const double tSat = liquid.saturationTemperature(film.p);
  // Surface enthalpy terms: h_O2(T_s) - L(T_s), with Cantera's h_O2.
  auto hSurface = [&](double ts) {
    std::vector<double> hrt(film.n);
    film.gas->setState_TP(ts, film.p);
    film.gas->getEnthalpy_RT(hrt.data());
    return hrt[film.o2] * Cantera::GasConstant * ts / film.gas->molecularWeight(film.o2) - liquid.at(ts).latent;
  };
  // g(T_s) = ln sum_c X_c(1) - ln(1 - p_sat(T_s) / p): positive near saturation, negative well below.
  auto g = [&](double ts, double& q, double& others) {
    const double hs = hSurface(ts);
    q = solveQ(film, far, ts, hs, steps);
    const auto s = film.shoot(far, q, hs, steps);
    others = 0;
    for (std::size_t i = 1; i < s.size(); ++i) others += std::exp(s[i]);
    const double room = 1 - liquid.at(ts).pressure / film.p;
    return std::log(others) - (room > 0 ? std::log(room) : -1e300);
  };
  Solution out;
  // Bisection in ln(T_sat - T_s) from 60 K below saturation to 1e-13 K below.
  double lo = std::log(1e-13), hi = std::log(tSat - liquid.minTemperature() - 1);
  double q = 0, others = 0;
  if (g(tSat - std::exp(lo), q, others) < 0) {
    // The root lies closer to saturation than a double resolves: the boiling limit.
    out.boiling = true;
    out.ts = tSat;
    out.q = solveQ(film, far, tSat, hSurface(tSat), steps);
    const auto s = film.shoot(far, out.q, hSurface(tSat), steps);
    for (std::size_t i = 1; i < s.size(); ++i) out.surfaceOthers += std::exp(s[i]);
    return out;
  }
  for (int it = 0; it < 60; ++it) {
    const double mid = 0.5 * (lo + hi);
    (g(tSat - std::exp(mid), q, others) > 0 ? lo : hi) = mid;
  }
  out.ts = tSat - std::exp(0.5 * (lo + hi));
  g(out.ts, out.q, out.surfaceOthers);
  return out;
}

// The model at Re = 0: its wet bulb by bisection on Q - mdot L (or saturation when it boils there),
// and mdot / d.
struct ModelResult {
  double ts{}, massPerD{};
  bool boiling{};
};
ModelResult model(const Spray& spray, const LiquidOxygen& liquid, const Spray::Gas& gas, std::vector<double>& work) {
  Parcel p;
  p.d = 100e-6;
  p.count = 1;
  auto at = [&](double t) {
    Parcel q = p;
    q.t = t;
    return spray.rates(q, gas, work);
  };
  const double tSat = liquid.saturationTemperature(gas.p);
  double lo = liquid.minTemperature(), hi = tSat;
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) break;
    const auto r = at(mid);
    if (r.boiling) { hi = mid; continue; }
    (r.heat - r.mass * liquid.at(r.surface).latent > 0 ? lo : hi) = mid;
  }
  ModelResult out;
  out.boiling = at(hi).boiling;
  out.ts = out.boiling ? tSat : hi;
  const auto r = at(out.ts);
  out.massPerD = r.mass / p.d;
  return out;
}
}  // namespace

int main(int argc, char** argv) {
  const int steps = argc > 1 ? std::atoi(argv[1]) : 2000;
  thermo::ReactionSource source("h2o2.yaml");
  auto solution = Cantera::newSolution("h2o2.yaml", "", "mixture-averaged");
  Film film;
  film.gas = solution->thermo();
  film.transport = solution->transport();
  film.n = film.gas->nSpecies();
  film.o2 = film.gas->speciesIndex("O2");
  film.x.resize(film.n);
  film.d.resize(film.n * film.n);
  film.hrt.resize(film.n);

  Definition def;
  def.experiment = Case::Chamber;
  def.nz = 4;
  def.nr = 2;
  def.contour = {{0, 0.005}, {0.01, 0.005}};
  def.species = source.medium().species();
  def.composition = source.massFractions("H2:1");
  def.backPressure = 1e6;
  def.ambientTemperature = 300;
  def.transport = thermo::transportFits("h2o2.yaml");
  Flow flow(def);
  const auto& medium = flow.medium();
  std::vector<double> yO2(medium.size(), 0.0);
  yO2[film.o2] = 1;
  const auto liquid = LiquidOxygen::load(kSaturation, medium.enthalpy(298.15, yO2.data()));
  Spray::Settings settings;
  settings.twoWay = false;
  Spray spray(flow, liquid, settings);
  std::vector<double> work(4 * medium.size() + 8);

  struct Case { const char* name; double p, t, yH2O; };
  const Case cases[] = {{"H2, 300 K, 1 MPa", 1e6, 300, 0},      {"H2, 300 K, 3 MPa", 3e6, 300, 0},
                        {"H2, 1,000 K, 1 MPa", 1e6, 1000, 0},   {"H2, 1,000 K, 3 MPa", 3e6, 1000, 0},
                        {"H2/H2O 0.4/0.6, 1,500 K, 1 MPa", 1e6, 1500, 0.6}, {"H2/H2O 0.4/0.6, 1,500 K, 3 MPa", 3e6, 1500, 0.6}};
  std::printf("S2. Wet-bulb drops in stagnant gas: resolved film (Cantera, %d and %d RK4 steps in a / r) against the model\n",
              steps, 2 * steps);
  int failures = 0;
  for (const auto& c : cases) {
    std::vector<double> y(film.n, 0.0);
    y[film.gas->speciesIndex("H2")] = 1 - c.yH2O;
    y[film.gas->speciesIndex("H2O")] = c.yH2O;
    film.p = c.p;
    film.gas->setState_TPY(c.t, c.p, y.data());
    std::vector<double> xFar(film.n);
    film.gas->getMoleFractions(xFar.data());
    film.carriers.clear();
    std::vector<double> far{c.t};
    for (std::size_t k = 0; k < film.n; ++k)
      if (k != film.o2 && xFar[k] > 0) {
        film.carriers.push_back(k);
        far.push_back(std::log(xFar[k]));
      }
    const auto ref = solve(film, liquid, far, steps);
    const auto fine = solve(film, liquid, far, 2 * steps);
    const Spray::Gas gas{c.t, c.p, c.p / (medium.gasConstant(y.data()) * c.t), 0, 0, y.data()};
    const auto mine = model(spray, liquid, gas, work);
    const double refPerD = 2 * kPi * fine.q;  // mdot / d = 2 pi mdot / (4 pi a)
    const double dm = mine.massPerD / refPerD - 1, dt = mine.ts - fine.ts;
    const bool ok = std::abs(dm) <= 0.10 && std::abs(dt) <= 2;
    if (!ok) ++failures;
    std::printf("  %s\n", c.name);
    std::printf("    reference: T_wb %.6f K%s, mdot/d %.6e kg/(m s), sum X_other at the surface %.3e; "
                "steps doubled: T_wb %+.2e K, mdot/d %+.2e\n",
                fine.ts, fine.boiling ? " (saturation: the root is within rounding of it)" : "", refPerD, fine.surfaceOthers,
                fine.ts - ref.ts, fine.q / ref.q - 1);
    std::printf("    model:     T_wb %.6f K%s, mdot/d %.6e kg/(m s)\n", mine.ts, mine.boiling ? " (boiling)" : "", mine.massPerD);
    std::printf("    mdot/d %+.2f%% (band 10%%), T_wb %+.3f K (band 2 K)  %s\n", 100 * dm, dt, ok ? "PASS" : "FAIL");
  }
  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
