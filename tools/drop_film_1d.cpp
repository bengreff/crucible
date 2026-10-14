// S2's reference (docs/evidence/SPRAY.md): the steady, spherically symmetric evaporation of a liquid
// oxygen drop at its wet bulb in stagnant gas, resolved in r with variable properties, against the
// drop model's rates (Spray::rates at Re = 0, with its film integrals at several point counts; 0 points is
// the 1/3 rule). Both use transport fitted down to 60 K (thermo::transportFits).
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
#include <cstdlib>
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
    // Mixture-averaged D_km = (1 - Y_k) / sum_{j != k} X_j / D_kj over the species present (the carriers
    // and O2), with 1 - Y_k summed from the others and their X floored at 1e-20 as Cantera floors them,
    // which gives the pure-gas limits (far away, and pure O2 at a boiling surface).
    auto xf = [&](std::size_t k) { return std::max(x[k], 1e-20); };
    auto yOf = [&](std::size_t k) { return x[k] * gas->molecularWeight(k) / w; };
    double restO2 = 0, sumO2 = 0;
    for (auto c : carriers) {
      restO2 += xf(c) * gas->molecularWeight(c) / w;
      sumO2 += xf(c) / d[n * c + o2];
    }
    const double cO2 = gas->molecularWeight(o2) / w * restO2 / sumO2;
    double a = 0;
    std::vector<double> cc(carriers.size());
    for (std::size_t i = 0; i < carriers.size(); ++i) {
      const auto c = carriers[i];
      double rest = xf(o2) * gas->molecularWeight(o2) / w, sum = xf(o2) / d[n * o2 + c];
      for (auto j : carriers)
        if (j != c) {
          rest += xf(j) * gas->molecularWeight(j) / w;
          sum += xf(j) / d[n * j + c];
        }
      cc[i] = gas->molecularWeight(c) / w * rest / sum;
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

// A root of f in the bracket [a, b] (f(a), f(b) of opposite signs) by the Illinois method, to |b - a| <= tol.
template <class F>
double illinois(F&& f, double a, double fa, double b, double fb, double tol) {
  for (int it = 0; it < 200 && std::abs(b - a) > tol; ++it) {
    const double c = b - fb * (b - a) / (fb - fa), fc = f(c);
    if (fc == 0) return c;
    if ((fc > 0) != (fb > 0)) {
      a = b;
      fa = fb;
    } else {
      fa /= 2;
    }
    b = c;
    fb = fc;
  }
  return std::abs(fa) < std::abs(fb) ? a : b;
}

// q for a given T_s (T(1) = T_s), on ln q. T(1) falls as q rises; a q that drives the film below 20 K
// fails the shot, so the bracket's upper end is moved back until its shot succeeds.
double solveQ(Film& film, const std::vector<double>& far, double ts, double hSurface, int steps) {
  auto miss = [&](double lq, bool& valid) {
    const auto s = film.shoot(far, std::exp(lq), hSurface, steps);
    valid = !s.empty();
    return valid ? s[0] - ts : -1.0;
  };
  bool valid = true;
  double lo = std::log(1e-9), flo = miss(lo, valid);
  for (int it = 0; flo < 0; ++it) {
    if (it == 40) throw std::runtime_error("drop_film_1d: no q leaves the film above T_s");
    flo = miss(lo -= 2, valid);
  }
  double hi = lo + 1, fhi = miss(hi, valid);
  for (int it = 0; fhi > 0; ++it) {
    if (it == 80) throw std::runtime_error("drop_film_1d: no q brings the film below T_s");
    lo = hi;
    flo = fhi;
    fhi = miss(hi += 1, valid);
  }
  for (int it = 0; !valid; ++it) {  // between a good positive miss and a failed shot: halve toward the good end
    if (it == 200)
      throw std::runtime_error("drop_film_1d: q bracket did not close at T_s " + std::to_string(ts) + ", ln q " +
                               std::to_string(lo) + " (miss " + std::to_string(flo) + ")");
    const double mid = 0.5 * (lo + hi);
    const double fm = miss(mid, valid);
    if (!valid) { hi = mid; continue; }
    if (fm > 0) { lo = mid; flo = fm; valid = false; continue; }
    hi = mid;
    fhi = fm;
  }
  return std::exp(illinois([&](double lq) { bool v = true; const double m = miss(lq, v); return v ? m : -1.0; }, lo, flo, hi, fhi, 1e-13));
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
  // On u = ln(T_sat - T_s), from 1e-13 K below saturation down by decades to the first sign change (the
  // film's transport fits are not asked about gas far colder than the wet bulb), at worst 1 K above the
  // table's foot.
  double q = 0, others = 0;
  double lo = std::log(1e-13);
  const double foot = std::log(tSat - liquid.minTemperature() - 1);
  double glo = g(tSat - std::exp(lo), q, others);
  if (glo < 0) {
    // The root lies closer to saturation than a double resolves: the boiling limit.
    out.boiling = true;
    out.ts = tSat;
    out.q = solveQ(film, far, tSat, hSurface(tSat), steps);
    const auto s = film.shoot(far, out.q, hSurface(tSat), steps);
    for (std::size_t i = 1; i < s.size(); ++i) out.surfaceOthers += std::exp(s[i]);
    return out;
  }
  double hi = lo, ghi = glo;
  while (ghi > 0) {
    if (hi >= foot) throw std::runtime_error("drop_film_1d: no wet bulb above the table's foot");
    lo = hi;
    glo = ghi;
    hi = std::min(hi + std::log(10.0), foot);
    ghi = g(tSat - std::exp(hi), q, others);
  }
  const double u = illinois([&](double v) { return g(tSat - std::exp(v), q, others); }, lo, glo, hi, ghi, 1e-10);
  out.ts = tSat - std::exp(u);
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
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  thermo::ReactionSource source("h2o2.yaml");
  // Cantera's phase with every species' range extended down to the fit floor, so its transport fits
  // are the engine's (thermo::transportFits with the same floor).
  constexpr double fitFloor = 60;
  auto root = Cantera::AnyMap::fromYamlFile("h2o2.yaml");
  for (auto& species : root["species"].asVector<Cantera::AnyMap>()) {
    auto& ranges = species["thermo"]["temperature-ranges"].asVector<double>();
    ranges.front() = std::min(ranges.front(), fitFloor);
  }
  auto solution = Cantera::newSolution(root["phases"].asVector<Cantera::AnyMap>().front(), root, "mixture-averaged");
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
  def.transport = thermo::transportFits("h2o2.yaml", fitFloor);
  Flow flow(def);
  const auto& medium = flow.medium();
  std::vector<double> yO2(medium.size(), 0.0);
  yO2[film.o2] = 1;
  const auto liquid = LiquidOxygen::load(kSaturation, medium.enthalpy(298.15, yO2.data()));
  // The engine's default film quadrature is judged; 0 (the 1/3 rule alone, as A-S) and finer ones are shown.
  const int pointsShown[] = {0, 2, 3, 4, 8, 16};
  std::vector<std::unique_ptr<Spray>> sprays;
  for (int points : pointsShown) {
    Spray::Settings settings;
    settings.twoWay = false;
    settings.filmPoints = points;
    sprays.push_back(std::make_unique<Spray>(flow, liquid, settings));
  }
  const int judged = Spray::Settings{}.filmPoints;
  std::vector<double> work(4 * medium.size() + 8);

  struct Case { const char* name; double p, t, yH2O; };
  const Case cases[] = {{"H2, 300 K, 1 MPa", 1e6, 300, 0},      {"H2, 300 K, 3 MPa", 3e6, 300, 0},
                        {"H2, 1,000 K, 1 MPa", 1e6, 1000, 0},   {"H2, 1,000 K, 3 MPa", 3e6, 1000, 0},
                        {"H2/H2O 0.4/0.6, 1,500 K, 1 MPa", 1e6, 1500, 0.6}, {"H2/H2O 0.4/0.6, 1,500 K, 3 MPa", 3e6, 1500, 0.6}};
  std::printf("S2. Wet-bulb drops in stagnant gas: resolved film (Cantera, %d and %d RK4 steps in a / r) against the model "
              "(judged at its default, %d film points)\n", steps, 2 * steps, judged);
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
    ModelResult mine;
    std::string shown;
    for (std::size_t i = 0; i < sprays.size(); ++i) {
      const auto r = model(*sprays[i], liquid, gas, work);
      char line[200];
      std::snprintf(line, sizeof line, "\n      %2d points: T_wb %.6f K%s, mdot/d %.6e kg/(m s) (%+.2f%%)", pointsShown[i],
                    r.ts, r.boiling ? " (boiling)" : "", r.massPerD, 100 * (r.massPerD / (2 * kPi * fine.q) - 1));
      shown += line;
      if (pointsShown[i] == judged) mine = r;
    }
    const double refPerD = 2 * kPi * fine.q;  // mdot / d = 2 pi mdot / (4 pi a)
    const double dm = mine.massPerD / refPerD - 1, dt = mine.ts - fine.ts;
    const bool ok = std::abs(dm) <= 0.10 && std::abs(dt) <= 2;
    if (!ok) ++failures;
    std::printf("  %s\n", c.name);
    std::printf("    reference: T_wb %.6f K%s, mdot/d %.6e kg/(m s), sum X_other at the surface %.3e; "
                "steps doubled: T_wb %+.2e K, mdot/d %+.2e\n",
                fine.ts, fine.boiling ? " (saturation: the root is within rounding of it)" : "", refPerD, fine.surfaceOthers,
                fine.ts - ref.ts, fine.q / ref.q - 1);
    std::printf("    model:%s\n", shown.c_str());
    std::printf("    mdot/d %+.2f%% (band 10%%), T_wb %+.3f K (band 2 K)  %s\n", 100 * dm, dt, ok ? "PASS" : "FAIL");
  }
  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
