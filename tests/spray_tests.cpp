// Verification of the liquid oxygen spray (core/spray.hpp) on the criteria of docs/evidence/SPRAY.md,
// stated before the code: S0 (properties), S1 (the drop model's algebra and integrator), S3 (budgets
// with parcels), S4 (threads), S5 (breakup) and S8 (turbulent dispersion). S2's reference is tools/drop_film_1d.cpp.
// With the argument "cost" it measures S6 instead: the S3 run with and without the spray; with "lafon",
// the model's side of S2(c) (outer step dt, default 2 us).
// Usage: crucible_spray_tests [cost [threads] | lafon [dt]]
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/spray.hpp"

using namespace crucible;

namespace {
constexpr double kPi = std::numbers::pi;
const std::string kSaturation = CRUCIBLE_SOURCE_DIR "/tools/o2_saturation_nist.tsv";
const std::string kLowPressure = CRUCIBLE_SOURCE_DIR "/tools/o2_lowp_nist.tsv";

int failures = 0;
double sq(double x) { return x * x; }
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-78s %.3e (limit %.0e)  %s\n", what, value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}
// Relative difference; infinite if either value is not finite, so a NaN cannot pass.
double relative(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b)) return std::numeric_limits<double>::infinity();
  return a == b ? 0.0 : std::abs(a - b) / std::max(std::abs(a), std::abs(b));
}

std::size_t oxygenIndex(const Medium& m) {
  for (std::size_t k = 0; k < m.size(); ++k)
    if (m.species()[k].name == "O2") return k;
  throw std::runtime_error("no O2");
}

// A round chamber with a straight wall, nz by nr, filled at rest with `fill` (mole fractions) at p and t,
// with h2o2.yaml's species and mixture-averaged transport (no-slip adiabatic walls).
Definition chamber(thermo::ReactionSource& source, int nz, int nr, double length, double radius, const std::string& fill,
                   double p, double t) {
  Definition d;
  d.experiment = Case::Chamber;
  d.nz = nz;
  d.nr = nr;
  d.contour = {{0, radius}, {length, radius}};
  d.species = source.medium().species();
  d.composition = source.massFractions(fill);
  d.backPressure = p;
  d.ambientTemperature = t;
  d.transport = thermo::transportFits("h2o2.yaml", 60);  // fits valid down to cold O2 (SPRAY.md, Results: S2)
  return d;
}

LiquidOxygen liquidFor(const Medium& m) {
  std::vector<double> y(m.size(), 0.0);
  y[oxygenIndex(m)] = 1;
  return LiquidOxygen::load(kSaturation, m.enthalpy(298.15, y.data()));
}

// The gas a parcel in cell q sees without loads, computed as Spray::gatherGas does.
struct CellGas {
  std::vector<double> y;
  Spray::Gas gas;
};
CellGas cellGas(const Flow& flow, std::size_t q) {
  const auto& medium = flow.medium();
  const auto ns = medium.size();
  CellGas out;
  out.y.resize(ns);
  const auto& partial = flow.partialDensities();
  double sum = 0;
  for (std::size_t k = 0; k < ns; ++k) sum += std::max(partial[q * ns + k], 0.0);
  for (std::size_t k = 0; k < ns; ++k) out.y[k] = std::max(partial[q * ns + k], 0.0) / sum;
  const auto& u = flow.state()[q];
  const auto& turbulence = flow.turbulence();
  const double k = turbulence.empty() ? 0 : turbulence[q * 2] / u[0];
  const double e = (u[3] - (u[1] * u[1] + u[2] * u[2]) / (2 * u[0])) / u[0] - k;
  const double t = medium.temperature(e, out.y.data(), 300);
  out.gas = {t, u[0] * medium.gasConstant(out.y.data()) * t, u[0], u[1] / u[0], u[2] / u[0], out.y.data()};
  return out;
}

// ---------------------------------------------------------------------------------------------------
// S0

void properties(const Medium& medium, const LiquidOxygen& liquid) {
  std::printf("S0. Properties\n");
  double worst = 0, worstSat = 0;
  const int n = 10000;
  const double t0 = liquid.minTemperature(), t1 = liquid.maxTemperature();
  for (int i = 0; i < n; ++i) {
    const double t = t0 + (t1 - t0) * (i + 0.5) / n;
    worst = std::max(worst, relative(liquid.temperature(liquid.enthalpy(t)), t));
    worstSat = std::max(worstSat, relative(liquid.saturationTemperature(liquid.at(t).pressure), t));
  }
  std::printf("  table %.4f to %.4f K\n", t0, t1);
  check(worst < 1e-12, "(a) T(h(T)) at 10,000 temperatures across the table", worst, 1e-12);
  check(worstSat < 1e-12, "(a) T_sat(p_sat(T)) likewise (the vapour-pressure inverse)", worstSat, 1e-12);

  std::vector<double> y(medium.size(), 0.0);
  y[oxygenIndex(medium)] = 1;
  const double h298 = medium.enthalpy(298.15, y.data());
  std::ifstream in(kLowPressure);
  if (!in) throw std::runtime_error("cannot open " + kLowPressure);
  double worstH = 0, at = 0, low = 1e9, high = 0;
  for (std::string line; std::getline(in, line);) {
    if (line.empty() || line[0] == '#' || line.rfind("Temperature", 0) == 0) continue;
    std::istringstream fields(line);
    std::vector<double> v;
    for (std::string f; std::getline(fields, f, '\t');) {
      try { v.push_back(std::stod(f)); } catch (const std::exception&) { v.push_back(std::nan("")); }
    }
    const double t = v[0], nist = h298 + v[5] * 1e3 - LiquidOxygen::idealGasEnthalpy298;
    const double diff = std::abs(medium.enthalpy(t, y.data()) - nist);
    if (diff > worstH) { worstH = diff; at = t; }
    for (double show : {75.0, 90.0, 100.0, 120.0, 150.0, 200.0, 250.0, 300.0})
      if (t == show) std::printf("  %5.0f K: h_medium - h_NIST %+8.1f J/kg\n", t, medium.enthalpy(t, y.data()) - nist);
    low = std::min(low, t);
    high = std::max(high, t);
  }
  std::printf("  medium O2 against NIST's 1e-6 MPa isobar, %.0f to %.0f K (the table's range): largest at %.0f K\n", low, high, at);
  check(worstH < 1000, "(b) |h_medium - h_NIST| on the engine's zero [J/kg]", worstH, 1000);
}

// ---------------------------------------------------------------------------------------------------
// S1

// Abramzon and Sirignano (1989) written out from the paper (SPRAY.md, Literature), independently of
// Spray::rates: the drop at tDrop (or at saturation when boiling) and slip speed `slip` in `gas`.
struct Expected {
  double reynolds{}, spalding{}, sherwood{}, nusselt{}, transfer{}, mass{}, heat{}, viscosity{};
};
double filmF(double b) { return std::pow(1 + b, 0.7) * std::log1p(b) / b; }
Expected abramzonSirignano(const Medium& m, const LiquidOxygen& liquid, double d, double tDrop, double slip,
                           const Spray::Gas& gas, bool boiling) {
  const auto ns = m.size();
  const auto o2 = oxygenIndex(m);
  const double ts = boiling ? liquid.saturationTemperature(gas.p) : tDrop;
  double ys = 1, wOther = 0;
  if (!boiling) {
    const double x = liquid.at(tDrop).pressure / gas.p;
    double num = 0, den = 0;
    for (std::size_t k = 0; k < ns; ++k)
      if (k != o2) { num += gas.y[k]; den += gas.y[k] / m.species()[k].molarMass; }
    wOther = num / den;
    const double wO2 = m.species()[o2].molarMass;
    ys = x * wO2 / (x * wO2 + (1 - x) * wOther);
  }
  std::vector<double> surface(ns), film(ns), diffusion(ns), work;
  for (std::size_t k = 0; k < ns; ++k) surface[k] = k == o2 ? ys : (1 - ys) * gas.y[k] / (1 - gas.y[o2]);
  const double tr = ts + (gas.t - ts) / 3;
  for (std::size_t k = 0; k < ns; ++k) film[k] = surface[k] + (gas.y[k] - surface[k]) / 3;
  std::vector<double> pure(ns, 0.0);
  pure[o2] = 1;
  const auto tp = m.transport(tr, gas.p, film.data(), diffusion.data(), work);
  const auto fp = m.properties(tr, film.data()), vp = m.properties(tr, pure.data());
  const double cpg = fp.cv + fp.r, cpf = vp.cv + vp.r;
  // Fick's law on Y_O2: the Blanc's-law D, W / W_c times Medium's D_km (SPRAY.md, Results: S2).
  const double rhoD = gas.p / (fp.r * tr) * diffusion[o2] * (wOther > 0 ? Medium::universalGasConstant / (fp.r * wOther) : 1);
  Expected e;
  e.viscosity = tp.viscosity;
  e.reynolds = gas.rho * slip * d / tp.viscosity;
  const double sc = tp.viscosity / rhoD, pr = cpg * tp.viscosity / tp.conductivity;
  const double le = tp.conductivity / (rhoD * cpg);
  const double f = std::max(1.0, std::pow(e.reynolds, 0.077));
  const double sh0 = 1 + std::cbrt(1 + e.reynolds * sc) * f, nu0 = 1 + std::cbrt(1 + e.reynolds * pr) * f;
  if (boiling) {
    const double latent = liquid.at(ts).latent;
    e.transfer = cpf * (gas.t - ts) / latent;
    e.nusselt = 2 + (nu0 - 2) / filmF(e.transfer);
    e.sherwood = e.nusselt;
    e.mass = kPi * d * tp.conductivity / cpf * e.nusselt * std::log1p(e.transfer);
    e.heat = e.mass * latent;
    return e;
  }
  e.spalding = (ys - gas.y[o2]) / (1 - ys);
  e.sherwood = 2 + (sh0 - 2) / filmF(e.spalding);
  // B_T = (1 + B_M)^phi - 1 with phi = (c_pF / c_pg)(Sh* / Nu*(B_T)) / Le, by bisection.
  auto g = [&](double bt) {
    const double nu = 2 + (nu0 - 2) / filmF(bt);
    return std::expm1(cpf / cpg * e.sherwood / nu / le * std::log1p(e.spalding)) - bt;
  };
  double lo = 0, hi = 1e6;
  for (int it = 0; it < 400; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) break;
    (g(mid) > 0 ? lo : hi) = mid;
  }
  e.transfer = std::abs(g(lo)) < std::abs(g(hi)) ? lo : hi;
  e.nusselt = 2 + (nu0 - 2) / filmF(e.transfer);
  e.mass = kPi * d * rhoD * e.sherwood * std::log1p(e.spalding);
  e.heat = kPi * d * tp.conductivity * e.nusselt * (gas.t - ts) * std::log1p(e.transfer) / e.transfer;
  return e;
}

// (a) One 100 um drop at its wet bulb in `fill` at rest at 1 MPa and tGas, one-way: d^2 against the
// d^2 law at 1,000 times, and its end against (d0^2 - d_min^2) / K.
void dSquared(thermo::ReactionSource& source, const char* fill, double tGas, int filmPoints = Spray::Settings{}.filmPoints) {
  Flow flow(chamber(source, 4, 2, 0.01, 0.005, fill, 1e6, tGas));
  const auto& medium = flow.medium();
  const auto liquid = liquidFor(medium);
  Spray::Settings settings;
  settings.twoWay = false;
  settings.breakup.on = false;
  settings.filmPoints = filmPoints;
  Spray probe(flow, liquid, settings);
  Parcel p;
  p.z = 0.005;
  p.r = 0.001;
  p.d = 100e-6;
  p.count = 1;
  p.t = 90;
  probe.add(p);
  const auto cg = cellGas(flow, static_cast<std::size_t>(probe.parcels()[0].cell));
  std::vector<double> work(4 * medium.size() + 8);
  auto balance = [&](double t) {
    Parcel q = p;
    q.t = t;
    const auto r = probe.rates(q, cg.gas, work);
    return r.heat - r.mass * liquid.at(r.surface).latent;
  };
  // Bisection on Q - mdot L. In hot hydrogen the wet bulb lies closer to saturation than a double
  // resolves (Y_s rounds to 1 first), and the drop is then held at saturation (boiling).
  const double tSat = liquid.saturationTemperature(cg.gas.p);
  double lo = liquid.minTemperature(), hi = tSat;
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (mid <= lo || mid >= hi) break;
    Parcel q = p;
    q.t = mid;
    if (probe.rates(q, cg.gas, work).boiling) { hi = mid; continue; }
    (balance(mid) > 0 ? lo : hi) = mid;
  }
  Parcel q = p;
  q.t = hi;
  const bool boiling = probe.rates(q, cg.gas, work).boiling;
  const double wetBulb = boiling ? tSat : (std::abs(balance(lo)) < std::abs(balance(hi)) ? lo : hi);
  Spray d2(flow, liquid, settings);
  p.t = wetBulb;
  d2.add(p);
  const auto r0 = d2.rates(d2.parcels()[0], cg.gas, work);
  const double rho = liquid.at(wetBulb).density, d0 = p.d, dMin = settings.minDiameter;
  const double k = 4 * r0.mass / (rho * kPi * d0);  // -d(d^2)/dt
  const double life = (d0 * d0 - dMin * dMin) / k;
  std::printf("  %s at %.0f K, %d film points: wet bulb %.9f K (T_sat %.9f K, %s), B_M %.4e, K %.6e m^2/s, time to d_min %.6e s\n", fill, tGas,
              filmPoints, wetBulb, tSat, r0.boiling ? "boiling" : "evaporating", r0.spalding, k, life);
  const int n = 1000;
  const double end = life * (1 - 1e-10), dt = end / n;
  double t = 0, worst = 0;
  bool alive = true;
  for (int i = 0; i < n; ++i) {
    d2.advance(dt);
    t += dt;
    if (d2.parcels().size() != 1) { alive = false; break; }
    const double d = d2.parcels()[0].d;
    worst = std::max(worst, std::abs(d * d - (d0 * d0 - k * t)) / (d0 * d0));
  }
  std::printf("  %ld sub-steps; drop temperature at the end %.9f K\n", d2.ledger().substeps,
              alive ? d2.parcels()[0].t : std::nan(""));
  check(alive && worst < 1e-10, "(a) d^2 against d0^2 - K t at 1,000 times (relative to d0^2)", alive ? worst : 1.0, 1e-10);
  d2.advance(2e-10 * life);
  const bool gone = d2.parcels().empty() && d2.ledger().vaporised == 1;
  check(alive && gone, "(a) alive at (1 - 1e-10) and vaporised by (1 + 1e-10) of (d0^2 - d_min^2)/K", alive && gone ? 0.0 : 1.0, 1e-10);
}

void dropModel(thermo::ReactionSource& source) {
  std::printf("S1. The drop model's algebra and integrator\n");
  dSquared(source, "N2:1", 300);
  // 1,200 K, not 1,000: at NASA-7's 1,000 K junction the inverse T(e) has two roots 2.6e-6 K apart and
  // which one it returns depends on the guess, so the drop's seen saturation point would too.
  dSquared(source, "H2:1", 1200);
  Flow flow(chamber(source, 4, 2, 0.01, 0.005, "H2:1", 1e6, 1000));
  const auto& medium = flow.medium();
  const auto ns = medium.size();
  const auto liquid = liquidFor(medium);
  Spray::Settings forms;
  forms.filmPoints = 0;  // the 1/3 rule's conductances, as the paper (the film integrals are S2's)
  Spray spray(flow, liquid, forms);
  std::vector<double> work(4 * ns + 8);

  // (b) The A-S forms: Re 0 to 400, B_M 1e-6 to 20, in O2/H2 at 1 MPa, 300 and 1,000 K; and boiling.
  const auto o2 = oxygenIndex(medium);
  std::size_t h2 = 0;
  for (std::size_t q = 0; q < ns; ++q)
    if (medium.species()[q].name == "H2") h2 = q;
  double worstB = 0, bMin = 1e9, bMax = 0, reMax = 0;
  int cases = 0;
  const char* worstWhat = "";
  for (double tGas : {300.0, 1000.0})
    for (double tDrop : {90.0, 115.0})
      for (double b : {1e-6, 1e-4, 1e-2, 0.1, 1.0, 3.0, 10.0, 20.0})
        for (double re : {0.0, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 400.0}) {
          // Y_O2 in the free stream for the B_M wanted (Raoult's Y_s at tDrop over an H2 background).
          const double x = liquid.at(tDrop).pressure / 1e6;
          const double wO2 = medium.species()[o2].molarMass, wH2 = medium.species()[h2].molarMass;
          const double ys = x * wO2 / (x * wO2 + (1 - x) * wH2);
          const double yInf = ys - b * (1 - ys);
          if (yInf < 0) continue;
          std::vector<double> y(ns, 0.0);
          y[o2] = yInf;
          y[h2] = 1 - yInf;
          Spray::Gas gas{tGas, 1e6, 1e6 / (medium.gasConstant(y.data()) * tGas), 0, 0, y.data()};
          const double d = 50e-6;
          const auto still = abramzonSirignano(medium, liquid, d, tDrop, 0, gas, false);
          gas.uz = re * still.viscosity / (gas.rho * d);
          const auto e = abramzonSirignano(medium, liquid, d, tDrop, gas.uz, gas, false);
          Parcel q;
          q.d = d;
          q.t = tDrop;
          const auto r = spray.rates(q, gas, work);
          const std::array<std::pair<double, double>, 7> pairs{{{r.reynolds, e.reynolds}, {r.spalding, e.spalding},
                                                                {r.sherwood, e.sherwood}, {r.nusselt, e.nusselt},
                                                                {r.transfer, e.transfer}, {r.mass, e.mass}, {r.heat, e.heat}}};
          const char* names[] = {"Re", "B_M", "Sh*", "Nu*", "B_T", "mdot", "Q"};
          for (std::size_t c = 0; c < pairs.size(); ++c) {
            const double err = re == 0 && c == 0 ? std::abs(pairs[c].first) : relative(pairs[c].first, pairs[c].second);
            if (err > worstB) { worstB = err; worstWhat = names[c]; }
          }
          if (r.boiling) worstB = std::max(worstB, 1.0);
          bMin = std::min(bMin, e.spalding);
          bMax = std::max(bMax, e.spalding);
          reMax = std::max(reMax, e.reynolds);
          ++cases;
        }
  std::printf("  %d evaporating cases, B_M %.3e to %.3e, Re to %.1f; largest difference in %s\n", cases, bMin, bMax, reMax, worstWhat);
  check(worstB < 1e-12, "(b) Re, B_M, Sh*, Nu*, B_T, mdot, Q against the paper's forms (evaporating)", worstB, 1e-12);
  double worstBoil = 0;
  for (double tGas : {300.0, 1000.0})
    for (double re : {0.0, 1.0, 30.0, 400.0}) {
      std::vector<double> y(ns, 0.0);
      y[h2] = 1;
      Spray::Gas gas{tGas, 1e6, 1e6 / (medium.gasConstant(y.data()) * tGas), 0, 0, y.data()};
      const double d = 50e-6, tSat = liquid.saturationTemperature(1e6);
      const auto still = abramzonSirignano(medium, liquid, d, tSat, 0, gas, true);
      gas.uz = re * still.viscosity / (gas.rho * d);
      const auto e = abramzonSirignano(medium, liquid, d, tSat, gas.uz, gas, true);
      Parcel q;
      q.d = d;
      q.t = tSat;
      const auto r = spray.rates(q, gas, work);
      if (!r.boiling) worstBoil = 1;
      for (auto [a, b] : {std::pair{r.nusselt, e.nusselt}, {r.transfer, e.transfer}, {r.mass, e.mass}, {r.heat, e.heat}})
        worstBoil = std::max(worstBoil, relative(a, b));
    }
  check(worstBoil < 1e-12, "(b) boiling at saturation: Nu*, B_T, mdot, Q against the paper's forms", worstBoil, 1e-12);

  // (c) The film's lean transport against Medium::transport: 10,000 states, T 60 to 3,500 K, p 0.1 to
  // 5 MPa, O2 with random amounts of the other species (each absent with probability 1/2), and pure O2.
  {
    std::uint64_t state = 7;
    auto uniform = [&] {  // splitmix64 on [0, 1)
      std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
      z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
      z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
      return static_cast<double>((z ^ (z >> 31)) >> 11) * 0x1.0p-53;
    };
    std::vector<double> y(ns), diffusion(ns), scratch;
    double worstT = 0;
    for (int i = 0; i < 10000; ++i) {
      double total = 0;
      for (std::size_t k = 0; k < ns; ++k) {
        y[k] = k == o2 ? uniform() : (i % 100 == 0 || uniform() < 0.5 ? 0 : uniform());
        total += y[k];
      }
      if (!(total > 0)) y[o2] = total = 1;
      for (auto& v : y) v /= total;
      const double t = 60 * std::pow(3500.0 / 60, uniform()), p = 1e5 * std::pow(50.0, uniform());
      const auto full = medium.transport(t, p, y.data(), diffusion.data(), scratch);
      const auto lean = spray.filmTransport(t, p, y.data());
      worstT = std::max({worstT, relative(lean.conductivity, full.conductivity), relative(lean.r, medium.gasConstant(y.data())),
                         diffusion[o2] == 0 ? std::abs(lean.diffusion) : relative(lean.diffusion, diffusion[o2])});
    }
    check(worstT < 1e-13, "(c) the film's lean transport against Medium::transport (lambda, D_O2, R)", worstT, 1e-13);
  }
}

// ---------------------------------------------------------------------------------------------------
// S5

// KH-RT from the published forms (Reitz 1987; Patterson and Reitz 1998; SPRAY.md, Literature), written
// out independently of Spray::waves.
Spray::Waves publishedWaves(double d, double u, double rhoG, double rhoL, double muL, double sigma, double g, double cRt) {
  const double a = 0.5 * d;
  const double weG = rhoG * u * u * a / sigma, weL = rhoL * u * u * a / sigma, reL = u * a / (muL / rhoL);
  const double ohnesorge = std::sqrt(weL) / reL, taylor = ohnesorge * std::sqrt(weG);
  Spray::Waves w;
  w.khLength = a * 9.02 * (1 + 0.45 * std::pow(ohnesorge, 0.5)) * (1 + 0.4 * std::pow(taylor, 0.7)) /
               std::pow(1 + 0.87 * std::pow(weG, 1.67), 0.6);
  w.khRate = std::sqrt(sigma / (rhoL * std::pow(a, 3))) * (0.34 + 0.38 * std::pow(weG, 1.5)) /
             ((1 + ohnesorge) * (1 + 1.4 * std::pow(taylor, 0.6)));
  w.rtLength = cRt * 2 * kPi * std::sqrt(3 * sigma / (g * (rhoL - rhoG)));
  w.rtRate = std::sqrt(2 / (3 * std::sqrt(3 * sigma)) * std::pow(g * (rhoL - rhoG), 1.5) / (rhoL + rhoG));
  return w;
}

void breakup(thermo::ReactionSource& source) {
  std::printf("S5. Breakup\n");
  double worst = 0;
  for (double d : {1e-6, 1e-5, 1e-4, 1e-3, 3e-3})
    for (double u : {1.0, 10.0, 100.0, 400.0})
      for (double rhoG : {0.5, 10.0, 60.0})
        for (double g : {10.0, 1e4, 1e8}) {
          const auto a = Spray::waves(d, u, rhoG, 1141, 1.9e-4, 0.0132, g, 0.1);
          const auto b = publishedWaves(d, u, rhoG, 1141, 1.9e-4, 0.0132, g, 0.1);
          for (auto [x, y] : {std::pair{a.khLength, b.khLength}, {a.khRate, b.khRate}, {a.rtLength, b.rtLength}, {a.rtRate, b.rtRate}})
            worst = std::max(worst, relative(x, y));
        }
  check(worst < 1e-12, "(a) KH and RT wavelengths and growth rates against the published forms", worst, 1e-12);

  // (b) A 1 mm blob at 90 K in H2 at 1 MPa and 300 K streaming at 250 m/s, one-way, steps of 0.1 us
  // (short enough that the blob takes one sub-step each, so its state at a child's creation is known).
  Definition def = chamber(source, 8, 4, 0.04, 0.01, "H2:1", 1e6, 300);
  Flow flow(def);
  const auto& medium = flow.medium();
  const auto ns = medium.size();
  std::vector<double> yH2(ns, 0.0);
  for (std::size_t k = 0; k < ns; ++k)
    if (medium.species()[k].name == "H2") yH2[k] = 1;
  const double rhoGas = 1e6 / (medium.gasConstant(yH2.data()) * 300);
  flow.setUniform({rhoGas, 250, 0, 1e6});
  const auto liquid = liquidFor(medium);
  for (int pass = 0; pass < 2; ++pass) {
    const bool khOnly = pass == 0;
    Spray::Settings settings;
    settings.twoWay = false;
    if (khOnly) settings.breakup.cTau = 1e30;  // RT never faster than KH
    Spray spray(flow, liquid, settings);
    Parcel blob;
    blob.z = 0.002;
    blob.r = 0.005;
    blob.d = 1e-3;
    blob.t = 90;
    blob.count = 1;
    spray.add(blob);
    const auto cg = cellGas(flow, static_cast<std::size_t>(spray.parcels()[0].cell));
    const double b0 = settings.breakup.b0;
    double worstChild = 0, worstMass = 0;
    long checked = 0;
    long rtEvents = 0;
    const int steps = khOnly ? 2000 : 5000;
    for (int step = 0; step < steps; ++step) {
      const auto before = spray.parcels();
      spray.advance(1e-7);
      worstMass = std::max(worstMass, std::abs(spray.balanceErrors()[0]));
      const auto& after = spray.parcels();
      if (after.empty() || before.empty()) continue;
      // The blob stays first while it lives (Spray::advance keeps the order), followed by the children
      // it released in this step. An RT event multiplies its drop count.
      if (after[0].count > before[0].count) ++rtEvents;
      if (!khOnly) continue;
      const Parcel& a = before[0];
      const Parcel& parent = after[0];
      for (std::size_t j = 1; j < after.size(); ++j) {
        const Parcel& c = after[j];
        if (!(c.z == parent.z && c.r == parent.r && c.uz == parent.uz && c.ur == parent.ur && c.ut == parent.ut &&
              c.distance == parent.distance))
          break;
        // The blob's drop at the breakup, before it shed this step: its mass now plus what this step
        // shed (the child's mass per blob drop less the shed it held before the step).
        const auto l = liquid.at(parent.t);
        const double shed = c.count * c.dropMass(l.density) / parent.count - a.shed;
        const double mass = parent.dropMass(l.density) + shed;
        const double dBefore = std::cbrt(6 * mass / (kPi * l.density));
        const double slip = std::hypot(cg.gas.uz - parent.uz, cg.gas.ur - parent.ur, parent.ut);
        const auto w = publishedWaves(dBefore, slip, cg.gas.rho, l.density, l.viscosity, l.tension, 1, 0.1);
        worstChild = std::max(worstChild, relative(c.d, 2 * b0 * w.khLength));
        ++checked;
      }
    }
    const auto& ledger = spray.ledger();
    std::printf("  %s: %zu parcels at %.1f ms, %ld children, %ld RT events of the blob, liquid %.6e of %.6e kg (exchanged %.3e kg)\n",
                khOnly ? "KH alone" : "KH-RT", spray.parcels().size(), steps * 1e-4, ledger.children, rtEvents, spray.liquidMass(),
                ledger.injectedMass, ledger.exchangedMass);
    check(worstMass < 1e-13, khOnly ? "(b) KH alone: liquid mass through every child (relative)" : "(b) KH-RT: liquid mass through every child and RT event",
          worstMass, 1e-13);
    check(khOnly ? rtEvents == 0 : rtEvents > 0, khOnly ? "(b) KH alone: no RT event" : "(b) KH-RT: the blob breaks by RT at least once",
          khOnly ? static_cast<double>(rtEvents) : (rtEvents > 0 ? 0.0 : 1.0), 0.5);
    if (khOnly) {
      std::printf("  %ld children of the blob checked against it at creation\n", checked);
      check(checked > 0 && worstChild < 1e-12, "(b) each child's diameter against 2 B0 Lambda at its creation", worstChild, 1e-12);
    }
  }
}

// ---------------------------------------------------------------------------------------------------
// S3, S4 and S6: a 16 by 6 cold chamber (radius 10 mm, length 16 mm) of H2 at 1 MPa and 300 K. H2 enters
// through the ring R/6 to R/2 at 3 g/s and 300 K; LOX at 90 K through r < 1.5 mm at 2 g/s, 50 m/s, in a
// 35 degree cone, Rosin-Rammler SMD 30 um (q 3), 2e5 parcels per second; walls restitution 0.5.

struct Run {
  std::vector<Conserved> state;
  std::vector<double> partial;
  std::vector<Parcel> parcels;
  Measurements start, end;
  Spray::Ledger ledger;
  std::array<double, 3> liquidErrors{};
  double liquidMass{}, liquidEnergy{}, liquidMomentum{}, wall{};
  long steps{};
};

// turbulent: SST with S8(c)'s inflow (intensity 0.05, viscosity ratio 10) and the Spalart-Rumsey ambient
// of tests/chamber_study.cpp, and the spray's dispersion.
Run chamberRun(thermo::ReactionSource& source, int threads, bool withSpray, double endTime, bool turbulent = false) {
  constexpr double radius = 0.01;
  Definition def = chamber(source, 16, 6, 0.016, radius, "H2:1", 1e6, 300);
  Supply supply;
  supply.innerRadius = radius / 6;
  supply.outerRadius = radius / 2;
  supply.massFlow = 3e-3;
  supply.totalTemperature = 300;
  supply.composition = source.massFractions("H2:1");
  if (turbulent) {
    supply.turbulenceIntensity = 0.05;
    supply.viscosityRatio = 10;
    def.turbulence.enabled = true;
    const Medium medium = def.medium();
    std::vector<double> diffusion(def.species.size()), work;
    const double a0 = medium.soundSpeed(300, def.composition.data());
    const double rho0 = 1e6 / (medium.gasConstant(def.composition.data()) * 300);
    const double mu0 = medium.transport(300, 1e6, def.composition.data(), diffusion.data(), work).viscosity;
    def.turbulence.ambientK = 9e-9 * a0 * a0;
    def.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
  }
  def.supplies = {supply};
  Flow flow(def);
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", threads, 1e-6, 1e-12, thermo::Chemistry::Frozen);
  Spray::Settings settings;
  Spray::Injector j;
  j.rOuter = 1.5e-3;
  j.massFlow = 2e-3;
  j.temperature = 90;
  j.velocity = 50;
  j.coneHalfAngle = 35 * kPi / 180;
  j.size = Spray::Injector::Size::RosinRammler;
  j.sauterDiameter = 30e-6;
  j.exponent = 3;
  j.parcelRate = 2e5;
  settings.injectors = {j};
  settings.restitution = 0.5;
  Spray spray(flow, liquidFor(flow.medium()), settings);
  if (withSpray) reacting.setSpray(&spray);
  Run run;
  run.start = flow.measurements();
  const auto t0 = std::chrono::steady_clock::now();
  while (flow.time() < endTime) {
    reacting.step(endTime - flow.time());
    ++run.steps;
  }
  run.wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  run.end = flow.measurements();
  run.state = flow.state();
  run.partial = flow.partialDensities();
  run.parcels = spray.parcels();
  run.ledger = spray.ledger();
  run.liquidErrors = spray.balanceErrors();
  run.liquidMass = spray.liquidMass();
  run.liquidEnergy = spray.liquidEnergy();
  run.liquidMomentum = spray.liquidMomentum();
  return run;
}

bool sameBits(const void* a, const void* b, std::size_t bytes) { return std::memcmp(a, b, bytes) == 0; }
bool identical(const Run& a, const Run& b) {
  if (a.state.size() != b.state.size() || a.partial.size() != b.partial.size() || a.parcels.size() != b.parcels.size())
    return false;
  if (!sameBits(a.state.data(), b.state.data(), a.state.size() * sizeof(Conserved))) return false;
  if (!sameBits(a.partial.data(), b.partial.data(), a.partial.size() * sizeof(double))) return false;
  for (std::size_t i = 0; i < a.parcels.size(); ++i) {
    const Parcel &p = a.parcels[i], &q = b.parcels[i];
    const std::array<double, 15> x{p.z, p.r, p.uz, p.ur, p.ut, p.d, p.t, p.count, p.shed, p.rtTime, p.distance,
                                   p.eddyZ, p.eddyR, p.eddyT, p.eddyTime};
    const std::array<double, 15> y{q.z, q.r, q.uz, q.ur, q.ut, q.d, q.t, q.count, q.shed, q.rtTime, q.distance,
                                   q.eddyZ, q.eddyR, q.eddyT, q.eddyTime};
    if (!sameBits(x.data(), y.data(), sizeof x) || p.random != q.random || p.cell != q.cell || p.injector != q.injector)
      return false;
  }
  return true;
}

// S3's checks on a run; a and b label the gas-and-liquid and the whole-system checks.
void closes(const Run& run, const char* a, const char* b) {
  const auto& g = run.end;
  const auto& l = run.ledger;
  std::printf("  %ld steps in %.1f s; %ld parcels injected, %ld children, %zu held at the end; %ld wall hits, %ld flashes, "
              "%ld vaporised, %ld clamped condensation sub-steps, %ld sub-steps\n",
              run.steps, run.wall, l.parcelsInjected, l.children, run.parcels.size(), l.wallHits, l.flashes, l.vaporised,
              l.clampedCondensation, l.substeps);
  std::printf("  liquid injected %.6e kg, evaporated %.6e kg, out %.6e kg, held %.6e kg; wall impulse %.4e N s, wall loss %.4e J\n",
              l.injectedMass, l.exchangedMass, l.outflowMass, run.liquidMass, l.wallImpulse, l.wallLoss);
  auto within = [&](double value, const char* label, const char* what) {
    char line[96];
    std::snprintf(line, sizeof line, "%s %s", label, what);
    check(std::abs(value) < 1e-11, line, std::abs(value), 1e-11);
  };
  within(g.massBalanceError, a, "gas mass balance with the exchange");
  within(g.energyBalanceError, a, "gas energy balance with the exchange");
  within(g.momentumBalanceError, a, "gas axial momentum balance with the exchange");
  within(run.liquidErrors[0], a, "liquid mass ledger");
  within(run.liquidErrors[1], a, "liquid energy ledger");
  within(run.liquidErrors[2], a, "liquid axial momentum ledger");
  // (b) Gas and liquid together, no exchange terms: (gas + liquid now) - (gas at the start + the gas's
  // boundary flux + injected - out [+ wall impulse, - wall loss]). The gas's boundary terms are its
  // balance residual's: flux = now - start - exchanged(gas) - error * scale, so the whole residual is
  // error * scale + liquid now - (injected - exchanged(gas) - out ...); it differs from the sum of the two
  // ledgers by exchanged(gas) - exchanged(liquid).
  const double m0 = run.start.mass, e0 = std::abs(run.start.energy);
  const double wm = g.massBalanceError * m0 + run.liquidMass - (l.injectedMass - g.exchangedMass - l.outflowMass);
  const double we = g.energyBalanceError * run.start.energy + run.liquidEnergy -
                    (l.injectedEnergy - g.exchangedEnergy - l.outflowEnergy - l.wallLoss);
  const double wp = g.momentumBalanceError * g.momentumBalanceScale + run.liquidMomentum -
                    (l.injectedMomentum - g.exchangedMomentum - l.outflowMomentum + l.wallImpulse);
  const double sm = m0 + l.injectedMass, se = e0 + std::abs(l.injectedEnergy) + std::abs(l.outflowEnergy) + l.wallLoss;
  const double sp = g.momentumBalanceScale + std::abs(l.injectedMomentum) + std::abs(l.outflowMomentum) + std::abs(l.wallImpulse);
  within(wm / sm, b, "whole system mass");
  within(we / se, b, "whole system energy");
  within(wp / sp, b, "whole system axial momentum");
}

void budgets(thermo::ReactionSource& source) {
  std::printf("S3. Budgets with parcels (16 by 6 cold chamber to 0.5 ms)\n");
  const Run run = chamberRun(source, 1, true, 0.5e-3);
  const auto& l = run.ledger;
  const bool happened = l.parcelsInjected > 0 && l.exchangedMass > 0 && l.wallHits > 0 && l.outflowMass > 0;
  closes(run, "(a)", "(b)");
  check(happened, "injection, evaporation, wall rebounds and parcels leaving all occurred", happened ? 0 : 1, 0.5);

  std::printf("S4. Threads\n");
  for (int threads : {4, 6}) {
    const Run other = chamberRun(source, threads, true, 0.5e-3);
    char what[96];
    std::snprintf(what, sizeof what, "final gas state and parcel list on %d threads against 1, bit for bit", threads);
    const bool same = identical(run, other);
    check(same, what, same ? 0 : 1, 0.5);
  }
}

// ---------------------------------------------------------------------------------------------------
// S8

void dispersion(thermo::ReactionSource& source) {
  std::printf("S8. Turbulent dispersion\n");
  // (b) The interaction time against Gosman and Ioannides' forms with the eddy length of SPRAY.md.
  const double sc = 0.7;
  double worstI = 0;
  int crossing = 0, lifetimes = 0;
  for (double k : {1e-6, 1e-2, 1.0, 1e2, 1e4})
    for (double omega : {1.0, 1e2, 1e4, 1e6})
      for (double tau : {1e-7, 1e-5, 1e-3, 1e-1})
        for (double slip : {0.0, 1e-3, 0.1, 10.0, 1e3}) {
          const double sigma = std::sqrt(2 * k / 3), length = 3 / sc * std::sqrt(2.0 / 3) * std::sqrt(k) / omega;
          const double life = length / sigma;
          const double cross = tau * slip > length ? -tau * std::log1p(-length / (tau * slip)) : 1e300;
          (cross < life ? crossing : lifetimes) += 1;
          worstI = std::max(worstI, relative(Spray::interaction(k, omega, tau, slip, sc), std::min(life, cross)));
        }
  std::printf("  %d cases bound by the crossing time, %d by the eddy lifetime\n", crossing, lifetimes);
  check(worstI < 1e-14 && crossing > 0 && lifetimes > 0, "(b) interaction time against the formulas (both branches)", worstI, 1e-14);

  // (a) Spreading from a point on the axis in frozen, uniform turbulence, one-way, no mass or heat passing.
  const double k = 1e-7, omega = 100, tauE = 3 / (sc * omega);
  const int n = 10000, eddies = 10;
  Definition def = chamber(source, 4, 2, 0.004, 0.002, "O2:1", 101325, 85);
  def.turbulence.enabled = true;
  def.turbulence.schmidt = sc;
  def.turbulence.ambientK = k;
  def.turbulence.ambientOmega = omega;
  Flow flow(def);
  const auto liquid = liquidFor(flow.medium());
  Spray::Settings settings;
  settings.twoWay = false;
  settings.breakup.on = false;
  Parcel p;
  p.z = 0.0015;
  p.r = 0;
  p.d = 10e-6;
  p.count = 1;
  p.t = 85;
  std::vector<double> work(4 * flow.medium().size() + 8);
  double tauP = 0;
  for (int pass = 0; pass < 3; ++pass) {  // tau_p = tau_e / 20 at zero slip (tau_p goes as d^2 there)
    Spray probe(flow, liquid, settings);
    probe.add(p);
    const auto cg = cellGas(flow, static_cast<std::size_t>(probe.parcels()[0].cell));
    tauP = probe.rates(probe.parcels()[0], cg.gas, work).tau;
    if (pass < 2) p.d *= std::sqrt(tauE / 20 / tauP);
  }
  Spray spray(flow, liquid, settings);
  for (int i = 0; i < n; ++i) spray.add(p);
  const int steps = 8 * eddies;
  for (int i = 0; i < steps; ++i) spray.advance(tauE / 8);
  double sz = 0, sz2 = 0, sr2 = 0, reMax = 0;
  for (const auto& q : spray.parcels()) {
    sz += q.z - p.z;
    sz2 += sq(q.z - p.z);
    sr2 += sq(q.r);
  }
  const auto& led = spray.ledger();
  const double mz = sz / n, mz2 = sz2 / n, mr2 = sr2 / n;
  const double sigma2 = 2 * k / 3, qd = std::exp(-tauE / tauP), N = eddies;
  const double expect = sigma2 * (N * tauE * tauE - 2 * tauP * tauE * (1 - std::pow(qd, N)) +
                                  tauP * tauP * (1 - qd) * (1 - std::pow(qd, 2 * N)) / (1 + qd));
  {
    Parcel fast = p;
    fast.uz = 3 * std::sqrt(sigma2);
    const auto cg = cellGas(flow, static_cast<std::size_t>(spray.parcels()[0].cell));
    reMax = spray.rates(fast, cg.gas, work).reynolds;
  }
  std::printf("  d %.3f um, tau_p %.4e s, tau_e %.4e s (ratio %.2f), Re at 3 sigma slip %.2e; %zu parcels left, %ld sub-steps, "
              "%ld clamped, %ld vaporised\n",
              p.d * 1e6, tauP, tauE, tauE / tauP, reMax, spray.parcels().size(), led.substeps, led.clampedCondensation,
              led.vaporised);
  std::printf("  <dz> %.4e m (3 se %.4e); <dz^2> %.6e m^2 against %.6e (%+.3f%%); <r^2> %.6e against %.6e (%+.3f%%)\n", mz,
              3 * std::sqrt(mz2 / n), mz2, expect, 100 * (mz2 / expect - 1), mr2, 2 * expect, 100 * (mr2 / (2 * expect) - 1));
  const bool all = spray.parcels().size() == static_cast<std::size_t>(n) && led.exchangedMass == 0;
  check(all, "(a) every parcel kept, no mass passed", all ? 0 : 1, 0.5);
  check(std::abs(mz) < 3 * std::sqrt(mz2 / n), "(a) <dz> against 0, in standard errors", std::abs(mz) / std::sqrt(mz2 / n), 3);
  check(std::abs(mz2 / expect - 1) < 3 * std::sqrt(2.0 / n), "(a) <dz^2> against the held-velocity expectation, in s.e.",
        std::abs(mz2 / expect - 1) / std::sqrt(2.0 / n), 3);
  check(std::abs(mr2 / (2 * expect) - 1) < 3 / std::sqrt(n), "(a) <r^2> against twice it, in standard errors",
        std::abs(mr2 / (2 * expect) - 1) * std::sqrt(n), 3);

  // (c) S3's chamber with SST and dispersion: budgets and threads.
  const Run run = chamberRun(source, 1, true, 0.2e-3, true);
  long eddying = 0;
  for (const auto& q : run.parcels) eddying += q.eddyTime > 0;
  std::printf("  (c) %ld of the parcels held are inside an eddy at the end\n", eddying);
  check(eddying > 0, "(c) dispersion acted (parcels inside an eddy)", eddying > 0 ? 0 : 1, 0.5);
  closes(run, "(c)", "(c)");
  const Run other = chamberRun(source, 4, true, 0.2e-3, true);
  const bool same = identical(run, other);
  check(same, "(c) final gas state and parcel list on 4 threads against 1, bit for bit", same ? 0 : 1, 0.5);
}

void cost(thermo::ReactionSource& source, int threads) {
  std::printf("S6. Cost of the S3 run on %d threads (measured)\n", threads);
  const Run gas = chamberRun(source, threads, false, 0.5e-3);
  const Run both = chamberRun(source, threads, true, 0.5e-3);
  const double spray = both.wall - gas.wall;
  std::printf("  gas alone %.3f s (%ld steps), with the spray %.3f s (%ld steps): spray %.3f s, %.1f%% of the step\n", gas.wall,
              gas.steps, both.wall, both.steps, spray, 100 * spray / both.wall);
  std::printf("  %ld parcel sub-steps: %.3e per second of spray time\n", both.ledger.substeps, both.ledger.substeps / spray);
}
}  // namespace

// S2(c), the model's side: Lafon et al.'s (2008) drop, 100 um at 90 K in quiescent H2 at 1,000 K and
// 10 atm, one-way, to vaporisation. The lifetime, and -d(d^2)/dt fitted by least squares over the
// samples with (d/d0)^2 from 0.6 down to 0.1 (after heat-up), as stated in SPRAY.md before their
// figures were read.
void lafon(thermo::ReactionSource& source, double dt) {
  Flow flow(chamber(source, 4, 2, 0.01, 0.005, "H2:1", 101325 * 10, 1000));
  const auto liquid = liquidFor(flow.medium());
  Spray::Settings settings;
  settings.twoWay = false;
  settings.breakup.on = false;
  Spray spray(flow, liquid, settings);
  Parcel p;
  p.z = 0.005;
  p.r = 0.001;
  p.d = 100e-6;
  p.count = 1;
  p.t = 90;
  spray.add(p);
  const double d0 = p.d;
  double t = 0, sT = 0, sY = 0, sTT = 0, sTY = 0, peak = 1, prev = 1, tPrev = 0;
  int samples = 0;
  std::vector<std::array<double, 2>> history;  // t, T
  std::printf("  t [ms]  (d/d0)^2  T [K]\n");
  while (!spray.parcels().empty()) {
    const auto& q = spray.parcels()[0];
    const double s = q.d * q.d / (d0 * d0);
    peak = std::max(peak, s);
    history.push_back({t, q.t});
    if (s <= 0.6 && s >= 0.1) { sT += t; sY += s; sTT += t * t; sTY += t * s; ++samples; }
    if (static_cast<long>(std::lround(t / dt)) % std::lround(5e-4 / dt) == 0) std::printf("  %6.2f  %.6f  %.4f\n", t * 1e3, s, q.t);
    prev = s;
    tPrev = t;
    spray.advance(dt);
    t += dt;
  }
  // The end within one step: d^2 extrapolated linearly from the last sample to d_min.
  const double dMin2 = settings.minDiameter * settings.minDiameter / (d0 * d0);
  const double slope = (samples * sTY - sT * sY) / (samples * sTT - sT * sT);  // of (d/d0)^2 per s
  const double life = tPrev + (prev - dMin2) / -slope;
  // Heat-up: the time T first reaches 99% of its rise to the final (wet-bulb) temperature.
  const double tFinal = history.back()[1];
  double t99 = 0;
  for (const auto& h : history)
    if (h[1] >= 90 + 0.99 * (tFinal - 90)) { t99 = h[0]; break; }
  std::printf("  lifetime %.4f ms (last step %.4f to %.4f ms); -d(d^2)/dt %.5f mm^2/s over (d/d0)^2 0.6 to 0.1 (%d samples);"
              " peak (d/d0)^2 %.5f; T reaches 99%% of its rise to %.4f K at %.4f ms (T_sat %.4f K)\n",
              life * 1e3, tPrev * 1e3, t * 1e3, -slope * d0 * d0 * 1e6, samples, peak, tFinal, t99 * 1e3,
              liquid.saturationTemperature(101325 * 10));
}

int main(int argc, char** argv) {
  try {
    thermo::ReactionSource source("h2o2.yaml");
    if (argc > 1 && std::string(argv[1]) == "cost") {
      cost(source, argc > 2 ? std::atoi(argv[2]) : 1);
      return 0;
    }
    if (argc > 1 && std::string(argv[1]) == "lafon") {
      lafon(source, argc > 2 ? std::atof(argv[2]) : 2e-6);
      return 0;
    }
    {
      Flow flow(chamber(source, 4, 2, 0.01, 0.005, "H2:1", 1e6, 300));
      properties(flow.medium(), liquidFor(flow.medium()));
    }
    dropModel(source);
    breakup(source);
    budgets(source);
    dispersion(source);
  } catch (const std::exception& e) {
    std::printf("error: %s\n", e.what());
    return 2;
  }
  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
