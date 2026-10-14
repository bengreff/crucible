// The medium's property functions against their previous implementation (docs/evidence/TRANSPORT_COST.md,
// criteria T1 and T4). The previous code (core/medium.cpp as of 9809e26) is reproduced here verbatim as the
// reference: the same formulas with the divisions inside the loops. On random states it checks
//   T1: mu, lambda and every D_km within 1e-13 relative; e and cv within 1e-14 of cv T and cv; the
//       enthalpies within 1e-14 relative; the e-to-T inversion within 1e-12 relative;
//   T4: the cost per call of the old and new functions on the same states (reported, no threshold).
// Usage: crucible_medium_cost [states]
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <random>
#include <vector>

#include "adapters/reaction.hpp"
#include "core/medium.hpp"

using namespace crucible;

namespace {
// ---- The previous implementation (reference) ----
const std::array<double, 7>& range(const Species& s, double t) { return t < s.tMid ? s.low : s.high; }
double cpOverR(const std::array<double, 7>& a, double t) { return a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4]))); }
double hOverR(const std::array<double, 7>& a, double t) {
  return t * (a[0] + t * (a[1] / 2 + t * (a[2] / 3 + t * (a[3] / 4 + t * a[4] / 5)))) + a[5];
}
struct Reference {
  std::vector<Species> species;
  std::vector<double> specificR;
  TransportFits fits;
  std::vector<double> wilkeMass, wilkeRoot;
  Reference(const Medium& m, TransportFits f) : species(m.species()), fits(std::move(f)) {
    const std::size_t n = species.size();
    for (const auto& s : species) specificR.push_back(Medium::universalGasConstant / s.molarMass);
    wilkeMass.resize(n * n);
    wilkeRoot.resize(n * n);
    for (std::size_t k = 0; k < n; ++k)
      for (std::size_t j = 0; j < n; ++j) {
        double wk = species[k].molarMass, wj = species[j].molarMass;
        wilkeMass[k * n + j] = std::sqrt(std::sqrt(wj / wk));
        wilkeRoot[k * n + j] = std::sqrt(8.0) * std::sqrt(1 + wk / wj);
      }
  }
  Medium::Properties properties(double t, const double* y) const {
    Medium::Properties out;
    for (std::size_t k = 0; k < species.size(); ++k) {
      const auto& a = range(species[k], t);
      double w = y[k] * specificR[k];
      out.r += w;
      out.e += w * (hOverR(a, t) - t);
      out.cv += w * (cpOverR(a, t) - 1);
    }
    return out;
  }
  double enthalpy(double t, const double* y) const {
    double h = 0;
    for (std::size_t k = 0; k < species.size(); ++k) h += y[k] * specificR[k] * hOverR(range(species[k], t), t);
    return h;
  }
  double energyFloor(const double* y) const {
    double e = 0;
    for (std::size_t k = 0; k < species.size(); ++k) e += y[k] * specificR[k] * species[k].low[5];
    return e;
  }
  double temperature(double e, const double* y, double guess) const {
    if (!(e > energyFloor(y)) || !std::isfinite(e)) return std::numeric_limits<double>::quiet_NaN();
    double lo = 0, hi = std::numeric_limits<double>::infinity();
    double t = guess > 0 && std::isfinite(guess) ? guess : 300;
    for (int n = 0; n < 100; ++n) {
      auto q = properties(t, y);
      double f = q.e - e;
      double next = t - f / q.cv;
      if (std::abs(next - t) <= 1e-13 * t) return next;
      if (f > 0) hi = t; else lo = t;
      if (!(next > lo && next < hi)) next = std::isfinite(hi) ? 0.5 * (lo + hi) : 2 * t;
      t = next;
    }
    return std::numeric_limits<double>::quiet_NaN();
  }
  Medium::Transport transport(double t, double p, const double* y, double* diffusion, std::vector<double>& work) const {
    const std::size_t n = species.size();
    work.resize(n * (4 + n));
    double *x = work.data(), *mu = x + n, *root = mu + n, *lambda = root + n, *binary = lambda + n;
    const double logt = std::log(t), sqrtT = std::sqrt(t), t14 = std::sqrt(sqrtT);
    auto poly = [&](const std::array<double, 5>& c) { return c[0] + logt * (c[1] + logt * (c[2] + logt * (c[3] + logt * c[4]))); };
    double moles = 0;
    for (std::size_t k = 0; k < n; ++k) moles += y[k] / species[k].molarMass;
    const double meanMass = 1 / moles;
    for (std::size_t k = 0; k < n; ++k) {
      x[k] = std::max(1e-20, y[k] / species[k].molarMass * meanMass);
      root[k] = t14 * poly(fits.viscosity[k]);
      mu[k] = root[k] * root[k];
      lambda[k] = sqrtT * poly(fits.conductivity[k]);
    }
    Medium::Transport out;
    double series = 0, parallel = 0;
    for (std::size_t k = 0; k < n; ++k) {
      double weight = 0;
      for (std::size_t j = 0; j < n; ++j) {
        double f = 1 + root[k] / root[j] * wilkeMass[k * n + j];
        weight += f * f / wilkeRoot[k * n + j] * x[j];
      }
      out.viscosity += x[k] * mu[k] / weight;
      series += x[k] * lambda[k];
      parallel += x[k] / lambda[k];
    }
    out.conductivity = 0.5 * (series + 1 / parallel);
    for (std::size_t k = 0, c = 0; k < n; ++k)
      for (std::size_t j = k; j < n; ++j, ++c) binary[k * n + j] = binary[j * n + k] = t * sqrtT * poly(fits.diffusion[c]);
    if (n == 1) { diffusion[0] = binary[0] / p; return out; }
    for (std::size_t k = 0; k < n; ++k) {
      double sum = 0, others = 0;
      for (std::size_t j = 0; j < n; ++j)
        if (j != k) { sum += x[j] / binary[j * n + k]; others += std::max(y[j], 0.0); }
      diffusion[k] = sum <= 0 ? binary[k * n + k] / p : others / (p * sum);
    }
    return out;
  }
};
double relative(double a, double b) { return a == b ? 0 : std::abs(a - b) / std::max(std::abs(a), std::abs(b)); }
}  // namespace

int main(int argc, char** argv) {
  const int states = argc > 1 ? std::atoi(argv[1]) : 100000;
  thermo::ReactionSource source("h2o2.yaml");
  Medium medium = source.medium();
  const auto fits = thermo::transportFits("h2o2.yaml");
  medium.setTransport(fits);
  const Reference ref(medium, fits);
  const std::size_t n = medium.size();

  // States: T uniform in [200, 6000] K, ln p uniform in [ln 1e3, ln 2e7]; compositions a third each
  // uniform-normalised, near-pure (one species 1 - 1e-12, the rest sharing 1e-12) and with a random
  // half of the species exactly zero.
  std::mt19937_64 rng(20261013);
  std::uniform_real_distribution<double> u(0, 1);
  std::vector<double> T(states), P(states), Y(states * n);
  for (int s = 0; s < states; ++s) {
    T[s] = 200 + 5800 * u(rng);
    P[s] = std::exp(std::log(1e3) + (std::log(2e7) - std::log(1e3)) * u(rng));
    double* y = &Y[s * n];
    double sum = 0;
    const int kind = s % 3;
    const std::size_t pure = static_cast<std::size_t>(u(rng) * n) % n;
    for (std::size_t k = 0; k < n; ++k) {
      y[k] = kind == 1 ? (k == pure ? 0 : u(rng)) : (kind == 2 && u(rng) < 0.5 ? 0 : u(rng));
      sum += y[k];
    }
    if (kind == 1) {
      for (std::size_t k = 0; k < n; ++k) y[k] = k == pure ? 1 - 1e-12 : 1e-12 * y[k] / sum;
    } else if (sum > 0) {
      for (std::size_t k = 0; k < n; ++k) y[k] /= sum;
    } else {
      y[pure] = 1;
    }
  }

  // T1
  double muErr = 0, lambdaErr = 0, dErr = 0, eErr = 0, cvErr = 0, hErr = 0, hkErr = 0, tErr = 0;
  std::vector<double> d0(n), d1(n), h1(n), w0, w1;
  for (int s = 0; s < states; ++s) {
    const double* y = &Y[s * n];
    auto a = ref.transport(T[s], P[s], y, d0.data(), w0);
    auto b = medium.transport(T[s], P[s], y, d1.data(), w1);
    muErr = std::max(muErr, relative(a.viscosity, b.viscosity));
    lambdaErr = std::max(lambdaErr, relative(a.conductivity, b.conductivity));
    for (std::size_t k = 0; k < n; ++k) dErr = std::max(dErr, relative(d0[k], d1[k]));
    auto pa = ref.properties(T[s], y);
    auto pb = medium.properties(T[s], y);
    eErr = std::max(eErr, std::abs(pa.e - pb.e) / (pa.cv * T[s]));
    cvErr = std::max(cvErr, relative(pa.cv, pb.cv));
    hErr = std::max(hErr, std::abs(ref.enthalpy(T[s], y) - medium.enthalpy(T[s], y)) / ((pa.cv + pa.r) * T[s]));
    medium.speciesEnthalpies(T[s], h1.data());
    for (std::size_t k = 0; k < n; ++k) {
      const double hk = ref.specificR[k] * hOverR(range(ref.species[k], T[s]), T[s]);
      const double scale = ref.specificR[k] * cpOverR(range(ref.species[k], T[s]), T[s]) * T[s];
      hkErr = std::max(hkErr, std::abs(hk - h1[k]) / scale);
    }
    // Inversion from a guess 30% off (several Newton steps).
    const double ta = ref.temperature(pa.e, y, 0.7 * T[s]), tb = medium.temperature(pa.e, y, 0.7 * T[s]);
    tErr = std::max(tErr, relative(ta, tb));
  }
  std::printf("T1 over %d states (%zu species):\n", states, n);
  std::printf("  mu %.2e  lambda %.2e  D_km %.2e   (relative, bound 1e-13)  %s\n", muErr, lambdaErr, dErr,
              muErr <= 1e-13 && lambdaErr <= 1e-13 && dErr <= 1e-13 ? "PASS" : "FAIL");
  std::printf("  e %.2e (of cv T)  cv %.2e  h %.2e (of cp T)  h_k %.2e (of cp_k T)   (bound 1e-14)  %s\n", eErr, cvErr, hErr,
              hkErr, eErr <= 1e-14 && cvErr <= 1e-14 && hErr <= 1e-14 && hkErr <= 1e-14 ? "PASS" : "FAIL");
  std::printf("  T from e %.2e   (relative, bound 1e-12)  %s\n", tErr, tErr <= 1e-12 ? "PASS" : "FAIL");
  const bool pass = muErr <= 1e-13 && lambdaErr <= 1e-13 && dErr <= 1e-13 && eErr <= 1e-14 && cvErr <= 1e-14 &&
                    hErr <= 1e-14 && hkErr <= 1e-14 && tErr <= 1e-12;

  // T4: cost per call, best of 5 sweeps over the states.
  auto timeIt = [&](auto&& body) {
    double best = 1e30;
    for (int rep = 0; rep < 5; ++rep) {
      auto t0 = std::chrono::steady_clock::now();
      double sink = 0;
      for (int s = 0; s < states; ++s) sink += body(s);
      double ns = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() * 1e9 / states;
      if (!(sink == 12345.678)) best = std::min(best, ns);
    }
    return best;
  };
  const double trOld = timeIt([&](int s) { return ref.transport(T[s], P[s], &Y[s * n], d0.data(), w0).viscosity + d0[0]; });
  const double trNew = timeIt([&](int s) { return medium.transport(T[s], P[s], &Y[s * n], d1.data(), w1).viscosity + d1[0]; });
  const double prOld = timeIt([&](int s) { return ref.properties(T[s], &Y[s * n]).e; });
  const double prNew = timeIt([&](int s) { return medium.properties(T[s], &Y[s * n]).e; });
  std::vector<double> E(states);
  for (int s = 0; s < states; ++s) E[s] = ref.properties(T[s], &Y[s * n]).e;
  const double tOld = timeIt([&](int s) { return ref.temperature(E[s], &Y[s * n], 0.99 * T[s]); });
  const double tNew = timeIt([&](int s) { return medium.temperature(E[s], &Y[s * n], 0.99 * T[s]); });
  std::printf("T4 cost per call (ns, best of 5 sweeps): transport %.1f -> %.1f (%.2fx), properties %.1f -> %.1f (%.2fx), "
              "temperature from a 1%% guess %.1f -> %.1f (%.2fx)\n",
              trOld, trNew, trOld / trNew, prOld, prNew, prOld / prNew, tOld, tNew, tOld / tNew);
  return pass ? 0 : 1;
}
