// Diagnostics for the PaSR test failures (not a test): criterion 2's reference gaps, and 3(d) field by field.
#define main pasr_main
#include "tests/pasr_tests.cpp"
#undef main
#include <chrono>

// RK4 as in the test, with Kahan-compensated accumulation of z.
std::vector<std::vector<double>> rk4k(thermo::ReactionSource& source, double rho, std::vector<double> z, double dt,
                                      long perSample, const std::vector<std::size_t>& set, double tauMix, double s) {
  const std::size_t n = z.size();
  std::vector<double> k1(n), k2(n), k3(n), k4(n), w(n), c(n, 0.0);
  auto f = [&](const std::vector<double>& x, std::vector<double>& out) {
    source.rates(rho, x.data(), out.data());
    const double factor = kappaEff(x.data(), out.data(), set, tauMix, s);
    for (double& v : out) v *= factor;
  };
  const double h = dt / (20.0 * double(perSample));
  std::vector<std::vector<double>> samples;
  try {
  for (int sample = 0; sample < 20; ++sample) {
    for (long step = 0; step < perSample; ++step) {
      f(z, k1);
      for (std::size_t i = 0; i < n; ++i) w[i] = z[i] + 0.5 * h * k1[i];
      f(w, k2);
      for (std::size_t i = 0; i < n; ++i) w[i] = z[i] + 0.5 * h * k2[i];
      f(w, k3);
      for (std::size_t i = 0; i < n; ++i) w[i] = z[i] + h * k3[i];
      f(w, k4);
      for (std::size_t i = 0; i < n; ++i) {
        const double inc = h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) - c[i];
        const double t = z[i] + inc;
        c[i] = (t - z[i]) - inc;
        z[i] = t;
      }
      if (!std::isfinite(z[0]) || z[0] <= 0 || z[0] > 1e4) { std::printf("    (failed in sample %d)\n", sample); return {}; }
    }
    samples.push_back(z);
  }
  } catch (const std::exception&) { std::printf("    (threw)\n"); return {}; }
  return samples;
}

int main(int argc, char** argv) {
  const std::string mechanism = "h2o2.yaml";
  thermo::ReactionSource source(mechanism);
  Medium medium = source.medium();
  const auto fits = thermo::transportFits(mechanism);
  medium.setTransport(fits);
  const std::size_t ns = medium.size();
  auto index = [&](const char* name) {
    for (std::size_t k = 0; k < ns; ++k) if (medium.species()[k].name == name) return k;
    throw std::runtime_error(std::string("missing species ") + name);
  };
  const std::size_t iH2 = index("H2"), iO2 = index("O2"), iH2O = index("H2O"), iN2 = index("N2"), iAR = index("AR");
  const std::vector<std::size_t> set = {iH2, iO2, iH2O};
  const std::vector<std::string> setNames = {"H2", "O2", "H2O"};
  const std::string which = argc > 1 ? argv[1] : "both";

  if (which != "ref") {
    std::printf("3(d) diagnosis\n");
    auto inertMix = [&](double z, double* x) { x[iAR] = 0.1 + 20 * z; x[iN2] = 0.9 - 20 * z; };
    for (int variant = 0; variant < 3; ++variant) {
      Duct d = duct(medium, fits, true, inertMix);
      Flow a(d.definition), b(d.definition);
      for (Flow* f : {&a, &b}) {
        f->setInitialState(d.cells, d.y);
        f->setTurbulence(uniformTurbulence(d.cells.size(), 10, 1e3));
      }
      double taken = 0;
      if (variant < 2) {
        thermo::ReactingFlow reacting(a, mechanism, 2, 1e-10, 1e-16, thermo::Chemistry::FiniteRate);
        if (variant == 0) reacting.setMixingClosure(1, setNames);
        taken = reacting.step();
      } else {
        std::vector<double> t, s;
        a.mixingInputs(1, set, t, s);
        taken = a.step();
      }
      const double tb = b.step(taken);
      std::printf(" variant %d (%s): taken %.6e, b took %.6e, times %.6e %.6e\n", variant,
                  variant == 0 ? "ReactingFlow + closure" : variant == 1 ? "ReactingFlow, no closure" : "Flow after mixingInputs",
                  taken, tb, a.time(), b.time());
      auto report = [&](const char* name, const std::vector<double>& x, const std::vector<double>& y, std::size_t stride) {
        for (std::size_t f = 0; f < stride; ++f) {
          double scale = 0, gap = 0; std::size_t at = 0;
          for (std::size_t q = 0; q < x.size() / stride; ++q) {
            scale = std::max(scale, std::abs(y[q * stride + f]));
            const double g = std::abs(x[q * stride + f] - y[q * stride + f]);
            if (g > gap) { gap = g; at = q; }
          }
          const double rel = scale > 0 ? gap / scale : (gap > 0 ? 1.0 : 0.0);
          if (rel > 1e-14) std::printf("   %s[%zu] rel %.3e at cell %zu: a %.10e b %.10e (scale %.3e)\n", name, f, rel, at,
                                       x[at * stride + f], y[at * stride + f], scale);
        }
      };
      std::vector<double> sa, sb;
      for (const auto& u : a.state()) sa.insert(sa.end(), u.begin(), u.end());
      for (const auto& u : b.state()) sb.insert(sb.end(), u.begin(), u.end());
      report("state", sa, sb, a.state()[0].size());
      report("partial", a.partialDensities(), b.partialDensities(), ns);
      report("turb", a.turbulence(), b.turbulence(), 2);
    }
  }
  if (which != "duct") {
    std::printf("criterion 2 reference gaps, laminar (s = 0), plain and Kahan RK4\n");
    std::vector<double> unburnt(ns + 1);
    const auto y = source.massFractions("H2:2, O2:1");
    unburnt[0] = 1200;
    std::copy(y.begin(), y.end(), unburnt.begin() + 1);
    const double rho = source.density(1200, p0, unburnt.data() + 1);
    const int kahan = argc > 2 ? std::atoi(argv[2]) : 0;
    const long lo = argc > 3 ? std::atol(argv[3]) : 256, hi = argc > 4 ? std::atol(argv[4]) : (1L << 15);
    {
      std::vector<std::vector<double>> coarse;
      for (long per = lo; per <= hi; per *= 2) {
        const auto t0 = std::chrono::steady_clock::now();
        auto fine = kahan ? rk4k(source, rho, unburnt, 3e-4, per, set, 1, 0) : rk4(source, rho, unburnt, 3e-4, per, set, 1, 0);
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        if (fine.empty()) { std::printf("  %s per %6ld: failed\n", kahan ? "kahan" : "plain", per); coarse.clear(); continue; }
        if (!coarse.empty()) {
          const auto d = difference(coarse, fine);
          std::size_t worst = 0; double wt = 0;
          for (std::size_t s = 0; s < 20; ++s) { const double g = std::abs(coarse[s][0] - fine[s][0]) / fine[s][0]; if (g > wt) { wt = g; worst = s; } }
          std::printf("  %s per %6ld: dT %.3e (sample %zu, T %.2f) dY %.3e; T end %.6f; %.1f s\n", kahan ? "kahan" : "plain", per,
                      d.t, worst, fine[worst][0], d.y, fine.back()[0], wall);
        } else std::printf("  %s per %6ld: T end %.6f; %.1f s\n", kahan ? "kahan" : "plain", per, fine.back()[0], wall);
        coarse = std::move(fine);
      }
    }
  }
}
