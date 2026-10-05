// Heat-release duration and timing of the criterion-2 substep for several closures (diagnostic).
#define main pasr_main
#include "tests/pasr_tests.cpp"
#undef main
int main() {
  thermo::ReactionSource source("h2o2.yaml");
  Medium medium = source.medium();
  const std::size_t ns = medium.size();
  auto index = [&](const char* name) { for (std::size_t k = 0; k < ns; ++k) if (medium.species()[k].name == name) return k; return ns; };
  const std::vector<std::size_t> set = {index("H2"), index("O2"), index("H2O")};
  std::vector<double> z(ns + 1);
  const auto y = source.massFractions("H2:2, O2:1");
  z[0] = 1200; std::copy(y.begin(), y.end(), z.begin() + 1);
  const double rho = source.density(1200, p0, z.data() + 1);
  std::vector<double> rate(ns + 1);
  for (const auto& [tm, s] : std::vector<std::pair<double, double>>{{1e-5, 0}, {1e-5, 0.5}, {1e-4, 0.5}, {1e-5, 0.9}, {1e-5, 0.99}, {1e-5, 1}}) {
    const int n = 30000;  // samples 1e-8 s apart
    std::vector<double> t(n), T(n), kap(n);
    int i = 0;
    thermo::ReactionStep step(source, 1e-10, 1e-16);
    std::vector<double> w = z;
    step.advance(rho, w.data(), 3e-4, {tm, s}, set, n, [&](double time, const double* x) {
      t[i] = time; T[i] = x[0];
      source.rates(rho, x, rate.data());
      kap[i] = thermo::reactingFraction(x, rate.data(), set, {tm, s});
      ++i;
    });
    const double Tend = T[n - 1], lo = 1200 + 0.1 * (Tend - 1200), hi = 1200 + 0.9 * (Tend - 1200);
    double t10 = -1, t90 = -1, kmin = 1;
    for (int j = 0; j < n; ++j) { if (t10 < 0 && T[j] >= lo) t10 = t[j]; if (t90 < 0 && T[j] >= hi) t90 = t[j]; kmin = std::min(kmin, kap[j]); }
    std::printf("tau_mix %.0e s, s %.2f: T 10%% at %.4e s, 90%% at %.4e s, 10-90%% rise %.3e s; min kappa_eff %.3e\n", tm, s, t10, t90, t90 - t10, kmin);
  }
}
