// Item-3 evidence (docs/LOW_MACH.md): steady subsonic venturi, outlet mass flow against the isentropic
// quasi-1D value, for each low-Mach treatment, four unchoked exit Mach numbers and the given grids.
// Started from the quasi-1D solution so the comparison measures steady discretization error, not the
// startup transient (from-rest startup is covered by the core venturi test). Drift between 0.75 t and t
// shows whether the run has settled.
// Usage: crucible_low_mach_study [tEnd=0.04] [grids=40,80]
#include <cmath>
#include <cstdio>
#include <chrono>
#include <cstdlib>
#include <future>
#include <numbers>
#include <string>
#include <vector>

#include "core/flow.hpp"

using namespace crucible;

int main(int argc, char** argv) {
  const double tEnd = argc > 1 ? std::atof(argv[1]) : 0.04;
  std::vector<int> grids = {40, 80};
  if (argc > 2) { grids.clear(); for (char* s = argv[2]; *s;) { grids.push_back(std::strtol(s, &s, 10)); if (*s) ++s; } }
  const std::vector<std::pair<LowMach, const char*>> schemes = {
      {LowMach::None, "hllc"}, {LowMach::Thornber, "thornber"}, {LowMach::HllcLm, "hllc-lm"}};
  const std::vector<double> backRatios = {0.998, 0.99, 0.985, 0.98};
  struct Row { std::string scheme; double back, me; int nz; double err1, err2, imbalance; long steps; double seconds; };
  std::vector<std::future<Row>> jobs;
  for (auto [scheme, name] : schemes) for (double back : backRatios) for (int nz : grids)
    jobs.push_back(std::async(std::launch::deferred, [=] {
      Definition d; d.nz = nz; d.nr = nz * 3 / 20; d.backPressure = back * d.totalPressure; d.lowMach = scheme;
      const double g = d.gas.gamma;
      const double me = std::sqrt(2 / (g - 1) * (std::pow(1 / back, (g - 1) / g) - 1));
      const double te = d.totalTemperature / (1 + (g - 1) / 2 * me * me);
      const double ideal = d.backPressure / (d.gas.specificR * te) * me * std::sqrt(g * d.gas.specificR * te) *
                           std::numbers::pi * d.exitRadius * d.exitRadius;
      Flow f(d);
      auto areaMach = [g](double m) {
        return std::pow(2 / (g + 1) * (1 + (g - 1) / 2 * m * m), (g + 1) / (2 * (g - 1))) / m;
      };
      const double aStar = std::numbers::pi * d.exitRadius * d.exitRadius / areaMach(me);
      std::vector<Primitive> init;
      for (int i = 0; i < nz; ++i) {
        const double rc = 0.5 * (f.mesh().radius[i] + f.mesh().radius[i + 1]);
        const double ratio = std::numbers::pi * rc * rc / aStar;
        double lo = 1e-6, hi = 1;
        for (int n = 0; n < 100; ++n) { const double mm = 0.5 * (lo + hi); (areaMach(mm) > ratio ? lo : hi) = mm; }
        const double mach = 0.5 * (lo + hi), fac = 1 + (g - 1) / 2 * mach * mach;
        const double t = d.totalTemperature / fac, pr = d.totalPressure * std::pow(fac, -g / (g - 1));
        for (int j = 0; j < d.nr; ++j)
          init.push_back({pr / (d.gas.specificR * t), mach * std::sqrt(g * d.gas.specificR * t), 0, pr});
      }
      f.setInitialState(init);
      const auto t0 = std::chrono::steady_clock::now();
      f.advanceTo(0.75 * tEnd);
      const double err1 = f.measurements().outletMassFlow / ideal - 1;
      f.advanceTo(tEnd);
      const auto m = f.measurements();
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      return Row{name, back, me, nz, err1, m.outletMassFlow / ideal - 1,
                 m.inletMassFlow / m.outletMassFlow - 1, static_cast<long>(m.steps), secs};
    }));
  // Run three at a time (shared machine).
  std::vector<std::future<Row>> running;
  std::printf("scheme    p_b/p0  M_exit  nz  err(0.75t)  err(t)    inlet/outlet-1  steps  s\n");
  for (std::size_t k = 0; k < jobs.size(); k += 3) {
    std::vector<std::future<Row>> batch;
    for (std::size_t q = k; q < std::min(k + 3, jobs.size()); ++q)
      batch.push_back(std::async(std::launch::async, [&jobs, q] { return jobs[q].get(); }));
    for (auto& b : batch) {
      const Row r = b.get();
      std::printf("%-9s %.3f   %.4f  %3d  %+9.5f  %+9.5f  %+11.2e  %6ld  %.0f\n", r.scheme.c_str(), r.back, r.me,
                  r.nz, r.err1, r.err2, r.imbalance, r.steps, r.seconds);
      std::fflush(stdout);
    }
  }
}
