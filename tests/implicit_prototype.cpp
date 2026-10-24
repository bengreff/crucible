// The implicit prototype's case P (docs/evidence/IMPLICIT_PROTOTYPE.md, criteria written before this file): the
// simple-wave pulse of tests/amr_subcycling_tests.cpp 3(a) on one level, one thread.
// Criterion 1: the L2 error of p against the explicit reference at CFL 0.1 (itself checked against CFL 0.05), at
// three physical steps halving from an acoustic CFL of 4, gives an observed order in [1.8, 2.2] on both pairs.
// Criterion 2 (printed here for P; judged on C): the mass and energy budgets close within 1e-11.
// Usage: crucible_implicit_prototype [columns (192)] [inner tolerance (1e-11)] [pseudo CFL (1e6)]
#include "core/flow.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace crucible;

namespace {
constexpr double gaussX[6] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
                              0.2386191860831969,  0.6612093864662645,  0.9324695142031521};
constexpr double gaussW[6] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
                              0.4679139345726910, 0.3607615730481386, 0.1713244923791704};
// As amr_subcycling_tests.cpp 3(a): isentropic, u - 2a / (gamma - 1) = -2 a0 / (gamma - 1).
struct Pulse {
  double p0, rho0, a0, gamma, amplitude, centre, halfWidth;
  std::array<double, 3> at(double z) const {
    const double s = (z - centre) / halfWidth, p = p0 * (1 + amplitude * std::exp(-std::log(2.0) * s * s));
    const double rho = rho0 * std::pow(p / p0, 1 / gamma), a = std::sqrt(gamma * p / rho), u = 2 * (a - a0) / (gamma - 1);
    return {rho, rho * u, p / (gamma - 1) + 0.5 * rho * u * u};
  }
  void fill(Flow& f) const {
    const Mesh& m = f.mesh();
    std::vector<Primitive> cells(m.cells.size());
    for (int i = 0; i < m.nz; ++i) {
      const double z0 = i * m.dz, z1 = (i + 1) * m.dz;
      std::array<double, 3> mean{};
      for (int k = 0; k < 6; ++k) {
        const auto u = at(0.5 * (z0 + z1) + 0.5 * (z1 - z0) * gaussX[k]);
        for (int c = 0; c < 3; ++c) mean[c] += 0.5 * gaussW[k] * u[c];
      }
      const Primitive w{mean[0], mean[1] / mean[0], 0, (gamma - 1) * (mean[2] - 0.5 * mean[1] * mean[1] / mean[0])};
      for (int j = 0; j < m.nr; ++j) cells[m.index(i, j)] = w;
    }
    f.setInitialState(cells);
  }
};
double seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
// The volume-weighted RMS of p - p_ref over the pulse's amplitude.
double pressureError(const Flow& f, const std::vector<double>& reference, double scale) {
  double s = 0, v = 0;
  for (std::size_t q = 0; q < f.state().size(); ++q) {
    const double vol = f.mesh().cells[q].volume;
    s += vol * std::pow(f.cellPrimitive(q).p - reference[q], 2);
    v += vol;
  }
  return std::sqrt(s / v) / scale;
}
std::vector<double> pressures(const Flow& f) {
  std::vector<double> p(f.state().size());
  for (std::size_t q = 0; q < p.size(); ++q) p[q] = f.cellPrimitive(q).p;
  return p;
}
}  // namespace

int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 192;
  Flow::ImplicitControl control;
  control.tolerance = argc > 2 ? std::atof(argv[2]) : 1e-11;
  control.maxIterations = 400;
  control.pseudoCfl = argc > 3 ? std::atof(argv[3]) : 1e6;
  constexpr double p0 = 1e5, t0 = 300, length = 1.8;
  Definition d;
  d.experiment = Case::UniformDuct;
  d.length = length;
  d.inletRadius = 0.035;
  d.nz = n;
  d.nr = 2;
  d.totalPressure = p0;
  d.totalTemperature = t0;
  const double gamma = d.gas.gamma, a0 = std::sqrt(gamma * d.gas.specificR * t0), rho0 = p0 / (d.gas.specificR * t0);
  const Pulse pulse{p0, rho0, a0, gamma, 1e-3, length / 6, 8 * length / 192};
  const double end = 2.0 / 3 * length / a0, dz = length / n;
  std::printf("case P: simple-wave pulse (1e-3 p0), %d x %d, to t = %.6e s (2/3 L / a0); one thread; inner tolerance %.1e, "
              "pseudo CFL %.3g\n", n, d.nr, end, control.tolerance, control.pseudoCfl);
  // The explicit references.
  std::vector<double> reference;
  double explicitWallPerStep = 0;
  for (double cfl : {0.1, 0.05}) {
    Definition e = d;
    e.cfl = cfl;
    Flow f(e);
    pulse.fill(f);
    const auto clock = std::chrono::steady_clock::now();
    f.advanceTo(end);
    const double wall = seconds(clock);
    const auto m = f.measurements();
    if (reference.empty()) {
      reference = pressures(f);
      explicitWallPerStep = wall / static_cast<double>(m.steps);
      std::printf("explicit CFL %.2f: %llu steps, wall %.2f s (%.3e s per step); budgets mass %.2e, energy %.2e\n", cfl,
                  static_cast<unsigned long long>(m.steps), wall, explicitWallPerStep, m.massBalanceError, m.energyBalanceError);
    } else {
      std::printf("explicit CFL %.2f: %llu steps; its difference from CFL 0.1, the reference's own error bound: %.3e\n", cfl,
                  static_cast<unsigned long long>(m.steps), pressureError(f, reference, pulse.amplitude * p0));
    }
  }
  // The implicit runs at acoustic CFL 4, 2, 1 (and 8, 0.5 and 0.25, printed): the step count of CFL 4 rounded up, then doubled.
  const long steps4 = static_cast<long>(std::ceil(end * a0 / (4 * dz)));
  std::vector<double> errors;
  for (long steps : {steps4 / 2, steps4, 2 * steps4, 4 * steps4, 8 * steps4, 16 * steps4}) {
    if (steps < 2) continue;
    Flow f(d);
    pulse.fill(f);
    const double dt = end / static_cast<double>(steps);
    const auto clock = std::chrono::steady_clock::now();
    for (long k = 0; k < steps; ++k) f.implicitStep(k + 1 == steps ? end - f.time() : dt, control);
    const double wall = seconds(clock);
    const auto& s = f.implicitStats();
    const auto m = f.measurements();
    const double error = pressureError(f, reference, pulse.amplitude * p0);
    errors.push_back(error);
    std::printf("implicit acoustic CFL %.3f (%ld steps): error %.4e; inner iterations %.2f per step (most %d), %llu "
                "unconverged, %llu rejected; wall %.2f s, %.3e s per inner iteration = %.2f explicit steps; "
                "budgets mass %.2e, energy %.2e\n",
                dt * a0 / dz, steps, error, static_cast<double>(s.iterations) / static_cast<double>(s.steps), s.mostIterations,
                static_cast<unsigned long long>(s.unconverged), static_cast<unsigned long long>(s.rejected), wall,
                wall / static_cast<double>(s.iterations), wall / static_cast<double>(s.iterations) / explicitWallPerStep,
                m.massBalanceError, m.energyBalanceError);
  }
  int failures = 0;
  // errors[0] is CFL 8 when present (steps4 / 2 >= 2), then 4, 2, 1 (judged), 0.5 and 0.25 (printed).
  const std::size_t c4 = errors.size() - 5;
  for (std::size_t k = c4; k + 1 < errors.size(); ++k) {
    if (k >= c4 + 2) {
      std::printf("printed, not judged: order between acoustic CFL %g and %g: %.3f\n", 4.0 / (1 << (k - c4)),
                  2.0 / (1 << (k - c4)), std::log2(errors[k] / errors[k + 1]));
      continue;
    }
    const double order = std::log2(errors[k] / errors[k + 1]);
    const bool ok = order >= 1.8 && order <= 2.2;
    failures += !ok;
    std::printf("criterion 1: observed order between acoustic CFL %g and %g: %.3f (bar [1.8, 2.2])  %s\n",
                4.0 / (1 << (k - c4)), 2.0 / (1 << (k - c4)), order, ok ? "PASS" : "FAIL");
  }
  return failures ? 1 : 0;
}
