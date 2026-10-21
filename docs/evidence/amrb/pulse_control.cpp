// Scratch control for 3(a): the same simple-wave pulse on uniform N against uniform 2N, no AMR, at the
// time 3(a) judges (2L / (3 a0)), with 3(a)'s norm: the base scheme's own order on this pulse.
// Built from the repository root after the library: c++ -std=c++20 -O2 -I. docs/evidence/amrb/pulse_control.cpp build/libcrucible_core.a -o pulse_control
#include "core/amr.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <vector>
using namespace crucible;
constexpr double gaussX[6] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
                              0.2386191860831969,  0.6612093864662645,  0.9324695142031521};
constexpr double gaussW[6] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
                              0.4679139345726910, 0.3607615730481386, 0.1713244923791704};
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
int main() {
  constexpr double p0 = 1e5, t0 = 300, length = 1.8;
  double previous = 0;
  for (int n : {192, 384, 768}) {
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
    Definition fine = d;
    fine.nz = 2 * n;
    fine.nr = 4;
    Amr amr(d, 1);
    Flow coarse(d, Mesh(amr.level(0).mesh())), f(fine, Mesh(amr.level(1).mesh()));
    coarse.setThreads(4);
    f.setThreads(4);
    pulse.fill(coarse);
    pulse.fill(f);
    const double t = 2 * length / (3 * a0), amp = rho0 * (std::pow(1 + pulse.amplitude, 1 / gamma) - 1);
    coarse.advanceTo(t);
    f.advanceTo(t);
    const Mesh& m = coarse.mesh();
    const Mesh& fm = f.mesh();
    double s = 0, v = 0, peakC = 0, peakF = 0;
    for (int i = 0; i < m.nz; ++i)
      for (int j = 0; j < m.nr; ++j) {
        double fs = 0, fv = 0;
        for (int a = 0; a < 2; ++a)
          for (int b = 0; b < 2; ++b) {
            const auto q = fm.index(2 * i + a, 2 * j + b);
            fs += fm.cells[q].volume * f.state()[q][0];
            fv += fm.cells[q].volume;
          }
        const auto q = m.index(i, j);
        s += m.cells[q].volume * std::abs(coarse.state()[q][0] - fs / fv);
        v += m.cells[q].volume;
        peakC = std::max(peakC, coarse.state()[q][0] - rho0);
        peakF = std::max(peakF, fs / fv - rho0);
      }
    const double e = s / v / amp;
    std::printf("%d columns: uniform %d against uniform %d, %.4e of the density amplitude; peak %.5f and %.5f of it", n, n, 2 * n, e,
                peakC / amp, peakF / amp);
    if (previous > 0) std::printf("; order %.3f", std::log2(previous / e));
    std::printf("\n");
    std::fflush(stdout);
    previous = e;
  }
}
