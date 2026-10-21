// AMR, step B1 (docs/evidence/AMR.md): two levels or three in a fixed pattern, subcycled, refluxed and
// averaged down, every level's kernels over the whole level. Criteria 1 and 3, as the section *Step B1*
// reads them (recorded before this file was written):
//   1(a) quiescent gas in the closed RL10 development chamber on the base plus two levels, 1,000 base
//        steps: |u| < 1e-10 a0 and |p / p0 - 1| < 1e-12 on each level's active cells;
//   1(b) uniform Mach 0.3 flow in a stretched-ring cylinder, refined in both directions, 400 base steps:
//        the same bounds on |u - u0| where no wave from an end can have arrived;
//   3(a) a simple-wave pulse through a refined middle third on bases of 192, 384 and 768 columns: the
//        difference from each run's uniform finest mesh converges at order 1.8 or better on both pairs,
//        and the reflection is below 1e-3 of the pulse on 768 columns and falls with resolution;
//   3(b) Sod through a refined middle third on 64, 128 and 256 columns: the L1 density error between
//        the uniform coarse and fine runs' and converging at order 0.8 or better on both pairs, the
//        shock within one base cell, and the pressure behind it within 1e-3 of its jump of the uniform
//        coarse run's.
// Step C's check 3 (AMR.md, Step C: setup), case c3: 1(b)'s case with its two levels moved by a caller's tag one
// base column per regrid (every 2 base steps), 400 base steps: check accepts every pattern (each regrid changes it),
// and 1(b)'s bounds hold on every level's active cells.
// Arguments: any of 1a 1b 3a 3b c3 (default all) and a thread count (default 4). Step B2's checks (AMR.md,
// Step B2: setup): b1 runs B1 (with the temperature-guess reset) instead of B2, ghost=N sets the ghost
// width, and dump=PREFIX writes each run's active cells, steps and budget to PREFIX_<case>.bin.
#include "core/amr.hpp"
#include "core/flow.hpp"
#include "tests/amr_dump.hpp"
#include "tests/rl10_wall.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <set>
#include <string>
#include <vector>

using namespace crucible;

namespace {
int failures = 0, threads = 4, ghost = 0;
bool b1 = false;
std::string dumpPrefix;
// The B1 or B2 choice and the ghost width (0: the default) for each AMR run, and its dump.
void configure(Amr& amr) {
  amr.setRestricted(!b1);
  if (ghost > 0) amr.setGhostWidth(ghost);
}
void dump(const Amr& amr, const std::string& name) {
  if (!dumpPrefix.empty()) writeAmrDump(amr, dumpPrefix + "_" + name + ".bin");
}
void check(bool ok, const std::string& what, double value, double limit) {
  std::printf("  %-74s %.3e (limit %.3g)  %s\n", what.c_str(), value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}
double seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
// Level l's active mask from a rule on the level's (i, j).
std::vector<unsigned char> pattern(const Flow& f, const std::function<bool(int, int)>& rule) {
  const Mesh& m = f.mesh();
  std::vector<unsigned char> a(m.cells.size(), 0);
  for (int i = 0; i < m.nz; ++i)
    for (int j = 0; j < m.nr; ++j) a[m.index(i, j)] = rule(i, j) ? 1 : 0;
  return a;
}
// The largest |u - (u0, 0)| / a0 and |p / p0 - 1| over the cells a mask admits (an empty mask admits
// all) and a filter on the centroid's z.
struct Departure { double speed = 0, pressure = 0; };
Departure departure(const Flow& f, const std::vector<unsigned char>& mask, double u0, double a0, double p0,
                    const std::function<bool(double)>& where = [](double) { return true; }) {
  Departure d;
  for (std::size_t q = 0; q < f.mesh().cells.size(); ++q) {
    if ((!mask.empty() && !mask[q]) || !where(f.mesh().cells[q].z)) continue;
    const Primitive w = f.cellPrimitive(q);
    d.speed = std::max(d.speed, std::hypot(w.uz - u0, w.ur) / a0);
    d.pressure = std::max(d.pressure, std::abs(w.p / p0 - 1));
  }
  return d;
}

// 1(a): the closed RL10 development chamber at rest.
void quiescentChamber() {
  std::printf("1(a) quiescent RL10 development chamber, base 150 x 16 plus two levels, 1,000 base steps\n");
  constexpr double p0 = 1e6, t0 = 300;
  Definition d = rl10Wall(150, 150, 16, 1.7347);
  d.experiment = Case::Chamber;
  d.backPressure = d.totalPressure = p0;
  d.ambientTemperature = d.totalTemperature = t0;
  d.outletRelaxation = 0;
  Supply core;
  core.innerRadius = 0;
  core.outerRadius = 0.05;
  core.massFlow = 1;
  core.totalTemperature = t0;
  core.composition = {1.0};
  core.opens = 1.0;
  core.ramp = 1e-3;
  Supply annulus = core;
  annulus.innerRadius = 0.05;
  annulus.outerRadius = 0.1;
  d.supplies = {core, annulus};
  const double a0 = std::sqrt(d.gas.gamma * d.gas.specificR * t0);

  Amr amr(d, 2);
  configure(amr);
  amr.setThreads(threads);
  const Mesh& base = amr.level(0).mesh();
  int throat = 0;
  for (int i = 1; i < base.nz; ++i)
    if (base.radius[i] + base.radius[i + 1] < base.radius[throat] + base.radius[throat + 1]) throat = i;
  std::printf("  throat column %d of %d (R %.4f to %.4f m)\n", throat, base.nz, base.radius[throat], base.radius[throat + 1]);
  amr.setActive(1, pattern(amr.level(1), [&](int i, int j) { return (j >> 1) >= 12 || ((i >> 1) >= throat - 10 && (i >> 1) < throat + 10); }));
  amr.setActive(2, pattern(amr.level(2), [&](int i, int j) {
    const int i1 = i >> 1, j1 = j >> 1;
    return j1 >= 26 || (i1 >= 2 * throat - 18 && i1 < 2 * throat + 18);
  }));
  const auto start = std::chrono::steady_clock::now();
  for (int n = 0; n < 1000; ++n) amr.step();
  dump(amr, "1a");
  std::printf("  AMR: t = %.4e s after %llu base steps (%llu redone), %.1f s of wall\n", amr.time(),
              static_cast<unsigned long long>(amr.steps()), static_cast<unsigned long long>(amr.rejectedSteps()), seconds(start));
  for (int l = 0; l < amr.levels(); ++l) {
    const Flow& f = amr.level(l);
    std::size_t active = 0;
    for (auto a : amr.active(l)) active += a;
    const auto e = departure(f, amr.active(l), 0, a0, p0);
    const std::string name = "level " + std::to_string(l) + " (" + std::to_string(f.mesh().nz) + " x " +
                             std::to_string(f.mesh().nr) + ", " + std::to_string(l ? active : f.mesh().cells.size()) + " active)";
    check(e.speed < 1e-10, name + ": largest |u| / a0", e.speed, 1e-10);
    check(e.pressure < 1e-12, name + ": largest |p / p0 - 1|", e.pressure, 1e-12);
  }
  // Printed beside: each of the base's and the finest level's meshes alone, 1,000 of its own steps.
  for (int l : {0, 2}) {
    Definition s = d;
    s.nz = 150 << l;
    s.nr = 16 << l;
    Flow f(s, Mesh(amr.level(l).mesh()));
    f.setThreads(threads);
    for (int n = 0; n < 1000; ++n) f.step();
    const auto e = departure(f, {}, 0, a0, p0);
    std::printf("      printed, not judged: %d x %d alone after 1,000 steps (t = %.4e s): |u| / a0 %.3e, |p / p0 - 1| %.3e\n",
                s.nz, s.nr, f.measurements().time, e.speed, e.pressure);
  }
}

// 1(b): uniform axial flow in a stretched-ring cylinder.
void uniformFlow() {
  std::printf("1(b) uniform Mach 0.3 flow, stretched-ring cylinder 400 x 8 plus two levels, 400 base steps\n");
  constexpr double p0 = 1e5, t0 = 300;
  Definition d;
  d.experiment = Case::UniformDuct;
  d.length = 0.6;
  d.inletRadius = 0.035;
  d.nz = 400;
  d.nr = 8;
  d.radialStretching = 1.7347;
  d.totalPressure = p0;
  d.totalTemperature = t0;
  const double a0 = std::sqrt(d.gas.gamma * d.gas.specificR * t0), u0 = 0.3 * a0, rho0 = p0 / (d.gas.specificR * t0);
  Amr amr(d, 2);
  configure(amr);
  amr.setThreads(threads);
  for (int l = 0; l < amr.levels(); ++l) amr.level(l).setUniform({rho0, u0, 0, p0});
  amr.setActive(1, pattern(amr.level(1), [](int i, int j) { return (i >> 1) >= 150 && (i >> 1) < 250 && (j >> 1) >= 2 && (j >> 1) < 6; }));
  amr.setActive(2, pattern(amr.level(2), [](int i, int j) { return (i >> 2) >= 170 && (i >> 2) < 230 && (j >> 2) >= 3 && (j >> 2) < 5; }));
  amr.averageDown();
  const auto start = std::chrono::steady_clock::now();
  for (int n = 0; n < 400; ++n) amr.step();
  dump(amr, "1b");
  const double t = amr.time(), lo = (u0 + a0) * t, hi = d.length - (a0 - u0) * t;
  std::printf("  t = %.4e s after %llu base steps (%llu redone), %.1f s of wall; judged for z in (%.4f, %.4f) m\n", t,
              static_cast<unsigned long long>(amr.steps()), static_cast<unsigned long long>(amr.rejectedSteps()), seconds(start), lo, hi);
  for (int l = 0; l < amr.levels(); ++l) {
    const Flow& f = amr.level(l);
    const auto e = departure(f, amr.active(l), u0, a0, p0, [&](double z) { return z > lo && z < hi; });
    const auto all = departure(f, amr.active(l), u0, a0, p0);
    const std::string name = "level " + std::to_string(l) + " (" + std::to_string(f.mesh().nz) + " x " + std::to_string(f.mesh().nr) + ")";
    check(e.speed < 1e-10, name + ": largest |u - u0| / a0", e.speed, 1e-10);
    check(e.pressure < 1e-12, name + ": largest |p / p0 - 1|", e.pressure, 1e-12);
    std::printf("      printed, not judged: over all its active cells %.3e and %.3e\n", all.speed, all.pressure);
  }
}

// Step C's check 3: 1(b)'s flow, the levels moved by a tag one base column per regrid.
void movingPattern() {
  std::printf("Step C check 3: 1(b)'s flow, two levels moved one base column per regrid, 400 base steps\n");
  constexpr double p0 = 1e5, t0 = 300;
  Definition d;
  d.experiment = Case::UniformDuct;
  d.length = 0.6;
  d.inletRadius = 0.035;
  d.nz = 400;
  d.nr = 8;
  d.radialStretching = 1.7347;
  d.totalPressure = p0;
  d.totalTemperature = t0;
  const double a0 = std::sqrt(d.gas.gamma * d.gas.specificR * t0), u0 = 0.3 * a0, rho0 = p0 / (d.gas.specificR * t0);
  Amr amr(d, 2);
  configure(amr);
  amr.setThreads(threads);
  for (int l = 0; l < amr.levels(); ++l) amr.level(l).setUniform({rho0, u0, 0, p0});
  // 1(b)'s pattern moved n base columns at the regrid after base step 2n: level 0's tag sets level 1, level 1's level 2.
  Amr::Tags tags;
  tags.throatHalfWidth = 0;
  tags.wallYPlus = 0;
  tags.mixture = tags.temperature = tags.pressure = tags.vorticity = 0;
  tags.custom = [&amr](int l, int i, int j) {
    const int n = static_cast<int>(amr.steps() / 2), column = (i >> l) - n, row = j >> l;
    return l == 0 ? column >= 50 && column < 150 && row >= 2 && row < 6 : column >= 70 && column < 130 && row >= 3 && row < 5;
  };
  amr.setTags(tags);
  const auto start = std::chrono::steady_clock::now();
  for (int n = 0; n < 400; ++n) amr.step();
  dump(amr, "c3");
  const double t = amr.time(), lo = (u0 + a0) * t, hi = d.length - (a0 - u0) * t;
  std::printf("  t = %.4e s after %llu base steps (%llu redone), %.1f s of wall; judged for z in (%.4f, %.4f) m\n", t,
              static_cast<unsigned long long>(amr.steps()), static_cast<unsigned long long>(amr.rejectedSteps()), seconds(start), lo, hi);
  // Two initial passes, then one regrid after every second base step, each moving the pattern.
  check(amr.regrids() == 202, "regrids that changed the pattern (202 expected)", static_cast<double>(amr.regrids()), 202);
  for (int l = 1; l < amr.levels(); ++l) {
    int first = amr.level(l).mesh().nz, last = -1;
    const Mesh& m = amr.level(l).mesh();
    for (int i = 0; i < m.nz; ++i)
      for (int j = 0; j < m.nr; ++j)
        if (amr.active(l)[m.index(i, j)]) first = std::min(first, i >> l), last = std::max(last, i >> l);
    std::printf("      printed, not judged: level %d's active cells span base columns %d to %d\n", l, first, last);
  }
  for (int l = 0; l < amr.levels(); ++l) {
    const Flow& f = amr.level(l);
    const auto e = departure(f, amr.active(l), u0, a0, p0, [&](double z) { return z > lo && z < hi; });
    const auto all = departure(f, amr.active(l), u0, a0, p0);
    const std::string name = "level " + std::to_string(l) + " (" + std::to_string(f.mesh().nz) + " x " + std::to_string(f.mesh().nr) + ")";
    check(e.speed < 1e-10, name + ": largest |u - u0| / a0", e.speed, 1e-10);
    check(e.pressure < 1e-12, name + ": largest |p / p0 - 1|", e.pressure, 1e-12);
    std::printf("      printed, not judged: over all its active cells %.3e and %.3e\n", all.speed, all.pressure);
  }
}

// 3(a): a simple-wave pulse through the refined middle third.
constexpr double gaussX[6] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
                              0.2386191860831969,  0.6612093864662645,  0.9324695142031521};
constexpr double gaussW[6] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910,
                              0.4679139345726910, 0.3607615730481386, 0.1713244923791704};
struct Pulse {
  double p0, rho0, a0, gamma, amplitude, centre, halfWidth;
  // Conserved (rho, rho u, E) of the simple wave at z: isentropic, u - 2a / (gamma - 1) = -2 a0 / (gamma - 1).
  std::array<double, 3> at(double z) const {
    const double s = (z - centre) / halfWidth, p = p0 * (1 + amplitude * std::exp(-std::log(2.0) * s * s));
    const double rho = rho0 * std::pow(p / p0, 1 / gamma), a = std::sqrt(gamma * p / rho), u = 2 * (a - a0) / (gamma - 1);
    return {rho, rho * u, p / (gamma - 1) + 0.5 * rho * u * u};
  }
  // Every cell of a uniform-radius flow: the cell average over its z extent (6-point Gauss).
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
// Volume-weighted mean density of a uniform run on (2n x 2m) cells over each cell of an n x m mesh.
std::vector<double> densityOnto(const Flow& fine, const Mesh& coarse) {
  const Mesh& m = fine.mesh();
  std::vector<double> rho(coarse.cells.size());
  for (int i = 0; i < coarse.nz; ++i)
    for (int j = 0; j < coarse.nr; ++j) {
      double s = 0, v = 0;
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
          const auto q = m.index(2 * i + a, 2 * j + b);
          s += m.cells[q].volume * fine.state()[q][0];
          v += m.cells[q].volume;
        }
      rho[coarse.index(i, j)] = s / v;
    }
  return rho;
}
void acousticPulse() {
  std::printf("3(a) simple-wave pulse (1e-3 p0, half-width at half maximum 8 cells of the 192-column base) through a refined middle third\n");
  constexpr double p0 = 1e5, t0 = 300, length = 1.8;
  std::vector<double> errors, reflections;
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
    configure(amr);
    // The uniform runs on the AMR levels' own meshes.
    Flow coarseRun(d, Mesh(amr.level(0).mesh())), fineRun(fine, Mesh(amr.level(1).mesh()));
    for (Flow* f : {&amr.level(0), &amr.level(1), &coarseRun, &fineRun}) {
      f->setThreads(threads);
      pulse.fill(*f);
    }
    amr.setThreads(threads);
    amr.setActive(1, pattern(amr.level(1), [&](int i, int) { return (i >> 1) >= n / 3 && (i >> 1) < 2 * n / 3; }));
    amr.averageDown();
    const auto start = std::chrono::steady_clock::now();
    const double end = length / a0;
    double reflection = 0, error = 0;
    for (int k = 1; k <= 300; ++k) {
      const double t = k * end / 300;
      amr.advanceTo(t);
      coarseRun.advanceTo(t);
      for (int i = 0; i < n / 3; ++i)
        for (int j = 0; j < d.nr; ++j) {
          const auto q = coarseRun.mesh().index(i, j);
          reflection = std::max(reflection, std::abs(amr.level(0).cellPrimitive(q).p - coarseRun.cellPrimitive(q).p) / (pulse.amplitude * p0));
        }
      if (k == 200) {
        fineRun.advanceTo(t);
        const auto reference = densityOnto(fineRun, coarseRun.mesh());
        const Mesh& m = coarseRun.mesh();
        double s = 0, v = 0;
        for (std::size_t q = 0; q < m.cells.size(); ++q) {
          s += m.cells[q].volume * std::abs(amr.level(0).state()[q][0] - reference[q]);
          v += m.cells[q].volume;
        }
        error = s / v / (rho0 * (std::pow(1 + pulse.amplitude, 1 / gamma) - 1));
      }
    }
    dump(amr, "3a_" + std::to_string(n));
    std::printf("  %d columns: %llu base steps, %.1f s of wall; difference from the uniform %d-column run %.4e of the "
                "pulse's density amplitude; reflection %.4e of its pressure amplitude\n",
                n, static_cast<unsigned long long>(amr.steps()), seconds(start), 2 * n, error, reflection);
    errors.push_back(error);
    reflections.push_back(reflection);
  }
  for (std::size_t k = 1; k < errors.size(); ++k) {
    const double order = std::log2(errors[k - 1] / errors[k]);
    check(order >= 1.8, "observed order, " + std::to_string(192 << (k - 1)) + " to " + std::to_string(192 << k) + " columns (at least)",
          order, 1.8);
  }
  check(reflections[2] < 1e-3, "reflection on 768 columns, relative to the pulse", reflections[2], 1e-3);
  for (std::size_t k = 1; k < reflections.size(); ++k)
    check(reflections[k] < reflections[k - 1], "reflection falls, " + std::to_string(192 << (k - 1)) + " to " + std::to_string(192 << k) +
          " columns (ratio below 1)", reflections[k] / reflections[k - 1], 1);
}

// 3(b): Sod through the refined middle third. The exact solution as core_tests' (gamma 1.4).
constexpr double sodPstar = 0.303130178050647, sodUstar = 0.92745262004895, sodG = 1.4;
double sodShockSpeed() {
  const double ar = std::sqrt(sodG * 0.1 / 0.125);
  return ar * std::sqrt((sodG + 1) / (2 * sodG) * sodPstar / 0.1 + (sodG - 1) / (2 * sodG));
}
double sodPostShockDensity() {
  const double g6 = (sodG - 1) / (sodG + 1);
  return 0.125 * (sodPstar / 0.1 + g6) / (g6 * sodPstar / 0.1 + 1);
}
double sodDensity(double x, double t) {
  const double xi = (x - 0.5) / t, al = std::sqrt(sodG), astar = al * std::pow(sodPstar, (sodG - 1) / (2 * sodG));
  if (xi < -al) return 1;
  if (xi < sodUstar - astar) return std::pow(2 / (sodG + 1) + (sodG - 1) / (sodG + 1) * (-xi / al), 2 / (sodG - 1));
  if (xi < sodUstar) return std::pow(sodPstar, 1 / sodG);
  if (xi < sodShockSpeed()) return sodPostShockDensity();
  return 0.125;
}
// Density L1 error of row 0 against the exact solution, per unit length.
double sodError(const Flow& f, double t, const std::function<bool(int)>& skip = [](int) { return false; }) {
  const Mesh& m = f.mesh();
  double e = 0;
  for (int i = 0; i < m.nz; ++i)
    if (!skip(i)) e += std::abs(f.state()[m.index(i, 0)][0] - sodDensity((i + 0.5) * m.dz, t)) * m.dz;
  return e;
}
void sod() {
  std::printf("3(b) Sod shock tube, level 1 over the middle third, to t = 0.25\n");
  constexpr double t = 0.25;
  const double jump = sodPstar - 0.1, middle = 0.5 * (sodPostShockDensity() + 0.125), shock = 0.5 + sodShockSpeed() * t;
  const double al = std::sqrt(sodG), tail = 0.5 + (sodUstar - al * std::pow(sodPstar, (sodG - 1) / (2 * sodG))) * t;
  std::vector<double> errors;
  for (int n : {64, 128, 256}) {
    Definition d;
    d.experiment = Case::ShockTube;
    d.length = 1;
    d.nz = n;
    d.nr = 2;
    Definition fine = d;
    fine.nz = 2 * n;
    fine.nr = 4;
    const int first = static_cast<int>(std::lround(n / 3.0)), last = static_cast<int>(std::lround(2 * n / 3.0));
    Amr amr(d, 1);
    configure(amr);
    amr.setThreads(threads);
    amr.setActive(1, pattern(amr.level(1), [&](int i, int) { return (i >> 1) >= first && (i >> 1) < last; }));
    amr.averageDown();
    Flow coarseRun(d, Mesh(amr.level(0).mesh())), fineRun(fine, Mesh(amr.level(1).mesh()));
    coarseRun.setThreads(threads);
    fineRun.setThreads(threads);
    amr.advanceTo(t);
    dump(amr, "3b_" + std::to_string(n));
    coarseRun.advanceTo(t);
    fineRun.advanceTo(t);
    const Flow& base = amr.level(0);
    const Mesh& m = base.mesh();
    const double h = m.dz;
    const double composite = sodError(base, t, [&](int i) { return i >= first && i < last; }) +
                             sodError(amr.level(1), t, [&](int i) { return (i >> 1) < first || (i >> 1) >= last; });
    const double coarseError = sodError(coarseRun, t), fineError = sodError(fineRun, t);
    std::printf("  %d columns (level 1 over base columns %d to %d, the interfaces at %.4f and %.4f): L1 %.4e; uniform %d: %.4e, uniform %d: %.4e\n",
                n, first, last - 1, first * h, last * h, composite, n, coarseError, 2 * n, fineError);
    check(composite <= coarseError && composite >= fineError,
          std::to_string(n) + " columns: L1 between the uniform runs' (its place, 0 fine to 1 coarse)",
          (composite - fineError) / (coarseError - fineError), 1);
    // The shock: the last place, from the outlet, where the density crosses the middle of its two states.
    double place = 0;
    for (int i = m.nz - 2; i >= last; --i) {
      const double r0 = base.state()[m.index(i, 0)][0], r1 = base.state()[m.index(i + 1, 0)][0];
      if (r0 >= middle && r1 < middle) { place = (i + 0.5 + (r0 - middle) / (r0 - r1)) * h; break; }
    }
    check(std::abs(place - shock) <= h, std::to_string(n) + " columns: shock place against exact " + std::to_string(shock) + ", in base cells",
          std::abs(place - shock) / h, 1);
    double behind = 0, other = 0;
    for (int i = 0; i < m.nz; ++i) {
      const double x = (i + 0.5) * h, dp = std::abs(base.cellPrimitive(m.index(i, 0)).p - coarseRun.cellPrimitive(m.index(i, 0)).p) / jump;
      if (i >= last && x <= shock - 5 * h) behind = std::max(behind, dp);
      if (i < last && i >= first && x >= tail) other = std::max(other, dp);
    }
    check(behind <= 1e-3, std::to_string(n) + " columns: |p - p_uniform coarse| / jump, interface to 5 cells behind the shock", behind, 1e-3);
    std::printf("      printed, not judged: the same from the rarefaction's tail (%.4f) to the interface %.3e\n", tail, other);
    errors.push_back(composite);
  }
  for (std::size_t k = 1; k < errors.size(); ++k) {
    const double order = std::log2(errors[k - 1] / errors[k]);
    check(order >= 0.8, "observed order, " + std::to_string(64 << (k - 1)) + " to " + std::to_string(64 << k) + " columns (at least)", order, 0.8);
  }
}
}  // namespace

int main(int argc, char** argv) {
  std::set<std::string> cases;
  for (int k = 1; k < argc; ++k) {
    const std::string a = argv[k];
    if (a == "1a" || a == "1b" || a == "3a" || a == "3b" || a == "c3") cases.insert(a);
    else if (a == "b1") b1 = true;
    else if (a.rfind("ghost=", 0) == 0) ghost = std::atoi(a.c_str() + 6);
    else if (a.rfind("dump=", 0) == 0) dumpPrefix = a.substr(5);
    else threads = std::max(1, std::atoi(a.c_str()));
  }
  if (cases.empty()) cases = {"1a", "1b", "3a", "3b", "c3"};
  std::printf("AMR step %s, criteria 1 and 3 (%d threads%s)\n", b1 ? "B1" : "B2", threads,
              ghost > 0 ? (", ghost width " + std::to_string(ghost)).c_str() : "");
  if (cases.count("1b")) uniformFlow();
  if (cases.count("c3")) movingPattern();
  if (cases.count("3b")) sod();
  if (cases.count("3a")) acousticPulse();
  if (cases.count("1a")) quiescentChamber();
  std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
