// Item-4 reacting-flow verification: 1-D detonations in 2H2 + O2 + 7Ar (h2o2.yaml) at 298 K and
// 6.67 kPa, Strang-split reacting Euler core. References come from Cantera equilibrium on the
// Hugoniot, independent of the flow solver.
//
// Piston mode (u_p > 0), the check: reactants at +u_p on the left half and -u_p on the right half
// collide at the midplane, which acts as a wall (a piston at u_p in the burned-gas frame). The
// reflected shock ignites the gas; the burned gas comes to rest in uniform equilibrium. For
// u_p > u_CJ the exact front speed in the lab is D - u_p with u_p = D (1 - v2/v1) on the
// equilibrium Hugoniot, and the plateau is that state (T2, p2). Mass conservation fixes the front
// speed once the reaction zone has stopped growing, so the slow recombination tail of this mixture
// cannot bias it. Tolerance stated before the run: front speed within 1% at the finest grid.
//
// Free mode (u_p = 0), reported only: a driver-initiated unsupported wave compared with D_CJ. Its
// heat release behind the sonic point (cm-scale recombination tail) cannot support the front,
// so at 0.45 m it runs below D_CJ by an amount that mixes grid and physics; not a pass/fail check.
//
// Usage: crucible_detonation_study <cells> <cfl> <threads> <u_p m/s, 0 = free> <length m, 0.5>
// Environment: TRACE=1 prints the front every 10 us; PROFILE=<file> writes the final field.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "cantera/core.h"

using namespace crucible;

namespace {
const char* kMoles = "H2:2, O2:1, AR:7";
constexpr double kT1 = 298.0, kP1 = 6670.0;

// Equilibrium Hugoniot of the reactants: at specific volume v2, the equilibrium state with
// e2 - e1 = (p1 + p2)(v1 - v2)/2 and the Rayleigh-line wave speed D = v1 sqrt((p2 - p1)/(v1 - v2)).
struct Hugoniot {
  std::shared_ptr<Cantera::Solution> sol = Cantera::newSolution("h2o2.yaml", "", "none");
  Cantera::ThermoPhase& gas = *sol->thermo();
  double v1, e1;
  std::vector<double> y0;
  Hugoniot() {
    gas.setState_TPX(kT1, kP1, kMoles);
    v1 = 1 / gas.density();
    e1 = gas.intEnergy_mass();
    y0.resize(gas.nSpecies());
    gas.getMassFractions(y0.data());
  }
  double speed(double v2) {  // leaves gas at the equilibrium Hugoniot state
    auto residual = [&](double t) {
      gas.setMassFractions(y0.data());
      gas.setState_TD(t, 1 / v2);
      gas.equilibrate("TV");
      return gas.intEnergy_mass() - e1 - 0.5 * (gas.pressure() + kP1) * (v1 - v2);
    };
    double lo = 500, hi = 8000;
    for (int n = 0; n < 80; ++n) {
      double m = 0.5 * (lo + hi);
      (residual(m) > 0 ? hi : lo) = m;
    }
    residual(0.5 * (lo + hi));
    return v1 * std::sqrt((gas.pressure() - kP1) / (v1 - v2));
  }
  double cjVolume() {  // golden-section minimum of D(v2)
    double a = 0.3 * v1, b = 0.9 * v1;
    const double r = 0.5 * (std::sqrt(5.0) - 1);
    double c = b - r * (b - a), d = a + r * (b - a), fc = speed(c), fd = speed(d);
    while (b - a > 1e-9 * v1) {
      if (fc < fd) { b = d; d = c; fd = fc; c = b - r * (b - a); fc = speed(c); }
      else { a = c; c = d; fc = fd; d = a + r * (b - a); fd = speed(d); }
    }
    return 0.5 * (a + b);
  }
  double pistonVolume(double up, double vcj) {  // strong branch: u2 = D (1 - v2/v1) rises as v2 falls
    double lo = 0.2 * v1, hi = vcj;
    for (int n = 0; n < 80; ++n) {
      double m = 0.5 * (lo + hi);
      (speed(m) * (1 - m / v1) > up ? lo : hi) = m;
    }
    return 0.5 * (lo + hi);
  }
};
}  // namespace

int main(int argc, char** argv) {
  const int cells = argc > 1 ? std::atoi(argv[1]) : 1000;
  const double cfl = argc > 2 ? std::atof(argv[2]) : 0.4;
  const int threads = argc > 3 ? std::atoi(argv[3]) : 4;
  const double up = argc > 4 ? std::atof(argv[4]) : 0;
  const bool piston = up > 0;
  const double length = argc > 5 ? std::atof(argv[5]) : 0.5;

  Hugoniot hugoniot;
  const double vcj = hugoniot.cjVolume(), dCJ = hugoniot.speed(vcj);
  std::printf("CJ: D %.3f m/s, T %.1f K, p/p1 %.4f, u_CJ %.2f m/s\n", dCJ, hugoniot.gas.temperature(),
              hugoniot.gas.pressure() / kP1, dCJ * (1 - vcj / hugoniot.v1));
  double target = dCJ, t2 = 0, p2 = 0;
  if (piston) {
    if (!(up > dCJ * (1 - vcj / hugoniot.v1))) { std::printf("u_p must exceed u_CJ\n"); return 1; }
    double v2 = hugoniot.pistonVolume(up, vcj), dp = hugoniot.speed(v2);
    target = dp - up;
    t2 = hugoniot.gas.temperature();
    p2 = hugoniot.gas.pressure();
    std::printf("Piston u_p %.1f m/s: D %.3f m/s (overdrive (D/D_CJ)^2 %.4f), lab front speed %.3f m/s, "
                "T2 %.1f K, p2/p1 %.4f\n", up, dp, dp * dp / (dCJ * dCJ), target, t2, p2 / kP1);
  }

  thermo::ReactionSource source("h2o2.yaml");
  Definition d;
  d.nz = cells; d.nr = 1; d.length = length; d.inletRadius = d.throatRadius = d.exitRadius = 0.01;
  d.experiment = Case::ShockTube; d.cfl = cfl;
  d.species = source.medium().species();
  d.composition = source.massFractions(kMoles);
  Flow flow(d);
  const double dx = d.length / cells;
  // Free mode initiation: the first 2 cm at 3000 K and 40 p1 (same reactants). A 1 cm, 2500 K,
  // 30 p1 driver failed to initiate on the 1 mm grid (decoupled shock at ~1100 m/s that re-ignited
  // after 0.15 m). Piston mode needs no driver: the collision shock ignites the gas.
  std::vector<Primitive> init;
  std::vector<double> y;
  const double rho1 = source.density(kT1, kP1, d.composition.data());
  for (int i = 0; i < cells; ++i) {
    double x = (i + 0.5) * dx;
    if (piston) {
      init.push_back({rho1, x < 0.5 * d.length ? up : -up, 0, kP1});
    } else {
      bool driver = x < 0.02;
      double t = driver ? 3000 : kT1, p = driver ? 40 * kP1 : kP1;
      init.push_back({source.density(t, p, d.composition.data()), 0, 0, p});
    }
    y.insert(y.end(), d.composition.begin(), d.composition.end());
  }
  flow.setInitialState(init, y);
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", threads);

  // Trace window: rightmost cell with p > 2 p1 (the right-hand inflow stays at p1). Speed: least-squares
  // slope of x(t) over the fit window, also split into its two halves to show steadiness.
  auto frontCell = [&] {
    for (int i = cells - 1; i >= 0; --i)
      if (flow.cellPrimitive(i).p > 2 * kP1) return i;
    return -1;
  };
  // Front position: the rightmost crossing of p = 10 p1 (mid-shock), linear between cell centres.
  // Using a whole cell index instead made x(t) a staircase whose fit aliased by a few 0.1%.
  auto front = [&] {
    const double level = 10 * kP1;
    for (int i = cells - 2; i >= 0; --i) {
      double pa = flow.cellPrimitive(i).p, pb = flow.cellPrimitive(i + 1).p;
      if (pa >= level && pb < level) return (i + 0.5 + (pa - level) / (pa - pb)) * dx;
    }
    return 0.0;
  };
  // Fit windows as fractions of the length: piston [0.70, 0.94] (front 0.4-0.88 of the half from
  // the midplane), free [0.60, 0.90].
  const double x0 = (piston ? 0.70 : 0.60) * length, x1 = (piston ? 0.94 : 0.90) * length,
               xm = 0.5 * (x0 + x1);
  const bool trace = std::getenv("TRACE") != nullptr;
  double nextTrace = 0;
  auto report = [&] {
    int f = frontCell(), back = std::max(0, f - int(0.02 / dx));
    double pmax = 0, tmax = 0;
    int ipmax = f, itmax = f;
    for (int i = back; i <= f; ++i) {
      if (flow.cellPrimitive(i).p > pmax) { pmax = flow.cellPrimitive(i).p; ipmax = i; }
      if (flow.temperature(i) > tmax) { tmax = flow.temperature(i); itmax = i; }
    }
    std::printf("  t %.1f us  x %.4f m  pmax/p1 %.2f at -%d cells  Tmax %.0f K at -%d cells\n",
                1e6 * flow.time(), front(), pmax / kP1, f - ipmax, tmax, f - itmax);
  };
  struct Fit {
    double st = 0, sx = 0, stt = 0, stx = 0;
    long n = 0;
    void add(double t, double x) { st += t; sx += x; stt += t * t; stx += t * x; ++n; }
    double slope() const { return (n * stx - st * sx) / (n * stt - st * st); }
  } all, first, second;
  auto clock = std::chrono::steady_clock::now();
  double x = 0;
  while ((x = front()) < x1) {
    if (x >= x0) {
      all.add(flow.time(), x);
      (x < xm ? first : second).add(flow.time(), x);
    }
    reacting.step();
    if (trace && flow.time() >= nextTrace) { report(); nextTrace += 10e-6; }
  }
  const double speed = all.slope();
  double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock).count();
  const auto& s = reacting.stats();
  std::printf("cells %d (dx %.4g mm) CFL %.2f: front %.3f m/s vs %.3f, error %+.4f%% (halves %+.4f%%, %+.4f%%)\n",
              cells, 1e3 * dx, cfl, speed, target, 100 * (speed / target - 1),
              100 * (first.slope() / target - 1), 100 * (second.slope() / target - 1));
  // Plateaus (burned gas at rest; exact state T2, p2, u = 0). Each parcel keeps the entropy of the
  // shock that burned it, so T records history. "early": 6-12 cm from the midplane, burned while
  // the wave was still settling, clear of the start-up entropy layer at the midplane (first gas
  // shocked inert then burned at constant volume, ~70 K cooler, within ~2 cm). "late": 0.40-0.60
  // of the half length, burned inside the fit window and behind the ~8 cm reaction zone at the end.
  // Windows chosen after viewing the 2 mm field.
  auto plateau = [&](const char* name, double from, double to) {
    double pm = 0, tm = 0, um = 0;
    int n = 0;
    for (int i = 0; i < cells; ++i) {
      double xc = (i + 0.5) * dx - 0.5 * length;
      if (xc < from || xc > to) continue;
      pm += flow.cellPrimitive(i).p; tm += flow.temperature(i); um += flow.cellPrimitive(i).uz; ++n;
    }
    std::printf("plateau %s [%.2f, %.2f] m: T %.1f K (exact %.1f, %+.3f%%), p %+.3f%%, u %.3f m/s\n", name, from,
                to, tm / n, t2, 100 * (tm / n / t2 - 1), 100 * (pm / n / p2 - 1), um / n);
  };
  if (piston) {
    plateau("early", 0.06, 0.12);
    plateau("late", 0.20 * length, 0.30 * length);
  }
  std::printf("steps %ld, replans %ld, asymmetric %ld (flow rejections %llu), max |T_chem - T_core| %.3e K, "
              "wall %.1f s (%d threads)\n", s.steps, s.replans, s.asymmetricSteps,
              static_cast<unsigned long long>(flow.measurements().rejectedSteps), s.maxTemperatureMismatch,
              wall, threads);
  if (const char* path = std::getenv("PROFILE")) {  // x, p, T, u, rho, Y(H2, O2, OH, H2O)
    std::FILE* out = std::fopen(path, "w");
    for (int i = 0; i < cells; ++i) {
      auto w = flow.cellPrimitive(i);
      auto yk = flow.massFractions(i);
      std::fprintf(out, "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", (i + 0.5) * dx, w.p,
                   flow.temperature(i), w.uz, w.rho, yk[0], yk[3], yk[4], yk[5]);
    }
    std::fclose(out);
  }
  return 0;
}
