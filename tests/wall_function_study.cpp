// Wall functions in the engine (docs/evidence/WALL_FUNCTIONS.md; criteria stated 5 October 2026 before
// any wall-function code).
//
// threads: criterion 0(d). A turbulent N2 chamber with wall functions on the side wall (contoured, so
//   its normals are not radial) and on plate rings, including the corner ring, at an isothermal 400 K
//   wall: 200 steps on 1 and on 4 threads. memcmp of the state, partial densities, rho k and rho omega
//   and the budgets; the mass and energy budgets within 1e-11. Reported: the same run without wall
//   functions, the wall heat flow and the rho k V the prescription added.
//
// apriori: the law on the reference profile, without a march (below).
//
// pipe: criteria 1 and 2. A fully developed N2 pipe (UniformDuct, 4 columns, axial body force, and for
//   criterion 2 uniform heating) started from a resolved reference profile of tools/sst_pipe_1d.py
//   interpolated (linear in r) to the engine's ring centroids, the density scaled so the bulk density
//   is the reference's. It marches until u_b changes by less than 1e-5 per ms (samples every 0.5 ms),
//   or to the end time.
//   Usage: crucible_wall_function_study pipe <profile.csv> <nr> <stretching> <force N/m^3> <T_wall K>
//          <heating W/m^3> <end ms> <threads> <wall: resolved|law|printed>
//   Prints u_b, c_f = 2 tau_w / (rho_b u_b^2) with tau_w from the wall force, T_axis (the axis ring's
//   mean), the mixing-cup temperature T_b (weighted by rho u, as the reference's), the Stanton number
//   St = q_w / (rho_b cp_b u_b (T_b - T_w)) with q_w from the wall heat flow and cp_b at T_b, the
//   first-cell y+ from the steady tau_w = f R / 2, and the budgets.
#include "adapters/reaction.hpp"
#include "core/flow.hpp"
#include "core/walls.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace crucible;

namespace {
int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-66s %.3e (limit %.0e)  %s\n", what, value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}
double sq(double x) { return x * x; }
struct Nitrogen {
  std::vector<Species> species;
  TransportFits fits;
  Medium medium;
};
Nitrogen nitrogen() {
  const std::string mechanism = "h2o2.yaml";
  thermo::ReactionSource source(mechanism);
  auto full = source.medium();
  const auto fits = thermo::transportFits(mechanism);
  const std::size_t n = full.size();
  std::size_t i = n;
  for (std::size_t k = 0; k < n; ++k) if (full.species()[k].name == "N2") i = k;
  if (i == n) throw std::runtime_error("no N2 in the mechanism");
  const std::size_t pair = i * n - i * (i - 1) / 2;
  Nitrogen out{{full.species()[i]}, {{fits.viscosity[i]}, {fits.conductivity[i]}, {fits.diffusion[pair]}}, {}};
  out.medium = Medium(out.species);
  out.medium.setTransport(out.fits);
  return out;
}

int threads(const Nitrogen& n2) {
  std::printf("criterion 0(d): a turbulent N2 chamber, wall functions on the side wall and the plate, 200 steps\n");
  auto run = [&](int count, bool law) {
    Definition d;
    d.experiment = Case::Chamber;
    d.nz = 16; d.nr = 6;
    d.contour = {{0, 5e-3}, {0.02, 5e-3}, {0.03, 3.5e-3}};
    d.species = n2.species; d.composition = {1};
    d.transport = n2.fits;
    d.backPressure = 101325; d.ambientTemperature = 300;
    d.wallTemperature = 400;
    d.turbulence.enabled = true; d.turbulence.ambientK = 1; d.turbulence.ambientOmega = 1e3;
    d.turbulence.wallFunctions = law;
    Supply s;
    s.innerRadius = 0; s.outerRadius = 3e-3; s.massFlow = 2e-3; s.totalTemperature = 300;
    s.composition = {1}; s.opens = 0; s.ramp = 2e-5;
    s.turbulenceIntensity = 0.05; s.viscosityRatio = 10;
    d.supplies = {s};
    Flow flow(d);
    flow.setThreads(count);
    for (int n = 0; n < 200; ++n) flow.step();
    return flow;
  };
  auto one = run(1, true), four = run(4, true);
  const auto a = one.measurements(), b = four.measurements();
  const bool same = one.state() == four.state() && one.partialDensities() == four.partialDensities() &&
                    one.turbulence() == four.turbulence() &&
                    std::memcmp(one.state().data(), four.state().data(), one.state().size() * sizeof(Conserved)) == 0 &&
                    std::memcmp(one.turbulence().data(), four.turbulence().data(), one.turbulence().size() * sizeof(double)) == 0;
  const std::array<double, 6> la{a.massBalanceError, a.energyBalanceError, a.momentumBalanceError, a.wallHeatFlow, a.prescribedTurbulentEnergy, a.clippedTurbulentEnergy};
  const std::array<double, 6> lb{b.massBalanceError, b.energyBalanceError, b.momentumBalanceError, b.wallHeatFlow, b.prescribedTurbulentEnergy, b.clippedTurbulentEnergy};
  const bool sameBudgets = std::memcmp(la.data(), lb.data(), sizeof la) == 0;
  std::printf("  t %.6e s after %llu steps; wall heat flow %.6e W; prescribed rho k V %.6e J; clipped %.3e J\n", a.time,
              static_cast<unsigned long long>(a.steps), a.wallHeatFlow, a.prescribedTurbulentEnergy, a.clippedTurbulentEnergy);
  std::printf("  budgets: mass %.3e, energy %.3e, axial momentum %.3e\n", a.massBalanceError, a.energyBalanceError, a.momentumBalanceError);
  check(same, "(d) state, partial densities, rho k and rho omega: 1 and 4 threads differ (1 = yes)", same ? 0 : 1, 0);
  check(sameBudgets, "(d) budgets, wall heat and turbulent ledgers identical (1 = differ)", sameBudgets ? 0 : 1, 0);
  check(std::abs(a.massBalanceError) < 1e-11, "(d) mass budget", std::abs(a.massBalanceError), 1e-11);
  check(std::abs(a.energyBalanceError) < 1e-11, "(d) energy budget", std::abs(a.energyBalanceError), 1e-11);
  auto resolved = run(4, false);
  const auto r = resolved.measurements();
  std::printf("  reported, the same run with the resolved-wall rule: wall heat flow %.6e W, wall axial force %.6e N (law %.6e N), budgets mass %.3e energy %.3e\n",
              r.wallHeatFlow, r.wallAxialForce, a.wallAxialForce, r.massBalanceError, r.energyBalanceError);
  // The first cells' prescribed values satisfy k = omega mu_t / rho with the engine's eddy viscosity.
  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}

struct Profile { std::vector<double> r, u, k, omega, t, rho; };
Profile readProfile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  Profile p;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    std::stringstream s(line);
    std::string cell;
    std::vector<double> v;
    while (std::getline(s, cell, ',')) v.push_back(std::stod(cell));
    if (v.size() < 6) throw std::runtime_error("short row in " + path);
    p.r.push_back(v[0]); p.u.push_back(v[1]); p.k.push_back(v[2]); p.omega.push_back(v[3]); p.t.push_back(v[4]); p.rho.push_back(v[5]);
  }
  // Ascending r.
  if (p.r.front() > p.r.back())
    for (auto* v : {&p.r, &p.u, &p.k, &p.omega, &p.t, &p.rho}) std::reverse(v->begin(), v->end());
  return p;
}
double interpolate(const std::vector<double>& x, const std::vector<double>& y, double at) {
  if (at <= x.front()) return y.front();
  if (at >= x.back()) return y.back();
  const auto i = static_cast<std::size_t>(std::upper_bound(x.begin(), x.end(), at) - x.begin());
  const double w = (at - x[i - 1]) / (x[i] - x[i - 1]);
  return y[i - 1] + w * (y[i] - y[i - 1]);
}

int pipe(const Nitrogen& n2, int argc, char** argv) {
  if (argc < 11) throw std::invalid_argument("pipe <profile.csv> <nr> <stretching> <force> <T_wall> <heating> <end ms> <threads> <resolved|law|printed>");
  const std::string path = argv[2], wall = argv[10];
  const int nr = std::atoi(argv[3]), threadCount = std::atoi(argv[9]);
  const double stretching = std::atof(argv[4]), force = std::atof(argv[5]), tWall = std::atof(argv[6]), heating = std::atof(argv[7]);
  const double endMs = std::atof(argv[8]);
  if (wall != "resolved" && wall != "law" && wall != "printed") throw std::invalid_argument("wall: resolved, law or printed");
  const auto profile = readProfile(path);
  const double radius = 5e-3, length = 2 * radius, y[1] = {1}, gas = n2.medium.gasConstant(y);
  Definition d;
  d.experiment = Case::UniformDuct;
  d.nz = 4; d.nr = nr; d.length = length; d.inletRadius = d.exitRadius = d.throatRadius = radius;
  d.radialStretching = stretching;
  d.species = n2.species; d.composition = {1};
  d.transport = n2.fits;
  d.turbulence.enabled = true; d.turbulence.ambientK = 0; d.turbulence.ambientOmega = 1;
  d.turbulence.wallFunctions = wall != "resolved"; d.turbulence.printedDerivative = wall == "printed";
  d.wallTemperature = tWall;
  // The reference's bulk density (its rings are the profile's points; the mass per length is the
  // integral of rho 2 pi r dr by the midpoint rule on its own rings).
  double refMass = 0, refArea = 0;
  {
    const auto& r = profile.r;
    for (std::size_t j = 0; j < r.size(); ++j) {
      const double lo = j == 0 ? 0 : 0.5 * (r[j - 1] + r[j]), hi = j + 1 == r.size() ? radius : 0.5 * (r[j] + r[j + 1]);
      const double a = sq(hi) - sq(lo);
      refMass += profile.rho[j] * a; refArea += a;
    }
  }
  const double refRhoB = refMass / refArea;
  d.totalPressure = d.backPressure = refRhoB * gas * tWall; d.totalTemperature = tWall;
  Flow flow(d);
  flow.setThreads(threadCount);
  const auto& m = flow.mesh();
  std::vector<Primitive> cells(m.cells.size());
  std::vector<double> kOmega(2 * m.cells.size());
  double mass = 0, volume = 0;
  for (std::size_t q = 0; q < m.cells.size(); ++q) {
    const double r = m.cells[q].r;
    mass += interpolate(profile.r, profile.rho, r) * m.cells[q].volume; volume += m.cells[q].volume;
  }
  const double scale = refRhoB / (mass / volume);
  for (std::size_t q = 0; q < m.cells.size(); ++q) {
    const double r = m.cells[q].r, rho = scale * interpolate(profile.r, profile.rho, r), t = interpolate(profile.r, profile.t, r);
    cells[q] = {rho, interpolate(profile.r, profile.u, r), 0, rho * gas * t};
    kOmega[2 * q] = interpolate(profile.r, profile.k, r); kOmega[2 * q + 1] = interpolate(profile.r, profile.omega, r);
  }
  flow.setInitialState(cells);
  flow.setTurbulence(kOmega);
  flow.setBodyForce(std::vector<std::array<double, 2>>(m.cells.size(), {force, 0}));
  if (heating != 0) flow.setHeating(std::vector<double>(m.cells.size(), heating));
  const double pi = std::numbers::pi, wallArea = 2 * pi * radius * length;
  std::printf("pipe: %s, nr %d, stretching %g, force %g N/m^3, wall %g K, heating %g W/m^3, %d threads; reference %s (bulk density %.9f)\n",
              wall.c_str(), nr, stretching, force, tWall, heating, threadCount, path.c_str(), refRhoB);
  struct Sample { double ub, cf, tAxis, tb, st, balance; };
  auto sample = [&]() {
    const auto meas = flow.measurements();
    Sample s{};
    s.ub = meas.axialMomentum / meas.mass;
    double tAxis = 0, heat = 0, mw = 0;
    for (int i = 0; i < m.nz; ++i) tAxis += flow.temperature(m.index(i, 0)) / m.nz;
    for (std::size_t q = 0; q < m.cells.size(); ++q) { const double w = flow.state()[q][1] * m.cells[q].volume; heat += w * flow.temperature(q); mw += w; }
    s.tAxis = tAxis; s.tb = heat / mw;
    const double rhoB = meas.mass / (pi * radius * radius * length), tauW = -meas.wallAxialForce / wallArea;
    s.cf = 2 * tauW / (rhoB * s.ub * s.ub);
    const double qw = -meas.wallHeatFlow / wallArea;
    const double cpB = n2.medium.cv(s.tb, y) + gas;
    s.st = s.tb != tWall ? qw / (rhoB * cpB * s.ub * (s.tb - tWall)) : 0;
    s.balance = (meas.wallAxialForce + meas.bodyAxialForce) / meas.bodyAxialForce;
    return s;
  };
  const auto wall0 = std::chrono::steady_clock::now();
  std::vector<double> history;
  Sample s{};
  bool steady = false;
  for (int half = 1; 0.5 * half <= endMs + 1e-9; ++half) {
    flow.advanceTo(0.5e-3 * half);
    s = sample();
    history.push_back(s.ub);
    const double change = history.size() > 2 ? std::abs(s.ub / history[history.size() - 3] - 1) : 1;
    std::printf("  t %6.1f ms  u_b %.8f  c_f %.8f  T_axis %.5f  T_b %.5f  St %.6e  force balance %+.2e  du_b/ms %.2e  steps %llu  wall %.0f s\n",
                flow.time() * 1e3, s.ub, s.cf, s.tAxis, s.tb, s.st, s.balance, change,
                static_cast<unsigned long long>(flow.measurements().steps),
                std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count());
    std::fflush(stdout);
    if (change < 1e-5) { steady = true; break; }
  }
  const auto meas = flow.measurements();
  // First-cell y+ from the steady wall shear f R / 2 (wall units at the wall: rho_w = p_1 / (R T_w)).
  std::vector<double> diffusion(1), work;
  const double muW = n2.medium.transport(tWall, 1e5, y, diffusion.data(), work).viscosity;
  double yPlusLo = 1e300, yPlusHi = 0;
  for (int i = 0; i < m.nz; ++i) {
    const auto q = m.index(i, m.nr - 1);
    const double rhoW = flow.cellPrimitive(q).p / (gas * tWall), uTau = std::sqrt(force * radius / 2 / rhoW);
    const double yp = rhoW * uTau * flow.wallDistances()[q] / muW;
    yPlusLo = std::min(yPlusLo, yp); yPlusHi = std::max(yPlusHi, yp);
  }
  std::printf("result: %s nr %d stretching %g: %s at %.1f ms; u_b %.8f m/s, c_f %.8f, T_axis %.6f K, T_axis - T_wall %.6f K, T_b %.6f K, St %.8e;"
              " first-cell y+ %.3f to %.3f; budgets mass %.2e momentum %.2e energy %.2e; prescribed rho k V %.3e J; %llu steps, %.0f s\n",
              wall.c_str(), nr, stretching, steady ? "steady" : "NOT steady", flow.time() * 1e3, s.ub, s.cf, s.tAxis, s.tAxis - tWall, s.tb, s.st,
              yPlusLo, yPlusHi, meas.massBalanceError, meas.momentumBalanceError, meas.energyBalanceError, meas.prescribedTurbulentEnergy,
              static_cast<unsigned long long>(meas.steps), std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count());
  return 0;
}
// A priori: the law, given the reference profile's (u, y, T, p) at a point, against the reference's own
// wall shear (f R / 2, the force balance) and wall heat flux ((Q pi R^2 + f int u dA) / (2 pi R), the
// energy balance), with the engine's property choices (Flow::wallSolve). No march: this separates the
// law's form from the engine's discretisation. Reported beside it, the incompressible law (Gamma = beta
// = 0: Spalding's formula with the wall's density and viscosity).
int apriori(const Nitrogen& n2, int argc, char** argv) {
  if (argc < 6) throw std::invalid_argument("apriori <profile.csv> <T_wall> <force N/m^3> <heating W/m^3>");
  const auto profile = readProfile(argv[2]);
  const double tWall = std::atof(argv[3]), force = std::atof(argv[4]), heating = std::atof(argv[5]);
  const double radius = 5e-3, y[1] = {1}, gas = n2.medium.gasConstant(y), pi = std::numbers::pi;
  const auto& r = profile.r;
  const std::size_t n = r.size();
  double flowArea = 0;  // int u dA by the midpoint rule on the profile's own rings
  for (std::size_t j = 0; j < n; ++j) {
    const double lo = j == 0 ? 0 : 0.5 * (r[j - 1] + r[j]), hi = j + 1 == n ? radius : 0.5 * (r[j] + r[j + 1]);
    flowArea += profile.u[j] * pi * (sq(hi) - sq(lo));
  }
  const double tauRef = force * radius / 2, qRef = (heating * pi * radius * radius + force * flowArea) / (2 * pi * radius);
  // Wall pressure: p + 2/3 rho k is uniform across the reference, and k = 0 at the wall.
  const double pWall = profile.rho[n - 1] * gas * profile.t[n - 1] + 2.0 / 3 * profile.rho[n - 1] * profile.k[n - 1];
  std::vector<double> diffusion(1), work;
  const auto wall = n2.medium.transport(tWall, pWall, y, diffusion.data(), work);
  const auto atWall = n2.medium.properties(tWall, y);
  const double cpW = atWall.cv + atWall.r, rhoW = pWall / (gas * tWall), uTau = std::sqrt(tauRef / rhoW);
  const double prW = cpW * wall.viscosity / wall.conductivity;
  const wallLaw::Constants c{0.3697, 3.752};
  std::printf("a priori: %s; T_wall %g K, p_wall %.6e Pa, rho_w %.6f, mu_w %.6e, Pr_w %.4f; tau_w %.6f Pa, q_w %.6e W/m^2, Re_tau %.1f\n",
              argv[2], tWall, pWall, rhoW, wall.viscosity, prW, tauRef, qRef, rhoW * uTau * radius / wall.viscosity);
  std::printf("%9s %10s %8s %9s %12s %12s %14s %12s %12s\n", "y+", "u+", "T1/Tw", "rho1/rhow", "tau law-1", "q law-1",
              "tau Spalding-1", "law Gamma", "law beta");
  for (double target : {1.0, 5.0, 11.0, 30.0, 100.0, 300.0, 1000.0, 3000.0}) {
    std::size_t best = 0;
    for (std::size_t j = 0; j < n; ++j)
      if (std::abs((radius - r[j]) * rhoW * uTau / wall.viscosity - target) < std::abs((radius - r[best]) * rhoW * uTau / wall.viscosity - target)) best = j;
    const double y1 = radius - r[best], u1 = profile.u[best], t1 = profile.t[best], p1 = profile.rho[best] * gas * t1;
    const double rho1W = p1 / (gas * tWall);  // the engine's rho_w: the first cell's p at T_w
    const double cp = std::abs(t1 - tWall) > 1e-3 * tWall ? (n2.medium.enthalpy(t1, y) - n2.medium.enthalpy(tWall, y)) / (t1 - tWall) : cpW;
    const auto w1 = n2.medium.transport(tWall, p1, y, diffusion.data(), work);
    const auto law = wallLaw::solveIsothermal({u1, y1, t1, tWall, rho1W, w1.viscosity, w1.conductivity, cp, std::cbrt(cpW * w1.viscosity / w1.conductivity)}, c);
    const auto flat = wallLaw::solveIsothermal({u1, y1, tWall, tWall, rho1W, w1.viscosity, w1.conductivity, cp, 1e-300}, c);
    std::printf("%9.2f %10.4f %8.4f %9.4f %+12.4e %+12.4e %+14.4e %12.4e %12.4e\n", y1 * rhoW * uTau / wall.viscosity, u1 / uTau, t1 / tWall,
                profile.rho[best] / rhoW, law.shear / tauRef - 1, law.heat / qRef - 1, flat.shear / tauRef - 1, law.gamma, law.beta);
  }
  return 0;
}
}  // namespace

int main(int argc, char** argv) {
  try {
    const auto n2 = nitrogen();
    const std::string mode = argc > 1 ? argv[1] : "threads";
    if (mode == "threads") return threads(n2);
    if (mode == "pipe") return pipe(n2, argc, argv);
    if (mode == "apriori") return apriori(n2, argc, argv);
    std::fprintf(stderr, "usage: %s threads | pipe ...\n", argv[0]);
    return 2;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
