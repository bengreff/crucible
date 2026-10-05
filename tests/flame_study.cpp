// C2 flame-speed study: a freely propagating premixed laminar flame computed by the full engine
// (compressible flow, mixture-averaged transport, CVODES finite-rate chemistry with Strang
// splitting) against Cantera's steady free flame (Flow1D, same mechanism, same transport model).
// Mechanism: h2o2.yaml (verification use only). Mixture: H2/air at equivalence ratio 1
// ("H2:2, O2:1, N2:3.76"), 300 K, 1 atm.
//
// Case: a straight duct of radius 1 mm and length 6 mm, one cell high (nr = 1), closed at the
// injector end by an adiabatic slip plate, with slip adiabatic walls and the subsonic outlet at
// 101325 Pa. At t = 0 the gas is at rest at 1 atm: premix at 300 K for z < 4 mm, the premix's
// constant-(h, p) equilibrium for z >= 4 mm. The flame propagates toward the closed end into gas
// at rest; the burned gas leaves through the outlet. There is no igniter.
//
// Measurements (stated 4 October 2026, before the first run):
//  - Flame position x_f: the first crossing, from the closed end, of T = (300 + T_ad) / 2, by
//    linear interpolation between cell centres.
//  - Consumption speed S_c = (-int_0^{x*} W w_H2 dz + j_H2(x*)) / (rho_u (Y_H2,u - Y_H2(x*))),
//    with x* = x_f + 1 mm, the rates evaluated by Cantera on the cell states, and j_H2 the
//    mixture-averaged diffusive flux (with the correction velocity) between the two points that
//    bracket x*. This is the H2 balance of the region ahead of x*, so for a steadily propagating
//    planar flame it equals the laminar flame speed.
//  - Displacement speed S_d = -dx_f/dt: least-squares slope of x_f over the last 0.2 ms (the gas
//    ahead of the flame is at rest against the closed end).
// Reference: Cantera Flow1D free flame, mixture-averaged transport, energy on, 3 cm domain,
// refined to (ratio 3, slope 0.06, curve 0.12), then (2, 0.02, 0.04), then (2, 0.01, 0.02);
// S_L = inlet velocity of the finest solution.
//
// Criteria (stated before the first run):
//  1. On the finest grid (dz = 10 um) S_c is within 1% of S_L.
//  2. |S_c / S_L - 1| falls at each refinement, dz = 40, 20, 10 um.
//  3. On the finest grid S_d is within 1% of S_L.
//  4. Settled: the means of S_c over the two halves of the last 0.2 ms differ by less than 0.2%.
//  Preconditions on the reference (if they fail, the comparison is not meaningful as stated):
//  the two finest refinement levels give S_L within 0.2% of each other, and the free flame's own
//  S_c (same balance, same x* offset) is within 0.2% of its S_L.
//  Mass and energy budgets are reported, not judged (C1 and C2 step 1 verify them).
//
// Change after the smoke test (4 October 2026, before any engine flame result; criteria
// unchanged): the first definition of S_c had no j_H2 term and the reference had two levels. On
// the reference alone, S_c without the flux came out 0.305% below S_L (2.325432 against
// 2.332536 m/s, 373 points), failing the precondition: H2 is still recombining 1 mm behind the
// flame and its diffusive flux there is not negligible. The flux term and the third level were
// added. The smoke test ran the engine to 0.02 ms only, before any flame had formed.
//
// Change after the first 40 um run (4 October 2026; criteria unchanged). With the fixed-pressure
// outlet the duct is a closed-open acoustic resonator (closed plate, pressure node at the outlet,
// about 20 kHz), and the flame drove it: pressure 92 to 114 kPa and S_c swinging 2.0 to 3.0 m/s
// through 1 ms, S_c +20.4% and drift 3.1% at the end (docs/evidence/FLAME_C2.md keeps that run).
// The reference is a flame in open space, so the outlet is now partially non-reflecting
// (Definition::outletRelaxation = 0.25, Poinsot and Lele; core test outletChecks): derived
// reflection about 0.13 at 20 kHz with burned gas at the outlet, steady pressure still 1 atm.
// The ambient is now the burned gas (composition and T_ad), so backflow cannot draw cold gas.
// The finer grids had not been run.
// Second change (same day; criteria unchanged). The first non-reflecting 40 um run blew up at
// 55 us: that outlet used isentropic Riemann invariants, which turn an entropy gradient at the
// outlet into an acoustic wave (measured in a perfect-gas duct: a 1 K hot spot grew to about
// 870 kPa). The outlet now works in W = p +- rho a u_z (core check: a 100 K hot spot leaves with
// a pressure change of 5.8e-4 of ambient). The blown-up run is kept in FLAME_C2.md.
//
// Usage: crucible_flame_study <dz um> <end time s> <threads> <output prefix>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "adapters/thermo.hpp"
#include "cantera/core.h"
#include "cantera/oneD/DomainFactory.h"
#include "cantera/onedim.h"
#include "core/flow.hpp"

using namespace crucible;

namespace {

constexpr double kPi = std::numbers::pi;
constexpr double kP = 101325.0, kTu = 300.0, kRadius = 1e-3, kLength = 6e-3, kStep = 4e-3, kOffset = 1e-3;
constexpr const char* kMixture = "H2:2, O2:1, N2:3.76";
constexpr const char* kMechanism = "h2o2.yaml";

double crossing(const std::vector<double>& z, const std::vector<double>& t, double level) {
  for (std::size_t i = 1; i < z.size(); ++i)
    if (t[i - 1] < level && t[i] >= level) return z[i - 1] + (level - t[i - 1]) / (t[i] - t[i - 1]) * (z[i] - z[i - 1]);
  return NAN;
}

double interpolate(const std::vector<double>& z, const std::vector<double>& v, double x) {
  for (std::size_t i = 1; i < z.size(); ++i)
    if (z[i] >= x) return v[i - 1] + (x - z[i - 1]) / (z[i] - z[i - 1]) * (v[i] - v[i - 1]);
  return v.back();
}

// Mixture-averaged diffusive mass flux of species k between two states a and b at positions za,
// zb (face averages of rho, W, D_km and Y; mole-fraction differences; correction velocity).
double diffusiveFlux(std::size_t k, const std::vector<double>& molarMass, double za, double zb, double rhoA, double rhoB,
                     const std::vector<double>& ya, const std::vector<double>& yb, const std::vector<double>& da,
                     const std::vector<double>& db) {
  const std::size_t ns = molarMass.size();
  double ma = 0, mb = 0, mf = 0;
  for (std::size_t l = 0; l < ns; ++l) {
    ma += ya[l] / molarMass[l];
    mb += yb[l] / molarMass[l];
    mf += 0.5 * (ya[l] + yb[l]) / molarMass[l];
  }
  const double rho = 0.5 * (rhoA + rhoB);
  double jk = 0, sum = 0;
  for (std::size_t l = 0; l < ns; ++l) {
    const double grad = (yb[l] / molarMass[l] / mb - ya[l] / molarMass[l] / ma) / (zb - za);
    const double j = -rho * molarMass[l] * mf * 0.5 * (da[l] + db[l]) * grad;
    sum += j;
    if (l == k) jk = j;
  }
  return jk - 0.5 * (ya[k] + yb[k]) * sum;
}

struct Reference {
  double speed, consumption, tad;
  std::size_t points;
  std::vector<double> z, t, yH2, yH, yOH;
};

// Steady free flame; the flame position is pinned by the fixed-temperature point.
std::vector<Reference> freeFlame(double tad, std::size_t iH2, double& rhoU, double& yH2u) {
  auto sol = Cantera::newSolution(kMechanism, "", "mixture-averaged");
  auto gas = sol->thermo();
  const std::size_t ns = gas->nSpecies();
  gas->setState_TPX(kTu, kP, kMixture);
  rhoU = gas->density();
  std::vector<double> x(ns), yin(ns), yout(ns);
  gas->getMoleFractions(x.data());
  gas->getMassFractions(yin.data());
  yH2u = yin[iH2];
  gas->equilibrate("HP");
  gas->getMassFractions(yout.data());
  const double rhoOut = gas->density();

  auto flow = Cantera::newDomain<Cantera::Flow1D>("gas-flow", sol, "flow");
  flow->setFreeFlow();
  const double width = 0.03;
  std::vector<double> grid(8);
  for (std::size_t i = 0; i < grid.size(); ++i) grid[i] = width * i / (grid.size() - 1);
  flow->setupGrid(grid.size(), grid.data());
  flow->setSteadyTolerances(1e-6, 1e-12);
  flow->setTransientTolerances(1e-6, 1e-12);
  auto inlet = Cantera::newDomain<Cantera::Inlet1D>("inlet", sol);
  const double uin = 2.0;
  inlet->setMoleFractions(x.data());
  inlet->setMdot(uin * rhoU);
  inlet->setTemperature(kTu);
  auto outlet = Cantera::newDomain<Cantera::Outlet1D>("outlet", sol);
  std::vector<std::shared_ptr<Cantera::Domain1D>> domains{inlet, flow, outlet};
  Cantera::Sim1D flame(domains);
  const std::vector<double> locs{0.0, 0.3, 0.5, 1.0};
  flow->setProfile("velocity", locs, {uin, uin, uin * rhoU / rhoOut, uin * rhoU / rhoOut});
  flow->setProfile("T", locs, {kTu, kTu, tad, tad});
  for (std::size_t k = 0; k < ns; ++k)
    flow->setProfile(gas->speciesName(k), locs, {yin[k], yin[k], yout[k], yout[k]});
  flame.setFixedTemperature(0.5 * (kTu + tad));
  flow->setEnergyEnabled(true);

  std::vector<Reference> out;
  auto rates = Cantera::newSolution(kMechanism, "", "none");
  std::vector<double> molarMass(ns);
  for (std::size_t k = 0; k < ns; ++k) molarMass[k] = gas->molecularWeight(k);
  for (auto [ratio, slope, curve] : {std::array<double, 3>{3, 0.06, 0.12}, std::array<double, 3>{2, 0.02, 0.04},
                                     std::array<double, 3>{2, 0.01, 0.02}}) {
    flow->setRefineCriteria(ratio, slope, curve);
    flame.solve(0, true);
    Reference r;
    r.tad = tad;
    r.z = flow->grid();
    r.t = flow->values("T");
    r.yH2 = flow->values("H2");
    r.yH = flow->values("H");
    r.yOH = flow->values("OH");
    r.points = r.z.size();
    r.speed = flow->values("velocity")[0];
    // Same consumption integral as the engine, trapezoid rule up to x* (linear in the last interval).
    std::vector<std::vector<double>> species(ns);
    for (std::size_t k = 0; k < ns; ++k) species[k] = flow->values(gas->speciesName(k));
    std::vector<double> c(r.points), y(ns), wdot(ns);
    for (std::size_t i = 0; i < r.points; ++i) {
      for (std::size_t k = 0; k < ns; ++k) y[k] = species[k][i];
      rates->thermo()->setState_TPY(r.t[i], kP, y.data());
      rates->kinetics()->getNetProductionRates(wdot.data());
      c[i] = -wdot[iH2] * rates->thermo()->molecularWeight(iH2);
    }
    const double xf = crossing(r.z, r.t, 0.5 * (kTu + tad)), xs = xf + kOffset;
    double integral = 0;
    for (std::size_t i = 1; i < r.points && r.z[i - 1] < xs; ++i) {
      const double b = std::min(r.z[i], xs), cb = interpolate(r.z, c, b);
      integral += 0.5 * (c[i - 1] + cb) * (b - r.z[i - 1]);
    }
    std::size_t a = 0;
    while (r.z[a + 1] < xs) ++a;
    std::vector<double> ya(ns), yb(ns), da(ns), db(ns);
    double rho[2];
    for (int side = 0; side < 2; ++side) {
      auto& yy = side ? yb : ya;
      for (std::size_t k = 0; k < ns; ++k) yy[k] = species[k][a + side];
      sol->thermo()->setState_TPY(r.t[a + side], kP, yy.data());
      rho[side] = sol->thermo()->density();
      sol->transport()->getMixDiffCoeffs((side ? db : da).data());
    }
    const double flux = diffusiveFlux(iH2, molarMass, r.z[a], r.z[a + 1], rho[0], rho[1], ya, yb, da, db);
    r.consumption = (integral + flux) / (rhoU * (yH2u - interpolate(r.z, r.yH2, xs)));
    out.push_back(std::move(r));
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 5) {
    std::fprintf(stderr, "usage: %s <dz um> <end s> <threads> <prefix>\n", argv[0]);
    return 2;
  }
  const double dz = std::atof(argv[1]) * 1e-6, end = std::atof(argv[2]);
  const int threads = std::atoi(argv[3]);
  const std::string prefix = argv[4];
  const int nz = static_cast<int>(std::lround(kLength / dz));
  auto clock = std::chrono::steady_clock::now();
  auto wall = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - clock).count(); };

  thermo::ReactionSource source(kMechanism);
  const auto medium = source.medium();
  const std::size_t ns = medium.size();
  std::size_t iH2 = ns;
  for (std::size_t k = 0; k < ns; ++k)
    if (medium.species()[k].name == "H2") iH2 = k;
  const auto yu = source.massFractions(kMixture);

  // Products: the premix's constant-(h, p) equilibrium.
  auto eq = Cantera::newSolution(kMechanism, "", "none");
  eq->thermo()->setState_TPY(kTu, kP, yu.data());
  eq->thermo()->equilibrate("HP");
  const double tad = eq->thermo()->temperature();
  std::vector<double> yb(ns);
  eq->thermo()->getMassFractions(yb.data());

  double rhoU = 0, yH2u = 0;
  const auto reference = freeFlame(tad, iH2, rhoU, yH2u);
  const Reference& ref = reference.back();
  std::printf("T_ad %.3f K; reference wall %.1f s\n", tad, wall());
  for (const auto& r : reference)
    std::printf("free flame: S_L %.6f m/s, S_c %.6f m/s (%+.3f%%), %zu points\n", r.speed, r.consumption,
                100 * (r.consumption / r.speed - 1), r.points);

  Definition d;
  d.experiment = Case::Chamber;
  d.nz = nz;
  d.nr = 1;
  d.contour = {{0.0, kRadius}, {kLength, kRadius}};
  d.species = medium.species();
  // Ambient = burned gas, so any backflow at the outlet draws products; the outlet is partially
  // non-reflecting (see the header). The initial fill is set per cell below.
  d.composition = yb;
  d.backPressure = kP;
  d.ambientTemperature = tad;
  d.outletRelaxation = 0.25;
  d.transport = thermo::transportFits(kMechanism);
  d.wallSlip = true;
  d.wallTemperature = 0;
  Flow flow(d);
  const auto& mesh = flow.mesh();
  std::vector<Primitive> cells(mesh.cells.size());
  std::vector<double> fractions(mesh.cells.size() * ns);
  for (int i = 0; i < nz; ++i) {
    const auto q = mesh.index(i, 0);
    const bool burnt = mesh.cells[q].z >= kStep;
    const auto& y = burnt ? yb : yu;
    const double t = burnt ? tad : kTu;
    cells[q] = {kP / (medium.gasConstant(y.data()) * t), 0, 0, kP};
    std::copy(y.begin(), y.end(), fractions.begin() + q * ns);
  }
  flow.setInitialState(cells, fractions);
  thermo::ReactingFlow reacting(flow, kMechanism, threads, 1e-6, 1e-12, thermo::Chemistry::FiniteRate);
  double area = 0;
  for (const auto& c : mesh.cells) area += c.volume;
  area /= kLength;
  std::printf("mesh %d x 1: dz %.2f um, area %.6e m^2 (pi R^2 %.6e)\n", nz, mesh.dz * 1e6, area, kPi * kRadius * kRadius);

  std::vector<double> zc(nz), temp(nz), yH2(nz), cons(nz), z(ns + 1), dzdt(ns + 1), molarMass(ns), work;
  for (std::size_t k = 0; k < ns; ++k) molarMass[k] = medium.species()[k].molarMass;
  auto transportMedium = medium;
  transportMedium.setTransport(thermo::transportFits(kMechanism));
  for (int i = 0; i < nz; ++i) zc[i] = mesh.cells[mesh.index(i, 0)].z;
  auto measure = [&](double& xf, double& sc) {
    for (int i = 0; i < nz; ++i) {
      const auto q = mesh.index(i, 0);
      const auto y = flow.massFractions(q);
      temp[i] = flow.temperature(q);
      yH2[i] = y[iH2];
      z[0] = temp[i];
      std::copy(y.begin(), y.end(), z.begin() + 1);
      const double rho = flow.state()[q][0];
      source.rates(rho, z.data(), dzdt.data());
      cons[i] = -rho * dzdt[iH2 + 1];
    }
    xf = crossing(zc, temp, 0.5 * (kTu + tad));
    const double xs = xf + kOffset;
    double integral = 0;
    for (int i = 0; i < nz; ++i) {
      const double lo = zc[i] - 0.5 * mesh.dz;
      integral += cons[i] * std::clamp((xs - lo) / mesh.dz, 0.0, 1.0) * mesh.dz;
    }
    int a = 0;
    while (a + 2 < nz && zc[a + 1] < xs) ++a;
    std::vector<double> ya = flow.massFractions(mesh.index(a, 0)), yb = flow.massFractions(mesh.index(a + 1, 0));
    std::vector<double> da(ns), db(ns);
    transportMedium.transport(temp[a], flow.cellPrimitive(mesh.index(a, 0)).p, ya.data(), da.data(), work);
    transportMedium.transport(temp[a + 1], flow.cellPrimitive(mesh.index(a + 1, 0)).p, yb.data(), db.data(), work);
    const double flux = diffusiveFlux(iH2, molarMass, zc[a], zc[a + 1], flow.state()[mesh.index(a, 0)][0],
                                      flow.state()[mesh.index(a + 1, 0)][0], ya, yb, da, db);
    sc = (integral + flux) / (rhoU * (yH2u - interpolate(zc, yH2, xs)));
  };

  FILE* history = std::fopen((prefix + "_history.csv").c_str(), "w");
  std::fprintf(history, "t,x_f,S_c,T_max,p_min,p_max,u_ahead,mass_budget,energy_budget,steps,wall\n");
  std::vector<double> times, fronts, speeds;
  const double every = 5e-6;
  for (int n = 1; n * every <= end + 1e-12; ++n) {
    reacting.advanceTo(n * every);
    double xf, sc;
    measure(xf, sc);
    double tmax = 0, pmin = 1e30, pmax = 0;
    for (int i = 0; i < nz; ++i) {
      const auto w = flow.cellPrimitive(mesh.index(i, 0));
      tmax = std::max(tmax, temp[i]);
      pmin = std::min(pmin, w.p);
      pmax = std::max(pmax, w.p);
    }
    const double uAhead = std::isfinite(xf) ? flow.cellPrimitive(mesh.index(std::max(0, static_cast<int>((xf - 5e-4) / mesh.dz)), 0)).uz : NAN;
    const auto m = flow.measurements();
    times.push_back(flow.time());
    fronts.push_back(xf);
    speeds.push_back(sc);
    std::fprintf(history, "%.9e,%.9e,%.9e,%.3f,%.3f,%.3f,%.6e,%.3e,%.3e,%ld,%.1f\n", flow.time(), xf, sc, tmax, pmin,
                 pmax, uAhead, m.massBalanceError, m.energyBalanceError, reacting.stats().steps, wall());
    std::fflush(history);
    if (n % 20 == 0)
      std::printf("t %.3f ms  x_f %.4f mm  S_c %.5f m/s  T_max %.1f K  steps %ld  wall %.0f s\n", flow.time() * 1e3,
                  xf * 1e3, sc, tmax, reacting.stats().steps, wall());
    std::fflush(stdout);
  }
  std::fclose(history);

  // Settling and displacement speed over the last 0.2 ms.
  const double window = 2e-4, tEnd = times.back();
  double first = 0, second = 0, st = 0, sx = 0, stt = 0, stx = 0;
  int nFirst = 0, nSecond = 0, nFit = 0;
  for (std::size_t k = 0; k < times.size(); ++k) {
    if (times[k] < tEnd - window - 1e-12) continue;
    if (times[k] < tEnd - window / 2) { first += speeds[k]; ++nFirst; } else { second += speeds[k]; ++nSecond; }
    st += times[k]; sx += fronts[k]; stt += times[k] * times[k]; stx += times[k] * fronts[k]; ++nFit;
  }
  first /= nFirst;
  second /= nSecond;
  const double sd = -(nFit * stx - st * sx) / (nFit * stt - st * st);
  const double sc = second, drift = std::abs(second / first - 1);

  {
    FILE* f = std::fopen((prefix + "_profile.csv").c_str(), "w");
    std::fprintf(f, "source,z_minus_xf,T,Y_H2,Y_H,Y_OH,p,u_z\n");
    const double xf = fronts.back();
    std::size_t iH = ns, iOH = ns;
    for (std::size_t k = 0; k < ns; ++k) {
      if (medium.species()[k].name == "H") iH = k;
      if (medium.species()[k].name == "OH") iOH = k;
    }
    for (int i = 0; i < nz; ++i) {
      const auto y = flow.massFractions(mesh.index(i, 0));
      const auto w = flow.cellPrimitive(mesh.index(i, 0));
      std::fprintf(f, "engine,%.6e,%.6e,%.6e,%.6e,%.6e,%.6f,%.6e\n", zc[i] - xf, temp[i], y[iH2], y[iH], y[iOH], w.p, w.uz);
    }
    const double xr = crossing(ref.z, ref.t, 0.5 * (kTu + tad));
    for (std::size_t i = 0; i < ref.points; ++i)
      std::fprintf(f, "freeflame,%.6e,%.6e,%.6e,%.6e,%.6e,%.6f,nan\n", ref.z[i] - xr, ref.t[i], ref.yH2[i], ref.yH[i], ref.yOH[i], kP);
    std::fclose(f);
  }

  const auto m = flow.measurements();
  std::printf("end %.3f ms after %ld steps, wall %.0f s; mass budget %.3e, energy budget %.3e\n", tEnd * 1e3,
              reacting.stats().steps, wall(), m.massBalanceError, m.energyBalanceError);
  std::printf("dz %.2f um: S_c %.6f m/s (%+.3f%%), S_d %.6f m/s (%+.3f%%), S_c drift over the last 0.2 ms %.3f%%; "
              "S_L %.6f m/s\n",
              mesh.dz * 1e6, sc, 100 * (sc / ref.speed - 1), sd, 100 * (sd / ref.speed - 1), 100 * drift, ref.speed);
  std::printf("reference: levels differ %.3f%%, free-flame S_c vs S_L %+.3f%%\n",
              100 * std::abs(reference[reference.size() - 2].speed / ref.speed - 1), 100 * (ref.consumption / ref.speed - 1));
  return 0;
}
