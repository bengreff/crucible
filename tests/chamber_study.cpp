// C1 chamber study: a premixed gaseous H2/O2 thrust chamber started "valves open, ignite"
// (TECHNICAL_PLAN, Dated plan to January). The chamber starts full of ambient N2 at rest, the
// valve ramps the supply open, the igniter fires, and the run marches in physical time until the
// chamber has settled. Mechanism: h2o2.yaml (verification use only, not a rocket-pressure model).
//
// Geometry (built here as a contour table): throat radius 10 mm, chamber radius 25 mm, circular
// throat arcs of radius 20 mm on both sides (R = r_c / r_t = 2), 30 degree convergent cone,
// 20 mm cylinder-to-cone blend, 15 degree conical exit to area ratio 10. The cylinder length puts
// the throat at mid-length, so an even nz places a mesh station exactly on the throat.
//
// Comparison with the one-dimensional ideal rocket (CEA problem, Mixture::idealRocket, same
// thermo) in LocalEquilibrium mode, with corrections stated before any run
// (docs/evidence/CHAMBER_C1.md):
//  - Chamber pressure: mass-flux-averaged equilibrium stagnation pressure p0 of the last cylinder
//    column (equilibrium isentropic stagnation from each cell's state).
//  - c*: two-dimensional throat flow passes Cd times the one-dimensional flow, so
//    c*_2D = p0 A_t / mdot = c*_1D(p0) / Cd. Cd from the Kliegel-Levine series (AIAA J 7(7), 1969),
//    coefficients as tabulated by Johnson and Wright (J. Fluids Eng. 130, 071202, 2008, Table 4):
//    Cd = 1 - a2/L^2 + a3/L^3 - a4/L^4, L = 1 + R, a2 = (g+1)/96, a3 = (g+1)(8g-27)/2304,
//    a4 = (g+1)(754g^2 - 757g + 3633)/276480.
//  - Vacuum Isp: conical divergence factor lambda = (1 + cos 15 deg)/2 on the momentum part,
//    Isp = lambda v_e + p_e eps c*/p0 with v_e = Isp_1D - p_e eps c*/p0.
//  - Vacuum thrust: device thrust plus ambient pressure times exit area (valid with a supersonic
//    exit, checked).
//
// Usage: crucible_chamber_study <eq|fr|frozen> <nz> <nr> <end time s> <threads> <output prefix>
// Writes <prefix>_history.csv (every 2 us) and <prefix>_field_<us>.csv snapshots.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "adapters/thermo.hpp"
#include "cantera/core.h"

using namespace crucible;

namespace {
constexpr double kPi = std::numbers::pi, kG0 = 9.80665;
constexpr double kRt = 0.010, kRc = 0.025, kArc = 0.020, kBlend = 0.020;
const double kAlpha = 30 * kPi / 180, kTheta = 15 * kPi / 180, kEps = 10;
constexpr double kMixtureRatio = 5, kMassFlow = 0.4, kAmbientP = 1000, kAmbientT = 300;

struct Geometry {
  std::vector<std::array<double, 2>> contour;
  double cylinderEnd = 0;  // z where the cylinder meets the blend arc
  double throat = 0;
};

// Contour with the throat at z = 0, then shifted so the injector face is z = 0.
Geometry geometry() {
  const double re = kRt * std::sqrt(kEps);
  const double z1 = -kArc * std::sin(kAlpha), r1 = kRt + kArc * (1 - std::cos(kAlpha));
  const double rb = kRc - kBlend * (1 - std::cos(kAlpha)), zb = z1 - (rb - r1) / std::tan(kAlpha);
  const double zc = zb - kBlend * std::sin(kAlpha);
  const double z2 = kArc * std::sin(kTheta), r2 = kRt + kArc * (1 - std::cos(kTheta));
  const double ze = z2 + (re - r2) / std::tan(kTheta);
  const double z0 = -ze;  // throat at mid-length
  std::vector<std::array<double, 2>> c;
  auto add = [&](double z, double r) {
    if (c.empty() || z > c.back()[0] + 1e-12) c.push_back({z, r});
  };
  const int n = 400;
  add(z0, kRc);
  add(zc, kRc);
  for (int k = 1; k <= n; ++k) {  // blend arc, centre (zc, kRc - kBlend)
    double a = kAlpha * k / n;
    add(zc + kBlend * std::sin(a), kRc - kBlend + kBlend * std::cos(a));
  }
  add(z1, r1);  // convergent cone (straight)
  for (int k = 1; k <= 2 * n; ++k) {  // throat arc, centre (0, kRt + kArc)
    double a = -kAlpha + (kAlpha + kTheta) * k / (2 * n);
    add(kArc * std::sin(a), kRt + kArc - kArc * std::cos(a));
  }
  add(ze, re);  // conical exit
  Geometry g;
  for (auto& p : c) g.contour.push_back({p[0] - z0, p[1]});
  g.cylinderEnd = zc - z0;
  g.throat = -z0;
  return g;
}

double dischargeCoefficient(double R, double g) {
  const double L = 1 + R;
  const double a2 = (g + 1) / 96, a3 = (g + 1) * (8 * g - 27) / 2304;
  const double a4 = (g + 1) * (754 * g * g - 757 * g + 3633) / 276480;
  return 1 - a2 / (L * L) + a3 / (L * L * L) - a4 / (L * L * L * L);
}

// Equilibrium isentropic stagnation pressure of a cell: hold s, raise p until h = h + u^2/2.
double stagnationPressure(Cantera::ThermoPhase& gas, double t, double p, const double* y, double kinetic) {
  gas.setState_TPY(t, p, y);
  gas.equilibrate("TP");
  const double s = gas.entropy_mass(), h0 = gas.enthalpy_mass() + kinetic;
  double p0 = p;
  for (int n = 0; n < 50; ++n) {
    gas.setState_SP(s, p0);
    gas.equilibrate("SP");
    double step = (h0 - gas.enthalpy_mass()) / (p0 / gas.density());
    p0 *= std::exp(step);
    if (std::abs(step) < 1e-13) break;
  }
  return p0;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 7) {
    std::fprintf(stderr, "usage: %s <eq|fr|frozen> <nz> <nr> <end s> <threads> <prefix>\n", argv[0]);
    return 2;
  }
  const std::string mode = argv[1], prefix = argv[6];
  const int nz = std::atoi(argv[2]), nr = std::atoi(argv[3]), threads = std::atoi(argv[5]);
  const double end = std::atof(argv[4]);
  const auto chemistry = mode == "eq" ? thermo::Chemistry::LocalEquilibrium
                         : mode == "fr" ? thermo::Chemistry::FiniteRate
                                        : thermo::Chemistry::Frozen;
  thermo::ReactionSource source("h2o2.yaml");
  const auto geo = geometry();

  Definition d;
  d.experiment = Case::Chamber;
  d.nz = nz;
  d.nr = nr;
  d.contour = geo.contour;
  d.species = source.medium().species();
  d.composition = source.massFractions("N2:1");
  d.backPressure = kAmbientP;
  d.ambientTemperature = kAmbientT;
  Supply s;
  s.outerRadius = kRc;
  s.massFlow = kMassFlow;
  s.totalTemperature = 300;
  auto h2 = source.massFractions("H2:1"), o2 = source.massFractions("O2:1");
  s.composition.resize(h2.size());
  for (std::size_t k = 0; k < h2.size(); ++k)
    s.composition[k] = (h2[k] + kMixtureRatio * o2[k]) / (1 + kMixtureRatio);
  s.opens = 0;
  s.ramp = 5e-4;
  d.supplies = {s};
  // Igniter: 0.5 J over 0.2 ms from 0.2 ms, in a 10 mm by 10 mm core 5 mm off the injector face.
  // Local equilibrium burns any premixed gas at once, so in "eq" mode it has no role.
  d.igniter = {0.005, 0.015, 0.010, 0.5, 2e-4, 2e-4};

  Flow flow(d);
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", threads, 1e-6, 1e-12, chemistry);
  const auto& mesh = flow.mesh();
  double rMin = 1e9;
  int throatStation = 0;
  for (int i = 0; i <= nz; ++i)
    if (mesh.radius[i] < rMin) { rMin = mesh.radius[i]; throatStation = i; }
  const double at = kPi * rMin * rMin, ae = kPi * mesh.radius[nz] * mesh.radius[nz], eps = ae / at;
  int endColumn = 0;
  for (int i = 0; i < nz; ++i)
    if (mesh.cells[mesh.index(i, 0)].z < geo.cylinderEnd) endColumn = i;
  std::printf("mesh %dx%d: dz %.3f mm, throat station %d at z %.4f m (design %.4f), r_t %.5f m, eps %.4f, "
              "chamber-end column %d at z %.4f m\n",
              nz, nr, mesh.dz * 1e3, throatStation, throatStation * mesh.dz, geo.throat, rMin, eps, endColumn,
              mesh.cells[mesh.index(endColumn, 0)].z);

  const std::size_t ns = flow.medium().size();
  auto speciesIndex = [&](const char* name) {
    for (std::size_t k = 0; k < ns; ++k)
      if (flow.medium().species()[k].name == name) return k;
    return ns;
  };
  const std::size_t iH2O = speciesIndex("H2O"), iN2 = speciesIndex("N2"), iO2 = speciesIndex("O2");

  auto writeField = [&](double t) {
    char name[512];
    std::snprintf(name, sizeof name, "%s_field_%06.0f.csv", prefix.c_str(), t * 1e6);
    FILE* f = std::fopen(name, "w");
    std::fprintf(f, "i,j,z,r,T,p,mach,uz,ur,Y_H2O,Y_O2,Y_N2\n");
    for (int i = 0; i < nz; ++i)
      for (int j = 0; j < nr; ++j) {
        auto q = mesh.index(i, j);
        auto w = flow.cellPrimitive(q);
        auto y = flow.massFractions(q);
        double temp = flow.temperature(q);
        double a = flow.medium().soundSpeed(temp, y.data());
        std::fprintf(f, "%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", i, j, mesh.cells[q].z,
                     mesh.cells[q].r, temp, w.p, std::hypot(w.uz, w.ur) / a, w.uz, w.ur, y[iH2O], y[iO2], y[iN2]);
      }
    std::fclose(f);
  };

  FILE* history = std::fopen((prefix + "_history.csv").c_str(), "w");
  std::fprintf(history,
               "t,supply,inlet,outlet,p_injector,p_chamber_end,F_vac,isp_vac_supply,exit_mach_min,T_max,"
               "igniter_energy,mass_budget,energy_budget,steps\n");
  std::vector<double> snapshots = {2e-5, 5e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 7.5e-4, 1e-3, 1.5e-3, 2e-3, 3e-3, 4e-3, 6e-3, 8e-3};
  std::size_t nextSnapshot = 0;
  const auto clock0 = std::chrono::steady_clock::now();
  double nextSample = 0;
  struct Row { double t, outlet, pInj, fVac; };
  std::vector<Row> rows;
  auto sample = [&]() {
    auto m = flow.measurements();
    double exitMin = 1e9, tMax = 0, pEnd = 0, area = 0;
    for (int j = 0; j < nr; ++j) {
      auto q = mesh.index(nz - 1, j);
      auto w = flow.cellPrimitive(q);
      auto y = flow.massFractions(q);
      exitMin = std::min(exitMin, std::hypot(w.uz, w.ur) / flow.medium().soundSpeed(flow.temperature(q), y.data()));
      double a = mesh.axialArea(endColumn + 1, j);
      pEnd += a * flow.cellPrimitive(mesh.index(endColumn, j)).p;
      area += a;
    }
    for (std::size_t q = 0; q < flow.state().size(); ++q) tMax = std::max(tMax, flow.temperature(q));
    const double fVac = m.deviceThrust + m.ambientAxialForce;
    std::fprintf(history, "%.9e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.5f,%.2f,%.6e,%.3e,%.3e,%llu\n", m.time,
                 m.supplyMassFlow, m.inletMassFlow, m.outletMassFlow, m.injectorPressure, pEnd / area, fVac,
                 m.supplyMassFlow > 0 ? fVac / (m.supplyMassFlow * kG0) : 0.0, exitMin, tMax, m.igniterEnergy,
                 m.massBalanceError, m.energyBalanceError, static_cast<unsigned long long>(m.steps));
    rows.push_back({m.time, m.outletMassFlow, m.injectorPressure, fVac});
  };
  writeField(0);
  while (flow.time() < end) {
    double target = std::min({end, nextSample, nextSnapshot < snapshots.size() ? snapshots[nextSnapshot] : end});
    if (target <= flow.time()) target = std::min(end, flow.time() + 2e-6);
    reacting.advanceTo(target);
    if (flow.time() >= nextSample) {
      sample();
      nextSample += 2e-6;
      if (rows.size() % 250 == 0) {
        double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();
        std::printf("t %.3f ms  p_inj %.4f MPa  outlet %.4f kg/s  F_vac %.1f N  steps %ld  wall %.0f s\n",
                    flow.time() * 1e3, rows.back().pInj / 1e6, rows.back().outlet, rows.back().fVac,
                    reacting.stats().steps, wall);
        std::fflush(stdout);
        std::fflush(history);
      }
    }
    if (nextSnapshot < snapshots.size() && flow.time() >= snapshots[nextSnapshot]) {
      writeField(flow.time());
      ++nextSnapshot;
    }
  }
  writeField(flow.time());
  std::fclose(history);
  if (std::getenv("CHEMPROFILE")) {
    // Time one equilibrium call per cell on the final state and list the slowest cells.
    struct Slow { double us, t, p, yH2O, yO2, yN2; std::size_t q; };
    std::vector<Slow> slow;
    std::vector<double> z(ns + 1);
    for (std::size_t q = 0; q < flow.state().size(); ++q) {
      auto y = flow.massFractions(q);
      z[0] = flow.temperature(q);
      std::copy(y.begin(), y.end(), z.begin() + 1);
      auto t0 = std::chrono::steady_clock::now();
      try { source.equilibrateUV(flow.state()[q][0], z.data()); } catch (const std::exception&) { z[0] = -1; }
      double us = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() * 1e6;
      slow.push_back({us, flow.temperature(q), flow.cellPrimitive(q).p, y[iH2O], y[iO2], y[iN2], q});
    }
    std::sort(slow.begin(), slow.end(), [](const Slow& a, const Slow& b) { return a.us > b.us; });
    double total = 0;
    for (const auto& c : slow) total += c.us;
    std::printf("equilibrium cost on final state: total %.0f ms over %zu cells\n", total / 1e3, slow.size());
    for (std::size_t k = 0; k < std::min<std::size_t>(12, slow.size()); ++k)
      std::printf("  %.0f us  cell (%zu,%zu)  T %.1f K  p %.0f Pa  Y_H2O %.3e  Y_O2 %.3e  Y_N2 %.3e\n", slow[k].us,
                  slow[k].q / nr, slow[k].q % nr, slow[k].t, slow[k].p, slow[k].yH2O, slow[k].yO2, slow[k].yN2);
  }
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();

  // Settling over the last millisecond (declared criterion: relative drift below 1e-3).
  auto drift = [&](auto get) {
    double lo = 1e300, hi = -1e300, last = get(rows.back());
    for (const auto& r : rows)
      if (r.t >= rows.back().t - 1e-3) { lo = std::min(lo, get(r)); hi = std::max(hi, get(r)); }
    return (hi - lo) / std::abs(last);
  };
  const auto m = flow.measurements();
  std::printf("\nend %.3f ms after %ld steps (%ld replans), wall %.0f s\n", m.time * 1e3, reacting.stats().steps,
              reacting.stats().replans, wall);
  std::printf("budgets: mass %.2e  energy %.2e\n", m.massBalanceError, m.energyBalanceError);
  std::printf("drift over last 1 ms: outlet mass flow %.2e, injector pressure %.2e, vacuum thrust %.2e\n",
              drift([](const Row& r) { return r.outlet; }), drift([](const Row& r) { return r.pInj; }),
              drift([](const Row& r) { return r.fVac; }));

  // Chamber-end stagnation pressure (mass-flux averaged) and the comparison.
  auto sol = Cantera::newSolution("h2o2.yaml", "", "none");
  auto& gas = *sol->thermo();
  double p0 = 0, flux = 0, pStatic = 0;
  for (int j = 0; j < nr; ++j) {
    auto q = mesh.index(endColumn, j);
    auto w = flow.cellPrimitive(q);
    auto y = flow.massFractions(q);
    double g = w.rho * w.uz * mesh.axialArea(endColumn + 1, j);
    p0 += g * stagnationPressure(gas, flow.temperature(q), w.p, y.data(), 0.5 * (w.uz * w.uz + w.ur * w.ur));
    pStatic += g * w.p;
    flux += g;
  }
  p0 /= flux;
  pStatic /= flux;
  const double fVac = m.deviceThrust + m.ambientAxialForce, mdot = m.outletMassFlow;
  const double cstarSim = p0 * at / mdot, ispSim = fVac / mdot / kG0;
  std::printf("chamber end: static %.5f MPa, equilibrium stagnation p0 %.5f MPa; injector face %.5f MPa\n",
              pStatic / 1e6, p0 / 1e6, m.injectorPressure / 1e6);
  std::printf("simulated: mdot out %.6f kg/s (supply %.6f), c* %.2f m/s, vacuum thrust %.2f N, vacuum Isp %.2f s\n",
              mdot, m.supplyMassFlow, cstarSim, fVac, ispSim);
  try {
    thermo::Mixture mix("h2o2.yaml");
    thermo::Reactant fuel{"H2:1", 300.0, std::nullopt}, ox{"O2:1", 300.0, std::nullopt};
    auto shifting = mix.idealRocket(fuel, ox, kMixtureRatio, p0, eps, thermo::Expansion::Equilibrium);
    auto frozen = mix.idealRocket(fuel, ox, kMixtureRatio, p0, eps, thermo::Expansion::FrozenAtChamber);
    const double cd = dischargeCoefficient(kArc / kRt, shifting.gammaChamber);
    const double lambda = (1 + std::cos(kTheta)) / 2;
    auto corrected = [&](const thermo::RocketPoint& r) {
      double pressurePart = r.pExit * eps * r.cstar / p0;
      return (lambda * (r.ispVac * kG0 - pressurePart) + pressurePart) / kG0;
    };
    std::printf("1-D ideal rocket at p0 (h2o2.yaml thermo): Tc %.1f K, c* %.2f m/s, Isp_vac shifting %.2f s, "
                "frozen %.2f s, gamma %.4f\n",
                shifting.tc, shifting.cstar, shifting.ispVac, frozen.ispVac, shifting.gammaChamber);
    std::printf("corrections: Cd %.5f (R = %.1f), lambda %.5f\n", cd, kArc / kRt, lambda);
    std::printf("predicted 2-D: c* %.2f m/s, Isp_vac shifting %.2f s, frozen %.2f s, F_vac %.2f N\n",
                shifting.cstar / cd, corrected(shifting), corrected(frozen), corrected(shifting) * kG0 * mdot);
    std::printf("simulated/predicted - 1: c* %+.4f%%, Isp_vac (shifting) %+.4f%%\n",
                100 * (cstarSim / (shifting.cstar / cd) - 1), 100 * (ispSim / corrected(shifting) - 1));
  } catch (const std::exception& e) {
    std::printf("ideal-rocket comparison failed at p0 = %.0f Pa: %s\n", p0, e.what());
  }
  return 0;
}
