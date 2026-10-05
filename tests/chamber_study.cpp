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
// Usage: crucible_chamber_study <eq|fr|frozen|frt|frp> <nz> <nr> <end time s> <threads> <output prefix>
//        [igniter energy J] [igniter duration s]
// Writes <prefix>_history.csv (every 2 us), <prefix>_mesh.csv (stations) and <prefix>_field_<us>.csv
// snapshots; tools/chamber_plots.py renders them.
//
// PaSR criterion 4 (stated in tests/pasr_tests.cpp): "frt" and "frp" are FiniteRate made viscous
// and turbulent as in tests/step_cost.cpp, "frp" with the closure (C_mix 1, S = {H2, O2, H2O}) and
// "frt" the control without it. Their default igniter is the one C1's FiniteRate run used, 300 J
// over 1 ms from 0.2 ms (CHAMBER_C1.md); the others default to 0.5 J over 0.2 ms. They judge (a)
// to (d) of criterion 4 and report its unjudged items. Light-off is the first history sample at
// which some cell has Y_H2O > 0.5 (a half-burnt cell; declared with the harness, before any run).
// The chamber volume for the s > 0.01 fraction is the cells upstream of the throat. The first-cell
// y+ is the laminar estimate sqrt(rho |u| d / mu) at the wall row, d the wall distance.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
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
    std::fprintf(stderr, "usage: %s <eq|fr|frozen|frt|frp> <nz> <nr> <end s> <threads> <prefix> [igniter J] [igniter s]\n",
                 argv[0]);
    return 2;
  }
  const std::string mode = argv[1], prefix = argv[6];
  const int nz = std::atoi(argv[2]), nr = std::atoi(argv[3]), threads = std::atoi(argv[5]);
  const double end = std::atof(argv[4]);
  const bool turbulent = mode == "frt" || mode == "frp", closure = mode == "frp";
  const double igniterEnergy = argc > 7 ? std::atof(argv[7]) : turbulent ? 300 : 0.5;
  const double igniterDuration = argc > 8 ? std::atof(argv[8]) : turbulent ? 1e-3 : 2e-4;
  const auto chemistry = mode == "eq" ? thermo::Chemistry::LocalEquilibrium
                         : mode == "fr" || turbulent ? thermo::Chemistry::FiniteRate
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
  // Igniter: from 0.2 ms in a 10 mm by 10 mm core 5 mm off the injector face; 0.5 J over 0.2 ms
  // unless given. Local equilibrium burns any premixed gas at once, so in "eq" mode it has no role.
  d.igniter = {0.005, 0.015, 0.010, igniterEnergy, 2e-4, igniterDuration};
  if (turbulent) {  // as tests/step_cost.cpp: Spalart-Rumsey ambient, supply I 0.05 and ratio 10
    d.transport = thermo::transportFits("h2o2.yaml");
    d.wallTemperature = 600;
    d.turbulence.enabled = true;
    d.supplies[0].turbulenceIntensity = 0.05;
    d.supplies[0].viscosityRatio = 10;
    const Medium medium = d.medium();
    std::vector<double> diffusion(d.species.size()), work;
    const double a0 = medium.soundSpeed(kAmbientT, d.composition.data());
    const double rho0 = kAmbientP / (medium.gasConstant(d.composition.data()) * kAmbientT);
    const double mu0 = medium.transport(kAmbientT, kAmbientP, d.composition.data(), diffusion.data(), work).viscosity;
    d.turbulence.ambientK = 9e-9 * a0 * a0;
    d.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
  }

  Flow flow(d);
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", threads, 1e-6, 1e-12, chemistry);
  if (closure) reacting.setMixingClosure(1.0, {"H2", "O2", "H2O"});
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
  std::printf("chemistry %s; igniter %.4g J over %.4g ms from 0.2 ms%s\n", mode.c_str(), igniterEnergy,
              igniterDuration * 1e3,
              turbulent ? (closure ? "; viscous, SST, PaSR closure on" : "; viscous, SST, closure off (control)") : "");

  {
    FILE* f = std::fopen((prefix + "_mesh.csv").c_str(), "w");
    std::fprintf(f, "i,z,radius\n");
    for (int i = 0; i <= nz; ++i) std::fprintf(f, "%d,%.9e,%.9e\n", i, i * mesh.dz, mesh.radius[i]);
    std::fclose(f);
  }
  const std::size_t ns = flow.medium().size();
  auto speciesIndex = [&](const char* name) {
    for (std::size_t k = 0; k < ns; ++k)
      if (flow.medium().species()[k].name == name) return k;
    return ns;
  };
  const std::size_t iH2O = speciesIndex("H2O"), iN2 = speciesIndex("N2"), iO2 = speciesIndex("O2");
  const std::vector<std::size_t> closureSpecies = {speciesIndex("H2"), iO2, iH2O};

  // Criterion 4 (d): kappa_eff in (0, 1] and s in [0, 1] at every field snapshot; the closure's
  // inputs are evaluated on the snapshot's state (Flow::mixingInputs) and kappa_eff from the
  // laminar rates there (reactingFraction). The control reports s too; its kappa_eff is 1.
  std::vector<double> mixTime, segregation, kappa;
  double worstKappa = 1, worstSegregation = 0, minKappa = 1, maxSegregation = 0;
  bool admissible = true;
  auto closureFields = [&]() {
    flow.mixingInputs(1.0, closureSpecies, mixTime, segregation);
    kappa.assign(flow.state().size(), 1.0);
    std::vector<double> z(ns + 1), dzdt(ns + 1);
    for (std::size_t q = 0; q < flow.state().size(); ++q) {
      if (closure) {
        auto y = flow.massFractions(q);
        z[0] = flow.temperature(q);
        std::copy(y.begin(), y.end(), z.begin() + 1);
        source.rates(flow.state()[q][0], z.data(), dzdt.data());
        kappa[q] = thermo::reactingFraction(z.data(), dzdt.data(), closureSpecies, {mixTime[q], segregation[q]});
      }
      if (!(kappa[q] > 0 && kappa[q] <= 1) || !(segregation[q] >= 0 && segregation[q] <= 1)) {
        if (admissible) { worstKappa = kappa[q]; worstSegregation = segregation[q]; }
        admissible = false;
      }
      minKappa = std::min(minKappa, kappa[q]);
      maxSegregation = std::max(maxSegregation, segregation[q]);
    }
  };
  // Fraction of the chamber volume (cells upstream of the throat) with s > 0.01.
  auto segregatedFraction = [&]() {
    double v = 0, total = 0;
    for (std::size_t q = 0; q < flow.state().size(); ++q)
      if (mesh.cells[q].z < geo.throat) {
        total += mesh.cells[q].volume;
        if (segregation[q] > 0.01) v += mesh.cells[q].volume;
      }
    return v / total;
  };

  auto writeField = [&](double t) {
    char name[512];
    std::snprintf(name, sizeof name, "%s_field_%06.0f.csv", prefix.c_str(), t * 1e6);
    FILE* f = std::fopen(name, "w");
    if (turbulent) closureFields();
    std::fprintf(f, "i,j,z,r,T,p,mach,uz,ur,Y_H2O,Y_O2,Y_N2%s\n", turbulent ? ",k,omega,mu_t,tau_mix,s,kappa_eff" : "");
    for (int i = 0; i < nz; ++i)
      for (int j = 0; j < nr; ++j) {
        auto q = mesh.index(i, j);
        auto w = flow.cellPrimitive(q);
        auto y = flow.massFractions(q);
        double temp = flow.temperature(q);
        double a = flow.medium().soundSpeed(temp, y.data());
        std::fprintf(f, "%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e", i, j, mesh.cells[q].z,
                     mesh.cells[q].r, temp, w.p, std::hypot(w.uz, w.ur) / a, w.uz, w.ur, y[iH2O], y[iO2], y[iN2]);
        if (turbulent)
          std::fprintf(f, ",%.6e,%.6e,%.6e,%.6e,%.6e,%.6e", flow.turbulence()[2 * q] / w.rho,
                       flow.turbulence()[2 * q + 1] / w.rho, flow.eddyViscosities()[q], mixTime[q], segregation[q],
                       kappa[q]);
        std::fprintf(f, "\n");
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
  // Criterion 4 (a) over the whole run and (b) at every sample.
  const double initialEnergy = flow.measurements().energy;
  double worstMass = 0, worstEnergyResidual = 0, negativeAt = -1, lightOff = -1, segregatedAtLightOff = 0;
  auto sample = [&]() {
    auto m = flow.measurements();
    if (turbulent) {
      worstMass = std::max(worstMass, std::abs(m.massBalanceError));
      worstEnergyResidual = std::max(worstEnergyResidual, std::abs(m.energyBalanceError * initialEnergy));
      double yMax = 0;
      for (std::size_t q = 0; q < flow.state().size(); ++q) {
        auto w = flow.cellPrimitive(q);
        if (negativeAt < 0 && !(w.rho > 0 && w.p > 0 && flow.temperature(q) > 0 && flow.turbulence()[2 * q] > 0 &&
                                flow.turbulence()[2 * q + 1] > 0))
          negativeAt = m.time;
        yMax = std::max(yMax, flow.massFractions(q)[iH2O]);
      }
      if (lightOff < 0 && yMax > 0.5) {
        lightOff = m.time;
        writeField(m.time);
        segregatedAtLightOff = segregatedFraction();
        std::printf("light-off at %.4f ms (largest Y_H2O %.3f); chamber volume with s > 0.01: %.4f\n", m.time * 1e3,
                    yMax, segregatedAtLightOff);
      }
    }
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
  auto march = [&](double end) {
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
  };
  // Settling over the last millisecond (declared criterion: relative drift below 1e-3).
  auto drift = [&](auto get) {
    double lo = 1e300, hi = -1e300, last = get(rows.back());
    for (const auto& r : rows)
      if (r.t >= rows.back().t - 1e-3) { lo = std::min(lo, get(r)); hi = std::max(hi, get(r)); }
    return (hi - lo) / std::abs(last);
  };
  auto settled = [&]() {
    const auto m = flow.measurements();
    return std::max({drift([](const Row& r) { return r.outlet; }), drift([](const Row& r) { return r.pInj; }),
                     drift([](const Row& r) { return r.fVac; })}) < 1e-3 &&
           std::abs(m.outletMassFlow / m.supplyMassFlow - 1) < 1e-3;
  };
  try {
    march(end);
    writeField(flow.time());
    // Criterion 4 (c): a run that fails only the settling check is extended once to 12 ms.
    const double energyEnd = std::abs(flow.measurements().energy);
    if (turbulent && end >= 8e-3 && end < 12e-3 && !settled() && worstMass < 1e-11 && worstEnergyResidual / energyEnd < 1e-11 &&
        negativeAt < 0 && admissible) {
      std::printf("criterion 4 (c) fails alone at %.3f ms: extending once to 12 ms\n", flow.time() * 1e3);
      march(12e-3);
      writeField(flow.time());
    }
  } catch (const std::exception& e) {
    std::fclose(history);
    std::printf("\nexception at %.6f ms after %ld steps: %s\n", flow.time() * 1e3, reacting.stats().steps, e.what());
    writeField(flow.time());  // the last state the run reached, for the diagnosis
    if (turbulent) std::printf("criterion 4 (b): the run did not reach its end  FAIL\n");
    return 1;
  }
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
  if (!turbulent) return 0;

  // PaSR criterion 4: judged (a) to (d), then the reported items.
  int failed = 0;
  auto judge = [&](bool ok, const char* what, double value, double limit) {
    std::printf("  %-62s %.3e (limit %.1e)  %s\n", what, value, limit, ok ? "ok" : "FAIL");
    if (!ok) ++failed;
  };
  const double energyEnd = std::abs(m.energy);
  std::printf("\nPaSR criterion 4, %s, %dx%d, to %.3f ms\n", closure ? "closure on" : "control (closure off)", nz, nr,
              m.time * 1e3);
  judge(worstMass < 1e-11, "(a) mass budget, largest over the run (initial fill)", worstMass, 1e-11);
  judge(worstEnergyResidual / energyEnd < 1e-11, "(a) energy budget, largest over the run (|E| at the end)",
        worstEnergyResidual / energyEnd, 1e-11);
  std::printf("  reported: energy budget against the initial N2 fill, largest %.3e, at the end %.3e\n",
              worstEnergyResidual / std::abs(initialEnergy), m.energyBalanceError);
  std::printf("  (b) positivity of rho, p, T, k, omega at every sample (%zu samples)  %s\n", rows.size(),
              negativeAt < 0 ? "ok" : "FAIL");
  if (negativeAt >= 0) { std::printf("      first non-positive value at %.6f ms\n", negativeAt * 1e3); ++failed; }
  const double settle = std::max({drift([](const Row& r) { return r.outlet; }), drift([](const Row& r) { return r.pInj; }),
                                  drift([](const Row& r) { return r.fVac; })});
  judge(settle < 1e-3, "(c) largest drift over the last 1 ms", settle, 1e-3);
  judge(std::abs(mdot / m.supplyMassFlow - 1) < 1e-3, "(c) outlet mass flow against the supply",
        std::abs(mdot / m.supplyMassFlow - 1), 1e-3);
  std::printf("  (d) kappa_eff in (0, 1] and s in [0, 1] at every field snapshot  %s (smallest kappa_eff %.4e, "
              "largest s %.4e)\n", admissible ? "ok" : "FAIL", minKappa, maxSegregation);
  if (!admissible) { std::printf("      first bad cell: kappa_eff %.6e, s %.6e\n", worstKappa, worstSegregation); ++failed; }
  std::printf("reported: light-off at %s ms; chamber volume with s > 0.01: %.4f at light-off, %.4f at the end\n",
              lightOff < 0 ? "never" : std::to_string(lightOff * 1e3).c_str(), segregatedAtLightOff, segregatedFraction());
  std::printf("reported: largest mass fraction clipped to zero after a reaction substep %.3e\n",
              reacting.stats().maxClippedFraction);
  {
    std::vector<double> yPlus, diffusion(ns), work;
    double chamberMax = 0;
    for (int i = 0; i < nz; ++i) {
      auto q = mesh.index(i, nr - 1);
      auto w = flow.cellPrimitive(q);
      auto y = flow.massFractions(q);
      const double mu = flow.medium().transport(flow.temperature(q), w.p, y.data(), diffusion.data(), work).viscosity;
      yPlus.push_back(std::sqrt(w.rho * std::hypot(w.uz, w.ur) * flow.wallDistances()[q] / mu));
      if (mesh.cells[q].z < geo.throat) chamberMax = std::max(chamberMax, yPlus.back());
    }
    std::vector<double> sorted = yPlus;
    std::sort(sorted.begin(), sorted.end());
    std::printf("reported: first-cell y+ at the wall row (laminar estimate): median %.1f, largest %.1f, largest "
                "upstream of the throat %.1f\n", sorted[sorted.size() / 2], sorted.back(), chamberMax);
  }
  std::printf("reported against the control: light-off, injector pressure %.6f MPa, c* %.2f m/s, vacuum Isp %.2f s\n",
              m.injectorPressure / 1e6, cstarSim, ispSim);
  std::printf("%d failures\n", failed);
  return failed ? 1 : 0;
}
