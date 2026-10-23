// The face inlet (a) and the plate's wall temperature: gates F2, F4 and F5 of docs/evidence/FACE_INLET.md
// (registered before the code). F1 is the fast ctest and rl10dev's short run; F3 is tests/rl10_tier1.cpp
// with CRUCIBLE_T1_OPEN_AREA set.
//
// f2: switch on with every supply on ring nodes (each phi 0 or 1) against the ring-centre code, laminar
//   and SST with the wall functions, 200 steps: the state's largest relative difference (per field, over
//   the field's largest magnitude) under 1e-10, the inlet mass flow and momentum flux and the wall's axial
//   force equal to 1e-12 relative. Reported: whether the runs are bit-identical.
// f4: AMR (16 x 6 with 2 finer levels, all cells active), supplies cutting face rings (a closed port among
//   them): over each coarse face ring, its children's summed area x flux against its own, for mass, axial
//   momentum and energy at the uniform initial state (1e-12 relative), and for mass and energy after 100
//   base steps (momentum reported).
// f5: (a) plateTemperature empty and equal to wallTemperature, wall functions and resolved SST, 200 steps:
//   state, partial densities, rho k and rho omega and the budgets bit-identical. (b) plateTemperature 0 with
//   the wall functions, ring-centre and open-area supplies (a supply ring with a wall share): after 200
//   steps and a transport evaluation at the final state, every plate face's q_w is zero to rounding
//   (|q_w| below 1e-9 of the side wall's largest) and its T_w lies in [T_1, T_1 + u_1^2 / (2 c_p(T_1))] to 1e-9 K.
//   Reported: (T_w - T_1) / (u_1^2 / (2 c_p)), near Pr^(1/3), and the open-area run's budgets.
// f6 (FACE_INLET.md, *The gas energy leak*): F4's layout on one level (its closed port a supply ring with a
//   wall share), each model set and plate temperature, and a ring-centre run on ring nodes: the mass, energy
//   and axial momentum budgets after one step and after 200 steps, each below 1e-13.
#include "adapters/reaction.hpp"
#include "core/amr.hpp"
#include "core/flow.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

using namespace crucible;

namespace {
int failures = 0;
void check(bool ok, const std::string& what, double value, double limit) {
  std::printf("  %-78s %.3e (limit %.0e)  %s\n", what.c_str(), value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}
double rel(double a, double b) {
  const double s = std::max(std::abs(a), std::abs(b));
  return s > 0 ? std::abs(a - b) / s : 0.0;
}
struct Mixture {
  std::vector<Species> species;
  TransportFits fits;
  std::vector<double> nitrogen, mixture;
};
Mixture gas() {
  const std::string mechanism = "h2o2.yaml";
  thermo::ReactionSource source(mechanism);
  return {source.medium().species(), thermo::transportFits(mechanism), source.massFractions("N2:1"),
          source.massFractions("H2:1, N2:4")};
}
// wall_function_study's chamber: contoured side wall at 400 K, N2 at 1 atm and 300 K.
Definition chamber(const Mixture& g, bool turbulent, bool law) {
  Definition d;
  d.experiment = Case::Chamber;
  d.nz = 16; d.nr = 6;
  d.contour = {{0, 5e-3}, {0.02, 5e-3}, {0.03, 3.5e-3}};
  d.species = g.species; d.composition = g.nitrogen;
  d.transport = g.fits;
  d.backPressure = 101325; d.ambientTemperature = 300;
  d.wallTemperature = 400;
  if (turbulent) {
    d.turbulence.enabled = true; d.turbulence.ambientK = 1; d.turbulence.ambientOmega = 1e3;
    d.turbulence.wallFunctions = law;
  }
  return d;
}
Supply supply(double inner, double outer, double massFlow, const std::vector<double>& y, double ramp = 2e-5) {
  Supply s;
  s.innerRadius = inner; s.outerRadius = outer; s.massFlow = massFlow; s.totalTemperature = 300;
  s.composition = y; s.opens = 0; s.ramp = ramp;
  s.turbulenceIntensity = 0.05; s.viscosityRatio = 10;
  return s;
}
Flow run(Definition d, int steps) {
  Flow flow(std::move(d));
  for (int n = 0; n < steps; ++n) flow.step();
  return flow;
}
// The largest over a field-major array's fields of max |a - b| / max |a| (stride fields per cell).
double fieldDifference(const double* a, const double* b, std::size_t count, std::size_t stride) {
  double worst = 0;
  for (std::size_t f = 0; f < stride; ++f) {
    double scale = 0, diff = 0;
    for (std::size_t q = f; q < count; q += stride) {
      scale = std::max(scale, std::abs(a[q]));
      diff = std::max(diff, std::abs(a[q] - b[q]));
    }
    if (scale > 0) worst = std::max(worst, diff / scale);
  }
  return worst;
}
bool identical(const Flow& a, const Flow& b) {
  const auto ma = a.measurements(), mb = b.measurements();
  const std::array<double, 7> la{ma.massBalanceError, ma.energyBalanceError, ma.momentumBalanceError, ma.wallHeatFlow,
                                 ma.wallAxialForce, ma.inletMomentumFlux, ma.prescribedTurbulentEnergy};
  const std::array<double, 7> lb{mb.massBalanceError, mb.energyBalanceError, mb.momentumBalanceError, mb.wallHeatFlow,
                                 mb.wallAxialForce, mb.inletMomentumFlux, mb.prescribedTurbulentEnergy};
  return a.state().size() == b.state().size() &&
         std::memcmp(a.state().data(), b.state().data(), a.state().size() * sizeof(Conserved)) == 0 &&
         a.partialDensities() == b.partialDensities() && a.turbulence().size() == b.turbulence().size() &&
         std::memcmp(a.turbulence().data(), b.turbulence().data(), a.turbulence().size() * sizeof(double)) == 0 &&
         std::memcmp(la.data(), lb.data(), sizeof la) == 0 && ma.steps == mb.steps;
}

void f2(const Mixture& g) {
  std::printf("F2: open-area supplies on ring nodes against the ring-centre code, 200 steps\n");
  for (const auto& [turbulent, label] : {std::pair{false, "laminar"}, std::pair{true, "SST, wall functions"}}) {
    Definition d = chamber(g, turbulent, turbulent);
    const Mesh mesh(d);
    auto node = [&](int j) { return mesh.radius[0] * mesh.fraction[static_cast<std::size_t>(j)]; };
    d.supplies = {supply(0, node(2), 5e-4, g.mixture), supply(node(3), node(4), 1e-3, g.nitrogen)};
    Definition o = d;
    o.openAreaSupplies = true;
    const Flow a = run(d, 200), b = run(o, 200);
    const auto ma = a.measurements(), mb = b.measurements();
    const double state = fieldDifference(reinterpret_cast<const double*>(a.state().data()),
                                         reinterpret_cast<const double*>(b.state().data()), a.state().size() * 4, 4);
    const double species = fieldDifference(a.partialDensities().data(), b.partialDensities().data(),
                                           a.partialDensities().size(), g.species.size());
    const double turbulence = turbulent ? fieldDifference(a.turbulence().data(), b.turbulence().data(), a.turbulence().size(), 2) : 0.0;
    std::printf("  %s: t %.6e s (%llu steps); bit-identical %s; budgets (open) mass %.3e energy %.3e\n", label, mb.time,
                static_cast<unsigned long long>(mb.steps), identical(a, b) ? "yes" : "no", mb.massBalanceError, mb.energyBalanceError);
    const std::string p = std::string(label) + ": ";
    check(std::max({state, species, turbulence}) < 1e-10, p + "state, partial densities, rho k, rho omega: largest relative difference",
          std::max({state, species, turbulence}), 1e-10);
    check(rel(ma.inletMassFlow, mb.inletMassFlow) < 1e-12, p + "inlet mass flow", rel(ma.inletMassFlow, mb.inletMassFlow), 1e-12);
    check(rel(ma.inletMomentumFlux, mb.inletMomentumFlux) < 1e-12, p + "inlet momentum flux",
          rel(ma.inletMomentumFlux, mb.inletMomentumFlux), 1e-12);
    check(rel(ma.wallAxialForce, mb.wallAxialForce) < 1e-12, p + "wall axial force", rel(ma.wallAxialForce, mb.wallAxialForce), 1e-12);
  }
}

void f4(const Mixture& g) {
  std::printf("F4: AMR, each coarse face ring's flux against its children's\n");
  Definition d = chamber(g, false, false);
  d.openAreaSupplies = true;
  d.supplies = {supply(0, 1.3e-3, 3e-4, g.mixture, 0), supply(2.1e-3, 2.9e-3, 0, g.nitrogen, 0),
                supply(3.3e-3, 4.6e-3, 2e-3, g.nitrogen, 0)};
  Amr amr(d, 2);
  for (int l = 1; l < amr.levels(); ++l) amr.setActive(l, std::vector<unsigned char>(amr.level(l).mesh().cells.size(), 1));
  auto compare = [&](const char* when, bool momentum) {
    std::array<double, 3> worst{};
    for (int l = 0; l + 1 < amr.levels(); ++l) {
      Flow& coarse = amr.level(l);
      Flow& fine = amr.level(l + 1);
      const auto c = coarse.injectorFluxes(), f = fine.injectorFluxes();
      const Mesh& cm = coarse.mesh();
      const Mesh& fm = fine.mesh();
      for (int j = 0; j < cm.nr; ++j) {
        int n = 0;
        for (int k : {0, 1, 3}) {
          const double own = cm.axialArea(0, j) * c[static_cast<std::size_t>(j)][k];
          const double children = fm.axialArea(0, 2 * j) * f[static_cast<std::size_t>(2 * j)][k] +
                                  fm.axialArea(0, 2 * j + 1) * f[static_cast<std::size_t>(2 * j + 1)][k];
          worst[static_cast<std::size_t>(n)] = std::max(worst[static_cast<std::size_t>(n)], rel(own, children));
          ++n;
        }
      }
    }
    std::printf("  %s (t %.6e s): largest relative difference mass %.3e, axial momentum %.3e, energy %.3e\n", when, amr.time(),
                worst[0], worst[1], worst[2]);
    check(worst[0] < 1e-12, std::string(when) + ": mass", worst[0], 1e-12);
    if (momentum) check(worst[1] < 1e-12, std::string(when) + ": axial momentum", worst[1], 1e-12);
    check(worst[2] < 1e-12, std::string(when) + ": energy", worst[2], 1e-12);
  };
  compare("uniform initial state", true);
  for (int n = 0; n < 100; ++n) amr.step();
  compare("after 100 base steps", false);
}

void f5(const Mixture& g) {
  std::printf("F5(a): plateTemperature empty against equal to wallTemperature, 200 steps\n");
  for (const auto& [law, label] : {std::pair{true, "wall functions"}, std::pair{false, "resolved SST"}}) {
    Definition d = chamber(g, true, law);
    d.supplies = {supply(0, 3e-3, 2e-3, g.nitrogen)};
    Definition e = d;
    e.plateTemperature = d.wallTemperature;
    const bool same = identical(run(d, 200), run(e, 200));
    check(same, std::string(label) + ": state, species, rho k, rho omega, budgets differ (1 = yes)", same ? 0 : 1, 0);
  }
  std::printf("F5(b): an adiabatic plate under the wall functions, 200 steps\n");
  for (bool open : {false, true}) {
    Definition d = chamber(g, true, true);
    d.plateTemperature = 0;
    d.openAreaSupplies = open;
    d.supplies = {supply(0, 3e-3, 2e-3, g.nitrogen)};
    Flow flow = run(d, 200);
    std::vector<Conserved> derivative;
    std::vector<double> species, turbulence;
    flow.transportDerivative(derivative, species, turbulence);  // the wall solutions at the final state
    double side = 0, plate = 0, below = 0, above = 0, lo = 1e300, hi = -1e300;
    int faces = 0;
    std::vector<double> diffusion(g.species.size()), work;
    for (const auto& w : flow.wallYPlus()) {
      if (!w.plate) { side = std::max(side, std::abs(w.heat)); continue; }
      ++faces;
      plate = std::max(plate, std::abs(w.heat));
      const double t1 = flow.temperature(w.cell), u1 = std::abs(flow.cellPrimitive(w.cell).ur);
      const auto y = flow.massFractions(w.cell);
      const auto pr = flow.medium().properties(t1, y.data());
      const double rise = u1 * u1 / (2 * (pr.cv + pr.r));
      below = std::max(below, t1 - w.temperature);
      above = std::max(above, w.temperature - (t1 + rise));
      if (rise > 1e-9 * t1) { lo = std::min(lo, (w.temperature - t1) / rise); hi = std::max(hi, (w.temperature - t1) / rise); }
    }
    const auto m = flow.measurements();
    const std::string p = open ? "open-area: " : "ring-centre: ";
    std::printf("  %s%d plate faces; largest |q_w| plate %.3e W/m^2, side wall %.3e W/m^2; (T_w - T_1) / (u_1^2 / 2 c_p) in [%.6f, %.6f]; "
                "budgets mass %.3e energy %.3e\n", p.c_str(), faces, plate, side, lo, hi, m.massBalanceError, m.energyBalanceError);
    check(faces > 0 && plate <= 1e-9 * side, p + "plate |q_w| over the side wall's largest", side > 0 ? plate / side : 1.0, 1e-9);
    check(below <= 1e-9, p + "T_1 - T_w at its largest (K)", below, 1e-9);
    check(above <= 1e-9, p + "T_w - (T_1 + u_1^2 / (2 c_p)) at its largest (K)", above, 1e-9);
  }
}

void f6(const Mixture& g) {
  std::printf("F6: the budgets with a supply ring that has a wall share, one step and 200 steps\n");
  struct Variant {
    const char* label;
    bool open, turbulent, law, corrections, viscous;
    double plate;
  };
  const Variant variants[] = {
      {"open-area, SST, wall functions and corrections, plate 400 K", true, true, true, true, true, 400},
      {"open-area, SST, wall functions, plate 400 K", true, true, true, false, true, 400},
      {"open-area, SST, wall functions, adiabatic plate", true, true, true, false, true, 0},
      {"open-area, resolved SST, plate 400 K", true, true, false, false, true, 400},
      {"open-area, laminar, plate 400 K", true, false, false, false, true, 400},
      {"open-area, inviscid", true, false, false, false, false, 0},
      {"ring-centre on ring nodes, SST, wall functions and corrections, plate 400 K", false, true, true, true, true, 400}};
  for (const auto& v : variants) {
    Definition d = chamber(g, v.turbulent, v.law);
    d.turbulence.propertyCorrections = v.corrections;
    if (v.viscous) d.plateTemperature = v.plate;
    else d.transport = {};
    const Mesh mesh(d);
    if (v.open) {
      d.openAreaSupplies = true;
      d.supplies = {supply(0, 1.3e-3, 3e-4, g.mixture, 0), supply(2.1e-3, 2.9e-3, 0, g.nitrogen, 0),
                    supply(3.3e-3, 4.6e-3, 2e-3, g.nitrogen, 0)};
    } else {
      auto node = [&](int j) { return mesh.radius[0] * mesh.fraction[static_cast<std::size_t>(j)]; };
      d.supplies = {supply(0, node(2), 5e-4, g.mixture, 0), supply(node(3), node(4), 1e-3, g.nitrogen, 0)};
    }
    // Each face ring's wall share w_j as Flow::openAreaShares sets it; the supply rings with one (0 < w_j < 1/2).
    std::string shared;
    if (v.open)
      for (int j = 0; j < mesh.nr; ++j) {
        const double lo = mesh.radius[0] * mesh.fraction[static_cast<std::size_t>(j)];
        const double hi = mesh.radius[0] * mesh.fraction[static_cast<std::size_t>(j) + 1];
        double covered = 0;
        for (const auto& s : d.supplies) {
          const double a = std::max(0.0, std::min(s.outerRadius, hi) * std::min(s.outerRadius, hi) -
                                             std::max(s.innerRadius, lo) * std::max(s.innerRadius, lo));
          if (s.outerRadius > lo && s.innerRadius < hi) covered += a / (hi * hi - lo * lo);
        }
        const double w = 1 - covered;
        if (w > 1e-12 && w < 0.5) shared += " j " + std::to_string(j) + " (w " + std::to_string(w).substr(0, 5) + ")";
      }
    Flow flow(d);
    flow.step();
    const auto one = flow.measurements();
    for (int n = 1; n < 200; ++n) flow.step();
    const auto end = flow.measurements();
    std::printf("  %s; supply rings with a wall share:%s\n    one step: mass %.3e energy %.3e momentum %.3e, plate heat "
                "%.4e W (booked, into the gas); 200 steps (t %.4e s): mass %.3e energy %.3e momentum %.3e\n",
                v.label, shared.empty() ? " none" : shared.c_str(), one.massBalanceError, one.energyBalanceError,
                one.momentumBalanceError, one.plateHeatFlow, end.time, end.massBalanceError, end.energyBalanceError,
                end.momentumBalanceError);
    const double worst = std::max({std::abs(one.massBalanceError), std::abs(one.energyBalanceError), std::abs(one.momentumBalanceError),
                                   std::abs(end.massBalanceError), std::abs(end.energyBalanceError), std::abs(end.momentumBalanceError)});
    check(worst < 1e-13, std::string(v.label) + ": largest |budget|", worst, 1e-13);
  }
}
}  // namespace

int main(int argc, char** argv) {
  try {
    const std::string mode = argc > 1 ? argv[1] : "all";
    const Mixture g = gas();
    if (mode == "all" || mode == "f2") f2(g);
    if (mode == "all" || mode == "f4") f4(g);
    if (mode == "all" || mode == "f5") f5(g);
    if (mode == "all" || mode == "f6") f6(g);
    std::printf("%d failures\n", failures);
    return failures ? 1 : 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
