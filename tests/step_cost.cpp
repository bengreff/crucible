// Cost of the reacting multi-species step per cell, for the compute budget in TECHNICAL_PLAN
// (section "Compute budget"). Not a verification test: it measures cell updates per second on the
// machine that runs it.
//
// The C1 chamber of tests/chamber_study.cpp (same contour and supply: premixed H2/O2 at O/F 5,
// 0.4 kg/s, h2o2.yaml, ambient N2 at 1000 Pa) is first marched in LocalEquilibrium mode to the setup
// time, so the chamber holds burned gas and the nozzle flows. From that one state each
// configuration takes a few untimed steps and then the timed steps:
//   inviscid: slip walls, no transport (the C1 configuration);
//   viscous:  mixture-averaged transport, SST-2003, no-slip walls at 600 K, supply turbulence
//             I 0.05 and viscosity ratio 10, Spalart-Rumsey ambient; k = 3/2 (0.05 |u|)^2 and
//             omega = rho k / (10 mu) in every cell at the start;
// each with Frozen chemistry (the flow step alone) and FiniteRate (CVODES, rtol 1e-6, atol 1e-12,
// as C1), on 1 thread and on the given thread count. The flow step is serial; only the chemistry
// is threaded. The valve is fully open from the start of each timed configuration.
// Last, the FiniteRate reaction substep alone on 1 thread (inviscid, after 20 FiniteRate steps from
// the setup state, so the cells are off equilibrium as in a march), 20 calls at each of three
// substep lengths, because the CVODES cost of a substep may depend on its length: half the
// measured flow step, and the halves of the explicit steps estimated for a wall-resolved RL10
// grid (TECHNICAL_PLAN, "Compute budget").
// With the optional last argument "equilibrium" (added 5 October 2026), the configurations are
// Frozen and LocalEquilibrium (a constant-(u, v) equilibrium per cell after each flow step, no
// CVODES) instead of Frozen and FiniteRate, and the substep sweep is skipped.
//
// Usage: crucible_step_cost <nz> <nr> <setup time s> <timed steps> <threads> [equilibrium]
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"

using namespace crucible;

// Chamber geometry and supply: copied from tests/chamber_study.cpp (keep them equal).
namespace {
constexpr double kPi = std::numbers::pi;
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
}  // namespace

int main(int argc, char** argv) {
  if (argc < 6) {
    std::fprintf(stderr, "usage: %s <nz> <nr> <setup time s> <timed steps> <threads> [equilibrium]\n", argv[0]);
    return 2;
  }
  const int nz = std::atoi(argv[1]), nr = std::atoi(argv[2]), steps = std::atoi(argv[4]), threads = std::atoi(argv[5]);
  const double setup = std::atof(argv[3]);
  const bool equilibrium = argc > 6 && std::string(argv[6]) == "equilibrium";
  const auto reactingMode = equilibrium ? thermo::Chemistry::LocalEquilibrium : thermo::Chemistry::FiniteRate;
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

  using Clock = std::chrono::steady_clock;
  auto since = [](Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); };
  Flow start(d);
  {
    thermo::ReactingFlow eq(start, "h2o2.yaml", threads, 1e-6, 1e-12, thermo::Chemistry::LocalEquilibrium);
    const auto t0 = Clock::now();
    eq.advanceTo(setup);
    const auto m = start.measurements();
    std::printf("setup %dx%d (%zu cells): LocalEquilibrium to %.3f ms, %lld steps, %.0f s wall; injector p %.4f MPa, outlet %.4f kg/s\n",
                nz, nr, start.state().size(), setup * 1e3, static_cast<long long>(m.steps), since(t0),
                m.injectorPressure / 1e6, m.outletMassFlow);
  }
  const std::size_t cells = start.state().size(), ns = start.medium().size();
  std::vector<Primitive> primitives(cells);
  std::vector<double> fractions(cells * ns), temperatures(cells);
  for (std::size_t q = 0; q < cells; ++q) {
    primitives[q] = start.cellPrimitive(q);
    temperatures[q] = start.temperature(q);
    auto y = start.massFractions(q);
    std::copy(y.begin(), y.end(), fractions.begin() + static_cast<std::ptrdiff_t>(q * ns));
  }

  const auto fits = thermo::transportFits("h2o2.yaml");
  std::printf("%-9s %-10s %7s %12s %12s %14s %16s %12s\n", "config", "chemistry", "threads", "dt start [s]", "mean dt [s]",
              "wall/step [s]", "cell updates/s", "replans/step");
  for (bool viscous : {false, true})
    for (auto chemistry : {thermo::Chemistry::Frozen, reactingMode})
      for (int t : {1, threads}) {
        if (t != 1 && (chemistry == thermo::Chemistry::Frozen || threads == 1)) continue;
        Definition c = d;
        c.supplies[0].ramp = 1e-9;
        if (viscous) {
          c.transport = fits;
          c.wallTemperature = 600;
          c.turbulence.enabled = true;
          c.supplies[0].turbulenceIntensity = 0.05;
          c.supplies[0].viscosityRatio = 10;
          const Medium medium = c.medium();
          const auto y = c.composition;
          std::vector<double> diffusion(ns), work;
          const double a0 = medium.soundSpeed(kAmbientT, y.data()), rho0 = kAmbientP / (medium.gasConstant(y.data()) * kAmbientT);
          const double mu0 = medium.transport(kAmbientT, kAmbientP, y.data(), diffusion.data(), work).viscosity;
          c.turbulence.ambientK = 9e-9 * a0 * a0;
          c.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
        }
        Flow flow(c);
        flow.setInitialState(primitives, fractions);
        if (viscous) {
          const Medium medium = c.medium();
          std::vector<double> kOmega(2 * cells), diffusion(ns), work;
          for (std::size_t q = 0; q < cells; ++q) {
            const auto& w = primitives[q];
            const double mu = medium.transport(temperatures[q], w.p, fractions.data() + q * ns, diffusion.data(), work).viscosity;
            const double k = 1.5 * std::pow(0.05 * std::hypot(w.uz, w.ur), 2) + c.turbulence.ambientK;
            kOmega[2 * q] = k;
            kOmega[2 * q + 1] = w.rho * k / (10 * mu);
          }
          flow.setTurbulence(kOmega);
        }
        thermo::ReactingFlow reacting(flow, "h2o2.yaml", t, 1e-6, 1e-12, chemistry);
        const double dt0 = flow.stableDt();
        for (int i = 0; i < 5; ++i) reacting.step();
        const double sim0 = flow.time();
        const long replans0 = reacting.stats().replans;
        const auto t0 = Clock::now();
        for (int i = 0; i < steps; ++i) reacting.step();
        const double wall = since(t0);
        std::printf("%-9s %-10s %7d %12.4e %12.4e %14.4e %16.4e %12.2f\n", viscous ? "viscous" : "inviscid",
                    chemistry == thermo::Chemistry::Frozen ? "frozen" : equilibrium ? "equil" : "finite", t, dt0, (flow.time() - sim0) / steps,
                    wall / steps, double(cells) * steps / wall, double(reacting.stats().replans - replans0) / steps);
        std::fflush(stdout);
      }
  if (equilibrium) return 0;
  std::printf("%-26s %16s\n", "FiniteRate substep [s]", "s per cell");
  Definition open = d;
  open.supplies[0].ramp = 1e-9;
  Flow flow(open);
  flow.setInitialState(primitives, fractions);
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", 1, 1e-6, 1e-12, thermo::Chemistry::FiniteRate);
  for (int i = 0; i < 20; ++i) reacting.step();
  const double half = 0.5 * flow.stableDt();
  for (double h : {half, 5e-10, 5e-11}) {
    reacting.react(h);
    const auto t0 = Clock::now();
    for (int i = 0; i < 20; ++i) reacting.react(h);
    std::printf("%-26.4e %16.4e\n", h, since(t0) / (20.0 * double(cells)));
  }
  return 0;
}
