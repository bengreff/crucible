// The spray on AMR's levels: criteria A1 to A3 of docs/evidence/SPRAY.md, Spray on AMR (stated before the
// code; A0, the uniform path unchanged, is docs/evidence/spray_amr/a0_dump.cpp against the build before).
// usage: crucible_spray_amr_tests [a1 | a2]   (no argument: all; a2 includes A3)
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/amr.hpp"
#include "core/spray.hpp"

using namespace crucible;

namespace {
const std::string kSaturation = CRUCIBLE_SOURCE_DIR "/tools/o2_saturation_nist.tsv";
int failures = 0;
void check(bool ok, const std::string& what, double value, double limit) {
  std::printf("  %-86s %.3e (limit %.0e)  %s\n", what.c_str(), value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}
double seconds(std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

constexpr double kRadius = 0.01, kLength = 0.016;

// S3's chamber (tests/spray_tests.cpp, chamberRun), and S8(c)'s with SST.
Definition chamber(thermo::ReactionSource& source, bool turbulent) {
  Definition d;
  d.experiment = Case::Chamber;
  d.nz = 16;
  d.nr = 6;
  d.contour = {{0, kRadius}, {kLength, kRadius}};
  d.species = source.medium().species();
  d.composition = source.massFractions("H2:1");
  d.backPressure = 1e6;
  d.ambientTemperature = 300;
  d.transport = thermo::transportFits("h2o2.yaml", 60);
  Supply supply;
  supply.innerRadius = kRadius / 6;
  supply.outerRadius = kRadius / 2;
  supply.massFlow = 3e-3;
  supply.totalTemperature = 300;
  supply.composition = source.massFractions("H2:1");
  if (turbulent) {
    supply.turbulenceIntensity = 0.05;
    supply.viscosityRatio = 10;
    d.turbulence.enabled = true;
    const Medium medium = d.medium();
    std::vector<double> diffusion(d.species.size()), work;
    const double a0 = medium.soundSpeed(300, d.composition.data());
    const double rho0 = 1e6 / (medium.gasConstant(d.composition.data()) * 300);
    const double mu0 = medium.transport(300, 1e6, d.composition.data(), diffusion.data(), work).viscosity;
    d.turbulence.ambientK = 9e-9 * a0 * a0;
    d.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
  }
  d.supplies = {supply};
  return d;
}
Spray::Settings settings() {
  Spray::Settings s;
  Spray::Injector j;
  j.rOuter = 1.5e-3;
  j.massFlow = 2e-3;
  j.temperature = 90;
  j.velocity = 50;
  j.coneHalfAngle = 35 * std::numbers::pi / 180;
  j.size = Spray::Injector::Size::RosinRammler;
  j.sauterDiameter = 30e-6;
  j.exponent = 3;
  j.parcelRate = 2e5;
  s.injectors = {j};
  s.restitution = 0.5;
  return s;
}
LiquidOxygen oxygen(const Flow& flow) {
  std::vector<double> y(flow.medium().size(), 0.0);
  for (std::size_t k = 0; k < y.size(); ++k)
    if (flow.medium().species()[k].name == "O2") y[k] = 1;
  return LiquidOxygen::load(kSaturation, flow.medium().enthalpy(298.15, y.data()));
}

// Every field of a parcel but the AMR phase.
std::vector<unsigned char> bytes(const std::vector<Parcel>& parcels) {
  std::vector<unsigned char> out;
  auto put = [&](const void* p, std::size_t n) {
    const auto* c = static_cast<const unsigned char*>(p);
    out.insert(out.end(), c, c + n);
  };
  for (const auto& p : parcels) {
    const double x[16]{p.z, p.r, p.uz, p.ur, p.ut, p.d, p.t, p.count, p.shed, p.rtTime, p.distance,
                       p.injectedMass, p.eddyZ, p.eddyR, p.eddyT, p.eddyTime};
    put(x, sizeof x);
    put(&p.random, sizeof p.random);
    put(&p.cell, sizeof p.cell);
    put(&p.injector, sizeof p.injector);
  }
  return out;
}
template <class T>
bool same(const std::vector<T>& a, const std::vector<T>& b) {
  return a.size() == b.size() && (a.empty() || std::memcmp(a.data(), b.data(), a.size() * sizeof(T)) == 0);
}
// The largest difference of the conserved state over each component's largest magnitude.
double stateDifference(const Flow& a, const Flow& b) {
  double out = 0;
  for (int c = 0; c < 4; ++c) {
    double scale = 0, diff = 0;
    for (std::size_t q = 0; q < a.state().size(); ++q) {
      scale = std::max(scale, std::abs(b.state()[q][c]));
      diff = std::max(diff, std::abs(a.state()[q][c] - b.state()[q][c]));
    }
    out = std::max(out, diff / scale);
  }
  return out;
}

// A1: single-level Amr (level 1 inactive) against a Flow on the AMR base's mesh, bit for bit. The reference
// drives the Flow as the precondition did (Flow::step, then Spray::advance): ReactingFlow's frozen step, the
// S3 run's own driver, plans with an extra Flow::stableDt, whose refresh moves the temperature guesses and so
// the gas by rounding (8e-15 relative at the second step with no spray, measured), so it is printed beside.
void a1(thermo::ReactionSource& source, bool turbulent, double end) {
  const char* name = turbulent ? "A1(b) S8(c), SST and dispersion, 0.2 ms" : "A1(a) S3, frozen, 0.5 ms";
  std::printf("%s\n", name);
  const Definition d = chamber(source, turbulent);
  Definition fine = d;
  fine.nz <<= 1;
  fine.nr <<= 1;
  Flow flow(d, Mesh::coarsened(Mesh(fine)));
  Spray uniform(flow, oxygen(flow), settings());
  Amr amr(d, 1);
  Spray spray(amr, oxygen(amr.level(0)), settings());
  amr.setAfterStep([&](int l, double dt) { spray.afterStep(l, dt); });
  amr.setRedo([&] { spray.redo(); });
  long steps = 0, differ = 0;
  const auto start = std::chrono::steady_clock::now();
  double uniformWall = 0;
  while (amr.time() < end) {
    const double t = amr.time();
    const double a = amr.step(std::min(end - t, spray.nextEvent(t) - t));
    const auto u0 = std::chrono::steady_clock::now();
    const double b = flow.step(std::min(end - flow.time(), uniform.nextEvent(flow.time()) - flow.time()));
    uniform.advance(b);
    uniformWall += seconds(u0);
    if (a != b) ++differ;
    ++steps;
  }
  const Flow& base = amr.level(0);
  std::printf("  %ld base steps (%llu redone), t %.17g against %.17g; %.1f s of wall\n", steps,
              static_cast<unsigned long long>(amr.rejectedSteps()), amr.time(), flow.time(), seconds(start));
  check(differ == 0 && amr.time() == flow.time(), "every step's dt and the end time equal (steps that differ)", static_cast<double>(differ), 0);
  const bool gas = same(base.state(), flow.state()) && same(base.partialDensities(), flow.partialDensities()) &&
                   same(base.turbulence(), flow.turbulence());
  check(gas, "final gas state, partial densities and turbulence bit-identical (largest relative difference)",
        stateDifference(base, flow), 0);
  check(same(bytes(spray.parcels()), bytes(uniform.parcels())),
        "parcel list bit-identical, every field (" + std::to_string(spray.parcels().size()) + " against " +
            std::to_string(uniform.parcels().size()) + " parcels)",
        same(bytes(spray.parcels()), bytes(uniform.parcels())) ? 0.0 : 1.0, 0);
  const auto& l = spray.ledger();
  const auto& r = uniform.ledger();
  const bool ledger = l.injectedMass == r.injectedMass && l.injectedEnergy == r.injectedEnergy && l.exchangedMass == r.exchangedMass &&
                      l.exchangedEnergy == r.exchangedEnergy && l.exchangedMomentum == r.exchangedMomentum &&
                      l.outflowMass == r.outflowMass && l.wallImpulse == r.wallImpulse && l.wallLoss == r.wallLoss &&
                      l.substeps == r.substeps && l.children == r.children && l.wallThrust == r.wallThrust;
  std::printf("      printed: ledgers equal (injected, exchanged, outflow, walls, sub-steps, children, wall thrust): %s;"
              " evaporated %.17g kg; level moves %ld finer, %ld coarser\n",
              ledger ? "yes" : "no", l.exchangedMass, l.finerMoves, l.coarserMoves);
  // Printed beside, not judged: the S3 run's own driver (ReactingFlow, frozen) on the same mesh, and on the
  // directly built 16 x 6 mesh.
  for (const bool coarsened : {true, false}) {
    const auto made = coarsened ? std::make_unique<Flow>(d, Mesh::coarsened(Mesh(fine))) : std::make_unique<Flow>(d);
    Flow& other = *made;
    thermo::ReactingFlow reacting(other, "h2o2.yaml", 1, 1e-6, 1e-12, thermo::Chemistry::Frozen);
    Spray sprayOther(other, oxygen(other), settings());
    reacting.setSpray(&sprayOther);
    while (other.time() < end) reacting.step(end - other.time());
    std::printf("      printed, not judged: ReactingFlow's frozen driver on the %s mesh: largest relative difference of the"
                " conserved state %.3e; evaporated %.17g kg, %zu parcels\n",
                coarsened ? "same (coarsened)" : "directly built", stateDifference(base, other), sprayOther.ledger().exchangedMass,
                sprayOther.parcels().size());
  }
  std::printf("      printed: the uniform run's spray %.2f s of %.2f s\n", uniform.wallTime(), uniformWall);
}

// A2: two levels over S3's chamber, a fixed injector-corner tag and a band sweeping down and back.
struct Run {
  std::vector<std::vector<Conserved>> state;
  std::vector<std::vector<double>> partial;
  std::vector<unsigned char> parcels;
};
Run a2(thermo::ReactionSource& source, int threads, double end) {
  std::printf("A2 on %d thread%s: 16 x 6 base plus 2 levels (64 x 24 finest), to %.2g ms\n", threads, threads > 1 ? "s" : "",
              end * 1e3);
  const Definition d = chamber(source, false);
  Amr amr(d, 2);
  amr.setThreads(threads);
  Amr::Tags tags;
  tags.interval = 2;
  tags.throatHalfWidth = 0;
  tags.wallYPlus = 0;
  tags.mixture = tags.temperature = tags.pressure = tags.vorticity = 0;
  const double period = 0.5e-3;
  tags.custom = [&amr, period](int l, int i, int j) {
    const Mesh& m = amr.level(l).mesh();
    const double z = (i + 0.5) * m.dz, r = m.cells[m.index(i, j)].r;
    if (z < 2e-3 && r < 3e-3) return true;
    const double x = std::fmod(amr.time(), period) / period;   // down to the outlet by half the period, then back
    const double centre = kLength * (1 - std::abs(1 - 2 * x));
    return std::abs(z - centre) < (l == 0 ? 3e-3 : 1.5e-3);
  };
  amr.setTags(tags);
  Spray spray(amr, oxygen(amr.level(0)), settings());
  amr.setAfterStep([&](int l, double dt) { spray.afterStep(l, dt); });
  amr.setRedo([&] { spray.redo(); });
  // The gas's initial content, for the scales (the base holds the composite; the state is uniform).
  const Flow& base = amr.level(0);
  double mass0 = 0, energy0 = 0;
  for (std::size_t q = 0; q < base.state().size(); ++q) {
    mass0 += base.state()[q][0] * base.mesh().cells[q].volume;
    energy0 += base.state()[q][3] * base.mesh().cells[q].volume;
  }
  const auto start = std::chrono::steady_clock::now();
  while (amr.time() < end) {
    const double t = amr.time();
    amr.step(std::min(end - t, spray.nextEvent(t) - t));
  }
  const double wall = seconds(start);
  const auto& l = spray.ledger();
  std::printf("  %llu base steps (%llu redone), %llu regrids that changed the pattern, %.1f s of wall\n",
              static_cast<unsigned long long>(amr.steps()), static_cast<unsigned long long>(amr.rejectedSteps()),
              static_cast<unsigned long long>(amr.regrids()), wall);
  check(amr.regrids() >= 50, "it exercises regrids: pattern-changing regrids (at least 50)", static_cast<double>(amr.regrids()), 50);
  check(l.finerMoves >= 100, "it exercises moves: parcel moves to a finer owner within a step (at least 100)",
        static_cast<double>(l.finerMoves), 100);
  check(l.coarserMoves >= 100, "it exercises moves: parcel moves to a coarser owner within a step (at least 100)",
        static_cast<double>(l.coarserMoves), 100);
  const auto b = amr.budget();
  double exchangedMass = 0, exchangedEnergy = 0, exchangedMomentum = 0;
  for (int k = 0; k < amr.levels(); ++k) {
    const auto& m = amr.level(k).measurements();
    exchangedMass += m.exchangedMass;
    exchangedEnergy += m.exchangedEnergy;
    exchangedMomentum += m.exchangedMomentum;
  }
  const double massScale = mass0 + l.injectedMass;
  const double energyScale = std::abs(energy0) + std::abs(l.injectedEnergy) + std::abs(l.outflowEnergy) + l.wallLoss;
  const double momentumScale = b.momentumGross + std::abs(l.injectedMomentum) + std::abs(l.outflowMomentum) + std::abs(l.wallImpulse);
  check(std::abs(b.massResidual) / massScale <= 1e-11, "(a) composite gas budget with the exchange: mass", std::abs(b.massResidual) / massScale, 1e-11);
  check(std::abs(b.energyResidual) / energyScale <= 1e-11, "(a) composite gas budget with the exchange: energy",
        std::abs(b.energyResidual) / energyScale, 1e-11);
  check(std::abs(b.momentumResidual) / momentumScale <= 1e-11, "(a) composite gas budget with the exchange: axial momentum",
        std::abs(b.momentumResidual) / momentumScale, 1e-11);
  const auto e = spray.balanceErrors();
  check(std::abs(e[0]) <= 1e-11, "(b) liquid ledger: mass", std::abs(e[0]), 1e-11);
  check(std::abs(e[1]) <= 1e-11, "(b) liquid ledger: energy", std::abs(e[1]), 1e-11);
  check(std::abs(e[2]) <= 1e-11, "(b) liquid ledger: axial momentum", std::abs(e[2]), 1e-11);
  // (c) gas plus liquid with no exchange term: the gas's residual with its exchange added back, plus the liquid's.
  const double systemMass = b.massResidual + exchangedMass + spray.liquidMass() - (l.injectedMass - l.outflowMass);
  const double systemEnergy = b.energyResidual + exchangedEnergy + spray.liquidEnergy() - (l.injectedEnergy - l.outflowEnergy - l.wallLoss);
  const double systemMomentum =
      b.momentumResidual + exchangedMomentum + spray.liquidMomentum() - (l.injectedMomentum - l.outflowMomentum + l.wallImpulse);
  check(std::abs(systemMass) / massScale <= 1e-11, "(c) whole system: mass", std::abs(systemMass) / massScale, 1e-11);
  check(std::abs(systemEnergy) / energyScale <= 1e-11, "(c) whole system: energy", std::abs(systemEnergy) / energyScale, 1e-11);
  check(std::abs(systemMomentum) / momentumScale <= 1e-11, "(c) whole system: axial momentum", std::abs(systemMomentum) / momentumScale, 1e-11);
  std::printf("      printed: exchanged mass by the gas %.17g kg, by the liquid ledger %.17g kg\n", exchangedMass, l.exchangedMass);
  std::printf("      printed, not judged: evaporated %.6e kg, %zu parcels held, %ld injected, %ld children, %ld sub-steps;"
              " level moves %ld finer, %ld coarser; spray %.2f s of %.2f s (%.1f%%)\n",
              l.exchangedMass, spray.parcels().size(), l.parcelsInjected, l.children, l.substeps, l.finerMoves, l.coarserMoves,
              spray.wallTime(), wall, 100 * spray.wallTime() / wall);
  for (int k = 1; k < amr.levels(); ++k) {
    long n = 0;
    for (const auto a : amr.active(k)) n += a;
    std::printf("      printed: level %d ends with %ld active cells of %zu\n", k, n, amr.active(k).size());
  }
  Run out;
  for (int k = 0; k < amr.levels(); ++k) {
    out.state.push_back(amr.level(k).state());
    out.partial.push_back(amr.level(k).partialDensities());
  }
  out.parcels = bytes(spray.parcels());
  return out;
}
}  // namespace

int main(int argc, char** argv) {
  try {
    thermo::ReactionSource source("h2o2.yaml");
    const std::string part = argc > 1 ? argv[1] : "";
    if (part.empty() || part == "a1") {
      a1(source, false, 0.5e-3);
      a1(source, true, 0.2e-3);
    }
    if (part.empty() || part == "a2") {
      const Run one = a2(source, 1, 0.5e-3);
      const Run four = a2(source, 4, 0.5e-3);
      std::printf("A3: A2 on 1 and 4 threads\n");
      bool gas = true;
      for (std::size_t k = 0; k < one.state.size(); ++k) gas = gas && same(one.state[k], four.state[k]) && same(one.partial[k], four.partial[k]);
      check(gas, "every level's state and partial densities bit-identical", gas ? 0.0 : 1.0, 0);
      check(one.parcels == four.parcels, "parcel list bit-identical, every field", one.parcels == four.parcels ? 0.0 : 1.0, 0);
    }
  } catch (const std::exception& e) {
    std::printf("exception: %s\n", e.what());
    return 2;
  }
  std::printf("%s (%d failed)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
