// SPRAY.md, Spray on AMR, criterion A0: the uniform path unchanged. Writes the raw end state of S3's run (0.5 ms)
// and S8(c)'s (0.2 ms, SST and dispersion), each the uniform Flow + Spray of tests/spray_tests.cpp's chamberRun,
// so the builds before and after the change can be compared with cmp.
// Build against the tree's libraries (docs/evidence/spray_amr/README.md); usage: a0_dump FILE
#include <cstdio>
#include <fstream>
#include <numbers>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/spray.hpp"

using namespace crucible;

namespace {
const std::string kSaturation = CRUCIBLE_SOURCE_DIR "/tools/o2_saturation_nist.tsv";

template <class T>
void put(std::ofstream& out, const std::vector<T>& v) {
  out.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
}

void run(thermo::ReactionSource& source, bool turbulent, double endTime, std::ofstream& out) {
  constexpr double radius = 0.01;
  Definition d;
  d.experiment = Case::Chamber;
  d.nz = 16;
  d.nr = 6;
  d.contour = {{0, radius}, {0.016, radius}};
  d.species = source.medium().species();
  d.composition = source.massFractions("H2:1");
  d.backPressure = 1e6;
  d.ambientTemperature = 300;
  d.transport = thermo::transportFits("h2o2.yaml", 60);
  Supply supply;
  supply.innerRadius = radius / 6;
  supply.outerRadius = radius / 2;
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
  Flow flow(d);
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", 1, 1e-6, 1e-12, thermo::Chemistry::Frozen);
  Spray::Settings settings;
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
  settings.injectors = {j};
  settings.restitution = 0.5;
  std::vector<double> y(flow.medium().size(), 0.0);
  for (std::size_t k = 0; k < y.size(); ++k)
    if (flow.medium().species()[k].name == "O2") y[k] = 1;
  Spray spray(flow, LiquidOxygen::load(kSaturation, flow.medium().enthalpy(298.15, y.data())), settings);
  reacting.setSpray(&spray);
  long steps = 0;
  while (flow.time() < endTime) {
    reacting.step(endTime - flow.time());
    ++steps;
  }
  put(out, flow.state());
  put(out, flow.partialDensities());
  put(out, flow.turbulence());
  for (const auto& p : spray.parcels()) {  // field by field, so a field added later does not change the bytes
    const std::vector<double> x{p.z, p.r, p.uz, p.ur, p.ut, p.d, p.t, p.count, p.shed, p.rtTime, p.distance,
                                p.injectedMass, p.eddyZ, p.eddyR, p.eddyT, p.eddyTime};
    put(out, x);
    put(out, std::vector<std::uint64_t>{p.random});
    put(out, std::vector<int>{p.cell, p.injector});
  }
  const auto& l = spray.ledger();
  const std::vector<double> ledger{l.injectedMass, l.injectedEnergy, l.injectedMomentum, l.exchangedMass, l.exchangedEnergy,
                                   l.exchangedMomentum, l.outflowMass, l.outflowEnergy, l.outflowMomentum, l.wallImpulse,
                                   l.wallLoss, l.droppedSwirl};
  put(out, ledger);
  std::printf("%s: %ld steps, %zu parcels, evaporated %.17g kg\n", turbulent ? "S8(c)" : "S3", steps, spray.parcels().size(),
              l.exchangedMass);
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) { std::printf("usage: a0_dump FILE\n"); return 2; }
  thermo::ReactionSource source("h2o2.yaml");
  std::ofstream out(argv[1], std::ios::binary);
  run(source, false, 0.5e-3, out);
  run(source, true, 0.2e-3, out);
  return 0;
}
