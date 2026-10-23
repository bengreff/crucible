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

