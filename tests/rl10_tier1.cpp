// RL10A-3-3A tier 1, the central case's dry-run build (docs/validation/rl10a33a/tier1_setup.md, its checks G1
// to G7 stated before this code). It builds the case's Definition: the PCHIP wall shared with rl10dev
// (tests/rl10_geometry.hpp), the face of choice 1 (8 rings, each a LOX slot between two post-tip walls and two
// H2 slots, the slots at the elements' own flow areas) and of choice 3 (10% of the H2 transpired at a uniform
// mass flux over the face outside the rings' cups), the supply turbulence at 5% and a pitch length scale. It
// reports the layout, the flows and the injection numbers, maps each supply onto candidate meshes by Flow's
// ring-centre rule, constructs Flow and Spray on each, and exits before any step. No measured value is read.
//
// Usage: crucible_rl10_tier1 [NZxNR ...]   (default 76x8 149x16 1192x128)
// A member's element dimensions replace the central ones through the environment (the skeleton's brackets):
//   CRUCIBLE_T1_POST_ID      d_L [mm], central 3.54
//   CRUCIBLE_T1_TIP          t / d_L, central 0.272
//   CRUCIBLE_T1_GAP          h / d_L, central 0.268
//   CRUCIBLE_T1_RECESS       recess / d_L, central 0.5
//   CRUCIBLE_T1_OUTER_RING   (R_c - r_8) / p, central 0.5 (bracket 0.5 to 1)
//   CRUCIBLE_T1_INNER_RING   pitch (r_1 = r_8 - 7p, central) or clearance (r_1 = 27 D_a / 2 pi)
//   CRUCIBLE_T1_TRANSPIRED   transpired share of the fuel, central 0.10
//   CRUCIBLE_T1_FACE_PRESSURE  G4's nominal face pressure [Pa], 2.8e6 (our development runs', not a measurement)
//   CRUCIBLE_RADIAL_STRETCHING  Definition::radialStretching, 0
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <numbers>
#include <optional>
#include <string>
#include <vector>

#include "adapters/reaction.hpp"
#include "core/flow.hpp"
#include "core/spray.hpp"
#include "tests/rl10_geometry.hpp"

using namespace crucible;

namespace {

constexpr double kPi = std::numbers::pi;
const std::string kSaturation = CRUCIBLE_SOURCE_DIR "/tools/o2_saturation_nist.tsv";
// Protocol section 2 (P1) and its addendum's plenum temperatures; the wall at the bracket's central 567 K.
constexpr double kOxidizer = 13.9525, kMixtureRatio = 5.055, kFuelT = 199.67, kOxidizerT = 101.17, kWallT = 567;
constexpr int kRings = 8, kElements = 216, kPerRing = 27;
constexpr double kAmbientP = 1000, kAmbientT = 300, kRamp = 5e-4;
constexpr double kIntensity = 0.05, kCmu = 0.09;

double env(const char* name, double fallback) {
  const char* v = std::getenv(name);
  return v ? std::atof(v) : fallback;
}

int failures = 0;
void check(const char* name, bool ok, const std::string& detail) {
  if (!ok) ++failures;
  std::printf("%-4s %s  %s\n", name, ok ? "pass" : "FAIL", detail.c_str());
}
std::string fmt(const char* f, double a, double b = 0, double c = 0, double d = 0) {
  char s[256];
  std::snprintf(s, sizeof s, f, a, b, c, d);
  return s;
}
double rel(double a, double b) { return std::abs(a - b) / std::abs(b); }

struct Ring {
  double r{}, w{};
  double h2InInner{}, h2InOuter{}, loxInner{}, loxOuter{}, h2OutInner{}, h2OutOuter{};
};
struct Stream {
  std::string name;
  double inner{}, outer{}, massFlow{}, velocity{}, ratio{};
  bool porous{false};
  [[nodiscard]] double area() const { return kPi * (outer * outer - inner * inner); }
};

}  // namespace

int main(int argc, char** argv) {
  try {
    std::vector<std::pair<int, int>> meshes;
    for (int a = 1; a < argc; ++a) {
      int nz = 0, nr = 0;
      if (std::sscanf(argv[a], "%dx%d", &nz, &nr) != 2) throw std::invalid_argument("mesh: NZxNR");
      meshes.push_back({nz, nr});
    }
    if (meshes.empty()) meshes = {{76, 8}, {149, 16}, {1192, 128}};

    // The element dimensions (tier1_skeleton.md, *Declared ranges*).
    const double dL = env("CRUCIBLE_T1_POST_ID", 3.54) * 1e-3;
    const double t = env("CRUCIBLE_T1_TIP", 0.272) * dL, h = env("CRUCIBLE_T1_GAP", 0.268) * dL;
    const double recess = env("CRUCIBLE_T1_RECESS", 0.5) * dL;
    const double outerOffset = env("CRUCIBLE_T1_OUTER_RING", 0.5);
    const std::string innerRule = std::getenv("CRUCIBLE_T1_INNER_RING") ? std::getenv("CRUCIBLE_T1_INNER_RING") : "pitch";
    if (innerRule != "pitch" && innerRule != "clearance") throw std::invalid_argument("CRUCIBLE_T1_INNER_RING: pitch or clearance");
    const double transpired = env("CRUCIBLE_T1_TRANSPIRED", 0.10);
    const double faceP = env("CRUCIBLE_T1_FACE_PRESSURE", 2.8e6);
    const double stretching = env("CRUCIBLE_RADIAL_STRETCHING", 0);
    const double dO = dL + 2 * t, dA = dO + 2 * h;
    const double aL = kPerRing * kPi * dL * dL / 4, aH = kPerRing * kPi * (dA * dA - dO * dO) / 4;

    // The wall at the first mesh's nz gives R_c (the same at every nz).
    const double rc = rl10::wall(meshes.front().first).rc;
    const double pitch = std::sqrt(kPi * rc * rc / kElements);
    const double r8 = rc - outerOffset * pitch;
    const double r1 = innerRule == "pitch" ? r8 - (kRings - 1) * pitch : kPerRing * dA / (2 * kPi);
    const double spacing = (r8 - r1) / (kRings - 1);
    std::printf("RL10A-3-3A tier 1, central case dry run (tier1_setup.md): no step is taken, no measured value is read\n");
    std::printf("element: d_L %.4f mm, t %.4f mm (%.3f d_L), h %.4f mm (%.3f d_L), post OD %.4f mm, annulus OD %.4f mm, "
                "recess %.4f mm (%.3f d_L)\n", dL * 1e3, t * 1e3, t / dL, h * 1e3, h / dL, dO * 1e3, dA * 1e3,
                recess * 1e3, recess / dL);
    std::printf("face: R_c %.4f mm, pitch %.4f mm, r_8 %.4f mm (R_c - %.3f p), r_1 %.4f mm (%s), ring spacing %.4f mm; "
                "per ring LOX area %.4f mm2, H2 area %.4f mm2\n", rc * 1e3, pitch * 1e3, r8 * 1e3, outerOffset,
                r1 * 1e3, innerRule.c_str(), spacing * 1e3, aL * 1e6, aH * 1e6);

    // The rings (choice 1).
    std::vector<Ring> rings(kRings);
    for (int i = 0; i < kRings; ++i) {
      Ring& g = rings[i];
      g.r = r1 + i * spacing;
      g.w = aL / (2 * kPi * g.r);
      g.loxInner = g.r - g.w / 2;
      g.loxOuter = g.r + g.w / 2;
      g.h2InOuter = g.loxInner - t;
      const double inner2 = g.h2InOuter * g.h2InOuter - aH / (2 * kPi);
      g.h2InInner = inner2 > 0 ? std::sqrt(inner2) : -1;
      g.h2OutInner = g.loxOuter + t;
      g.h2OutOuter = std::sqrt(g.h2OutInner * g.h2OutInner + aH / (2 * kPi));
    }
    std::printf("\nring  r [mm]    H2 inner slot [mm]   tip   LOX slot [mm] (w)          tip   H2 outer slot [mm]   cup width [mm]\n");
    for (int i = 0; i < kRings; ++i) {
      const Ring& g = rings[i];
      std::printf("%d  %8.3f   [%7.3f, %7.3f]   %.3f  [%7.3f, %7.3f] (%.4f)   %.3f  [%7.3f, %7.3f]   %.3f\n", i + 1,
                  g.r * 1e3, g.h2InInner * 1e3, g.h2InOuter * 1e3, t * 1e3, g.loxInner * 1e3, g.loxOuter * 1e3,
                  g.w * 1e3, t * 1e3, g.h2OutInner * 1e3, g.h2OutOuter * 1e3, (g.h2OutOuter - g.h2InInner) * 1e3);
    }

    // G2: the layout.
    std::printf("\n");
    bool ordered = rings.front().h2InInner > 0 && rings.back().h2OutOuter < rc;
    for (int i = 0; i < kRings; ++i) {
      const Ring& g = rings[i];
      ordered = ordered && g.h2InInner < g.h2InOuter && g.h2InOuter < g.loxInner && g.loxInner < g.loxOuter &&
                g.loxOuter < g.h2OutInner && g.h2OutInner < g.h2OutOuter;
      if (i + 1 < kRings) ordered = ordered && g.h2OutOuter < rings[i + 1].h2InInner;
    }
    check("G2", ordered, fmt("no overlap: ring 1's inner H2 slot from %.4f mm (axis clear), ring 8's outer to %.4f mm "
                             "(wall at %.4f mm)", rings.front().h2InInner * 1e3, rings.back().h2OutOuter * 1e3, rc * 1e3));
    double areaError = 0, cups = 0;
    for (const Ring& g : rings) {
      areaError = std::max({areaError, rel(kPi * (g.loxOuter * g.loxOuter - g.loxInner * g.loxInner), aL),
                            rel(kPi * (g.h2InOuter * g.h2InOuter - g.h2InInner * g.h2InInner), aH / 2),
                            rel(kPi * (g.h2OutOuter * g.h2OutOuter - g.h2OutInner * g.h2OutInner), aH / 2)});
      cups += kPi * (g.h2OutOuter * g.h2OutOuter - g.h2InInner * g.h2InInner);
    }
    check("G2", areaError < 1e-12, fmt("slot areas: largest relative error %.3g against A_L and A_H / 2", areaError));

    // The streams: 16 H2 slots, 9 porous bands (choice 3).
    const double fuel = kOxidizer / kMixtureRatio, sweep = fuel * transpired;
    std::vector<Stream> streams;
    std::vector<std::array<double, 2>> gaps;
    gaps.push_back({0, rings.front().h2InInner});
    for (int i = 0; i + 1 < kRings; ++i) gaps.push_back({rings[i].h2OutOuter, rings[i + 1].h2InInner});
    gaps.push_back({rings.back().h2OutOuter, rc});
    double porous = 0;
    for (const auto& g : gaps) porous += kPi * (g[1] * g[1] - g[0] * g[0]);
    const double flux = sweep / porous;
    check("G2", rel(cups + porous, kPi * rc * rc) < 1e-12,
          fmt("cups %.2f mm2 + porous %.2f mm2 = pi R_c^2 to %.3g", cups * 1e6, porous * 1e6,
              rel(cups + porous, kPi * rc * rc)));
    for (std::size_t g = 0; g < gaps.size(); ++g) {
      Stream s{"porous " + std::to_string(g), gaps[g][0], gaps[g][1]};
      s.massFlow = flux * s.area();
      s.porous = true;
      streams.push_back(s);
    }
    for (int i = 0; i < kRings; ++i) {
      const double each = (fuel - sweep) / (2 * kRings);
      streams.push_back({"H2 in " + std::to_string(i + 1), rings[i].h2InInner, rings[i].h2InOuter, each});
      streams.push_back({"H2 out " + std::to_string(i + 1), rings[i].h2OutInner, rings[i].h2OutOuter, each});
    }
    std::sort(streams.begin(), streams.end(), [](const Stream& a, const Stream& b) { return a.inner < b.inner; });

    // The medium, the liquid, and G3.
    const std::string mechanism = std::getenv("CRUCIBLE_MECHANISM") ? std::getenv("CRUCIBLE_MECHANISM")
                                                                     : CRUCIBLE_SOURCE_DIR "/tools/mechanisms/h2o2_highT.yaml";
    thermo::ReactionSource source(mechanism);
    const TransportFits fits = thermo::transportFits(mechanism);
    Definition base;  // the medium with its transport data
    base.species = source.medium().species();
    base.transport = fits;
    const Medium medium = base.medium();
    const auto h2 = source.massFractions("H2:1"), o2 = source.massFractions("O2:1");
    const LiquidOxygen liquid = LiquidOxygen::load(kSaturation, medium.enthalpy(298.15, o2.data()));
    const double rhoL = liquid.at(kOxidizerT).density;
    double fuelSum = 0, sweepSum = 0, fluxSpread = 0;
    for (const Stream& s : streams) {
      fuelSum += s.massFlow;
      if (s.porous) {
        sweepSum += s.massFlow;
        fluxSpread = std::max(fluxSpread, rel(s.massFlow / s.area(), flux));
      }
    }
    const double oxSum = kRings * (kOxidizer / kRings);
    check("G3", rel(oxSum, kOxidizer) < 1e-12 && rel(fuelSum, fuel) < 1e-12,
          fmt("oxidizer %.6f kg/s, fuel %.6f kg/s, O/F %.6f, chamber flow %.6f kg/s", oxSum, fuelSum, oxSum / fuelSum,
              oxSum + fuelSum));
    check("G3", rel(sweepSum / fuel, transpired) < 1e-12 && fluxSpread < 1e-12,
          fmt("transpired %.6f of the fuel at %.6f kg/(m2 s) on every porous band (spread %.3g); ring O/F %.5f",
              sweepSum / fuel, flux, fluxSpread, (kOxidizer / kRings) / ((fuel - sweep) / kRings)));

    // G4: the injection numbers at the nominal face pressure, and each stream's viscosity ratio from L = p.
    const double rhoG = faceP / (medium.gasConstant(h2.data()) * kFuelT);
    std::vector<double> diffusion(medium.size()), work;
    const double muG = medium.transport(kFuelT, faceP, h2.data(), diffusion.data(), work).viscosity;
    const double ringOF = (kOxidizer / kRings) / ((fuel - sweep) / kRings);
    std::printf("\nG4 (reported) at a nominal %.3g MPa face pressure, H2 at %.2f K: rho_H2 %.4f kg/m3, mu %.4g Pa s; "
                "LOX at %.2f K: rho_L %.2f kg/m3\n", faceP * 1e-6, kFuelT, rhoG, muG, kOxidizerT, rhoL);
    for (Stream& s : streams) {
      s.velocity = s.massFlow / (rhoG * s.area());
      const double k = 1.5 * std::pow(kIntensity * s.velocity, 2);
      const double omega = std::sqrt(k) / (std::pow(kCmu, 0.25) * pitch);
      s.ratio = rhoG * k / (omega * muG);
    }
    for (int i = 0; i < kRings; ++i) {
      const double vL = (kOxidizer / kRings) / (rhoL * aL);
      const double vG = ((fuel - sweep) / kRings) / (rhoG * aH);
      const double j = rhoG * vG * vG / (rhoL * vL * vL);
      std::printf("  ring %d: V_L %.4f m/s, V_H2 %.3f m/s, J %.3f, (V_g - V_L) / (V_L MR) %.3f\n", i + 1, vL, vG, j,
                  (vG - vL) / (vL * ringOF));
    }
    const auto slot = std::find_if(streams.begin(), streams.end(), [](const Stream& s) { return !s.porous; });
    std::printf("  transpired: %.4f m/s; viscosity ratios from I %.2f and L = p: slots %.4g, porous %.4g\n",
                flux / rhoG, kIntensity, slot->ratio, streams[0].ratio);

    // The Definition (tier1_skeleton.md, *Tier 1's settings*).
    auto define = [&](int nz, int nr, bool merged) {
      Definition d;
      d.experiment = Case::Chamber;
      d.nz = nz;
      d.nr = nr;
      d.contour = rl10::wall(nz).contour;
      d.species = medium.species();
      d.composition = source.massFractions("N2:1");
      d.backPressure = kAmbientP;
      d.ambientTemperature = kAmbientT;
      d.radialStretching = stretching;
      d.transport = fits;
      d.wallTemperature = kWallT;
      d.turbulence.enabled = true;
      d.turbulence.wallFunctions = true;
      d.turbulence.propertyCorrections = true;
      std::vector<double> dd(medium.size()), w;
      const double a0 = medium.soundSpeed(kAmbientT, d.composition.data());
      const double rho0 = kAmbientP / (medium.gasConstant(d.composition.data()) * kAmbientT);
      const double mu0 = medium.transport(kAmbientT, kAmbientP, d.composition.data(), dd.data(), w).viscosity;
      d.turbulence.ambientK = 9e-9 * a0 * a0;
      d.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
      auto supply = [&](const Stream& s) {
        Supply p;
        p.innerRadius = s.inner;
        p.outerRadius = s.outer;
        p.massFlow = s.massFlow;
        p.totalTemperature = kFuelT;
        p.composition = h2;
        p.opens = 0;
        p.ramp = kRamp;
        p.turbulenceIntensity = kIntensity;
        p.viscosityRatio = s.ratio;
        return p;
      };
      for (const Stream& s : streams)
        if (s.porous || !merged) d.supplies.push_back(supply(s));
      if (merged) {  // each ring's two H2 slots as one band across its cup, at the slots' flow
        for (int i = 0; i < kRings; ++i) {
          Stream s{"H2 band", rings[i].h2InInner, rings[i].h2OutOuter, (fuel - sweep) / kRings};
          s.velocity = s.massFlow / (rhoG * s.area());
          const double k = 1.5 * std::pow(kIntensity * s.velocity, 2);
          s.ratio = rhoG * k / (std::sqrt(k) / (std::pow(kCmu, 0.25) * pitch) * muG);
          d.supplies.push_back(supply(s));
        }
        std::sort(d.supplies.begin(), d.supplies.end(),
                  [](const Supply& a, const Supply& b) { return a.innerRadius < b.innerRadius; });
      }
      // rl10dev's igniter, scaled from C1's (tier 1's igniter is open, tier1_skeleton.md *Open*).
      const double scale = rc / 0.025;
      d.igniter = {0.005 * scale, 0.015 * scale, 0.010 * scale, 0.5, 2e-4, 2e-4};
      return d;
    };

    // G1, G5, G6 per mesh; G7 at the finest.
    for (const auto& [nz, nr] : meshes) {
      std::printf("\n=== %d x %d (radial stretching %g)\n", nz, nr, stretching);
      std::optional<rl10::Wall> found;
      try {
        found = rl10::wall(nz);
      } catch (const std::exception& e) {
        check("G1", false, std::string("the contour: ") + e.what());
        continue;
      }
      const rl10::Wall& wall = *found;
      const Definition d = define(nz, nr, false);
      const Mesh mesh(d);
      const double throatR = mesh.radius[wall.throatStation];
      check("G1", rel(throatR, wall.rt) < 1e-12 && std::abs(wall.face - rl10::kNominalFace) <= 0.03 &&
                      rel(mesh.radius[0], rc) < 1e-12 && rel(mesh.radius[nz], wall.rt * std::sqrt(61.0)) < 1e-12,
            fmt("throat station %.0f radius to %.3g; face at %.4f in (%+.4f from -13.1)", wall.throatStation,
                rel(throatR, wall.rt), wall.face, wall.face - rl10::kNominalFace) +
                fmt("; face radius to %.3g; exit radius to %.3g", rel(mesh.radius[0], rc), rel(mesh.radius[nz], wall.rt * std::sqrt(61.0))) +
                fmt("; dz %.4f mm, face ring heights %.4f to %.4f mm",
                    mesh.dz * 1e3, (mesh.fraction[1] - mesh.fraction[0]) * rc * 1e3,
                    (mesh.fraction[nr] - mesh.fraction[nr - 1]) * rc * 1e3));
      // Flow's rule: a face ring belongs to the supply that holds its centre.
      int resolved = 0, empty = 0, worstIndex = -1;
      double worst = 0;
      std::printf("  supply         [inner, outer] mm        true mm2   rings  mesh mm2   A_mesh/A_true  J_mesh/J\n");
      for (std::size_t s = 0; s < streams.size(); ++s) {
        const Stream& st = streams[s];
        double area = 0;
        int count = 0;
        for (int j = 0; j < nr; ++j) {
          const double centre = mesh.ringMiddle(0, j);
          if (centre >= st.inner && centre < st.outer) { area += mesh.axialArea(0, j); ++count; }
        }
        const double ratio = area / st.area();
        if (count == 0) ++empty;
        if (std::abs(ratio - 1) <= 0.05) ++resolved;
        if (count > 0 && std::abs(ratio - 1) > worst) { worst = std::abs(ratio - 1); worstIndex = static_cast<int>(s); }
        std::printf("  %-12s [%8.3f, %8.3f]  %9.2f   %4d  %9.2f   %8.4f   %s\n", st.name.c_str(), st.inner * 1e3,
                    st.outer * 1e3, st.area() * 1e6, count, area * 1e6, ratio,
                    count ? fmt("%.4f", 1 / (ratio * ratio)).c_str() : "-");
      }
      check("G5", resolved == static_cast<int>(streams.size()),
            fmt("%.0f of %.0f supplies within 5%% of their true area; %.0f hold no ring centre", resolved, streams.size(), empty) +
                (worstIndex >= 0 ? " (worst held: " + streams[worstIndex].name + fmt(", %.3f)", 1 + worst) : std::string()));
      // The merged layout's bands: each carries its ring's two slots' flow over the cup, so its J against the slots'
      // is (A_H / A_mesh)^2.
      std::printf("  merged band    [inner, outer] mm      slots mm2   cup mm2   rings  mesh mm2   J_mesh/J\n");
      for (int i = 0; i < kRings; ++i) {
        const Ring& g = rings[i];
        double area = 0;
        int count = 0;
        for (int j = 0; j < nr; ++j) {
          const double centre = mesh.ringMiddle(0, j);
          if (centre >= g.h2InInner && centre < g.h2OutOuter) { area += mesh.axialArea(0, j); ++count; }
        }
        std::printf("  ring %d       [%8.3f, %8.3f]  %9.2f  %8.2f   %4d  %9.2f   %s\n", i + 1, g.h2InInner * 1e3,
                    g.h2OutOuter * 1e3, aH * 1e6, kPi * (g.h2OutOuter * g.h2OutOuter - g.h2InInner * g.h2InInner) * 1e6,
                    count, area * 1e6, count ? fmt("%.4f", std::pow(aH / area, 2)).c_str() : "-");
      }
      for (bool merged : {false, true}) {
        const char* name = merged ? "merged (each ring's H2 slots as one band across its cup)" : "true slots";
        try {
          Definition dm = define(nz, nr, merged);
          Flow flow(dm);
          Spray::Settings settings;
          const Medium& fm = flow.medium();
          for (int i = 0; i < kRings; ++i) {
            Spray::Injector j;
            j.rInner = rings[i].loxInner;
            j.rOuter = rings[i].loxOuter;
            j.massFlow = kOxidizer / kRings;
            j.opens = 0;
            j.ramp = kRamp;
            j.temperature = kOxidizerT;
            settings.injectors.push_back(j);
          }
          settings.restitution = 0.5;
          [[maybe_unused]] Spray spray(flow, LiquidOxygen::load(kSaturation, fm.enthalpy(298.15, o2.data())), settings);
          check("G6", true, std::string(name) + fmt(": Flow and Spray constructed (%.0f supplies, %.0f injectors, %.0f cells)",
                                                    dm.supplies.size(), settings.injectors.size(), flow.mesh().cells.size()));
        } catch (const std::exception& e) {
          check("G6", false, std::string(name) + ": " + e.what());
        }
      }
      if (nz == meshes.back().first && nr == meshes.back().second) {
        // G7 at the finest mesh given: the cup's criterion (tier1_skeleton.md).
        double minSlot = 1e9, minTip = 1e9;
        for (const Ring& g : rings) {
          auto across = [&](double a, double b) {  // the face rings' heights over [a, b]: cells across it
            double n = 0;
            for (int j = 0; j < nr; ++j) {
              const double lo = mesh.fraction[j] * rc, hi = mesh.fraction[j + 1] * rc;
              n += std::max(0.0, std::min(hi, b) - std::max(lo, a)) / (hi - lo);
            }
            return n;
          };
          minSlot = std::min({minSlot, across(g.h2InInner, g.h2InOuter), across(g.loxInner, g.loxOuter),
                              across(g.h2OutInner, g.h2OutOuter)});
          minTip = std::min({minTip, across(g.h2InOuter, g.loxInner), across(g.loxOuter, g.h2OutInner)});
        }
        check("G7", false,
              fmt("the cup: %.2f cells across its depth (need 4), fewest across a slot %.2f and a tip %.2f (need 2); ",
                  recess / mesh.dz, minSlot, minTip) +
                  "Mesh is one body-fitted map of r(z) from the face plane, with no cells behind it, so it cannot carry the cup's walls");
      }
    }
    std::printf("\n%d check(s) failed\n", failures);
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
}
