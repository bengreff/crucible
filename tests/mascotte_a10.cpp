// Mascotte A-10 (docs/evidence/SPRAY.md, S7): the subcritical LOX/GH2 single-element spray flame at
// 1 MPa, run blind against the modelling choices pre-registered in
// docs/validation/MASCOTTE_A10_PREDICTION.md (inputs: docs/validation/MASCOTTE_A10_INPUTS.md).
//
// Geometry: the square 50 by 50 mm chamber as a round one of radius R (default 28.0 mm: 70 rings of
// 0.4 mm put ring faces on the annulus's 2.8 and 6.0 mm), from the injector exit plane (z = 0) to
// 478 mm, constant area, ending in a partially non-reflecting subsonic outlet at 1.0 MPa (the
// convergent and the throat are not modelled). H2 enters through 2.8 <= r < 6.0 mm at 23.7 g/s and
// 287 K total; LOX at 50 g/s and 85 K from the 5 mm post (r < 2.5 mm), the face around it a wall.
// SST, adiabatic no-slip walls, the chamber filled at rest with the overall mixture's complete-
// combustion products at 1 MPa and 2,000 K; supply and LOX valves ramp open over 1 ms.
//
// CRUCIBLE_CFL overrides Definition::cfl.
// Usage: crucible_mascotte_a10 <size> <nz> <nr> <end s> <average from s> <threads> <prefix> <table|none>
//                             [parcels per second] [table audit interval]
//   size: blob (Reitz's blob, KH-RT atomises), or rr:<SMD um>:<q>:<core length mm> (Rosin-Rammler
//   at that Sauter mean diameter and exponent, released from the liquid core's surface over its length).
//   table: a Table A file for LocalEquilibrium (none: the direct Cantera call per cell).
// Writes <prefix>_history.csv (every 0.1 ms), <prefix>_stations.csv (T at the CARS stations,
// time-averaged over [average from, end] and over its two halves), <prefix>_field.csv (time-averaged
// T, Y_OH and rho Y_OH per cell) and <prefix>_oh.csv (the radially integrated OH per unit length).
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/equilibrium_table.hpp"
#include "core/spray.hpp"
using namespace crucible;

namespace {
constexpr double kPi = std::numbers::pi;
const std::string kSaturation = CRUCIBLE_SOURCE_DIR "/tools/o2_saturation_nist.tsv";
constexpr double kLength = 0.478, kPressure = 1.0e6;
constexpr double kPostRadius = 2.5e-3, kAnnulusInner = 2.8e-3, kAnnulusOuter = 6.0e-3;
constexpr double kOxygenFlow = 0.050, kHydrogenFlow = 0.0237, kOxygenT = 85, kHydrogenT = 287;
constexpr double kFillT = 2000, kRamp = 1e-3, kSample = 5e-5, kPrint = 1e-4;
constexpr std::array<double, 3> kStations{0.010, 0.250, 0.410};

std::size_t indexOf(const Medium& m, const char* name) {
  const auto& s = m.species();
  for (std::size_t k = 0; k < s.size(); ++k)
    if (s[k].name == name) return k;
  throw std::runtime_error(std::string("species not in the medium: ") + name);
}

struct Averages {
  double weight{};
  std::vector<double> t, yoh, rhoyoh;
  void add(const Flow& flow, std::size_t oh, double w) {
    const std::size_t n = flow.mesh().cells.size();
    if (t.empty()) t.assign(n, 0), yoh.assign(n, 0), rhoyoh.assign(n, 0);
    const auto& partial = flow.partialDensities();
    const std::size_t ns = flow.medium().size();
    for (std::size_t q = 0; q < n; ++q) {
      const double rho = flow.state()[q][0];
      t[q] += w * flow.temperature(q);
      yoh[q] += w * partial[q * ns + oh] / rho;
      rhoyoh[q] += w * partial[q * ns + oh];
    }
    weight += w;
  }
  [[nodiscard]] double at(const std::vector<double>& v, std::size_t q) const { return v[q] / weight; }
};

// T at axial station z in ring j, linear between the neighbouring columns' centres.
double station(const Mesh& mesh, const Averages& a, double z, int j) {
  const double s = z / mesh.dz - 0.5;
  const int i = std::clamp(static_cast<int>(std::floor(s)), 0, mesh.nz - 2);
  const double w = std::clamp(s - i, 0.0, 1.0);
  return (1 - w) * a.at(a.t, mesh.index(i, j)) + w * a.at(a.t, mesh.index(i + 1, j));
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 9) {
    std::fprintf(stderr,
                 "usage: %s <blob|rr:SMD_um:q:core_mm> <nz> <nr> <end s> <average from s> <threads> <prefix> "
                 "<table|none> [parcels/s] [audit]\n",
                 argv[0]);
    return 2;
  }
  try {
    const std::string size = argv[1], prefix = argv[7], tablePath = argv[8];
    const int nz = std::atoi(argv[2]), nr = std::atoi(argv[3]), threads = std::atoi(argv[6]);
    const double endTime = std::atof(argv[4]), averageFrom = std::atof(argv[5]);
    const double parcelRate = argc > 9 ? std::atof(argv[9]) : 5e5;
    const long audit = argc > 10 ? std::atol(argv[10]) : 0;
    const double radius = std::getenv("CRUCIBLE_MASCOTTE_RADIUS") ? std::atof(std::getenv("CRUCIBLE_MASCOTTE_RADIUS")) : 0.028;

    thermo::ReactionSource source("h2o2.yaml");
    Definition d;
    d.experiment = Case::Chamber;
    d.nz = nz;
    d.nr = nr;
    d.contour = {{0, radius}, {kLength, radius}};
    d.species = source.medium().species();
    const auto h2 = source.massFractions("H2:1"), h2o = source.massFractions("H2O:1");
    // Complete-combustion products of the overall mixture: per kg of H2, MR kg of O2 burn MR/8 kg of H2.
    const double mr = kOxygenFlow / kHydrogenFlow, water = 9 * mr / 8 / (1 + mr);
    d.composition.resize(h2.size());
    for (std::size_t k = 0; k < h2.size(); ++k) d.composition[k] = water * h2o[k] + (1 - water) * h2[k];
    d.backPressure = kPressure;
    d.ambientTemperature = kFillT;
    d.outletRelaxation = 0.5;
    d.transport = thermo::transportFits("h2o2.yaml", 60);  // fits valid down to cold O2 (SPRAY.md, Results: S2)
    d.wallTemperature = 0;
    if (const char* cfl = std::getenv("CRUCIBLE_CFL")) d.cfl = std::atof(cfl);
    Supply s;
    s.innerRadius = kAnnulusInner;
    s.outerRadius = kAnnulusOuter;
    s.massFlow = kHydrogenFlow;
    s.totalTemperature = kHydrogenT;
    s.composition = h2;
    s.opens = 0;
    s.ramp = kRamp;
    s.turbulenceIntensity = 0.05;
    s.viscosityRatio = 10;
    d.supplies = {s};
    d.turbulence.enabled = true;
    {  // the Spalart-Rumsey ambient of tests/chamber_study.cpp, on the fill
      const Medium medium = d.medium();
      std::vector<double> diffusion(d.species.size()), work;
      const double a0 = medium.soundSpeed(kFillT, d.composition.data());
      const double rho0 = kPressure / (medium.gasConstant(d.composition.data()) * kFillT);
      const double mu0 = medium.transport(kFillT, kPressure, d.composition.data(), diffusion.data(), work).viscosity;
      d.turbulence.ambientK = 9e-9 * a0 * a0;
      d.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
    }
    Flow flow(d);
    const bool table = tablePath != "none";
    EquilibriumTable equilibrium;  // outlives the reacting flow
    if (table) equilibrium = EquilibriumTable::load(tablePath);
    thermo::ReactingFlow reacting(flow, "h2o2.yaml", threads, 1e-6, 1e-12, thermo::Chemistry::LocalEquilibrium);
    if (table) reacting.useEquilibriumTable(equilibrium, audit);

    const Medium& medium = flow.medium();
    std::vector<double> oxygen(medium.size(), 0.0);
    oxygen[indexOf(medium, "O2")] = 1;
    Spray::Settings settings;
    Spray::Injector j;
    j.rOuter = kPostRadius;
    j.massFlow = kOxygenFlow;
    j.opens = 0;
    j.ramp = kRamp;
    j.temperature = kOxygenT;
    j.parcelRate = parcelRate;
    if (size == "blob") {
      j.size = Spray::Injector::Size::Blob;
    } else if (size.rfind("rr:", 0) == 0) {
      std::stringstream in(size.substr(3));
      std::string a, b, c;
      std::getline(in, a, ':');
      std::getline(in, b, ':');
      std::getline(in, c, ':');
      j.size = Spray::Injector::Size::RosinRammler;
      j.sauterDiameter = std::stod(a) * 1e-6;
      j.exponent = std::stod(b);
      j.coreLength = std::stod(c) * 1e-3;
    } else {
      throw std::invalid_argument("size: blob or rr:<SMD um>:<q>:<core length mm>");
    }
    settings.injectors = {j};
    settings.restitution = 0.5;
    Spray spray(flow, LiquidOxygen::load(kSaturation, medium.enthalpy(298.15, oxygen.data())), settings);
    reacting.setSpray(&spray);

    const Mesh& mesh = flow.mesh();
    const std::size_t oh = indexOf(medium, "OH");
    const auto liquid = spray.liquid().at(kOxygenT);
    const double area = kPi * (kAnnulusOuter * kAnnulusOuter - kAnnulusInner * kAnnulusInner);
    std::printf("Mascotte A-10: %s, %dx%d (dz %.4g mm, dr %.4g mm), R %.4g mm, L %.4g mm, cfl %.3g, %d threads, %.3g parcels/s\n",
                size.c_str(), nz, nr, mesh.dz * 1e3, radius / nr * 1e3, radius * 1e3, kLength * 1e3, d.cfl, threads, parcelRate);
    std::printf("chemistry LocalEquilibrium, %s; LOX %.4g kg/m3 at %.4g K, %.4g m/s; H2 annulus %.4g mm2\n",
                table ? tablePath.c_str() : "direct Cantera", liquid.density, kOxygenT,
                kOxygenFlow / (liquid.density * kPi * kPostRadius * kPostRadius), area * 1e6);

    FILE* history = std::fopen((prefix + "_history.csv").c_str(), "w");
    if (!history) throw std::runtime_error("cannot write " + prefix + "_history.csv");
    std::fprintf(history,
                 "t,steps,dt,wall_s,inlet,outlet,p_injector,T_max,parcels,liquid_mass,evaporated,liquid_out,"
                 "mass_budget,energy_budget,liquid_mass_err,liquid_energy_err,table_cells,table_clamped,"
                 "table_fallbacks,audit_T\n");
    Averages all, first, second;
    const double middle = 0.5 * (averageFrom + endTime);
    const auto t0 = std::chrono::steady_clock::now();
    double nextPrint = 0, nextSample = averageFrom;
    long steps = 0;
    while (flow.time() < endTime) {
      reacting.step(endTime - flow.time());
      ++steps;
      const double now = flow.time();
      if (now >= nextSample - 1e-12 && now >= averageFrom) {
        all.add(flow, oh, 1);
        (now < middle ? first : second).add(flow, oh, 1);
        nextSample += kSample;
      }
      if (now >= nextPrint - 1e-12 || now >= endTime) {
        const auto m = flow.measurements();
        double tMax = 0;
        for (std::size_t q = 0; q < mesh.cells.size(); ++q) tMax = std::max(tMax, flow.temperature(q));
        const auto& st = reacting.stats();
        const auto err = spray.balanceErrors();
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::fprintf(history, "%.6e,%ld,%.4e,%.2f,%.6e,%.6e,%.6e,%.2f,%zu,%.6e,%.6e,%.6e,%.3e,%.3e,%.3e,%.3e,%ld,%ld,%ld,%.3g\n",
                     now, steps, m.dt, wall, m.inletMassFlow, m.outletMassFlow, m.injectorPressure, tMax,
                     spray.parcels().size(), spray.liquidMass(), spray.ledger().exchangedMass,
                     spray.ledger().outflowMass, m.massBalanceError, m.energyBalanceError, err[0], err[1],
                     st.tableCells, st.tableClamped, st.tableFallbacks, st.auditTemperature);
        std::fflush(history);
        std::printf("t %.3f ms, %ld steps, dt %.3g s, wall %.0f s, outlet %.4g kg/s, T_max %.0f K, %zu parcels, "
                    "liquid %.4g g, fallbacks %ld\n",
                    now * 1e3, steps, m.dt, wall, m.outletMassFlow, tMax, spray.parcels().size(),
                    spray.liquidMass() * 1e3, st.tableFallbacks);
        std::fflush(stdout);
        nextPrint += kPrint;
      }
    }
    std::fclose(history);
    if (all.weight == 0) return 0;

    FILE* f = std::fopen((prefix + "_stations.csv").c_str(), "w");
    std::fprintf(f, "z_mm,r_mm,T,T_first_half,T_second_half\n");
    double drift = 0;
    for (double z : kStations)
      for (int jr = 0; jr < nr; ++jr) {
        const double r = mesh.cells[mesh.index(0, jr)].r;
        const double a = station(mesh, all, z, jr), b = station(mesh, first, z, jr), c = station(mesh, second, z, jr);
        if (r <= 0.020) drift = std::max(drift, std::abs(b - c));
        std::fprintf(f, "%.1f,%.4f,%.2f,%.2f,%.2f\n", z * 1e3, r * 1e3, a, b, c);
      }
    std::fclose(f);
    f = std::fopen((prefix + "_field.csv").c_str(), "w");
    std::fprintf(f, "i,j,z_mm,r_mm,T,Y_OH,rhoY_OH\n");
    for (int i = 0; i < nz; ++i)
      for (int jr = 0; jr < nr; ++jr) {
        const std::size_t q = mesh.index(i, jr);
        std::fprintf(f, "%d,%d,%.4f,%.4f,%.2f,%.6e,%.6e\n", i, jr, mesh.cells[q].z * 1e3, mesh.cells[q].r * 1e3,
                     all.at(all.t, q), all.at(all.yoh, q), all.at(all.rhoyoh, q));
      }
    std::fclose(f);
    // OH per unit length, its peak and the downstream half-maximum.
    f = std::fopen((prefix + "_oh.csv").c_str(), "w");
    std::fprintf(f, "z_mm,OH_kg_per_m\n");
    std::vector<double> line(nz, 0.0);
    for (int i = 0; i < nz; ++i) {
      for (int jr = 0; jr < nr; ++jr) {
        const std::size_t q = mesh.index(i, jr);
        line[i] += all.at(all.rhoyoh, q) * mesh.cells[q].volume / mesh.dz;
      }
      std::fprintf(f, "%.4f,%.6e\n", (i + 0.5) * mesh.dz * 1e3, line[i]);
    }
    std::fclose(f);
    const int peak = static_cast<int>(std::max_element(line.begin(), line.end()) - line.begin());
    int half = peak;
    while (half + 1 < nz && line[half + 1] >= 0.5 * line[peak]) ++half;
    std::printf("averaged %.4g to %.4g ms (%.0f samples); station T drift between halves (r <= 20 mm): %.1f K\n",
                averageFrom * 1e3, endTime * 1e3, all.weight, drift);
    std::printf("OH per unit length: peak at z %.1f mm, half maximum downstream at z %.1f mm\n",
                (peak + 0.5) * mesh.dz * 1e3, (half + 1) * mesh.dz * 1e3);
    const auto ledger = spray.ledger();
    std::printf("liquid: injected %.6e kg, evaporated %.6e kg, out %.6e kg, held %.6e kg; balance errors %.2e %.2e %.2e\n",
                ledger.injectedMass, ledger.exchangedMass, ledger.outflowMass, spray.liquidMass(),
                spray.balanceErrors()[0], spray.balanceErrors()[1], spray.balanceErrors()[2]);
  } catch (const std::exception& e) {
    std::printf("error: %s\n", e.what());
    return 1;
  }
  return 0;
}
