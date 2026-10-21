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
// CRUCIBLE_CFL overrides Definition::cfl; CRUCIBLE_RECONSTRUCT_T sets Definition::temperatureReconstruction,
// CRUCIBLE_BOUND_H Definition::boundedFaceEnthalpy, and CRUCIBLE_BOUND_B the same judged on the excess above
// the mixing line (Flow::setBoundedExcess); CRUCIBLE_PRANDTL_T overrides the turbulent Prandtl
// number (a diagnostic: Pr_t = Sc_t makes turbulent transport carry the excess enthalpy as one scalar);
// CRUCIBLE_DELTA_H sets Definition::Turbulence::consistentFlux (the turbulent fluxes from the two cells'
// differences of Y_k and h).
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
#include <cctype>
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
    if (std::getenv("CRUCIBLE_RECONSTRUCT_T")) d.temperatureReconstruction = true;
    if (std::getenv("CRUCIBLE_BOUND_H") || std::getenv("CRUCIBLE_BOUND_B")) d.boundedFaceEnthalpy = true;
    if (const char* prt = std::getenv("CRUCIBLE_PRANDTL_T")) d.turbulence.prandtl = std::atof(prt);
    if (std::getenv("CRUCIBLE_DELTA_H")) d.turbulence.consistentFlux = true;
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
                 "table_fallbacks,audit_T,spray_wall_s\n");
    Averages all, first, second;
    const double middle = 0.5 * (averageFrom + endTime);
    const auto t0 = std::chrono::steady_clock::now();
    double nextPrint = 0, nextSample = averageFrom;
    long steps = 0;
    // CRUCIBLE_HOT_TRACE=<K>: the hottest cell, every step, kept in a ring of the last kRing steps; when it
    // first passes <K>, the ring, the next kRing steps and its 9x9 neighbourhood are written to
    // <prefix>_hot_trace.csv and <prefix>_hot_cells.csv (the over-temperature in low, docs/evidence/COST.md).
    // excess = h + |u|^2/2 - (f_H h_H2(287 K) + (1 - f_H) h_LOX): the specific total enthalpy above the
    // propellants' mixing line at the cell's hydrogen element fraction (the fill sits slightly off it).
    const char* hotEnv = std::getenv("CRUCIBLE_HOT_TRACE");
    const double hotT = hotEnv ? std::atof(hotEnv) : 0;
    constexpr std::size_t kRing = 2000;
    std::vector<std::string> ring;
    std::size_t ringNext = 0;
    long hotAfter = -1;
    double hotLast = 0;
    FILE* hotTrace = nullptr;
    const std::size_t ns = medium.size();
    std::vector<double> hydrogenShare(ns, 0);  // kg of H per kg of species k
    for (std::size_t k = 0; k < ns; ++k) {
      const std::string& name = medium.species()[k].name;
      for (std::size_t c = 0; c < name.size(); ++c) {
        if (name[c] != 'H' || (c + 1 < name.size() && std::islower(static_cast<unsigned char>(name[c + 1])))) continue;
        int n = 0;
        while (c + 1 < name.size() && std::isdigit(static_cast<unsigned char>(name[c + 1]))) n = 10 * n + (name[++c] - '0');
        hydrogenShare[k] += (n ? n : 1) * 1.00794 / medium.species()[k].molarMass;
      }
    }
    const double hHydrogen = medium.enthalpy(kHydrogenT, h2.data());
    const double cellDr = radius / nr;
    auto cellLine = [&](std::size_t q) {
      const auto y = flow.massFractions(q);
      const auto w = flow.cellPrimitive(q);
      const double t = flow.temperature(q);
      double fh = 0;
      for (std::size_t k = 0; k < ns; ++k) fh += y[k] * hydrogenShare[k];
      const auto& lg = spray.ledger();
      const double hLox = lg.injectedMass > 0 ? lg.injectedEnergy / lg.injectedMass : 0;
      const double excess = medium.enthalpy(t, y.data()) + 0.5 * (w.uz * w.uz + w.ur * w.ur) - (fh * hHydrogen + (1 - fh) * hLox);
      const double* ex = spray.exchanged().data() + q * (ns + 3);
      double depositMass = 0;
      for (std::size_t k = 0; k < ns; ++k) depositMass += ex[k];
      int count = 0;
      double held = 0;
      for (const auto& pc : spray.parcels())
        if (pc.cell == static_cast<int>(q)) {
          ++count;
          held += pc.count * pc.dropMass(spray.liquid().at(pc.t).density);
        }
      char b[512];
      std::snprintf(b, sizeof b, "%d,%d,%.3f,%.3f,%.2f,%.6e,%.6e,%.3f,%.3f,%.5e,%.5e,%.5e,%.5e,%.5f,%.6e,%.6e,%.6e,%d,%.4e",
                    static_cast<int>(mesh.cells[q].z / mesh.dz), static_cast<int>(mesh.cells[q].r / cellDr),
                    mesh.cells[q].z * 1e3, mesh.cells[q].r * 1e3, t, w.p, w.rho, w.uz, w.ur, y[indexOf(medium, "H2")],
                    y[indexOf(medium, "O2")], y[indexOf(medium, "H2O")], y[oh], fh, excess,
                    depositMass / (w.rho * mesh.cells[q].volume), ex[ns + 2] / (w.rho * mesh.cells[q].volume), count, held);
      return std::string(b);
    };
    const char* hotHeader = "i,j,z_mm,r_mm,T,p,rho,uz,ur,Y_H2,Y_O2,Y_H2O,Y_OH,f_H,excess_J_per_kg,"
                            "spray_mass_per_gas_mass,spray_energy_J_per_kg,parcels,liquid_kg";
    // CRUCIBLE_BUDGET=<t0 s>:<i0>:<i1>:<j0>:<j1>: from the first step past t0, the one-cell budgets of the
    // cells i0..i1 x j0..j1 (Flow::setBudgetProbe, the tracked quantity rho f_H), written every 10 steps to
    // <prefix>_budget.csv: each cell's V, rho V, rho E V, rho f_H V, p V, rho k V, T, f_H and excess with the
    // probe's running sums (face 0 to 3 x Flow::FaceTerm, then Flow::SourceTerm), and h_H2(287 K), h_LOX.
    const char* budgetEnv = std::getenv("CRUCIBLE_BUDGET");
    double budgetFrom = 0;
    int budgetBox[4] = {0, 0, 0, 0};
    if (budgetEnv && std::sscanf(budgetEnv, "%lf:%d:%d:%d:%d", &budgetFrom, &budgetBox[0], &budgetBox[1], &budgetBox[2], &budgetBox[3]) != 5)
      throw std::runtime_error("CRUCIBLE_BUDGET: <t0>:<i0>:<i1>:<j0>:<j1>");
    FILE* budget = nullptr;
    auto budgetLines = [&] {
      const auto& probe = flow.budgetProbe();
      const auto& lg = spray.ledger();
      const double hLox = lg.injectedMass > 0 ? lg.injectedEnergy / lg.injectedMass : 0;
      for (std::size_t c = 0; c < probe.cells.size(); ++c) {
        const std::size_t q = probe.cells[c];
        const double v = mesh.cells[q].volume;
        const auto& u = flow.state()[q];
        const auto y = flow.massFractions(q);
        const auto w = flow.cellPrimitive(q);
        const double t = flow.temperature(q);
        double fh = 0;
        for (std::size_t k = 0; k < ns; ++k) fh += y[k] * hydrogenShare[k];
        double tracked = 0;
        for (std::size_t k = 0; k < ns; ++k) tracked += hydrogenShare[k] * flow.partialDensities()[q * ns + k];
        const double rhoK = flow.turbulence().empty() ? 0 : flow.turbulence()[q * 2];
        const double excess = medium.enthalpy(t, y.data()) + 0.5 * (w.uz * w.uz + w.ur * w.ur) - (fh * hHydrogen + (1 - fh) * hLox);
        std::fprintf(budget, "%.12e,%ld,%d,%d,%.12e,%.15e,%.15e,%.15e,%.15e,%.15e,%.6f,%.9e,%.9e,%.9e,%.9e", flow.time(), steps,
                     static_cast<int>(q / mesh.nr), static_cast<int>(q % mesh.nr), v, u[0] * v, u[3] * v, tracked * v, w.p * v,
                     rhoK * v, t, fh, excess, hHydrogen, hLox);
        for (std::size_t n = 0; n < Flow::kProbeStride; ++n) std::fprintf(budget, ",%.15e", probe.sums[c * Flow::kProbeStride + n]);
        std::fprintf(budget, "\n");
      }
      std::fflush(budget);
    };
    // CRUCIBLE_BOUND_B: the bound judges the excess above the mixing line (Flow::setBoundedExcess), from the
    // first step with injected LOX (h_LOX is the spray's injected energy per kilogram); h before.
    bool excessBound = !std::getenv("CRUCIBLE_BOUND_B");
    while (flow.time() < endTime) {
      if (!excessBound && spray.ledger().injectedMass > 0) {
        const double hLox = spray.ledger().injectedEnergy / spray.ledger().injectedMass;
        std::vector<double> c(ns);
        for (std::size_t k = 0; k < ns; ++k) c[k] = hLox + (hHydrogen - hLox) * hydrogenShare[k];
        flow.setBoundedExcess(c);
        excessBound = true;
      }
      reacting.step(endTime - flow.time());
      ++steps;
      if (hotT > 0 && hotAfter != 0) {
        std::size_t hot = 0;
        for (std::size_t q = 1; q < mesh.cells.size(); ++q)
          if (flow.temperature(q) > flow.temperature(hot)) hot = q;
        char head[96];
        std::snprintf(head, sizeof head, "%.9e,%ld,%.4e,", flow.time(), steps, flow.time() - hotLast);
        hotLast = flow.time();
        std::string line = head + cellLine(hot);
        if (hotTrace) {
          std::fprintf(hotTrace, "%s\n", line.c_str());
          if (--hotAfter == 0) std::fclose(hotTrace);
        } else if (flow.temperature(hot) > hotT) {
          hotTrace = std::fopen((prefix + "_hot_trace.csv").c_str(), "w");
          std::fprintf(hotTrace, "t,step,dt,%s\n", hotHeader);
          for (std::size_t n = 0; n < ring.size(); ++n) std::fprintf(hotTrace, "%s\n", ring[(ringNext + n) % ring.size()].c_str());
          std::fprintf(hotTrace, "%s\n", line.c_str());
          FILE* cellsFile = std::fopen((prefix + "_hot_cells.csv").c_str(), "w");
          std::fprintf(cellsFile, "%s\n", hotHeader);
          const int i0 = static_cast<int>(mesh.cells[hot].z / mesh.dz), j0 = static_cast<int>(mesh.cells[hot].r / cellDr);
          for (int i = std::max(0, i0 - 4); i <= std::min(nz - 1, i0 + 4); ++i)
            for (int jr = std::max(0, j0 - 4); jr <= std::min(nr - 1, j0 + 4); ++jr)
              std::fprintf(cellsFile, "%s\n", cellLine(mesh.index(i, jr)).c_str());
          std::fclose(cellsFile);
          std::printf("hot trace: T %.0f K at step %ld (t %.6f ms), cell %s\n", flow.temperature(hot), steps, flow.time() * 1e3,
                      line.c_str());
          std::fflush(stdout);
          hotAfter = static_cast<long>(kRing);
        } else {
          if (ring.size() < kRing) ring.push_back(line);
          else ring[ringNext] = line;
          ringNext = (ringNext + 1) % kRing;
        }
      }
      if (budgetEnv && !budget && flow.time() >= budgetFrom) {
        std::vector<std::size_t> cells;
        for (int i = budgetBox[0]; i <= budgetBox[1]; ++i)
          for (int jr = budgetBox[2]; jr <= budgetBox[3]; ++jr) cells.push_back(mesh.index(i, jr));
        flow.setBudgetProbe(cells, hydrogenShare);
        budget = std::fopen((prefix + "_budget.csv").c_str(), "w");
        if (!budget) throw std::runtime_error("cannot write " + prefix + "_budget.csv");
        std::fprintf(budget, "t,step,i,j,V,rhoV,rhoEV,rhofHV,pV,rhokV,T,f_H,excess,h_H2,h_LOX");
        for (int face = 0; face < 4; ++face)
          for (int term = 0; term < Flow::kFaceTerms; ++term) std::fprintf(budget, ",f%d_%d", face, term);
        for (int term = 0; term < Flow::kSourceTerms; ++term) std::fprintf(budget, ",s_%d", term);
        std::fprintf(budget, "\n");
        budgetLines();
      } else if (budget && steps % 10 == 0) budgetLines();
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
        std::fprintf(history, "%.6e,%ld,%.4e,%.2f,%.6e,%.6e,%.6e,%.2f,%zu,%.6e,%.6e,%.6e,%.3e,%.3e,%.3e,%.3e,%ld,%ld,%ld,%.3g,%.2f\n",
                     now, steps, m.dt, wall, m.inletMassFlow, m.outletMassFlow, m.injectorPressure, tMax,
                     spray.parcels().size(), spray.liquidMass(), spray.ledger().exchangedMass,
                     spray.ledger().outflowMass, m.massBalanceError, m.energyBalanceError, err[0], err[1],
                     st.tableCells, st.tableClamped, st.tableFallbacks, st.auditTemperature, spray.wallTime());
        std::fflush(history);
        std::printf("t %.3f ms, %ld steps, dt %.3g s, wall %.0f s, outlet %.4g kg/s, T_max %.0f K, %zu parcels, "
                    "liquid %.4g g, spray %.0f s, fallbacks %ld\n",
                    now * 1e3, steps, m.dt, wall, m.outletMassFlow, tMax, spray.parcels().size(),
                    spray.liquidMass() * 1e3, spray.wallTime(), st.tableFallbacks);
        std::fflush(stdout);
        nextPrint += kPrint;
      }
    }
    std::fclose(history);
    if (budget) {
      budgetLines();
      std::fclose(budget);
    }
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
