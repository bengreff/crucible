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
// Usage: crucible_chamber_study <eq|eqt|eqtt|fr|frozen|frt|frp> <nz> <nr> <end time s> <threads> <output prefix>
//        [igniter energy J] [igniter duration s]
// "eqt" is "eq" with the equilibrium from Table A (docs/evidence/TABLE_A.md): the table file is
// CRUCIBLE_EQ_TABLE, and CRUCIBLE_EQ_AUDIT (default 0, none) sets the audit's interval in reaction calls.
// "eqtt" (added 5 October 2026) is "eqt" made viscous and turbulent as "frt" (no PaSR closure, and
// eqt's igniter): the turbulent C1 of the wall-function criteria (docs/evidence/WALL_FUNCTIONS.md).
// It reports the first-cell y+ and the wall heat flow and axial force at the end; criterion 4 is not judged.
// CRUCIBLE_WALL_FUNCTIONS (law or printed) and CRUCIBLE_RADIAL_STRETCHING set the wall law and the ring
// clustering for the wall-function criteria 3 and 4. CRUCIBLE_STATE_START = <dump>,<nz>,<nr>,<stretching>
// (criterion 3's shorter reference) starts from a CRUCIBLE_STATE_DUMP of another grid, interpolated, with
// the supply fully open and no igniter; time restarts at zero and the history adds c* every 0.05 ms.
// The history carries the wall heat flow and the wall axial force in every run.
// Writes <prefix>_history.csv (every 2 us), <prefix>_mesh.csv (stations) and <prefix>_field_<us>.csv
// snapshots (<prefix>_field_<us>_failed.csv for the last state of a run that stopped); tools/chamber_plots.py
// renders them.
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
#include <optional>
#include <stdexcept>
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

// Mass-flux-averaged equilibrium stagnation pressure (and static pressure) of a column.
double columnStagnation(const Flow& flow, Cantera::ThermoPhase& gas, int column, double& pStatic) {
  const auto& mesh = flow.mesh();
  double p0 = 0, flux = 0;
  pStatic = 0;
  for (int j = 0; j < mesh.nr; ++j) {
    auto q = mesh.index(column, j);
    auto w = flow.cellPrimitive(q);
    auto y = flow.massFractions(q);
    double g = w.rho * w.uz * mesh.axialArea(column + 1, j);
    p0 += g * stagnationPressure(gas, flow.temperature(q), w.p, y.data(), 0.5 * (w.uz * w.uz + w.ur * w.ur));
    pStatic += g * w.p;
    flux += g;
  }
  pStatic /= flux;
  return p0 / flux;
}

// Criterion 3's start (docs/evidence/WALL_FUNCTIONS.md): a CRUCIBLE_STATE_DUMP of the same contour on another
// grid, mapped column by column in eta = r / R_wall. p, T, u_z, u_r, the mass fractions, k and omega are linear
// in z between the old columns and in eta between the old centroids (the first centroid's values inside it).
// Between the last old centroid and the wall they run to the wall values: u = 0, T = T_w and k = 0, with p and
// the mass fractions held; there omega is the larger of the last centroid's value and the engine's wall rule
// 6 nu_w / (beta1 d^2) at the cell's own wall distance d. With equilibrium chemistry (`gas` given) the
// interpolated composition is then equilibrated at the cell's T and p, elements kept: the hot core's radicals
// held at the wall's 600 K would otherwise recombine in the first steps and double the pressure there.
// The density follows from p and T.
void restoreInterpolated(Flow& flow, const Definition& d, const std::string& path, int nzOld, int nrOld,
                         double stretchingOld, Cantera::ThermoPhase* gas) {
  Definition old = d;
  old.nz = nzOld;
  old.nr = nrOld;
  old.radialStretching = stretchingOld;
  const Mesh from(old);
  const Mesh& to = flow.mesh();
  const Medium& medium = flow.medium();
  const std::size_t ns = medium.size(), nOld = from.cells.size(), nNew = to.cells.size();
  std::vector<Conserved> state(nOld);
  std::vector<double> partial(nOld * ns), turbulence(2 * nOld);
  FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) throw std::runtime_error("cannot read " + path);
  const bool ok = std::fread(state.data(), sizeof(Conserved), nOld, f) == nOld &&
                  std::fread(partial.data(), sizeof(double), partial.size(), f) == partial.size() &&
                  std::fread(turbulence.data(), sizeof(double), turbulence.size(), f) == turbulence.size() &&
                  std::fgetc(f) == EOF;
  std::fclose(f);
  if (!ok) throw std::runtime_error(path + " does not match the stated grid");
  // Old cell values: p, T, u_z, u_r, k, omega, then the mass fractions.
  const std::size_t nv = 6 + ns;
  std::vector<double> value(nOld * nv), eta(nOld);
  for (std::size_t q = 0; q < nOld; ++q) {
    const auto& u = state[q];
    double* v = value.data() + q * nv;
    double sum = 0;
    for (std::size_t k = 0; k < ns; ++k) sum += partial[q * ns + k];
    for (std::size_t k = 0; k < ns; ++k) v[6 + k] = partial[q * ns + k] / sum;
    const double uz = u[1] / u[0], ur = u[2] / u[0], k = turbulence[2 * q] / u[0];
    const double t = medium.temperature(u[3] / u[0] - (uz * uz + ur * ur) / 2 - k, v + 6, 1000.0);
    v[0] = u[0] * medium.gasConstant(v + 6) * t;
    v[1] = t;
    v[2] = uz;
    v[3] = ur;
    v[4] = k;
    v[5] = turbulence[2 * q + 1] / u[0];
  }
  // Cells are placed by the mesh's own coordinates: the column, and eta the ring's mid fraction of the
  // local radius.
  for (int i = 0; i < nzOld; ++i)
    for (int j = 0; j < nrOld; ++j) eta[from.index(i, j)] = 0.5 * (from.fraction[j] + from.fraction[j + 1]);
  // One old column at eta; returns how far into the wall zone eta lies (0 inside the last centroid).
  auto column = [&](int i, double e, double* out) {
    const std::size_t first = from.index(i, 0), last = from.index(i, nrOld - 1);
    if (e >= eta[last]) {
      const double s = std::min(1.0, (e - eta[last]) / (1 - eta[last]));
      const double* v = value.data() + last * nv;
      std::copy(v, v + nv, out);
      out[1] = v[1] + s * (d.wallTemperature - v[1]);
      out[2] = (1 - s) * v[2];
      out[3] = (1 - s) * v[3];
      out[4] = (1 - s) * v[4];
      return s;
    }
    int j = 0;
    while (j + 1 < nrOld && eta[from.index(i, j + 1)] <= e) ++j;
    const std::size_t a = from.index(i, j), b = from.index(i, std::min(j + 1, nrOld - 1));
    const double w = e <= eta[first] ? 0.0 : std::clamp((e - eta[a]) / (eta[b] - eta[a]), 0.0, 1.0);
    for (std::size_t n = 0; n < nv; ++n) out[n] = (1 - w) * value[a * nv + n] + w * value[b * nv + n];
    return 0.0;
  };
  std::vector<Primitive> cells(nNew);
  std::vector<double> fractions(nNew * ns), kOmega(2 * nNew), left(nv), right(nv), diffusion(ns), work;
  for (int i = 0; i < to.nz; ++i)
    for (int j = 0; j < to.nr; ++j) {
      const std::size_t q = to.index(i, j);
      const double e = 0.5 * (to.fraction[j] + to.fraction[j + 1]), x = (i + 0.5) * to.dz / from.dz - 0.5;
      const int i0 = std::clamp(static_cast<int>(std::floor(x)), 0, nzOld - 1), i1 = std::min(i0 + 1, nzOld - 1);
      const double w = std::clamp(x - i0, 0.0, 1.0);
      const double s = std::max(column(i0, e, left.data()), column(i1, e, right.data()));
      for (std::size_t n = 0; n < nv; ++n) left[n] = (1 - w) * left[n] + w * right[n];
      double* y = left.data() + 6;
      const double p = left[0], t = left[1];
      if (gas) {
        gas->setState_TPY(t, p, y);
        gas->equilibrate("TP");
        gas->getMassFractions(y);
      }
      double omega = left[5];
      if (s > 0) {
        const double muW = medium.transport(d.wallTemperature, p, y, diffusion.data(), work).viscosity;
        const double nuW = muW * medium.gasConstant(y) * d.wallTemperature / p, dist = flow.wallDistances()[q];
        omega = std::max(omega, 6 * nuW / (0.075 * dist * dist));  // beta1 0.075, as core/transport.cpp
      }
      cells[q] = {p / (medium.gasConstant(y) * t), left[2], left[3], p};
      std::copy(y, y + ns, fractions.begin() + q * ns);
      kOmega[2 * q] = left[4];
      kOmega[2 * q + 1] = omega;
    }
  flow.setInitialState(cells, fractions);
  flow.setTurbulence(kOmega);
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 7) {
    std::fprintf(stderr, "usage: %s <eq|eqt|eqtt|fr|frozen|frt|frp> <nz> <nr> <end s> <threads> <prefix> [igniter J] [igniter s]\n",
                 argv[0]);
    return 2;
  }
  const std::string mode = argv[1], prefix = argv[6];
  const int nz = std::atoi(argv[2]), nr = std::atoi(argv[3]), threads = std::atoi(argv[5]);
  const double end = std::atof(argv[4]);
  // pasr: the criterion-4 runs (judged); turbulent: those and eqtt.
  const bool pasr = mode == "frt" || mode == "frp", turbulent = pasr || mode == "eqtt", closure = mode == "frp";
  const double igniterEnergy = argc > 7 ? std::atof(argv[7]) : pasr ? 300 : 0.5;
  const double igniterDuration = argc > 8 ? std::atof(argv[8]) : pasr ? 1e-3 : 2e-4;
  const bool tabulated = mode == "eqt" || mode == "eqtt";
  struct Start { std::string path; int nz{}, nr{}; double stretching{}; };
  std::optional<Start> start;
  if (const char* spec = std::getenv("CRUCIBLE_STATE_START")) {
    const std::string text = spec;
    const auto a = text.find(','), b = text.find(',', a + 1), c = text.find(',', b + 1);
    if (a == std::string::npos || b == std::string::npos || c == std::string::npos)
      throw std::invalid_argument("CRUCIBLE_STATE_START: <dump>,<nz>,<nr>,<radial stretching>");
    start = Start{text.substr(0, a), std::stoi(text.substr(a + 1, b - a - 1)), std::stoi(text.substr(b + 1, c - b - 1)),
                  std::stod(text.substr(c + 1))};
  }
  const auto chemistry = mode == "eq" || tabulated ? thermo::Chemistry::LocalEquilibrium
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
  s.ramp = start ? 0.0 : 5e-4;  // a restart has the valve fully open
  d.supplies = {s};
  // Igniter: from 0.2 ms in a 10 mm by 10 mm core 5 mm off the injector face; 0.5 J over 0.2 ms
  // unless given. Local equilibrium burns any premixed gas at once, so in "eq" mode it has no role.
  d.igniter = {0.005, 0.015, 0.010, start ? 0.0 : igniterEnergy, 2e-4, igniterDuration};
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
  // Wall-function criteria 3 and 4 (docs/evidence/WALL_FUNCTIONS.md): CRUCIBLE_WALL_FUNCTIONS = law or
  // printed puts the wall law on every no-slip wall face; CRUCIBLE_RADIAL_STRETCHING clusters the rings
  // at the wall (Definition::radialStretching). Neither set: the run is unchanged.
  if (const char* law = std::getenv("CRUCIBLE_WALL_FUNCTIONS")) {
    const std::string mode = law;
    if (mode != "law" && mode != "printed") throw std::invalid_argument("CRUCIBLE_WALL_FUNCTIONS: law or printed");
    d.turbulence.wallFunctions = true;
    d.turbulence.printedDerivative = mode == "printed";
    std::printf("wall functions: %s (kappa %.4f, B %.3f)\n", mode.c_str(), d.turbulence.wallKappa, d.turbulence.wallB);
  }
  // CRUCIBLE_SST_CORRECTION = hp: the SST property corrections (Definition::Turbulence::propertyCorrections).
  if (const char* correction = std::getenv("CRUCIBLE_SST_CORRECTION")) {
    if (std::string(correction) != "hp") throw std::invalid_argument("CRUCIBLE_SST_CORRECTION: hp");
    d.turbulence.propertyCorrections = true;
    std::printf("SST property corrections: Hasan, Elias, Menter and Pecnik (2025), full form\n");
  }
  if (const char* stretching = std::getenv("CRUCIBLE_RADIAL_STRETCHING")) {
    d.radialStretching = std::atof(stretching);
    std::printf("radial stretching %g\n", d.radialStretching);
  }

  Flow flow(d);
  if (start) {
    auto startSolution = Cantera::newSolution("h2o2.yaml", "", "none");
    restoreInterpolated(flow, d, start->path, start->nz, start->nr, start->stretching,
                        chemistry == thermo::Chemistry::LocalEquilibrium ? startSolution->thermo().get() : nullptr);
    std::printf("start: %s (%dx%d, radial stretching %g) interpolated to %dx%d; supply fully open, no igniter\n",
                start->path.c_str(), start->nz, start->nr, start->stretching, nz, nr);
  }
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", threads, 1e-6, 1e-12, chemistry);
  if (closure) reacting.setMixingClosure(1.0, {"H2", "O2", "H2O"});
  EquilibriumTable table;
  long audit = 0;
  if (tabulated) {
    const char* path = std::getenv("CRUCIBLE_EQ_TABLE");
    if (!path) throw std::runtime_error("eqt needs CRUCIBLE_EQ_TABLE");
    if (const char* a = std::getenv("CRUCIBLE_EQ_AUDIT")) audit = std::atol(a);
    table = EquilibriumTable::load(path);
    reacting.useEquilibriumTable(table, audit);
    std::printf("table %s: %s\naudit every %ld reaction calls\n", path, table.provenance.c_str(), audit);
  }
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
  std::printf("chemistry %s; igniter %.4g J over %.4g ms from 0.2 ms%s\n", mode.c_str(), d.igniter.energy,
              igniterDuration * 1e3,
              pasr ? (closure ? "; viscous, SST, PaSR closure on" : "; viscous, SST, closure off (control)")
              : turbulent ? "; viscous, SST" : "");

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

  auto writeField = [&](double t, const char* tag = "") {
    char name[512];
    std::snprintf(name, sizeof name, "%s_field_%06.0f%s.csv", prefix.c_str(), t * 1e6, tag);
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
               "igniter_energy,mass_budget,energy_budget,steps,wall_heat_flow,wall_axial_force,cstar\n");
  auto historySolution = Cantera::newSolution("h2o2.yaml", "", "none");
  std::vector<double> snapshots = {2e-5, 5e-5, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 7.5e-4, 1e-3, 1.5e-3, 2e-3, 3e-3, 4e-3, 6e-3, 8e-3};
  std::size_t nextSnapshot = 0;
  const auto clock0 = std::chrono::steady_clock::now();
  double nextSample = 0;
  struct Row { double t, outlet, pInj, fVac, heat, force, cstar; };
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
    double cstar = NAN, pStatic = 0;  // a restart's settling test reads c* at each 0.05 ms window's end
    if (start && (rows.size() + 1) % 25 == 0)
      cstar = columnStagnation(flow, *historySolution->thermo(), endColumn, pStatic) * at / m.outletMassFlow;
    std::fprintf(history, "%.9e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.5f,%.2f,%.6e,%.3e,%.3e,%llu,%.9e,%.9e,%.9e\n",
                 m.time, m.supplyMassFlow, m.inletMassFlow, m.outletMassFlow, m.injectorPressure, pEnd / area, fVac,
                 m.supplyMassFlow > 0 ? fVac / (m.supplyMassFlow * kG0) : 0.0, exitMin, tMax, m.igniterEnergy,
                 m.massBalanceError, m.energyBalanceError, static_cast<unsigned long long>(m.steps), m.wallHeatFlow,
                 m.wallAxialForce, cstar);
    rows.push_back({m.time, m.outletMassFlow, m.injectorPressure, fVac, m.wallHeatFlow, m.wallAxialForce, cstar});
  };
  writeField(0);
  // Criterion 3's settling test for a restart (docs/evidence/WALL_FUNCTIONS.md, declared before the runs): 0.05 ms
  // windows (25 samples); per window the means of the wall heat flow and wall axial force, and c*, vacuum Isp and
  // vacuum thrust at its end. A quantity is settled when its last three window values contract, rho = (x3 - x2) /
  // (x2 - x1) in [0, 0.9], with the geometric remainder |x3 - x2| rho / (1 - rho) at most a quarter of its tightest
  // band: 0.005% for c*, Isp and thrust, 0.5% for the heat flow and force. Corrected 6 October 02:50, before the
  // reference ran: a quantity whose last three window values span at most half of that allowance is also settled
  // (window-to-window noise makes rho meaningless there). The run stops when all five are settled at two consecutive
  // windows (02:51, the same restart: one window passed F_vac on a noisy rho).
  struct Window { double cstar, isp, fVac, heat, force; };
  std::vector<Window> windows;
  bool settledStop = false;
  int settledWindows = 0;
  auto settlingWindow = [&]() {
    Window w{rows.back().cstar, rows.back().fVac / (rows.back().outlet * kG0), rows.back().fVac, 0, 0};
    for (std::size_t i = rows.size() - 25; i < rows.size(); ++i) { w.heat += rows[i].heat / 25; w.force += rows[i].force / 25; }
    windows.push_back(w);
    const char* names[] = {"c*", "Isp_vac", "F_vac", "wall heat", "wall force"};
    const double bands[] = {5e-5, 5e-5, 5e-5, 5e-3, 5e-3};
    auto value = [](const Window& v, int k) { return k == 0 ? v.cstar : k == 1 ? v.isp : k == 2 ? v.fVac : k == 3 ? v.heat : v.force; };
    std::printf("window %zu to %.3f ms: c* %.3f m/s, Isp_vac %.4f s, F_vac %.3f N, wall heat %.6e W, wall force %.6e N", windows.size(),
                rows.back().t * 1e3, w.cstar, w.isp, w.fVac, w.heat, w.force);
    if (windows.size() >= 3) {
      int count = 0;
      std::printf("\n   ");
      for (int k = 0; k < 5; ++k) {
        double x1 = value(windows[windows.size() - 3], k), x2 = value(windows[windows.size() - 2], k), x3 = value(w, k);
        double rho = x2 != x1 ? (x3 - x2) / (x2 - x1) : (x3 == x2 ? 0.0 : INFINITY);
        bool contracts = rho >= 0 && rho <= 0.9;
        double remainder = contracts ? std::abs(x3 - x2) * rho / (1 - rho) / std::abs(x3) : INFINITY;
        double span = (std::max({x1, x2, x3}) - std::min({x1, x2, x3})) / std::abs(x3);
        bool ok = (contracts && remainder <= bands[k]) || span <= 0.5 * bands[k];
        count += ok;
        std::printf(" %s rho %+.3f remainder %.2e span %.2e %s;", names[k], rho, remainder, span, ok ? "settled" : "not settled");
      }
      settledWindows = count == 5 ? settledWindows + 1 : 0;
      if (settledWindows == 2) {
        settledStop = true;
        double n = 0, heat = 0, force = 0, fv = 0, isp = 0;
        for (std::size_t i = rows.size() - 50; i < rows.size(); ++i, ++n) {
          heat += rows[i].heat; force += rows[i].force; fv += rows[i].fVac; isp += rows[i].fVac / (rows[i].outlet * kG0);
        }
        std::printf("\nsettling test met at %.3f ms. Means over the last 0.1 ms: c* %.3f m/s (two window ends), Isp_vac %.4f s, "
                    "F_vac %.3f N, wall heat flow %.6e W, wall axial force %.6e N",
                    rows.back().t * 1e3, 0.5 * (windows[windows.size() - 2].cstar + w.cstar), isp / n, fv / n, heat / n, force / n);
      }
    }
    std::printf("\n");
    std::fflush(stdout);
    std::fflush(history);
  };
  auto march = [&](double end) {
  while (flow.time() < end && !settledStop) {
    double target = std::min({end, nextSample, nextSnapshot < snapshots.size() ? snapshots[nextSnapshot] : end});
    if (target <= flow.time()) target = std::min(end, flow.time() + 2e-6);
    reacting.advanceTo(target);
    if (flow.time() >= nextSample) {
      sample();
      nextSample += 2e-6;
      if (start && rows.size() % 25 == 0) settlingWindow();
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
    if (pasr && end >= 8e-3 && end < 12e-3 && !settled() && worstMass < 1e-11 && worstEnergyResidual / energyEnd < 1e-11 &&
        negativeAt < 0 && admissible) {
      std::printf("criterion 4 (c) fails alone at %.3f ms: extending once to 12 ms\n", flow.time() * 1e3);
      march(12e-3);
      writeField(flow.time());
    }
  } catch (const std::exception& e) {
    std::fclose(history);
    std::printf("\nexception at %.6f ms after %ld steps: %s\n", flow.time() * 1e3, reacting.stats().steps, e.what());
    writeField(flow.time(), "_failed");  // the last state the run reached, for the diagnosis
    if (pasr) std::printf("criterion 4 (b): the run did not reach its end  FAIL\n");
    return 1;
  }
  std::fclose(history);
  if (const char* path = std::getenv("CRUCIBLE_STATE_DUMP")) {
    // The final state's raw doubles (conserved, partial densities, turbulence), to compare builds bit for bit.
    FILE* f = std::fopen(path, "wb");
    if (!f) {
      std::fprintf(stderr, "cannot write %s\n", path);
      return 1;
    }
    std::fwrite(flow.state().data(), sizeof(Conserved), flow.state().size(), f);
    std::fwrite(flow.partialDensities().data(), sizeof(double), flow.partialDensities().size(), f);
    std::fwrite(flow.turbulence().data(), sizeof(double), flow.turbulence().size(), f);
    std::fclose(f);
  }
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
  if (d.turbulence.propertyCorrections)
    std::printf("SST property corrections: integral of Phi_k booked %.6e J; S_n floored (F1-weighted) %.3e\n", m.correctionEnergy, m.correctionFloored);
  {
    const auto& st = reacting.stats();
    std::printf("reaction wall %.1f s (%.1f%% of the run), %.3f us per cell update\n", st.reactWall,
                100 * st.reactWall / wall, 1e6 * st.reactWall / (static_cast<double>(st.steps) * flow.state().size()));
    if (tabulated) {
      const double updates = static_cast<double>(st.tableCells + st.tableFallbacks);
      std::printf("table: %ld cell updates from it (%ld clamped at an axis end), %ld outside it took the direct call (%.3e)\n",
                  st.tableCells, st.tableClamped, st.tableFallbacks, st.tableFallbacks / updates);
      if (audit)
        std::printf("audit: %ld cells, max |T_core(table Y) - T_Cantera| %.3f K, max |dY| %.3e; criterion 2(b) "
                    "(max |dT| <= 3 K): %s\n", st.auditedCells, st.auditTemperature, st.auditMassFraction,
                    st.auditTemperature <= 3 ? "pass" : "FAIL");
    }
  }
  std::printf("drift over last 1 ms: outlet mass flow %.2e, injector pressure %.2e, vacuum thrust %.2e\n",
              drift([](const Row& r) { return r.outlet; }), drift([](const Row& r) { return r.pInj; }),
              drift([](const Row& r) { return r.fVac; }));

  // Chamber-end stagnation pressure (mass-flux averaged) and the comparison.
  double pStatic = 0;
  const double p0 = columnStagnation(flow, *historySolution->thermo(), endColumn, pStatic);
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
  // The first-cell y+ at the wall row (laminar estimate) and, for eqtt, the wall heat flow and axial force.
  auto reportWall = [&]() {
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
  };
  if (!pasr) {
    if (turbulent) {
      reportWall();
      std::printf("reported: wall heat flow into the gas %.6e W, wall axial force on the gas %.6e N\n", m.wallHeatFlow,
                  m.wallAxialForce);
    }
    return 0;
  }

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
  reportWall();
  std::printf("reported against the control: light-off, injector pressure %.6f MPa, c* %.2f m/s, vacuum Isp %.2f s\n",
              m.injectorPressure / 1e6, cstarSim, ispSim);
  std::printf("%d failures\n", failed);
  return failed ? 1 : 0;
}
