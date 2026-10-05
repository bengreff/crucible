// Verification of the SST-2003 machinery (TECHNICAL_PLAN step 7). Criteria stated before the first
// run (4 October 2026):
//   1. Wall distance (core/walls.cpp). Against references computed without the projection formula:
//      (a) straight duct, equal and wall-clustered rings (b = 2): d = R - r for every centroid;
//      (b) the default converging-diverging nozzle on wall-clustered rings (b = 1.5): the minimum
//          over the side-wall segments of a ternary search of the distance along each segment;
//      (c) the same duct with plate rings 0, 3, 4 and 7 marked as walls: the minimum over the side
//          wall and those plate segments, by the same search;
//      every centroid within 1e-13 of the largest wall radius. (d) No side wall and no plate: every
//      distance infinite.
// Criteria 2 to 4 stated 4 October 2026, before the first run of each. Media: N2 alone (one
// species with the h2o2.yaml thermo and transport fits) for 2 and 4; h2o2.yaml for 3.
//   2. Decaying homogeneous turbulence (the source split and the energy accounting). N2 at rest,
//      300 K and 1e5 Pa, in a straight duct 1 m long and 1 m in radius on 4 x 2 cells with adiabatic
//      slip walls (wall distance infinite, so F1 = 0 and beta = beta2), k0 = 5e4 m^2/s^2 and
//      omega0 = 1e5 1/s. Exact: omega = omega0 / (1 + beta2 omega0 t), k = k0 (1 + beta2 omega0 t)^(-beta*/beta2),
//      at constant rho and E (the decayed k becomes heat). Every step applies the source twice over
//      dt / 2 (Strang). After N = 4, 8, 16, 32 equal steps to beta2 omega0 t = 0.8:
//      (a) the larger relative error of k and omega converges at an observed order of at least 1.9
//          between N = 16 and 32 and is below 1e-4 at N = 32 (the scheme's own value is 2.2e-5,
//          derived from the formula below, not from the engine);
//      (b) k and omega equal the MPRK22 formula of core/transport.cpp evaluated directly here, two
//          half steps per step, within 1e-12 relative;
//      (c) in every cell rho and E are unchanged within 1e-13 relative, the temperature is within
//          1e-8 K of the one with e(T) = e(T0) + k0 - k (the cell's own k), k and omega are uniform
//          within 1e-13 relative, and |u| < 1e-8 m/s;
//      (d) stiff: omega0 = 1e8 1/s, 10 steps of beta2 omega0 dt = 300: k and omega stay positive and
//          finite and decrease at every step.
//   3. Turbulent transport operator. A straight duct (R 1 cm, L 4 cm) with adiabatic slip walls, so
//      the wall distance is infinite and F1 = F2 = 0: mu_t = rho k / omega, sigma_k = 1,
//      sigma_omega = 0.856. Smooth fields that meet the wall conditions the operator imposes (u_r = 0,
//      zero normal derivative of u_z, T, composition, k and omega), with mu_t about 8 mu.
//      (a) The engine's mu_t equals rho k / omega within 1e-12 relative in every cell.
//      (b) The discrete transport divergence of axial momentum, energy, the species, rho k and
//          rho omega, against the exact cell average of the continuum operator (stress
//          (mu + mu_t)(grad u + grad u^T - 2/3 div u I) - 2/3 rho k I with its hoop part; conduction
//          lambda + cp mu_t / Pr_t; every D_km plus mu_t / (rho Sc_t); the k flux in the energy flux;
//          k and omega diffusing with mu + sigma mu_t), by the Gauss quadrature of test 1 of
//          tests/transport_tests.cpp, converges at an observed order of at least 1.8 between the two
//          finest grids (32 x 8 to 256 x 64), volume-weighted L1 off the wall row and the end columns.
//      (c) Radial momentum. The laminar operator already fails this near the axis (order 1.19 on the
//          contoured duct, 1.52 on a straight one; docs/evidence/TRANSPORT_C2.md). The same field
//          without turbulence, against the laminar exact operator, gives the laminar order p_lam;
//          the turbulent radial-momentum order must be at least min(1.8, p_lam - 0.1): turbulence
//          adds no failure of its own. A pass below 1.8 does not clear the laminar failure.
//   4. Fully developed pipe against the independent reference tools/sst_pipe_1d.py (argument --pipe;
//      long, run through the slot). N2, R 0.4 mm, no-slip wall at 300 K, fill 101325 Pa at 300 K,
//      axial body force 2.95e5 N/m^3 (Re_tau about 182), rings clustered with b = 2, nz 4, nr 16, 32
//      and 64 (wall y+ about 0.9, 0.4, 0.2, derived). (First stated as nz 2; the mesh needs at least
//      4 axial cells, so the first launch stopped before any output. The flow is axially uniform.) Each run starts from the reference solution on
//      the same rings (tests/data/sst_pipe_N.csv, written by the tool) and marches 4 ms, about 14
//      times the slowest momentum relaxation time R^2 / (j01^2 nu_t) (estimated). On every grid:
//      (a) steady: u_b and T_axis change by less than 1e-5 relative over the last 1 ms, and the wall
//          shear force equals the body force within 1e-5;
//      (b) budgets: mass below 1e-12; axial momentum and energy (wall shear, body force and its work,
//          wall heat) below 1e-9, the round-off allowance of about 4e6 steps;
//      (c) against the reference on the same rings: u_b within 1e-3 relative, c_f = 2 tau_w / (rho_b u_b^2)
//          within 2e-3, and the axis temperature rise T_axis - T_wall within 1e-2. The two codes share
//          the equations, the rings, the Menter wall omega rule and, in this one-dimensional flow,
//          the same centroid-difference gradients, so the difference measures the implementation,
//          not the grid. The grid error is the model's (first order from the wall omega; reference
//          Richardson limits c_f 0.0104865, u_b 99.439 m/s, T_axis 304.2755 K): reported, not judged.
#include "adapters/reaction.hpp"
#include "core/flow.hpp"
#include "core/walls.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace crucible;

namespace {
int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-58s %.3e (limit %.1e)  %s\n", what, value, limit, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

// Distance from (z, r) to the segment by ternary search on the (convex) squared distance.
double searchDistance(double z, double r, double z0, double r0, double z1, double r1) {
  auto dist2 = [&](double t) { double a = z0 + t * (z1 - z0) - z, b = r0 + t * (r1 - r0) - r; return a * a + b * b; };
  double lo = 0, hi = 1;
  for (int it = 0; it < 200; ++it) {
    double m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3;
    if (dist2(m1) < dist2(m2)) hi = m2; else lo = m1;
  }
  return std::sqrt(std::min({dist2(0.5 * (lo + hi)), dist2(0), dist2(1)}));
}

double referenceWorst(const Mesh& m, const std::vector<double>& d, const std::vector<bool>& plate) {
  double worst = 0, rmax = *std::max_element(m.radius.begin(), m.radius.end());
  for (std::size_t q = 0; q < m.cells.size(); ++q) {
    double z = m.cells[q].z, r = m.cells[q].r, best = std::numeric_limits<double>::infinity();
    for (int i = 0; i < m.nz; ++i) best = std::min(best, searchDistance(z, r, i * m.dz, m.radius[i], (i + 1) * m.dz, m.radius[i + 1]));
    for (std::size_t j = 0; j < plate.size(); ++j)
      if (plate[j]) best = std::min(best, searchDistance(z, r, 0, m.fraction[j] * m.radius[0], 0, m.fraction[j + 1] * m.radius[0]));
    worst = std::max(worst, std::abs(d[q] - best) / rmax);
  }
  return worst;
}

constexpr double pi = std::numbers::pi;
constexpr double betaStar = 0.09, beta2 = 0.0828, sigmaK2 = 1.0, sigmaW2 = 0.856;
double sq(double x) { return x * x; }

// N2 alone: one species with the mechanism's thermo and transport fits (the pair fits are the
// upper triangle, row-major).
struct Nitrogen {
  std::vector<Species> species;
  TransportFits fits;
  Medium medium;
};
Nitrogen nitrogen(const Medium& full, const TransportFits& fits, std::size_t i) {
  const std::size_t n = full.size(), pair = i * n - i * (i - 1) / 2;
  Nitrogen out{{full.species()[i]}, {{fits.viscosity[i]}, {fits.conductivity[i]}, {fits.diffusion[pair]}}, {}};
  out.medium = Medium(out.species);
  out.medium.setTransport(out.fits);
  return out;
}

// One MPRK22 step of the pure decay (no production), as in Flow::turbulenceSource.
void decayStep(double& k, double& w, double tau) {
  const double dk0 = betaStar * w, dw0 = beta2 * w;
  const double k1 = k / (1 + tau * dk0), w1 = w / (1 + tau * dw0);
  const double dk1 = betaStar * w1, dw1 = beta2 * w1;
  k = k / (1 + 0.5 * tau * (dk0 * k + dk1 * k1) / k1);
  w = w / (1 + 0.5 * tau * (dw0 * w + dw1 * w1) / w1);
}

// Test-3 field on a straight duct: every field but u_r has zero normal derivative at the wall,
// and u_r vanishes there (the conditions an adiabatic slip wall imposes).
struct TurbulentField {
  double length = 0.04, radius = 0.01;
  double uz0 = 50, ur0 = 10, t0 = 600, dT = 900, p = 2e5, k0 = 10, w0 = 2e4;
  std::size_t iH2{}, iO2{}, iH2O{}, iOH{}, iN2{};
  double a() const { return 2 * pi / length; }
  double flat(double r) const { return sq(1 - sq(r / radius)); }  // (1 - s^2)^2: 1 on the axis, 0 with zero slope at the wall
  double uz(double z, double r) const { return uz0 * (1 + 0.3 * std::sin(a() * z)) * (1 + 0.5 * flat(r)); }
  double ur(double z, double r) const { double s = r / radius; return ur0 * s * (1 - s * s) * std::cos(a() * z); }
  double t(double z, double r) const { return t0 + dT * flat(r) * (1 + 0.2 * std::cos(a() * z)); }
  double k(double z, double r) const { return k0 * (1 + 0.3 * std::sin(a() * z + 1)) * (1 + 0.5 * (1 - flat(r))); }
  double omega(double z, double r) const { return w0 * (1 + 0.2 * std::cos(a() * z + 2)) * (1 + 0.4 * flat(r)); }
  void y(double z, double r, std::size_t n, double* out) const {
    double s2 = 1 - flat(r), x = a() * z;
    std::fill(out, out + n, 0.0);
    out[iH2] = 0.03 + 0.02 * s2 * std::sin(x);
    out[iO2] = 0.2 + 0.05 * s2 * std::cos(x);
    out[iH2O] = 0.1 + 0.05 * (1 - s2);
    out[iOH] = 0.005 * (1 + s2 * std::sin(x));
    out[iN2] = 1 - out[iH2] - out[iO2] - out[iH2O] - out[iOH];
  }
};

// Exact transport flux per unit area through the unit normal (nz, nr) at a point, with and without
// turbulence: {-(tau.n)_z, -(tau.n)_r, energy, species..., k flux, omega flux}; hoop = tau_thth.
struct ExactTurbulent {
  const TurbulentField& f;
  const Medium& medium;
  double prandtl, schmidt;
  std::size_t n;
  double step = 1e-6;
  mutable std::vector<double> y, x, d, h, work;
  ExactTurbulent(const TurbulentField& field, const Medium& m, double pr, double sc)
      : f(field), medium(m), prandtl(pr), schmidt(sc), n(m.size()), y(n), x(n), d(n), h(n) {}
  void moles(double z, double r, double* out) const {
    f.y(z, r, n, out);
    double sum = 0;
    for (std::size_t k = 0; k < n; ++k) sum += out[k] / medium.species()[k].molarMass;
    for (std::size_t k = 0; k < n; ++k) out[k] = out[k] / medium.species()[k].molarMass / sum;
  }
  template <class F> std::array<double, 2> gradient(F g, double z, double r) const {
    return {(g(z + step, r) - g(z - step, r)) / (2 * step), (g(z, r + step) - g(z, r - step)) / (2 * step)};
  }
  void flux(double z, double r, double nz, double nr, std::vector<double>& turbulent, std::vector<double>& laminar,
            double* hoopTurbulent, double* hoopLaminar) const {
    auto uzg = gradient([&](double a, double b) { return f.uz(a, b); }, z, r);
    auto urg = gradient([&](double a, double b) { return f.ur(a, b); }, z, r);
    auto tg = gradient([&](double a, double b) { return f.t(a, b); }, z, r);
    auto kg = gradient([&](double a, double b) { return f.k(a, b); }, z, r);
    auto wg = gradient([&](double a, double b) { return f.omega(a, b); }, z, r);
    const double t = f.t(z, r), uz = f.uz(z, r), ur = f.ur(z, r), k = f.k(z, r), omega = f.omega(z, r);
    f.y(z, r, n, y.data());
    auto props = medium.transport(t, f.p, y.data(), d.data(), work);
    auto thermo = medium.properties(t, y.data());
    const double rho = f.p / (thermo.r * t), mut = rho * k / omega, cp = thermo.cv + thermo.r;
    double mean = 0;
    for (std::size_t j = 0; j < n; ++j) mean += y[j] / medium.species()[j].molarMass;
    std::vector<std::array<double, 2>> xg(n);
    for (std::size_t j = 0; j < n; ++j) xg[j] = gradient([&](double a, double b) { moles(a, b, x.data()); return x[j]; }, z, r);
    medium.speciesEnthalpies(t, h.data());
    const double div = uzg[0] + urg[1] + ur / r;
    for (int pass = 0; pass < 2; ++pass) {
      const bool turb = pass == 0;
      auto& out = turb ? turbulent : laminar;
      out.assign(3 + n + (turb ? 2 : 0), 0.0);
      const double mu = props.viscosity + (turb ? mut : 0), normal = turb ? 2.0 / 3 * rho * k : 0;
      const double tzz = mu * (2 * uzg[0] - 2.0 / 3 * div) - normal, trr = mu * (2 * urg[1] - 2.0 / 3 * div) - normal;
      const double tzr = mu * (uzg[1] + urg[0]);
      double& hoop = *(turb ? hoopTurbulent : hoopLaminar);
      hoop = mu * (2 * ur / r - 2.0 / 3 * div) - normal;
      const double fz = tzz * nz + tzr * nr, fr = tzr * nz + trr * nr;
      double heat = -(props.conductivity + (turb ? cp * mut / prandtl : 0)) * (tg[0] * nz + tg[1] * nr), sum = 0;
      for (std::size_t j = 0; j < n; ++j) {
        const double diffusion = d[j] + (turb ? mut / (rho * schmidt) : 0);
        out[3 + j] = -rho * medium.species()[j].molarMass * mean * diffusion * (xg[j][0] * nz + xg[j][1] * nr);
        sum += out[3 + j];
      }
      for (std::size_t j = 0; j < n; ++j) { out[3 + j] -= y[j] * sum; heat += h[j] * out[3 + j]; }
      if (turb) {
        out[3 + n] = -(props.viscosity + sigmaK2 * mut) * (kg[0] * nz + kg[1] * nr);
        out[4 + n] = -(props.viscosity + sigmaW2 * mut) * (wg[0] * nz + wg[1] * nr);
        heat += out[3 + n];
      }
      out[0] = -fz; out[1] = -fr; out[2] = -(fz * uz + fr * ur) + heat;
    }
  }
};

const double gaussX[5] = {0.04691007703066800, 0.2307653449471585, 0.5, 0.7692346550528415, 0.9530899229693320};
const double gaussW[5] = {0.1184634425280945, 0.2393143352496832, 0.2844444444444444, 0.2393143352496832, 0.1184634425280945};

// Reference solution of tools/sst_pipe_1d.py (--profile): one row per ring, axis first.
struct Profile { std::vector<double> r, u, k, omega, t, rho; };
Profile readProfile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  Profile p;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line)) {
    std::stringstream s(line);
    std::string cell;
    std::vector<double> v;
    while (std::getline(s, cell, ',')) v.push_back(std::stod(cell));
    if (v.size() < 6) throw std::runtime_error("short row in " + path);
    p.r.push_back(v[0]); p.u.push_back(v[1]); p.k.push_back(v[2]); p.omega.push_back(v[3]); p.t.push_back(v[4]); p.rho.push_back(v[5]);
  }
  return p;
}
}  // namespace

void decayTests(const Nitrogen& n2);
void operatorTests(const Medium& medium, const TransportFits& fits, const TurbulentField& field);
void pipeTests(const Nitrogen& n2, const std::vector<int>& grids);

int main(int argc, char** argv) {
  const std::string mechanism = "h2o2.yaml";
  thermo::ReactionSource source(mechanism);
  auto medium = source.medium();
  const auto fits = thermo::transportFits(mechanism);
  medium.setTransport(fits);
  auto index = [&](const char* name) {
    for (std::size_t k = 0; k < medium.size(); ++k) if (medium.species()[k].name == name) return k;
    throw std::runtime_error(std::string("missing species ") + name);
  };
  const Nitrogen n2 = nitrogen(medium, fits, index("N2"));
  if (argc > 1 && std::strcmp(argv[1], "--pipe") == 0) {
    std::vector<int> grids;
    for (int a = 2; a < argc; ++a) grids.push_back(std::atoi(argv[a]));
    if (grids.empty()) grids = {16, 32, 64};
    pipeTests(n2, grids);
    std::printf("%d failures\n", failures);
    return failures == 0 ? 0 : 1;
  }
  std::printf("1. wall distance\n");
  for (double b : {0.0, 2.0}) {
    Definition def;
    def.experiment = Case::UniformDuct;
    def.nz = 12; def.nr = 10; def.radialStretching = b; def.length = 0.05; def.inletRadius = def.exitRadius = def.throatRadius = 0.01;
    Mesh m(def);
    auto d = wallDistance(m, true, {});
    double worst = 0;
    for (const auto& c : m.cells) worst = std::max(worst, std::abs(d[&c - m.cells.data()] - (0.01 - c.r)) / 0.01);
    check(worst < 1e-13, b > 0 ? "(a) straight duct, b = 2: |d - (R - r)| / R" : "(a) straight duct, equal rings: |d - (R - r)| / R", worst, 1e-13);
  }
  Definition nozzle;
  nozzle.nz = 40; nozzle.nr = 12; nozzle.radialStretching = 1.5;
  Mesh m(nozzle);
  double worst = referenceWorst(m, wallDistance(m, true, {}), {});
  check(worst < 1e-13, "(b) nozzle, b = 1.5: against the segment search", worst, 1e-13);
  std::vector<bool> plate(m.nr, false);
  for (int j : {0, 3, 4, 7}) plate[j] = true;
  worst = referenceWorst(m, wallDistance(m, true, plate), plate);
  check(worst < 1e-13, "(c) nozzle with plate rings 0, 3, 4, 7: segment search", worst, 1e-13);
  auto none = wallDistance(m, false, std::vector<bool>(m.nr, false));
  bool allInfinite = std::all_of(none.begin(), none.end(), [](double v) { return std::isinf(v); });
  check(allInfinite, "(d) no walls: every distance infinite (1 = yes)", allInfinite ? 1 : 0, 1);
  decayTests(n2);
  TurbulentField field;
  field.iH2 = index("H2"); field.iO2 = index("O2"); field.iH2O = index("H2O"); field.iOH = index("OH"); field.iN2 = index("N2");
  operatorTests(medium, fits, field);
  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}

Definition n2Duct(const Nitrogen& n2, double length, double radius, int nz, int nr) {
  Definition d;
  d.experiment = Case::UniformDuct;
  d.nz = nz; d.nr = nr; d.length = length; d.inletRadius = d.exitRadius = d.throatRadius = radius;
  d.species = n2.species; d.composition = {1};
  d.transport = n2.fits;
  d.turbulence.enabled = true;
  return d;
}

void decayTests(const Nitrogen& n2) {
  std::printf("2. decaying homogeneous turbulence, N2 at 300 K and 1e5 Pa, slip walls\n");
  const double t0 = 300, p0 = 1e5, k0 = 5e4, w0 = 1e5, y[1] = {1};
  const double rho0 = p0 / (n2.medium.gasConstant(y) * t0), e0 = n2.medium.internalEnergy(t0, y);
  auto run = [&](double omega0, int steps, double dt, std::vector<std::array<double, 2>>* history) {
    Definition d = n2Duct(n2, 1, 1, 4, 2);
    d.totalPressure = d.backPressure = p0; d.totalTemperature = t0;
    d.wallSlip = true; d.wallTemperature = 0;
    d.turbulence.ambientK = k0; d.turbulence.ambientOmega = omega0;
    auto flow = std::make_unique<Flow>(d);
    flow->setInitialState(std::vector<Primitive>(flow->mesh().cells.size(), Primitive{rho0, 0, 0, p0}));
    for (int s = 0; s < steps; ++s) {
      flow->step(dt);
      if (history) history->push_back({flow->turbulence()[0] / flow->state()[0][0], flow->turbulence()[1] / flow->state()[0][0]});
    }
    return flow;
  };
  const double tEnd = 0.8 / (beta2 * w0);
  std::vector<double> errors;
  double formulaWorst = 0, rhoWorst = 0, energyWorst = 0, tWorst = 0, uniformWorst = 0, speedWorst = 0;
  for (int steps : {4, 8, 16, 32}) {
    std::vector<Conserved> start;
    auto flow = run(w0, 0, 1, nullptr);
    start = flow->state();
    flow = run(w0, steps, tEnd / steps, nullptr);
    const double t = flow->time(), growth = 1 + beta2 * w0 * t;
    const double wExact = w0 / growth, kExact = k0 * std::pow(growth, -betaStar / beta2);
    double kFormula = k0, wFormula = w0;
    for (int s = 0; s < steps; ++s) { decayStep(kFormula, wFormula, 0.5 * tEnd / steps); decayStep(kFormula, wFormula, 0.5 * tEnd / steps); }
    const auto& tu = flow->turbulence();
    double error = 0;
    for (std::size_t q = 0; q < start.size(); ++q) {
      const auto& u = flow->state()[q];
      const double k = tu[2 * q] / u[0], w = tu[2 * q + 1] / u[0];
      error = std::max({error, std::abs(k / kExact - 1), std::abs(w / wExact - 1)});
      formulaWorst = std::max({formulaWorst, std::abs(k / kFormula - 1), std::abs(w / wFormula - 1)});
      rhoWorst = std::max(rhoWorst, std::abs(u[0] / start[q][0] - 1));
      energyWorst = std::max(energyWorst, std::abs(u[3] / start[q][3] - 1));
      tWorst = std::max(tWorst, std::abs(flow->temperature(q) - n2.medium.temperature(e0 + k0 - k, y, t0)));
      const double k0cell = tu[0] / flow->state()[0][0], w0cell = tu[1] / flow->state()[0][0];
      uniformWorst = std::max({uniformWorst, std::abs(k / k0cell - 1), std::abs(w / w0cell - 1)});
      speedWorst = std::max(speedWorst, std::hypot(u[1], u[2]) / u[0]);
    }
    std::printf("  N %-3d t %.6e s (target %.6e)  k %.6e  omega %.6e  error %.3e  T %.4f K\n", steps, t, tEnd,
                tu[0] / flow->state()[0][0], tu[1] / flow->state()[0][0], error, flow->temperature(0));
    errors.push_back(error);
  }
  const double order = std::log2(errors[2] / errors[3]);
  check(order >= 1.9, "2a. observed order, N 16 -> 32", order, 1.9);
  check(errors[3] < 1e-4, "2a. relative error at N 32", errors[3], 1e-4);
  check(formulaWorst < 1e-12, "2b. against the MPRK22 formula, worst relative", formulaWorst, 1e-12);
  check(rhoWorst < 1e-13, "2c. rho unchanged, worst relative", rhoWorst, 1e-13);
  check(energyWorst < 1e-13, "2c. E unchanged, worst relative", energyWorst, 1e-13);
  check(tWorst < 1e-8, "2c. T against e(T0) + k0 - k [K]", tWorst, 1e-8);
  check(uniformWorst < 1e-13, "2c. k and omega uniform, worst relative", uniformWorst, 1e-13);
  check(speedWorst < 1e-8, "2c. largest |u| [m/s]", speedWorst, 1e-8);
  // Stiff: beta2 omega0 dt = 300.
  const double stiff = 1e8;
  std::vector<std::array<double, 2>> history;
  auto flow = run(stiff, 10, 300 / (beta2 * stiff), &history);
  bool ok = true;
  double previousK = k0, previousW = stiff;
  for (const auto& hw : history) {
    ok = ok && std::isfinite(hw[0]) && std::isfinite(hw[1]) && hw[0] > 0 && hw[1] > 0 && hw[0] < previousK && hw[1] < previousW;
    previousK = hw[0]; previousW = hw[1];
  }
  const double growth = 1 + beta2 * stiff * flow->time();
  std::printf("  stiff: after 10 steps (t %.3e s) k %.4e (exact %.4e), omega %.4e (exact %.4e), %llu rejected\n", flow->time(),
              history.back()[0], k0 * std::pow(growth, -betaStar / beta2), history.back()[1], stiff / growth,
              static_cast<unsigned long long>(flow->measurements().rejectedSteps));
  check(ok, "2d. stiff: positive, finite, decreasing every step (1 = yes)", ok ? 1 : 0, 1);
}

void operatorTests(const Medium& medium, const TransportFits& fits, const TurbulentField& field) {
  std::printf("3. turbulent transport operator on a straight slip-wall duct (h2o2.yaml)\n");
  const std::size_t n = medium.size();
  Definition::Turbulence settings;
  ExactTurbulent exact(field, medium, settings.prandtl, settings.schmidt);
  const char* names[6] = {"axial momentum", "radial momentum", "energy", "species", "rho k", "rho omega"};
  std::vector<std::array<double, 6>> errors;
  std::vector<std::array<double, 4>> laminarErrors;
  double eddyWorst = 0;
  for (int nz : {32, 64, 128, 256}) {
    const int nr = nz / 4;
    Definition d;
    d.experiment = Case::UniformDuct;
    d.nz = nz; d.nr = nr; d.length = field.length; d.inletRadius = d.exitRadius = d.throatRadius = field.radius;
    d.species = medium.species();
    d.composition.assign(n, 0.0); d.composition[field.iN2] = 1;
    d.totalPressure = d.backPressure = field.p;
    d.transport = fits;
    d.wallSlip = true; d.wallTemperature = 0;
    Flow laminar(d);
    d.turbulence.enabled = true; d.turbulence.ambientK = field.k0; d.turbulence.ambientOmega = field.w0;
    Flow turbulent(d);
    const auto& m = turbulent.mesh();
    std::vector<Primitive> cells(m.cells.size());
    std::vector<double> y(m.cells.size() * n), kOmega(2 * m.cells.size());
    for (std::size_t q = 0; q < m.cells.size(); ++q) {
      const auto& c = m.cells[q];
      field.y(c.z, c.r, n, y.data() + q * n);
      double t = field.t(c.z, c.r);
      cells[q] = {field.p / (medium.gasConstant(y.data() + q * n) * t), field.uz(c.z, c.r), field.ur(c.z, c.r), field.p};
      kOmega[2 * q] = field.k(c.z, c.r); kOmega[2 * q + 1] = field.omega(c.z, c.r);
    }
    laminar.setInitialState(cells, y);
    turbulent.setInitialState(cells, y);
    turbulent.setTurbulence(kOmega);
    std::vector<Conserved> rate, rateLaminar;
    std::vector<double> speciesRate, speciesLaminar, turbulenceRate;
    turbulent.transportDerivative(rate, speciesRate, turbulenceRate);
    laminar.transportDerivative(rateLaminar, speciesLaminar);
    for (std::size_t q = 0; q < m.cells.size(); ++q) {
      const double rho = turbulent.state()[q][0], k = turbulent.turbulence()[2 * q] / rho, w = turbulent.turbulence()[2 * q + 1] / rho;
      eddyWorst = std::max(eddyWorst, std::abs(turbulent.eddyViscosities()[q] / (rho * k / w) - 1));
    }
    std::array<double, 6> num{}, den{};
    std::array<double, 4> numL{}, denL{};
    std::vector<double> fluxT, fluxL;
    for (int i = 1; i < nz - 1; ++i)
      for (int j = 0; j < nr - 1; ++j) {
        const auto q = m.index(i, j);
        const auto& c = m.cells[q];
        std::vector<double> total(5 + n, 0.0), totalL(3 + n, 0.0);
        auto face = [&](double z0, double r0, double z1, double r1, double nzOut, double nrOut) {
          double len = std::hypot(z1 - z0, r1 - r0);
          if (!(r0 + r1 > 0)) return;
          for (int g = 0; g < 5; ++g) {
            double z = z0 + gaussX[g] * (z1 - z0), r = r0 + gaussX[g] * (r1 - r0), hoopT, hoopL;
            exact.flux(z, r, nzOut, nrOut, fluxT, fluxL, &hoopT, &hoopL);
            for (std::size_t k = 0; k < fluxT.size(); ++k) total[k] -= gaussW[g] * 2 * pi * r * len * fluxT[k];
            for (std::size_t k = 0; k < fluxL.size(); ++k) totalL[k] -= gaussW[g] * 2 * pi * r * len * fluxL[k];
          }
        };
        const double a = m.fraction[j], b = m.fraction[j + 1], z0 = i * m.dz, z1 = (i + 1) * m.dz, R = field.radius;
        face(z0, a * R, z0, b * R, -1, 0);
        face(z1, a * R, z1, b * R, 1, 0);
        face(z0, a * R, z1, a * R, 0, -1);
        face(z0, b * R, z1, b * R, 0, 1);
        for (int g = 0; g < 5; ++g)
          for (int h = 0; h < 5; ++h) {
            double z = z0 + gaussX[g] * m.dz, r = (a + gaussX[h] * (b - a)) * R, hoopT, hoopL;
            exact.flux(z, r, 1, 0, fluxT, fluxL, &hoopT, &hoopL);
            total[1] -= gaussW[g] * gaussW[h] * m.dz * (b - a) * R * 2 * pi * hoopT;
            totalL[1] -= gaussW[g] * gaussW[h] * m.dz * (b - a) * R * 2 * pi * hoopL;
          }
        for (int k = 0; k < 3; ++k) {
          double e = total[k] / c.volume, eL = totalL[k] / c.volume;
          num[k] += c.volume * std::abs(rate[q][k + 1] - e); den[k] += c.volume * std::abs(e);
          numL[k] += c.volume * std::abs(rateLaminar[q][k + 1] - eL); denL[k] += c.volume * std::abs(eL);
        }
        for (std::size_t k = 0; k < n; ++k) {
          double e = total[3 + k] / c.volume, eL = totalL[3 + k] / c.volume;
          num[3] += c.volume * std::abs(speciesRate[q * n + k] - e); den[3] += c.volume * std::abs(e);
          numL[3] += c.volume * std::abs(speciesLaminar[q * n + k] - eL); denL[3] += c.volume * std::abs(eL);
        }
        for (int t = 0; t < 2; ++t) {
          double e = total[3 + n + t] / c.volume;
          num[4 + t] += c.volume * std::abs(turbulenceRate[2 * q + t] - e); den[4 + t] += c.volume * std::abs(e);
        }
      }
    std::array<double, 6> e{};
    std::array<double, 4> eL{};
    std::printf("  %4dx%-3d", nz, nr);
    for (int k = 0; k < 6; ++k) { e[k] = num[k] / den[k]; std::printf("  %s %.3e", names[k], e[k]); }
    std::printf("\n           laminar:");
    for (int k = 0; k < 4; ++k) { eL[k] = numL[k] / denL[k]; std::printf("  %s %.3e", names[k], eL[k]); }
    std::printf("\n");
    errors.push_back(e);
    laminarErrors.push_back(eL);
  }
  check(eddyWorst < 1e-12, "3a. mu_t against rho k / omega, worst relative", eddyWorst, 1e-12);
  double laminarOrder[4];
  for (int k = 0; k < 4; ++k) laminarOrder[k] = std::log2(laminarErrors[2][k] / laminarErrors[3][k]);
  std::printf("  laminar orders (two finest): axial %.3f  radial %.3f  energy %.3f  species %.3f\n", laminarOrder[0],
              laminarOrder[1], laminarOrder[2], laminarOrder[3]);
  for (int k = 0; k < 6; ++k) {
    const double order = std::log2(errors[2][k] / errors[3][k]);
    const double limit = k == 1 ? std::min(1.8, laminarOrder[1] - 0.1) : 1.8;
    char what[128];
    std::snprintf(what, sizeof what, "%s %s: observed order (two finest grids)", k == 1 ? "3c." : "3b.", names[k]);
    check(order >= limit, what, order, limit);
  }
}

void pipeTests(const Nitrogen& n2, const std::vector<int>& grids) {
  std::printf("4. fully developed pipe against tools/sst_pipe_1d.py (N2, R 0.4 mm, wall 300 K, f 2.95e5 N/m^3, b = 2)\n");
  const double radius = 4e-4, tWall = 300, force = 2.95e5, length = 2 * radius, y[1] = {1};
  const double gas = n2.medium.gasConstant(y);
  for (int nr : grids) {
    const auto profile = readProfile(std::string(CRUCIBLE_SOURCE_DIR) + "/tests/data/sst_pipe_" + std::to_string(nr) + ".csv");
    if (profile.r.size() != static_cast<std::size_t>(nr)) throw std::runtime_error("profile ring count");
    Definition d = n2Duct(n2, length, radius, 4, nr);
    d.radialStretching = 2;
    d.totalPressure = d.backPressure = 101325; d.totalTemperature = tWall;
    d.wallTemperature = tWall;
    d.turbulence.ambientK = 0; d.turbulence.ambientOmega = 1;
    Flow flow(d);
    const auto& m = flow.mesh();
    std::vector<Primitive> cells(m.cells.size());
    std::vector<double> kOmega(2 * m.cells.size());
    double centroidWorst = 0;
    for (int i = 0; i < m.nz; ++i)
      for (int j = 0; j < m.nr; ++j) {
        const auto q = m.index(i, j);
        centroidWorst = std::max(centroidWorst, std::abs(m.cells[q].r - profile.r[j]) / radius);
        cells[q] = {profile.rho[j], profile.u[j], 0, profile.rho[j] * gas * profile.t[j]};
        kOmega[2 * q] = profile.k[j]; kOmega[2 * q + 1] = profile.omega[j];
      }
    if (centroidWorst > 1e-12) throw std::runtime_error("reference rings differ from the engine's");
    flow.setInitialState(cells);
    flow.setTurbulence(kOmega);
    flow.setBodyForce(std::vector<std::array<double, 2>>(m.cells.size(), {force, 0}));
    // Reference measures from the profile (its mass is the fill's; its wall shear is f R / 2).
    double refMass = 0, refMomentum = 0;
    for (int j = 0; j < nr; ++j) {
      const double vol = 0.5 * radius * radius * (sq(m.fraction[j + 1]) - sq(m.fraction[j]));
      refMass += profile.rho[j] * vol; refMomentum += profile.rho[j] * profile.u[j] * vol;
    }
    const double refUb = refMomentum / refMass, refRhoB = refMass / (0.5 * radius * radius);
    const double refCf = 2 * (force * radius / 2) / (refRhoB * refUb * refUb), refRise = profile.t[0] - tWall;
    auto sample = [&](double& ub, double& tAxis, double& cf, double& forceBalance) {
      const auto meas = flow.measurements();
      ub = meas.axialMomentum / meas.mass;
      tAxis = 0;
      for (int i = 0; i < m.nz; ++i) tAxis += flow.temperature(m.index(i, 0)) / m.nz;
      const double area = 2 * pi * radius * length, tauW = -meas.wallAxialForce / area;
      cf = 2 * tauW / (meas.mass / (pi * radius * radius * length) * ub * ub);
      forceBalance = (meas.wallAxialForce + meas.bodyAxialForce) / meas.bodyAxialForce;
    };
    const auto wall0 = std::chrono::steady_clock::now();
    double ub = 0, tAxis = 0, cf = 0, balance = 0, ub3 = 0, tAxis3 = 0;
    for (int half = 1; half <= 8; ++half) {
      flow.advanceTo(0.5e-3 * half);
      sample(ub, tAxis, cf, balance);
      if (half == 6) { ub3 = ub; tAxis3 = tAxis; }
      std::printf("  nr %-3d t %.1f ms  u_b %.8f m/s  c_f %.10f  T_axis %.6f K  force balance %+.2e  steps %llu  wall %.0f s\n", nr,
                  flow.time() * 1e3, ub, cf, tAxis, balance, static_cast<unsigned long long>(flow.measurements().steps),
                  std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count());
      std::fflush(stdout);
    }
    const auto meas = flow.measurements();
    // Axial uniformity (reported): largest difference of a cell's u_z from its ring's mean.
    double axialWorst = 0;
    for (int j = 0; j < m.nr; ++j) {
      double mean = 0;
      for (int i = 0; i < m.nz; ++i) mean += flow.cellPrimitive(m.index(i, j)).uz / m.nz;
      for (int i = 0; i < m.nz; ++i) axialWorst = std::max(axialWorst, std::abs(flow.cellPrimitive(m.index(i, j)).uz - mean) / ub);
    }
    std::printf("  nr %-3d reference: u_b %.8f m/s  c_f %.10f  T_axis %.6f K;  axial non-uniformity %.2e;  budgets mass %.2e momentum %.2e energy %.2e\n",
                nr, refUb, refCf, profile.t[0], axialWorst, meas.massBalanceError, meas.momentumBalanceError, meas.energyBalanceError);
    std::printf("  nr %-3d against the reference Richardson limits (reported): u_b %+.3e  c_f %+.3e  T_axis rise %+.3e\n", nr,
                ub / 99.439 - 1, cf / 0.0104865 - 1, (tAxis - tWall) / (304.2755 - tWall) - 1);
    char what[128];
    auto label = [&](const char* text) { std::snprintf(what, sizeof what, "nr %d: %s", nr, text); return what; };
    check(std::abs(ub / ub3 - 1) < 1e-5, label("4a. u_b change over the last 1 ms"), std::abs(ub / ub3 - 1), 1e-5);
    check(std::abs(tAxis / tAxis3 - 1) < 1e-5, label("4a. T_axis change over the last 1 ms"), std::abs(tAxis / tAxis3 - 1), 1e-5);
    check(std::abs(balance) < 1e-5, label("4a. wall shear force against the body force"), std::abs(balance), 1e-5);
    check(std::abs(meas.massBalanceError) < 1e-12, label("4b. mass budget"), std::abs(meas.massBalanceError), 1e-12);
    check(std::abs(meas.momentumBalanceError) < 1e-9, label("4b. axial momentum budget"), std::abs(meas.momentumBalanceError), 1e-9);
    check(std::abs(meas.energyBalanceError) < 1e-9, label("4b. energy budget"), std::abs(meas.energyBalanceError), 1e-9);
    check(std::abs(ub / refUb - 1) < 1e-3, label("4c. u_b against the reference"), std::abs(ub / refUb - 1), 1e-3);
    check(std::abs(cf / refCf - 1) < 2e-3, label("4c. c_f against the reference"), std::abs(cf / refCf - 1), 2e-3);
    check(std::abs((tAxis - tWall) / refRise - 1) < 1e-2, label("4c. T_axis - T_wall against the reference"),
          std::abs((tAxis - tWall) / refRise - 1), 1e-2);
  }
}
