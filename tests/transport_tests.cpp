// C2 verification of the molecular-transport operator (core/transport.cpp) with h2o2.yaml
// mixture-averaged transport. Criteria stated before the first run (4 October 2026):
//   1. Operator order. On a contoured duct (wall slope up to 0.2) carrying smooth fields that satisfy
//      the no-slip isothermal impermeable wall, the discrete transport divergence of axial momentum, radial
//      momentum, energy and the species, against the exact cell average of the continuum operator
//      (Gauss quadrature of the analytic fluxes over the actual cell faces and of the hoop source
//      over the cell, same transport properties), converges at an observed order of at least 1.8
//      between the two finest grids, in the volume-weighted L1 norm over the cells off the wall row
//      and the two end columns. The wall row is reported, not judged (a two-point wall gradient has
//      an O(1) local truncation there; the solution still converges, see 2).
//   2. Pipe decay. In an axially uniform N2 pipe with a no-slip isothermal wall the profile
//      u_z = U J0(j01 r / R) exp(-nu j01^2 t / R^2) is exact. After one e-fold the volume-weighted
//      L1 error converges at an observed order of at least 1.8 between the two finest grids and
//      is below 1e-3 on the finest. The wall shear enters the momentum budget: balance below 1e-12.
//      Test correction (4 October 2026, after the first runs, criterion unchanged): the first
//      composition varied as s^2 and so carried a species flux through the wall, which the operator
//      forbids. The reference integrated that flux and the wall-row gradient fit assumed none, so
//      the energy and species errors grew toward the wall (observed orders 0.50 and 0.28). The
//      composition now varies as 1 - (1 - s^2)^2. Results of both are in docs/evidence/TRANSPORT_C2.md.
//   3. Budgets with transport. A Chamber (closed plate, isothermal side wall and plate, ambient
//      exit) started from the field of test 1: mass, energy (with the wall heat) and axial momentum
//      (with the wall shear) budgets below 1e-12 after 200 steps.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

#include "adapters/reaction.hpp"
#include "core/flow.hpp"

using namespace crucible;

namespace {
int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-58s %.3e (limit %.1e)  %s\n", what, value, limit, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}
constexpr double pi = std::numbers::pi;

// The test-1 field: contour R(z), velocity, temperature and composition at (z, r).
struct Field {
  double length = 0.04, r0 = 0.01, bump = 0.25;
  double uz0 = 50, ur0 = 10, tWall = 600, dT = 900, p = 2e5;
  std::size_t iH2, iO2, iH2O, iOH, iN2;
  double radius(double z) const { return r0 * (1 - bump * std::pow(std::sin(pi * z / length), 2)); }
  double uz(double z, double r) const { double s = r / radius(z); return uz0 * (1 - s * s) * (1 + 0.3 * std::sin(2 * pi * z / length)); }
  double ur(double z, double r) const { double s = r / radius(z); return ur0 * s * (1 - s * s) * std::cos(2 * pi * z / length); }
  double t(double z, double r) const { double s = r / radius(z); return tWall + dT * (1 - s * s) * (1 + 0.2 * std::cos(2 * pi * z / length)); }
  void y(double z, double r, std::size_t n, double* out) const {
    // The composition varies as 1 - (1 - s^2)^2, flat at the wall in both directions: the wall is
    // impermeable (no species flux), as the operator imposes.
    double s2 = 1 - std::pow(1 - std::pow(r / radius(z), 2), 2), a = 2 * pi * z / length;
    std::fill(out, out + n, 0.0);
    out[iH2] = 0.03 + 0.02 * s2 * std::sin(a);
    out[iO2] = 0.2 + 0.05 * s2 * std::cos(a);
    out[iH2O] = 0.1 + 0.05 * (1 - s2);
    out[iOH] = 0.005 * (1 + s2 * std::sin(a));
    out[iN2] = 1 - out[iH2] - out[iO2] - out[iH2O] - out[iOH];
  }
};

// Exact transport flux per unit area through unit normal n at a point: {-(tau.n)_z, -(tau.n)_r,
// -(tau.n).u + q.n} and the species fluxes; also returns tau_thth through hoop.
struct Exact {
  const Field& f;
  const Medium& medium;
  std::size_t n;
  double step = 1e-6;
  mutable std::vector<double> y, x, d, h, work, yp, ym;
  Exact(const Field& field, const Medium& m) : f(field), medium(m), n(m.size()), y(n), x(n), d(n), h(n), yp(n), ym(n) {}
  void moles(double z, double r, double* out) const {
    f.y(z, r, n, out);
    double sum = 0;
    for (std::size_t k = 0; k < n; ++k) sum += out[k] / medium.species()[k].molarMass;
    for (std::size_t k = 0; k < n; ++k) out[k] = out[k] / medium.species()[k].molarMass / sum;
  }
  template <class F> std::array<double, 2> gradient(F g, double z, double r) const {
    return {(g(z + step, r) - g(z - step, r)) / (2 * step), (g(z, r + step) - g(z, r - step)) / (2 * step)};
  }
  std::vector<double> flux(double z, double r, double nz, double nr, double* hoop) const {
    auto uzg = gradient([&](double a, double b) { return f.uz(a, b); }, z, r);
    auto urg = gradient([&](double a, double b) { return f.ur(a, b); }, z, r);
    auto tg = gradient([&](double a, double b) { return f.t(a, b); }, z, r);
    double t = f.t(z, r), uz = f.uz(z, r), ur = f.ur(z, r);
    f.y(z, r, n, y.data());
    auto props = medium.transport(t, f.p, y.data(), d.data(), work);
    double mu = props.viscosity, div = uzg[0] + urg[1] + ur / r;
    double tzz = mu * (2 * uzg[0] - 2.0 / 3 * div), trr = mu * (2 * urg[1] - 2.0 / 3 * div), tzr = mu * (uzg[1] + urg[0]);
    if (hoop) *hoop = mu * (2 * ur / r - 2.0 / 3 * div);
    double fz = tzz * nz + tzr * nr, fr = tzr * nz + trr * nr;
    double heat = -props.conductivity * (tg[0] * nz + tg[1] * nr);
    double rho = f.p / (medium.gasConstant(y.data()) * t), mean = 0;
    for (std::size_t k = 0; k < n; ++k) mean += y[k] / medium.species()[k].molarMass;
    std::vector<double> out(3 + n);
    double sum = 0;
    for (std::size_t k = 0; k < n; ++k) {
      auto xk = [&](double a, double b) { moles(a, b, x.data()); return x[k]; };
      auto g = gradient(xk, z, r);
      out[3 + k] = -rho * medium.species()[k].molarMass * mean * d[k] * (g[0] * nz + g[1] * nr);
      sum += out[3 + k];
    }
    medium.speciesEnthalpies(t, h.data());
    for (std::size_t k = 0; k < n; ++k) { out[3 + k] -= y[k] * sum; heat += h[k] * out[3 + k]; }
    out[0] = -fz; out[1] = -fr; out[2] = -(fz * uz + fr * ur) + heat;
    return out;
  }
};

const double gaussX[5] = {0.04691007703066800, 0.2307653449471585, 0.5, 0.7692346550528415, 0.9530899229693320};
const double gaussW[5] = {0.1184634425280945, 0.2393143352496832, 0.2844444444444444, 0.2393143352496832, 0.1184634425280945};

Definition fieldDefinition(const Field& f, const Medium& medium, const TransportFits& fits, int nz, int nr) {
  Definition d;
  d.experiment = Case::UniformDuct;
  d.nz = nz; d.nr = nr;
  d.species = medium.species();
  d.composition.assign(medium.size(), 0.0); d.composition[f.iN2] = 1;
  for (int i = 0; i <= nz; ++i) { double z = f.length * i / nz; d.contour.push_back({z, f.radius(z)}); }
  d.totalPressure = d.backPressure = f.p;
  d.transport = fits;
  d.wallTemperature = f.tWall;
  return d;
}
void setField(Flow& flow, const Field& f) {
  const auto& m = flow.mesh();
  const std::size_t n = flow.medium().size();
  std::vector<Primitive> cells(m.cells.size());
  std::vector<double> y(m.cells.size() * n);
  for (std::size_t q = 0; q < m.cells.size(); ++q) {
    const auto& c = m.cells[q];
    f.y(c.z, c.r, n, y.data() + q * n);
    double t = f.t(c.z, c.r);
    cells[q] = {f.p / (flow.medium().gasConstant(y.data() + q * n) * t), f.uz(c.z, c.r), f.ur(c.z, c.r), f.p};
  }
  flow.setInitialState(cells, y);
}

double besselJ0(double x) {
  double term = 1, sum = 1;
  for (int m = 1; m < 40; ++m) { term *= -(x * x / 4) / (m * m); sum += term; }
  return sum;
}
}  // namespace

int main() {
  const std::string mechanism = "h2o2.yaml";
  thermo::ReactionSource source(mechanism);
  auto medium = source.medium();
  const auto fits = thermo::transportFits(mechanism);
  medium.setTransport(fits);
  const std::size_t n = medium.size();
  auto index = [&](const char* name) {
    for (std::size_t k = 0; k < n; ++k) if (medium.species()[k].name == name) return k;
    throw std::runtime_error(std::string("missing species ") + name);
  };
  Field field;
  field.iH2 = index("H2"); field.iO2 = index("O2"); field.iH2O = index("H2O"); field.iOH = index("OH"); field.iN2 = index("N2");
  Exact exact(field, medium);

  std::printf("1. transport operator order on a contoured duct (h2o2.yaml, no-slip wall at %.0f K)\n", field.tWall);
  const char* names[4] = {"axial momentum", "radial momentum", "energy", "species"};
  std::vector<std::array<double, 4>> errors;
  std::vector<int> grids = {32, 64, 128, 256};
  for (int nz : grids) {
    int nr = nz / 4;
    Flow flow(fieldDefinition(field, medium, fits, nz, nr));
    setField(flow, field);
    std::vector<Conserved> rate;
    std::vector<double> speciesRate;
    flow.transportDerivative(rate, speciesRate);
    const auto& m = flow.mesh();
    std::array<double, 4> num{}, den{}, wallNum{}, wallDen{};
    for (int i = 1; i < nz - 1; ++i)
      for (int j = 0; j < nr; ++j) {
        auto q = m.index(i, j);
        const auto& c = m.cells[q];
        // Exact cell average: -(1/V) (outward face integrals) + (1/V) integral of -tau_thth / r dV.
        std::vector<double> total(3 + n, 0.0);
        auto face = [&](double z0, double r0, double z1, double r1, double nzOut, double nrOut) {
          double len = std::hypot(z1 - z0, r1 - r0);
          if (!(r0 + r1 > 0)) return;
          for (int g = 0; g < 5; ++g) {
            double z = z0 + gaussX[g] * (z1 - z0), r = r0 + gaussX[g] * (r1 - r0);
            auto fl = exact.flux(z, r, nzOut, nrOut, nullptr);
            for (std::size_t k = 0; k < fl.size(); ++k) total[k] -= gaussW[g] * 2 * pi * r * len * fl[k];
          }
        };
        double a = double(j) / nr, b = double(j + 1) / nr, z0 = i * m.dz, z1 = (i + 1) * m.dz;
        double ra = m.radius[i], rb = m.radius[i + 1], dr = rb - ra;
        face(z0, a * ra, z0, b * ra, -1, 0);
        face(z1, a * rb, z1, b * rb, 1, 0);
        // Radial faces: outward unit normals (a dr, -dz)/L below and (-b dr, dz)/L above.
        face(z0, a * ra, z1, a * rb, a * dr / std::hypot(m.dz, a * dr), -m.dz / std::hypot(m.dz, a * dr));
        face(z0, b * ra, z1, b * rb, -b * dr / std::hypot(m.dz, b * dr), m.dz / std::hypot(m.dz, b * dr));
        for (int g = 0; g < 5; ++g)
          for (int h = 0; h < 5; ++h) {
            double z = z0 + gaussX[g] * m.dz, wallR = ra + gaussX[g] * dr, sigma = a + gaussX[h] * (b - a);
            double hoop = 0;
            exact.flux(z, sigma * wallR, 1, 0, &hoop);
            total[1] -= gaussW[g] * gaussW[h] * m.dz * (b - a) * 2 * pi * wallR * hoop;
          }
        bool wallRow = j == nr - 1;
        auto& en = wallRow ? wallNum : num;
        auto& ed = wallRow ? wallDen : den;
        for (int k = 0; k < 3; ++k) {
          double e = total[k] / c.volume;
          en[k] += c.volume * std::abs(rate[q][k + 1] - e);
          ed[k] += c.volume * std::abs(e);
        }
        for (std::size_t k = 0; k < n; ++k) {
          double e = total[3 + k] / c.volume;
          en[3] += c.volume * std::abs(speciesRate[q * n + k] - e);
          ed[3] += c.volume * std::abs(e);
        }
      }
    std::array<double, 4> e{};
    std::printf("  %4dx%-3d", nz, nr);
    for (int k = 0; k < 4; ++k) { e[k] = num[k] / den[k]; std::printf("  %s %.3e", names[k], e[k]); }
    std::printf("\n           wall row:");
    for (int k = 0; k < 4; ++k) std::printf("  %.3e", wallNum[k] / wallDen[k]);
    std::printf("\n");
    errors.push_back(e);
  }
  for (int k = 0; k < 4; ++k) {
    double order = std::log2(errors[errors.size() - 2][k] / errors.back()[k]);
    char what[96];
    std::snprintf(what, sizeof what, "%s: observed order (two finest grids)", names[k]);
    check(order >= 1.8, what, order, 1.8);
  }

  std::printf("2. pipe decay, u_z = U J0(j01 r/R) exp(-nu j01^2 t/R^2), N2 at 1 kPa, R = 1 mm\n");
  const double j01 = 2.404825557695773, radius = 1e-3, pressure = 1000, temperature = 300, speed = 1;
  std::vector<double> decay;
  double momentumWorst = 0;
  for (int nr : {8, 16, 32}) {
    Definition d;
    d.experiment = Case::UniformDuct;
    d.nz = 4; d.nr = nr; d.length = 1e-3; d.inletRadius = d.exitRadius = d.throatRadius = radius;
    d.species = medium.species();
    d.composition.assign(n, 0.0); d.composition[field.iN2] = 1;
    d.totalPressure = d.backPressure = pressure; d.totalTemperature = temperature;
    d.transport = fits; d.wallTemperature = temperature;
    Flow flow(d);
    std::vector<double> diffusion(n), work;
    auto props = medium.transport(temperature, pressure, d.composition.data(), diffusion.data(), work);
    const double rho = pressure / (medium.gasConstant(d.composition.data()) * temperature), nu = props.viscosity / rho;
    const auto& m = flow.mesh();
    std::vector<Primitive> cells(m.cells.size());
    for (std::size_t q = 0; q < cells.size(); ++q) cells[q] = {rho, speed * besselJ0(j01 * m.cells[q].r / radius), 0, pressure};
    flow.setInitialState(cells);
    const double tEnd = radius * radius / (nu * j01 * j01);
    flow.advanceTo(tEnd);
    double num = 0, den = 0, factor = std::exp(-1.0);
    for (std::size_t q = 0; q < cells.size(); ++q) {
      double ex = speed * besselJ0(j01 * m.cells[q].r / radius) * factor;
      num += m.cells[q].volume * std::abs(flow.cellPrimitive(q).uz - ex);
      den += m.cells[q].volume * std::abs(ex);
    }
    auto meas = flow.measurements();
    momentumWorst = std::max(momentumWorst, std::abs(meas.momentumBalanceError));
    std::printf("  nr %-3d  t %.4e s  steps %llu  L1 error %.3e  momentum balance %.2e  energy balance %.2e\n", nr, tEnd,
                static_cast<unsigned long long>(meas.steps), num / den, meas.momentumBalanceError, meas.energyBalanceError);
    decay.push_back(num / den);
  }
  check(std::log2(decay[1] / decay[2]) >= 1.8, "pipe decay: observed order (nr 16 -> 32)", std::log2(decay[1] / decay[2]), 1.8);
  check(decay.back() < 1e-3, "pipe decay: L1 error at nr 32", decay.back(), 1e-3);
  check(momentumWorst < 1e-12, "pipe decay: momentum balance (wall shear)", momentumWorst, 1e-12);

  std::printf("3. budgets with transport: closed chamber, isothermal wall and plate, ambient exit\n");
  {
    Definition d = fieldDefinition(field, medium, fits, 32, 8);
    d.experiment = Case::Chamber;
    d.ambientTemperature = field.tWall;
    Flow flow(d);
    setField(flow, field);
    for (int s = 0; s < 200; ++s) flow.step();
    auto meas = flow.measurements();
    std::printf("  t %.3e s, wall heat flow %.4g W, wall axial force %.4g N\n", meas.time, meas.wallHeatFlow, meas.wallAxialForce);
    check(std::abs(meas.massBalanceError) < 1e-12, "mass balance", std::abs(meas.massBalanceError), 1e-12);
    check(std::abs(meas.energyBalanceError) < 1e-12, "energy balance (with wall heat)", std::abs(meas.energyBalanceError), 1e-12);
    check(std::abs(meas.momentumBalanceError) < 1e-12, "axial momentum balance (with wall shear)", std::abs(meas.momentumBalanceError), 1e-12);
  }
  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
