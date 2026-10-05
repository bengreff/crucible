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
//      axial body force 2.95e5 N/m^3 (Re_tau about 182), rings clustered with b = 2, nz 2, nr 16, 32
//      and 64 (wall y+ about 0.9, 0.4, 0.2, derived). Each run starts from the reference solution on
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
#include "core/flow.hpp"
#include "core/walls.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
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
}  // namespace

int main() {
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
  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
