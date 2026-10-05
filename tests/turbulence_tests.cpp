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
