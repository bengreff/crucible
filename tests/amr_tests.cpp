// AMR, criterion 0 (docs/evidence/AMR.md, criteria stated 13 October 2026 before any AMR code): the
// explicit-geometry Mesh and its coarsening constructor.
//   (c) For a cylinder and for a cone, at uniform and at stretched rings, Mesh::coarsened of a 2n x 2m
//       mesh against the directly built n x m mesh: every array within 1e-14, each array's difference
//       taken relative to its largest magnitude (the cells' five fields, the axial faces' areas and ring
//       middles, the radial faces' two area components, middles, reference radii and r^2 moments, the
//       stations, the fractions, the column integrals of R^2 and R^3, the wall radius at five points per
//       column, and the side-wall distance of every centroid). The one stated exception: the radial
//       face's r^2 moment f^2 (r0^2 + r1^2) / 2 is a trapezoid mean, so the coarse value (the mean of the
//       two finer faces') differs on a cone from the direct one by exactly -f^2 (r1 - r0)^2 / 8; that
//       difference is checked within 1e-14 of the array's largest magnitude.
//   (d) The RL10 development contour (as tests/chamber_study.cpp rl10dev: TM-107318 Table E1 by PCHIP,
//       the face placed for the 150-column base), built directly at level 3 (1200 x 128, stretching
//       1.7347) and coarsened to levels 2, 1 and 0: every cell's outward face area vectors sum to
//       (0, its radialPressureMeasure) within 1e-14 of the largest component among its four faces.
//       Printed beside, not judged: the same sums relative to the pressure measure and to pi R^2 at the
//       cell's wider station (the size of the terms the face formulas difference), and level 0 coarsened
//       against level 0 built directly (the chords' error the polyline removes).
// The denominators above were chosen with the code, before its first run.
#include "core/contour.hpp"
#include "core/flow.hpp"
#include "core/walls.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace crucible;

namespace {
int failures = 0;
void check(bool ok, const std::string& what, double value, double limit) {
  std::printf("  %-74s %.3e (limit %.0e)  %s\n", what.c_str(), value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}
// Largest |a - b - shift| over the arrays, over the largest |b|; infinite if any value is not finite.
template <class Get, class Shift>
double arrayError(std::size_t n, Get get, Shift shift) {
  double worst = 0, scale = 0;
  for (std::size_t k = 0; k < n; ++k) {
    const auto [a, b] = get(k);
    if (!std::isfinite(a) || !std::isfinite(b)) return std::numeric_limits<double>::infinity();
    worst = std::max(worst, std::abs(a - b - shift(k)));
    scale = std::max(scale, std::abs(b));
  }
  return scale > 0 ? worst / scale : worst;
}
template <class Get>
double arrayError(std::size_t n, Get get) {
  return arrayError(n, get, [](std::size_t) { return 0.0; });
}

void compare(const std::string& name, const Definition& fine, const Definition& coarse) {
  const Mesh m = Mesh::coarsened(Mesh(fine)), d(coarse);
  std::printf("(c) %s: %dx%d coarsened against %dx%d built directly\n", name.c_str(), fine.nz, fine.nr, coarse.nz,
              coarse.nr);
  if (m.nz != d.nz || m.nr != d.nr || std::abs(m.dz - d.dz) > 1e-14 * d.dz) {
    check(false, "shape and dz", 1, 0);
    return;
  }
  const double limit = 1e-14;
  const std::size_t cells = d.cells.size(), axial = static_cast<std::size_t>(d.nz + 1) * d.nr,
                    radial = static_cast<std::size_t>(d.nz) * (d.nr + 1);
  auto axialPair = [&](auto f) {
    return [&, f](std::size_t k) {
      const int face = static_cast<int>(k / d.nr), j = static_cast<int>(k % d.nr);
      return std::array<double, 2>{f(m, face, j), f(d, face, j)};
    };
  };
  auto radialPair = [&](auto f) {
    return [&, f](std::size_t k) {
      const int i = static_cast<int>(k / (d.nr + 1)), j = static_cast<int>(k % (d.nr + 1));
      return std::array<double, 2>{f(m, i, j), f(d, i, j)};
    };
  };
  const char* fields[] = {"volume", "z", "r", "radialPressureMeasure", "radialSecondMoment"};
  for (int field = 0; field < 5; ++field) {
    auto get = [&](std::size_t q) {
      auto pick = [field](const Cell& c) {
        return field == 0 ? c.volume : field == 1 ? c.z : field == 2 ? c.r : field == 3 ? c.radialPressureMeasure
                                                                                       : c.radialSecondMoment;
      };
      return std::array<double, 2>{pick(m.cells[q]), pick(d.cells[q])};
    };
    const double e = arrayError(cells, get);
    check(e < limit, std::string("cell ") + fields[field], e, limit);
  }
  double e = arrayError(axial, axialPair([](const Mesh& x, int f, int j) { return x.axialArea(f, j); }));
  check(e < limit, "axial face area", e, limit);
  e = arrayError(axial, axialPair([](const Mesh& x, int f, int j) { return x.ringMiddle(f, j); }));
  check(e < limit, "axial face ring middle", e, limit);
  e = arrayError(radial, radialPair([](const Mesh& x, int i, int j) { return x.radialAreaVector(i, j)[0]; }));
  check(e < limit, "radial face area vector, z component", e, limit);
  e = arrayError(radial, radialPair([](const Mesh& x, int i, int j) { return x.radialAreaVector(i, j)[1]; }));
  check(e < limit, "radial face area vector, r component", e, limit);
  e = arrayError(radial, radialPair([](const Mesh& x, int i, int j) { return x.radialFaceMiddle(i, j); }));
  check(e < limit, "radial face middle", e, limit);
  e = arrayError(radial, radialPair([](const Mesh& x, int i, int j) { return x.radialFaceRadius(i, j); }));
  check(e < limit, "radial face reference radius r_f", e, limit);
  e = arrayError(radial, radialPair([](const Mesh& x, int i, int j) { return x.radialFaceSecondMoment(i, j); }),
                 [&](std::size_t k) {
                   const int i = static_cast<int>(k / (d.nr + 1)), j = static_cast<int>(k % (d.nr + 1));
                   const double f = d.fraction[j];
                   return -f * f * (d.radius[i + 1] - d.radius[i]) * (d.radius[i + 1] - d.radius[i]) / 8;
                 });
  check(e < limit, "radial face r^2 moment, less -f^2 (r1 - r0)^2 / 8", e, limit);
  e = arrayError(d.radius.size(), [&](std::size_t k) { return std::array<double, 2>{m.radius[k], d.radius[k]}; });
  check(e < limit, "station radius", e, limit);
  e = arrayError(d.fraction.size(),
                 [&](std::size_t k) { return std::array<double, 2>{m.fraction[k], d.fraction[k]}; });
  check(e < limit, "ring fraction", e, limit);
  e = arrayError(static_cast<std::size_t>(d.nz),
                 [&](std::size_t k) { return std::array<double, 2>{m.columnSquares()[k], d.columnSquares()[k]}; });
  check(e < limit, "column integral of R^2", e, limit);
  e = arrayError(static_cast<std::size_t>(d.nz),
                 [&](std::size_t k) { return std::array<double, 2>{m.columnCubes()[k], d.columnCubes()[k]}; });
  check(e < limit, "column integral of R^3", e, limit);
  e = arrayError(static_cast<std::size_t>(d.nz) * 5, [&](std::size_t k) {
    const int i = static_cast<int>(k / 5);
    const double z = (i + 0.25 * static_cast<double>(k % 5)) * d.dz;
    return std::array<double, 2>{m.wallRadiusAt(i, z), d.wallRadiusAt(i, z)};
  });
  check(e < limit, "wall radius at five points per column", e, limit);
  const auto wm = wallDistance(m, true, {}), wd = wallDistance(d, true, {});
  e = arrayError(cells, [&](std::size_t q) { return std::array<double, 2>{wm[q], wd[q]}; });
  check(e < limit, "side-wall distance of each centroid", e, limit);
}

// Outward area vectors of cell (i, j): axial faces normal to z, radial faces as Mesh::radialAreaVector
// (pointing to larger j). Largest closure error relative to the largest component among the four faces,
// to the pressure measure, and to pi R^2 at the wider station.
std::array<double, 3> closure(const Mesh& m) {
  std::array<double, 3> worst{0, 0, 0};
  for (int i = 0; i < m.nz; ++i)
    for (int j = 0; j < m.nr; ++j) {
      const auto lo = m.radialAreaVector(i, j), hi = m.radialAreaVector(i, j + 1);
      const double left = m.axialArea(i, j), right = m.axialArea(i + 1, j);
      const double sz = right - left + hi[0] - lo[0], sr = hi[1] - lo[1];
      const double measure = m.cells[m.index(i, j)].radialPressureMeasure;
      const double scale = std::max({left, right, std::abs(lo[0]), std::abs(lo[1]), std::abs(hi[0]), std::abs(hi[1])});
      const double e = std::max(std::abs(sz), std::abs(sr - measure));
      worst[0] = std::max(worst[0], std::isfinite(e) ? e / scale : std::numeric_limits<double>::infinity());
      worst[1] = std::max(worst[1], e / measure);
      worst[2] = std::max(worst[2], e / (M_PI * std::max(m.radius[i], m.radius[i + 1]) * std::max(m.radius[i], m.radius[i + 1])));
    }
  return worst;
}

// The RL10 development wall as tests/chamber_study.cpp rl10dev builds it, with the injector face placed
// for a base of `baseColumns` and the mesh at nz x nr.
Definition rl10Wall(int baseColumns, int nz, int nr, double stretching) {
  constexpr double inch = 0.0254, exitStation = 46.0, nominalFace = -13.1, rt = 0.062738;
  std::vector<std::array<double, 2>> nodes;
  const std::string path = std::string(CRUCIBLE_SOURCE_DIR) + "/docs/validation/rl10a33a/tm107318_table_E1_si.csv";
  FILE* f = std::fopen(path.c_str(), "r");
  if (!f) throw std::runtime_error("cannot read " + path);
  char line[512];
  bool header = true;
  while (std::fgets(line, sizeof line, f)) {
    if (line[0] == '#') continue;
    if (header) { header = false; continue; }
    double zIn, zM, ratio;
    if (std::sscanf(line, "%lf,%lf,%lf", &zIn, &zM, &ratio) == 3) nodes.push_back({zM, rt * std::sqrt(ratio)});
  }
  std::fclose(f);
  if (nodes.size() != 16) throw std::runtime_error("Table E1: expected 16 nodes in " + path);
  const int k = static_cast<int>(std::lround(baseColumns * -nominalFace / (exitStation - nominalFace)));
  const double face = -exitStation * k / (baseColumns - k);
  nodes.push_back({0.0, rt});
  nodes.push_back({face * inch, rt * 5.13 / 2.47});
  nodes.push_back({exitStation * inch, rt * std::sqrt(61.0)});
  std::sort(nodes.begin(), nodes.end());
  Definition d;
  d.experiment = Case::UniformDuct;
  for (const auto& p : Pchip(nodes).sample(200)) d.contour.push_back({p[0] - face * inch, p[1]});
  d.nz = nz;
  d.nr = nr;
  d.radialStretching = stretching;
  return d;
}
}  // namespace

int main() {
  std::printf("AMR criterion 0, step A\n");
  for (double stretching : {0.0, 1.7347}) {
    Definition cylinder;
    cylinder.experiment = Case::UniformDuct;
    cylinder.length = 0.1;
    cylinder.inletRadius = 0.03;
    cylinder.radialStretching = stretching;
    Definition cone = cylinder;
    cone.contour = {{0.0, 0.05}, {0.1, 0.02}};
    for (auto* d : {&cylinder, &cone}) {
      Definition fine = *d, coarse = *d;
      fine.nz = 32; fine.nr = 16; coarse.nz = 16; coarse.nr = 8;
      compare(std::string(d == &cylinder ? "cylinder" : "cone") + (stretching > 0 ? ", stretched rings" : ", uniform rings"),
              fine, coarse);
    }
  }
  std::printf("(d) RL10 development contour, base 150 x 16, stretching 1.7347\n");
  Mesh level = Mesh(rl10Wall(150, 1200, 128, 1.7347));
  for (int l = 3; l >= 0; --l) {
    const auto e = closure(level);
    check(e[0] < 1e-14, "level " + std::to_string(l) + " (" + std::to_string(level.nz) + " x " +
                            std::to_string(level.nr) + ") closure, relative to the faces",
          e[0], 1e-14);
    std::printf("      printed, not judged: relative to the pressure measure %.3e, to pi R^2 %.3e\n", e[1], e[2]);
    if (l > 0) level = Mesh::coarsened(level);
  }
  {
    const Mesh direct(rl10Wall(150, 150, 16, 1.7347));
    double volume = 0, throat = 0;
    int throatColumn = 0;
    for (int i = 1; i < direct.nz; ++i)
      if (direct.radius[i] < direct.radius[throatColumn]) throatColumn = i;
    for (std::size_t q = 0; q < direct.cells.size(); ++q)
      volume = std::max(volume, std::abs(level.cells[q].volume / direct.cells[q].volume - 1));
    {
      double a = 0, b = 0;
      for (int j = 0; j < direct.nr; ++j) { a += level.axialArea(throatColumn, j); b += direct.axialArea(throatColumn, j); }
      throat = a / b - 1;
    }
    double side = 0;
    for (int i = 0; i < direct.nz; ++i) {
      const auto a = level.radialAreaVector(i, direct.nr), b = direct.radialAreaVector(i, direct.nr);
      side = std::max(side, std::abs(std::hypot(a[0], a[1]) / std::hypot(b[0], b[1]) - 1));
    }
    std::printf("      level 0 coarsened against built directly (printed, not judged): cell volume %.3e largest, "
                "throat station area %+.3e, side-wall face area %.3e largest\n", volume, throat, side);
  }
  std::printf("%s (%d failure%s)\n", failures ? "FAILED" : "PASSED", failures, failures == 1 ? "" : "s");
  return failures ? 1 : 0;
}
