#pragma once
// C1's contour, shared by tests/chamber_study.cpp and the AMR budget run (tests/amr_budget.cpp).
// Throat radius 10 mm, chamber radius 25 mm, circular throat arcs of radius 20 mm on both sides
// (R = r_c / r_t = 2), 30 degree convergent cone, 20 mm cylinder-to-cone blend, 15 degree conical exit
// to area ratio 10. The cylinder length puts the throat at mid-length, so an even nz places a mesh
// station exactly on the throat.
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

namespace crucible::c1_geometry {

inline constexpr double kRt = 0.010, kRc = 0.025, kArc = 0.020, kBlend = 0.020;
inline const double kAlpha = 30 * std::numbers::pi / 180, kTheta = 15 * std::numbers::pi / 180, kEps = 10;

struct Geometry {
  std::vector<std::array<double, 2>> contour;
  double cylinderEnd = 0;  // z where the cylinder meets the blend arc
  double throat = 0;
};

// Contour with the throat at z = 0, then shifted so the injector face is z = 0.
inline Geometry geometry() {
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

}  // namespace crucible::c1_geometry
