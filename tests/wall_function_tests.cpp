// Wall functions, criterion 0: the algebra (docs/evidence/WALL_FUNCTIONS.md, criteria stated 5 October
// 2026 before any wall-function code). Nichols and Nelson's law in core/walls.cpp (crucible::wallLaw)
// against exact or analytic references; kappa 0.4 and B 5.5 here (the algebra does not depend on them).
//   (a) beta = 0, Gamma = 0: eq. 10.6 equals Spalding's eq. 10.3 and eq. 10.12 (mu_1 = mu_w) equals
//       eq. 10.11 within 1e-12 relative for u+ in [0, 40]. y+(u+) is continuous as Gamma goes from 0
//       to 1e-2: Gamma enters eq. 10.6 only through u_eq+, which must equal an independent composite
//       Gauss-Legendre quadrature of its integral within 1e-13 relative for Gamma = 0, 1e-300, 1e-100
//       and 10^(e / 100) from 1e-30 to 1e-2, beta in {-0.05, 0, 0.05} and u+ in [0.5, 40] (where
//       1 + beta s - Gamma s^2 >= 0.05 on [0, u+]); and y+ must not decrease as Gamma grows on that
//       ladder by more than its rounding, 1e-15 (1 + kappa u_eq+) relative (exp(kappa u_eq+)
//       amplifies a relative rounding of u_eq+ by kappa u_eq+).
//   (b) Inversion. For first-cell y+ in {0.5, 1, 5, 11, 30, 100, 300, 1000, 3000}, Gamma in {0, 1e-5,
//       1e-4, 1e-3} and beta in {-0.05, 0, 0.02, 0.05, 0.1}: u+ is found from eq. 10.6 at that
//       (Gamma, beta) by bisection here, a layer is built with that u_tau (cp chosen to give Gamma,
//       T_1 from eq. 10.8), and solveIsothermal must return u_tau within 1e-12 relative, satisfy
//       eq. 10.6 at its own (u+, Gamma, beta) within 1e-12, and its q_w must give T_1 back through
//       eq. 10.8 within 1e-12. A case whose eq. 10.8 profile reaches T <= 0 at u+ has no layer and is
//       listed, not solved.
//   (c) dy+/du+ from the corrected eq. 10.13 against a five-point central difference of eq. 10.6
//       (step 1e-3 of the smaller of u+ and its distance to the edge of the profile's range) within
//       1e-7 relative over the cases of (b); the printed eq. 10.13's error is reported.
#include "core/walls.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using namespace crucible;

namespace {
int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-78s %.3e (limit %.0e)  %s\n", what, value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}
// Relative difference; infinite if either value is not finite, so a NaN cannot pass.
double relative(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b)) return std::numeric_limits<double>::infinity();
  return a == b ? 0.0 : std::abs(a - b) / std::max(std::abs(a), std::abs(b));
}
// u_eq+ = integral from 0 to u of ds / sqrt(1 + beta s - Gamma s^2): 64 panels of 20-point
// Gauss-Legendre, nodes by Newton's method on P_20.
double quadrature(double u, double gamma, double beta) {
  constexpr int n = 20, panels = 64;
  static std::array<double, n> x{}, w{};
  if (w[0] == 0)
    for (int i = 0; i < n / 2; ++i) {
      double z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), pp = 0, previous;
      do {
        double p1 = 1, p2 = 0;
        for (int j = 1; j <= n; ++j) { const double p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
        pp = n * (z * p1 - p2) / (z * z - 1);
        previous = z;
        z -= p1 / pp;
      } while (std::abs(z - previous) > 1e-16);
      x[i] = -z; x[n - 1 - i] = z;
      w[i] = w[n - 1 - i] = 2 / ((1 - z * z) * pp * pp);
    }
  double sum = 0;
  for (int k = 0; k < panels; ++k) {
    const double a = u * k / panels, b = u * (k + 1) / panels, half = 0.5 * (b - a), mid = 0.5 * (a + b);
    double panel = 0;
    for (int i = 0; i < n; ++i) { const double s = mid + half * x[i]; panel += w[i] / std::sqrt(1 + beta * s - gamma * s * s); }
    sum += half * panel;
  }
  return sum;
}
}  // namespace

int main() {
  const wallLaw::Constants c{0.4, 5.5};
  std::printf("criterion 0, kappa %.2f, B %.2f\n(a) incompressible adiabatic limit\n", c.kappa, c.b);
  {
    double worstY = 0, worstEddy = 0;
    for (int n = 0; n <= 4000; ++n) {
      const double u = 0.01 * n;
      worstY = std::max(worstY, relative(wallLaw::yPlus(u, 0, 0, c), wallLaw::spalding(u, c)));
      worstEddy = std::max(worstEddy, relative(wallLaw::eddyRatio(u, 0, 0, 1, c), wallLaw::spaldingEddy(u, c)));
    }
    check(worstY < 1e-12, "eq. 10.6 against Spalding's eq. 10.3, u+ 0 to 40", worstY, 1e-12);
    check(worstEddy < 1e-12, "eq. 10.12 against eq. 10.11, u+ 0 to 40", worstEddy, 1e-12);
    double worstEq = 0, worstDrop = 0;
    int points = 0;
    for (double u : {0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0}) {
      for (double beta : {-0.05, 0.0, 0.05}) {
        std::vector<double> ladder{0.0, 1e-300, 1e-100};
        for (int e = -3000; e <= -200; ++e) ladder.push_back(std::pow(10.0, e / 100.0));
        double previous = 0;
        for (double gamma : ladder) {
          // The profile's smallest 1 + beta s - Gamma s^2 on [0, u] is at an end (the quadratic is concave).
          if (std::min(1.0, 1 + beta * u - gamma * u * u) < 0.05) break;
          const double ue = wallLaw::equivalentVelocity(u, gamma, beta);
          worstEq = std::max(worstEq, relative(ue, quadrature(u, gamma, beta)));
          const double y = wallLaw::yPlus(u, gamma, beta, c);
          const double allowed = 1e-15 * (1 + c.kappa * ue);
          if (y < previous) worstDrop = std::max(worstDrop, (previous - y) / previous / allowed);
          if (!std::isfinite(y)) worstDrop = std::numeric_limits<double>::infinity();
          previous = y;
          ++points;
        }
      }
    }
    std::printf("  %d (u+, Gamma, beta) points\n", points);
    check(worstEq < 1e-13, "continuity: u_eq+ against Gauss-Legendre quadrature, Gamma 0 to 1e-2", worstEq, 1e-13);
    check(worstDrop <= 1, "continuity: largest decrease of y+ as Gamma grows, over its rounding allowance", worstDrop, 1);
  }

  std::printf("(b) inversion and (c) derivative\n");
  double worstTau = 0, worstResidual = 0, worstT1 = 0, worstDerivative = 0, worstPrinted = 0, printedAt[3] = {};
  int solved = 0, maxIterations = 0;
  std::vector<std::array<double, 3>> skipped;
  for (double target : {0.5, 1.0, 5.0, 11.0, 30.0, 100.0, 300.0, 1000.0, 3000.0})
    for (double gamma : {0.0, 1e-5, 1e-4, 1e-3})
      for (double beta : {-0.05, 0.0, 0.02, 0.05, 0.1}) {
        // u+ at this (Gamma, beta) by bisection on eq. 10.6, inside the range where 1 + beta s - Gamma s^2 > 0.
        double lo = 0, hi = 200;
        auto feasible = [&](double s) { return 1 + beta * s - gamma * s * s > 0; };
        while (!feasible(hi)) hi *= 0.999;
        const double hi0 = hi;
        if (wallLaw::yPlus(hi, gamma, beta, c) < target) { skipped.push_back({target, gamma, beta}); continue; }
        for (int it = 0; it < 200; ++it) {
          const double mid = 0.5 * (lo + hi);
          (wallLaw::yPlus(mid, gamma, beta, c) < target ? lo : hi) = mid;
        }
        const double u = 0.5 * (lo + hi);
        // Ratio 1 + beta u+ - Gamma u+^2 = T_1 / T_w; skip where it is not positive (checked above by hi).
        const double tw = 500, rhoW = 2.0, muW = 3e-5, kW = 0.04, recovery = 0.9, y1 = 2e-4;
        const double uTau = wallLaw::yPlus(u, gamma, beta, c) * muW / (rhoW * y1);
        const double cp = gamma > 0 ? recovery * uTau * uTau / (2 * gamma * tw) : std::numeric_limits<double>::infinity();
        const double t1 = tw * (1 + beta * u - gamma * u * u);
        const wallLaw::Layer layer{u * uTau, y1, t1, tw, rhoW, muW, kW, cp, recovery};
        const auto s = wallLaw::solveIsothermal(layer, c);
        ++solved;
        maxIterations = std::max(maxIterations, s.iterations);
        worstTau = std::max(worstTau, relative(s.uTau, uTau));
        const double yp = rhoW * s.uTau * y1 / muW, up = layer.u1 / s.uTau;
        const double g = recovery * s.uTau * s.uTau / (2 * cp * tw), b = s.heat * muW / (rhoW * tw * kW * s.uTau);
        worstResidual = std::max(worstResidual, relative(wallLaw::yPlus(up, g, b, c), yp));
        worstT1 = std::max(worstT1, relative(tw * (1 + b * up - g * up * up), t1));
        // Five-point central difference of eq. 10.6 at fixed (Gamma, beta), inside the feasible range.
        const double h = 1e-3 * std::min(u, hi0 - u);
        auto y = [&](double v) { return wallLaw::yPlus(v, gamma, beta, c); };
        const double fd = (-y(u + 2 * h) + 8 * y(u + h) - 8 * y(u - h) + y(u - 2 * h)) / (12 * h);
        worstDerivative = std::max(worstDerivative, relative(wallLaw::yPlusDerivative(u, gamma, beta, c), fd));
        const double printed = relative(wallLaw::yPlusDerivativePrinted(u, gamma, beta, c), fd);
        if (printed > worstPrinted) { worstPrinted = printed; printedAt[0] = target; printedAt[1] = gamma; printedAt[2] = beta; }
      }
  std::printf("  %d cases solved, at most %d iterations; %zu without a layer (T <= 0 before eq. 10.6 reaches y+):\n",
              solved, maxIterations, skipped.size());
  for (const auto& k : skipped) std::printf("    y+ %g, Gamma %g, beta %g\n", k[0], k[1], k[2]);
  check(worstTau < 1e-12, "(b) solved u_tau against the constructed one", worstTau, 1e-12);
  check(worstResidual < 1e-12, "(b) eq. 10.6 at the solved u+, Gamma, beta", worstResidual, 1e-12);
  check(worstT1 < 1e-12, "(b) T_1 from the solved q_w through eq. 10.8", worstT1, 1e-12);
  check(worstDerivative < 1e-7, "(c) corrected eq. 10.13 against the five-point difference", worstDerivative, 1e-7);
  std::printf("  reported: printed eq. 10.13 against the same difference, largest %.3e at y+ %g, Gamma %g, beta %g\n",
              worstPrinted, printedAt[0], printedAt[1], printedAt[2]);
  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
