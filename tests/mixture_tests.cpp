// Item-4 verification: the core's thermally perfect species mixture.
// Tolerances stated before the run:
//   1. Core NASA-7 thermo against Cantera on the same mechanism (gri30): gas constant, internal
//      energy, cv and frozen sound speed within 1e-12 relative; temperature inversion of
//      Cantera's internal energy within 1e-9 K. Five compositions, T = 300..3500 K.
//   2. Two-gas shock tube (air-like gamma 1.4 | helium-like gamma 5/3, constant cp each) against
//      the exact two-gamma Riemann solution: density L1(400 cells) <= 0.5 L1(100 cells); star
//      pressure midway between contact and shock within 1% at 400 cells; species masses sum to
//      the total mass within 1e-12 relative.
//   3. Same two gases advected at uniform pressure and velocity: measured and reported only
//      (spurious pressure at a variable-gamma contact is a known property of fully conservative
//      schemes; its size decides whether a correction is needed).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "adapters/reaction.hpp"
#include "core/flow.hpp"

using namespace crucible;

namespace {
int failures = 0;
void check(bool ok, const char* what, double value, double tol) {
  std::printf("  %-48s %12.4e  (tol %.1e) %s\n", what, value, tol, ok ? "ok" : "FAIL");
  failures += !ok;
}

Species constantCp(const char* name, double gamma, double molarMass) {
  Species s{name, molarMass, 1000, {}, {}};
  s.low[0] = s.high[0] = gamma / (gamma - 1);
  return s;
}

// Exact Riemann solution for two calorically perfect gases (Toro ch. 4, gamma per side).
struct Side { double rho, u, p, g; };
struct Exact {
  Side l, r;
  double pStar{}, uStar{};
  double f(const Side& s, double p) const {
    double a = std::sqrt(s.g * s.p / s.rho);
    if (p > s.p) {
      double A = 2 / ((s.g + 1) * s.rho), B = (s.g - 1) / (s.g + 1) * s.p;
      return (p - s.p) * std::sqrt(A / (p + B));
    }
    return 2 * a / (s.g - 1) * (std::pow(p / s.p, (s.g - 1) / (2 * s.g)) - 1);
  }
  Exact(Side left, Side right) : l(left), r(right) {
    double lo = 1e-12 * std::min(l.p, r.p), hi = 100 * std::max(l.p, r.p);
    for (int n = 0; n < 200; ++n) {
      double m = 0.5 * (lo + hi);
      if (f(l, m) + f(r, m) + r.u - l.u > 0) hi = m; else lo = m;
    }
    pStar = 0.5 * (lo + hi);
    uStar = 0.5 * (l.u + r.u) + 0.5 * (f(r, pStar) - f(l, pStar));
  }
  // Density at similarity coordinate s = x/t.
  double density(double s) const {
    bool left = s < uStar;
    const Side& k = left ? l : r;
    double sign = left ? -1 : 1, a = std::sqrt(k.g * k.p / k.rho), g = k.g;
    if (pStar > k.p) {
      double shock = k.u + sign * a * std::sqrt((g + 1) / (2 * g) * pStar / k.p + (g - 1) / (2 * g));
      double behind = k.rho * (pStar / k.p + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * pStar / k.p + 1);
      return (left ? s < shock : s > shock) ? k.rho : behind;
    }
    double aStar = a * std::pow(pStar / k.p, (g - 1) / (2 * g));
    double head = k.u + sign * a, tail = uStar + sign * aStar;
    if (left ? s < head : s > head) return k.rho;
    if (left ? s > tail : s < tail) return k.rho * std::pow(pStar / k.p, 1 / g);
    double c = 2 / (g + 1) - sign * (g - 1) / ((g + 1) * a) * (k.u - s);
    return k.rho * std::pow(c, 2 / (g - 1));
  }
  double shockSpeed() const {  // right-moving shock (right side compressed)
    double a = std::sqrt(r.g * r.p / r.rho), g = r.g;
    return r.u + a * std::sqrt((g + 1) / (2 * g) * pStar / r.p + (g - 1) / (2 * g));
  }
};

Definition tube(int n) {
  Definition d;
  d.nz = n; d.nr = 1; d.length = 1; d.inletRadius = d.throatRadius = d.exitRadius = 0.05;
  d.experiment = Case::ShockTube;
  d.species = {constantCp("air", 1.4, 28.96), constantCp("helium", 5.0 / 3.0, 4.0026)};
  d.composition = {1, 0};
  return d;
}

void initialise(Flow& f, int n, const Side& l, const Side& r) {
  std::vector<Primitive> cells;
  std::vector<double> y;
  for (int i = 0; i < n; ++i) {
    bool left = (i + 0.5) / n < 0.5;
    const Side& s = left ? l : r;
    cells.push_back({s.rho, s.u, 0, s.p});
    y.insert(y.end(), {left ? 1.0 : 0.0, left ? 0.0 : 1.0});
  }
  f.setInitialState(cells, y);
}
}  // namespace

int main() {
  std::printf("1. Core NASA-7 mixture thermo against Cantera (gri30)\n");
  {
    thermo::ReactionSource cantera("gri30.yaml");
    const Medium medium = cantera.medium();
    double worstProp = 0, worstT = 0;
    for (const char* moles : {"CH4:1, O2:2", "O2:1", "CH4:1", "H2O:2, CO2:1, CO:0.3, OH:0.2, H:0.1, O:0.1",
                              "N2:0.79, O2:0.21"}) {
      auto y = cantera.massFractions(moles);
      for (double t : {300.0, 800.0, 999.0, 1001.0, 1500.0, 2500.0, 3500.0}) {
        double rho = 1.7;
        double p = cantera.pressure(t, rho, y.data());
        double e = cantera.internalEnergy(t, rho, y.data()), cv = cantera.cv(t, rho, y.data());
        double r = p / (rho * t);
        double a = std::sqrt((cv + r) / cv * r * t);
        worstProp = std::max({worstProp, std::abs(medium.gasConstant(y.data()) / r - 1),
                              std::abs(medium.internalEnergy(t, y.data()) - e) / (cv * t),
                              std::abs(medium.cv(t, y.data()) / cv - 1),
                              std::abs(medium.soundSpeed(t, y.data()) / a - 1)});
        worstT = std::max(worstT, std::abs(medium.temperature(e, y.data(), 1000) - t));
      }
    }
    check(worstProp <= 1e-12, "worst relative difference R, e/(cv T), cv, a", worstProp, 1e-12);
    check(worstT <= 1e-9, "worst temperature inversion error [K]", worstT, 1e-9);
  }

  std::printf("2. Two-gas shock tube against the exact two-gamma solution\n");
  {
    const double gl = 1.4, gr = 5.0 / 3.0;
    Side l{1.0, 0, 1e5, gl}, r{0.125, 0, 1e4, gr};
    Exact exact(l, r);
    const double t = 4e-4;
    std::printf("  exact p* %.6g Pa, u* %.6g m/s, shock %.6g m/s\n", exact.pStar, exact.uStar, exact.shockSpeed());
    double l1[2], pMid = 0, massClosure = 0, oscillation = 0;
    int idx = 0;
    for (int n : {100, 400}) {
      Flow f(tube(n));
      initialise(f, n, l, r);
      f.advanceTo(t);
      double err = 0, total = 0, partial = 0;
      for (int i = 0; i < n; ++i) {
        double x = (i + 0.5) / n - 0.5;
        err += std::abs(f.state()[i][0] - exact.density(x / t)) / n;
        total += f.state()[i][0];
        partial += f.partialDensities()[2 * i] + f.partialDensities()[2 * i + 1];
      }
      l1[idx++] = err;
      massClosure = std::abs(partial / total - 1);
      if (n == 400) {
        double contact = 0.5 + exact.uStar * t, shock = 0.5 + exact.shockSpeed() * t;
        pMid = f.cellPrimitive(static_cast<std::size_t>(0.5 * (contact + shock) * n)).p;
        // Spread of pressure across the whole star region, 4 cells clear of the rarefaction tail and the shock.
        double rhoStarL = l.rho * std::pow(exact.pStar / l.p, 1 / gl);
        double tailSpeed = exact.uStar - std::sqrt(gl * exact.pStar / rhoStarL);
        int i0 = static_cast<int>((0.5 + tailSpeed * t) * n) + 4, i1 = static_cast<int>(shock * n) - 4;
        for (int i = i0; i <= i1; ++i)
          oscillation = std::max(oscillation, std::abs(f.cellPrimitive(i).p / exact.pStar - 1));
      }
    }
    std::printf("  density L1: 100=%.6g 400=%.6g\n", l1[0], l1[1]);
    check(l1[1] <= 0.5 * l1[0], "L1(400)/L1(100)", l1[1] / l1[0], 0.5);
    check(std::abs(pMid / exact.pStar - 1) <= 0.01, "star pressure, contact-shock midpoint", std::abs(pMid / exact.pStar - 1), 0.01);
    check(massClosure <= 1e-12, "sum of species masses vs total mass", massClosure, 1e-12);
    std::printf("  max |p/p*-1| across the star region (400 cells): %.4e (reported)\n", oscillation);
  }

  std::printf("3. Variable-gamma contact advected at uniform p, u (reported)\n");
  for (int n : {100, 400}) {
    Flow f(tube(n));
    const double u = 100, p = 1e5;
    Side l{p / (287.1 * 300), u, p, 1.4}, r{p / (2077.3 * 300), u, p, 5.0 / 3.0};
    initialise(f, n, l, r);
    f.advanceTo(0.25 / u);
    double dp = 0, du = 0;
    for (int i = 0; i < n; ++i) {
      auto w = f.cellPrimitive(i);
      dp = std::max(dp, std::abs(w.p / p - 1));
      du = std::max(du, std::abs(w.uz / u - 1));
    }
    std::printf("  %d cells, contact moved 0.25 m: max |p/p0-1| %.4e, max |u/u0-1| %.4e\n", n, dp, du);
  }

  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
