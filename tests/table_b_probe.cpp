// Table B probe (6 October 2026): the design's two derived claims, tested on the engine's own
// CVODES reactor (h2o2.yaml, H2/O2 at O/F 6 by mass, 3 MPa chamber). Exploratory research, not a
// criterion (docs/evidence/TABLE_B.md).
//   (1) c = Y_H2O moves monotonically from the reactants to equilibrium in a constant-(e, rho)
//       reactor.
//   (2) Partial equilibrium: a state that has left the manifold (the chamber's equilibrium,
//       frozen, then expanded isentropically or cooled at constant density) relaxes onto the lower
//       branch of its new (e, rho) node. At equal c the probe reports the differences in T, in the
//       mass fractions and in lambda = (dc/dt) / (c_eq - c). At the start this is the snap error.
//
// Usage: crucible_table_b_probe
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

#include "adapters/reaction.hpp"
#include "cantera/core.h"

using namespace crucible::thermo;

namespace {

struct Trajectory {
  std::vector<double> t, c, temperature, lambda;
  std::vector<std::vector<double>> y;
  double ceq = 0;
  double largestTurnBack = 0;  // largest fall in c against the approach, relative to |c_eq - c_0|
};

double bisect(const std::function<double(double)>& g, double lo, double hi) {
  double glo = g(lo);
  for (int i = 0; i < 200 && hi - lo > 1e-10 * hi; ++i) {
    double mid = 0.5 * (lo + hi), gm = g(mid);
    if ((gm < 0) == (glo < 0)) { lo = mid; glo = gm; } else { hi = mid; }
  }
  return 0.5 * (lo + hi);
}

}  // namespace

int main() {
  ReactionSource source("h2o2.yaml");
  const std::size_t n = source.nSpecies();
  std::size_t iw = 0;
  for (std::size_t k = 0; k < n; ++k) if (source.speciesName(k) == "H2O") iw = k;
  ReactionStep step(source, 1e-10, 1e-16);

  auto yH2 = source.massFractions("H2:1"), yO2 = source.massFractions("O2:1");
  std::vector<double> reactants(n);
  for (std::size_t k = 0; k < n; ++k) reactants[k] = yH2[k] / 7 + yO2[k] * 6 / 7;
  // Complete-combustion products of the same elements (fuel rich): H2O and the excess H2.
  std::vector<double> complete(n, 0.0);
  for (std::size_t k = 0; k < n; ++k) if (source.speciesName(k) == "H2O") complete[k] = (6.0 / 7) * 18.01528 / 15.999;
  for (std::size_t k = 0; k < n; ++k) if (source.speciesName(k) == "H2") complete[k] = 1 - (6.0 / 7) * 18.01528 / 15.999;

  auto energy = [&](double t, double rho, const std::vector<double>& y) { return source.internalEnergy(t, rho, y.data()); };
  auto temperatureAt = [&](double e, double rho, const std::vector<double>& y) {
    return bisect([&](double t) { return energy(t, rho, y) - e; }, 100, 6000);
  };
  auto rate = [&](double rho, const std::vector<double>& z) {
    std::vector<double> d(n + 1);
    source.rates(rho, z.data(), d.data());
    return d[1 + iw];
  };

  auto integrate = [&](double rho, std::vector<double> z) {
    Trajectory tr;
    std::vector<double> zeq = z;
    source.equilibrateUV(rho, zeq.data());
    tr.ceq = zeq[1 + iw];
    const double c0 = z[1 + iw], sign = tr.ceq >= c0 ? 1.0 : -1.0, span = std::fabs(tr.ceq - c0);
    double best = c0, t = 0;
    auto record = [&](double time) {
      double c = z[1 + iw];
      tr.t.push_back(time); tr.c.push_back(c); tr.temperature.push_back(z[0]);
      tr.lambda.push_back(rate(rho, z) / (tr.ceq - c));
      tr.y.emplace_back(z.begin() + 1, z.end());
      if (sign * (c - best) > 0) best = c;
      tr.largestTurnBack = std::max(tr.largestTurnBack, sign * (best - c) / span);
    };
    record(0);
    // Geometric output times, 40 per decade from 1e-10 s to 1 s, one CVODES run per interval.
    for (int i = 0; i <= 400; ++i) {
      double next = 1e-10 * std::pow(10.0, i / 40.0);
      step.advance(rho, z.data(), next - t);
      t = next;
      record(t);
      if (std::fabs(z[1 + iw] - tr.ceq) < 1e-7 * tr.ceq && i > 40) break;
    }
    return tr;
  };

  // The lower branch's start at a node, as the design states it.
  auto lowerStart = [&](double e, double rho, bool& fromReactants) {
    std::vector<double> z(n + 1);
    if (energy(200, rho, reactants) <= e) {
      fromReactants = true;
      z[0] = temperatureAt(e, rho, reactants);
      std::copy(reactants.begin(), reactants.end(), z.begin() + 1);
      return z;
    }
    fromReactants = false;
    auto frozenAt = [&](double tStar) {
      std::vector<double> w(n + 1);
      w[0] = tStar;
      std::copy(reactants.begin(), reactants.end(), w.begin() + 1);
      source.equilibrateTV(rho, w.data());
      return std::vector<double>(w.begin() + 1, w.end());
    };
    double tStar = bisect([&](double ts) { return energy(200, rho, frozenAt(ts)) - e; }, 300, 6000);
    auto y = frozenAt(tStar * (1 - 1e-9));
    z[0] = temperatureAt(e, rho, y);
    std::copy(y.begin(), y.end(), z.begin() + 1);
    return z;
  };

  auto atC = [&](const Trajectory& tr, double c, double& t, std::vector<double>& y, double& lambda) {
    for (std::size_t i = 1; i < tr.c.size(); ++i) {
      double a = tr.c[i - 1], b = tr.c[i];
      if ((c - a) * (c - b) <= 0 && a != b) {
        double w = (c - a) / (b - a);
        t = tr.temperature[i - 1] + w * (tr.temperature[i] - tr.temperature[i - 1]);
        lambda = tr.lambda[i - 1] + w * (tr.lambda[i] - tr.lambda[i - 1]);
        y.resize(n);
        for (std::size_t k = 0; k < n; ++k) y[k] = tr.y[i - 1][k] + w * (tr.y[i][k] - tr.y[i - 1][k]);
        return true;
      }
    }
    return false;
  };

  std::vector<double> wk(n);
  {
    auto s0 = Cantera::newSolution("h2o2.yaml", "", "none");
    s0->thermo()->getMolecularWeights(wk.data());
  }
  auto molesOf = [&](const std::vector<double>& yy) {
    double m = 0;
    for (std::size_t k = 0; k < n; ++k) m += yy[k] / wk[k];
    return m;
  };
  std::printf("Table B probe, h2o2.yaml, O/F 6, CVODES rtol 1e-10\n\n(1) Ignition from the reactants at 3 MPa\n");
  std::printf("%8s %12s %12s %16s %10s %26s %16s\n", "T0 K", "c_eq", "T_eq K", "largest turn-back", "t_99 s",
              "moles: largest rise / fall", "c where N falls");
  for (double t0 : {900.0, 1000.0, 1200.0, 1500.0}) {
    std::vector<double> z(n + 1);
    z[0] = t0;
    std::copy(reactants.begin(), reactants.end(), z.begin() + 1);
    double rho = source.density(t0, 3e6, reactants.data());
    auto tr = integrate(rho, z);
    double t99 = 0;
    for (std::size_t i = 0; i < tr.c.size(); ++i) if (tr.c[i] >= 0.99 * tr.ceq) { t99 = tr.t[i]; break; }
    // Moles per kilogram along the branch: the largest rise against its total fall, and the c by
    // which N has fallen 1% of its total.
    double n0 = molesOf(tr.y.front()), n1 = molesOf(tr.y.back()), lowest = n0, rise = 0, cFall = -1;
    for (std::size_t i = 0; i < tr.y.size(); ++i) {
      double m = molesOf(tr.y[i]);
      rise = std::max(rise, m - lowest);
      lowest = std::min(lowest, m);
      if (cFall < 0 && n0 - m > 0.01 * (n0 - n1)) cFall = tr.c[i];
    }
    std::printf("%8.0f %12.6f %12.1f %16.3e %10.3e %26.3e %16.3e\n", t0, tr.ceq, tr.temperature.back(), tr.largestTurnBack, t99,
                rise / (n0 - n1), cFall / tr.ceq);
  }

  // The chamber: constant-(h, p) equilibrium of the reactants at 300 K, 3 MPa.
  auto sol = Cantera::newSolution("h2o2.yaml", "", "none");
  auto& gas = *sol->thermo();
  gas.setMassFractions(reactants.data());
  gas.setState_TP(300, 3e6);
  gas.equilibrate("HP");
  const double tc = gas.temperature(), rhoc = gas.density(), sc = gas.entropy_mass();
  std::vector<double> chamber(n);
  gas.getMassFractions(chamber.data());
  std::printf("\nchamber equilibrium: T %.1f K, rho %.4f kg/m3, c %.6f\n", tc, rhoc, chamber[iw]);

  // Departure from partial equilibrium: the largest |ln(forward/reverse rate of progress)| over
  // the three bimolecular chain reactions, which is zero when they are equilibrated.
  auto& kin = *sol->kinetics();
  std::vector<std::size_t> chain;
  for (std::size_t i = 0; i < kin.nReactions(); ++i) {
    auto eq = kin.reaction(i)->equation();
    if (eq == "H + O2 <=> O + OH" || eq == "H2 + O <=> H + OH" || eq == "H2 + OH <=> H + H2O") chain.push_back(i);
  }
  std::printf("chain reactions found: %zu\n", chain.size());
  auto departure = [&](double t, double rho, const std::vector<double>& y) {
    gas.setMassFractions_NoNorm(y.data());
    gas.setState_TD(t, rho);
    std::vector<double> f(kin.nReactions()), r(kin.nReactions());
    kin.getFwdRatesOfProgress(f.data());
    kin.getRevRatesOfProgress(r.data());
    double d = 0;
    for (auto i : chain) d = std::max(d, std::fabs(std::log(f[i] / r[i])));
    return d;
  };

  std::printf("\n(2) Off-manifold starts: the chamber's composition, frozen\n");
  std::printf("%-22s %8s %8s %9s %6s | %-34s | %-34s | %-34s | %-34s\n", "case", "T K", "c0/ceq", "branch", "turn",
              "phi 0 (snap): dT K, max|dY|, lam", "phi 0.1", "phi 0.5", "phi 0.9");
  struct Case { const char* name; double rhoRatio; double tCooled; };
  std::vector<Case> cases = {{"expanded rho/rho_c 0.3", 0.3, 0}, {"expanded rho/rho_c 0.1", 0.1, 0},
                             {"expanded rho/rho_c 0.03", 0.03, 0}, {"cooled to 2500 K", 1.0, 2500},
                             {"cooled to 2000 K", 1.0, 2000}, {"cooled to 1500 K", 1.0, 1500}};
  for (const auto& cs : cases) {
    double rho = rhoc * cs.rhoRatio, t;
    if (cs.tCooled > 0) {
      t = cs.tCooled;
    } else {
      t = bisect([&](double tt) {
        gas.setMassFractions_NoNorm(chamber.data());
        gas.setState_TD(tt, rho);
        return gas.entropy_mass() - sc;
      }, 200, tc);
    }
    std::vector<double> z(n + 1);
    z[0] = t;
    std::copy(chamber.begin(), chamber.end(), z.begin() + 1);
    const double e = energy(t, rho, chamber);
    auto off = integrate(rho, z);
    bool fromReactants = false;
    auto branch = integrate(rho, lowerStart(e, rho, fromReactants));
    std::printf("%-22s %8.1f %8.4f %9s %6.0e |", cs.name, t, chamber[iw] / off.ceq, fromReactants ? "reactants" : "frozen eq",
                branch.largestTurnBack);
    const double c0 = chamber[iw];
    std::vector<std::pair<double, double>> pe;
    for (double phi : {0.0, 0.1, 0.5, 0.9}) {
      double c = c0 + phi * (off.ceq - c0), tb, to, lb, lo;
      std::vector<double> yb, yo;
      if (!atC(branch, c, tb, yb, lb) || !atC(off, c, to, yo, lo)) {
        std::printf(" %-34s |", "branch does not reach this c");
        continue;
      }
      double dy = 0;
      for (std::size_t k = 0; k < n; ++k) dy = std::max(dy, std::fabs(yb[k] - yo[k]));
      // At equal c, elements, e and rho, the manifold state's T differs from the trajectory's.
      std::printf(" %+9.2f %9.2e %12.4e |", tb - to, dy, lb / lo - 1);
      pe.push_back({departure(tb, rho, yb), departure(to, rho, yo)});
    }
    std::printf("\n%-22s   chain-reaction departure |ln(fwd/rev)|, branch and off-manifold, at phi 0, 0.1, 0.5, 0.9:", "");
    for (auto& [b, o] : pe) std::printf("  %.2e %.2e", b, o);
    double tHalf = 0;
    for (std::size_t i = 0; i < off.c.size(); ++i) if ((off.c[i] - c0) >= 0.5 * (off.ceq - c0)) { tHalf = off.t[i]; break; }
    std::printf("\n%-22s   off-manifold time to phi 0.5: %.2e s\n", "", tHalf);
  }
  std::printf("\ncolumns at each phi: T(branch) - T(off-manifold), the largest |Y_k| difference, lambda(branch)/lambda(off) - 1\n");

  // (3) A parcel expanding gradually from the chamber, rho = rho_c exp(-t / tau), split as the engine
  // splits it: a frozen adiabatic expansion (de = -p d(1/rho)) then reaction at fixed (e, rho), 200
  // splits per e-folding. At rho/rho_c 0.3, 0.1, 0.03 and 0.01 the parcel is compared with the lower
  // branch of its node at its own c: the snap the engine would make there.
  std::printf("\n(3) Gradual expansion from the chamber, the parcel against its node's lower branch at equal c\n");
  // The slow variable: moles per kilogram, sum Y_k / W_k. The chain reactions conserve moles, so
  // only recombination and dissociation change it.
  std::vector<double> w(n);
  gas.getMolecularWeights(w.data());
  auto moleRate = [&](double rho, const std::vector<double>& zz) {
    std::vector<double> d(n + 1);
    source.rates(rho, zz.data(), d.data());
    double r = 0;
    for (std::size_t k = 0; k < n; ++k) r += d[1 + k] / w[k];
    return r;
  };
  std::size_t ih = 0, ioh = 0;
  for (std::size_t k = 0; k < n; ++k) {
    if (source.speciesName(k) == "H") ih = k;
    if (source.speciesName(k) == "OH") ioh = k;
  }
  std::printf("%10s %9s %9s %9s %9s %10s %10s %11s %11s %11s %21s %14s\n", "tau s", "rho/rho_c", "T K", "T_eq K", "phi", "dT snap K",
              "max|dY|", "lam ratio-1", "Y_H ratio-1", "Y_OH ratio-1", "departure parcel, branch", "dN/dt ratio-1");
  for (double tau : {3e-6, 3e-5, 3e-4}) {
    std::vector<double> z(n + 1);
    z[0] = tc;
    std::copy(chamber.begin(), chamber.end(), z.begin() + 1);
    double rho = rhoc;
    const double dt = tau / 200;
    std::size_t target = 0;
    const double targets[] = {0.3, 0.1, 0.03, 0.01};
    while (target < 4) {
      double rhoNew = rho * std::exp(-dt / tau);
      std::vector<double> y(z.begin() + 1, z.end());
      double p = source.pressure(z[0], rho, y.data());
      double e = energy(z[0], rho, y) - p * (1 / rhoNew - 1 / rho);
      rho = rhoNew;
      z[0] = temperatureAt(e, rho, y);
      step.advance(rho, z.data(), dt);
      if (rho / rhoc <= targets[target]) {
        y.assign(z.begin() + 1, z.end());
        e = energy(z[0], rho, y);
        std::vector<double> zeq = z;
        source.equilibrateUV(rho, zeq.data());
        bool fromReactants = false;
        auto branch = integrate(rho, lowerStart(e, rho, fromReactants));
        double c = z[1 + iw], tb, lb;
        std::vector<double> yb;
        double phi = (c - branch.c.front()) / (branch.ceq - branch.c.front());
        if (atC(branch, c, tb, yb, lb)) {
          double dy = 0;
          for (std::size_t k = 0; k < n; ++k) dy = std::max(dy, std::fabs(yb[k] - y[k]));
          double lo = rate(rho, z) / (branch.ceq - c);
          std::vector<double> zb(n + 1);
          zb[0] = tb;
          std::copy(yb.begin(), yb.end(), zb.begin() + 1);
          std::printf("%10.1e %9.3f %9.1f %9.1f %9.4f %+10.2f %10.2e %+11.3e %+11.3e %+11.3e %10.2e %10.2e %+14.3e\n", tau, rho / rhoc,
                      z[0], zeq[0], phi, tb - z[0], dy, lb / lo - 1, yb[ih] / y[ih] - 1, yb[ioh] / y[ioh] - 1,
                      departure(z[0], rho, y), departure(tb, rho, yb), moleRate(rho, zb) / moleRate(rho, z) - 1);
          // The same snap at equal moles per kilogram (the slow variable) instead of equal c.
          auto moles = [&](const std::vector<double>& yy) {
            double m = 0;
            for (std::size_t k = 0; k < n; ++k) m += yy[k] / w[k];
            return m;
          };
          const double np = moles(y);
          bool found = false;
          for (std::size_t i = 1; i < branch.c.size() && !found; ++i) {
            double a = moles(branch.y[i - 1]), b = moles(branch.y[i]);
            if ((np - a) * (np - b) <= 0 && a != b) {
              double f = (np - a) / (b - a), tn = branch.temperature[i - 1] + f * (branch.temperature[i] - branch.temperature[i - 1]);
              std::vector<double> yn(n);
              for (std::size_t k = 0; k < n; ++k) yn[k] = branch.y[i - 1][k] + f * (branch.y[i][k] - branch.y[i - 1][k]);
              double dyn = 0;
              for (std::size_t k = 0; k < n; ++k) dyn = std::max(dyn, std::fabs(yn[k] - y[k]));
              std::vector<double> zn(n + 1);
              zn[0] = tn;
              std::copy(yn.begin(), yn.end(), zn.begin() + 1);
              std::printf("%10s %9s   snap at equal moles: dT %+.2f K, dc %+.3e, max|dY| %.2e, Y_H ratio-1 %+.3e, dN/dt ratio-1 %+.3e\n",
                          "", "", tn - z[0], yn[iw] - c, dyn, yn[ih] / y[ih] - 1, moleRate(rho, zn) / moleRate(rho, z) - 1);
              found = true;
            }
          }
          std::printf("%10s %9s   c %.6f, c_eq %.6f, chamber c %.6f\n", "", "", c, branch.ceq, chamber[iw]);
          if (!found) std::printf("%10s %9s   snap at equal moles: N %.6f not on the branch\n", "", "", np);
        } else {
          std::printf("%10.1e %9.3f %9.1f %9.1f   c %.6f outside the branch [%.6f, %.6f]\n", tau, rho / rhoc, z[0], zeq[0], c,
                      branch.c.front(), branch.ceq);
        }
        ++target;
      }
    }
  }
  std::printf("\nphi is the parcel's place on its node's lower branch, from its start (0) to equilibrium (1)\n");

  // (4) Burnt gas mixed into fresh gas (6 October 2026, 03:01 CDT; TABLE_B.md, the handover). A mass
  // fraction m of the chamber's equilibrium products is mixed into the fresh reactants at 300 K and
  // 3 MPa as a finite-volume cell mixes them: Y, e and 1/rho are mass-weighted. The mixture then
  // reacts at constant (e, rho) (CVODES: the truth). The table would snap it onto its node's lower
  // branch and advance it along the branch, whose own clock is the table's (lambda is the branch's
  // dc/dt). For the snap at equal c, at equal moles and by the stated rule, the probe reports the time
  // to reach phi 0.5 and 0.9 (from the mixture's c toward c_eq) against the truth's.
  // Trajectories here are sampled densely in c (steps halved until c moves by at most 0.5% of its
  // span), with an absolute tolerance of 1e-30 and no time limit short of 1e15 s, as the builder will
  // integrate cold nodes.
  ReactionStep coldStep(source, 1e-10, 1e-30);
  auto integrateDense = [&](double rho, std::vector<double> z) {
    Trajectory tr;
    std::vector<double> zeq = z;
    source.equilibrateUV(rho, zeq.data());
    tr.ceq = zeq[1 + iw];
    const double c0 = z[1 + iw], dcMax = 0.005 * std::fabs(tr.ceq - c0);
    auto record = [&](double time) {
      tr.t.push_back(time); tr.c.push_back(z[1 + iw]); tr.temperature.push_back(z[0]);
      tr.lambda.push_back(rate(rho, z) / (tr.ceq - z[1 + iw]));
      tr.y.emplace_back(z.begin() + 1, z.end());
    };
    record(0);
    double t = 0, dt = 1e-10;
    while (t < 1e15 && std::fabs(z[1 + iw] - tr.ceq) > 1e-6 * tr.ceq) {
      auto saved = z;
      coldStep.advance(rho, z.data(), dt);
      double dc = std::fabs(z[1 + iw] - saved[1 + iw]);
      if (dc > dcMax && dt > 1e-13) { z = saved; dt *= 0.5; continue; }
      t += dt;
      record(t);
      if (dc < 0.25 * dcMax) dt *= 2;
    }
    return tr;
  };
  auto atCTime = [&](const Trajectory& tr, double c) {
    for (std::size_t i = 1; i < tr.c.size(); ++i) {
      double a = tr.c[i - 1], b = tr.c[i];
      if ((c - a) * (c - b) <= 0 && a != b) return tr.t[i - 1] + (c - a) / (b - a) * (tr.t[i] - tr.t[i - 1]);
    }
    return -1.0;
  };
  auto molesAtC = [&](const Trajectory& tr, double c) {
    double tt, ll;
    std::vector<double> yy;
    return atC(tr, c, tt, yy, ll) ? molesOf(yy) : -1.0;
  };
  std::printf("\n(4) Burnt gas mixed into fresh gas at 300 K, then reacting at constant (e, rho)\n");
  const double rhoFresh = source.density(300, 3e6, reactants.data());
  const double eFresh = energy(300, rhoFresh, reactants), eBurnt = energy(tc, rhoc, chamber);
  std::printf("fresh: e %.4e J/kg, rho %.4f; burnt: e %.4e J/kg, rho %.4f\n", eFresh, rhoFresh, eBurnt, rhoc);
  for (double m : {0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9}) {
    std::vector<double> y(n);
    for (std::size_t k = 0; k < n; ++k) y[k] = m * chamber[k] + (1 - m) * reactants[k];
    const double e = m * eBurnt + (1 - m) * eFresh, rho = 1 / (m / rhoc + (1 - m) / rhoFresh);
    std::vector<double> z(n + 1);
    z[0] = temperatureAt(e, rho, y);
    std::copy(y.begin(), y.end(), z.begin() + 1);
    const double cMix = y[iw], nMix = molesOf(y);
    auto truth = integrateDense(rho, z);
    bool fromReactants = false;
    auto start = lowerStart(e, rho, fromReactants);
    auto branch = integrateDense(rho, start);
    const double c0 = branch.c.front(), ceq = branch.ceq;
    const double n0 = molesOf(branch.y.front()), neq = molesOf(branch.y.back());
    std::printf("m %.2f: T_mix %.1f K, rho %.4f, c_mix %.5f, N_mix %.6f | truth: c_eq %.5f, T_eq %.1f K, %zu samples | "
                "branch from %s at %.1f K: c_0 %.5f, N_0 %.6f, N_eq %.6f, %zu samples\n",
                m, z[0], rho, cMix, nMix, truth.ceq, truth.temperature.back(), truth.c.size(),
                fromReactants ? "reactants" : "frozen eq", start[0], c0, n0, neq, branch.c.size());
    // The rule: w from the branch's share of its fall in N at the cell's c; c_N where the branch's N
    // equals the cell's, sought beyond c_a (0.5% of the fall) and clamped there.
    double phiN = -1, w = -1, ca = -1, cN = -1;
    if (cMix >= c0 && cMix <= ceq) {
      phiN = (n0 - molesAtC(branch, cMix)) / (n0 - neq);
      w = std::clamp((phiN - 0.005) / 0.005, 0.0, 1.0);
    }
    for (std::size_t i = 0; i < branch.c.size(); ++i) {
      double fall = (n0 - molesOf(branch.y[i])) / (n0 - neq);
      if (ca < 0 && fall >= 0.005) ca = branch.c[i];
      if (ca >= 0 && i > 0 && cN < 0) {
        double a = molesOf(branch.y[i - 1]), b = molesOf(branch.y[i]);
        if ((nMix - a) * (nMix - b) <= 0 && a != b) cN = branch.c[i - 1] + (nMix - a) / (b - a) * (branch.c[i] - branch.c[i - 1]);
      }
    }
    if (cN < 0 && ca >= 0) cN = nMix > n0 ? ca : ceq;  // clamped
    std::printf("        rule: phi_N(c_mix) %.4f, w %.3f, c_a %.5f, c_N %.5f (c_N/c_eq %.4f)\n", phiN, w, ca, cN, cN / ceq);
    struct Snap { const char* name; double c; };
    std::vector<Snap> snaps = {{"equal c", cMix}, {"equal N", cN}};
    if (w >= 0) snaps.push_back({"rule", (1 - w) * cMix + w * cN});
    for (double phi : {0.5, 0.9}) {
      const double cTarget = cMix + phi * (truth.ceq - cMix);
      const double tTruth = atCTime(truth, cTarget);
      double tt, ll;
      std::vector<double> yt, yb;
      bool okT = atC(truth, cTarget, tt, yt, ll);
      double tb = 0;
      bool okB = atC(branch, cTarget, tb, yb, ll);
      std::printf("        phi %.1f: truth t %.3e s, N moved %+.3e of N_mix; T(branch) - T(truth) at this c %s", phi, tTruth,
                  okT ? molesOf(yt) / nMix - 1 : 0.0, okB && okT ? "" : "n/a");
      if (okB && okT) std::printf("%+.2f K", tb - tt);
      for (auto& s : snaps) {
        double from = atCTime(branch, s.c), to = atCTime(branch, cTarget);
        if (s.c < c0 || s.c > ceq || from < 0 || to < 0) { std::printf(" | %s: off the branch", s.name); continue; }
        std::printf(" | %s: t %.3e s (ratio %.3f)", s.name, std::max(0.0, to - from), std::max(0.0, to - from) / tTruth);
      }
      std::printf("\n");
    }
    // The snap's own jump in T at each target.
    for (auto& s : snaps) {
      double tb, ll;
      std::vector<double> yb;
      if (s.c >= c0 && s.c <= ceq && atC(branch, s.c, tb, yb, ll)) {
        double dy = 0;
        for (std::size_t k = 0; k < n; ++k) dy = std::max(dy, std::fabs(yb[k] - y[k]));
        // The truth where it passes the snap's c: how far the snapped state is from it there.
        double tt2, ll2;
        std::vector<double> yt2;
        bool ok = atC(truth, s.c, tt2, yt2, ll2);
        std::printf("        snap %s: dT %+.2f K, max|dY| %.2e; truth at the same c: T(snap) - T(truth) %s", s.name, tb - z[0], dy, ok ? "" : "n/a\n");
        if (ok) std::printf("%+.2f K\n", tb - tt2);
      }
    }
  }
  return 0;
}
