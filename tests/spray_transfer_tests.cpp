// Verification of the heat-transfer number's solve near saturation (Spray::rates; docs/evidence/SPEED.md, section 6).
// Criteria stated 24 October 2026, before the fix was written, from the spray NaN of tier 1 development (1 level,
// the lead shortened to 1 ms, at about 1.35 ms; speed_2026-10-24/nan_diag.log): a 0.70 mm drop at 97.43 K with
// 1 - Y_s 7.1e-9 (B_M 7.7e7) in gas at 187.12 K and 2.034e5 Pa with Y_O2 0.450. There B_T's secant met a residual
// expm1(phi ln(1 + B_M)) - B_T past the double range, its fallback step B_T + r was infinite, and the step test
// |next - b| <= 1e-11 next held for next = inf; B_T came back infinite and the film temperature NaN.
// Mechanism tools/mechanisms/h2o2_highT.yaml, transport fits from 60 K, Spray::Settings' film points, as tier 1.
// The sweep: gas at 187.12 K and 2.034e5 Pa, Y_O2 0.45 with the rest H2, or the rest N2; a 0.70 mm drop; slip 0, 3
// and 30 m/s; drop temperatures T_sat - delta for delta from 3 K to 1e-12 K in 121 log steps, and 97.427679170131341 K.
//   (a) No call throws, and mass, heat, conductance, stiffness, Nusselt and B_T are all finite.
//   (b) Where the drop evaporates (neither boiling nor clamped), B_T solves its equation: ln(1 + B_T) =
//       phi ln(1 + B_M) is c_pF mdot = conductance B_T (Rates' definitions), within 1e-9 relative, with c_pF
//       O2's c_p at the 1/3-rule temperature as rates() takes it.
// Outside this file: (c) spray_verification prints the same lines before and after the fix (diffed), and (d) the
// tier 1 case runs past 1.4 ms (SPEED.md, section 6).
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/spray.hpp"

using namespace crucible;

namespace {

int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-74s %.3e (limit %.0e)  %s\n", what, value, limit, ok ? "PASS" : "FAIL");
  if (!ok) ++failures;
}

std::size_t index(const Medium& m, const char* name) {
  for (std::size_t k = 0; k < m.size(); ++k)
    if (m.species()[k].name == name) return k;
  throw std::runtime_error(std::string("no ") + name);
}

}  // namespace

int main() {
  try {
    const std::string mechanism = CRUCIBLE_SOURCE_DIR "/tools/mechanisms/h2o2_highT.yaml";
    thermo::ReactionSource source(mechanism);
    Definition d;
    d.experiment = Case::Chamber;
    d.nz = 4;
    d.nr = 2;
    d.contour = {{0, 0.005}, {0.01, 0.005}};
    d.species = source.medium().species();
    d.composition = source.massFractions("N2:1");
    d.backPressure = 2.034e5;
    d.ambientTemperature = 187.12;
    d.transport = thermo::transportFits(mechanism, 60);
    Flow flow(d);
    const auto& medium = flow.medium();
    const auto ns = medium.size(), o2 = index(medium, "O2");
    std::vector<double> unit(ns, 0.0);
    unit[o2] = 1;
    const auto liquid = LiquidOxygen::load(CRUCIBLE_SOURCE_DIR "/tools/o2_saturation_nist.tsv", medium.enthalpy(298.15, unit.data()));
    Spray::Settings settings;
    settings.twoWay = false;
    settings.breakup.on = false;
    const Spray spray(flow, liquid, settings);
    std::vector<double> work(4 * ns + 8);
    const double tGas = 187.12094790942112, p = 203443.90669341254, tSat = liquid.saturationTemperature(p);
    std::printf("Spray: B_T near saturation (T_sat %.9f K at %.0f Pa)\n", tSat, p);
    for (const char* carrier : {"H2", "N2"}) {
      std::vector<double> y(ns, 0.0);
      y[o2] = 0.45005158130273049;
      y[index(medium, carrier)] = 1 - y[o2];
      const double rho = p / (medium.gasConstant(y.data()) * tGas);
      for (double slip : {0.0, 3.0, 30.0}) {
        int calls = 0, thrown = 0, nonFinite = 0, evaporating = 0;
        double worst = 0, largest = 0;
        std::string first;
        for (int i = 0; i <= 121; ++i) {
          Parcel q;
          q.d = 0.70274543041639649e-3;
          q.count = 1;
          q.t = i < 121 ? tSat - 3 * std::pow(10.0, -12.477121254719663 * i / 120) : 97.427679170131341;
          if (!(q.t < tSat)) continue;
          const Spray::Gas gas{tGas, p, rho, slip, 0, y.data()};
          ++calls;
          try {
            const auto r = spray.rates(q, gas, work);
            if (!std::isfinite(r.mass) || !std::isfinite(r.heat) || !std::isfinite(r.conductance) || !std::isfinite(r.stiffness) ||
                !std::isfinite(r.nusselt) || !std::isfinite(r.transfer)) {
              ++nonFinite;
              continue;
            }
            if (r.boiling || r.clamped) continue;
            ++evaporating;
            const double tr = r.surface + (tGas - r.surface) / 3;
            const auto vapour = medium.properties(tr, unit.data());
            const double cpf = vapour.cv + vapour.r;
            const double error = std::abs(cpf * r.mass - r.conductance * r.transfer) / (cpf * r.mass);
            worst = std::max(worst, std::isfinite(error) ? error : 1.0);
            largest = std::max(largest, r.transfer);
          } catch (const std::exception& e) {
            if (thrown++ == 0) first = e.what();
          }
        }
        std::printf("  rest %s, slip %g m/s: %d calls, %d evaporating, largest B_T %.4e%s%s\n", carrier, slip, calls, evaporating, largest,
                    thrown ? "; first throw: " : "", first.c_str());
        const std::string label = std::string("rest ") + carrier + ", slip " + std::to_string(static_cast<int>(slip)) + " m/s: ";
        check(thrown == 0 && nonFinite == 0, (label + "(a) calls that threw or gave a non-finite rate").c_str(), thrown + nonFinite, 0);
        check(worst <= 1e-9 && evaporating > 0, (label + "(b) |c_pF mdot - conductance B_T| / c_pF mdot").c_str(), worst, 1e-9);
      }
    }
  } catch (const std::exception& e) {
    std::printf("error: %s\n", e.what());
    return 2;
  }
  std::printf("%s: %d failure(s)\n", failures ? "FAIL" : "PASS", failures);
  return failures ? 1 : 0;
}
