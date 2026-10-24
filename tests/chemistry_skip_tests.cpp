// Verification of the chemistry skip (ReactingFlow::setChemistrySkip; docs/evidence/SPEED.md, "The chemistry
// skip"). Criteria stated 24 October 2026, before the skip was first built or run. Mechanisms h2o2.yaml and the
// engine's tools/mechanisms/h2o2_highT.yaml; the tolerances of tier 1 (rtol 1e-6, atol 1e-12).
// The skip, per FiniteRate cell and reaction substep: a species is reactive if some reaction has a nonzero net
// coefficient for it. A cell with no reactive species present (Y_k != 0) is left unchanged; so is a cell at or
// below 500 K whose one reactive species is H2 (or O2), if every reaction direction that needs no reactive
// species but H2 (O2) changes no reactive species but H2 and H (O2 and O).
//   1. The mechanisms. Written here independently of the engine: the reactive species of each mechanism are all
//      but N2 and AR (and HE where present); the reaction directions H2 alone can start number at least one and
//      change only H2 and H, and those O2 alone can start at least one and change only O2 and O.
//   2. The bound. Pure H2 (O2) can only dissociate, and from no atoms the composition moves toward the
//      dissociation equilibrium, so the atom's equilibrium mole fraction bounds the change at any step length.
//      x_a <= sqrt(Kp p_ref / p) with Kp = exp(-(2 g_a - g_m) / (R T)), g from the medium's NASA-7 data (valid
//      from 200 K in both mechanisms), p_ref 1e5 Pa. At 500 K and 1 Pa: below 1e-16 for H and O.
//   3. The integrator agrees. ReactionStep over 1e-3 s (about 1e4 tier-1 steps) from X_H2 1, X_H2 0.5 with N2,
//      X_O2 1 and X_O2 0.3 with N2, at 200, 300 and 500 K and 1e2, 1e5 and 1e7 Pa: every Y_k within 1e-16 of its
//      start and T within 1e-14 relative.
//   4. The flow. A duct at rest (R 5 mm, L 40 mm, nz 16, nr 4) at 101325 Pa, in zones of axial columns: 0-3
//      X_H2 0.5 and N2 at 300 K; 4-7 N2 at 1500 K; 8-11 X_O2 0.3 and N2 at 300 K; 12-13 X_H2 0.5 and N2 at
//      600 K; 14-15 X_H2 0.2, X_O2 0.1 and N2 at 300 K. One reaction substep of 1e-6 s with the skip on and
//      one with it off, each on its own flow from the same state.
//      (a) Skip on: skippedCells 48 (columns 0-11) and integratedCells 16; off: 0 and 64.
//      (b) Skip on: every cell of columns 0-11 keeps its partial densities and conserved state bit for bit.
//      (c) On against off: every partial density within 1e-15 of the cell's density; columns 12-15 bit for bit.
//   5. Trace atoms (stated 24 October 2026 at about 01:10 CDT, after tier 1's lead showed 7% of level 3's updates
//      skipped, CVODES leaving traces of H in the H2 it integrates, and before the extension was run): the
//      molecule may carry its atom at 0 < Y <= 1e-14 (ReactingFlow::kAtomTrace).
//      (a) The mechanisms: every direction whose reactants lie within {H2, H} and the inert species changes only H2
//          and H, at least one of them (the recombinations); likewise {O2, O}.
//      (b) A duct as in 4, columns: 0-3 X_H2 0.5 and N2 with Y_H 1e-15 at 300 K; 4-7 the same with Y_H 1e-10;
//          8-11 X_O2 0.3 and N2 with Y_O 1e-15 at 300 K; 12-13 as 0-3 at 600 K; 14-15 as 0-3 with Y_H -1e-15
//          (Y_H2 making up the sum). Skip on: skippedCells 32 (columns 0-3 and 8-11), integratedCells 32; off: 0, 64.
//      (c) Skip on: columns 0-3 and 8-11 keep their partial densities and conserved state bit for bit. On against
//          off: every partial density within 1e-14 of the cell's density; columns 4-7 and 12-15 bit for bit.
//   Run 1 (24 October 2026): 4(a) failed, integratedCells 0 with the skip off, since the engine counted the
//   integrations only with the skip on, against the stat's own definition (every CVODES integration). The counter
//   fixed; the criteria unchanged. Run 2: every check passes.
//   Criterion 5, run 1 (24 October 2026, the trace extension): 5(a) passes; 5(b) and 5(c) fail in both mechanisms,
//   skipped 40 and integrated 24, and columns 14-15 not bit for bit (gap 1.0e-15 of the density). The engine built
//   the skip's mask from Flow::massFractions, which clips a negative partial density to 0, so the cells with Y_H
//   -1e-15 were skipped and kept the negative partial density the integrator's path clips (4a7628f too, against its
//   own comment). The mask now reads the partial densities; the criteria unchanged. Run 2: every check passes.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/flow.hpp"

using namespace crucible;

namespace {

int failures = 0;
void check(bool ok, const char* what, double value, double limit) {
  std::printf("  %-74s %.3e (limit %.1e)  %s\n", what, value, limit, ok ? "ok" : "FAIL");
  if (!ok) ++failures;
}

constexpr double p0 = 101325;

// g / (R T) of a species from its NASA-7 data.
double gibbs(const Species& s, double t) {
  const auto& a = t < s.tMid ? s.low : s.high;
  const double h = a[0] + a[1] * t / 2 + a[2] * t * t / 3 + a[3] * t * t * t / 4 + a[4] * t * t * t * t / 5 + a[5] / t;
  const double sr = a[0] * std::log(t) + a[1] * t + a[2] * t * t / 2 + a[3] * t * t * t / 3 + a[4] * t * t * t * t / 4 + a[6];
  return h - sr;
}

void mechanism(const std::string& name) {
  std::printf("%s\n", name.c_str());
  thermo::ReactionSource source(name);
  Medium medium = source.medium();
  const std::size_t ns = medium.size();
  auto index = [&](const char* s) {
    for (std::size_t k = 0; k < ns; ++k) if (medium.species()[k].name == s) return k;
    throw std::runtime_error(std::string("missing species ") + s);
  };
  const auto reactions = source.reactionSpecies();

  // 1. Reactive species, and the directions one molecule alone can start.
  std::vector<bool> reactive(ns, false);
  for (const auto& r : reactions)
    for (auto k : r.changed) reactive[k] = true;
  std::string inert;
  bool inertOk = true;
  for (std::size_t k = 0; k < ns; ++k) {
    const auto& n = medium.species()[k].name;
    const bool expected = n == "N2" || n == "AR" || n == "HE";
    if (!reactive[k]) inert += n + " ";
    if (reactive[k] == expected) inertOk = false;
  }
  std::printf("  inert species: %s\n", inert.c_str());
  check(inertOk, "1. reactive species are all but N2, AR, HE (1 = differs)", inertOk ? 0 : 1, 0);
  for (const auto& [m, a] : {std::pair{"H2", "H"}, std::pair{"O2", "O"}}) {
    const std::size_t im = index(m), ia = index(a);
    int started = 0;
    bool only = true;
    auto alone = [&](const std::vector<std::size_t>& set) {
      for (auto k : set) if (k != im && reactive[k]) return false;
      return true;
    };
    auto direction = [&](const std::vector<std::size_t>& from, const std::vector<std::size_t>& to, std::size_t i,
                         const char* way) {
      if (!alone(from)) return;
      bool changesReactive = false;
      for (auto k : to) if (reactive[k] && k != im) changesReactive = true;
      for (auto k : from) if (reactive[k]) changesReactive = true;
      if (!changesReactive) return;  // changes nothing
      ++started;
      for (auto k : to) if (reactive[k] && k != im && k != ia) only = false;
      std::printf("  %s alone starts reaction %zu %s\n", m, i, way);
    };
    for (std::size_t i = 0; i < reactions.size(); ++i) {
      direction(reactions[i].reactants, reactions[i].products, i, "forward");
      if (reactions[i].reversible) direction(reactions[i].products, reactions[i].reactants, i, "in reverse");
    }
    char what[96];
    std::snprintf(what, sizeof what, "1. %s alone: directions started (at least 1), all to %s and %s only", m, m, a);
    check(started >= 1 && only, what, started, 1);
  }

  // 2. The equilibrium bound at 500 K and 1 Pa.
  for (const auto& [m, a] : {std::pair{"H2", "H"}, std::pair{"O2", "O"}}) {
    const auto& sm = medium.species()[index(m)];
    const auto& sa = medium.species()[index(a)];
    const double t = 500, p = 1;
    const double kp = std::exp(-(2 * gibbs(sa, t) - gibbs(sm, t)));
    const double x = std::sqrt(kp * 1e5 / p);
    char what[96];
    std::snprintf(what, sizeof what, "2. equilibrium x_%s in pure %s at 500 K and 1 Pa (Kp %.2e)", a, m, kp);
    check(x < 1e-16, what, x, 1e-16);
  }

  // 3. The integrator over 1e-3 s.
  double worstY = 0, worstT = 0;
  for (const char* moles : {"H2:1", "H2:0.5, N2:0.5", "O2:1", "O2:0.3, N2:0.7"})
    for (double t : {200.0, 300.0, 500.0})
      for (double p : {1e2, 1e5, 1e7}) {
        std::vector<double> z(ns + 1);
        const auto y = source.massFractions(moles);
        z[0] = t;
        std::copy(y.begin(), y.end(), z.begin() + 1);
        const double rho = source.density(t, p, y.data());
        const auto start = z;
        thermo::ReactionStep step(source, 1e-6, 1e-12);
        step.advance(rho, z.data(), 1e-3);
        for (std::size_t k = 0; k < ns; ++k) worstY = std::max(worstY, std::abs(z[k + 1] - start[k + 1]));
        worstT = std::max(worstT, std::abs(z[0] - start[0]) / start[0]);
      }
  check(worstY <= 1e-16, "3. ReactionStep over 1e-3 s: largest |dY_k|", worstY, 1e-16);
  check(worstT <= 1e-14, "3. ReactionStep over 1e-3 s: largest |dT| / T", worstT, 1e-14);

  // 5(a). The directions within {molecule, atom}.
  for (const auto& [m, a] : {std::pair{"H2", "H"}, std::pair{"O2", "O"}}) {
    const std::size_t im = index(m), ia = index(a);
    int started = 0;
    bool only = true;
    auto within = [&](const std::vector<std::size_t>& set) {
      for (auto k : set) if (k != im && k != ia && reactive[k]) return false;
      return true;
    };
    auto direction = [&](const std::vector<std::size_t>& from, const std::vector<std::size_t>& to) {
      if (!within(from)) return;
      ++started;
      if (!within(to)) only = false;
    };
    for (const auto& r : reactions) {
      direction(r.reactants, r.products);
      if (r.reversible) direction(r.products, r.reactants);
    }
    char what[96];
    std::snprintf(what, sizeof what, "5(a) within %s and %s: directions (at least 1), all to %s and %s only", m, a, m, a);
    check(started >= 1 && only, what, started, 1);
  }
}

// The duct of 4.
struct Duct {
  Definition definition;
  std::vector<Primitive> cells;
  std::vector<double> y;
  std::vector<int> column;
};

Duct duct(thermo::ReactionSource& source, const Medium& medium, const TransportFits& fits) {
  Duct d;
  auto& def = d.definition;
  def.experiment = Case::UniformDuct;
  def.nz = 16; def.nr = 4; def.length = 0.04; def.inletRadius = def.exitRadius = def.throatRadius = 0.005;
  def.species = medium.species();
  const std::size_t n = medium.size();
  def.composition.assign(n, 0.0);
  for (std::size_t k = 0; k < n; ++k) def.composition[k] = medium.species()[k].name == "N2";
  def.totalPressure = def.backPressure = p0;
  def.transport = fits;
  def.wallSlip = true; def.wallTemperature = 0;
  Mesh m(def);
  d.cells.resize(m.cells.size());
  d.y.resize(m.cells.size() * n);
  d.column.resize(m.cells.size());
  for (std::size_t q = 0; q < m.cells.size(); ++q) {
    const int c = std::min(15, static_cast<int>(m.cells[q].z / (def.length / def.nz)));
    d.column[q] = c;
    const char* moles = c < 4 ? "H2:0.5, N2:0.5" : c < 8 ? "N2:1" : c < 12 ? "O2:0.3, N2:0.7" : c < 14 ? "H2:0.5, N2:0.5"
                                                                                                       : "H2:0.2, O2:0.1, N2:0.7";
    const double t = c < 4 ? 300 : c < 8 ? 1500 : c < 12 ? 300 : c < 14 ? 600 : 300;
    const auto y = source.massFractions(moles);
    std::copy(y.begin(), y.end(), d.y.begin() + static_cast<std::ptrdiff_t>(q * n));
    d.cells[q] = {p0 / (medium.gasConstant(y.data()) * t), 0, 0, p0};
  }
  return d;
}

void flow(const std::string& name) {
  std::printf("%s, the duct\n", name.c_str());
  thermo::ReactionSource source(name);
  Medium medium = source.medium();
  const auto fits = thermo::transportFits(name);
  medium.setTransport(fits);
  const std::size_t ns = medium.size();
  Duct d = duct(source, medium, fits);
  Flow on(d.definition), off(d.definition);
  for (Flow* f : {&on, &off}) f->setInitialState(d.cells, d.y);
  const auto before = on.partialDensities();
  const auto beforeState = on.state();
  thermo::ReactingFlow a(on, name, 2, 1e-6, 1e-12, thermo::Chemistry::FiniteRate);
  thermo::ReactingFlow b(off, name, 2, 1e-6, 1e-12, thermo::Chemistry::FiniteRate);
  b.setChemistrySkip(false);
  a.react(1e-6);
  b.react(1e-6);
  const auto& sa = a.stats();
  const auto& sb = b.stats();
  std::printf("  skip on: skipped %ld, integrated %ld; off: skipped %ld, integrated %ld\n", sa.skippedCells,
              sa.integratedCells, sb.skippedCells, sb.integratedCells);
  const bool counts = sa.skippedCells == 48 && sa.integratedCells == 16 && sb.skippedCells == 0 && sb.integratedCells == 64;
  check(counts, "4(a) the counts (1 = differ)", counts ? 0 : 1, 0);
  const auto& pa = on.partialDensities();
  const auto& pb = off.partialDensities();
  bool kept = true, same = true;
  double gap = 0;
  for (std::size_t q = 0; q < d.cells.size(); ++q) {
    const bool skipped = d.column[q] < 12;
    const std::size_t o = q * ns;
    if (skipped && (std::memcmp(pa.data() + o, before.data() + o, ns * sizeof(double)) != 0 ||
                    std::memcmp(on.state()[q].data(), beforeState[q].data(), sizeof(beforeState[q])) != 0))
      kept = false;
    if (!skipped && std::memcmp(pa.data() + o, pb.data() + o, ns * sizeof(double)) != 0) same = false;
    for (std::size_t k = 0; k < ns; ++k) gap = std::max(gap, std::abs(pa[o + k] - pb[o + k]) / on.state()[q][0]);
  }
  check(kept, "4(b) skipped cells unchanged bit for bit (1 = changed)", kept ? 0 : 1, 0);
  check(gap <= 1e-15, "4(c) on against off: largest partial density gap / rho", gap, 1e-15);
  check(same, "4(c) columns 12-15 on against off bit for bit (1 = differ)", same ? 0 : 1, 0);
}

// 5(b) and (c): the duct with trace atoms.
void traces(const std::string& name) {
  std::printf("%s, the duct with trace atoms\n", name.c_str());
  thermo::ReactionSource source(name);
  Medium medium = source.medium();
  const auto fits = thermo::transportFits(name);
  medium.setTransport(fits);
  const std::size_t ns = medium.size();
  auto index = [&](const char* s) {
    for (std::size_t k = 0; k < ns; ++k) if (medium.species()[k].name == s) return k;
    throw std::runtime_error(std::string("missing species ") + s);
  };
  const std::size_t h2 = index("H2"), h = index("H"), o2 = index("O2"), o = index("O");
  Duct d = duct(source, medium, fits);
  for (std::size_t q = 0; q < d.cells.size(); ++q) {
    const int c = d.column[q];
    const bool oxygen = c >= 8 && c < 12;
    auto y = source.massFractions(oxygen ? "O2:0.3, N2:0.7" : "H2:0.5, N2:0.5");
    const double trace = c >= 4 && c < 8 ? 1e-10 : c >= 14 ? -1e-15 : 1e-15;
    y[oxygen ? o : h] = trace;
    y[oxygen ? o2 : h2] -= trace;
    const double t = c >= 12 && c < 14 ? 600 : 300;
    std::copy(y.begin(), y.end(), d.y.begin() + static_cast<std::ptrdiff_t>(q * ns));
    d.cells[q] = {p0 / (medium.gasConstant(y.data()) * t), 0, 0, p0};
  }
  Flow on(d.definition), off(d.definition);
  for (Flow* f : {&on, &off}) f->setInitialState(d.cells, d.y);
  const auto before = on.partialDensities();
  const auto beforeState = on.state();
  thermo::ReactingFlow a(on, name, 2, 1e-6, 1e-12, thermo::Chemistry::FiniteRate);
  thermo::ReactingFlow b(off, name, 2, 1e-6, 1e-12, thermo::Chemistry::FiniteRate);
  b.setChemistrySkip(false);
  a.react(1e-6);
  b.react(1e-6);
  const auto& sa = a.stats();
  const auto& sb = b.stats();
  std::printf("  skip on: skipped %ld, integrated %ld; off: skipped %ld, integrated %ld\n", sa.skippedCells,
              sa.integratedCells, sb.skippedCells, sb.integratedCells);
  const bool counts = sa.skippedCells == 32 && sa.integratedCells == 32 && sb.skippedCells == 0 && sb.integratedCells == 64;
  check(counts, "5(b) the counts (1 = differ)", counts ? 0 : 1, 0);
  const auto& pa = on.partialDensities();
  const auto& pb = off.partialDensities();
  bool kept = true, same = true;
  double gap = 0;
  for (std::size_t q = 0; q < d.cells.size(); ++q) {
    const bool skipped = d.column[q] < 4 || (d.column[q] >= 8 && d.column[q] < 12);
    const std::size_t off0 = q * ns;
    if (skipped && (std::memcmp(pa.data() + off0, before.data() + off0, ns * sizeof(double)) != 0 ||
                    std::memcmp(on.state()[q].data(), beforeState[q].data(), sizeof(beforeState[q])) != 0))
      kept = false;
    if (!skipped && std::memcmp(pa.data() + off0, pb.data() + off0, ns * sizeof(double)) != 0) same = false;
    for (std::size_t k = 0; k < ns; ++k) gap = std::max(gap, std::abs(pa[off0 + k] - pb[off0 + k]) / on.state()[q][0]);
  }
  check(kept, "5(c) skipped cells unchanged bit for bit (1 = changed)", kept ? 0 : 1, 0);
  check(gap <= 1e-14, "5(c) on against off: largest partial density gap / rho", gap, 1e-14);
  check(same, "5(c) columns 4-7 and 12-15 on against off bit for bit (1 = differ)", same ? 0 : 1, 0);
}

}  // namespace

int main() {
  try {
    for (const std::string name : {std::string("h2o2.yaml"), std::string(CRUCIBLE_SOURCE_DIR "/tools/mechanisms/h2o2_highT.yaml")}) {
      mechanism(name);
      flow(name);
      traces(name);
    }
  } catch (const std::exception& e) {
    std::printf("error: %s\n", e.what());
    return 2;
  }
  std::printf("%d failures\n", failures);
  return failures == 0 ? 0 : 1;
}
