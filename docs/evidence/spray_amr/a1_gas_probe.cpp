// Probe: is a single-level Amr (level 1 inactive) bit-identical to a Flow on the AMR base's mesh, gas alone?
#include <cstdio>
#include <cstring>
#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/amr.hpp"
using namespace crucible;
Definition chamber(thermo::ReactionSource& source, int nz, int nr, double length, double radius) {
  Definition d;
  d.experiment = Case::Chamber; d.nz = nz; d.nr = nr;
  d.contour = {{0, radius}, {length, radius}};
  d.species = source.medium().species();
  d.composition = source.massFractions("H2:1");
  d.backPressure = 1e6; d.ambientTemperature = 300;
  d.transport = thermo::transportFits("h2o2.yaml", 60);
  Supply s; s.innerRadius = radius / 6; s.outerRadius = radius / 2; s.massFlow = 3e-3; s.totalTemperature = 300;
  s.composition = source.massFractions("H2:1");
  d.supplies = {s};
  return d;
}
int main() {
  thermo::ReactionSource source("h2o2.yaml");
  Definition def = chamber(source, 16, 6, 0.016, 0.01);
  Amr amr(def, 1);
  Definition fine = def; fine.nz <<= 1; fine.nr <<= 1;
  Flow flow(def, Mesh::coarsened(Mesh(fine)));
  long n = 0;
  while (flow.time() < 0.1e-3) { double a = amr.step(0.1e-3 - amr.time()); double b = flow.step(0.1e-3 - flow.time()); if (a != b) { std::printf("dt differs at %ld: %.17g %.17g\n", n, a, b); return 1; } ++n; }
  const auto& x = amr.level(0).state(); const auto& y = flow.state();
  bool same = std::memcmp(x.data(), y.data(), x.size() * sizeof x[0]) == 0 &&
              std::memcmp(amr.level(0).partialDensities().data(), flow.partialDensities().data(), flow.partialDensities().size() * 8) == 0;
  std::printf("%ld steps, t %.17g vs %.17g; state and partials bit-identical: %s\n", n, amr.time(), flow.time(), same ? "yes" : "no");
  return same ? 0 : 1;
}
