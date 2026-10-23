// A1(b) probe (SPRAY.md, Spray on AMR: results): S8(c) chamber, SST, single-level Amr against a Flow on the same mesh, to 0.2 ms.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <memory>
#include <numbers>
#include <string>
#include <vector>
#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/amr.hpp"
#include "core/spray.hpp"
using namespace crucible;
#include "common.inc"
}
// mode 0: gas alone, Flow::step; mode 1: gas alone, Flow with stableDt() before each step; mode 2: mode 1 with the spray;
// mode 3: mode 1 (gas alone) with Flow::step limited to the first stableDt()'s value, as Amr::step takes it (the second rounding source)
int main(int argc, char** argv) {
  const int mode = std::atoi(argv[1]);
  thermo::ReactionSource source("h2o2.yaml");
  const Definition d = chamber(source, true);
  Definition fine = d; fine.nz <<= 1; fine.nr <<= 1;
  Flow flow(d, Mesh::coarsened(Mesh(fine)));
  auto s = settings(); if (mode < 2 || mode == 3) s.injectors.clear();
  Spray uniform(flow, oxygen(flow), s);
  Amr amr(d, 1);
  Spray spray(amr, oxygen(amr.level(0)), s);
  amr.setAfterStep([&](int l, double dt) { spray.afterStep(l, dt); });
  amr.setRedo([&] { spray.redo(); });
  const double end = 0.2e-3;
  long n = 0, first = -1, shorter = 0;
  while (amr.time() < end) {
    const double t = amr.time();
    const double a = amr.step(std::min(end - t, spray.nextEvent(t) - t));
    const double dt1 = mode >= 1 ? flow.stableDt() : 0;
    const double b = flow.step(mode == 3 ? std::min(end - flow.time(), dt1) : end - flow.time());
    if (mode == 3 && b != dt1 && b != end - (flow.time() - b)) ++shorter;
    uniform.advance(b);
    const bool same_ = a == b && same(amr.level(0).state(), flow.state()) && same(amr.level(0).partialDensities(), flow.partialDensities()) &&
                       same(amr.level(0).turbulence(), flow.turbulence()) && same(bytes(spray.parcels()), bytes(uniform.parcels()));
    if (!same_ && first < 0) { first = n; std::printf("  first difference: step %ld, dt %s (%.17g vs %.17g)%s\n", n, a == b ? "equal" : "differs", a, b, mode == 3 && b != dt1 ? ", Flow::step's own stableDt() smaller than the first call's" : ""); }
    ++n;
  }
  if (mode == 3) std::printf("  steps where Flow::step's own stableDt() was smaller than the first call's: %ld\n", shorter);
  std::printf("mode %d: %ld steps; first difference at step %ld; end gas difference %.3e; parcels %zu vs %zu, identical %d\n", mode, n, first,
              stateDifference(amr.level(0), flow), spray.parcels().size(), uniform.parcels().size(), same(bytes(spray.parcels()), bytes(uniform.parcels())));
}
