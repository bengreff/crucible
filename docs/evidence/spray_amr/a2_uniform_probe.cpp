// A2's printed comparison (SPRAY.md, Spray on AMR: results): S3's chamber with the spray on uniform meshes, the 16 x 6 base
// and the 64 x 24 finest, to 0.5 ms, frozen, on 1 thread. Printed, not judged.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
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
int main() {
  thermo::ReactionSource source("h2o2.yaml");
  for (int shift : {0, 2}) {
    Definition d = chamber(source, false);
    d.nz <<= shift;
    d.nr <<= shift;
    Flow flow(d);
    Spray spray(flow, oxygen(flow), settings());
    const double end = 0.5e-3;
    const auto start = std::chrono::steady_clock::now();
    long n = 0;
    while (flow.time() < end) {
      const double t = flow.time();
      spray.advance(flow.step(std::min(end - t, spray.nextEvent(t) - t)));
      ++n;
    }
    const auto& l = spray.ledger();
    std::printf("uniform %d x %d: %ld steps, %.1f s; evaporated %.6e kg, %zu parcels held, %ld injected, %ld children, %ld sub-steps\n",
                d.nz, d.nr, n, seconds(start), l.exchangedMass, spray.parcels().size(), l.parcelsInjected, l.children, l.substeps);
  }
}
