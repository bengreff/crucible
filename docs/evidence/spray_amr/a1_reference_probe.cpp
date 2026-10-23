// A1 reference-driver probe (SPRAY.md, Spray on AMR: results): S3 chamber, frozen, against a single-level Amr with the spray.
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
int main(int argc, char** argv) {
  // 0: ReactingFlow (frozen) with the spray; 1: ReactingFlow, no injectors; 2: Flow::step + Spray::advance; 3: mode 2 with a
  // stableDt() before each Flow::step. Each mode stops at the first step that differs from the single-level Amr.
  const int mode = std::atoi(argv[1]);
  thermo::ReactionSource source("h2o2.yaml");
  const Definition d = chamber(source, false);
  Definition fine = d; fine.nz <<= 1; fine.nr <<= 1;
  Flow flow(d, Mesh::coarsened(Mesh(fine)));
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", 1, 1e-6, 1e-12, thermo::Chemistry::Frozen);
  auto s = settings(); if (mode == 1) s.injectors.clear();
  Spray uniform(flow, oxygen(flow), s);
  if (mode < 2) reacting.setSpray(&uniform);
  Amr amr(d, 1);
  Spray spray(amr, oxygen(amr.level(0)), s);
  amr.setAfterStep([&](int l, double dt) { spray.afterStep(l, dt); });
  amr.setRedo([&] { spray.redo(); });
  const double end = 0.5e-3;
  for (long n = 0; n < 6000 && amr.time() < end; ++n) {
    const double t = amr.time();
    const double a = amr.step(std::min(end - t, spray.nextEvent(t) - t));
    double b;
    if (mode >= 2) { if (mode == 3) (void)flow.stableDt(); b = flow.step(end - flow.time()); uniform.advance(b); } else b = reacting.step(end - flow.time());
    const bool gas = same(amr.level(0).state(), flow.state()) && same(amr.level(0).partialDensities(), flow.partialDensities());
    const bool par = same(bytes(spray.parcels()), bytes(uniform.parcels()));
    if (a != b || !gas || !par) {
      std::printf("step %ld: dt %.17g vs %.17g; gas same %d, parcels same %d (%zu vs %zu)\n", n, a, b, gas, par, spray.parcels().size(), uniform.parcels().size());
      if (!par) {
        for (std::size_t i = 0; i < std::min(spray.parcels().size(), uniform.parcels().size()); ++i) {
          const auto& p = spray.parcels()[i]; const auto& q = uniform.parcels()[i];
          if (std::memcmp(&p, &q, offsetof(Parcel, phase)) != 0) { std::printf("  parcel %zu: z %.17g %.17g r %.17g %.17g uz %.17g %.17g d %.17g %.17g t %.17g %.17g cell %d %d\n", i, p.z, q.z, p.r, q.r, p.uz, q.uz, p.d, q.d, p.t, q.t, p.cell, q.cell); break; }
        }
      }
      if (!gas) std::printf("  gas difference %.3e\n", stateDifference(amr.level(0), flow));
      return 0;
    }
  }
  std::printf("no divergence\n");
}
