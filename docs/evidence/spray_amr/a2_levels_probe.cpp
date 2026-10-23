// A2's printed comparison, continued (SPRAY.md, Spray on AMR: results): S3's chamber with the spray on Amr(def, 2), 16 x 6
// plus 2 levels, to 0.5 ms, frozen, 1 thread. Mode 0: every level active everywhere, no regrid (the 64 x 24 finest owns
// every parcel). Mode 1: the injector-corner tag alone (z < 2 mm, r < 3 mm), regridded every 2 base steps. Mode 2: A2's
// run (the corner and the sweeping band). Modes 3 and 4: modes 1 and 0 with no injectors (the gas alone). Modes 5 and 6:
// mode 1 with one-way coupling, and with the point-implicit loads off. Every 500 base steps it prints each level's largest
// owned u_z and where. Printed, not judged.
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
  const int given = argc > 1 ? std::atoi(argv[1]) : 0;
  const bool gasAlone = given == 3 || given == 4;
  const int mode = given == 3 || given >= 5 ? 1 : given == 4 ? 0 : given;
  thermo::ReactionSource source("h2o2.yaml");
  const Definition d = chamber(source, false);
  Amr amr(d, 2);
  amr.setThreads(1);
  if (mode == 0) {
    for (int l = 1; l < amr.levels(); ++l) amr.setActive(l, std::vector<unsigned char>(amr.level(l).mesh().cells.size(), 1));
  } else {
    Amr::Tags tags;
    tags.interval = 2;
    tags.throatHalfWidth = 0;
    tags.wallYPlus = 0;
    tags.mixture = tags.temperature = tags.pressure = tags.vorticity = 0;
    const double period = 0.5e-3;
    tags.custom = [&amr, period, mode](int l, int i, int j) {
      const Mesh& m = amr.level(l).mesh();
      const double z = (i + 0.5) * m.dz, r = m.cells[m.index(i, j)].r;
      if (z < 2e-3 && r < 3e-3) return true;
      if (mode == 1) return false;
      const double x = std::fmod(amr.time(), period) / period;
      const double centre = kLength * (1 - std::abs(1 - 2 * x));
      return std::abs(z - centre) < (l == 0 ? 3e-3 : 1.5e-3);
    };
    amr.setTags(tags);
  }
  auto s = settings();
  if (gasAlone && given < 5) s.injectors.clear();
  if (given == 5) s.twoWay = false;
  if (given == 6) s.implicitCoupling = false;
  Spray spray(amr, oxygen(amr.level(0)), s);
  amr.setAfterStep([&](int l, double dt) { spray.afterStep(l, dt); });
  amr.setRedo([&] { spray.redo(); });
  const double end = 0.5e-3;
  const auto start = std::chrono::steady_clock::now();
  while (amr.time() < end) {
    const double t = amr.time();
    amr.step(std::min(end - t, spray.nextEvent(t) - t));
    if (amr.steps() % 500 == 0) {
      std::printf("  t %.3f ms, %zu parcels:", amr.time() * 1e3, spray.parcels().size());
      for (int k = 0; k <= amr.top(); ++k) {
        const Flow& f = amr.level(k);
        double best = 0, z = 0, r = 0;
        for (std::size_t q = 0; q < f.state().size(); ++q) {
          if (k > 0 && !amr.active(k)[q]) continue;
          if (k < amr.top() && amr.active(k + 1)[amr.level(k + 1).mesh().index(2 * (q / f.mesh().nr), 2 * (q % f.mesh().nr))]) continue;
          const double u = f.state()[q][1] / f.state()[q][0];
          if (std::abs(u) > std::abs(best)) { best = u; z = f.mesh().cells[q].z; r = f.mesh().cells[q].r; }
        }
        std::printf(" L%d u_z %.1f at (%.2f, %.2f) mm;", k, best, z * 1e3, r * 1e3);
      }
      std::printf("\n");
    }
  }
  const auto& l = spray.ledger();
  // The end state near the corner: the owner cell's u_z and T along finest rings 0 and 4, and the parcels per 4 finest columns.
  const Mesh& fm = amr.level(amr.levels() - 1).mesh();
  const int f = amr.levels() - 1;
  for (int j : {0, 4}) {
    std::printf("  ring %d (r %.2f mm):", j, fm.cells[fm.index(0, j)].r * 1e3);
    for (int i = 0; i < 24; ++i) {
      int lv = 0;
      std::size_t q = amr.level(0).mesh().index(i >> f, j >> f);
      for (int m = amr.top(); m >= 1; --m) {
        const auto qm = amr.level(m).mesh().index(i >> (f - m), j >> (f - m));
        if (amr.active(m)[qm]) { lv = m; q = qm; break; }
      }
      const auto& u = amr.level(lv).state()[q];
      std::printf(" %d:%.1f", lv, u[1] / u[0]);
    }
    std::printf("  (level:u_z m/s, finest columns 0 to 23, dz %.3f mm)\n", fm.dz * 1e3);
  }
  std::vector<int> bins(16, 0);
  double dsum = 0;
  for (const auto& p : spray.parcels()) { bins[std::min(15, static_cast<int>(p.z / fm.dz) / 4)]++; dsum += p.d; }
  std::printf("  parcels per 4 finest columns:");
  for (int b : bins) std::printf(" %d", b);
  std::printf("; mean d %.1f um\n", spray.parcels().empty() ? 0 : dsum / spray.parcels().size() * 1e6);
  std::printf("mode %d: %llu base steps, %llu regrids, %.1f s; evaporated %.6e kg, %zu parcels held, %ld injected, %ld children, "
              "%ld sub-steps; level moves %ld finer, %ld coarser\n",
              given, static_cast<unsigned long long>(amr.steps()), static_cast<unsigned long long>(amr.regrids()), seconds(start),
              l.exchangedMass, spray.parcels().size(), l.parcelsInjected, l.children, l.substeps, l.finerMoves, l.coarserMoves);
}
