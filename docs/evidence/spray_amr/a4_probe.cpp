// The mixed-level jet's diagnosis (SPRAY.md, A defect the printed comparison found): S3's chamber, frozen, 1 thread.
// usage: a4_probe KIND [MASK] [END_MS]
//   KIND u0, u1, u2: uniform 16 x 6, 32 x 12, 64 x 24 (Flow::step, then Spray::advance)
//        cell: Amr(def, 1), level 1 static over base cell (2, 1) alone (z 2 to 3 mm, r 1.67 to 3.33 mm)
//        block: Amr(def, 1), level 1 static over base cells i 0 to 2, j 0 to 1 (z < 3 mm, r < 3.33 mm)
//        axis: Amr(def, 1), level 1 static over base cell (2, 0) alone (z 2 to 3 mm, r 0 to 1.67 mm)
//        all1: Amr(def, 1), level 1 active everywhere
//   A KIND ending in "+o2" has no injectors and a second supply of gaseous O2 at 300 K, 2e-3 kg/s over r < 1.5 mm:
//   the gas alone with a composition front across the patch.
//        corner2: a2_levels_probe mode 1 (Amr(def, 2), the corner tag, regridded every 2 base steps)
//   MASK: sprayDevelopmentMask (bits 1, 2, 4: the mass, momentum and energy deposits kept; 8: a deposit to a cell
//         that is not owned throws). Default 15.
// Every 0.02 ms it prints each level's largest owned |u_z| and where. Printed, not judged.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <numbers>
#include <string>
#include <vector>
#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"
#include "core/amr.hpp"
#include "core/spray.hpp"
namespace crucible { extern int sprayDevelopmentMask; extern int amrDevelopmentSpread; }
using namespace crucible;
#include "common.inc"
}
namespace {
void trace(const std::function<int()>& top, const std::function<const Flow&(int)>& level,
           const std::function<bool(int, std::size_t)>& owned, double t, std::size_t parcels, double evaporated) {
  std::printf("  t %.3f ms, %zu parcels, evaporated %.4e kg:", t * 1e3, parcels, evaporated);
  for (int k = 0; k <= top(); ++k) {
    const Flow& f = level(k);
    double best = 0, z = 0, r = 0;
    for (std::size_t q = 0; q < f.state().size(); ++q) {
      if (!owned(k, q)) continue;
      const double u = f.state()[q][1] / f.state()[q][0];
      if (std::abs(u) > std::abs(best)) { best = u; z = f.mesh().cells[q].z; r = f.mesh().cells[q].r; }
    }
    std::printf(" L%d u_z %.1f at (%.2f, %.2f) mm;", k, best, z * 1e3, r * 1e3);
  }
  std::printf("\n");
}
}  // namespace
int main(int argc, char** argv) {
  std::string kind = argc > 1 ? argv[1] : "cell";
  const bool o2 = kind.size() > 3 && kind.substr(kind.size() - 3) == "+o2";
  if (o2) kind.resize(kind.size() - 3);
  auto define = [o2](thermo::ReactionSource& source) {
    Definition d = chamber(source, false);
    if (o2) {
      Supply s = d.supplies[0];
      s.innerRadius = 0;
      s.outerRadius = 1.5e-3;
      s.massFlow = 2e-3;
      s.composition = source.massFractions("O2:1");
      d.supplies.push_back(s);
    }
    return d;
  };
  auto spraySettings = [o2] {
    auto s = settings();
    if (o2) s.injectors.clear();
    return s;
  };
  sprayDevelopmentMask = argc > 2 ? std::atoi(argv[2]) : 15;
  amrDevelopmentSpread = argc > 4 ? std::atoi(argv[4]) : 0;
  const double end = (argc > 3 ? std::atof(argv[3]) : 0.2) * 1e-3;
  thermo::ReactionSource source("h2o2.yaml");
  const auto start = std::chrono::steady_clock::now();
  double next = 0.02e-3;
  if (kind[0] == 'u') {
    const int shift = kind[1] - '0';
    Definition d = define(source);
    d.nz <<= shift;
    d.nr <<= shift;
    Flow flow(d);
    Spray spray(flow, oxygen(flow), spraySettings());
    while (flow.time() < end) {
      const double t = flow.time();
      spray.advance(flow.step(std::min(end - t, spray.nextEvent(t) - t)));
      if (flow.time() >= next * (1 - 1e-9)) {
        trace([] { return 0; }, [&](int) -> const Flow& { return flow; }, [](int, std::size_t) { return true; }, flow.time(),
              spray.parcels().size(), spray.ledger().exchangedMass);
        next += 0.02e-3;
      }
    }
    const auto& l = spray.ledger();
    std::printf("%s%s mask %d: %d x %d, %.1f s; evaporated %.6e kg, %zu parcels held, %ld children\n", kind.c_str(), o2 ? "+o2" : "", sprayDevelopmentMask,
                d.nz, d.nr, seconds(start), l.exchangedMass, spray.parcels().size(), l.children);
    return 0;
  }
  const Definition d = define(source);
  const int finer = kind == "corner2" ? 2 : 1;
  Amr amr(d, finer);
  amr.setThreads(1);
  if (kind == "corner2") {
    Amr::Tags tags;
    tags.interval = 2;
    tags.throatHalfWidth = 0;
    tags.wallYPlus = 0;
    tags.mixture = tags.temperature = tags.pressure = tags.vorticity = 0;
    tags.custom = [&amr](int l, int i, int j) {
      const Mesh& m = amr.level(l).mesh();
      return (i + 0.5) * m.dz < 2e-3 && m.cells[m.index(i, j)].r < 3e-3;
    };
    amr.setTags(tags);
  } else {
    const Mesh& fm = amr.level(1).mesh();
    std::vector<unsigned char> active(fm.cells.size(), 0);
    for (int i = 0; i < fm.nz; ++i)
      for (int j = 0; j < fm.nr; ++j) {
        const int I = i >> 1, J = j >> 1;
        const bool on = kind == "all1" || (kind == "cell" && I == 2 && J == 1) || (kind == "axis" && I == 2 && J == 0) || (kind == "block" && I <= 2 && J <= 1);
        active[fm.index(i, j)] = on;
      }
    amr.setActive(1, active);
  }
  Spray spray(amr, oxygen(amr.level(0)), spraySettings());
  amr.setAfterStep([&](int l, double dt) { spray.afterStep(l, dt); });
  amr.setRedo([&] { spray.redo(); });
  auto owned = [&amr](int k, std::size_t q) {
    const Flow& f = amr.level(k);
    if (k > 0 && !amr.active(k)[q]) return false;
    return !(k < amr.top() && amr.active(k + 1)[amr.level(k + 1).mesh().index(2 * static_cast<int>(q / f.mesh().nr), 2 * static_cast<int>(q % f.mesh().nr))]);
  };
  while (amr.time() < end) {
    const double t = amr.time();
    amr.step(std::min(end - t, spray.nextEvent(t) - t));
    if (amr.time() >= next * (1 - 1e-9)) {
      trace([&] { return amr.top(); }, [&](int k) -> const Flow& { return amr.level(k); }, owned, amr.time(), spray.parcels().size(),
            spray.ledger().exchangedMass);
      next += 0.02e-3;
    }
  }
  const auto& l = spray.ledger();
  std::printf("%s%s mask %d: %llu base steps, %llu regrids, %.1f s; evaporated %.6e kg, %zu parcels held, %ld children; moves %ld finer, %ld coarser\n",
              kind.c_str(), o2 ? "+o2" : "", sprayDevelopmentMask, static_cast<unsigned long long>(amr.steps()), static_cast<unsigned long long>(amr.regrids()),
              seconds(start), l.exchangedMass, spray.parcels().size(), l.children, l.finerMoves, l.coarserMoves);
}
