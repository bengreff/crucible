#pragma once
// AMR step B2's comparison dumps (docs/evidence/AMR.md, Step B2: setup): per level its active cells (every
// cell of the base) as the four conserved values, the partial densities, and rho k and rho omega with
// turbulence, then the base step count and the composite budget, raw, for cmp between B1 and B2.
#include "core/amr.hpp"
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace crucible {
inline void writeAmrDump(const Amr& amr, const std::string& path) {
  std::FILE* out = std::fopen(path.c_str(), "wb");
  if (!out) throw std::runtime_error("Cannot write " + path);
  for (int l = 0; l < amr.levels(); ++l) {
    const Flow& f = amr.level(l);
    const std::size_t cells = f.state().size();
    const std::size_t ns = f.partialDensities().size() / cells, nt = f.turbulence().size() / cells;
    for (std::size_t q = 0; q < cells; ++q) {
      if (l > 0 && !amr.active(l)[q]) continue;
      std::fwrite(f.state()[q].data(), sizeof(double), 4, out);
      std::fwrite(f.partialDensities().data() + q * ns, sizeof(double), ns, out);
      if (nt) std::fwrite(f.turbulence().data() + q * nt, sizeof(double), nt, out);
    }
  }
  const std::uint64_t steps = amr.steps();
  std::fwrite(&steps, sizeof steps, 1, out);
  const Amr::Budget b = amr.budget();
  const double budget[] = {b.mass, b.energy, b.momentum, b.massResidual, b.energyResidual, b.momentumResidual, b.momentumGross, b.igniterEnergy};
  std::fwrite(budget, sizeof(double), 8, out);
  std::fclose(out);
}
}  // namespace crucible
