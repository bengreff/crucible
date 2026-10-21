#pragma once
// The RL10 development wall, shared by the AMR tests.
#include "core/contour.hpp"
#include "core/flow.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace crucible {

// The RL10 development wall as tests/chamber_study.cpp rl10dev builds it, with the injector face placed
// for a base of `baseColumns` and the mesh at nz x nr.
inline Definition rl10Wall(int baseColumns, int nz, int nr, double stretching) {
  constexpr double inch = 0.0254, exitStation = 46.0, nominalFace = -13.1, rt = 0.062738;
  std::vector<std::array<double, 2>> nodes;
  const std::string path = std::string(CRUCIBLE_SOURCE_DIR) + "/docs/validation/rl10a33a/tm107318_table_E1_si.csv";
  FILE* f = std::fopen(path.c_str(), "r");
  if (!f) throw std::runtime_error("cannot read " + path);
  char line[512];
  bool header = true;
  while (std::fgets(line, sizeof line, f)) {
    if (line[0] == '#') continue;
    if (header) { header = false; continue; }
    double zIn, zM, ratio;
    if (std::sscanf(line, "%lf,%lf,%lf", &zIn, &zM, &ratio) == 3) nodes.push_back({zM, rt * std::sqrt(ratio)});
  }
  std::fclose(f);
  if (nodes.size() != 16) throw std::runtime_error("Table E1: expected 16 nodes in " + path);
  const int k = static_cast<int>(std::lround(baseColumns * -nominalFace / (exitStation - nominalFace)));
  const double face = -exitStation * k / (baseColumns - k);
  nodes.push_back({0.0, rt});
  nodes.push_back({face * inch, rt * 5.13 / 2.47});
  nodes.push_back({exitStation * inch, rt * std::sqrt(61.0)});
  std::sort(nodes.begin(), nodes.end());
  Definition d;
  d.experiment = Case::UniformDuct;
  for (const auto& p : Pchip(nodes).sample(200)) d.contour.push_back({p[0] - face * inch, p[1]});
  d.nz = nz;
  d.nr = nr;
  d.radialStretching = stretching;
  return d;
}
}  // namespace crucible
