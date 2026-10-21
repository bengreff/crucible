#pragma once
// The RL10A-3-3A thrust chamber's wall (docs/validation/RL10A33A_PROTOCOL.md, sections 1 and 2), shared by
// chamber_study's rl10dev and the tier-1 drivers (factored from rl10dev, 21 October 2026, its output unchanged).
// r = r_t sqrt(A / A_t) through TM-107318 Table E1's 16 nodes, the throat (0, r_t), a face node at the chamber
// radius r_t 5.13 / 2.47 and the exit node (Figure E1's 46.0 in, r_t sqrt 61), by PCHIP (core/contour). The
// injector face is the one nearest Figure E1's -13.1 in that puts a mesh station on the throat at uniform axial
// spacing; it must fall in the protocol's -13.4 to -12.8 in bracket. r_t is CRUCIBLE_RL10_THROAT_RADIUS in
// metres, default 0.062738 (the radius reading, protocol section 1).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/contour.hpp"

namespace crucible::rl10 {

constexpr double kInch = 0.0254, kExitStation = 46.0, kNominalFace = -13.1;

struct Wall {
  double rt{}, rc{};
  double face{};               // the injector face's station [in], throat at 0
  int throatStation{};         // the throat's mesh station, of nz
  std::vector<std::array<double, 2>> nodes;  // sorted, z [m] from the throat
  std::size_t throatNode{};
  Pchip pchip;
  std::vector<std::array<double, 2>> contour;  // (z from the face, r) [m], 200 samples per interval
};

inline Wall wall(int nz) {
  double rt = 0.062738;
  if (const char* v = std::getenv("CRUCIBLE_RL10_THROAT_RADIUS")) rt = std::atof(v);
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
  const int k = static_cast<int>(std::lround(nz * -kNominalFace / (kExitStation - kNominalFace)));
  const double face = -kExitStation * k / (nz - k);  // inches
  if (face > -12.8 || face < -13.4)
    throw std::invalid_argument("rl10dev: the throat station puts the face at " + std::to_string(face) +
                                " in, outside -13.4 to -12.8 in; choose another nz");
  const double rc = rt * 5.13 / 2.47;
  nodes.push_back({0.0, rt});
  nodes.push_back({face * kInch, rc});
  nodes.push_back({kExitStation * kInch, rt * std::sqrt(61.0)});
  std::sort(nodes.begin(), nodes.end());
  Pchip pchip(nodes);
  std::size_t throatNode = 0;
  while (nodes[throatNode][0] != 0) ++throatNode;
  std::vector<std::array<double, 2>> contour;
  for (const auto& p : pchip.sample(200)) contour.push_back({p[0] - face * kInch, p[1]});
  return {rt, rc, face, k, std::move(nodes), throatNode, std::move(pchip), std::move(contour)};
}

}  // namespace crucible::rl10
