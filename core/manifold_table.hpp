#pragma once
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace crucible {

// Table B's data layout: a tabulated finite-rate progress manifold (TECHNICAL_PLAN, *Lightweight
// engine*; docs/evidence/TABLE_B.md). This holds the layout only, built offline by
// crucible_build_manifold_table from Cantera and CVODES on the full mechanism; no lookup yet (that
// is stated with the runtime, before its own code). No Cantera here.
//
// Each node of (f_H, Z_N, eta, ln rho) is a constant-(e, rho) homogeneous reactor, built as two
// branches meeting at equilibrium (TABLE_B.md, "The manifold"): a lower branch approaching
// equilibrium from below in the progress variable c, and an upper branch from above. The two are
// joined into one curve of 40 stored points per node (TABLE_B.md, "The builder, stated before its
// code"): the lower branch's own start (1 point), 12 points geometric in s = (c - c_0)/(c_eq - c_0)
// from 1e-9 to 0.1, 18 uniform in s from 0.1 to 1 (c_eq); then the upper branch's 9 points uniform
// from c_eq to c_cp (complete combustion, no dissociation). A frozen branch (one that has not
// reached equilibrium by 1e15 s) places the same s fractions on the span it actually reached.
//
// Node elements: Y_H = f_H (1 - Z_N), Y_O = (1 - f_H)(1 - Z_N), Y_N2 = Z_N. The reacting species
// (made of H and O alone) are stored per point; N2 and Ar are fixed by their elements (Z_N and 0)
// and are not stored. e is normalised per (f_H, Z_N) node as e_lo + eta (e_hi - e_lo), e_lo the
// complete-combustion products' energy at 200 K, e_hi the same products frozen at 6000 K.
struct ManifoldTable {
  // Axes. f_H and Z_N are explicit node lists (clustered at stoichiometric and the ends, as Table A's
  // f); eta is uniform on [0, 1]; ln rho is uniform over [lnRhoMin, lnRhoMax].
  std::vector<double> fH;  // increasing, 0 to 1 (both endpoints are nodes: the degenerate ends)
  std::vector<double> zN;  // increasing, the given Z_N list
  int nEta{};
  int nLnRho{};
  double lnRhoMin{}, lnRhoMax{};

  // The reacting (H/O) species' names, in the order their Y is stored per point.
  std::vector<std::string> reactingSpecies;

  // The curve's point layout (TABLE_B.md, "The builder"): 1 + lowerGeometric + lowerUniform points
  // on the lower branch (the branch's own start, then geometric, then uniform in s), then
  // upperUniform points on the upper branch, uniform in c from c_eq to c_cp.
  int lowerGeometric{12}, lowerUniform{18}, upperUniform{9};
  [[nodiscard]] int lowerPoints() const { return 1 + lowerGeometric + lowerUniform; }
  [[nodiscard]] int pointsPerNode() const { return lowerPoints() + upperUniform; }

  enum Flag : std::uint32_t {
    kFrozen = 1u << 0,      // a branch had not reached equilibrium by 1e15 s
    kDegenerate = 1u << 1,  // f_H = 0 or 1: c cannot describe the node (TABLE_B.md, "Nodes c cannot describe")
    kZeroSpan = 1u << 2,    // a branch's |c_eq - c_start| < 1e-12
    kOutside = 1u << 3,     // a point on the node's branches is above the mechanism's thermo range
  };

  // Per node (f_H, Z_N, eta, ln rho): the branch endpoints and the gates' bookkeeping.
  struct Node {
    double c0{}, cEq{}, cCp{};
    std::uint32_t flags{};
    std::int32_t progressVariable{};  // the candidate used (0: Y_H2O, 1: Y_H2O+HO2+H2O2, 2: O2 consumed)
  };
  std::vector<Node> nodes;  // size fH.size() * zN.size() * nEta * nLnRho, ln rho fastest

  // The per-(f_H, Z_N) energy normalisation: complete-combustion products at 200 K (e_lo) and the
  // same products frozen at 6000 K (e_hi). Size fH.size() * zN.size().
  std::vector<double> eLo, eHi;

  // Per point: the reacting species' Y (reactingSpecies.size() values), then lambda =
  // (dc/dt)/(c_eq - c). Flat: node index * pointsPerNode() * (reactingSpecies.size() + 1) + point *
  // (reactingSpecies.size() + 1) + [0 .. size-1] = Y, [size] = lambda.
  std::vector<double> points;

  std::string provenance;

  [[nodiscard]] std::size_t fzIndex(std::size_t i, std::size_t j) const { return i * zN.size() + j; }
  [[nodiscard]] std::size_t nodeIndex(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
    return ((i * zN.size() + j) * static_cast<std::size_t>(nEta) + k) * static_cast<std::size_t>(nLnRho) + l;
  }
  [[nodiscard]] double eta(std::size_t k) const { return nEta > 1 ? static_cast<double>(k) / (nEta - 1) : 0.0; }
  [[nodiscard]] double lnRho(std::size_t l) const {
    return nLnRho > 1 ? lnRhoMin + (lnRhoMax - lnRhoMin) * static_cast<double>(l) / (nLnRho - 1) : lnRhoMin;
  }
  [[nodiscard]] std::size_t pointWidth() const { return reactingSpecies.size() + 1; }
  [[nodiscard]] std::size_t pointOffset(std::size_t node, int point) const {
    return (node * static_cast<std::size_t>(pointsPerNode()) + static_cast<std::size_t>(point)) * pointWidth();
  }
  [[nodiscard]] double* pointY(std::size_t node, int point) { return points.data() + pointOffset(node, point); }
  [[nodiscard]] const double* pointY(std::size_t node, int point) const { return points.data() + pointOffset(node, point); }
  [[nodiscard]] double& lambda(std::size_t node, int point) { return points[pointOffset(node, point) + reactingSpecies.size()]; }
  [[nodiscard]] double lambda(std::size_t node, int point) const { return points[pointOffset(node, point) + reactingSpecies.size()]; }

  void save(const std::string& path) const;
  static ManifoldTable load(const std::string& path);
};

}  // namespace crucible
