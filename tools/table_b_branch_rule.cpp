// E0's replay (docs/evidence/TABLE_B.md, "The entropy candidate in the builder and the E0 replay", and the
// runtime's xi in "per-node refinement"): the branch rule of a Table B built with the specific entropy,
// applied to a chamber study's state samples (CRUCIBLE_STATE_SAMPLES with T and Y). Per sample:
//   - s_cell from the cell's T, rho and Y on the mechanism's NASA-7 data (core Species), with the table's p°;
//     above a species' top temperature, s° continues as cp ln(T / T_top) (the constant-cp continuation);
//   - per branch, the 16 corners' blended s_eq and span (multilinear weights in f_H, Z_N, eta, ln rho), and
//     xi = 1 - sqrt((s_eq - s_cell) / span), exact for the blend since every stored point satisfies
//     s = s_ref - span (1 - xi)^2 on its own corner;
//   - each corner's Y at its own s at that xi, by a search of its points, linear in d = sqrt(s_ref - s);
//     the blend of the 16, and max_k |dY_k| against the cell's reacting species;
//   - the nearer branch by that max|dY|; CVODES, by cause, if s_cell is on neither curve, on an exact tie,
//     or if the nearer branch's max|dY| is above dY_max = 0.01.
// And the rule that replaces it (TABLE_B.md, "The cell placed by its moles per kilogram"): N = sum Y_k / W_k
// of the cell against the blended curve's N, tau from 0 (start) to 1 (equilibrium or reach) on the lower
// branch and 1 to 2 (complete combustion) on the upper; phi_N, the 02:59 weight w, the branch from N, tau_N by
// Illinois, s* = (1 - w) s_cell + w s(tau_N) for w < 1, and the same dY_max on the placed point's Y.
// One output row per sample, in the samples' order, for tools/fallback_rate.py to join (it applies the
// domain: the axes and the gates; E0_COLUMN picks the rule). Also E4: s per cell and each rule per cell,
// timed on one thread.
// R1 and R2 (TABLE_B.md, "The runtime's step on the placed point"): the runtime's own placement
// (core/manifold_runtime.*) on every sample, against this replay's (columns n_runtime, tau_runtime); and its step at
// the cell on every CRUCIBLE_R2_EVERY-th sample it places (default 50).
// CRUCIBLE_E0_OWN_EQ_EVERY = N (default 20; 0 for none): every N-th sample inside the axes also gets its own
// s_eq from Cantera's equilibrateUV, for E0's stated measurement of the blend's s_eq against the cell's.
//
// Usage: crucible_table_b_branch_rule <table> <mechanism> <samples CSV> <out CSV>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "adapters/reaction.hpp"
#include "cantera/core.h"
#include "core/manifold_runtime.hpp"
#include "core/manifold_table.hpp"
#include "core/medium.hpp"

using crucible::ManifoldTable;
using crucible::Medium;
using crucible::Species;

namespace {
constexpr double kDYMax = 0.01;

// s°/R of one species (NASA-7), with the constant-cp continuation above tHigh.
double entropyOverR(const Species& s, double t) {
  const double c = std::min(t, s.tHigh);
  const auto& a = c < s.tMid ? s.low : s.high;
  double v = a[0] * std::log(c) + c * (a[1] + c * (a[2] / 2 + c * (a[3] / 3 + c * a[4] / 4))) + a[6];
  if (t > s.tHigh) v += (a[0] + c * (a[1] + c * (a[2] + c * (a[3] + c * a[4])))) * std::log(t / s.tHigh);
  return v;
}

// s [J/(kg K)] = sum_k Y_k / W_k (R s°_k - R ln(X_k p / p°)), with X_k p = rho Y_k R T / W_k.
double cellEntropy(const std::vector<Species>& sp, double t, double rho, const double* y, double pRef) {
  constexpr double r = Medium::universalGasConstant;
  double s = 0;
  for (std::size_t k = 0; k < sp.size(); ++k) {
    if (!(y[k] > 0)) continue;
    const double n = y[k] / sp[k].molarMass;
    s += n * r * (entropyOverR(sp[k], t) - std::log(rho * n * r * t / pRef));
  }
  return s;
}

struct Corner { std::size_t node; double w; };

// One branch of one corner: its stored points [first, first + count), c monotone, s_ref the end at
// equilibrium (or a frozen branch's reach) and span |s_ref - s_far|.
struct Branch { std::size_t first, count; bool ascending; double sRef, span; };

Branch branchOf(const ManifoldTable& t, std::size_t node, int b) {
  const auto& n = t.nodes[node];
  Branch br{};
  if (b == 0) {
    br.first = 0; br.count = n.lowerCount; br.ascending = true;
    br.sRef = t.c(node, static_cast<int>(n.lowerCount) - 1);
    br.span = br.sRef - t.c(node, 0);
  } else {
    br.first = n.lowerCount; br.count = n.pointCount - n.lowerCount; br.ascending = false;
    br.sRef = t.c(node, static_cast<int>(n.lowerCount));
    br.span = br.sRef - t.c(node, static_cast<int>(n.pointCount) - 1);
  }
  return br;
}

// The bracket of u = 1 - xi on branch br of one corner: its s is s_ref - span u^2, bracketed by a search in c
// (monotone), and interpolated linearly in d = sqrt(s_ref - s) between the bracket's points a and b.
struct Bracket { int a, b; double th; };
Bracket cornerBracket(const ManifoldTable& t, std::size_t node, const Branch& br, double u) {
  const int p0 = static_cast<int>(br.first);
  if (br.span <= 0 || br.count < 2) return {p0, p0, 0.0};
  const double sTarget = br.sRef - br.span * u * u;
  // Points ordered by increasing d: the lower branch is stored from its start (largest d) to s_ref.
  auto at = [&](std::size_t q) { return br.ascending ? p0 + static_cast<int>(br.count - 1 - q) : p0 + static_cast<int>(q); };
  // In d order, s falls; find q with s(q) >= sTarget >= s(q + 1).
  std::size_t lo = 0, hi = br.count - 1;
  while (hi - lo > 1) {
    const std::size_t mid = (lo + hi) / 2;
    (t.c(node, at(mid)) >= sTarget ? lo : hi) = mid;
  }
  const double dLo = std::sqrt(std::max(br.sRef - t.c(node, at(lo)), 0.0)), dHi = std::sqrt(std::max(br.sRef - t.c(node, at(hi)), 0.0));
  const double d = std::sqrt(br.span) * u;
  return {at(lo), at(hi), dHi > dLo ? std::clamp((d - dLo) / (dHi - dLo), 0.0, 1.0) : 0.0};
}

// The corner's Y (reacting species) at u on branch br.
void cornerY(const ManifoldTable& t, std::size_t node, const Branch& br, double u, double* out) {
  const std::size_t ns = t.reactingSpecies.size();
  const Bracket k = cornerBracket(t, node, br, u);
  const double* a = t.pointY(node, k.a);
  const double* b = t.pointY(node, k.b);
  for (std::size_t q = 0; q < ns; ++q) out[q] = a[q] + k.th * (b[q] - a[q]);
}

struct BranchResult { bool on{}; double sEq{}, span{}, xi{}, dY{}; const char* off{""}; };
struct Rule { BranchResult b[2]; int branch{-1}; std::string cause; };

// Brackets x on a sorted axis; false if outside.
bool bracket(const std::vector<double>& a, double x, std::size_t& i, double& f) {
  if (!(x >= a.front() && x <= a.back())) return false;
  i = std::min<std::size_t>(std::upper_bound(a.begin(), a.end(), x) - a.begin(), a.size() - 1);
  i = i > 0 ? i - 1 : 0;
  if (i + 1 >= a.size()) i = a.size() - 2;
  f = (x - a[i]) / (a[i + 1] - a[i]);
  return true;
}
bool bracketUniform(double lo, double hi, int n, double x, std::size_t& i, double& f) {
  if (!(x >= lo - 1e-12 * std::fabs(hi - lo) && x <= hi + 1e-12 * std::fabs(hi - lo))) return false;
  const double g = std::clamp((x - lo) / (hi - lo) * (n - 1), 0.0, static_cast<double>(n - 1));
  i = std::min<std::size_t>(static_cast<std::size_t>(g), static_cast<std::size_t>(n - 2));
  f = g - static_cast<double>(i);
  return true;
}

// The 16 corners and their weights, or a cause if the cell is outside the axes or a corner is unusable.
std::string cornersOf(const ManifoldTable& t, double fH, double zN, double eta, double rho, std::array<Corner, 16>& c) {
  std::size_t ii, jj, kk, ll;
  double fi, fj, fk, fl;
  if (!std::isfinite(fH)) return "no H or O";
  if (!bracket(t.zN, zN, jj, fj)) return "outside the axes";
  if (!bracket(t.fH, fH, ii, fi)) return "outside the axes";
  if (!bracketUniform(0.0, 1.0, t.nEta, eta, kk, fk)) return "outside the axes";
  if (!bracketUniform(t.lnRhoMin, t.lnRhoMax, t.nLnRho, std::log(rho), ll, fl)) return "outside the axes";
  for (int q = 0; q < 16; ++q) {
    const int a = q >> 3 & 1, b = q >> 2 & 1, e = q >> 1 & 1, r = q & 1;
    const double w = (a ? fi : 1 - fi) * (b ? fj : 1 - fj) * (e ? fk : 1 - fk) * (r ? fl : 1 - fl);
    c[q] = {t.nodeIndex(ii + a, jj + b, kk + e, ll + r), w};
    const auto& n = t.nodes[c[q].node];
    if (n.pointCount == 0 || n.lowerCount == 0 || n.lowerCount >= n.pointCount) return "a corner was refused";
    if (n.flags & (ManifoldTable::kOutside | ManifoldTable::kDegenerate)) return "a corner flagged outside or degenerate";
  }
  return "";
}

Rule applyRule(const ManifoldTable& t, const std::array<Corner, 16>& corners, double sCell, const double* yCell,
               std::vector<double>& yc, std::vector<double>& yb) {
  const std::size_t ns = t.reactingSpecies.size();
  Rule rule;
  for (int b = 0; b < 2; ++b) {
    auto& r = rule.b[b];
    std::array<Branch, 16> br;
    for (int q = 0; q < 16; ++q) {
      br[q] = branchOf(t, corners[q].node, b);
      r.sEq += corners[q].w * br[q].sRef;
      r.span += corners[q].w * br[q].span;
    }
    if (sCell > r.sEq) { r.off = "above the blended s_eq"; continue; }
    if (!(r.span > 0) || sCell < r.sEq - r.span) { r.off = "below the blended start"; continue; }
    const double u = std::sqrt((r.sEq - sCell) / r.span);
    r.xi = 1 - u;
    std::fill(yb.begin(), yb.end(), 0.0);
    for (int q = 0; q < 16; ++q) {
      if (corners[q].w == 0) continue;
      cornerY(t, corners[q].node, br[q], u, yc.data());
      for (std::size_t k = 0; k < ns; ++k) yb[k] += corners[q].w * yc[k];
    }
    for (std::size_t k = 0; k < ns; ++k) r.dY = std::max(r.dY, std::fabs(yb[k] - yCell[k]));
    r.on = true;
  }
  if (!rule.b[0].on && !rule.b[1].on) {
    rule.cause = std::string("off the curve (lower ") + rule.b[0].off + "; upper " + rule.b[1].off + ")";
    return rule;
  }
  if (rule.b[0].on && rule.b[1].on && rule.b[0].dY == rule.b[1].dY) { rule.cause = "a tie"; return rule; }
  const int best = !rule.b[1].on || (rule.b[0].on && rule.b[0].dY < rule.b[1].dY) ? 0 : 1;
  if (rule.b[best].dY > kDYMax) { rule.cause = "dY above dY_max"; return rule; }
  rule.branch = best;
  rule.cause = best == 0 ? "table (lower branch)" : "table (upper branch)";
  return rule;
}

// The N placement (TABLE_B.md, "The cell placed by its moles per kilogram"). nPt: N of every stored point's
// reacting species, by global point index; N is linear in Y, so a corner's N at u is the bracket's
// interpolation of it, and the blend's N is the weights' sum.
constexpr double kNTol = 1e-12;
struct NRule { std::string cause; int branch{-1}; double phi{NAN}, w{NAN}, tau{NAN}, dY{NAN}, gap{NAN}; int iterations{}; };

struct Blend {
  const ManifoldTable& t;
  const std::vector<double>& nPt;
  const std::array<Corner, 16>& corners;
  std::array<Branch, 16> br[2];
  double sEq[2]{}, span[2]{};
  Blend(const ManifoldTable& tab, const std::vector<double>& n, const std::array<Corner, 16>& c) : t(tab), nPt(n), corners(c) {
    for (int b = 0; b < 2; ++b)
      for (int q = 0; q < 16; ++q) {
        br[b][q] = branchOf(t, corners[q].node, b);
        sEq[b] += corners[q].w * br[b][q].sRef;
        span[b] += corners[q].w * br[b][q].span;
      }
  }
  // u on branch b: the lower branch's tau is 1 - u, the upper's 1 + u.
  [[nodiscard]] double n(int b, double u) const {
    double v = 0;
    for (int q = 0; q < 16; ++q) {
      if (corners[q].w == 0) continue;
      const Bracket k = cornerBracket(t, corners[q].node, br[b][q], u);
      const std::size_t f = t.nodes[corners[q].node].firstPoint;
      v += corners[q].w * (nPt[f + k.a] + k.th * (nPt[f + k.b] - nPt[f + k.a]));
    }
    return v;
  }
  // N at a branch's end, read from its stored end point: u = 0 is s_ref, u = 1 the far end.
  [[nodiscard]] double nEnd(int b, int u) const {
    double v = 0;
    for (int q = 0; q < 16; ++q) {
      const auto& n = t.nodes[corners[q].node];
      const std::uint32_t p = b == 0 ? (u == 1 ? 0 : n.lowerCount - 1) : (u == 0 ? n.lowerCount : n.pointCount - 1);
      v += corners[q].w * nPt[n.firstPoint + p];
    }
    return v;
  }
  void y(int b, double u, std::vector<double>& yc, std::vector<double>& yb) const {
    std::fill(yb.begin(), yb.end(), 0.0);
    for (int q = 0; q < 16; ++q) {
      if (corners[q].w == 0) continue;
      cornerY(t, corners[q].node, br[b][q], u, yc.data());
      for (std::size_t k = 0; k < yb.size(); ++k) yb[k] += corners[q].w * yc[k];
    }
  }
};

// The base tau values (the builder's slots): lower, xi = 0, 12 geometric 1e-9..0.1, 18 uniform to 1; upper,
// 1 + j/8 for its 9 points.
std::vector<double> baseTau(const ManifoldTable& t, int b) {
  std::vector<double> v;
  if (b == 0) {
    v.push_back(0.0);
    for (int g = 1; g <= t.lowerGeometric; ++g)
      v.push_back(1e-9 * std::pow(0.1 / 1e-9, t.lowerGeometric > 1 ? static_cast<double>(g - 1) / (t.lowerGeometric - 1) : 1.0));
    for (int u = 1; u <= t.lowerUniform; ++u) v.push_back(0.1 + 0.9 * static_cast<double>(u) / t.lowerUniform);
  } else {
    for (int j = 0; j < t.upperUniform; ++j) v.push_back(1.0 + static_cast<double>(j) / (t.upperUniform - 1));
  }
  return v;
}

// check: sample the blend's N at the base tau values, report a rise as its own cause and bracket tau_N there
// (the replay); without it, Illinois runs on the whole branch (the runtime's cost, for E4).
NRule placeByN(const Blend& bl, const std::vector<double>& wR, double sCell, const double* yCell, bool check,
               const std::vector<double>* tauBase, std::vector<double>& yc, std::vector<double>& yb) {
  const std::size_t ns = wR.size();
  NRule r;
  double nCell = 0;
  for (std::size_t k = 0; k < ns; ++k) nCell += yCell[k] / wR[k];
  const double n0 = bl.nEnd(0, 1), nLoEnd = bl.nEnd(0, 0), nUpEnd = bl.nEnd(1, 0), nCp = bl.nEnd(1, 1);
  r.gap = (nLoEnd - nUpEnd) / nLoEnd;
  if (!(n0 > nLoEnd)) { r.cause = "no fall in N"; return r; }
  r.phi = (n0 - nCell) / (n0 - nLoEnd);
  r.w = std::clamp((r.phi - 0.005) / 0.005, 0.0, 1.0);
  const char* wClass = r.w == 0 ? "w 0" : (r.w < 1 ? "w between" : "w 1");
  int b = 0;
  double tauN = NAN;
  if (r.w > 0) {
    if (nCell >= nLoEnd) b = 0;
    else if (nCell <= nUpEnd) b = 1;
    else { r.cause = "frozen reach"; return r; }
    if (b == 1 && nCell < nCp) { r.cause = "N below the blended N_cp"; return r; }
    auto tauOf = [&](double u) { return b == 0 ? 1 - u : 1 + u; };
    auto uOf = [&](double tau) { return b == 0 ? 1 - tau : tau - 1; };
    auto f = [&](double tau) { return bl.n(b, uOf(tau)) - nCell; };
    double a = tauOf(1.0 * (b == 0)), c = tauOf(1.0 * (b == 1)), fa = (b == 0 ? n0 : nUpEnd) - nCell, fc = (b == 0 ? nLoEnd : nCp) - nCell;
    if (check) {
      const auto& g = tauBase[b];
      std::vector<double> ng(g.size());
      for (std::size_t j = 0; j < g.size(); ++j) ng[j] = f(g[j]) + nCell;
      std::size_t start = 0;
      if (b == 0) {
        while (start < g.size() && n0 - ng[start] < 0.005 * (n0 - nLoEnd)) ++start;
        start = start > 0 ? start - 1 : 0;
      }
      for (std::size_t j = start + (b == 0 ? 1 : 0); j + 1 < g.size(); ++j)
        if (ng[j + 1] > ng[j] + kNTol * ng[j]) { r.cause = "blended N not monotone"; return r; }
      std::size_t j = start;
      while (j + 1 < g.size() - 1 && ng[j + 1] > nCell) ++j;
      a = g[j]; c = g[j + 1]; fa = ng[j] - nCell; fc = ng[j + 1] - nCell;
    }
    // Illinois on [a, c], fa >= 0 >= fc.
    const double tol = kNTol * nCell;
    if (std::fabs(fa) <= tol) tauN = a;
    else if (std::fabs(fc) <= tol) tauN = c;
    else {
      for (r.iterations = 1; r.iterations <= 60; ++r.iterations) {
        tauN = c - fc * (c - a) / (fc - fa);
        const double fn = f(tauN);
        if (std::fabs(fn) <= tol) break;
        if (fn * fc < 0) { a = c; fa = fc; } else fa /= 2;
        c = tauN; fc = fn;
      }
    }
  }
  double u;  // the placed point's u on branch b
  if (r.w == 1) {
    u = b == 0 ? 1 - tauN : tauN - 1;
  } else {
    const double sStar = r.w == 0 ? sCell : (1 - r.w) * sCell + r.w * (bl.sEq[0] - bl.span[0] * (1 - tauN) * (1 - tauN));
    if (sStar > bl.sEq[0]) { r.cause = std::string("off the curve in s (above the blended s_eq; ") + wClass + ")"; return r; }
    if (!(bl.span[0] > 0) || sStar < bl.sEq[0] - bl.span[0]) { r.cause = std::string("off the curve in s (below the blended start; ") + wClass + ")"; return r; }
    u = std::sqrt((bl.sEq[0] - sStar) / bl.span[0]);
  }
  r.branch = b;
  r.tau = b == 0 ? 1 - u : 1 + u;
  bl.y(b, u, yc, yb);
  r.dY = 0;
  for (std::size_t k = 0; k < ns; ++k) r.dY = std::max(r.dY, std::fabs(yb[k] - yCell[k]));
  const std::string where = std::string(b == 0 ? "lower; " : "upper; ") + wClass;
  r.cause = r.dY > kDYMax ? "dY above dY_max (" + where + ")" : "table (" + where + ")";
  return r;
}

std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> c;
  std::stringstream s(line);
  std::string f;
  while (std::getline(s, f, ',')) c.push_back(f);
  return c;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc != 5) {
    std::fprintf(stderr, "usage: %s <table> <mechanism> <samples CSV> <out CSV>\n", argv[0]);
    return 2;
  }
  const ManifoldTable table = ManifoldTable::load(argv[1]);
  for (const auto& n : table.nodes)
    if (n.pointCount > 0 && !(n.flags & ManifoldTable::kDegenerate) && n.progressVariable != 3)
      throw std::runtime_error("The table's progress variable is not the specific entropy (candidate 3).");
  crucible::thermo::ReactionSource source(argv[2]);
  const Medium medium = source.medium();
  const auto& species = medium.species();
  const std::size_t nsAll = species.size(), ns = table.reactingSpecies.size();
  auto solution = Cantera::newSolution(argv[2], "", "none");
  auto gas = solution->thermo();
  if (gas->nSpecies() != nsAll) throw std::runtime_error("Mechanism species count differs from the medium's.");
  if (std::fabs(gas->refPressure() - table.referencePressure) > 1e-9 * table.referencePressure)
    throw std::runtime_error("The mechanism's reference pressure differs from the table's.");

  std::ifstream in(argv[3]);
  if (!in) throw std::runtime_error(std::string("Cannot read ") + argv[3]);
  std::string line;
  std::getline(in, line);
  const auto header = split(line);
  std::map<std::string, std::size_t> col;
  for (std::size_t q = 0; q < header.size(); ++q) col[header[q]] = q;
  for (const char* need : {"f_H", "Z_N", "eta", "rho", "T"})
    if (!col.count(need)) throw std::runtime_error(std::string("The samples have no column ") + need);
  std::vector<std::size_t> yCol(nsAll), reactingCol(ns);
  for (std::size_t k = 0; k < nsAll; ++k) {
    const auto it = col.find("Y_" + species[k].name);
    if (it == col.end()) throw std::runtime_error("The samples have no column Y_" + species[k].name);
    yCol[k] = it->second;
  }
  for (std::size_t k = 0; k < ns; ++k) {
    std::size_t q = 0;
    while (q < nsAll && species[q].name != table.reactingSpecies[k]) ++q;
    if (q == nsAll) throw std::runtime_error("Reacting species " + table.reactingSpecies[k] + " is not in the mechanism.");
    reactingCol[k] = q;
  }
  // N per stored point (reacting species), and the reacting species' molar masses.
  std::vector<double> wR(ns);
  for (std::size_t k = 0; k < ns; ++k) wR[k] = species[reactingCol[k]].molarMass;
  std::vector<double> nPt(table.points.size() / table.pointWidth());
  for (std::size_t p = 0; p < nPt.size(); ++p) {
    const double* y = table.points.data() + p * table.pointWidth();
    for (std::size_t k = 0; k < ns; ++k) nPt[p] += y[k] / wR[k];
  }
  const std::vector<double> tauBase[2] = {baseTau(table, 0), baseTau(table, 1)};
  // The runtime, with the species' element fractions from Cantera.
  std::vector<crucible::ManifoldRuntime::Elements> elements(nsAll);
  {
    const std::size_t eH = gas->elementIndex("H"), eO = gas->elementIndex("O"), eN = gas->elementIndex("N");
    for (std::size_t k = 0; k < nsAll; ++k) {
      const double w = gas->molecularWeight(k);
      elements[k] = {gas->nAtoms(k, eH) * gas->atomicWeight(eH) / w, gas->nAtoms(k, eO) * gas->atomicWeight(eO) / w,
                     eN != Cantera::npos ? gas->nAtoms(k, eN) * gas->atomicWeight(eN) / w : 0.0};
    }
  }
  const crucible::ManifoldRuntime runtime(table, medium, elements);
  auto rwork = runtime.work();
  using RCause = crucible::ManifoldRuntime::Cause;
  long r2Every = 50;
  if (const char* e = std::getenv("CRUCIBLE_R2_EVERY")) r2Every = std::atol(e);
  // R1's tallies: cause disagreements by pair, and the worst |tau| difference where both place.
  std::map<std::string, long> r1Mismatch;
  double r1Tau = 0;
  long r1Both = 0;
  // R2's tallies.
  struct R2 { long n = 0; long aFail = 0; double b = 0, c = 0, d = 0; long dN = 0; double dWLess = 0; long dWLessN = 0; };
  R2 r2;
  const double dts[] = {1e-8, 1e-7, 1e-6, 1e-5};
  double r2e[4] = {0, 0, 0, 0}, r2eY[4] = {0, 0, 0, 0};
  std::array<double, 7> r2eWorst[4] = {};
  long r2ePartial = 0;
  std::map<std::string, long> r2Causes;
  long r2Placed = 0;
  long ownEvery = 20;
  if (const char* e = std::getenv("CRUCIBLE_E0_OWN_EQ_EVERY")) ownEvery = std::atol(e);

  struct Sample { double fH, zN, eta, rho, T; std::vector<double> y; };
  std::vector<Sample> samples;
  while (std::getline(in, line)) {
    const auto c = split(line);
    Sample s;
    auto num = [&](const char* name) { return std::strtod(c[col[name]].c_str(), nullptr); };
    s.fH = num("f_H"); s.zN = num("Z_N"); s.eta = num("eta"); s.rho = num("rho"); s.T = num("T");
    s.y.resize(nsAll);
    for (std::size_t k = 0; k < nsAll; ++k) s.y[k] = std::strtod(c[yCol[k]].c_str(), nullptr);
    samples.push_back(std::move(s));
  }

  std::FILE* out = std::fopen(argv[4], "w");
  if (!out) throw std::runtime_error(std::string("Cannot write ") + argv[4]);
  std::fprintf(out, "e0,branch,s_cell,s_eq_lower,span_lower,xi_lower,dY_lower,s_eq_upper,span_upper,xi_upper,dY_upper,s_eq_own,n_rule,n_branch,phi_N,w,tau,dY_n,n_runtime,tau_runtime\n");
  std::vector<double> yc(ns), yb(ns), yr(ns);
  std::map<std::string, long> count, countN;
  long unconverged = 0;
  std::map<int, long> iterationCount;
  std::vector<double> reachGap;
  long reachWithFrozen = 0;
  std::map<std::string, std::vector<double>> dYByClass;
  double worstSCheck = 0;
  long checked = 0, inside = 0;
  std::vector<std::array<double, 3>> ownEq;  // (s_eq_own - s_eq_lower blend, s_eq_own - s_cell, xi_lower)
  for (std::size_t n = 0; n < samples.size(); ++n) {
    const auto& s = samples[n];
    std::array<Corner, 16> corners;
    const std::string domain = cornersOf(table, s.fH, s.zN, s.eta, s.rho, corners);
    const double eCell = medium.internalEnergy(s.T, s.y.data());
    const auto rp = runtime.place(eCell, s.rho, s.T, s.y.data(), rwork);
    if (!domain.empty()) {
      ++count["(domain) " + domain];
      std::fprintf(out, "%s,-1,,,,,,,,,,,%s,-1,,,,,%s,\n", domain.c_str(), domain.c_str(), crucible::ManifoldRuntime::causeName(rp.cause));
      if (rp.cause != RCause::OutsideAxes && rp.cause != RCause::NoHydrogenOrOxygen && rp.cause != RCause::CornerRefused &&
          rp.cause != RCause::CornerFlagged)
        ++r1Mismatch[std::string("replay domain '") + domain + "', runtime '" + crucible::ManifoldRuntime::causeName(rp.cause) + "'"];
      continue;
    }
    ++inside;
    const double sCell = cellEntropy(species, s.T, s.rho, s.y.data(), table.referencePressure);
    gas->setMassFractions_NoNorm(s.y.data());
    gas->setState_TD(s.T, s.rho);
    if (n % 97 == 0) {
      worstSCheck = std::max(worstSCheck, std::fabs(sCell - gas->entropy_mass()));
      ++checked;
    }
    for (std::size_t k = 0; k < ns; ++k) yr[k] = s.y[reactingCol[k]];
    const Rule rule = applyRule(table, corners, sCell, yr.data(), yc, yb);
    ++count[rule.cause];
    double sOwn = NAN;
    if (ownEvery > 0 && inside % ownEvery == 0) {
      gas->equilibrate("UV");
      sOwn = gas->entropy_mass();
      ownEq.push_back({sOwn - rule.b[0].sEq, sOwn - sCell, rule.b[0].xi});
    }
    std::fprintf(out, "%s,%d,%.12e,%.12e,%.6e,%.6e,%.6e,%.12e,%.6e,%.6e,%.6e,", rule.cause.c_str(), rule.branch, sCell,
                 rule.b[0].sEq, rule.b[0].span, rule.b[0].xi, rule.b[0].dY, rule.b[1].sEq, rule.b[1].span, rule.b[1].xi,
                 rule.b[1].dY);
    if (std::isfinite(sOwn)) std::fprintf(out, "%.12e,", sOwn); else std::fprintf(out, ",");
    const Blend blend(table, nPt, corners);
    const NRule nr = placeByN(blend, wR, sCell, yr.data(), true, tauBase, yc, yb);
    ++countN[nr.cause];
    if (nr.iterations > 60) ++unconverged;
    if (nr.w > 0 && nr.branch >= 0) ++iterationCount[nr.iterations];
    if (nr.cause == "frozen reach") {
      reachGap.push_back(nr.gap);
      bool frozen = false;
      for (const auto& c : corners) frozen |= c.w > 0 && (table.nodes[c.node].flags & ManifoldTable::kFrozen);
      reachWithFrozen += frozen;
    }
    if (nr.branch >= 0) dYByClass[nr.cause.substr(nr.cause.find('(') + 1)].push_back(nr.dY);
    std::fprintf(out, "%s,%d,%.6e,%.6e,%.9e,%.6e,%s,%.12e\n", nr.cause.c_str(), nr.branch, nr.phi, nr.w, nr.tau, nr.dY,
                 crucible::ManifoldRuntime::causeName(rp.cause), rp.tau);
    // R1: the replay's cause names the branch and w class; the runtime's names neither, and adds the frozen corner.
    {
      const bool replayTable = nr.cause.rfind("table", 0) == 0, runtimeTable = rp.cause == RCause::Table;
      std::string replayBase = nr.cause.substr(0, nr.cause.find(" ("));
      std::string runtimeBase = crucible::ManifoldRuntime::causeName(rp.cause);
      runtimeBase = runtimeBase.substr(0, runtimeBase.find(" ("));
      const bool frozenCorner = rp.cause == RCause::FrozenCorner && replayTable;
      if (!frozenCorner && replayBase != runtimeBase) ++r1Mismatch["replay '" + nr.cause + "', runtime '" + crucible::ManifoldRuntime::causeName(rp.cause) + "'"];
      if (frozenCorner) ++r1Mismatch["(new cause) frozen corner on the branch, replay '" + nr.cause + "'"];
      if ((replayTable || frozenCorner) && (runtimeTable || frozenCorner) && rp.branch == nr.branch) {
        r1Tau = std::max(r1Tau, std::fabs(rp.tau - nr.tau));
        ++r1Both;
      }
      (void)runtimeTable;
    }
    // R2: the step at the cell.
    if (rp.cause == RCause::Table && r2Every > 0 && (r2Placed++ % r2Every) == 0) {
      ++r2.n;
      // The samples' Y are printed to 10 digits (sums off 1 by up to 1e-10): R2 works on the cell renormalised.
      std::vector<double> y0 = s.y, yq(ns);
      {
        double sum = 0;
        for (double v : y0) sum += v;
        for (double& v : y0) v /= sum;
      }
      std::vector<double> y1 = y0;
      const auto rp0 = runtime.place(eCell, s.rho, s.T, y0.data(), rwork);
      if (rp0.cause != RCause::Table) continue;
      // (a) dt = 0 gives the placed point bit for bit.
      runtime.advance(eCell, s.rho, s.T, y1.data(), 0.0, rwork);
      runtime.curve(eCell, s.rho, y0.data(), rp0.branch, rp0.u, yq.data(), rwork);
      for (std::size_t k = 0; k < ns; ++k) if (y1[reactingCol[k]] != yq[k]) { ++r2.aFail; break; }
      // (b) dt = 1 s gives the blended equilibrium.
      y1 = y0;
      runtime.advance(eCell, s.rho, s.T, y1.data(), 1.0, rwork);
      runtime.curve(eCell, s.rho, y0.data(), rp0.branch, 0.0, yq.data(), rwork);
      for (std::size_t k = 0; k < ns; ++k) r2.b = std::max(r2.b, std::fabs(y1[reactingCol[k]] - yq[k]));
      // (c) the elements, and (d) the re-placement, after dt = 1e-6 s.
      y1 = y0;
      const auto step = runtime.advance(eCell, s.rho, s.T, y1.data(), 1e-6, rwork);
      ++r2Causes[crucible::ManifoldRuntime::causeName(step.cause)];
      const auto z0 = source.elementMassFractions(y0.data()), z1 = source.elementMassFractions(y1.data());
      for (std::size_t m = 0; m < z0.size(); ++m)
        if (z0[m] > 1e-12) r2.c = std::max(r2.c, std::fabs(z1[m] - z0[m]) / z0[m]);
      const double t1 = medium.temperature(eCell, y1.data(), s.T);
      const auto again = runtime.place(eCell, s.rho, t1, y1.data(), rwork);
      if (rp0.w == 1) {
        if (again.cause == RCause::Table && again.branch == step.branch) r2.d = std::max(r2.d, std::fabs(again.tau - step.tau));
        else ++r2.dN;
      } else if (again.cause == RCause::Table) {
        r2.dWLess = std::max(r2.dWLess, std::fabs(again.tau - step.tau));
        ++r2.dWLessN;
      }
      // (e) two halves against one whole step.
      for (int q = 0; q < 4; ++q) {
        if (rp0.w < 1) { if (q == 0) ++r2ePartial; continue; }
        std::vector<double> ya = y0, yb2 = y0;
        const auto whole = runtime.advance(eCell, s.rho, s.T, ya.data(), dts[q], rwork);
        const auto half = runtime.advance(eCell, s.rho, s.T, yb2.data(), dts[q] / 2, rwork);
        if (whole.cause != RCause::Table || half.cause != RCause::Table) continue;
        const double th = medium.temperature(eCell, yb2.data(), s.T);
        const auto half2 = runtime.advance(eCell, s.rho, th, yb2.data(), dts[q] / 2, rwork);
        if (half2.cause != RCause::Table || half2.branch != whole.branch) continue;
        const double uw = whole.u, uh = half2.u;
        if (uw > 0) r2e[q] = std::max(r2e[q], std::fabs(uh - uw) / std::max(uw, uh));
        double dYq = 0;
        for (std::size_t k = 0; k < ns; ++k) dYq = std::max(dYq, std::fabs(ya[reactingCol[k]] - yb2[reactingCol[k]]));
        if (dYq > r2eY[q]) {
          r2eY[q] = dYq;
          r2eWorst[q] = {static_cast<double>(whole.branch), rp0.u, uw, uh, static_cast<double>(whole.substeps), s.eta, s.T};
        }
      }
    }
  }
  std::fclose(out);

  std::printf("E0 replay: %zu samples, %ld inside the axes with usable corners (unweighted; fallback_rate.py weights by steps and applies the gates)\n",
              samples.size(), inside);
  for (const auto& [cause, k] : count) std::printf("  %s: %ld\n", cause.c_str(), k);
  std::printf("The N placement (unweighted): %ld placements with Illinois unconverged after 60 iterations\n", unconverged);
  for (const auto& [cause, k] : countN) std::printf("  %s: %ld\n", cause.c_str(), k);
  if (!reachGap.empty()) {
    std::sort(reachGap.begin(), reachGap.end());
    std::printf("  frozen reach: (N_lo(1) - N_up(1)) / N_lo(1) min %.3e, median %.3e, max %.3e; %ld of %zu with a frozen corner of nonzero weight\n",
                reachGap.front(), reachGap[reachGap.size() / 2], reachGap.back(), reachWithFrozen, reachGap.size());
  }
  std::printf("  Illinois iterations of the placed cells with w > 0:");
  for (const auto& [it, k] : iterationCount) std::printf(" %d: %ld;", it, k);
  std::printf("\n");
  for (auto& [cls, v] : dYByClass) {
    std::sort(v.begin(), v.end());
    auto pct = [&](double p) { return v[static_cast<std::size_t>(p * (v.size() - 1))]; };
    std::printf("  dY of the placed point, (%s: %zu cells, median %.3e, 95%% %.3e, 99%% %.3e, max %.3e\n", cls.c_str(), v.size(),
                pct(0.5), pct(0.95), pct(0.99), v.back());
  }
  std::printf("R1 (the runtime's placement against the replay's): %zu kinds of disagreement; worst |tau| difference %.3e over %ld samples placed by both on the same branch\n",
              r1Mismatch.size(), r1Tau, r1Both);
  for (const auto& [what, k] : r1Mismatch) std::printf("  %s: %ld\n", what.c_str(), k);
  std::printf("R2 (the step at the cell) on %ld placed samples: (a) dt 0 not bit for bit in %ld; (b) dt 1 s against the blended equilibrium, max |dY| %.3e; "
              "(c) elements after 1e-6 s, worst relative %.3e; (d) w 1, re-placed after 1e-6 s, worst |dtau| %.3e (%ld not re-placed on the same branch); "
              "w < 1: worst |dtau| %.3e over %ld\n", r2.n, r2.aFail, r2.b, r2.c, r2.d, r2.dN, r2.dWLess, r2.dWLessN);
  for (const auto& [what, k] : r2Causes) std::printf("  the 1e-6 s step's cause: %s: %ld\n", what.c_str(), k);
  std::printf("  (e) on the cells placed at w 1 (%ld at w < 1 left out: (d) measures their re-placement)\n", r2ePartial);
  for (int q = 0; q < 4; ++q)
    std::printf("  (e) dt %.0e s: two halves against one step, worst |du| / max(u) %.3e, max |dY| %.3e (at branch %.0f, u %.6e to %.6e in one step and %.6e in two, %.0f sub-steps in one, sample eta %.4f, T %.1f K)\n",
                dts[q], r2e[q], r2eY[q], r2eWorst[q][0], r2eWorst[q][1], r2eWorst[q][2], r2eWorst[q][3], r2eWorst[q][4], r2eWorst[q][5], r2eWorst[q][6]);
  std::printf("s_cell against Cantera's entropy_mass on %ld samples: worst |difference| %.3e J/(kg K)\n", checked, worstSCheck);
  if (!ownEq.empty()) {
    std::vector<double> a, b;
    for (const auto& v : ownEq) { a.push_back(v[0]); b.push_back(v[1]); }
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    auto pct = [](const std::vector<double>& v, double p) { return v[static_cast<std::size_t>(p * (v.size() - 1))]; };
    std::printf("own s_eq (equilibrateUV) on %zu samples: s_eq_own - blended lower s_eq: min %.3f, 5%% %.3f, median %.3f, 95%% %.3f, max %.3f J/(kg K); "
                "s_eq_own - s_cell: min %.3e, median %.3f, max %.3f\n",
                ownEq.size(), a.front(), pct(a, 0.05), pct(a, 0.5), pct(a, 0.95), a.back(), b.front(), pct(b, 0.5), b.back());
  }

  // E4: s per cell, and the whole rule per cell (s, the 16 corners, both branches), on one thread, over the
  // samples inside the axes, repeated to at least 1 s.
  std::vector<std::size_t> insideIdx;
  std::vector<std::array<Corner, 16>> insideCorners;
  for (std::size_t n = 0; n < samples.size(); ++n) {
    std::array<Corner, 16> corners;
    if (cornersOf(table, samples[n].fH, samples[n].zN, samples[n].eta, samples[n].rho, corners).empty()) {
      insideIdx.push_back(n);
      insideCorners.push_back(corners);
    }
  }
  if (!insideIdx.empty()) {
    double sink = 0;
    auto timeIt = [&](auto&& body) {
      long reps = 0;
      const auto t0 = std::chrono::steady_clock::now();
      double sec = 0;
      do {
        for (std::size_t q = 0; q < insideIdx.size(); ++q) sink += body(q);
        ++reps;
        sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      } while (sec < 1.0);
      return 1e9 * sec / (static_cast<double>(reps) * static_cast<double>(insideIdx.size()));
    };
    const double nsS = timeIt([&](std::size_t q) {
      const auto& s = samples[insideIdx[q]];
      return cellEntropy(species, s.T, s.rho, s.y.data(), table.referencePressure);
    });
    const double nsRule = timeIt([&](std::size_t q) {
      const auto& s = samples[insideIdx[q]];
      const double sc = cellEntropy(species, s.T, s.rho, s.y.data(), table.referencePressure);
      for (std::size_t k = 0; k < ns; ++k) yr[k] = s.y[reactingCol[k]];
      return static_cast<double>(applyRule(table, insideCorners[q], sc, yr.data(), yc, yb).branch);
    });
    auto nTime = [&](bool check) {
      return timeIt([&](std::size_t q) {
        const auto& s = samples[insideIdx[q]];
        const double sc = cellEntropy(species, s.T, s.rho, s.y.data(), table.referencePressure);
        for (std::size_t k = 0; k < ns; ++k) yr[k] = s.y[reactingCol[k]];
        const Blend blend(table, nPt, insideCorners[q]);
        return static_cast<double>(placeByN(blend, wR, sc, yr.data(), check, tauBase, yc, yb).branch);
      });
    };
    const double nsN = nTime(false), nsNCheck = nTime(true);
    std::vector<double> ytmp;
    const double nsPlace = timeIt([&](std::size_t q) {
      const auto& s = samples[insideIdx[q]];
      return static_cast<double>(runtime.place(medium.internalEnergy(s.T, s.y.data()), s.rho, s.T, s.y.data(), rwork).branch);
    });
    const double nsAdvance = timeIt([&](std::size_t q) {
      const auto& s = samples[insideIdx[q]];
      ytmp = s.y;
      return static_cast<double>(runtime.advance(medium.internalEnergy(s.T, s.y.data()), s.rho, s.T, ytmp.data(), 1e-7, rwork).substeps);
    });
    const double nsEnergy = timeIt([&](std::size_t q) { const auto& s = samples[insideIdx[q]]; return medium.internalEnergy(s.T, s.y.data()); });
    std::printf("E4, the runtime (one thread): e from T (the timing's own input, not the runtime's) %.1f ns; place %.1f ns; advance over 1e-7 s %.1f ns per cell, each including e\n",
                nsEnergy, nsPlace, nsAdvance);
    std::printf("E4 (one thread, %zu cells): s_cell %.1f ns per cell; s_cell with the branch rule (16 corners, both branches) %.1f ns per cell; "
                "s_cell with the N placement %.1f ns per cell (Illinois on the whole branch), %.1f with the replay's monotone check (sink %g)\n",
                insideIdx.size(), nsS, nsRule, nsN, nsNCheck, sink > 0 ? 1.0 : 0.0);
  }
  return 0;
}
