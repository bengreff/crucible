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
// One output row per sample, in the samples' order, for tools/fallback_rate.py to join (it applies the
// domain: the axes and the gates). Also E4: s per cell and the whole rule per cell, timed on one thread.
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

// The corner's Y (reacting species) at u = 1 - xi on branch br: its s is s_ref - span u^2, bracketed by a
// search in c (monotone), and interpolated linearly in d = sqrt(s_ref - s) between the bracket.
void cornerY(const ManifoldTable& t, std::size_t node, const Branch& br, double u, double* out) {
  const std::size_t ns = t.reactingSpecies.size();
  const int p0 = static_cast<int>(br.first);
  if (br.span <= 0 || br.count < 2) {
    const double* y = t.pointY(node, p0);
    std::copy(y, y + ns, out);
    return;
  }
  const double sTarget = br.sRef - br.span * u * u;
  // Points ordered by increasing d: the lower branch is stored from its start (largest d) to s_ref.
  auto cAt = [&](std::size_t q) { return t.c(node, br.ascending ? p0 + static_cast<int>(br.count - 1 - q) : p0 + static_cast<int>(q)); };
  auto yAt = [&](std::size_t q) { return t.pointY(node, br.ascending ? p0 + static_cast<int>(br.count - 1 - q) : p0 + static_cast<int>(q)); };
  // In d order, s falls; find q with s(q) >= sTarget >= s(q + 1).
  std::size_t lo = 0, hi = br.count - 1;
  while (hi - lo > 1) {
    const std::size_t mid = (lo + hi) / 2;
    (cAt(mid) >= sTarget ? lo : hi) = mid;
  }
  const double dLo = std::sqrt(std::max(br.sRef - cAt(lo), 0.0)), dHi = std::sqrt(std::max(br.sRef - cAt(hi), 0.0));
  const double d = std::sqrt(br.span) * u;
  const double th = dHi > dLo ? std::clamp((d - dLo) / (dHi - dLo), 0.0, 1.0) : 0.0;
  const double* a = yAt(lo);
  const double* b = yAt(hi);
  for (std::size_t k = 0; k < ns; ++k) out[k] = a[k] + th * (b[k] - a[k]);
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
  std::fprintf(out, "e0,branch,s_cell,s_eq_lower,span_lower,xi_lower,dY_lower,s_eq_upper,span_upper,xi_upper,dY_upper,s_eq_own\n");
  std::vector<double> yc(ns), yb(ns), yr(ns);
  std::map<std::string, long> count;
  double worstSCheck = 0;
  long checked = 0, inside = 0;
  std::vector<std::array<double, 3>> ownEq;  // (s_eq_own - s_eq_lower blend, s_eq_own - s_cell, xi_lower)
  for (std::size_t n = 0; n < samples.size(); ++n) {
    const auto& s = samples[n];
    std::array<Corner, 16> corners;
    const std::string domain = cornersOf(table, s.fH, s.zN, s.eta, s.rho, corners);
    if (!domain.empty()) {
      ++count["(domain) " + domain];
      std::fprintf(out, "%s,-1,,,,,,,,,,\n", domain.c_str());
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
    if (std::isfinite(sOwn)) std::fprintf(out, "%.12e\n", sOwn); else std::fprintf(out, "\n");
  }
  std::fclose(out);

  std::printf("E0 replay: %zu samples, %ld inside the axes with usable corners (unweighted; fallback_rate.py weights by steps and applies the gates)\n",
              samples.size(), inside);
  for (const auto& [cause, k] : count) std::printf("  %s: %ld\n", cause.c_str(), k);
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
    std::printf("E4 (one thread, %zu cells): s_cell %.1f ns per cell; s_cell with the branch rule (16 corners, both branches) %.1f ns per cell (sink %g)\n",
                insideIdx.size(), nsS, nsRule, sink > 0 ? 1.0 : 0.0);
  }
  return 0;
}
