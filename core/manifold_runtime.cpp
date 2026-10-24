#include "core/manifold_runtime.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace crucible {

namespace {
constexpr double kDYMax = 0.01;
// The lower branch's induction part: below this fraction of its fall in N, a cell is not on the table
// (TABLE_B.md, "the induction cut"; measured from the table's march against CVODES on C1's front cells).
constexpr double kPhiInduction = 0.3;
constexpr double kNTol = 1e-12;
constexpr double kArgonMax = 1e-12;
constexpr double kCornerDrop = 1e-12;
constexpr std::size_t kNone = static_cast<std::size_t>(-1);

// s°/R of one species (NASA-7), with the constant-cp continuation above tHigh and the lowest range below tLow.
double entropyOverR(const Species& s, double t) {
  const double c = std::min(t, s.tHigh);
  const auto& a = c < s.tLow ? s.lowest : c < s.tMid ? s.low : s.high;
  double v = a[0] * std::log(c) + c * (a[1] + c * (a[2] / 2 + c * (a[3] / 3 + c * a[4] / 4))) + a[6];
  if (t > s.tHigh) v += (a[0] + c * (a[1] + c * (a[2] + c * (a[3] + c * a[4])))) * std::log(t / s.tHigh);
  return v;
}

struct Corner { std::size_t node; double w; };

// One branch of one corner: its stored points [first, first + count), s monotone, s_ref the end at equilibrium
// (or a frozen branch's reach), span = s_ref - s_far, frozen if s_ref is a reach.
struct Branch { std::size_t first, count; bool ascending; double sRef, span; bool frozen; };

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
  br.frozen = (n.flags & ManifoldTable::kFrozen) && br.span > 0 && br.sRef != n.cEq;
  return br;
}

// The bracket of u on branch br of one corner: its s is s_ref - span u^2, bracketed by a search in s
// (monotone), and interpolated linearly in d = sqrt(s_ref - s) between the bracket's points a and b.
struct Bracket { int a, b; double th; };
Bracket cornerBracket(const ManifoldTable& t, std::size_t node, const Branch& br, double u) {
  const int p0 = static_cast<int>(br.first);
  if (br.span <= 0 || br.count < 2) return {p0, p0, 0.0};
  const double sTarget = br.sRef - br.span * u * u;
  // Points ordered by increasing d: the lower branch is stored from its start (largest d) to s_ref.
  auto at = [&](std::size_t q) { return br.ascending ? p0 + static_cast<int>(br.count - 1 - q) : p0 + static_cast<int>(q); };
  std::size_t lo = 0, hi = br.count - 1;
  while (hi - lo > 1) {
    const std::size_t mid = (lo + hi) / 2;
    (t.c(node, at(mid)) >= sTarget ? lo : hi) = mid;
  }
  const double dLo = std::sqrt(std::max(br.sRef - t.c(node, at(lo)), 0.0)), dHi = std::sqrt(std::max(br.sRef - t.c(node, at(hi)), 0.0));
  const double d = std::sqrt(br.span) * u;
  return {at(lo), at(hi), dHi > dLo ? std::clamp((d - dLo) / (dHi - dLo), 0.0, 1.0) : 0.0};
}

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

using Cause = ManifoldRuntime::Cause;

Cause cornersOf(const ManifoldTable& t, double fH, double zN, double eta, double rho, std::array<Corner, 16>& c) {
  std::size_t ii, jj, kk, ll;
  double fi, fj, fk, fl;
  if (!std::isfinite(fH)) return Cause::NoHydrogenOrOxygen;
  if (!bracket(t.zN, zN, jj, fj) || !bracket(t.fH, fH, ii, fi) || !bracketUniform(0.0, 1.0, t.nEta, eta, kk, fk) ||
      !bracketUniform(t.lnRhoMin, t.lnRhoMax, t.nLnRho, std::log(rho), ll, fl))
    return Cause::OutsideAxes;
  // A corner weighted below kCornerDrop is dropped: not checked, not blended, the kept weights renormalised
  // (TABLE_B.md, "The bar set, the face mechanism measured, and the corner drop").
  double kept = 0;
  for (int q = 0; q < 16; ++q) {
    const int a = q >> 3 & 1, b = q >> 2 & 1, e = q >> 1 & 1, r = q & 1;
    const double w = (a ? fi : 1 - fi) * (b ? fj : 1 - fj) * (e ? fk : 1 - fk) * (r ? fl : 1 - fl);
    c[q] = {t.nodeIndex(ii + a, jj + b, kk + e, ll + r), w < kCornerDrop ? 0.0 : w};
    kept += c[q].w;
  }
  for (int q = 0; q < 16; ++q) {
    c[q].w /= kept;
    if (c[q].w == 0) continue;
    const auto& n = t.nodes[c[q].node];
    if (n.pointCount == 0 || n.lowerCount == 0 || n.lowerCount >= n.pointCount) return Cause::CornerRefused;
    if (n.flags & (ManifoldTable::kOutside | ManifoldTable::kDegenerate)) return Cause::CornerFlagged;
  }
  return Cause::Table;
}

// The 16 corners' blended curve on both branches. u on branch b: the lower branch's tau is 1 - u, the upper's 1 + u.
struct Blend {
  const ManifoldTable& t;
  const std::vector<double>& nPt;
  const std::array<Corner, 16>& corners;
  std::array<Branch, 16> br[2]{};
  double sEq[2]{}, span[2]{};
  Blend(const ManifoldTable& tab, const std::vector<double>& n, const std::array<Corner, 16>& c) : t(tab), nPt(n), corners(c) {
    for (int b = 0; b < 2; ++b)
      for (int q = 0; q < 16; ++q) {
        if (corners[q].w == 0) continue;
        br[b][q] = branchOf(t, corners[q].node, b);
        sEq[b] += corners[q].w * br[b][q].sRef;
        span[b] += corners[q].w * br[b][q].span;
      }
  }
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
      if (corners[q].w == 0) continue;
      const auto& n = t.nodes[corners[q].node];
      const std::uint32_t p = b == 0 ? (u == 1 ? 0 : n.lowerCount - 1) : (u == 0 ? n.lowerCount : n.pointCount - 1);
      v += corners[q].w * nPt[n.firstPoint + p];
    }
    return v;
  }
  void y(int b, double u, std::vector<double>& yb) const {
    const std::size_t ns = yb.size();
    std::fill(yb.begin(), yb.end(), 0.0);
    for (int q = 0; q < 16; ++q) {
      if (corners[q].w == 0) continue;
      const Bracket k = cornerBracket(t, corners[q].node, br[b][q], u);
      const double* a = t.pointY(corners[q].node, k.a);
      const double* c = t.pointY(corners[q].node, k.b);
      for (std::size_t s = 0; s < ns; ++s) yb[s] += corners[q].w * (a[s] + k.th * (c[s] - a[s]));
    }
  }
  // The blended rate: exp of the weighted mean of ln lambda over the corners of nonzero weight and span;
  // NaN if one of them has no positive finite lambda at u; +0 if none has span.
  [[nodiscard]] double lambda(int b, double u) const {
    double sum = 0, weight = 0;
    for (int q = 0; q < 16; ++q) {
      if (corners[q].w == 0 || !(br[b][q].span > 0)) continue;
      const Bracket k = cornerBracket(t, corners[q].node, br[b][q], u);
      const double la = t.lambda(corners[q].node, k.a), lb = t.lambda(corners[q].node, k.b);
      const double l = la + k.th * (lb - la);
      if (!(l > 0) || !std::isfinite(l)) return NAN;
      sum += corners[q].w * std::log(l);
      weight += corners[q].w;
    }
    return weight > 0 ? std::exp(sum / weight) : 0.0;
  }
};
}  // namespace

const char* ManifoldRuntime::causeName(Cause c) {
  switch (c) {
    case Cause::Table: return "table";
    case Cause::NoHydrogenOrOxygen: return "no H or O";
    case Cause::OutsideAxes: return "outside the axes";
    case Cause::Argon: return "AR above 1e-12";
    case Cause::CornerRefused: return "a corner was refused";
    case Cause::CornerFlagged: return "a corner flagged outside or degenerate";
    case Cause::NoFallInN: return "no fall in N";
    case Cause::FrozenReach: return "frozen reach";
    case Cause::BelowNCp: return "N below the blended N_cp";
    case Cause::OffCurveAboveSEq: return "off the curve in s (above the blended s_eq)";
    case Cause::OffCurveBelowStart: return "off the curve in s (below the blended start)";
    case Cause::DYAboveMax: return "dY above dY_max";
    case Cause::Induction: return "the lower branch's induction part (phi_N below 0.3)";
    case Cause::FrozenCorner: return "frozen corner on the branch";
    case Cause::NoRate: return "no rate at a corner";
    case Cause::Count: break;
  }
  return "?";
}

ManifoldRuntime::ManifoldRuntime(const ManifoldTable& table, const Medium& medium, std::vector<Elements> elements)
    : table_(table), medium_(medium), elements_(std::move(elements)), n2_(kNone), ar_(kNone), h2o_(kNone), h2_(kNone), o2_(kNone) {
  const auto& sp = medium_.species();
  if (elements_.size() != sp.size()) throw std::invalid_argument("Table B runtime: one element entry per medium species.");
  for (const auto& n : table_.nodes)
    if (n.pointCount > 0 && !(n.flags & ManifoldTable::kDegenerate) && n.progressVariable != 3)
      throw std::invalid_argument("Table B runtime: the table's progress variable is not the specific entropy (candidate 3).");
  if (!(table_.referencePressure > 0)) throw std::invalid_argument("Table B runtime: the table has no reference pressure.");
  reacting_.assign(table_.reactingSpecies.size(), kNone);
  for (std::size_t k = 0; k < sp.size(); ++k) {
    const auto& name = sp[k].name;
    if (name == "N2") { n2_ = k; continue; }
    if (name == "AR" || name == "Ar") { ar_ = k; continue; }
    std::size_t r = 0;
    while (r < table_.reactingSpecies.size() && table_.reactingSpecies[r] != name) ++r;
    if (r == table_.reactingSpecies.size()) throw std::invalid_argument("Table B runtime: species " + name + " is not in the table.");
    reacting_[r] = k;
    if (name == "H2O") h2o_ = k;
    if (name == "H2") h2_ = k;
    if (name == "O2") o2_ = k;
  }
  for (std::size_t r = 0; r < reacting_.size(); ++r)
    if (reacting_[r] == kNone) throw std::invalid_argument("Table B runtime: reacting species " + table_.reactingSpecies[r] + " is not in the medium.");
  if (h2o_ == kNone || h2_ == kNone || o2_ == kNone) throw std::invalid_argument("Table B runtime: H2, O2 and H2O are needed.");
  hInWater_ = elements_[h2o_].h;
  reactingMass_.resize(reacting_.size());
  for (std::size_t r = 0; r < reacting_.size(); ++r) reactingMass_[r] = sp[reacting_[r]].molarMass;
  const std::size_t ns = reacting_.size();
  pointN_.assign(table_.points.size() / table_.pointWidth(), 0.0);
  for (std::size_t p = 0; p < pointN_.size(); ++p) {
    const double* y = table_.points.data() + p * table_.pointWidth();
    for (std::size_t k = 0; k < ns; ++k) pointN_[p] += y[k] / reactingMass_[k];
  }
  // The builder's base slots (TABLE_B.md, "The entropy candidate in the builder"): lower xi = 0, 12 geometric
  // from 1e-9 to 0.1, 18 uniform to 1, so u = 1 - xi; upper u = j / 8.
  auto& lo = baseU_[0];
  lo.push_back(1.0);
  for (int g = 1; g <= table_.lowerGeometric; ++g)
    lo.push_back(1 - 1e-9 * std::pow(0.1 / 1e-9, table_.lowerGeometric > 1 ? static_cast<double>(g - 1) / (table_.lowerGeometric - 1) : 1.0));
  for (int u = 1; u <= table_.lowerUniform; ++u) lo.push_back(1 - (0.1 + 0.9 * static_cast<double>(u) / table_.lowerUniform));
  lo.back() = 0.0;
  for (int j = table_.upperUniform - 1; j >= 0; --j) baseU_[1].push_back(static_cast<double>(j) / (table_.upperUniform - 1));
}

ManifoldRuntime::Work ManifoldRuntime::work() const {
  Work w;
  w.yc.resize(reacting_.size());
  w.yb.resize(reacting_.size());
  w.yr.resize(reacting_.size());
  w.products.resize(medium_.size());
  return w;
}

double ManifoldRuntime::entropy(const Medium& medium, double t, double rho, const double* y, double pRef) {
  constexpr double r = Medium::universalGasConstant;
  const auto& sp = medium.species();
  double s = 0;
  for (std::size_t k = 0; k < sp.size(); ++k) {
    if (!(y[k] > 0)) continue;
    const double n = y[k] / sp[k].molarMass;
    s += n * r * (entropyOverR(sp[k], t) - std::log(rho * n * r * t / pRef));
  }
  return s;
}

ManifoldRuntime::Coordinates ManifoldRuntime::coordinates(double e, const double* y, Work& work) const {
  double zH = 0, zO = 0, zN = 0;
  for (std::size_t k = 0; k < elements_.size(); ++k) {
    zH += y[k] * elements_[k].h;
    zO += y[k] * elements_[k].o;
    zN += y[k] * elements_[k].n;
  }
  Coordinates c;
  c.zN = zN;
  c.fH = zH + zO > 0 ? zH / (zH + zO) : NAN;
  // Complete-combustion products of the cell's elements, as the builder's productsOf and the chamber study's samples.
  auto& p = work.products;
  std::fill(p.begin(), p.end(), 0.0);
  if (zH * (1 - hInWater_) >= zO * hInWater_) {
    p[h2o_] = zO / (1 - hInWater_);
    p[h2_] = std::max(zH - p[h2o_] * hInWater_, 0.0);
  } else {
    p[h2o_] = zH / hInWater_;
    p[o2_] = std::max(zO - p[h2o_] * (1 - hInWater_), 0.0);
  }
  if (n2_ != kNone) p[n2_] = zN;
  const double eLo = medium_.internalEnergy(200.0, p.data()), eHi = medium_.internalEnergy(6000.0, p.data());
  c.eta = (e - eLo) / (eHi - eLo);
  return c;
}

ManifoldRuntime::Placement ManifoldRuntime::place(double e, double rho, double t, const double* y, Work& work) const {
  Placement r;
  if (ar_ != kNone && y[ar_] > kArgonMax) { r.cause = Cause::Argon; return r; }
  const Coordinates co = coordinates(e, y, work);
  std::array<Corner, 16> corners;
  r.cause = cornersOf(table_, co.fH, co.zN, co.eta, rho, corners);
  if (r.cause != Cause::Table) return r;
  const Blend bl(table_, pointN_, corners);
  const std::size_t ns = reacting_.size();
  for (std::size_t k = 0; k < ns; ++k) work.yr[k] = y[reacting_[k]];
  const double sCell = entropy(medium_, t, rho, y, table_.referencePressure);
  double nCell = 0;
  for (std::size_t k = 0; k < ns; ++k) nCell += work.yr[k] / reactingMass_[k];
  const double n0 = bl.nEnd(0, 1), nLoEnd = bl.nEnd(0, 0), nUpEnd = bl.nEnd(1, 0), nCp = bl.nEnd(1, 1);
  if (!(n0 > nLoEnd)) { r.cause = Cause::NoFallInN; return r; }
  r.phi = (n0 - nCell) / (n0 - nLoEnd);
  r.w = std::clamp((r.phi - 0.005) / 0.005, 0.0, 1.0);
  int b = 0;
  double tauN = NAN;
  if (r.w > 0) {
    if (nCell >= nLoEnd) b = 0;
    else if (nCell <= nUpEnd) b = 1;
    else { r.cause = Cause::FrozenReach; return r; }
    if (b == 1 && nCell < nCp) { r.cause = Cause::BelowNCp; return r; }
    auto uOf = [&](double tau) { return b == 0 ? 1 - tau : tau - 1; };
    auto f = [&](double tau) { return bl.n(b, uOf(tau)) - nCell; };
    // Illinois on the whole branch: tau from a (N above the cell's) to c (below).
    double a = b == 0 ? 0.0 : 1.0, c = b == 0 ? 1.0 : 2.0;
    double fa = (b == 0 ? n0 : nUpEnd) - nCell, fc = (b == 0 ? nLoEnd : nCp) - nCell;
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
  double u;
  if (r.w == 1) {
    u = b == 0 ? 1 - tauN : tauN - 1;
  } else {
    const double sStar = r.w == 0 ? sCell : (1 - r.w) * sCell + r.w * (bl.sEq[0] - bl.span[0] * (1 - tauN) * (1 - tauN));
    if (sStar > bl.sEq[0]) { r.cause = Cause::OffCurveAboveSEq; return r; }
    if (!(bl.span[0] > 0) || sStar < bl.sEq[0] - bl.span[0]) { r.cause = Cause::OffCurveBelowStart; return r; }
    u = std::sqrt((bl.sEq[0] - sStar) / bl.span[0]);
  }
  r.branch = b;
  r.u = u;
  r.tau = b == 0 ? 1 - u : 1 + u;
  bl.y(b, u, work.yb);
  r.dY = 0;
  for (std::size_t k = 0; k < ns; ++k) r.dY = std::max(r.dY, std::fabs(work.yb[k] - work.yr[k]));
  if (r.dY > kDYMax) { r.cause = Cause::DYAboveMax; return r; }
  if (b == 0 && r.phi < kPhiInduction) { r.cause = Cause::Induction; return r; }
  for (int q = 0; q < 16; ++q)
    if (corners[q].w > 0 && bl.br[b][q].frozen && corners[q].w > r.frozenW) {
      r.cause = Cause::FrozenCorner;
      r.frozenW = corners[q].w;
      r.frozenNode = corners[q].node;
    }
  return r;
}

ManifoldRuntime::Placement ManifoldRuntime::advance(double e, double rho, double t, double* y, double dt, Work& work) const {
  Placement r = place(e, rho, t, y, work);
  if (r.cause != Cause::Table) return r;
  // The corners again (place() keeps none): the same coordinates give the same corners and weights.
  const Coordinates co = coordinates(e, y, work);
  std::array<Corner, 16> corners;
  cornersOf(table_, co.fH, co.zN, co.eta, rho, corners);
  const Blend bl(table_, pointN_, corners);
  const int b = r.branch;
  const auto& base = baseU_[b];
  double u = r.u, remaining = dt;
  while (remaining > 0 && u > 0) {
    // The next base u strictly below u (base is decreasing and ends at 0).
    std::size_t j = 0;
    while (j + 1 < base.size() && base[j] >= u) ++j;
    const double ub = base[j] < u ? base[j] : 0.0;
    const double l0 = bl.lambda(b, u);
    if (std::isnan(l0)) { r.cause = Cause::NoRate; return r; }
    if (l0 == 0) break;  // no corner with span: the branch is its equilibrium at every u
    const double up = std::max(u * std::exp(-l0 * remaining / 2), ub);
    const double lm = bl.lambda(b, 0.5 * (u + up));
    if (!(lm > 0)) { r.cause = Cause::NoRate; return r; }
    ++r.substeps;
    if (ub > 0) {
      const double tb = 2 * std::log(u / ub) / lm;
      if (tb < remaining) { u = ub; remaining -= tb; continue; }
    }
    u *= std::exp(-lm * remaining / 2);
    remaining = 0;
  }
  r.u = u;
  r.tau = b == 0 ? 1 - u : 1 + u;
  bl.y(b, u, work.yb);
  for (std::size_t k = 0; k < reacting_.size(); ++k) y[reacting_[k]] = work.yb[k];
  return r;
}

bool ManifoldRuntime::curve(double e, double rho, const double* y, int branch, double u, double* out, Work& work) const {
  const Coordinates co = coordinates(e, y, work);
  std::array<Corner, 16> corners;
  if (cornersOf(table_, co.fH, co.zN, co.eta, rho, corners) != Cause::Table) return false;
  const Blend bl(table_, pointN_, corners);
  bl.y(branch, u, work.yb);
  std::copy(work.yb.begin(), work.yb.end(), out);
  return true;
}

}  // namespace crucible
