// C++ port of tools/sst_wall_layer.py (and the parts of tools/sst_pipe_1d.py and tools/eq_wall_model.py
// it exercises with planar=True, constant_pressure=True, correction="hp", dic=True, hp_form="full",
// their_closure=False): the corrected-SST planar constant-stress, constant-heat-flux inner layer used
// as a wall model, and its a priori check against the pipe references. Standalone: reads
// tools/n2_properties.csv directly (no Cantera, no engine dependency).
//
// Branches of sst_pipe_1d.py's Pipe not taken under these settings (powerlaw, non-planar geometry, the
// segregated solver, check_operator, Pipe's own start/measures, hp_form "inner", their_closure) are not
// ported; see sst_wall_layer.py and sst_pipe_1d.py for them.
//
// Build: target crucible_wall_layer (CMakeLists.txt), unconditional in the default configuration.
// Usage: crucible_wall_layer apriori <heated_hp|cold_hp|heated|cold> [y+ targets, comma-separated]
//        crucible_wall_layer point <u1> <T1> <y1> <p> <Tw> <tau0> <q0>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// --- SST-2003 constants and the pipe reference's radius (sst_pipe_1d.py, eq_wall_model.py) ---
constexpr std::array<double, 2> SIGMA_K{0.85, 1.0};
constexpr std::array<double, 2> SIGMA_W{0.5, 0.856};
constexpr std::array<double, 2> BETA{0.075, 0.0828};
constexpr std::array<double, 2> GAMMA{5.0 / 9.0, 0.44};
constexpr double BETA_STAR = 0.09, A1 = 0.31;
constexpr double RADIUS = 5e-3;

// D^ic = D(Rt, Mt) / D(Rt, 0), D = [1 - exp(-Rt / (3.5 + 0.39 Mt^0.77))]^2 (Hasan, Elias, Menter and
// Pecnik, J. Fluid Mech. 1019, A8, 2025, eq. 4.3-4.5).
double hpDamping(double rt, double mt) {
  rt = std::max(rt, 1e-12);
  const double r = std::expm1(-rt / (3.5 + 0.39 * std::pow(mt, 0.77))) / std::expm1(-rt / 3.5);
  return r * r;
}

// np.geomspace(a, b, count): count log-spaced points from a to b inclusive, endpoints exact.
std::vector<double> geomspace(double a, double b, int count) {
  std::vector<double> v(count);
  if (count == 1) {
    v[0] = a;
    return v;
  }
  const double ratio = b / a;
  for (int i = 0; i < count; ++i) v[i] = a * std::pow(ratio, static_cast<double>(i) / (count - 1));
  v.front() = a;
  v.back() = b;
  return v;
}

// --- N2 properties (tools/n2_properties.csv): mu, lambda, cp linear in T (np.interp, clamped at the
// ends), R constant (the file's first row).
struct Properties {
  std::vector<double> t, mu, lam, cp;
  double gasR{};

  static Properties load(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("cannot open " + path);
    std::string line;
    std::getline(f, line);  // comment
    std::getline(f, line);  // header
    Properties p;
    while (std::getline(f, line)) {
      if (line.empty()) continue;
      std::stringstream ss(line);
      std::string tok;
      std::vector<double> v;
      while (std::getline(ss, tok, ',')) v.push_back(std::stod(tok));
      p.t.push_back(v[0]);
      p.mu.push_back(v[1]);
      p.lam.push_back(v[2]);
      p.cp.push_back(v[3]);
      if (p.gasR == 0.0) p.gasR = v[4];
    }
    return p;
  }

  static double interp(double x, const std::vector<double>& xs, const std::vector<double>& ys) {
    if (x <= xs.front()) return ys.front();
    if (x >= xs.back()) return ys.back();
    const auto it = std::upper_bound(xs.begin(), xs.end(), x);
    const std::size_t i = static_cast<std::size_t>(it - xs.begin());  // xs[i-1] <= x < xs[i]
    const double f = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
    return ys[i - 1] + f * (ys[i] - ys[i - 1]);
  }
  double muAt(double T) const { return interp(T, t, mu); }
  double lamAt(double T) const { return interp(T, t, lam); }
  double cpAt(double T) const { return interp(T, t, cp); }
};

const Properties& props() {
  static const Properties p = Properties::load(std::string(CRUCIBLE_SOURCE_DIR) + "/tools/n2_properties.csv");
  return p;
}

// Closure fields of Pipe.closure/hp_fields for one cell array (sst_pipe_1d.py).
struct Closure {
  std::vector<double> rho, mu, lam, cp, muT, sk, sw, beta, gamma, prod, omegaProd, cross, f1, sn;
  double rhoWall{}, wWall{};
};

struct Diffusion {
  std::vector<double> k, w;  // volume-integrated k and omega diffusion (corrected, hp full form)
};

// A planar constant-stress, constant-heat-flux layer (Pipe with planar=True, constant_pressure=True;
// Layer in sst_wall_layer.py), cells ordered top (j = 0) to wall (j = n - 1).
class Layer {
 public:
  int n{}, nInner{};
  double y1{}, top{}, p0{}, tWall{}, prT{}, omegaFactor{};
  double gasR{}, muWall{}, lamWall{};
  std::vector<double> d, vol, dc;  // size n, n, n - 1
  double dw{};
  double tau{0.0}, flux{0.0}, layerForce{0.0}, layerHeat{0.0};

  Layer(double y1_, double top_, double p, double tw, double first, int nInner_, int nOuter_,
        double prT_ = 0.9, double omegaFactor_ = 10.0)
      : n(nInner_ + nOuter_), nInner(nInner_), y1(y1_), top(top_), p0(p), tWall(tw), prT(prT_),
        omegaFactor(omegaFactor_) {
    const auto& P = props();
    gasR = P.gasR;
    muWall = P.muAt(tw);
    lamWall = P.lamAt(tw);
    std::vector<double> yf(n + 1);
    yf[0] = 0.0;
    const auto inner = geomspace(first, y1, nInner_);
    for (int i = 0; i < nInner_; ++i) yf[1 + i] = inner[i];
    const auto outer = geomspace(y1, top, nOuter_ + 1);
    for (int i = 0; i < nOuter_; ++i) yf[1 + nInner_ + i] = outer[1 + i];
    d.resize(n);
    vol.resize(n);
    for (int i = 0; i < n; ++i) {
      const int j = n - 1 - i;  // the pipe's order (top first) from yf's wall-first order
      d[i] = 0.5 * (yf[j] + yf[j + 1]);
      vol[i] = yf[j + 1] - yf[j];
    }
    dc.resize(n - 1);
    for (int i = 0; i < n - 1; ++i) dc[i] = d[i] - d[i + 1];
    dw = d[n - 1];
  }

  // Layer.grad: zero-gradient continuation at the top (face[0] = face[1]), not the pipe's axis face[0] = 0.
  std::vector<double> grad(const std::vector<double>& phi, double wall) const {
    std::vector<double> face(n + 1);
    for (int j = 1; j < n; ++j) face[j] = (phi[j] - phi[j - 1]) / dc[j - 1];
    face[n] = (wall - phi[n - 1]) / dw;
    face[0] = face[1];
    std::vector<double> g(n);
    for (int j = 0; j < n; ++j) g[j] = 0.5 * (face[j] + face[j + 1]);
    return g;
  }

  std::vector<double> conductance(const std::vector<double>& gammaC, double gammaWall) const {
    std::vector<double> c(n + 1);
    c[0] = 0.0;
    for (int j = 1; j < n; ++j) c[j] = 0.5 * (gammaC[j - 1] + gammaC[j]) / dc[j - 1];
    c[n] = gammaWall / dw;
    return c;
  }

  std::vector<double> faceFlux(const std::vector<double>& cond, const std::vector<double>& phi, double wall) const {
    std::vector<double> flux(n + 1, 0.0);
    for (int j = 1; j < n; ++j) flux[j] = cond[j] * (phi[j] - phi[j - 1]);
    flux[n] = cond[n] * (wall - phi[n - 1]);
    return flux;
  }

  // Pipe.closure (SST blending F1/F2, the limiter, production) plus the HP compressibility damping of
  // mu_t (correction == "hp" and dic, always here) and Pipe.hp_fields (hp_form == "full" always here).
  Closure closure(const double* u, const double* k, const double* w, const double* t, double pu) const {
    const auto& P = props();
    Closure s;
    s.rho.resize(n);
    s.mu.resize(n);
    s.lam.resize(n);
    s.cp.resize(n);
    for (int j = 0; j < n; ++j) {
      s.rho[j] = pu / (gasR * t[j] + 2.0 / 3.0 * k[j]);
      s.mu[j] = P.muAt(t[j]);
      s.lam[j] = P.lamAt(t[j]);
      s.cp[j] = P.cpAt(t[j]);
    }
    s.rhoWall = pu / (gasR * tWall);
    s.wWall = omegaFactor * 6.0 * (muWall / s.rhoWall) / (BETA[0] * dw * dw);
    const std::vector<double> uV(u, u + n), kV(k, k + n), wV(w, w + n);
    const auto dudr = grad(uV, 0.0);
    const auto dkdr = grad(kV, 0.0);
    const auto dwdr = grad(wV, s.wWall);
    s.muT.resize(n);
    s.sk.resize(n);
    s.sw.resize(n);
    s.beta.resize(n);
    s.gamma.resize(n);
    s.prod.resize(n);
    s.omegaProd.resize(n);
    s.cross.resize(n);
    s.f1.resize(n);
    for (int j = 0; j < n; ++j) {
      const double nu = s.mu[j] / s.rho[j];
      const double strain = std::abs(dudr[j]);
      const double cd = std::max(2.0 * s.rho[j] * SIGMA_W[1] / w[j] * dkdr[j] * dwdr[j], 1e-10);
      const double arg1 =
          std::min(std::max(std::sqrt(k[j]) / (BETA_STAR * w[j] * d[j]), 500.0 * nu / (d[j] * d[j] * w[j])),
                    4.0 * s.rho[j] * SIGMA_W[1] * k[j] / (cd * d[j] * d[j]));
      const double f1 = std::tanh(arg1 * arg1 * arg1 * arg1);
      const double arg2 = std::max(2.0 * std::sqrt(k[j]) / (BETA_STAR * w[j] * d[j]), 500.0 * nu / (d[j] * d[j] * w[j]));
      const double f2 = std::tanh(arg2 * arg2);
      double muT = s.rho[j] * A1 * k[j] / std::max(A1 * w[j], strain * f2);
      const double sound = std::sqrt(s.cp[j] / (s.cp[j] - gasR) * gasR * t[j]);
      muT *= hpDamping(s.rho[j] * k[j] / (s.mu[j] * w[j]), std::sqrt(2.0 * std::max(k[j], 0.0)) / sound);
      s.muT[j] = muT;
      s.f1[j] = f1;
      s.sk[j] = f1 * SIGMA_K[0] + (1 - f1) * SIGMA_K[1];
      s.sw[j] = f1 * SIGMA_W[0] + (1 - f1) * SIGMA_W[1];
      s.beta[j] = f1 * BETA[0] + (1 - f1) * BETA[1];
      s.gamma[j] = f1 * GAMMA[0] + (1 - f1) * GAMMA[1];
      s.prod[j] = std::min(muT * strain * strain, 10.0 * BETA_STAR * s.rho[j] * w[j] * k[j]);
      s.omegaProd[j] = s.gamma[j] * s.rho[j] * (s.prod[j] / muT);
      s.cross[j] = 2.0 * (1 - f1) * s.rho[j] * SIGMA_W[1] * dkdr[j] * dwdr[j];
    }
    hpFields(s, kV, wV);
    return s;
  }

  // S_n (floored at psi / 10) and, in the full form, the corrected cross-diffusion.
  void hpFields(Closure& s, const std::vector<double>& k, const std::vector<double>& w) const {
    std::vector<double> psi(n);
    for (int j = 0; j < n; ++j) psi[j] = std::sqrt(s.rho[j]) / s.mu[j];
    const double sqW = std::sqrt(s.rhoWall);
    const auto gradPsi = grad(psi, sqW / muWall);
    s.sn.resize(n);
    for (int j = 0; j < n; ++j) {
      const double den = psi[j] - d[j] * gradPsi[j];
      s.sn[j] = 1.0 / std::max(den, 0.1 * psi[j]);
    }
    std::vector<double> rk(n), sqw(n);
    for (int j = 0; j < n; ++j) {
      rk[j] = s.rho[j] * k[j];
      sqw[j] = std::sqrt(s.rho[j]) * w[j];
    }
    const auto gradRk = grad(rk, 0.0);
    const auto gradSqw = grad(sqw, sqW * s.wWall);
    for (int j = 0; j < n; ++j) {
      const double sq = std::sqrt(s.rho[j]);
      s.cross[j] = 2.0 * (1 - s.f1[j]) * SIGMA_W[1] / sq * gradRk[j] * gradSqw[j];
    }
  }

  // Pipe.diffusion, hp full form: f1 * inner + (1 - f1) * outer for both k and omega.
  Diffusion diffusion(const Closure& s, const std::vector<double>& k, const std::vector<double>& w) const {
    std::vector<double> muk(n), muw(n), snMu(n), rk(n), muW(n);
    for (int j = 0; j < n; ++j) {
      muk[j] = s.mu[j] + s.sk[j] * s.muT[j];
      muw[j] = s.mu[j] + s.sw[j] * s.muT[j];
      snMu[j] = s.sn[j] / s.mu[j];
      rk[j] = s.rho[j] * k[j];
      muW[j] = s.mu[j] * w[j];
    }
    const double sqW = std::sqrt(s.rhoWall), wWall = s.wWall;
    std::vector<double> mukSn(n), muwSn(n);
    for (int j = 0; j < n; ++j) {
      mukSn[j] = muk[j] * snMu[j];
      muwSn[j] = muw[j] * snMu[j];
    }
    const auto fluxInnerK = faceFlux(conductance(mukSn, muWall / sqW), rk, 0.0);
    const auto fluxInnerW = faceFlux(conductance(muwSn, muWall / sqW), muW, muWall * wWall);
    std::vector<double> innerK(n), innerW(n);
    for (int j = 0; j < n; ++j) {
      innerK[j] = snMu[j] * (fluxInnerK[j + 1] - fluxInnerK[j]);
      innerW[j] = s.rho[j] / s.mu[j] * snMu[j] * (fluxInnerW[j + 1] - fluxInnerW[j]);
    }
    std::vector<double> mukSq(n), muwSq(n), sqw(n);
    for (int j = 0; j < n; ++j) {
      const double sq = std::sqrt(s.rho[j]);
      mukSq[j] = muk[j] / sq;
      muwSq[j] = muw[j] / sq;
      sqw[j] = sq * w[j];
    }
    const auto fluxOuterK = faceFlux(conductance(mukSq, muWall / sqW), rk, 0.0);
    const auto fluxOuterW = faceFlux(conductance(muwSq, muWall / sqW), sqw, sqW * wWall);
    Diffusion out;
    out.k.resize(n);
    out.w.resize(n);
    for (int j = 0; j < n; ++j) {
      const double sq = std::sqrt(s.rho[j]);
      const double outerK = (fluxOuterK[j + 1] - fluxOuterK[j]) / sq;
      const double outerW = fluxOuterW[j + 1] - fluxOuterW[j];
      out.k[j] = s.f1[j] * innerK[j] + (1 - s.f1[j]) * outerK;
      out.w[j] = s.f1[j] * innerW[j] + (1 - s.f1[j]) * outerW;
    }
    return out;
  }

  // Layer.residual: the top face's shear and total energy flux are prescribed (self.tau, self.flux); k
  // and omega have zero flux there (the default of face_flux, not overridden).
  std::vector<double> residual(const std::vector<double>& x) const {
    const double* u = x.data();
    const double* k = x.data() + n;
    const double* w = x.data() + 2 * n;
    const double* t = x.data() + 3 * n;
    const double pu = x[4 * n];
    const auto s = closure(u, k, w, t, pu);
    const std::vector<double> uV(u, u + n), kV(k, k + n), wV(w, w + n), tV(t, t + n);
    std::vector<double> muPlusMuT(n);
    for (int j = 0; j < n; ++j) muPlusMuT[j] = s.mu[j] + s.muT[j];
    auto fu = faceFlux(conductance(muPlusMuT, muWall), uV, 0.0);
    fu[0] = -(tau - layerForce * top);
    const auto diff = diffusion(s, kV, wV);
    std::vector<double> condT(n);
    for (int j = 0; j < n; ++j) condT[j] = s.lam[j] + s.cp[j] * s.muT[j] / prT;
    auto ft = faceFlux(conductance(condT, lamWall), tV, tWall);
    std::vector<double> uf(n + 1, 0.0);
    for (int j = 1; j < n; ++j) uf[j] = 0.5 * (u[j - 1] + u[j]);
    uf[0] = u[0] + (u[0] - u[1]) * (top - d[0]) / dc[0];
    ft[0] = -flux - uf[0] * fu[0];
    std::vector<double> res(4 * n + 1);
    for (int j = 0; j < n; ++j) res[j] = (fu[j + 1] - fu[j]) + layerForce * vol[j];
    for (int j = 0; j < n; ++j) res[n + j] = diff.k[j] + (s.prod[j] - BETA_STAR * s.rho[j] * w[j] * k[j]) * vol[j];
    for (int j = 0; j < n; ++j)
      res[2 * n + j] = diff.w[j] + (s.omegaProd[j] - s.beta[j] * s.rho[j] * w[j] * w[j] + s.cross[j] / w[j]) * vol[j];
    for (int j = 0; j < n; ++j) {
      const double netFt = ft[j + 1] - ft[j];
      const double netUfFu = uf[j + 1] * fu[j + 1] - uf[j] * fu[j];
      res[3 * n + j] = netFt + netUfFu + diff.k[j] + (layerForce * u[j] + layerHeat) * vol[j];
    }
    res[4 * n] = pu - p0;
    return res;
  }

  // Pipe.jacobian: a 5-colour, bandwidth-2 (per variable block) finite-difference Jacobian, plus the
  // dense last column (d residual / d pu) and the trivial last row (pu - p0, constant_pressure).
  std::vector<Eigen::Triplet<double>> jacobianTriplets(const std::vector<double>& x, const std::vector<double>& r0) const {
    std::vector<Eigen::Triplet<double>> trip;
    trip.reserve(static_cast<std::size_t>(n) * 90 + 4 * n + 1);
    double maxAbsU = 0.0, maxK = 0.0;
    for (int j = 0; j < n; ++j) {
      maxAbsU = std::max(maxAbsU, std::abs(x[j]));
      maxK = std::max(maxK, x[n + j]);
    }
    const double floorArr[4] = {1e-10 * maxAbsU, 1e-14 * maxK, 0.0, 1e-9 * tWall};
    for (int v = 0; v < 4; ++v) {
      for (int color = 0; color < 5; ++color) {
        std::vector<int> cells;
        for (int c = color; c < n; c += 5) cells.push_back(c);
        std::vector<double> h(cells.size());
        std::vector<double> xp = x;
        for (std::size_t m = 0; m < cells.size(); ++m) {
          h[m] = 1e-7 * std::abs(x[v * n + cells[m]]) + floorArr[v];
          xp[v * n + cells[m]] += h[m];
        }
        const auto rp = residual(xp);
        for (int off = -2; off <= 2; ++off) {
          for (std::size_t m = 0; m < cells.size(); ++m) {
            const int rcell = cells[m] + off;
            if (rcell < 0 || rcell >= n) continue;
            for (int eq = 0; eq < 4; ++eq) {
              const double dr = rp[eq * n + rcell] - r0[eq * n + rcell];
              trip.emplace_back(eq * n + rcell, v * n + cells[m], dr / h[m]);
            }
          }
        }
      }
    }
    const double hp = 1e-7 * x[4 * n];
    std::vector<double> xp = x;
    xp[4 * n] += hp;
    const auto rp = residual(xp);
    for (int row = 0; row < 4 * n; ++row) trip.emplace_back(row, 4 * n, (rp[row] - r0[row]) / hp);
    trip.emplace_back(4 * n, 4 * n, 1.0);  // constant_pressure: d(pu - p0)/d(pu) = 1
    return trip;
  }

  // Layer.newton: pseudo-time Newton with a positivity step limit on k and omega.
  std::vector<double> newton(std::vector<double> x, double tol = 1e-12, int maxSteps = 400) {
    const double rho0 = p0 / (gasR * tWall);
    const double uTau = std::sqrt(tau / rho0);
    double dtau = 1e-3 * top / uTau;
    double change = 0.0;
    for (int step = 0; step < maxSteps; ++step) {
      const auto r0 = residual(x);
      const auto jTrip = jacobianTriplets(x, r0);
      const auto s = closure(x.data(), x.data() + n, x.data() + 2 * n, x.data() + 3 * n, x[4 * n]);
      std::vector<Eigen::Triplet<double>> aTrip;
      aTrip.reserve(jTrip.size() + static_cast<std::size_t>(4 * n));
      for (const auto& t : jTrip) aTrip.emplace_back(t.row(), t.col(), -t.value());
      for (int j = 0; j < n; ++j) {
        aTrip.emplace_back(j, j, s.rho[j] * vol[j] / dtau);
        aTrip.emplace_back(n + j, n + j, s.rho[j] * vol[j] / dtau);
        aTrip.emplace_back(2 * n + j, 2 * n + j, s.rho[j] * vol[j] / dtau);
        aTrip.emplace_back(3 * n + j, 3 * n + j, s.rho[j] * s.cp[j] * vol[j] / dtau);
      }
      Eigen::SparseMatrix<double> A(4 * n + 1, 4 * n + 1);
      A.setFromTriplets(aTrip.begin(), aTrip.end());
      A.makeCompressed();
      if (!patternReady_) {
        solver_.analyzePattern(A);
        patternReady_ = true;
      }
      solver_.factorize(A);
      if (solver_.info() != Eigen::Success) throw std::runtime_error("layer Newton: sparse factorisation failed");
      const Eigen::VectorXd rhs = Eigen::Map<const Eigen::VectorXd>(r0.data(), static_cast<Eigen::Index>(r0.size()));
      const Eigen::VectorXd dx = solver_.solve(rhs);
      double alpha = 1.0;
      for (int v : {1, 2}) {
        for (int j = 0; j < n; ++j) {
          const double dv = dx[v * n + j];
          if (dv < 0) alpha = std::min(alpha, 0.5 * x[v * n + j] / (-dv));
        }
      }
      for (int i = 0; i < 4 * n + 1; ++i) x[i] += alpha * dx[i];
      double mdxU = 0, mxU = 0, mdxK = 0, mxK = 0, mw = 0, mdxT = 0;
      for (int j = 0; j < n; ++j) {
        mdxU = std::max(mdxU, std::abs(dx[j]));
        mxU = std::max(mxU, std::abs(x[j]));
        mdxK = std::max(mdxK, std::abs(dx[n + j]));
        mxK = std::max(mxK, x[n + j]);
        mw = std::max(mw, std::abs(dx[2 * n + j]) / x[2 * n + j]);
        mdxT = std::max(mdxT, std::abs(dx[3 * n + j]));
      }
      change = std::max({mdxU / mxU, mdxK / mxK, mw, mdxT / tWall});
      dtau = (alpha == 1.0) ? dtau * 3.0 : dtau * 0.5;
      if (alpha == 1.0 && change < tol && dtau * uTau / top > 1e6) return x;
    }
    throw std::runtime_error("layer Newton did not converge (last update " + std::to_string(change) + ")");
  }

  // Layer.start: the Kawai-Larsson mixing-length layer (kappa 0.41, A+ 17) for u and T; k and omega from it.
  std::vector<double> start(double tauIn, double q) const {
    const auto& P = props();
    std::vector<double> u(n), t(n);
    double uu = 0.0, tt = tWall, ya = 0.0;
    for (int idx = n - 1; idx >= 0; --idx) {
      const double yb = d[idx];
      for (int i = 0; i < 8; ++i) {
        const double a = ya + (yb - ya) * i / 8.0, b = ya + (yb - ya) * (i + 1) / 8.0;
        const double mu = P.muAt(tt), lam = P.lamAt(tt), cp = P.cpAt(tt);
        const double rho = p0 / (gasR * tt);
        const double ym = 0.5 * (a + b);
        const double sq = std::sqrt(rho * tauIn);
        const double damp = 1.0 - std::exp(-ym * sq / mu / 17.0);
        const double mut = 0.41 * ym * sq * damp * damp;
        uu += (b - a) * tauIn / (mu + mut);
        tt += (b - a) * (q - tauIn * uu) / (cp * (lam / cp + mut / prT));
      }
      u[idx] = uu;
      t[idx] = tt;
      ya = yb;
    }
    std::vector<double> k(n), w(n);
    for (int j = 0; j < n; ++j) {
      const double rho = p0 / (gasR * t[j]);
      const double mu = P.muAt(t[j]);
      const double ystar = d[j] * std::sqrt(rho * tauIn) / mu;
      k[j] = tauIn / (rho * std::sqrt(BETA_STAR)) * std::min(1.0, (ystar / 10.0) * (ystar / 10.0));
      w[j] = std::max(std::sqrt(tauIn / rho) / (std::sqrt(BETA_STAR) * 0.41 * d[j]), 6.0 * mu / rho / (BETA[0] * d[j] * d[j]));
    }
    std::vector<double> x(4 * n + 1);
    for (int j = 0; j < n; ++j) {
      x[j] = u[j];
      x[n + j] = k[j];
      x[2 * n + j] = w[j];
      x[3 * n + j] = t[j];
    }
    x[4 * n] = p0;
    return x;
  }

  std::pair<double, double> atY1(const std::vector<double>& x) const {
    const int below = n - nInner, above = below - 1;
    const double f = (y1 - d[below]) / (d[above] - d[below]);
    const double u = x[below] + f * (x[above] - x[below]);
    const double t = x[3 * n + below] + f * (x[3 * n + above] - x[3 * n + below]);
    return {u, t};
  }

  std::pair<double, double> wall(const std::vector<double>& x) const {
    return {muWall * x[n - 1] / dw, lamWall * (x[3 * n + (n - 1)] - tWall) / dw};
  }

 private:
  mutable Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_;
  mutable bool patternReady_{false};
};

struct ModelResult {
  double tauW{}, qW{};
};

// sst_wall_layer.py's model(): the top at c * max(y1, 30 wall units), the 2-D Newton on (ln tau_w, the
// scaled top flux) with a finite-difference Jacobian, warm-started from the previous inversion's state.
ModelResult model(double u1, double t1, double y1, double p, double tw, double tau0, double q0, double c = 10.0,
                   int nInner = 200, int nOuter = 60, double sourceForce = 0.0, double sourceHeat = 0.0) {
  const auto& P = props();
  const double rhoW = p / (P.gasR * tw);
  const double unit = P.muAt(tw) / (rhoW * std::sqrt(tau0 / rhoW));
  Layer layer(y1, c * std::max(y1, 30.0 * unit), p, tw, std::min(0.05 * unit, 0.01 * y1), nInner, nOuter);
  layer.layerForce = sourceForce;
  layer.layerHeat = sourceHeat;
  const double qScale = std::max(std::abs(q0), 1e-4 * rhoW * P.cpAt(tw) * std::sqrt(tau0 / rhoW) * tw);
  const double tScale = std::max(std::abs(t1 - tw), 1e-3 * tw);
  std::vector<double> state = layer.start(tau0, q0);

  auto miss = [&](const std::array<double, 2>& v) -> std::array<double, 2> {
    layer.tau = std::exp(v[0]);
    layer.flux = v[1] * qScale;
    state = layer.newton(state);
    const auto [ua, ta] = layer.atY1(state);
    return {(ua - u1) / u1, (ta - t1) / tScale};
  };

  std::array<double, 2> v{std::log(tau0), q0 / qScale};
  bool converged = false;
  for (int it = 0; it < 40 && !converged; ++it) {
    const auto f = miss(v);
    const auto base = state;
    if (std::max(std::abs(f[0]), std::abs(f[1])) < 1e-10) {
      converged = true;
      break;
    }
    double jac[2][2];
    const double h[2] = {1e-6, 1e-6 * std::max(std::abs(v[1]), 1.0)};
    for (int j = 0; j < 2; ++j) {
      auto vp = v;
      vp[j] += h[j];
      state = base;
      const auto fp = miss(vp);
      jac[0][j] = (fp[0] - f[0]) / h[j];
      jac[1][j] = (fp[1] - f[1]) / h[j];
    }
    state = base;
    const double det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
    const double rhs0 = -f[0], rhs1 = -f[1];
    double dv0 = (rhs0 * jac[1][1] - jac[0][1] * rhs1) / det;
    const double dv1 = (jac[0][0] * rhs1 - rhs0 * jac[1][0]) / det;
    dv0 = std::max(-0.3, std::min(0.3, dv0));
    v[0] += dv0;
    v[1] += dv1;
  }
  if (!converged) throw std::runtime_error("inversion did not converge");
  const auto [tauW, qW] = layer.wall(state);
  return {tauW, qW};
}

// --- eq_wall_model.py: CASES and references() (the reference profile CSV, r, u, T, rho, tau, q, p_wall). ---
struct CaseDef {
  std::string profile;
  double twall, force, heat;
};

const std::map<std::string, CaseDef>& cases() {
  static const std::map<std::string, CaseDef> c = {
      {"cold", {"reference_re_tau_1e4_profile_6400.csv", 300.0, 45690.0, 0.0}},
      {"heated", {"reference_heated_profile_6400.csv", 600.0, 19800.0, 1.95e9}},
      {"cold_hp", {"reference_hp_cold_profile_6400.csv", 300.0, 45690.0, 0.0}},
      {"heated_hp", {"reference_hp_heated_profile_6400.csv", 600.0, 19800.0, 1.95e9}},
  };
  return c;
}

struct Reference {
  std::vector<double> r, u, t, rho;
  double tau{}, q{}, pWall{}, tWall{};
};

Reference references(const std::string& caseName) {
  const auto it = cases().find(caseName);
  if (it == cases().end()) throw std::runtime_error("unknown case " + caseName);
  const auto& c = it->second;
  const std::string path = std::string(CRUCIBLE_SOURCE_DIR) + "/docs/evidence/wall_functions/" + c.profile;
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open " + path);
  std::string line;
  std::getline(f, line);  // header
  std::vector<double> r, u, k, t, rho;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::string tok;
    std::vector<double> v;
    while (std::getline(ss, tok, ',')) v.push_back(std::stod(tok));
    r.push_back(v[0]);
    u.push_back(v[1]);
    k.push_back(v[2]);
    t.push_back(v[4]);
    rho.push_back(v[5]);
  }
  std::vector<std::size_t> order(r.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return r[a] < r[b]; });
  Reference out;
  const std::size_t N = r.size();
  out.r.resize(N);
  out.u.resize(N);
  out.t.resize(N);
  out.rho.resize(N);
  std::vector<double> kSorted(N);
  for (std::size_t i = 0; i < N; ++i) {
    out.r[i] = r[order[i]];
    out.u[i] = u[order[i]];
    out.t[i] = t[order[i]];
    out.rho[i] = rho[order[i]];
    kSorted[i] = k[order[i]];
  }
  std::vector<double> edges(N + 1);
  edges[0] = 0.0;
  edges[N] = RADIUS;
  for (std::size_t i = 1; i < N; ++i) edges[i] = 0.5 * (out.r[i - 1] + out.r[i]);
  constexpr double kPi = 3.14159265358979323846;
  double flow = 0.0;
  for (std::size_t i = 0; i < N; ++i) flow += out.u[i] * kPi * (edges[i + 1] * edges[i + 1] - edges[i] * edges[i]);
  out.tau = c.force * RADIUS / 2.0;
  out.q = (c.heat * kPi * RADIUS * RADIUS + c.force * flow) / (2.0 * kPi * RADIUS);
  out.pWall = out.rho[N - 1] * props().gasR * out.t[N - 1] + 2.0 / 3.0 * out.rho[N - 1] * kSorted[N - 1];
  out.tWall = c.twall;
  return out;
}

std::vector<double> parseTargets(const std::string& s) {
  std::vector<double> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) out.push_back(std::stod(tok));
  return out;
}

void apriori(const std::string& caseName, const std::vector<double>& targets) {
  const auto ref = references(caseName);
  const auto& P = props();
  const double rhoW = ref.pWall / (P.gasR * ref.tWall);
  const double uTau = std::sqrt(ref.tau / rhoW);
  const double muW = P.muAt(ref.tWall);
  std::vector<double> yplus(ref.r.size());
  for (std::size_t i = 0; i < ref.r.size(); ++i) yplus[i] = (RADIUS - ref.r[i]) * rhoW * uTau / muW;
  std::printf("%s: tau_w %.6f Pa, q_w %.6e W/m^2, p_wall %.1f Pa, Re_tau %.1f\n", caseName.c_str(), ref.tau, ref.q,
              ref.pWall, RADIUS * rhoW * uTau / muW);
  std::printf("%9s %8s %12s %12s %10s\n", "y+", "T1/Tw", "tau model-1", "q model-1", "ms/inv");
  for (const double target : targets) {
    std::size_t i = 0;
    double best = std::abs(yplus[0] - target);
    for (std::size_t j = 1; j < yplus.size(); ++j) {
      const double d = std::abs(yplus[j] - target);
      if (d < best) {
        best = d;
        i = j;
      }
    }
    const double y1 = RADIUS - ref.r[i];
    const auto t0 = std::chrono::steady_clock::now();
    try {
      const auto m = model(ref.u[i], ref.t[i], y1, ref.pWall, ref.tWall, ref.tau, ref.q);
      const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
      std::printf("%9.2f %8.4f %+12.4e %+12.4e %10.2f\n", yplus[i], ref.t[i] / ref.tWall, m.tauW / ref.tau - 1.0,
                  m.qW / ref.q - 1.0, ms);
    } catch (const std::exception& err) {
      std::printf("%9.2f %8.4f  no solution: %s\n", yplus[i], ref.t[i] / ref.tWall, err.what());
    }
  }
}

void pointCommand(int argc, char** argv) {
  if (argc < 9) throw std::runtime_error("point needs u1 T1 y1 p Tw tau0 q0");
  const double u1 = std::stod(argv[2]), t1 = std::stod(argv[3]), y1 = std::stod(argv[4]), p = std::stod(argv[5]),
               tw = std::stod(argv[6]), tau0 = std::stod(argv[7]), q0 = std::stod(argv[8]);
  const auto t0 = std::chrono::steady_clock::now();
  const auto m = model(u1, t1, y1, p, tw, tau0, q0);
  const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  std::printf("tau_w %.10g\n", m.tauW);
  std::printf("q_w %.10g\n", m.qW);
  std::printf("time_ms %.3f\n", ms);
}

}  // namespace

int main(int argc, char** argv) {
  std::setvbuf(stdout, nullptr, _IOLBF, 0);
  try {
    if (argc < 2) throw std::runtime_error("usage: crucible_wall_layer apriori|point ...");
    const std::string mode = argv[1];
    if (mode == "apriori") {
      const std::string caseName = argc > 2 ? argv[2] : "heated_hp";
      const std::string targetsStr = argc > 3 ? argv[3] : "1,5,30,100,300,1000,3000";
      apriori(caseName, parseTargets(targetsStr));
    } else if (mode == "point") {
      pointCommand(argc, argv);
    } else {
      throw std::runtime_error("unknown mode " + mode);
    }
  } catch (const std::exception& err) {
    std::fprintf(stderr, "crucible_wall_layer: %s\n", err.what());
    return 1;
  }
  return 0;
}
