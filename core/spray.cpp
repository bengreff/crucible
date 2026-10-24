#include "core/spray.hpp"
#include "core/amr.hpp"
#include "core/checkpoint.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <limits>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace crucible {
namespace {
constexpr double pi = std::numbers::pi;
constexpr int filmPasses = 3;  // of the film integrals against B_T: within 0.1% of converged in S2 (rates())
constexpr double tableTop = 5.01e6;  // Pa: the last subcritical row of the NIST table, rounded down
double sq(double x) { return x * x; }
// splitmix64 (Steele, Lea and Flood 2014): the stream's next 64 bits, and a number uniform in (0, 1).
std::uint64_t next(std::uint64_t& state) {
  std::uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}
double uniform(std::uint64_t& state) { return (static_cast<double>(next(state) >> 11) + 0.5) * 0x1.0p-53; }
// Two independent standard normals (Box-Muller).
std::array<double, 2> normals(std::uint64_t& state) {
  const double radius = std::sqrt(-2 * std::log(uniform(state))), angle = 2 * pi * uniform(state);
  return {radius * std::cos(angle), radius * std::sin(angle)};
}
// ln(1 + b) / b, 1 at b = 0.
double logRatio(double b) { return b > 1e-8 ? std::log1p(b) / b : 1 - b / 2; }
// Abramzon and Sirignano's film correction F(B) = (1 + B)^0.7 ln(1 + B) / B.
double filmCorrection(double b) { return std::pow(1 + b, 0.7) * logRatio(b); }
}  // namespace

// ---------------------------------------------------------------------------------------------------
// Liquid oxygen

LiquidOxygen LiquidOxygen::load(const std::string& path, double gasEnthalpy298) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Cannot open the liquid oxygen table " + path);
  LiquidOxygen out;
  const double offset = gasEnthalpy298 - idealGasEnthalpy298;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line.rfind("Temperature", 0) == 0) continue;
    if (line.find("undefined") != std::string::npos) continue;  // the critical point itself
    std::istringstream fields(line);
    std::vector<double> v;
    for (std::string f; std::getline(fields, f, '\t');) {
      if (f.empty()) continue;
      try { v.push_back(std::stod(f)); } catch (const std::exception&) { v.push_back(std::nan("")); }
    }
    if (v.size() < 25) throw std::runtime_error("Liquid oxygen table: short row in " + path);
    // Columns: T, P [MPa], rho_l, v_l, u_l, h_l [kJ/kg], s_l, cv_l, cp_l [J/(g K)], a_l, JT_l, mu_l,
    // k_l, sigma, then the vapour's rho, v, u, h, ...
    out.rows_.push_back({v[0], std::log(v[1] * 1e6), v[2], v[5] * 1e3 + offset, v[8] * 1e3, v[11], v[12], v[13],
                         v[17] * 1e3 + offset});
  }
  const auto n = out.rows_.size();
  if (n < 2) throw std::runtime_error("Liquid oxygen table: too few rows in " + path);
  out.t0_ = out.rows_.front()[0];
  out.step_ = (out.rows_.back()[0] - out.t0_) / static_cast<double>(n - 1);
  for (std::size_t i = 0; i < n; ++i) {
    for (double x : out.rows_[i])
      if (!std::isfinite(x)) throw std::runtime_error("Liquid oxygen table: non-finite value");
    if (std::abs(out.rows_[i][0] - (out.t0_ + out.step_ * static_cast<double>(i))) > 1e-6)
      throw std::runtime_error("Liquid oxygen table: rows are not uniform in T");
    if (i > 0 && !(out.rows_[i][3] > out.rows_[i - 1][3] && out.rows_[i][1] > out.rows_[i - 1][1]))
      throw std::runtime_error("Liquid oxygen table: enthalpy or vapour pressure not increasing");
  }
  return out;
}
std::array<double, 9> LiquidOxygen::row(double t) const {
  const double x = std::clamp((t - t0_) / step_, 0.0, static_cast<double>(rows_.size() - 1));
  const auto i = std::min(static_cast<std::size_t>(x), rows_.size() - 2);
  const double f = x - static_cast<double>(i);
  std::array<double, 9> out;
  for (std::size_t c = 0; c < 9; ++c) out[c] = rows_[i][c] + f * (rows_[i + 1][c] - rows_[i][c]);
  return out;
}
LiquidOxygen::State LiquidOxygen::at(double t) const {
  const auto r = row(t);
  return {r[2], r[3], r[4], r[5], r[6], r[7], std::exp(r[1]), r[8] - r[3]};
}
double LiquidOxygen::enthalpy(double t) const { return row(t)[3]; }
double LiquidOxygen::vapourPressureSlope(double t) const {
  const double x = std::clamp((t - t0_) / step_, 0.0, static_cast<double>(rows_.size() - 1));
  const auto i = std::min(static_cast<std::size_t>(x), rows_.size() - 2);
  return (rows_[i + 1][1] - rows_[i][1]) / step_;
}
namespace {
// T at which column c (increasing) of a uniform-T table takes value x, inverting the linear interpolant.
double invert(const std::vector<std::array<double, 9>>& rows, std::size_t c, double x, double t0, double step) {
  if (x <= rows.front()[c]) return t0;
  if (x >= rows.back()[c]) return t0 + step * static_cast<double>(rows.size() - 1);
  const auto it = std::upper_bound(rows.begin(), rows.end(), x,
                                   [c](double v, const std::array<double, 9>& r) { return v < r[c]; });
  const auto i = static_cast<std::size_t>(it - rows.begin()) - 1;
  const double f = (x - rows[i][c]) / (rows[i + 1][c] - rows[i][c]);
  return t0 + step * (static_cast<double>(i) + f);
}
}  // namespace
double LiquidOxygen::temperature(double h) const { return invert(rows_, 3, h, t0_, step_); }
double LiquidOxygen::saturationTemperature(double p) const { return invert(rows_, 1, std::log(p), t0_, step_); }

double Parcel::dropMass(double density) const { return density * pi / 6 * d * d * d; }

// ---------------------------------------------------------------------------------------------------
// Spray

double Spray::Injector::opening(double time) const {
  if (time < opens) return 0;
  return ramp > 0 ? std::min(1.0, (time - opens) / ramp) : 1.0;
}

Spray::Spray(Flow& flow, LiquidOxygen liquid, Settings settings)
    : flow_(flow), liquid_(std::move(liquid)), settings_(std::move(settings)) {
  const auto& medium = flow_.medium();
  const auto& species = medium.species();
  const auto o2 = std::find_if(species.begin(), species.end(), [](const Species& s) { return s.name == "O2"; });
  if (o2 == species.end()) throw std::invalid_argument("The spray needs O2 in the medium.");
  oxygen_ = static_cast<std::size_t>(o2 - species.begin());
  if (!medium.hasTransport()) throw std::invalid_argument("The spray's film model needs the medium's transport.");
  const auto& fits = flow_.definition().transport;
  const auto ns = species.size();
  conductivityFit_ = fits.conductivity;
  for (std::size_t k = 0; k < ns; ++k) {  // pairs k <= j in row order: (k, j) at k ns - k (k - 1) / 2 + j - k
    const auto lo = std::min(k, oxygen_), hi = std::max(k, oxygen_);
    oxygenDiffusionFit_.push_back(fits.diffusion[lo * ns - lo * (lo - 1) / 2 + hi - lo]);
  }
  if (settings_.filmPoints < 0 || settings_.filmPoints > 32) throw std::invalid_argument("Spray: filmPoints in 0 to 32.");
  // Gauss-Legendre nodes by Newton's method on P_n, mapped to [0, 1].
  const int points = settings_.filmPoints;
  for (int i = 0; i < points; ++i) {
    double x = std::cos(pi * (i + 0.75) / (points + 0.5)), dp = 1;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1, p1 = x;
      for (int k = 2; k <= points; ++k) {
        const double p2 = ((2 * k - 1) * x * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      if (points == 1) p0 = 1;
      dp = points * (x * p1 - p0) / (x * x - 1);
      const double dx = p1 / dp;
      x -= dx;
      if (std::abs(dx) <= 1e-16) break;
    }
    gaussNode_.push_back(0.5 * (1 - x));
    gaussWeight_.push_back(1 / ((1 - x * x) * dp * dp));
  }
  // gaussIntegral_[i n + j] = int_0^{x_i} l_j, l_j the Lagrange basis on the nodes (exact by the same rule).
  const auto n = gaussNode_.size();
  gaussIntegral_.assign(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      for (std::size_t m = 0; m < n; ++m) {
        const double x = gaussNode_[i] * gaussNode_[m];
        double l = 1;
        for (std::size_t q = 0; q < n; ++q)
          if (q != j) l *= (x - gaussNode_[q]) / (gaussNode_[j] - gaussNode_[q]);
        gaussIntegral_[i * n + j] += gaussNode_[i] * gaussWeight_[m] * l;
      }
  for (const auto& j : settings_.injectors) {
    if (!(j.rOuter > j.rInner) || !(j.rInner >= 0) || !(j.massFlow >= 0) || !(j.parcelRate > 0) || !(j.z >= 0))
      throw std::invalid_argument("Spray injector: needs rOuter > rInner >= 0, massFlow >= 0, parcelRate > 0, z >= 0.");
    if (!(j.temperature >= liquid_.minTemperature() && j.temperature <= liquid_.maxTemperature()))
      throw std::invalid_argument("Spray injector: temperature outside the liquid table.");
    if (j.size == Injector::Size::RosinRammler && !(j.sauterDiameter > 0 && j.exponent > 1))
      throw std::invalid_argument("Spray injector: Rosin-Rammler needs a positive SMD and exponent above 1.");
  }
  if (!(settings_.restitution >= 0 && settings_.restitution <= 1)) throw std::invalid_argument("Restitution outside [0, 1].");
  levels_ = {&flow_};
  setCells();
  work_.resize(static_cast<std::size_t>(flow_.threads()));
  for (auto& w : work_) w.resize(4 * ns + 8);
  std::uint64_t seed = settings_.seed;
  for (std::size_t k = 0; k < settings_.injectors.size(); ++k) {
    injectorRandom_.push_back(seed + 0x632be59bd9b4e019ULL * (k + 1));
    pending_.push_back({});
  }
  unit_.assign(ns, 0);
  unit_[oxygen_] = 1;
}

Spray::Spray(Amr& amr, LiquidOxygen liquid, Settings settings) : Spray(amr.level(0), std::move(liquid), std::move(settings)) {
  amr_ = &amr;
  levels_.clear();
  for (int l = 0; l < amr.levels(); ++l) levels_.push_back(&amr.level(l));
  setCells();
}

// The composite per-cell arrays over levels_.
void Spray::setCells() {
  const auto ns = flow_.medium().size();
  offset_.assign(1, 0);
  for (const Flow* f : levels_) offset_.push_back(offset_.back() + f->mesh().cells.size());
  const auto cells = offset_.back();
  cellGas_.assign(cells * (6 + ns), 0);
  for (std::size_t q = 0; q < cells; ++q) cellGas_[q * (6 + ns)] = 300;  // temperature guess
  cellLoad_.assign(cells * 5, 0);
  occupied_.assign(cells, 0);
  if (!flow_.turbulence().empty()) cellTurbulence_.assign(cells * 2, 0);
  cellSize_.resize(cells);
  cellVolume_.resize(cells);
  amounts_.clear();
  for (std::size_t l = 0; l < levels_.size(); ++l) {
    const auto& mesh = levels_[l]->mesh();
    const auto at = offset_[l];
    for (int i = 0; i < mesh.nz; ++i)
      for (int jr = 0; jr < mesh.nr; ++jr) {
        const double width = 0.5 * (mesh.radius[i] + mesh.radius[i + 1]) * (mesh.fraction[jr + 1] - mesh.fraction[jr]);
        cellSize_[at + mesh.index(i, jr)] = std::min(mesh.dz, width);
      }
    for (std::size_t q = 0; q < mesh.cells.size(); ++q) cellVolume_[at + q] = mesh.cells[q].volume;
    amounts_.emplace_back(mesh.cells.size() * (ns + 3), 0.0);
  }
}
int Spray::levelOf(std::size_t cell) const {
  int l = 0;
  while (cell >= offset_[static_cast<std::size_t>(l) + 1]) ++l;
  return l;
}
// A coarse cell is the exact union of its children (Mesh::coarsened), so level l's cell over the finest
// cell (i, j) is (i, j) shifted right by the levels between; the owner is the highest stepping level active there.
int Spray::owner(int i, int j) const {
  if (!amr_) return static_cast<int>(finest().index(i, j));
  const int f = static_cast<int>(levels_.size()) - 1;
  for (int l = amr_->top(); l > 0; --l) {
    const auto q = levels_[static_cast<std::size_t>(l)]->mesh().index(i >> (f - l), j >> (f - l));
    if (amr_->active(l)[q]) return static_cast<int>(offset_[static_cast<std::size_t>(l)] + q);
  }
  return static_cast<int>(levels_[0]->mesh().index(i >> f, j >> f));
}
// A covered cell resolves into its owned children, a cell that is not active into the owned cell below holding it.
void Spray::mark(int l, int i, int j) {
  const auto& mesh = levels_[static_cast<std::size_t>(l)]->mesh();
  const auto q = mesh.index(i, j);
  if (amr_) {
    if (l < amr_->top() && amr_->active(l + 1)[levels_[static_cast<std::size_t>(l) + 1]->mesh().index(2 * i, 2 * j)]) {
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) mark(l + 1, 2 * i + a, 2 * j + b);
      return;
    }
    if (l > 0 && !amr_->active(l)[q]) { mark(l - 1, i >> 1, j >> 1); return; }
  }
  occupied_[offset_[static_cast<std::size_t>(l)] + q] = 1;
}

// The wall radius at z, linear in each column, and its slope dR/dz.
double Spray::wallRadius(double z, double* slope) const {
  const auto& mesh = finest();
  const int i = std::clamp(static_cast<int>(z / mesh.dz), 0, mesh.nz - 1);
  const double s = (mesh.radius[i + 1] - mesh.radius[i]) / mesh.dz;
  if (slope) *slope = s;
  return mesh.radius[i] + (z - mesh.dz * i) * s;
}
int Spray::locate(double z, double r) const {
  const auto& mesh = finest();
  const double span = mesh.dz * mesh.nz;
  if (!(z >= 0 && z < span && r >= 0)) return -1;
  const int i = std::min(static_cast<int>(z / mesh.dz), mesh.nz - 1);
  const double wall = wallRadius(z, nullptr);
  if (r > wall) return -1;
  const double s = r / wall;
  const auto it = std::upper_bound(mesh.fraction.begin(), mesh.fraction.end(), s);
  const int j = std::clamp(static_cast<int>(it - mesh.fraction.begin()) - 1, 0, mesh.nr - 1);
  return owner(i, j);
}

void Spray::add(Parcel parcel) {
  parcel.cell = locate(parcel.z, parcel.r);
  if (parcel.cell < 0) throw std::invalid_argument("Parcel outside the domain.");
  if (parcel.random == 0) parcel.random = settings_.seed ^ (0xd1b54a32d192ed03ULL * ++nextId_);
  const auto s = liquid_.at(parcel.t);
  const double m = parcel.count * parcel.dropMass(s.density);
  if (!(parcel.injectedMass > 0)) parcel.injectedMass = m;
  parcels_.push_back(parcel);
  ledger_.injectedMass += m;
  ledger_.injectedEnergy += m * (s.enthalpy + 0.5 * (sq(parcel.uz) + sq(parcel.ur) + sq(parcel.ut)));
  ledger_.injectedMomentum += m * parcel.uz;
}

// The gas each occupied cell offers its parcels: T, p, rho, u_z, u_r, c_v, then the mass fractions.
void Spray::gatherGas() {
  const auto& medium = flow_.medium();
  const auto ns = medium.size();
  const std::size_t nt = flow_.turbulence().empty() ? 0 : 2;
  // The moving parcels' cells and their neighbours (on AMR, the owned cells over the 3 x 3 block of
  // the owner level's cells): in one step a parcel moves less than a cell (its speed is below the gas's
  // u + c, which bounds the step), so it stays among these; moveParcel checks.
  std::fill(occupied_.begin(), occupied_.end(), 0);
  for (const auto n : moved_) {
    const auto cell = static_cast<std::size_t>(parcels_[n].cell);
    const int l = levelOf(cell);
    const auto& mesh = levels_[static_cast<std::size_t>(l)]->mesh();
    const auto q = static_cast<int>(cell - offset_[static_cast<std::size_t>(l)]);
    const int i = q / mesh.nr, j = q % mesh.nr;
    for (int a = std::max(i - 1, 0); a <= std::min(i + 1, mesh.nz - 1); ++a)
      for (int b = std::max(j - 1, 0); b <= std::min(j + 1, mesh.nr - 1); ++b) mark(l, a, b);
  }
  flow_.pool().blocks(occupied_.size(), 256, [&](std::size_t begin, std::size_t end, int) {
    for (std::size_t c = begin; c < end; ++c) {
      if (!occupied_[c]) continue;
      const int l = levelOf(c);
      const Flow& level = *levels_[static_cast<std::size_t>(l)];
      const auto& state = level.state();
      const auto& partial = level.partialDensities();
      const auto& turbulence = level.turbulence();
      const std::size_t q = c - offset_[static_cast<std::size_t>(l)];
      double* g = cellGas_.data() + c * (6 + ns);
      double* y = g + 6;
      double sum = 0;
      for (std::size_t k = 0; k < ns; ++k) sum += std::max(partial[q * ns + k], 0.0);
      for (std::size_t k = 0; k < ns; ++k) y[k] = std::max(partial[q * ns + k], 0.0) / sum;
      const auto& u = state[q];
      const double k = nt ? turbulence[q * nt] / u[0] : 0;
      const double e = (u[3] - (sq(u[1]) + sq(u[2])) / (2 * u[0])) / u[0] - k;
      const double t = medium.temperature(e, y, g[0]);
      if (!std::isfinite(t)) throw std::runtime_error("Spray: non-admissible gas in an occupied cell.");
      g[0] = t;
      g[1] = u[0] * medium.gasConstant(y) * t;
      g[2] = u[0];
      g[3] = u[1] / u[0];
      g[4] = u[2] / u[0];
      g[5] = medium.cv(t, y);
      if (nt) {
        cellTurbulence_[c * 2] = k;
        cellTurbulence_[c * 2 + 1] = turbulence[q * nt + 1] / u[0];
      }
      if (g[1] > tableTop) throw std::runtime_error("Spray: a parcel's cell is above oxygen's critical pressure.");
    }
  });
}

Spray::FilmTransport Spray::filmTransport(double t, double p, const double* y) const {
  const auto& species = flow_.medium().species();
  const auto n = species.size();
  const double logt = std::log(t), sqrtT = std::sqrt(t);
  auto poly = [&](const std::array<double, 5>& c) {
    return c[0] + logt * (c[1] + logt * (c[2] + logt * (c[3] + logt * c[4])));
  };
  double moles = 0;
  for (std::size_t k = 0; k < n; ++k) moles += y[k] / species[k].molarMass;
  const double meanMass = 1 / moles;
  double series = 0, parallel = 0, sum = 0, others = 0;
  for (std::size_t k = 0; k < n; ++k) {
    if (!(y[k] > 0)) continue;
    const double x = std::max(1e-20, y[k] / species[k].molarMass * meanMass);
    const double lambda = sqrtT * poly(conductivityFit_[k]);
    series += x * lambda;
    parallel += x / lambda;
    if (k != oxygen_) {
      sum += x / (t * sqrtT * poly(oxygenDiffusionFit_[k]));
      others += y[k];
    }
  }
  FilmTransport out;
  out.conductivity = 0.5 * (series + 1 / parallel);
  out.diffusion = n == 1 ? t * sqrtT * poly(oxygenDiffusionFit_[oxygen_]) / p : sum > 0 ? others / (p * sum) : 0;
  out.r = Medium::universalGasConstant * moles;
  // A fit taken below its range turns negative (H2's conductivity below about 97 K on fits from 200 K or more,
  // thermo::transportFits); the film then has a negative phi and its shape leaves [0, 1].
  if (!(out.conductivity > 0) || !(out.diffusion >= 0) || !std::isfinite(out.diffusion))
    throw std::runtime_error("Spray: the film's transport is not positive at " + std::to_string(t) +
                             " K: the transport fits do not reach the film's temperature (thermo::transportFits' minTemperature)");
  return out;
}

// The molar mass of the gas's non-O2 part [kg/kmol], 0 if it has none.
double Spray::carrierMass(const double* y) const {
  const auto& species = flow_.medium().species();
  double mass = 0, moles = 0;
  for (std::size_t k = 0; k < species.size(); ++k)
    if (k != oxygen_) {
      mass += y[k];
      moles += y[k] / species[k].molarMass;
    }
  return moles > 0 ? mass / moles : 0;
}

// The stagnant film (Re = 0) between the surface and the free stream, with variable properties. With
// eta = -ln(1 - Y_O2), O2 the only species with a net flux and Fick's law on its Y with the Blanc's-law
// D, mdot / (4 pi a) = int rho D d eta and ln(1 + B_T) = int phi d eta along the film's (T, Y) path, where
// phi = c_pF rho D / lambda, and theta = T - T_s + L_eff / c_pF falls as e^(-int phi d eta). With
// constant properties these are A-S's rho D ln(1 + B_M) and phi ln(1 + B_M). So the film's rho D is the
// mean of rho D over eta, and its lambda is c_pF times that over the mean of phi. On s = eta / aM in
// [0, 1], aM = ln(1 + B_M), Y is exact by construction, 1 - Y_O2 = (1 - Y_inf) e^(-aM s), and T follows
// T - T_s = (T_inf - T_s) expm1(aT (1 - F(s))) / expm1(aT), aT = ln(1 + B_T), with F = int_0^s phi / int_0^1 phi
// (s for constant phi). shape[] holds F at the nodes on entry and this pass's F on return, so passes
// converge on the path's own shape. The other species keep the free stream's proportions. Gauss-Legendre
// in s, on [0, 8 / aM] and [8 / aM, 1] when aM > 8 (the carrier gas then lies in a layer near the free
// stream). Boiling (aM infinite): a pure O2 film, where the exact form is lambda = 1 / int dxi / lambda
// over xi = a / r, with T on the constant-phi shape. SPRAY.md, Results: S2.
std::array<double, 2> Spray::filmConductances(const Gas& gas, double ts, double aM, double aT, double cpf, double* y,
                                              double* shape) const {
  const auto ns = flow_.medium().size(), n = gaussNode_.size();
  const double yInf = gas.y[oxygen_];
  const bool pure = !std::isfinite(aM) || !(yInf < 1);
  const double wOther = carrierMass(gas.y);
  thread_local std::vector<double> phi;
  phi.resize(2 * n);
  double sumRhoD = 0, sumPhi = 0, resistance = 0;
  auto piece = [&](double x0, double x1, std::size_t first) {
    for (std::size_t i = 0; i < n; ++i) {
      const double s = x0 + (x1 - x0) * gaussNode_[i], w = (x1 - x0) * gaussWeight_[i];
      const double carrier = pure ? 0 : (1 - yInf) * std::exp(-aM * s);  // 1 - Y_O2
      for (std::size_t k = 0; k < ns; ++k)
        y[k] = k == oxygen_ ? 1 - carrier : (pure ? 0 : carrier * gas.y[k] / (1 - yInf));
      const double f = pure ? s : shape[first + i];
      const double t = std::abs(aT) > 1e-12 ? ts + (gas.t - ts) * std::expm1(aT * (1 - f)) / std::expm1(aT)
                                            : ts + (gas.t - ts) * (1 - f);
      const auto film = filmTransport(t, gas.p, y);
      const double rhoD = gas.p / (film.r * t) * film.diffusion * (wOther > 0 ? Medium::universalGasConstant / (film.r * wOther) : 1);
      phi[first + i] = cpf * rhoD / film.conductivity;
      sumRhoD += w * rhoD;
      sumPhi += w * phi[first + i];
      resistance += w / film.conductivity;
    }
  };
  const double split = !pure && aM > 8 ? 8 / aM : 1;
  piece(0, split, 0);
  if (split < 1) piece(split, 1, n);
  if (pure) return {1 / resistance, sumRhoD};
  // F at the nodes from this pass's phi: the interpolant's integral within each piece, plus the first
  // piece's whole integral for the second.
  double before = 0;
  for (std::size_t p = 0, first = 0; p < (split < 1 ? 2u : 1u); ++p, first += n) {
    const double width = p == 0 ? split : 1 - split;
    double whole = 0;
    for (std::size_t i = 0; i < n; ++i) {
      double integral = 0;
      for (std::size_t j = 0; j < n; ++j) integral += gaussIntegral_[i * n + j] * phi[first + j];
      shape[first + i] = (before + width * integral) / sumPhi;
      whole += gaussWeight_[i] * phi[first + i];
    }
    before += width * whole;
  }
  return {cpf * sumRhoD / sumPhi, sumRhoD};
}

Spray::Rates Spray::rates(const Parcel& parcel, const Gas& gas, std::vector<double>& work) const {
  const auto& medium = flow_.medium();
  const auto ns = medium.size();
  if (gas.p > tableTop) throw std::runtime_error("Spray: gas above oxygen's critical pressure.");
  Rates out;
  const auto liquid = liquid_.at(parcel.t);
  const double tSat = liquid_.saturationTemperature(gas.p);
  out.boiling = parcel.t >= tSat || liquid.pressure >= gas.p;
  // Surface composition: O2 at Raoult's mole fraction, the other species in the free stream's proportions.
  double* ys = work.data();
  double* yr = ys + ns;
  double* diffusion = yr + ns;
  const double yInf = gas.y[oxygen_];
  double ySurface = 1, slope = 0;  // slope: dY_s/dT
  const double wOther = carrierMass(gas.y);
  if (!out.boiling) {
    const double xs = liquid.pressure / gas.p;
    if (wOther > 0) {
      const double wO2 = medium.species()[oxygen_].molarMass;
      const double mix = xs * wO2 + (1 - xs) * wOther;
      ySurface = xs * wO2 / mix;
      slope = wO2 * wOther / sq(mix) * xs * liquid_.vapourPressureSlope(parcel.t);
    }
    // Within rounding of saturation (Y_s rounds to 1, B_M is infinite): boiling.
    if (!(ySurface < 1)) {
      out.boiling = true;
      ySurface = 1;
      slope = 0;
    }
  }
  const double ts = out.boiling ? tSat : parcel.t;
  for (std::size_t k = 0; k < ns; ++k)
    ys[k] = k == oxygen_ ? ySurface : (yInf < 1 ? (1 - ySurface) * gas.y[k] / (1 - yInf) : 0);
  // The 1/3 rule.
  const double tr = ts + (gas.t - ts) / 3;
  for (std::size_t k = 0; k < ns; ++k) yr[k] = ys[k] + (gas.y[k] - ys[k]) / 3;
  thread_local std::vector<double> scratch;
  const auto transport = medium.transport(tr, gas.p, yr, diffusion, scratch);
  const auto film = medium.properties(tr, yr);
  const auto vapour = medium.properties(tr, unit_.data());
  const double cpg = film.cv + film.r, cpf = vapour.cv + vapour.r;
  const double rhoFilm = gas.p / (film.r * tr);
  // The 1/3 rule's properties set the convective groups and the drag; the film integrals replace its
  // conductances lambda and rho D (filmConductances) unless filmPoints is 0.
  // A-S's rho D belongs to Fick's law on Y_O2. Medium's D_km belongs to gradients of X_k; for O2 against
  // carriers in fixed proportions (molar mass W_c) the two differ by W / W_c (the Blanc's-law binary D).
  double rhoD = rhoFilm * diffusion[oxygen_] * (wOther > 0 ? Medium::universalGasConstant / (film.r * wOther) : 1);
  const double mu = transport.viscosity;
  double lambda = transport.conductivity;
  const double slip = std::sqrt(sq(gas.uz - parcel.uz) + sq(gas.ur - parcel.ur) + sq(gas.ut - parcel.ut));
  out.reynolds = gas.rho * slip * parcel.d / mu;
  const double pr = cpg * mu / lambda, sc = mu / rhoD;
  const double f = out.reynolds > 1 ? std::pow(out.reynolds, 0.077) : 1.0;
  const double sh0 = 1 + std::cbrt(1 + out.reynolds * sc) * f;
  const double nu0 = 1 + std::cbrt(1 + out.reynolds * pr) * f;
  const bool integrated = !gaussNode_.empty();
  double bt = 0;
  if (out.boiling) {
    // Held at saturation: all the heat that arrives evaporates liquid.
    bt = gas.t > ts ? cpf * (gas.t - ts) / liquid_.at(ts).latent : 0;
    if (integrated) {
      const auto c = filmConductances(gas, ts, std::numeric_limits<double>::infinity(), std::log1p(bt), cpf, yr, nullptr);
      lambda = c[0];
      rhoD = c[1];
    }
    out.nusselt = 2 + (nu0 - 2) / filmCorrection(bt);
    out.sherwood = out.nusselt;
    out.mass = pi * parcel.d * lambda / cpf * out.nusselt * std::log1p(bt);
  } else {
    if (ySurface - yInf < 0) out.clamped = true;
    const double bm = std::max(0.0, (ySurface - yInf) / (1 - ySurface));
    out.spalding = bm;
    out.sherwood = 2 + (sh0 - 2) / filmCorrection(bm);
    // B_T = (1 + B_M)^phi - 1 with phi = c_pF rho D Sh* / (lambda Nu*(B_T)): secant steps on the residual
    // from a fixed-point step off `start` (A-S's own start is B_M; a fixed-point step whenever the secant
    // leaves B > 0), and bisection on the residual when that does not settle within 30 iterations.
    const double lnBm = std::log1p(bm);
    auto solve = [&](double start) {
      auto residual = [&](double b) {
        const double phi = cpf * rhoD * out.sherwood / (lambda * (2 + (nu0 - 2) / filmCorrection(b)));
        const double r = std::expm1(phi * lnBm) - b;
        if (std::isnan(r)) throw std::runtime_error("Spray: the heat-transfer number's residual is NaN.");
        return r;
      };
      double b = start;
      double rt = residual(b);
      bool settled = rt == 0;
      double bPrev = b, rPrev = rt;
      if (!settled) b += rt;
      for (int it = 0; it < 30 && !settled && std::isfinite(b); ++it) {
        rt = residual(b);
        if (rt == 0) {
          settled = true;
          break;
        }
        double next = b - rt * (b - bPrev) / (rt - rPrev);
        if (!(next > 0) || !std::isfinite(next)) next = b + rt;
        bPrev = b;
        rPrev = rt;
        // The step after it lands at rounding (superlinear). Near saturation (B_M about 1e8) the residual of an
        // iterate below the root can pass the double range, so the step is infinite, or be so large that the next
        // secant step is tiny far from the root; neither settles (SPEED.md, section 6: the spray NaN).
        settled = std::isfinite(next) && std::abs(next - b) <= 1e-11 * next && std::abs(rt) <= 1e-6 * next;
        b = next;
      }
      if (!settled) {
        double lo = 0, hi = std::max(1.0, 2 * (std::isfinite(b) ? b : bPrev));
        while (residual(hi) > 0) hi *= 2;
        for (int it = 0; it < 2000; ++it) {
          const double mid = 0.5 * (lo + hi);
          if (mid <= lo || mid >= hi) break;
          (residual(mid) > 0 ? lo : hi) = mid;
        }
        b = std::abs(residual(lo)) < std::abs(residual(hi)) ? lo : hi;
        if (!std::isfinite(b)) throw std::runtime_error("Spray: the heat-transfer number did not converge.");
      }
      return b;
    };
    bt = solve(bm);
    if (integrated) {
      // The film's path takes B_T and its shape from the last pass, from the constant-phi shape. A fixed
      // number of passes keeps the rates a smooth function of the drop and the gas.
      const double aM = std::log((1 - yInf) / (1 - ySurface));  // ln(1 + B_M) before the clamp
      const double split = aM > 8 ? 8 / aM : 1;
      const auto n = gaussNode_.size();
      thread_local std::vector<double> shape;
      shape.resize(2 * n);
      for (std::size_t i = 0; i < n; ++i) {
        shape[i] = split * gaussNode_[i];
        shape[n + i] = split + (1 - split) * gaussNode_[i];
      }
      for (int pass = 0; pass < filmPasses; ++pass) {
        const auto c = filmConductances(gas, ts, aM, std::log1p(bt), cpf, yr, shape.data());
        lambda = c[0];
        rhoD = c[1];
        bt = solve(bt);
      }
    }
    out.nusselt = 2 + (nu0 - 2) / filmCorrection(bt);
    out.mass = pi * parcel.d * rhoD * out.sherwood * std::log1p(bm);
    // d(mdot)/dT through Y_s alone: pi d rho D Sh* (dY_s/dT) / (1 - Y_s); it sets the stiffness of the
    // drop's temperature near its wet bulb, where Y_s nears 1.
    if (ySurface > yInf) out.stiffness = liquid.latent * pi * parcel.d * rhoD * out.sherwood * slope / (1 - ySurface);
  }
  // Conduction to the surface, Q = mdot c_pF (T_inf - T_s) / B_T = pi d lambda Nu* (T_inf - T_s) ln(1 + B_T) / B_T.
  out.conductance = pi * parcel.d * lambda * out.nusselt * logRatio(bt);
  out.stiffness += out.conductance;
  out.surface = ts;
  out.heat = out.conductance * (gas.t - ts);
  // Boiling: everything that arrives evaporates (the same value as the line above, written exactly).
  if (out.boiling && gas.t > ts) out.heat = out.mass * liquid_.at(ts).latent;
  out.transfer = bt;
  // Drag: the momentum relaxation time tau = rho_l d^2 / (18 mu g(Re, B)), C_D = (24 / Re) g.
  out.tau = liquid.density * sq(parcel.d) / (18 * mu * dragFactor(out.reynolds, bt));
  out.dragCoefficient = 24 * dragFactor(out.reynolds, bt) / std::max(out.reynolds, 1e-300);
  return out;
}

// Renksizbulut and Yuen (1983): C_D = (24 / Re)(1 + 0.2 Re^0.63) / (1 + B)^0.2 (see SPRAY.md, Literature).
double Spray::dragFactor(double reynolds, double transfer) {
  return (1 + 0.2 * std::pow(reynolds, 0.63)) / std::pow(1 + transfer, 0.2);
}

Spray::Waves Spray::waves(double d, double u, double rhoGas, double rhoLiquid, double muLiquid, double tension,
                          double deceleration, double cRt) {
  const double a = d / 2;
  const double weGas = rhoGas * u * u * a / tension, weLiquid = rhoLiquid * u * u * a / tension;
  const double reLiquid = rhoLiquid * u * a / muLiquid;
  const double z = std::sqrt(weLiquid) / reLiquid, t = z * std::sqrt(weGas);
  Waves w;
  w.khLength = 9.02 * a * (1 + 0.45 * std::sqrt(z)) * (1 + 0.4 * std::pow(t, 0.7)) /
               std::pow(1 + 0.87 * std::pow(weGas, 1.67), 0.6);
  w.khRate = (0.34 + 0.38 * std::pow(weGas, 1.5)) / ((1 + z) * (1 + 1.4 * std::pow(t, 0.6))) /
             std::sqrt(rhoLiquid * a * a * a / tension);
  const double drive = deceleration * (rhoLiquid - rhoGas);
  w.rtLength = 2 * pi * cRt * std::sqrt(3 * tension / drive);
  w.rtRate = std::sqrt(2 / (3 * std::sqrt(3 * tension)) * std::pow(drive, 1.5) / (rhoLiquid + rhoGas));
  return w;
}

Spray::Gas Spray::seen(int cell) const {
  const double* g = cellGas_.data() + static_cast<std::size_t>(cell) * (6 + flow_.medium().size());
  return {g[0], g[1], g[2], g[3], g[4], g + 6};
}
Spray::Gas Spray::seenBy(const Parcel& parcel) const {
  Gas gas = seen(parcel.cell);
  gas.uz += parcel.eddyZ;
  gas.ur += parcel.eddyR;
  gas.ut = parcel.eddyT;
  return gas;
}

// Eddy lifetime l_e / sqrt(2k/3) = 3 / (Sc_t omega) and crossing time -tau ln(1 - l_e / (tau slip)), the
// shorter (SPRAY.md, Turbulent dispersion).
double Spray::interaction(double k, double omega, double tau, double slip, double schmidt) {
  const double length = 3 / schmidt * std::sqrt(2 * k / 3) / omega, lifetime = 3 / (schmidt * omega);
  const double reach = tau * slip;
  return length >= reach ? lifetime : std::min(lifetime, -tau * std::log1p(-length / reach));
}
double Spray::blobDiameter(int injector) const {
  const auto& j = settings_.injectors[static_cast<std::size_t>(injector)];
  return 2 * (j.rOuter - j.rInner);
}

// The point-implicit seen state (SPRAY.md, Coupling): each parcel's drag and heat conductances at the
// cell's own gas, summed into cells in parcel order; then each loaded cell's velocity and temperature
// become the backward-Euler end state under those loads.
void Spray::applyLoads(double dt) {
  const auto ns = flow_.medium().size();
  flow_.pool().blocks(moved_.size(), 64, [&](std::size_t begin, std::size_t end, int worker) {
    for (std::size_t n = begin; n < end; ++n) {
      const auto i = moved_[n];
      const auto& p = parcels_[i];
      const auto r = rates(p, seenBy(p), work_[static_cast<std::size_t>(worker)]);
      slots_[i].drag = p.count * p.dropMass(liquid_.at(p.t).density) / r.tau;
      slots_[i].conductance = p.count * r.conductance;
      slots_[i].surface = r.surface;
    }
  });
  std::fill(cellLoad_.begin(), cellLoad_.end(), 0.0);
  for (const auto i : moved_) {
    const auto& p = parcels_[i];
    double* l = cellLoad_.data() + static_cast<std::size_t>(p.cell) * 5;
    l[0] += slots_[i].drag;
    l[1] += slots_[i].drag * (p.uz - p.eddyZ);  // the gas feels the drag toward its velocity plus the eddy's
    l[2] += slots_[i].drag * (p.ur - p.eddyR);
    l[3] += slots_[i].conductance;
    l[4] += slots_[i].conductance * slots_[i].surface;
  }
  for (std::size_t q = 0; q < cellVolume_.size(); ++q) {
    const double* l = cellLoad_.data() + q * 5;
    if (l[0] == 0 && l[3] == 0) continue;
    double* g = cellGas_.data() + q * (6 + ns);
    const double gasMass = g[2] * cellVolume_[q];
    g[3] = (gasMass * g[3] + dt * l[1]) / (gasMass + dt * l[0]);
    g[4] = (gasMass * g[4] + dt * l[2]) / (gasMass + dt * l[0]);
    const double capacity = gasMass * g[5];
    const double t = (capacity * g[0] + dt * l[4]) / (capacity + dt * l[3]);
    g[1] *= t / g[0];
    g[0] = t;
  }
}

// One parcel across dt in sub-steps (SPRAY.md, Integration). Each sub-step's change of the parcel's
// mass, momentum and energy from evaporation, heating and drag goes to the gas of the cell it started
// in, as minus that change; walls, the outlet and breakup are booked separately.
void Spray::moveParcel(std::size_t index, double dt, int worker) {
  Slot& slot = slots_[index];
  slot.deposits.clear();
  slot.children.clear();
  slot.wallImpulse = slot.wallLoss = slot.swirl = slot.flashed = slot.outMass = slot.outEnergy = slot.outMomentum = 0;
  slot.alive = true;
  slot.substeps = slot.wallHits = slot.flashes = slot.vaporised = slot.clamped = slot.finer = slot.coarser = 0;
  Parcel& p = parcels_[index];
  auto& work = work_[static_cast<std::size_t>(worker)];
  const auto& mesh = finest();
  const double span = mesh.dz * mesh.nz;
  const auto& breakup = settings_.breakup;
  const double e = settings_.restitution;
  const double minMass = liquid_.at(p.t).density * pi / 6 * std::pow(settings_.minDiameter, 3);
  auto deposit = [&](int cell, double mass, double pz, double pr, double energy) {
    if (!slot.deposits.empty() && slot.deposits.back().cell == cell) {
      auto& d = slot.deposits.back();
      d.mass += mass; d.pz += pz; d.pr += pr; d.energy += energy;
    } else {
      slot.deposits.push_back({cell, mass, pz, pr, energy});
    }
  };
  // The parcel's own state across the sub-steps: mass and enthalpy per drop.
  double m = p.dropMass(liquid_.at(p.t).density), h = liquid_.enthalpy(p.t);
  const bool disperse = settings_.dispersion && !cellTurbulence_.empty();
  const double schmidt = flow_.definition().turbulence.schmidt;
  double remaining = dt;
  while (remaining > 0) {
    if (!occupied_[static_cast<std::size_t>(p.cell)])
      throw std::runtime_error("Spray: a parcel crossed more than one cell in a step.");
    const int cell = p.cell;
    // A new eddy where the last interaction ended: u' from the cell's k, its time set below.
    bool newEddy = false;
    if (disperse && p.eddyTime <= 0) {
      const double k = cellTurbulence_[static_cast<std::size_t>(cell) * 2];
      const double sigma = k > 0 ? std::sqrt(2 * k / 3) : 0;
      const auto a = normals(p.random), b = normals(p.random);
      p.eddyZ = sigma * a[0];
      p.eddyR = sigma * a[1];
      p.eddyT = sigma * b[0];
      newEddy = k > 0 && cellTurbulence_[static_cast<std::size_t>(cell) * 2 + 1] > 0;
      if (!newEddy) p.eddyZ = p.eddyR = p.eddyT = 0;
    }
    const Gas gas = seenBy(p);
    ++slot.substeps;
    // Start-of-sub-step totals (shed mass travels with the drop at its temperature and velocity).
    const double ke0 = 0.5 * (sq(p.uz) + sq(p.ur) + sq(p.ut));
    const double mass0 = p.count * (m + p.shed);
    const double pz0 = mass0 * p.uz, px0 = mass0 * p.ur, py0 = mass0 * p.ut, e0 = mass0 * (h + ke0);
    // Below the smallest followed size, or flashed whole: the liquid becomes vapour at once.
    auto vaporise = [&]() {
      deposit(cell, mass0, pz0, px0, e0);
      slot.swirl += std::abs(py0);
      ++slot.vaporised;
      slot.alive = false;
    };
    if (m < minMass) { vaporise(); return; }
    // A drop above saturation (it reached lower pressure) flashes the enthalpy above saturation.
    const double tSat = liquid_.saturationTemperature(gas.p);
    if (p.t > tSat) {
      const double hSat = liquid_.enthalpy(tSat);
      const double f = (h - hSat) / liquid_.at(tSat).latent;
      ++slot.flashes;
      slot.flashed += p.count * std::min(f, 1.0) * (m + p.shed);
      if (f >= 1) { vaporise(); return; }
      m *= 1 - f;
      p.shed *= 1 - f;
      h = hSat;
      p.t = tSat;
      p.d = std::cbrt(6 * m / (pi * liquid_.at(p.t).density));
      if (m < minMass) { vaporise(); return; }
    }
    // Stage 1 and the sub-step's length.
    const auto r1 = rates(p, gas, work);
    if (r1.clamped) ++slot.clamped;
    if (newEddy) {
      const double slip = std::sqrt(sq(gas.uz - p.uz) + sq(gas.ur - p.ur) + sq(gas.ut - p.ut));
      const auto c = static_cast<std::size_t>(cell) * 2;
      p.eddyTime = interaction(cellTurbulence_[c], cellTurbulence_[c + 1], r1.tau, slip, schmidt);
    }
    const auto l1 = liquid_.at(p.t);
    double step = std::min(remaining, 0.25 * r1.tau);
    if (r1.mass > 0) step = std::min(step, 0.25 * m / r1.mass);
    if (!r1.boiling && r1.stiffness > 0) step = std::min(step, 0.25 * m * l1.cp / r1.stiffness);
    const double speed = std::sqrt(sq(p.uz) + sq(p.ur) + sq(p.ut));
    if (speed > 0) step = std::min(step, 0.5 * cellSize_[static_cast<std::size_t>(cell)] / speed);
    if (p.eddyTime > 0) step = std::min(step, p.eddyTime);  // the sub-step ends with the interaction
    if (step >= remaining) step = remaining;
    remaining = step == remaining ? 0 : remaining - step;
    // Mass and enthalpy: Heun on s = m^(2/3) and h.
    const double s0 = std::cbrt(m) * std::cbrt(m);
    const double fs1 = -2.0 / 3 * r1.mass / std::cbrt(m);
    // The vapour leaves at the surface temperature (saturation when boiling) with that latent heat.
    const double fh1 = (r1.heat - r1.mass * liquid_.at(r1.surface).latent) / m;
    double s1 = s0 + step * fs1, h1 = h + step * fh1;
    double m1 = m, t1 = p.t;
    if (s1 * std::sqrt(s1) < minMass || !(s1 > 0)) {
      m1 = 0;
    } else {
      Parcel q = p;
      const double hSat = liquid_.enthalpy(tSat);
      // Boiling: held at saturation while heat arrives (by the sign of the heat, not of fh1, which is
      // zero to within the rounding of a contracted multiply-add).
      const bool held = r1.boiling && r1.heat > 0;
      if (!r1.boiling && h1 > hSat) h1 = hSat;
      q.t = held ? p.t : h1 >= hSat ? tSat : liquid_.temperature(h1);
      const double mq = s1 * std::sqrt(s1);
      q.d = std::cbrt(6 * mq / (pi * liquid_.at(q.t).density));
      const auto r2 = rates(q, gas, work);
      const double fs2 = -2.0 / 3 * r2.mass / std::cbrt(mq);
      const double fh2 = (r2.heat - r2.mass * liquid_.at(r2.surface).latent) / mq;
      s1 = s0 + 0.5 * step * (fs1 + fs2);
      h1 = held ? h : h + 0.5 * step * (fh1 + fh2);
      m1 = r1.mass == 0 && r2.mass == 0 ? m : s1 * std::sqrt(s1);
      if (held) {
        t1 = p.t;
      } else if (h1 >= hSat) {  // reached saturation within the step: boiling from the next
        h1 = hSat;
        t1 = tSat;
      } else {
        t1 = liquid_.temperature(h1);
      }
    }
    // Velocity: exact relaxation toward the seen gas over the step at stage 1's drag; position exact.
    // Frame: x radial in the parcel's start plane, y azimuthal, z axial; the seen gas has (u_r, u_t, u_z),
    // u_t the eddy's alone.
    const double decay = std::exp(-step / r1.tau), lag = -std::expm1(-step / r1.tau) * r1.tau;
    const double vx = gas.ur + (p.ur - gas.ur) * decay, vy = gas.ut + (p.ut - gas.ut) * decay;
    const double vz = gas.uz + (p.uz - gas.uz) * decay;
    const double x = p.r + gas.ur * step + (p.ur - gas.ur) * lag, y = gas.ut * step + (p.ut - gas.ut) * lag;
    double z = p.z + gas.uz * step + (p.uz - gas.uz) * lag;
    p.distance += std::sqrt(sq(x - p.r) + sq(y) + sq(z - p.z));
    if (m1 < minMass) {  // vaporised within the sub-step
      vaporise();
      return;
    }
    // The exchange: minus the parcel's change over the sub-step.
    const double ke1 = 0.5 * (vx * vx + vy * vy + vz * vz);
    const double mass1 = p.count * (m1 + p.shed);
    const double e1 = mass1 * (h1 + ke1);
    deposit(cell, mass0 - mass1, pz0 - mass1 * vz, px0 - mass1 * vx, e0 - e1);
    slot.swirl += std::abs(mass1 * vy - py0);
    m = m1;
    h = h1;
    p.t = t1;
    p.d = std::cbrt(6 * m / (pi * liquid_.at(p.t).density));
    // Back to the meridional plane through the new position.
    double rr = std::hypot(x, y), ur = vx, ut = vy, uz = vz;
    if (p.eddyTime > 0) p.eddyTime = step >= p.eddyTime ? 0 : p.eddyTime - step;
    if (p.eddyTime <= 0) p.eddyZ = p.eddyR = p.eddyT = 0;
    if (rr > 0) {
      const double c = x / rr, sn = y / rr;
      ur = c * vx + sn * vy;
      ut = -sn * vx + c * vy;
      const double er = p.eddyR;  // the eddy's velocity is fixed in space: it turns with the frame
      p.eddyR = c * er + sn * p.eddyT;
      p.eddyT = -sn * er + c * p.eddyT;
    }
    // Walls and the outlet. Injector plate:
    const double mass = p.count * (m + p.shed);
    if (z < 0) {
      z = -e * z;
      slot.wallImpulse += -(1 + e) * mass * uz;
      slot.wallLoss += 0.5 * mass * (1 - e * e) * uz * uz;
      uz = -e * uz;
      ++slot.wallHits;
    }
    // Side wall: mirror about the column's wall line.
    for (int bounce = 0; bounce < 4; ++bounce) {
      double slope = 0;
      const double wall = wallRadius(z, &slope);
      if (rr <= wall) break;
      const double norm = std::sqrt(1 + slope * slope), nz = -slope / norm, nr = 1 / norm;
      const double depth = (rr - wall) / norm;
      z -= (1 + e) * depth * nz;
      rr -= (1 + e) * depth * nr;
      const double vn = uz * nz + ur * nr;
      if (vn > 0) {
        slot.wallImpulse += -(1 + e) * vn * nz * mass;
        slot.wallLoss += 0.5 * mass * (1 - e * e) * vn * vn;
        uz -= (1 + e) * vn * nz;
        ur -= (1 + e) * vn * nr;
      }
      ++slot.wallHits;
      if (rr < 0) { rr = -rr; ur = -ur; }
    }
    if (z >= span) {
      slot.outMass += mass;
      slot.outEnergy += mass * (h + 0.5 * (uz * uz + ur * ur + ut * ut));
      slot.outMomentum += mass * uz;
      slot.alive = false;
      return;
    }
    p.z = z;
    p.r = rr;
    p.uz = uz;
    p.ur = ur;
    p.ut = ut;
    p.cell = locate(p.z, p.r);
    if (p.cell < 0) {  // on a boundary to rounding (a corner): pulled just inside
      p.z = std::clamp(p.z, 0.0, std::nextafter(span, 0.0));
      p.r = std::min(p.r, wallRadius(p.z, nullptr) * (1 - 1e-12));
      p.cell = locate(p.z, p.r);
      if (p.cell < 0) throw std::runtime_error("Spray: a parcel left the domain through a wall.");
    }
    if (amr_) {
      const int from = levelOf(static_cast<std::size_t>(cell)), to = levelOf(static_cast<std::size_t>(p.cell));
      if (to > from) ++slot.finer;
      else if (to < from) ++slot.coarser;
    }
    // Breakup, at the sub-step's end state inside the domain; it keeps every total of the liquid.
    if (breakup.on) {
      const double u = std::sqrt(sq(gas.uz - p.uz) + sq(gas.ur - p.ur) + sq(gas.ut - p.ut));
      if (u > 0) {
        const auto l = liquid_.at(p.t);
        const double a = p.d / 2;
        const auto w = waves(p.d, u, gas.rho, l.density, l.viscosity, l.tension, u / r1.tau, breakup.cRt);
        const double stable = breakup.b0 * w.khLength;
        const double tauKh = 3.726 * breakup.b1 * a / (w.khLength * w.khRate);
        const double tauRt = breakup.cTau / w.rtRate;
        const double length = p.injector >= 0 ? breakup.cBl * std::sqrt(l.density / gas.rho) * blobDiameter(p.injector) : 0;
        if (p.distance >= length && w.rtLength < p.d && tauRt < tauKh) {
          p.rtTime += step;
          if (p.rtTime >= tauRt) {
            const double ratio = std::pow(p.d / w.rtLength, 3);
            p.count *= ratio;
            p.shed /= ratio;
            p.d = w.rtLength;
            m = p.dropMass(l.density);
            p.rtTime = 0;
          }
        } else {
          p.rtTime = 0;
          if (stable < a) {
            const double aNew = stable + (a - stable) * std::exp(-step / tauKh);
            const double mNew = l.density * pi / 6 * 8 * aNew * aNew * aNew;
            p.shed += m - mNew;
            m = mNew;
            p.d = 2 * aNew;
            if (p.count * p.shed >= breakup.shedFraction * p.injectedMass) {
              Parcel child = p;
              child.d = 2 * stable;
              child.count = p.count * p.shed / child.dropMass(l.density);
              child.shed = 0;
              child.random = next(p.random);
              slot.children.push_back(child);
              p.shed = 0;
            }
          }
        }
      }
    }
  }
}

// The opening's time integral from the valve's opening to t: the delivered mass is massFlow times it.
namespace {
double openedTime(const Spray::Injector& j, double t) {
  if (t <= j.opens) return 0;
  if (!(j.ramp > 0)) return t - j.opens;
  const double x = t - j.opens;
  return x < j.ramp ? x * x / (2 * j.ramp) : j.ramp / 2 + (x - j.ramp);
}
}  // namespace

void Spray::setParcelRate(double rate) {
  if (!(rate > 0)) throw std::invalid_argument("Spray: the parcel rate must be positive.");
  for (auto& j : settings_.injectors) j.parcelRate = rate;
}

// The injectors deliver [t0, t1]'s liquid into their pending stores, and release a parcel every
// 1 / parcelRate of open time, each carrying an equal share of what is pending.
void Spray::release(double t0, double t1) {
  ledger_.injectionThrust = 0;
  for (std::size_t k = 0; k < settings_.injectors.size(); ++k) {
    const auto& j = settings_.injectors[k];
    auto& pend = pending_[k];
    const auto l = liquid_.at(j.temperature);
    const double area = pi * (sq(j.rOuter) - sq(j.rInner));
    const double speed = j.velocity > 0 ? j.velocity : j.massFlow / (l.density * area);
    const double specific = l.enthalpy + 0.5 * speed * speed;
    const double delivered = j.massFlow * (openedTime(j, t1) - openedTime(j, t0));
    const double cosMax = std::cos(j.coneHalfAngle);
    // The injector's reaction, at the cone's mean axial speed (the parcels' own momentum is booked exactly).
    ledger_.injectionThrust += j.massFlow * j.opening(t1) * speed * 0.5 * (1 + cosMax);
    if (!(delivered > 0)) continue;
    pend.mass += delivered;
    pend.energy += delivered * specific;
    ledger_.injectedMass += delivered;
    ledger_.injectedEnergy += delivered * specific;
    pend.time += t1 - t0;
    const double interval = 1 / j.parcelRate;
    const auto n = static_cast<long>(std::floor(pend.time / interval));
    if (n <= 0) continue;
    pend.time -= static_cast<double>(n) * interval;
    const double share = pend.mass / static_cast<double>(n);
    auto& random = injectorRandom_[k];
    for (long c = 0; c < n; ++c) {
      Parcel p;
      p.injector = static_cast<int>(k);
      p.t = j.temperature;
      p.random = next(random);
      if (j.size == Injector::Size::Blob) {
        p.d = blobDiameter(static_cast<int>(k));
      } else {
        const double x = j.sauterDiameter * std::tgamma(1 - 1 / j.exponent);
        const double top = -std::expm1(-std::pow(3 * j.sauterDiameter / x, j.exponent));
        p.d = x * std::pow(-std::log1p(-uniform(random) * top), 1 / j.exponent);
      }
      if (j.size == Injector::Size::RosinRammler && j.coreLength > 0) {
        p.z = j.z + uniform(random) * j.coreLength;
        p.r = j.rOuter;
      } else {
        p.z = j.z;
        p.r = std::sqrt(sq(j.rInner) + uniform(random) * (sq(j.rOuter) - sq(j.rInner)));
      }
      const double cosTheta = 1 - uniform(random) * (1 - cosMax), sinTheta = std::sqrt(1 - cosTheta * cosTheta);
      const double phi = 2 * pi * uniform(random);
      p.uz = speed * cosTheta;
      p.ur = speed * sinTheta * std::cos(phi);
      p.ut = speed * sinTheta * std::sin(phi);
      p.count = share / p.dropMass(l.density);
      p.cell = locate(p.z, p.r);
      if (p.cell < 0) throw std::runtime_error("Spray: an injector releases outside the domain.");
      const double mass = p.count * p.dropMass(l.density);
      pend.mass -= mass;
      pend.energy -= mass * specific;
      ledger_.injectedMomentum += mass * p.uz;
      ++ledger_.parcelsInjected;
      p.injectedMass = mass;
      parcels_.push_back(p);
    }
  }
}

namespace {
struct Timer {
  double& total;
  std::chrono::steady_clock::time_point start{std::chrono::steady_clock::now()};
  ~Timer() { total += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); }
};
}  // namespace

void Spray::move(double dt, int phase) {
  const auto ns = flow_.medium().size();
  if (work_.size() != static_cast<std::size_t>(flow_.threads())) {
    work_.resize(static_cast<std::size_t>(flow_.threads()));
    for (auto& w : work_) w.resize(4 * ns + 8);
  }
  gatherGas();
  slots_.resize(parcels_.size());
  if (settings_.implicitCoupling && settings_.twoWay) applyLoads(dt);
  flow_.pool().blocks(moved_.size(), 16, [&](std::size_t begin, std::size_t end, int worker) {
    for (std::size_t n = begin; n < end; ++n) moveParcel(moved_[n], dt, worker);
  });
  for (auto& a : amounts_) std::fill(a.begin(), a.end(), 0.0);
  std::vector<char> deposited(levels_.size(), 0);
  std::vector<Parcel> kept;
  kept.reserve(parcels_.size());
  std::size_t next = 0;  // moved_ is in list order
  for (std::size_t i = 0; i < parcels_.size(); ++i) {
    if (next == moved_.size() || moved_[next] != i) {
      kept.push_back(parcels_[i]);
      continue;
    }
    ++next;
    const auto& s = slots_[i];
    for (const auto& d : s.deposits) {
      const auto cell = static_cast<std::size_t>(d.cell);
      const auto l = static_cast<std::size_t>(levelOf(cell));
      double* a = amounts_[l].data() + (cell - offset_[l]) * (ns + 3);
      deposited[l] = 1;
      a[oxygen_] += d.mass;
      a[ns] += d.pz;
      a[ns + 1] += d.pr;
      a[ns + 2] += d.energy;
      ledger_.exchangedMass += d.mass;
      ledger_.exchangedMomentum += d.pz;
      ledger_.exchangedEnergy += d.energy;
    }
    ledger_.wallImpulse += s.wallImpulse;
    ledger_.wallLoss += s.wallLoss;
    ledger_.droppedSwirl += s.swirl;
    ledger_.outflowMass += s.outMass;
    ledger_.outflowEnergy += s.outEnergy;
    ledger_.outflowMomentum += s.outMomentum;
    ledger_.substeps += s.substeps;
    ledger_.wallHits += s.wallHits;
    ledger_.flashes += s.flashes;
    ledger_.flashedMass += s.flashed;
    ledger_.vaporised += s.vaporised;
    ledger_.clampedCondensation += s.clamped;
    ledger_.finerMoves += s.finer;
    ledger_.coarserMoves += s.coarser;
    if (s.alive) {
      kept.push_back(parcels_[i]);
      kept.back().phase = phase;
    }
    for (const auto& c : s.children) {
      kept.push_back(c);
      kept.back().phase = phase;
      ++ledger_.children;
    }
  }
  parcels_.swap(kept);
  if (settings_.twoWay)
    for (std::size_t l = 0; l < levels_.size(); ++l)
      if (!amr_ || deposited[l]) levels_[l]->exchange(amounts_[l]);
}

void Spray::advance(double dt) {
  if (amr_) throw std::logic_error("Spray on AMR: step it with afterStep, not advance.");
  if (!(dt > 0)) return;
  Timer timer{wall_};
  const double t1 = flow_.time(), t0 = t1 - dt;
  const double impulse = ledger_.wallImpulse;
  if (!parcels_.empty()) {
    moved_.resize(parcels_.size());
    for (std::size_t i = 0; i < moved_.size(); ++i) moved_[i] = i;
    move(dt, 0);
  }
  ledger_.wallThrust = (ledger_.wallImpulse - impulse) / dt;
  release(t0, t1);
}

// SPRAY.md, Spray on AMR, The clock: the k-th step of level m within the base step ends at phase
// k 2^(top - m); it moves the parcels owned at or below m that are short of that phase, which are all at
// (k - 1) 2^(top - m), so all by the level's dt.
void Spray::afterStep(int level, double dt) {
  if (!amr_) throw std::logic_error("Spray: afterStep is for a spray on AMR.");
  if (!(dt > 0)) return;
  Timer timer{wall_};
  const int top = amr_->top();
  if (level == 0) {
    snapshot_ = {parcels_, ledger_, pending_, injectorRandom_, nextId_};
    snapshotStep_ = amr_->steps();
    calls_.assign(static_cast<std::size_t>(top) + 1, 0);
    baseDt_ = dt;
    baseImpulse_ = ledger_.wallImpulse;
    // The pattern can have changed since the last base step (a regrid); within one it is fixed and each
    // move finds its parcels' cells again.
    for (auto& p : parcels_) {
      p.phase = 0;
      p.cell = locate(p.z, p.r);
      if (p.cell < 0) throw std::runtime_error("Spray on AMR: a parcel is outside the domain.");
    }
  }
  if (level > top || calls_.size() != static_cast<std::size_t>(top) + 1)
    throw std::logic_error("Spray on AMR: a level stepped outside its base step.");
  const int width = 1 << (top - level);
  const int end = ++calls_[static_cast<std::size_t>(level)] * width;
  moved_.clear();
  for (std::size_t i = 0; i < parcels_.size(); ++i) {
    const auto& p = parcels_[i];
    if (levelOf(static_cast<std::size_t>(p.cell)) > level || p.phase >= end) continue;
    if (p.phase != end - width) throw std::logic_error("Spray on AMR: a parcel to move is not at the step's start phase.");
    moved_.push_back(i);
  }
  if (!moved_.empty()) move(dt, end);
  if (level == 0) {
    const double t1 = flow_.time();
    const auto before = parcels_.size();
    release(t1 - dt, t1);
    for (std::size_t i = before; i < parcels_.size(); ++i) parcels_[i].phase = 1 << top;
  }
  if (level == top && end == 1 << top) ledger_.wallThrust = (ledger_.wallImpulse - baseImpulse_) / baseDt_;
}

void Spray::redo() {
  if (!amr_ || snapshotStep_ != amr_->steps()) return;  // the failed attempt did not reach the spray
  parcels_ = snapshot_.parcels;
  ledger_ = snapshot_.ledger;
  pending_ = snapshot_.pending;
  injectorRandom_ = snapshot_.injectorRandom;
  nextId_ = snapshot_.nextId;
}

double Spray::nextEvent(double time) const {
  double out = std::numeric_limits<double>::infinity();
  for (const auto& j : settings_.injectors) {
    if (j.opens > time) out = std::min(out, j.opens);
    if (j.ramp > 0 && j.opens + j.ramp > time) out = std::min(out, j.opens + j.ramp);
  }
  return out;
}

double Spray::exchangeDt(double fraction) const {
  if (parcels_.empty()) return std::numeric_limits<double>::infinity();
  const auto& medium = flow_.medium();
  const auto ns = medium.size();
  const std::size_t n = parcels_.size(), threads = static_cast<std::size_t>(flow_.threads());
  // The gas of level l's cell q as it stands (its own temperature guess; the spray's cache is not touched).
  auto gasOf = [&](std::size_t c, double* y) {
    const int l = levelOf(c);
    const Flow& level = *levels_[static_cast<std::size_t>(l)];
    const std::size_t q = c - offset_[static_cast<std::size_t>(l)];
    const auto& partial = level.partialDensities();
    double sum = 0;
    for (std::size_t k = 0; k < ns; ++k) sum += std::max(partial[q * ns + k], 0.0);
    for (std::size_t k = 0; k < ns; ++k) y[k] = std::max(partial[q * ns + k], 0.0) / sum;
    const auto& u = level.state()[q];
    const double t = level.temperature(q);
    return Gas{t, u[0] * medium.gasConstant(y) * t, u[0], u[1] / u[0], u[2] / u[0], y};
  };
  std::vector<int> cell(n);
  std::vector<double> conductance(n), vapour(n);
  std::vector<std::vector<double>> work(threads, std::vector<double>(4 * ns + 8)), fractions(threads, std::vector<double>(ns));
  flow_.pool().blocks(n, 64, [&](std::size_t begin, std::size_t end, int worker) {
    auto& w = work[static_cast<std::size_t>(worker)];
    double* y = fractions[static_cast<std::size_t>(worker)].data();
    for (std::size_t i = begin; i < end; ++i) {
      const Parcel& p = parcels_[i];
      cell[i] = locate(p.z, p.r);
      if (cell[i] < 0) continue;
      Gas gas = gasOf(static_cast<std::size_t>(cell[i]), y);
      gas.uz += p.eddyZ;
      gas.ur += p.eddyR;
      gas.ut = p.eddyT;
      const auto r = rates(p, gas, w);
      conductance[i] = p.count * r.conductance;
      vapour[i] = p.count * r.mass;
    }
  });
  // Summed per cell in parcel order, so the bound does not depend on the threads.
  std::vector<double> sumG(offset_.back(), 0.0), sumM(offset_.back(), 0.0);
  std::vector<std::size_t> occupied;
  for (std::size_t i = 0; i < n; ++i) {
    if (cell[i] < 0) continue;
    const auto c = static_cast<std::size_t>(cell[i]);
    if (sumG[c] == 0 && sumM[c] == 0) occupied.push_back(c);
    sumG[c] += conductance[i];
    sumM[c] += vapour[i];
  }
  double out = std::numeric_limits<double>::infinity();
  std::vector<double> y(ns);
  for (const auto c : occupied) {
    const int l = levelOf(c);
    const Gas gas = gasOf(c, y.data());
    const double m = gas.rho * levels_[static_cast<std::size_t>(l)]->mesh().cells[c - offset_[static_cast<std::size_t>(l)]].volume;
    double tau = std::numeric_limits<double>::infinity();
    if (sumG[c] > 0) tau = std::min(tau, m * medium.cv(gas.t, y.data()) / sumG[c]);
    if (sumM[c] > 0) tau = std::min(tau, m / sumM[c]);
    out = std::min(out, static_cast<double>(1 << l) * fraction * tau);
  }
  return out;
}

double Spray::liquidMass() const {
  double out = 0;
  for (const auto& p : parcels_) out += p.count * (p.dropMass(liquid_.at(p.t).density) + p.shed);
  for (const auto& q : pending_) out += q.mass;
  return out;
}
double Spray::liquidEnergy() const {
  double out = 0;
  for (const auto& p : parcels_) {
    const auto l = liquid_.at(p.t);
    out += p.count * (p.dropMass(l.density) + p.shed) * (l.enthalpy + 0.5 * (sq(p.uz) + sq(p.ur) + sq(p.ut)));
  }
  for (const auto& q : pending_) out += q.energy;
  return out;
}
double Spray::liquidMomentum() const {
  double out = 0;
  for (const auto& p : parcels_) out += p.count * (p.dropMass(liquid_.at(p.t).density) + p.shed) * p.uz;
  return out;
}
std::array<double, 3> Spray::balanceErrors() const {
  const auto& l = ledger_;
  const double tiny = std::numeric_limits<double>::min();
  const double mass = (liquidMass() - (l.injectedMass - l.exchangedMass - l.outflowMass)) / std::max(l.injectedMass, tiny);
  const double energy = (liquidEnergy() - (l.injectedEnergy - l.exchangedEnergy - l.outflowEnergy - l.wallLoss)) /
                        std::max(std::abs(l.injectedEnergy) + std::abs(l.exchangedEnergy) + std::abs(l.outflowEnergy) + l.wallLoss, tiny);
  const double momentum = (liquidMomentum() - (l.injectedMomentum - l.exchangedMomentum - l.outflowMomentum + l.wallImpulse)) /
                          std::max(std::abs(l.injectedMomentum) + std::abs(l.exchangedMomentum) + std::abs(l.outflowMomentum) + std::abs(l.wallImpulse), tiny);
  return {mass, energy, momentum};
}
// The checkpoint's fields (docs/evidence/RESTART.md). Parcels go field by field, so Parcel's padding never
// reaches the file and equal states give equal bytes.
namespace {
template <class P, class Io> void parcelFields(P& p, Io& io) {
  io.field(p.z), io.field(p.r), io.field(p.uz), io.field(p.ur), io.field(p.ut), io.field(p.d), io.field(p.t);
  io.field(p.count), io.field(p.shed), io.field(p.rtTime), io.field(p.distance), io.field(p.injectedMass);
  io.field(p.eddyZ), io.field(p.eddyR), io.field(p.eddyT), io.field(p.eddyTime), io.field(p.random);
  io.field(p.cell), io.field(p.injector), io.field(p.phase);
}
template <class V, class Io> void parcelList(V& list, Io& io) {
  std::uint64_t n = list.size();
  io.field(n);
  if constexpr (Io::reading) list.resize(static_cast<std::size_t>(n));
  for (auto& p : list) parcelFields(p, io);
}
}  // namespace
template <class Self, class Io> void Spray::checkpoint(Self& self, Io& io) {
  io.mark("spray");
  io.same(static_cast<std::uint64_t>(self.cellGas_.size()), "the spray's composite cells");
  io.same(static_cast<std::uint64_t>(self.settings_.injectors.size()), "the injector count");
  for (auto& j : self.settings_.injectors) io.field(j.parcelRate);
  parcelList(self.parcels_, io);
  io.field(self.ledger_), io.field(self.pending_), io.field(self.injectorRandom_), io.field(self.nextId_);
  io.field(self.calls_), io.field(self.baseDt_), io.field(self.baseImpulse_);
  parcelList(self.snapshot_.parcels, io);
  io.field(self.snapshot_.ledger), io.field(self.snapshot_.pending), io.field(self.snapshot_.injectorRandom);
  io.field(self.snapshot_.nextId), io.field(self.snapshotStep_);
  io.field(self.cellGas_);
}
void Spray::writeCheckpoint(CheckpointWriter& out) const { checkpoint(*this, out); }
void Spray::readCheckpoint(CheckpointReader& in) { checkpoint(*this, in); }
}  // namespace crucible
