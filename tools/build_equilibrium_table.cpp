// Builds Table A, the tabulated shifting equilibrium (core/equilibrium_table.hpp; TECHNICAL_PLAN,
// *Lightweight engine*; docs/evidence/TABLE_A.md), offline from Cantera's isothermal equilibrium
// (ReactionSource::equilibrateTV, the same solver as the march's equilibrateUV) of the reacting
// (H/O) species alone; the builder checks that every other species is the only carrier of its
// elements. Along each (f, rho_r) line the first node starts from the complete-combustion products and each later node
// from the previous node's equilibrium (equilibrium is unique, so the start only changes the cost).
// Cantera conserves the elements to its own tolerance (about 1e-10); each node is then projected onto
// its exact elements (smallest change weighted by Y_k), and its e is that of the projected Y. Prints
// the provenance and the conservation errors before and after the projection (TABLE_A criterion 1).
// The table file is not kept in git.
//
// The T axis is the mechanism's thermo range: outside it Cantera's equilibrium falls back to a
// solver costing about 70 ms per call against 30 to 60 us inside (measured, 5 October), and the
// NASA polynomials are extrapolated. States hotter than the range take the direct call in the run.
//
// Usage: crucible_build_equilibrium_table <out> <f nodes, lean side> <f nodes, rich side> <f spacing at f_st>
//                                         <T nodes> <rho_r nodes> <rho_r min> <rho_r max> <threads>
// The f axis is uniform on [0, f_st] and on [f_st, 1], with f_st that of H2O, plus nodes at
// f_st +- d, 0 + d and 1 - d, d growing from the given spacing by 1.5 per node while below the
// uniform spacing (0: none): near f_st the equilibrium has a near-kink (excess O2 on one side, H2 on
// the other) whose width is set by dissociation, and near the ends one element is a trace. T and ln rho_r are uniform.
//
//        crucible_build_equilibrium_table check <table> <samples> <threads> [rho min] [rho max] [seed]
// checks a table against direct Cantera (TABLE_A criteria 1 for lookups and 2a). States are uniform
// in (f, z = Y_N2, eta, ln rho), rho over [rho min, rho max] (default 1e-4 to 30 kg/m3), eta = (e - e_lo) / (e_hi - e_lo) between the query's complete-combustion
// products at the range's lowest T and the same products, frozen, at 5000 K: a superset of the
// table's energies (states above it are counted outside). |dT| is the core's temperature
// (Medium::temperature) from the table's Y at the query's e, against the march's direct call
// (equilibrateUV).
#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "adapters/reaction.hpp"
#include "cantera/base/global.h"
#include "cantera/core.h"
#include "core/equilibrium_table.hpp"

using namespace crucible;

namespace {
constexpr double kHighT = 5000;  // check: frozen products at the top of the sampled energies

struct Products {
  std::size_t h2, o2, h2o, n2;
  double hInWater;  // mass fraction of H in H2O, the stoichiometric f
};

// Complete-combustion products of the elements (Z_H, Z_O, Z_N): all of the scarcer element in H2O.
std::vector<double> products(const Products& s, std::size_t ns, double zh, double zo, double zn) {
  std::vector<double> y(ns, 0.0);
  if (zh * (1 - s.hInWater) >= zo * s.hInWater) {  // H-rich: all O in H2O
    y[s.h2o] = zo / (1 - s.hInWater);
    y[s.h2] = std::max(zh - y[s.h2o] * s.hInWater, 0.0);
  } else {
    y[s.h2o] = zh / s.hInWater;
    y[s.o2] = std::max(zo - y[s.h2o] * (1 - s.hInWater), 0.0);
  }
  y[s.n2] = zn;
  return y;
}

}  // namespace

// Indices of the species complete combustion uses.
Products productIndices(thermo::ReactionSource& source) {
  auto index = [&](const char* name) {
    for (std::size_t k = 0; k < source.nSpecies(); ++k)
      if (source.speciesName(k) == name) return k;
    throw std::runtime_error(std::string("species missing: ") + name);
  };
  return {index("H2"), index("O2"), index("H2O"), index("N2"), 0};
}

// Uniform on [0, 1) from (seed, sample, draw) by SplitMix64, so a sample does not depend on threads.
double uniform(std::uint64_t seed, std::uint64_t sample, std::uint64_t draw) {
  std::uint64_t x = seed * 0x9e3779b97f4a7c15ULL + sample * 0xbf58476d1ce4e5b9ULL + draw * 0x94d049bb133111ebULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x ^= x >> 31;
  return static_cast<double>(x >> 11) * 0x1.0p-53;
}

int check(const std::string& path, std::size_t samples, int threads, double rhoMin, double rhoMax, std::uint64_t seed) {
  const auto table = EquilibriumTable::load(path);
  std::printf("%s\n", table.provenance.c_str());
  const std::string mechanism = "h2o2.yaml";
  thermo::ReactionSource probe(mechanism);
  const std::size_t ns = probe.nSpecies();
  if (table.species.size() != ns) throw std::runtime_error("table species differ from the mechanism");
  for (std::size_t k = 0; k < ns; ++k)
    if (table.species[k] != probe.speciesName(k)) throw std::runtime_error("table species differ from the mechanism");
  Products s = productIndices(probe);
  s.hInWater = table.elements[s.h2o][0];
  double lowT;
  {
    auto solution = Cantera::newSolution(mechanism, "", "none");
    lowT = solution->thermo()->minTemp();
  }
  const Medium medium = probe.medium();

  // Per sample: |dT| (NaN when not compared), Cantera's T, the interpolated T's error, and the result.
  struct Sample {
    double f, z, eta, lnRho, tRef = std::nan(""), dT = std::nan(""), dTable = std::nan("");
    int result = -1;
  };
  std::vector<Sample> out(samples);
  std::vector<std::array<double, 16>> worstY(threads);
  std::atomic<std::size_t> next{0};
  std::mutex report;
  double worstElement = 0, worstSum = 0;
  std::size_t failures = 0;
  const auto clock0 = std::chrono::steady_clock::now();
  auto run = [&](int w) {
    thermo::ReactionSource source(mechanism);
    std::vector<double> z(ns + 1), y, yt(ns);
    double element = 0, sum = 0;
    std::size_t failed = 0;
    auto& dy = worstY[w];
    dy.fill(0);
    for (std::size_t i; (i = next.fetch_add(1)) < samples;) {
      Sample& c = out[i];
      c.f = uniform(seed, i, 0);
      c.z = uniform(seed, i, 1);
      c.eta = uniform(seed, i, 2);
      c.lnRho = std::log(rhoMin) + std::log(rhoMax / rhoMin) * uniform(seed, i, 3);
      const double target[3] = {c.f * (1 - c.z), (1 - c.f) * (1 - c.z), c.z};
      const double rho = std::exp(c.lnRho);
      y = products(s, ns, target[0], target[1], target[2]);
      const double eLow = source.internalEnergy(lowT, 1, y.data()), eHigh = source.internalEnergy(kHighT, 1, y.data());
      const double e = eLow + c.eta * (eHigh - eLow);
      double t = 0;
      const auto result = table.lookup(y.data(), e, rho, yt.data(), t);
      c.result = static_cast<int>(result);
      if (result == EquilibriumTable::Result::Outside) continue;  // the run takes the direct call
      // Criterion 1 for the lookup: every element of the mechanism.
      const auto want = source.elementMassFractions(y.data()), got = source.elementMassFractions(yt.data());
      double total = -1;
      for (std::size_t q = 0; q < ns; ++q) total += yt[q];
      for (std::size_t m = 0; m < want.size(); ++m) element = std::max(element, std::abs(got[m] - want[m]));
      sum = std::max(sum, std::abs(total));
      try {
        z[0] = medium.temperature(e, y.data(), 1000);  // the frozen start, as the march's
        std::copy(y.begin(), y.end(), z.begin() + 1);
        source.equilibrateUV(rho, z.data());
      } catch (const std::exception& ex) {
        ++failed;
        std::lock_guard<std::mutex> lock(report);
        std::fprintf(stderr, "reference f %.5f z %.4f eta %.5f rho %.3e failed: %s\n", c.f, c.z, c.eta, rho, ex.what());
        continue;
      }
      c.tRef = z[0];
      c.dT = std::abs(medium.temperature(e, yt.data(), t) - z[0]);
      c.dTable = std::abs(t - z[0]);
      for (std::size_t q = 0; q < ns && q < dy.size(); ++q) dy[q] = std::max(dy[q], std::abs(yt[q] - z[q + 1]));
    }
    std::lock_guard<std::mutex> lock(report);
    worstElement = std::max(worstElement, element);
    worstSum = std::max(worstSum, sum);
    failures += failed;
  };
  std::vector<std::thread> pool;
  for (int w = 0; w < threads; ++w) pool.emplace_back(run, w);
  for (auto& t : pool) t.join();
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();

  std::size_t outside = 0, clamped = 0;
  std::vector<double> band, all;
  const Sample* worst = nullptr;
  for (const auto& c : out) {
    if (c.result == static_cast<int>(EquilibriumTable::Result::Outside)) ++outside;
    if (c.result == static_cast<int>(EquilibriumTable::Result::Clamped)) ++clamped;
    if (std::isnan(c.dT)) continue;
    all.push_back(c.dT);
    if (c.tRef >= 200 && c.tRef <= 4000) {
      band.push_back(c.dT);
      if (!worst || c.dT > worst->dT) worst = &c;
    }
  }
  auto quantile = [](std::vector<double> v, double p) {
    if (v.empty()) return std::nan("");
    const auto k = static_cast<std::size_t>(p * static_cast<double>(v.size() - 1));
    std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(k), v.end());
    return v[k];
  };
  std::printf("%zu samples (seed %llu, rho %.3g to %.3g kg/m3) in %.0f s on %d threads: %zu outside the table (direct "
              "call in the run), %zu clamped, %zu Cantera failures\n", samples, static_cast<unsigned long long>(seed),
              rhoMin, rhoMax, wall, threads, outside, clamped, failures);
  std::printf("criterion 1 (lookups): element error %.2e, |sum Y - 1| %.2e (limit 1e-12): %s\n", worstElement, worstSum,
              worstElement <= 1e-12 && worstSum <= 1e-12 ? "pass" : "FAIL");
  // Bands of Cantera's T, for the resolution choice.
  const double edges[] = {200, 1000, 2000, 2500, 3000, 3500, 4000};
  for (int b = 0; b + 1 < 7; ++b) {
    std::vector<double> v;
    for (const auto& c : out)
      if (!std::isnan(c.dT) && c.tRef >= edges[b] && c.tRef < edges[b + 1]) v.push_back(c.dT);
    std::printf("  T_eq %4.0f-%4.0f K: %6zu states, |dT| median %.3f, 99th %.3f, max %.3f K\n", edges[b], edges[b + 1],
                v.size(), quantile(v, 0.5), quantile(v, 0.99), v.empty() ? std::nan("") : *std::max_element(v.begin(), v.end()));
  }
  std::vector<double> tableT;
  for (const auto& c : out)
    if (!std::isnan(c.dTable) && c.tRef >= 200 && c.tRef <= 4000) tableT.push_back(c.dTable);
  std::printf("interpolated table T against Cantera's (consistency only): 99th %.3f K\n", quantile(tableT, 0.99));
  std::printf("max |dY_k|:");
  for (std::size_t q = 0; q < ns && q < 16; ++q) {
    double m = 0;
    for (const auto& dy : worstY) m = std::max(m, dy[q]);
    std::printf(" %s %.2e", table.species[q].c_str(), m);
  }
  std::printf("\n");
  if (worst)
    std::printf("worst state in 200-4000 K: f %.5f z %.4f eta %.5f rho %.4e T_eq %.1f K |dT| %.3f K\n", worst->f,
                worst->z, worst->eta, std::exp(worst->lnRho), worst->tRef, worst->dT);
  const double p99 = quantile(band, 0.99);
  std::printf("criterion 2(a): T_eq 200-4000 K, %zu states: |dT| median %.3f K, 99th percentile %.3f K, max %.3f K "
              "(limit 3 K on the 99th): %s\n", band.size(), quantile(band, 0.5), p99,
              band.empty() ? std::nan("") : *std::max_element(band.begin(), band.end()), p99 <= 3 ? "pass" : "FAIL");
  return worstElement <= 1e-12 && worstSum <= 1e-12 && p99 <= 3 && failures == 0 ? 0 : 1;
}

int main(int argc, char** argv) {
  if (argc >= 5 && std::string(argv[1]) == "check")
    return check(argv[2], std::strtoull(argv[3], nullptr, 10), std::max(1, std::atoi(argv[4])),
                 argc > 5 ? std::atof(argv[5]) : 1e-4, argc > 6 ? std::atof(argv[6]) : 30,
                 argc > 7 ? std::strtoull(argv[7], nullptr, 10) : 1);
  if (argc != 10) {
    std::fprintf(stderr, "usage: %s <out> <f lean> <f rich> <f spacing at f_st> <T> <rho_r> <rho_r min> <rho_r max> <threads>\n", argv[0]);
    return 2;
  }
  const std::string out = argv[1];
  const int fLean = std::atoi(argv[2]), fRich = std::atoi(argv[3]);
  const double fCluster = std::atof(argv[4]);
  const int nT = std::atoi(argv[5]), nRho = std::atoi(argv[6]), threads = std::max(1, std::atoi(argv[9]));
  const double rhoMin = std::atof(argv[7]), rhoMax = std::atof(argv[8]);
  const std::string mechanism = "h2o2.yaml";

  thermo::ReactionSource probe(mechanism);
  double minT, maxT;
  {
    auto solution = Cantera::newSolution(mechanism, "", "none");
    minT = solution->thermo()->minTemp();
    maxT = solution->thermo()->maxTemp();
  }
  const std::size_t ns = probe.nSpecies(), ne = probe.nElements();
  EquilibriumTable table;
  Products s = productIndices(probe);
  // Element columns from pure species: element m's mass fraction in species k.
  std::vector<std::vector<double>> columns;
  std::vector<double> unit(ns, 0.0);
  for (std::size_t k = 0; k < ns; ++k) {
    unit[k] = 1;
    columns.push_back(probe.elementMassFractions(unit.data()));
    unit[k] = 0;
  }
  const auto mH = std::max_element(columns[s.h2].begin(), columns[s.h2].end()) - columns[s.h2].begin();
  const auto mO = std::max_element(columns[s.o2].begin(), columns[s.o2].end()) - columns[s.o2].begin();
  for (std::size_t k = 0; k < ns; ++k) {
    table.species.push_back(probe.speciesName(k));
    const double h = columns[k][mH], o = columns[k][mO];
    if (h + o > 1 - 1e-12) {
      table.elements.push_back({h, o});
      continue;
    }
    if (h + o > 1e-12) throw std::runtime_error("species " + table.species[k] + " mixes H/O with another element");
    // Inert: the only carrier of each of its elements.
    for (std::size_t m = 0; m < ne; ++m)
      for (std::size_t r = 0; r < ns; ++r)
        if (r != k && columns[k][m] > 0 && columns[r][m] > 0)
          throw std::runtime_error("species " + table.species[k] + " shares an element with " + probe.speciesName(r));
    table.elements.push_back({0, 0});
    table.inert.push_back(k);
  }
  s.hInWater = table.elements[s.h2o][0];

  for (int i = 0; i <= fLean; ++i) table.f.push_back(s.hInWater * i / fLean);
  for (int i = 1; i <= fRich; ++i) table.f.push_back(s.hInWater + (1 - s.hInWater) * i / fRich);
  std::size_t clustered = 0;
  if (fCluster > 0) {
    for (double d = fCluster; d < s.hInWater / fLean; d *= 1.5, ++clustered) table.f.push_back(s.hInWater - d);
    for (double d = fCluster; d < (1 - s.hInWater) / fRich; d *= 1.5, ++clustered) table.f.push_back(s.hInWater + d);
    // Likewise at the ends, where one element is a trace.
    for (double d = fCluster; d < s.hInWater / fLean; d *= 1.5, ++clustered) table.f.push_back(d);
    for (double d = fCluster; d < (1 - s.hInWater) / fRich; d *= 1.5, ++clustered) table.f.push_back(1 - d);
    std::sort(table.f.begin(), table.f.end());
  }
  table.nT = nT;
  table.tMin = minT;
  table.tMax = maxT;
  table.nRho = nRho;
  table.lnRhoMin = std::log(rhoMin);
  table.lnRhoMax = std::log(rhoMax);
  const std::size_t nf = table.f.size(), nv = table.width();
  table.values.resize(nf * nT * nRho * nv);
  for (std::size_t q : table.inert)
    for (int k = 0; k < nT; ++k) {
      unit[q] = 1;
      table.inertEnergy.push_back(probe.internalEnergy(table.temperature(k), 1, unit.data()));
      unit[q] = 0;
    }

  // Work items: one (f, rho_r) line of T nodes each.
  std::atomic<std::size_t> next{0};
  const std::size_t lines = nf * nRho;
  std::mutex report;
  double worstElement = 0, worstSum = 0, worstRaw = 0;
  std::size_t failures = 0;
  const auto clock0 = std::chrono::steady_clock::now();
  auto run = [&]() {
    thermo::ReactionSource source(mechanism);
    std::vector<double> zstate(ns + 1), y(ns);
    double element = 0, sum = 0, raw = 0;
    std::size_t failed = 0;
    // Element errors of y against the node's (H, O), and sum Y - 1.
    auto errors = [&](const double* y, const double target[2], double err[2]) {
      double total = 0;
      for (int m = 0; m < 2; ++m) err[m] = -target[m];
      for (std::size_t q = 0; q < ns; ++q) {
        total += y[q];
        for (int m = 0; m < 2; ++m) err[m] += y[q] * table.elements[q][m];
      }
      return total - 1;
    };
    for (std::size_t line; (line = next.fetch_add(1)) < lines;) {
      const std::size_t i = line / nRho, l = line % nRho;
      const double f = table.f[i];
      const double target[2] = {f, 1 - f};
      const double rho = std::exp(table.lnRho(l));
      y = products(s, ns, target[0], target[1], 0);
      for (int k = 0; k < nT; ++k) {
        const double t = table.temperature(k);
        double* v = table.values.data() + table.node(i, k, l) * nv;
        try {
          zstate[0] = t;
          std::copy(y.begin(), y.end(), zstate.begin() + 1);
          source.equilibrateTV(rho, zstate.data());
        } catch (const std::exception& ex) {
          ++failed;
          std::lock_guard<std::mutex> lock(report);
          std::fprintf(stderr, "node f %.4f T %.1f rho_r %.3e failed: %s\n", f, t, rho, ex.what());
          y = products(s, ns, target[0], target[1], 0);
          continue;
        }
        for (std::size_t q = 0; q < ns; ++q) y[q] = std::max(zstate[q + 1], 0.0);
        // Projection: y_q += y_q sum_m a_qm lambda_m with (sum_q y_q a_qm a_qn) lambda = -err.
        double err[2];
        errors(y.data(), target, err);
        raw = std::max({raw, std::abs(err[0]), std::abs(err[1])});
        for (int pass = 0; pass < 2; ++pass) {
          double m[2][2] = {}, lambda[2] = {-err[0], -err[1]};
          for (std::size_t q = 0; q < ns; ++q)
            for (int a = 0; a < 2; ++a)
              for (int b = 0; b < 2; ++b) m[a][b] += y[q] * table.elements[q][a] * table.elements[q][b];
          // An element absent from the node (f = 0 or 1) has a zero row: leave it out.
          for (int a = 0; a < 2; ++a)
            if (m[a][a] == 0) { m[a][a] = 1; m[a][1 - a] = m[1 - a][a] = 0; lambda[a] = 0; }
          const double det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
          const double l0 = (lambda[0] * m[1][1] - m[0][1] * lambda[1]) / det;
          const double l1 = (m[0][0] * lambda[1] - m[1][0] * lambda[0]) / det;
          for (std::size_t q = 0; q < ns; ++q) y[q] += y[q] * (table.elements[q][0] * l0 + table.elements[q][1] * l1);
          errors(y.data(), target, err);
        }
        sum = std::max(sum, std::abs(errors(y.data(), target, err)));
        element = std::max({element, std::abs(err[0]), std::abs(err[1])});
        v[0] = source.internalEnergy(t, rho, y.data());
        std::copy(y.begin(), y.end(), v + 1);
      }
    }
    std::lock_guard<std::mutex> lock(report);
    worstElement = std::max(worstElement, element);
    worstSum = std::max(worstSum, sum);
    worstRaw = std::max(worstRaw, raw);
    failures += failed;
  };
  std::vector<std::thread> pool;
  for (int w = 0; w < threads; ++w) pool.emplace_back(run);
  for (auto& t : pool) t.join();
  const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - clock0).count();

  char text[1024];
  const std::time_t now = std::time(nullptr);
  char date[64];
  std::strftime(date, sizeof date, "%Y-%m-%d %H:%M %Z", std::localtime(&now));
  std::snprintf(text, sizeof text,
                "Table A built %s by crucible_build_equilibrium_table from %s, Cantera %s; f %zu nodes (%d lean, %d "
                "rich uniform, node at f_st %.6f, %zu clustered at it and the ends from spacing %.2g), T %d over [%.0f, %.0f] K, ln rho_r %d over "
                "[%.3g, %.3g] kg/m3; %zu inert species; %zu nodes in %.0f s on %d threads; node element error %.2e from Cantera, %.2e after the "
                "projection, |sum Y - 1| %.2e, failures %zu",
                date, mechanism.c_str(), Cantera::version().c_str(), table.f.size(), fLean, fRich, s.hInWater, clustered, fCluster,
                nT, minT, maxT, nRho, rhoMin, rhoMax, table.inert.size(), nf * nT * nRho, wall, threads,
                worstRaw, worstElement, worstSum, failures);
  table.provenance = text;
  std::printf("%s\n", text);
  std::printf("criterion 1 (nodes): element error %.2e, |sum Y - 1| %.2e (limit 1e-12): %s\n", worstElement, worstSum,
              worstElement <= 1e-12 && worstSum <= 1e-12 && failures == 0 ? "pass" : "FAIL");
  table.save(out);
  std::printf("wrote %s (%.1f MB)\n", out.c_str(), table.values.size() * 8.0 / 1e6);
  return failures == 0 ? 0 : 1;
}
