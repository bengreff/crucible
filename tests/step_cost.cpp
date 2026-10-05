// Cost of the reacting multi-species step per cell, for the compute budget in TECHNICAL_PLAN
// (section "Compute budget"). Not a verification test: it measures cell updates per second on the
// machine that runs it.
//
// The C1 chamber of tests/chamber_study.cpp (same contour and supply: premixed H2/O2 at O/F 5,
// 0.4 kg/s, h2o2.yaml, ambient N2 at 1000 Pa) is first marched in LocalEquilibrium mode to the setup
// time, so the chamber holds burned gas and the nozzle flows. From that one state each
// configuration takes a few untimed steps and then the timed steps:
//   inviscid: slip walls, no transport (the C1 configuration);
//   viscous:  mixture-averaged transport, SST-2003, no-slip walls at 600 K, supply turbulence
//             I 0.05 and viscosity ratio 10, Spalart-Rumsey ambient; k = 3/2 (0.05 |u|)^2 and
//             omega = rho k / (10 mu) in every cell at the start;
// each with Frozen chemistry (the flow step alone) and FiniteRate (CVODES, rtol 1e-6, atol 1e-12,
// as C1), on 1 thread and on the given thread count. The flow step is serial; only the chemistry
// is threaded. The valve is fully open from the start of each timed configuration.
// Last, the FiniteRate reaction substep alone on 1 thread (inviscid, after 20 FiniteRate steps from
// the setup state, so the cells are off equilibrium as in a march), 20 calls at each of three
// substep lengths, because the CVODES cost of a substep may depend on its length: half the
// measured flow step, and the halves of the explicit steps estimated for a wall-resolved RL10
// grid (TECHNICAL_PLAN, "Compute budget").
// With the optional last argument "equilibrium" (added 5 October 2026), the configurations are
// Frozen and LocalEquilibrium (a constant-(u, v) equilibrium per cell after each flow step, no
// CVODES) instead of Frozen and FiniteRate, and the substep sweep is skipped.
// With the optional last argument "cache" (added 5 October 2026, before its first run), nothing is
// timed as above. Instead it measures what a chemistry cache could reuse, on the reaction calls of
// the viscous FiniteRate march (5 untimed steps, then <timed steps> steps on the given threads, each
// step as ReactingFlow::step, replans included, with every cell of every reaction call recorded:
// the input T, Y, rho and substep length, and the output T, Y after the call). A query's key is its
// bin in x = (T, Y_1..Y_K, ln rho, ln dt), with bin widths (1000 K, 1, 1) times w, at
// w = 1e-2, 1e-3, 1e-4, 1e-5, and ln dt binned at a fixed 1e-2 (a march can take its steps from a
// fixed ladder of lengths at no cost to the physics, so the step length should not split bins).
// Calls are taken in march order. A query whose bin already holds a record is a hit, and the record
// supplies the increment dz = z_out - z_in by two retrievals:
//   constant: the record's dz times dt / dt_record;
//   linear:   the record's dz + J (x - x_record), with J = d(dz)/dx at the record by forward
//             differences of a CVODES call at rtol 1e-10, atol 1e-14 (ISAT's retrieval, Pope 1997).
// Any other query is a miss, and the first query in a bin becomes its record. Two tables are
// reported: "any" (the record may come from any cell) and "other" (a hit needs a record from a
// different cell, so a cell's own earlier calls do not count; this removes the reuse that comes from
// a cell changing little between steps).
// The retrieval error against the call's own result is measured in two forms:
//   the integrator's weighted norm, max_k |error_k| / (1e-6 |z_k| + 1e-12), so 1 is CVODES's own
//   local tolerance in C1;
//   plainly, |error in T| in K and the largest |error in Y_k|.
// Reported per w and table: queries, hits, records, the 50%, 99% and largest weighted error of each
// retrieval, the fraction of all queries that are hits within weighted error 1, 10 and 100, the
// largest plain errors, and the cost of one retrieval (key, lookup and the linear sum; measured).
// The integration error itself is reported on the same scale: every 97th call is repeated at rtol
// 1e-10, atol 1e-14 and compared with the march's result.
//
// Usage: crucible_step_cost <nz> <nr> <setup time s> <timed steps> <threads> [equilibrium | cache]
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "adapters/reacting_flow.hpp"
#include "adapters/reaction.hpp"

using namespace crucible;

// Chamber geometry and supply: copied from tests/chamber_study.cpp (keep them equal).
namespace {
constexpr double kPi = std::numbers::pi;
constexpr double kRt = 0.010, kRc = 0.025, kArc = 0.020, kBlend = 0.020;
const double kAlpha = 30 * kPi / 180, kTheta = 15 * kPi / 180, kEps = 10;
constexpr double kMixtureRatio = 5, kMassFlow = 0.4, kAmbientP = 1000, kAmbientT = 300;

struct Geometry {
  std::vector<std::array<double, 2>> contour;
  double cylinderEnd = 0;  // z where the cylinder meets the blend arc
  double throat = 0;
};

// Contour with the throat at z = 0, then shifted so the injector face is z = 0.
Geometry geometry() {
  const double re = kRt * std::sqrt(kEps);
  const double z1 = -kArc * std::sin(kAlpha), r1 = kRt + kArc * (1 - std::cos(kAlpha));
  const double rb = kRc - kBlend * (1 - std::cos(kAlpha)), zb = z1 - (rb - r1) / std::tan(kAlpha);
  const double zc = zb - kBlend * std::sin(kAlpha);
  const double z2 = kArc * std::sin(kTheta), r2 = kRt + kArc * (1 - std::cos(kTheta));
  const double ze = z2 + (re - r2) / std::tan(kTheta);
  const double z0 = -ze;  // throat at mid-length
  std::vector<std::array<double, 2>> c;
  auto add = [&](double z, double r) {
    if (c.empty() || z > c.back()[0] + 1e-12) c.push_back({z, r});
  };
  const int n = 400;
  add(z0, kRc);
  add(zc, kRc);
  for (int k = 1; k <= n; ++k) {  // blend arc, centre (zc, kRc - kBlend)
    double a = kAlpha * k / n;
    add(zc + kBlend * std::sin(a), kRc - kBlend + kBlend * std::cos(a));
  }
  add(z1, r1);  // convergent cone (straight)
  for (int k = 1; k <= 2 * n; ++k) {  // throat arc, centre (0, kRt + kArc)
    double a = -kAlpha + (kAlpha + kTheta) * k / (2 * n);
    add(kArc * std::sin(a), kRt + kArc - kArc * std::cos(a));
  }
  add(ze, re);  // conical exit
  Geometry g;
  for (auto& p : c) g.contour.push_back({p[0] - z0, p[1]});
  g.cylinderEnd = zc - z0;
  g.throat = -z0;
  return g;
}

double quantile(std::vector<double> v, double f) {
  if (v.empty()) return 0;
  const auto i = static_cast<std::size_t>(f * double(v.size() - 1));
  std::nth_element(v.begin(), v.begin() + static_cast<std::ptrdiff_t>(i), v.end());
  return v[i];
}

// The "cache" study (header).
int cacheStudy(thermo::ReactionSource& source, const Definition& c, const std::function<void(Flow&, const Definition&)>& initialise,
               int steps, int threads) {
  using Clock = std::chrono::steady_clock;
  auto since = [](Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); };
  Flow flow(c);
  initialise(flow, c);
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", threads, 1e-6, 1e-12, thermo::Chemistry::FiniteRate);
  for (int i = 0; i < 5; ++i) reacting.step();
  const std::size_t cells = flow.state().size(), ns = flow.medium().size(), nzv = ns + 1, nx = ns + 3;
  // The query stream, in march order: cell, rho, dt, and z = (T, Y) before and after.
  std::vector<std::size_t> cellOf;
  std::vector<double> rhoOf, dtOf, in, out;
  auto record = [&](double h) {
    for (std::size_t q = 0; q < cells; ++q) {
      cellOf.push_back(q);
      rhoOf.push_back(flow.state()[q][0]);
      dtOf.push_back(h);
      in.push_back(flow.temperature(q));
      for (double y : flow.massFractions(q)) in.push_back(y);
    }
    reacting.react(h);
    for (std::size_t q = 0; q < cells; ++q) {
      out.push_back(flow.temperature(q));
      for (double y : flow.massFractions(q)) out.push_back(y);
    }
  };
  const double t0sim = flow.time();
  for (int i = 0; i < steps; ++i) {
    // As ReactingFlow::step for FiniteRate, recording every reaction call.
    double planned = std::min(flow.stableDt(), flow.nextEvent(flow.time()) - flow.time());
    const auto saved = flow.partialDensities();
    for (int attempt = 0;; ++attempt) {
      record(0.5 * planned);
      const double limit = flow.stableDt();
      if (limit >= planned) break;
      if (attempt == 20) throw std::runtime_error("Reaction step keeps lowering the flow time step.");
      flow.setPartialDensities(saved);
      planned = limit;
    }
    record(0.5 * flow.step(planned));
  }
  const std::size_t n = cellOf.size();
  std::printf("cache study: viscous FiniteRate, %d steps (%.3e s simulated), %zu reaction calls of %zu cells, %zu queries\n", steps,
              flow.time() - t0sim, n / cells, cells, n);

  auto weighted = [&](const double* err, const double* z) {
    double e = 0;
    for (std::size_t k = 0; k < nzv; ++k) e = std::max(e, std::abs(err[k]) / (1e-6 * std::abs(z[k]) + 1e-12));
    return e;
  };
  auto x = [&](std::size_t i, std::size_t k) {
    if (k < nzv) return in[i * nzv + k];
    return k == nzv ? std::log(rhoOf[i]) : std::log(dtOf[i]);
  };

  // The integration error on the same scale.
  thermo::ReactionStep fine(source, 1e-10, 1e-14);
  std::vector<double> z(nzv), err(nzv), integration;
  double intT = 0, intY = 0;
  for (std::size_t i = 0; i < n; i += 97) {
    std::copy(in.begin() + std::ptrdiff_t(i * nzv), in.begin() + std::ptrdiff_t((i + 1) * nzv), z.begin());
    fine.advance(rhoOf[i], z.data(), dtOf[i]);
    for (std::size_t k = 0; k < nzv; ++k) err[k] = out[i * nzv + k] - z[k];
    integration.push_back(weighted(err.data(), z.data()));
    intT = std::max(intT, std::abs(err[0]));
    for (std::size_t k = 1; k < nzv; ++k) intY = std::max(intY, std::abs(err[k]));
  }
  std::printf("integration error (rtol 1e-6 march against rtol 1e-10, %zu calls): weighted 50%% %.3e, 99%% %.3e, largest %.3e; "
              "largest |dT| %.3e K, largest |dY| %.3e\n",
              integration.size(), quantile(integration, 0.5), quantile(integration, 0.99), quantile(integration, 1.0), intT, intY);

  // J = d(dz)/dx at a record, forward differences at rtol 1e-10, memoised over the levels.
  std::unordered_map<std::size_t, std::vector<double>> jacobians;
  auto jacobian = [&](std::size_t j) -> const std::vector<double>& {
    auto it = jacobians.find(j);
    if (it != jacobians.end()) return it->second;
    std::vector<double> jac(nzv * nx), base(nzv), pert(nzv);
    auto dz = [&](std::vector<double>& result, std::size_t k, double h) {
      std::copy(in.begin() + std::ptrdiff_t(j * nzv), in.begin() + std::ptrdiff_t((j + 1) * nzv), result.begin());
      double rho = rhoOf[j], dt = dtOf[j];
      if (k < nzv) result[k] += h;
      else if (k == nzv) rho *= std::exp(h);
      else if (k == nzv + 1) dt *= std::exp(h);
      const std::vector<double> start = result;
      fine.advance(rho, result.data(), dt);
      for (std::size_t m = 0; m < nzv; ++m) result[m] -= start[m];
    };
    dz(base, nx, 0);
    for (std::size_t k = 0; k < nx; ++k) {
      const double h = k == 0 ? 1e-2 : k < nzv ? 1e-6 : 1e-5;
      dz(pert, k, h);
      for (std::size_t m = 0; m < nzv; ++m) jac[m * nx + k] = (pert[m] - base[m]) / h;
    }
    return jacobians.emplace(j, std::move(jac)).first->second;
  };

  std::printf("%-7s %-6s %9s %9s %9s %11s %11s %11s %11s %11s %11s %7s %7s %7s %7s %7s %7s %12s %11s\n", "w", "table", "queries", "hits",
              "records", "const 50%", "const 99%", "const max", "lin 50%", "lin 99%", "lin max", "c<=1", "c<=10", "c<=100",
              "l<=1", "l<=10", "l<=100", "lin max dT K", "lin max dY");
  std::vector<double> width(nx, 1.0);
  width[0] = 1000;
  std::vector<std::int64_t> key(nx);
  std::string keyBytes(nx * sizeof(std::int64_t), '\0');
  auto keyOf = [&](std::size_t i, double w) -> const std::string& {
    for (std::size_t k = 0; k < nx; ++k) key[k] = static_cast<std::int64_t>(std::floor(x(i, k) / (k == nx - 1 ? 1e-2 : width[k] * w)));
    std::memcpy(keyBytes.data(), key.data(), keyBytes.size());
    return keyBytes;
  };
  std::vector<double> estimate(nzv);
  auto linear = [&](std::size_t i, std::size_t j, const std::vector<double>& jac) {
    for (std::size_t m = 0; m < nzv; ++m) {
      double v = out[j * nzv + m] - in[j * nzv + m];
      for (std::size_t k = 0; k < nx; ++k) v += jac[m * nx + k] * (x(i, k) - x(j, k));
      estimate[m] = v;
    }
  };
  for (double w : {1e-2, 1e-3, 1e-4, 1e-5})
    for (bool other : {false, true}) {
      // One record per bin: the bin's first query. Under "other", a query whose bin's record came
      // from its own cell is a miss (and the record stays).
      std::unordered_map<std::string, std::size_t> table;
      std::vector<double> constErr, linErr;
      std::size_t records = 0, hits = 0, c1 = 0, c10 = 0, c100 = 0, l1 = 0, l10 = 0, l100 = 0;
      double maxT = 0, maxY = 0;
      for (std::size_t i = 0; i < n; ++i) {
        auto [it, fresh] = table.try_emplace(keyOf(i, w), i);
        if (fresh) { ++records; continue; }
        const std::size_t j = it->second;
        if (other && cellOf[j] == cellOf[i]) continue;
        ++hits;
        const double* truth = &out[i * nzv];
        const double scale = dtOf[i] / dtOf[j];
        for (std::size_t m = 0; m < nzv; ++m) err[m] = scale * (out[j * nzv + m] - in[j * nzv + m]) - (truth[m] - in[i * nzv + m]);
        const double ec = weighted(err.data(), truth);
        linear(i, j, jacobian(j));
        for (std::size_t m = 0; m < nzv; ++m) err[m] = estimate[m] - (truth[m] - in[i * nzv + m]);
        const double el = weighted(err.data(), truth);
        maxT = std::max(maxT, std::abs(err[0]));
        for (std::size_t m = 1; m < nzv; ++m) maxY = std::max(maxY, std::abs(err[m]));
        constErr.push_back(ec);
        linErr.push_back(el);
        c1 += ec <= 1;
        c10 += ec <= 10;
        c100 += ec <= 100;
        l1 += el <= 1;
        l10 += el <= 10;
        l100 += el <= 100;
      }
      const double q = double(n);
      std::printf("%-7.0e %-6s %9zu %9zu %9zu %11.3e %11.3e %11.3e %11.3e %11.3e %11.3e %7.4f %7.4f %7.4f %7.4f %7.4f %7.4f %12.3e %11.3e\n", w,
                  other ? "other" : "any", n, hits, records, quantile(constErr, 0.5), quantile(constErr, 0.99),
                  quantile(constErr, 1.0), quantile(linErr, 0.5), quantile(linErr, 0.99), quantile(linErr, 1.0), c1 / q, c10 / q,
                  c100 / q, l1 / q, l10 / q, l100 / q, maxT, maxY);
      if (w == 1e-3 && !other) {
        // The cost of one retrieval: key, lookup and the linear sum, over every query on the built table.
        const auto t0 = Clock::now();
        double sink = 0;
        for (std::size_t i = 0; i < n; ++i) {
          const auto it = table.find(keyOf(i, w));
          const std::size_t j = it->second;
          const auto found = jacobians.find(j);
          if (found != jacobians.end()) linear(i, j, found->second);
          sink += estimate[0];
        }
        std::printf("retrieval: %.3e s per query (key, lookup, linear sum where the record has J; checksum %.3e)\n", since(t0) / q,
                    sink);
      }
      std::fflush(stdout);
    }
  std::printf("Jacobians computed: %zu\n", jacobians.size());
  return 0;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc < 6) {
    std::fprintf(stderr, "usage: %s <nz> <nr> <setup time s> <timed steps> <threads> [equilibrium | cache]\n", argv[0]);
    return 2;
  }
  const int nz = std::atoi(argv[1]), nr = std::atoi(argv[2]), steps = std::atoi(argv[4]), threads = std::atoi(argv[5]);
  const double setup = std::atof(argv[3]);
  const bool equilibrium = argc > 6 && std::string(argv[6]) == "equilibrium";
  const bool cache = argc > 6 && std::string(argv[6]) == "cache";
  const auto reactingMode = equilibrium ? thermo::Chemistry::LocalEquilibrium : thermo::Chemistry::FiniteRate;
  thermo::ReactionSource source("h2o2.yaml");
  const auto geo = geometry();
  Definition d;
  d.experiment = Case::Chamber;
  d.nz = nz;
  d.nr = nr;
  d.contour = geo.contour;
  d.species = source.medium().species();
  d.composition = source.massFractions("N2:1");
  d.backPressure = kAmbientP;
  d.ambientTemperature = kAmbientT;
  Supply s;
  s.outerRadius = kRc;
  s.massFlow = kMassFlow;
  s.totalTemperature = 300;
  auto h2 = source.massFractions("H2:1"), o2 = source.massFractions("O2:1");
  s.composition.resize(h2.size());
  for (std::size_t k = 0; k < h2.size(); ++k)
    s.composition[k] = (h2[k] + kMixtureRatio * o2[k]) / (1 + kMixtureRatio);
  s.opens = 0;
  s.ramp = 5e-4;
  d.supplies = {s};

  using Clock = std::chrono::steady_clock;
  auto since = [](Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); };
  Flow start(d);
  {
    thermo::ReactingFlow eq(start, "h2o2.yaml", threads, 1e-6, 1e-12, thermo::Chemistry::LocalEquilibrium);
    const auto t0 = Clock::now();
    eq.advanceTo(setup);
    const auto m = start.measurements();
    std::printf("setup %dx%d (%zu cells): LocalEquilibrium to %.3f ms, %lld steps, %.0f s wall; injector p %.4f MPa, outlet %.4f kg/s\n",
                nz, nr, start.state().size(), setup * 1e3, static_cast<long long>(m.steps), since(t0),
                m.injectorPressure / 1e6, m.outletMassFlow);
  }
  const std::size_t cells = start.state().size(), ns = start.medium().size();
  std::vector<Primitive> primitives(cells);
  std::vector<double> fractions(cells * ns), temperatures(cells);
  for (std::size_t q = 0; q < cells; ++q) {
    primitives[q] = start.cellPrimitive(q);
    temperatures[q] = start.temperature(q);
    auto y = start.massFractions(q);
    std::copy(y.begin(), y.end(), fractions.begin() + static_cast<std::ptrdiff_t>(q * ns));
  }

  const auto fits = thermo::transportFits("h2o2.yaml");
  // The timed configurations: inviscid as C1, or viscous (header).
  auto configure = [&](bool viscous) {
    Definition c = d;
    c.supplies[0].ramp = 1e-9;
    if (viscous) {
      c.transport = fits;
      c.wallTemperature = 600;
      c.turbulence.enabled = true;
      c.supplies[0].turbulenceIntensity = 0.05;
      c.supplies[0].viscosityRatio = 10;
      const Medium medium = c.medium();
      const auto y = c.composition;
      std::vector<double> diffusion(ns), work;
      const double a0 = medium.soundSpeed(kAmbientT, y.data()), rho0 = kAmbientP / (medium.gasConstant(y.data()) * kAmbientT);
      const double mu0 = medium.transport(kAmbientT, kAmbientP, y.data(), diffusion.data(), work).viscosity;
      c.turbulence.ambientK = 9e-9 * a0 * a0;
      c.turbulence.ambientOmega = 1e-6 * rho0 * a0 * a0 / mu0;
    }
    return c;
  };
  auto initialise = [&](Flow& flow, const Definition& c) {
    flow.setInitialState(primitives, fractions);
    if (!c.turbulence.enabled) return;
    const Medium medium = c.medium();
    std::vector<double> kOmega(2 * cells), diffusion(ns), work;
    for (std::size_t q = 0; q < cells; ++q) {
      const auto& w = primitives[q];
      const double mu = medium.transport(temperatures[q], w.p, fractions.data() + q * ns, diffusion.data(), work).viscosity;
      const double k = 1.5 * std::pow(0.05 * std::hypot(w.uz, w.ur), 2) + c.turbulence.ambientK;
      kOmega[2 * q] = k;
      kOmega[2 * q + 1] = w.rho * k / (10 * mu);
    }
    flow.setTurbulence(kOmega);
  };
  if (cache) return cacheStudy(source, configure(true), initialise, steps, threads);

  std::printf("%-9s %-10s %7s %12s %12s %14s %16s %12s\n", "config", "chemistry", "threads", "dt start [s]", "mean dt [s]",
              "wall/step [s]", "cell updates/s", "replans/step");
  for (bool viscous : {false, true})
    for (auto chemistry : {thermo::Chemistry::Frozen, reactingMode})
      for (int t : {1, threads}) {
        if (t != 1 && (chemistry == thermo::Chemistry::Frozen || threads == 1)) continue;
        const Definition c = configure(viscous);
        Flow flow(c);
        initialise(flow, c);
        thermo::ReactingFlow reacting(flow, "h2o2.yaml", t, 1e-6, 1e-12, chemistry);
        const double dt0 = flow.stableDt();
        for (int i = 0; i < 5; ++i) reacting.step();
        const double sim0 = flow.time();
        const long replans0 = reacting.stats().replans;
        const auto t0 = Clock::now();
        for (int i = 0; i < steps; ++i) reacting.step();
        const double wall = since(t0);
        std::printf("%-9s %-10s %7d %12.4e %12.4e %14.4e %16.4e %12.2f\n", viscous ? "viscous" : "inviscid",
                    chemistry == thermo::Chemistry::Frozen ? "frozen" : equilibrium ? "equil" : "finite", t, dt0, (flow.time() - sim0) / steps,
                    wall / steps, double(cells) * steps / wall, double(reacting.stats().replans - replans0) / steps);
        std::fflush(stdout);
      }
  if (equilibrium) return 0;
  std::printf("%-26s %16s\n", "FiniteRate substep [s]", "s per cell");
  Definition open = d;
  open.supplies[0].ramp = 1e-9;
  Flow flow(open);
  flow.setInitialState(primitives, fractions);
  thermo::ReactingFlow reacting(flow, "h2o2.yaml", 1, 1e-6, 1e-12, thermo::Chemistry::FiniteRate);
  for (int i = 0; i < 20; ++i) reacting.step();
  const double half = 0.5 * flow.stableDt();
  for (double h : {half, 5e-10, 5e-11}) {
    reacting.react(h);
    const auto t0 = Clock::now();
    for (int i = 0; i < 20; ++i) reacting.react(h);
    std::printf("%-26.4e %16.4e\n", h, since(t0) / (20.0 * double(cells)));
  }
  return 0;
}
