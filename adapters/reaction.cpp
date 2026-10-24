#include "adapters/reaction.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "adapters/thermo.hpp"
#include "cantera/core.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/ext/cvodes/cvodes.h"
#include "cantera/ext/nvector/nvector_serial.h"
#include "cantera/ext/sunlinsol/sunlinsol_dense.h"
#include "cantera/ext/sunmatrix/sunmatrix_dense.h"

namespace crucible::thermo {

struct ReactionSource::Impl {
  std::shared_ptr<Cantera::Solution> solution;
  std::shared_ptr<Cantera::ThermoPhase> gas;
  std::shared_ptr<Cantera::Kinetics> kinetics;
  std::vector<double> wdot, uk;
  // jacobian(): Cantera's derivatives, molar cv_k, and each species' NASA-7 data as reportParams gives it
  // (c[0] the mid temperature, c[1..7] the high range, c[8..14] the low range); empty until first used.
  std::vector<double> dwdT, dwdP, cvk;
  std::vector<std::array<double, 15>> nasa;

  // Put the phase at (T, rho, Y) without renormalising Y (integrator iterates may not sum to 1).
  void set(double t, double rho, const double* y) {
    gas->setMassFractions_NoNorm(y);
    gas->setState_TD(t, rho);
  }
};

ReactionSource::ReactionSource(const std::string& mechanism) : impl_(std::make_unique<Impl>()) {
  Cantera::addDataDirectory(dataDirectory());
  impl_->solution = Cantera::newSolution(mechanism, "", "none");
  impl_->gas = impl_->solution->thermo();
  impl_->kinetics = impl_->solution->kinetics();
  impl_->wdot.resize(nSpecies());
  impl_->uk.resize(nSpecies());
}

ReactionSource::~ReactionSource() = default;

std::size_t ReactionSource::nSpecies() const { return impl_->gas->nSpecies(); }
std::size_t ReactionSource::nElements() const { return impl_->gas->nElements(); }
std::string ReactionSource::elementName(std::size_t m) const { return impl_->gas->elementName(m); }
std::string ReactionSource::speciesName(std::size_t k) const { return impl_->gas->speciesName(k); }

std::vector<ReactionSource::ReactionSpecies> ReactionSource::reactionSpecies() const {
  auto& kinetics = *impl_->kinetics;
  std::vector<ReactionSpecies> out(kinetics.nReactions());
  for (std::size_t i = 0; i < out.size(); ++i) {
    for (std::size_t k = 0; k < nSpecies(); ++k) {
      const double nuR = kinetics.reactantStoichCoeff(k, i), nuP = kinetics.productStoichCoeff(k, i);
      if (nuR != 0) out[i].reactants.push_back(k);
      if (nuP != 0) out[i].products.push_back(k);
      if (nuP != nuR) out[i].changed.push_back(k);
    }
    out[i].reversible = kinetics.reaction(i)->reversible;
  }
  return out;
}

Medium ReactionSource::medium() const {
  auto& gas = *impl_->gas;
  std::vector<Species> species;
  for (std::size_t k = 0; k < nSpecies(); ++k) {
    int type;
    double tlow, thigh, pref, c[15];
    gas.speciesThermo().reportParams(k, type, c, tlow, thigh, pref);
    if (type != NASA2) throw std::runtime_error("Species " + speciesName(k) + " is not NASA-7 data");
    Species s{speciesName(k), gas.molecularWeight(k), c[0], {}, {}};
    std::copy(c + 8, c + 15, s.low.begin());
    std::copy(c + 1, c + 8, s.high.begin());
    // Fits that reach 6,000 K (h2o2_highT) get the constant-cp continuation above their top, labelled
    // outside the model: there, extrapolated, they turn cv negative from 10,630 K (OH) and e(T) cannot
    // be inverted. Fits to 3,500 K (GRI h2o2.yaml) keep their extrapolation, which is monotone to
    // 40,000 K, so every result on them is unchanged.
    if (thigh >= 6000) s.tHigh = thigh;
    species.push_back(s);
  }
  return Medium(species);
}

TransportFits transportFits(const std::string& mechanism, double minTemperature) {
  std::shared_ptr<Cantera::Solution> solution;
  if (minTemperature > 0) {
    // Cantera fits over the phase's temperature range, from the highest of its species' lower limits.
    // Each species' NASA-7 range is extended down to minTemperature (its low polynomial extrapolated, as
    // Medium evaluates it), so the kinetic theory behind the fits is sampled there too.
    auto root = Cantera::AnyMap::fromYamlFile(mechanism);
    for (auto& species : root["species"].asVector<Cantera::AnyMap>()) {
      auto& ranges = species["thermo"]["temperature-ranges"].asVector<double>();
      ranges.front() = std::min(ranges.front(), minTemperature);
    }
    solution = Cantera::newSolution(root["phases"].asVector<Cantera::AnyMap>().front(), root, "mixture-averaged");
  } else {
    solution = Cantera::newSolution(mechanism, "", "mixture-averaged");
  }
  auto& transport = *solution->transport();
  if (transport.CKMode()) throw std::runtime_error("CHEMKIN-mode transport fits are not supported");
  const std::size_t n = solution->thermo()->nSpecies();
  TransportFits fits;
  fits.viscosity.resize(n);
  fits.conductivity.resize(n);
  for (std::size_t k = 0; k < n; ++k) {
    transport.getViscosityPolynomial(k, fits.viscosity[k].data());
    transport.getConductivityPolynomial(k, fits.conductivity[k].data());
    for (std::size_t j = k; j < n; ++j) {
      fits.diffusion.emplace_back();
      transport.getBinDiffusivityPolynomial(k, j, fits.diffusion.back().data());
    }
  }
  return fits;
}

std::vector<double> ReactionSource::massFractions(const std::string& moles) {
  impl_->gas->setMoleFractionsByName(moles);
  std::vector<double> y(nSpecies());
  impl_->gas->getMassFractions(y.data());
  return y;
}

double ReactionSource::density(double t, double p, const double* y) {
  impl_->gas->setMassFractions_NoNorm(y);
  impl_->gas->setState_TP(t, p);
  return impl_->gas->density();
}

double ReactionSource::pressure(double t, double rho, const double* y) {
  impl_->set(t, rho, y);
  return impl_->gas->pressure();
}

double ReactionSource::internalEnergy(double t, double rho, const double* y) {
  impl_->set(t, rho, y);
  return impl_->gas->intEnergy_mass();
}

double ReactionSource::cv(double t, double rho, const double* y) {
  impl_->set(t, rho, y);
  return impl_->gas->cv_mass();
}

std::vector<double> ReactionSource::elementMassFractions(const double* y) {
  auto& gas = *impl_->gas;
  std::vector<double> e(nElements(), 0.0);
  for (std::size_t m = 0; m < e.size(); ++m)
    for (std::size_t k = 0; k < nSpecies(); ++k)
      e[m] += y[k] * gas.nAtoms(k, m) * gas.atomicWeight(m) / gas.molecularWeight(k);
  return e;
}

void ReactionSource::rates(double rho, const double* z, double* dzdt) {
  auto& gas = *impl_->gas;
  impl_->set(z[0], rho, z + 1);
  impl_->kinetics->getNetProductionRates(impl_->wdot.data());  // [kmol/(m^3 s)]
  gas.getPartialMolarIntEnergies(impl_->uk.data());            // [J/kmol]
  double heat = 0.0;
  for (std::size_t k = 0; k < nSpecies(); ++k) {
    dzdt[k + 1] = impl_->wdot[k] * gas.molecularWeight(k) / rho;
    heat += impl_->uk[k] * impl_->wdot[k];
  }
  dzdt[0] = -heat / (rho * gas.cv_mass());
}

void ReactionSource::jacobian(double rho, const double* z, double* jac, double* dzdt) {
  auto& gas = *impl_->gas;
  auto& kinetics = *impl_->kinetics;
  const std::size_t ns = nSpecies(), n = ns + 1;
  if (impl_->nasa.empty()) {
    // Exact derivatives: third-body and falloff terms in. Cantera takes one-sided differences in T and [M]
    // for falloff and the equilibrium constants; the rest is analytic. Their step is 1e-7 relative, not the
    // default 1e-8: against converged differences, 1e-8 leaves round-off of 1.0e-6 in dT/dT at 2,500 K,
    // and the worst column is 3.1e-7 (h2o2) and 4.4e-7 (h2o2_highT) at 1e-7 (measured, TABLE_B.md).
    Cantera::AnyMap settings;
    settings["skip-third-bodies"] = false;
    settings["skip-falloff"] = false;
    settings["rtol-delta"] = 1e-7;
    kinetics.setDerivativeSettings(settings);
    impl_->dwdT.resize(ns);
    impl_->dwdP.resize(ns);
    impl_->cvk.resize(ns);
    impl_->nasa.resize(ns);
    for (std::size_t k = 0; k < ns; ++k) {
      int type;
      double tlow, thigh, pref;
      gas.speciesThermo().reportParams(k, type, impl_->nasa[k].data(), tlow, thigh, pref);
      if (type != NASA2) throw std::runtime_error("Species " + speciesName(k) + " is not NASA-7 data");
    }
  }
  const double t = z[0];
  const double* y = z + 1;
  impl_->set(t, rho, y);
  // Cantera 3.2's getNetProductionRates_ddT reads a stale state after a composition change at unchanged
  // (T, rho): the falloff reactions' species rows come out wrong by up to 5e4 times (measured, TABLE_B.md,
  // "The analytic Jacobian"). Invalidating the kinetics cache gives the fresh values exactly.
  kinetics.invalidateCache();
  auto& wdot = impl_->wdot;
  auto& uk = impl_->uk;
  auto& cvk = impl_->cvk;
  kinetics.getNetProductionRates(wdot.data());       // [kmol/(m^3 s)]
  kinetics.getNetProductionRates_ddT(impl_->dwdT.data());  // concentrations held
  kinetics.getNetProductionRates_ddP(impl_->dwdP.data());  // T and concentrations held
  const Eigen::SparseMatrix<double> dwdc = kinetics.netProductionRates_ddCi();  // (k, j): d w_k / d C_j
  gas.getPartialMolarIntEnergies(uk.data());         // [J/kmol]
  gas.getPartialMolarCp(cvk.data());
  const double gasConstant = Cantera::GasConstant, p = gas.pressure();
  for (double& c : cvk) c -= gasConstant;            // molar cv_k of the ideal gas
  const double cv = gas.cv_mass();
  double heat = 0.0, heatT = 0.0, dcvdT = 0.0;
  for (std::size_t k = 0; k < ns; ++k) {
    const double wk = gas.molecularWeight(k);
    const double dwdT = impl_->dwdT[k] + impl_->dwdP[k] * p / t;  // at fixed volume, P = R T sum C
    heat += uk[k] * wdot[k];
    heatT += cvk[k] * wdot[k] + uk[k] * dwdT;
    const auto& c = impl_->nasa[k];
    const double* a = t <= c[0] ? c.data() + 8 : c.data() + 1;
    dcvdT += y[k] / wk * gasConstant * (a[1] + t * (2 * a[2] + t * (3 * a[3] + t * 4 * a[4])));
    jac[k + 1] = wk / rho * dwdT;  // column 0: d/dT
  }
  const double fT = -heat / (rho * cv);
  jac[0] = -heatT / (rho * cv) - fT * dcvdT / cv;
  // Columns 1..K: d/dY_j = (rho / W_j) d/dC_j, plus the pressure term R T dw/dP.
  for (std::size_t j = 0; j < ns; ++j) {
    double* column = jac + (j + 1) * n;
    const double wj = gas.molecularWeight(j);
    double sum = 0.0;
    for (std::size_t k = 0; k < ns; ++k) {
      const double d = gasConstant * t * impl_->dwdP[k];
      column[k + 1] = gas.molecularWeight(k) / wj * d;
      sum += uk[k] * d;
    }
    for (Eigen::SparseMatrix<double>::InnerIterator it(dwdc, static_cast<Eigen::Index>(j)); it; ++it) {
      const auto k = static_cast<std::size_t>(it.row());
      column[k + 1] += gas.molecularWeight(k) / wj * it.value();
      sum += uk[k] * it.value();
    }
    column[0] = -(sum + fT * cvk[j]) / (wj * cv);
  }
  if (dzdt) {
    for (std::size_t k = 0; k < ns; ++k) dzdt[k + 1] = wdot[k] * gas.molecularWeight(k) / rho;
    dzdt[0] = fT;
  }
}

void ReactionSource::equilibrateUV(double rho, double* z) { equilibrate(rho, z, "UV"); }
void ReactionSource::equilibrateTV(double rho, double* z) { equilibrate(rho, z, "TV"); }

void ReactionSource::equilibrate(double rho, double* z, const char* constraints) {
  // Equilibrium depends only on the elements; integrator round-off negatives (~1e-20)
  // are clipped so the solver starts from an admissible composition.
  auto& gas = *impl_->gas;
  const std::size_t ns = nSpecies(), ne = nElements();
  std::vector<double> y(z + 1, z + 1 + ns);
  for (double& v : y) v = std::max(v, 0.0);
  // Elements below kTraceElement (mass fraction) are left out of the problem: the species that
  // carry them keep their amounts. Advection leaves such traces (1e-50 and below) in every cell,
  // and the equilibrium solvers can take ~0.5 s on them. The energy they could release is below
  // 1e-13 of the cell's internal energy, far inside the equilibrium tolerance.
  constexpr double kTraceElement = 1e-14;
  auto elements = elementMassFractions(y.data());
  std::vector<bool> frozen(ns, false);
  for (std::size_t m = 0; m < ne; ++m)
    if (elements[m] < kTraceElement)
      for (std::size_t k = 0; k < ns; ++k)
        if (gas.nAtoms(k, m) > 0) frozen[k] = true;
  // If every remaining element is carried by exactly one remaining species the composition is
  // already the equilibrium one (for example N2 alone), and so is the temperature.
  bool trivial = true;
  for (std::size_t m = 0; m < ne && trivial; ++m) {
    if (elements[m] < kTraceElement) continue;
    int carriers = 0;
    for (std::size_t k = 0; k < ns; ++k) carriers += !frozen[k] && gas.nAtoms(k, m) > 0;
    trivial = carriers == 1;
  }
  if (trivial) return;
  double frozenMass = 0;
  std::vector<double> active(y);
  for (std::size_t k = 0; k < ns; ++k)
    if (frozen[k]) { frozenMass += y[k]; active[k] = 0; }
  impl_->set(z[0], rho, active.data());
  gas.equilibrate(constraints);
  z[0] = gas.temperature();
  gas.getMassFractions(active.data());
  for (std::size_t k = 0; k < ns; ++k) z[k + 1] = frozen[k] ? y[k] : active[k] * (1 - frozenMass);
}

namespace {

struct Context {
  ReactionSource* source;
  double rho;
  Mixing mixing;  // segregation 0: no closure
  const std::vector<std::size_t>* species = nullptr;
  std::size_t n = 0;
};

int rhs(realtype, N_Vector z, N_Vector dz, void* data) {
  auto* c = static_cast<Context*>(data);
  try {
    c->source->rates(c->rho, NV_DATA_S(z), NV_DATA_S(dz));
  } catch (const Cantera::CanteraError&) {
    return 1;  // recoverable: CVODES retries with a smaller step
  }
  if (c->mixing.segregation > 0) {
    const double f = reactingFraction(NV_DATA_S(z), NV_DATA_S(dz), *c->species, c->mixing);
    for (std::size_t i = 0; i < c->n; ++i) NV_Ith_S(dz, i) *= f;
  }
  return 0;
}

// The Newton matrix: ReactionSource::jacobian, times kappa_eff (held) under the PaSR closure.
int jacobian(realtype, N_Vector z, N_Vector, SUNMatrix jac, void* data, N_Vector work, N_Vector, N_Vector) {
  auto* c = static_cast<Context*>(data);
  try {
    c->source->jacobian(c->rho, NV_DATA_S(z), SM_DATA_D(jac), NV_DATA_S(work));
  } catch (const Cantera::CanteraError&) {
    return 1;  // recoverable: CVODES retries with a smaller step
  }
  if (c->mixing.segregation > 0) {
    const double f = reactingFraction(NV_DATA_S(z), NV_DATA_S(work), *c->species, c->mixing);
    const std::size_t entries = c->n * c->n;
    for (std::size_t i = 0; i < entries; ++i) SM_DATA_D(jac)[i] *= f;
  }
  return 0;
}

}  // namespace

double reactingFraction(const double* z, const double* dzdt, const std::vector<std::size_t>& species,
                        const Mixing& mixing) {
  if (mixing.segregation == 0) return 1;
  double chemical = 0;
  bool active = false;
  for (std::size_t i : species) {
    const double y = z[i + 1], rate = std::abs(dzdt[i + 1]);
    if (!(y > 0) || rate == 0) continue;
    chemical = std::max(chemical, y / rate);
    active = true;
  }
  if (!active) return 1;
  const double kappa = chemical / (chemical + mixing.time);
  return 1 - mixing.segregation * (1 - kappa);
}

struct ReactionStep::Impl {
  ReactionSource* source;
  Context context;
  sunindextype n;
  N_Vector z = nullptr;
  SUNMatrix jac = nullptr;
  SUNLinearSolver linear = nullptr;
  void* cvode = nullptr;
  long steps = 0;

  ~Impl() {
    if (cvode) CVodeFree(&cvode);
    if (linear) SUNLinSolFree(linear);
    if (jac) SUNMatDestroy(jac);
    if (z) N_VDestroy_Serial(z);
  }
};

ReactionStep::ReactionStep(ReactionSource& source, double rtol, double atol, Jacobian jacobianKind)
    : impl_(std::make_unique<Impl>()) {
  auto& s = *impl_;
  s.source = &source;
  s.context = {&source, 0.0, {}, nullptr, source.nSpecies() + 1};
  s.n = static_cast<sunindextype>(source.nSpecies() + 1);
  s.z = N_VNew_Serial(s.n);
  N_VConst(300.0, s.z);
  s.cvode = CVodeCreate(CV_BDF);
  if (!s.cvode || CVodeInit(s.cvode, rhs, 0.0, s.z) != CV_SUCCESS)
    throw std::runtime_error("CVODES initialisation failed");
  // Species use the given absolute tolerance; temperature [K] gets a floor of rtol kelvin.
  N_Vector abs = N_VNew_Serial(s.n);
  N_VConst(atol, abs);
  NV_Ith_S(abs, 0) = rtol;
  CVodeSVtolerances(s.cvode, rtol, abs);
  N_VDestroy_Serial(abs);
  CVodeSetUserData(s.cvode, &s.context);
  CVodeSetMaxNumSteps(s.cvode, 100000);
  // CVODES's first step from a cold start in a dense radical pool (200 K, recombining in 5e-14 to 2e-10 s)
  // is 3e-7 to 2e-6 s, and its default 10 corrector failures, each cutting h by 4, stop 10 to 100 times
  // short of a step Newton can take (measured, TABLE_B.md). 40 reach 1e-24 of the first step. A run that
  // never fails 10 times in one step is unchanged.
  CVodeSetMaxConvFails(s.cvode, 40);
  s.jac = SUNDenseMatrix(s.n, s.n);
  s.linear = SUNLinSol_Dense(s.z, s.jac);
  CVodeSetLinearSolver(s.cvode, s.linear, s.jac);
  if (jacobianKind == Jacobian::Analytic) CVodeSetJacFn(s.cvode, jacobian);
}

ReactionStep::~ReactionStep() = default;

void ReactionStep::advance(double rho, double* z, double dt, int samples, const Observer& observe) {
  advance(rho, z, dt, Mixing{}, {}, samples, observe);
}

void ReactionStep::advance(double rho, double* z, double dt, const Mixing& mixing,
                           const std::vector<std::size_t>& species, int samples, const Observer& observe) {
  auto& s = *impl_;
  if (!(mixing.segregation >= 0 && mixing.segregation <= 1) || !(mixing.time >= 0))
    throw std::invalid_argument("PaSR mixing inputs out of range.");
  for (std::size_t i : species)
    if (i + 1 >= static_cast<std::size_t>(s.n)) throw std::invalid_argument("PaSR species index out of range.");
  s.context.rho = rho;
  s.context.mixing = mixing;
  s.context.species = &species;
  for (sunindextype i = 0; i < s.n; ++i) NV_Ith_S(s.z, i) = z[i];
  if (CVodeReInit(s.cvode, 0.0, s.z) != CV_SUCCESS) throw std::runtime_error("CVODES reinit failed");
  CVodeSetStopTime(s.cvode, dt);
  realtype t = 0.0;
  for (int i = 1; i <= samples; ++i) {
    const double target = dt * i / samples;
    const int flag = CVode(s.cvode, target, s.z, &t, CV_NORMAL);
    if (flag < 0) throw std::runtime_error("CVODES failed with flag " + std::to_string(flag));
    if (observe) observe(t, NV_DATA_S(s.z));
  }
  for (sunindextype i = 0; i < s.n; ++i) z[i] = NV_Ith_S(s.z, i);
  CVodeGetNumSteps(s.cvode, &s.steps);
}

long ReactionStep::steps() const { return impl_->steps; }

double ReactionStep::follow(double rho, double* z, double tEnd, const StepObserver& observe) {
  auto& s = *impl_;
  s.context.rho = rho;
  s.context.mixing = Mixing{};  // no PaSR closure: the fallback's own rate, as the table's build wants
  s.context.species = nullptr;
  for (sunindextype i = 0; i < s.n; ++i) NV_Ith_S(s.z, i) = z[i];
  if (CVodeReInit(s.cvode, 0.0, s.z) != CV_SUCCESS) throw std::runtime_error("CVODES reinit failed");
  CVodeSetStopTime(s.cvode, tEnd);
  realtype t = 0.0;
  constexpr long kMaxInternalSteps = 2000000;
  long taken = 0;
  for (;;) {
    const int flag = CVode(s.cvode, tEnd, s.z, &t, CV_ONE_STEP);
    if (flag < 0) throw std::runtime_error("CVODES failed with flag " + std::to_string(flag));
    if (++taken > kMaxInternalSteps) throw std::runtime_error("ReactionStep::follow exceeded its internal step cap");
    const bool keepGoing = observe ? observe(t, NV_DATA_S(s.z)) : true;
    if (!keepGoing || flag == CV_TSTOP_RETURN || t >= tEnd) break;
  }
  for (sunindextype i = 0; i < s.n; ++i) z[i] = NV_Ith_S(s.z, i);
  CVodeGetNumSteps(s.cvode, &s.steps);
  return t;
}

}  // namespace crucible::thermo
