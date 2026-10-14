#pragma once
#include "core/flow.hpp"
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace crucible {
// Liquid oxygen spray: Lagrangian parcels in the axisymmetric gas of a Flow (TECHNICAL_PLAN step 5;
// design, limits and criteria in docs/evidence/SPRAY.md, stated before this code). Subcritical only:
// every cell a parcel meets must lie below oxygen's critical pressure.

// Saturated liquid oxygen from the NIST WebBook table (tools/o2_saturation_nist.tsv, rows uniform in
// T from the triple point to just below the critical point), linear in T between rows (ln p for the
// vapour pressure). Enthalpies are on the engine's scale: h = h_medium,O2(298.15 K) + h_NIST - h_ig,NIST(298.15 K),
// so liquid and gas share one zero. The latent heat is NIST's own, h_v - h_l on saturation: the vapour
// leaves the drop with the real saturated-vapour enthalpy and the ideal-gas medium carries it from there.
class LiquidOxygen {
 public:
  static constexpr double idealGasEnthalpy298 = 271260.1698;  // J/kg, NIST, 298.15 K at 1e-6 MPa (DEF)
  static constexpr double criticalPressure = 5.043e6;          // Pa
  // gasEnthalpy298: the medium's O2 enthalpy at 298.15 K [J/kg].
  static LiquidOxygen load(const std::string& path, double gasEnthalpy298);
  struct State {
    double density{}, enthalpy{}, cp{}, viscosity{}, conductivity{}, tension{}, pressure{}, latent{};
  };
  // At t, clamped to the table.
  [[nodiscard]] State at(double t) const;
  [[nodiscard]] double enthalpy(double t) const;
  [[nodiscard]] double vapourPressureSlope(double t) const;     // d ln p_sat / dT [1/K]
  [[nodiscard]] double temperature(double h) const;               // exact inverse of enthalpy()
  [[nodiscard]] double saturationTemperature(double p) const;     // exact inverse of the ln p interpolant
  [[nodiscard]] double minTemperature() const { return t0_; }
  [[nodiscard]] double maxTemperature() const { return t0_ + step_ * static_cast<double>(rows_.size() - 1); }

 private:
  // T, ln p, rho_l, h_l, cp_l, mu_l, k_l, sigma, h_v (enthalpies on the engine scale, SI).
  std::vector<std::array<double, 9>> rows_;
  double t0_{}, step_{};
  [[nodiscard]] std::array<double, 9> row(double t) const;
};

// One parcel: `count` identical drops. Position and velocity in the meridional half-plane through the
// parcel (z axial, r >= 0 radial, ut azimuthal); the gas has no swirl.
struct Parcel {
  double z{}, r{}, uz{}, ur{}, ut{};
  double d{}, t{};       // drop diameter [m], uniform drop temperature [K]
  double count{};        // drops in the parcel
  double shed{};         // KH: mass shed per drop and not yet released as a child parcel [kg]
  double rtTime{};       // RT: time its waves have grown at the present diameter [s]
  double distance{};     // path length travelled since injection [m]
  double injectedMass{}; // the parcel's liquid mass when injected or added [kg]; children inherit it
  // Turbulent dispersion: the eddy's velocity fluctuation (axial, radial, azimuthal in the parcel's
  // meridional plane) and the time left of its interaction [s]; 0 when it has none.
  double eddyZ{}, eddyR{}, eddyT{}, eddyTime{};
  std::uint64_t random{};  // the parcel's own random stream (splitmix64 state)
  int cell{-1};
  int injector{-1};      // its injector (sets the RT breakup length), -1 for a parcel placed by add()
  [[nodiscard]] double dropMass(double density) const;
};

class Spray {
 public:
  // A liquid stream released at the injector face plane z (0 for the face itself) across rInner <= r <= rOuter,
  // with Supply's valve schedule (opens, ramp). velocity 0: the mean speed massFlow / (rho_l area).
  // Blob: drops of the release annulus's hydraulic diameter 2 (rOuter - rInner) and
  // breakup does the atomisation. RosinRammler: diameters from a Rosin-Rammler distribution with the
  // given Sauter mean diameter and exponent q, cut at 3 times the SMD; with coreLength > 0 they are
  // released from the liquid core's surface (r = rOuter) uniformly in z over [z, z + coreLength].
  struct Injector {
    enum class Size { Blob, RosinRammler };
    double z{}, rInner{}, rOuter{};
    double massFlow{}, opens{}, ramp{};
    double temperature{90}, velocity{}, coneHalfAngle{};
    Size size{Size::Blob};
    double sauterDiameter{}, exponent{3}, coreLength{};
    double parcelRate{2e6};  // parcels released per second [1/s]
    [[nodiscard]] double opening(double time) const;
  };
  struct Breakup {
    bool on{true};
    // KH releases a child each time a parcel's shed mass reaches shedFraction of its injected mass (Reitz 1987).
    double b0{0.61}, b1{10}, cRt{0.1}, cTau{1}, cBl{1}, shedFraction{0.03};
  };
  struct Settings {
    std::vector<Injector> injectors;
    Breakup breakup;
    double restitution{1};      // normal coefficient of restitution at walls (1: specular)
    double minDiameter{1e-7};   // below it a parcel's liquid becomes vapour at once [m]
    // Gauss points per piece of the film integrals that give the conductances (SPRAY.md, Results: S2);
    // 0: the 1/3 rule's properties alone, as Abramzon and Sirignano.
    int filmPoints{3};
    // Turbulent dispersion (SPRAY.md, Turbulent dispersion), when the flow carries SST's k and omega.
    bool dispersion{true};
    bool implicitCoupling{true};
    // false: one-way, the drop model alone in a fixed gas (tests): the exchange is booked in the ledger
    // but not given to the gas, and the seen gas takes no implicit loads.
    bool twoWay{true};
    std::uint64_t seed{1};
  };
  // Running totals since construction [kg, J, kg m/s]. Liquid energy is enthalpy plus kinetic energy.
  // exchanged*: to the gas (Flow::exchange). wallImpulse: the axial momentum the walls gave the parcels,
  // which is the device's thrust impulse from them; wallLoss: kinetic energy left in the walls (restitution < 1).
  struct Ledger {
    double injectedMass{}, injectedEnergy{}, injectedMomentum{};
    double exchangedMass{}, exchangedEnergy{}, exchangedMomentum{};
    double outflowMass{}, outflowEnergy{}, outflowMomentum{};
    double wallImpulse{}, wallLoss{};
    double droppedSwirl{};       // |azimuthal momentum| passed to the gas as heat (no swirl in the gas)
    long parcelsInjected{}, children{}, wallHits{}, flashes{}, vaporised{}, clampedCondensation{};
    long substeps{};
    // The last advance(): injector thrust and wall thrust over its dt [N].
    double injectionThrust{}, wallThrust{};
  };

  Spray(Flow& flow, LiquidOxygen liquid, Settings settings);
  // Moves the parcels across the flow's last step dt (call after Flow::step), exchanges mass, momentum
  // and energy with the gas through Flow::exchange, then releases the step's new parcels.
  void advance(double dt);
  // Next injector schedule event after `time` (infinity if none).
  [[nodiscard]] double nextEvent(double time) const;
  [[nodiscard]] const std::vector<Parcel>& parcels() const { return parcels_; }
  [[nodiscard]] const Ledger& ledger() const { return ledger_; }
  [[nodiscard]] double wallTime() const { return wall_; }  // wall time spent in advance() [s]
  // The last advance()'s exchange per cell, Flow::exchange's input: ns species masses, axial and radial
  // momentum, energy [kg, kg m/s, J] (ns + 3 per cell).
  [[nodiscard]] const std::vector<double>& exchanged() const { return amounts_; }
  [[nodiscard]] const LiquidOxygen& liquid() const { return liquid_; }
  // Liquid held now (parcels plus mass the injectors have delivered but not yet released).
  [[nodiscard]] double liquidMass() const;
  [[nodiscard]] double liquidEnergy() const;
  [[nodiscard]] double liquidMomentum() const;
  // Relative closure of the liquid's mass, energy and axial momentum ledgers.
  [[nodiscard]] std::array<double, 3> balanceErrors() const;
  // For tests: a parcel placed directly (its cell is found here).
  void add(Parcel parcel);

  // The drop model alone: one drop in a fixed gas state (cell-centred values, no coupling).
  struct Gas {
    double t{}, p{}, rho{}, uz{}, ur{};
    const double* y{};  // mass fractions (medium order)
    double ut{};        // azimuthal velocity (an eddy's; the mean gas has none)
  };
  struct Rates {
    double mass{};      // evaporation rate per drop [kg/s] (>= 0)
    double heat{};      // heat conducted from the gas to the surface per drop [W]
    double conductance{};  // heat / (T_gas - T_surface) [W/K] (implicit coupling)
    double surface{};   // surface temperature [K]: the drop's, or saturation when boiling
    double stiffness{}; // d(Q - mdot L)/dT magnitude per drop [W/K], for the sub-step bound
    double tau{};       // momentum relaxation time [s]
    double reynolds{}, spalding{}, transfer{}, sherwood{}, nusselt{}, dragCoefficient{};
    bool boiling{}, clamped{};  // clamped: the gas holds more O2 than the surface (B_M set to 0)
  };
  // work: at least 4 ns + 8 doubles.
  [[nodiscard]] Rates rates(const Parcel& parcel, const Gas& gas, std::vector<double>& work) const;
  // Medium::transport's conductivity and O2's mixture diffusion coefficient alone, by the same formulas
  // (all the film integrals need), with the mixture gas constant; absent species are skipped (Medium
  // floors their mole fractions at 1e-20). S1(c) holds it to Medium::transport.
  struct FilmTransport { double conductivity{}, diffusion{}, r{}; };
  [[nodiscard]] FilmTransport filmTransport(double t, double p, const double* y) const;
  // Gosman and Ioannides' interaction time with an eddy of the gas's k and omega, for a parcel of
  // momentum relaxation time tau at slip speed `slip` against the gas plus the eddy, with the eddy
  // length that gives tracers the gas's nu_t / Sc_t (SPRAY.md, Turbulent dispersion).
  [[nodiscard]] static double interaction(double k, double omega, double tau, double slip, double schmidt);
  // g(Re, B) in C_D = (24 / Re) g (SPRAY.md, Literature).
  [[nodiscard]] static double dragFactor(double reynolds, double transfer);
  // KH-RT for a drop of diameter d at slip speed u: the KH wavelength and growth rate, and the RT
  // wavelength (with C_RT) and growth rate at the drop's drag deceleration g_t (SPRAY.md, Literature).
  struct Waves { double khLength{}, khRate{}, rtLength{}, rtRate{}; };
  [[nodiscard]] static Waves waves(double d, double u, double rhoGas, double rhoLiquid, double muLiquid,
                                   double tension, double deceleration, double cRt);

 private:
  struct Deposit { int cell; double mass, pz, pr, energy; };
  struct Slot {  // per-parcel results of the parallel pass
    std::vector<Deposit> deposits;
    std::vector<Parcel> children;
    double wallImpulse{}, wallLoss{}, swirl{};
    double outMass{}, outEnergy{}, outMomentum{};
    bool alive{true};
    long substeps{}, wallHits{}, flashes{}, vaporised{}, clamped{};
    double conductance{}, drag{}, surface{};  // implicit-coupling loads of the parcel (all its drops)
  };
  struct Pending { double mass{}, energy{}, time{}; };  // delivered, not yet released
  Flow& flow_;
  LiquidOxygen liquid_;
  Settings settings_;
  std::size_t oxygen_{};
  std::vector<Parcel> parcels_;
  std::vector<Slot> slots_;
  Ledger ledger_;
  std::vector<Pending> pending_;  // per injector
  std::vector<std::uint64_t> injectorRandom_;
  std::uint64_t nextId_{};
  // Per cell: the gas the parcels see this step (T, p, rho, u_z, u_r, c_v, then mass fractions; T and
  // the velocity are the point-implicit values once the loads are in), the implicit-coupling sums
  // (sum K, sum K v_z, sum K v_r, sum H, sum H T_s), whether a parcel is in it, and its smaller extent.
  std::vector<double> cellGas_, cellLoad_;
  std::vector<double> cellTurbulence_;  // per occupied cell with SST: k and omega
  double wall_{};
  std::vector<char> occupied_;
  std::vector<double> cellSize_;
  std::vector<double> amounts_;   // Flow::exchange input
  std::vector<double> unit_;      // pure O2 mass fractions
  std::vector<std::vector<double>> work_;  // per worker scratch
  std::vector<double> gaussNode_, gaussWeight_;  // Gauss-Legendre on [0, 1], settings_.filmPoints of them
  std::vector<double> gaussIntegral_;             // their integration matrix (filmConductances)
  // filmTransport's fits: each species' conductivity, and the binary D of O2 with each species.
  std::vector<std::array<double, 5>> conductivityFit_, oxygenDiffusionFit_;
  [[nodiscard]] int locate(double z, double r) const;
  [[nodiscard]] double wallRadius(double z, double* slope) const;
  [[nodiscard]] Gas seen(int cell) const;
  [[nodiscard]] Gas seenBy(const Parcel& parcel) const;  // the cell's, plus the parcel's eddy
  [[nodiscard]] double carrierMass(const double* y) const;
  // The stagnant film's lambda and rho D_O2 (in rates()).
  [[nodiscard]] std::array<double, 2> filmConductances(const Gas& gas, double ts, double aM, double aT, double cpf,
                                                       double* y, double* shape) const;
  [[nodiscard]] double blobDiameter(int injector) const;
  void gatherGas();
  void applyLoads(double dt);
  void moveParcel(std::size_t index, double dt, int worker);
  void release(double t0, double t1);
};
}  // namespace crucible
