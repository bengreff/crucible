#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>
#include "core/medium.hpp"

namespace crucible {
// SI throughout. z is axial, r is radial; no swirl in this first gas model.
using Conserved = std::array<double, 4>; // rho, rho*u_z, rho*u_r, total energy density
struct Primitive { double rho{}, uz{}, ur{}, p{}; };
// Thermodynamic face/cell data the flux needs besides (rho, u, p): specific internal energy, frozen
// sound speed, and the internal-energy floor below which the composition has no admissible state.
struct Thermal { double e{}, a{}, floor{}; };
// Nozzle: reservoir inlet and a prepared 1-D expansion. Chamber: closed injector plate carrying
// supplies, an initial fill of ambient gas at rest, and an ambient exit (see Supply, Igniter).
enum class Case { Nozzle, UniformDuct, ShockTube, Chamber };
// Low-Mach treatment of the HLLC dissipation (docs/evidence/LOW_MACH.md). Thornber: velocity jumps at interior
// faces scaled by min(1, local Mach) before the flux (Thornber et al., JCP 227, 2008). HllcLm: acoustic
// wave terms of the HLLC dissipation scaled by sin(pi/2 min(1, M/0.1)) (Fleischmann et al., JCP 423, 2020).
enum class LowMach { None, Thornber, HllcLm };
// A propellant stream through the injector face (z = 0) of a Chamber. The face rings whose centre
// radius lies in [innerRadius, outerRadius) carry it with a uniform mass flux, so its delivered
// mass flow is exact and its area is the mesh's ring area. Imposed: mass flow, frozen total enthalpy
// (stagnation temperature) and composition; the static pressure follows from the interior along the
// outgoing characteristic. A face that would need supersonic inflow delivers the same flow as a
// sonic (choked) stream. The valve raises the flow linearly from zero at `opens` to full at
// `opens + ramp`; valve travel is a declared input, not modelled equipment.
struct Supply {
    double innerRadius{}, outerRadius{};
    double massFlow{};           // kg/s at full opening, over the whole stream
    double totalTemperature{300};
    std::vector<double> composition;
    double opens{}, ramp{};
    [[nodiscard]] double opening(double time) const;
};
// A bounded energy deposit: `energy` joules at constant power over [start, start + duration] into
// the cells whose centroid lies in zMin <= z <= zMax, r <= rMax, in proportion to their volume.
struct Igniter {
    double zMin{}, zMax{}, rMax{}, energy{}, start{}, duration{};
};
struct Definition {
    int nz{160}, nr{24};
    double length{0.6}, inletRadius{0.035}, throatRadius{0.020}, exitRadius{0.035};
    double throatFraction{0.36};
    double totalPressure{300000}, totalTemperature{300}, backPressure{15000};
    double cfl{0.4};
    bool secondOrder{true};
    LowMach lowMach{LowMach::None};
    Case experiment{Case::Nozzle};
    Gas gas{};
    // Species of the medium (empty: the calorically perfect `gas` as one species) and the mass
    // fractions of the initial fill and the reservoir (empty: pure first species).
    std::vector<Species> species;
    std::vector<double> composition;
    // Wall contour as data: (z, r) points with z increasing from the injector face, sampled at the
    // mesh stations by linear interpolation. When given it replaces the cosine nozzle shape and
    // sets the length.
    std::vector<std::array<double, 2>> contour;
    // Chamber only: supplies, the igniter (none when energy is zero), and the ambient temperature.
    // The ambient gas (composition, backPressure, ambientTemperature) fills the chamber at rest at
    // t = 0 and is what any backflow at the exit draws in.
    std::vector<Supply> supplies;
    Igniter igniter{};
    double ambientTemperature{300};
    // Molecular transport (empty: inviscid), in the medium's species order. With transport the walls
    // (the side wall and the injector plate) are no-slip unless wallSlip, and isothermal at
    // wallTemperature, or adiabatic when it is zero. Supply rings exchange no diffusive flux; the
    // nozzle inlet and the outlet take a zero normal gradient (only tangential-gradient stress).
    TransportFits transport;
    bool wallSlip{false};
    double wallTemperature{0};
    // Subsonic outlet (Nozzle and Chamber). Zero: the boundary state takes backPressure, a pressure
    // node that reflects every acoustic wave. Positive: partially non-reflecting (Poinsot and Lele
    // 1992); the incoming characteristic relaxes the outlet pressure toward backPressure at the rate
    // K = sigma (1 - M^2) a / span, and outgoing waves leave (reflection -1 / (1 + 2 i omega / K)).
    double outletRelaxation{0};
    // Radial rings. Zero: equal heights, ring j spanning the fractions j / nr to (j + 1) / nr of the
    // local radius. Positive (b): clustered toward the wall, the fraction at s = j / nr being
    // tanh(b s) / tanh(b); the wall ring is then about 2b / sinh(2b) of the equal height (0.15 at
    // b = 2) and the axis ring b / tanh(b) of it.
    double radialStretching{0};
    [[nodiscard]] double span() const;  // axial length of the domain
    [[nodiscard]] Medium medium() const;
    [[nodiscard]] std::vector<double> massFractions() const;
    void validate() const;
};
// r is the exact volume centroid radius (integral of r dV / V) where r-weighted cell averages sit;
// radialSecondMoment is the integral of r^2 dV / V; radialPressureMeasure is the integral of dV / r.
struct Cell { double volume{}, z{}, r{}, radialPressureMeasure{}, radialSecondMoment{}; };
struct Mesh {
    explicit Mesh(const Definition& definition);
    int nz{}, nr{};
    double dz{};
    std::vector<double> radius;
    // Radial position of ring boundary j (0..nr) as a fraction of the local wall radius.
    std::vector<double> fraction;
    bool uniform{true};
    std::vector<Cell> cells;
    std::size_t index(int i, int j) const { return static_cast<std::size_t>(i)*nr+j; }
    double axialArea(int face, int j) const;
    // Radius of the middle of ring j on axial face `face`, and of radial face j at the middle of column i.
    double ringMiddle(int face, int j) const;
    double radialFaceMiddle(int i, int j) const;
    std::array<double, 2> radialAreaVector(int i, int face) const;
    // Reference radius r_f and r^2 moment of a radial face: a quantity linear in r (or in r^2)
    // integrates over the face's radial area component exactly at these values.
    double radialFaceRadius(int i, int face) const;
    double radialFaceSecondMoment(int i, int face) const;
};
// Axial forces use +z = exhaust direction for forces on the gas. Thrusts are forces on the
// device (supply/reservoir, walls, and later coils) positive against the exhaust (-z).
// Gas budget: d(P_z)/dt = inletMomentumFlux - outletMomentumFlux + wallAxialForce + bodyAxialForce.
// deviceThrust = inletMomentumFlux + wallAxialForce + bodyAxialForce - ambientAxialForce
//             = exitPlaneThrust + d(P_z)/dt  (reaction of every force the device applies to the gas).
struct Measurements {
    double time{}, dt{}, inletMassFlow{}, outletMassFlow{};
    double inletMomentumFlux{}, outletMomentumFlux{}, wallAxialForce{}, bodyAxialForce{}, ambientAxialForce{};
    double deviceThrust{}, exitPlaneThrust{}, axialMomentum{}, momentumBalanceError{};
    double mass{}, energy{}, massBalanceError{}, energyBalanceError{};
    double exitMach{}, minPressure{}, maxPressure{}, maxMach{};
    // Chamber: scheduled supply flow (the Isp denominator), igniter power and energy delivered so
    // far, and the area-averaged pressure of the cells on the injector face.
    double supplyMassFlow{}, igniterPower{}, igniterEnergy{}, injectorPressure{};
    // Heat conducted into the gas through the walls [W] (molecular transport only).
    double wallHeatFlow{};
    std::uint64_t steps{}, rejectedSteps{};
};
struct FieldSnapshot {
    Definition definition;
    std::vector<Primitive> cells;
    std::vector<double> radius, fraction;
    Measurements measurements;
    double appliedTotalPressure{};
    std::uint64_t appliedControlSequence{}, generation{};
};
Conserved conservative(Primitive w, double internalEnergy);
// Calorically perfect forms (one species, constant cp).
Conserved conservative(Primitive w, Gas gas);
Primitive primitive(const Conserved& u, Gas gas);
Thermal thermal(Primitive w, Gas gas);
// Species fluxes are not part of these: they follow the mass flux (see Flow::rhs).
Conserved hllc(Primitive left, Thermal tl, Primitive right, Thermal tr, double nz, double nr);
Conserved hllcLm(Primitive left, Thermal tl, Primitive right, Thermal tr, double nz, double nr);
Conserved hllc(Primitive left, Primitive right, double nz, double nr, Gas gas);
Conserved hllcLm(Primitive left, Primitive right, double nz, double nr, Gas gas);
// Thornber low-Mach reconstruction correction applied to a face's left/right states.
void thornberScale(Primitive& left, Primitive& right, double soundLeft, double soundRight);
void thornberScale(Primitive& left, Primitive& right, Gas gas);
double areaMach(double mach, double gamma);
double machFromArea(double areaRatio, bool supersonic, double gamma);
double chokedMassFlow(const Definition& definition);

class Flow {
public:
    explicit Flow(Definition definition);
    double step(double maxDt=std::numeric_limits<double>::infinity());
    void advanceTo(double time);
    void setTotalPressure(double pressure);
    void setUniform(Primitive state);
    void setInitialState(const std::vector<Primitive>& cells);
    // Per-cell mass fractions, cell-major (cells.size() x species).
    void setInitialState(const std::vector<Primitive>& cells, const std::vector<double>& massFractions);
    // Replace a cell's composition at fixed density and total energy (a constant-volume
    // adiabatic reaction substep leaves exactly these unchanged).
    void setMassFractions(std::size_t cell, const double* y);
    // Restore partial densities saved from partialDensities() (bulk state unchanged), e.g. to undo a split reaction substep.
    void setPartialDensities(const std::vector<double>& partial);
    // Constant volumetric force density (N/m^3, {z, r}) per cell, e.g. a Lorentz force. Its axial
    // integral enters bodyAxialForce (reaction on the equipment) and its work the energy budget.
    void setBodyForce(std::vector<std::array<double, 2>> forcePerVolume);
    [[nodiscard]] Measurements measurements() const;
    [[nodiscard]] FieldSnapshot snapshot() const;
    [[nodiscard]] const Definition& definition() const { return definition_; }
    [[nodiscard]] const Mesh& mesh() const { return mesh_; }
    [[nodiscard]] const std::vector<Conserved>& state() const { return state_; }
    [[nodiscard]] const Medium& medium() const { return medium_; }
    // Partial densities rho*Y_k, cell-major.
    [[nodiscard]] const std::vector<double>& partialDensities() const { return species_; }
    [[nodiscard]] std::vector<double> massFractions(std::size_t cell) const;
    [[nodiscard]] Primitive cellPrimitive(std::size_t cell) const;
    [[nodiscard]] double temperature(std::size_t cell) const;
    // Largest stable step for the current state (CFL bound, with the explicit diffusion bound when
    // the medium has transport).
    double stableDt();
    // The molecular-transport part of the right-hand side alone for the current state, per unit
    // volume (verification of the operator).
    void transportDerivative(std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative);
    [[nodiscard]] double time() const { return time_; }
    // Next schedule discontinuity (valve opens or finishes opening, igniter on or off) after `time`,
    // infinity if none. A step that would cross one ends on it instead.
    [[nodiscard]] double nextEvent(double time) const;
    [[nodiscard]] double totalPressure() const { return totalPressure_; }
private:
    // Exchange rates of one right-hand-side evaluation, combined with the RK weights.
    struct BoundaryRates {
        double mass{}, energy{}, inlet{}, outlet{}, inletMomentum{}, outletMomentum{}, wallAxial{}, bodyAxial{}, heat{}, wallHeat{};
        [[nodiscard]] double netMomentum() const { return inletMomentum-outletMomentum+wallAxial+bodyAxial; }
        [[nodiscard]] double grossMomentum() const;
        static BoundaryRates average(const BoundaryRates& a, const BoundaryRates& b);
    };
    Definition definition_;
    Mesh mesh_;
    Medium medium_;
    std::size_t ns_{};
    std::vector<double> inletComposition_;
    std::vector<Conserved> state_, stage_, next_, rhs_, slopesZ_;
    std::vector<Primitive> primitives_, radialLow_, radialHigh_;
    // Species: partial densities (state, stages, derivative), cell mass fractions and their
    // reconstructions, all cell-major with stride ns_. temperature_ warm-starts the energy inversion.
    std::vector<double> species_, speciesStage_, speciesNext_, speciesRhs_, fractions_, fractionSlopesZ_;
    std::vector<double> fractionsLow_, fractionsHigh_, temperature_, sound_;
    std::vector<std::array<double, 2>> bodyForce_;
    std::vector<double> pressureSource_;
    double time_{}, dt_{}, totalPressure_{}, initialMass_{}, initialEnergy_{};
    double integratedMassFlux_{}, integratedEnergyFlux_{}, initialMomentum_{};
    double integratedMomentumSource_{}, integratedMomentumGross_{}, integratedHeat_{};
    // Chamber: supply index of each injector-face ring (-1: plate), each supply's ring area, the
    // igniter's cells and their total volume.
    std::vector<int> faceSupply_;
    std::vector<double> supplyArea_;
    std::vector<std::size_t> igniterCells_;
    double igniterVolume_{}, igniterPower_{};
    std::uint64_t steps_{}, rejectedSteps_{};
    BoundaryRates lastRates_{};
    void resetAccounting();
    void radialProfiles();
    void refresh(const std::vector<Conserved>& state, const std::vector<double>& species);
    BoundaryRates rhs(const std::vector<Conserved>& state, const std::vector<double>& species,
                      std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative, double time);
    // Face flux per unit area of a supply delivering mass flux g into the interior face state.
    Conserved supplyFlux(const Supply& supply, double g, const Primitive& inside, const double* yInside) const;
    bool admissible(const Conserved& u, const double* partial) const;
    Thermal faceThermal(const Primitive& w, double* y) const;
    Primitive inlet(Primitive inside) const;
    // ambientInflow is set when a Chamber exit draws ambient gas in (subsonic backflow).
    // inside: the reconstructed face state of the last cell; cell: that cell's own state.
    Primitive outlet(Primitive inside, const double* y, const Primitive& cell, bool& ambientInflow) const;
    std::vector<double> ambient_;  // Chamber: ambient composition
    // Molecular transport (core/transport.cpp): cell viscosity, conductivity, mixture diffusion
    // coefficients and mole fractions; least-squares gradients of u_z, u_r, T and the mole fractions
    // (cell-major, 3 + ns_ fields of {d/dz, d/dr}); each cell's inverse least-squares matrix.
    std::vector<double> viscosity_, conductivity_, diffusion_, moles_, gradients_, transportWork_, faceEnthalpy_;
    std::vector<std::array<double, 3>> leastSquares_;
    void prepareTransport();
    void transportProperties();
    void transportFluxes(std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative, BoundaryRates& rates);
};
} // namespace crucible
