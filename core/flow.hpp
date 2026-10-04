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
enum class Case { Nozzle, UniformDuct, ShockTube };
// Low-Mach treatment of the HLLC dissipation (docs/evidence/LOW_MACH.md). Thornber: velocity jumps at interior
// faces scaled by min(1, local Mach) before the flux (Thornber et al., JCP 227, 2008). HllcLm: acoustic
// wave terms of the HLLC dissipation scaled by sin(pi/2 min(1, M/0.1)) (Fleischmann et al., JCP 423, 2020).
enum class LowMach { None, Thornber, HllcLm };
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
    std::vector<Cell> cells;
    std::size_t index(int i, int j) const { return static_cast<std::size_t>(i)*nr+j; }
    double axialArea(int face, int j) const;
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
    std::uint64_t steps{}, rejectedSteps{};
};
struct FieldSnapshot {
    Definition definition;
    std::vector<Primitive> cells;
    std::vector<double> radius;
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
    // Largest stable step for the current state (CFL bound).
    double stableDt();
    [[nodiscard]] double time() const { return time_; }
    [[nodiscard]] double totalPressure() const { return totalPressure_; }
private:
    // Exchange rates of one right-hand-side evaluation, combined with the RK weights.
    struct BoundaryRates {
        double mass{}, energy{}, inlet{}, outlet{}, inletMomentum{}, outletMomentum{}, wallAxial{}, bodyAxial{};
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
    double integratedMomentumSource_{}, integratedMomentumGross_{};
    std::uint64_t steps_{}, rejectedSteps_{};
    BoundaryRates lastRates_{};
    void resetAccounting();
    void radialProfiles();
    void refresh(const std::vector<Conserved>& state, const std::vector<double>& species);
    BoundaryRates rhs(const std::vector<Conserved>& state, const std::vector<double>& species,
                      std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative);
    bool admissible(const Conserved& u, const double* partial) const;
    Thermal faceThermal(const Primitive& w, double* y) const;
    Primitive inlet(Primitive inside) const;
    Primitive outlet(Primitive inside, const double* y) const;
};
} // namespace crucible
