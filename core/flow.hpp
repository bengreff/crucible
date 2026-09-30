#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace crucible {
// SI throughout. z is axial, r is radial; no swirl in this first gas model.
using Conserved = std::array<double, 4>; // rho, rho*u_z, rho*u_r, total energy density
struct Primitive { double rho{}, uz{}, ur{}, p{}; };
struct Gas { double gamma{1.4}, specificR{287.05}; };
enum class Case { Nozzle, UniformDuct, ShockTube };
struct Definition {
    int nz{160}, nr{24};
    double length{0.6}, inletRadius{0.035}, throatRadius{0.020}, exitRadius{0.035};
    double throatFraction{0.36};
    double totalPressure{300000}, totalTemperature{300}, backPressure{15000};
    double cfl{0.4};
    bool secondOrder{true};
    Case experiment{Case::Nozzle};
    Gas gas{};
    void validate() const;
};
struct Cell { double volume{}, z{}, r{}, radialPressureMeasure{}; };
struct Mesh {
    explicit Mesh(const Definition& definition);
    int nz{}, nr{};
    double dz{};
    std::vector<double> radius;
    std::vector<Cell> cells;
    std::size_t index(int i, int j) const { return static_cast<std::size_t>(i)*nr+j; }
    double axialArea(int face, int j) const;
    std::array<double, 2> radialAreaVector(int i, int face) const;
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
Conserved conservative(Primitive w, Gas gas);
Primitive primitive(const Conserved& u, Gas gas);
Conserved hllc(Primitive left, Primitive right, double nz, double nr, Gas gas);
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
    [[nodiscard]] Measurements measurements() const;
    [[nodiscard]] FieldSnapshot snapshot() const;
    [[nodiscard]] const Definition& definition() const { return definition_; }
    [[nodiscard]] const Mesh& mesh() const { return mesh_; }
    [[nodiscard]] const std::vector<Conserved>& state() const { return state_; }
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
    std::vector<Conserved> state_, stage_, next_, rhs_, slopesZ_, slopesR_;
    std::vector<Primitive> primitives_;
    double time_{}, dt_{}, totalPressure_{}, initialMass_{}, initialEnergy_{};
    double integratedMassFlux_{}, integratedEnergyFlux_{}, initialMomentum_{};
    double integratedMomentumSource_{}, integratedMomentumGross_{};
    std::uint64_t steps_{}, rejectedSteps_{};
    BoundaryRates lastRates_{};
    void resetAccounting();
    double stableDt();
    BoundaryRates rhs(const std::vector<Conserved>& state, std::vector<Conserved>& derivative);
    Primitive inlet(Primitive inside) const;
    Primitive outlet(Primitive inside) const;
};
} // namespace crucible
