#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <vector>
#include "core/medium.hpp"
#include "core/pool.hpp"

namespace crucible {
class CheckpointWriter;
class CheckpointReader;
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
// mass flow is exact and its area is the mesh's ring area (under Definition::openAreaSupplies, its
// true annulus, shared with the rings it overlaps). Imposed: mass flow, frozen total enthalpy
// (stagnation temperature) and composition; the static pressure follows from the interior along the
// outgoing characteristic. A face that would need supersonic inflow delivers the same flow as a
// sonic (choked) stream. The valve raises the flow linearly from zero at `opens` to full at
// `opens + ramp`; valve travel is a declared input, not modelled equipment.
// With turbulence the stream carries k = 3/2 (I u)^2 and omega = rho k / (viscosityRatio mu) at the
// face (I the turbulenceIntensity, both declared inputs, required then), and its total enthalpy
// includes the turbulent part: h(T) + u^2 / 2 + 5/3 k = h(T0), k plus the work of the Reynolds
// normal stress 2/3 rho k, which the momentum flux also carries.
struct Supply {
    double innerRadius{}, outerRadius{};
    double massFlow{};           // kg/s at full opening, over the whole stream
    double totalTemperature{300};
    std::vector<double> composition;
    double opens{}, ramp{};
    double turbulenceIntensity{}, viscosityRatio{};
    [[nodiscard]] double opening(double time) const;
};
// The face of a supply: per unit area, the mass, momentum and energy flux and the rho k and
// rho omega flux, and the stream's state (k and omega zero without turbulence).
struct SupplyFace {
    Conserved flux{};
    std::array<double, 2> turbulenceFlux{};
    double p{}, t{}, u{}, k{}, omega{};
    bool choked{false};
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
    // Reconstruct T in place of rho (face rho = p / (R(Y_face) T_face)): the face temperature stays
    // within the limited range of T where the composition jumps, which rho, p and Y reconstructed
    // separately do not guarantee (COST.md, lever 3).
    bool temperatureReconstruction{false};
    // At interior faces, a side whose reconstructed enthalpy h = e + k + p/rho lies outside the two
    // cells' range takes its cell's own state (first order on that side): reconstructing T and Y,
    // or rho and Y, separately can pair one cell's temperature with the other's composition
    // (COST.md, lever 3).
    bool boundedFaceEnthalpy{false};
    // At interior faces, each side's static enthalpy h = e + p/rho is the cell's plus its limited slope (as the
    // other reconstructed variables', one-sided at the plate, the exit and the wall row, r^2 on the axis row),
    // with zero slope where the cell's excess b (Flow::setBoundedExcess; h with k when none is set) is above or
    // below that of every face neighbour; the face density is p / (R T) with T solved from that h at the face
    // composition. One-sided slopes are not zero at an extremum, and h of separately reconstructed rho (or T),
    // p and Y need not be: outflow below a maximum's own b raises it (COST.md, lever 3, (x)).
    bool extremumBoundedEnthalpy{false};
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
    // A discriminator, not a model (COST.md, lever 3's backups): on interior faces every species'
    // molecular flux -(lambda / cp) dY_k/dn and the molecular energy flux -(lambda / cp) dh/dn, both
    // built as the turbulent ones under Turbulence::consistentFlux, which it requires.
    bool unityLewis{false};
    bool wallSlip{false};
    double wallTemperature{0};
    // The injector plate's wall temperature (Chamber): empty, wallTemperature; zero, adiabatic.
    std::optional<double> plateTemperature;
    // Chamber: false, each face ring carries the supply holding its centre radius (above). True, every
    // supply enters on its true annulus (docs/evidence/FACE_INLET.md): each face ring takes each
    // supply's flux in proportion to the share phi of its area the supply covers, the rest of the ring
    // (w = 1 - sum phi) being the plate; a ring counts as plate for the stencil, the wall distance and
    // the wall functions' first cell if w >= 1/2, else as its largest supply. A supply's mass flow and
    // face velocity are then its own at any mesh. A zero-flow supply is a closed port (pressure only).
    bool openAreaSupplies{false};
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
    // SST-2003 URANS (TECHNICAL_PLAN step 7; equations of the NASA Turbulence Modeling Resource
    // sst.html). rho k and rho omega are carried by the mass flux like the species; k is part of the
    // total energy (e = E / rho - |u|^2 / 2 - k); the eddy viscosity enters the stress (with the
    // Reynolds normal stress -2/3 rho k), the heat flux (cp mu_t / prandtl) and the species flux
    // (mu_t / (rho schmidt) added to every D_km); the sources are Strang-split around each step.
    // Needs transport. No-slip walls take k = 0 and omega = wallOmegaFactor 6 nu_w / (beta1 d1^2) on
    // the wall face, d1 the wall distance of the adjacent centroid and nu_w the kinematic viscosity
    // at the wall temperature (the cell's for an adiabatic wall); slip walls pass no k or omega flux
    // and are not walls for the wall distance. The initial fill carries ambientK [m^2/s^2] and
    // ambientOmega [1/s], as does ambient gas drawn in at a Chamber exit. Chamber supplies carry
    // their declared inflow turbulence (Supply). The nozzle inlet's is not yet modelled, so Nozzle
    // does not accept turbulence.
    // Wall functions (docs/evidence/WALL_FUNCTIONS.md): Nichols and Nelson's compressible, heated law
    // (crucible::wallLaw) on every no-slip wall face in place of the resolved wall, on isothermal walls.
    // The face passes tau_w against the first cell's tangential velocity and q_w, and no species, k or
    // omega; the first cell's k and omega are prescribed after every stage (eqs. 10.12, 10.17 to 10.20)
    // and its turbulence equations are not solved. wallKappa and wallB are the SST model's own log layer
    // (the resolved pipe at Re_tau 9,941, WALL_FUNCTIONS.md); printedDerivative takes mu_t from the
    // printed eq. 10.13 (to report its effect) instead of the derived one.
    struct Turbulence {
        bool enabled{false};
        double prandtl{0.9}, schmidt{0.7}, wallOmegaFactor{10};
        double ambientK{0}, ambientOmega{0};
        bool wallFunctions{false}, printedDerivative{false};
        double wallKappa{0.3697}, wallB{3.752};
        // The variable-property and intrinsic-compressibility corrections of Hasan, Elias, Menter and
        // Pecnik (J. Fluid Mech. 1019, A8, 2025), full form, as tools/sst_pipe_1d.py --correction hp
        // (docs/evidence/WALL_FUNCTIONS.md, 6 October 2026): k and omega diffuse by F1 times the
        // semi-local inner form plus (1 - F1) times the outer form, the cross-diffusion source is the
        // corrected one, mu_t carries D^ic, and the corrected k diffusion's non-divergence part enters E
        // (Measurements::correctionEnergy).
        bool propertyCorrections{false};
        // The turbulent fluxes from the two cells' differences of Y_k and h along the line of centres:
        // j_t = -(mu_t / Sc_t) dY_k/dn and the turbulent energy flux -(mu_t / Pr_t) dh/dn +
        // (1 - Sc_t / Pr_t) sum_k h_k j_t,k on interior faces, so that the species and conduction parts
        // agree by construction (COST.md, lever 3; off: mu_t / Sc_t added to D_km, cp mu_t / Pr_t grad T).
        bool consistentFlux{false};
    };
    Turbulence turbulence;
    [[nodiscard]] double span() const;  // axial length of the domain
    [[nodiscard]] Medium medium() const;
    [[nodiscard]] std::vector<double> massFractions() const;
    void validate() const;
};
// r is the exact volume centroid radius (integral of r dV / V) where r-weighted cell averages sit;
// radialSecondMoment is the integral of r^2 dV / V; radialPressureMeasure is the integral of dV / r.
struct Cell { double volume{}, z{}, r{}, radialPressureMeasure{}, radialSecondMoment{}; };
// The geometry the kernels read, held as explicit arrays (docs/evidence/AMR.md, step A): the direct
// constructor fills them from the wall radius at each axial station and the ring fractions, and
// coarsened() builds the level below a finer mesh from its arrays, so a coarse cell's sloped faces
// are the polyline of its children's.
struct Mesh {
    explicit Mesh(const Definition& definition);
    // The mesh with the columns and the rings of `finer` paired (nz and nr even): volumes, the
    // pressure measure, the r and r^2 integrals and the face area vectors summed; the radial face's
    // reference radius f R_col with R_col the ratio of the summed integrals of R^3 and R^2 over z; its
    // r^2 moment the mean of the two finer faces' (a trapezoid over the finer nodes); the side wall
    // the finest polyline.
    static Mesh coarsened(const Mesh& finer);
    int nz{}, nr{};
    double dz{};
    // The wall radius at axial station i (0..nz) and the ring boundaries' fractions of it (0..nr); on
    // a coarsened mesh, the finer mesh's at every other station and boundary.
    std::vector<double> radius;
    std::vector<double> fraction;
    bool uniform{true};
    std::vector<Cell> cells;
    // The side wall as the finest polyline: wallRadius at stations wallDz apart from z = 0, wallStride
    // of its segments per column (1 on a directly built mesh, where wallRadius is radius).
    std::vector<double> wallRadius;
    double wallDz{};
    int wallStride{1};
    std::size_t index(int i, int j) const { return static_cast<std::size_t>(i)*nr+j; }
    double axialArea(int face, int j) const { return axialAreas_[static_cast<std::size_t>(face)*nr+j]; }
    // Radius of the middle of ring j on axial face `face`, and of radial face j at the middle of column i.
    double ringMiddle(int face, int j) const { return ringMiddles_[static_cast<std::size_t>(face)*nr+j]; }
    double radialFaceMiddle(int i, int j) const { return radialMiddles_[radialIndex(i,j)]; }
    std::array<double, 2> radialAreaVector(int i, int face) const { return radialAreaVectors_[radialIndex(i,face)]; }
    // Reference radius r_f and r^2 moment of a radial face: a quantity linear in r (or in r^2)
    // integrates over the face's radial area component exactly at these values.
    double radialFaceRadius(int i, int face) const { return radialRadii_[radialIndex(i,face)]; }
    double radialFaceSecondMoment(int i, int face) const { return radialSecondMoments_[radialIndex(i,face)]; }
    // The side wall's radius at axial position z inside column i, on the finest polyline.
    double wallRadiusAt(int i, double z) const;
    // Per column the integrals of R^2 and of R^3 over z, R(z) the wall radius: a cell between ring
    // fractions a and b has the volume pi (b^2 - a^2) times the first and the centroid radius
    // 2 (b^3 - a^3) / (3 (b^2 - a^2)) times R_col, the integral of R^3 over that of R^2.
    const std::vector<double>& columnSquares() const { return columnSquares_; }
    const std::vector<double>& columnCubes() const { return columnCubes_; }
private:
    Mesh() = default;
    std::size_t radialIndex(int i, int j) const { return static_cast<std::size_t>(i)*(nr+1)+j; }
    std::vector<double> axialAreas_, ringMiddles_;
    std::vector<std::array<double, 2>> radialAreaVectors_;
    std::vector<double> radialMiddles_, radialRadii_, radialSecondMoments_;
    std::vector<double> columnSquares_, columnCubes_;
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
    // The same three balance errors from the compensated exchange sums (AMR.md, the tagged run's mass error).
    double compensatedMassBalanceError{}, compensatedEnergyBalanceError{}, compensatedMomentumBalanceError{};
    double exitMach{}, minPressure{}, maxPressure{}, maxMach{};
    // Chamber: scheduled supply flow (the Isp denominator), igniter power and energy delivered so
    // far, and the area-averaged pressure of the cells on the injector face.
    double supplyMassFlow{}, igniterPower{}, igniterEnergy{}, injectorPressure{};
    // Heat conducted into the gas through the walls [W] (molecular transport only), and the part of it through
    // the injector face (i = 0: its walls, or the plate's share under open-area supplies).
    double wallHeatFlow{}, plateHeatFlow{};
    // Turbulence: rho k V added so far where a stage left rho k below zero [J] (Flow::step), and by
    // the wall functions' prescription of the first cells' k (at fixed total energy).
    double clippedTurbulentEnergy{}, prescribedTurbulentEnergy{};
    // Turbulence::propertyCorrections: the energy the corrected k diffusion's non-divergence part has
    // added so far [J] (in the energy balance), and the F1-weighted fraction of cells whose S_n
    // denominator is at its floor psi / 10.
    double correctionEnergy{}, correctionFloored{};
    // Flow::exchange: mass, total energy and axial momentum added so far [kg, J, kg m/s] (in the balances),
    // and the momentum balance error's denominator [kg m/s].
    double exchangedMass{}, exchangedEnergy{}, exchangedMomentum{}, momentumBalanceScale{};
    std::uint64_t steps{}, rejectedSteps{};
    // The booked boundary rates above are the last step's. These are their time integrals since t = 0 [kg, N s, J]:
    // a history that samples every few steps reports their change over its period (periodMean), not the last step's.
    struct Integrals {
        double inletMass{}, outletMass{}, inletMomentum{}, outletMomentum{}, wallAxial{}, bodyAxial{}, wallHeat{}, plateHeat{};
    } integrals;
};
// The boundary rates (mass flows, momentum fluxes, axial forces, heat flows) and the thrusts derived from them,
// averaged over (previous.time, now.time]; the rest of now as it is. The last step of a period is shortened to land
// on the sample, and with a split after the flow stages (Table A's projection) its rates differ from the period's:
// on C1 at CFL 0.4 by up to 0.4% of the vacuum thrust (docs/evidence/IMPLICIT_PROTOTYPE.md, criterion 5's cause).
// With now.time <= previous.time, now unchanged.
Measurements periodMean(const Measurements& now, const Measurements& previous);
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
    // A flow on a given mesh with definition's nz x nr cells: a level of the AMR hierarchy
    // (core/amr.hpp), the finest built directly and the others by Mesh::coarsened.
    Flow(Definition definition, Mesh mesh);
    // Threads for the loops of the step (1 by default). The result is the same bit for bit on any
    // number: every cell and face is computed alone, and the sums over them keep the serial order.
    void setThreads(int threads);
    [[nodiscard]] int threads() const { return pool_->threads(); }
    // The step's threads, for per-cell work between steps (ReactingFlow). Not reentrant.
    [[nodiscard]] Pool& pool() { return *pool_; }
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
    // Turbulence: specific k [m^2/s^2] and omega [1/s] per cell, cell-major pairs. The temperature is
    // kept: the total energy takes the change of rho k. Resets the budgets like setInitialState.
    void setTurbulence(const std::vector<double>& kOmega);
    // Constant volumetric force density (N/m^3, {z, r}) per cell, e.g. a Lorentz force. Its axial
    // integral enters bodyAxialForce (reaction on the equipment) and its work the energy budget.
    void setBodyForce(std::vector<std::array<double, 2>> forcePerVolume);
    // Constant volumetric heating (W/m^3) per cell, e.g. to hold a heated pipe fully developed. It
    // enters the energy budget as an external exchange.
    void setHeating(std::vector<double> powerPerVolume);
    // Add to each cell, from outside the gas (a spray, core/spray.hpp), per cell in stride ns + 3:
    // each species' mass [kg], axial and radial momentum [kg m/s] and total energy [J]. Booked in the
    // budgets as an external exchange (Measurements::exchanged*); throws if a cell is left non-admissible.
    void exchange(const std::vector<double>& amounts);
    [[nodiscard]] Measurements measurements() const;
    [[nodiscard]] FieldSnapshot snapshot() const;
    [[nodiscard]] const Definition& definition() const { return definition_; }
    [[nodiscard]] const Mesh& mesh() const { return mesh_; }
    [[nodiscard]] const std::vector<Conserved>& state() const { return state_; }
    [[nodiscard]] const Medium& medium() const { return medium_; }
    // Partial densities rho*Y_k, cell-major.
    [[nodiscard]] const std::vector<double>& partialDensities() const { return species_; }
    // Turbulence: rho k and rho omega, cell-major pairs (empty without it); the exact wall distance
    // of each centroid (infinite without no-slip walls); the eddy viscosity of the last stableDt()
    // or mixingInputs().
    [[nodiscard]] const std::vector<double>& turbulence() const { return turbulence_; }
    [[nodiscard]] const std::vector<double>& wallDistances() const { return wallDistance_; }
    [[nodiscard]] const std::vector<double>& eddyViscosities() const { return eddy_; }
    // AMR step B2 (docs/evidence/AMR.md, Step B2: setup): the cells whose state at the end of a step is
    // read (N_1), which an after-step hook such as a reaction step must update; null: every cell.
    [[nodiscard]] const std::vector<std::size_t>* stepCells() const { return region(1); }
    // Wall functions: per no-slip wall face (the side wall by column, then the plate rings), its cell,
    // whether it is on the plate, and the first cell's y+ = rho_w u_tau d / mu_w from the last
    // transport evaluation, the second stage of the last step, with its q_w and T_w (the recovery
    // temperature on an adiabatic wall; empty without wall functions).
    struct WallYPlus { std::size_t cell; bool plate; double yPlus, heat, temperature; };
    [[nodiscard]] std::vector<WallYPlus> wallYPlus() const;
    // The PaSR closure's inputs per cell at the current state (TECHNICAL_PLAN step 7, "PaSR"; needs
    // turbulence): the mixing time tau_mix = cmix sqrt(nu_eff / epsilon), nu_eff = (mu + mu_t) / rho,
    // epsilon = beta* k omega (infinite where k = 0), and the segregation s, the largest over the
    // given species (indices into the medium) with 0 < X_i < 1 of
    // min(1, (mu_t / (rho schmidt)) |grad X_i|^2 / (beta* omega X_i (1 - X_i))), 0 if there is none.
    void mixingInputs(double cmix, const std::vector<std::size_t>& species, std::vector<double>& time,
                      std::vector<double>& segregation);
    [[nodiscard]] std::vector<double> massFractions(std::size_t cell) const;
    [[nodiscard]] Primitive cellPrimitive(std::size_t cell) const;
    [[nodiscard]] double temperature(std::size_t cell) const;
    // Largest stable step for the current state (CFL bound, with the explicit diffusion bound when
    // the medium has transport).
    double stableDt();
    // The cell that set the last stableDt() (the lowest index among equal steps, so independent of the
    // thread count) and its convective and diffusive rates per unit volume: the step is cfl / (sum).
    struct StepLimit { std::size_t cell = 0; double convective = 0, diffusive = 0; };
    [[nodiscard]] const StepLimit& stepLimit() const { return stepLimit_; }
    // The molecular-transport part of the right-hand side alone for the current state, per unit
    // volume (verification of the operator).
    void transportDerivative(std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative);
    void transportDerivative(std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative,
                             std::vector<double>& turbulenceDerivative);
    [[nodiscard]] double time() const { return time_; }
    // A Chamber's injector-face convective flux per unit area, per ring, at the current state and time
    // (verification; docs/evidence/FACE_INLET.md, F4).
    [[nodiscard]] std::vector<Conserved> injectorFluxes();
    // The number of owned wall-function cells whose conserved state, mass fractions, rho k and rho omega
    // are not bitwise what their last prescription left (verification; AMR.md, step D, D2 as amended);
    // 0 without wall functions.
    [[nodiscard]] std::size_t prescriptionMismatches() const;
    // The igniter's mean power over the last step, before this level's share of it (AMR.md, D1's cause).
    [[nodiscard]] double igniterPower() const { return igniterPower_; }
    // Each supply's area on the injector face (its rings' areas, or under openAreaSupplies its shares').
    [[nodiscard]] const std::vector<double>& supplyAreas() const { return supplyArea_; }
    // Next schedule discontinuity (valve opens or finishes opening, igniter on or off) after `time`,
    // infinity if none. A step that would cross one ends on it instead.
    [[nodiscard]] double nextEvent(double time) const;
    [[nodiscard]] double totalPressure() const { return totalPressure_; }
    // Face of a supply delivering mass flux g [kg/(m^2 s)] against the interior face state inside
    // (composition yInside, mass fractions); g = 0 is a closed valve.
    [[nodiscard]] SupplyFace supplyFace(const Supply& supply, double g, const Primitive& inside, const double* yInside) const;
    // A read-only budget probe (docs/evidence/COST.md, lever 3; single level). For each listed cell, the
    // integrals over the accepted steps, since setBudgetProbe, of what entered it through each of its
    // faces by mechanism, and of its sources. Faces: 0 the axial face at its low z (the injector plate or
    // a supply ring in the first column of a chamber), 1 the one at its high z, 2 the radial face below
    // it (zero at the axis), 3 the one above it. Per face, as gains of the cell in kg, J and kg of the
    // tracked combination sum_k weights[k] rho Y_k (e.g. the hydrogen element): the convective mass,
    // total energy and tracked fluxes; the viscous work, the molecular and turbulent conduction, the
    // diffusion of k (part of E), the species-enthalpy diffusion sum_k h_k j_k and the tracked diffusion
    // sum_k w_k j_k, each split into its molecular part and its turbulent part (each with its share of
    // the correction flux); and, as the cells used them, the viscous energy flux (the parts' sum to
    // rounding) and the tracked diffusion. A wall-function face's heat flux is booked as molecular
    // conduction, an open end's viscous energy flux as work. After the four faces: the cell's sources.
    // The parts are recomputed from the stage's own gradients and properties and do not touch the step.
    enum FaceTerm { kInviscidMass, kInviscidEnergy, kInviscidTracked, kViscousWork, kConduction, kTurbulentConduction,
                    kTurbulentEnergyDiffusion, kSpeciesEnthalpy, kTurbulentSpeciesEnthalpy, kTrackedDiffusion,
                    kTurbulentTrackedDiffusion, kViscousEnergy, kTrackedTransport, kFaceTerms };
    enum SourceTerm { kIgniter, kHeating, kBodyWork, kCorrection, kExchangedMass, kExchangedEnergy, kExchangedTracked,
                      kSourceTerms };
    static constexpr std::size_t kProbeStride=4*kFaceTerms+kSourceTerms;
    struct BudgetProbe {
        std::vector<std::size_t> cells;
        std::vector<double> weights;  // per species
        std::vector<double> sums;     // per cell in stride kProbeStride: face f's term t at f kFaceTerms + t, then the sources
    };
    void setBudgetProbe(std::vector<std::size_t> cells, std::vector<double> weights);
    // Under boundedFaceEnthalpy, the bound judges the excess b = h + |u|^2/2 - sum_k c_k Y_k (h without k) in
    // place of h when c (per species) is given, e.g. c_k = h_ox + (h_fuel - h_ox) times the fuel element's share
    // of species k: the total enthalpy above the propellants' mixing line, which h and Y bounded separately
    // do not bound (COST.md, lever 3). Empty: h. Definition::extremumBoundedEnthalpy's extremum test uses the same measure.
    void setBoundedExcess(std::vector<double> c) { boundedExcess_ = std::move(c); }
    [[nodiscard]] const BudgetProbe& budgetProbe() const { return probe_; }
    // A restart checkpoint (core/checkpoint.hpp, docs/evidence/RESTART.md): the state, the temperature guess,
    // the primitives and mass fractions the last refresh left, the active cells, the budgets' sums and what
    // the measurements read. Read into a Flow built from the same Definition (on AMR, the same level), the
    // continued run is the uninterrupted one bit for bit.
    void writeCheckpoint(CheckpointWriter& out) const;
    void readCheckpoint(CheckpointReader& in);
    // The implicit prototype (docs/evidence/IMPLICIT_PROTOTYPE.md; one level, inviscid, no turbulence). One physical
    // step of dt by BDF2, solved by dual time: LU-SGS sweeps on the Rusanov split of the convective Jacobian, the
    // scalar diagonal plus V / (a dt), each cell's pseudo-time step pseudoCfl times its convective limit. The history
    // is the state at the start of the last implicit step; the first step, or one after the state was stepped by
    // anything else, is backward Euler, and a step of another size than the last takes the variable-step coefficients.
    // The inner loop ends when the defect d = U - U_known - a dt R(U), as an RMS over the volume of each bulk component
    // over its scale (rho, rho a, rho a, E of U^n), is at most tolerance; the step is then U_known + a dt R of the last
    // iterate, so the budgets, booked with the same coefficients, close to rounding. An iterate's cell update is damped
    // until the cell is admissible; a final state that is not admissible, or a failed evaluation, halves the step (at
    // most 14 times). Returns the step taken.
    // extrapolate: the first iterate U^n + w (U^n - U^(n-1)) in each cell where it is admissible (U^n elsewhere and
    // without a history); trace: each inner iteration's defect to stdout. project, when set, is called before each
    // residual with the iterate swapped into the state (a closure's projection of the composition, such as Table A's,
    // inside the step instead of split after it).
    struct ImplicitControl {
        double tolerance{1e-6}; int maxIterations{100}; double pseudoCfl{1e6}; bool extrapolate{false}, trace{false};
        std::function<void()> project;
    };
    struct ImplicitStats {
        std::uint64_t steps{}, iterations{}, rejected{}, unconverged{}, rejectedIterations{};  // the last: spent on rejected attempts
        int lastIterations{}, mostIterations{};
        double lastDefect{}, lastFirstDefect{};
    };
    double implicitStep(double dt, const ImplicitControl& control);
    [[nodiscard]] const ImplicitStats& implicitStats() const { return implicitStats_; }
private:
    // The implicit prototype's history (the state at the start of the last implicit step and that step), its time and
    // step count at the end of that step (a mismatch means something else stepped), the budgets' increments of that
    // step, and the sweeps' scratch.
    std::vector<Conserved> implicitPrevious_;
    std::vector<double> implicitPreviousSpecies_;
    double implicitPreviousDt_{}, implicitEndTime_{-1};
    std::uint64_t implicitEndSteps_{};
    std::array<double, 6> implicitIncrement_{};
    std::vector<double> implicitKnown_, implicitDelta_, implicitDiagonal_, implicitAxialRate_, implicitRadialRate_;
    ImplicitStats implicitStats_;
    void implicitSweeps(double diagonalTime, double pseudoCfl);
    friend class Amr;
    template<class Self, class Io> static void checkpoint(Self& self, Io& io);
    // The budget probe: its sums, the current right-hand side's terms and the attempt's (dt/2 per stage).
    BudgetProbe probe_;
    std::vector<double> boundedExcess_;
    std::vector<double> probeStage_, probeAttempt_;
    void probeConvective();
    // Exchange rates of one right-hand-side evaluation, combined with the RK weights.
    struct BoundaryRates {
        double mass{}, energy{}, inlet{}, outlet{}, inletMomentum{}, outletMomentum{}, wallAxial{}, bodyAxial{}, heat{}, wallHeat{};
        double correction{};  // the integral of Phi_k dV (Turbulence::propertyCorrections)
        double plateHeat{};   // wallHeat's part through the injector face (i = 0)
        [[nodiscard]] double netMomentum() const { return inletMomentum-outletMomentum+wallAxial+bodyAxial; }
        [[nodiscard]] double grossMomentum() const;
        static BoundaryRates average(const BoundaryRates& a, const BoundaryRates& b);
        void accumulate(const BoundaryRates& rates, double dt);  // this += dt rates, field by field
    };
    Definition definition_;
    Mesh mesh_;
    Medium medium_;
    // Species count, turbulence fields (2 with turbulence, else 0) and the stride of the face
    // fractions (ns_ mass fractions, then specific k and omega).
    std::size_t ns_{}, nt_{}, nw_{};
    std::vector<double> inletComposition_;
    std::vector<Conserved> state_, stage_, next_, rhs_, slopesZ_;
    std::vector<Primitive> primitives_, radialLow_, radialHigh_;
    // Species: partial densities (state, stages, derivative), cell-major with stride ns_. Turbulence:
    // rho k and rho omega likewise with stride nt_. Cell mass fractions with specific k and omega and
    // their reconstructions, stride nw_. temperature_ warm-starts the energy inversion.
    std::vector<double> species_, speciesStage_, speciesNext_, speciesRhs_, fractions_, fractionSlopesZ_;
    std::vector<double> turbulence_, turbulenceStage_, turbulenceNext_, turbulenceRhs_, turbulenceStart_;
    std::vector<double> fractionsLow_, fractionsHigh_, temperature_, sound_;
    std::vector<std::array<double, 2>> bodyForce_;
    std::vector<double> heating_;
    std::vector<double> pressureSource_;
    std::unique_ptr<Pool> pool_;
    // Pool blocks: cells, and faces. Face blocks are small so the dozen costly supply faces of the
    // first row (supplyFace) spread over the workers.
    static constexpr std::size_t kCells=32,kFaces=4;
    // Face fluxes of one right-hand side: bulk per unit area, and the transported mass fractions
    // with k and omega times the area (stride nw_; at a supply face, its k and omega fluxes per unit
    // area). Axial face (i, j) at i nr + j for i = 0..nz; radial face j of column i at i (nr + 1) + j,
    // j = 1..nr (j = 0 is the axis); the face areas likewise. The cells sum them afterwards with the
    // arithmetic and in the order of the serial face loops (each area product inside the sum, where
    // the compiler may fuse it), so the result does not depend on the thread count.
    std::vector<Conserved> axialFlux_, radialFlux_;
    std::vector<double> axialTransported_, radialTransported_, axialArea_, radialArea_;
    // Per worker: face mass fractions (left, right), the stage's admissibility and smallest stable step,
    // and the transport scratch (Medium::transport's work, the wall's diffusion coefficients, the face
    // mass fractions and species enthalpies).
    struct TransportScratch { std::vector<double> work, diffusion, face, enthalpy, gradient, turbulent; };
    std::vector<std::vector<double>> faceFractions_;
    std::vector<TransportScratch> transportScratch_;
    std::vector<char> workerOk_;
    std::vector<double> workerDt_;
    std::vector<StepLimit> workerLimit_;
    StepLimit stepLimit_;
    // Per cell: inside the igniter, and its volume inside it (its own volume on a single level); rho k V
    // added by the last stage's clip (Flow::step).
    std::vector<char> igniterCell_;
    std::vector<double> igniterCellVolume_, clippedCell_;
    double time_{}, dt_{}, totalPressure_{}, initialMass_{}, initialEnergy_{};
    double integratedMassFlux_{}, integratedEnergyFlux_{}, initialMomentum_{};
    double integratedMomentumSource_{}, integratedMomentumGross_{}, integratedHeat_{}, clippedTurbulentEnergy_{};
    double prescribedTurbulentEnergy_{};
    BoundaryRates integratedRates_{};  // the booked rates' time integrals since t = 0 (Measurements::integrals)
    // The delivered igniter energy is integratedHeat_ plus heatCompensation_, a Neumaier compensated sum
    // (docs/evidence/AMR.md, Before step C).
    double heatCompensation_{};
    [[nodiscard]] double deliveredHeat() const { return integratedHeat_+heatCompensation_; }
    // The exchange sums carry the same compensation beside their plain values, which stay as they were
    // (AMR.md, the tagged run's mass error): a level's sum can run to the whole throughput.
    double massCompensation_{}, energyCompensation_{}, correctionCompensation_{}, momentumCompensation_{};
    [[nodiscard]] double compensatedMassFlux() const { return integratedMassFlux_+massCompensation_; }
    [[nodiscard]] double compensatedEnergyFlux() const { return integratedEnergyFlux_+energyCompensation_; }
    [[nodiscard]] double compensatedCorrection() const { return integratedCorrection_+correctionCompensation_; }
    [[nodiscard]] double compensatedMomentumSource() const { return integratedMomentumSource_+momentumCompensation_; }
    // Chamber: supply index of each injector-face ring (-1: plate), each supply's ring area, the
    // igniter's cells and their total volume.
    std::vector<int> faceSupply_;
    std::vector<double> supplyArea_;
    // Under openAreaSupplies: each face ring's supply shares (ring, supply, phi) in ring order from
    // shareBegin_[j], its wall share w (plateShare_), and from the last right-hand side the plate
    // share's axial momentum flux per unit area (plateMomentum_).
    struct FaceShare { int ring, supply; double phi; };
    std::vector<FaceShare> faceShares_;
    std::vector<std::size_t> shareBegin_;
    std::vector<double> plateShare_, plateMomentum_;
    void openAreaShares();
    [[nodiscard]] double plateTemperature() const { return definition_.plateTemperature.value_or(definition_.wallTemperature); }
    std::vector<std::size_t> igniterCells_;
    double igniterVolume_{}, igniterPower_{};
    std::uint64_t steps_{}, rejectedSteps_{};
    BoundaryRates lastRates_{};
    // The step in pieces, so the AMR driver can fill a level's non-active cells between the stages:
    // beginStep saves rho k, rho omega and the frozen source coefficients; firstStage (igniter power,
    // the first half source, u^(1)) and secondStage (the wall functions' prescription of u^(1), u^(n+1),
    // the second half source) return false on a non-admissible state; accept swaps in u^(n+1) and books
    // the budgets. step() is these in this order, with the same arithmetic.
    struct Attempt { double clipped{}, prescribed{}; BoundaryRates first{}, second{}; };
    void beginStep();
    bool firstStage(double dt, Attempt& attempt);
    bool secondStage(double dt, Attempt& attempt);
    void accept(double dt, double newTime, const Attempt& attempt);
    // A stage's clip of negative rho k in cell q (kept per cell), the stage's clips summed in cell
    // order, and the wall functions' prescription of a stage with its rho k V summed (false if not admissible).
    void clipTurbulence(std::vector<double>& turbulence, std::size_t q);
    void sumClipped(Attempt& attempt) const;
    bool prescribe(std::vector<Conserved>& state, std::vector<double>& species, std::vector<double>& turbulence, Attempt& attempt);
    // AMR: per cell whether this level updates it (empty: all). A cell that is not active is carried
    // through the stages unchanged and left out of the stable step and the admissibility checks.
    std::vector<unsigned char> active_;
    [[nodiscard]] bool isActive(std::size_t q) const { return active_.empty() || active_[q]; }
    // AMR: per cell whether this level books the boundary exchange beside it and its sources (igniter,
    // heating, body force) in the budgets (empty: all): active and not covered by the level above.
    std::vector<unsigned char> owned_;
    [[nodiscard]] bool isOwned(std::size_t q) const { return owned_.empty() || owned_[q]; }
    // AMR: with recordFaces_, each attempt sums at every face the flux that moved the cells,
    // dt (F(u^n) + F(u^(1))) / 2 times the area, convective and viscous, per face in stride 4 + nw_ (bulk,
    // then the transported mass fractions, k and omega), the axial faces first and then the radial ones.
    bool recordFaces_{false};
    std::vector<double> faceIntegral_;
    void recordFaceIntegrals(double weight);
    // AMR step B2 (docs/evidence/AMR.md, Step B2: setup): the kernels run over lists of cells and faces
    // instead of the whole level. within_[k] lists the cells within Chebyshev distance k (index space) of
    // an active cell, N_k, k = 0..g for the ghost width g, in ascending order, and distance_ holds each
    // cell's distance (g + 1 beyond N_g). filled_ lists the non-active cells the driver fills (all four
    // children of each parent of N_g \ A, so they contain N_g \ A), carried through the stages as B1
    // carries its non-active cells. faces_ lists the faces with an active cell on a side in
    // faceIntegral_'s numbering, the axialFaces_ axial ones first; columns_ the rows of N_1 in each column
    // that has any, as [lo, hi) intervals. With no regions (a single level, the base, or B1) every loop
    // runs over the whole level.
    struct Column { int i; std::vector<std::array<int, 2>> rows; };
    std::vector<std::vector<std::size_t>> within_;
    std::vector<unsigned char> distance_;
    std::vector<std::size_t> filled_, faces_;
    std::size_t axialFaces_{};
    std::vector<Column> columns_;
    void setRegions(int width);
    void clearRegions();
    // The cells of N_k (null: every cell, with no regions or for k = kWhole), and whether cell q is among them.
    static constexpr int kWhole=-1;
    [[nodiscard]] const std::vector<std::size_t>* region(int k) const {
        return within_.empty() || k<0?nullptr:&within_[static_cast<std::size_t>(k)];
    }
    [[nodiscard]] bool inRegion(std::size_t q,int k) const { return distance_.empty() || k<0 || distance_[q]<=k; }
    // The ghost width g of the regions (0 with none, where the depths are not used).
    [[nodiscard]] int ghostWidth() const { return within_.empty()?0:static_cast<int>(within_.size())-1; }
    // Runs f(q, worker) on the pool for each cell of a list (null: every cell).
    template<class F> void eachCell(const std::vector<std::size_t>* list,F f) const;
    // Applies f(q) on the pool to each cell of a list (null: every cell); true if it returned true for all.
    template<class F> bool allOf(const std::vector<std::size_t>* list,F f);
    void resetAccounting();
    // Column i's limited axial slopes and radial face states in rows lo..hi - 1.
    void axialSlopes(int i,int lo,int hi);
    void radialProfiles(int i,int lo,int hi);
    // One axial or radial face's area-weighted flux into axialFlux_ / radialFlux_ (y: the worker's scratch).
    void axialFace(int i, int j, double time, std::vector<double>& yl, std::vector<double>& yr);
    void radialFace(int i, int j, std::vector<double>& yl, std::vector<double>& yr);
    // Cell q's sum of its face fluxes, igniter and body force (before transport, the pressure source and the volume).
    void gather(std::size_t q, std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative,
                std::vector<double>& turbulenceDerivative) const;
    // Over the cells of N_depth (region).
    void refresh(const std::vector<Conserved>& state, const std::vector<double>& species, const std::vector<double>& turbulence, int depth);
    BoundaryRates rhs(const std::vector<Conserved>& state, const std::vector<double>& species, const std::vector<double>& turbulence,
                      std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative,
                      std::vector<double>& turbulenceDerivative, double time);
    // turbulence: the cell's rho k and rho omega (unused without turbulence).
    bool admissible(const Conserved& u, const double* partial, const double* turbulence) const;
    // Specific k of a cell of a state (zero without turbulence).
    [[nodiscard]] double turbulentEnergy(const Conserved& u, const std::vector<double>& turbulence, std::size_t q) const {
        return nt_?turbulence[q*nt_]/u[0]:0.0;
    }
    Thermal faceThermal(const Primitive& w, double* y) const;
    // A cell's reconstructed variables: (rho, u_z, u_r, p), or (T, u_z, u_r, p) under temperatureReconstruction;
    // toFace turns a reconstructed face state back into (rho, u_z, u_r, p), normalising y.
    [[nodiscard]] Conserved reconstructionValues(std::size_t q) const;
    void toFace(Primitive& w, double* y) const;
    // boundedFaceEnthalpy: the bounded measure of a state (w, e with k, y), cell q's, and the bound applied to
    // one side's face state (w, t, y) between cells q (its own) and other.
    [[nodiscard]] double boundMeasure(const Primitive& w, double e, const double* y) const;
    [[nodiscard]] double cellEnthalpy(std::size_t q) const;
    void boundEnthalpy(Primitive& w, Thermal& t, double* y, std::size_t q, std::size_t other) const;
    // extremumBoundedEnthalpy: whether cell q's excess is above or below every face neighbour's, and a face
    // state's density from static enthalpy h at its composition y (normalised in place; unchanged if h has no T).
    [[nodiscard]] bool excessExtremum(std::size_t q) const;
    void enthalpyFace(Primitive& w, double h, double* y) const;
    Primitive inlet(Primitive inside) const;
    // ambientInflow is set when a Chamber exit draws ambient gas in (subsonic backflow).
    // inside: the reconstructed face state of the last cell; cell: that cell's own state.
    Primitive outlet(Primitive inside, const double* y, const Primitive& cell, bool& ambientInflow) const;
    std::vector<double> ambient_;  // Chamber: ambient composition (then ambient k and omega)
    // Molecular transport (core/transport.cpp): cell viscosity, conductivity, mixture diffusion
    // coefficients and mole fractions; least-squares gradients of u_z, u_r, T and the mole fractions
    // (cell-major, 3 + ns_ fields of {d/dz, d/dr}); each cell's inverse least-squares matrix. The
    // transport face fluxes, indexed like axialFlux_ and radialFlux_: bulk per unit area, then the
    // species and k, omega fluxes per unit area (stride ns_ + nt_), summed per cell like the convective ones.
    // fieldValues_ holds each cell's gradient fields (transportValue, cell-major), gathered by transportGradients.
    std::vector<double> viscosity_, conductivity_, diffusion_, moles_, gradients_, fieldValues_;
    std::vector<std::array<double, 3>> leastSquares_;
    std::vector<Conserved> axialViscous_, radialViscous_;
    std::vector<double> axialViscousTransported_, radialViscousTransported_;
    // Turbulence: eddy viscosity mu_t, turbulent conductivity cp mu_t / Pr_t and the turbulent parts
    // sigma_k mu_t, sigma_omega mu_t of the k and omega diffusivities (stride 2) per cell; the wall
    // distance, and the wall-face omega of the cells next to a no-slip wall (wallCells_, each once). sources_ holds what
    // the source split freezes per half step: S^2 and S of the mean flow, the cross-diffusion
    // product grad k . grad omega, the molecular nu and the inverse wall distance.
    // With the property corrections: crossSource is the corrected cross-diffusion's gradient product
    // grad(rho k) . grad(sqrt(rho) omega) / rho^(3/2) (crossGradient, the conventional one, stays in F1),
    // and sound the speed of sound for D^ic's M_t.
    struct SourceCoefficients { double strain2{}, strain{}, crossGradient{}, nu{}, inverseDistance{}, crossSource{}, sound{1}; };
    std::vector<double> eddy_, eddyConductivity_, eddyDiffusion_, wallDistance_, wallOmega_;
    std::vector<double> staticEnthalpy_, heatCapacity_;  // per cell, h and cp at its T and Y (consistentFlux)
    // extremumBoundedEnthalpy: per cell, h = e + p/rho and the excess (boundMeasure); h's axial slope and its
    // radial face values (low, high).
    std::vector<double> cellStaticEnthalpy_, cellExcess_, enthalpySlopesZ_, enthalpyLow_, enthalpyHigh_;
    // Property corrections: the unit vector away from the nearest wall; per cell the fields psi =
    // sqrt(rho) / mu, rho k, sqrt(rho) omega and mu omega (stride 4, gradient fields 3 + ns_ + nt_ on);
    // rho_w and mu_w at the no-slip wall cells; S_n / mu, F1 and whether S_n is floored; the inner and
    // outer k and omega fluxes per unit area of each face (stride 4, indexed like axialFlux_ and radialFlux_).
    std::vector<std::array<double, 2>> wallDirection_;
    std::vector<double> correctionFields_, wallDensity_, wallViscosity_, snOverMu_, blendF1_;
    std::vector<char> snFloored_;
    std::vector<double> axialCorrection_, radialCorrection_, correctionCell_;
    double integratedCorrection_{};
    double exchangedMass_{}, exchangedEnergy_{}, exchangedMomentum_{};
    [[nodiscard]] std::size_t gradientFields() const { return 3+ns_+nt_+(corrections_?4:0); }
    bool corrections_{false};
    std::vector<std::size_t> wallCells_;
    std::vector<SourceCoefficients> sources_, sourcesStart_;
    // Wall functions. Each no-slip wall face (the side wall by column, then the plate rings): its
    // cell, whether it is on the plate, the unit normal out of the gas and the centroid's normal
    // distance to it. Its solution at the current state: the wall's traction on the gas, q_w (gas to
    // wall), the first cell's mu_t and omega, and its y+. Per entry of wallCells_ the face that sets that
    // cell's k and omega (the nearer, the side wall on a tie); per ring the plate face (-1: none); per
    // cell whether its turbulence is prescribed, and the rho k V the last prescription added.
    struct WallFace { std::size_t cell; bool plate; double nz, nr, distance; };
    struct WallSolution { double tz{}, tr{}, heat{}, eddy{}, omega{}, yPlus{}, temperature{}; };
    std::vector<WallFace> wallFaces_;
    std::vector<WallSolution> wallSolutions_;
    std::vector<std::size_t> wallCellFace_;
    std::vector<int> plateWallFace_;
    std::vector<char> prescribedCell_;
    std::vector<double> prescribedChange_;
    // What each wall-function cell's last prescription left, per wallCells_ entry: its state, mass
    // fractions and turbulence (D2).
    std::vector<double> prescribedValues_;
    // The law at one wall face for the cell state (w, t1, mass fractions y, molecular viscosity mu1).
    WallSolution wallSolve(const WallFace& face, const Primitive& w, double t1, const double* y, double mu1,
                           TransportScratch& scratch) const;
    // Sets rho k and rho omega of the wall-function cells of a stage from the law at that stage's
    // state (total energy fixed) and records each change in prescribedChange_; false if a cell is
    // then not admissible.
    bool prescribeWallTurbulence(const std::vector<Conserved>& state, const std::vector<double>& species,
                                 std::vector<double>& turbulence);
    // The same for wallCells_[n] alone, with one worker's scratch and mass-fraction buffer y.
    bool prescribeCell(std::size_t n, const std::vector<Conserved>& state, const std::vector<double>& species,
                       std::vector<double>& turbulence, TransportScratch& scratch, double* y);
    // Under AMR, after a reflux: the wall-function cells marked in touched (per cell) prescribed again
    // at the current state, the owned cells' change booked in prescribedTurbulentEnergy_ (AMR.md,
    // step D); false if one is then not admissible.
    bool represcribe(const std::vector<char>& touched);
    void prepareTransport();
    // Over N_depth, with the wall-face omega of the wall cells among them.
    void transportProperties(int depth);
    // Gradient field f of a cell (u_z, u_r, T, the mole fractions, then k and omega) and its value
    // on a wall face with unit normal (nz, nr).
    [[nodiscard]] double transportValue(std::size_t q, std::size_t f) const;
    [[nodiscard]] double wallValue(std::size_t q, std::size_t f, double nz, double nr, bool plate) const;
    // The gradients over N_depth, from the field values gathered over N_(depth + 1).
    void transportGradients(int depth);
    // Eddy viscosity and diffusivities from the current gradients over N_depth; with sources, also sources_.
    void eddyViscosity(bool sources, int depth);
    // Advances rho k and rho omega over tau under the SST sources at fixed rho and E, over N_depth.
    void turbulenceSource(const std::vector<Conserved>& state, std::vector<double>& turbulence,
                          const std::vector<SourceCoefficients>& coefficients, double tau, int depth) const;
    void transportFluxes(std::vector<Conserved>& derivative, std::vector<double>& speciesDerivative,
                         std::vector<double>& turbulenceDerivative, BoundaryRates& rates);
};
template<class F> void Flow::eachCell(const std::vector<std::size_t>* list, F f) const {
    pool_->blocks(list ? list->size() : state_.size(), kCells, [&](std::size_t begin, std::size_t end, int worker) {
        if(list) for(std::size_t n = begin; n < end; ++n) f((*list)[n], worker);
        else for(std::size_t q = begin; q < end; ++q) f(q, worker);
    });
}
} // namespace crucible
