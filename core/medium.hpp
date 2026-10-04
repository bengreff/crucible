#pragma once
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace crucible {
struct Gas { double gamma{1.4}, specificR{287.05}; };
// Thermally perfect species: NASA 7-coefficient polynomials, low range below tMid, high above.
// cp/R = a0 + a1 T + a2 T^2 + a3 T^3 + a4 T^4,  h/(R T) = a0 + a1 T/2 + ... + a4 T^4/5 + a5/T.
struct Species {
    std::string name;
    double molarMass{};  // kg/kmol
    double tMid{1000};
    std::array<double, 7> low{}, high{};
};
// Molecular transport fits in the form of Cantera's GasTransport (degree 4 in ln T):
// sqrt(mu_k / sqrt(T)) [Pa s], lambda_k / sqrt(T) [W/(m K)], and D_kj p / T^1.5 [m^2 Pa/s] for the
// species pairs k <= j in row order (00, 01, ..., 0n, 11, ...).
struct TransportFits {
    std::vector<std::array<double, 5>> viscosity, conductivity, diffusion;
};
// Ideal-gas mixture of thermally perfect species. Every gas in the core is one of these; a
// calorically perfect gas is a single species with constant cp (Medium::perfectGas).
class Medium {
public:
    // Mixture gas constant [J/(kg K)], internal energy [J/kg] and cv [J/(kg K)] at one temperature.
    struct Properties { double r{}, e{}, cv{}; };
    static constexpr double universalGasConstant=8314.46261815324;  // J/(kmol K), exact SI (N_A k_B)
    Medium()=default;
    explicit Medium(std::vector<Species> species);
    static Medium perfectGas(Gas gas);
    [[nodiscard]] std::size_t size() const { return species_.size(); }
    [[nodiscard]] const std::vector<Species>& species() const { return species_; }
    [[nodiscard]] Properties properties(double t, const double* y) const;  // one pass over the species
    [[nodiscard]] double gasConstant(const double* y) const;              // J/(kg K)
    [[nodiscard]] double internalEnergy(double t, const double* y) const; // J/kg, NASA reference
    [[nodiscard]] double enthalpy(double t, const double* y) const;       // J/kg
    void speciesEnthalpies(double t, double* h) const;                    // J/kg, per species
    [[nodiscard]] double cv(double t, const double* y) const;             // J/(kg K)
    // Frozen sound speed sqrt(cp/cv R T).
    [[nodiscard]] double soundSpeed(double t, const double* y) const;
    // Internal energy as T -> 0 on the low-range polynomials: admissible states lie above it.
    [[nodiscard]] double energyFloor(const double* y) const;
    // sum_k max(w_k, 0) R_k a5_k for unnormalised weights w (e.g. partial densities).
    [[nodiscard]] double energyFloorOfClipped(const double* w) const;
    // Temperature with internal energy e (safeguarded Newton from guess); NaN if e is not above the floor.
    [[nodiscard]] double temperature(double e, const double* y, double guess) const;

    // Mixture-averaged molecular transport, the same formulas as Cantera's "mixture-averaged"
    // model: Wilke viscosity, conductivity as the mean of the series and parallel averages, and
    // mass-based mixture diffusion coefficients D_km = (1 - Y_k) / sum_{j != k} X_j / D_jk, with
    // mole fractions floored at 1e-20 as Cantera does (1 - Y_k is summed from the other species).
    struct Transport { double viscosity{}, conductivity{}; };
    void setTransport(TransportFits fits);
    [[nodiscard]] bool hasTransport() const { return !fits_.viscosity.empty(); }
    // Writes D_km [m^2/s] to diffusion (size()); work is scratch reused between calls.
    Transport transport(double t, double p, const double* y, double* diffusion, std::vector<double>& work) const;
private:
    std::vector<Species> species_;
    std::vector<double> specificR_;
    TransportFits fits_;
    std::vector<double> wilkeMass_, wilkeRoot_;  // (W_j/W_k)^(1/4) and sqrt(8 (1 + W_k/W_j)), row k
};
} // namespace crucible
