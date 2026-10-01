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
    [[nodiscard]] double cv(double t, const double* y) const;             // J/(kg K)
    // Frozen sound speed sqrt(cp/cv R T).
    [[nodiscard]] double soundSpeed(double t, const double* y) const;
    // Internal energy as T -> 0 on the low-range polynomials: admissible states lie above it.
    [[nodiscard]] double energyFloor(const double* y) const;
    // sum_k max(w_k, 0) R_k a5_k for unnormalised weights w (e.g. partial densities).
    [[nodiscard]] double energyFloorOfClipped(const double* w) const;
    // Temperature with internal energy e (safeguarded Newton from guess); NaN if e is not above the floor.
    [[nodiscard]] double temperature(double e, const double* y, double guess) const;
private:
    std::vector<Species> species_;
    std::vector<double> specificR_;
};
} // namespace crucible
