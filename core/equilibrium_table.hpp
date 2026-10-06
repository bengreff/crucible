#pragma once
#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace crucible {
// Tabulated constant-(u, v) shifting equilibrium of H/O mixtures with inert diluents (TECHNICAL_PLAN,
// *Lightweight engine*, Table A; criteria in docs/evidence/TABLE_A.md). Built offline from Cantera
// by crucible_build_equilibrium_table; no Cantera here.
//
// Species made of H and O alone react; every other species is the only carrier of its elements (N2,
// Ar in h2o2.yaml; the builder checks), so element conservation fixes its amount. In an ideal gas
// each species' chemical potential depends on T and its own partial pressure, so the reacting
// species' equilibrium at a given T depends only on their own elements and partial density, not on
// the diluent. Axes: f = Z_H / (Z_H + Z_O); the temperature T; ln rho_r, rho_r = rho (Z_H + Z_O) the
// reacting species' partial density. Each node holds the reacting mixture's equilibrium at its
// (f, T, rho_r): [e, Y_1..Y_ns] per unit reacting mass, e the specific internal energy. A lookup
// finds the T at which the reacting mass's interpolated e plus the diluents' e equals the query's e
// (increasing in T), then interpolates Y there. Z_H and Z_O are linear in f, so the returned Y
// carries the query's elements to rounding.
struct EquilibriumTable {
    std::vector<std::string> species;
    std::vector<std::array<double,2>> elements;  // per species: mass fractions of H and O in it (0, 0 if inert)
    std::vector<std::size_t> inert;              // the inert species
    std::vector<double> f;                       // node coordinates, increasing, from 0 to 1
    int nT{}, nRho{};
    double tMin{}, tMax{}, lnRhoMin{}, lnRhoMax{};
    std::vector<double> inertEnergy;             // per inert species and T node: e_k(T) [J/kg]
    std::vector<double> values;                  // per (f, T, rho_r) node, rho_r fastest: e [J/kg], Y_1..Y_ns
    std::string provenance;

    [[nodiscard]] std::size_t width() const { return species.size()+1; }
    [[nodiscard]] std::size_t node(std::size_t i, std::size_t k, std::size_t l) const {
        return (i*static_cast<std::size_t>(nT)+k)*static_cast<std::size_t>(nRho)+l;
    }
    [[nodiscard]] double temperature(std::size_t k) const { return tMin+(tMax-tMin)*static_cast<double>(k)/(nT-1); }
    [[nodiscard]] double lnRho(std::size_t l) const { return lnRhoMin+(lnRhoMax-lnRhoMin)*static_cast<double>(l)/(nRho-1); }

    enum class Result { Table, Clamped, Outside };
    // Equilibrium mass fractions (out, size ns) and temperature t at the elements of y (mass
    // fractions summing to 1), specific internal energy e and density rho. Clamped: e lies below the
    // tMin equilibrium's, which is returned (the low-temperature limit: complete combustion), or rho_r
    // lies below the axis and is read at its end (the reacting mass is then below exp(lnRhoMin) per
    // m^3). Outside: e above the tMax equilibrium's, or rho_r above the axis; out is then unusable.
    // t is a consistency check only: the caller's temperature comes from e and the returned Y.
    Result lookup(const double* y, double e, double rho, double* out, double& t) const;

    void save(const std::string& path) const;
    static EquilibriumTable load(const std::string& path);
};
} // namespace crucible
