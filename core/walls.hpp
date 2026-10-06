#pragma once
#include "core/flow.hpp"
#include <array>
#include <vector>

namespace crucible {
// Exact distance from each cell centroid to the no-slip walls. For a surface of revolution the
// nearest point lies in the point's own meridional half-plane, so the distance is the plane
// distance to the wall segments: the side wall, the polyline through (i dz, radius[i]) when
// sideWall, and ring j of the injector face (z = 0, from fraction[j] to fraction[j + 1] of
// radius[0]) where plate[j]. Infinity where there is no wall at all.
// direction, if given, receives each centroid's unit vector away from its nearest wall point (the
// gradient of the distance; zero where there is no wall).
std::vector<double> wallDistance(const Mesh& mesh, bool sideWall, const std::vector<bool>& plate,
                                 std::vector<std::array<double, 2>>* direction = nullptr);

// Wall functions: Nichols and Nelson's compressible, heated law (R. H. Nichols, Turbulence Models
// and Their Application to Complex Flows, rev. 4.01, ch. 10; docs/evidence/WALL_FUNCTIONS.md).
// Spalding's single formula over sublayer, buffer and log layer with its log term replaced by White
// and Christoph's law (eq. 10.6), and Crocco-Busemann (eq. 10.8), T = T_w (1 + beta u+ - Gamma u+^2):
//   y+ = u+ + exp(kappa u_eq+ - kappa B) - exp(-kappa B) (1 + x + x^2 / 2 + x^3 / 6),  x = kappa u+,
//   u_eq+ = integral from 0 to u+ of ds / sqrt(1 + beta s - Gamma s^2)  (the van Driest velocity),
//   Gamma = r u_tau^2 / (2 cp T_w),  beta = q_w mu_w / (rho_w T_w k_w u_tau),  r = Pr_w^(1/3).
// The printed form (1 / sqrt(Gamma)) [asin((2 Gamma u+ - beta) / Q) - phi] is this integral; it is
// evaluated here as an angle difference through atan2 so the limit Gamma -> 0 has no cancellation.
// The derivative is dy+/du+ = 1 + kappa y+_White / sqrt(1 + beta u+ - Gamma u+^2) - kappa exp(-kappa B)
// (1 + x + x^2 / 2): eq. 10.13 with the exponent -1/2 that differentiating eq. 10.7 gives (the
// printed +1/2 is kept only to report its error). q_w is the heat flux from the gas into the wall.
namespace wallLaw {
struct Constants { double kappa{0.4}, b{5.5}; };
// u_eq+ above, for 1 + beta s - Gamma s^2 > 0 on [0, u].
double equivalentVelocity(double u, double gamma, double beta);
// Eq. 10.6, and dy+/du+ at fixed Gamma and beta (corrected eq. 10.13, and the printed one).
double yPlus(double u, double gamma, double beta, Constants c);
double yPlusDerivative(double u, double gamma, double beta, Constants c);
double yPlusDerivativePrinted(double u, double gamma, double beta, Constants c);
// Eq. 10.12, mu_t / mu_w = dy+/du+ - mu_1 / mu_w (mu_1 the first cell's viscosity), corrected.
double eddyRatio(double u, double gamma, double beta, double viscosityRatio, Constants c);
// Spalding's eq. 10.3 and the eddy viscosity eq. 10.11, written out directly (criterion 0(a)).
double spalding(double u, Constants c);
double spaldingEddy(double u, Constants c);
// One wall face: the first cell's tangential speed, wall distance and temperature; the wall
// temperature; rho_w (the cell's pressure at T_w), mu_w and k_w at T_w with the cell's composition;
// cp the layer mean (h(T_1) - h(T_w)) / (T_1 - T_w); r = Pr_w^(1/3).
struct Layer { double u1, y1, t1, tw, rhoW, muW, kW, cp, recovery; };
struct Solution { double uTau, uPlus, yPlus, gamma, beta, shear, heat; int iterations; };
// Isothermal wall: with u+ = u_1 / u_tau, Crocco-Busemann at T_1 fixes beta u+ = T_1 / T_w - 1 + G and
// Gamma u+^2 = G = r u_1^2 / (2 cp T_w), so u_eq+ = u+ I with I = u_eq+(1; G, beta u+) and eq. 10.6
// becomes one equation in u+, u+ y+(u+) = rho_w u_1 y_1 / mu_w, solved by bracketed Newton-bisection
// to 1e-13 relative. shear = rho_w u_tau^2 and heat = q_w, both written so u_1 -> 0 is regular.
Solution solveIsothermal(const Layer& layer, Constants c);
}
}
