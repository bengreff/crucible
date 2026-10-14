#pragma once
#include <array>
#include <cstddef>
#include <vector>

namespace crucible {
// A wall contour r(z) through given points (z strictly increasing): the monotone piecewise cubic
// Hermite interpolant (PCHIP; Fritsch and Carlson, SIAM J. Numer. Anal. 17, 1980) with the node
// derivatives of SciPy's PchipInterpolator. Inside, the derivative is zero where the neighbouring
// secant slopes differ in sign or one is zero, else their weighted harmonic mean (Fritsch and
// Butland); at the ends, the shape-preserving three-point formula. It passes through every point and
// is monotone between neighbours, so a local minimum of r (a throat) is one of the given points.
// Used for area-against-station tables, r = r_t sqrt(A / A_t) (TECHNICAL_PLAN, *Geometry as data*).
class Pchip {
public:
    explicit Pchip(std::vector<std::array<double,2>> points);
    [[nodiscard]] double operator()(double z) const;
    [[nodiscard]] double derivative(double z) const;
    // The second derivative at point k from the interval ending there (left) and the one starting there (right).
    [[nodiscard]] double secondLeft(std::size_t k) const;
    [[nodiscard]] double secondRight(std::size_t k) const;
    // Samples for Definition::contour: every point, and perInterval - 1 equally spaced ones inside each interval.
    [[nodiscard]] std::vector<std::array<double,2>> sample(int perInterval) const;
    [[nodiscard]] const std::vector<std::array<double,2>>& points() const { return points_; }
    [[nodiscard]] const std::vector<double>& slopes() const { return slopes_; }
private:
    [[nodiscard]] std::size_t interval(double z) const;
    std::vector<std::array<double,2>> points_;
    std::vector<double> slopes_;
};
}  // namespace crucible
