#pragma once
#include "core/flow.hpp"
#include <vector>

namespace crucible {
// Exact distance from each cell centroid to the no-slip walls. For a surface of revolution the
// nearest point lies in the point's own meridional half-plane, so the distance is the plane
// distance to the wall segments: the side wall, the polyline through (i dz, radius[i]) when
// sideWall, and ring j of the injector face (z = 0, from fraction[j] to fraction[j + 1] of
// radius[0]) where plate[j]. Infinity where there is no wall at all.
std::vector<double> wallDistance(const Mesh& mesh, bool sideWall, const std::vector<bool>& plate);
}
