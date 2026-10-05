#include "core/walls.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace crucible {
std::vector<double> wallDistance(const Mesh& m, bool sideWall, const std::vector<bool>& plate) {
    if(!plate.empty() && plate.size()!=static_cast<std::size_t>(m.nr))
        throw std::invalid_argument("Plate wall flags must be empty or one per ring.");
    struct Segment { double z0, r0, z1, r1; };
    std::vector<Segment> walls;
    if(sideWall) for(int i=0;i<m.nz;++i) walls.push_back({i*m.dz,m.radius[i],(i+1)*m.dz,m.radius[i+1]});
    for(std::size_t j=0;j<plate.size();++j)
        if(plate[j]) walls.push_back({0,m.fraction[j]*m.radius[0],0,m.fraction[j+1]*m.radius[0]});
    std::vector<double> out(m.cells.size(),std::numeric_limits<double>::infinity());
    for(std::size_t q=0;q<m.cells.size();++q) {
        const double z=m.cells[q].z,r=m.cells[q].r;
        for(const auto& s:walls) {
            // Foot of the perpendicular, clamped to the segment.
            double ez=s.z1-s.z0,er=s.r1-s.r0;
            double t=std::clamp(((z-s.z0)*ez+(r-s.r0)*er)/(ez*ez+er*er),0.0,1.0);
            out[q]=std::min(out[q],std::hypot(z-s.z0-t*ez,r-s.r0-t*er));
        }
    }
    return out;
}
}
