#include "core/contour.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace crucible {
namespace {
double sign(double x) { return x>0?1.0:x<0?-1.0:0.0; }
// SciPy's PchipInterpolator._edge_case: the three-point end derivative, kept shape-preserving.
double edge(double h0,double h1,double m0,double m1) {
    double d=((2*h0+h1)*m0-h0*m1)/(h0+h1);
    if(sign(d)!=sign(m0)) return 0;
    if(sign(m0)!=sign(m1) && std::abs(d)>3*std::abs(m0)) return 3*m0;
    return d;
}
}  // namespace

Pchip::Pchip(std::vector<std::array<double,2>> points):points_(std::move(points)) {
    const std::size_t n=points_.size();
    if(n<2) throw std::invalid_argument("Pchip: at least two points");
    for(std::size_t k=1;k<n;++k)
        if(!(points_[k][0]>points_[k-1][0])) throw std::invalid_argument("Pchip: z must increase strictly");
    std::vector<double> h(n-1),m(n-1);
    for(std::size_t k=0;k+1<n;++k) { h[k]=points_[k+1][0]-points_[k][0];m[k]=(points_[k+1][1]-points_[k][1])/h[k]; }
    slopes_.assign(n,0.0);
    if(n==2) { slopes_[0]=slopes_[1]=m[0];return; }
    for(std::size_t k=1;k+1<n;++k) {
        if(sign(m[k-1])!=sign(m[k]) || m[k-1]==0 || m[k]==0) continue;
        const double w1=2*h[k]+h[k-1],w2=h[k]+2*h[k-1];
        slopes_[k]=(w1+w2)/(w1/m[k-1]+w2/m[k]);
    }
    slopes_[0]=edge(h[0],h[1],m[0],m[1]);
    slopes_[n-1]=edge(h[n-2],h[n-3],m[n-2],m[n-3]);
}
std::size_t Pchip::interval(double z) const {
    auto it=std::upper_bound(points_.begin(),points_.end(),z,[](double v,const std::array<double,2>& p){return v<p[0];});
    const auto k=static_cast<std::size_t>(std::max<std::ptrdiff_t>(it-points_.begin()-1,0));
    return std::min(k,points_.size()-2);
}
double Pchip::operator()(double z) const {
    const std::size_t k=interval(z);
    const double h=points_[k+1][0]-points_[k][0],t=(z-points_[k][0])/h;
    const double r0=points_[k][1],r1=points_[k+1][1];
    return (2*t*t*t-3*t*t+1)*r0+(t*t*t-2*t*t+t)*h*slopes_[k]+(-2*t*t*t+3*t*t)*r1+(t*t*t-t*t)*h*slopes_[k+1];
}
double Pchip::derivative(double z) const {
    const std::size_t k=interval(z);
    const double h=points_[k+1][0]-points_[k][0],t=(z-points_[k][0])/h;
    const double r0=points_[k][1],r1=points_[k+1][1];
    return (6*t*t-6*t)*(r0-r1)/h+(3*t*t-4*t+1)*slopes_[k]+(3*t*t-2*t)*slopes_[k+1];
}
double Pchip::secondLeft(std::size_t k) const {
    if(k==0) throw std::out_of_range("Pchip::secondLeft: no interval ends at the first point");
    const double h=points_[k][0]-points_[k-1][0],secant=(points_[k][1]-points_[k-1][1])/h;
    return (-6*secant+2*slopes_[k-1]+4*slopes_[k])/h;
}
double Pchip::secondRight(std::size_t k) const {
    if(k+1>=points_.size()) throw std::out_of_range("Pchip::secondRight: no interval starts at the last point");
    const double h=points_[k+1][0]-points_[k][0],secant=(points_[k+1][1]-points_[k][1])/h;
    return (6*secant-4*slopes_[k]-2*slopes_[k+1])/h;
}
std::vector<std::array<double,2>> Pchip::sample(int perInterval) const {
    std::vector<std::array<double,2>> out;
    for(std::size_t k=0;k+1<points_.size();++k) {
        out.push_back(points_[k]);
        const double h=points_[k+1][0]-points_[k][0];
        for(int s=1;s<perInterval;++s) {
            const double z=points_[k][0]+h*s/perInterval;
            out.push_back({z,(*this)(z)});
        }
    }
    out.push_back(points_.back());
    return out;
}
}  // namespace crucible
