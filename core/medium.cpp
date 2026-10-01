#include "core/medium.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace crucible {
namespace {
const std::array<double, 7>& range(const Species& s, double t) { return t<s.tMid?s.low:s.high; }
double cpOverR(const std::array<double, 7>& a, double t) {
    return a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4])));
}
double hOverR(const std::array<double, 7>& a, double t) {
    return t*(a[0]+t*(a[1]/2+t*(a[2]/3+t*(a[3]/4+t*a[4]/5))))+a[5];
}
}

Medium::Medium(std::vector<Species> species):species_(std::move(species)) {
    if(species_.empty()) throw std::invalid_argument("A medium needs at least one species.");
    for(const auto& s:species_) {
        if(!(s.molarMass>0) || !(s.tMid>0)) throw std::invalid_argument("Species molar mass and mid temperature must be positive.");
        specificR_.push_back(universalGasConstant/s.molarMass);
    }
}
Medium Medium::perfectGas(Gas gas) {
    Species s{"gas",universalGasConstant/gas.specificR,1000,{},{}};
    s.low[0]=s.high[0]=gas.gamma/(gas.gamma-1);
    return Medium({s});
}
Medium::Properties Medium::properties(double t,const double* y) const {
    Properties out;
    for(std::size_t k=0;k<size();++k) {
        const auto& a=range(species_[k],t);double w=y[k]*specificR_[k];
        out.r+=w;out.e+=w*(hOverR(a,t)-t);out.cv+=w*(cpOverR(a,t)-1);
    }
    return out;
}
double Medium::gasConstant(const double* y) const {
    double r=0;
    for(std::size_t k=0;k<size();++k) r+=y[k]*specificR_[k];
    return r;
}
double Medium::enthalpy(double t,const double* y) const {
    double h=0;
    for(std::size_t k=0;k<size();++k) h+=y[k]*specificR_[k]*hOverR(range(species_[k],t),t);
    return h;
}
double Medium::internalEnergy(double t,const double* y) const {
    return enthalpy(t,y)-gasConstant(y)*t;
}
double Medium::cv(double t,const double* y) const {
    double c=0;
    for(std::size_t k=0;k<size();++k) c+=y[k]*specificR_[k]*(cpOverR(range(species_[k],t),t)-1);
    return c;
}
double Medium::soundSpeed(double t,const double* y) const {
    double c=cv(t,y),r=gasConstant(y);
    return std::sqrt((c+r)/c*r*t);
}
double Medium::energyFloor(const double* y) const {
    double e=0;
    for(std::size_t k=0;k<size();++k) e+=y[k]*specificR_[k]*species_[k].low[5];
    return e;
}
double Medium::energyFloorOfClipped(const double* w) const {
    double e=0;
    for(std::size_t k=0;k<size();++k) e+=std::max(w[k],0.0)*specificR_[k]*species_[k].low[5];
    return e;
}
double Medium::temperature(double e,const double* y,double guess) const {
    if(!(e>energyFloor(y)) || !std::isfinite(e)) return std::numeric_limits<double>::quiet_NaN();
    // e(T) is increasing (cv > 0): keep a bracket [lo, hi] and fall back to bisection when a
    // Newton step leaves it.
    double lo=0,hi=std::numeric_limits<double>::infinity();
    double t=guess>0 && std::isfinite(guess)?guess:300;
    for(int n=0;n<100;++n) {
        auto q=properties(t,y);
        double f=q.e-e;
        double next=t-f/q.cv;
        if(std::abs(next-t)<=1e-13*t) return next;
        if(f>0) hi=t; else lo=t;
        if(!(next>lo && next<hi)) next=std::isfinite(hi)?0.5*(lo+hi):2*t;
        t=next;
    }
    throw std::runtime_error("Temperature iteration did not converge.");
}
} // namespace crucible
