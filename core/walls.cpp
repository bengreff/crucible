#include "core/walls.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace crucible {
std::vector<double> wallDistance(const Mesh& m, bool sideWall, const std::vector<bool>& plate,
                                 std::vector<std::array<double,2>>* direction) {
    if(!plate.empty() && plate.size()!=static_cast<std::size_t>(m.nr))
        throw std::invalid_argument("Plate wall flags must be empty or one per ring.");
    struct Segment { double z0, r0, z1, r1; };
    std::vector<Segment> walls;
    if(sideWall) for(int i=0;i<m.nz;++i) walls.push_back({i*m.dz,m.radius[i],(i+1)*m.dz,m.radius[i+1]});
    for(std::size_t j=0;j<plate.size();++j)
        if(plate[j]) walls.push_back({0,m.fraction[j]*m.radius[0],0,m.fraction[j+1]*m.radius[0]});
    std::vector<double> out(m.cells.size(),std::numeric_limits<double>::infinity());
    if(direction) direction->assign(m.cells.size(),{0.0,0.0});
    for(std::size_t q=0;q<m.cells.size();++q) {
        const double z=m.cells[q].z,r=m.cells[q].r;
        for(const auto& s:walls) {
            // Foot of the perpendicular, clamped to the segment.
            double ez=s.z1-s.z0,er=s.r1-s.r0;
            double t=std::clamp(((z-s.z0)*ez+(r-s.r0)*er)/(ez*ez+er*er),0.0,1.0);
            const double dz=z-s.z0-t*ez,dr=r-s.r0-t*er,d=std::hypot(dz,dr);
            if(d<out[q]) { out[q]=d;if(direction) (*direction)[q]={dz/d,dr/d}; }
        }
    }
    return out;
}

namespace wallLaw {
namespace {
// e^x minus its Taylor polynomial of degree n - 1, for x >= 0: the series where it would cancel.
double expTail(double x,int n) {
    if(x>2) {
        double poly=0,term=1;
        for(int k=0;k<n;++k) { poly+=term;term*=x/(k+1); }
        return std::exp(x)-poly;
    }
    double term=1,sum=0;
    for(int k=1;k<=n;++k) term*=x/k;
    for(int k=n;term>1e-18*sum || k==n;++k) { sum+=term;term*=x/(k+1); }
    return sum;
}
double atanOverT(double t) { return std::abs(t)<1e-4?1-t*t/3+t*t*t*t/5:std::atan(t)/t; }
// exp(e) - (1 + x + ... + x^(n-1) / (n-1)!), x >= 0, in whichever of two forms rounds less:
// exp(x) expm1(e - x) + tail(x), exact in sign and free of cancellation when e >= x (both terms
// non-negative) and near the wall; or exp(e) - polynomial(x) when e < x and x is large, where the
// first form subtracts two terms of size exp(x) (with heating u_eq+ < u+, and at u+ 300 the first
// form would lose every digit).
double excess(double x,double e,int n) {
    const double a1=std::exp(x)*std::expm1(e-x),a2=expTail(x,n);
    if(e>=x) return a1+a2;
    double poly=0,term=1;
    for(int k=0;k<n;++k) { poly+=term;term*=x/(k+1); }
    const double b1=std::exp(e);
    return std::abs(a1)+a2<=b1+poly?a1+a2:b1-poly;
}
// dy+/du+ - 1 = kappa exp(-kappa B) [exp(kappa u_eq) / D - (1 + x + x^2 / 2)], D = sqrt(1 + beta u - Gamma u^2).
double derivativeExcess(double u,double gamma,double beta,Constants c) {
    const double x=c.kappa*u,ue=equivalentVelocity(u,gamma,beta);
    return c.kappa*std::exp(-c.kappa*c.b)*excess(x,c.kappa*ue-0.5*std::log1p(beta*u-gamma*u*u),3);
}
}
// With D = sqrt(1 + beta u - Gamma u^2), the angle difference asin(x_u) - asin(x_0) of the printed
// form has sine 2 sqrt(Gamma) u N / ((1 + D) Q^2) and cosine C / Q^2, where
// N = Gamma (2 + 2D - beta u) + beta^2 and C = Gamma (4D - 2 beta u) + beta^2. For C > 0 the integral is
// atan(sqrt(Gamma) g) / sqrt(Gamma) with g = 2 u (N / C) / (1 + D), which is g at Gamma = 0
// (2u / (1 + D) = (2 / beta)(D - 1), exact) and u at Gamma = beta = 0.
double equivalentVelocity(double u,double gamma,double beta) {
    if(gamma==0 && beta==0) return u;
    const double d=std::sqrt(1+beta*u-gamma*u*u);
    const double n=gamma*(2+2*d-beta*u)+beta*beta,cosine=gamma*(4*d-2*beta*u)+beta*beta;
    if(cosine>0) { const double g=2*u*(n/cosine)/(1+d);return g*atanOverT(std::sqrt(gamma)*g); }
    return std::atan2(2*std::sqrt(gamma)*u*n/(1+d),cosine)/std::sqrt(gamma);
}
// y+ = u + exp(-kappa B) [exp(kappa u_eq) - (1 + x + x^2 / 2 + x^3 / 6)].
double yPlus(double u,double gamma,double beta,Constants c) {
    return u+std::exp(-c.kappa*c.b)*excess(c.kappa*u,c.kappa*equivalentVelocity(u,gamma,beta),4);
}
double yPlusDerivative(double u,double gamma,double beta,Constants c) { return 1+derivativeExcess(u,gamma,beta,c); }
// Printed eq. 10.13: 2 y+_White (kappa sqrt(Gamma) / Q) [1 - (2 Gamma u - beta)^2 / Q^2]^(1/2) = 4 kappa Gamma D y+_White / Q^2.
double yPlusDerivativePrinted(double u,double gamma,double beta,Constants c) {
    const double d=std::sqrt(1+beta*u-gamma*u*u),q2=beta*beta+4*gamma;
    const double white=std::exp(c.kappa*(equivalentVelocity(u,gamma,beta)-c.b));
    const double x=c.kappa*u;
    return 1+c.kappa*white*(q2>0?4*gamma*d/q2:1.0)-c.kappa*std::exp(-c.kappa*c.b)*(1+x+x*x/2);
}
double eddyRatio(double u,double gamma,double beta,double viscosityRatio,Constants c) {
    return (1-viscosityRatio)+derivativeExcess(u,gamma,beta,c);
}
// The exponential less its polynomial is summed as its series below x = 1, where the difference
// would cancel (eq. 10.11 at u+ 0.01 has x^3 / 6 = 1e-8 of the terms subtracted).
namespace {
double referenceTail(double x,int n) {
    double term=1,poly=0;
    for(int k=0;k<n;++k) { poly+=term;term*=x/(k+1); }
    if(x>1) return std::exp(x)-poly;
    double sum=0;
    for(int k=n;k<n+40;++k) { sum+=term;term*=x/(k+1); }
    return sum;
}
}
double spalding(double u,Constants c) { return u+std::exp(-c.kappa*c.b)*referenceTail(c.kappa*u,4); }
double spaldingEddy(double u,Constants c) { return c.kappa*std::exp(-c.kappa*c.b)*referenceTail(c.kappa*u,3); }
// u+ y+(u+) - Re_1 = H(u+) with y+(u+) = u+ + exp(-kappa B) [exp(x I) - (1 + x + x^2 / 2 + x^3 / 6)] and
// H' = y+ + u+ dy+/du+ along this path, dy+/du+ = 1 + kappa exp(-kappa B) [I exp(x I) - (1 + x + x^2 / 2)].
// The bracket doubles up from u+ = 1, and Newton starts at its top, where H >= 0.
Solution solveIsothermal(const Layer& l,Constants c) {
    if(!(l.u1>=0 && l.y1>0 && l.t1>0 && l.tw>0 && l.rhoW>0 && l.muW>0 && l.kW>0 && l.cp>0 && l.recovery>0))
        throw std::invalid_argument("Invalid wall-function layer.");
    const double g=l.recovery*l.u1*l.u1/(2*l.cp*l.tw),theta=l.t1/l.tw-1+g;
    const double reynolds=l.rhoW*l.u1*l.y1/l.muW,scale=std::exp(-c.kappa*c.b);
    if(reynolds==0) return {0,0,0,0,0,0,theta*l.tw*l.kW/l.y1,0};
    const double ratio=equivalentVelocity(1,g,theta);
    const double logRatio=std::log(ratio);
    auto law=[&](double u,double& y,double& dy) {
        const double x=c.kappa*u;
        y=u+scale*excess(x,x*ratio,4);
        dy=1+c.kappa*scale*excess(x,x*ratio+logRatio,3);
    };
    double lo=0,hi=1,y,dy;
    for(law(hi,y,dy);hi*y<reynolds;law(hi,y,dy)) { lo=hi;hi*=2;if(hi>1e4) throw std::runtime_error("Wall function: no bracket."); }
    double u=hi;
    // A Newton step below 1e-13 u+ is taken and ends the solve (it may land on the bracket's end);
    // a larger one outside the bracket is replaced by bisection.
    int it=0;
    for(;it<200;++it) {
        law(u,y,dy);
        const double h=u*y-reynolds,dh=y+u*dy;
        if(h==0) break;
        if(h<0) lo=u; else hi=u;
        double next=u-h/dh;
        if(dh>0 && std::abs(next-u)<=1e-13*u) { u=next;break; }
        if(!(next>lo && next<hi) || !(dh>0)) next=0.5*(lo+hi);
        u=next;
        if(hi-lo<=1e-15*hi) break;
    }
    if(it==200) throw std::runtime_error("Wall function: no convergence.");
    law(u,y,dy);
    // y = Re_1 / u+ at the root, so u_tau = u_1 / u+ and y+ / u+ = Re_1 / u+^2 stay finite as u_1 -> 0.
    const double over=y/u;
    return {l.u1/u,u,y,g/(u*u),theta/u,l.muW*l.u1*over/l.y1,theta*l.tw*l.kW*over/l.y1,it+1};
}
}
}
