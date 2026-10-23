#include "core/flow.hpp"
#include "core/walls.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>

namespace crucible {
namespace {
constexpr double pi=std::numbers::pi;
double sq(double x) { return x*x; }
Conserved values(Primitive w) { return {w.rho,w.uz,w.ur,w.p}; }
Primitive unpack(Conserved w) { return {w[0],w[1],w[2],w[3]}; }
double minmod(double a, double b) {
    return a*b<=0 ? 0 : std::copysign(std::min(std::abs(a),std::abs(b)),a);
}
Conserved physicalFlux(Primitive w, double e, double nz, double nr) {
    auto u=conservative(w,e);
    double un=w.uz*nz+w.ur*nr;
    return {w.rho*un,u[1]*un+w.p*nz,u[2]*un+w.p*nr,(u[3]+w.p)*un};
}
Primitive reflect(Primitive w, double nz, double nr) {
    double un=w.uz*nz+w.ur*nr;
    w.uz-=2*un*nz; w.ur-=2*un*nr;
    return w;
}
// Positive density, finite values, and internal energy above the composition's floor.
bool admissibleBulk(const Conserved& u, double floor) {
    if (!(u[0]>0)) return false;
    for(double x:u) if(!std::isfinite(x)) return false;
    return u[3]-(sq(u[1])+sq(u[2]))/(2*u[0])>floor*u[0];
}
// Mass fractions from partial densities: negatives clipped, then normalised.
void normalise(const double* partial, std::size_t n, double* y) {
    double sum=0;
    for(std::size_t k=0;k<n;++k) { y[k]=std::max(partial[k],0.0); sum+=y[k]; }
    if(!(sum>0)) throw std::runtime_error("Non-admissible gas state: no species mass in a cell.");
    for(std::size_t k=0;k<n;++k) y[k]/=sum;
}
// Low-Mach corrections act on interior faces and the mirror (slip) wall; inlet/outlet models keep their own states.
Conserved faceFlux(LowMach lowMach,Primitive l,Thermal tl,Primitive r,Thermal tr,double nz,double nr,bool interior) {
    if(lowMach==LowMach::HllcLm && interior) return hllcLm(l,tl,r,tr,nz,nr);
    if(lowMach==LowMach::Thornber && interior) thornberScale(l,r,tl.a,tr.a);
    return hllc(l,tl,r,tr,nz,nr);
}
}

double Supply::opening(double time) const {
    if(time<opens) return 0;
    return ramp>0?std::min(1.0,(time-opens)/ramp):1.0;
}
double Definition::span() const { return contour.empty()?length:contour.back()[0]-contour.front()[0]; }
Medium Definition::medium() const {
    Medium m=species.empty()?Medium::perfectGas(gas):Medium(species);
    if(!transport.viscosity.empty()) m.setTransport(transport);
    return m;
}
std::vector<double> Definition::massFractions() const {
    std::size_t n=species.empty()?1:species.size();
    if(composition.empty()) { std::vector<double> y(n,0.0); y[0]=1; return y; }
    return composition;
}

void Definition::validate() const {
    if(unityLewis && !(turbulence.enabled && turbulence.consistentFlux && !transport.viscosity.empty()))
        throw std::invalid_argument("Unity Lewis number needs transport and turbulence with the consistent turbulent fluxes.");
    if((openAreaSupplies || plateTemperature) && experiment!=Case::Chamber)
        throw std::invalid_argument("Open-area supplies and the plate temperature are a Chamber's.");
    if(plateTemperature && (!(*plateTemperature>=0) || !std::isfinite(*plateTemperature)))
        throw std::invalid_argument("The plate temperature must be finite and non-negative.");
    if(nz<4 || nr<1 || static_cast<std::uint64_t>(nz)*nr>16000000)
        throw std::invalid_argument("Mesh must have at least 4 axial / 1 radial cells and at most 16 million cells.");
    for(double x:{length,inletRadius,throatRadius,exitRadius,totalPressure,totalTemperature,
                  backPressure,cfl,gas.specificR})
        if(!(x>0) || !std::isfinite(x)) throw std::invalid_argument("Lengths, thermodynamic inputs and CFL must be finite and positive.");
    if(!(gas.gamma>1 && gas.gamma<2) || !(cfl<=0.8) ||
       !(throatFraction>0.05 && throatFraction<0.95))
        throw std::invalid_argument("Invalid gamma, CFL, or throat location.");
    if(experiment==Case::Nozzle && contour.empty() && (throatRadius>inletRadius || throatRadius>exitRadius))
        throw std::invalid_argument("The throat must be no wider than inlet and exit.");
    const std::size_t n=species.empty()?1:species.size();
    auto checkComposition=[&](const std::vector<double>& y) {
        double sum=0;
        for(double v:y) { if(!(v>=0) || !std::isfinite(v)) throw std::invalid_argument("Mass fractions must be finite and non-negative."); sum+=v; }
        if(y.size()!=n || std::abs(sum-1)>1e-12)
            throw std::invalid_argument("Composition must give one mass fraction per species, summing to one.");
    };
    checkComposition(massFractions());
    if(!(wallTemperature>=0) || !std::isfinite(wallTemperature))
        throw std::invalid_argument("Wall temperature must be finite and non-negative (zero: adiabatic).");
    if(!(radialStretching>=0 && radialStretching<=20))
        throw std::invalid_argument("Radial stretching must lie in [0, 20] (zero: rings of equal height).");
    if(!(outletRelaxation>=0) || !std::isfinite(outletRelaxation))
        throw std::invalid_argument("Outlet relaxation must be finite and non-negative (zero: fixed back pressure).");
    if(turbulence.wallFunctions && !turbulence.enabled) throw std::invalid_argument("Wall functions need turbulence.");
    if(turbulence.propertyCorrections && (!turbulence.enabled || wallSlip))
        throw std::invalid_argument("The SST property corrections need turbulence and no-slip walls (S_n is defined from the wall).");
    if(turbulence.enabled) {
        const auto& t=turbulence;
        if(transport.viscosity.empty()) throw std::invalid_argument("Turbulence needs molecular transport.");
        if(experiment==Case::Nozzle)
            throw std::invalid_argument("Turbulent inflow at the nozzle inlet is not yet implemented: turbulence runs on Chamber, UniformDuct and ShockTube.");
        if(experiment==Case::Chamber)
            for(const auto& s:supplies)
                if(!(s.turbulenceIntensity>0 && s.viscosityRatio>0) || !std::isfinite(s.turbulenceIntensity+s.viscosityRatio))
                    throw std::invalid_argument("With turbulence every supply needs a positive turbulence intensity and viscosity ratio.");
        for(double x:{t.prandtl,t.schmidt,t.wallOmegaFactor,t.ambientOmega})
            if(!(x>0) || !std::isfinite(x)) throw std::invalid_argument("Turbulent Prandtl and Schmidt numbers, the wall omega factor and the ambient omega must be positive.");
        if(!(t.ambientK>=0) || !std::isfinite(t.ambientK)) throw std::invalid_argument("Ambient k must be finite and non-negative.");
        if(t.wallFunctions) {
            if(wallSlip) throw std::invalid_argument("Wall functions need no-slip walls.");
            if(!(t.wallKappa>0) || !std::isfinite(t.wallKappa) || !std::isfinite(t.wallB))
                throw std::invalid_argument("Wall-function kappa must be positive and B finite.");
        }
    }
    if(!contour.empty()) {
        if(contour.size()<2) throw std::invalid_argument("A wall contour needs at least two points.");
        for(std::size_t k=0;k<contour.size();++k) {
            if(!(contour[k][1]>0) || !std::isfinite(contour[k][0]) || !std::isfinite(contour[k][1]))
                throw std::invalid_argument("Contour radii must be positive and finite.");
            if(k>0 && !(contour[k][0]>contour[k-1][0])) throw std::invalid_argument("Contour z must increase strictly.");
        }
    }
    if(experiment==Case::Chamber) {
        if(!(ambientTemperature>0) || !std::isfinite(ambientTemperature)) throw std::invalid_argument("Ambient temperature must be positive.");
        for(const auto& s:supplies) {
            if(!(s.innerRadius>=0 && s.outerRadius>s.innerRadius) || !std::isfinite(s.outerRadius))
                throw std::invalid_argument("A supply needs 0 <= inner radius < outer radius.");
            if(!(s.massFlow>=0) || !std::isfinite(s.massFlow) || !(s.totalTemperature>0) || !std::isfinite(s.totalTemperature))
                throw std::invalid_argument("A supply needs a finite non-negative mass flow and a positive temperature.");
            if(!(s.opens>=0) || !(s.ramp>=0) || !std::isfinite(s.opens+s.ramp))
                throw std::invalid_argument("Valve opening time and ramp must be finite and non-negative.");
            checkComposition(s.composition);
        }
        const auto& g=igniter;
        if(!(g.energy>=0) || !std::isfinite(g.energy)) throw std::invalid_argument("Igniter energy must be finite and non-negative.");
        if(g.energy>0 && (!(g.duration>0) || !(g.start>=0) || !(g.zMax>=g.zMin) || !(g.rMax>0) ||
                          !std::isfinite(g.duration+g.start+g.zMax+g.zMin+g.rMax)))
            throw std::invalid_argument("Igniter needs a positive duration and radius, a non-negative start and zMin <= zMax.");
    }
}

Mesh::Mesh(const Definition& d):nz(d.nz),nr(d.nr),dz(d.span()/d.nz) {
    d.validate();
    radius.resize(nz+1);
    std::size_t segment=0;
    for(int i=0;i<=nz;++i) {
        double x=static_cast<double>(i)/nz;
        if(!d.contour.empty()) {
            // Linear interpolation of the contour table at the station.
            const auto& c=d.contour;
            double z=i==nz?c.back()[0]:c.front()[0]+i*dz;
            while(segment+2<c.size() && z>c[segment+1][0]) ++segment;
            double t=std::clamp((z-c[segment][0])/(c[segment+1][0]-c[segment][0]),0.0,1.0);
            radius[i]=c[segment][1]+t*(c[segment+1][1]-c[segment][1]);
        }
        else if(d.experiment!=Case::Nozzle) radius[i]=d.inletRadius;
        else if(x<d.throatFraction) {
            double t=x/d.throatFraction;
            radius[i]=d.throatRadius+(d.inletRadius-d.throatRadius)*(1+std::cos(pi*t))/2;
        } else {
            double t=(x-d.throatFraction)/(1-d.throatFraction);
            radius[i]=d.throatRadius+(d.exitRadius-d.throatRadius)*(1-std::cos(pi*t))/2;
        }
    }
    fraction.resize(nr+1);
    uniform=!(d.radialStretching>0);
    for(int j=0;j<=nr;++j) {
        double s=static_cast<double>(j)/nr;
        fraction[j]=uniform || j==0 || j==nr?s:std::tanh(d.radialStretching*s)/std::tanh(d.radialStretching);
    }
    cells.resize(static_cast<std::size_t>(nz)*nr);
    for(int i=0;i<nz;++i) for(int j=0;j<nr;++j) {
        double a=fraction[j],b=fraction[j+1];
        double r0=radius[i],r1=radius[i+1];
        double v=pi*dz/3*(r0*r0+r0*r1+r1*r1)*(b*b-a*a);
        // Exact radial moments of the frustum ring (wall radius linear in z across the cell).
        double first=2*pi/3*(b*b*b-a*a*a)*dz*(r0+r1)*(r0*r0+r1*r1)/4;
        double second=pi/2*(b*b*b*b-a*a*a*a)*dz*(r0*r0*r0*r0+r0*r0*r0*r1+r0*r0*r1*r1+r0*r1*r1*r1+r1*r1*r1*r1)/5;
        // The axial centroid lies toward the wider end where the wall slopes: (z, r) is then the volume
        // centroid, a point inside the cell. With the column centre and the radial centroid instead, a
        // thin ring at a sloping wall would sit outside the wall line.
        double centroid=(i+0.5)*dz+dz*(r1-r0)*(r1+r0)/(4*(r0*r0+r0*r1+r1*r1));
        cells[index(i,j)]={v,centroid,first/v,pi*dz*(r0+r1)*(b-a),second/v};
    }
    wallRadius=radius;wallDz=dz;wallStride=1;
    columnSquares_.resize(nz);columnCubes_.resize(nz);
    for(int i=0;i<nz;++i) {
        double r0=radius[i],r1=radius[i+1];
        columnSquares_[i]=dz*(r0*r0+r0*r1+r1*r1)/3;columnCubes_[i]=dz*(r0+r1)*(r0*r0+r1*r1)/4;
    }
    axialAreas_.resize(static_cast<std::size_t>(nz+1)*nr);ringMiddles_.resize(axialAreas_.size());
    for(int face=0;face<=nz;++face) for(int j=0;j<nr;++j) {
        const auto k=static_cast<std::size_t>(face)*nr+j;
        if(uniform) {
            axialAreas_[k]=pi*sq(radius[face])*(2*j+1)/sq(nr);
            ringMiddles_[k]=radius[face]*(j+0.5)/nr;
        } else {
            axialAreas_[k]=pi*sq(radius[face])*(sq(fraction[j+1])-sq(fraction[j]));
            ringMiddles_[k]=radius[face]*0.5*(fraction[j]+fraction[j+1]);
        }
    }
    const auto radialCount=static_cast<std::size_t>(nz)*(nr+1);
    radialAreaVectors_.resize(radialCount);radialMiddles_.resize(radialCount);
    radialRadii_.resize(radialCount);radialSecondMoments_.resize(radialCount);
    for(int i=0;i<nz;++i) for(int face=0;face<=nr;++face) {
        const auto k=radialIndex(i,face);
        double f=fraction[face];
        radialMiddles_[k]=uniform?0.5*(radius[i]+radius[i+1])*face/nr:0.5*(radius[i]+radius[i+1])*f;
        double r0=f*radius[i],r1=f*radius[i+1];
        radialAreaVectors_[k]={-pi*(r0+r1)*(r1-r0),pi*(r0+r1)*dz};
        // On a sloping column the reconstruction's faces sit on the same column radius as the centroids
        // (cell.r is fraction-moment times R_col = 3 (r0 + r1)(r0^2 + r1^2) / (4 (r0^2 + r0 r1 + r1^2))),
        // so each centroid lies between its faces. With the face's area-weighted radius instead, a thin
        // ring at a sloping wall had both faces on one side of its centroid and its slopes extrapolated
        // outside it.
        r0=radius[i];r1=radius[i+1];
        if(face==0) radialRadii_[k]=0;
        else if(r0==r1) radialRadii_[k]=2*f*(r0*r0+r0*r1+r1*r1)/(3*(r0+r1));
        else radialRadii_[k]=f*3*(r0+r1)*(r0*r0+r1*r1)/(4*(r0*r0+r0*r1+r1*r1));
        radialSecondMoments_[k]=f*f*(sq(radius[i])+sq(radius[i+1]))/2;
    }
}
Mesh Mesh::coarsened(const Mesh& finer) {
    if(finer.nz%2 || finer.nr%2 || finer.nz<2 || finer.nr<2)
        throw std::invalid_argument("A mesh coarsens only with an even number of columns and of rings.");
    Mesh m;
    m.nz=finer.nz/2;m.nr=finer.nr/2;m.dz=2*finer.dz;m.uniform=finer.uniform;
    m.radius.resize(m.nz+1);m.fraction.resize(m.nr+1);
    for(int i=0;i<=m.nz;++i) m.radius[i]=finer.radius[2*i];
    for(int j=0;j<=m.nr;++j) m.fraction[j]=finer.fraction[2*j];
    m.wallRadius=finer.wallRadius;m.wallDz=finer.wallDz;m.wallStride=2*finer.wallStride;
    m.cells.resize(static_cast<std::size_t>(m.nz)*m.nr);
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
        double v=0,z=0,r=0,pressure=0,second=0;
        for(int a=2*i;a<2*i+2;++a) for(int b=2*j;b<2*j+2;++b) {
            const auto& c=finer.cells[finer.index(a,b)];
            v+=c.volume;z+=c.volume*c.z;r+=c.volume*c.r;pressure+=c.radialPressureMeasure;second+=c.volume*c.radialSecondMoment;
        }
        m.cells[m.index(i,j)]={v,z/v,r/v,pressure,second/v};
    }
    m.columnSquares_.resize(m.nz);m.columnCubes_.resize(m.nz);
    for(int i=0;i<m.nz;++i) {
        m.columnSquares_[i]=finer.columnSquares_[2*i]+finer.columnSquares_[2*i+1];
        m.columnCubes_[i]=finer.columnCubes_[2*i]+finer.columnCubes_[2*i+1];
    }
    // An axial face is two finer faces stacked on one line of constant z.
    m.axialAreas_.resize(static_cast<std::size_t>(m.nz+1)*m.nr);m.ringMiddles_.resize(m.axialAreas_.size());
    for(int face=0;face<=m.nz;++face) for(int j=0;j<m.nr;++j) {
        const auto k=static_cast<std::size_t>(face)*m.nr+j;
        m.axialAreas_[k]=finer.axialArea(2*face,2*j)+finer.axialArea(2*face,2*j+1);
        m.ringMiddles_[k]=m.radius[face]*0.5*(m.fraction[j]+m.fraction[j+1]);
    }
    // A radial face is the polyline of two finer faces; its middle is their shared node.
    const auto radialCount=static_cast<std::size_t>(m.nz)*(m.nr+1);
    m.radialAreaVectors_.resize(radialCount);m.radialMiddles_.resize(radialCount);
    m.radialRadii_.resize(radialCount);m.radialSecondMoments_.resize(radialCount);
    for(int i=0;i<m.nz;++i) for(int face=0;face<=m.nr;++face) {
        const auto k=m.radialIndex(i,face);
        const auto a=finer.radialAreaVector(2*i,2*face),b=finer.radialAreaVector(2*i+1,2*face);
        m.radialAreaVectors_[k]={a[0]+b[0],a[1]+b[1]};
        m.radialMiddles_[k]=finer.radius[2*i+1]*m.fraction[face];
        m.radialRadii_[k]=face==0?0:m.fraction[face]*(m.columnCubes_[i]/m.columnSquares_[i]);
        m.radialSecondMoments_[k]=0.5*(finer.radialFaceSecondMoment(2*i,2*face)+finer.radialFaceSecondMoment(2*i+1,2*face));
    }
    return m;
}
double Mesh::wallRadiusAt(int i,double z) const {
    // The finest segment of column i holding z (the column's first or last at its ends).
    const int k=std::clamp(i*wallStride+static_cast<int>(std::floor((z-i*dz)/wallDz)),i*wallStride,(i+1)*wallStride-1);
    return wallRadius[k]+(wallRadius[k+1]-wallRadius[k])*(z-k*wallDz)/wallDz;
}
Conserved conservative(Primitive w,double e) {
    return {w.rho,w.rho*w.uz,w.rho*w.ur,w.rho*e+0.5*w.rho*(sq(w.uz)+sq(w.ur))};
}
Thermal thermal(Primitive w,Gas gas) {
    return {w.p/((gas.gamma-1)*w.rho),std::sqrt(gas.gamma*w.p/w.rho),0};
}
Conserved conservative(Primitive w,Gas gas) { return conservative(w,thermal(w,gas).e); }
Primitive primitive(const Conserved& u,Gas gas) {
    if(!admissibleBulk(u,0)) throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
    return {u[0],u[1]/u[0],u[2]/u[0],(gas.gamma-1)*(u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))};
}
Conserved hllc(Primitive l,Thermal tl,Primitive r,Thermal tr,double nz,double nr) {
    const auto ul=conservative(l,tl.e),ur=conservative(r,tr.e);
    const auto fl=physicalFlux(l,tl.e,nz,nr),fr=physicalFlux(r,tr.e,nz,nr);
    double vl=l.uz*nz+l.ur*nr,vr=r.uz*nz+r.ur*nr;
    double al=tl.a,ar=tr.a;
    double sl=std::min(vl-al,vr-ar),sr=std::max(vl+al,vr+ar);
    if(sl>=0) return fl;
    if(sr<=0) return fr;
    auto hlle=[&] {
        Conserved f{};
        for(int k=0;k<4;++k) f[k]=(sr*fl[k]-sl*fr[k]+sl*sr*(ur[k]-ul[k]))/(sr-sl);
        return f;
    };
    double denom=l.rho*(sl-vl)-r.rho*(sr-vr);
    if(std::abs(denom)<1e-30) return hlle();
    double sm=(r.p-l.p+l.rho*vl*(sl-vl)-r.rho*vr*(sr-vr))/denom;
    Primitive w=sm>=0?l:r;
    const auto& u=sm>=0?ul:ur;
    const auto& f=sm>=0?fl:fr;
    double s=sm>=0?sl:sr,v=sm>=0?vl:vr,floor=sm>=0?tl.floor:tr.floor;
    double pstar=w.p+w.rho*(s-v)*(sm-v);
    if(!(pstar>0) || std::abs(s-sm)<1e-20) return hlle();
    double density=w.rho*(s-v)/(s-sm);
    Conserved star{density,density*(w.uz+(sm-v)*nz),density*(w.ur+(sm-v)*nr),
        ((s-v)*u[3]-w.p*v+pstar*sm)/(s-sm)};
    if(!admissibleBulk(star,floor)) return hlle();
    Conserved result{};
    for(int k=0;k<4;++k) result[k]=f[k]+s*(star[k]-u[k]);
    return result;
}
// HLLC written as central flux plus dissipation, with the two acoustic wave terms scaled by phi.
// phi = 1 recovers HLLC exactly; Ma_limit = 0.1 is the published value, not tuned here.
Conserved hllcLm(Primitive l,Thermal tl,Primitive r,Thermal tr,double nz,double nr) {
    const auto ul=conservative(l,tl.e),ur=conservative(r,tr.e);
    const auto fl=physicalFlux(l,tl.e,nz,nr),fr=physicalFlux(r,tr.e,nz,nr);
    double vl=l.uz*nz+l.ur*nr,vr=r.uz*nz+r.ur*nr;
    double al=tl.a,ar=tr.a;
    double sl=std::min(vl-al,vr-ar),sr=std::max(vl+al,vr+ar);
    if(sl>=0) return fl;
    if(sr<=0) return fr;
    double denom=l.rho*(sl-vl)-r.rho*(sr-vr);
    if(std::abs(denom)<1e-30) return hllc(l,tl,r,tr,nz,nr);
    double sm=(r.p-l.p+l.rho*vl*(sl-vl)-r.rho*vr*(sr-vr))/denom;
    auto star=[&](const Primitive& w,const Conserved& u,double s,double v) {
        double pstar=w.p+w.rho*(s-v)*(sm-v),density=w.rho*(s-v)/(s-sm);
        return Conserved{density,density*(w.uz+(sm-v)*nz),density*(w.ur+(sm-v)*nr),
            ((s-v)*u[3]-w.p*v+pstar*sm)/(s-sm)};
    };
    if(std::abs(sl-sm)<1e-20 || std::abs(sr-sm)<1e-20) return hllc(l,tl,r,tr,nz,nr);
    auto sL=star(l,ul,sl,vl),sR=star(r,ur,sr,vr);
    if(!admissibleBulk(sL,tl.floor) || !admissibleBulk(sR,tr.floor)) return hllc(l,tl,r,tr,nz,nr);
    double mach=std::max(std::hypot(l.uz,l.ur)/al,std::hypot(r.uz,r.ur)/ar);
    double phi=std::sin(0.5*pi*std::min(1.0,mach/0.1));
    Conserved result{};
    for(int k=0;k<4;++k)
        result[k]=0.5*(fl[k]+fr[k])+0.5*(phi*sl*(sL[k]-ul[k])+std::abs(sm)*(sL[k]-sR[k])+phi*sr*(sR[k]-ur[k]));
    return result;
}
Conserved hllc(Primitive l,Primitive r,double nz,double nr,Gas gas) {
    return hllc(l,thermal(l,gas),r,thermal(r,gas),nz,nr);
}
Conserved hllcLm(Primitive l,Primitive r,double nz,double nr,Gas gas) {
    return hllcLm(l,thermal(l,gas),r,thermal(r,gas),nz,nr);
}
void thornberScale(Primitive& l,Primitive& r,double al,double ar) {
    double mach=std::max(std::hypot(l.uz,l.ur)/al,std::hypot(r.uz,r.ur)/ar);
    double z=std::min(1.0,mach);
    for(double Primitive::*c:{&Primitive::uz,&Primitive::ur}) {
        double mean=0.5*(l.*c+r.*c),half=0.5*(l.*c-r.*c);
        l.*c=mean+z*half;r.*c=mean-z*half;
    }
}
void thornberScale(Primitive& l,Primitive& r,Gas gas) {
    thornberScale(l,r,thermal(l,gas).a,thermal(r,gas).a);
}
double areaMach(double m,double g) {
    return std::pow(2/(g+1)*(1+(g-1)*m*m/2),(g+1)/(2*(g-1)))/m;
}
double machFromArea(double ratio,bool supersonic,double g) {
    if(ratio<1-1e-12) throw std::invalid_argument("Area ratio is below one.");
    double lo=supersonic?1:1e-8,hi=supersonic?30:1;
    for(int n=0;n<80;++n) {
        double m=(lo+hi)/2;
        if((areaMach(m,g)>ratio)==supersonic) hi=m; else lo=m;
    }
    return (lo+hi)/2;
}
double chokedMassFlow(const Definition& d) {
    double g=d.gas.gamma;
    return pi*sq(d.throatRadius)*d.totalPressure/std::sqrt(d.totalTemperature)*
        std::sqrt(g/d.gas.specificR)*std::pow(2/(g+1),(g+1)/(2*(g-1)));
}
Flow::Flow(Definition d):Flow(d,Mesh(d)) {}
Flow::Flow(Definition d,Mesh mesh):definition_(d),mesh_(std::move(mesh)),medium_(d.medium()),ns_(medium_.size()),nt_(d.turbulence.enabled?2:0),
    nw_(ns_+nt_),inletComposition_(d.massFractions()),totalPressure_(d.totalPressure),corrections_(nt_ && d.turbulence.propertyCorrections) {
    if(mesh_.nz!=d.nz || mesh_.nr!=d.nr) throw std::invalid_argument("The mesh does not have the definition's nz x nr cells.");
    auto count=mesh_.cells.size();
    state_.resize(count); stage_.resize(count); next_.resize(count); rhs_.resize(count);
    slopesZ_.resize(count); primitives_.resize(count); radialLow_.resize(count); radialHigh_.resize(count); pressureSource_.resize(count);
    for(auto* v:{&species_,&speciesStage_,&speciesNext_,&speciesRhs_}) v->resize(count*ns_);
    for(auto* v:{&fractions_,&fractionSlopesZ_,&fractionsLow_,&fractionsHigh_}) v->resize(count*nw_);
    if(d.extremumBoundedEnthalpy)
        for(auto* v:{&cellStaticEnthalpy_,&cellExcess_,&enthalpySlopesZ_,&enthalpyLow_,&enthalpyHigh_}) v->resize(count);
    for(auto* v:{&turbulence_,&turbulenceStage_,&turbulenceNext_,&turbulenceRhs_,&turbulenceStart_}) v->resize(count*nt_);
    temperature_.assign(count,d.totalTemperature); sound_.resize(count);
    axialFlux_.resize(static_cast<std::size_t>(d.nz+1)*d.nr);axialTransported_.resize(axialFlux_.size()*nw_);
    radialFlux_.resize(static_cast<std::size_t>(d.nz)*(d.nr+1));radialTransported_.resize(radialFlux_.size()*nw_);
    axialArea_.resize(axialFlux_.size());radialArea_.resize(radialFlux_.size());
    for(int i=0;i<=d.nz;++i) for(int j=0;j<d.nr;++j) axialArea_[static_cast<std::size_t>(i)*d.nr+j]=mesh_.axialArea(i,j);
    for(int i=0;i<d.nz;++i) for(int j=1;j<=d.nr;++j) {
        auto ar=mesh_.radialAreaVector(i,j);radialArea_[static_cast<std::size_t>(i)*(d.nr+1)+j]=std::hypot(ar[0],ar[1]);
    }
    igniterCell_.assign(count,0);igniterCellVolume_.assign(count,0.0);clippedCell_.assign(count,0.0);
    setThreads(1);
    // Initial fill: the reservoir composition, prepared with its frozen properties at the total temperature.
    const double* y0=inletComposition_.data();
    const double gasR=medium_.gasConstant(y0),cv0=medium_.cv(d.totalTemperature,y0),g=(cv0+gasR)/cv0;
    std::vector<Primitive> cells(count);
    for(int i=0;i<d.nz;++i) {
        Primitive base{d.totalPressure/(gasR*d.totalTemperature),0,0,d.totalPressure};
        double radius=(mesh_.radius[i]+mesh_.radius[i+1])/2;
        if(d.experiment==Case::ShockTube)
            base=i<d.nz/2?Primitive{1,0,0,1}:Primitive{0.125,0,0,0.1};
        if(d.experiment==Case::Chamber) base={d.backPressure/(gasR*d.ambientTemperature),0,0,d.backPressure};
        if(d.experiment==Case::Nozzle) {
            // The 1D preparation is shared by a whole axial column.
            double m=machFromArea(sq(radius/d.throatRadius),(i+0.5)*mesh_.dz>d.length*d.throatFraction,g);
            double factor=1+(g-1)*m*m/2;
            double t=d.totalTemperature/factor;
            base.p=d.totalPressure/std::pow(factor,g/(g-1));
            base.rho=base.p/(gasR*t);
            base.uz=m*std::sqrt(g*gasR*t);
        }
        for(int j=0;j<d.nr;++j) {
            auto q=mesh_.index(i,j);auto w=base;
            if(d.experiment==Case::Nozzle)
                w.ur=w.uz*mesh_.cells[q].r/radius*(mesh_.radius[i+1]-mesh_.radius[i])/mesh_.dz;
            cells[q]=w;
        }
    }
    if(d.experiment==Case::Chamber) {
        ambient_=inletComposition_;
        if(nt_) { ambient_.push_back(d.turbulence.ambientK);ambient_.push_back(d.turbulence.ambientOmega); }
        faceSupply_.assign(d.nr,-1);supplyArea_.assign(d.supplies.size(),0);
        if(d.openAreaSupplies) openAreaShares();
        else {
            for(int j=0;j<d.nr;++j) {
                double centre=mesh_.ringMiddle(0,j);
                for(std::size_t s=0;s<d.supplies.size();++s)
                    if(centre>=d.supplies[s].innerRadius && centre<d.supplies[s].outerRadius) {
                        if(faceSupply_[j]>=0) throw std::invalid_argument("Supplies overlap at the injector face.");
                        faceSupply_[j]=static_cast<int>(s);supplyArea_[s]+=mesh_.axialArea(0,j);
                    }
            }
            for(double a:supplyArea_)
                if(!(a>0)) throw std::invalid_argument("A supply contains no injector-face ring centre; widen it or refine the radial mesh.");
        }
        const auto& g=d.igniter;
        if(g.energy>0) {
            for(std::size_t q=0;q<count;++q) {
                const auto& c=mesh_.cells[q];
                if(c.z>=g.zMin && c.z<=g.zMax && c.r<=g.rMax) {
                    igniterCells_.push_back(q);igniterCell_[q]=1;igniterCellVolume_[q]=c.volume;igniterVolume_+=c.volume;
                }
            }
            if(igniterCells_.empty()) throw std::invalid_argument("The igniter region contains no cell centroid.");
        }
    }
    if(nt_) {
        std::vector<bool> plate;
        if(d.experiment==Case::Chamber && !d.wallSlip) for(int s:faceSupply_) plate.push_back(s<0);
        wallDistance_=wallDistance(mesh_,!d.wallSlip,plate,corrections_?&wallDirection_:nullptr);
    }
    if(medium_.hasTransport()) prepareTransport();
    std::vector<double> fractions;
    for(std::size_t q=0;q<count;++q) fractions.insert(fractions.end(),inletComposition_.begin(),inletComposition_.end());
    setInitialState(cells,fractions);
}
// Each face ring's shares of the supplies' true annuli (docs/evidence/FACE_INLET.md): a_sj the overlap
// of supply s with ring j (the ring's own area when the supply covers it whole, so an aligned layout
// keeps the ring areas bit for bit), phi_sj = a_sj / the ring's area, the supply's area the sum of its
// a_sj, the wall share w_j = 1 - sum phi (0 below 1e-12), and the ring's kind: the plate if w_j >= 1/2,
// else its largest supply (the lower index on a tie).
void Flow::openAreaShares() {
    const auto& d=definition_;const auto& m=mesh_;
    const double face=m.radius[0];
    std::vector<std::size_t> order(d.supplies.size());
    for(std::size_t s=0;s<order.size();++s) order[s]=s;
    std::sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b) { return d.supplies[a].innerRadius<d.supplies[b].innerRadius; });
    for(std::size_t n=1;n<order.size();++n)
        if(d.supplies[order[n]].innerRadius<d.supplies[order[n-1]].outerRadius) throw std::invalid_argument("Supplies overlap at the injector face.");
    for(const auto& s:d.supplies)
        if(s.outerRadius>face*(1+1e-12)) throw std::invalid_argument("A supply reaches past the injector face's wall radius.");
    faceShares_.clear();shareBegin_.assign(static_cast<std::size_t>(m.nr)+1,0);
    plateShare_.assign(m.nr,1.0);plateMomentum_.assign(m.nr,0.0);
    for(int j=0;j<m.nr;++j) {
        shareBegin_[j]=faceShares_.size();
        const double lo=face*m.fraction[j],hi=face*m.fraction[j+1];
        const double area=m.axialArea(0,j);
        double covered=0,largest=0;
        for(std::size_t s=0;s<d.supplies.size();++s) {
            const auto& u=d.supplies[s];
            double a=u.innerRadius<=lo && u.outerRadius>=hi?area:pi*(sq(std::min(u.outerRadius,hi))-sq(std::max(u.innerRadius,lo)));
            if(!(a>0)) continue;
            const double phi=a/area;
            faceShares_.push_back({j,static_cast<int>(s),phi});supplyArea_[s]+=a;covered+=phi;
            if(phi>largest) { largest=phi;faceSupply_[j]=static_cast<int>(s); }
        }
        double w=1-covered;
        if(w<1e-12) w=0;
        plateShare_[j]=w;
        if(w>=0.5) faceSupply_[j]=-1;
    }
    shareBegin_[m.nr]=faceShares_.size();
    for(double a:supplyArea_)
        if(!(a>0)) throw std::invalid_argument("A supply has no area on the injector face.");
}
void Flow::setThreads(int threads) {
    pool_=std::make_unique<Pool>(threads);
    const auto n=static_cast<std::size_t>(pool_->threads());
    faceFractions_.assign(2*n,std::vector<double>(nw_));
    transportScratch_.assign(n,TransportScratch{{},std::vector<double>(ns_),std::vector<double>(ns_),std::vector<double>(ns_),{},std::vector<double>(ns_)});
    workerOk_.assign(n,1);workerDt_.assign(n,0.0);workerLimit_.assign(n,{});
}
void Flow::setBudgetProbe(std::vector<std::size_t> cells,std::vector<double> weights) {
    if(!active_.empty()) throw std::logic_error("Flow::setBudgetProbe: single level only.");
    if(definition_.openAreaSupplies) throw std::logic_error("Flow::setBudgetProbe: the probe books the ring-centre supplies' faces only.");
    if(weights.size()!=ns_) throw std::invalid_argument("Flow::setBudgetProbe: one weight per species.");
    for(std::size_t q:cells) if(q>=state_.size()) throw std::out_of_range("Flow::setBudgetProbe: no such cell.");
    probe_.cells=std::move(cells);probe_.weights=std::move(weights);
    probe_.sums.assign(probe_.cells.size()*kProbeStride,0.0);
    probeStage_.assign(probe_.sums.size(),0.0);probeAttempt_.assign(probe_.sums.size(),0.0);
}
// The probe cells' convective face terms and sources of the current right-hand side, with gather's
// faces, signs and supply composition.
void Flow::probeConvective() {
    const auto& m=mesh_;const auto& d=definition_;
    std::fill(probeStage_.begin(),probeStage_.end(),0.0);
    auto tracked=[&](const double* t) { double s=0;for(std::size_t k=0;k<ns_;++k) s+=probe_.weights[k]*t[k];return s; };
    for(std::size_t c=0;c<probe_.cells.size();++c) {
        const std::size_t q=probe_.cells[c];
        const int i=static_cast<int>(q/m.nr),j=static_cast<int>(q%m.nr);
        const std::size_t a0=q,a1=q+m.nr,r0=static_cast<std::size_t>(i)*(m.nr+1)+j,r1=r0+1;
        double* out=probeStage_.data()+c*kProbeStride;
        auto put=[&](int face,double sign,double area,const Conserved& flux,double moved) {
            double* o=out+face*kFaceTerms;
            o[kInviscidMass]=sign*area*flux[0];o[kInviscidEnergy]=sign*area*flux[3];o[kInviscidTracked]=sign*moved;
        };
        const bool injector=i==0 && d.experiment==Case::Chamber;
        const int supply=injector?faceSupply_[j]:-1;
        double moved=0;
        if(supply>=0) for(std::size_t k=0;k<ns_;++k) moved+=probe_.weights[k]*axialArea_[a0]*axialFlux_[a0][0]*d.supplies[supply].composition[k];
        else if(!injector) moved=tracked(&axialTransported_[a0*nw_]);
        put(0,1,axialArea_[a0],axialFlux_[a0],moved);
        put(1,-1,axialArea_[a1],axialFlux_[a1],tracked(&axialTransported_[a1*nw_]));
        if(j>0) put(2,1,radialArea_[r0],radialFlux_[r0],tracked(&radialTransported_[r0*nw_]));
        put(3,-1,radialArea_[r1],radialFlux_[r1],tracked(&radialTransported_[r1*nw_]));
        double* source=out+4*kFaceTerms;
        if(igniterPower_>0 && igniterCell_[q]) source[kIgniter]=igniterPower_*igniterCellVolume_[q]/igniterVolume_;
        if(!heating_.empty()) source[kHeating]=heating_[q]*m.cells[q].volume;
        if(!bodyForce_.empty()) {
            const auto& f=bodyForce_[q];
            source[kBodyWork]=(f[0]*primitives_[q].uz+f[1]*primitives_[q].ur)*m.cells[q].volume;
        }
    }
}
template<class F> bool Flow::allOf(const std::vector<std::size_t>* list,F f) {
    std::fill(workerOk_.begin(),workerOk_.end(),char{1});
    eachCell(list,[&](std::size_t q,int worker) { if(!f(q)) workerOk_[worker]=0; });
    return std::find(workerOk_.begin(),workerOk_.end(),char{0})==workerOk_.end();
}
// AMR step B2: the regions of the active cells for ghost width `width` (docs/evidence/AMR.md, Step B2:
// setup). The driver sets filled_.
void Flow::setRegions(int width) {
    const auto& m=mesh_;
    const std::size_t count=m.cells.size();
    distance_.assign(count,static_cast<unsigned char>(width+1));
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
        if(!active_[m.index(i,j)]) continue;
        for(int ii=std::max(0,i-width);ii<=std::min(m.nz-1,i+width);++ii)
            for(int jj=std::max(0,j-width);jj<=std::min(m.nr-1,j+width);++jj) {
                auto& d=distance_[m.index(ii,jj)];
                d=std::min(d,static_cast<unsigned char>(std::max(std::abs(ii-i),std::abs(jj-j))));
            }
    }
    within_.assign(static_cast<std::size_t>(width+1),{});
    for(std::size_t q=0;q<count;++q) for(int k=distance_[q];k<=width;++k) within_[static_cast<std::size_t>(k)].push_back(q);
    faces_.clear();
    const std::size_t na=axialFlux_.size();
    for(std::size_t f=0;f<na;++f) {
        const int i=static_cast<int>(f/m.nr),j=static_cast<int>(f%m.nr);
        if((i>0 && active_[m.index(i-1,j)]) || (i<m.nz && active_[m.index(i,j)])) faces_.push_back(f);
    }
    axialFaces_=faces_.size();
    for(std::size_t e=0;e<radialFlux_.size();++e) {
        const int i=static_cast<int>(e/(m.nr+1)),j=static_cast<int>(e%(m.nr+1));
        if(j>0 && (active_[m.index(i,j-1)] || (j<m.nr && active_[m.index(i,j)]))) faces_.push_back(na+e);
    }
    columns_.clear();
    for(int i=0;i<m.nz;++i) {
        Column c{i,{}};
        for(int j=0;j<m.nr;) {
            if(distance_[m.index(i,j)]>1) { ++j;continue; }
            const int lo=j;
            while(j<m.nr && distance_[m.index(i,j)]<=1) ++j;
            c.rows.push_back({lo,j});
        }
        if(!c.rows.empty()) columns_.push_back(std::move(c));
    }
}
void Flow::clearRegions() {
    within_.clear();distance_.clear();filled_.clear();faces_.clear();columns_.clear();axialFaces_=0;
}
double Flow::BoundaryRates::grossMomentum() const {
    return std::abs(inletMomentum)+std::abs(outletMomentum)+std::abs(wallAxial)+std::abs(bodyAxial);
}
Flow::BoundaryRates Flow::BoundaryRates::average(const BoundaryRates& a,const BoundaryRates& b) {
    auto mean=[](double x,double y){return 0.5*(x+y);};
    return {mean(a.mass,b.mass),mean(a.energy,b.energy),mean(a.inlet,b.inlet),mean(a.outlet,b.outlet),
        mean(a.inletMomentum,b.inletMomentum),mean(a.outletMomentum,b.outletMomentum),
        mean(a.wallAxial,b.wallAxial),mean(a.bodyAxial,b.bodyAxial),mean(a.heat,b.heat),mean(a.wallHeat,b.wallHeat),
        mean(a.correction,b.correction),mean(a.plateHeat,b.plateHeat)};
}
void Flow::resetAccounting() {
    initialMass_=initialEnergy_=initialMomentum_=0;
    for(std::size_t q=0;q<state_.size();++q) {
        initialMass_+=state_[q][0]*mesh_.cells[q].volume;
        initialMomentum_+=state_[q][1]*mesh_.cells[q].volume;
        initialEnergy_+=state_[q][3]*mesh_.cells[q].volume;
    }
    integratedMassFlux_=integratedEnergyFlux_=integratedMomentumSource_=integratedMomentumGross_=integratedHeat_=heatCompensation_=clippedTurbulentEnergy_=prescribedTurbulentEnergy_=integratedCorrection_=0;
    massCompensation_=energyCompensation_=correctionCompensation_=momentumCompensation_=0;
    exchangedMass_=exchangedEnergy_=exchangedMomentum_=0;
    lastRates_={};time_=dt_=0;steps_=rejectedSteps_=0;
}
void Flow::setUniform(Primitive w) {
    setInitialState(std::vector<Primitive>(state_.size(),w));
}
void Flow::setInitialState(const std::vector<Primitive>& cells) {
    std::vector<double> fractions;
    for(std::size_t q=0;q<cells.size();++q) fractions.insert(fractions.end(),inletComposition_.begin(),inletComposition_.end());
    setInitialState(cells,fractions);
}
void Flow::setInitialState(const std::vector<Primitive>& cells,const std::vector<double>& y) {
    if(cells.size()!=state_.size() || y.size()!=cells.size()*ns_) throw std::invalid_argument("Initial field does not match the mesh.");
    std::vector<Conserved> initialized(cells.size());
    std::vector<double> partial(y.size()),turbulence(cells.size()*nt_);
    const auto& tu=definition_.turbulence;
    for(std::size_t q=0;q<cells.size();++q) {
        const auto& w=cells[q];const double* yq=y.data()+q*ns_;
        if(!(w.rho>0) || !(w.p>0)) throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
        double t=w.p/(w.rho*medium_.gasConstant(yq));
        initialized[q]=conservative(w,medium_.internalEnergy(t,yq));
        temperature_[q]=t;
        for(std::size_t k=0;k<ns_;++k) partial[q*ns_+k]=w.rho*yq[k];
        // The fill carries the ambient turbulence, its k on top of the thermal energy.
        if(nt_) { turbulence[q*nt_]=w.rho*tu.ambientK;turbulence[q*nt_+1]=w.rho*tu.ambientOmega;initialized[q][3]+=w.rho*tu.ambientK; }
        if(!admissible(initialized[q],partial.data()+q*ns_,turbulence.data()+q*nt_))
            throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
    }
    state_.swap(initialized);species_.swap(partial);turbulence_.swap(turbulence);resetAccounting();
}
void Flow::setMassFractions(std::size_t q,const double* y) {
    if(q>=state_.size()) throw std::out_of_range("Cell index out of range.");
    std::vector<double> partial(ns_);
    for(std::size_t k=0;k<ns_;++k) partial[k]=state_[q][0]*y[k];
    if(!admissible(state_[q],partial.data(),turbulence_.data()+q*nt_)) throw std::runtime_error("Composition gives a non-admissible state.");
    std::copy(partial.begin(),partial.end(),species_.begin()+static_cast<std::ptrdiff_t>(q*ns_));
}
void Flow::setPartialDensities(const std::vector<double>& partial) {
    if(partial.size()!=species_.size()) throw std::invalid_argument("Partial densities do not match the mesh.");
    for(std::size_t q=0;q<state_.size();++q)
        if(!admissible(state_[q],partial.data()+q*ns_,turbulence_.data()+q*nt_)) throw std::runtime_error("Partial densities give a non-admissible state.");
    species_=partial;
}
void Flow::setTurbulence(const std::vector<double>& kOmega) {
    if(!nt_) throw std::logic_error("The definition has no turbulence.");
    if(kOmega.size()!=turbulence_.size()) throw std::invalid_argument("Turbulence field does not match the mesh.");
    auto state=state_;auto turbulence=turbulence_;
    for(std::size_t q=0;q<state.size();++q) {
        double rho=state[q][0],k=kOmega[q*nt_],omega=kOmega[q*nt_+1];
        if(!(k>=0) || !(omega>0) || !std::isfinite(k) || !std::isfinite(omega))
            throw std::invalid_argument("Turbulence needs finite k >= 0 and omega > 0.");
        state[q][3]+=rho*k-turbulence[q*nt_];
        turbulence[q*nt_]=rho*k;turbulence[q*nt_+1]=rho*omega;
        if(!admissible(state[q],species_.data()+q*ns_,turbulence.data()+q*nt_))
            throw std::runtime_error("Turbulence gives a non-admissible state.");
    }
    state_.swap(state);turbulence_.swap(turbulence);resetAccounting();
}
void Flow::setBodyForce(std::vector<std::array<double,2>> force) {
    if(!force.empty() && force.size()!=state_.size()) throw std::invalid_argument("Body force does not match the mesh.");
    for(const auto& f:force) if(!std::isfinite(f[0]) || !std::isfinite(f[1])) throw std::invalid_argument("Body force must be finite.");
    bodyForce_=std::move(force);
}
void Flow::setHeating(std::vector<double> power) {
    if(!power.empty() && power.size()!=state_.size()) throw std::invalid_argument("Heating does not match the mesh.");
    for(double h:power) if(!std::isfinite(h)) throw std::invalid_argument("Heating must be finite.");
    heating_=std::move(power);
}
void Flow::setTotalPressure(double p) {
    if(!(p>0) || !std::isfinite(p)) throw std::invalid_argument("Reservoir pressure must be positive and finite.");
    totalPressure_=p;
}
void Flow::exchange(const std::vector<double>& amounts) {
    const std::size_t stride=ns_+3;
    if(amounts.size()!=state_.size()*stride) throw std::invalid_argument("Flow::exchange: amounts must hold cells x (species + 3).");
    for(std::size_t q=0;q<state_.size();++q) {
        const double* a=amounts.data()+q*stride;
        if(std::all_of(a,a+stride,[](double x) { return x==0; })) continue;
        const double v=mesh_.cells[q].volume;
        double mass=0;
        for(std::size_t k=0;k<ns_;++k) { species_[q*ns_+k]+=a[k]/v;mass+=a[k]; }
        state_[q][0]+=mass/v;state_[q][1]+=a[ns_]/v;state_[q][2]+=a[ns_+1]/v;state_[q][3]+=a[ns_+2]/v;
        exchangedMass_+=mass;exchangedMomentum_+=a[ns_];exchangedEnergy_+=a[ns_+2];
        integratedMomentumGross_+=std::abs(a[ns_]);
        for(std::size_t c=0;c<probe_.cells.size();++c) {
            if(probe_.cells[c]!=q) continue;
            double* source=probe_.sums.data()+c*kProbeStride+4*kFaceTerms;
            source[kExchangedMass]+=mass;source[kExchangedEnergy]+=a[ns_+2];
            for(std::size_t k=0;k<ns_;++k) source[kExchangedTracked]+=probe_.weights[k]*a[k];
        }
        if(!admissible(state_[q],species_.data()+q*ns_,nt_?turbulence_.data()+q*nt_:nullptr))
            throw std::runtime_error("Flow::exchange left a cell non-admissible.");
    }
}
bool Flow::admissible(const Conserved& u,const double* partial,const double* turbulence) const {
    // energyFloor is linear in the mass fractions, so it is evaluated on the clipped partial densities.
    double sum=0;
    for(std::size_t k=0;k<ns_;++k) { if(!std::isfinite(partial[k])) return false; sum+=std::max(partial[k],0.0); }
    if(!(sum>0)) return false;
    double k=0;
    if(nt_) {
        // rho k >= 0 and rho omega > 0; the thermal energy lies above the floor once k is taken out.
        if(!(u[0]>0) || !(turbulence[0]>=0) || !(turbulence[1]>0) || !std::isfinite(turbulence[0]) || !std::isfinite(turbulence[1])) return false;
        k=turbulence[0]/u[0];
    }
    return admissibleBulk(u,medium_.energyFloorOfClipped(partial)/sum+k);
}
std::vector<double> Flow::massFractions(std::size_t q) const {
    std::vector<double> y(ns_);
    normalise(species_.data()+q*ns_,ns_,y.data());
    return y;
}
double Flow::temperature(std::size_t q) const {
    auto y=massFractions(q);const auto& u=state_[q];
    return medium_.temperature((u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))/u[0]-turbulentEnergy(u,turbulence_,q),y.data(),temperature_[q]);
}
Primitive Flow::cellPrimitive(std::size_t q) const {
    auto y=massFractions(q);const auto& u=state_[q];const double k=turbulentEnergy(u,turbulence_,q);
    double t=medium_.temperature((u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))/u[0]-k,y.data(),temperature_[q]);
    if(!std::isfinite(t) || !admissibleBulk(u,medium_.energyFloor(y.data())+k))
        throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
    return {u[0],u[1]/u[0],u[2]/u[0],u[0]*medium_.gasConstant(y.data())*t};
}
// Cell primitives, mass fractions (then specific k and omega), temperatures and frozen sound speeds of a state.
void Flow::refresh(const std::vector<Conserved>& state,const std::vector<double>& species,const std::vector<double>& turbulence,int depth) {
    eachCell(region(depth),[&](std::size_t q,int) {
        const auto& u=state[q];double* y=fractions_.data()+q*nw_;
        normalise(species.data()+q*ns_,ns_,y);
        double k=0;
        if(nt_) { k=y[ns_]=turbulence[q*nt_]/u[0];y[ns_+1]=turbulence[q*nt_+1]/u[0]; }
        if(!admissibleBulk(u,medium_.energyFloor(y)+k))
            throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
        double t=medium_.temperature((u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))/u[0]-k,y,temperature_[q]);
        auto props=medium_.properties(t,y);
        temperature_[q]=t;sound_[q]=std::sqrt((props.cv+props.r)/props.cv*props.r*t);
        primitives_[q]={u[0],u[1]/u[0],u[2]/u[0],u[0]*props.r*t};
        if(definition_.extremumBoundedEnthalpy) {
            cellStaticEnthalpy_[q]=props.e+props.r*t;
            cellExcess_[q]=boundMeasure(primitives_[q],props.e+k,y);
        }
    });
}
Conserved Flow::reconstructionValues(std::size_t q) const {
    auto w=values(primitives_[q]);
    if(definition_.temperatureReconstruction) w[0]=temperature_[q];
    return w;
}
void Flow::toFace(Primitive& w,double* y) const {
    if(!definition_.temperatureReconstruction) return;
    normalise(y,ns_,y);
    w.rho=w.p/(medium_.gasConstant(y)*w.rho);
}
double Flow::boundMeasure(const Primitive& w,double e,const double* y) const {
    if(boundedExcess_.empty()) return e+w.p/w.rho;
    double b=e-(nt_?y[ns_]:0.0)+w.p/w.rho+0.5*(w.uz*w.uz+w.ur*w.ur);
    for(std::size_t k=0;k<ns_;++k) b-=boundedExcess_[k]*y[k];
    return b;
}
double Flow::cellEnthalpy(std::size_t q) const {
    const double* y=fractions_.data()+q*nw_;
    return boundMeasure(primitives_[q],medium_.properties(temperature_[q],y).e+(nt_?y[ns_]:0.0),y);
}
void Flow::boundEnthalpy(Primitive& w,Thermal& t,double* y,std::size_t q,std::size_t other) const {
    if(!definition_.boundedFaceEnthalpy) return;
    const double h=boundMeasure(w,t.e,y),a=cellEnthalpy(q),b=cellEnthalpy(other);
    if(h>=std::min(a,b) && h<=std::max(a,b)) return;
    w=primitives_[q];
    std::copy_n(fractions_.data()+q*nw_,nw_,y);
    t=faceThermal(w,y);
}
bool Flow::excessExtremum(std::size_t q) const {
    const auto& m=mesh_;
    const int i=static_cast<int>(q/m.nr),j=static_cast<int>(q%m.nr);
    const double b=cellExcess_[q];
    bool above=true,below=true;
    auto compare=[&](int a,int c){const double n=cellExcess_[m.index(a,c)];above=above && b>n;below=below && b<n;};
    if(i>0) compare(i-1,j);
    if(i<m.nz-1) compare(i+1,j);
    if(j>0) compare(i,j-1);
    if(j<m.nr-1) compare(i,j+1);
    return above || below;
}
// Newton on h(T) = e(T, y) + R T with the bracket of Medium::temperature, from the reconstructed T.
void Flow::enthalpyFace(Primitive& w,double h,double* y) const {
    normalise(y,ns_,y);
    const double gasR=medium_.gasConstant(y);
    double lo=0,hi=std::numeric_limits<double>::infinity(),t=w.p/(w.rho*gasR);
    if(!(t>0) || !std::isfinite(t)) t=300;
    for(int n=0;n<100;++n) {
        auto props=medium_.properties(t,y);
        const double f=props.e+gasR*t-h;
        double next=t-f/(props.cv+gasR);
        if(std::abs(next-t)<=1e-13*t) { if(next>0) w.rho=w.p/(gasR*next); return; }
        if(f>0) hi=t; else lo=t;
        if(!(next>lo && next<hi)) next=std::isfinite(hi)?0.5*(lo+hi):2*t;
        t=next;
    }
}
// Face thermodynamics from reconstructed (rho, p), mass fractions and (with turbulence) k; y is
// clipped and normalised in place. The face energy includes k, so the flux carries rho k in E.
Thermal Flow::faceThermal(const Primitive& w,double* y) const {
    normalise(y,ns_,y);
    double t=w.p/(w.rho*medium_.gasConstant(y));
    auto props=medium_.properties(t,y);
    const double k=nt_?y[ns_]:0.0;
    return {props.e+k,std::sqrt((props.cv+props.r)/props.cv*props.r*t),medium_.energyFloor(y)+k};
}
// Subsonic reservoir boundary at the reservoir composition. Exact for a calorically perfect gas;
// for a thermally perfect mixture it uses the frozen ratio of specific heats at the total temperature.
Primitive Flow::inlet(Primitive w) const {
    if(definition_.experiment!=Case::Nozzle) return w;
    const double* y=inletComposition_.data();
    const double gasR=medium_.gasConstant(y),cv0=medium_.cv(definition_.totalTemperature,y);
    const double g=(cv0+gasR)/cv0,a0=std::sqrt(g*gasR*definition_.totalTemperature);
    const double invariant=w.uz-2*std::sqrt(g*w.p/w.rho)/(g-1);
    // Preserve the outgoing acoustic characteristic.
    double lo=0,hi=a0/std::sqrt(1+(g-1)/2);
    const double lowInvariant=-2*a0/(g-1);
    const double highInvariant=hi-2*hi/(g-1);
    if(invariant<lowInvariant-1e-10*a0 || invariant>highInvariant+1e-10*a0)
        throw std::runtime_error("Inlet requires reverse or supersonic flow, outside this reservoir boundary model.");
    for(int n=0;n<40;++n) {
        double v=(lo+hi)/2,a=std::sqrt(a0*a0-(g-1)*v*v/2);
        if(v-2*a/(g-1)>invariant) hi=v; else lo=v;
    }
    double v=(lo+hi)/2,t=definition_.totalTemperature-v*v*(g-1)/(2*g*gasR);
    double p=totalPressure_*std::pow(t/definition_.totalTemperature,g/(g-1));
    return {p/(gasR*t),v,0,p};
}
// Outlet with the interior face's frozen ratio of specific heats (exact for a calorically perfect gas).
// A Chamber exit in subsonic backflow draws ambient gas at rest: ambient pressure, temperature and
// composition, with the interior's axial velocity (a declared simplification for start-up transients).
Primitive Flow::outlet(Primitive w,const double* y,const Primitive& cell,bool& ambientInflow) const {
    ambientInflow=false;
    const bool chamber=definition_.experiment==Case::Chamber;
    if(definition_.experiment!=Case::Nozzle && !chamber) return w;
    double t=w.p/(w.rho*medium_.gasConstant(y)),a=medium_.soundSpeed(t,y),g=a*a*w.rho/w.p;
    auto ambient=[&](double u) {
        ambientInflow=true;double p=definition_.backPressure;
        return Primitive{p/(medium_.gasConstant(ambient_.data())*definition_.ambientTemperature),u,0,p};
    };
    if(w.uz>=a) return w; // All characteristics leave a supersonic outlet.
    // Non-reflecting option, in the acoustic characteristic variables W = p +- rho a u_z with the
    // face's impedance. The outgoing W+ is the face's; the incoming W- moves the fraction beta from
    // the last cell's own value toward the one that gives backPressure. The cell lags the face by
    // half a cell (dz / 2a), so the face's incoming wave relaxes at 2 beta a / dz, which is K / 2
    // in the Poinsot-Lele form for the beta below (measured: 3.9 beta a / dz on 800 and 1600
    // cells). Measured failures of the alternatives: taking W- from the extrapolated face feeds the
    // interior slope back every step and grows without bound; isentropic Riemann invariants
    // u -+ 2a/(gamma-1) read a temperature difference between face and cell as an acoustic wave,
    // and a 1 K hot spot leaving the duct grows into a blow-up. Pressure and velocity follow the
    // same rule in backflow, where the gas drawn in takes only the ambient entropy and
    // composition, so the outlet does not change its acoustic behaviour when u_z crosses zero.
    if(definition_.outletRelaxation>0) {
        double mach=w.uz/a,beta=std::min(1.0,0.25*definition_.outletRelaxation*(1-mach*mach)*mesh_.dz/definition_.span());
        double impedance=w.rho*a,outgoing=w.p+impedance*w.uz,incoming=cell.p-impedance*cell.uz;
        incoming+=beta*(2*definition_.backPressure-outgoing-incoming);
        double p=0.5*(outgoing+incoming),u=0.5*(outgoing-incoming)/impedance;
        if(!(p>0)) throw std::runtime_error("Non-reflecting outlet produced a non-positive pressure.");
        if(u < -1e-10*a) {
            if(!chamber) throw std::runtime_error("Outlet backflow is outside this nozzle prototype's boundary model.");
            ambientInflow=true;
            return {p/(medium_.gasConstant(ambient_.data())*definition_.ambientTemperature),std::max(u,-a),0,p};
        }
        return {w.rho*std::pow(p/w.p,1/g),u,w.ur,p};
    }
    if(w.uz < -1e-10*a) {
        if(chamber) return ambient(std::max(w.uz,-a));
        throw std::runtime_error("Outlet backflow is outside this nozzle prototype's boundary model.");
    }
    double p=definition_.backPressure,rho=w.rho*std::pow(p/w.p,1/g);
    double u=w.uz+2*(a-std::sqrt(g*p/rho))/(g-1);
    if(u < -1e-10*a) {
        if(chamber) return ambient(std::max(u,-a));
        throw std::runtime_error("Imposed outlet pressure requires unsupported backflow.");
    }
    return {rho,u,w.ur,p};
}
double Flow::nextEvent(double time) const {
    double next=std::numeric_limits<double>::infinity();
    auto consider=[&](double t){ if(t>time) next=std::min(next,t); };
    if(definition_.experiment!=Case::Chamber) return next;
    for(const auto& s:definition_.supplies) { consider(s.opens);consider(s.opens+s.ramp); }
    if(definition_.igniter.energy>0) { consider(definition_.igniter.start);consider(definition_.igniter.start+definition_.igniter.duration); }
    return next;
}
// Supply face: mass flux g, total enthalpy h0 and composition are imposed. The face pressure p and
// velocity u = g/rho(p, T) satisfy the outgoing (left-running) characteristic from the interior
// state, p = p_in (1 + (gamma-1)/(2 a_in) (u - u_in))^(2 gamma/(gamma-1)), and T solves
// h(T) + c u^2 = h0, c = 1/2 without turbulence. If that needs u above the face sound speed the face
// is choked and carries the sonic state of the same stream. g = 0 is a closed valve (a wall).
// With turbulence k = 3/2 (I u)^2, so h + u^2/2 + 5/3 k = h0 gives c = (1 + 5 I^2) / 2, and the
// momentum flux g u + p + 2/3 rho k is g u (1 + I^2) + p.
SupplyFace Flow::supplyFace(const Supply& s,double g,const Primitive& in,const double* yIn) const {
    const double* y=s.composition.data();
    const double gasR=medium_.gasConstant(y),h0=medium_.enthalpy(s.totalTemperature,y);
    const double tIn=in.p/(in.rho*medium_.gasConstant(yIn)),aIn=medium_.soundSpeed(tIn,yIn),gIn=aIn*aIn*in.rho/in.p;
    const double i2=nt_?sq(s.turbulenceIntensity):0.0,c=0.5*(1+5*i2);
    auto characteristic=[&](double u) {
        double base=1+(gIn-1)/(2*aIn)*(u-in.uz);
        if(!(base>0)) throw std::runtime_error("Supply face: the interior gas recedes faster than the supply can follow.");
        return in.p*std::pow(base,2*gIn/(gIn-1));
    };
    SupplyFace face;
    if(!(g>0)) { face.p=characteristic(0);face.flux={0,face.p,0,0};return face; }
    // Static temperature of the stream at face pressure p (Newton from the stagnation temperature).
    auto staticTemperature=[&](double p) {
        double t=s.totalTemperature;
        for(int n=0;n<60;++n) {
            double u=g*gasR*t/p,f=medium_.enthalpy(t,y)+c*u*u-h0;
            double step=f/(medium_.cv(t,y)+gasR+2*c*u*u/t);
            t=std::max(0.5*t,t-step);
            if(std::abs(step)<1e-13*t) break;
        }
        return t;
    };
    // Sonic state: h(T) + c a(T)^2 = h0.
    double tSonic=s.totalTemperature;
    for(int n=0;n<60;++n) {
        double cv=medium_.cv(tSonic,y),gamma=(cv+gasR)/cv;
        double f=medium_.enthalpy(tSonic,y)+c*gamma*gasR*tSonic-h0;
        double step=f/(cv+gasR+c*gamma*gasR);
        tSonic=std::max(0.5*tSonic,tSonic-step);
        if(std::abs(step)<1e-13*tSonic) break;
    }
    const double aSonic=medium_.soundSpeed(tSonic,y),pSonic=g*gasR*tSonic/aSonic;
    auto residual=[&](double p){ double t=staticTemperature(p);return p-characteristic(g*gasR*t/p); };
    double p=pSonic,t=tSonic,u=aSonic;
    face.choked=!(residual(pSonic)<0);
    if(!face.choked) {
        // Subsonic face: the residual increases with p; bracket and bisect.
        double lo=pSonic,hi=2*std::max(pSonic,in.p);
        for(int n=0;residual(hi)<0;++n) { if(n==60) throw std::runtime_error("Supply face pressure not bracketed."); lo=hi;hi*=2; }
        for(int n=0;n<200 && hi-lo>1e-14*hi;++n) { double mid=0.5*(lo+hi);(residual(mid)<0?lo:hi)=mid; }
        p=0.5*(lo+hi);t=staticTemperature(p);u=g*gasR*t/p;
    }
    face.p=p;face.t=t;face.u=u;
    face.flux={g,g*u*(1+i2)+p,0,g*h0};
    if(nt_) {
        face.k=1.5*i2*u*u;
        std::vector<double> diffusion(ns_),work;
        const double mu=medium_.transport(t,p,y,diffusion.data(),work).viscosity;
        face.omega=p/(gasR*t)*face.k/(s.viscosityRatio*mu);
        face.turbulenceFlux={g*face.k,g*face.omega};
    }
    return face;
}
// With transport the diffusive rate 2 nu sum(A^2) / V is added to the convective one, with nu the
// largest of 4/3 mu / rho, lambda / (rho cv) and the mixture diffusion coefficients: pure diffusion
// then steps at cfl/2 of the explicit Euler bound V^2 / (nu sum(A^2)). With turbulence the eddy
// parts are added (mu + mu_t, lambda + cp mu_t / Pr_t, D_km + mu_t / (rho Sc_t); the k and omega
// diffusivities mu + sigma mu_t lie below 4/3 (mu + mu_t)), and the source coefficients of the
// first half of the split are frozen here.
double Flow::stableDt() {
    // AMR step B2: with turbulence the source coefficients over N_(g-1), which read the gradients'
    // neighbours in N_g; the step over the active cells.
    const int g=ghostWidth();
    refresh(state_,species_,turbulence_,nt_?g:0);
    const bool transport=medium_.hasTransport();
    if(transport) transportProperties(nt_?g:0);
    if(nt_) { transportGradients(g-1);eddyViscosity(true,g-1); }
    // Each worker's smallest step, then the smallest of those (a minimum does not depend on the order).
    std::fill(workerDt_.begin(),workerDt_.end(),std::numeric_limits<double>::infinity());
    const auto* active=region(0);
    pool_->blocks(active?active->size():mesh_.cells.size(),kCells,[&](std::size_t begin,std::size_t end,int worker) {
      double dt=std::numeric_limits<double>::infinity();
      StepLimit limit;
      for(std::size_t n=begin;n<end;++n) {
        const std::size_t q=active?(*active)[n]:n;
        const int i=static_cast<int>(q/mesh_.nr),j=static_cast<int>(q%mesh_.nr);
        const auto& w=primitives_[q];
        double a=sound_[q],volume=mesh_.cells[q].volume;
        double sumArea2=sq(mesh_.axialArea(i,j))+sq(mesh_.axialArea(i+1,j));
        double rate=(mesh_.axialArea(i,j)+mesh_.axialArea(i+1,j))*(std::abs(w.uz)+a);
        for(int f:{j,j+1}) {
            auto ar=mesh_.radialAreaVector(i,f);
            rate+=std::abs(w.uz*ar[0]+w.ur*ar[1])+a*std::hypot(ar[0],ar[1]);
            sumArea2+=sq(ar[0])+sq(ar[1]);
        }
        const double convective=rate;
        if(transport) {
            const double* y=fractions_.data()+q*nw_;
            double nu=std::max(4.0/3*viscosity_[q],conductivity_[q]/medium_.cv(temperature_[q],y))/w.rho;
            for(std::size_t k=0;k<ns_;++k) nu=std::max(nu,diffusion_[q*ns_+k]);
            if(nt_) {
                const auto& tu=definition_.turbulence;
                nu=std::max({nu,4.0/3*(viscosity_[q]+eddy_[q])/w.rho,(conductivity_[q]+eddyConductivity_[q])/(medium_.cv(temperature_[q],y)*w.rho)});
                for(std::size_t k=0;k<ns_;++k) nu=std::max(nu,diffusion_[q*ns_+k]+eddy_[q]/(w.rho*tu.schmidt));
            }
            rate+=2*nu*sumArea2/volume;
        }
        const double cell=definition_.cfl*volume/rate;
        if(cell<dt && isActive(q)) { dt=cell;limit={q,convective/volume,(rate-convective)/volume}; }
      }
      if(dt<workerDt_[worker] || (dt==workerDt_[worker] && limit.cell<workerLimit_[worker].cell)) {
          workerDt_[worker]=dt;workerLimit_[worker]=limit;
      }
    });
    std::size_t best=0;
    for(std::size_t w=1;w<workerDt_.size();++w)
        if(workerDt_[w]<workerDt_[best] || (workerDt_[w]==workerDt_[best] && workerLimit_[w].cell<workerLimit_[best].cell)) best=w;
    stepLimit_=workerLimit_[best];
    return workerDt_[best];
}
// Radial reconstruction. Cell averages are r-weighted, so they are point values at the volume
// centroid r_c, not at the cell midpoint; near the axis the two differ by O(dr), which made the
// first-face states first-order. Slopes are therefore taken in physical r between centroids and
// evaluated at each face's area-weighted radius. The axis row uses the parity of the field:
// u_r is odd (u_r = s r, s limited against the slope to the next centroid) and rho, u_z, p are
// even (linear in r^2 about <r^2>). pressureSource_ is the exact integral of p/r dV over the
// reconstructed profile, so a pressure gradient balanced by a body force stays at rest.
void Flow::radialProfiles(int i,int lo,int hi) {
    const auto& m=mesh_;
    for(int j=lo;j<hi;++j) {
        auto q=m.index(i,j);const auto& cell=m.cells[q];auto c=reconstructionValues(q);
        double low=m.radialFaceRadius(i,j)-cell.r,high=m.radialFaceRadius(i,j+1)-cell.r;
        Conserved slope{};
        auto centroidSlope=[&](int a,int b,int k){
            return (reconstructionValues(m.index(i,b))[k]-reconstructionValues(m.index(i,a))[k])/(m.cells[m.index(i,b)].r-m.cells[m.index(i,a)].r);};
        if(definition_.secondOrder && j==0) {
            Conserved curvature{};
            if(m.nr>=3) {
                auto moment=[&](int b){return m.cells[m.index(i,b)].radialSecondMoment;};
                auto w1=reconstructionValues(m.index(i,1)),w2=reconstructionValues(m.index(i,2));
                for(int k:{0,1,3}) curvature[k]=minmod((w1[k]-c[k])/(moment(1)-moment(0)),(w2[k]-w1[k])/(moment(2)-moment(1)));
            }
            if(m.nr>=2) slope[2]=minmod(c[2]/cell.r,centroidSlope(0,1,2));
            Primitive top=unpack(c);
            double face=m.radialFaceSecondMoment(i,1)-cell.radialSecondMoment;
            top.rho+=curvature[0]*face;top.uz+=curvature[1]*face;top.ur+=slope[2]*high;top.p+=curvature[3]*face;
            radialLow_[q]=unpack(c);radialHigh_[q]=top;
            pressureSource_[q]=c[3]*cell.radialPressureMeasure+curvature[3]*(cell.r*cell.volume-cell.radialSecondMoment*cell.radialPressureMeasure);
            continue;
        }
        if(definition_.secondOrder && j<m.nr-1)
            for(int k=0;k<4;++k) slope[k]=minmod(centroidSlope(j-1,j,k),centroidSlope(j,j+1,k));
        // Wall row: limited one-sided slope (see Flow::axialSlopes), with density and pressure face
        // values kept at least half the cell value.
        if(definition_.secondOrder && j==m.nr-1 && m.nr>=3) {
            double reach=std::max(-low,high);
            for(int k=0;k<4;++k) slope[k]=minmod(centroidSlope(j-1,j,k),centroidSlope(j-2,j-1,k));
            for(int k:{0,3}) slope[k]=std::clamp(slope[k],-0.5*c[k]/reach,0.5*c[k]/reach);
        }
        auto lowState=c,highState=c;
        for(int k=0;k<4;++k){lowState[k]+=slope[k]*low;highState[k]+=slope[k]*high;}
        radialLow_[q]=unpack(lowState);radialHigh_[q]=unpack(highState);
        pressureSource_[q]=c[3]*cell.radialPressureMeasure+slope[3]*(cell.volume-cell.r*cell.radialPressureMeasure);
    }
    // Mass fractions (and specific k and omega) are even in r and reconstructed like density: r^2
    // curvature on the axis row, limited centroid slopes inside, the limited one-sided slope (face >=
    // half the cell value) at the wall. The axis-row curvature is bounded like the wall slope: the
    // profile stays within half the cell value of it at the first face and on the axis. Unbounded, a
    // cell near zero beside large off-axis values gets a large face value, and outward flow there
    // removes more than the cell holds at any step size (omega past a throat, 5 October 2026).
    for(int j=lo;j<hi;++j) {
        auto q=m.index(i,j);const auto& cell=m.cells[q];
        double low=m.radialFaceRadius(i,j)-cell.r,high=m.radialFaceRadius(i,j+1)-cell.r;
        auto y=[&](int b,std::size_t k){return fractions_[m.index(i,b)*nw_+k];};
        auto centroidSlope=[&](int a,int b,std::size_t k){return (y(b,k)-y(a,k))/(m.cells[m.index(i,b)].r-m.cells[m.index(i,a)].r);};
        for(std::size_t k=0;k<nw_;++k) {
            double c=y(j,k),lowValue=c,highValue=c;
            if(definition_.secondOrder && j==0) {
                if(m.nr>=3) {
                    auto moment=[&](int b){return m.cells[m.index(i,b)].radialSecondMoment;};
                    double curvature=minmod((y(1,k)-c)/(moment(1)-moment(0)),(y(2,k)-y(1,k))/(moment(2)-moment(1)));
                    const double reach=std::max(m.radialFaceSecondMoment(i,1)-cell.radialSecondMoment,cell.radialSecondMoment);
                    curvature=std::clamp(curvature,-0.5*c/reach,0.5*c/reach);
                    highValue+=curvature*(m.radialFaceSecondMoment(i,1)-cell.radialSecondMoment);
                }
            } else if(definition_.secondOrder) {
                double slope=0;
                if(j<m.nr-1) slope=minmod(centroidSlope(j-1,j,k),centroidSlope(j,j+1,k));
                else if(m.nr>=3) {
                    double reach=std::max(-low,high);
                    slope=std::clamp(minmod(centroidSlope(j-1,j,k),centroidSlope(j-2,j-1,k)),-0.5*c/reach,0.5*c/reach);
                }
                lowValue+=slope*low;highValue+=slope*high;
            }
            fractionsLow_[q*nw_+k]=lowValue;fractionsHigh_[q*nw_+k]=highValue;
        }
    }
    // extremumBoundedEnthalpy: h's radial face values, profiled as density (r^2 curvature on the axis row, limited
    // centroid slopes inside, the limited one-sided slope at the wall), flat where the cell's excess is an extremum.
    if(definition_.extremumBoundedEnthalpy)
        for(int j=lo;j<hi;++j) {
            auto q=m.index(i,j);const auto& cell=m.cells[q];
            auto h=[&](int b){return cellStaticEnthalpy_[m.index(i,b)];};
            auto centroidSlope=[&](int a,int b){return (h(b)-h(a))/(m.cells[m.index(i,b)].r-m.cells[m.index(i,a)].r);};
            double c=h(j),lowValue=c,highValue=c;
            if(definition_.secondOrder && !excessExtremum(q)) {
                if(j==0) {
                    if(m.nr>=3) {
                        auto moment=[&](int b){return m.cells[m.index(i,b)].radialSecondMoment;};
                        highValue+=minmod((h(1)-c)/(moment(1)-moment(0)),(h(2)-h(1))/(moment(2)-moment(1)))
                                   *(m.radialFaceSecondMoment(i,1)-cell.radialSecondMoment);
                    }
                } else {
                    double low=m.radialFaceRadius(i,j)-cell.r,high=m.radialFaceRadius(i,j+1)-cell.r,slope=0;
                    if(j<m.nr-1) slope=minmod(centroidSlope(j-1,j),centroidSlope(j,j+1));
                    else if(m.nr>=3) slope=minmod(centroidSlope(j-1,j),centroidSlope(j-2,j-1));
                    lowValue+=slope*low;highValue+=slope*high;
                }
            }
            enthalpyLow_[q]=lowValue;enthalpyHigh_[q]=highValue;
        }
}
// Limited axial slopes of column i. Boundary columns have no outer neighbour. A zero slope there
// makes wall, inlet and outlet face states first-order (measured: wall-row entropy error order
// ~1.1). Use the limited one-sided slope minmod(c-n1, n1-n2), which is second-order for smooth data
// and falls back towards zero across a jump. Density/pressure/mass-fraction face values stay >= half
// the cell value.
void Flow::axialSlopes(int i,int lo,int hi) {
    const auto& m=mesh_;
    if(!definition_.secondOrder) {
        for(int j=lo;j<hi;++j) {
            auto q=m.index(i,j);slopesZ_[q]={};
            std::fill_n(fractionSlopesZ_.begin()+static_cast<std::ptrdiff_t>(q*nw_),nw_,0.0);
            if(definition_.extremumBoundedEnthalpy) enthalpySlopesZ_[q]=0;
        }
        return;
    }
    auto oneSided=[&](std::size_t q,std::size_t n1,std::size_t n2,double sign) {
        auto c=reconstructionValues(q),a=reconstructionValues(n1),b=reconstructionValues(n2);
        auto& slope=slopesZ_[q];
        for(int k=0;k<4;++k) slope[k]=sign*minmod(c[k]-a[k],a[k]-b[k]);
        for(int k:{0,3}) slope[k]=std::clamp(slope[k],-c[k],c[k]);
        for(std::size_t k=0;k<nw_;++k) {
            double yc=fractions_[q*nw_+k];
            fractionSlopesZ_[q*nw_+k]=std::clamp(sign*minmod(yc-fractions_[n1*nw_+k],fractions_[n1*nw_+k]-fractions_[n2*nw_+k]),-yc,yc);
        }
        if(definition_.extremumBoundedEnthalpy) {
            const auto& h=cellStaticEnthalpy_;
            enthalpySlopesZ_[q]=excessExtremum(q)?0.0:sign*minmod(h[q]-h[n1],h[n1]-h[n2]);
        }
    };
    for(int j=lo;j<hi;++j) {
        if(i==0) { oneSided(m.index(0,j),m.index(1,j),m.index(2,j),-1);continue; }
        if(i==m.nz-1) { oneSided(m.index(m.nz-1,j),m.index(m.nz-2,j),m.index(m.nz-3,j),1);continue; }
        auto q=m.index(i,j); auto c=reconstructionValues(q);
        auto ql=m.index(i-1,j),qr=m.index(i+1,j);
        auto l=reconstructionValues(ql),r=reconstructionValues(qr);
        for(int k=0;k<4;++k) slopesZ_[q][k]=minmod(c[k]-l[k],r[k]-c[k]);
        for(std::size_t k=0;k<nw_;++k)
            fractionSlopesZ_[q*nw_+k]=minmod(fractions_[q*nw_+k]-fractions_[ql*nw_+k],fractions_[qr*nw_+k]-fractions_[q*nw_+k]);
        if(definition_.extremumBoundedEnthalpy) {
            const auto& h=cellStaticEnthalpy_;
            enthalpySlopesZ_[q]=excessExtremum(q)?0.0:minmod(h[q]-h[ql],h[qr]-h[q]);
        }
    }
}
// Species follow the mass flux with the upwind face composition (Larrouturou, JCP 95, 1991):
// positivity-preserving, and the species fluxes sum to the mass flux. rho k and rho omega follow it
// likewise with the upwind face k and omega.
void Flow::axialFace(int i,int j,double time,std::vector<double>& yl,std::vector<double>& yr) {
    const auto& d=definition_;const auto& m=mesh_;
    const auto f=static_cast<std::size_t>(i)*m.nr+j;
    auto& stored=axialFlux_[f];double* transported=axialTransported_.data()+f*nw_;
    auto reconstructed=[&](std::size_t q,double direction,std::vector<double>& y) {
        auto w=reconstructionValues(q);
        for(int k=0;k<4;++k) w[k]+=direction*0.5*slopesZ_[q][k];
        for(std::size_t k=0;k<nw_;++k) y[k]=fractions_[q*nw_+k]+direction*0.5*fractionSlopesZ_[q*nw_+k];
        auto face=unpack(w);
        toFace(face,y.data());
        return face;
    };
    auto il=m.index(std::max(0,i-1),j),ir=m.index(std::min(m.nz-1,i),j);
    // Boundary models receive the reconstructed interior face state.
    Primitive l=reconstructed(il,1,yl),r=reconstructed(ir,-1,yr);
    double area=axialArea_[f];
    if(i==0 && d.experiment==Case::Chamber) {
        // Injector face: a supply ring (its species flux is g Y_supply, g = flux[0]) or the closed
        // plate (a slip wall, no species flux). Open-area supplies: phi_sj of each supply's face and w_j
        // of the plate's, the species and rho k, rho omega fluxes stored per unit area.
        if(d.openAreaSupplies) {
            Conserved sum{};
            std::fill(transported,transported+nw_,0.0);
            for(std::size_t n=shareBegin_[j];n<shareBegin_[j+1];++n) {
                const auto& share=faceShares_[n];const auto& supply=d.supplies[static_cast<std::size_t>(share.supply)];
                const double g=supply.massFlow*supply.opening(time)/supplyArea_[static_cast<std::size_t>(share.supply)];
                const auto face=supplyFace(supply,g,r,yr.data());
                for(int k=0;k<4;++k) sum[k]+=share.phi*face.flux[k];
                for(std::size_t k=0;k<ns_;++k) transported[k]+=share.phi*face.flux[0]*supply.composition[k];
                for(std::size_t k=0;k<nt_;++k) transported[ns_+k]+=share.phi*face.turbulenceFlux[k];
            }
            plateMomentum_[j]=0;
            if(const double w=plateShare_[j];w>0) {
                Thermal t=faceThermal(r,yr.data());
                const auto plate=faceFlux(d.lowMach,reflect(r,1,0),t,r,t,1,0,true);
                plateMomentum_[j]=w*plate[1];sum[1]+=plateMomentum_[j];sum[2]+=w*plate[2];
            }
            stored=sum;
            return;
        }
        if(int s=faceSupply_[j];s>=0) {
            const auto& supply=d.supplies[s];
            double g=supply.massFlow*supply.opening(time)/supplyArea_[s];
            const auto face=supplyFace(supply,g,r,yr.data());
            stored=face.flux;
            for(std::size_t k=0;k<nt_;++k) transported[ns_+k]=face.turbulenceFlux[k];
        } else {
            Thermal t=faceThermal(r,yr.data());
            stored=faceFlux(d.lowMach,reflect(r,1,0),t,r,t,1,0,true);
            stored[0]=0;stored[3]=0;
        }
        return;
    }
    // Inflow at the reservoir composition; a transmissive inlet (UniformDuct, ShockTube) passes the
    // interior k and omega.
    if(i==0) { l=inlet(r);std::copy(inletComposition_.begin(),inletComposition_.end(),yl.begin());std::copy(yr.begin()+ns_,yr.end(),yl.begin()+ns_); }
    if(d.extremumBoundedEnthalpy && i>0 && i<m.nz) {
        enthalpyFace(l,cellStaticEnthalpy_[il]+0.5*enthalpySlopesZ_[il],yl.data());
        enthalpyFace(r,cellStaticEnthalpy_[ir]-0.5*enthalpySlopesZ_[ir],yr.data());
    }
    Thermal tl=faceThermal(l,yl.data());
    if(i==m.nz) { bool ambientInflow=false;r=outlet(l,yl.data(),primitives_[il],ambientInflow);yr=ambientInflow?ambient_:yl; }
    Thermal tr=faceThermal(r,yr.data());
    if(i>0 && i<m.nz) { boundEnthalpy(l,tl,yl.data(),il,ir);boundEnthalpy(r,tr,yr.data(),ir,il); }
    stored=faceFlux(d.lowMach,l,tl,r,tr,1,0,i>0 && i<m.nz);
    const auto& y=stored[0]>=0?yl:yr;
    for(std::size_t k=0;k<nw_;++k) transported[k]=area*stored[0]*y[k];
}
std::vector<Conserved> Flow::injectorFluxes() {
    if(definition_.experiment!=Case::Chamber) throw std::logic_error("Flow::injectorFluxes: a Chamber's.");
    refresh(state_,species_,turbulence_,kWhole);
    axialSlopes(0,0,mesh_.nr);
    std::vector<Conserved> out(static_cast<std::size_t>(mesh_.nr));
    for(int j=0;j<mesh_.nr;++j) { axialFace(0,j,time_,faceFractions_[0],faceFractions_[1]);out[j]=axialFlux_[j]; }
    return out;
}
void Flow::radialFace(int i,int j,std::vector<double>& yl,std::vector<double>& yr) {
    const auto& m=mesh_;
    const auto f=static_cast<std::size_t>(i)*(m.nr+1)+j;
    auto il=m.index(i,j-1),ir=m.index(i,std::min(m.nr-1,j));
    auto ar=m.radialAreaVector(i,j);double area=radialArea_[f];
    double nz=ar[0]/area,nr=ar[1]/area;
    Primitive l=radialHigh_[il];
    std::copy_n(fractionsHigh_.begin()+static_cast<std::ptrdiff_t>(il*nw_),nw_,yl.begin());
    toFace(l,yl.data());
    Primitive r=l;
    if(j==m.nr) { r=reflect(l,nz,nr);yr=yl; }
    else { r=radialLow_[ir];std::copy_n(fractionsLow_.begin()+static_cast<std::ptrdiff_t>(ir*nw_),nw_,yr.begin());toFace(r,yr.data()); }
    if(definition_.extremumBoundedEnthalpy && j<m.nr) { enthalpyFace(l,enthalpyHigh_[il],yl.data());enthalpyFace(r,enthalpyLow_[ir],yr.data()); }
    Thermal tl=faceThermal(l,yl.data()),tr=faceThermal(r,yr.data());
    if(j<m.nr) { boundEnthalpy(l,tl,yl.data(),il,ir);boundEnthalpy(r,tr,yr.data(),ir,il); }
    auto& flux=radialFlux_[f];
    flux=faceFlux(definition_.lowMach,l,tl,r,tr,nz,nr,true);
    // A stationary slip wall has exactly zero mass and energy exchange.
    if(j==m.nr) { flux[0]=0;flux[3]=0; }
    const auto& y=flux[0]>=0?yl:yr;
    for(std::size_t k=0;k<nw_;++k) radialTransported_[f*nw_+k]=area*flux[0]*y[k];
}
// The igniter's share, the faces in the order the serial face loops reached the cell (axial faces
// i then i + 1, radial faces j then j + 1), then the body force.
void Flow::gather(std::size_t q,std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative,
                  std::vector<double>& turbulenceDerivative) const {
    const auto& m=mesh_;
    const int j=static_cast<int>(q%m.nr);
    const std::size_t i=q/m.nr,a0=q,a1=q+m.nr,r0=i*(m.nr+1)+j,r1=r0+1;
    const bool injector=i==0 && definition_.experiment==Case::Chamber;
    const int supply=injector?faceSupply_[j]:-1;
    Conserved u{};
    if(igniterPower_>0 && igniterCell_[q]) u[3]+=igniterPower_*igniterCellVolume_[q]/igniterVolume_;
    for(int k=0;k<4;++k) u[k]+=axialArea_[a0]*axialFlux_[a0][k];
    for(int k=0;k<4;++k) u[k]-=axialArea_[a1]*axialFlux_[a1][k];
    if(j>0) for(int k=0;k<4;++k) u[k]+=radialArea_[r0]*radialFlux_[r0][k];
    for(int k=0;k<4;++k) u[k]-=radialArea_[r1]*radialFlux_[r1][k];
    if(!bodyForce_.empty()) {
        const auto& f=bodyForce_[q];double v=m.cells[q].volume;
        double work=(f[0]*primitives_[q].uz+f[1]*primitives_[q].ur)*v;
        u[1]+=f[0]*v;u[2]+=f[1]*v;u[3]+=work;
    }
    if(!heating_.empty()) u[3]+=heating_[q]*m.cells[q].volume;
    derivative[q]=u;
    for(std::size_t k=0;k<nw_;++k) {
        double s=0;
        if(injector && definition_.openAreaSupplies) s+=axialArea_[a0]*axialTransported_[a0*nw_+k];
        else if(supply>=0) {
            const double area=axialArea_[a0],g=axialFlux_[a0][0];
            if(k<ns_) s+=area*g*definition_.supplies[supply].composition[k];
            else s+=area*axialTransported_[a0*nw_+k];
        }
        else if(!injector) s+=axialTransported_[a0*nw_+k];
        s-=axialTransported_[a1*nw_+k];
        if(j>0) s+=radialTransported_[r0*nw_+k];
        s-=radialTransported_[r1*nw_+k];
        if(k<ns_) speciesDerivative[q*ns_+k]=s; else turbulenceDerivative[q*nt_+k-ns_]=s;
    }
}
// The faces are computed alone into face arrays and summed per cell in the order of the serial face
// loops, and the boundary exchange rates are summed from them in that order, so the result is the
// same bit for bit on any number of threads.
Flow::BoundaryRates Flow::rhs(const std::vector<Conserved>& state,const std::vector<double>& species,const std::vector<double>& turbulence,
                              std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative,
                              std::vector<double>& turbulenceDerivative,double time) {
    const auto& d=definition_;const auto& m=mesh_;
    // AMR step B2: the cells within 2 of an active cell, the slopes and profiles over N_1 (their stencils
    // reach N_2), and the faces with an active cell on a side.
    refresh(state,species,turbulence,2);
    const bool restricted=!within_.empty();
    // Per column: axial slopes, radial profiles and the radial faces (which use only the column's profiles).
    pool_->blocks(restricted?columns_.size():static_cast<std::size_t>(m.nz),1,[&](std::size_t begin,std::size_t end,int worker) {
        auto& yl=faceFractions_[2*worker];auto& yr=faceFractions_[2*worker+1];
        for(std::size_t c=begin;c<end;++c) {
            if(!restricted) {
                const auto i=static_cast<int>(c);
                axialSlopes(i,0,m.nr);radialProfiles(i,0,m.nr);
                for(int j=1;j<=m.nr;++j) radialFace(i,j,yl,yr);
                continue;
            }
            const int i=columns_[c].i;
            for(const auto& rows:columns_[c].rows) { axialSlopes(i,rows[0],rows[1]);radialProfiles(i,rows[0],rows[1]); }
            for(int j=1;j<=m.nr;++j)
                if(isActive(m.index(i,j-1)) || (j<m.nr && isActive(m.index(i,j)))) radialFace(i,j,yl,yr);
        }
    });
    pool_->blocks(restricted?axialFaces_:axialFlux_.size(),kFaces,[&](std::size_t begin,std::size_t end,int worker) {
        auto& yl=faceFractions_[2*worker];auto& yr=faceFractions_[2*worker+1];
        for(std::size_t n=begin;n<end;++n) {
            const std::size_t f=restricted?faces_[n]:n;
            axialFace(static_cast<int>(f/m.nr),static_cast<int>(f%m.nr),time,yl,yr);
        }
    });
    const bool transport=medium_.hasTransport();
    // A wall-function cell's k and omega are prescribed after the stage, not advanced.
    const auto& prescribed=prescribedCell_;
    auto finish=[&](std::size_t q) {
        derivative[q][2]+=pressureSource_[q];
        for(double& v:derivative[q]) v/=m.cells[q].volume;
        for(std::size_t k=0;k<ns_;++k) speciesDerivative[q*ns_+k]/=m.cells[q].volume;
        for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]=prescribed[q]?0.0:turbulenceDerivative[q*nt_+k]/m.cells[q].volume;
    };
    eachCell(region(0),[&](std::size_t q,int) { gather(q,derivative,speciesDerivative,turbulenceDerivative);if(!transport) finish(q); });
    if(!probe_.cells.empty()) probeConvective();
    // Boundary exchange rates: the igniter (constant power over the step, set by step()), the inlet
    // and outlet faces, the side wall, and the volumetric forces (e.g. a Lorentz force), whose axial
    // integral is rates.bodyAxial, reacting on the equipment (coils) in deviceThrust, and whose work
    // is an external energy exchange. Under AMR only the owned cells' (isOwned) are booked, the
    // igniter as its owned cells' share of the power.
    BoundaryRates rates{};
    if(owned_.empty()) { rates.energy+=igniterPower_;rates.heat+=igniterPower_; }
    else if(igniterPower_>0) {
        double share=0;
        for(auto q:igniterCells_) if(owned_[q]) share+=igniterCellVolume_[q];
        const double power=igniterPower_*share/igniterVolume_;
        rates.energy+=power;rates.heat+=power;
    }
    for(int j=0;j<m.nr;++j) {
        if(!isOwned(m.index(0,j))) continue;
        const auto& flux=axialFlux_[j];const double area=axialArea_[j];
        if(d.experiment==Case::Chamber && d.openAreaSupplies) {
            // The plate share's pressure force is the wall's; the supplies' shares (closed ports too) the inlet's.
            rates.wallAxial+=area*plateMomentum_[j];
            rates.mass+=area*flux[0];rates.energy+=area*flux[3];rates.inlet+=area*flux[0];rates.inletMomentum+=area*(flux[1]-plateMomentum_[j]);
            continue;
        }
        if(d.experiment==Case::Chamber && faceSupply_[j]<0) { rates.wallAxial+=area*flux[1];continue; }
        rates.mass+=area*flux[0];rates.energy+=area*flux[3];rates.inlet+=area*flux[0];rates.inletMomentum+=area*flux[1];
    }
    for(int j=0;j<m.nr;++j) {
        if(!isOwned(m.index(m.nz-1,j))) continue;
        const auto f=static_cast<std::size_t>(m.nz)*m.nr+j;
        const auto& flux=axialFlux_[f];const double area=axialArea_[f];
        rates.mass-=area*flux[0];rates.energy-=area*flux[3];rates.outlet+=area*flux[0];rates.outletMomentum+=area*flux[1];
    }
    for(int i=0;i<m.nz;++i) {
        if(!isOwned(m.index(i,m.nr-1))) continue;
        const auto f=static_cast<std::size_t>(i)*(m.nr+1)+m.nr;
        rates.wallAxial-=radialArea_[f]*radialFlux_[f][1];
    }
    if(!bodyForce_.empty()) for(std::size_t q=0;q<state.size();++q) {
        if(!isOwned(q)) continue;
        const auto& f=bodyForce_[q];double v=m.cells[q].volume;
        double work=(f[0]*primitives_[q].uz+f[1]*primitives_[q].ur)*v;
        rates.bodyAxial+=f[0]*v;rates.energy+=work;
    }
    if(!heating_.empty()) for(std::size_t q=0;q<state.size();++q) if(isOwned(q)) rates.energy+=heating_[q]*m.cells[q].volume;
    if(transport) {
        transportFluxes(derivative,speciesDerivative,turbulenceDerivative,rates);
        eachCell(region(0),[&](std::size_t q,int) { finish(q); });
    }
    return rates;
}
// The wall functions' first cells (core/transport.cpp): each cell's own state gives T_1 and mu_1; the
// face that governs the cell gives k and omega; E is unchanged.
bool Flow::prescribeWallTurbulence(const std::vector<Conserved>& state,const std::vector<double>& species,
                                   std::vector<double>& turbulence) {
    std::fill(workerOk_.begin(),workerOk_.end(),char{1});
    pool_->blocks(wallCells_.size(),kFaces,[&](std::size_t begin,std::size_t end,int worker) {
        auto& s=transportScratch_[worker];double* y=faceFractions_[2*worker].data();
        for(std::size_t n=begin;n<end;++n)
            if(!prescribeCell(n,state,species,turbulence,s,y)) workerOk_[worker]=0;
    });
    return std::find(workerOk_.begin(),workerOk_.end(),char{0})==workerOk_.end();
}
bool Flow::prescribeCell(std::size_t n,const std::vector<Conserved>& state,const std::vector<double>& species,
                         std::vector<double>& turbulence,TransportScratch& s,double* y) {
    const std::size_t q=wallCells_[n];const auto& u=state[q];
    prescribedChange_[q]=0;
    normalise(species.data()+q*ns_,ns_,y);
    const double k=turbulence[q*nt_]/u[0];
    if(!admissibleBulk(u,medium_.energyFloor(y)+k)) return false;
    const double t=medium_.temperature((u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))/u[0]-k,y,temperature_[q]);
    const Primitive w{u[0],u[1]/u[0],u[2]/u[0],u[0]*medium_.gasConstant(y)*t};
    const double mu1=medium_.transport(t,w.p,y,s.diffusion.data(),s.work).viscosity;
    const auto law=wallSolve(wallFaces_[wallCellFace_[n]],w,t,y,mu1,s);
    // rho k = omega mu_t.
    const double rhoK=law.omega*law.eddy;
    prescribedChange_[q]=(rhoK-turbulence[q*nt_])*mesh_.cells[q].volume;
    turbulence[q*nt_]=rhoK;turbulence[q*nt_+1]=u[0]*law.omega;
    return admissible(u,species.data()+q*ns_,turbulence.data()+q*nt_);
}
double Flow::prescriptionResidual() {
    if(!(nt_ && definition_.turbulence.wallFunctions)) return 0;
    std::vector<double> law=turbulence_;
    const auto change=prescribedChange_;
    (void)prescribeWallTurbulence(state_,species_,law);
    prescribedChange_=change;
    double worst=0;
    for(std::size_t q:wallCells_) if(isOwned(q)) for(std::size_t k=0;k<nt_;++k) {
        const double a=turbulence_[q*nt_+k],b=law[q*nt_+k];
        worst=std::max(worst,std::abs(a-b)/std::max(std::abs(b),std::numeric_limits<double>::min()));
    }
    return worst;
}
bool Flow::represcribe(const std::vector<char>& touched) {
    bool admitted=true;
    for(std::size_t n=0;n<wallCells_.size();++n) {
        const std::size_t q=wallCells_[n];
        if(!touched[q]) continue;
        if(!prescribeCell(n,state_,species_,turbulence_,transportScratch_[0],faceFractions_[0].data())) admitted=false;
        if(isOwned(q)) prescribedTurbulentEnergy_+=prescribedChange_[q];
    }
    return admitted;
}
// One step: SSPRK2 on the fluxes. With turbulence the SST sources are Strang-split around it: half a
// step of source (coefficients frozen at the start, from stableDt), the flux step, and half a step
// of source with coefficients at the new state. A rejected attempt restarts from the saved rho k
// and rho omega with half the step.
// The k and omega diffusion is not sign-preserving where the line of centres is not along the face
// normal (contoured rings): the tangential part of the face gradient can draw k out of a cell at
// k = 0. Next to quiescent gas (k = 0) the numerical domain of dependence carries values like
// 1e-185 ahead of the turbulent front, and the update can leave rho k at -1e-198, which no step
// reduction cures. After each stage rho k is therefore set to max(rho k, 0); E is not touched (k is
// part of it), so the energy budget stays exact. The rho k V added is summed over stages in
// clippedTurbulentEnergy, so a clip that matters shows up in the measurements.
void Flow::beginStep() {
    if(nt_) { turbulenceStart_=turbulence_;sourcesStart_=sources_; }
}
// Each cell's clip is kept, and the stage's are summed in cell order; likewise the wall functions'
// prescription, in the order of wallCells_.
void Flow::clipTurbulence(std::vector<double>& turbulence,std::size_t q) {
    if(!nt_) return;
    clippedCell_[q]=0;
    if(turbulence[q*nt_]<0) { clippedCell_[q]=turbulence[q*nt_]*mesh_.cells[q].volume;turbulence[q*nt_]=0; }
}
void Flow::sumClipped(Attempt& attempt) const {
    if(!nt_) return;
    // Under AMR the owned cells' only (AMR.md, step D).
    if(const auto* active=region(0)) { for(std::size_t q:*active) if(isOwned(q)) attempt.clipped-=clippedCell_[q]; }
    else for(double c:clippedCell_) attempt.clipped-=c;
}
bool Flow::prescribe(std::vector<Conserved>& state,std::vector<double>& species,std::vector<double>& turbulence,Attempt& attempt) {
    if(!(nt_ && definition_.turbulence.wallFunctions)) return true;
    const bool admitted=prescribeWallTurbulence(state,species,turbulence);
    for(std::size_t q:wallCells_) if(isOwned(q)) attempt.prescribed+=prescribedChange_[q];
    return admitted;
}
bool Flow::firstStage(double dt,Attempt& attempt) {
    attempt=Attempt{};
    // Mean igniter power over [t, t + dt]; steps end on its switching times, so this is exact.
    const auto& ig=definition_.igniter;
    igniterPower_=0;
    if(ig.energy>0) {
        double overlap=std::min(time_+dt,ig.start+ig.duration)-std::max(time_,ig.start);
        if(overlap>0) igniterPower_=ig.energy/ig.duration*overlap/dt;
    }
    if(recordFaces_) {
        if(within_.empty()) std::fill(faceIntegral_.begin(),faceIntegral_.end(),0.0);
        else for(std::size_t f:faces_) std::fill_n(faceIntegral_.begin()+static_cast<std::ptrdiff_t>(f*(4+nw_)),4+nw_,0.0);
    }
    if(nt_) {
        turbulence_=turbulenceStart_;
        turbulenceSource(state_,turbulence_,sourcesStart_,0.5*dt,ghostWidth()-1);
        if(!allOf(region(0),[&](std::size_t q) { return !isActive(q) || admissible(state_[q],species_.data()+q*ns_,turbulence_.data()+q*nt_); }))
            return false;
    }
    attempt.first=rhs(state_,species_,turbulence_,rhs_,speciesRhs_,turbulenceRhs_,time_);
    if(recordFaces_) recordFaceIntegrals(0.5*dt);
    for(std::size_t n=0;n<probeStage_.size();++n) probeAttempt_[n]=0.5*dt*probeStage_[n];
    // A cell that is not active keeps its value through both stages; B2 keeps only the filled ones, the
    // others (off cells) having been set to their state when the regions were set.
    const auto keep=[&](std::size_t q) {
        stage_[q]=state_[q];
        std::copy_n(species_.begin()+static_cast<std::ptrdiff_t>(q*ns_),ns_,speciesStage_.begin()+static_cast<std::ptrdiff_t>(q*ns_));
        std::copy_n(turbulence_.begin()+static_cast<std::ptrdiff_t>(q*nt_),nt_,turbulenceStage_.begin()+static_cast<std::ptrdiff_t>(q*nt_));
        if(nt_) clippedCell_[q]=0;
        return true;
    };
    const bool ok=allOf(region(0),[&](std::size_t q) {
        if(!isActive(q)) return keep(q);
        for(int k=0;k<4;++k) stage_[q][k]=state_[q][k]+dt*rhs_[q][k];
        for(std::size_t k=q*ns_;k<(q+1)*ns_;++k) speciesStage_[k]=species_[k]+dt*speciesRhs_[k];
        for(std::size_t k=q*nt_;k<(q+1)*nt_;++k) turbulenceStage_[k]=turbulence_[k]+dt*turbulenceRhs_[k];
        clipTurbulence(turbulenceStage_,q);
        return admissible(stage_[q],speciesStage_.data()+q*ns_,turbulenceStage_.data()+q*nt_);
    });
    if(!within_.empty()) allOf(&filled_,keep);
    sumClipped(attempt);
    return ok;
}
bool Flow::secondStage(double dt,Attempt& attempt) {
    if(!prescribe(stage_,speciesStage_,turbulenceStage_,attempt)) return false;
    attempt.second=rhs(stage_,speciesStage_,turbulenceStage_,rhs_,speciesRhs_,turbulenceRhs_,time_+dt);
    if(recordFaces_) recordFaceIntegrals(0.5*dt);
    for(std::size_t n=0;n<probeStage_.size();++n) probeAttempt_[n]+=0.5*dt*probeStage_[n];
    const auto keep=[&](std::size_t q) {
        next_[q]=stage_[q];
        std::copy_n(speciesStage_.begin()+static_cast<std::ptrdiff_t>(q*ns_),ns_,speciesNext_.begin()+static_cast<std::ptrdiff_t>(q*ns_));
        std::copy_n(turbulenceStage_.begin()+static_cast<std::ptrdiff_t>(q*nt_),nt_,turbulenceNext_.begin()+static_cast<std::ptrdiff_t>(q*nt_));
        if(nt_) clippedCell_[q]=0;
        return true;
    };
    const bool ok=allOf(region(0),[&](std::size_t q) {
        if(!isActive(q)) return keep(q);
        for(int k=0;k<4;++k) next_[q][k]=0.5*(state_[q][k]+stage_[q][k]+dt*rhs_[q][k]);
        for(std::size_t k=q*ns_;k<(q+1)*ns_;++k) speciesNext_[k]=0.5*(species_[k]+speciesStage_[k]+dt*speciesRhs_[k]);
        for(std::size_t k=q*nt_;k<(q+1)*nt_;++k) turbulenceNext_[k]=0.5*(turbulence_[k]+turbulenceStage_[k]+dt*turbulenceRhs_[k]);
        clipTurbulence(turbulenceNext_,q);
        return admissible(next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_);
    });
    if(!within_.empty()) allOf(&filled_,keep);
    sumClipped(attempt);
    if(!(ok && prescribe(next_,speciesNext_,turbulenceNext_,attempt))) return false;
    if(nt_) {
        // B2: the cells within distance 1 of the active ones, which the level above's fill reads.
        refresh(next_,speciesNext_,turbulenceNext_,2);
        transportProperties(2);transportGradients(1);eddyViscosity(true,1);
        turbulenceSource(next_,turbulenceNext_,sources_,0.5*dt,1);
        return allOf(region(0),[&](std::size_t q) { return !isActive(q) || admissible(next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_); });
    }
    return true;
}
void Flow::accept(double dt,double newTime,const Attempt& attempt) {
    state_.swap(next_);species_.swap(speciesNext_);turbulence_.swap(turbulenceNext_);
    time_=newTime;dt_=dt;++steps_;
    for(std::size_t n=0;n<probeAttempt_.size();++n) probe_.sums[n]+=probeAttempt_[n];
    lastRates_=BoundaryRates::average(attempt.first,attempt.second);
    // Neumaier: the plain sum as before, and its rounding gathered beside it.
    auto add=[](double& plain,double& compensation,double term) {
        const double sum=plain+term;
        compensation+=std::abs(plain)>=std::abs(term)?(plain-sum)+term:(term-sum)+plain;
        plain=sum;
    };
    add(integratedHeat_,heatCompensation_,dt*lastRates_.heat);
    clippedTurbulentEnergy_+=attempt.clipped;prescribedTurbulentEnergy_+=attempt.prescribed;
    add(integratedMassFlux_,massCompensation_,dt*lastRates_.mass);
    add(integratedEnergyFlux_,energyCompensation_,dt*lastRates_.energy);
    add(integratedCorrection_,correctionCompensation_,dt*lastRates_.correction);
    add(integratedMomentumSource_,momentumCompensation_,dt*lastRates_.netMomentum());
    integratedMomentumGross_+=dt*lastRates_.grossMomentum();
}
// The convective face fluxes, and with transport the viscous ones (stored alone by transportFluxes in
// the same orientation; they carry no mass, and their transported entries are per unit area). At a
// supply or boundary face the entries are not those the cells use; no coarse-fine face lies on a boundary.
void Flow::recordFaceIntegrals(double weight) {
    const std::size_t stride=4+nw_,na=axialFlux_.size();
    const bool viscous=medium_.hasTransport();
    const bool restricted=!within_.empty();
    pool_->blocks(restricted?faces_.size():na+radialFlux_.size(),kCells,[&](std::size_t begin,std::size_t end,int) {
        for(std::size_t n=begin;n<end;++n) {
            const std::size_t f=restricted?faces_[n]:n;
            const bool axial=f<na;const std::size_t e=axial?f:f-na;
            const double area=axial?axialArea_[e]:radialArea_[e];
            const auto& flux=axial?axialFlux_[e]:radialFlux_[e];
            const double* moved=(axial?axialTransported_:radialTransported_).data()+e*nw_;
            double* g=faceIntegral_.data()+f*stride;
            if(!viscous) {
                for(int k=0;k<4;++k) g[k]+=weight*(area*flux[k]);
                for(std::size_t k=0;k<nw_;++k) g[4+k]+=weight*moved[k];
                continue;
            }
            const auto& v=axial?axialViscous_[e]:radialViscous_[e];
            const double* t=(axial?axialViscousTransported_:radialViscousTransported_).data()+e*nw_;
            g[0]+=weight*(area*flux[0]);
            for(int k=1;k<4;++k) g[k]+=weight*(area*flux[k]+area*v[k]);
            for(std::size_t k=0;k<nw_;++k) g[4+k]+=weight*(moved[k]+area*t[k]);
        }
    });
}
double Flow::step(double maxDt) {
    if(!(maxDt>0)) throw std::invalid_argument("Step interval must be positive.");
    double dt=std::min(stableDt(),maxDt);
    const double event=nextEvent(time_);
    bool landing=false;
    if(event-time_<=dt) { dt=event-time_;landing=true; }
    beginStep();
    for(int attempt=0;attempt<14;++attempt) {
        Attempt a;
        if(firstStage(dt,a) && secondStage(dt,a)) { accept(dt,landing?event:time_+dt,a);return dt; }
        ++rejectedSteps_;dt*=0.5;landing=false;
    }
    if(nt_) turbulence_=turbulenceStart_;
    throw std::runtime_error("Could not advance an admissible gas state after 14 timestep reductions.");
}
void Flow::advanceTo(double target) {
    if(!std::isfinite(target) || target<time_) throw std::invalid_argument("Invalid target time.");
    while(time_<target) step(target-time_);
}
Measurements Flow::measurements() const {
    Measurements result{};
    result.time=time_;result.dt=dt_;result.steps=steps_;result.rejectedSteps=rejectedSteps_;
    result.clippedTurbulentEnergy=clippedTurbulentEnergy_;result.prescribedTurbulentEnergy=prescribedTurbulentEnergy_;
    result.inletMassFlow=lastRates_.inlet;result.outletMassFlow=lastRates_.outlet;
    result.inletMomentumFlux=lastRates_.inletMomentum;result.outletMomentumFlux=lastRates_.outletMomentum;
    result.wallAxialForce=lastRates_.wallAxial;result.bodyAxialForce=lastRates_.bodyAxial;result.wallHeatFlow=lastRates_.wallHeat;result.plateHeatFlow=lastRates_.plateHeat;
    result.minPressure=std::numeric_limits<double>::infinity();
    double exitArea=0;
    for(int i=0;i<mesh_.nz;++i) for(int j=0;j<mesh_.nr;++j) {
        auto q=mesh_.index(i,j);auto w=cellPrimitive(q);auto y=massFractions(q);
        double mach=std::hypot(w.uz,w.ur)/medium_.soundSpeed(w.p/(w.rho*medium_.gasConstant(y.data())),y.data());
        result.mass+=state_[q][0]*mesh_.cells[q].volume;
        result.energy+=state_[q][3]*mesh_.cells[q].volume;
        result.axialMomentum+=state_[q][1]*mesh_.cells[q].volume;
        result.minPressure=std::min(result.minPressure,w.p);result.maxPressure=std::max(result.maxPressure,w.p);
        result.maxMach=std::max(result.maxMach,mach);
        if(i==mesh_.nz-1) { double a=mesh_.axialArea(mesh_.nz,j);exitArea+=a;result.exitMach+=a*mach; }
    }
    result.exitMach/=exitArea;
    if(definition_.experiment==Case::Chamber) {
        for(const auto& s:definition_.supplies) result.supplyMassFlow+=s.massFlow*s.opening(time_);
        result.igniterPower=lastRates_.heat;result.igniterEnergy=deliveredHeat();
        double area=0;
        for(int j=0;j<mesh_.nr;++j) { double a=mesh_.axialArea(0,j);area+=a;result.injectorPressure+=a*cellPrimitive(mesh_.index(0,j)).p; }
        result.injectorPressure/=area;
    }
    // Ambient pressure acts on the closed exterior of the device except the exit opening.
    result.ambientAxialForce=definition_.backPressure*exitArea;
    if(steps_>0) {
        result.exitPlaneThrust=result.outletMomentumFlux-result.ambientAxialForce;
        result.deviceThrust=result.inletMomentumFlux+result.wallAxialForce+result.bodyAxialForce-result.ambientAxialForce;
    }
    // Normalized by initial |P_z| plus the integrated gross momentum exchange (P_z can start at zero).
    result.exchangedMass=exchangedMass_;result.exchangedEnergy=exchangedEnergy_;result.exchangedMomentum=exchangedMomentum_;
    result.momentumBalanceScale=std::max(std::abs(initialMomentum_)+integratedMomentumGross_,std::numeric_limits<double>::min());
    result.momentumBalanceError=(result.axialMomentum-initialMomentum_-integratedMomentumSource_-exchangedMomentum_)/result.momentumBalanceScale;
    result.massBalanceError=(result.mass-initialMass_-integratedMassFlux_-exchangedMass_)/initialMass_;
    result.energyBalanceError=(result.energy-initialEnergy_-integratedEnergyFlux_-integratedCorrection_-exchangedEnergy_)/initialEnergy_;
    result.compensatedMomentumBalanceError=(result.axialMomentum-initialMomentum_-compensatedMomentumSource()-exchangedMomentum_)/result.momentumBalanceScale;
    result.compensatedMassBalanceError=(result.mass-initialMass_-compensatedMassFlux()-exchangedMass_)/initialMass_;
    result.compensatedEnergyBalanceError=(result.energy-initialEnergy_-compensatedEnergyFlux()-compensatedCorrection()-exchangedEnergy_)/initialEnergy_;
    result.correctionEnergy=integratedCorrection_;
    if(corrections_) {
        double weight=0,floored=0;
        for(std::size_t q=0;q<state_.size();++q) { weight+=blendF1_[q];if(snFloored_[q]) floored+=blendF1_[q]; }
        result.correctionFloored=weight>0?floored/weight:0;
    }
    return result;
}
FieldSnapshot Flow::snapshot() const {
    FieldSnapshot out;
    out.definition=definition_;out.radius=mesh_.radius;out.fraction=mesh_.fraction;out.measurements=measurements();out.appliedTotalPressure=totalPressure_;
    out.cells.reserve(state_.size());for(std::size_t q=0;q<state_.size();++q) out.cells.push_back(cellPrimitive(q));
    return out;
}
} // namespace crucible
