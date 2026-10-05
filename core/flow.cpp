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
    if(turbulence.enabled) {
        const auto& t=turbulence;
        if(transport.viscosity.empty()) throw std::invalid_argument("Turbulence needs molecular transport.");
        if(experiment!=Case::UniformDuct && experiment!=Case::ShockTube)
            throw std::invalid_argument("Turbulent inflow conditions are not yet implemented: turbulence runs on UniformDuct and ShockTube only.");
        for(double x:{t.prandtl,t.schmidt,t.wallOmegaFactor,t.ambientOmega})
            if(!(x>0) || !std::isfinite(x)) throw std::invalid_argument("Turbulent Prandtl and Schmidt numbers, the wall omega factor and the ambient omega must be positive.");
        if(!(t.ambientK>=0) || !std::isfinite(t.ambientK)) throw std::invalid_argument("Ambient k must be finite and non-negative.");
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
        cells[index(i,j)]={v,(i+0.5)*dz,first/v,pi*dz*(r0+r1)*(b-a),second/v};
    }
}
double Mesh::axialArea(int face,int j) const {
    if(uniform) return pi*sq(radius[face])*(2*j+1)/sq(nr);
    return pi*sq(radius[face])*(sq(fraction[j+1])-sq(fraction[j]));
}
double Mesh::ringMiddle(int face,int j) const {
    if(uniform) return radius[face]*(j+0.5)/nr;
    return radius[face]*0.5*(fraction[j]+fraction[j+1]);
}
double Mesh::radialFaceMiddle(int i,int j) const {
    if(uniform) return 0.5*(radius[i]+radius[i+1])*j/nr;
    return 0.5*(radius[i]+radius[i+1])*fraction[j];
}
std::array<double,2> Mesh::radialAreaVector(int i,int face) const {
    double f=fraction[face];
    double r0=f*radius[i],r1=f*radius[i+1];
    return {-pi*(r0+r1)*(r1-r0),pi*(r0+r1)*dz};
}
double Mesh::radialFaceRadius(int i,int face) const {
    double f=fraction[face],r0=radius[i],r1=radius[i+1];
    return face==0?0:2*f*(r0*r0+r0*r1+r1*r1)/(3*(r0+r1));
}
double Mesh::radialFaceSecondMoment(int i,int face) const {
    double f=fraction[face];
    return f*f*(sq(radius[i])+sq(radius[i+1]))/2;
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
Flow::Flow(Definition d):definition_(d),mesh_(d),medium_(d.medium()),ns_(medium_.size()),nt_(d.turbulence.enabled?2:0),
    nw_(ns_+nt_),inletComposition_(d.massFractions()),totalPressure_(d.totalPressure) {
    auto count=mesh_.cells.size();
    state_.resize(count); stage_.resize(count); next_.resize(count); rhs_.resize(count);
    slopesZ_.resize(count); primitives_.resize(count); radialLow_.resize(count); radialHigh_.resize(count); pressureSource_.resize(count);
    for(auto* v:{&species_,&speciesStage_,&speciesNext_,&speciesRhs_}) v->resize(count*ns_);
    for(auto* v:{&fractions_,&fractionSlopesZ_,&fractionsLow_,&fractionsHigh_}) v->resize(count*nw_);
    for(auto* v:{&turbulence_,&turbulenceStage_,&turbulenceNext_,&turbulenceRhs_,&turbulenceStart_}) v->resize(count*nt_);
    temperature_.assign(count,d.totalTemperature); sound_.resize(count);
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
        const auto& g=d.igniter;
        if(g.energy>0) {
            for(std::size_t q=0;q<count;++q) {
                const auto& c=mesh_.cells[q];
                if(c.z>=g.zMin && c.z<=g.zMax && c.r<=g.rMax) { igniterCells_.push_back(q);igniterVolume_+=c.volume; }
            }
            if(igniterCells_.empty()) throw std::invalid_argument("The igniter region contains no cell centroid.");
        }
    }
    if(nt_) {
        std::vector<bool> plate;
        if(d.experiment==Case::Chamber && !d.wallSlip) for(int s:faceSupply_) plate.push_back(s<0);
        wallDistance_=wallDistance(mesh_,!d.wallSlip,plate);
    }
    if(medium_.hasTransport()) prepareTransport();
    std::vector<double> fractions;
    for(std::size_t q=0;q<count;++q) fractions.insert(fractions.end(),inletComposition_.begin(),inletComposition_.end());
    setInitialState(cells,fractions);
}
double Flow::BoundaryRates::grossMomentum() const {
    return std::abs(inletMomentum)+std::abs(outletMomentum)+std::abs(wallAxial)+std::abs(bodyAxial);
}
Flow::BoundaryRates Flow::BoundaryRates::average(const BoundaryRates& a,const BoundaryRates& b) {
    auto mean=[](double x,double y){return 0.5*(x+y);};
    return {mean(a.mass,b.mass),mean(a.energy,b.energy),mean(a.inlet,b.inlet),mean(a.outlet,b.outlet),
        mean(a.inletMomentum,b.inletMomentum),mean(a.outletMomentum,b.outletMomentum),
        mean(a.wallAxial,b.wallAxial),mean(a.bodyAxial,b.bodyAxial),mean(a.heat,b.heat),mean(a.wallHeat,b.wallHeat)};
}
void Flow::resetAccounting() {
    initialMass_=initialEnergy_=initialMomentum_=0;
    for(std::size_t q=0;q<state_.size();++q) {
        initialMass_+=state_[q][0]*mesh_.cells[q].volume;
        initialMomentum_+=state_[q][1]*mesh_.cells[q].volume;
        initialEnergy_+=state_[q][3]*mesh_.cells[q].volume;
    }
    integratedMassFlux_=integratedEnergyFlux_=integratedMomentumSource_=integratedMomentumGross_=integratedHeat_=0;
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
void Flow::setTotalPressure(double p) {
    if(!(p>0) || !std::isfinite(p)) throw std::invalid_argument("Reservoir pressure must be positive and finite.");
    totalPressure_=p;
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
void Flow::refresh(const std::vector<Conserved>& state,const std::vector<double>& species,const std::vector<double>& turbulence) {
    for(std::size_t q=0;q<state.size();++q) {
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
// h(T) + u^2/2 = h0. If that needs u above the face sound speed the face is choked and carries the
// sonic state of the same stream. g = 0 is a closed valve (a wall).
Conserved Flow::supplyFlux(const Supply& s,double g,const Primitive& in,const double* yIn) const {
    const double* y=s.composition.data();
    const double gasR=medium_.gasConstant(y),h0=medium_.enthalpy(s.totalTemperature,y);
    const double tIn=in.p/(in.rho*medium_.gasConstant(yIn)),aIn=medium_.soundSpeed(tIn,yIn),gIn=aIn*aIn*in.rho/in.p;
    auto characteristic=[&](double u) {
        double base=1+(gIn-1)/(2*aIn)*(u-in.uz);
        if(!(base>0)) throw std::runtime_error("Supply face: the interior gas recedes faster than the supply can follow.");
        return in.p*std::pow(base,2*gIn/(gIn-1));
    };
    if(!(g>0)) return {0,characteristic(0),0,0};
    // Static temperature of the stream at face pressure p (Newton from the stagnation temperature).
    auto staticTemperature=[&](double p) {
        double t=s.totalTemperature;
        for(int n=0;n<60;++n) {
            double u=g*gasR*t/p,f=medium_.enthalpy(t,y)+0.5*u*u-h0;
            double step=f/(medium_.cv(t,y)+gasR+u*u/t);
            t=std::max(0.5*t,t-step);
            if(std::abs(step)<1e-13*t) break;
        }
        return t;
    };
    // Sonic state: h(T) + a(T)^2/2 = h0.
    double tSonic=s.totalTemperature;
    for(int n=0;n<60;++n) {
        double cv=medium_.cv(tSonic,y),gamma=(cv+gasR)/cv;
        double f=medium_.enthalpy(tSonic,y)+0.5*gamma*gasR*tSonic-h0;
        double step=f/(cv+gasR+0.5*gamma*gasR);
        tSonic=std::max(0.5*tSonic,tSonic-step);
        if(std::abs(step)<1e-13*tSonic) break;
    }
    const double aSonic=medium_.soundSpeed(tSonic,y),pSonic=g*gasR*tSonic/aSonic;
    auto residual=[&](double p){ double t=staticTemperature(p);return p-characteristic(g*gasR*t/p); };
    double p=pSonic,u=aSonic;
    if(residual(pSonic)<0) {
        // Subsonic face: the residual increases with p; bracket and bisect.
        double lo=pSonic,hi=2*std::max(pSonic,in.p);
        for(int n=0;residual(hi)<0;++n) { if(n==60) throw std::runtime_error("Supply face pressure not bracketed."); lo=hi;hi*=2; }
        for(int n=0;n<200 && hi-lo>1e-14*hi;++n) { double mid=0.5*(lo+hi);(residual(mid)<0?lo:hi)=mid; }
        p=0.5*(lo+hi);u=g*gasR*staticTemperature(p)/p;
    }
    return {g,g*u+p,0,g*h0};
}
// With transport the diffusive rate 2 nu sum(A^2) / V is added to the convective one, with nu the
// largest of 4/3 mu / rho, lambda / (rho cv) and the mixture diffusion coefficients: pure diffusion
// then steps at cfl/2 of the explicit Euler bound V^2 / (nu sum(A^2)). With turbulence the eddy
// parts are added (mu + mu_t, lambda + cp mu_t / Pr_t, D_km + mu_t / (rho Sc_t); the k and omega
// diffusivities mu + sigma mu_t lie below 4/3 (mu + mu_t)), and the source coefficients of the
// first half of the split are frozen here.
double Flow::stableDt() {
    refresh(state_,species_,turbulence_);
    const bool transport=medium_.hasTransport();
    if(transport) transportProperties();
    if(nt_) { transportGradients();eddyViscosity(true); }
    double dt=std::numeric_limits<double>::infinity();
    for(int i=0;i<mesh_.nz;++i) for(int j=0;j<mesh_.nr;++j) {
        auto q=mesh_.index(i,j);const auto& w=primitives_[q];
        double a=sound_[q],volume=mesh_.cells[q].volume;
        double sumArea2=sq(mesh_.axialArea(i,j))+sq(mesh_.axialArea(i+1,j));
        double rate=(mesh_.axialArea(i,j)+mesh_.axialArea(i+1,j))*(std::abs(w.uz)+a);
        for(int f:{j,j+1}) {
            auto ar=mesh_.radialAreaVector(i,f);
            rate+=std::abs(w.uz*ar[0]+w.ur*ar[1])+a*std::hypot(ar[0],ar[1]);
            sumArea2+=sq(ar[0])+sq(ar[1]);
        }
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
        dt=std::min(dt,definition_.cfl*volume/rate);
    }
    return dt;
}
// Radial reconstruction. Cell averages are r-weighted, so they are point values at the volume
// centroid r_c, not at the cell midpoint; near the axis the two differ by O(dr), which made the
// first-face states first-order. Slopes are therefore taken in physical r between centroids and
// evaluated at each face's area-weighted radius. The axis row uses the parity of the field:
// u_r is odd (u_r = s r, s limited against the slope to the next centroid) and rho, u_z, p are
// even (linear in r^2 about <r^2>). pressureSource_ is the exact integral of p/r dV over the
// reconstructed profile, so a pressure gradient balanced by a body force stays at rest.
void Flow::radialProfiles() {
    const auto& m=mesh_;
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
        auto q=m.index(i,j);const auto& cell=m.cells[q];auto c=values(primitives_[q]);
        double low=m.radialFaceRadius(i,j)-cell.r,high=m.radialFaceRadius(i,j+1)-cell.r;
        Conserved slope{};
        auto centroidSlope=[&](int a,int b,int k){
            return (values(primitives_[m.index(i,b)])[k]-values(primitives_[m.index(i,a)])[k])/(m.cells[m.index(i,b)].r-m.cells[m.index(i,a)].r);};
        if(definition_.secondOrder && j==0) {
            Conserved curvature{};
            if(m.nr>=3) {
                auto moment=[&](int b){return m.cells[m.index(i,b)].radialSecondMoment;};
                auto w1=values(primitives_[m.index(i,1)]),w2=values(primitives_[m.index(i,2)]);
                for(int k:{0,1,3}) curvature[k]=minmod((w1[k]-c[k])/(moment(1)-moment(0)),(w2[k]-w1[k])/(moment(2)-moment(1)));
            }
            if(m.nr>=2) slope[2]=minmod(c[2]/cell.r,centroidSlope(0,1,2));
            Primitive top=unpack(c);
            double face=m.radialFaceSecondMoment(i,1)-cell.radialSecondMoment;
            top.rho+=curvature[0]*face;top.uz+=curvature[1]*face;top.ur+=slope[2]*high;top.p+=curvature[3]*face;
            radialLow_[q]=primitives_[q];radialHigh_[q]=top;
            pressureSource_[q]=c[3]*cell.radialPressureMeasure+curvature[3]*(cell.r*cell.volume-cell.radialSecondMoment*cell.radialPressureMeasure);
            continue;
        }
        if(definition_.secondOrder && j<m.nr-1)
            for(int k=0;k<4;++k) slope[k]=minmod(centroidSlope(j-1,j,k),centroidSlope(j,j+1,k));
        // Wall row: limited one-sided slope (see oneSided above), with density and pressure face
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
    // half the cell value) at the wall.
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
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
}
Flow::BoundaryRates Flow::rhs(const std::vector<Conserved>& state,const std::vector<double>& species,const std::vector<double>& turbulence,
                              std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative,
                              std::vector<double>& turbulenceDerivative,double time) {
    const auto& d=definition_;const auto& m=mesh_;
    std::fill(derivative.begin(),derivative.end(),Conserved{});
    std::fill(speciesDerivative.begin(),speciesDerivative.end(),0.0);
    std::fill(turbulenceDerivative.begin(),turbulenceDerivative.end(),0.0);
    refresh(state,species,turbulence);
    std::fill(slopesZ_.begin(),slopesZ_.end(),Conserved{});
    std::fill(fractionSlopesZ_.begin(),fractionSlopesZ_.end(),0.0);
    if(d.secondOrder) for(int i=1;i<m.nz-1;++i) for(int j=0;j<m.nr;++j) {
        auto q=m.index(i,j); auto c=values(primitives_[q]);
        auto ql=m.index(i-1,j),qr=m.index(i+1,j);
        auto l=values(primitives_[ql]),r=values(primitives_[qr]);
        for(int k=0;k<4;++k) slopesZ_[q][k]=minmod(c[k]-l[k],r[k]-c[k]);
        for(std::size_t k=0;k<nw_;++k)
            fractionSlopesZ_[q*nw_+k]=minmod(fractions_[q*nw_+k]-fractions_[ql*nw_+k],fractions_[qr*nw_+k]-fractions_[q*nw_+k]);
    }
    // Boundary rows/columns have no outer neighbour. A zero slope there makes wall, inlet and
    // outlet face states first-order (measured: wall-row entropy error order ~1.1). Use the
    // limited one-sided slope minmod(c-n1, n1-n2), which is second-order for smooth data and
    // falls back towards zero across a jump. Density/pressure/mass-fraction face values stay >= half the cell value.
    auto oneSided=[&](std::size_t q,std::size_t n1,std::size_t n2,double sign) {
        auto c=values(primitives_[q]),a=values(primitives_[n1]),b=values(primitives_[n2]);
        auto& slope=slopesZ_[q];
        for(int k=0;k<4;++k) slope[k]=sign*minmod(c[k]-a[k],a[k]-b[k]);
        for(int k:{0,3}) slope[k]=std::clamp(slope[k],-c[k],c[k]);
        for(std::size_t k=0;k<nw_;++k) {
            double yc=fractions_[q*nw_+k];
            fractionSlopesZ_[q*nw_+k]=std::clamp(sign*minmod(yc-fractions_[n1*nw_+k],fractions_[n1*nw_+k]-fractions_[n2*nw_+k]),-yc,yc);
        }
    };
    if(d.secondOrder) for(int j=0;j<m.nr;++j) {
        oneSided(m.index(0,j),m.index(1,j),m.index(2,j),-1);
        oneSided(m.index(m.nz-1,j),m.index(m.nz-2,j),m.index(m.nz-3,j),1);
    }
    radialProfiles();
    std::vector<double> yl(nw_),yr(nw_);
    auto reconstructed=[&](std::size_t q,double direction,std::vector<double>& y) {
        auto w=values(primitives_[q]);
        for(int k=0;k<4;++k) w[k]+=direction*0.5*slopesZ_[q][k];
        for(std::size_t k=0;k<nw_;++k) y[k]=fractions_[q*nw_+k]+direction*0.5*fractionSlopesZ_[q*nw_+k];
        return unpack(w);
    };
    // Low-Mach corrections act on interior faces and the mirror (slip) wall; inlet/outlet models keep their own states.
    auto faceFlux=[&](Primitive l,Thermal tl,Primitive r,Thermal tr,double nz,double nr,bool interior) {
        if(d.lowMach==LowMach::HllcLm && interior) return hllcLm(l,tl,r,tr,nz,nr);
        if(d.lowMach==LowMach::Thornber && interior) thornberScale(l,r,tl.a,tr.a);
        return hllc(l,tl,r,tr,nz,nr);
    };
    // Species follow the mass flux with the upwind face composition (Larrouturou, JCP 95, 1991):
    // positivity-preserving, and the species fluxes sum to the mass flux. rho k and rho omega
    // follow it likewise with the upwind face k and omega.
    auto speciesFlux=[&](double massFlux,std::size_t into,std::size_t from,double area,bool addInto,bool addFrom) {
        const auto& y=massFlux>=0?yl:yr;
        for(std::size_t k=0;k<ns_;++k) {
            double f=area*massFlux*y[k];
            if(addFrom) speciesDerivative[from*ns_+k]-=f;
            if(addInto) speciesDerivative[into*ns_+k]+=f;
        }
        for(std::size_t k=0;k<nt_;++k) {
            double f=area*massFlux*y[ns_+k];
            if(addFrom) turbulenceDerivative[from*nt_+k]-=f;
            if(addInto) turbulenceDerivative[into*nt_+k]+=f;
        }
    };
    BoundaryRates rates{};
    // Igniter: constant power over the step (set by step()), shared by the cells by volume.
    if(igniterPower_>0) for(auto q:igniterCells_) derivative[q][3]+=igniterPower_*m.cells[q].volume/igniterVolume_;
    rates.energy+=igniterPower_;rates.heat+=igniterPower_;
    for(int i=0;i<=m.nz;++i) for(int j=0;j<m.nr;++j) {
        auto il=m.index(std::max(0,i-1),j),ir=m.index(std::min(m.nz-1,i),j);
        // Boundary models receive the reconstructed interior face state.
        Primitive l=reconstructed(il,1,yl),r=reconstructed(ir,-1,yr);
        double area=m.axialArea(i,j);
        if(i==0 && d.experiment==Case::Chamber) {
            // Injector face: a supply ring or the closed plate (a slip wall).
            Conserved flux{};
            if(int s=faceSupply_[j];s>=0) {
                const auto& supply=d.supplies[s];
                double g=supply.massFlow*supply.opening(time)/supplyArea_[s];
                flux=supplyFlux(supply,g,r,yr.data());
                for(std::size_t k=0;k<ns_;++k) speciesDerivative[ir*ns_+k]+=area*g*supply.composition[k];
                rates.mass+=area*g;rates.energy+=area*flux[3];rates.inlet+=area*g;rates.inletMomentum+=area*flux[1];
            } else {
                Thermal t=faceThermal(r,yr.data());
                flux=faceFlux(reflect(r,1,0),t,r,t,1,0,true);
                flux[0]=0;flux[3]=0;
                rates.wallAxial+=area*flux[1];
            }
            for(int k=0;k<4;++k) derivative[ir][k]+=area*flux[k];
            continue;
        }
        // Inflow at the reservoir composition; a transmissive inlet (UniformDuct, ShockTube) passes the
        // interior k and omega.
        if(i==0) { l=inlet(r);std::copy(inletComposition_.begin(),inletComposition_.end(),yl.begin());std::copy(yr.begin()+ns_,yr.end(),yl.begin()+ns_); }
        Thermal tl=faceThermal(l,yl.data());
        if(i==m.nz) { bool ambientInflow=false;r=outlet(l,yl.data(),primitives_[il],ambientInflow);yr=ambientInflow?ambient_:yl; }
        Thermal tr=faceThermal(r,yr.data());
        auto flux=faceFlux(l,tl,r,tr,1,0,i>0 && i<m.nz);
        for(int k=0;k<4;++k) {
            if(i>0) derivative[il][k]-=area*flux[k];
            if(i<m.nz) derivative[ir][k]+=area*flux[k];
        }
        speciesFlux(flux[0],ir,il,area,i<m.nz,i>0);
        if(i==0) {
            rates.mass+=area*flux[0];rates.energy+=area*flux[3];rates.inlet+=area*flux[0];rates.inletMomentum+=area*flux[1];
        }
        if(i==m.nz) {
            rates.mass-=area*flux[0];rates.energy-=area*flux[3];rates.outlet+=area*flux[0];rates.outletMomentum+=area*flux[1];
        }
    }
    for(int i=0;i<m.nz;++i) for(int j=1;j<=m.nr;++j) {
        auto il=m.index(i,j-1),ir=m.index(i,std::min(m.nr-1,j));
        auto ar=m.radialAreaVector(i,j);double area=std::hypot(ar[0],ar[1]);
        double nz=ar[0]/area,nr=ar[1]/area;
        Primitive l=radialHigh_[il];
        std::copy_n(fractionsHigh_.begin()+static_cast<std::ptrdiff_t>(il*nw_),nw_,yl.begin());
        Primitive r=l;
        if(j==m.nr) { r=reflect(l,nz,nr);yr=yl; }
        else { r=radialLow_[ir];std::copy_n(fractionsLow_.begin()+static_cast<std::ptrdiff_t>(ir*nw_),nw_,yr.begin()); }
        Thermal tl=faceThermal(l,yl.data()),tr=faceThermal(r,yr.data());
        auto flux=faceFlux(l,tl,r,tr,nz,nr,true);
        // A stationary slip wall has exactly zero mass and energy exchange.
        if(j==m.nr) { flux[0]=0;flux[3]=0; }
        for(int k=0;k<4;++k) {
            derivative[il][k]-=area*flux[k];
            if(j<m.nr) derivative[ir][k]+=area*flux[k];
        }
        speciesFlux(flux[0],ir,il,area,j<m.nr,true);
        if(j==m.nr) rates.wallAxial-=area*flux[1];
    }
    // Volumetric forces (e.g. a Lorentz force): the axial integral is rates.bodyAxial, whose reaction
    // acts on the equipment (coils) in deviceThrust; their work is an external energy exchange.
    if(!bodyForce_.empty()) for(std::size_t q=0;q<state.size();++q) {
        const auto& f=bodyForce_[q];double v=m.cells[q].volume;
        double work=(f[0]*primitives_[q].uz+f[1]*primitives_[q].ur)*v;
        derivative[q][1]+=f[0]*v;derivative[q][2]+=f[1]*v;derivative[q][3]+=work;
        rates.bodyAxial+=f[0]*v;rates.energy+=work;
    }
    if(medium_.hasTransport()) transportFluxes(derivative,speciesDerivative,turbulenceDerivative,rates);
    for(std::size_t q=0;q<state.size();++q) {
        derivative[q][2]+=pressureSource_[q];
        for(double& v:derivative[q]) v/=m.cells[q].volume;
        for(std::size_t k=0;k<ns_;++k) speciesDerivative[q*ns_+k]/=m.cells[q].volume;
        for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]/=m.cells[q].volume;
    }
    return rates;
}
// One step: SSPRK2 on the fluxes. With turbulence the SST sources are Strang-split around it: half a
// step of source (coefficients frozen at the start, from stableDt), the flux step, and half a step
// of source with coefficients at the new state. A rejected attempt restarts from the saved rho k
// and rho omega with half the step.
double Flow::step(double maxDt) {
    if(!(maxDt>0)) throw std::invalid_argument("Step interval must be positive.");
    double dt=std::min(stableDt(),maxDt);
    const double event=nextEvent(time_);
    bool landing=false;
    if(event-time_<=dt) { dt=event-time_;landing=true; }
    const std::size_t n=state_.size();
    const auto& ig=definition_.igniter;
    if(nt_) { turbulenceStart_=turbulence_;sourcesStart_=sources_; }
    for(int attempt=0;attempt<14;++attempt) {
        // Mean igniter power over [t, t + dt]; steps end on its switching times, so this is exact.
        igniterPower_=0;
        if(ig.energy>0) {
            double overlap=std::min(time_+dt,ig.start+ig.duration)-std::max(time_,ig.start);
            if(overlap>0) igniterPower_=ig.energy/ig.duration*overlap/dt;
        }
        bool ok=true;
        if(nt_) {
            turbulence_=turbulenceStart_;
            turbulenceSource(state_,turbulence_,sourcesStart_,0.5*dt);
            for(std::size_t q=0;q<n;++q) ok=ok&&admissible(state_[q],species_.data()+q*ns_,turbulence_.data()+q*nt_);
        }
        BoundaryRates first{},second{};
        if(ok) {
            first=rhs(state_,species_,turbulence_,rhs_,speciesRhs_,turbulenceRhs_,time_);
            for(std::size_t q=0;q<n;++q) {
                for(int k=0;k<4;++k) stage_[q][k]=state_[q][k]+dt*rhs_[q][k];
                for(std::size_t k=q*ns_;k<(q+1)*ns_;++k) speciesStage_[k]=species_[k]+dt*speciesRhs_[k];
                for(std::size_t k=q*nt_;k<(q+1)*nt_;++k) turbulenceStage_[k]=turbulence_[k]+dt*turbulenceRhs_[k];
                ok=ok&&admissible(stage_[q],speciesStage_.data()+q*ns_,turbulenceStage_.data()+q*nt_);
            }
        }
        if(ok) {
            second=rhs(stage_,speciesStage_,turbulenceStage_,rhs_,speciesRhs_,turbulenceRhs_,time_+dt);
            for(std::size_t q=0;q<n;++q) {
                for(int k=0;k<4;++k) next_[q][k]=0.5*(state_[q][k]+stage_[q][k]+dt*rhs_[q][k]);
                for(std::size_t k=q*ns_;k<(q+1)*ns_;++k) speciesNext_[k]=0.5*(species_[k]+speciesStage_[k]+dt*speciesRhs_[k]);
                for(std::size_t k=q*nt_;k<(q+1)*nt_;++k) turbulenceNext_[k]=0.5*(turbulence_[k]+turbulenceStage_[k]+dt*turbulenceRhs_[k]);
                ok=ok&&admissible(next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_);
            }
        }
        if(ok && nt_) {
            refresh(next_,speciesNext_,turbulenceNext_);
            transportProperties();transportGradients();eddyViscosity(true);
            turbulenceSource(next_,turbulenceNext_,sources_,0.5*dt);
            for(std::size_t q=0;q<n;++q) ok=ok&&admissible(next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_);
        }
        if(ok) {
            state_.swap(next_);species_.swap(speciesNext_);turbulence_.swap(turbulenceNext_);
            time_=landing?event:time_+dt;dt_=dt;++steps_;
            lastRates_=BoundaryRates::average(first,second);
            integratedHeat_+=dt*lastRates_.heat;
            integratedMassFlux_+=dt*lastRates_.mass;
            integratedEnergyFlux_+=dt*lastRates_.energy;
            integratedMomentumSource_+=dt*lastRates_.netMomentum();
            integratedMomentumGross_+=dt*lastRates_.grossMomentum();
            return dt;
        }
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
    result.inletMassFlow=lastRates_.inlet;result.outletMassFlow=lastRates_.outlet;
    result.inletMomentumFlux=lastRates_.inletMomentum;result.outletMomentumFlux=lastRates_.outletMomentum;
    result.wallAxialForce=lastRates_.wallAxial;result.bodyAxialForce=lastRates_.bodyAxial;result.wallHeatFlow=lastRates_.wallHeat;
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
        result.igniterPower=lastRates_.heat;result.igniterEnergy=integratedHeat_;
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
    result.momentumBalanceError=(result.axialMomentum-initialMomentum_-integratedMomentumSource_)/
        std::max(std::abs(initialMomentum_)+integratedMomentumGross_,std::numeric_limits<double>::min());
    result.massBalanceError=(result.mass-initialMass_-integratedMassFlux_)/initialMass_;
    result.energyBalanceError=(result.energy-initialEnergy_-integratedEnergyFlux_)/initialEnergy_;
    return result;
}
FieldSnapshot Flow::snapshot() const {
    FieldSnapshot out;
    out.definition=definition_;out.radius=mesh_.radius;out.fraction=mesh_.fraction;out.measurements=measurements();out.appliedTotalPressure=totalPressure_;
    out.cells.reserve(state_.size());for(std::size_t q=0;q<state_.size();++q) out.cells.push_back(cellPrimitive(q));
    return out;
}
} // namespace crucible
