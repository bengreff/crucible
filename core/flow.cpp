#include "core/flow.hpp"
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

Medium Definition::medium() const { return species.empty()?Medium::perfectGas(gas):Medium(species); }
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
    if(experiment==Case::Nozzle && (throatRadius>inletRadius || throatRadius>exitRadius))
        throw std::invalid_argument("The throat must be no wider than inlet and exit.");
    auto y=massFractions();
    double sum=0;
    for(double v:y) { if(!(v>=0) || !std::isfinite(v)) throw std::invalid_argument("Mass fractions must be finite and non-negative."); sum+=v; }
    if(y.size()!=(species.empty()?1:species.size()) || std::abs(sum-1)>1e-12)
        throw std::invalid_argument("Composition must give one mass fraction per species, summing to one.");
}

Mesh::Mesh(const Definition& d):nz(d.nz),nr(d.nr),dz(d.length/d.nz) {
    d.validate();
    radius.resize(nz+1);
    for(int i=0;i<=nz;++i) {
        double x=static_cast<double>(i)/nz;
        if(d.experiment!=Case::Nozzle) radius[i]=d.inletRadius;
        else if(x<d.throatFraction) {
            double t=x/d.throatFraction;
            radius[i]=d.throatRadius+(d.inletRadius-d.throatRadius)*(1+std::cos(pi*t))/2;
        } else {
            double t=(x-d.throatFraction)/(1-d.throatFraction);
            radius[i]=d.throatRadius+(d.exitRadius-d.throatRadius)*(1-std::cos(pi*t))/2;
        }
    }
    cells.resize(static_cast<std::size_t>(nz)*nr);
    for(int i=0;i<nz;++i) for(int j=0;j<nr;++j) {
        double a=static_cast<double>(j)/nr,b=static_cast<double>(j+1)/nr;
        double r0=radius[i],r1=radius[i+1];
        double v=pi*dz/3*(r0*r0+r0*r1+r1*r1)*(b*b-a*a);
        // Exact radial moments of the frustum ring (wall radius linear in z across the cell).
        double first=2*pi/3*(b*b*b-a*a*a)*dz*(r0+r1)*(r0*r0+r1*r1)/4;
        double second=pi/2*(b*b*b*b-a*a*a*a)*dz*(r0*r0*r0*r0+r0*r0*r0*r1+r0*r0*r1*r1+r0*r1*r1*r1+r1*r1*r1*r1)/5;
        cells[index(i,j)]={v,(i+0.5)*dz,first/v,pi*dz*(r0+r1)*(b-a),second/v};
    }
}
double Mesh::axialArea(int face,int j) const {
    return pi*sq(radius[face])*(2*j+1)/sq(nr);
}
std::array<double,2> Mesh::radialAreaVector(int i,int face) const {
    double f=static_cast<double>(face)/nr;
    double r0=f*radius[i],r1=f*radius[i+1];
    return {-pi*(r0+r1)*(r1-r0),pi*(r0+r1)*dz};
}
double Mesh::radialFaceRadius(int i,int face) const {
    double f=static_cast<double>(face)/nr,r0=radius[i],r1=radius[i+1];
    return face==0?0:2*f*(r0*r0+r0*r1+r1*r1)/(3*(r0+r1));
}
double Mesh::radialFaceSecondMoment(int i,int face) const {
    double f=static_cast<double>(face)/nr;
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
Flow::Flow(Definition d):definition_(d),mesh_(d),medium_(d.medium()),ns_(medium_.size()),
    inletComposition_(d.massFractions()),totalPressure_(d.totalPressure) {
    auto count=mesh_.cells.size();
    state_.resize(count); stage_.resize(count); next_.resize(count); rhs_.resize(count);
    slopesZ_.resize(count); primitives_.resize(count); radialLow_.resize(count); radialHigh_.resize(count); pressureSource_.resize(count);
    for(auto* v:{&species_,&speciesStage_,&speciesNext_,&speciesRhs_,&fractions_,&fractionSlopesZ_,&fractionsLow_,&fractionsHigh_})
        v->resize(count*ns_);
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
        mean(a.wallAxial,b.wallAxial),mean(a.bodyAxial,b.bodyAxial)};
}
void Flow::resetAccounting() {
    initialMass_=initialEnergy_=initialMomentum_=0;
    for(std::size_t q=0;q<state_.size();++q) {
        initialMass_+=state_[q][0]*mesh_.cells[q].volume;
        initialMomentum_+=state_[q][1]*mesh_.cells[q].volume;
        initialEnergy_+=state_[q][3]*mesh_.cells[q].volume;
    }
    integratedMassFlux_=integratedEnergyFlux_=integratedMomentumSource_=integratedMomentumGross_=0;
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
    std::vector<double> partial(y.size());
    for(std::size_t q=0;q<cells.size();++q) {
        const auto& w=cells[q];const double* yq=y.data()+q*ns_;
        if(!(w.rho>0) || !(w.p>0)) throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
        double t=w.p/(w.rho*medium_.gasConstant(yq));
        initialized[q]=conservative(w,medium_.internalEnergy(t,yq));
        temperature_[q]=t;
        for(std::size_t k=0;k<ns_;++k) partial[q*ns_+k]=w.rho*yq[k];
        if(!admissible(initialized[q],partial.data()+q*ns_))
            throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
    }
    state_.swap(initialized);species_.swap(partial);resetAccounting();
}
void Flow::setMassFractions(std::size_t q,const double* y) {
    if(q>=state_.size()) throw std::out_of_range("Cell index out of range.");
    std::vector<double> partial(ns_);
    for(std::size_t k=0;k<ns_;++k) partial[k]=state_[q][0]*y[k];
    if(!admissible(state_[q],partial.data())) throw std::runtime_error("Composition gives a non-admissible state.");
    std::copy(partial.begin(),partial.end(),species_.begin()+static_cast<std::ptrdiff_t>(q*ns_));
}
void Flow::setPartialDensities(const std::vector<double>& partial) {
    if(partial.size()!=species_.size()) throw std::invalid_argument("Partial densities do not match the mesh.");
    for(std::size_t q=0;q<state_.size();++q)
        if(!admissible(state_[q],partial.data()+q*ns_)) throw std::runtime_error("Partial densities give a non-admissible state.");
    species_=partial;
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
bool Flow::admissible(const Conserved& u,const double* partial) const {
    // energyFloor is linear in the mass fractions, so it is evaluated on the clipped partial densities.
    double sum=0;
    for(std::size_t k=0;k<ns_;++k) { if(!std::isfinite(partial[k])) return false; sum+=std::max(partial[k],0.0); }
    if(!(sum>0)) return false;
    return admissibleBulk(u,medium_.energyFloorOfClipped(partial)/sum);
}
std::vector<double> Flow::massFractions(std::size_t q) const {
    std::vector<double> y(ns_);
    normalise(species_.data()+q*ns_,ns_,y.data());
    return y;
}
double Flow::temperature(std::size_t q) const {
    auto y=massFractions(q);const auto& u=state_[q];
    return medium_.temperature((u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))/u[0],y.data(),temperature_[q]);
}
Primitive Flow::cellPrimitive(std::size_t q) const {
    auto y=massFractions(q);const auto& u=state_[q];
    double t=medium_.temperature((u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))/u[0],y.data(),temperature_[q]);
    if(!std::isfinite(t) || !admissibleBulk(u,medium_.energyFloor(y.data())))
        throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
    return {u[0],u[1]/u[0],u[2]/u[0],u[0]*medium_.gasConstant(y.data())*t};
}
// Cell primitives, mass fractions, temperatures and frozen sound speeds of a state.
void Flow::refresh(const std::vector<Conserved>& state,const std::vector<double>& species) {
    for(std::size_t q=0;q<state.size();++q) {
        const auto& u=state[q];double* y=fractions_.data()+q*ns_;
        normalise(species.data()+q*ns_,ns_,y);
        if(!admissibleBulk(u,medium_.energyFloor(y)))
            throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
        double t=medium_.temperature((u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))/u[0],y,temperature_[q]);
        auto props=medium_.properties(t,y);
        temperature_[q]=t;sound_[q]=std::sqrt((props.cv+props.r)/props.cv*props.r*t);
        primitives_[q]={u[0],u[1]/u[0],u[2]/u[0],u[0]*props.r*t};
    }
}
// Face thermodynamics from reconstructed (rho, p) and mass fractions; y is clipped and normalised in place.
Thermal Flow::faceThermal(const Primitive& w,double* y) const {
    normalise(y,ns_,y);
    double t=w.p/(w.rho*medium_.gasConstant(y));
    auto props=medium_.properties(t,y);
    return {props.e,std::sqrt((props.cv+props.r)/props.cv*props.r*t),medium_.energyFloor(y)};
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
Primitive Flow::outlet(Primitive w,const double* y) const {
    if(definition_.experiment!=Case::Nozzle) return w;
    double t=w.p/(w.rho*medium_.gasConstant(y)),a=medium_.soundSpeed(t,y),g=a*a*w.rho/w.p;
    if(w.uz>=a) return w; // All characteristics leave a supersonic outlet.
    if(w.uz < -1e-10*a) throw std::runtime_error("Outlet backflow is outside this nozzle prototype's boundary model.");
    double p=definition_.backPressure,rho=w.rho*std::pow(p/w.p,1/g);
    double u=w.uz+2*(a-std::sqrt(g*p/rho))/(g-1);
    if(u < -1e-10*a) throw std::runtime_error("Imposed outlet pressure requires unsupported backflow.");
    return {rho,u,w.ur,p};
}
double Flow::stableDt() {
    refresh(state_,species_);
    double dt=std::numeric_limits<double>::infinity();
    for(int i=0;i<mesh_.nz;++i) for(int j=0;j<mesh_.nr;++j) {
        auto q=mesh_.index(i,j);const auto& w=primitives_[q];
        double a=sound_[q];
        double rate=(mesh_.axialArea(i,j)+mesh_.axialArea(i+1,j))*(std::abs(w.uz)+a);
        for(int f:{j,j+1}) {
            auto ar=mesh_.radialAreaVector(i,f);
            rate+=std::abs(w.uz*ar[0]+w.ur*ar[1])+a*std::hypot(ar[0],ar[1]);
        }
        dt=std::min(dt,definition_.cfl*mesh_.cells[q].volume/rate);
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
    // Mass fractions are even in r and reconstructed like density: r^2 curvature on the axis row,
    // limited centroid slopes inside, the limited one-sided slope (face >= half the cell value) at the wall.
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
        auto q=m.index(i,j);const auto& cell=m.cells[q];
        double low=m.radialFaceRadius(i,j)-cell.r,high=m.radialFaceRadius(i,j+1)-cell.r;
        auto y=[&](int b,std::size_t k){return fractions_[m.index(i,b)*ns_+k];};
        auto centroidSlope=[&](int a,int b,std::size_t k){return (y(b,k)-y(a,k))/(m.cells[m.index(i,b)].r-m.cells[m.index(i,a)].r);};
        for(std::size_t k=0;k<ns_;++k) {
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
            fractionsLow_[q*ns_+k]=lowValue;fractionsHigh_[q*ns_+k]=highValue;
        }
    }
}
Flow::BoundaryRates Flow::rhs(const std::vector<Conserved>& state,const std::vector<double>& species,
                              std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative) {
    const auto& d=definition_;const auto& m=mesh_;
    std::fill(derivative.begin(),derivative.end(),Conserved{});
    std::fill(speciesDerivative.begin(),speciesDerivative.end(),0.0);
    refresh(state,species);
    std::fill(slopesZ_.begin(),slopesZ_.end(),Conserved{});
    std::fill(fractionSlopesZ_.begin(),fractionSlopesZ_.end(),0.0);
    if(d.secondOrder) for(int i=1;i<m.nz-1;++i) for(int j=0;j<m.nr;++j) {
        auto q=m.index(i,j); auto c=values(primitives_[q]);
        auto ql=m.index(i-1,j),qr=m.index(i+1,j);
        auto l=values(primitives_[ql]),r=values(primitives_[qr]);
        for(int k=0;k<4;++k) slopesZ_[q][k]=minmod(c[k]-l[k],r[k]-c[k]);
        for(std::size_t k=0;k<ns_;++k)
            fractionSlopesZ_[q*ns_+k]=minmod(fractions_[q*ns_+k]-fractions_[ql*ns_+k],fractions_[qr*ns_+k]-fractions_[q*ns_+k]);
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
        for(std::size_t k=0;k<ns_;++k) {
            double yc=fractions_[q*ns_+k];
            fractionSlopesZ_[q*ns_+k]=std::clamp(sign*minmod(yc-fractions_[n1*ns_+k],fractions_[n1*ns_+k]-fractions_[n2*ns_+k]),-yc,yc);
        }
    };
    if(d.secondOrder) for(int j=0;j<m.nr;++j) {
        oneSided(m.index(0,j),m.index(1,j),m.index(2,j),-1);
        oneSided(m.index(m.nz-1,j),m.index(m.nz-2,j),m.index(m.nz-3,j),1);
    }
    radialProfiles();
    std::vector<double> yl(ns_),yr(ns_);
    auto reconstructed=[&](std::size_t q,double direction,std::vector<double>& y) {
        auto w=values(primitives_[q]);
        for(int k=0;k<4;++k) w[k]+=direction*0.5*slopesZ_[q][k];
        for(std::size_t k=0;k<ns_;++k) y[k]=fractions_[q*ns_+k]+direction*0.5*fractionSlopesZ_[q*ns_+k];
        return unpack(w);
    };
    // Low-Mach corrections act on interior faces and the mirror (slip) wall; inlet/outlet models keep their own states.
    auto faceFlux=[&](Primitive l,Thermal tl,Primitive r,Thermal tr,double nz,double nr,bool interior) {
        if(d.lowMach==LowMach::HllcLm && interior) return hllcLm(l,tl,r,tr,nz,nr);
        if(d.lowMach==LowMach::Thornber && interior) thornberScale(l,r,tl.a,tr.a);
        return hllc(l,tl,r,tr,nz,nr);
    };
    // Species follow the mass flux with the upwind face composition (Larrouturou, JCP 95, 1991):
    // positivity-preserving, and the species fluxes sum to the mass flux.
    auto speciesFlux=[&](double massFlux,std::size_t into,std::size_t from,double area,bool addInto,bool addFrom) {
        const auto& y=massFlux>=0?yl:yr;
        for(std::size_t k=0;k<ns_;++k) {
            double f=area*massFlux*y[k];
            if(addFrom) speciesDerivative[from*ns_+k]-=f;
            if(addInto) speciesDerivative[into*ns_+k]+=f;
        }
    };
    BoundaryRates rates{};
    for(int i=0;i<=m.nz;++i) for(int j=0;j<m.nr;++j) {
        auto il=m.index(std::max(0,i-1),j),ir=m.index(std::min(m.nz-1,i),j);
        // Boundary models receive the reconstructed interior face state.
        Primitive l=reconstructed(il,1,yl),r=reconstructed(ir,-1,yr);
        if(i==0) { l=inlet(r);yl=inletComposition_; }
        Thermal tl=faceThermal(l,yl.data());
        if(i==m.nz) { r=outlet(l,yl.data());yr=yl; }
        Thermal tr=faceThermal(r,yr.data());
        auto flux=faceFlux(l,tl,r,tr,1,0,i>0 && i<m.nz);double area=m.axialArea(i,j);
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
        std::copy_n(fractionsHigh_.begin()+static_cast<std::ptrdiff_t>(il*ns_),ns_,yl.begin());
        Primitive r=l;
        if(j==m.nr) { r=reflect(l,nz,nr);yr=yl; }
        else { r=radialLow_[ir];std::copy_n(fractionsLow_.begin()+static_cast<std::ptrdiff_t>(ir*ns_),ns_,yr.begin()); }
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
    for(std::size_t q=0;q<state.size();++q) {
        derivative[q][2]+=pressureSource_[q];
        for(double& v:derivative[q]) v/=m.cells[q].volume;
        for(std::size_t k=0;k<ns_;++k) speciesDerivative[q*ns_+k]/=m.cells[q].volume;
    }
    return rates;
}
double Flow::step(double maxDt) {
    if(!(maxDt>0)) throw std::invalid_argument("Step interval must be positive.");
    double dt=std::min(stableDt(),maxDt);
    const std::size_t n=state_.size();
    for(int attempt=0;attempt<14;++attempt) {
        auto first=rhs(state_,species_,rhs_,speciesRhs_);
        bool ok=true;
        for(std::size_t q=0;q<n;++q) {
            for(int k=0;k<4;++k) stage_[q][k]=state_[q][k]+dt*rhs_[q][k];
            for(std::size_t k=q*ns_;k<(q+1)*ns_;++k) speciesStage_[k]=species_[k]+dt*speciesRhs_[k];
            ok=ok&&admissible(stage_[q],speciesStage_.data()+q*ns_);
        }
        BoundaryRates second{};
        if(ok) {
            second=rhs(stage_,speciesStage_,rhs_,speciesRhs_);
            for(std::size_t q=0;q<n;++q) {
                for(int k=0;k<4;++k) next_[q][k]=0.5*(state_[q][k]+stage_[q][k]+dt*rhs_[q][k]);
                for(std::size_t k=q*ns_;k<(q+1)*ns_;++k) speciesNext_[k]=0.5*(species_[k]+speciesStage_[k]+dt*speciesRhs_[k]);
                ok=ok&&admissible(next_[q],speciesNext_.data()+q*ns_);
            }
        }
        if(ok) {
            state_.swap(next_);species_.swap(speciesNext_);time_+=dt;dt_=dt;++steps_;
            lastRates_=BoundaryRates::average(first,second);
            integratedMassFlux_+=dt*lastRates_.mass;
            integratedEnergyFlux_+=dt*lastRates_.energy;
            integratedMomentumSource_+=dt*lastRates_.netMomentum();
            integratedMomentumGross_+=dt*lastRates_.grossMomentum();
            return dt;
        }
        ++rejectedSteps_;dt*=0.5;
    }
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
    result.wallAxialForce=lastRates_.wallAxial;result.bodyAxialForce=lastRates_.bodyAxial;
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
    out.definition=definition_;out.radius=mesh_.radius;out.measurements=measurements();out.appliedTotalPressure=totalPressure_;
    out.cells.reserve(state_.size());for(std::size_t q=0;q<state_.size();++q) out.cells.push_back(cellPrimitive(q));
    return out;
}
} // namespace crucible
