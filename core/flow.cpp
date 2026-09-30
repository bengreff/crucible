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
Conserved physicalFlux(Primitive w, double nz, double nr, Gas gas) {
    auto u=conservative(w,gas);
    double un=w.uz*nz+w.ur*nr;
    return {w.rho*un,u[1]*un+w.p*nz,u[2]*un+w.p*nr,(u[3]+w.p)*un};
}
Primitive reflect(Primitive w, double nz, double nr) {
    double un=w.uz*nz+w.ur*nr;
    w.uz-=2*un*nz; w.ur-=2*un*nr;
    return w;
}
bool admissible(const Conserved& u, Gas gas) {
    if (!(u[0]>0)) return false;
    for(double x:u) if(!std::isfinite(x)) return false;
    return (gas.gamma-1)*(u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))>0;
}
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
        // Reference center for reconstruction; all conservation uses exact ring measures.
        cells[index(i,j)]={v,(i+0.5)*dz,(a+b)*(r0+r1)/4,pi*dz*(r0+r1)*(b-a)};
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
Conserved conservative(Primitive w,Gas gas) {
    return {w.rho,w.rho*w.uz,w.rho*w.ur,w.p/(gas.gamma-1)+0.5*w.rho*(sq(w.uz)+sq(w.ur))};
}
Primitive primitive(const Conserved& u,Gas gas) {
    if(!admissible(u,gas)) throw std::runtime_error("Non-admissible gas state: density or internal energy is not positive.");
    return {u[0],u[1]/u[0],u[2]/u[0],(gas.gamma-1)*(u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))};
}
Conserved hllc(Primitive l,Primitive r,double nz,double nr,Gas gas) {
    const auto ul=conservative(l,gas),ur=conservative(r,gas);
    const auto fl=physicalFlux(l,nz,nr,gas),fr=physicalFlux(r,nz,nr,gas);
    double vl=l.uz*nz+l.ur*nr,vr=r.uz*nz+r.ur*nr;
    double al=std::sqrt(gas.gamma*l.p/l.rho),ar=std::sqrt(gas.gamma*r.p/r.rho);
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
    double s=sm>=0?sl:sr,v=sm>=0?vl:vr;
    double pstar=w.p+w.rho*(s-v)*(sm-v);
    if(!(pstar>0) || std::abs(s-sm)<1e-20) return hlle();
    double density=w.rho*(s-v)/(s-sm);
    Conserved star{density,density*(w.uz+(sm-v)*nz),density*(w.ur+(sm-v)*nr),
        ((s-v)*u[3]-w.p*v+pstar*sm)/(s-sm)};
    if(!admissible(star,gas)) return hlle();
    Conserved result{};
    for(int k=0;k<4;++k) result[k]=f[k]+s*(star[k]-u[k]);
    return result;
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
Flow::Flow(Definition d):definition_(d),mesh_(d),totalPressure_(d.totalPressure) {
    auto count=mesh_.cells.size();
    state_.resize(count); stage_.resize(count); next_.resize(count); rhs_.resize(count);
    slopesZ_.resize(count); slopesR_.resize(count); primitives_.resize(count);
    for(int i=0;i<d.nz;++i) {
        Primitive base{d.totalPressure/(d.gas.specificR*d.totalTemperature),0,0,d.totalPressure};
        double radius=(mesh_.radius[i]+mesh_.radius[i+1])/2;
        if(d.experiment==Case::ShockTube)
            base=i<d.nz/2?Primitive{1,0,0,1}:Primitive{0.125,0,0,0.1};
        if(d.experiment==Case::Nozzle) {
            // The 1D preparation is shared by a whole axial column.
            double m=machFromArea(sq(radius/d.throatRadius),(i+0.5)*mesh_.dz>d.length*d.throatFraction,d.gas.gamma);
            double factor=1+(d.gas.gamma-1)*m*m/2;
            double t=d.totalTemperature/factor;
            base.p=d.totalPressure/std::pow(factor,d.gas.gamma/(d.gas.gamma-1));
            base.rho=base.p/(d.gas.specificR*t);
            base.uz=m*std::sqrt(d.gas.gamma*d.gas.specificR*t);
        }
        for(int j=0;j<d.nr;++j) {
            auto q=mesh_.index(i,j);auto w=base;
            if(d.experiment==Case::Nozzle)
                w.ur=w.uz*mesh_.cells[q].r/radius*(mesh_.radius[i+1]-mesh_.radius[i])/mesh_.dz;
            state_[q]=conservative(w,d.gas);
        }
    }
    resetAccounting();
}
void Flow::resetAccounting() {
    initialMass_=initialEnergy_=0;
    for(std::size_t q=0;q<state_.size();++q) {
        initialMass_+=state_[q][0]*mesh_.cells[q].volume;
        initialEnergy_+=state_[q][3]*mesh_.cells[q].volume;
    }
    integratedMassFlux_=integratedEnergyFlux_=0;
    lastRates_={};time_=dt_=0;steps_=rejectedSteps_=0;
}
void Flow::setUniform(Primitive w) {
    auto u=conservative(w,definition_.gas); (void)primitive(u,definition_.gas);
    std::fill(state_.begin(),state_.end(),u);resetAccounting();
}
void Flow::setInitialState(const std::vector<Primitive>& cells) {
    if(cells.size()!=state_.size()) throw std::invalid_argument("Initial field does not match the mesh.");
    std::vector<Conserved> initialized;initialized.reserve(cells.size());
    for(auto w:cells) {
        auto u=conservative(w,definition_.gas);(void)primitive(u,definition_.gas);
        initialized.push_back(u);
    }
    state_.swap(initialized);resetAccounting();
}
void Flow::setTotalPressure(double p) {
    if(!(p>0) || !std::isfinite(p)) throw std::invalid_argument("Reservoir pressure must be positive and finite.");
    totalPressure_=p;
}
Primitive Flow::inlet(Primitive w) const {
    if(definition_.experiment!=Case::Nozzle) return w;
    const double g=definition_.gas.gamma,a0=std::sqrt(g*definition_.gas.specificR*definition_.totalTemperature);
    const double invariant=w.uz-2*std::sqrt(g*w.p/w.rho)/(g-1);
    // Subsonic reservoir boundary: preserve the outgoing acoustic characteristic.
    double lo=0,hi=a0/std::sqrt(1+(g-1)/2);
    const double lowInvariant=-2*a0/(g-1);
    const double highInvariant=hi-2*hi/(g-1);
    if(invariant<lowInvariant-1e-10*a0 || invariant>highInvariant+1e-10*a0)
        throw std::runtime_error("Inlet requires reverse or supersonic flow, outside this reservoir boundary model.");
    for(int n=0;n<40;++n) {
        double v=(lo+hi)/2,a=std::sqrt(a0*a0-(g-1)*v*v/2);
        if(v-2*a/(g-1)>invariant) hi=v; else lo=v;
    }
    double v=(lo+hi)/2,t=definition_.totalTemperature-v*v*(g-1)/(2*g*definition_.gas.specificR);
    double p=totalPressure_*std::pow(t/definition_.totalTemperature,g/(g-1));
    return {p/(definition_.gas.specificR*t),v,0,p};
}
Primitive Flow::outlet(Primitive w) const {
    if(definition_.experiment!=Case::Nozzle) return w;
    double g=definition_.gas.gamma,a=std::sqrt(g*w.p/w.rho);
    if(w.uz>=a) return w; // All characteristics leave a supersonic outlet.
    if(w.uz < -1e-10*a) throw std::runtime_error("Outlet backflow is outside this nozzle prototype's boundary model.");
    double p=definition_.backPressure,rho=w.rho*std::pow(p/w.p,1/g);
    double u=w.uz+2*(a-std::sqrt(g*p/rho))/(g-1);
    if(u < -1e-10*a) throw std::runtime_error("Imposed outlet pressure requires unsupported backflow.");
    return {rho,u,w.ur,p};
}
double Flow::stableDt() {
    double dt=std::numeric_limits<double>::infinity();
    for(int i=0;i<mesh_.nz;++i) for(int j=0;j<mesh_.nr;++j) {
        auto q=mesh_.index(i,j);auto w=primitive(state_[q],definition_.gas);
        double a=std::sqrt(definition_.gas.gamma*w.p/w.rho);
        double rate=(mesh_.axialArea(i,j)+mesh_.axialArea(i+1,j))*(std::abs(w.uz)+a);
        for(int f:{j,j+1}) {
            auto ar=mesh_.radialAreaVector(i,f);
            rate+=std::abs(w.uz*ar[0]+w.ur*ar[1])+a*std::hypot(ar[0],ar[1]);
        }
        dt=std::min(dt,definition_.cfl*mesh_.cells[q].volume/rate);
    }
    return dt;
}
Flow::BoundaryRates Flow::rhs(const std::vector<Conserved>& state,std::vector<Conserved>& derivative) {
    const auto& d=definition_;const auto& m=mesh_;
    std::fill(derivative.begin(),derivative.end(),Conserved{});
    for(std::size_t q=0;q<state.size();++q) primitives_[q]=primitive(state[q],d.gas);
    std::fill(slopesZ_.begin(),slopesZ_.end(),Conserved{});
    std::fill(slopesR_.begin(),slopesR_.end(),Conserved{});
    if(d.secondOrder) for(int i=1;i<m.nz-1;++i) for(int j=0;j<m.nr;++j) {
        auto q=m.index(i,j); auto c=values(primitives_[q]);
        auto l=values(primitives_[m.index(i-1,j)]),r=values(primitives_[m.index(i+1,j)]);
        for(int k=0;k<4;++k) slopesZ_[q][k]=minmod(c[k]-l[k],r[k]-c[k]);
    }
    if(d.secondOrder) for(int i=0;i<m.nz;++i) for(int j=1;j<m.nr-1;++j) {
        auto q=m.index(i,j); auto c=values(primitives_[q]);
        auto l=values(primitives_[m.index(i,j-1)]),r=values(primitives_[m.index(i,j+1)]);
        for(int k=0;k<4;++k) slopesR_[q][k]=minmod(c[k]-l[k],r[k]-c[k]);
    }
    // Boundary rows/columns have no outer neighbour. A zero slope there makes wall, inlet and
    // outlet face states first-order (measured: wall-row entropy error order ~1.1). Use the
    // limited one-sided slope minmod(c-n1, n1-n2), which is second-order for smooth data and
    // falls back towards zero across a jump. Density/pressure face values stay >= half the cell value.
    auto oneSided=[&](std::size_t q,std::size_t n1,std::size_t n2,double sign,Conserved& slope) {
        auto c=values(primitives_[q]),a=values(primitives_[n1]),b=values(primitives_[n2]);
        for(int k=0;k<4;++k) slope[k]=sign*minmod(c[k]-a[k],a[k]-b[k]);
        for(int k:{0,3}) slope[k]=std::clamp(slope[k],-c[k],c[k]);
    };
    if(d.secondOrder) for(int j=0;j<m.nr;++j) {
        oneSided(m.index(0,j),m.index(1,j),m.index(2,j),-1,slopesZ_[m.index(0,j)]);
        oneSided(m.index(m.nz-1,j),m.index(m.nz-2,j),m.index(m.nz-3,j),1,slopesZ_[m.index(m.nz-1,j)]);
    }
    if(d.secondOrder && m.nr>=3) for(int i=0;i<m.nz;++i)
        oneSided(m.index(i,m.nr-1),m.index(i,m.nr-2),m.index(i,m.nr-3),1,slopesR_[m.index(i,m.nr-1)]);
    auto reconstructed=[&](std::size_t q,bool axial,double direction) {
        auto w=values(primitives_[q]); const auto& slope=axial?slopesZ_[q]:slopesR_[q];
        for(int k=0;k<4;++k) w[k]+=direction*0.5*slope[k];
        return unpack(w);
    };
    BoundaryRates rates{};
    for(int i=0;i<=m.nz;++i) for(int j=0;j<m.nr;++j) {
        auto il=m.index(std::max(0,i-1),j),ir=m.index(std::min(m.nz-1,i),j);
        // Boundary models receive the reconstructed interior face state.
        Primitive l=reconstructed(il,true,1),r=reconstructed(ir,true,-1);
        if(i==0) l=inlet(r);
        if(i==m.nz) r=outlet(l);
        auto flux=hllc(l,r,1,0,d.gas);double area=m.axialArea(i,j);
        for(int k=0;k<4;++k) {
            if(i>0) derivative[il][k]-=area*flux[k];
            if(i<m.nz) derivative[ir][k]+=area*flux[k];
        }
        if(i==0) { rates.mass+=area*flux[0];rates.energy+=area*flux[3];rates.inlet+=area*flux[0]; }
        if(i==m.nz) {
            rates.mass-=area*flux[0];rates.energy-=area*flux[3];rates.outlet+=area*flux[0];
            rates.outletForce+=area*(flux[1]-d.backPressure);
        }
    }
    for(int i=0;i<m.nz;++i) for(int j=1;j<=m.nr;++j) {
        auto il=m.index(i,j-1),ir=m.index(i,std::min(m.nr-1,j));
        auto ar=m.radialAreaVector(i,j);double area=std::hypot(ar[0],ar[1]);
        double nz=ar[0]/area,nr=ar[1]/area;
        Primitive l=reconstructed(il,false,1),r=j==m.nr?reflect(l,nz,nr):reconstructed(ir,false,-1);
        auto flux=hllc(l,r,nz,nr,d.gas);
        // A stationary slip wall has exactly zero mass and energy exchange.
        if(j==m.nr) { flux[0]=0;flux[3]=0; }
        for(int k=0;k<4;++k) {
            derivative[il][k]-=area*flux[k];
            if(j<m.nr) derivative[ir][k]+=area*flux[k];
        }
    }
    for(std::size_t q=0;q<state.size();++q) {
        derivative[q][2]+=primitives_[q].p*m.cells[q].radialPressureMeasure;
        for(double& v:derivative[q]) v/=m.cells[q].volume;
    }
    return rates;
}
double Flow::step(double maxDt) {
    if(!(maxDt>0)) throw std::invalid_argument("Step interval must be positive.");
    double dt=std::min(stableDt(),maxDt);
    for(int attempt=0;attempt<14;++attempt) {
        auto first=rhs(state_,rhs_);
        bool ok=true;
        for(std::size_t q=0;q<state_.size();++q) {
            for(int k=0;k<4;++k) stage_[q][k]=state_[q][k]+dt*rhs_[q][k];
            ok=ok&&admissible(stage_[q],definition_.gas);
        }
        BoundaryRates second{};
        if(ok) {
            second=rhs(stage_,rhs_);
            for(std::size_t q=0;q<state_.size();++q) {
                for(int k=0;k<4;++k) next_[q][k]=0.5*(state_[q][k]+stage_[q][k]+dt*rhs_[q][k]);
                ok=ok&&admissible(next_[q],definition_.gas);
            }
        }
        if(ok) {
            state_.swap(next_);time_+=dt;dt_=dt;++steps_;
            integratedMassFlux_+=0.5*dt*(first.mass+second.mass);
            integratedEnergyFlux_+=0.5*dt*(first.energy+second.energy);
            lastRates_={0.5*(first.mass+second.mass),0.5*(first.energy+second.energy),
                0.5*(first.inlet+second.inlet),0.5*(first.outlet+second.outlet),0.5*(first.outletForce+second.outletForce)};
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
    result.inletMassFlow=lastRates_.inlet;result.outletMassFlow=lastRates_.outlet;result.outletForce=lastRates_.outletForce;
    result.minPressure=std::numeric_limits<double>::infinity();
    double exitArea=0;
    for(int i=0;i<mesh_.nz;++i) for(int j=0;j<mesh_.nr;++j) {
        auto q=mesh_.index(i,j);auto w=primitive(state_[q],definition_.gas);
        double mach=std::hypot(w.uz,w.ur)/std::sqrt(definition_.gas.gamma*w.p/w.rho);
        result.mass+=state_[q][0]*mesh_.cells[q].volume;
        result.energy+=state_[q][3]*mesh_.cells[q].volume;
        result.minPressure=std::min(result.minPressure,w.p);result.maxPressure=std::max(result.maxPressure,w.p);
        result.maxMach=std::max(result.maxMach,mach);
        if(i==mesh_.nz-1) { double a=mesh_.axialArea(mesh_.nz,j);exitArea+=a;result.exitMach+=a*mach; }
    }
    result.exitMach/=exitArea;
    result.massBalanceError=(result.mass-initialMass_-integratedMassFlux_)/initialMass_;
    result.energyBalanceError=(result.energy-initialEnergy_-integratedEnergyFlux_)/initialEnergy_;
    return result;
}
FieldSnapshot Flow::snapshot() const {
    FieldSnapshot out;
    out.definition=definition_;out.radius=mesh_.radius;out.measurements=measurements();out.appliedTotalPressure=totalPressure_;
    out.cells.reserve(state_.size());for(const auto& u:state_) out.cells.push_back(primitive(u,definition_.gas));
    return out;
}
} // namespace crucible
