#include "core/flow.hpp"
#include "core/session.hpp"
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <thread>

using namespace crucible;
void require(bool condition,const char* message) {if(!condition) throw std::runtime_error(message);}
void near(double actual,double expected,double tolerance,const char* message) {
    if(std::abs(actual-expected)>tolerance) {
        std::cerr<<message<<": actual="<<actual<<" expected="<<expected<<" tolerance="<<tolerance<<'\n';
        throw std::runtime_error(message);
    }
}
void waitFor(const std::function<bool()>& predicate) {
    auto deadline=std::chrono::steady_clock::now()+std::chrono::seconds(8);
    while(!predicate()) {
        if(std::chrono::steady_clock::now()>deadline) throw std::runtime_error("Worker timeout");
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
}
// Independent analytic Sod solution: gamma=1.4, left (rho,p)=(1,1), right=(.125,.1).
double sodDensity(double x,double t) {
    constexpr double g=1.4,pstar=0.303130178050647,ustar=0.92745262004895;
    double xi=(x-0.5)/t,al=std::sqrt(g),ar=std::sqrt(g*0.1/0.125);
    double astar=al*std::pow(pstar,(g-1)/(2*g));
    double shock=ar*std::sqrt((g+1)/(2*g)*pstar/0.1+(g-1)/(2*g));
    if(xi<-al) return 1;
    if(xi<ustar-astar) return std::pow(2/(g+1)+(g-1)/(g+1)*(-xi/al),2/(g-1));
    if(xi<ustar) return std::pow(pstar,1/g);
    if(xi<shock) return 0.125*(pstar/0.1+(g-1)/(g+1))/((g-1)/(g+1)*pstar/0.1+1);
    return .125;
}
double shockError(int n) {
    Definition d;d.nz=n;d.nr=2;d.length=1;d.experiment=Case::ShockTube;
    Flow f(d);f.advanceTo(.15);double error=0;
    for(int i=0;i<n;++i) error+=std::abs(f.state()[f.mesh().index(i,0)][0]-sodDensity((i+.5)/n,.15))/n;
    auto m=f.measurements();require(std::abs(m.massBalanceError)<1e-12,"Sod mass conservation");
    require(std::abs(m.energyBalanceError)<1e-12,"Sod energy conservation");
    return error;
}
// Small-amplitude cylindrical acoustic eigenmode. J1(kR)=0 imposes a rigid wall.
// This independent reference tests radial transport and the cylindrical source together.
double bessel(int order,double x) {
    double term=order==0?1:x/2,sum=term;
    for(int k=1;k<30;++k){term*=-(x*x/4)/(k*(k+order));sum+=term;}
    return sum;
}
double radialWaveError(int nr) {
    Definition d;d.experiment=Case::UniformDuct;d.nz=4;d.nr=nr;
    Flow flow(d);std::vector<Primitive> initial;
    constexpr double epsilon=1e-5,p0=100000,rho0=1,root=3.8317059702075125;
    double sound=std::sqrt(d.gas.gamma*p0/rho0),k=root/d.inletRadius;
    for(auto c:flow.mesh().cells) {
        double perturbation=epsilon*bessel(0,k*c.r);
        initial.push_back({rho0*(1+perturbation/d.gas.gamma),0,0,p0*(1+perturbation)});
    }
    flow.setInitialState(initial);flow.advanceTo(std::numbers::pi/(2*sound*k));
    double amplitude=epsilon*p0/(rho0*sound),error=0;
    for(int j=0;j<nr;++j){auto q=flow.mesh().index(1,j);auto w=primitive(flow.state()[q],d.gas);
        double expected=amplitude*bessel(1,k*flow.mesh().cells[q].r);
        error+=std::abs(w.ur-expected)/amplitude*(2*j+1)/(nr*nr);}
    auto m=flow.measurements();require(std::abs(m.massBalanceError)<1e-12 && std::abs(m.energyBalanceError)<1e-12,"Radial acoustic budgets");
    return error;
}
int main() {
    try {
        {Definition d;d.experiment=Case::UniformDuct;d.nz=12;d.nr=5;Mesh m(d);double v=0;
         for(auto c:m.cells) v+=c.volume;
         near(v,std::numbers::pi*d.inletRadius*d.inletRadius*d.length,1e-16,"Cylinder volume");}
        {Primitive w{1.3,20,-2,90000};Gas g;auto q=primitive(conservative(w,g),g);
         near(q.p,w.p,1e-9,"State round trip");auto f=hllc(w,w,1,0,g);
         near(f[0],w.rho*w.uz,1e-12,"Identical-state mass flux");}
        {Definition d;d.nz=20;d.nr=6;d.experiment=Case::UniformDuct;Flow f(d);
         Primitive w{1,100,0,100000};f.setUniform(w);for(int n=0;n<30;++n) f.step();
         for(auto u:f.state()) {auto p=primitive(u,d.gas);near(p.rho,1,1e-13,"Uniform density");near(p.ur,0,1e-10,"Axis pressure cancellation");}}
        {Definition d;d.nz=24;d.nr=8;d.backPressure=d.totalPressure;Flow f(d);
         f.setUniform({d.totalPressure/(d.gas.specificR*d.totalTemperature),0,0,d.totalPressure});
         for(int n=0;n<30;++n) f.step();
         for(auto u:f.state()) {auto p=primitive(u,d.gas);near(p.uz,0,1e-7,"Curved-wall rest axial balance");near(p.ur,0,1e-7,"Curved-wall rest radial balance");}}
        double radialCoarse=radialWaveError(16),radialFine=radialWaveError(48);
        std::cout<<"Radial acoustic normalized velocity error: 16="<<radialCoarse<<" 48="<<radialFine<<'\n';
        require(radialFine<radialCoarse*.6 && radialFine<.02,"Radial acoustic mode must converge");
        double coarse=shockError(80),fine=shockError(240);
        std::cout<<"Sod density L1: 80="<<coarse<<" 240="<<fine<<'\n';
        require(fine<coarse*0.65 && fine<0.015,"Sod must converge to independent analytic solution");
        {Definition d;d.nz=64;d.nr=8;Flow f(d);f.advanceTo(.004);auto m=f.measurements();
         std::cout<<"Nozzle mdot="<<m.outletMassFlow<<" reference="<<chokedMassFlow(d)<<" exit M="<<m.exitMach<<'\n';
         require(std::abs(m.massBalanceError)<1e-11 && std::abs(m.energyBalanceError)<1e-11,"Nozzle budgets");
         require(m.exitMach>1 && std::abs(m.outletMassFlow/chokedMassFlow(d)-1)<0.02,"Nozzle ideal mass-flow limit");
         f.setTotalPressure(330000);f.advanceTo(.012);auto changed=f.measurements();
         require(std::abs(changed.outletMassFlow/(1.1*chokedMassFlow(d))-1)<0.02,"Pressure change must alter physical outlet flow");
         require(std::abs(changed.massBalanceError)<1e-11 && std::abs(changed.energyBalanceError)<1e-11,"Live-control budgets");}
        {Definition d;d.nz=32;d.nr=6;Flow f(d);
         for(double pressure:{240000.0,360000.0,300000.0}) {
             f.setTotalPressure(pressure);f.advanceTo(f.time()+.006);
             auto m=f.measurements();
             require(std::abs(m.outletMassFlow/(chokedMassFlow(d)*pressure/d.totalPressure)-1)<.025,"UI pressure range must settle consistently");
         }}
        {Definition d;d.nz=8;d.nr=2;Flow f(d);f.setUniform({1,-10,0,100000});bool rejected=false;
         try{f.step();}catch(const std::runtime_error&){rejected=true;}
         require(rejected,"Unsupported reverse boundary flow must fail explicitly");
         near(f.time(),0,0,"Rejected boundary must not advance time");}
        {Definition d;d.nz=24;d.nr=4;Session session(d);
         waitFor([&]{return session.status()==RunState::Paused;});auto original=session.latest();
         session.run();waitFor([&]{return session.latest()->measurements.steps>10;});
         auto sequence=session.setTotalPressure(330000);
         waitFor([&]{return session.latest()->appliedControlSequence==sequence;});
         session.pause();waitFor([&]{return session.status()==RunState::Paused;});
         auto held=session.latest();std::this_thread::sleep_for(std::chrono::milliseconds(25));
         near(session.latest()->measurements.time,held->measurements.time,0,"Pause stability");
         near(original->measurements.time,0,0,"Snapshot immutability");
         require(session.controls().size()==1,"Control applied once");
         Flow replay(d);
         for(auto event:session.controls()){replay.advanceTo(event.time);replay.setTotalPressure(event.totalPressure);}
         replay.advanceTo(held->measurements.time);auto repeated=replay.snapshot();
         for(std::size_t k=0;k<held->cells.size();++k)
             near(repeated.cells[k].p,held->cells[k].p,held->cells[k].p*1e-9,"Accepted control history must replay independently of wall time");}
        std::cout<<"All core verification checks passed.\n";
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
