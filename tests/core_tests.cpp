#include "core/flow.hpp"
#include "core/session.hpp"
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <thread>
#include <string>
#include <utility>
#include <vector>

using namespace crucible;
// --low-mach=thornber|hllclm reruns every check with that flux variant (docs/LOW_MACH.md).
LowMach testScheme=LowMach::None;
Definition base() {Definition d;d.lowMach=testScheme;return d;}
double sq(double x) {return x*x;}
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
    Definition d=base();d.nz=n;d.nr=2;d.length=1;d.experiment=Case::ShockTube;
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
    Definition d=base();d.experiment=Case::UniformDuct;d.nz=4;d.nr=nr;
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

// Axis-row (max-norm) error of the same mode: u_r is odd in r, so a midpoint/zero-slope axis
// reconstruction gives an O(dr) face error in the mass flux, an O(1) divergence error in that row.
double radialWaveAxisError(int nr) {
    Definition d=base();d.experiment=Case::UniformDuct;d.nz=4;d.nr=nr;
    Flow flow(d);std::vector<Primitive> initial;
    constexpr double epsilon=1e-5,p0=100000,rho0=1,root=3.8317059702075125;
    double sound=std::sqrt(d.gas.gamma*p0/rho0),k=root/d.inletRadius;
    for(auto c:flow.mesh().cells) {double s=epsilon*bessel(0,k*c.r);initial.push_back({rho0*(1+s/d.gas.gamma),0,0,p0*(1+s)});}
    flow.setInitialState(initial);flow.advanceTo(std::numbers::pi/(2*sound*k));
    auto q=flow.mesh().index(1,0);auto w=primitive(flow.state()[q],d.gas);
    // The r-weighted average of J1(kr) over the axis cell (independent quadrature).
    double rt=d.inletRadius/nr,num=0,den=0;
    for(int n=0;n<400;++n){double r=(n+.5)*rt/400;num+=bessel(1,k*r)*r;den+=r;}
    return std::abs(w.ur-epsilon*p0/(rho0*sound)*num/den)/(epsilon*p0/(rho0*sound));
}
// Current loop (radius a, current I) field at (z,r) relative to the loop plane; AGM elliptic integrals.
std::array<double,2> loopField(double z,double r,double a,double current) {
    constexpr double mu0=1.25663706212e-6;double m=4*a*r/(sq(a+r)+z*z);
    double x=1,y=std::sqrt(1-m),sum=0,power=.5,big=1,e=1-m/2;
    for(int n=0;n<40;++n){double an=(x+y)/2,cn=(x-y)/2;y=std::sqrt(x*y);x=an;power*=2;sum+=power*cn*cn;}
    big=std::numbers::pi/(2*x);e=big*(1-m/2-sum);
    double base=mu0*current/(2*std::numbers::pi*std::sqrt(sq(a+r)+z*z)),d2=sq(a-r)+z*z;
    double bz=base*(big+(a*a-r*r-z*z)/d2*e),br=r>0?base*z/r*(-big+(a*a+r*r+z*z)/d2*e):0;
    return {bz,br};
}
// Gas at rest in a coil-like body force: F = -grad(B^2/2mu0) of a 0.3 T current loop (radius
// 0.05 m) around a straight duct; exact equilibrium p = p0 - B^2/2mu0, uniform density. Returns
// the initial acceleration residual {max all, max axis row, volume L1}, normalized by max|F|/rho.
std::array<double,5> coilRestResidual(int nz) {
    Definition d=base();d.experiment=Case::UniformDuct;d.nz=nz;d.nr=nz*3/20;
    constexpr double mu0=1.25663706212e-6,a=.05,zc=.3,rho=1.16,p0=1e6;
    double current=2*a*.3/mu0;
    auto magneticPressure=[&](double z,double r){auto b=loopField(z-zc,r,a,current);return (sq(b[0])+sq(b[1]))/(2*mu0);};
    Flow flow(d);const auto& m=flow.mesh();std::vector<Primitive> initial;std::vector<std::array<double,2>> force;
    double scale=0;
    for(auto c:m.cells) {
        // Cell averages: equilibrium pressure by quadrature over the ring; force at the centroid
        // by central differences (h = 1e-7 m; error ~1e-9 of the force).
        initial.push_back({rho,0,0,0});
        double h=1e-7;
        std::array<double,2> f{-(magneticPressure(c.z+h,c.r)-magneticPressure(c.z-h,c.r))/(2*h),
                               -(magneticPressure(c.z,c.r+h)-magneticPressure(c.z,c.r-h))/(2*h)};
        force.push_back(f);scale=std::max(scale,std::hypot(f[0],f[1])/rho);
    }
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
        double z0=i*m.dz,ra=m.radius[i]*j/m.nr,rb=m.radius[i]*(j+1)/m.nr,sum=0,weight=0;
        for(int u=0;u<8;++u) for(int v=0;v<8;++v){double z=z0+(u+.5)*m.dz/8,r=ra+(v+.5)*(rb-ra)/8;
            sum+=r*magneticPressure(z,r);weight+=r;}
        initial[m.index(i,j)].p=p0-sum/weight;
    }
    flow.setInitialState(initial);flow.setBodyForce(force);
    constexpr double dt=1e-9;flow.step(dt);
    // Residual acceleration in units of the largest applied force per mass: [0] max |a|,
    // [1] max |a| on the axis row, [2] volume L1 of |a|, [3] max |a_r|, [4] max |a_r| on the axis row.
    std::array<double,5> e{};double volume=0;
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j){auto q=m.index(i,j);auto w=primitive(flow.state()[q],d.gas);
        double accel=std::hypot(w.uz,w.ur)/dt/scale,radial=std::abs(w.ur)/dt/scale;
        e[0]=std::max(e[0],accel);e[3]=std::max(e[3],radial);
        if(j==0){e[1]=std::max(e[1],accel);e[4]=std::max(e[4],radial);}
        e[2]+=accel*m.cells[q].volume;volume+=m.cells[q].volume;}
    e[2]/=volume;
    return e;
}

// Independent exact Riemann solver for gamma-law gas (Toro, ch. 4); state is (rho,u,p).
struct Exact { double rho,u,p; };
Exact exactRiemann(Exact l,Exact r,double xi,double g=1.4) {
    double al=std::sqrt(g*l.p/l.rho),ar=std::sqrt(g*r.p/r.rho),g6=(g-1)/(g+1);
    auto f=[&](double p,const Exact& k,double a) {
        if(p>k.p) return (p-k.p)*std::sqrt(2/((g+1)*k.rho)/(p+g6*k.p));
        return 2*a/(g-1)*(std::pow(p/k.p,(g-1)/(2*g))-1);
    };
    double lo=1e-14,hi=1e7;
    for(int n=0;n<300;++n){double p=std::sqrt(lo*hi);if(f(p,l,al)+f(p,r,ar)+r.u-l.u>0)hi=p;else lo=p;}
    double ps=std::sqrt(lo*hi),us=0.5*(l.u+r.u)+0.5*(f(ps,r,ar)-f(ps,l,al));
    if(xi<=us) {
        if(ps>l.p) {
            double s=l.u-al*std::sqrt((g+1)/(2*g)*ps/l.p+(g-1)/(2*g));
            return xi<=s?l:Exact{l.rho*(ps/l.p+g6)/(g6*ps/l.p+1),us,ps};
        }
        double head=l.u-al,tail=us-al*std::pow(ps/l.p,(g-1)/(2*g));
        if(xi<=head) return l;
        if(xi>tail) return {l.rho*std::pow(ps/l.p,1/g),us,ps};
        double c=2/(g+1)+g6/al*(l.u-xi);
        return {l.rho*std::pow(c,2/(g-1)),2/(g+1)*(al+(g-1)/2*l.u+xi),l.p*std::pow(c,2*g/(g-1))};
    }
    if(ps>r.p) {
        double s=r.u+ar*std::sqrt((g+1)/(2*g)*ps/r.p+(g-1)/(2*g));
        return xi>=s?r:Exact{r.rho*(ps/r.p+g6)/(g6*ps/r.p+1),us,ps};
    }
    double head=r.u+ar,tail=us+ar*std::pow(ps/r.p,(g-1)/(2*g));
    if(xi>=head) return r;
    if(xi<tail) return {r.rho*std::pow(ps/r.p,1/g),us,ps};
    double c=2/(g+1)-g6/ar*(r.u-xi);
    return {r.rho*std::pow(c,2/(g-1)),2/(g+1)*(-ar+(g-1)/2*r.u+xi),r.p*std::pow(c,2*g/(g-1))};
}
// Axial Riemann problem in a straight duct with nr radial rows. Returns the density L1 error
// against the exact solution; also requires identical rows and exactly zero radial velocity.
double riemannError(Exact l,Exact r,int n,double t,int nr=2) {
    Definition d=base();d.nz=n;d.nr=nr;d.length=1;d.experiment=Case::ShockTube;Flow f(d);
    std::vector<Primitive> initial;
    for(int i=0;i<n;++i) for(int j=0;j<nr;++j){auto w=i<n/2?l:r;initial.push_back({w.rho,w.u,0,w.p});}
    f.setInitialState(initial);f.advanceTo(t);double error=0;
    double rowDeviation=0,radialSpeed=0;
    for(int i=0;i<n;++i) {
        auto first=primitive(f.state()[f.mesh().index(i,0)],d.gas);
        error+=std::abs(first.rho-exactRiemann(l,r,((i+.5)/n-.5)/t).rho)/n;
        for(int j=0;j<nr;++j) {
            auto w=primitive(f.state()[f.mesh().index(i,j)],d.gas);
            rowDeviation=std::max(rowDeviation,std::abs(w.rho/first.rho-1));
            radialSpeed=std::max(radialSpeed,std::abs(w.ur));
        }
    }
    std::cout<<"Planar rows: max density deviation "<<rowDeviation<<", max |u_r| "<<radialSpeed<<'\n';
    require(radialSpeed<1e-9,"Planar axial flow must not create radial velocity");
    if(testScheme==LowMach::None) require(rowDeviation<1e-12,"Planar axial flow must be identical in every radial row");
    auto m=f.measurements();
    require(std::abs(m.massBalanceError)<1e-12 && std::abs(m.energyBalanceError)<1e-12,"Riemann-problem budgets");
    return error;
}
// Volume-weighted entropy deviation from the reservoir, overall and in the wall row.
std::pair<double,double> nozzleEntropyError(int nz,int nr,double t) {
    Definition d=base();d.nz=nz;d.nr=nr;Flow f(d);f.advanceTo(t);
    double g=d.gas.gamma,rho0=d.totalPressure/(d.gas.specificR*d.totalTemperature),all=0,volume=0,wall=0,wallVolume=0;
    for(int i=0;i<nz;++i) for(int j=0;j<nr;++j){auto q=f.mesh().index(i,j);auto w=primitive(f.state()[q],d.gas);
        double v=f.mesh().cells[q].volume,s=std::abs(std::log(w.p/d.totalPressure)-g*std::log(w.rho/rho0));
        all+=v*s;volume+=v;if(j==nr-1){wall+=v*s;wallVolume+=v;}}
    return {all/volume,wall/wallVolume};
}
// Steady subsonic quasi-1D exit Mach number for exit static pressure pe and mass flow mdot.
double exitMachFor(double mdot,double pe,double area,const Definition& d) {
    double g=d.gas.gamma,lo=1e-6,hi=1;
    auto flow=[&](double m){return pe*area*m*std::sqrt(g/(d.gas.specificR*d.totalTemperature)*(1+(g-1)/2*m*m));};
    for(int n=0;n<100;++n){double m=(lo+hi)/2;if(flow(m)>mdot)hi=m;else lo=m;}
    return (lo+hi)/2;
}
// Nozzle started from reservoir rest with a subsonic back pressure; the outlet is a static-pressure boundary.
Measurements fromRest(int nz,double backRatio,double t) {
    Definition d=base();d.nz=nz;d.nr=nz*3/20;d.backPressure=backRatio*d.totalPressure;Flow f(d);
    f.setUniform({d.totalPressure/(d.gas.specificR*d.totalTemperature),0,0,d.totalPressure});f.advanceTo(t);
    auto m=f.measurements();
    require(std::abs(m.massBalanceError)<1e-12 && std::abs(m.energyBalanceError)<1e-12,"Subsonic-outlet budgets");
    return m;
}
// Fully subsonic venturi: isentropic quasi-1D exit state at the imposed back pressure.
// At exit Mach ~0.15 mass flow amplifies total-pressure error by ~1/(gamma M^2) ~ 33.
void venturi() {
    Definition d=base();d.backPressure=0.985*d.totalPressure;double g=d.gas.gamma;
    double me=std::sqrt(2/(g-1)*(std::pow(1/0.985,(g-1)/g)-1)),te=d.totalTemperature/(1+(g-1)/2*me*me);
    double ideal=d.backPressure/(d.gas.specificR*te)*me*std::sqrt(g*d.gas.specificR*te)*std::numbers::pi*d.exitRadius*d.exitRadius;
    double coarse=fromRest(40,.985,.16).outletMassFlow/ideal-1,fine=fromRest(80,.985,.16).outletMassFlow/ideal-1;
    std::cout<<"Subsonic venturi mass-flow error: 40="<<coarse<<" 80="<<fine<<" order="<<std::log2(coarse/fine)<<'\n';
    require(std::log2(coarse/fine)>1.7 && std::abs(fine)<0.03,"Subsonic venturi must converge to isentropic quasi-1D flow");
}
int main(int argc,char** argv) {
    try {
        for(int a=1;a<argc;++a) {
            std::string arg=argv[a];
            if(arg=="--low-mach=thornber") testScheme=LowMach::Thornber;
            else if(arg=="--low-mach=hllclm") testScheme=LowMach::HllcLm;
        }
        bool slow=false;for(int a=1;a<argc;++a) slow|=std::string(argv[a])=="--slow";
        if(slow) {venturi();std::cout<<"Slow verification checks passed.\n";return 0;}
        {Definition d=base();d.experiment=Case::UniformDuct;d.nz=12;d.nr=5;Mesh m(d);double v=0;
         for(auto c:m.cells) v+=c.volume;
         near(v,std::numbers::pi*d.inletRadius*d.inletRadius*d.length,1e-16,"Cylinder volume");}
        {Primitive w{1.3,20,-2,90000};Gas g;auto q=primitive(conservative(w,g),g);
         near(q.p,w.p,1e-9,"State round trip");auto f=hllc(w,w,1,0,g);
         near(f[0],w.rho*w.uz,1e-12,"Identical-state mass flux");}
        {Definition d=base();d.nz=20;d.nr=6;d.experiment=Case::UniformDuct;Flow f(d);
         Primitive w{1,100,0,100000};f.setUniform(w);for(int n=0;n<30;++n) f.step();
         for(auto u:f.state()) {auto p=primitive(u,d.gas);near(p.rho,1,1e-13,"Uniform density");near(p.ur,0,1e-10,"Axis pressure cancellation");}}
        {Definition d=base();d.nz=24;d.nr=8;d.backPressure=d.totalPressure;Flow f(d);
         f.setUniform({d.totalPressure/(d.gas.specificR*d.totalTemperature),0,0,d.totalPressure});
         for(int n=0;n<30;++n) f.step();
         for(auto u:f.state()) {auto p=primitive(u,d.gas);near(p.uz,0,1e-7,"Curved-wall rest axial balance");near(p.ur,0,1e-7,"Curved-wall rest radial balance");}}
        // The exact solver reproduces the independently hard-coded Sod star state.
        {auto star=exactRiemann({1,0,1},{.125,0,.1},0.5);near(star.p,0.303130178050647,1e-12,"Exact solver Sod p*");
         near(star.u,0.92745262004895,1e-12,"Exact solver Sod u*");
         near(exactRiemann({1,0,1},{.125,0,.1},1.5).rho,sodDensity(.5+.15*1.5,.15),1e-12,"Exact solver matches Sod");}
        // HLLC: conservation symmetry, rotation invariance, supersonic upwinding, exact stationary contact.
        {Gas g;Primitive a{1.2,40,-15,90000},b{0.7,-20,30,40000};double c=0.6,s=0.8;
         auto f=hllc(a,b,c,s,g),back=hllc(b,a,-c,-s,g);
         for(int k=0;k<4;++k) near(f[k],-back[k],1e-9*std::abs(f[k])+1e-9,"HLLC must be antisymmetric under side/normal exchange");
         auto rot=[&](Primitive w){return Primitive{w.rho,c*w.uz+s*w.ur,-s*w.uz+c*w.ur,w.p};};
         auto aligned=hllc(rot(a),rot(b),1,0,g);
         near(f[0],aligned[0],1e-9*std::abs(f[0]),"HLLC rotation: mass");
         near(f[1],c*aligned[1]-s*aligned[2],1e-9*std::abs(f[1]),"HLLC rotation: axial momentum");
         near(f[2],s*aligned[1]+c*aligned[2],1e-9*std::abs(f[2]),"HLLC rotation: radial momentum");
         near(f[3],aligned[3],1e-9*std::abs(f[3]),"HLLC rotation: energy");
         Primitive fast{1,900,0,100000};auto up=hllc(fast,{0.5,850,5,80000},1,0,g);
         near(up[0],900,1e-12,"Supersonic HLLC takes the upstream flux");
         auto contact=hllc({1,0,7,1e5},{0.1,0,-3,1e5},1,0,g);
         near(contact[0],0,1e-10,"HLLC stationary contact: no mass flux");near(contact[1],1e5,1e-9,"HLLC stationary contact: pressure flux");
         near(contact[3],0,1e-6,"HLLC stationary contact: no energy flux (roundoff |S|E eps)");}
        {Definition d=base();d.nz=40;d.nr=3;d.length=1;d.experiment=Case::ShockTube;Flow f(d);std::vector<Primitive> initial;
         for(int i=0;i<40;++i) for(int j=0;j<3;++j) initial.push_back({i<20?1.0:0.1,0,0,1e5});
         f.setInitialState(initial);for(int n=0;n<200;++n) f.step();
         for(int i=0;i<40;++i) for(int j=0;j<3;++j){auto w=primitive(f.state()[f.mesh().index(i,j)],d.gas);
             near(w.rho,i<20?1.0:0.1,1e-12,"Stationary contact must be preserved");near(w.uz,0,1e-9,"Stationary contact velocity");}}
        // Exact-solution Riemann problems in a multi-row duct (Toro 2009, tests 1-3).
        {double c=riemannError({1,0,1},{.125,0,.1},100,.15,5),fi=riemannError({1,0,1},{.125,0,.1},400,.15,5);
         std::cout<<"Planar Sod (5 rows) density L1: 100="<<c<<" 400="<<fi<<'\n';
         require(fi<c*0.45,"Planar Sod in radial rows must converge");}
        {double c=riemannError({1,-2,.4},{1,2,.4},100,.15),fi=riemannError({1,-2,.4},{1,2,.4},400,.15);
         std::cout<<"Toro 123 strong rarefaction density L1: 100="<<c<<" 400="<<fi<<'\n';
         require(fi<c*0.6 && fi<0.03,"Strong rarefaction must stay admissible and converge");}
        {double c=riemannError({1,0,1000},{1,0,.01},100,.012),fi=riemannError({1,0,1000},{1,0,.01},400,.012);
         std::cout<<"Toro strong shock (p ratio 1e5) density L1: 100="<<c<<" 400="<<fi<<'\n';
         require(fi<c*0.6,"Strong shock must converge");}
        // Nozzle entropy error (exact steady solution is isentropic): second order including the wall row.
        {auto coarse=nozzleEntropyError(40,6,.004),fine=nozzleEntropyError(80,12,.004);
         double all=std::log2(coarse.first/fine.first),wall=std::log2(coarse.second/fine.second);
         std::cout<<"Nozzle entropy error order: volume="<<all<<" wall row="<<wall<<'\n';
         require(all>1.8 && wall>1.8,"Nozzle entropy error must converge at second order, including the wall");}
        // Choked nozzle with a normal shock in the divergent section (back pressure 0.7 p0).
        // Quasi-1D: choked mass flow and subsonic exit at the back pressure fix the exit Mach number.
        {Definition d=base();auto m=fromRest(40,.7,.02);double ideal=chokedMassFlow(d);
         double me=exitMachFor(ideal,.7*d.totalPressure,std::numbers::pi*d.exitRadius*d.exitRadius,d);
         std::cout<<"Internal shock: mdot/choked-1="<<m.outletMassFlow/ideal-1<<" exit M="<<m.exitMach<<" quasi-1D="<<me<<" max M="<<m.maxMach<<'\n';
         require(std::abs(m.outletMassFlow/ideal-1)<.003 && std::abs(m.exitMach/me-1)<.005 && m.maxMach>1.5,"Internal normal shock must match quasi-1D exit state");}
        // Device thrust: the axial momentum ledger closes, and the gradual nozzle approaches the
        // quasi-1D ideal thrust coefficient at matched p0, pa and area ratio.
        {Definition d=base();d.nz=80;d.nr=12;Flow f(d);f.advanceTo(.008);auto m=f.measurements();double g=d.gas.gamma;
         require(std::abs(m.momentumBalanceError)<1e-13,"Axial momentum budget must close");
         near(m.deviceThrust,m.inletMomentumFlux+m.wallAxialForce-m.ambientAxialForce,1e-9,"Device thrust composition");
         double me=machFromArea(sq(d.exitRadius/d.throatRadius),true,g),pe=d.totalPressure*std::pow(1+(g-1)/2*me*me,-g/(g-1));
         double at=std::numbers::pi*sq(d.throatRadius),ae=std::numbers::pi*sq(d.exitRadius);
         double ideal=std::sqrt(2*g*g/(g-1)*std::pow(2/(g+1),(g+1)/(g-1))*(1-std::pow(pe/d.totalPressure,(g-1)/g)))+(pe-d.backPressure)*ae/(d.totalPressure*at);
         double cf=m.deviceThrust/(d.totalPressure*at),steady=(m.deviceThrust-m.exitPlaneThrust)/m.deviceThrust;
         std::cout<<"Thrust coefficient: device="<<cf<<" quasi-1D ideal="<<ideal<<" ratio-1="<<cf/ideal-1
                  <<" (device-exit)/device="<<steady<<" momentum residual="<<m.momentumBalanceError<<'\n';
         require(std::abs(cf/ideal-1)<.01 && std::abs(steady)<1e-3,"Gradual-nozzle device thrust must approach quasi-1D ideal");}
        {Definition d=base();d.nz=40;d.nr=3;d.length=1;d.experiment=Case::ShockTube;Flow f(d);f.advanceTo(.1);
         require(std::abs(f.measurements().momentumBalanceError)<1e-13,"Momentum budget from rest (Sod)");}
        // Exact frustum-ring moments against independent midpoint quadrature on a nozzle cell.
        {Definition d=base();d.nz=20;d.nr=5;Mesh m(d);auto c=m.cells[m.index(7,2)];double r0=m.radius[7],r1=m.radius[8],v=0,first=0,second=0;
         for(int u=0;u<400;++u){double z=(u+.5)/400,R=r0+(r1-r0)*z;for(int w=0;w<400;++w){double r=R*(2+(w+.5)/400)/5,dv=r*R;
             v+=dv;first+=dv*r;second+=dv*r*r;}}
         near(c.r,first/v,1e-7*c.r,"Cell centroid radius");near(c.radialSecondMoment,second/v,1e-6*c.radialSecondMoment,"Cell r^2 moment");}
        // Gas at rest under a coil-like body force F = -grad(B^2/2mu0) balanced by p = p0 - B^2/2mu0
        // must stay at rest; the residual must vanish under refinement everywhere, including the axis.
        {std::array<double,5> coarse{},fine{};
         for(int n:{40,80,160,320}){auto e=coilRestResidual(n);coarse=fine;fine=e;
             std::cout<<"Coil-force rest residual nz="<<n<<": max="<<e[0]<<" axis row="<<e[1]<<" L1="<<e[2]<<" radial max="<<e[3]<<" radial axis row="<<e[4]<<'\n';}
         require(fine[0]<coarse[0]*.7 && fine[1]<coarse[1]*.6 && fine[4]<coarse[4]*.6 && fine[2]<coarse[2]*.3,
                 "Body-force equilibrium residual must converge, including its maximum and the axis row");}
        // The axis row of the radial acoustic mode exposes reconstruction about the cell centroid:
        // midpoint-referenced slopes made it first order (measured 0.038, 0.013, 0.0044).
        {double coarse=radialWaveAxisError(16),fine=radialWaveAxisError(48);
         std::cout<<"Radial acoustic axis-row error: 16="<<coarse<<" 48="<<fine<<" 144="<<radialWaveAxisError(144)<<'\n';
         require(fine<coarse/9 && fine<5e-4,"Axis-row radial acoustic error must converge at least at second order");}
        double radialCoarse=radialWaveError(16),radialFine=radialWaveError(48);
        std::cout<<"Radial acoustic normalized velocity error: 16="<<radialCoarse<<" 48="<<radialFine<<'\n';
        require(radialFine<radialCoarse*.6 && radialFine<.02,"Radial acoustic mode must converge");
        double coarse=shockError(80),fine=shockError(240);
        std::cout<<"Sod density L1: 80="<<coarse<<" 240="<<fine<<'\n';
        require(fine<coarse*0.65 && fine<0.015,"Sod must converge to independent analytic solution");
        {Definition d=base();d.nz=64;d.nr=8;Flow f(d);f.advanceTo(.004);auto m=f.measurements();
         std::cout<<"Nozzle mdot="<<m.outletMassFlow<<" reference="<<chokedMassFlow(d)<<" exit M="<<m.exitMach<<'\n';
         require(std::abs(m.massBalanceError)<1e-11 && std::abs(m.energyBalanceError)<1e-11,"Nozzle budgets");
         require(m.exitMach>1 && std::abs(m.outletMassFlow/chokedMassFlow(d)-1)<0.02,"Nozzle ideal mass-flow limit");
         f.setTotalPressure(330000);f.advanceTo(.012);auto changed=f.measurements();
         require(std::abs(changed.outletMassFlow/(1.1*chokedMassFlow(d))-1)<0.02,"Pressure change must alter physical outlet flow");
         require(std::abs(changed.massBalanceError)<1e-11 && std::abs(changed.energyBalanceError)<1e-11,"Live-control budgets");}
        {Definition d=base();d.nz=32;d.nr=6;Flow f(d);
         for(double pressure:{240000.0,360000.0,300000.0}) {
             f.setTotalPressure(pressure);f.advanceTo(f.time()+.006);
             auto m=f.measurements();
             require(std::abs(m.outletMassFlow/(chokedMassFlow(d)*pressure/d.totalPressure)-1)<.025,"UI pressure range must settle consistently");
         }}
        {Definition d=base();d.nz=8;d.nr=2;Flow f(d);f.setUniform({1,-10,0,100000});bool rejected=false;
         try{f.step();}catch(const std::runtime_error&){rejected=true;}
         require(rejected,"Unsupported reverse boundary flow must fail explicitly");
         near(f.time(),0,0,"Rejected boundary must not advance time");}
        {Definition d=base();d.nz=24;d.nr=4;Session session(d);
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
