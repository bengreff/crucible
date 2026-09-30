// Grid-convergence study of the default nozzle. The exact steady inviscid solution is
// isentropic with uniform total enthalpy, so entropy and total-enthalpy deviations have an
// exact zero reference. The quasi-1D choked mass flow is exact to Hall's 2D discharge
// correction, (gamma+1)/(96 (Rc/rt)^2) < 3e-5 for this throat (Rc = 0.63 m, rt = 0.02 m).
#include "core/flow.hpp"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

using namespace crucible;
int main(int argc,char** argv) {
    int base=40,levels=4;std::vector<double> times{0.008};bool first=false;
    for(int a=1;a<argc;++a) {
        std::string arg=argv[a];
        if(arg=="--base" && a+1<argc) base=std::stoi(argv[++a]);
        else if(arg=="--levels" && a+1<argc) levels=std::stoi(argv[++a]);
        else if(arg=="--time" && a+1<argc) times={std::stod(argv[++a])};
        else if(arg=="--times" && a+1<argc) {times.clear();std::stringstream list(argv[++a]);std::string t;
            while(std::getline(list,t,',')) times.push_back(std::stod(t));}
        else if(arg=="--first-order") first=true;
        else {std::fprintf(stderr,"crucible_convergence [--base 40] [--levels 4] [--time 0.008 | --times 0.004,0.008] [--first-order]\n");return 1;}
    }
    std::printf("nz,nr,time_s,steps,compute_s,mdot_out_rel_err,mdot_in_rel_err,entropy_L1,entropy_wall_L1,entropy_axis_L1,enthalpy_L1,exit_mach,outlet_force_N,mass_res,energy_res\n");
    for(int level=0;level<levels;++level) {
        Definition d;d.nz=base<<level;d.nr=d.nz*3/20;d.secondOrder=!first;
        Flow flow(d);auto start=std::chrono::steady_clock::now();
        for(double end:times) {flow.advanceTo(end);
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const double g=d.gas.gamma,rho0=d.totalPressure/(d.gas.specificR*d.totalTemperature);
        const double h0=g/(g-1)*d.gas.specificR*d.totalTemperature;
        double volume=0,entropy=0,enthalpy=0,wallVolume=0,wall=0,axisVolume=0,axis=0;
        const auto& m=flow.mesh();
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
            auto q=m.index(i,j);auto w=primitive(flow.state()[q],d.gas);double v=m.cells[q].volume;
            double s=std::abs(std::log(w.p/d.totalPressure)-g*std::log(w.rho/rho0));
            double h=std::abs((g/(g-1)*w.p/w.rho+0.5*(w.uz*w.uz+w.ur*w.ur))/h0-1);
            volume+=v;entropy+=v*s;enthalpy+=v*h;
            if(j==m.nr-1){wallVolume+=v;wall+=v*s;}
            if(j==0){axisVolume+=v;axis+=v*s;}
        }
        auto r=flow.measurements();double ideal=chokedMassFlow(d);
        std::printf("%d,%d,%.6g,%llu,%.3f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.8f,%.6f,%.2e,%.2e\n",d.nz,d.nr,r.time,
            static_cast<unsigned long long>(r.steps),seconds,r.outletMassFlow/ideal-1,r.inletMassFlow/ideal-1,
            entropy/volume,wall/wallVolume,axis/axisVolume,enthalpy/volume,r.exitMach,r.outletForce,
            r.massBalanceError,r.energyBalanceError);
        std::fflush(stdout);
    }}
}
