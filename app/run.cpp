#include "core/flow.hpp"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char** argv) {
    try {
        crucible::Definition d;double end=0.01;int steps=0;
        for(int a=1;a<argc;++a) {
            std::string arg=argv[a];
            auto value=[&]() -> const char* {if(++a>=argc) throw std::invalid_argument("Missing option value");return argv[a];};
            if(arg=="--nz") d.nz=std::stoi(value());
            else if(arg=="--nr") d.nr=std::stoi(value());
            else if(arg=="--time") end=std::stod(value());
            else if(arg=="--steps") steps=std::stoi(value());
            else if(arg=="--first-order") d.secondOrder=false;
            else if(arg=="--help") {std::cout<<"crucible_run [--nz 160] [--nr 24] [--time 0.01] [--steps N] [--first-order]\n";return 0;}
            else throw std::invalid_argument("Unknown argument: "+arg);
        }
        auto start=std::chrono::steady_clock::now();
        crucible::Flow flow(d);
        auto initialized=std::chrono::steady_clock::now();
        if(steps>0) for(int n=0;n<steps;++n) flow.step();else flow.advanceTo(end);
        auto done=std::chrono::steady_clock::now();
        auto m=flow.measurements();double seconds=std::chrono::duration<double>(done-initialized).count();
        std::cout<<std::setprecision(10)<<"model=axisymmetric_ideal_gas_euler\n"
          <<"cells="<<d.nz*d.nr<<"\nsteps="<<m.steps<<"\nphysical_time_s="<<m.time
          <<"\ninitialize_s="<<std::chrono::duration<double>(initialized-start).count()<<"\ncompute_s="<<seconds
          <<"\ncell_updates_per_s="<<d.nz*static_cast<double>(d.nr)*m.steps/seconds
          <<"\ninlet_mass_flow_kg_s="<<m.inletMassFlow<<"\noutlet_mass_flow_kg_s="<<m.outletMassFlow
          <<"\nideal_choked_mass_flow_kg_s="<<crucible::chokedMassFlow(d)
          <<"\noutlet_force_N="<<m.outletForce<<"\nexit_mach="<<m.exitMach
          <<"\nmass_budget_relative_error="<<m.massBalanceError
          <<"\nenergy_budget_relative_error="<<m.energyBalanceError<<"\nrejected_steps="<<m.rejectedSteps<<'\n';
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
