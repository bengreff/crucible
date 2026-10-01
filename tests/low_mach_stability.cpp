// Item-3 stability probe (docs/LOW_MACH.md): planar Sod in a 5-row duct; rows must stay identical.
// Prints the max relative row density deviation and |u_r| over time for HLLC and Thornber.
#include "core/flow.hpp"
#include <cstdio>
#include <cmath>
using namespace crucible;
int main(){
  for(auto s:{LowMach::None,LowMach::Thornber}) for(int n:{400,800}){
    Definition d;d.nz=n;d.nr=5;d.length=1;d.experiment=Case::ShockTube;d.lowMach=s;Flow f(d);
    std::vector<Primitive> init;
    for(int i=0;i<n;++i)for(int j=0;j<5;++j)init.push_back(i<n/2?Primitive{1,0,0,1}:Primitive{.125,0,0,.1});
    f.setInitialState(init);
    for(double t:{0.05,0.1,0.15,0.2,0.3,0.45}){f.advanceTo(t);
      double dev=0,ur=0;int at=0;
      for(int i=0;i<n;++i)for(int j=0;j<5;++j){auto a=primitive(f.state()[f.mesh().index(i,0)],d.gas),b=primitive(f.state()[f.mesh().index(i,j)],d.gas);
        double e=std::abs(b.rho/a.rho-1);if(e>dev){dev=e;at=i;} ur=std::max(ur,std::abs(b.ur));}
      printf("scheme %d n=%d t=%.2f rowdev %.2e at x=%.3f  max|ur| %.2e\n",(int)s,n,t,dev,(at+.5)/n,ur);}
  }
}
