// Scratch diagnostic for 3(b): pressure profiles of AMR, uniform N and uniform 2N at t (default 0.25).
// Built from the repository root after the library: c++ -std=c++20 -O2 -I. docs/evidence/amrb/sod_profile.cpp build/libcrucible_core.a -o sod_profile
// Run: ./sod_profile N [t] > sod_profile_N.txt
#include "core/amr.hpp"
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace crucible;
int main(int argc, char** argv) {
  const int n = argc > 1 ? std::atoi(argv[1]) : 64;
  const double t = argc > 2 ? std::atof(argv[2]) : 0.25;
  Definition d;
  d.experiment = Case::ShockTube;
  d.length = 1;
  d.nz = n;
  d.nr = 2;
  Definition fine = d;
  fine.nz = 2 * n;
  fine.nr = 4;
  const int first = static_cast<int>(std::lround(n / 3.0)), last = static_cast<int>(std::lround(2 * n / 3.0));
  Amr amr(d, 1);
  std::vector<unsigned char> a(amr.level(1).mesh().cells.size(), 0);
  for (int i = 0; i < 2 * n; ++i)
    for (int j = 0; j < 4; ++j) a[amr.level(1).mesh().index(i, j)] = (i >> 1) >= first && (i >> 1) < last;
  amr.setActive(1, a);
  amr.averageDown();
  Flow c(d, Mesh(amr.level(0).mesh())), f(fine, Mesh(amr.level(1).mesh()));
  auto content=[&](const Flow& g){ double m0=0,e0=0; for(std::size_t q=0;q<g.mesh().cells.size();++q){m0+=g.mesh().cells[q].volume*g.state()[q][0];e0+=g.mesh().cells[q].volume*g.state()[q][3];} return std::array<double,2>{m0,e0}; };
  const auto before=content(amr.level(0));
  amr.advanceTo(t);
  const auto after=content(amr.level(0));
  std::printf("# composite mass change %.3e relative, energy change %.3e relative (energy also moves by the end pressures' work: zero, u = 0 at both ends)\n",(after[0]-before[0])/before[0],(after[1]-before[1])/before[1]);
  const auto cb=content(c);
  c.advanceTo(t);
  const auto ca=content(c);
  std::printf("# uniform N: mass change %.3e, energy change %.3e; AMR steps %llu redone %llu\n",(ca[0]-cb[0])/cb[0],(ca[1]-cb[1])/cb[1],(unsigned long long)amr.steps(),(unsigned long long)amr.rejectedSteps());
  f.advanceTo(t);
  { const auto ma=amr.level(0).measurements(), mc=c.measurements(); std::printf("# balance errors: AMR base (composite) mass %.3e energy %.3e; uniform N mass %.3e energy %.3e\n",ma.massBalanceError,ma.energyBalanceError,mc.massBalanceError,mc.energyBalanceError); }
  const Mesh& m = c.mesh();
  const Mesh& fm = f.mesh();
  std::printf("# i x rho_amr rho_N rho_2N p_amr p_N p_2N (p_amr-p_N)/jump (p_2N-p_N)/jump\n");
  for (int i = 0; i < n; ++i) {
    const auto pa = amr.level(0).cellPrimitive(m.index(i, 0)), pc = c.cellPrimitive(m.index(i, 0));
    const auto f0 = f.cellPrimitive(fm.index(2 * i, 0)), f1 = f.cellPrimitive(fm.index(2 * i + 1, 0));
    const double p2 = 0.5 * (f0.p + f1.p), r2 = 0.5 * (f0.rho + f1.rho), jump = 0.203130178050647;
    std::printf("%3d %.4f %.5f %.5f %.5f %.6f %.6f %.6f %+.3e %+.3e\n", i, (i + 0.5) / n, pa.rho, pc.rho, r2, pa.p, pc.p, p2,
                (pa.p - pc.p) / jump, (p2 - pc.p) / jump);
  }
}
