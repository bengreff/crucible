// HP equilibrium of gaseous H2 and O2 on a Cantera mechanism, for the check against CEA's chamber
// problem (docs/evidence/TABLE_B.md, *The 3,500 K ceiling*; driven by tools/cea_equilibrium_check.py).
// Reads lines "T_reactants_K O/F p_Pa" on stdin; both reactants enter at T_reactants, and the
// mixture's enthalpy is their mass-weighted sum, as in CEA. Writes, per line, the adiabatic
// equilibrium T and the mole fractions of H2, H2O, OH, H, O and O2.
// Usage: crucible_equilibrium_hp <mechanism.yaml> < cases
#include <cstdio>
#include <iostream>

#include "cantera/core.h"

int main(int argc, char** argv) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s <mechanism.yaml> < cases\n", argv[0]);
    return 2;
  }
  auto sol = Cantera::newSolution(argv[1], "", "none");
  auto gas = sol->thermo();
  const char* out[] = {"H2", "H2O", "OH", "H", "O", "O2"};
  std::printf("# %s: T0_K of p_Pa T_K", argv[1]);
  for (const char* s : out) std::printf(" X_%s", s);
  std::printf("\n");
  double T0, of, p;
  while (std::cin >> T0 >> of >> p) {
    gas->setState_TPX(T0, p, "H2:1");
    const double hf = gas->enthalpy_mass();
    gas->setState_TPX(T0, p, "O2:1");
    const double ho = gas->enthalpy_mass();
    const double yo = of / (1.0 + of);
    gas->setState_TPY(T0, p, Cantera::Composition{{"H2", 1.0 - yo}, {"O2", yo}});
    gas->setState_HP(yo * ho + (1.0 - yo) * hf, p);
    gas->equilibrate("HP");
    std::printf("%.2f %.4f %.6e %.4f", T0, of, p, gas->temperature());
    for (const char* s : out) std::printf(" %.6e", gas->moleFraction(s));
    std::printf("\n");
  }
  return 0;
}
