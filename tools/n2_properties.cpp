// Pure N2 properties from Cantera (h2o2.yaml, mixture-averaged transport) on a temperature table, for
// the independent 1-D SST pipe reference (tools/sst_pipe_1d.py). The engine uses its own fits of the
// same Cantera data; the reference reads Cantera directly.
// Build: c++ -std=c++20 -O2 -I$HOME/src/cantera/include -I/opt/homebrew/include tools/n2_properties.cpp
//        $HOME/src/cantera/build/lib/libcantera.a -framework Accelerate -o /tmp/n2_properties
// Run:   /tmp/n2_properties > tools/n2_properties.csv
#include <cstdio>
#include "cantera/core.h"
#include "cantera/transport/Transport.h"
int main() {
  auto sol = Cantera::newSolution("h2o2.yaml", "", "mixture-averaged");
  auto gas = sol->thermo();
  auto tr = sol->transport();
  std::printf("# Cantera %s, h2o2.yaml, mixture-averaged, pure N2 at 101325 Pa (mu, lambda and cp do not depend on p here)\n",
              CANTERA_VERSION);
  std::printf("T_K,mu_Pa_s,lambda_W_m_K,cp_J_kg_K,R_J_kg_K\n");
  for (int i = 0; i <= 3300; ++i) {  // 200 to 3,500 K: N2 is inert in h2o2.yaml (no N atoms in the mechanism)
    double t = 200 + i;
    gas->setState_TPX(t, 101325, "N2:1");
    std::printf("%.1f,%.12e,%.12e,%.12e,%.12e\n", t, tr->viscosity(), tr->thermalConductivity(), gas->cp_mass(),
                Cantera::GasConstant / gas->meanMolecularWeight());
  }
}
