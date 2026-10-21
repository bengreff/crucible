// Throwaway check (not committed as a test): the C1 eqtt content at a short time on uniform 32 x 6 and 64 x 12,
// to compare with the AMR run's composite content.
#define main amr_budget_main
#include "tests/amr_budget.cpp"
#undef main

int main(int argc, char** argv) {
  const double end = argc > 1 ? std::atof(argv[1]) : 5e-5;
  const EquilibriumTable table = EquilibriumTable::load(std::getenv("CRUCIBLE_EQ_TABLE"));
  thermo::ReactionSource source("h2o2.yaml");
  for (int n : {1, 2}) {
    const Definition d = c1eqtt(source, 32 * n, 6 * n);
    Flow flow(d);
    thermo::ReactingFlow reacting(flow, "h2o2.yaml", 1, 1e-6, 1e-12, thermo::Chemistry::LocalEquilibrium);
    reacting.useEquilibriumTable(table, 0);
    const auto m0 = flow.measurements();
    while (flow.time() < end) reacting.step(end - flow.time());
    const auto m = flow.measurements();
    std::printf("%d x %d at %.4g ms: mass %.9e to %.9e kg, energy %.9e to %.9e J, wall heat flow %.4e W, inlet %.4e kg/s, steps %llu\n",
                32 * n, 6 * n, end * 1e3, m0.mass, m.mass, m0.energy, m.energy, m.wallHeatFlow, m.inletMassFlow,
                static_cast<unsigned long long>(m.steps));
  }
}
