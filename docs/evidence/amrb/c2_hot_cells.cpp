// Throwaway diagnostic (not committed as a test): in the criterion 2 AMR run, which cells pass 3500 K (Table A's
// top), by level and kind: owned (active, not covered), covered, or filled (not active).
#define main amr_budget_main
#include "tests/amr_budget.cpp"
#undef main

int main(int argc, char** argv) {
  const double end = argc > 1 ? std::atof(argv[1]) : 5e-4;
  const EquilibriumTable table = EquilibriumTable::load(std::getenv("CRUCIBLE_EQ_TABLE"));
  thermo::ReactionSource source("h2o2.yaml");
  const Definition d = c1eqtt(source, 32, 6);
  Amr amr(d, 1);
  std::vector<std::unique_ptr<thermo::ReactingFlow>> reacting;
  for (int l = 0; l < amr.levels(); ++l) {
    reacting.push_back(std::make_unique<thermo::ReactingFlow>(amr.level(l), "h2o2.yaml", 1, 1e-6, 1e-12,
                                                              thermo::Chemistry::LocalEquilibrium));
    reacting.back()->useEquilibriumTable(table, 0);
  }
  amr.setAfterStep([&](int l, double dt) { reacting[static_cast<std::size_t>(l)]->react(dt); });
  const Mesh& fine = amr.level(1).mesh();
  std::vector<unsigned char> active(fine.cells.size(), 0);
  for (int i = 0; i < fine.nz; ++i)
    for (int j = 0; j < fine.nr; ++j) active[fine.index(i, j)] = ((i >> 1) < 16 && (j >> 1) <= 1) || (j >> 1) == 5;
  amr.setActive(1, active);
  amr.averageDown();
  const char* kinds[3] = {"owned", "covered", "filled"};
  double worst[2][3] = {};
  long hot[2][3] = {};
  int reported = 0;
  while (amr.time() < end) {
    amr.step(end - amr.time());
    for (int l = 0; l < amr.levels(); ++l) {
      const Flow& f = amr.level(l);
      const Mesh& m = f.mesh();
      const auto& act = amr.active(l);
      for (int i = 0; i < m.nz; ++i)
        for (int j = 0; j < m.nr; ++j) {
          const auto q = m.index(i, j);
          const bool isActive = act.empty() || act[q];
          const int kind = !isActive ? 2 : amr.covered(l, q) ? 1 : 0;
          const double t = f.temperature(q);
          worst[l][kind] = std::max(worst[l][kind], t);
          if (t > 3500) {
            ++hot[l][kind];
            if (reported < 12) {
              ++reported;
              std::printf("t %.4f ms: level %d %s cell (%d, %d) at z %.4f r %.4f: T %.1f K\n", amr.time() * 1e3, l, kinds[kind], i, j,
                          m.cells[q].z, m.cells[q].r, t);
            }
          }
        }
    }
  }
  for (int l = 0; l < 2; ++l)
    for (int k = 0; k < 3; ++k)
      std::printf("level %d %-8s largest T %.1f K, cell-steps above 3500 K %ld\n", l, kinds[k], worst[l][k], hot[l][k]);
}
