import sys
p = sys.argv[1] + "/tests/rl10_tier1.cpp"
s = open(p).read()
anchor = "                                       : amr.step(std::min(other, bound));\n"
assert s.count(anchor) == 1, "anchor"
probe = anchor + r'''    if (std::getenv("LBPROBE_ON")) {  // print-only probe (5 November): the liquid budget's three terms, per step
      static std::array<double, 3> last{};
      static Spray::Ledger before{};
      static bool first = true;
      const auto e = spray.balanceErrors();
      const auto& L = spray.ledger();
      bool moved = first;
      for (int c = 0; c < 3; ++c) moved |= std::abs(e[c] - last[c]) > 1e-13;
      if (moved)
        std::printf("LBPROBE t %.9e dt %.4e mass %.4e energy %.4e momentum %.4e; liquid %.9e kg, %zu parcels; injected +%.4e "
                    "exchanged +%.4e out +%.4e; vaporised +%ld walls +%ld children +%ld flashes +%ld clamped +%ld substeps +%ld\n",
                    amr.time(), dt, e[0], e[1], e[2], spray.liquidMass(), spray.parcels().size(), L.injectedMass - before.injectedMass,
                    L.exchangedMass - before.exchangedMass, L.outflowMass - before.outflowMass, L.vaporised - before.vaporised,
                    L.wallHits - before.wallHits, L.children - before.children, L.flashes - before.flashes,
                    L.clampedCondensation - before.clampedCondensation, L.substeps - before.substeps);
      last = e;
      before = L;
      first = false;
    }
'''
s = s.replace(anchor, probe)
open(p, "w").write(s)
print("lbprobe applied")
