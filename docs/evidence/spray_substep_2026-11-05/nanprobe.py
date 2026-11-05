import sys
p = sys.argv[1] + "/core/spray.cpp"
s = open(p).read()
old = '        if (std::isnan(r)) throw std::runtime_error("Spray: the heat-transfer number\'s residual is NaN.");\n'
assert s.count(old) == 1, "anchor"
new = r'''        if (std::isnan(r)) {  // print-only probe (5 November): the residual's terms at the throw
          std::printf("NANPROBE b %.6e phi %.6e lnBm %.6e bm %.6e cpf %.6e rhoD %.6e Sh %.6e lambda %.6e nu0 %.6e Fb %.6e; "
                      "drop d %.6e T %.6f z %.6e r %.6e count %.6e boiling %d; Tsat %.6f ts %.6f tr %.6f ySurface %.9e yInf %.9e "
                      "wOther %.6e Re %.6e Pr %.6e Sc %.6e; gas T %.6f p %.6f rho %.6e u %.4e %.4e %.4e; y",
                      b, phi, lnBm, bm, cpf, rhoD, out.sherwood, lambda, nu0, filmCorrection(b), parcel.d, parcel.t, parcel.z,
                      parcel.r, parcel.count, static_cast<int>(out.boiling), tSat, ts, tr, ySurface, yInf, wOther, out.reynolds, pr,
                      sc, gas.t, gas.p, gas.rho, gas.uz, gas.ur, gas.ut);
          for (std::size_t k = 0; k < ns; ++k) std::printf(" %.3e", gas.y[k]);
          std::printf("\n");
          std::fflush(stdout);
          throw std::runtime_error("Spray: the heat-transfer number's residual is NaN.");
        }
'''
s = s.replace(old, new)
open(p, "w").write(s)
print("nanprobe applied")
s = open(p).read()
old2 = "    // Velocity: exact relaxation toward the seen gas over the step at stage 1's drag; position exact.\n"
assert s.count(old2) == 1, "anchor 2"
new2 = r'''    {  // print-only probe (5 November): the stored temperature against the integrated enthalpy
      static const bool hprobe = std::getenv("LBPROBE_ON") != nullptr;
      if (hprobe && m1 >= minMass && std::abs(liquid_.enthalpy(t1) - h1) > 1e-8 * std::abs(h1))
        std::printf("HPROBE h1 %.9e h(t1) %.9e t1 %.6f p.t %.6f table %.3f to %.3f (h %.6e to %.6e) held %d boiling %d clampedS %d "
                    "d %.6e m1 %.6e count %.6e step %.4e; gas T %.4f p %.2f; Tsat %.4f heat %.6e mass %.6e\n",
                    h1, liquid_.enthalpy(t1), t1, p.t, liquid_.minTemperature(), liquid_.maxTemperature(),
                    liquid_.enthalpy(liquid_.minTemperature()), liquid_.enthalpy(liquid_.maxTemperature()),
                    static_cast<int>(r1.boiling && r1.heat > 0), static_cast<int>(r1.boiling), static_cast<int>(r1.clamped), p.d, m1,
                    p.count, step, gas.t, gas.p, tSat, r1.heat, r1.mass);
    }
''' + old2
s = s.replace(old2, new2)
open(p, "w").write(s)
print("hprobe applied")
