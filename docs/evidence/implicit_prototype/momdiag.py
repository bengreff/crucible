#!/usr/bin/env python3
"""Print-only diagnostics for criterion 5's thrust offset (not for main). Applied to a source tree:
python3 momdiag.py <tree>. With CRUCIBLE_MOMDIAG=<t> set, each history sample at or after time t prints the booked
axial momentum rates (inlet, outlet, wall), the momentum budget's error and the gas's axial momentum; and each
explicit step ending at or after t + 0.016 ms
prints its two stages' rates and the booked mean; and each implicit step ending at or after t prints its booked rates (at the last inner iterate) beside the same rates
evaluated at the accepted state, and the accepted state's injector pressure beside the last iterate's. The march
is unchanged: the extra residual is evaluated into scratch arrays after the step is accepted, and every cache it
refreshes is recomputed by the next step's first residual."""
import sys
tree = sys.argv[1]

def edit(path, pairs):
    s = open(f"{tree}/{path}").read()
    for old, new in pairs:
        if old not in s: sys.exit(f"anchor not found in {path}: {old[:70]}")
        assert s.count(old) == 1, f"anchor not unique in {path}: {old[:70]}"
        s = s.replace(old, new)
    open(f"{tree}/{path}", "w").write(s)

edit("tests/chamber_study.cpp", [(
"    rows.push_back({m.time, m.outletMassFlow, m.injectorPressure, fVac, m.wallHeatFlow, m.wallAxialForce, cstar});\n",
"    rows.push_back({m.time, m.outletMassFlow, m.injectorPressure, fVac, m.wallHeatFlow, m.wallAxialForce, cstar});\n"
"    if (const char* from = std::getenv(\"CRUCIBLE_MOMDIAG\"); from && m.time >= std::atof(from))\n"
"      std::printf(\"MOMDIAG sample t %.9e inlet %.9e outlet %.9e wall %.9e body %.9e ambient %.9e F_vac %.9e exit-plane %.9e \"\n"
"                  \"momentum error %.3e axial momentum %.9e p_inj %.9e\\n\", m.time, m.inletMomentumFlux, m.outletMomentumFlux,\n"
"                  m.wallAxialForce, m.bodyAxialForce, m.ambientAxialForce, fVac, m.exitPlaneThrust, m.momentumBalanceError,\n"
"                  m.axialMomentum, m.injectorPressure);\n")])

edit("core/flow_implicit.cpp", [(
"                if(defect>control.tolerance) ++s.unconverged;\n                return dt;\n",
"                if(defect>control.tolerance) ++s.unconverged;\n"
"                if(const char* from=std::getenv(\"CRUCIBLE_MOMDIAG\"); from && time_>=std::atof(from)) {\n"
"                    // The last iterate's injector pressure (stage_ holds it, projected), then the accepted state's rates.\n"
"                    double pm=0,ps=0,area=0;\n"
"                    for(int j=0;j<m.nr;++j) {\n"
"                        const double ar=m.axialArea(0,j);area+=ar;\n"
"                        state_.swap(stage_);species_.swap(speciesStage_);pm+=ar*cellPrimitive(m.index(0,j)).p;\n"
"                        state_.swap(stage_);species_.swap(speciesStage_);ps+=ar*cellPrimitive(m.index(0,j)).p;\n"
"                    }\n"
"                    std::vector<Conserved> d(n);std::vector<double> sd(n*ns_),td;\n"
"                    const auto at=rhs(state_,species_,turbulence_,d,sd,td,time_);\n"
"                    double net=0;\n"
"                    for(std::size_t q=0;q<n;++q) net+=d[q][1]*m.cells[q].volume;\n"
"                    std::printf(\"MOMDIAG implicit t %.9e iterations %d defect %.3e booked: inlet %.9e outlet %.9e wall %.9e net %.6e; \"\n"
"                                \"at the accepted state: inlet %.9e outlet %.9e wall %.9e net %.6e (sum of cells %.6e); \"\n"
"                                \"p_inj iterate %.9e accepted %.9e\\n\",time_,it,defect,rates.inletMomentum,rates.outletMomentum,\n"
"                                rates.wallAxial,rates.netMomentum(),at.inletMomentum,at.outletMomentum,at.wallAxial,at.netMomentum(),\n"
"                                net,pm/area,ps/area);\n"
"                }\n"
"                return dt;\n")])
edit("core/flow.cpp", [(
"    lastRates_=BoundaryRates::average(attempt.first,attempt.second);\n",
"    lastRates_=BoundaryRates::average(attempt.first,attempt.second);\n"
"    if(const char* from=std::getenv(\"CRUCIBLE_MOMDIAG\"); from && time_>=std::atof(from)+0.016e-3) {\n"
"        // The explicit step's two stages: the rates at U and at U* = U + dt R(U), and their mean (booked).\n"
"        const auto& a=attempt.first;const auto& b=attempt.second;\n"
"        std::printf(\"MOMDIAG explicit t %.12e dt %.6e first: inlet %.9e outlet %.9e wall %.9e net %.6e; second: inlet %.9e \"\n"
"                    \"outlet %.9e wall %.9e net %.6e; booked net %.6e\\n\",time_,dt,a.inletMomentum,a.outletMomentum,a.wallAxial,\n"
"                    a.netMomentum(),b.inletMomentum,b.outletMomentum,b.wallAxial,b.netMomentum(),lastRates_.netMomentum());\n"
"    }\n")])
s = open(f"{tree}/core/flow.cpp").read()
open(f"{tree}/core/flow.cpp", "w").write("#include <cstdio>\n#include <cstdlib>\n" + s)
s = open(f"{tree}/core/flow_implicit.cpp").read()
open(f"{tree}/core/flow_implicit.cpp", "w").write("#include <cstdlib>\n" + s)
print("momdiag applied to", tree)
