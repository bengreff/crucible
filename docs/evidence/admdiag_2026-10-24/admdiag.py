#!/usr/bin/env python3
"""Print-only diagnostics for "AMR: could not advance an admissible state" (not for main). Applied to a source tree:
python3 admdiag.py <tree>. With CRUCIBLE_ADMDIAG set, each rejected base step prints which check failed (the level,
the exit in Amr::advance, the stage) and the first failing cells' state (z, r, rho, internal energy against its floor,
rho k, rho omega, the smallest partial density), capped at 600 lines. On a tree with the reflux clip, each clip also
prints its level, cell and rho k. Nothing is changed but the printing."""
import sys
tree = sys.argv[1]

def edit(path, pairs):
    s = open(f"{tree}/{path}").read()
    for old, new, required in pairs:
        if old not in s:
            if required: sys.exit(f"anchor not found in {path}: {old[:70]}")
            print(f"{path}: optional anchor absent ({old[:50]})"); continue
        assert s.count(old) == 1, f"anchor not unique in {path}: {old[:70]}"
        s = s.replace(old, new)
    open(f"{tree}/{path}", "w").write(s)

edit("core/flow.hpp", [(
"    bool admissible(const Conserved& u, const double* partial, const double* turbulence) const;\n",
"    bool admissible(const Conserved& u, const double* partial, const double* turbulence) const;\n"
"    bool admdiag(const char* where, std::size_t q, const Conserved& u, const double* partial, const double* turbulence) const;\n", True)])

flow_fn = r'''
bool admdiagOn() { static const bool on=std::getenv("CRUCIBLE_ADMDIAG")!=nullptr; return on; }
std::atomic<int> admdiagLines{0};
}  // namespace (admdiag)
bool Flow::admdiag(const char* where,std::size_t q,const Conserved& u,const double* partial,const double* turbulence) const {
    if(admissible(u,partial,turbulence)) return true;
    if(admdiagOn() && admdiagLines.fetch_add(1)<600) {
        double sum=0,least=1e300;
        for(std::size_t k=0;k<ns_;++k) { sum+=std::max(partial[k],0.0);least=std::min(least,partial[k]); }
        const double k=nt_ && u[0]>0?turbulence[0]/u[0]:0;
        const double e=(u[3]-(u[1]*u[1]+u[2]*u[2])/(2*u[0]))/u[0]-k;
        const double floor=sum>0?medium_.energyFloorOfClipped(partial)/sum:0;
        std::printf("ADMDIAG %s t %.9e cell %zu z %.6e r %.6e rho %.6e e %.6e floor %.6e rhok %.6e rhoomega %.6e partial min %.3e sum %.6e\n",
                    where,time_,q,mesh_.cells[q].z,mesh_.cells[q].r,u[0],e,floor,nt_?turbulence[0]:0.0,nt_?turbulence[1]:0.0,least,sum);
    }
    return false;
}
namespace {  // (admdiag)
'''
edit("core/flow.cpp", [
("bool admissibleBulk(const Conserved& u, double floor) {", flow_fn + "bool admissibleBulk(const Conserved& u, double floor) {", True),
("        if(!allOf(region(0),[&](std::size_t q) { return !isActive(q) || admissible(state_[q],species_.data()+q*ns_,turbulence_.data()+q*nt_); }))\n            return false;",
 "        if(!allOf(region(0),[&](std::size_t q) { return !isActive(q) || admdiag(\"first stage, start turbulence source\",q,state_[q],species_.data()+q*ns_,turbulence_.data()+q*nt_); }))\n            return false;", True),
("        clipTurbulence(turbulenceStage_,q);\n        return admissible(stage_[q],speciesStage_.data()+q*ns_,turbulenceStage_.data()+q*nt_);",
 "        clipTurbulence(turbulenceStage_,q);\n        return admdiag(\"first stage\",q,stage_[q],speciesStage_.data()+q*ns_,turbulenceStage_.data()+q*nt_);", True),
("        clipTurbulence(turbulenceNext_,q);\n        return admissible(next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_);",
 "        clipTurbulence(turbulenceNext_,q);\n        return admdiag(\"second stage\",q,next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_);", True),
("        return allOf(region(0),[&](std::size_t q) { return !isActive(q) || admissible(next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_); });",
 "        return allOf(region(0),[&](std::size_t q) { return !isActive(q) || admdiag(\"second stage, end turbulence source\",q,next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_); });", True),
("    const bool admitted=prescribeWallTurbulence(state,species,turbulence);\n",
 "    const bool admitted=prescribeWallTurbulence(state,species,turbulence);\n    if(!admitted && std::getenv(\"CRUCIBLE_ADMDIAG\")) std::printf(\"ADMDIAG prescribe (wall functions) not admitted, t %.9e\\n\",time_);\n", True),
])
# flow.cpp needs <atomic>, <cstdio>, <cstdlib>, and the helper must sit in an anonymous namespace before admissibleBulk.
s = open(f"{tree}/core/flow.cpp").read()
s = "#include <atomic>\n#include <cstdio>\n#include <cstdlib>\n" + s
open(f"{tree}/core/flow.cpp", "w").write(s)

edit("core/amr.cpp", [
("    if(!f.firstStage(dt,attempt)) return false;\n    if(l>0) fill(l,f.stage_,f.speciesStage_,f.turbulenceStage_,0.5*(half+1));\n    if(!f.secondStage(dt,attempt)) return false;",
 "    const bool diag=std::getenv(\"CRUCIBLE_ADMDIAG\")!=nullptr;\n"
 "    if(!f.firstStage(dt,attempt)) { if(diag) std::printf(\"ADMDIAG exit: level %d first stage, t %.9e dt %.6e half %d\\n\",l,t,dt,half); return false; }\n"
 "    if(l>0) fill(l,f.stage_,f.speciesStage_,f.turbulenceStage_,0.5*(half+1));\n"
 "    if(!f.secondStage(dt,attempt)) { if(diag) std::printf(\"ADMDIAG exit: level %d second stage, t %.9e dt %.6e half %d\\n\",l,t,dt,half); return false; }", True),
("        averageDown(l);\n        if(!reflux(l)) return false;",
 "        averageDown(l);\n        if(!reflux(l)) { if(diag) std::printf(\"ADMDIAG exit: level %d reflux from level %d, t %.9e dt %.6e\\n\",l,l+1,t,dt); return false; }", True),
("        if(!c.admissible(c.state_[u],c.species_.data()+u*ns,nt?c.turbulence_.data()+u*nt:nullptr)) return false;",
 "        if(!c.admdiag(\"reflux\",u,c.state_[u],c.species_.data()+u*ns,nt?c.turbulence_.data()+u*nt:nullptr)) return false;", True),
("        const double added=-c.turbulence_[u*nt]*cm.cells[u].volume;\n",
 "        const double added=-c.turbulence_[u*nt]*cm.cells[u].volume;\n"
 "        if(std::getenv(\"CRUCIBLE_ADMDIAG\")) std::printf(\"ADMDIAG reflux clip t %.9e level %d cell %zu z %.6e r %.6e rhok %.6e rho %.6e\\n\",c.time_,l,u,cm.cells[u].z,cm.cells[u].r,c.turbulence_[u*nt],c.state_[u][0]);\n", False),
("        restore();++rejected_;dt*=0.5;landing=false;",
 "        if(std::getenv(\"CRUCIBLE_ADMDIAG\")) std::printf(\"ADMDIAG rejected base step: t %.9e dt %.6e attempt %d\\n\",time_,dt,attempt);\n        restore();++rejected_;dt*=0.5;landing=false;", True),
])
s = open(f"{tree}/core/amr.cpp").read()
s = "#include <cstdio>\n#include <cstdlib>\n" + s
open(f"{tree}/core/amr.cpp", "w").write(s)
print("admdiag applied to", tree)
