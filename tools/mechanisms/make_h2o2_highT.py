"""Writes tools/mechanisms/h2o2_highT.yaml: Cantera's h2o2.yaml (reactions, transport, the
ideal-gas phase) with every species' thermo replaced by the NASA 7-coefficient fits to 6,000 K of
McBride, Gordon and Reno, NASA TM-4513 (1993), as shipped in Cantera's nasa_gas.yaml. The engine
reads NASA-7 only (adapters/reaction.cpp), so these fits enter with no change to the core.
Ar's single 200 to 6,000 K range is written as two identical ranges split at 1,000 K (the same
function). TM-4513's standard state is 1 bar, which nasa_gas.yaml does not declare (Cantera then
reads the fits at 1 atm, shifting every species' s by R ln(1.01325) and each Kp by (1.01325)^dn), so
each species here declares reference-pressure 1 bar. OH's heat of formation is then moved to Ruscic
et al. (2002)'s, the value CEA's thermo.inp carries (the Director's decision, 20 October 2026;
TM-4513's predates the revision): a6 of both ranges is shifted so h(298.15 K) = HF298_OH, which
leaves cp and s unchanged.
Checked against CEA's own NASA-9 data by check_h2o2_highT.py (docs/evidence/TABLE_B.md).
Usage: python3 tools/mechanisms/make_h2o2_highT.py [cantera data dir]"""
import os, sys, yaml

R = 8.31446261815324  # J/(mol K)
# Ruscic, Wagner, Harding et al., J. Phys. Chem. A 106 (2002) 2727, "On the enthalpy of formation of
# hydroxyl radical and gas-phase bond dissociation energies of water and hydroxyl"; the value in CEA's
# thermo.inp (McBride, Zehe and Gordon, NASA TP-2002-211556; RocketCEA 1.2.3), OH entry, J/mol.
HF298_OH = 37278.206

def h298(a):  # NASA-7 enthalpy at 298.15 K, J/mol
    T = 298.15
    return R * T * (a[0] + a[1] * T / 2 + a[2] * T**2 / 3 + a[3] * T**3 / 4 + a[4] * T**4 / 5 + a[5] / T)

data = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/src/cantera/data")
mech = yaml.safe_load(open(os.path.join(data, "h2o2.yaml")))
nasa = {str(s["name"]).upper(): s for s in yaml.safe_load(open(os.path.join(data, "nasa_gas.yaml")))["species"]}

for sp in mech["species"]:
    src = nasa[sp["name"].upper()]["thermo"]
    ranges, coeffs = list(src["temperature-ranges"]), [list(c) for c in src["data"]]
    if len(coeffs) == 1:
        ranges, coeffs = [ranges[0], 1000.0, ranges[1]], [coeffs[0], list(coeffs[0])]
    sp["thermo"] = {"model": "NASA7", "temperature-ranges": ranges, "data": coeffs,
                    "reference-pressure": 1.0e5,
                    "note": f"NASA TM-4513 (McBride, Gordon and Reno 1993) via Cantera nasa_gas.yaml; source {src.get('note', '')}"}
    if sp["name"].upper() == "OH":
        shift = (HF298_OH - h298(coeffs[0])) / R
        for c in coeffs:
            c[5] += shift
        sp["thermo"]["note"] += (f"; Hf298 moved from TM-4513's {HF298_OH - R * shift:.3f} to {HF298_OH:.3f} J/mol, "
                                 "Ruscic et al. 2002 as in CEA thermo.inp (a6 of both ranges shifted)")
mech["phases"] = [p for p in mech["phases"] if p["name"] == "ohmech"]
mech["description"] = ("h2o2.yaml (Cantera) with NASA TM-4513 NASA-7 thermo to 6,000 K for every species (standard state "
                       "1 bar, as TM-4513), OH's heat of "
                       "formation from Ruscic et al. 2002 (as in CEA); "
                       "written by tools/mechanisms/make_h2o2_highT.py.\n" + mech.get("description", ""))
for k in ("generator", "cantera-version", "date", "input-files"):
    mech.pop(k, None)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "h2o2_highT.yaml")
with open(out, "w") as f:
    yaml.safe_dump(mech, f, sort_keys=False, width=110)
print("wrote", out)
