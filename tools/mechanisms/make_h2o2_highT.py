"""Writes tools/mechanisms/h2o2_highT.yaml: Cantera's h2o2.yaml (reactions, transport, the
ideal-gas phase) with every species' thermo replaced by the NASA 7-coefficient fits to 6,000 K of
McBride, Gordon and Reno, NASA TM-4513 (1993), as shipped in Cantera's nasa_gas.yaml. The engine
reads NASA-7 only (adapters/reaction.cpp), so these fits enter with no change to the core.
Ar's single 200 to 6,000 K range is written as two identical ranges split at 1,000 K (the same
function). Checked against CEA's own NASA-9 data by check_h2o2_highT.py (docs/evidence/TABLE_B.md).
Usage: python3 tools/mechanisms/make_h2o2_highT.py [cantera data dir]"""
import os, sys, yaml

data = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/src/cantera/data")
mech = yaml.safe_load(open(os.path.join(data, "h2o2.yaml")))
nasa = {str(s["name"]).upper(): s for s in yaml.safe_load(open(os.path.join(data, "nasa_gas.yaml")))["species"]}

for sp in mech["species"]:
    src = nasa[sp["name"].upper()]["thermo"]
    ranges, coeffs = list(src["temperature-ranges"]), [list(c) for c in src["data"]]
    if len(coeffs) == 1:
        ranges, coeffs = [ranges[0], 1000.0, ranges[1]], [coeffs[0], list(coeffs[0])]
    sp["thermo"] = {"model": "NASA7", "temperature-ranges": ranges, "data": coeffs,
                    "note": f"NASA TM-4513 (McBride, Gordon and Reno 1993) via Cantera nasa_gas.yaml; source {src.get('note', '')}"}
mech["phases"] = [p for p in mech["phases"] if p["name"] == "ohmech"]
mech["description"] = ("h2o2.yaml (Cantera) with NASA TM-4513 NASA-7 thermo to 6,000 K for every species; "
                       "written by tools/mechanisms/make_h2o2_highT.py.\n" + mech.get("description", ""))
for k in ("generator", "cantera-version", "date", "input-files"):
    mech.pop(k, None)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "h2o2_highT.yaml")
with open(out, "w") as f:
    yaml.safe_dump(mech, f, sort_keys=False, width=110)
print("wrote", out)
