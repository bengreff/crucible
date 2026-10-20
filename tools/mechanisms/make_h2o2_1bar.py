"""Writes tools/mechanisms/h2o2_1bar.yaml: Cantera's h2o2.yaml (GRI-Mech 3.0 H/O submechanism),
byte-for-byte identical except that every species' thermo block gains reference-pressure: 1.0e5
(1 bar), the same declaration make_h2o2_highT.py adds for the TM-4513 fits. check_h2o2_refp.py found
every reacting species' (H2, H, O, O2, OH, H2O, HO2, H2O2) raw NASA-7 s and g already consistent with
CEA's 1-bar NASA-9 data (|d(s/R)| a few e-3 at most), so without this declaration Cantera reads them
at its hardcoded default of 1 atm, a reference-pressure mismatch, not a refit: no species' NASA-7
coefficients (cp, h) are touched, only the reference-pressure tag is added. AR and N2 (inert in C1,
not checked for this reason in the gate) measured closer to the 1-atm signature or in between; they
get the same declaration for uniformity, as make_h2o2_highT.py did.

Usage: python3 tools/mechanisms/make_h2o2_1bar.py [cantera data dir]
"""
import os, sys, yaml

data = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/src/cantera/data")
mech = yaml.safe_load(open(os.path.join(data, "h2o2.yaml")))

for sp in mech["species"]:
    sp["thermo"]["reference-pressure"] = 1.0e5

mech["description"] = ("h2o2.yaml (Cantera, GRI-Mech 3.0 H/O submechanism) with reference-pressure 1 bar "
                        "declared on every species' thermo (no other change: same NASA-7 coefficients, same "
                        "temperature ranges, same reactions); written by tools/mechanisms/make_h2o2_1bar.py.\n"
                        + mech.get("description", ""))

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "h2o2_1bar.yaml")
with open(out, "w") as f:
    yaml.safe_dump(mech, f, sort_keys=False, width=110)
print("wrote", out)
