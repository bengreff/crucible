"""Diagnostic (TABLE_B.md, *The equilibrium check*): writes h2o2_highT.yaml with CEA's own NASA-9 thermo
(RocketCEA 1.2.3's thermo.inp, ranges to 6,000 K) for the listed species, to separate data from
everything else in the check against CEA. Usage (repo root): python3 diag_nasa9_swap.py ALL|H2O,OH,... out.yaml
Reference pressure is whatever h2o2_highT.yaml declares."""
import sys, yaml
inp = open("/tmp/ceaenv_new/lib/python3.14/site-packages/rocketcea/thermo.inp").read().splitlines()
def nasa9(name):
    for i, line in enumerate(inp):
        if line[:18].split() and line[:18].split()[0].upper() == name and inp[i + 1][:2].strip().isdigit():
            n, ivals, j = int(inp[i + 1][:2]), [], i + 2
            f = lambda s: float(s.replace("D", "E"))
            for _ in range(n):
                lo, hi = float(inp[j][:11]), float(inp[j][11:22])
                a = [f(inp[j + 1][16 * k:16 * k + 16]) for k in range(5)] + [f(inp[j + 2][16 * k:16 * k + 16]) for k in range(2)]
                b = [f(inp[j + 2][48:64]), f(inp[j + 2][64:80])]
                ivals.append((lo, hi, a, b)); j += 3
            return ivals
which, out = sys.argv[1].split(","), sys.argv[2]
m = yaml.safe_load(open("tools/mechanisms/h2o2_highT.yaml"))
for sp in m["species"]:
    n = sp["name"].upper()
    if which == ["ALL"] or n in which:
        iv = [v for v in nasa9(n) if v[0] < 5999.0]
        sp["thermo"] = {"model": "NASA9", "temperature-ranges": [iv[0][0]] + [v[1] for v in iv],
                        "data": [a + b for (lo, hi, a, b) in iv]}
yaml.safe_dump(m, open(out, "w"), sort_keys=False)
