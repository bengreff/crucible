"""Species thermo of h2o2_highT.yaml (NASA TM-4513 NASA-7) and Cantera's h2o2.yaml (GRI NASA-7,
to 3,500 K) against CEA's own NASA-9 data (thermo.inp as shipped with RocketCEA; McBride, Zehe and
Gordon, NASA TP-2002-211556). Prints, per species and temperature band, the largest |cp| error
relative to CEA, |h| error (kJ/mol) and |g/RT| error (the quantity equilibrium sees).
Usage: python3 check_h2o2_highT.py <thermo.inp> [cantera data dir]"""
import math, os, sys, yaml

R = 8.31446261815324
here = os.path.dirname(os.path.abspath(__file__))
inp = open(sys.argv[1]).read().splitlines()
data = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser("~/src/cantera/data")
names = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"]

def nasa9(name):
    for i, line in enumerate(inp):
        if line[:18].split() and line[:18].split()[0].upper() == name and inp[i + 1][:2].strip().isdigit():
            n, ivals, j = int(inp[i + 1][:2]), [], i + 2
            for _ in range(n):
                lo, hi = float(inp[j][:11]), float(inp[j][11:22])
                f = lambda s: float(s.replace("D", "E"))
                a = [f(inp[j + 1][16 * k:16 * k + 16]) for k in range(5)] + [f(inp[j + 2][16 * k:16 * k + 16]) for k in range(2)]
                b = [f(inp[j + 2][48:64]), f(inp[j + 2][64:80])]
                ivals.append((lo, hi, a, b)); j += 3
            return ivals
    raise KeyError(name)

def eval9(iv, T):
    lo, hi, a, b = next(v for v in iv if v[0] <= T <= v[1])
    cp = a[0] / T**2 + a[1] / T + a[2] + a[3] * T + a[4] * T**2 + a[5] * T**3 + a[6] * T**4
    h = -a[0] / T**2 + a[1] * math.log(T) / T + a[2] + a[3] * T / 2 + a[4] * T**2 / 3 + a[5] * T**3 / 4 + a[6] * T**4 / 5 + b[0] / T
    s = -a[0] / T**2 / 2 - a[1] / T + a[2] * math.log(T) + a[3] * T + a[4] * T**2 / 2 + a[5] * T**3 / 3 + a[6] * T**4 / 4 + b[1]
    return cp, h, s  # cp/R, h/RT, s/R

def nasa7(file, name):
    for s in yaml.safe_load(open(file))["species"]:
        if str(s["name"]).upper() == name:
            return s["thermo"]["temperature-ranges"], s["thermo"]["data"]

def eval7(rd, T):
    r, d = rd
    a = d[0] if T <= r[1] else d[-1]
    cp = a[0] + a[1] * T + a[2] * T**2 + a[3] * T**3 + a[4] * T**4
    h = a[0] + a[1] * T / 2 + a[2] * T**2 / 3 + a[3] * T**3 / 4 + a[4] * T**4 / 5 + a[5] / T
    s = a[0] * math.log(T) + a[1] * T + a[2] * T**2 / 2 + a[3] * T**3 / 3 + a[4] * T**4 / 4 + a[6]
    return cp, h, s

bands = [(300, 3000), (3000, 3500), (3500, 5000), (5000, 6000)]
for label, file in [("h2o2_highT.yaml (TM-4513)", os.path.join(here, "h2o2_highT.yaml")), ("h2o2.yaml (GRI, extrapolated above 3,500 K)", os.path.join(data, "h2o2.yaml"))]:
    print(f"\n{label} against CEA thermo.inp: max |dcp|/cp, max |dh| kJ/mol, max |d(g/RT)| per band")
    print("species " + "  ".join(f"{lo}-{hi} K".center(30) for lo, hi in bands))
    worst = [[0, 0, 0] for _ in bands]
    for n in names:
        iv, rd = nasa9(n), nasa7(file, n)
        row = []
        for bi, (lo, hi) in enumerate(bands):
            dcp = dh = dg = 0.0
            for T in range(lo, hi + 1, 10):
                c9, h9, s9 = eval9(iv, T); c7, h7, s7 = eval7(rd, T)
                dcp = max(dcp, abs(c7 - c9) / c9); dh = max(dh, abs(h7 - h9) * R * T / 1000)
                dg = max(dg, abs((h7 - s7) - (h9 - s9)))
            row.append(f"{dcp:9.2e} {dh:8.3f} {dg:9.2e}")
            worst[bi] = [max(worst[bi][0], dcp), max(worst[bi][1], dh), max(worst[bi][2], dg)]
        print(f"{n:7s} " + "  ".join(row))
    print("worst   " + "  ".join(f"{w[0]:9.2e} {w[1]:8.3f} {w[2]:9.2e}" for w in worst))
