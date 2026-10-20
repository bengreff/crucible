"""Per-species standard state of Cantera's h2o2.yaml (GRI-Mech 3.0 H/O submechanism, NASA-7, no
reference-pressure declared) against CEA's own NASA-9 data (thermo.inp as shipped with RocketCEA;
McBride, Zehe and Gordon, NASA TP-2002-211556), read the same way as check_h2o2_highT.py (raw
polynomial evaluation; no Cantera python import).

Prints, per species and temperature (298.15, 1000, 3000 K), d(s/R) = s/R(h2o2.yaml) - s/R(CEA) and
d(g/RT) = g/RT(h2o2.yaml) - g/RT(CEA), g/RT = h/RT - s/R. A species whose d(s/R) is about
-ln(1.01325) = -0.013163 is consistent with h2o2.yaml's raw NASA-7 s being read as if it were at a
1-bar reference while actually a 1-atm fit would read LOWER than a 1-bar fit by that amount (S at the
lower reference pressure, 1 bar, is the larger one); a species whose d(s/R) is within about 1e-3 of
zero is consistent with h2o2.yaml's fit already matching CEA's 1-bar-based value. The sign actually
measured is printed and not forced either way. OH's h2o2.yaml enthalpy is known to differ from CEA's
(pre-Ruscic et al. 2002) and will move d(g/RT) without moving d(s/R).

Usage: python3 check_h2o2_refp.py <thermo.inp> [cantera data dir]
"""
import math, os, sys, yaml

R = 8.31446261815324
here = os.path.dirname(os.path.abspath(__file__))
inp = open(sys.argv[1]).read().splitlines()
data = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser("~/src/cantera/data")
names = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"]
TS = [298.15, 1000.0, 3000.0]
LN_1ATM_PER_BAR = math.log(1.01325)


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


mech_file = os.path.join(data, "h2o2.yaml")
print(f"# h2o2.yaml ({mech_file}) NASA-7 (no reference-pressure declared) raw s/R, g/RT against CEA thermo.inp NASA-9")
print(f"# d(s/R) = s/R(h2o2.yaml) - s/R(CEA); d(g/RT) = g/RT(h2o2.yaml) - g/RT(CEA), g/RT = h/RT - s/R")
print(f"# -ln(1.01325) = {-LN_1ATM_PER_BAR:.6f} (h2o2.yaml 1-atm-fit signature, see docstring); ~0 (1-bar-fit signature)")
header = "species " + "  ".join(f"T={t:<9.2f}".rjust(24) for t in TS)
print(header)
print("        " + "  ".join(f"{'d(s/R)':>10} {'d(g/RT)':>10}" for _ in TS))
rows = {}
for n in names:
    iv, rd = nasa9(n), nasa7(mech_file, n)
    row = []
    for T in TS:
        cp9, h9, s9 = eval9(iv, T)
        cp7, h7, s7 = eval7(rd, T)
        ds = s7 - s9
        dg = (h7 - s7) - (h9 - s9)
        row.append((ds, dg))
    rows[n] = row
    print(f"{n:7s} " + "  ".join(f"{ds:+10.6f} {dg:+10.6f}" for ds, dg in row))

print()
print("# reading: |d(s/R)| within 2e-3 of 0 -> 1-bar-consistent; within 2e-3 of ln(1.01325) in magnitude -> 1-atm-consistent (sign as measured)")
for n in names:
    ds_vals = [ds for ds, dg in rows[n]]
    ds_mean = sum(ds_vals) / len(ds_vals)
    ds_spread = max(ds_vals) - min(ds_vals)
    if abs(ds_mean) < 2e-3:
        verdict = "1-bar-consistent (s matches CEA)"
    elif abs(abs(ds_mean) - LN_1ATM_PER_BAR) < 2e-3:
        verdict = f"1-atm-consistent ({'GRI lower' if ds_mean < 0 else 'GRI higher'} by ~ln(1.01325))"
    else:
        verdict = "neither band (see raw d(s/R))"
    print(f"{n:7s} mean d(s/R) {ds_mean:+10.6f} (spread over T {ds_spread:.2e}) -> {verdict}")
