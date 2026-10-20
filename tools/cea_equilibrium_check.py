"""The equilibrium check stated in docs/evidence/TABLE_B.md (*The 3,500 K ceiling*): CEA's chamber
(HP) problem, through RocketCEA, against Cantera's HP equilibrium (build/crucible_equilibrium_hp) on
the same reactants: gaseous H2 and O2, both at 300, 1,500 or 3,000 K, O/F 1 to 12, 0.1 to 20 MPa.
Judged on the cases whose CEA chamber temperature is 3,000 to 5,000 K: pass if |dT| <= 3 K on the
first mechanism. The others are printed beside, not judged. Major mole fractions are reported.
Usage: python cea_equilibrium_check.py <crucible_equilibrium_hp> <judged.yaml> [other.yaml ...]"""
import subprocess, sys
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer

PSIA = 6894.757293168  # Pa
T0S = [300.0, 1500.0, 3000.0]
OFS = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0]
PS = [0.1e6, 0.3e6, 1e6, 3e6, 10e6, 20e6]
SPECIES = ["H2", "H2O", "OH", "H", "O", "O2"]
TLO, THI, BAR = 3000.0, 5000.0, 3.0

exe, mechs = sys.argv[1], sys.argv[2:]
cases = [(t0, of, p) for t0 in T0S for of in OFS for p in PS]

cea = {}
for t0 in T0S:
    tag = f"{int(t0)}"
    add_new_fuel(f"H2_G{tag}", f"fuel H2 H 2 wt%=100. t(k)={t0}")
    add_new_oxidizer(f"O2_G{tag}", f"oxid O2 O 2 wt%=100. t(k)={t0}")
    c = CEA_Obj(oxName=f"O2_G{tag}", fuelName=f"H2_G{tag}")
    for of in OFS:
        for p in PS:
            pc = p / PSIA
            tc = c.get_Tcomb(Pc=pc, MR=of) / 1.8
            _, x = c.get_SpeciesMoleFractions(Pc=pc, MR=of, eps=2.0, frozen=0, frozenAtThroat=0, min_fraction=1e-10)
            x = {k.lstrip("*"): v[1] for k, v in x.items()}  # index 1 is the chamber (0 is the injector)
            cea[(t0, of, p)] = (tc, [x.get(s, 0.0) for s in SPECIES])

stdin = "".join(f"{t0} {of} {p}\n" for t0, of, p in cases)
can = []
for m in mechs:
    rows = [l.split() for l in subprocess.run([exe, m], input=stdin, capture_output=True, text=True, check=True).stdout.splitlines()
            if not l.startswith("#")]
    can.append({(t0, of, p): (float(r[3]), [float(v) for v in r[4:]]) for (t0, of, p), r in zip(cases, rows)})

names = [m.split("/")[-1] for m in mechs]
print(f"# CEA (RocketCEA, HP chamber problem) against Cantera HP equilibrium; judged: {names[0]}; beside: {', '.join(names[1:]) or '-'}")
print(f"# dT = T_Cantera - T_CEA (K); dX = largest |X_Cantera - X_CEA| over {', '.join(SPECIES)}")
head = f"{'T0':>6} {'O/F':>5} {'p_MPa':>6} {'T_CEA':>8}" + "".join(f" {'dT':>7} {'dX':>8} {'arg':>4}" for _ in names)
print(head + "   (" + " | ".join(names) + ")")
worst = [(0.0, None) for _ in names]
worstx = [(0.0, None) for _ in names]
n = 0
for k in cases:
    tc, xc = cea[k]
    inrange = TLO <= tc <= THI
    n += inrange
    line = f"{k[0]:6.0f} {k[1]:5.1f} {k[2]/1e6:6.1f} {tc:8.2f}"
    for i, d in enumerate(can):
        t, x = d[k]
        dx = [abs(a - b) for a, b in zip(x, xc)]
        j = max(range(len(SPECIES)), key=lambda q: dx[q])
        line += f" {t - tc:7.2f} {dx[j]:8.2e} {SPECIES[j]:>4}"
        if inrange and abs(t - tc) > abs(worst[i][0]):
            worst[i] = (t - tc, k)
        if inrange and dx[j] > worstx[i][0]:
            worstx[i] = (dx[j], (k, SPECIES[j]))
    print(line + ("" if inrange else "   (outside 3,000-5,000 K)"))
print(f"\n# {n} of {len(cases)} cases have CEA Tc in {TLO:.0f} to {THI:.0f} K")
for i, m in enumerate(names):
    (dt, k), (dx, kx) = worst[i], worstx[i]
    verdict = ("PASS" if abs(dt) <= BAR else "FAIL") if i == 0 else "(not judged)"
    print(f"# {m}: worst dT {dt:+.2f} K at T0 {k[0]:.0f} K, O/F {k[1]}, {k[2]/1e6} MPa; worst dX {dx:.2e} ({kx[1]}) at "
          f"T0 {kx[0][0]:.0f} K, O/F {kx[0][1]}, {kx[0][2]/1e6} MPa  {verdict}")
