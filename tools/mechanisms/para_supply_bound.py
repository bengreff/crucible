"""The bound on the supply's para-rich hydrogen before any conversion model (docs/evidence/LOW_T_THERMO.md, (d)).
CEA-style bookkeeping on h2o2_highT.yaml's species: an HP equilibrium chamber at p_c (an infinite-area combustor),
then a shifting-equilibrium isentropic one-dimensional expansion, u = sqrt(2 (h_c - h(s_c, p))), G = rho u; the
throat is G's maximum over p (golden section), c* = p_c / G_t, and the vacuum Isp at area ratio eps is
u_e + p_e / (rho_e u_e) on the supersonic branch. Each case runs twice: the H2 supply with its equilibrium
enthalpy h_eq(T_s), and with the frozen para-rich enthalpy h_eq(T_s) - dh(T_s, x_p). dh comes from the
rotational partition functions on make_low_t.py's v = 0 levels (Huber and Herzberg 1979), with a common zero at
para's J = 0; translation and vibration are the same in both forms and cancel. Only relative changes are printed.
Usage: ~/venvs/crucible/bin/python tools/mechanisms/para_supply_bound.py"""
import math, os
import cantera as ct

here = os.path.dirname(os.path.abspath(__file__))
gas = ct.Solution(os.path.join(here, "h2o2_highT.yaml"))
R = ct.gas_constant  # J/(kmol K)
M = dict(zip(gas.species_names, gas.molecular_weights))  # kg/kmol
c2 = 1.438776877  # second radiation constant hc/k, cm K
E = [c2 * (59.3392 * j * (j + 1) - 0.04599 * (j * (j + 1)) ** 2 + 5.2e-5 * (j * (j + 1)) ** 3) for j in range(20)]

def rotation(T, js, weight):
    """Partition function and mean rotational energy (K) over the levels js."""
    z = m1 = 0.0
    for j in js:
        g = weight * (2 * j + 1) * math.exp(-E[j] / T)
        z += g; m1 += g * E[j]
    return z, m1 / z

def e_rot(T, xp=None):
    """Mean rotational energy (K): equilibrium if xp is None, else ortho and para frozen at para fraction xp."""
    zp, ep = rotation(T, range(0, 20, 2), 1)
    zo, eo = rotation(T, range(1, 20, 2), 3)
    if xp is None:
        return (zp * ep + zo * eo) / (zp + zo)
    return xp * ep + (1 - xp) * eo

def dh(T, xp):
    """h_eq(T) - h_frozen(T, xp), J/kg of H2."""
    return (e_rot(T) - e_rot(T, xp)) * R / M["H2"]

def h_h2_eq(T):
    """The equilibrium H2 supply's enthalpy, J/kg: the mechanism's (TM-4513) at and above 200 K; below, its
    200 K value plus the equilibrium partition function's change."""
    if T >= 200.0:
        gas.TPX = T, ct.one_atm, "H2:1"
        return gas.enthalpy_mass
    gas.TPX = 200.0, ct.one_atm, "H2:1"
    return gas.enthalpy_mass + (2.5 * (T - 200.0) + e_rot(T) - e_rot(200.0)) * R / M["H2"]

def h_o2(state):
    if state == "LOX 90.17 K":
        return -12979.0e3 / M["O2"]  # CEA's O2(L) at 90.17 K, -12.979 kJ/mol
    gas.TPX = 298.15, ct.one_atm, "O2:1"
    return gas.enthalpy_mass

def heat_release(h0, y_h2, y_o2):
    """h0 minus the enthalpy at 298.15 K of the complete-combustion products (H2O gas and the excess H2)."""
    n_o2 = y_o2 / M["O2"]
    n_h2 = y_h2 / M["H2"] - 2 * n_o2
    assert n_h2 > 0
    gas.TPY = 298.15, ct.one_atm, {"H2O": 2 * n_o2 * M["H2O"], "H2": n_h2 * M["H2"]}
    return h0 - gas.enthalpy_mass

def performance(h0, y_h2, y_o2, pc, eps_list):
    gas.HPY = h0, pc, {"H2": y_h2, "O2": y_o2}
    gas.equilibrate("HP")
    tc, yc, sc = gas.T, gas.Y, gas.entropy_mass
    def state(p):
        gas.TPY = tc, pc, yc
        gas.SP = sc, p
        gas.equilibrate("SP")
        u = math.sqrt(max(2.0 * (h0 - gas.enthalpy_mass), 0.0))
        return gas.density * u, u, gas.density
    a, b = 0.4 * pc, 0.75 * pc
    g = (math.sqrt(5) - 1) / 2
    x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = state(x1)[0], state(x2)[0]
    while b - a > 1e-10 * pc:
        if f1 > f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - g * (b - a); f1 = state(x1)[0]
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + g * (b - a); f2 = state(x2)[0]
    pt = 0.5 * (a + b)
    gt = state(pt)[0]
    isp = []
    for eps in eps_list:
        lo, hi = math.log(pt * 1e-7), math.log(pt)  # supersonic branch: G falls as p falls below p_t
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if gt / state(math.exp(mid))[0] > eps:
                lo = mid
            else:
                hi = mid
        pe = math.exp(0.5 * (lo + hi))
        ge, ue, rhoe = state(pe)
        isp.append(ue + pe / (rhoe * ue))
    return pc / gt, isp

print("check: dh(199.67 K, x_p 1) = %.1f kJ/kg (LOW_T_THERMO.md (b): 163)" % (dh(199.67, 1.0) / 1e3))
print("dh = h_eq - h_frozen, kJ/kg of H2, by T_s:")
for T in (100.0, 150.0, 200.0, 250.0, 300.0):
    print("  %5.0f K: x_p 0.998 %6.1f, x_p 0.95 %6.1f, equilibrium para fraction %.4f"
          % (T, dh(T, 0.998) / 1e3, dh(T, 0.95) / 1e3,
             rotation(T, range(0, 20, 2), 1)[0] / (rotation(T, range(0, 20, 2), 1)[0] + rotation(T, range(1, 20, 2), 3)[0])))
eps_list = (20.0, 60.0, 100.0)
rows = []
print("O/F  p_c MPa  T_s K  x_p    O2            dh kJ/kg  f=Y_H2 dh/q  dc*/c*     dIsp/Isp eps 20, 60, 100")
for o2 in ("LOX 90.17 K", "GOX 298.15 K"):
    for of in (4.0, 5.0, 6.0):
        y_h2, y_o2 = 1 / (1 + of), of / (1 + of)
        for pc in (2e6, 3e6, 4e6):
            for ts in (150.0, 200.0, 250.0, 300.0):
                h0 = y_h2 * h_h2_eq(ts) + y_o2 * h_o2(o2)
                q = heat_release(h0, y_h2, y_o2)
                cstar0, isp0 = performance(h0, y_h2, y_o2, pc, eps_list)
                for xp in (0.95, 0.998):
                    d = dh(ts, xp)
                    cstar1, isp1 = performance(h0 - y_h2 * d, y_h2, y_o2, pc, eps_list)
                    row = (of, pc / 1e6, ts, xp, o2, d / 1e3, y_h2 * d / q, cstar1 / cstar0 - 1,
                           [i1 / i0 - 1 for i0, i1 in zip(isp0, isp1)], q)
                    rows.append(row)
                    print("%3.0f  %5.0f  %5.0f  %.3f  %-12s  %7.1f  %10.3e  %10.3e  %s"
                          % (row[:7] + (row[7], ", ".join("%.3e" % x for x in row[8]))))
def worst(rs):
    return max(max(r[6], abs(r[7]), max(abs(x) for x in r[8])) for r in rs)
for ts in (150.0, 200.0, 250.0, 300.0):
    rs = [r for r in rows if r[2] == ts]
    print("T_s %3.0f K (%d cases): largest f %.3e, largest |dc*/c*| %.3e, largest |dIsp/Isp| %.3e, largest of the three %.3e%s"
          % (ts, len(rs), max(r[6] for r in rs), max(abs(r[7]) for r in rs), max(max(abs(x) for x in r[8]) for r in rs),
             worst(rs), "  <- judged against 2e-4" if ts == 200.0 else "  (reported, not judged)"))
rs = [r for r in rows if r[2] == 200.0]
print("judgement at T_s 200 K: largest %.3e against 2e-4: %s" % (worst(rs), "below" if worst(rs) < 2e-4 else "NOT below"))

# Beyond the registration (derived, not judged): what an in-chamber conversion rate can move.
print("not judged: h_eq - h_frozen(x_p 0.998) as the gas heats, kJ/kg of H2:",
      ", ".join("%.0f K %.3g" % (T, dh(T, 0.998) / 1e3) for T in (300.0, 400.0, 500.0, 700.0, 1000.0)))
for tstar in (500.0, 1000.0):
    print("not judged: equilibration taken as instantaneous above %.0f K: largest Y_H2 dh(%.0f K, x_p 0.998) / q over the grid %.2e"
          % (tstar, tstar, max(dh(tstar, 0.998) / (1 + r[0]) / r[9] for r in rows)))
for T in (20.0, 20.27):
    zp, zo = rotation(T, range(0, 20, 2), 1)[0], rotation(T, range(1, 20, 2), 3)[0]
    print("not judged: equilibrium para fraction at %.2f K %.5f" % (T, zp / (zp + zo)))
for ts in (150.0, 200.0, 250.0, 300.0):
    target = h_h2_eq(ts) - dh(ts, 0.998)  # the para-rich supply's enthalpy, J/kg
    lo, hi = 20.0, ts
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if h_h2_eq(mid) > target:
            hi = mid
        else:
            lo = mid
    print("not judged: T_s %3.0f K, x_p 0.998, equilibrated at the inlet at the same enthalpy: %.1f K (%.1f K colder, "
          "density %.2f%% higher at the same pressure)" % (ts, mid, ts - mid, 100 * (ts / mid - 1)))
