# RL10-range equilibrium samples for Table B's static interpolation check (TABLE_B.md, RL10 range).
# Per f_H: HP equilibrium at p_c from reactants (H2 + O2 gases at 300 K, as rl10dev; and the same lowered by
# dh_cold for the liquid LOX at 101 K and H2 at 199.67 K, derived below), then shifting-equilibrium isentropic
# expansion to p/p_c. eta as tests/chamber_study.cpp's state samples.
import cantera as ct, sys, numpy as np
mech = sys.argv[1]; out = sys.argv[2]
g = ct.Solution(mech)
names = g.species_names
iH, iO, iN = g.element_index('H'), g.element_index('O'), g.element_index('N')
hInWater = 2 * g.atomic_weight(iH) / g.molecular_weights[g.species_index('H2O')]
def elem_mass(y):
    g.TPY = 300, ct.one_atm, y
    return np.array([g.elemental_mass_fraction(e) for e in g.element_names])
def eta_of(state_y, T, rho):
    z = elem_mass(state_y); zH, zO, zN = z[iH], z[iO], z[iN]
    prod = np.zeros(g.n_species)
    if zH * (1 - hInWater) >= zO * hInWater:
        prod[g.species_index('H2O')] = zO / (1 - hInWater); prod[g.species_index('H2')] = max(zH - prod[g.species_index('H2O')] * hInWater, 0)
    else:
        prod[g.species_index('H2O')] = zH / hInWater; prod[g.species_index('O2')] = max(zO - prod[g.species_index('H2O')] * (1 - hInWater), 0)
    prod[g.species_index('N2')] = zN
    g.TDY = 200.0, rho, prod; eLo = g.int_energy_mass
    g.TDY = 6000.0, rho, prod; eHi = g.int_energy_mass
    g.TDY = T, rho, state_y; e = g.int_energy_mass
    return zH / (zH + zO), zN, (e - eLo) / (eHi - eLo)
# dh_cold (derived): LOX at 101.17 K against O2 gas at 300 K, about -12.6 kJ/mol - 0.054 kJ/mol (O2 gas 298 to 300 K)
# = -0.396 MJ/kg of O2; H2 at 199.67 K against 300 K, cp about 14.0 kJ/(kg K) x -100 K = -1.40 MJ/kg of H2.
def run(fH, pc, cold, ratios):
    of = (1 - fH) / fH
    y = np.zeros(g.n_species); y[g.species_index('H2')] = fH; y[g.species_index('O2')] = 1 - fH
    g.TPY = 300, pc, y
    if cold:
        h = g.enthalpy_mass + (-0.396e6) * (1 - fH) + (-1.40e6) * fH
        g.HPY = h, pc, y
    g.equilibrate('HP'); s0 = g.entropy_mass; yc = g.Y.copy()
    rows = []
    for r in ratios:
        g.SPY = s0, pc * r, g.Y; g.equilibrate('SP')
        T, rho, Y = g.T, g.density, g.Y.copy()
        f, zn, eta = eta_of(Y, T, rho)
        rows.append((f, zn, eta, rho, T, Y))
        g.TDY = T, rho, Y
    return rows
fHs = [float(v) for v in sys.argv[3].split(',')]
ratios = [float(v) for v in sys.argv[4].split(',')]
pcs = [float(v) for v in sys.argv[5].split(',')]
with open(out, 'w') as fo:
    fo.write('t,steps,i,j,f_H,Z_N,eta,rho,T,' + ','.join('Y_' + n for n in names) + '\n')
    k = 0
    for fH in fHs:
        for pc in pcs:
            for cold in (False, True):
                for (f, zn, eta, rho, T, Y) in run(fH, pc, cold, ratios):
                    fo.write('%.6e,1,%d,0,%.9e,%.9e,%.9e,%.9e,%.9e,' % (0.0, k, f, zn, eta, rho, T) + ','.join('%.9e' % v for v in Y) + '\n')
                    k += 1
print(out, k, 'samples')
