import csv, sys
rows = list(csv.DictReader(open(sys.argv[1] + '_hot_cells.csv')))
I = sorted({int(r['i']) for r in rows}); J = sorted({int(r['j']) for r in rows})
g = {(int(r['i']), int(r['j'])): r for r in rows}
for k, fmt in [('T', '%6.0f'), ('rho', '%6.2f'), ('f_H', '%6.3f'), ('excess_J_per_kg', '%6.0f'), ('Y_O2', '%6.3f'), ('Y_H2', '%6.3f'), ('Y_OH', '%6.3f'), ('spray_energy_J_per_kg', '%6.0f'), ('parcels', '%6s')]:
    sc = 1e-3 if k == 'excess_J_per_kg' or k == 'spray_energy_J_per_kg' else 1
    print('==', k + (' [kJ/kg]' if sc != 1 else ''), ' rows j (radial, down) x cols i (axial):', I)
    for j in reversed(J):
        print('%3d ' % j + ' '.join((fmt % (float(g[(i, j)][k]) * sc) if k != 'parcels' else fmt % g[(i, j)][k]) if (i, j) in g else '      ' for i in I))
