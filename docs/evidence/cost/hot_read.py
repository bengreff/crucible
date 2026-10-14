import csv, sys
pre = sys.argv[1]
rows = [r for r in csv.DictReader(open(pre + '_hot_trace.csv')) if None not in r.values()]
f = lambda r, k: float(r[k])
print('rows', len(rows))
# cells visited as hottest, with counts
from collections import Counter
c = Counter((r['i'], r['j']) for r in rows)
print('hottest cells (i,j):count', c.most_common(8))
cross = next(n for n, r in enumerate(rows) if f(r, 'T') > 3450)
print('crossing row', cross)
cols = ['t', 'step', 'i', 'j', 'z_mm', 'r_mm', 'T', 'rho', 'uz', 'ur', 'Y_H2', 'Y_O2', 'Y_H2O', 'Y_OH', 'f_H', 'excess_J_per_kg',
        'spray_mass_per_gas_mass', 'spray_energy_J_per_kg', 'parcels']
print(','.join(cols))
idx = sorted(set(list(range(0, len(rows), 200)) + list(range(max(0, cross - 20), min(len(rows), cross + 6))) + [len(rows) - 1]))
for n in idx:
    r = rows[n]
    print(','.join(r[k] if k not in ('t',) else '%.6f' % (f(r, k) * 1e3) for k in cols))
