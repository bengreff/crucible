import csv, sys
# Summary of a run's hot trace: the hottest cell every ~0.01 ms around the crossing, with T, f_H and b.
pre = sys.argv[1]
rows = [r for r in csv.DictReader(open('/tmp/inq2c/dh/%s_hot_trace.csv' % pre)) if None not in r.values()]
f = lambda r, k: float(r[k])
cross = next(n for n, r in enumerate(rows) if f(r, 'T') > 3450)
print(pre, 'rows', len(rows), 'crossing row', cross, 't0 %.4f t1 %.4f ms' % (f(rows[0], 't') * 1e3, f(rows[-1], 't') * 1e3))
last = -1
for n, r in enumerate(rows):
    t = f(r, 't') * 1e3
    if n == cross or t - last >= 0.01 or n == len(rows) - 1:
        last = t
        print('%.4f i %s j %s z %s r %s T %.1f f_H %s b %.3e rho %s ur %s' % (t, r['i'], r['j'], r['z_mm'], r['r_mm'], f(r, 'T'), r['f_H'], f(r, 'excess_J_per_kg'), r['rho'], r['ur']))
