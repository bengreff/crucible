#!/usr/bin/env python3
"""Diagnosis only (not a criterion's instrument): criterion 4's differences with the explicit's period-mean rates
(outlet flow, F_vac) stamped at their periods' midpoints instead of their ends. p_inj, a state value, keeps its time."""
import csv, gzip, math, sys
def load(p):
    with (gzip.open(p, "rt") if p.endswith(".gz") else open(p)) as f: rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}
def at(ts, ys, t):
    lo, hi = 0, len(ts) - 1
    if t <= ts[0]: return ys[0]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ts[mid] <= t: lo = mid
        else: hi = mid
    return ys[lo] + (t - ts[lo]) / (ts[hi] - ts[lo]) * (ys[hi] - ys[lo])
ex = load(sys.argv[1]); t = ex["t"]
mid = [t[0]] + [(a + b) / 2 for a, b in zip(t, t[1:])]
norm = {}
for c in ("p_injector", "outlet", "F_vac"):
    v = [y for tt, y in zip(t, ex[c]) if 6e-3 - 1e-12 <= tt <= 8e-3 + 1e-12]; norm[c] = sum(v) / len(v)
print("sample periods (s): min %.3e max %.3e" % (min(b - a for a, b in zip(t, t[1:])), max(b - a for a, b in zip(t, t[1:]))))
for p in sys.argv[2:]:
    h = load(p); print("==", p)
    for c in ("p_injector", "outlet", "F_vac"):
        for name, ts in (("end", t), ("midpoint", mid)):
            if c == "p_injector" and name == "midpoint": continue
            d = [(y - at(ts, ex[c], tt)) / norm[c] for tt, y in zip(h["t"], h[c]) if 0.5e-3 - 1e-12 <= tt <= 8e-3 + 1e-12]
            print("  %-10s stamped at %-8s RMS %.5f%%  largest %.5f%%  (%d times)" % (c, name, 100 * math.sqrt(sum(x * x for x in d) / len(d)), 100 * max(map(abs, d)), len(d)))
