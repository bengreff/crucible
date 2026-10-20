#!/usr/bin/env python3
"""The probe block's excess (MJ/kg) and T (K) at chosen times, from crucible_mascotte_a10's <prefix>_budget.csv.

usage: budget_block.py <budget.csv> <t us> [<t us> ...]   (each time: the nearest record)
"""
import csv
import sys

rows = {}
for r in csv.DictReader(open(sys.argv[1])):
    rows.setdefault(float(r["t"]) * 1e6, {})[(int(r["i"]), int(r["j"]))] = (float(r["excess"]), float(r["T"]))
times = sorted(rows)
cells = rows[times[0]]
ni, nj = max(i for i, _ in cells) + 1, max(j for _, j in cells) + 1
for want in map(float, sys.argv[2:]):
    t = min(times, key=lambda x: abs(x - want))
    print("t %.1f us: excess MJ/kg (T K), j 0..%d across, i 0..%d down" % (t, nj - 1, ni - 1))
    for i in range(ni):
        print("  i%d " % i + " ".join("%5.2f(%4.0f)" % (rows[t][(i, j)][0] / 1e6, rows[t][(i, j)][1]) for j in range(nj)))
