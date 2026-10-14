"""RL10A-3-3A size check (docs/validation/RL10A33A_PROTOCOL.md, section 1, item (c)): the smooth lateral area,
chamber volume and L* of the protocol's PCHIP wall under the two readings of Table 2.5.1's 2.47 in and 5.13 in, to set
beside Table 2.4.1's "Hot-gas-side Effective Surface Area" (7200 in^2). Geometry only; no performance value.

Usage: .venv-cea/bin/python docs/validation/rl10a33a/surface_area.py
"""
import csv
import os

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))
rows = [line for line in open(os.path.join(HERE, 'tm107318_table_E1.csv')) if not line.startswith('#')]
data = list(csv.reader(rows))[1:]
z = [float(x[0]) for x in data]
eps = [float(x[1]) for x in data]
for rt, rc in ((1.235, 2.565), (2.47, 5.13)):
    for face, exit_ in ((-13.1, 46.0), (-13.5, 45.64)):
        pts = sorted(list(zip(z, [rt * np.sqrt(e) for e in eps])) + [(0.0, rt), (face, rc), (exit_, rt * np.sqrt(61))])
        p = PchipInterpolator([q[0] for q in pts], [q[1] for q in pts])
        d = p.derivative()
        knots = [q[0] for q in pts]
        spans = list(zip(knots[:-1], knots[1:]))
        lateral = lambda s: 2 * np.pi * p(s) * np.sqrt(1 + d(s) ** 2)
        area = sum(quad(lateral, a, b, limit=200)[0] for a, b in spans)
        chamber = sum(quad(lateral, a, b, limit=200)[0] for a, b in spans if b <= 0)
        volume = sum(quad(lambda s: np.pi * p(s) ** 2, a, b, limit=200)[0] for a, b in spans if b <= 0)
        print(f"r_t {rt} in, r_c {rc} in, face {face} in, exit {exit_} in: wall area face to exit {area:.0f} in2 "
              f"(7200 / area = {7200 / area:.2f}), face to throat {chamber:.0f} in2, chamber volume to the throat "
              f"{volume:.0f} in3, L* {volume / (np.pi * rt * rt):.1f} in")
