"""T, the excess b and f_H around a hot-trace crossing, from a driver's <prefix>_hot_cells.csv (COST.md, lever 3).

Usage: python3 hot_cells_map.py <prefix>_hot_cells.csv [title]
"""
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1])))
grid = {(int(r["i"]), int(r["j"])): r for r in rows}
I = sorted({k[0] for k in grid})
J = sorted({k[1] for k in grid})
if len(sys.argv) > 2:
    print("# " + sys.argv[2])
for name, key, fmt, scale in [("T (K)", "T", "{:5.0f}", 1), ("excess b (MJ/kg)", "excess_J_per_kg", "{:5.2f}", 1e-6), ("f_H", "f_H", "{:5.3f}", 1)]:
    print(f"{name}; rows j, columns i {I}")
    for j in J:
        print(f"{j:2d}", " ".join(fmt.format(float(grid[(i, j)][key]) * scale) for i in I))
