#!/usr/bin/env python3
"""The S7 band (docs/validation/MASCOTTE_A10_PREDICTION.md, sections 5 and 6) from the four runs' outputs.

usage: mascotte_band.py <low prefix> <high prefix> <blob prefix> <blob-fine prefix> <out csv>

P1, at each station and ring (r <= 25 mm):
  band = [min, max] of the three production runs' averaged T, widened on each side by
  |T(blob-fine) - T(blob)| at that point and by the largest |T_first_half - T_second_half| of the
  three runs at that point.
P2, from the averaged OH:
  the peak of the radially integrated OH per unit length, the downstream half-maximum, and the radius
  of the largest OH concentration (rho Y_OH, linear in z between columns) at x = 50, 100 and 200 mm;
  the range is the three runs' min to max, widened by 20% (min * 0.8 to max * 1.2).
Prints the markdown tables for section 6 and writes the full band to <out csv>.
"""
import csv
import sys
from collections import defaultdict


def stations(prefix):
    rows = {}
    with open(prefix + "_stations.csv") as f:
        for r in csv.DictReader(f):
            rows[(float(r["z_mm"]), float(r["r_mm"]))] = (
                float(r["T"]), float(r["T_first_half"]), float(r["T_second_half"]))
    return rows


def oh_line(prefix):
    z, v = [], []
    with open(prefix + "_oh.csv") as f:
        for r in csv.DictReader(f):
            z.append(float(r["z_mm"]))
            v.append(float(r["OH_kg_per_m"]))
    i = max(range(len(v)), key=v.__getitem__)
    half = None
    for k in range(i, len(v) - 1):  # downstream half-maximum, linear between columns
        if v[k + 1] <= 0.5 * v[i]:
            half = z[k] + (0.5 * v[i] - v[k]) / (v[k + 1] - v[k]) * (z[k + 1] - z[k])
            break
    return z[i], half


def oh_radius(prefix, x):
    columns = defaultdict(dict)
    with open(prefix + "_field.csv") as f:
        for r in csv.DictReader(f):
            columns[float(r["z_mm"])][float(r["r_mm"])] = float(r["rhoY_OH"])
    zs = sorted(columns)
    k = max(j for j in range(len(zs) - 1) if zs[j] <= x)
    w = (x - zs[k]) / (zs[k + 1] - zs[k])
    a, b = columns[zs[k]], columns[zs[k + 1]]
    return max(a, key=lambda r: (1 - w) * a[r] + w * b[r])


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    low, high, blob, fine, out = sys.argv[1:]
    runs = {"low": stations(low), "high": stations(high), "blob": stations(blob)}
    finer = stations(fine)
    for name, s in runs.items():
        drift = max(abs(v[1] - v[2]) for (z, r), v in s.items() if r <= 20.0)
        print(f"{name}: largest half-to-half drift at r <= 20 mm {drift:.1f} K (settled if <= 50 K)")
    band = []
    for key in sorted(runs["blob"]):
        z, r = key
        if r > 25.0:
            continue
        t = [runs[n][key][0] for n in runs]
        grid = abs(finer[key][0] - runs["blob"][key][0])
        drift = max(abs(runs[n][key][1] - runs[n][key][2]) for n in runs)
        band.append((z, r, *t, grid, drift, min(t) - grid - drift, max(t) + grid + drift))
    with open(out, "w") as f:
        f.write("z_mm,r_mm,T_low,T_high,T_blob,grid,drift,band_lo,band_hi\n")
        for b in band:
            f.write(",".join(f"{v:.4f}" if i < 2 else f"{v:.1f}" for i, v in enumerate(b)) + "\n")
    print("\nP1 (K), every 2 mm (nearest ring):")
    print("| x (mm) | r (mm) | low | high | blob | grid | drift | band |")
    print("|---|---|---|---|---|---|---|---|")
    for z in sorted({b[0] for b in band}):
        rings = [b for b in band if b[0] == z]
        for target in range(0, 26, 2):
            b = min(rings, key=lambda b: abs(b[1] - target))
            print(f"| {z:.0f} | {b[1]:.1f} | {b[2]:.0f} | {b[3]:.0f} | {b[4]:.0f} | {b[5]:.0f} | {b[6]:.0f} "
                  f"| {b[7]:.0f} to {b[8]:.0f} |")
    print("\nP2 (mm):")
    prefixes = {"low": low, "high": high, "blob": blob, "blob-fine": fine}
    rows = {n: (*oh_line(p), *(oh_radius(p, x) for x in (50.0, 100.0, 200.0))) for n, p in prefixes.items()}
    names = ["OH peak x", "half-maximum x", "OH radius at 50", "OH radius at 100", "OH radius at 200"]
    print("| quantity | low | high | blob | blob-fine | range ±20% |")
    print("|---|---|---|---|---|---|")
    for q, name in enumerate(names):
        v = [rows[n][q] for n in ("low", "high", "blob")]
        cell = lambda x: "none" if x is None else f"{x:.1f}"
        rng = "n/a" if None in v else f"{0.8 * min(v):.1f} to {1.2 * max(v):.1f}"
        print(f"| {name} | {cell(v[0])} | {cell(v[1])} | {cell(v[2])} | {cell(rows['blob-fine'][q])} | {rng} |")


if __name__ == "__main__":
    main()
