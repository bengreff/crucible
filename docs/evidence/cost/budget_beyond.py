#!/usr/bin/env python3
"""The inviscid b supply beyond the two cells' range, over a whole probe block (COST.md, lever 3).

usage: budget_beyond.py <budget.csv> <t0 s> <t1 s>

For every probed cell and each face with a probed cell beside it, over 10-step intervals from the first record
at or after t0 to the last at or before t1: the face's b per unit mass flux against the range of the two cells'
excess at the interval's ends, as budget_read.py's (b) does. Prints the block's positive supply of b, the supply
beyond the range (its share), the faces' worst excursion, and the five cells with the largest beyond-range share.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from budget_read import FACES, b_of, f, load, s, window  # noqa: E402


def main():
    path, t0, t1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    rows = load(path)
    total_supply, total_beyond, worst_all, per_cell = 0.0, 0.0, 0.0, []
    for cell, series in sorted(rows.items()):
        a, z, inside = window(series, t0, t1)
        df = lambda face, term: f(z, face, term) - f(a, face, term)
        ds = lambda term: s(z, term) - s(a, term)
        dh = lambda r: r["h_H2"] - r["h_LOX"]
        hl, dhm = 0.5 * (a["h_LOX"] + z["h_LOX"]), 0.5 * (dh(a) + dh(z))
        items = []
        for k in range(4):
            items.append(df(k, "inviscid_energy") - hl * df(k, "inviscid_mass") - dhm * df(k, "inviscid_tracked"))
            items += [df(k, "conduction"), df(k, "turbulent_conduction"), df(k, "turbulent_k_diffusion"), df(k, "viscous_work"),
                      df(k, "species_enthalpy") - dhm * df(k, "tracked_diffusion"),
                      df(k, "turbulent_species_enthalpy") - dhm * df(k, "turbulent_tracked_diffusion"),
                      df(k, "viscous_energy") - sum(df(k, t) for t in ("viscous_work", "conduction", "turbulent_conduction",
                                                                        "turbulent_k_diffusion", "species_enthalpy",
                                                                        "turbulent_species_enthalpy"))]
        items += [z["pV"] - a["pV"], -(z["rhokV"] - a["rhokV"]),
                  ds("exchanged_energy") - hl * ds("exchanged_mass") - dhm * ds("exchanged_tracked")]
        items += [ds(t) for t in ("igniter", "heating", "body_work", "correction")]
        supply = sum(v for v in items if v > 0)
        neighbour = {0: (cell[0] - 1, cell[1]), 1: (cell[0] + 1, cell[1]), 2: (cell[0], cell[1] - 1), 3: (cell[0], cell[1] + 1)}
        beyond_cell, worst_cell, worst_face = 0.0, 0.0, ""
        for k in range(4):
            other = rows.get(neighbour[k])
            if other is None:
                continue
            by_step = {r["step"]: r for r in other}
            beyond = 0.0
            for p, q in zip(inside, inside[1:]):
                m = f(q, k, "inviscid_mass") - f(p, k, "inviscid_mass")
                if q["step"] not in by_step or p["step"] not in by_step or abs(m) < 1e-30:
                    continue
                e = f(q, k, "inviscid_energy") - f(p, k, "inviscid_energy")
                h = f(q, k, "inviscid_tracked") - f(p, k, "inviscid_tracked")
                bf = (e - hl * m - dhm * h) / m
                ends = [p["excess"], q["excess"], by_step[p["step"]]["excess"], by_step[q["step"]]["excess"]]
                lo, hi = min(ends), max(ends)
                if bf < lo or bf > hi:
                    excursion = (bf - hi) if bf > hi else (lo - bf)
                    if excursion > worst_cell:
                        worst_cell, worst_face = excursion, FACES[k]
                    beyond += max(0.0, m * (bf - (hi if bf > hi else lo)))
            beyond_cell += beyond
        total_supply += supply
        total_beyond += beyond_cell
        worst_all = max(worst_all, worst_cell)
        per_cell.append((beyond_cell / supply if supply > 0 else 0.0, cell, beyond_cell, supply, worst_cell, worst_face))
    first = min(min(r["t"] for r in v if r["t"] >= t0) for v in rows.values())
    last = max(max(r["t"] for r in v if r["t"] <= t1) for v in rows.values())
    print("%s: %d cells, t %.4f to %.4f ms" % (path, len(rows), first * 1e3, last * 1e3))
    print("  positive supply of b %.4e J; inviscid supply beyond the two cells' range %.4e J (%.1f%%); worst excursion %.3g MJ/kg"
          % (total_supply, total_beyond, 100 * total_beyond / total_supply if total_supply > 0 else 0, worst_all / 1e6))
    for share, cell, beyond, supply, worst, face in sorted(per_cell, reverse=True)[:5]:
        print("  cell i %d j %d: beyond %.3e J of supply %.3e J (%.1f%%), worst %.3g MJ/kg (face %s)"
              % (cell[0], cell[1], beyond, supply, 100 * share, worst / 1e6, face))


if __name__ == "__main__":
    main()
