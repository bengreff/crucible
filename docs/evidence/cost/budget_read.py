#!/usr/bin/env python3
"""The one-cell budget's readout (COST.md, lever 3), from crucible_mascotte_a10's <prefix>_budget.csv.

usage: budget_read.py <budget.csv> <t0 s> <t1 s> <i>:<j0>-<j1> [<i>:<j0>-<j1> ...]

For each named cell over the records from the first at or after t0 to the last at or before t1:
(a) closure: the change of rho V, rho E V and rho f_H V against the probe's face and source sums, relative
    to the largest term, and each face's parts (work, conduction, k diffusion, species enthalpy) against
    the energy the cell used (kViscousEnergy), and its tracked diffusion against kTrackedTransport;
(b) the excess b's budget, rho b V = rho E V - rho k V + p V - h_LOX rho V - (h_H2 - h_LOX) rho f_H V, split
    per face and mechanism, with V dp/dt, -d(rho k V)/dt and the sources, and the positive supply's shares
    by the pre-registered groups; and the inviscid faces' b per unit mass flux (each 10-step interval)
    against the range of the two cells beside them.
"""
import csv
import sys

FACE_TERMS = ["inviscid_mass", "inviscid_energy", "inviscid_tracked", "viscous_work", "conduction",
              "turbulent_conduction", "turbulent_k_diffusion", "species_enthalpy", "turbulent_species_enthalpy",
              "tracked_diffusion", "turbulent_tracked_diffusion", "viscous_energy", "tracked_transport"]
SOURCES = ["igniter", "heating", "body_work", "correction", "exchanged_mass", "exchanged_energy", "exchanged_tracked"]
FACES = ["low z", "high z", "below", "above"]


def load(path):
    rows = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            key = (int(r["i"]), int(r["j"]))
            rows.setdefault(key, []).append({k: float(v) for k, v in r.items()})
    return rows


def f(r, face, term):
    return r["f%d_%d" % (face, FACE_TERMS.index(term))]


def s(r, term):
    return r["s_%d" % SOURCES.index(term)]


def window(series, t0, t1):
    inside = [r for r in series if t0 <= r["t"] <= t1]
    return inside[0], inside[-1], inside


def b_of(r):
    return r["rhoEV"] - r["rhokV"] + r["pV"] - r["h_LOX"] * r["rhoV"] - (r["h_H2"] - r["h_LOX"]) * r["rhofHV"]


def main():
    path, t0, t1 = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    cells = []
    for spec in sys.argv[4:]:
        i, js = spec.split(":")
        j0, j1 = js.split("-")
        cells += [(int(i), j) for j in range(int(j0), int(j1) + 1)]
    rows = load(path)
    total_supply = {}
    for cell in cells:
        a, z, inside = window(rows[cell], t0, t1)
        d = lambda name: z[name] - a[name]
        df = lambda face, term: f(z, face, term) - f(a, face, term)
        ds = lambda term: s(z, term) - s(a, term)
        print("cell i %d j %d: t %.6f to %.6f ms, steps %d to %d, T %.0f to %.0f K, excess %.4g to %.4g J/kg"
              % (cell[0], cell[1], a["t"] * 1e3, z["t"] * 1e3, a["step"], z["step"], a["T"], z["T"], a["excess"], z["excess"]))
        # (a) closure
        mass = [df(k, "inviscid_mass") for k in range(4)] + [ds("exchanged_mass")]
        energy = ([df(k, "inviscid_energy") for k in range(4)] + [df(k, "viscous_energy") for k in range(4)]
                  + [ds(t) for t in ("igniter", "heating", "body_work", "correction", "exchanged_energy")])
        tracked = ([df(k, "inviscid_tracked") for k in range(4)] + [df(k, "tracked_transport") for k in range(4)]
                   + [ds("exchanged_tracked")])
        for name, change, terms in (("rho V", d("rhoV"), mass), ("rho E V", d("rhoEV"), energy), ("rho f_H V", d("rhofHV"), tracked)):
            scale = max(abs(x) for x in terms + [change]) or 1
            print("  (a) %-9s change %+.6e, sum of terms %+.6e, residual / largest term %.2e"
                  % (name, change, sum(terms), (change - sum(terms)) / scale))
        for k in range(4):
            parts = sum(df(k, t) for t in ("viscous_work", "conduction", "turbulent_conduction", "turbulent_k_diffusion",
                                           "species_enthalpy", "turbulent_species_enthalpy"))
            used = df(k, "viscous_energy")
            moved = df(k, "tracked_diffusion") + df(k, "turbulent_tracked_diffusion")
            print("      face %-6s viscous energy %+.6e, parts %+.6e (%.1e); tracked %+.6e, parts %+.6e"
                  % (FACES[k], used, parts, (parts - used) / (abs(used) or 1), df(k, "tracked_transport"), moved))
        # (b) the excess
        dh = lambda r: r["h_H2"] - r["h_LOX"]
        hl, dhm = 0.5 * (a["h_LOX"] + z["h_LOX"]), 0.5 * (dh(a) + dh(z))
        if abs(z["h_LOX"] - a["h_LOX"]) > 1e-6 * abs(hl):
            print("  note: h_LOX moved by %.3e J/kg over the window (means used)" % (z["h_LOX"] - a["h_LOX"]))
        items = {}
        for k in range(4):
            items[(FACES[k], "inviscid")] = df(k, "inviscid_energy") - hl * df(k, "inviscid_mass") - dhm * df(k, "inviscid_tracked")
            items[(FACES[k], "molecular conduction")] = df(k, "conduction")
            items[(FACES[k], "molecular species")] = df(k, "species_enthalpy") - dhm * df(k, "tracked_diffusion")
            items[(FACES[k], "turbulent conduction")] = df(k, "turbulent_conduction")
            items[(FACES[k], "turbulent species")] = df(k, "turbulent_species_enthalpy") - dhm * df(k, "turbulent_tracked_diffusion")
            items[(FACES[k], "turbulent k diffusion")] = df(k, "turbulent_k_diffusion")
            items[(FACES[k], "viscous work")] = df(k, "viscous_work")
            # the energy booked to this face by a wall function or an open end beyond its parts (zero at interior faces)
            items[(FACES[k], "unsplit viscous")] = df(k, "viscous_energy") - sum(
                df(k, t) for t in ("viscous_work", "conduction", "turbulent_conduction", "turbulent_k_diffusion",
                                   "species_enthalpy", "turbulent_species_enthalpy"))
        items[("cell", "V dp")] = d("pV")
        items[("cell", "-d(rho k V)")] = -d("rhokV")
        items[("cell", "spray exchange")] = ds("exchanged_energy") - hl * ds("exchanged_mass") - dhm * ds("exchanged_tracked")
        for t in ("igniter", "heating", "body_work", "correction"):
            items[("cell", t)] = ds(t)
        change = b_of(z) - b_of(a)
        # b_of uses each record's own h_LOX; the means above differ by (h_LOX drift) x (rho V and rho f_H V), shown as the residual.
        print("  (b) rho b V change %+.6e J; sum of items %+.6e; residual %+.3e" % (change, sum(items.values()), change - sum(items.values())))
        for key, value in sorted(items.items(), key=lambda kv: -abs(kv[1])):
            if value != 0:
                print("      %-7s %-28s %+.6e" % (key[0], key[1], value))
        positive = {k: v for k, v in items.items() if v > 0}
        supply = sum(positive.values())
        group = lambda pick: sum(v for k, v in positive.items() if pick(k)) / supply if supply > 0 else 0
        interior = lambda k: k[0] in ("high z", "above") or (k[0] == "below" and cell[1] > 0) or (k[0] == "low z" and cell[0] > 0)
        boundary = lambda k: k[0] in FACES and not interior(k)
        shares = {
            "molecular diffusion, interior faces": group(lambda k: interior(k) and k[1].startswith("molecular")),
            "turbulent diffusion, interior faces": group(lambda k: interior(k) and k[1].startswith("turbulent")),
            "plate or axis faces": group(boundary),
            "sources that should be zero": group(lambda k: k[0] == "cell" and k[1] in ("igniter", "heating", "body_work", "correction")),
            "inviscid, interior faces": group(lambda k: interior(k) and k[1] == "inviscid"),
            "V dp": group(lambda k: k == ("cell", "V dp")),
            "other": group(lambda k: not (k[0] in FACES) and k[1] not in ("igniter", "heating", "body_work", "correction", "V dp")
                           or (interior(k) and k[1] in ("viscous work", "unsplit viscous"))),
        }
        print("  positive supply %.6e J; shares: %s" % (supply, ", ".join("%s %.1f%%" % (k, 100 * v) for k, v in shares.items())))
        total_supply[cell] = (supply, shares)
        # inviscid faces' b per unit mass flux against the two cells beside them, each record interval
        neighbour = {0: (cell[0] - 1, cell[1]), 1: (cell[0] + 1, cell[1]), 2: (cell[0], cell[1] - 1), 3: (cell[0], cell[1] + 1)}
        for k in range(4):
            other = rows.get(neighbour[k])
            if other is None:
                continue
            by_step = {r["step"]: r for r in other}
            outside, counted, worst, carried, beyond = 0, 0, 0.0, 0.0, 0.0
            for p, q in zip(inside, inside[1:]):
                m = f(q, k, "inviscid_mass") - f(p, k, "inviscid_mass")
                if q["step"] not in by_step or p["step"] not in by_step or abs(m) < 1e-30:
                    continue
                e = f(q, k, "inviscid_energy") - f(p, k, "inviscid_energy")
                h = f(q, k, "inviscid_tracked") - f(p, k, "inviscid_tracked")
                bf = (e - hl * m - dhm * h) / m
                ends = [p["excess"], q["excess"], by_step[p["step"]]["excess"], by_step[q["step"]]["excess"]]
                lo, hi = min(ends), max(ends)
                counted += 1
                if bf < lo or bf > hi:
                    outside += 1
                    worst = max(worst, (bf - hi) if bf > hi else (lo - bf))
                    carried += max(0.0, m * bf)
                    # the supply beyond the nearer bound: an inflow above the range or an outflow below it adds b
                    beyond += max(0.0, m * (bf - (hi if bf > hi else lo)))
            if counted:
                print("      face %-6s inviscid b per unit mass outside the two cells' range in %d of %d intervals (worst by %.4g J/kg);"
                      " those intervals' positive b supply %.3e J (%.1f%% of the cell's), beyond the range %.3e J (%.1f%%)"
                      % (FACES[k], outside, counted, worst, carried, 100 * carried / supply if supply > 0 else 0,
                         beyond, 100 * beyond / supply if supply > 0 else 0))
    if len(cells) > 1:
        supply = sum(v[0] for v in total_supply.values())
        print("all named cells: positive supply %.6e J; shares: %s" % (supply, ", ".join(
            "%s %.1f%%" % (name, 100 * sum(v[0] * v[1][name] for v in total_supply.values()) / supply)
            for name in next(iter(total_supply.values()))[1])))


if __name__ == "__main__":
    main()
