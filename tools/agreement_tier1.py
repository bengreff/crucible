#!/usr/bin/env python3
"""SPEED_PLAN.md 7.2: the implicit step's agreement with the explicit reference on tier 1.

Usage: agreement_tier1.py DIR [implicit run names...]   (default i10 i5)
DIR holds e1/, e2/ and the implicit runs, each with <name>_DEVELOPMENT_series.csv and <name>.log.
Plain Python (no numpy). Prints the window L and window H tables, the budgets and the speed-ups.
"""
import math
import re
import sys

OBS = ["mdot_kg_s", "fvac_N", "p0_Pa", "pface_Pa", "qwall_W"]
REPORTED = ["qface_W"]
L = (0.5e-3, 20e-3)
H = (23e-3, 30e-3)


def series(path):
    rows = {}
    head = None
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cells = line.strip().split(",")
            if head is None:
                head = cells
                continue
            if len(cells) != len(head):
                continue  # a row still being written
            r = dict(zip(head, (float(c) for c in cells)))
            rows[round(r["t_s"] / 1e-5)] = r
    return rows


def window(rows, w):
    return {k: r for k, r in rows.items() if w[0] - 1e-12 <= r["t_s"] <= w[1] + 1e-12}


def log_facts(path):
    budgets, walls = [], []
    with open(path) as f:
        text = f.read()
    for m in re.finditer(r"--- t ([0-9.]+) ms.*?wall ([0-9.]+) s", text, re.S):
        walls.append((float(m.group(1)) * 1e-3, float(m.group(2))))
    for m in re.finditer(r"budgets: \(a\) gas ([-0-9.e+naif]+)", text):
        budgets.append(float(m.group(1)))
    end = None
    m = re.search(r"\nreal ([0-9.]+)", text)
    if m:
        end = float(m.group(1))
    err = re.findall(r"^error: (.*)$", text, re.M)
    return budgets, walls, end, err


def wall_at(walls, t):
    for tt, w in walls:
        if abs(tt - t) < 1e-9:
            return w
    return None


def main():
    root = sys.argv[1]
    names = sys.argv[2:] or ["i10", "i5"]
    ref = series(f"{root}/e1/e1_DEVELOPMENT_series.csv")
    tmax_ref = max(r["t_s"] for r in ref.values())
    print(f"e1 reaches {tmax_ref * 1e3:.2f} ms")

    def nonfinite(rows, w=None):
        rs = window(rows, w) if w else rows
        return sum(1 for r in rs.values() for o in OBS + REPORTED if not math.isfinite(r[o]))

    floor = {}
    try:
        e2 = series(f"{root}/e2/e2_DEVELOPMENT_series.csv")
        pre = [k for k in e2 if k in ref and e2[k]["t_s"] < 20e-3 - 1e-12]
        same = all(all((ref[k][o] == e2[k][o]) or (math.isnan(ref[k][o]) and math.isnan(e2[k][o]))
                       for o in ref[k]) for k in pre)
        print(f"e2 identical to e1 before t_L on {len(pre)} common samples: {same}")
        if max(r["t_s"] for r in e2.values()) >= H[1] - 1e-9 and tmax_ref >= H[1] - 1e-9:
            for o in OBS + REPORTED:
                a = [ref[k][o] for k in window(ref, H)]
                b = [e2[k][o] for k in window(e2, H)]
                ma, mb = sum(a) / len(a), sum(b) / len(b)
                floor[o] = abs(mb - ma) / abs(ma)
    except FileNotFoundError:
        print("e2: no series")

    for name in ["e1", "e2"] + names:
        try:
            budgets, walls, end, err = log_facts(f"{root}/{name}/{name}.log")
        except FileNotFoundError:
            continue
        rows = series(f"{root}/{name}/{name}_DEVELOPMENT_series.csv")
        worst = max((abs(b) for b in budgets), default=float("nan"))
        print(f"{name}: to {max(r['t_s'] for r in rows.values()) * 1e3:.2f} ms, largest |gas budget| {worst:.2e} "
              f"over {len(budgets)} lines, non-finite samples {nonfinite(rows)} (in L {nonfinite(rows, L)}, in H "
              f"{nonfinite(rows, H)}), real {end} s{', error: ' + err[-1] if err else ''}")

    _, ref_walls, _, _ = log_facts(f"{root}/e1/e1.log")
    for name in names:
        try:
            rows = series(f"{root}/{name}/{name}_DEVELOPMENT_series.csv")
        except FileNotFoundError:
            print(f"{name}: no series")
            continue
        _, walls, _, _ = log_facts(f"{root}/{name}/{name}.log")
        print(f"\n{name}, window L [0.5, 20] ms")
        common = [k for k in window(ref, L) if k in rows]
        span = window(ref, L)
        print(f"  {len(common)} of e1's {len(span)} samples in the window")
        ok = len(common) == len(span) and len(span) > 0
        for o in OBS + REPORTED:
            if not common:
                break
            mag = math.sqrt(sum(ref[k][o] ** 2 for k in common) / len(common))
            d = [rows[k][o] - ref[k][o] for k in common]
            rms = math.sqrt(sum(x * x for x in d) / len(d)) / mag
            big = max(abs(x) for x in d) / mag
            judged = o in OBS
            if judged and not (rms <= 5e-3 and big <= 2e-2):
                ok = False
            print(f"  {o:10s} RMS {rms * 100:8.4f}%  largest {big * 100:8.4f}%{'' if judged else '  (reported)'}")
        print(f"  window L: {'PASS' if ok else 'FAIL'}")
        if tmax_ref >= H[1] - 1e-9 and max(r["t_s"] for r in rows.values()) >= H[1] - 1e-9:
            okh = True
            print(f"{name}, window H [23, 30] ms (means)")
            for o in OBS + REPORTED:
                a = [ref[k][o] for k in window(ref, H)]
                b = [rows[k][o] for k in window(rows, H)]
                ma, mb = sum(a) / len(a), sum(b) / len(b)
                diff = (mb - ma) / abs(ma)
                limit = max(0.03 if o == "qwall_W" else 0.01, 2 * floor.get(o, 0))
                if o in OBS and abs(diff) > limit:
                    okh = False
                print(f"  {o:10s} {diff * 100:+8.4f}%  limit {limit * 100:.3f}%  seed floor "
                      f"{floor.get(o, float('nan')) * 100:.4f}%")
            print(f"  window H: {'PASS' if okh else 'FAIL'}")
        else:
            print(f"{name}, window H: not judged (e1 or {name} short of 30 ms)")
        for t0, t1 in [(0.5e-3, 20e-3), (23e-3, 30e-3)]:
            a0, a1 = wall_at(ref_walls, t0), wall_at(ref_walls, t1)
            b0, b1 = wall_at(walls, t0), wall_at(walls, t1)
            if None not in (a0, a1, b0, b1):
                print(f"  speed-up on [{t0 * 1e3:g}, {t1 * 1e3:g}] ms: {(a1 - a0) / (b1 - b0):.2f} "
                      f"(e1 {a1 - a0:.0f} s, {name} {b1 - b0:.0f} s)")


if __name__ == "__main__":
    main()
