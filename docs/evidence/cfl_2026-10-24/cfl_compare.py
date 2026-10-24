#!/usr/bin/env python3
"""Lever B's comparison (SPEED.md section 5c), run after mac_cfl.sh: python3 cfl_compare.py [/tmp/cfl].

For each pair, the difference of each quantity at each common 10 us row over the window, divided by the
reference run's mean |Q| over the window: its largest value, its RMS and its value at the window's end.
Criterion (a) from each run's log and done-file; criterion (b) from the table (F_vac, the outlet flow, p0)."""
import csv, math, re, sys

base = sys.argv[1] if len(sys.argv) > 1 else "/tmp/cfl"
cols = [("outlet flow", "mdot_kg_s"), ("F_vac", "fvac_N"), ("p0", "p0_Pa"), ("face pressure", "pface_Pa"),
        ("hot-gas wall heat", "qwall_W")]
judged = ("outlet flow", "F_vac", "p0")


def series(name):
    with open(f"{base}/{name}/{name}_DEVELOPMENT_series.csv") as f:
        rows = [r for r in f if not r.startswith("#")]
    return {round(float(r["t_s"]) * 1e5): r for r in csv.DictReader(rows)}


def run(name):
    log = open(f"{base}/{name}/{name}.log", errors="replace").read()
    done = open(f"{base}/{name}/{name}_DEVELOPMENT_done").read().strip()
    real = re.findall(r"^real ([\d.]+)", log, re.M)
    steps = re.findall(r"^  (\d+) base steps \((\d+) redone\)", log, re.M)
    end = re.search(r"^end \[.*?\]: .*?(\d+) non-finite values;.*?worst budgets \(a\) (\S+), \(b\) (\S+), \(c\) (\S+);", log, re.M)
    out = {"done": done, "wall": float(real[-1]) if real else math.nan,
           "steps": int(steps[-1][0]) if steps else -1, "redone": int(steps[-1][1]) if steps else -1}
    if end:
        out["nonfinite"] = int(end.group(1))
        out["budget"] = max(float(end.group(i).rstrip(",;")) for i in (2, 3, 4))
    return out


def compare(test, ref, lo, hi):
    a, b = series(test), series(ref)
    keys = sorted(k for k in b if k in a and lo * 1e5 - 0.5 <= k <= hi * 1e5 + 0.5)
    out = {}
    for label, c in cols:
        pairs = [(float(a[k][c]), float(b[k][c])) for k in keys]
        pairs = [(x, y) for x, y in pairs if math.isfinite(x) and math.isfinite(y)]
        scale = sum(abs(y) for _, y in pairs) / len(pairs)
        d = [(x - y) / scale for x, y in pairs]
        out[label] = (max(abs(v) for v in d), math.sqrt(sum(v * v for v in d) / len(d)), d[-1], len(d), len(keys) - len(d))
    return out


def identical_before(test, ref, t):
    a, b = series(test), series(ref)
    keys = [k for k in b if k in a and k < t * 1e5 - 0.5]
    same = all(a[k] == b[k] for k in keys)
    return same, len(keys)


def show(title, res):
    print(title)
    for label, (big, rms, last, n, dropped) in res.items():
        print(f"  {label:18s} largest {100 * big:.4f}%, RMS {100 * rms:.4f}%, at the end {100 * last:+.4f}%"
              f" ({n} rows{', %d undefined dropped' % dropped if dropped else ''})")


names = ["L4", "L8", "H4", "H8", "H4s"]
runs = {}
for n in names:
    try:
        runs[n] = run(n)
    except FileNotFoundError:
        pass
print("Runs:")
for n, r in runs.items():
    print(f"  {n:4s} {r['done']}; wall {r['wall']:.0f} s; {r['steps']} base steps ({r['redone']} redone);"
          f" non-finite {r.get('nonfinite', '?')}; worst budget {r.get('budget', math.nan):.2e}")
for n in ("L8", "H8"):
    if n in runs:
        r = runs[n]
        ok = r["done"] == "exit 0" and r.get("nonfinite", 1) == 0 and r.get("budget", 1) <= 1e-11
        print(f"(a) {n}: {'passes' if ok else 'FAILS'}")
if "L4" in runs and "L8" in runs:
    res = compare("L8", "L4", 0.1e-3, 1.0e-3)
    show("L8 against L4 over [0.1, 1.0] ms, over L4's mean |Q|:", res)
    worst = max(res[q][1] for q in judged)
    print(f"(b) L: the largest RMS of F_vac, the outlet flow and p0 is {100 * worst:.4f}%: {'passes' if worst <= 2e-3 else 'FAILS'} (bar 0.2%)")
if all(n in runs for n in ("H4", "H8", "H4s")):
    same, n = identical_before("H4s", "H4", 1.0e-3)
    print(f"H4s against H4 before LOX opens (t < 1 ms, {n} rows): {'bit for bit' if same else 'DIFFERENT'}")
    floor = compare("H4s", "H4", 1.0e-3, 2.0e-3)
    show("H4s against H4 over [1.0, 2.0] ms (the floor):", floor)
    res = compare("H8", "H4", 1.0e-3, 2.0e-3)
    show("H8 against H4 over [1.0, 2.0] ms:", res)
    for q in judged:
        bar = max(2e-3, floor[q][1])
        print(f"(b) H, {q}: RMS {100 * res[q][1]:.4f}% against the bar {100 * bar:.4f}%: {'passes' if res[q][1] <= bar else 'FAILS'}")
for a, b in (("L8", "L4"), ("H8", "H4")):
    if a in runs and b in runs and runs[a]["wall"] > 0:
        print(f"{b} / {a}: wall {runs[b]['wall'] / runs[a]['wall']:.2f}, base steps {runs[b]['steps'] / max(runs[a]['steps'], 1):.2f}")
