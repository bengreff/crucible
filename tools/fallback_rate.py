#!/usr/bin/env python3
"""The share of a run's cell updates that a Table B's domain sends to CVODES (docs/evidence/TABLE_B.md,
"The fallback rate in a C1-like run").

Usage: fallback_rate.py <state samples CSV> <node report CSV> [candidate]

The state samples are tests/chamber_study.cpp's CRUCIBLE_STATE_SAMPLES (each cell's f_H, Z_N, eta and rho
every 2 us, with the flow steps since the previous sample), with eta in the table's mechanism's energies
(crucible_state_samples_eta converts them if the run used another); the node report is the builder's
CRUCIBLE_TABLE_B_NODE_REPORT. A cell update counts as a fallback if its state is outside the axes or any
of its 16 corners failed a gate (2(a) turn-back above 1e-9 or an N-gate failure, 2(b) above 3 K, 2(c), or
the outside flag). The runtime's other causes (off the curve, Delta Y, frozen reach, AR) are not counted,
so the rate is a lower bound. Each sample is weighted by its steps, over the whole run and the last ms.
"""
import bisect
import csv
import itertools
import math
import sys
from collections import Counter


def main():
    samples_path, report_path = sys.argv[1], sys.argv[2]
    candidate = sys.argv[3] if len(sys.argv) > 3 else "0"
    failed = {}
    for row in csv.DictReader(open(report_path)):
        if row["candidate"] != candidate:
            continue
        key = tuple(float(row[c]) for c in ("f_H", "Z_N", "eta", "rho"))
        # 2(b) is blank where it was not evaluated (a zero-span branch), as the builder's gate reads it.
        dT = [float(row[c]) if row[c] else 0.0 for c in ("dT_2b_lower", "dT_2b_upper")]
        bad = (float(row["turn_back"]) > 1e-9 or row["n_gate_ok"] != row["n_gate_checked"] or max(dT) > 3
               or row["integrator_failure"] == "1" or row["outside"] == "1")
        failed[key] = bad
    axes = [sorted({k[a] for k in failed}) for a in range(4)]
    names = ("f_H", "Z_N", "eta", "rho")
    print(f"table: {len(failed)} nodes, {sum(failed.values())} failing a gate; axes", end="")
    for n, a in zip(names, axes):
        print(f"  {n} [{a[0]:.4g}, {a[-1]:.4g}] ({len(a)})", end="")
    print()

    def bracket(a, x):
        if not (a[0] <= x <= a[-1]):
            return None
        i = min(max(bisect.bisect_right(a, x) - 1, 0), len(a) - 2)
        return (a[i], a[i + 1])

    # The node report prints rounded coordinates; the cell is compared with them as printed.
    def classify(f, zn, eta, rho):
        if math.isnan(f):
            return "no H or O (pure N2)"
        brackets = [None] * 4
        for a in (1, 0, 2, 3):  # Z_N first: a cell of mostly ambient N2 is outside for that reason
            x = (f, zn, eta, rho)[a]
            brackets[a] = bracket(axes[a], x)
            if brackets[a] is None:
                return f"{names[a]} {'above' if x > axes[a][-1] else 'below'} the axis"
        if any(failed[c] for c in itertools.product(*brackets)):
            return "a corner failed a gate"
        return "table"

    end = 0.0
    rows = []
    for row in csv.DictReader(open(samples_path)):
        t = float(row["t"])
        end = max(end, t)
        f = float(row["f_H"]) if row["f_H"] not in ("nan", "-nan") else math.nan
        rows.append((t, int(row["steps"]), classify(f, float(row["Z_N"]), float(row["eta"]), float(row["rho"]))))
    for label, t0 in (("whole run", -1.0), ("last ms", end - 1e-3)):
        weight = Counter()
        for t, steps, cause in rows:
            if t > t0:
                weight[cause] += steps
        total = sum(weight.values())
        fallback = total - weight["table"]
        print(f"{label}: {total} cell updates, fallback {fallback} = {100 * fallback / total:.3f}% "
              f"({'under' if fallback < 0.01 * total else 'over'} 1%)")
        for cause, w in weight.most_common():
            if cause != "table":
                print(f"  {cause}: {100 * w / total:.3f}%")


if __name__ == "__main__":
    main()
