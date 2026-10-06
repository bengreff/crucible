#!/usr/bin/env python3
"""Compare two crucible_chamber_study histories at the printed times (every 0.5 ms from 0.5 ms).

Usage: compare_histories.py <reference_history.csv> <test_history.csv> [tolerance=1e-3]
Prints the relative difference (test / reference - 1) of injector pressure, outlet mass flow and
vacuum thrust at each printed time, and the largest over all of them (TABLE_A criterion 3).
"""
import csv
import sys


def rows(path):
    with open(path) as f:
        return {round(float(r["t"]) * 1e6): r for r in csv.DictReader(f)}


def main():
    ref, test = rows(sys.argv[1]), rows(sys.argv[2])
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
    keys = [("p_injector", "p_inj"), ("outlet", "outlet"), ("F_vac", "F_vac")]
    worst = {k: 0.0 for k, _ in keys}
    print("t [ms]  " + "  ".join(f"{label:>10}" for _, label in keys))
    for us in range(500, 8001, 500):
        if us not in ref or us not in test:
            continue
        d = {k: float(test[us][k]) / float(ref[us][k]) - 1 for k, _ in keys}
        for k in d:
            worst[k] = max(worst[k], abs(d[k]))
        print(f"{us / 1000:6.3f}  " + "  ".join(f"{d[k]:+10.3e}" for k, _ in keys))
    print("max     " + "  ".join(f"{worst[k]:10.3e}" for k, _ in keys))
    ok = worst["p_injector"] <= tol and worst["outlet"] <= tol
    print(f"injector pressure and outlet mass flow within {tol:g} at every printed time: {'pass' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
