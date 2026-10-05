#!/usr/bin/env python3
"""Headless plots for tests/flame_study.cpp output (matplotlib Agg, PNG only).

Usage: flame_plots.py <output png prefix> <S_L m/s> <run prefix> [<run prefix> ...]
Each run prefix has <prefix>_history.csv and <prefix>_profile.csv; the label is the text after
the last '/'. Writes <output>_history.png (consumption speed, flame position, pressure range,
gas velocity ahead of the flame) and <output>_profiles.png (T, Y_H2, Y_H, Y_OH against distance
from the flame, the engine's last state for every run and the free flame).
"""
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return rows


def columns(rows, keys):
    return {k: np.array([float(r[k]) for r in rows]) for k in keys}


def main():
    out, speed, runs = sys.argv[1], float(sys.argv[2]), sys.argv[3:]
    fig, ax = plt.subplots(4, 1, figsize=(8, 11), sharex=True)
    for run in runs:
        label = run.rsplit("/", 1)[-1]
        h = columns(read(run + "_history.csv"), ["t", "x_f", "S_c", "p_min", "p_max", "u_ahead"])
        t = h["t"] * 1e3
        ax[0].plot(t, h["S_c"], label=label)
        ax[1].plot(t, h["x_f"] * 1e3, label=label)
        ax[2].plot(t, (h["p_max"] - h["p_min"]) / 1e3, label=label)
        ax[3].plot(t, h["u_ahead"], label=label)
    ax[0].axhline(speed, color="k", ls="--", lw=1, label="free flame S_L")
    ax[0].set_ylim(0.8 * speed, 1.2 * speed)
    ax[0].set_ylabel("S_c [m/s]")
    ax[1].set_ylabel("flame position [mm]")
    ax[2].set_ylabel("p_max - p_min [kPa]")
    ax[3].set_ylabel("u_z 0.5 mm ahead [m/s]")
    ax[3].set_xlabel("t [ms]")
    for a in ax:
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out + "_history.png", dpi=110)

    fields = [("T", "T [K]"), ("Y_H2", "Y_H2"), ("Y_H", "Y_H"), ("Y_OH", "Y_OH")]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    for i, run in enumerate(runs):
        rows = read(run + "_profile.csv")
        label = run.rsplit("/", 1)[-1]
        for source in ("engine", "freeflame"):
            if source == "freeflame" and i > 0:
                continue
            sel = [r for r in rows if r["source"] == source]
            c = columns(sel, ["z_minus_xf"] + [f for f, _ in fields])
            z = c["z_minus_xf"] * 1e3
            for a, (f, name) in zip(ax.flat, fields):
                if source == "engine":
                    a.plot(z, c[f], ".-", ms=3, lw=0.8, label=label)
                else:
                    a.plot(z, c[f], "k-", lw=1.2, label="free flame")
                a.set_ylabel(name)
    for a in ax.flat:
        a.set_xlim(-1.0, 1.5)
        a.set_xlabel("z - x_f [mm]")
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out + "_profiles.png", dpi=110)


if __name__ == "__main__":
    main()
