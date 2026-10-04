#!/usr/bin/env python3
"""Headless plots for tests/chamber_study.cpp output (matplotlib Agg, PNG only).

Usage: chamber_plots.py <prefix> [title]
Reads <prefix>_history.csv, <prefix>_mesh.csv and <prefix>_field_*.csv; writes
<prefix>_history.png and <prefix>_frames.png (contact sheet: temperature above the axis,
Mach number below, one panel per snapshot).
"""
import csv
import glob
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def read(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def history(prefix, title):
    h = read(prefix + "_history.csv")
    t = h["t"] * 1e3
    fig, ax = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
    ax[0].plot(t, h["p_injector"] / 1e6, label="injector face")
    ax[0].plot(t, h["p_chamber_end"] / 1e6, label="chamber end (static)")
    ax[0].set_ylabel("pressure [MPa]")
    ax[0].legend()
    ax[1].plot(t, h["supply"], label="scheduled supply")
    ax[1].plot(t, h["inlet"], "--", label="delivered")
    ax[1].plot(t, h["outlet"], label="through exit")
    ax[1].set_ylabel("mass flow [kg/s]")
    ax[1].legend()
    ax[2].plot(t, h["F_vac"])
    ax[2].set_ylabel("vacuum thrust [N]")
    ax3 = ax[3]
    ok = h["supply"] > 0.05 * h["supply"].max()
    ax3.plot(t[ok], h["isp_vac_supply"][ok])
    ax3.set_ylabel("vacuum Isp on supply [s]")
    ax3.set_xlabel("time [ms]")
    for a in ax:
        a.grid(alpha=0.3)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(prefix + "_history.png", dpi=110)
    plt.close(fig)


def frames(prefix, title):
    mesh = read(prefix + "_mesh.csv")
    zs, rs = mesh["z"], mesh["radius"]
    files = sorted(glob.glob(prefix + "_field_*.csv"))
    if not files:
        return
    first = read(files[0])
    nz, nr = int(first["i"].max()) + 1, int(first["j"].max()) + 1
    # Cell corners: station z, wall radius times j/nr.
    zc = np.repeat(zs[:, None], nr + 1, axis=1)
    rc = rs[:, None] * (np.arange(nr + 1)[None, :] / nr)
    fields = [(f, read(f)) for f in files]
    tmax = max(d["T"].max() for _, d in fields)
    mmax = max(d["mach"].max() for _, d in fields)
    cols = 3
    rows = (len(fields) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.6 * cols + 1.6, 2.0 * rows + 0.6), squeeze=False,
                             layout="constrained")
    for ax, (f, d) in zip(axes.flat, fields):
        T = d["T"].reshape(nz, nr)
        M = d["mach"].reshape(nz, nr)
        pt = ax.pcolormesh(zc * 1e3, rc * 1e3, T, vmin=200, vmax=tmax, cmap="inferno", shading="flat")
        pm = ax.pcolormesh(zc * 1e3, -rc * 1e3, M, vmin=0, vmax=mmax, cmap="viridis", shading="flat")
        ax.plot(zs * 1e3, rs * 1e3, "k", lw=0.8)
        ax.plot(zs * 1e3, -rs * 1e3, "k", lw=0.8)
        us = int(f.rsplit("_", 1)[1].split(".")[0])
        ax.set_title(f"t = {us / 1000:.3f} ms", fontsize=9)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=7)
    for ax in list(axes.flat)[len(fields):]:
        ax.axis("off")
    fig.colorbar(pt, ax=axes, shrink=0.6, location="right", label="T [K] (upper half)")
    fig.colorbar(pm, ax=axes, shrink=0.6, location="right", label="Mach (lower half)")
    fig.suptitle(title + "  (z, r in mm)")
    fig.savefig(prefix + "_frames.png", dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    prefix = sys.argv[1]
    title = sys.argv[2] if len(sys.argv) > 2 else prefix
    history(prefix, title)
    frames(prefix, title)
