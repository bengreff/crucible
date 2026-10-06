#!/usr/bin/env python3
"""kappa and B of the SST model's own log layer, from a resolved profile of tools/sst_pipe_1d.py.

The wall functions' constants (docs/evidence/WALL_FUNCTIONS.md, design item 7): a least-squares line of
u+ against ln y+ over 50 < y+ < 0.1 Re_tau, in wall units of the wall (tau_w = f R / 2 from the force
balance; rho_w from the first cell's p + 2/3 rho k at T_w; mu_w from tools/n2_properties.csv). Also
prints the fit's sensitivity to its range, the local slope, and Spalding's formula with the fitted
constants and with Nichols' 0.4 and 5.5 against the profile.

Usage: log_law_fit.py --profile p.csv --radius R --force F --twall T
"""
import argparse
import math
import os

import numpy as np


def spalding(u, kappa, b):
    x = kappa * u
    return u + math.exp(-kappa * b) * (math.exp(x) - 1 - x - x * x / 2 - x ** 3 / 6)


def spalding_u(y, kappa, b):
    lo, hi = 0.0, 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if spalding(mid, kappa, b) < y else (lo, mid)
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("--radius", type=float, required=True)
    ap.add_argument("--force", type=float, required=True)
    ap.add_argument("--twall", type=float, required=True)
    a = ap.parse_args()
    props = np.genfromtxt(os.path.join(os.path.dirname(__file__), "n2_properties.csv"), delimiter=",", skip_header=2)
    mu_w, gas_r = float(np.interp(a.twall, props[:, 0], props[:, 1])), props[0, 4]
    r, u, k, w, t, rho, _ = np.genfromtxt(a.profile, delimiter=",", skip_header=1).T
    rho_w = (rho[-1] * gas_r * t[-1] + 2 / 3 * rho[-1] * k[-1]) / (gas_r * a.twall)
    u_tau = math.sqrt(a.force * a.radius / 2 / rho_w)
    nu = mu_w / rho_w
    order = np.argsort(a.radius - r)
    yp, up = ((a.radius - r) * u_tau / nu)[order], (u / u_tau)[order]
    re_tau = a.radius * u_tau / nu
    print(f"Re_tau {re_tau:.2f}, u_tau {u_tau:.6f} m/s, rho_w {rho_w:.6f} kg/m^3, mu_w {mu_w:.6e} Pa s")
    fits = {}
    for lo, hi in [(50, 0.1), (30, 0.1), (100, 0.1), (50, 0.05), (50, 0.2)]:
        m = (yp > lo) & (yp < hi * re_tau)
        slope, b = np.polyfit(np.log(yp[m]), up[m], 1)
        resid = np.abs(up[m] - (slope * np.log(yp[m]) + b)).max()
        fits[(lo, hi)] = (1 / slope, b)
        print(f"fit over {lo} < y+ < {hi} Re_tau ({m.sum()} points): kappa {1 / slope:.4f}  B {b:.4f}  "
              f"largest |u+ - line| {resid:.4f}" + ("   <- the stated rule" if (lo, hi) == (50, 0.1) else ""))
    kappa, b = fits[(50, 0.1)]
    slope = np.gradient(up, np.log(yp))
    print(f"{'y+':>8} {'model u+':>9} {'local kappa':>12} {'Spalding fitted':>18} {'Spalding 0.4/5.5':>18}")
    for target in [1, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 5000]:
        i = int(np.argmin(np.abs(yp - target)))
        if yp[i] > 0.5 * re_tau:
            break
        s1, s2 = spalding_u(yp[i], kappa, b), spalding_u(yp[i], 0.4, 5.5)
        print(f"{yp[i]:8.1f} {up[i]:9.3f} {1 / slope[i]:12.4f} {s1:9.3f} ({100 * (s1 / up[i] - 1):+5.1f}%) "
              f"{s2:9.3f} ({100 * (s2 / up[i] - 1):+5.1f}%)")


if __name__ == "__main__":
    main()
