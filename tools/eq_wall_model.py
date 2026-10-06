#!/usr/bin/env python3
"""A priori test of the equilibrium ODE wall model against the SST pipe references (tools/sst_pipe_1d.py).

The model (Kawai and Larsson, Phys. Fluids 24, 015105, 2012; Larsson et al., Mech. Eng. Rev. 3, 2016):
in a layer of constant shear stress and constant total heat flux between the wall (u = 0, T = T_w) and
the first cell's centroid (u_1, T_1 at distance y_1), with p uniform (p_1) and rho = p / (R T),
  (mu + mu_t) du/dy = tau_w
  tau_w u + c_p (mu / Pr + mu_t / Pr_t) dT/dy = q_w          (q_w: heat flux from the gas into the wall)
  mu_t = kappa rho y sqrt(tau_w / rho) [1 - exp(-y* / A+)]^2,  y* = y sqrt(rho tau_w) / mu  (semi-local)
mu, lambda and c_p from tools/n2_properties.csv (Cantera, the reference's own properties). Given
(u_1, T_1, y_1) the pair (tau_w, q_w) is found by shooting: RK4 from the wall on log-spaced points and a
2-D Newton on u(y_1) = u_1, T(y_1) = T_1.

kappa is the SST model's own (0.3697, WALL_FUNCTIONS.md design item 7); A+ is fitted here so the cold
profile's a priori shear error is smallest over y+ 5 to 1,000; Pr_t 0.9 is the reference's.

Usage: eq_wall_model.py [--reference sst|hp] [--aplus A] [--fit] | --transforms
"""
import argparse
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EVIDENCE = os.path.join(HERE, "..", "docs", "evidence", "wall_functions")
P = np.genfromtxt(os.path.join(HERE, "n2_properties.csv"), delimiter=",", skip_header=2)
GAS_R = P[0, 4]
RADIUS = 5e-3
CASES = {
    "cold": dict(profile="reference_re_tau_1e4_profile_6400.csv", twall=300.0, force=45690.0, heat=0.0),
    "heated": dict(profile="reference_heated_profile_6400.csv", twall=600.0, force=19800.0, heat=1.95e9),
    # SSTs with the Hasan, Elias, Menter and Pecnik corrections (sst_pipe_1d.py --correction hp), the same cases.
    "cold_hp": dict(profile="reference_hp_cold_profile_6400.csv", twall=300.0, force=45690.0, heat=0.0),
    "heated_hp": dict(profile="reference_hp_heated_profile_6400.csv", twall=600.0, force=19800.0, heat=1.95e9),
}


def props(t):
    return np.interp(t, P[:, 0], P[:, 1]), np.interp(t, P[:, 0], P[:, 2]), np.interp(t, P[:, 0], P[:, 3])


def shoot(tau, q, y1, tw, p, kappa, aplus, prt, n=400):
    """Integrate from the wall to y1; returns u(y1), T(y1)."""
    y = np.concatenate([[0.0], y1 * np.geomspace(1e-6, 1.0, n)])

    def rhs(yy, u, t):
        mu, lam, cp = props(t)
        rho = p / (GAS_R * t)
        ystar = yy * math.sqrt(rho * tau) / mu
        mut = kappa * yy * math.sqrt(rho * tau) * (1 - math.exp(-ystar / aplus)) ** 2
        return tau / (mu + mut), (q - tau * u) / (cp * (lam / cp + mut / prt))

    u, t = 0.0, tw
    for a, b in zip(y[:-1], y[1:]):
        h = b - a
        k1 = rhs(a, u, t)
        k2 = rhs(a + h / 2, u + h / 2 * k1[0], t + h / 2 * k1[1])
        k3 = rhs(a + h / 2, u + h / 2 * k2[0], t + h / 2 * k2[1])
        k4 = rhs(b, u + h * k3[0], t + h * k3[1])
        u += h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        t += h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    return u, t


def solve(u1, t1, y1, tw, p, kappa, aplus, prt, tau0, q0):
    """Newton on (ln tau_w, q_w) with a finite-difference Jacobian."""
    x = np.array([math.log(tau0), q0])
    scale = np.array([u1, max(abs(t1 - tw), 1e-3 * tw)])
    for _ in range(60):
        f = (np.array(shoot(math.exp(x[0]), x[1], y1, tw, p, kappa, aplus, prt)) - [u1, t1]) / scale
        if np.max(np.abs(f)) < 1e-10:
            break
        jac = np.empty((2, 2))
        steps = [1e-6, 1e-6 * max(abs(x[1]), 1.0)]
        for j in range(2):
            xp = x.copy()
            xp[j] += steps[j]
            fp = (np.array(shoot(math.exp(xp[0]), xp[1], y1, tw, p, kappa, aplus, prt)) - [u1, t1]) / scale
            jac[:, j] = (fp - f) / steps[j]
        dx = np.linalg.solve(jac, -f)
        dx[0] = max(-0.5, min(0.5, dx[0]))
        x += dx
    else:
        raise RuntimeError("no convergence")
    return math.exp(x[0]), x[1]


def references(case):
    c = CASES[case]
    r, u, k, w, t, rho, _ = np.genfromtxt(os.path.join(EVIDENCE, c["profile"]), delimiter=",", skip_header=1).T
    order = np.argsort(r)
    r, u, k, t, rho = r[order], u[order], k[order], t[order], rho[order]
    edges = np.concatenate([[0.0], 0.5 * (r[1:] + r[:-1]), [RADIUS]])
    flow = np.sum(u * math.pi * np.diff(edges ** 2))
    tau = c["force"] * RADIUS / 2
    q = (c["heat"] * math.pi * RADIUS ** 2 + c["force"] * flow) / (2 * math.pi * RADIUS)
    pw = rho[-1] * GAS_R * t[-1] + 2 / 3 * rho[-1] * k[-1]
    return r, u, t, rho, tau, q, pw, c["twall"]


def apriori(case, kappa, aplus, prt, targets=(1, 5, 11, 30, 100, 300, 1000, 3000), quiet=False):
    r, u, t, rho, tau, q, pw, tw = references(case)
    mu_w = props(tw)[0]
    rho_w = pw / (GAS_R * tw)
    u_tau = math.sqrt(tau / rho_w)
    yplus = (RADIUS - r) * rho_w * u_tau / mu_w
    out = []
    if not quiet:
        print(f"{case}: tau_w {tau:.6f} Pa, q_w {q:.6e} W/m^2, Re_tau {RADIUS * rho_w * u_tau / mu_w:.1f}; kappa {kappa}, A+ {aplus}, Pr_t {prt}")
        print(f"{'y+':>9} {'T1/Tw':>8} {'tau model-1':>12} {'q model-1':>12}")
    for target in targets:
        i = int(np.argmin(np.abs(yplus - target)))
        y1 = RADIUS - r[i]
        p1 = rho[i] * GAS_R * t[i]
        ts, qs = solve(u[i], t[i], y1, tw, p1, kappa, aplus, prt, tau, q if q != 0 else 1.0)
        out.append((yplus[i], ts / tau - 1, qs / q - 1))
        if not quiet:
            print(f"{yplus[i]:9.2f} {t[i] / tw:8.4f} {ts / tau - 1:+12.4e} {qs / q - 1:+12.4e}")
    return out


def transforms(cold="cold", heated="heated"):
    """The heated profile under van Driest (y+, wall units) and Trettel-Larsson (y*, semi-local) scaling
    against the cold profile's u+(y+), with the local slopes 1 / (du / d ln y)."""
    def one(case):
        r, u, t, rho, tau, q, pw, tw = references(case)
        mu, mu_w = props(t)[0], props(tw)[0]
        rho_w = pw / (GAS_R * tw)
        ut = math.sqrt(tau / rho_w)
        y = RADIUS - r
        o = np.argsort(y)
        y, u, rho, mu = y[o], u[o], rho[o], mu[o]
        yp, up, ystar = y * rho_w * ut / mu_w, u / ut, y * np.sqrt(rho * tau) / mu
        s = np.sqrt(rho / rho_w)
        g = s * (1 + 0.5 * y / rho * np.gradient(rho, y) - y / mu * np.gradient(mu, y))
        integral = lambda f: np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(up))]) + f[0] * up[0]
        return yp, up, integral(s), ystar, integral(g)
    yc, uc, _, _, _ = one(cold)
    yh, uh, uvd, ys, utl = one(heated)
    slope = lambda yy, uu: 1 / np.gradient(uu, np.log(yy))
    kc, kv, kt = slope(yc, uc), slope(yh, uvd), slope(ys, utl)
    print(f"{heated} profile (T_w 600 K, Re_tau {yh.max():.0f}) under compressible scalings, against the {cold} profile")
    print(f"{'y+ or y*':>9} {'cold u+':>8} {'heated u+':>10} {'u_vD+(y+)':>10} {'u_TL+(y*)':>10} {'vD-cold':>8} {'TL-cold':>8}"
          f" {'kappa cold':>10} {'kappa vD':>9} {'kappa TL':>9}")
    for target in (1, 5, 11, 30, 100, 300, 1000):
        i, j, k = (int(np.argmin(np.abs(a - target))) for a in (yc, yh, ys))
        print(f"{target:9d} {uc[i]:8.3f} {uh[j]:10.3f} {uvd[j]:10.3f} {utl[k]:10.3f} {uvd[j] - uc[i]:+8.3f} {utl[k] - uc[i]:+8.3f}"
              f" {kc[i]:10.3f} {kv[j]:9.3f} {kt[k]:9.3f}")
    print(f"largest y* in the heated pipe: {ys.max():.0f} (the axis); y* 1,000 and above do not exist there")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kappa", type=float, default=0.3697)
    ap.add_argument("--aplus", type=float, default=17.0)
    ap.add_argument("--prt", type=float, default=0.9)
    ap.add_argument("--fit", action="store_true", help="fit A+ to the cold profile first")
    ap.add_argument("--transforms", action="store_true", help="only the scaling diagnostic of the heated profile")
    ap.add_argument("--reference", choices=("sst", "hp"), default="sst", help="standard SST or the corrected SST")
    a = ap.parse_args()
    cold, heated = ("cold", "heated") if a.reference == "sst" else ("cold_hp", "heated_hp")
    if a.transforms:
        transforms(cold, heated)
        return
    aplus = a.aplus
    if a.fit:
        best = None
        for trial in np.arange(10.0, 30.01, 0.5):
            errs = [abs(e[1]) for e in apriori(cold, a.kappa, trial, a.prt, targets=(5, 11, 30, 100, 300, 1000), quiet=True)]
            rms = math.sqrt(sum(e * e for e in errs) / len(errs))
            if best is None or rms < best[0]:
                best = (rms, trial)
        aplus = best[1]
        print(f"A+ fitted on the cold profile (y+ 5 to 1,000, rms shear error {best[0]:.3e}): {aplus}")
    for case in (cold, heated):
        apriori(case, a.kappa, aplus, a.prt)


if __name__ == "__main__":
    main()
