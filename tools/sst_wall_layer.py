#!/usr/bin/env python3
"""The corrected SST's own inner layer as a wall model, and its a priori check against the pipe references.

The fallback of docs/evidence/WALL_FUNCTIONS.md (stated 6 October 2026, 02:21 CDT, before this code): the
1-D constant-stress, constant-heat-flux layer with the k and omega equations of tools/sst_pipe_1d.py (SSTs
with the Hasan, Elias, Menter and Pecnik corrections in full, D^ic), its closure, discretisation and wall
rule reused unchanged, on a planar layer instead of a pipe:
  from the wall (u = 0, T = T_w, k = 0, omega = 10 * 6 nu_w / (beta1 d1^2) on the layer's own first cell)
  to a top at y_top = c y1, p + 2/3 rho k uniform; no body force and no heating, so the shear flux is tau_w
  at every face and the total energy flux (conduction, turbulent transport, shear work and the corrected k
  diffusion, as in the reference's energy equation) is constant apart from int Phi_k. At the top the shear
  flux is tau_w and the total energy flux is prescribed; k and omega have zero flux there.
Geometric grid in y with y1 a face, the first face at y+ 0.05 or below. Given (u1, T1, y1, p, T_w), Newton on
(ln tau_w, the top energy flux) until u and T at y1 (linear between the two centroids beside the face) match;
q_w is the wall conduction.

The a priori check gives the layer the reference profile's (u, T) at y+ 1, 5, 30, 100, 300, 1,000 and 3,000
(pressure: the reference's uniform p + 2/3 rho k) and compares tau_w and q_w with the reference's own (the
force and energy balances). Checks of the layer itself: c = 20 against c = 10, and the cells doubled.

Usage: sst_wall_layer.py [--case heated_hp|cold_hp] [--targets 1,5,30,...] [--checks]
"""
import argparse
import math
import sys

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from sst_pipe_1d import BETA, BETA_STAR, Pipe
from eq_wall_model import CASES, GAS_R, RADIUS, props, references


class Layer(Pipe):
    """A planar constant-stress layer, cells ordered as the pipe's: the top (j = 0) to the wall (j = n - 1)."""

    def __init__(self, y1, top, p, t_wall, first, n_inner=200, n_outer=60, pr_t=0.9, omega_factor=10.0,
                 correction="hp", dic=True):
        n = n_inner + n_outer
        super().__init__(top, p, t_wall, 0.0, n, 2.0, pr_t, omega_factor, 0.0, correction=correction, dic=dic,
                         planar=True, constant_pressure=True)
        yf = np.concatenate([[0.0], np.geomspace(first, y1, n_inner), np.geomspace(y1, top, n_outer + 1)[1:]])
        self.yf = yf
        yc = 0.5 * (yf[1:] + yf[:-1])
        # The pipe's coordinate r runs from the top (0) to the wall (top); distances are taken from y directly.
        self.rf = (top - yf)[::-1]
        self.rf[0], self.rf[-1] = 0.0, top
        self.wf = np.ones(n + 1)
        self.vol = np.diff(yf)[::-1].copy()
        self.d = yc[::-1].copy()
        self.rc = top - self.d
        self.dc = -np.diff(self.d)
        self.dw = self.d[-1]
        self.y1, self.top, self.n_inner = y1, top, n_inner
        self.tau, self.flux = 0.0, 0.0
        self.layer_force, self.layer_heat = 0.0, 0.0      # the diagnostic's sources (--pipe-gradients)

    def grad(self, phi, wall):
        face = np.empty(self.n + 1)
        face[1:-1] = np.diff(phi) / self.dc
        face[-1] = (wall - phi[-1]) / self.dw
        face[0] = face[1]
        return 0.5 * (face[:-1] + face[1:])

    def residual(self, x):
        u, k, w, t, pu = self.unpack(x)
        s = self.closure(x)
        rho, mu, mu_t, vol = s["rho"], s["mu"], s["mu_t"], self.vol
        net = lambda flux: flux[1:] - flux[:-1]
        # Face fluxes are positive toward the top (the pipe's inward); the top face's are prescribed.
        fu = self.face_flux(self.conductance(mu + mu_t, self.mu_wall), u, 0.0)
        fu[0] = -(self.tau - self.layer_force * self.top)
        diff_k, diff_w, fk = self.diffusion(s, k, w)
        ft = self.face_flux(self.conductance(s["lam"] + s["cp"] * mu_t / self.pr_t, self.lam_wall), t, self.t_wall)
        uf = np.zeros(self.n + 1)
        uf[1:-1] = 0.5 * (u[:-1] + u[1:])
        uf[0] = u[0] + (u[0] - u[1]) * (self.top - self.d[0]) / self.dc[0]
        ft[0] = -self.flux - uf[0] * fu[0]          # total energy flux toward the wall at the top: self.flux
        ru = net(fu) + self.layer_force * vol
        rk = diff_k + (s["prod"] - BETA_STAR * rho * w * k) * vol
        rw = diff_w + (s["omega_prod"] - s["beta"] * rho * w ** 2 + s["cross"] / w) * vol
        rt = net(ft) + net(uf * fu) + (diff_k if self.hp else net(fk)) + (self.layer_force * u + self.layer_heat) * vol
        return np.concatenate([ru, rk, rw, rt, [pu - self.p0]])

    def newton(self, x, tol=1e-12, max_steps=400, log=False):
        n = self.n
        u_tau = math.sqrt(self.tau / (self.p0 / (self.gas_r * self.t_wall)))
        dtau = 1e-3 * self.top / u_tau
        for step in range(max_steps):
            r0 = self.residual(x)
            jac = self.jacobian(x, r0)
            s = self.closure(x)
            m = np.concatenate([s["rho"] * self.vol] * 3 + [s["rho"] * s["cp"] * self.vol, [0.0]])
            dx = spsolve((diags(m / dtau) - jac).tocsc(), r0)
            alpha = 1.0
            for v in (1, 2):
                xv, dv = x[v * n:(v + 1) * n], dx[v * n:(v + 1) * n]
                neg = dv < 0
                if np.any(neg):
                    alpha = min(alpha, float(np.min(0.5 * xv[neg] / -dv[neg])))
            x = x + alpha * dx
            change = max(np.max(np.abs(dx[:n])) / np.max(np.abs(x[:n])), np.max(np.abs(dx[n:2 * n])) / np.max(x[n:2 * n]),
                         np.max(np.abs(dx[2 * n:3 * n]) / x[2 * n:3 * n]), np.max(np.abs(dx[3 * n:4 * n])) / self.t_wall)
            dtau = dtau * 3 if alpha == 1.0 else dtau * 0.5
            if log:
                print(f"    step {step:3d} alpha {alpha:.3f} dtau {dtau:.2e} update {change:.3e}", file=sys.stderr)
            if alpha == 1.0 and change < tol and dtau * u_tau / self.top > 1e6:
                return x
        raise RuntimeError(f"layer Newton did not converge (last update {change:.3e})")

    def start(self, tau, q):
        """The Kawai-Larsson mixing-length layer (kappa 0.41, A+ 17) for u and T; k and omega from it."""
        y = self.d[::-1]
        u, t = np.empty(self.n), np.empty(self.n)
        uu, tt, ya = 0.0, self.t_wall, 0.0
        for j, yb in enumerate(y):
            for a, b in ((ya + (yb - ya) * i / 8, ya + (yb - ya) * (i + 1) / 8) for i in range(8)):
                mu, lam, cp = props(tt)
                rho = self.p0 / (GAS_R * tt)
                ym = 0.5 * (a + b)
                mut = 0.41 * ym * math.sqrt(rho * tau) * (1 - math.exp(-ym * math.sqrt(rho * tau) / mu / 17.0)) ** 2
                uu += (b - a) * tau / (mu + mut)
                tt += (b - a) * (q - tau * uu) / (cp * (lam / cp + mut / self.pr_t))
            u[j], t[j], ya = uu, tt, yb
        rho = self.p0 / (GAS_R * t)
        mu = props(t)[0]
        ystar = y * np.sqrt(rho * tau) / mu
        k = tau / (rho * math.sqrt(BETA_STAR)) * np.minimum(1.0, (ystar / 10.0) ** 2)
        w = np.maximum(np.sqrt(tau / rho) / (math.sqrt(BETA_STAR) * 0.41 * y), 6 * mu / rho / (BETA[0] * y ** 2))
        return np.concatenate([u[::-1], k[::-1], w[::-1], t[::-1], [self.p0]])

    def at_y1(self, x):
        u, k, w, t, pu = self.unpack(x)
        below = self.n - self.n_inner          # the first cell under y1 (toward the wall), in the pipe's order
        above = below - 1
        f = (self.y1 - self.d[below]) / (self.d[above] - self.d[below])
        return u[below] + f * (u[above] - u[below]), t[below] + f * (t[above] - t[below])

    def wall(self, x):
        u, k, w, t, pu = self.unpack(x)
        return self.mu_wall * u[-1] / self.dw, self.lam_wall * (t[-1] - self.t_wall) / self.dw


def model(u1, t1, y1, p, tw, tau0, q0, c=10.0, cells=(200, 60), correction="hp", log=False, sources=(0.0, 0.0)):
    """(tau_w, q_w) of the layer through (u1, T1) at y1."""
    rho_w = p / (GAS_R * tw)
    unit = props(tw)[0] / (rho_w * math.sqrt(tau0 / rho_w))      # one wall unit at the starting tau_w
    # The top at c max(y1, y+ 30): below y+ 30 a top at c y1 lies in the buffer layer, where k's flux is not zero
    # (amendment of 6 October, about 02:30 CDT).
    layer = Layer(y1, c * max(y1, 30 * unit), p, tw, min(0.05 * unit, 0.01 * y1), cells[0], cells[1], correction=correction)
    layer.layer_force, layer.layer_heat = sources
    q_scale = max(abs(q0), 1e-4 * rho_w * float(props(tw)[2]) * math.sqrt(tau0 / rho_w) * tw)
    t_scale = max(abs(t1 - tw), 1e-3 * tw)
    state = {"x": layer.start(tau0, q0)}

    def miss(v):
        layer.tau, layer.flux = math.exp(v[0]), v[1] * q_scale
        state["x"] = layer.newton(state["x"], log=log)
        ua, ta = layer.at_y1(state["x"])
        return np.array([(ua - u1) / u1, (ta - t1) / t_scale])

    v = np.array([math.log(tau0), q0 / q_scale])
    for it in range(40):
        f = miss(v)
        base = state["x"]
        if log:
            print(f"  inversion {it}: tau {math.exp(v[0]):.6e} flux {v[1] * q_scale:.6e} miss {f}", file=sys.stderr)
        if np.max(np.abs(f)) < 1e-10:
            break
        jac = np.empty((2, 2))
        for j, h in enumerate((1e-6, 1e-6 * max(abs(v[1]), 1.0))):
            vp = v.copy()
            vp[j] += h
            state["x"] = base
            jac[:, j] = (miss(vp) - f) / h
        state["x"] = base
        dv = np.linalg.solve(jac, -f)
        dv[0] = max(-0.3, min(0.3, dv[0]))
        v = v + dv
    else:
        raise RuntimeError("inversion did not converge")
    tau_w, q_w = layer.wall(state["x"])
    return tau_w, q_w, layer, state["x"]


def apriori(case, targets, checks, pipe_gradients=False):
    r, u, t, rho, tau, q, pw, tw = references(case)
    rho_w = pw / (GAS_R * tw)
    u_tau = math.sqrt(tau / rho_w)
    yplus = (RADIUS - r) * rho_w * u_tau / props(tw)[0]
    print(f"{case}: tau_w {tau:.6f} Pa, q_w {q:.6e} W/m^2, p_wall {pw:.1f} Pa, Re_tau {RADIUS * rho_w * u_tau / props(tw)[0]:.1f}",
          flush=True)
    print(f"{'y+':>9} {'T1/Tw':>8} {'tau model-1':>12} {'q model-1':>12}" + ("  c 20: tau, q change   cells x2: tau, q change" if checks else ""),
          flush=True)
    for target in targets:
        i = int(np.argmin(np.abs(yplus - target)))
        y1 = RADIUS - r[i]
        # Diagnostic: planar sources of half the pipe's force and heating give the pipe's linear fall of stress and
        # heat flux, tau_w (1 - y/R) and about q_w (1 - y/R), inside the layer.
        sources = (0.5 * CASES[case]["force"], 0.5 * CASES[case]["heat"]) if pipe_gradients else (0.0, 0.0)
        try:
            ts, qs, _, _ = model(u[i], t[i], y1, pw, tw, tau, q, sources=sources)
        except RuntimeError as err:
            print(f"{yplus[i]:9.2f} {t[i] / tw:8.4f}  no solution: {err}", flush=True)
            continue
        line = f"{yplus[i]:9.2f} {t[i] / tw:8.4f} {ts / tau - 1:+12.4e} {qs / q - 1:+12.4e}"
        if checks:
            tc, qc, _, _ = model(u[i], t[i], y1, pw, tw, ts, qs, c=20.0, sources=sources)
            tg, qg, _, _ = model(u[i], t[i], y1, pw, tw, ts, qs, cells=(400, 120), sources=sources)
            line += f"  {tc / ts - 1:+.2e} {qc / qs - 1:+.2e}   {tg / ts - 1:+.2e} {qg / qs - 1:+.2e}"
        print(line, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="heated_hp")
    ap.add_argument("--targets", default="1,5,30,100,300,1000,3000")
    ap.add_argument("--checks", action="store_true", help="also c = 20 and the cells doubled at every point")
    ap.add_argument("--pipe-gradients", action="store_true", help="diagnostic: the pipe's fall of stress and heat flux")
    a = ap.parse_args()
    apriori(a.case, [float(v) for v in a.targets.split(",")], a.checks, a.pipe_gradients)


if __name__ == "__main__":
    main()
