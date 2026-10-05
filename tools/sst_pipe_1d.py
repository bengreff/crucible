#!/usr/bin/env python3
"""Independent 1-D reference for fully developed turbulent pipe flow under SST-2003.

Reference for the engine's periodic-pipe verification (TECHNICAL_PLAN step 7, verification 4). It
shares no code with the engine: its own grid, discretisation and iteration, and it reads Cantera's
properties directly (tools/n2_properties.csv) where the engine uses its own fits.

Problem (steady, axially periodic, u_r = 0, pure N2 between the axis and a no-slip isothermal wall
at radius R, driven by a uniform axial body force f per unit volume):
  axial momentum   (1/r) d/dr [r (mu + mu_t) du/dr] + f = 0
  k                (1/r) d/dr [r (mu + sigma_k mu_t) dk/dr] + Pt - beta* rho omega k = 0
  omega            (1/r) d/dr [r (mu + sigma_w mu_t) dw/dr] + (gamma / nu_t) Pt - beta rho omega^2
                     + 2 (1 - F1) rho sigma_w2 / omega dk/dr dw/dr = 0
  total energy     (1/r) d/dr [r (q + u tau_rz + (mu + sigma_k mu_t) dk/dr)] + f u = 0,
                     q = (lambda + cp mu_t / Pr_t) dT/dr, tau_rz = (mu + mu_t) du/dr
  radial momentum  p + 2/3 rho k uniform (the Reynolds normal stresses -2/3 rho k)
  state            p = rho R_gas T, mass per unit length fixed by the fill (p0, T_wall)
with P = mu_t (du/dr)^2 (the -2/3 rho k div u part is zero here), Pt = min(P, 10 beta* rho omega k),
mu_t = rho a1 k / max(a1 omega, S F2), S = |du/dr|, and the SST-2003 constants and blending of the
NASA Turbulence Modeling Resource (sst.html). Wall: u = 0, T = T_wall, k = 0,
omega = 10 * 6 nu / (beta1 d1^2) with d1 the distance of the first cell centroid from the wall.

Finite volumes on rings r_j = R tanh(b j / N) / tanh(b), centroids in r dr, face gradients by
centroid differences, face coefficients the mean of the two cells (the wall face takes the wall's).
A segregated iteration on a coarse grid gives the start; each grid is then solved by Newton's
method (sparse finite-difference Jacobian, pseudo-time continuation), the unknowns u, k, omega, T
per cell and the uniform p + 2/3 rho k.

Usage: sst_pipe_1d.py --radius R --p0 P --twall T --force F --grids 100,200,400 [--stretch 2]
       [--prt 0.9] [--omega-factor 10] [--profile out.csv]
"""
import argparse
import math
import os
import sys

import numpy as np
from scipy.linalg import solve_banded
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

SIGMA_K = (0.85, 1.0)
SIGMA_W = (0.5, 0.856)
BETA = (0.075, 0.0828)
GAMMA = (5.0 / 9.0, 0.44)
BETA_STAR, A1 = 0.09, 0.31


def properties():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "n2_properties.csv")
    data = np.genfromtxt(path, delimiter=",", skip_header=1, names=True)
    t = data["T_K"]
    gas_r = float(data["R_J_kg_K"][0])
    return (lambda temp: np.interp(temp, t, data["mu_Pa_s"]),
            lambda temp: np.interp(temp, t, data["lambda_W_m_K"]),
            lambda temp: np.interp(temp, t, data["cp_J_kg_K"]), gas_r)


class Pipe:
    def __init__(self, radius, p0, t_wall, force, n, stretch=2.0, pr_t=0.9, omega_factor=10.0):
        self.mu_of, self.lam_of, self.cp_of, self.gas_r = properties()
        self.radius, self.t_wall, self.force, self.n, self.pr_t, self.omega_factor = radius, t_wall, force, n, pr_t, omega_factor
        rf = radius * np.tanh(stretch * np.arange(n + 1) / n) / np.tanh(stretch)
        rf[0], rf[-1] = 0.0, radius
        self.rf = rf
        self.vol = 0.5 * (rf[1:] ** 2 - rf[:-1] ** 2)              # per radian, per unit length
        self.rc = (2.0 / 3.0) * (rf[1:] ** 3 - rf[:-1] ** 3) / (rf[1:] ** 2 - rf[:-1] ** 2)
        self.d = radius - self.rc                                  # wall distance
        self.dc = np.diff(self.rc)
        self.dw = radius - self.rc[-1]
        self.mass = 0.5 * radius ** 2 * p0 / (self.gas_r * t_wall) # per radian, per unit length
        self.mu_wall, self.lam_wall = float(self.mu_of(t_wall)), float(self.lam_of(t_wall))

    def unpack(self, x):
        n = self.n
        return x[:n], x[n:2 * n], x[2 * n:3 * n], x[3 * n:4 * n], x[4 * n]

    def grad(self, phi, wall):
        face = np.empty(self.n + 1)
        face[0] = 0.0
        face[1:-1] = np.diff(phi) / self.dc
        face[-1] = (wall - phi[-1]) / self.dw
        return 0.5 * (face[:-1] + face[1:])

    def conductance(self, gamma_c, gamma_wall):
        c = np.empty(self.n + 1)
        c[0] = 0.0
        c[1:-1] = self.rf[1:-1] * 0.5 * (gamma_c[:-1] + gamma_c[1:]) / self.dc
        c[-1] = self.rf[-1] * gamma_wall / self.dw
        return c

    @staticmethod
    def face_flux(cond, phi, wall):
        """Inward-positive diffusive flux through each face (outer side of cell j is face j + 1)."""
        flux = np.zeros(len(phi) + 1)
        flux[1:-1] = cond[1:-1] * np.diff(phi)
        flux[-1] = cond[-1] * (wall - phi[-1])
        return flux

    def closure(self, x):
        u, k, w, t, pu = self.unpack(x)
        s = {}
        s["rho"] = rho = pu / (self.gas_r * t + 2.0 / 3.0 * k)
        s["mu"] = mu = self.mu_of(t)
        s["lam"], s["cp"] = self.lam_of(t), self.cp_of(t)
        rho_wall = pu / (self.gas_r * self.t_wall)
        s["w_wall"] = self.omega_factor * 6 * (self.mu_wall / rho_wall) / (BETA[0] * self.dw ** 2)
        dudr, dkdr, dwdr = self.grad(u, 0.0), self.grad(k, 0.0), self.grad(w, s["w_wall"])
        nu, d = mu / rho, self.d
        strain = np.abs(dudr)
        cd = np.maximum(2 * rho * SIGMA_W[1] / w * dkdr * dwdr, 1e-10)
        arg1 = np.minimum(np.maximum(np.sqrt(k) / (BETA_STAR * w * d), 500 * nu / (d ** 2 * w)),
                          4 * rho * SIGMA_W[1] * k / (cd * d ** 2))
        f1 = np.tanh(arg1 ** 4)
        arg2 = np.maximum(2 * np.sqrt(k) / (BETA_STAR * w * d), 500 * nu / (d ** 2 * w))
        f2 = np.tanh(arg2 ** 2)
        s["mu_t"] = mu_t = rho * A1 * k / np.maximum(A1 * w, strain * f2)
        blend = lambda pair: f1 * pair[0] + (1 - f1) * pair[1]
        s["sk"], s["sw"], s["beta"], s["gamma"] = blend(SIGMA_K), blend(SIGMA_W), blend(BETA), blend(GAMMA)
        s["prod"] = np.minimum(mu_t * strain ** 2, 10 * BETA_STAR * rho * w * k)
        s["cross"] = 2 * (1 - f1) * rho * SIGMA_W[1] * dkdr * dwdr  # C in C / omega
        s["f1"] = f1
        return s

    def residual(self, x):
        u, k, w, t, pu = self.unpack(x)
        s = self.closure(x)
        rho, mu, mu_t, vol = s["rho"], s["mu"], s["mu_t"], self.vol
        net = lambda flux: flux[1:] - flux[:-1]
        fu = self.face_flux(self.conductance(mu + mu_t, self.mu_wall), u, 0.0)
        fk = self.face_flux(self.conductance(mu + s["sk"] * mu_t, self.mu_wall), k, 0.0)
        fw = self.face_flux(self.conductance(mu + s["sw"] * mu_t, self.mu_wall), w, s["w_wall"])
        ft = self.face_flux(self.conductance(s["lam"] + s["cp"] * mu_t / self.pr_t, self.lam_wall), t, self.t_wall)
        uf = np.zeros(self.n + 1)
        uf[1:-1] = 0.5 * (u[:-1] + u[1:])
        work = uf * fu + fk                                        # shear work and k diffusion, inward
        ru = net(fu) + self.force * vol
        rk = net(fk) + (s["prod"] - BETA_STAR * rho * w * k) * vol
        rw = net(fw) + (s["gamma"] * rho * s["prod"] / mu_t - s["beta"] * rho * w ** 2 + s["cross"] / w) * vol
        rt = net(ft) + net(work) + self.force * u * vol
        rm = np.sum(rho * vol) - self.mass
        return np.concatenate([ru, rk, rw, rt, [rm]])

    def jacobian(self, x, r0):
        n = self.n
        rows, cols, vals = [], [], []
        floor = [1e-10 * np.max(np.abs(x[:n])), 1e-14 * np.max(x[n:2 * n]), 0.0, 1e-9 * self.t_wall]
        for v in range(4):
            for color in range(5):
                cells = np.arange(color, n, 5)
                h = 1e-7 * np.abs(x[v * n + cells]) + floor[v]
                xp = x.copy()
                xp[v * n + cells] += h
                dr = self.residual(xp) - r0
                for off in range(-2, 3):
                    rcell = cells + off
                    ok = (rcell >= 0) & (rcell < n)
                    for eq in range(4):
                        rows.append(eq * n + rcell[ok])
                        cols.append(v * n + cells[ok])
                        vals.append(dr[eq * n + rcell[ok]] / h[ok])
        hp = 1e-7 * x[4 * n]
        xp = x.copy()
        xp[4 * n] += hp
        dr = (self.residual(xp) - r0) / hp
        rows.append(np.arange(4 * n))
        cols.append(np.full(4 * n, 4 * n))
        vals.append(dr[:4 * n])
        # Mass row, analytic: sum rho V with rho = pu / (R T + 2/3 k).
        u, k, w, t, pu = self.unpack(x)
        den = self.gas_r * t + 2.0 / 3.0 * k
        rho = pu / den
        rows += [np.full(n, 4 * n), np.full(n, 4 * n), [4 * n]]
        cols += [3 * n + np.arange(n), n + np.arange(n), [4 * n]]
        vals += [-rho * self.gas_r * self.vol / den, -rho * (2.0 / 3.0) * self.vol / den, [np.sum(self.vol / den)]]
        return coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(4 * n + 1, 4 * n + 1)).tocsc()

    def newton(self, x, tol=1e-12, max_steps=400, log=True):
        n = self.n
        u_tau = math.sqrt(self.force * self.radius / (2 * self.mass / (0.5 * self.radius ** 2)))
        dtau = 1e-3 * self.radius / u_tau
        for step in range(max_steps):
            r0 = self.residual(x)
            jac = self.jacobian(x, r0)
            s = self.closure(x)
            m = np.concatenate([s["rho"] * self.vol] * 3 + [s["rho"] * s["cp"] * self.vol, [0.0]])
            from scipy.sparse import diags
            dx = spsolve((diags(m / dtau) - jac).tocsc(), r0)
            # Keep k and omega positive: at most halve them in one step.
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
                print(f"  N {n} step {step:3d}  alpha {alpha:.3f}  dtau {dtau:.2e}  largest relative update {change:.3e}",
                      file=sys.stderr)
            if alpha == 1.0 and change < tol and dtau * u_tau / self.radius > 1e6:
                return x, step + 1, change
        raise RuntimeError(f"Newton did not converge on N {n} (last update {change:.3e})")

    def start(self, coarse=None):
        """Initial state: interpolated from a coarser solution (in wall distance), else the log law."""
        n, d = self.n, self.d
        rho0 = self.mass / (0.5 * self.radius ** 2)
        nu = self.mu_wall / rho0
        if coarse is None:
            u_tau = math.sqrt(self.force * self.radius / (2 * rho0))
            yplus = d * u_tau / nu
            u = u_tau * np.minimum(yplus, 2.5 * np.log(np.maximum(yplus, 1e-30)) + 5.5)
            k = u_tau ** 2 / math.sqrt(BETA_STAR) * np.minimum(1.0, (yplus / 10.0) ** 2)
            w = np.maximum(u_tau / (math.sqrt(BETA_STAR) * 0.41 * d), 6 * nu / (BETA[0] * d ** 2))
            return np.concatenate([u, k, w, np.full(n, self.t_wall), [self.mass / np.sum(self.vol) * self.gas_r * self.t_wall]])
        pipe, xc = coarse
        uc, kc, wc, tc, pu = pipe.unpack(xc)
        dcoarse = np.concatenate([[0.0], pipe.d[::-1]])
        at = lambda phi, wall: np.interp(d, dcoarse, np.concatenate([[wall], phi[::-1]]))
        w = np.maximum(np.interp(d, pipe.d[::-1], wc[::-1]), 6 * nu / (BETA[0] * d ** 2))
        return np.concatenate([at(uc, 0.0), at(kc, 0.0), w, at(tc, self.t_wall), [pu]])

    def measures(self, x):
        u, k, w, t, pu = self.unpack(x)
        s = self.closure(x)
        tau_w = self.mu_wall * u[-1] / self.dw
        rho_wall = pu / (self.gas_r * self.t_wall)
        area = 0.5 * self.radius ** 2
        rho_b = self.mass / area
        u_b = np.sum(s["rho"] * u * self.vol) / (rho_b * area)
        u_tau = math.sqrt(tau_w / rho_wall)
        heat_in = self.lam_wall * (self.t_wall - t[-1]) / self.dw * self.radius
        k_in = self.mu_wall * (0.0 - k[-1]) / self.dw * self.radius
        return dict(tau_w=tau_w, force_balance=tau_w * self.radius / (self.force * area) - 1,
                    energy_balance=-(heat_in + k_in) / (self.force * np.sum(u * self.vol)) - 1,
                    u_b=u_b, rho_b=rho_b, c_f=2 * tau_w / (rho_b * u_b ** 2), u_tau=u_tau,
                    re_tau=rho_wall * u_tau * self.radius / self.mu_wall, re_b=rho_b * u_b * 2 * self.radius / self.mu_wall,
                    y1plus=rho_wall * u_tau * self.dw / self.mu_wall, t_axis=t[0], p_wall=pu,
                    profile=dict(r=self.rc, u=u, k=k, omega=w, T=t, rho=s["rho"], mu_t=s["mu_t"]))


def segregated(pipe, x, iterations=20000, tol=1e-8):
    """Under-relaxed segregated iteration (one tridiagonal solve per equation), for a starting state."""
    n = pipe.n
    for it in range(iterations):
        u, k, w, t, pu = pipe.unpack(x)
        s = pipe.closure(x)
        rho, mu, mu_t, vol = s["rho"], s["mu"], s["mu_t"], pipe.vol

        def solve(cond, sink, src, phi, wall, relax):
            ab = np.zeros((3, n))
            diag = cond[:-1] + cond[1:] + sink
            ab[1] = diag / relax
            ab[0, 1:] = -cond[1:-1]
            ab[2, :-1] = -cond[1:-1]
            b = src + (1 / relax - 1) * diag * phi
            b[-1] += cond[-1] * wall
            return solve_banded((1, 1), ab, b)

        u_new = solve(pipe.conductance(mu + mu_t, pipe.mu_wall), 0.0, pipe.force * vol, u, 0.0, 0.8)
        k_new = np.maximum(solve(pipe.conductance(mu + s["sk"] * mu_t, pipe.mu_wall), BETA_STAR * rho * w * vol,
                                 s["prod"] * vol, k, 0.0, 0.6), 1e-30)
        c = s["cross"]
        src = (s["gamma"] * rho * s["prod"] / mu_t + s["beta"] * rho * w ** 2 + np.where(c > 0, c / w, 2 * c / w)) * vol
        sink = (2 * s["beta"] * rho * w + np.where(c > 0, 0.0, -c / w ** 2)) * vol
        w_new = np.maximum(solve(pipe.conductance(mu + s["sw"] * mu_t, pipe.mu_wall), sink, src, w, s["w_wall"], 0.6), 1e-30)
        x_mid = np.concatenate([u_new, k_new, w_new, t, [pu]])
        # Energy with the updated velocity and k: everything but conduction as a source.
        r = pipe.residual(x_mid)[3 * n:4 * n]
        ct = pipe.conductance(s["lam"] + s["cp"] * mu_t / pipe.pr_t, pipe.lam_wall)
        conduction = pipe.face_flux(ct, t, pipe.t_wall)
        t_new = solve(ct, 0.0, r - (conduction[1:] - conduction[:-1]), t, pipe.t_wall, 0.8)
        shape = 1 / (pipe.gas_r * t_new + 2.0 / 3.0 * k_new)
        pu_new = pipe.mass / np.sum(shape * vol)
        x_new = np.concatenate([u_new, k_new, w_new, t_new, [pu_new]])
        change = max(np.max(np.abs(u_new - u)) / np.max(np.abs(u_new)), np.max(np.abs(k_new - k)) / np.max(k_new),
                     np.max(np.abs(w_new - w) / w_new), np.max(np.abs(t_new - t)) / pipe.t_wall)
        x = x_new
        if change < tol and it > 100:
            break
    return x


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--radius", type=float, required=True)
    ap.add_argument("--p0", type=float, required=True, help="fill pressure at the wall temperature (Pa)")
    ap.add_argument("--twall", type=float, required=True)
    ap.add_argument("--force", type=float, required=True, help="axial body force per unit volume (N/m^3)")
    ap.add_argument("--grids", required=True, help="comma-separated ring counts, coarse to fine")
    ap.add_argument("--stretch", type=float, default=2.0)
    ap.add_argument("--prt", type=float, default=0.9)
    ap.add_argument("--omega-factor", type=float, default=10.0)
    ap.add_argument("--profile", help="CSV of the finest profile")
    a = ap.parse_args()
    grids = [int(g) for g in a.grids.split(",")]
    first = Pipe(a.radius, a.p0, a.twall, a.force, min(grids[0], 50), a.stretch, a.prt, a.omega_factor)
    previous = (first, segregated(first, first.start()))
    rows = []
    for n in grids:
        pipe = Pipe(a.radius, a.p0, a.twall, a.force, n, a.stretch, a.prt, a.omega_factor)
        x, steps, change = pipe.newton(pipe.start(previous))
        m = pipe.measures(x)
        rows.append((n, m))
        print(f"N {n:5d}  Newton steps {steps:3d}  last update {change:.1e}  force balance {m['force_balance']:+.1e}  "
              f"energy balance {m['energy_balance']:+.1e}  y1+ {m['y1plus']:.4f}  Re_tau {m['re_tau']:.4f}  "
              f"u_b {m['u_b']:.8f} m/s  c_f {m['c_f']:.10f}  T_axis {m['t_axis']:.6f} K", flush=True)
        previous = (pipe, x)
    for key in ("c_f", "u_b", "t_axis"):
        vals = [m[key] for _, m in rows]
        for i in range(2, len(vals)):
            e1, e2 = vals[i - 1] - vals[i - 2], vals[i] - vals[i - 1]
            if e1 != 0 and e2 != 0 and e1 / e2 > 0:
                order = math.log(e1 / e2) / math.log(grids[i] / grids[i - 1])
                extrap = vals[i] + e2 / (2 ** order - 1)
                print(f"{key}: grids {grids[i - 2]}-{grids[i - 1]}-{grids[i]}  observed order {order:.3f}  "
                      f"Richardson {extrap:.10g}")
    if a.profile:
        prof = rows[-1][1]["profile"]
        with open(a.profile, "w") as fh:
            fh.write("r_m,u_m_s,k_m2_s2,omega_1_s,T_K,rho_kg_m3,mu_t_Pa_s\n")
            for row in zip(prof["r"], prof["u"], prof["k"], prof["omega"], prof["T"], prof["rho"], prof["mu_t"]):
                fh.write(",".join(f"{v:.12e}" for v in row) + "\n")


if __name__ == "__main__":
    main()
