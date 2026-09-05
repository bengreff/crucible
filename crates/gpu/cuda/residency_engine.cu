// crucible-gpu — S13c CLOSE: the COMPOSED full-physics 3-D resident ENGINE
// STEP — the ◆C3 world's whole SDC step on-device, state resident across the
// march behind a persistent handle. Composes every validated leg into ONE
// step, bit-for-formula from the CPU production path (`Euler::eval_rhs` +
// `Combustion::accumulate_inner` + `Sdc::step` with the flow + reaction classes
// + `Euler::srd`, as `crucible_engine::run` schedules them on
// `rl10_startup_3d`):
//   • the blend EOS `prim_checked_hinted` (unburnt ⊕ burnt (p,h,Z) surfaces,
//     the three-way branch select, warm-start hints on the pure limits, the
//     mass-weighted blended sound speed → the prim's e/Γ₁ aux slots);
//   • the general-EOS HLLC-Batten (Γ₁ aux slot, datum-free Roe c²);
//   • the r/θ/z sweeps on the real cut geometry (per-sector κ + apertures, the
//     axis parity pair, slip-wall ghosts) + the COUP-7 injector MASS-FLOW
//     INFLOW ghost (the unburnt TableEos fixed-count solve, sonic-capped) and
//     the PRESSURE-OUTFLOW ghost (the pump-down ambient);
//   • the SOLV-1 §3.3 geometric + wall-closure sources + the bounded spark
//     IGNITER deposit (the sector-gated ramped energy source);
//   • the SOLV-4.4 combustion source (Nagumo front + front-thickening
//     diffusion, θ-faces, per-sector κ, the loud [0,1] guard);
//   • the 2-node Lobatto IMEX-SDC composition with stagewise SRD and the
//     class-R implicit auto-ignition node solve after every accepted
//     composition (trapezoid-quadrature base, realized-rate roll-over);
//   • the COUP-2 ledger (port + source, net + gross) and the κV-weighted
//     stored totals as fixed-topology tree reductions — the audit's operands
//     (the audit CHECK itself runs host-side on these O(1) scalars, through
//     the CPU's own `Sdc` audit arithmetic);
//   • `stable_dt` with the θ-arc member + the front-carrier signal + the
//     scale-separation guard (max-tree).
//
// Determinism (META-1 §2.5): every physics kernel is a per-cell gather (one
// writer per cell, a face's flux identical from both sides); reductions are
// fixed-topology trees (shape a function of (N, TPB) only); the only atomics
// are the halt flag. Same-build reruns bit-identical; CPU↔GPU = FMA/libm ECT.
// f64 everywhere (ruling: never precision).
//
// Device NaN sentinels (the projection's tangency corner, a non-physical
// state) raise the halt flag — never a silent branch (META-1 P6).
#include <cuda_runtime.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NP 9
#define NC 7
#define NGH 3
#define I_RHO 0
#define I_MR 1
#define I_MT 2
#define I_MZ 3
#define I_EN 4
#define I_RC 5
#define I_RB 6
#define I_EI 7
#define I_G1 8

#define EPS_P_PROJECTION 1e-11
#define N_P_ITER_MAX 48
#define N_P_BISECT 64
#define N_P_SCAN 64
#define H_MARGIN_TABLE 1e-12
#define H_MARGIN_BLEND 1e-9
#define HINT_SPREAD 1.05
#define EPS_B_PURE_UNBURNT 1.0e-9
#define BURN_COMPLETE 1.0e-3
#define EPS_B_PURE_BURNT (2.0 * BURN_COMPLETE)
#define N_TAU_REFREEZE 2
#define N_INFLOW_ITER 8
#define IGN_A 0.2
#define C_NAGUMO_SLOPE 0.8
#define S_T_MACH_LIMIT (2.0 / 3.0)
#define EPS_BURN_BOUND 1.0e-4

#define GK_DOMAIN_REFLECT 0
#define GK_DOMAIN_TRANSMISSIVE 1
#define GK_WALL_MIRROR 2
#define GK_WALL_SLIP 3
#define GK_AXIS 4
#define GK_INFLOW 5     // COUP-7 mass-flow inflow (z_lo)
#define GK_OUTFLOW 6    // pressure outflow (z_hi)

#define RTPB 256
// ledger columns (column-major [4*NC][N]): port_net, port_abs, src_net, src_abs
#define LED_PN 0
#define LED_PA (NC)
#define LED_SN (2*NC)
#define LED_SA (3*NC)

__device__ __constant__ double E_TAU = 6.283185307179586;

// ---- tables ----------------------------------------------------------------------
struct Col3 {
    const double *pp, *hp, *zp; int np, nh, nz; int lp, lh, lz;
    int s0, s1, s2; const double* data; int vlog;
};
__device__ __forceinline__ int t_cell(const double* pts, int n, double q) {
    int cnt = 0; for (int i = 0; i < n; i++) cnt += (pts[i] <= q);
    int hi = cnt < n-1 ? cnt : n-1; return (hi > 1 ? hi : 1) - 1;
}
__device__ __forceinline__ double t_frac(const double* pts, int i, int is_log, double q) {
    if (is_log) return (log(q) - log(pts[i])) / (log(pts[i+1]) - log(pts[i]));
    return (q - pts[i]) / (pts[i+1] - pts[i]);
}
__device__ double t_interp(const Col3& c, double qp, double qh, double qz) {
    int ip = t_cell(c.pp, c.np, qp), ih = t_cell(c.hp, c.nh, qh), iz = t_cell(c.zp, c.nz, qz);
    double t0 = t_frac(c.pp, ip, c.lp, qp), t1 = t_frac(c.hp, ih, c.lh, qh), t2 = t_frac(c.zp, iz, c.lz, qz);
    int cell[3] = {ip, ih, iz}; double t[3] = {t0, t1, t2}; int str[3] = {c.s0, c.s1, c.s2};
    double acc = 0.0;
    for (int corner = 0; corner < 8; corner++) {
        double w = 1.0; int idx = 0;
        for (int d = 0; d < 3; d++) { int up = (corner >> d) & 1; w *= up ? t[d] : (1.0 - t[d]); idx += (cell[d] + up) * str[d]; }
        double v = c.vlog ? log(c.data[idx]) : c.data[idx];
        acc += w * v;
    }
    return c.vlog ? exp(acc) : acc;
}
struct BlendEnv {
    double pu_lo, pu_hi, hu_floor, hu_ceil;
    double pb_lo, pb_hi, hb_lo, hb_hi;
    double h_off;
};
struct Tables {
    Col3 urho, usnd, utmp;    // unburnt (p,h,Z): ρ, a, T
    Col3 brho, bsnd;          // burnt (p,h,Z): ρ, a
    Col3 flame, delay;        // ignition (p,T_u,Z): S_L, τ_ign
    BlendEnv v;
    double p_floor, tu_floor; // ignition surface non-reactive floor
};

// ---- blend partition + closure (blend_eos.rs) ---------------------------------------
__device__ __forceinline__ void b_partition(const BlendEnv& v, double h, double b, double* h_u, double* h_b) {
    if ((h >= v.hu_floor && h <= v.hu_ceil) || b <= EPS_B_PURE_UNBURNT) { *h_u = h; *h_b = h; return; }
    double hu_pin = fmin(fmax(h, v.hu_floor), v.hu_ceil);
    double balance = (h - (1.0 - b) * hu_pin) / b;
    *h_u = hu_pin; *h_b = fmin(fmax(balance, v.hb_lo - v.h_off), v.hb_hi - v.h_off);
}
__device__ double b_inv_rho(const Tables& T, double p, double h, double z, double b) {
    double h_u, h_b; b_partition(T.v, h, b, &h_u, &h_b);
    double vv = 0.0;
    if (b < 1.0 - EPS_B_PURE_BURNT) { double ru = t_interp(T.urho, p, h_u, z); vv += (1.0 - b) / ru; }
    if (b > EPS_B_PURE_UNBURNT)     { double rb = t_interp(T.brho, p, h_b + T.v.h_off, z); vv += b / rb; }
    return vv;
}
struct GTable { const Col3* rc; double e_q, inv, z, rho;
    __device__ double operator()(double p) const { return t_interp(*rc, p, e_q + p * inv, z) - rho; } };
struct GBlend { const Tables* T; double e, inv, z, b;
    __device__ double operator()(double p) const { return b_inv_rho(*T, p, e + p * inv, z, b) - inv; } };

template <typename G>
__device__ double illinois(double a, double b, double ga, double gb, const G& g) {
    for (int it = 0; it < N_P_ITER_MAX; it++) {
        if (fabs(b - a) <= EPS_P_PROJECTION * fmax(fabs(a), fabs(b))) return 0.5 * (a + b);
        double denom = gb - ga;
        double p = (denom != 0.0) ? b - gb * (b - a) / denom : 0.5 * (a + b);
        double loe = fmin(a, b), hie = fmax(a, b);
        if (!(p > loe && p < hie)) p = 0.5 * (a + b);
        double gp = g(p);
        if (gp == 0.0) return p;
        if (gp * gb < 0.0) { a = b; ga = gb; } else { ga *= 0.5; }
        b = p; gb = gp;
    }
    if (fabs(b - a) <= EPS_P_PROJECTION * fmax(fabs(a), fabs(b)) * 10.0) return 0.5 * (a + b);
    double aa, bb, gaa;
    if (a < b) { aa = a; bb = b; gaa = ga; } else { aa = b; bb = a; gaa = gb; }
    for (int it = 0; it < N_P_BISECT; it++) {
        double m = 0.5 * (aa + bb);
        if (fabs(bb - aa) <= EPS_P_PROJECTION * fmax(fabs(aa), fabs(bb))) return m;
        double gm = g(m);
        if (gm == 0.0) return m;
        if (gaa * gm < 0.0) bb = m; else { aa = m; gaa = gm; }
    }
    return 0.5 * (aa + bb);
}
template <typename G>
__device__ double scan_first(double lo, double hi, double ga0, const G& g) {
    double ratio = hi / lo, prev_p = lo, prev_g = ga0;
    for (int k = 1; k <= N_P_SCAN; k++) {
        double pk = lo * pow(ratio, (double)k / (double)N_P_SCAN);
        double gk = g(pk);
        if (gk == 0.0) return pk;
        if (prev_g * gk < 0.0) return illinois(prev_p, pk, prev_g, gk, g);
        prev_p = pk; prev_g = gk;
    }
    return nan("");   // tangency corner — halt (never a silent branch)
}
// TableEos::project_pressure (cold).
__device__ double project_table_cold(const Col3& rc, double rho, double e_q, double z, double inv,
                                     double p_lo_env, double p_hi_env, double h_lo_env, double h_hi_env) {
    double p_from_h_lo = rho * (h_lo_env - e_q), p_from_h_hi = rho * (h_hi_env - e_q);
    double lo = fmax(p_lo_env, p_from_h_lo * (1.0 + H_MARGIN_TABLE));
    double hi = fmin(p_hi_env, p_from_h_hi * (1.0 - H_MARGIN_TABLE));
    if (!(isfinite(lo) && isfinite(hi)) || lo >= hi || hi <= 0.0) return nan("");
    lo = fmax(lo, 0.0);
    GTable g{&rc, e_q, inv, z, rho};
    double ga0 = g(lo), gb0 = g(hi);
    if (ga0 == 0.0) return lo;
    if (gb0 == 0.0) return hi;
    if (ga0 * gb0 < 0.0) return illinois(lo, hi, ga0, gb0, g);
    return scan_first(lo, hi, ga0, g);
}
// TableEos::project_pressure_hinted (warm; falls through to cold).
__device__ double project_table_hinted(const Col3& rc, double rho, double e_q, double z, double inv,
                                       double p_hint, double p_lo_env, double p_hi_env,
                                       double h_lo_env, double h_hi_env) {
    double p_from_h_lo = rho * (h_lo_env - e_q), p_from_h_hi = rho * (h_hi_env - e_q);
    double lo_adm = fmax(fmax(p_lo_env, p_from_h_lo * (1.0 + H_MARGIN_TABLE)), 0.0);
    double hi_adm = fmin(p_hi_env, p_from_h_hi * (1.0 - H_MARGIN_TABLE));
    double a = fmax(p_hint / HINT_SPREAD, lo_adm);
    double b = fmin(p_hint * HINT_SPREAD, hi_adm);
    if (!(isfinite(a) && isfinite(b)) || a >= b)
        return project_table_cold(rc, rho, e_q, z, inv, p_lo_env, p_hi_env, h_lo_env, h_hi_env);
    GTable g{&rc, e_q, inv, z, rho};
    double ga_full = g(lo_adm), gb_full = g(hi_adm);
    if (!(ga_full * gb_full < 0.0))
        return project_table_cold(rc, rho, e_q, z, inv, p_lo_env, p_hi_env, h_lo_env, h_hi_env);
    double ga = g(a);
    if (ga == 0.0) return a;
    double gb = g(b);
    if (gb == 0.0) return b;
    if (ga * gb < 0.0) return illinois(a, b, ga, gb, g);
    return illinois(lo_adm, hi_adm, ga_full, gb_full, g);
}
// BurnBlendEos::project_pressure, the mid-b path: the SOLV-4 0.4.10 warm start
// (full-bracket straddle ⇒ Illinois on the tight bracket about the hint, or on
// the full bracket), else the cold first-crossing scan.
__device__ double project_mid(const Tables& T, double rho, double e, double z, double b, double inv,
                              double hint, int has_hint) {
    const BlendEnv& v = T.v;
    double h_lo = v.hb_lo - v.h_off, h_hi = v.hb_hi - v.h_off;
    double p_lo_env = fmax(v.pu_lo, v.pb_lo), p_hi_env = fmin(v.pu_hi, v.pb_hi);
    double p_from_h_lo = rho * (h_lo - e), p_from_h_hi = rho * (h_hi - e);
    double lo = fmax(fmax(p_lo_env, p_from_h_lo * (1.0 + H_MARGIN_BLEND)), 0.0);
    double hi = fmin(p_hi_env, p_from_h_hi * (1.0 - H_MARGIN_BLEND));
    if (!(isfinite(lo) && isfinite(hi)) || lo >= hi || hi <= 0.0) return nan("");
    lo = fmax(lo, 0.0);
    GBlend g{&T, e, inv, z, b};
    double ga0 = g(lo), gb0 = g(hi);
    if (ga0 == 0.0) return lo;
    if (gb0 == 0.0) return hi;
    if (has_hint && isfinite(hint) && hint > 0.0 && ga0 * gb0 < 0.0) {
        double a = fmax(hint / HINT_SPREAD, lo), bb = fmin(hint * HINT_SPREAD, hi);
        if (a < bb) {
            double ga = g(a);
            double gb = g(bb);
            if (ga == 0.0) return a;
            if (gb == 0.0) return bb;
            if (ga * gb < 0.0) return illinois(a, bb, ga, gb, g);
            return illinois(lo, hi, ga0, gb0, g);
        }
    }
    return scan_first(lo, hi, ga0, g);
}
// BurnBlendEos::prim_checked_hinted → w (NP). Returns 0 on non-physical/NaN.
__device__ int prim_blend(const Tables& T, const double* u, double hint, int has_hint, double* w) {
    double rho = u[I_RHO];
    if (!isfinite(rho) || rho <= 0.0) return 0;
    double inv = 1.0 / rho;
    double ur = u[I_MR] * inv, ut = u[I_MT] * inv, uz = u[I_MZ] * inv;
    double e = u[I_EN] * inv - 0.5 * (ur * ur + ut * ut + uz * uz);
    if (!isfinite(e)) return 0;
    double z = u[I_RC] * inv;
    if (!isfinite(z)) return 0;
    double b_raw = u[I_RB] * inv;
    if (!isfinite(b_raw)) return 0;
    double b = fmin(fmax(b_raw, 0.0), 1.0);
    bool need_u = b < 1.0 - EPS_B_PURE_BURNT, need_b = b > EPS_B_PURE_UNBURNT;
    bool hinted = has_hint && isfinite(hint) && hint > 0.0;
    const BlendEnv& v = T.v;
    double p;
    if (!need_b) {
        p = hinted ? project_table_hinted(T.urho, rho, e, z, inv, hint, v.pu_lo, v.pu_hi, v.hu_floor, v.hu_ceil)
                   : project_table_cold(T.urho, rho, e, z, inv, v.pu_lo, v.pu_hi, v.hu_floor, v.hu_ceil);
    } else if (!need_u) {
        double e_q = e + v.h_off;
        p = hinted ? project_table_hinted(T.brho, rho, e_q, z, inv, hint, v.pb_lo, v.pb_hi, v.hb_lo, v.hb_hi)
                   : project_table_cold(T.brho, rho, e_q, z, inv, v.pb_lo, v.pb_hi, v.hb_lo, v.hb_hi);
    } else {
        p = project_mid(T, rho, e, z, b, inv, hint, hinted ? 1 : 0);
    }
    if (!isfinite(p) || p <= 0.0) return 0;
    double h = e + p * inv;
    double h_u, h_b; b_partition(v, h, b, &h_u, &h_b);
    double a2 = 0.0;
    if (need_u) { double au = t_interp(T.usnd, p, h_u, z); a2 += (1.0 - b) * au * au; }
    if (need_b) { double ab = t_interp(T.bsnd, p, h_b + v.h_off, z); a2 += b * ab * ab; }
    double g1 = rho * a2 / p;
    if (!isfinite(g1)) return 0;
    w[0] = rho; w[1] = ur; w[2] = ut; w[3] = uz; w[4] = p; w[5] = z; w[6] = b_raw; w[7] = e; w[8] = g1;
    return 1;
}

// ---- general-EOS HLLC-Batten (hllc.rs on the aux slots) ------------------------------
__device__ __forceinline__ double e_total_energy(const double* w) {
    return w[0] * (w[7] + 0.5 * (w[1]*w[1] + w[2]*w[2] + w[3]*w[3]));
}
__device__ __forceinline__ double e_sound(const double* w) { return sqrt(w[8] * w[4] / w[0]); }
__device__ void e_physical_flux(const double* w, int n, double* f) {
    double rho = w[0], p = w[4], un = w[n], m = rho * un, e = e_total_energy(w);
    f[0] = m; f[1] = m * w[1]; f[2] = m * w[2]; f[3] = m * w[3]; f[n] += p;
    f[4] = un * (e + p); f[5] = m * w[5]; f[6] = m * w[6];
}
__device__ void e_hllc(const double* wl, const double* wr, int n, double* f) {
    double rho_l = wl[0], p_l = wl[4], rho_r = wr[0], p_r = wr[4];
    double un_l = wl[n], un_r = wr[n];
    double c_l = e_sound(wl), c_r = e_sound(wr);
    double sql = sqrt(rho_l), sqr = sqrt(rho_r), inv = 1.0 / (sql + sqr);
    double u1 = (sql*wl[1] + sqr*wr[1]) * inv;
    double u2 = (sql*wl[2] + sqr*wr[2]) * inv;
    double u3 = (sql*wl[3] + sqr*wr[3]) * inv;
    double urn = (n == 1 ? u1 : (n == 2 ? u2 : u3));
    // h_roe / q2_roe are formed (as the CPU does) but the datum-free Roe c²
    // does not read them.
    double cl2 = wl[8] * wl[4] / wl[0], cr2 = wr[8] * wr[4] / wr[0];
    double c_roe = sqrt((sql * cl2 + sqr * cr2) * inv);
    double s_l = fmin(un_l - c_l, urn - c_roe);
    double s_r = fmax(un_r + c_r, urn + c_roe);
    double ml = rho_l * (s_l - un_l), mr = rho_r * (s_r - un_r);
    double s_m = (mr * un_r - ml * un_l + p_l - p_r) / (mr - ml);
    for (int k = 0; k < NC; k++) f[k] = 0.0;
    if (s_l >= 0.0) { e_physical_flux(wl, n, f); return; }
    if (s_r <= 0.0) { e_physical_flux(wr, n, f); return; }
    const double* w = (s_m >= 0.0) ? wl : wr; double s_k = (s_m >= 0.0) ? s_l : s_r;
    double rho = w[0], p = w[4], un = w[n], e = e_total_energy(w);
    double p_star = rho * (un - s_k) * (un - s_m) + p;
    double u_k[NC] = {rho, rho*w[1], rho*w[2], rho*w[3], e, rho*w[5], rho*w[6]};
    double fac = (s_k - un) / (s_k - s_m), rho_s = rho * fac;
    double u_s[NC];
    u_s[0] = rho_s; u_s[1] = rho_s*w[1]; u_s[2] = rho_s*w[2]; u_s[3] = rho_s*w[3];
    u_s[n] = rho_s * s_m;
    u_s[4] = fac * e + (p_star * s_m - p * un) / (s_k - s_m);
    u_s[5] = rho_s*w[5]; u_s[6] = rho_s*w[6];
    double fk[NC]; for (int k = 0; k < NC; k++) fk[k] = 0.0; e_physical_flux(w, n, fk);
    for (int k = 0; k < NC; k++) f[k] = fk[k] + s_k * (u_s[k] - u_k[k]);
}

// ---- PPM on a 7-cell pencil (recon.rs) -----------------------------------------------
__device__ __forceinline__ double mc_slope(double wm, double w0, double wp) {
    double dl = w0 - wm, dr = wp - w0;
    if (dl * dr <= 0.0) return 0.0;
    double dc = 0.5 * (dl + dr);
    double a = fabs(dc); a = fmin(a, 2.0*fabs(dl)); a = fmin(a, 2.0*fabs(dr));
    return copysign(a, dc);
}
__device__ void p_slope(const double* wm, const double* w0, const double* wp, double* s) {
    for (int k = 0; k < NP; k++) s[k] = mc_slope(wm[k], w0[k], wp[k]);
}
__device__ void p_iface(const double* wa, const double* wb, const double* sa, const double* sb, double* out) {
    for (int k = 0; k < NP; k++) out[k] = 0.5*(wa[k]+wb[k]) - (sb[k]-sa[k])/6.0;
}
__device__ void p_edge(const double* lo_if, const double* hi_if, const double* c, double* elo, double* ehi) {
    for (int k = 0; k < NP; k++) {
        double lo = lo_if[k], hi = hi_if[k], cc = c[k];
        if ((hi-cc)*(cc-lo) <= 0.0) { lo = cc; hi = cc; }
        else {
            double d = hi-lo, six = 6.0*(cc-0.5*(lo+hi));
            if (d*six > d*d) lo = 3.0*cc-2.0*hi;
            else if (d*six < -(d*d)) hi = 3.0*cc-2.0*lo;
        }
        elo[k] = lo; ehi[k] = hi;
    }
}
__device__ void pencil_faces(const double* p, double* lfl, double* lfr, double* rfl, double* rfr) {
    double s1[NP],s2[NP],s3[NP],s4[NP],s5[NP];
    p_slope(p+0*NP, p+1*NP, p+2*NP, s1);
    p_slope(p+1*NP, p+2*NP, p+3*NP, s2);
    p_slope(p+2*NP, p+3*NP, p+4*NP, s3);
    p_slope(p+3*NP, p+4*NP, p+5*NP, s4);
    p_slope(p+4*NP, p+5*NP, p+6*NP, s5);
    double if1[NP],if2[NP],if3[NP],if4[NP];
    p_iface(p+1*NP, p+2*NP, s1, s2, if1);
    p_iface(p+2*NP, p+3*NP, s2, s3, if2);
    p_iface(p+3*NP, p+4*NP, s3, s4, if3);
    p_iface(p+4*NP, p+5*NP, s4, s5, if4);
    double e2lo[NP],e2hi[NP],e3lo[NP],e3hi[NP],e4lo[NP],e4hi[NP];
    p_edge(if1, if2, p+2*NP, e2lo, e2hi);
    p_edge(if2, if3, p+3*NP, e3lo, e3hi);
    p_edge(if3, if4, p+4*NP, e4lo, e4hi);
    for (int k = 0; k < NP; k++) { lfl[k] = e2hi[k]; lfr[k] = e3lo[k]; rfl[k] = e3hi[k]; rfr[k] = e4lo[k]; }
}

// ---- the world + per-step parameters --------------------------------------------------
struct World {
    int n_r, n_z, nt;
    double r_min, dr, z_min, dz;
    const int* act; const double* kappa; const double* ap; int has_geom;
    const int* rs_r; const int* rl_r; const int* klo_r; const int* khi_r; const double* nlo_r; const double* nhi_r;
    const int* rs_z; const int* rl_z; const int* klo_z; const int* khi_z; const double* nlo_z; const double* nhi_z;
};
struct Igniter {
    int armed; double r_lo, r_hi, z_lo, z_hi; int theta_gated; int nt; int jt;
    double t_on, t_off, t_ramp, q0;
};
struct StepParams {
    double t;              // the evaluation time (sources, BC schedules)
    double mdot_per_area;  // injector inflow
    double h_total, c_frac;
    double p_amb;          // pressure-outflow ambient at t
    double wrinkling, theta;
};

__device__ __forceinline__ double face_r(const World& W, int f) { return W.r_min + (double)f * W.dr; }
__device__ __forceinline__ double area_r(const World& W, int f) { return face_r(W, f) * (E_TAU / (double)W.nt) * W.dz; }
__device__ __forceinline__ double vol_of(const World& W, int i_r) {
    double ri = face_r(W, i_r), ro = face_r(W, i_r + 1);
    return 0.5 * (ro*ro - ri*ri) * (E_TAU / (double)W.nt) * W.dz;
}
__device__ __forceinline__ double area_z(const World& W, int i_r) {
    double ri = face_r(W, i_r), ro = face_r(W, i_r + 1);
    return 0.5 * (ro*ro - ri*ri) * (E_TAU / (double)W.nt);
}
__device__ __forceinline__ double area_th(const World& W) { return W.dr * W.dz; }
__device__ __forceinline__ double r_center(const World& W, int i_r) { return W.r_min + ((double)i_r + 0.5) * W.dr; }
__device__ __forceinline__ double z_center(const World& W, int i_z) { return W.z_min + ((double)i_z + 0.5) * W.dz; }
__device__ __forceinline__ double theta_center(int j, int nt) { return ((double)j + 0.5) * E_TAU / (double)nt; }
__device__ __forceinline__ void slip(double* m, double nr, double nz) {
    double vn = m[1]*nr + m[3]*nz; m[1] -= 2.0*vn*nr; m[3] -= 2.0*vn*nz;
}

// TableEos::mass_flow_inflow_ghost on the UNBURNT surface (h_offset 0), sign/normal.
__device__ int inflow_ghost(const Tables& T, double mdot_per_area, double h_total, double c_frac,
                            double p_int, double sign, int normal, double* m) {
    if (!(isfinite(p_int) && p_int > 0.0)) return 0;
    double h_s = h_total, rho = nan(""), u = 0.0;
    for (int it = 0; it < N_INFLOW_ITER; it++) {
        rho = t_interp(T.urho, p_int, h_s, c_frac);
        double a = t_interp(T.usnd, p_int, h_s, c_frac);
        u = fmin(mdot_per_area / rho, a);
        h_s = h_total - 0.5 * u * u;
    }
    double a = t_interp(T.usnd, p_int, h_s, c_frac);
    double e_true = h_s - p_int / rho;
    double g1 = rho * a * a / p_int;
    m[0] = rho; m[1] = 0.0; m[2] = 0.0; m[3] = 0.0; m[4] = p_int; m[5] = c_frac; m[6] = 0.0; m[7] = e_true; m[8] = g1;
    m[normal] = sign * u;
    return isfinite(rho) && isfinite(g1);
}

// Gather one pencil entry of a meridional sweep (dir 1 = r, 3 = z), ghosts by
// the run-end rule (fill_ghosts_low/high, wall_ghosts_low/high, the axis pair,
// the injector inflow, the pressure outflow). Returns 0 on a bad ghost.
__device__ int gather_rz(const World& W, const Tables& T, const StepParams& P, const double* __restrict__ prim,
                         int j, int i_r, int i_z, int dir, int m, double* out) {
    long NRZ = (long)W.n_r * W.n_z;
    int rz = i_r * W.n_z + i_z;
    int start, len, klo, khi; const double *nlo, *nhi; int normal, pos;
    if (dir == 1) { start = W.rs_r[rz]; len = W.rl_r[rz]; klo = W.klo_r[rz]; khi = W.khi_r[rz];
                    nlo = W.nlo_r + 2*rz; nhi = W.nhi_r + 2*rz; normal = 1; pos = i_r; }
    else          { start = W.rs_z[rz]; len = W.rl_z[rz]; klo = W.klo_z[rz]; khi = W.khi_z[rz];
                    nlo = W.nlo_z + 2*rz; nhi = W.nhi_z + 2*rz; normal = 3; pos = i_z; }
    int q = pos - start + m;
    int src; int jsrc = j; int kind = -1; const double* nh = nlo; double sign = 1.0;
    if (q >= 0 && q < len) { src = start + q; }
    else if (q < 0) {
        int k = -q; kind = klo; nh = nlo; sign = 1.0;
        src = (kind == GK_DOMAIN_TRANSMISSIVE || kind == GK_INFLOW || kind == GK_OUTFLOW)
              ? start : start + min(k - 1, len - 1);
        if (kind == GK_AXIS) jsrc = (j + W.nt / 2) % W.nt;
    } else {
        int k = q - len + 1; kind = khi; nh = nhi; sign = -1.0;
        src = (kind == GK_DOMAIN_TRANSMISSIVE || kind == GK_INFLOW || kind == GK_OUTFLOW)
              ? start + len - 1 : start + max(len - k, 0);
    }
    long cs = (dir == 1) ? ((long)jsrc * NRZ + (long)src * W.n_z + i_z)
                         : ((long)jsrc * NRZ + (long)i_r * W.n_z + src);
    for (int k = 0; k < NP; k++) out[k] = prim[cs*NP + k];
    if (kind == GK_DOMAIN_REFLECT || kind == GK_WALL_MIRROR) out[normal] = -out[normal];
    else if (kind == GK_WALL_SLIP) slip(out, nh[0], nh[1]);
    else if (kind == GK_AXIS) { out[1] = -out[1]; out[2] = -out[2]; }
    else if (kind == GK_INFLOW) {
        double p_int = out[4];
        if (!inflow_ghost(T, P.mdot_per_area, P.h_total, P.c_frac, p_int, sign, normal, out)) return 0;
    } else if (kind == GK_OUTFLOW) {
        // Outflow through a LOW face is velocity toward −normal; HIGH face +normal.
        double out_mach = (sign > 0.0 ? -out[normal] : out[normal]) / e_sound(out);
        if (!(out_mach >= 1.0)) out[4] = P.p_amb;
    }
    return 1;
}

// ---- combustion closures (combustion.rs) ------------------------------------------------
__device__ __forceinline__ void front_coeffs(double s_t, double delta, double theta, double* d_c, double* k) {
    if (s_t <= 0.0) { *d_c = 0.0; *k = 0.0; return; }
    double w = theta * delta, f = 1.0 - 2.0 * IGN_A;
    *d_c = w * s_t / f; *k = 2.0 * s_t / (f * w);
}
// (react, ρD_c) at a prim (accumulate_inner's per-cell body + rho_dc_of), plus
// T_u/τ operands for the class-R/front-carrier paths. `live` = reactive.
struct Closure { int live; double t_u, rho_u, s_t, d_c, k, tau; };
__device__ Closure closure_at(const Tables& T, const StepParams& P, const double* w, double delta, int want_tau) {
    Closure c; c.live = 0; c.t_u = 0; c.rho_u = 0; c.s_t = 0; c.d_c = 0; c.k = 0; c.tau = 0;
    double rho = w[0], p = w[4], z = w[5];
    double b_c = fmin(fmax(w[6], 0.0), 1.0);
    double h = w[7] + p / rho;
    if (b_c >= 1.0 - BURN_COMPLETE || h < T.v.hu_floor) return c;
    double h_u, h_b; b_partition(T.v, h, b_c, &h_u, &h_b);
    double t_u = t_interp(T.utmp, p, h_u, z);
    if (p < T.p_floor || t_u < T.tu_floor) return c;
    c.live = 1; c.t_u = t_u;
    c.rho_u = t_interp(T.urho, p, h_u, z);
    double s_l = t_interp(T.flame, p, t_u, z);
    c.s_t = fmax(s_l * P.wrinkling, 0.0);
    front_coeffs(c.s_t, delta, P.theta, &c.d_c, &c.k);
    if (want_tau) c.tau = t_interp(T.delay, p, t_u, z);
    return c;
}

// ---- S14 FLUX-BUFFER path: each face's PPM+HLLC ONCE, register-lean ---------------------------
// The fused rate kernels hold a 7-cell pencil + two full face reconstructions
// (5 slopes × 4 interface values × 6 edges × NP doubles live) — 230 registers +
// spill (measured). Here a cell computes only its LOW face from a 6-cell
// pencil, component by component with ~12 live scalars (the SAME formulas as
// `pencil_faces`: the low face of cell c+1 IS the high face of cell c, bit for
// bit), into a device flux buffer; a gather-divergence kernel then forms the
// well-balanced single difference. Run-end cells additionally compute their
// high (ghost-side) face. Selected by CRUCIBLE_GPU_FLUXBUF=1 (both paths kept
// for the A/B measurement; their outputs are bit-identical by construction).
__device__ __forceinline__ void edge_k(double lo_if, double hi_if, double cc, double* elo, double* ehi) {
    double lo = lo_if, hi = hi_if;
    if ((hi-cc)*(cc-lo) <= 0.0) { lo = cc; hi = cc; }
    else {
        double d = hi-lo, six = 6.0*(cc-0.5*(lo+hi));
        if (d*six > d*d) lo = 3.0*cc-2.0*hi;
        else if (d*six < -(d*d)) hi = 3.0*cc-2.0*lo;
    }
    *elo = lo; *ehi = hi;
}
// pen: 6 cells (pencil offsets −3..+2 about the face's high-side cell); the
// face lies between pen[2] and pen[3]. wl = hi edge of pen[2], wr = lo edge of pen[3].
__device__ void face_lean(const double* pen, double* wl, double* wr) {
    for (int k = 0; k < NP; k++) {
        double p0 = pen[0*NP+k], p1 = pen[1*NP+k], p2 = pen[2*NP+k], p3 = pen[3*NP+k], p4 = pen[4*NP+k], p5 = pen[5*NP+k];
        double s1 = mc_slope(p0, p1, p2), s2 = mc_slope(p1, p2, p3), s3 = mc_slope(p2, p3, p4), s4 = mc_slope(p3, p4, p5);
        double if1 = 0.5*(p1+p2) - (s2-s1)/6.0;
        double if2 = 0.5*(p2+p3) - (s3-s2)/6.0;
        double if3 = 0.5*(p3+p4) - (s4-s3)/6.0;
        double lo2, hi2, lo3, hi3;
        edge_k(if1, if2, p2, &lo2, &hi2);
        edge_k(if2, if3, p3, &lo3, &hi3);
        wl[k] = hi2; wr[k] = lo3;
    }
}
// Per cell: the LOW-face flux (dir 1 = r, 3 = z) into flux_lo; run-end cells
// also their HIGH-face flux into flux_hi. Inactive cells: untouched.
__global__ void k_face_rz(const double* __restrict__ prim, const World W, const Tables T, const StepParams P,
                          long N, int dir, double* __restrict__ flux_lo, double* __restrict__ flux_hi,
                          int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z, i_z = rz % W.n_z;
    if (!W.act[rz]) return;
    double pen[6*NP];
    for (int m = -3; m <= 2; m++)
        if (!gather_rz(W, T, P, prim, j, i_r, i_z, dir, m, pen + (m+3)*NP)) { atomicExch(bad, 1); return; }
    double wl[NP], wr[NP], f[NC];
    face_lean(pen, wl, wr);
    e_hllc(wl, wr, dir, f);
    for (int k = 0; k < NC; k++) flux_lo[c*NC+k] = f[k];
    int start = (dir == 1) ? W.rs_r[rz] : W.rs_z[rz];
    int len = (dir == 1) ? W.rl_r[rz] : W.rl_z[rz];
    int pos = (dir == 1) ? i_r : i_z;
    if (pos + 1 == start + len) {
        for (int m = -2; m <= 3; m++)
            if (!gather_rz(W, T, P, prim, j, i_r, i_z, dir, m, pen + (m+2)*NP)) { atomicExch(bad, 1); return; }
        face_lean(pen, wl, wr);
        e_hllc(wl, wr, dir, f);
        for (int k = 0; k < NC; k++) flux_hi[c*NC+k] = f[k];
    }
}
// The divergence: r WRITES the rate (+ zeroes the ledger), z ADDS — the same
// single-difference/aperture/ledger arithmetic as k_rate_r / k_rate_z.
__global__ void k_div_rz(const double* __restrict__ flux_lo, const double* __restrict__ flux_hi, const World W,
                         long N, int dir, double* __restrict__ rate, double* __restrict__ led) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int rz = (int)(c % NRZ);
    int i_r = rz / W.n_z, i_z = rz % W.n_z;
    if (dir == 1) for (int col = 0; col < 4*NC; col++) led[(long)col*N + c] = 0.0;
    if (!W.act[rz]) { if (dir == 1) for (int k = 0; k < NC; k++) rate[c*NC+k] = 0.0; return; }
    int start = (dir == 1) ? W.rs_r[rz] : W.rs_z[rz];
    int len = (dir == 1) ? W.rl_r[rz] : W.rl_z[rz];
    int pos = (dir == 1) ? i_r : i_z;
    long stride = (dir == 1) ? (long)W.n_z : 1L;
    int last = (pos + 1 == start + len);
    const double* fL = flux_lo + c*NC;
    const double* fR = last ? flux_hi + c*NC : flux_lo + (c + stride)*NC;
    if (dir == 1) {
        double ap_lo = W.ap[0*N + c];
        double ap_hi = last ? W.ap[1*N + c] : W.ap[0*N + (c + W.n_z)];
        double a_lo = area_r(W, i_r) * ap_lo, a_hi = area_r(W, i_r + 1) * ap_hi;
        double kv = W.kappa[c] * vol_of(W, i_r);
        for (int k = 0; k < NC; k++) {
            double afl = a_lo * fL[k], afh = a_hi * fR[k];
            rate[c*NC+k] = (afl - afh) / kv;
            if (i_r == start) { led[(long)(LED_PN+k)*N + c] += afl;  led[(long)(LED_PA+k)*N + c] += fabs(afl); }
            if (last)         { led[(long)(LED_PN+k)*N + c] -= afh;  led[(long)(LED_PA+k)*N + c] += fabs(afh); }
        }
    } else {
        double ap_lo = W.ap[2*N + c];
        double ap_hi = last ? W.ap[3*N + c] : W.ap[2*N + (c + 1)];
        double inv_dz = 1.0 / W.dz, kap = W.kappa[c], a_z = area_z(W, i_r);
        for (int k = 0; k < NC; k++) {
            double afl = fL[k] * ap_lo, afh = fR[k] * ap_hi;
            rate[c*NC+k] += (afl - afh) * inv_dz / kap;
            if (i_z == start) { led[(long)(LED_PN+k)*N + c] += a_z * afl;  led[(long)(LED_PA+k)*N + c] += a_z * fabs(afl); }
            if (last)         { led[(long)(LED_PN+k)*N + c] -= a_z * afh;  led[(long)(LED_PA+k)*N + c] += a_z * fabs(afh); }
        }
    }
}
// θ: the LOW face of ring cell j (between j−1 and j), periodic; then the divergence (ADDS).
__global__ void k_face_theta(const double* __restrict__ prim, const World W, long N, double* __restrict__ flux_lo) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    if (!W.act[rz]) return;
    int n = W.nt;
    double pen[6*NP];
    for (int m = -3; m <= 2; m++) {
        int jj = ((j + m) % n + n) % n;
        long cs = (long)jj * NRZ + rz;
        for (int k = 0; k < NP; k++) pen[(m+3)*NP+k] = prim[cs*NP+k];
    }
    double wl[NP], wr[NP], f[NC];
    face_lean(pen, wl, wr);
    e_hllc(wl, wr, 2, f);
    for (int k = 0; k < NC; k++) flux_lo[c*NC+k] = f[k];
}
__global__ void k_div_theta(const double* __restrict__ flux_lo, const World W, long N, double* __restrict__ rate) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z;
    if (!W.act[rz]) return;
    int n = W.nt;
    int jm = (j + n - 1) % n, jp = (j + 1) % n;
    const double* fL = flux_lo + c*NC;
    const double* fR = flux_lo + ((long)jp * NRZ + rz)*NC;
    double ap_lo = W.ap[5*N + ((long)jm * NRZ + rz)];
    double ap_hi = W.ap[5*N + c];
    double inv = area_th(W) / vol_of(W, i_r);
    double kap = W.kappa[c];
    for (int k = 0; k < NC; k++) rate[c*NC+k] += ((fL[k]*ap_lo) - (fR[k]*ap_hi)) * inv / kap;
}

// ---- kernels: prims -------------------------------------------------------------------------
__global__ void k_fill_prims(const double* __restrict__ cons, const World W, const Tables T, long N,
                             const double* __restrict__ hint, int has_hint,
                             double* __restrict__ prim, int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int rz = (int)(c % NRZ);
    if (!W.act[rz]) { for (int k = 0; k < NP; k++) prim[c*NP+k] = 0.0; return; }
    double w[NP];
    double h = has_hint ? hint[c*NP + 4] : 0.0;
    if (!prim_blend(T, cons + c*NC, h, has_hint, w)) { atomicExch(bad, 1); return; }
    for (int k = 0; k < NP; k++) prim[c*NP+k] = w[k];
}

// ---- kernels: the class-A rate + ledger --------------------------------------------------------
// led: column-major [4*NC][N]: 0..NC-1 port_net, NC..2NC-1 port_abs, 2NC.. src_net, 3NC.. src_abs
#define LED_PN 0
#define LED_PA (NC)
#define LED_SN (2*NC)
#define LED_SA (3*NC)
__global__ void k_rate_r(const double* __restrict__ prim, const World W, const Tables T, const StepParams P,
                         long N, double* __restrict__ rate, double* __restrict__ led, int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z, i_z = rz % W.n_z;
    for (int col = 0; col < 4*NC; col++) led[(long)col*N + c] = 0.0;
    if (!W.act[rz]) { for (int k = 0; k < NC; k++) rate[c*NC+k] = 0.0; return; }
    double p[7*NP];
    for (int m = -3; m <= 3; m++)
        if (!gather_rz(W, T, P, prim, j, i_r, i_z, 1, m, p + (m+3)*NP)) { atomicExch(bad, 1); return; }
    double lfl[NP],lfr[NP],rfl[NP],rfr[NP];
    pencil_faces(p, lfl, lfr, rfl, rfr);
    double fL[NC], fR[NC];
    e_hllc(lfl, lfr, 1, fL);
    e_hllc(rfl, rfr, 1, fR);
    int start = W.rs_r[rz], len = W.rl_r[rz];
    double ap_lo = W.ap[0*N + c];
    double ap_hi = (i_r + 1 < start + len) ? W.ap[0*N + (c + W.n_z)] : W.ap[1*N + c];
    double a_lo = area_r(W, i_r) * ap_lo;
    double a_hi = area_r(W, i_r + 1) * ap_hi;
    double kv = W.kappa[c] * vol_of(W, i_r);
    for (int k = 0; k < NC; k++) {
        double afl = a_lo * fL[k], afh = a_hi * fR[k];
        rate[c*NC+k] = (afl - afh) / kv;
        // COUP-2 ports: run-boundary faces only (interior faces telescope).
        if (i_r == start)           { led[(long)(LED_PN+k)*N + c] += afl;  led[(long)(LED_PA+k)*N + c] += fabs(afl); }
        if (i_r + 1 == start + len) { led[(long)(LED_PN+k)*N + c] -= afh;  led[(long)(LED_PA+k)*N + c] += fabs(afh); }
    }
}
__global__ void k_rate_theta(const double* __restrict__ prim, const World W, long N, double* __restrict__ rate) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z;
    if (!W.act[rz]) return;
    int n = W.nt;
    double p[7*NP];
    for (int m = -3; m <= 3; m++) {
        int jj = ((j + m) % n + n) % n;
        long cs = (long)jj * NRZ + rz;
        for (int k = 0; k < NP; k++) p[(m+3)*NP+k] = prim[cs*NP+k];
    }
    double lfl[NP],lfr[NP],rfl[NP],rfr[NP];
    pencil_faces(p, lfl, lfr, rfl, rfr);
    double fL[NC], fR[NC];
    e_hllc(lfl, lfr, 2, fL);
    e_hllc(rfl, rfr, 2, fR);
    int jm = (j + n - 1) % n;
    double ap_lo = W.ap[5*N + ((long)jm * NRZ + rz)];
    double ap_hi = W.ap[5*N + c];
    double inv = area_th(W) / vol_of(W, i_r);
    double kap = W.kappa[c];
    for (int k = 0; k < NC; k++) rate[c*NC+k] += ((fL[k]*ap_lo) - (fR[k]*ap_hi)) * inv / kap;
}
__global__ void k_rate_z(const double* __restrict__ prim, const World W, const Tables T, const StepParams P,
                         long N, double* __restrict__ rate, double* __restrict__ led, int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z, i_z = rz % W.n_z;
    if (!W.act[rz]) return;
    double p[7*NP];
    for (int m = -3; m <= 3; m++)
        if (!gather_rz(W, T, P, prim, j, i_r, i_z, 3, m, p + (m+3)*NP)) { atomicExch(bad, 1); return; }
    double lfl[NP],lfr[NP],rfl[NP],rfr[NP];
    pencil_faces(p, lfl, lfr, rfl, rfr);
    double fL[NC], fR[NC];
    e_hllc(lfl, lfr, 3, fL);
    e_hllc(rfl, rfr, 3, fR);
    int start = W.rs_z[rz], len = W.rl_z[rz];
    double ap_lo = W.ap[2*N + c];
    double ap_hi = (i_z + 1 < start + len) ? W.ap[2*N + (c + 1)] : W.ap[3*N + c];
    double inv_dz = 1.0 / W.dz;
    double kap = W.kappa[c];
    double a_z = area_z(W, i_r);
    for (int k = 0; k < NC; k++) {
        double afl = fL[k] * ap_lo, afh = fR[k] * ap_hi;
        rate[c*NC+k] += (afl - afh) * inv_dz / kap;
        if (i_z == start)           { led[(long)(LED_PN+k)*N + c] += a_z * afl;  led[(long)(LED_PA+k)*N + c] += a_z * fabs(afl); }
        if (i_z + 1 == start + len) { led[(long)(LED_PN+k)*N + c] -= a_z * afh;  led[(long)(LED_PA+k)*N + c] += a_z * fabs(afh); }
    }
}
// Geometric + wall-closure + igniter sources (add_sources) with the ledger.
__global__ void k_sources(const double* __restrict__ prim, const World W, const Igniter G, const StepParams P,
                          long N, double* __restrict__ rate, double* __restrict__ led) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z, i_z = rz % W.n_z;
    if (!W.act[rz]) return;
    double a_in = area_r(W, i_r), a_out = area_r(W, i_r + 1);
    double vol = vol_of(W, i_r);
    double geo = (a_out - a_in) / vol;
    const double* wc = prim + c*NP;
    double rho = wc[0], ur = wc[1], ut = wc[2], pp = wc[4];
    double kappa = W.kappa[c];
    double kv = kappa * vol;
    double s_mr = (a_out*pp - a_in*pp)/vol + rho*ut*ut*geo;
    double s_mt = rho*ur*ut*geo;
    rate[c*NC+I_MR] += s_mr;
    rate[c*NC+I_MT] -= s_mt;
    led[(long)(LED_SN+I_MR)*N + c] += kv * s_mr;  led[(long)(LED_SA+I_MR)*N + c] += fabs(kv * s_mr);
    led[(long)(LED_SN+I_MT)*N + c] -= kv * s_mt;  led[(long)(LED_SA+I_MT)*N + c] += fabs(kv * s_mt);
    if (W.has_geom) {
        double w_r = W.ap[1*N+c]*a_out - W.ap[0*N+c]*a_in - kappa*(a_out - a_in);
        double w_th = (W.ap[5*N+c] - W.ap[4*N+c]) * area_th(W);
        double w_z = (W.ap[3*N+c] - W.ap[2*N+c]) * area_z(W, i_r);
        double inv_kv = 1.0 / (kappa * vol);
        double wr_kv = w_r * inv_kv, wz_kv = w_z * inv_kv;
        rate[c*NC+I_MR] += pp * wr_kv;
        rate[c*NC+I_MZ] += pp * wz_kv;
        led[(long)(LED_SN+I_MR)*N + c] += kv * (pp * wr_kv);  led[(long)(LED_SA+I_MR)*N + c] += fabs(kv * (pp * wr_kv));
        led[(long)(LED_SN+I_MZ)*N + c] += kv * (pp * wz_kv);  led[(long)(LED_SA+I_MZ)*N + c] += fabs(kv * (pp * wz_kv));
        if (w_th != 0.0) {
            double wt_kv = w_th * inv_kv;
            rate[c*NC+I_MT] += pp * wt_kv;
            led[(long)(LED_SN+I_MT)*N + c] += kv * (pp * wt_kv);  led[(long)(LED_SA+I_MT)*N + c] += fabs(kv * (pp * wt_kv));
        }
    }
    // External intake: the igniter deposit (run.rs source_fn), every slot
    // added (the CPU adds src[k] = 0 to the others too).
    double src[NC]; for (int k = 0; k < NC; k++) src[k] = 0.0;
    if (G.armed) {
        double r = r_center(W, i_r), z = z_center(W, i_z), th = theta_center(j, W.nt), t = P.t;
        if (t >= G.t_on && t < G.t_off && r > G.r_lo && r < G.r_hi && z > G.z_lo && z < G.z_hi) {
            int in_theta = 1;
            if (G.theta_gated) {
                double dth = E_TAU / (double)G.nt;
                double jf = floor(th / dth);
                double jm = fmod(jf, (double)G.nt); if (jm < 0.0) jm += (double)G.nt;   // rem_euclid
                in_theta = ((int)jm == G.jt);
            }
            if (in_theta) {
                double ramp = (G.t_ramp > 0.0) ? fmin((t - G.t_on) / G.t_ramp, 1.0) : 1.0;
                src[I_EN] = G.q0 * ramp;
            }
        }
    }
    for (int k = 0; k < NC; k++) {
        rate[c*NC+k] += src[k];
        led[(long)(LED_SN+k)*N + c] += kv * src[k];
        led[(long)(LED_SA+k)*N + c] += fabs(kv * src[k]);
    }
}
// SOLV-4.4 combustion source (accumulate_inner): ADDS to rate[I_RB] + its ledger.
__global__ void k_combustion(const double* __restrict__ prim, const double* __restrict__ cons,
                             const World W, const Tables T, const StepParams P, long N,
                             double* __restrict__ rate, double* __restrict__ led, int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z, i_z = rz % W.n_z;
    if (!W.act[rz]) return;
    double delta = sqrt(W.dr * W.dz);
    const double* w = prim + c*NP;
    double rho = w[0], b = w[6];
    if (!(b >= -EPS_BURN_BOUND && b <= 1.0 + EPS_BURN_BOUND)) { atomicExch(bad, 1); return; }
    double b_c = fmin(fmax(b, 0.0), 1.0);
    Closure cl = closure_at(T, P, w, delta, 0);
    double react = cl.live ? cl.rho_u * cl.k * b_c * (1.0 - b_c) * (b_c - IGN_A) : 0.0;
    double rho_dc_here = cl.live ? rho * cl.d_c : 0.0;
    double vol = vol_of(W, i_r);
    double kv = W.kappa[c] * vol;
    double diff = 0.0;
    // meridional faces r−, r+, z−, z+ (accumulate_inner's order)
    for (int fc = 0; fc < 4; fc++) {
        int d = fc;   // aperture slot: 0 r−, 1 r+, 2 z−, 3 z+
        double ap = W.ap[(long)d*N + c];
        if (ap <= 0.0) continue;
        int nr = i_r, nz = i_z;
        if (fc == 0) { if (i_r == 0) continue; nr = i_r - 1; }
        else if (fc == 1) { nr = i_r + 1; }
        else if (fc == 2) { if (i_z == 0) continue; nz = i_z - 1; }
        else { nz = i_z + 1; }
        if (nr >= W.n_r || nz >= W.n_z) continue;
        int nrz = nr * W.n_z + nz;
        if (!W.act[nrz]) continue;
        long cn = (long)j * NRZ + nrz;
        double b_nbr = cons[cn*NC + I_RB] / cons[cn*NC + I_RHO];
        Closure cn_cl = closure_at(T, P, prim + cn*NP, delta, 0);
        double rho_dc_nbr = cn_cl.live ? prim[cn*NP] * cn_cl.d_c : 0.0;
        double rho_dc_face = 0.5 * (rho_dc_here + rho_dc_nbr);
        double area = (fc < 2) ? area_r(W, fc == 0 ? i_r : i_r + 1) : area_z(W, i_r);
        double dd = (fc < 2) ? W.dr : W.dz;
        diff += ap * area * rho_dc_face * (b_nbr - b) / dd;
    }
    if (W.nt > 1) {
        double a_th = area_th(W);
        double arc = r_center(W, i_r) * (E_TAU / (double)W.nt);
        int jns[2] = { (j + W.nt - 1) % W.nt, (j + 1) % W.nt };
        int dirs[2] = { 4, 5 };   // θ−, θ+
        for (int q = 0; q < 2; q++) {
            double ap_th = W.ap[(long)dirs[q]*N + c];
            if (ap_th <= 0.0) continue;
            long cn = (long)jns[q] * NRZ + rz;
            double b_nbr = cons[cn*NC + I_RB] / cons[cn*NC + I_RHO];
            Closure cn_cl = closure_at(T, P, prim + cn*NP, delta, 0);
            double rho_dc_nbr = cn_cl.live ? prim[cn*NP] * cn_cl.d_c : 0.0;
            double rho_dc_face = 0.5 * (rho_dc_here + rho_dc_nbr);
            diff += ap_th * a_th * rho_dc_face * (b_nbr - b) / arc;
        }
    }
    double diff_div = (kv > 0.0) ? diff / kv : 0.0;
    double src = react + diff_div;
    rate[c*NC + I_RB] += src;
    led[(long)(LED_SN+I_RB)*N + c] += kv * src;
    led[(long)(LED_SA+I_RB)*N + c] += fabs(kv * src);
}

// ---- kernels: class-R -------------------------------------------------------------------------
// Node-0 realizable rate (auto_rate_node0) from the node-0 prims.
__global__ void k_r0(const double* __restrict__ prim, const World W, const Tables T, const StepParams P,
                     long N, double dt, double* __restrict__ r0, int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int rz = (int)(c % NRZ);
    r0[c] = 0.0;
    if (!W.act[rz]) return;
    const double* w = prim + c*NP;
    double delta = sqrt(W.dr * W.dz);
    Closure cl = closure_at(T, P, w, delta, 1);
    if (!cl.live) return;
    double rho = w[0];
    double b_c = fmin(fmax(w[6], 0.0), 1.0);
    double r = rho * (1.0 - b_c) / cl.tau;
    if (!isfinite(r)) { atomicExch(bad, 1); return; }
    if (r <= 0.0 || dt <= 0.0) { r0[c] = r; return; }
    double cap = rho * (1.0 - BURN_COMPLETE);
    double rb = rho * b_c;
    r0[c] = fmin(r, fmax((cap - rb) / dt, 0.0));
}
// The implicit node solve on the composed state (apply_reaction): writes ρb in
// place (one writer per cell), r_trial, and the per-cell (net, gross) records.
__global__ void k_class_r(double* __restrict__ cons, const World W, const Tables T, const StepParams P, long N,
                          const double* __restrict__ r0, const double* __restrict__ r_prev,
                          double wr0, double wrprev, double wrnew,
                          double* __restrict__ r_trial, double* __restrict__ rec_net, double* __restrict__ rec_gross,
                          int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int rz = (int)(c % NRZ);
    int i_r = rz / W.n_z;
    r_trial[c] = 0.0; rec_net[c] = 0.0; rec_gross[c] = 0.0;
    if (!W.act[rz]) return;
    double kappa = W.kappa[c];
    if (kappa <= 0.0) return;
    double* uc = cons + c*NC;
    double u[NC]; for (int k = 0; k < NC; k++) u[k] = uc[k];
    double base_raw = u[I_RB] + wr0 * r0[c] + wrprev * r_prev[c];
    if (!isfinite(base_raw)) { atomicExch(bad, 1); return; }
    double rho = u[I_RHO], inv = 1.0 / rho;
    double cap = rho * (1.0 - BURN_COMPLETE);
    double x, r;
    if (u[I_RB] >= cap) { x = u[I_RB]; r = (u[I_RB] - base_raw) / wrnew; }
    else {
        double base = fmin(fmax(base_raw, 0.0), cap);
        if (base == cap) { x = cap; r = (cap - base_raw) / wrnew; }
        else {
            x = base;
            int done = 0;
            double delta = sqrt(W.dr * W.dz);
            for (int it = 0; it < N_TAU_REFREEZE && !done; it++) {
                double ut[NC]; for (int k = 0; k < NC; k++) ut[k] = u[k];
                ut[I_RB] = x;
                double w[NP];
                if (!prim_blend(T, ut, 0.0, 0, w)) { atomicExch(bad, 1); return; }
                Closure cl = closure_at(T, P, w, delta, 1);
                if (!cl.live) { x = base; done = 2; break; }   // floors: source zero
                x = fmin((base + wrnew * rho / cl.tau) / (1.0 + wrnew / cl.tau), cap);
                if (x == cap) break;
            }
            if (done != 2) x = fmin(fmax(x, 0.0), cap);
            r = (x - base_raw) / wrnew;
        }
    }
    (void)inv;
    uc[I_RB] = x;
    r_trial[c] = r;
    double kv = kappa * vol_of(W, i_r);
    rec_net[c] = kv * (x - u[I_RB]);
    rec_gross[c] = kv * (fabs(wr0 * r0[c]) + fabs(wrprev * r_prev[c]) + fabs(wrnew * r));
}

// ---- kernels: stored totals, SRD, composition, reductions, stable_dt --------------------------
// stored: column-major [2*NC][N]: κV·U_k then κV·|U_k|
__global__ void k_stored(const double* __restrict__ cons, const World W, long N, double* __restrict__ st) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int rz = (int)(c % NRZ);
    int i_r = rz / W.n_z;
    double kv = W.act[rz] ? W.kappa[c] * vol_of(W, i_r) : 0.0;
    for (int k = 0; k < NC; k++) {
        double v = W.act[rz] ? cons[c*NC+k] : 0.0;
        st[(long)k*N + c] = kv * v;
        st[(long)(NC+k)*N + c] = kv * fabs(v);
    }
}
struct SrdTab {
    int ns, na;
    const int* mem_off; const int* mem; const double* mem_kv; const int* cnt;
    const int* aff; const int* aff_owner; const int* inv_off; const int* inv;
};
__global__ void k_srd_q(const double* __restrict__ cons, const SrdTab S, double* __restrict__ q) {
    int s = blockIdx.x*blockDim.x + threadIdx.x;
    if (s >= S.ns) return;
    double num[NC]; for (int k = 0; k < NC; k++) num[k] = 0.0;
    double den = 0.0;
    for (int m = S.mem_off[s]; m < S.mem_off[s+1]; m++) {
        int c = S.mem[m];
        double w = S.mem_kv[m] / (double)S.cnt[c];
        for (int k = 0; k < NC; k++) num[k] += w * cons[(long)c*NC+k];
        den += w;
    }
    for (int k = 0; k < NC; k++) q[(long)s*NC+k] = num[k] / den;
}
__global__ void k_srd_apply(const SrdTab S, const double* __restrict__ q, double* __restrict__ cons) {
    int a = blockIdx.x*blockDim.x + threadIdx.x;
    if (a >= S.na) return;
    int c = S.aff[a];
    double acc[NC];
    for (int k = 0; k < NC; k++) acc[k] = S.aff_owner[a] ? 0.0 : cons[(long)c*NC+k];
    for (int i = S.inv_off[a]; i < S.inv_off[a+1]; i++) {
        int s = S.inv[i];
        for (int k = 0; k < NC; k++) acc[k] += q[(long)s*NC+k];
    }
    double n = (double)S.cnt[c];
    for (int k = 0; k < NC; k++) cons[(long)c*NC+k] = acc[k] / n;
}
__global__ void k_compose_pred(const double* __restrict__ u0, const double* __restrict__ r0, double dt, long n, double* __restrict__ cons) {
    long i = (long)blockIdx.x*blockDim.x + threadIdx.x; if (i >= n) return;
    cons[i] = u0[i] + dt * r0[i];
}
__global__ void k_compose_corr(const double* __restrict__ u0, const double* __restrict__ r0, const double* __restrict__ rl,
                               double half_dt, long n, double* __restrict__ cons) {
    long i = (long)blockIdx.x*blockDim.x + threadIdx.x; if (i >= n) return;
    cons[i] = u0[i] + half_dt * r0[i] + half_dt * rl[i];
}
// Multi-column fixed-topology sum: x column-major [ncol][n]; partial [ncol][nblk]; out[ncol].
__global__ void k_sum_partial(const double* __restrict__ x, long n, double* __restrict__ partial, int nblk) {
    __shared__ double sh[RTPB];
    int col = blockIdx.y;
    long i = (long)blockIdx.x*blockDim.x + threadIdx.x;
    double v = (i < n) ? x[(long)col*n + i] : 0.0;
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int s = RTPB/2; s > 0; s >>= 1) { if (threadIdx.x < s) sh[threadIdx.x] += sh[threadIdx.x+s]; __syncthreads(); }
    if (threadIdx.x == 0) partial[(long)col*nblk + blockIdx.x] = sh[0];
}
__global__ void k_sum_final(const double* __restrict__ partial, int nblk, double* __restrict__ out) {
    __shared__ double sh[RTPB];
    int col = blockIdx.x;
    double v = 0.0;
    for (int j = threadIdx.x; j < nblk; j += RTPB) v += partial[(long)col*nblk + j];
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int s = RTPB/2; s > 0; s >>= 1) { if (threadIdx.x < s) sh[threadIdx.x] += sh[threadIdx.x+s]; __syncthreads(); }
    if (threadIdx.x == 0) out[col] = sh[0];
}
// stable_dt: σ per cell (θ-arc + front carrier + the scale-separation guard) → max-tree.
__global__ void k_sig_partial(const double* __restrict__ prim, const World W, const Tables T, const StepParams P,
                              long N, double* __restrict__ partial, int* __restrict__ bad) {
    __shared__ double sh[RTPB];
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    double v = 0.0;
    if (c < N) {
        long NRZ = (long)W.n_r * W.n_z;
        int rz = (int)(c % NRZ);
        if (W.act[rz]) {
            const double* w = prim + c*NP;
            int i_r = rz / W.n_z;
            double cs = e_sound(w);
            double sig = (fabs(w[1]) + cs)/W.dr + (fabs(w[3]) + cs)/W.dz;
            double dth = E_TAU / (double)W.nt;
            double inv_sq = 1.0/(W.dr*W.dr) + 1.0/(W.dz*W.dz);
            if (W.nt > 1) {
                double arc = r_center(W, i_r) * dth;
                sig += (fabs(w[2]) + cs) / arc;
                inv_sq += 1.0 / (arc * arc);
            }
            double delta = sqrt(W.dr * W.dz);
            Closure cl = closure_at(T, P, w, delta, 0);
            if (cl.live && cl.s_t > 0.0) {
                if (cl.s_t > S_T_MACH_LIMIT * cs) atomicExch(bad, 1);   // scale-separation refusal
                sig += 2.0 * cl.d_c * inv_sq + C_NAGUMO_SLOPE * (cl.rho_u / w[0]) * cl.k;
            }
            v = sig;
            if (!isfinite(v)) atomicExch(bad, 1);
        }
    }
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int s = RTPB/2; s > 0; s >>= 1) { if (threadIdx.x < s) sh[threadIdx.x] = fmax(sh[threadIdx.x], sh[threadIdx.x+s]); __syncthreads(); }
    if (threadIdx.x == 0) partial[blockIdx.x] = sh[0];
}
__global__ void k_max_final(const double* __restrict__ partial, int nblk, double* __restrict__ out) {
    __shared__ double sh[RTPB];
    double v = 0.0;
    for (int j = threadIdx.x; j < nblk; j += RTPB) v = fmax(v, partial[j]);
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int s = RTPB/2; s > 0; s >>= 1) { if (threadIdx.x < s) sh[threadIdx.x] = fmax(sh[threadIdx.x], sh[threadIdx.x+s]); __syncthreads(); }
    if (threadIdx.x == 0) *out = sh[0];
}

// ---- host side: the persistent engine handle ----------------------------------------------------
struct HostCol { const double* p; int np; int lp; const double* h; int nh; int lh; const double* z; int nz; int lz;
                 const int* str; const double* data; int vlog; };
struct HostWorld {
    int n_r, n_z, nt; double r_min, dr, z_min, dz; int has_geom;
    const int* act; const double* kappa; const double* ap;
    const int* rs_r; const int* rl_r; const int* klo_r; const int* khi_r; const double* nlo_r; const double* nhi_r;
    const int* rs_z; const int* rl_z; const int* klo_z; const int* khi_z; const double* nlo_z; const double* nhi_z;
};
struct HostSrd {
    int ns, na; const int* mem_off; const int* mem; const double* mem_kv; const int* cnt;
    const int* aff; const int* aff_owner; const int* inv_off; const int* inv;
};
struct HostTables {
    HostCol urho, usnd, utmp, brho, bsnd, flame, delay;
    double pu_lo, pu_hi, hu_floor, hu_ceil, pb_lo, pb_hi, hb_lo, hb_hi, h_off;
    double p_floor, tu_floor;
};
struct HostIgniter { int armed; double r_lo, r_hi, z_lo, z_hi; int theta_gated; int nt; int jt; double t_on, t_off, t_ramp, q0; };

template <typename U> static U* up(const U* h, long n, void** keep, int* nk) {
    U* d; cudaMalloc(&d, n * sizeof(U));
    cudaMemcpy(d, h, n * sizeof(U), cudaMemcpyHostToDevice);
    keep[(*nk)++] = d; return d;
}
static Col3 upload_col(const HostCol& h, void** keep, int* nk) {
    long nd = (long)h.np * h.nh * h.nz;
    Col3 c;
    c.pp = up(h.p, h.np, keep, nk); c.hp = up(h.h, h.nh, keep, nk); c.zp = up(h.z, h.nz, keep, nk);
    c.np = h.np; c.nh = h.nh; c.nz = h.nz; c.lp = h.lp; c.lh = h.lh; c.lz = h.lz;
    c.s0 = h.str[0]; c.s1 = h.str[1]; c.s2 = h.str[2];
    c.data = up(h.data, nd, keep, nk); c.vlog = h.vlog;
    return c;
}
struct Engine {
    World W; Tables T; SrdTab S; Igniter G;
    double wrinkling, theta, h_total, c_frac;
    long N, NRZ, nscal; int nblk;
    double *cons, *u0, *prim, *prim_dt, *rate_e0, *rate_l, *led0, *ledl, *st, *q;
    double *r0, *r_prev, *r_trial, *rec_net, *rec_gross;
    double *partial, *scal;
    int* bad;
    int primed;
    void* keep[128]; int nk;
    // CRUCIBLE_GPU_PROFILE=1: per-kernel-group wall time (ms) accumulated over
    // the handle's life (cudaEvent pairs; a synchronizing measurement, so it
    // is a PROFILE, never left on in a production march).
    int profile; double prof_ms[10]; long prof_steps;
    int fluxbuf; double *flux_lo, *flux_hi;   // CRUCIBLE_GPU_FLUXBUF=1: the S14 path
};
enum ProfSlot { P_FILL = 0, P_RATE_R, P_RATE_TH, P_RATE_Z, P_SRC, P_COMB, P_CLASS_R, P_SRD, P_REDUCE, P_OTHER };
static const char* PROF_NAME[10] = {"fill_prims", "rate_r", "rate_theta", "rate_z", "sources", "combustion",
                                    "class_r", "srd", "reductions", "compose+copies"};
struct ProfScope {
    Engine* E; int slot; cudaEvent_t a, b;
    ProfScope(Engine* e, int s) : E(e), slot(s) {
        if (E->profile) { cudaEventCreate(&a); cudaEventCreate(&b); cudaEventRecord(a); }
    }
    ~ProfScope() {
        if (E->profile) {
            cudaEventRecord(b); cudaEventSynchronize(b);
            float ms = 0; cudaEventElapsedTime(&ms, a, b); E->prof_ms[slot] += ms;
            cudaEventDestroy(a); cudaEventDestroy(b);
        }
    }
};
static void reduce_cols(Engine* E, const double* x, int ncol, double* out_host) {
    ProfScope ps(E, P_REDUCE);
    dim3 grid((unsigned)E->nblk, (unsigned)ncol);
    k_sum_partial<<<grid, RTPB>>>(x, E->N, E->partial, E->nblk);
    k_sum_final<<<ncol, RTPB>>>(E->partial, E->nblk, E->scal);
    cudaMemcpy(out_host, E->scal, ncol * sizeof(double), cudaMemcpyDeviceToHost);
}
static void launch_srd(Engine* E) {
    if (E->S.ns == 0) return;
    ProfScope ps(E, P_SRD);
    int tpb = 128;
    k_srd_q<<<(E->S.ns + tpb - 1)/tpb, tpb>>>(E->cons, E->S, E->q);
    k_srd_apply<<<(E->S.na + tpb - 1)/tpb, tpb>>>(E->S, E->q, E->cons);
}
// One class-A + combustion RHS evaluation at (cons, t) into `rate`, with the
// ledger per-cell arrays `led` and the prim cache (hinted by the previous fill).
static void eval_rhs(Engine* E, const StepParams& P, double* rate, double* led) {
    int tpb = 128; long cblk = (E->N + tpb - 1)/tpb;
    { ProfScope ps(E, P_FILL);
      k_fill_prims<<<cblk,tpb>>>(E->cons, E->W, E->T, E->N, E->prim, E->primed, E->prim, E->bad); }
    E->primed = 1;
    if (E->fluxbuf) {
        { ProfScope ps(E, P_RATE_R);
          k_face_rz<<<cblk,tpb>>>(E->prim, E->W, E->T, P, E->N, 1, E->flux_lo, E->flux_hi, E->bad);
          k_div_rz<<<cblk,tpb>>>(E->flux_lo, E->flux_hi, E->W, E->N, 1, rate, led); }
        { ProfScope ps(E, P_RATE_TH);
          k_face_theta<<<cblk,tpb>>>(E->prim, E->W, E->N, E->flux_lo);
          k_div_theta<<<cblk,tpb>>>(E->flux_lo, E->W, E->N, rate); }
        { ProfScope ps(E, P_RATE_Z);
          k_face_rz<<<cblk,tpb>>>(E->prim, E->W, E->T, P, E->N, 3, E->flux_lo, E->flux_hi, E->bad);
          k_div_rz<<<cblk,tpb>>>(E->flux_lo, E->flux_hi, E->W, E->N, 3, rate, led); }
    } else {
        { ProfScope ps(E, P_RATE_R);  k_rate_r<<<cblk,tpb>>>(E->prim, E->W, E->T, P, E->N, rate, led, E->bad); }
        { ProfScope ps(E, P_RATE_TH); k_rate_theta<<<cblk,tpb>>>(E->prim, E->W, E->N, rate); }
        { ProfScope ps(E, P_RATE_Z);  k_rate_z<<<cblk,tpb>>>(E->prim, E->W, E->T, P, E->N, rate, led, E->bad); }
    }
    { ProfScope ps(E, P_SRC);     k_sources<<<cblk,tpb>>>(E->prim, E->W, E->G, P, E->N, rate, led); }
    { ProfScope ps(E, P_COMB);    k_combustion<<<cblk,tpb>>>(E->prim, E->cons, E->W, E->T, P, E->N, rate, led, E->bad); }
}

extern "C" void* gpu_engine_create(const HostWorld* hw, const HostSrd* hs, const HostTables* ht,
                                   const HostIgniter* hg, double wrinkling, double theta,
                                   double h_total, double c_frac) {
    Engine* E = (Engine*)calloc(1, sizeof(Engine));
    E->nk = 0;
    long NRZ = (long)hw->n_r * hw->n_z, N = NRZ * hw->nt;
    E->N = N; E->NRZ = NRZ; E->nscal = N * NC; E->nblk = (int)((N + RTPB - 1)/RTPB);
    World& W = E->W;
    W.n_r = hw->n_r; W.n_z = hw->n_z; W.nt = hw->nt; W.r_min = hw->r_min; W.dr = hw->dr; W.z_min = hw->z_min; W.dz = hw->dz;
    W.has_geom = hw->has_geom;
    W.act = up(hw->act, NRZ, E->keep, &E->nk);
    W.kappa = up(hw->kappa, N, E->keep, &E->nk);
    W.ap = up(hw->ap, 6*N, E->keep, &E->nk);
    W.rs_r = up(hw->rs_r, NRZ, E->keep, &E->nk); W.rl_r = up(hw->rl_r, NRZ, E->keep, &E->nk);
    W.klo_r = up(hw->klo_r, NRZ, E->keep, &E->nk); W.khi_r = up(hw->khi_r, NRZ, E->keep, &E->nk);
    W.nlo_r = up(hw->nlo_r, 2*NRZ, E->keep, &E->nk); W.nhi_r = up(hw->nhi_r, 2*NRZ, E->keep, &E->nk);
    W.rs_z = up(hw->rs_z, NRZ, E->keep, &E->nk); W.rl_z = up(hw->rl_z, NRZ, E->keep, &E->nk);
    W.klo_z = up(hw->klo_z, NRZ, E->keep, &E->nk); W.khi_z = up(hw->khi_z, NRZ, E->keep, &E->nk);
    W.nlo_z = up(hw->nlo_z, 2*NRZ, E->keep, &E->nk); W.nhi_z = up(hw->nhi_z, 2*NRZ, E->keep, &E->nk);
    SrdTab& S = E->S; S.ns = hs->ns; S.na = hs->na;
    if (hs->ns > 0) {
        long nmem = hs->mem_off[hs->ns], ninv = hs->inv_off[hs->na];
        S.mem_off = up(hs->mem_off, hs->ns + 1, E->keep, &E->nk);
        S.mem = up(hs->mem, nmem, E->keep, &E->nk);
        S.mem_kv = up(hs->mem_kv, nmem, E->keep, &E->nk);
        S.cnt = up(hs->cnt, N, E->keep, &E->nk);
        S.aff = up(hs->aff, hs->na, E->keep, &E->nk);
        S.aff_owner = up(hs->aff_owner, hs->na, E->keep, &E->nk);
        S.inv_off = up(hs->inv_off, hs->na + 1, E->keep, &E->nk);
        S.inv = up(hs->inv, ninv, E->keep, &E->nk);
    }
    Tables& T = E->T;
    T.urho = upload_col(ht->urho, E->keep, &E->nk); T.usnd = upload_col(ht->usnd, E->keep, &E->nk);
    T.utmp = upload_col(ht->utmp, E->keep, &E->nk); T.brho = upload_col(ht->brho, E->keep, &E->nk);
    T.bsnd = upload_col(ht->bsnd, E->keep, &E->nk); T.flame = upload_col(ht->flame, E->keep, &E->nk);
    T.delay = upload_col(ht->delay, E->keep, &E->nk);
    T.v.pu_lo = ht->pu_lo; T.v.pu_hi = ht->pu_hi; T.v.hu_floor = ht->hu_floor; T.v.hu_ceil = ht->hu_ceil;
    T.v.pb_lo = ht->pb_lo; T.v.pb_hi = ht->pb_hi; T.v.hb_lo = ht->hb_lo; T.v.hb_hi = ht->hb_hi; T.v.h_off = ht->h_off;
    T.p_floor = ht->p_floor; T.tu_floor = ht->tu_floor;
    Igniter& G = E->G;
    G.armed = hg->armed; G.r_lo = hg->r_lo; G.r_hi = hg->r_hi; G.z_lo = hg->z_lo; G.z_hi = hg->z_hi;
    G.theta_gated = hg->theta_gated; G.nt = hg->nt; G.jt = hg->jt; G.t_on = hg->t_on; G.t_off = hg->t_off;
    G.t_ramp = hg->t_ramp; G.q0 = hg->q0;
    E->wrinkling = wrinkling; E->theta = theta; E->h_total = h_total; E->c_frac = c_frac;
    size_t nb = (size_t)E->nscal * sizeof(double);
    cudaMalloc(&E->cons, nb); cudaMalloc(&E->u0, nb);
    cudaMalloc(&E->prim, N*NP*sizeof(double)); cudaMalloc(&E->prim_dt, N*NP*sizeof(double));
    cudaMalloc(&E->rate_e0, nb); cudaMalloc(&E->rate_l, nb);
    cudaMalloc(&E->led0, (size_t)4*NC*N*sizeof(double)); cudaMalloc(&E->ledl, (size_t)4*NC*N*sizeof(double));
    cudaMalloc(&E->st, (size_t)2*NC*N*sizeof(double));
    cudaMalloc(&E->q, (size_t)(S.ns > 0 ? S.ns : 1)*NC*sizeof(double));
    cudaMalloc(&E->r0, N*sizeof(double)); cudaMalloc(&E->r_prev, N*sizeof(double)); cudaMalloc(&E->r_trial, N*sizeof(double));
    cudaMalloc(&E->rec_net, N*sizeof(double)); cudaMalloc(&E->rec_gross, N*sizeof(double));
    cudaMalloc(&E->partial, (size_t)4*NC*E->nblk*sizeof(double));
    cudaMalloc(&E->scal, (size_t)4*NC*sizeof(double));
    cudaMalloc(&E->bad, sizeof(int));
    cudaMemset(E->bad, 0, sizeof(int));
    cudaMemset(E->prim, 0, N*NP*sizeof(double));
    E->primed = 0;
    const char* pf = getenv("CRUCIBLE_GPU_PROFILE");
    E->profile = (pf && pf[0] == '1') ? 1 : 0;
    const char* fb = getenv("CRUCIBLE_GPU_FLUXBUF");
    E->fluxbuf = (fb && fb[0] == '1') ? 1 : 0;
    cudaMalloc(&E->flux_lo, nb); cudaMalloc(&E->flux_hi, nb);
    cudaMemset(E->flux_lo, 0, nb); cudaMemset(E->flux_hi, 0, nb);
    for (int i = 0; i < 10; i++) E->prof_ms[i] = 0.0;
    E->prof_steps = 0;
    return E;
}
// Dump the profile (stderr) — called by destroy when profiling.
static void prof_dump(const Engine* E) {
    if (!E->profile) return;
    double tot = 0.0; for (int i = 0; i < 10; i++) tot += E->prof_ms[i];
    fprintf(stderr, "[gpu profile] %ld steps, %.1f ms total in kernels (%.4f s/step)\n",
            E->prof_steps, tot, tot / 1000.0 / (E->prof_steps > 0 ? E->prof_steps : 1));
    for (int i = 0; i < 10; i++)
        fprintf(stderr, "[gpu profile]   %-16s %9.1f ms  %5.1f%%\n", PROF_NAME[i], E->prof_ms[i],
                tot > 0 ? 100.0 * E->prof_ms[i] / tot : 0.0);
}
extern "C" void gpu_engine_destroy(void* h) {
    Engine* E = (Engine*)h;
    prof_dump(E);
    cudaFree(E->cons); cudaFree(E->u0); cudaFree(E->prim); cudaFree(E->prim_dt); cudaFree(E->rate_e0); cudaFree(E->rate_l);
    cudaFree(E->led0); cudaFree(E->ledl); cudaFree(E->st); cudaFree(E->q); cudaFree(E->r0); cudaFree(E->r_prev);
    cudaFree(E->r_trial); cudaFree(E->rec_net); cudaFree(E->rec_gross); cudaFree(E->partial); cudaFree(E->scal); cudaFree(E->bad);
    cudaFree(E->flux_lo); cudaFree(E->flux_hi);
    for (int i = 0; i < E->nk; i++) cudaFree(E->keep[i]);
    free(E);
}
// State I/O (the FND-6 checkpoint round-trip): cons and the prim cache (the
// warm-start hints — restored so a resumed march is bit-faithful).
extern "C" void gpu_engine_upload(void* h, const double* cons, const double* prim_hint, int primed) {
    Engine* E = (Engine*)h;
    cudaMemcpy(E->cons, cons, (size_t)E->nscal*sizeof(double), cudaMemcpyHostToDevice);
    if (primed && prim_hint) {
        cudaMemcpy(E->prim, prim_hint, (size_t)E->N*NP*sizeof(double), cudaMemcpyHostToDevice);
        E->primed = 1;
    } else {
        cudaMemset(E->prim, 0, (size_t)E->N*NP*sizeof(double));
        E->primed = 0;
    }
    cudaMemset(E->bad, 0, sizeof(int));
}
extern "C" int gpu_engine_download(void* h, double* cons, double* prim_hint) {
    Engine* E = (Engine*)h;
    cudaDeviceSynchronize();
    cudaMemcpy(cons, E->cons, (size_t)E->nscal*sizeof(double), cudaMemcpyDeviceToHost);
    if (prim_hint) cudaMemcpy(prim_hint, E->prim, (size_t)E->N*NP*sizeof(double), cudaMemcpyDeviceToHost);
    return E->primed;
}
// stable_dt at the current state (warm hints from the prim cache; the cache
// itself is NOT touched — Euler::stable_dt_ws reads hints only). Returns Δt
// or a negative sentinel; *bad set on a refusal.
extern "C" double gpu_engine_stable_dt(void* h, double t, double cfl, int* bad) {
    Engine* E = (Engine*)h;
    StepParams P; P.t = t; P.mdot_per_area = 0.0; P.h_total = E->h_total; P.c_frac = E->c_frac; P.p_amb = 0.0;
    P.wrinkling = E->wrinkling; P.theta = E->theta;
    int tpb = 128; long cblk = (E->N + tpb - 1)/tpb;
    cudaMemset(E->bad, 0, sizeof(int));
    { ProfScope ps(E, P_FILL); k_fill_prims<<<cblk,tpb>>>(E->cons, E->W, E->T, E->N, E->prim, E->primed, E->prim_dt, E->bad); }
    { ProfScope ps(E, P_REDUCE);
      k_sig_partial<<<E->nblk,RTPB>>>(E->prim_dt, E->W, E->T, P, E->N, E->partial, E->bad);
      k_max_final<<<1,RTPB>>>(E->partial, E->nblk, E->scal); }
    cudaDeviceSynchronize();
    double max_sig = 0.0; cudaMemcpy(&max_sig, E->scal, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(bad, E->bad, sizeof(int), cudaMemcpyDeviceToHost);
    if (*bad || !(max_sig > 0.0)) return -1.0;
    return cfl / max_sig;
}
// One full SDC step (flow + reaction, stagewise SRD). Writes the audit operands:
// out[0..7) stored_before, [7..14) stored_before_abs, [14..21) after, [21..28) after_abs,
// [28..56) l0 (pn, pa, sn, sa × NC), [56..84) l_last, [84] burn_applied, [85] burn_gross.
// Returns the halt flag (0 = ok).
extern "C" int gpu_engine_step(void* h, double t, double dt, double mdot_per_area, double p_amb_t,
                               double p_amb_t1, double* out) {
    Engine* E = (Engine*)h;
    int tpb = 128; long cblk = (E->N + tpb - 1)/tpb, sblk = (E->nscal + tpb - 1)/tpb;
    StepParams P0; P0.t = t; P0.mdot_per_area = mdot_per_area; P0.h_total = E->h_total; P0.c_frac = E->c_frac;
    P0.p_amb = p_amb_t; P0.wrinkling = E->wrinkling; P0.theta = E->theta;
    StepParams P1 = P0; P1.t = t + dt; P1.p_amb = p_amb_t1;
    cudaMemset(E->bad, 0, sizeof(int));
    // stored BEFORE
    k_stored<<<cblk,tpb>>>(E->cons, E->W, E->N, E->st);
    reduce_cols(E, E->st, 2*NC, out + 0);
    // node 0
    { ProfScope ps(E, P_OTHER); cudaMemcpy(E->u0, E->cons, (size_t)E->nscal*sizeof(double), cudaMemcpyDeviceToDevice); }
    eval_rhs(E, P0, E->rate_e0, E->led0);
    reduce_cols(E, E->led0, 4*NC, out + 28);
    { ProfScope ps(E, P_CLASS_R); k_r0<<<cblk,tpb>>>(E->prim, E->W, E->T, P0, E->N, dt, E->r0, E->bad); }
    { ProfScope ps(E, P_OTHER); cudaMemcpy(E->r_prev, E->r0, (size_t)E->N*sizeof(double), cudaMemcpyDeviceToDevice); }
    double burn_applied = 0.0, burn_gross = 0.0;
    double half = 0.5 * dt;
    for (int sweep = 0; sweep <= 2; sweep++) {
        int predictor = (sweep == 0);
        double wq0 = predictor ? 0.0 : half, wqprev = predictor ? 0.0 : -half, wqnew = dt;
        if (!predictor) {
            eval_rhs(E, P1, E->rate_l, E->ledl);
            ProfScope ps(E, P_OTHER);
            k_compose_corr<<<sblk,tpb>>>(E->u0, E->rate_e0, E->rate_l, half, E->nscal, E->cons);
        } else {
            ProfScope ps(E, P_OTHER);
            k_compose_pred<<<sblk,tpb>>>(E->u0, E->rate_e0, dt, E->nscal, E->cons);
        }
        launch_srd(E);
        { ProfScope ps(E, P_CLASS_R);
          k_class_r<<<cblk,tpb>>>(E->cons, E->W, E->T, P1, E->N, E->r0, E->r_prev, wq0, wqprev, wqnew,
                                  E->r_trial, E->rec_net, E->rec_gross, E->bad); }
        double rec[2];
        // two single-column reductions (net, gross)
        reduce_cols(E, E->rec_net, 1, rec + 0);
        reduce_cols(E, E->rec_gross, 1, rec + 1);
        burn_applied = rec[0]; burn_gross = rec[1];
        if (sweep < 2) { double* tmp = E->r_prev; E->r_prev = E->r_trial; E->r_trial = tmp; }
    }
    reduce_cols(E, E->ledl, 4*NC, out + 56);
    k_stored<<<cblk,tpb>>>(E->cons, E->W, E->N, E->st);
    reduce_cols(E, E->st, 2*NC, out + 14);
    out[84] = burn_applied; out[85] = burn_gross;
    E->prof_steps += 1;
    cudaDeviceSynchronize();
    int bad = 0; cudaMemcpy(&bad, E->bad, sizeof(int), cudaMemcpyDeviceToHost);
    return bad;
}
// Select the sweep path at runtime (0 = fused kernels, 1 = the S14 flux buffer).
extern "C" void gpu_engine_set_fluxbuf(void* h, int flag) { ((Engine*)h)->fluxbuf = flag ? 1 : 0; }
// Single-shot class-A + combustion RHS at the current state (cross-check oracle);
// rate dense [N*NC]; returns the halt flag.
extern "C" int gpu_engine_rhs(void* h, double t, double mdot_per_area, double p_amb, double* rate, double* led_out) {
    Engine* E = (Engine*)h;
    StepParams P; P.t = t; P.mdot_per_area = mdot_per_area; P.h_total = E->h_total; P.c_frac = E->c_frac;
    P.p_amb = p_amb; P.wrinkling = E->wrinkling; P.theta = E->theta;
    cudaMemset(E->bad, 0, sizeof(int));
    eval_rhs(E, P, E->rate_e0, E->led0);
    reduce_cols(E, E->led0, 4*NC, led_out);
    cudaDeviceSynchronize();
    cudaMemcpy(rate, E->rate_e0, (size_t)E->nscal*sizeof(double), cudaMemcpyDeviceToHost);
    int bad = 0; cudaMemcpy(&bad, E->bad, sizeof(int), cudaMemcpyDeviceToHost);
    return bad;
}
