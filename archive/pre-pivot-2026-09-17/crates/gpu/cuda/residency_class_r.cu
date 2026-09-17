// crucible-gpu — S13c residency: the CLASS-R implicit AUTO-IGNITION node solve
// on device (SOLV-4 §3.6 v0.4.4/0.4.8; COUP-3 §3.3's cell-local stiff slot).
// Bit-for-formula from crucible_solvers::euler::Combustion::implicit_auto_update
// (combustion.rs): backward-Euler in ρb at frozen τ_ign, a fixed
// N_TAU_REFREEZE(=2) re-evaluation of τ at the current iterate, the burnt
// fixed point ρ(1−BURN_COMPLETE) honored EXACTLY (a stiff update parks AT the
// cap — the exact integral of the declared discontinuous law), the symmetric
// base projection onto [0, cap], and the cold-side non-reactive floors.
//
// THE NEW SURFACE vs the S27 legs: the per-refreeze state query is the blend's
// FULL `prim_checked` — a three-way branch select on b = x/ρ:
//   b ≤ EPS_B_PURE_UNBURNT          → the unburnt TableEos projection (cold)
//   b ≥ 1 − EPS_B_PURE_BURNT        → the burnt TableEos projection at e + h_off
//   otherwise                       → the mid-b two-branch blend projection
// (blend_eos.rs project_pressure). A stiff node walks x: base → mid-b → cap,
// so one cell crosses all three within a solve; the select is exactly the
// CPU's need_u/need_b logic, no clamp, no smoothing. Both projections are the
// in-tree device pieces (residency_eos.cu / residency_blend_eos.cu) inlined
// with their OWN bracket margins (TableEos 1e-12, blend 1e-9 — table_eos.rs vs
// blend_eos.rs; a latent 1e-12 in the S27 blend port is corrected there too).
//
// Two-root → CONTINUITY (Ben ruling 2026-08-31; SOLV-4 0.4.9): the first scan
// crossing from the cold end, both projections. The near-vacuum golden-section
// tangency corner stays the shared follow-on (NaN sentinel → the node result
// is NaN, never a silent branch; the interior fixture never reaches it).
//
// Determinism (META-1 §2.5): a per-cell gather, no reduction — same-build
// reruns bit-identical; CPU↔GPU is FMA-order + libm(ln/exp/pow) ECT.
#include <cuda_runtime.h>
#include <math.h>

#define EPS_P_PROJECTION 1e-11
#define N_P_ITER_MAX 48
#define N_P_BISECT 64
#define N_P_SCAN 64
#define H_MARGIN_TABLE 1e-12     // table_eos.rs H_BRACKET_MARGIN
#define H_MARGIN_BLEND 1e-9      // blend_eos.rs H_BRACKET_MARGIN
#define EPS_B_PURE_UNBURNT 1.0e-9
#define BURN_COMPLETE 1.0e-3
#define EPS_B_PURE_BURNT (2.0 * BURN_COMPLETE)
#define N_TAU_REFREEZE 2
#define NC 7

struct Col3 {
    const double *pp, *hp, *zp; int np, nh, nz; int lp, lh, lz;
    int s0, s1, s2; const double* data; int vlog;
};
__device__ __forceinline__ int r_cell(const double* pts, int n, double q) {
    int cnt = 0; for (int i = 0; i < n; i++) cnt += (pts[i] <= q);
    int hi = cnt < n-1 ? cnt : n-1; return (hi > 1 ? hi : 1) - 1;
}
__device__ __forceinline__ double r_frac(const double* pts, int i, int is_log, double q) {
    if (is_log) return (log(q) - log(pts[i])) / (log(pts[i+1]) - log(pts[i]));
    return (q - pts[i]) / (pts[i+1] - pts[i]);
}
// Fixed 8-corner multilinear in interp_rule space (BoundColumn::interpolate order).
__device__ double r_interp(const Col3& c, double qp, double qh, double qz) {
    int ip = r_cell(c.pp, c.np, qp), ih = r_cell(c.hp, c.nh, qh), iz = r_cell(c.zp, c.nz, qz);
    double t0 = r_frac(c.pp, ip, c.lp, qp), t1 = r_frac(c.hp, ih, c.lh, qh), t2 = r_frac(c.zp, iz, c.lz, qz);
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

// ---- the blend's partition + closure (blend_eos.rs) ----------------------
struct BlendEnv {
    double pu_lo, pu_hi, hu_floor, hu_ceil;   // unburnt (p,h) envelope
    double pb_lo, pb_hi, hb_lo, hb_hi;        // burnt (p,h) envelope (RAW, un-shifted)
    double h_off;
};
// partition_h → (h_u, h_b); hb window = burnt env − h_off.
__device__ __forceinline__ void r_partition(const BlendEnv& v, double h, double b, double* h_u, double* h_b) {
    if ((h >= v.hu_floor && h <= v.hu_ceil) || b <= EPS_B_PURE_UNBURNT) { *h_u = h; *h_b = h; return; }
    double hu_pin = fmin(fmax(h, v.hu_floor), v.hu_ceil);
    double balance = (h - (1.0 - b) * hu_pin) / b;
    *h_u = hu_pin; *h_b = fmin(fmax(balance, v.hb_lo - v.h_off), v.hb_hi - v.h_off);
}
__device__ double r_blend_inv_rho(const BlendEnv& v, const Col3& urho, const Col3& brho,
                                  double p, double h, double z, double b) {
    double h_u, h_b; r_partition(v, h, b, &h_u, &h_b);
    double vv = 0.0;
    if (b < 1.0 - EPS_B_PURE_BURNT) { double ru = r_interp(urho, p, h_u, z); vv += (1.0 - b) / ru; }
    if (b > EPS_B_PURE_UNBURNT)     { double rb = r_interp(brho, p, h_b + v.h_off, z); vv += b / rb; }
    return vv;
}

// ---- the residual functors: g(p) for the table path and the blend path ----
struct GTable {   // ρ_tab(p, e_q + p/ρ, z) − ρ
    const Col3* rc; double e_q, inv, z, rho;
    __device__ double operator()(double p) const { return r_interp(*rc, p, e_q + p * inv, z) - rho; }
};
struct GBlend {   // v_blend(p, e + p/ρ, z, b) − 1/ρ
    const BlendEnv* v; const Col3* urho; const Col3* brho; double e, inv, z, b;
    __device__ double operator()(double p) const {
        return r_blend_inv_rho(*v, *urho, *brho, p, e + p * inv, z, b) - inv;
    }
};

// Deterministic Illinois regula-falsi + bisection backstop (TableEos::illinois_root).
template <typename G>
__device__ double r_illinois(double a, double b, double ga, double gb, const G& g) {
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

// The fixed log-scan for the FIRST sign change (continuity) → Illinois; NaN on
// no crossing (tangency corner, follow-on).
template <typename G>
__device__ double r_scan_first(double lo, double hi, double ga0, const G& g) {
    double ratio = hi / lo, prev_p = lo, prev_g = ga0;
    for (int k = 1; k <= N_P_SCAN; k++) {
        double pk = lo * pow(ratio, (double)k / (double)N_P_SCAN);
        double gk = g(pk);
        if (gk == 0.0) return pk;
        if (prev_g * gk < 0.0) return r_illinois(prev_p, pk, prev_g, gk, g);
        prev_p = pk; prev_g = gk;
    }
    return nan("");
}

// TableEos::project_pressure (cold): admissible bracket → straddle → Illinois,
// else the first-crossing scan. Margin 1e-12.
__device__ double r_project_table(const Col3& rc, double rho, double e_q, double z, double inv,
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
    if (ga0 * gb0 < 0.0) return r_illinois(lo, hi, ga0, gb0, g);
    return r_scan_first(lo, hi, ga0, g);
}

// BurnBlendEos::project_pressure, mid-b: products h-window ∩ p-window bracket,
// the scan ALWAYS first (first crossing → Illinois). Margin 1e-9.
__device__ double r_project_mid(const BlendEnv& v, const Col3& urho, const Col3& brho,
                                double rho, double e, double z, double b, double inv) {
    double h_lo = v.hb_lo - v.h_off, h_hi = v.hb_hi - v.h_off;
    double p_lo_env = fmax(v.pu_lo, v.pb_lo), p_hi_env = fmin(v.pu_hi, v.pb_hi);
    double p_from_h_lo = rho * (h_lo - e), p_from_h_hi = rho * (h_hi - e);
    double lo = fmax(fmax(p_lo_env, p_from_h_lo * (1.0 + H_MARGIN_BLEND)), 0.0);
    double hi = fmin(p_hi_env, p_from_h_hi * (1.0 - H_MARGIN_BLEND));
    if (!(isfinite(lo) && isfinite(hi)) || lo >= hi || hi <= 0.0) return nan("");
    lo = fmax(lo, 0.0);
    GBlend g{&v, &urho, &brho, e, inv, z, b};
    double ga0 = g(lo), gb0 = g(hi);
    if (ga0 == 0.0) return lo;
    if (gb0 == 0.0) return hi;
    return r_scan_first(lo, hi, ga0, g);
}

// The blend's prim_checked pressure: the THREE-WAY branch select on b.
__device__ __forceinline__ double r_blend_pressure(const BlendEnv& v, const Col3& urho, const Col3& brho,
                                                   double rho, double e, double z, double b, double inv) {
    bool need_u = b < 1.0 - EPS_B_PURE_BURNT;
    bool need_b = b > EPS_B_PURE_UNBURNT;
    if (!need_b) return r_project_table(urho, rho, e, z, inv, v.pu_lo, v.pu_hi, v.hu_floor, v.hu_ceil);
    if (!need_u) return r_project_table(brho, rho, e + v.h_off, z, inv, v.pb_lo, v.pb_hi, v.hb_lo, v.hb_hi);
    return r_project_mid(v, urho, brho, rho, e, z, b, inv);
}

// One class-R node solve per cell (implicit_auto_update). u = the composed
// conserved state (7 slots, interleaved), base = ρb + quadrature terms.
__global__ void k_class_r(const double* __restrict__ u, const double* __restrict__ base_in,
                          int n, double w_new, Col3 urho, Col3 utemp, Col3 brho, Col3 dly,
                          BlendEnv v, double p_floor, double tu_floor,
                          double* __restrict__ x_out, double* __restrict__ rate_out) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double* uc = u + i * NC;
    double rho = uc[0], inv = 1.0 / rho;
    double ur = uc[1] * inv, ut = uc[2] * inv, uz = uc[3] * inv;
    double e = uc[4] * inv - 0.5 * (ur * ur + ut * ut + uz * uz);
    double z = uc[5] * inv;
    double rb_adv = uc[6];
    double base_raw = base_in[i];
    double cap = rho * (1.0 - BURN_COMPLETE);
    // Advected ρb at/past the cap: zero source, state unchanged (never shaved).
    if (rb_adv >= cap) { x_out[i] = rb_adv; rate_out[i] = (rb_adv - base_raw) / w_new; return; }
    // Symmetric base projection onto the law's invariant set [0, cap].
    double base = fmin(fmax(base_raw, 0.0), cap);
    if (base == cap) { x_out[i] = cap; rate_out[i] = (cap - base_raw) / w_new; return; }
    double x = base;
    for (int it = 0; it < N_TAU_REFREEZE; it++) {
        double b_raw = x * inv;
        double b = fmin(fmax(b_raw, 0.0), 1.0);
        double p = r_blend_pressure(v, urho, brho, rho, e, z, b, inv);
        if (!isfinite(p)) { x_out[i] = nan(""); rate_out[i] = nan(""); return; }   // tangency corner sentinel
        double h = e + p * inv;
        // below_unburnt_floor: colder than any representable reactant → source zero.
        if (h < v.hu_floor) { x_out[i] = base; rate_out[i] = (base - base_raw) / w_new; return; }
        // unburnt_temperature: the pinned reactant sub-state's T_u.
        double h_u, h_b; r_partition(v, h, b, &h_u, &h_b);
        double t_u = r_interp(utemp, p, h_u, z);
        if (p < p_floor || t_u < tu_floor) { x_out[i] = base; rate_out[i] = (base - base_raw) / w_new; return; }
        double tau = r_interp(dly, p, t_u, z);
        // BE at frozen τ, parked at the cap.
        x = fmin((base + w_new * rho / tau) / (1.0 + w_new / tau), cap);
        if (x == cap) break;
    }
    x = fmin(fmax(x, 0.0), cap);
    x_out[i] = x; rate_out[i] = (x - base_raw) / w_new;
}

static Col3 upload_col(const double* p, int np, int lp, const double* h, int nh, int lh,
                       const double* z, int nz, int lz, const int* str, const double* data, int vlog,
                       double** keep, int* nkeep) {
    double *d_p, *d_h, *d_z, *d_d; long nd = (long)np * nh * nz;
    cudaMalloc(&d_p, np * sizeof(double)); cudaMalloc(&d_h, nh * sizeof(double));
    cudaMalloc(&d_z, nz * sizeof(double)); cudaMalloc(&d_d, nd * sizeof(double));
    cudaMemcpy(d_p, p, np * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_h, h, nh * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_z, z, nz * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_d, data, nd * sizeof(double), cudaMemcpyHostToDevice);
    keep[(*nkeep)++] = d_p; keep[(*nkeep)++] = d_h; keep[(*nkeep)++] = d_z; keep[(*nkeep)++] = d_d;
    Col3 c{d_p, d_h, d_z, np, nh, nz, lp, lh, lz, str[0], str[1], str[2], d_d, vlog};
    return c;
}

// FFI: the class-R node solve over N cells. Four columns: unburnt ρ + T
// (shared unburnt axes), burnt ρ (own axes), τ_ign (ignition (p,T_u,Z) axes).
extern "C" void gpu_class_r_update(
    const double* u, const double* base, int n, double w_new,
    // unburnt axes + rho + temp data
    const double* up, int unp, int ulp, const double* uh, int unh, int ulh,
    const double* uz, int unz, int ulz, const int* ustr,
    const double* urho_data, int urho_vlog, const double* utemp_data, int utemp_vlog,
    // burnt axes + rho data
    const double* bp, int bnp, int blp, const double* bh, int bnh, int blh,
    const double* bz, int bnz, int blz, const int* bstr, const double* brho_data, int brho_vlog,
    // ignition axes + delay data
    const double* ip, int inp, int ilp, const double* it, int inh, int ilh,
    const double* iz, int inz, int ilz, const int* istr, const double* dly_data, int dly_vlog,
    // envelopes + scalars
    double pu_lo, double pu_hi, double hu_floor, double hu_ceil,
    double pb_lo, double pb_hi, double hb_lo, double hb_hi, double h_off,
    double p_floor, double tu_floor,
    double* x_out, double* rate_out) {
    size_t nb = (size_t)n * sizeof(double);
    double *d_u, *d_base, *d_x, *d_rate;
    cudaMalloc(&d_u, nb * NC); cudaMalloc(&d_base, nb); cudaMalloc(&d_x, nb); cudaMalloc(&d_rate, nb);
    cudaMemcpy(d_u, u, nb * NC, cudaMemcpyHostToDevice);
    cudaMemcpy(d_base, base, nb, cudaMemcpyHostToDevice);
    double* keep[16]; int nkeep = 0;
    Col3 urho = upload_col(up, unp, ulp, uh, unh, ulh, uz, unz, ulz, ustr, urho_data, urho_vlog, keep, &nkeep);
    Col3 utemp = upload_col(up, unp, ulp, uh, unh, ulh, uz, unz, ulz, ustr, utemp_data, utemp_vlog, keep, &nkeep);
    Col3 brho = upload_col(bp, bnp, blp, bh, bnh, blh, bz, bnz, blz, bstr, brho_data, brho_vlog, keep, &nkeep);
    Col3 dly = upload_col(ip, inp, ilp, it, inh, ilh, iz, inz, ilz, istr, dly_data, dly_vlog, keep, &nkeep);
    BlendEnv v{pu_lo, pu_hi, hu_floor, hu_ceil, pb_lo, pb_hi, hb_lo, hb_hi, h_off};
    int tpb = 64; long blk = (n + tpb - 1) / tpb;
    k_class_r<<<blk, tpb>>>(d_u, d_base, n, w_new, urho, utemp, brho, dly, v, p_floor, tu_floor, d_x, d_rate);
    cudaDeviceSynchronize();
    cudaMemcpy(x_out, d_x, nb, cudaMemcpyDeviceToHost);
    cudaMemcpy(rate_out, d_rate, nb, cudaMemcpyDeviceToHost);
    cudaFree(d_u); cudaFree(d_base); cudaFree(d_x); cudaFree(d_rate);
    for (int k = 0; k < nkeep; k++) cudaFree(keep[k]);
}
