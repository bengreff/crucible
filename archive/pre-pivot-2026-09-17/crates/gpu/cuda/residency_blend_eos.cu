// crucible-gpu — S13c residency: the BLEND EOS (p,h,Z) projection on device —
// the two-branch (unburnt ⊕ burnt) shifting-equilibrium projection that a
// burn-progress run's conserved→primitive conversion performs, and the class-R
// auto-ignition node solve's per-refreeze state query. Bit-for-formula from
// crucible_solvers::euler::BurnBlendEos::{project_pressure, blend_inv_rho}
// (blend_eos.rs): the mass-weighted specific-volume density closure
// 1/ρ = (1−b)/ρ_u(p,h_u,Z) + b/ρ_b(p,h_b+off,Z) on the partitioned sub-state
// enthalpies, with the same deterministic Illinois + first-crossing scan the
// TableEos port uses. The multi-root case resolves by CONTINUITY (first
// crossing) — Ben ruling 2026-08-31 — NOT the CPU's >1-crossing refusal (the
// unique-root fixture matches the CPU bit-for-formula; the device carries the
// continuity policy the ruling prescribes).
//
// SCOPE: the MID-b path (both branches live) — the genuinely new blend physics;
// the pure limits (b≈0 / b≈1) delegate to the single-branch TableEos projection
// (already resident, `residency_eos.cu`). The near-vacuum golden-section
// tangency corner is the shared follow-on. Reuses the interp + Illinois pattern.
#include <cuda_runtime.h>
#include <math.h>

#define EPS_P_PROJECTION 1e-11
// blend_eos.rs's OWN margin (1e-9, wider than TableEos's 1e-12 — the partition
// arithmetic's ulp slack); S27 shipped 1e-12 here, latent because the fixture's
// bracket was envelope-dominated. Corrected S28 (the class-R port inlines the
// same projection with the same 1e-9).
#define H_BRACKET_MARGIN 1e-9
#define N_P_ITER_MAX 48
#define N_P_BISECT 64
#define N_P_SCAN 64
#define EPS_B_PURE_UNBURNT 1.0e-9
#define BURN_COMPLETE 1.0e-3
#define EPS_B_PURE_BURNT (2.0 * BURN_COMPLETE)

struct Col3 {
    const double *pp, *hp, *zp; int np, nh, nz; int lp, lh, lz;
    int s0, s1, s2; const double* data; int vlog;
};
__device__ __forceinline__ int b_cell(const double* pts, int n, double q) {
    int cnt = 0; for (int i = 0; i < n; i++) cnt += (pts[i] <= q);
    int hi = cnt < n-1 ? cnt : n-1; return (hi > 1 ? hi : 1) - 1;
}
__device__ __forceinline__ double b_frac(const double* pts, int i, int is_log, double q) {
    if (is_log) return (log(q) - log(pts[i])) / (log(pts[i+1]) - log(pts[i]));
    return (q - pts[i]) / (pts[i+1] - pts[i]);
}
__device__ double b_interp(const Col3& c, double qp, double qh, double qz) {
    int ip = b_cell(c.pp, c.np, qp), ih = b_cell(c.hp, c.nh, qh), iz = b_cell(c.zp, c.nz, qz);
    double t0 = b_frac(c.pp, ip, c.lp, qp), t1 = b_frac(c.hp, ih, c.lh, qh), t2 = b_frac(c.zp, iz, c.lz, qz);
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

// partition_h (blend_eos.rs) → (h_u, h_b). hu_floor/ceil = unburnt h-envelope;
// hb_lo/hi = burnt h-envelope shifted by −h_offset.
__device__ void b_partition(double h, double b, double hu_floor, double hu_ceil,
                            double hb_lo, double hb_hi, double* h_u, double* h_b) {
    if ((h >= hu_floor && h <= hu_ceil) || b <= EPS_B_PURE_UNBURNT) { *h_u = h; *h_b = h; return; }
    double hu_pin = fmin(fmax(h, hu_floor), hu_ceil);
    double balance = (h - (1.0 - b) * hu_pin) / b;
    *h_u = hu_pin; *h_b = fmin(fmax(balance, hb_lo), hb_hi);
}

// blend_inv_rho: mass-weighted specific volume on the partitioned sub-states.
__device__ double b_inv_rho(double p, double h, double z, double b,
                            const Col3& urho, const Col3& brho,
                            double hu_floor, double hu_ceil, double hb_lo, double hb_hi, double h_off) {
    double h_u, h_b; b_partition(h, b, hu_floor, hu_ceil, hb_lo, hb_hi, &h_u, &h_b);
    double v = 0.0;
    if (b < 1.0 - EPS_B_PURE_BURNT) { double ru = b_interp(urho, p, h_u, z); v += (1.0 - b) / ru; }
    if (b > EPS_B_PURE_UNBURNT)     { double rb = b_interp(brho, p, h_b + h_off, z); v += b / rb; }
    return v;
}
__device__ __forceinline__ double b_g(double p, double e, double inv, double z, double b,
                                      const Col3& urho, const Col3& brho,
                                      double hu_floor, double hu_ceil, double hb_lo, double hb_hi, double h_off) {
    return b_inv_rho(p, e + p * inv, z, b, urho, brho, hu_floor, hu_ceil, hb_lo, hb_hi, h_off) - inv;
}

// Illinois + bisection backstop on the blend g (same shape as the TableEos port).
__device__ double b_illinois(double a, double b_, double ga, double gb,
                             double e, double inv, double z, double bb, const Col3& urho, const Col3& brho,
                             double huf, double huc, double hbl, double hbh, double hoff) {
    for (int it = 0; it < N_P_ITER_MAX; it++) {
        if (fabs(b_ - a) <= EPS_P_PROJECTION * fmax(fabs(a), fabs(b_))) return 0.5 * (a + b_);
        double denom = gb - ga;
        double p = (denom != 0.0) ? b_ - gb * (b_ - a) / denom : 0.5 * (a + b_);
        double loe = fmin(a, b_), hie = fmax(a, b_);
        if (!(p > loe && p < hie)) p = 0.5 * (a + b_);
        double gp = b_g(p, e, inv, z, bb, urho, brho, huf, huc, hbl, hbh, hoff);
        if (gp == 0.0) return p;
        if (gp * gb < 0.0) { a = b_; ga = gb; } else { ga *= 0.5; }
        b_ = p; gb = gp;
    }
    if (fabs(b_ - a) <= EPS_P_PROJECTION * fmax(fabs(a), fabs(b_)) * 10.0) return 0.5 * (a + b_);
    double aa, bb2, gaa;
    if (a < b_) { aa = a; bb2 = b_; gaa = ga; } else { aa = b_; bb2 = a; gaa = gb; }
    for (int it = 0; it < N_P_BISECT; it++) {
        double m = 0.5 * (aa + bb2);
        if (fabs(bb2 - aa) <= EPS_P_PROJECTION * fmax(fabs(aa), fabs(bb2))) return m;
        double gm = b_g(m, e, inv, z, bb, urho, brho, huf, huc, hbl, hbh, hoff);
        if (gm == 0.0) return m;
        if (gaa * gm < 0.0) bb2 = m; else { aa = m; gaa = gm; }
    }
    return 0.5 * (aa + bb2);
}

// The MID-b blend projection: bracket = products h-window ∩ p-window; first
// crossing (continuity) → Illinois. Returns NaN on no-crossing (tangency
// corner, follow-on) or an invalid bracket.
__device__ double b_project_mid(double rho, double e, double z, double b, double inv,
                                const Col3& urho, const Col3& brho,
                                double pu_lo, double pu_hi, double huf, double huc,
                                double pb_lo, double pb_hi, double hbf, double hbc, double h_off) {
    // products window (blend_eos.rs mid-b): h in [hb_lo−off, hb_hi−off]; p in ∩.
    double h_lo = hbf - h_off, h_hi = hbc - h_off;
    double p_lo_env = fmax(pu_lo, pb_lo), p_hi_env = fmin(pu_hi, pb_hi);
    double p_from_h_lo = rho * (h_lo - e), p_from_h_hi = rho * (h_hi - e);
    double lo = fmax(fmax(p_lo_env, p_from_h_lo * (1.0 + H_BRACKET_MARGIN)), 0.0);
    double hi = fmin(p_hi_env, p_from_h_hi * (1.0 - H_BRACKET_MARGIN));
    if (!(isfinite(lo) && isfinite(hi)) || lo >= hi || hi <= 0.0) return nan("");
    // hb_lo/hi for partition = burnt h-envelope shifted by −off (same as window).
    double ga0 = b_g(lo, e, inv, z, b, urho, brho, huf, huc, hbf - h_off, hbc - h_off, h_off);
    double gb0 = b_g(hi, e, inv, z, b, urho, brho, huf, huc, hbf - h_off, hbc - h_off, h_off);
    if (ga0 == 0.0) return lo;
    if (gb0 == 0.0) return hi;
    // first-crossing scan (continuity — Ben ruling; no >1 refusal on device).
    double ratio = hi / lo, prev_p = lo, prev_g = ga0;
    for (int k = 1; k <= N_P_SCAN; k++) {
        double pk = lo * pow(ratio, (double)k / (double)N_P_SCAN);
        double gk = b_g(pk, e, inv, z, b, urho, brho, huf, huc, hbf - h_off, hbc - h_off, h_off);
        if (gk == 0.0) return pk;
        if (prev_g * gk < 0.0)
            return b_illinois(prev_p, pk, prev_g, gk, e, inv, z, b, urho, brho,
                              huf, huc, hbf - h_off, hbc - h_off, h_off);
        prev_p = pk; prev_g = gk;
    }
    return nan("");   // no crossing — tangency corner (follow-on)
}

// One mid-b projection per state.
__global__ void k_blend_project(const double* __restrict__ rho, const double* __restrict__ e,
                                const double* __restrict__ z, const double* __restrict__ b, int n,
                                Col3 urho, Col3 brho, double pu_lo, double pu_hi, double huf, double huc,
                                double pb_lo, double pb_hi, double hbf, double hbc, double h_off,
                                double* __restrict__ p_out) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double r = rho[i], inv = 1.0 / r;
    p_out[i] = b_project_mid(r, e[i], z[i], b[i], inv, urho, brho,
                             pu_lo, pu_hi, huf, huc, pb_lo, pb_hi, hbf, hbc, h_off);
}

// FFI: project N mid-b blend states. Two rho surfaces (unburnt, burnt) each with
// own axes; envelopes + h_offset scalars. Writes p_out per state.
extern "C" void gpu_blend_project(
    const double* rho, const double* e, const double* z, const double* b, int n,
    const double* up, int unp, int ulp, const double* uh, int unh, int ulh,
    const double* uz, int unz, int ulz, const int* ustr, const double* urho_data, int urho_vlog,
    const double* bp, int bnp, int blp, const double* bh, int bnh, int blh,
    const double* bz, int bnz, int blz, const int* bstr, const double* brho_data, int brho_vlog,
    double pu_lo, double pu_hi, double huf, double huc,
    double pb_lo, double pb_hi, double hbf, double hbc, double h_off,
    double* p_out) {
    long ncell = n; size_t nb = (size_t)n * sizeof(double);
    double *d_rho,*d_e,*d_z,*d_b,*d_p, *d_up,*d_uh,*d_uz,*d_ud, *d_bp,*d_bh,*d_bz,*d_bd;
    cudaMalloc(&d_rho,nb);cudaMalloc(&d_e,nb);cudaMalloc(&d_z,nb);cudaMalloc(&d_b,nb);cudaMalloc(&d_p,nb);
    cudaMalloc(&d_up,unp*sizeof(double));cudaMalloc(&d_uh,unh*sizeof(double));cudaMalloc(&d_uz,unz*sizeof(double));
    cudaMalloc(&d_bp,bnp*sizeof(double));cudaMalloc(&d_bh,bnh*sizeof(double));cudaMalloc(&d_bz,bnz*sizeof(double));
    long und = (long)unp*unh*unz, bnd = (long)bnp*bnh*bnz;
    cudaMalloc(&d_ud,und*sizeof(double));cudaMalloc(&d_bd,bnd*sizeof(double));
    cudaMemcpy(d_rho,rho,nb,cudaMemcpyHostToDevice);cudaMemcpy(d_e,e,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_z,z,nb,cudaMemcpyHostToDevice);cudaMemcpy(d_b,b,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_up,up,unp*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_uh,uh,unh*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_uz,uz,unz*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_ud,urho_data,und*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_bp,bp,bnp*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_bh,bh,bnh*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_bz,bz,bnz*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_bd,brho_data,bnd*sizeof(double),cudaMemcpyHostToDevice);
    Col3 urho{d_up,d_uh,d_uz,unp,unh,unz,ulp,ulh,ulz,ustr[0],ustr[1],ustr[2],d_ud,urho_vlog};
    Col3 brho{d_bp,d_bh,d_bz,bnp,bnh,bnz,blp,blh,blz,bstr[0],bstr[1],bstr[2],d_bd,brho_vlog};
    int tpb = 64; long blk = (n + tpb - 1) / tpb;
    k_blend_project<<<blk,tpb>>>(d_rho,d_e,d_z,d_b,n,urho,brho,
        pu_lo,pu_hi,huf,huc,pb_lo,pb_hi,hbf,hbc,h_off,d_p);
    cudaDeviceSynchronize();
    cudaMemcpy(p_out,d_p,nb,cudaMemcpyDeviceToHost);
    cudaFree(d_rho);cudaFree(d_e);cudaFree(d_z);cudaFree(d_b);cudaFree(d_p);
    cudaFree(d_up);cudaFree(d_uh);cudaFree(d_uz);cudaFree(d_ud);
    cudaFree(d_bp);cudaFree(d_bh);cudaFree(d_bz);cudaFree(d_bd);
}
