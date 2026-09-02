// crucible-gpu — S13c residency: the REAL HDF5 TableEos (p,h,Z) projection on
// device, replacing the class-A/D GammaLaw stand-in with the actual
// equilibrium surface. Bit-for-formula from crucible_solvers::euler::TableEos
// (table_eos.rs): the fixed 8-corner multilinear interp in interp_rule space +
// the deterministic Illinois regula-falsi pressure projection (warm-started
// hint + cold full-bracket), with the two-root behaviour RESOLVED BY
// CONTINUITY (Ben ruling 2026-08-31): the warm path keeps the root near the
// hint, the cold path takes the first scan crossing — neither halts.
//
// SCOPE (this milestone): the DOMINANT paths — warm-hinted fast + cold
// straddling-bracket Illinois + the near-vacuum first-crossing scan — which
// cover equilibrium chamber states (the ◆C3 regime). The rare near-vacuum
// GOLDEN-SECTION TANGENCY corner (no sign change over the whole scan) is NOT
// ported here — it returns a NaN sentinel so the harness fixture (interior
// states) never silently takes a wrong branch; that corner rides a follow-on.
//
// Determinism (META-1 §2.5): a per-cell gather (no reduction) — same-build
// reruns bit-identical; CPU↔GPU is FMA-order + libm(ln/exp) ECT (the interp's
// corner-reduction order is identical to the CPU; only ln/exp differ per host).
#include <cuda_runtime.h>
#include <math.h>

// projection constants (table_eos.rs)
#define EPS_P_PROJECTION 1e-11
#define H_BRACKET_MARGIN 1e-12
#define N_P_ITER_MAX 48
#define N_P_BISECT 64
#define N_P_SCAN 64

// A 3-D bound column marshaled flat (shared axes p,h,Z across the 3 columns;
// each column its own data block + value log-flag).
struct Col3 {
    const double *pp, *hp, *zp;   // axis grid points
    int np, nh, nz;               // axis lengths
    int lp, lh, lz;               // per-axis log flag
    int s0, s1, s2;               // row-major strides
    const double* data;           // value block
    int vlog;                     // value log flag
};

// partition_point(pts <= q) then hi=min(cnt,n-1), i=max(hi,1)-1 (table_eos.rs
// interpolate). Integer result is exact ⇒ identical cell index CPU↔GPU.
__device__ __forceinline__ int e_cell(const double* pts, int n, double q) {
    int cnt = 0;
    for (int i = 0; i < n; i++) cnt += (pts[i] <= q);
    int hi = cnt < n-1 ? cnt : n-1;
    return (hi > 1 ? hi : 1) - 1;
}
__device__ __forceinline__ double e_frac(const double* pts, int i, int is_log, double q) {
    if (is_log) return (log(q) - log(pts[i])) / (log(pts[i+1]) - log(pts[i]));
    return (q - pts[i]) / (pts[i+1] - pts[i]);
}
// Fixed 8-corner multilinear in interp_rule space (identical reduction order to
// BoundColumn::interpolate). Clamps rather than refuses (the projection bracket
// keeps the query in-domain; interior fixture stays inside).
__device__ double e_interp(const Col3& c, double qp, double qh, double qz) {
    int ip = e_cell(c.pp, c.np, qp), ih = e_cell(c.hp, c.nh, qh), iz = e_cell(c.zp, c.nz, qz);
    double t0 = e_frac(c.pp, ip, c.lp, qp);
    double t1 = e_frac(c.hp, ih, c.lh, qh);
    double t2 = e_frac(c.zp, iz, c.lz, qz);
    int cell[3] = {ip, ih, iz};
    double t[3] = {t0, t1, t2};
    int str[3] = {c.s0, c.s1, c.s2};
    double acc = 0.0;
    for (int corner = 0; corner < 8; corner++) {
        double w = 1.0; int idx = 0;
        for (int d = 0; d < 3; d++) {
            int up = (corner >> d) & 1;
            w *= up ? t[d] : (1.0 - t[d]);
            idx += (cell[d] + up) * str[d];
        }
        double v = c.vlog ? log(c.data[idx]) : c.data[idx];
        acc += w * v;
    }
    return c.vlog ? exp(acc) : acc;
}

// g(p) = ρ_tab(p, e_q + p/ρ, z) − ρ.
__device__ __forceinline__ double e_g(const Col3& rho_col, double p, double e_q, double inv,
                                      double z, double rho) {
    return e_interp(rho_col, p, e_q + p * inv, z) - rho;
}

// Deterministic Illinois regula-falsi + bisection backstop (illinois_root).
__device__ double e_illinois(double a, double b, double ga, double gb,
                             const Col3& rc, double e_q, double inv, double z, double rho) {
    for (int it = 0; it < N_P_ITER_MAX; it++) {
        if (fabs(b - a) <= EPS_P_PROJECTION * fmax(fabs(a), fabs(b))) return 0.5 * (a + b);
        double denom = gb - ga;
        double p = (denom != 0.0) ? b - gb * (b - a) / denom : 0.5 * (a + b);
        double loe = fmin(a, b), hie = fmax(a, b);
        if (!(p > loe && p < hie)) p = 0.5 * (a + b);
        double gp = e_g(rc, p, e_q, inv, z, rho);
        if (gp == 0.0) return p;
        if (gp * gb < 0.0) { a = b; ga = gb; }
        else { ga *= 0.5; }
        b = p; gb = gp;
    }
    if (fabs(b - a) <= EPS_P_PROJECTION * fmax(fabs(a), fabs(b)) * 10.0) return 0.5 * (a + b);
    // bisection backstop (bracket sign-changed by construction)
    double aa, bb, gaa;
    if (a < b) { aa = a; bb = b; gaa = ga; } else { aa = b; bb = a; gaa = gb; }
    for (int it = 0; it < N_P_BISECT; it++) {
        double m = 0.5 * (aa + bb);
        if (fabs(bb - aa) <= EPS_P_PROJECTION * fmax(fabs(aa), fabs(bb))) return m;
        double gm = e_g(rc, m, e_q, inv, z, rho);
        if (gm == 0.0) return m;
        if (gaa * gm < 0.0) bb = m; else { aa = m; gaa = gm; }
    }
    return 0.5 * (aa + bb);
}

// cold project_pressure: admissible bracket → straddle→Illinois; else scan for
// the FIRST sign change → Illinois; the golden-section tangency corner is not
// ported (returns NaN — interior fixture never reaches it).
__device__ double e_project_cold(const Col3& rc, double rho, double e_q, double z,
                                 double inv, double p_lo_env, double p_hi_env,
                                 double h_lo_env, double h_hi_env) {
    double p_from_h_lo = rho * (h_lo_env - e_q);
    double p_from_h_hi = rho * (h_hi_env - e_q);
    double lo = fmax(p_lo_env, p_from_h_lo * (1.0 + H_BRACKET_MARGIN));
    double hi = fmin(p_hi_env, p_from_h_hi * (1.0 - H_BRACKET_MARGIN));
    if (!(isfinite(lo) && isfinite(hi)) || lo >= hi || hi <= 0.0) return nan("");
    lo = fmax(lo, 0.0);
    double ga0 = e_g(rc, lo, e_q, inv, z, rho);
    double gb0 = e_g(rc, hi, e_q, inv, z, rho);
    if (ga0 == 0.0) return lo;
    if (gb0 == 0.0) return hi;
    if (ga0 * gb0 < 0.0) return e_illinois(lo, hi, ga0, gb0, rc, e_q, inv, z, rho);
    // near-vacuum: fixed log-scan for the FIRST sign change.
    double prev_p = lo, prev_g = ga0;
    double ratio = hi / lo;
    for (int k = 1; k <= N_P_SCAN; k++) {
        double pk = lo * pow(ratio, (double)k / (double)N_P_SCAN);
        double gk = e_g(rc, pk, e_q, inv, z, rho);
        if (gk == 0.0) return pk;
        if (prev_g * gk < 0.0)
            return e_illinois(prev_p, pk, prev_g, gk, rc, e_q, inv, z, rho);
        prev_p = pk; prev_g = gk;
    }
    return nan("");   // golden-section tangency corner — not ported (follow-on)
}

// warm project_pressure_hinted: tight bracket around the hint, gated on the
// full bracket straddling (continuity — same root as cold); else cold.
__device__ double e_project_hinted(const Col3& rc, double rho, double e_q, double z,
                                   double inv, double p_hint, double p_lo_env, double p_hi_env,
                                   double h_lo_env, double h_hi_env) {
    const double HINT_SPREAD = 1.05;
    double p_from_h_lo = rho * (h_lo_env - e_q);
    double p_from_h_hi = rho * (h_hi_env - e_q);
    double lo_adm = fmax(fmax(p_lo_env, p_from_h_lo * (1.0 + H_BRACKET_MARGIN)), 0.0);
    double hi_adm = fmin(p_hi_env, p_from_h_hi * (1.0 - H_BRACKET_MARGIN));
    double a = fmax(p_hint / HINT_SPREAD, lo_adm);
    double b = fmin(p_hint * HINT_SPREAD, hi_adm);
    if (!(isfinite(a) && isfinite(b)) || a >= b)
        return e_project_cold(rc, rho, e_q, z, inv, p_lo_env, p_hi_env, h_lo_env, h_hi_env);
    double ga_full = e_g(rc, lo_adm, e_q, inv, z, rho);
    double gb_full = e_g(rc, hi_adm, e_q, inv, z, rho);
    if (!(ga_full * gb_full < 0.0))
        return e_project_cold(rc, rho, e_q, z, inv, p_lo_env, p_hi_env, h_lo_env, h_hi_env);
    double ga = e_g(rc, a, e_q, inv, z, rho);
    if (ga == 0.0) return a;
    double gb = e_g(rc, b, e_q, inv, z, rho);
    if (gb == 0.0) return b;
    if (ga * gb < 0.0) return e_illinois(a, b, ga, gb, rc, e_q, inv, z, rho);
    return e_illinois(lo_adm, hi_adm, ga_full, gb_full, rc, e_q, inv, z, rho);
}

// One projection per state: (rho, e_q, z, hint) → p, then sound & temperature
// at (p, h=e_q+p/rho, z). g1 = rho*a*a/p (the GammaLaw-equivalent stiffness the
// prim carries). Writes p_out, a_out, t_out, g1_out per state.
__global__ void k_project(const double* __restrict__ rho, const double* __restrict__ e_q,
                          const double* __restrict__ z, const double* __restrict__ hint,
                          int n, Col3 rc, Col3 sc, Col3 tc,
                          double p_lo, double p_hi, double h_lo, double h_hi,
                          double* __restrict__ p_out, double* __restrict__ a_out,
                          double* __restrict__ t_out, double* __restrict__ g1_out) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double r = rho[i], eq = e_q[i], zz = z[i], inv = 1.0 / r;
    double ph = hint[i];
    double p = (isfinite(ph) && ph > 0.0)
        ? e_project_hinted(rc, r, eq, zz, inv, ph, p_lo, p_hi, h_lo, h_hi)
        : e_project_cold(rc, r, eq, zz, inv, p_lo, p_hi, h_lo, h_hi);
    double h = eq + p * inv;
    double a = e_interp(sc, p, h, zz);
    double tt = e_interp(tc, p, h, zz);
    p_out[i] = p; a_out[i] = a; t_out[i] = tt; g1_out[i] = r * a * a / p;
}

// FFI: project N states on the real equilibrium surface. Axes (pp,hp,zp) shared
// across the 3 columns; each column its own data + value-log flag + strides.
extern "C" void gpu_table_project(
    const double* rho, const double* e_q, const double* z, const double* hint, int n,
    const double* pp, int np, int lp, const double* hp, int nh, int lh,
    const double* zp, int nz, int lz, const int* strides,
    const double* rho_data, int rho_vlog, const double* snd_data, int snd_vlog,
    const double* tmp_data, int tmp_vlog,
    double p_lo, double p_hi, double h_lo, double h_hi,
    double* p_out, double* a_out, double* t_out, double* g1_out) {
    size_t nb = (size_t)n * sizeof(double);
    double *d_rho,*d_eq,*d_z,*d_hint,*d_p,*d_a,*d_t,*d_g1;
    double *d_pp,*d_hp,*d_zp,*d_rd,*d_sd,*d_td; int *d_str;
    cudaMalloc(&d_rho,nb);cudaMalloc(&d_eq,nb);cudaMalloc(&d_z,nb);cudaMalloc(&d_hint,nb);
    cudaMalloc(&d_p,nb);cudaMalloc(&d_a,nb);cudaMalloc(&d_t,nb);cudaMalloc(&d_g1,nb);
    cudaMalloc(&d_pp,np*sizeof(double));cudaMalloc(&d_hp,nh*sizeof(double));cudaMalloc(&d_zp,nz*sizeof(double));
    long ndata = (long)np*nh*nz;
    cudaMalloc(&d_rd,ndata*sizeof(double));cudaMalloc(&d_sd,ndata*sizeof(double));cudaMalloc(&d_td,ndata*sizeof(double));
    cudaMalloc(&d_str,3*sizeof(int));
    cudaMemcpy(d_rho,rho,nb,cudaMemcpyHostToDevice);cudaMemcpy(d_eq,e_q,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_z,z,nb,cudaMemcpyHostToDevice);cudaMemcpy(d_hint,hint,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_pp,pp,np*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_hp,hp,nh*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_zp,zp,nz*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_rd,rho_data,ndata*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_sd,snd_data,ndata*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_td,tmp_data,ndata*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_str,strides,3*sizeof(int),cudaMemcpyHostToDevice);
    Col3 rc{d_pp,d_hp,d_zp,np,nh,nz,lp,lh,lz,strides[0],strides[1],strides[2],d_rd,rho_vlog};
    Col3 sc{d_pp,d_hp,d_zp,np,nh,nz,lp,lh,lz,strides[0],strides[1],strides[2],d_sd,snd_vlog};
    Col3 tc{d_pp,d_hp,d_zp,np,nh,nz,lp,lh,lz,strides[0],strides[1],strides[2],d_td,tmp_vlog};
    int tpb = 128; long blk = (n + tpb - 1) / tpb;
    k_project<<<blk,tpb>>>(d_rho,d_eq,d_z,d_hint,n,rc,sc,tc,p_lo,p_hi,h_lo,h_hi,d_p,d_a,d_t,d_g1);
    cudaDeviceSynchronize();
    cudaMemcpy(p_out,d_p,nb,cudaMemcpyDeviceToHost);cudaMemcpy(a_out,d_a,nb,cudaMemcpyDeviceToHost);
    cudaMemcpy(t_out,d_t,nb,cudaMemcpyDeviceToHost);cudaMemcpy(g1_out,d_g1,nb,cudaMemcpyDeviceToHost);
    cudaFree(d_rho);cudaFree(d_eq);cudaFree(d_z);cudaFree(d_hint);cudaFree(d_p);cudaFree(d_a);
    cudaFree(d_t);cudaFree(d_g1);cudaFree(d_pp);cudaFree(d_hp);cudaFree(d_zp);cudaFree(d_rd);
    cudaFree(d_sd);cudaFree(d_td);cudaFree(d_str);
}
