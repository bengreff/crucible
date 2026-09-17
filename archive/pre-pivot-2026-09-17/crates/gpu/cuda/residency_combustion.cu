// crucible-gpu — S13c residency: the COMBUSTION SOURCE (SOLV-4.4) on device —
// the bistable-Nagumo propagation reaction + the matched front-thickening
// diffusion that ride eval_rhs as the class-A burn-progress source. Bit-for-
// formula from crucible_solvers::euler::Combustion::accumulate_inner
// (combustion.rs). This is the flame-propagation physics: ρ_u·K·b(1−b)(b−a)
// (the pushed front) + ∇·(ρD_c∇b) (the matched thickening) with the coefficients
// D_c = w·S_T/(1−2a), K = 2·S_T/((1−2a)·w), w = Θ·Δ, S_T = S_L·wrinkling — so the
// front SPEED is closure-set (grid-independent), the S13c ignition headline.
//
// SCOPE (this milestone): the class-A propagation source + front diffusion,
// N_θ=1 box interior, direct interps (no root-find — the cell's p is already in
// the prim). The class-R implicit auto-ignition node solve (spontaneous light)
// is the remaining combustion leg. Surfaces: the UNBURNT (p,h,Z) temperature +
// density columns (T_u, ρ_u) and the IGNITION (p,T_u,Z) flame-speed column S_L.
//
// Determinism (META-1 §2.5): a per-cell gather (no reduction) — same-build
// reruns bit-identical; CPU↔GPU is interp FMA + libm(ln/exp) ECT.
#include <cuda_runtime.h>
#include <math.h>

// SOLV-4.4 constants (combustion.rs)
#define IGN_A 0.2
#define BURN_COMPLETE 1.0e-3
#define EPS_B_PURE_UNBURNT 1.0e-9

struct Col3 {
    const double *pp, *hp, *zp; int np, nh, nz; int lp, lh, lz;
    int s0, s1, s2; const double* data; int vlog;
};
__device__ __forceinline__ int c_cell(const double* pts, int n, double q) {
    int cnt = 0; for (int i = 0; i < n; i++) cnt += (pts[i] <= q);
    int hi = cnt < n-1 ? cnt : n-1; return (hi > 1 ? hi : 1) - 1;
}
__device__ __forceinline__ double c_frac(const double* pts, int i, int is_log, double q) {
    if (is_log) return (log(q) - log(pts[i])) / (log(pts[i+1]) - log(pts[i]));
    return (q - pts[i]) / (pts[i+1] - pts[i]);
}
__device__ double c_interp(const Col3& c, double qp, double qh, double qz) {
    int ip = c_cell(c.pp, c.np, qp), ih = c_cell(c.hp, c.nh, qh), iz = c_cell(c.zp, c.nz, qz);
    double t0 = c_frac(c.pp, ip, c.lp, qp), t1 = c_frac(c.hp, ih, c.lh, qh), t2 = c_frac(c.zp, iz, c.lz, qz);
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

// partition_h → h_u (blend_eos.rs partition_h, first element).
__device__ __forceinline__ double c_hu(double h, double b, double hu_floor, double hu_ceil) {
    if ((h >= hu_floor && h <= hu_ceil) || b <= EPS_B_PURE_UNBURNT) return h;
    return fmin(fmax(h, hu_floor), hu_ceil);
}
// front_coeffs (combustion.rs): (D_c, K); (0,0) if S_T<=0.
__device__ __forceinline__ void c_front(double s_t, double delta, double theta, double* d_c, double* k) {
    if (s_t <= 0.0) { *d_c = 0.0; *k = 0.0; return; }
    double w = theta * delta, f = 1.0 - 2.0 * IGN_A;
    *d_c = w * s_t / f; *k = 2.0 * s_t / (f * w);
}

// One cell's (react into ρc, ρ·D_c) — accumulate_inner's per-cell body +
// rho_dc_of, sharing the reaction gate. b is the RAW prim b; b_c its clamp.
__device__ void c_terms(double p, double h, double z, double rho, double b,
                        const Col3& tempc, const Col3& rhoc, const Col3& flamec,
                        double hu_floor, double hu_ceil, double p_floor, double tu_floor,
                        double wrinkling, double theta, double delta,
                        double* react, double* rho_dc) {
    *react = 0.0; *rho_dc = 0.0;
    double b_c = fmin(fmax(b, 0.0), 1.0);
    // burnt fixed point, or colder than any representable reactant → non-reactive
    if (b_c >= 1.0 - BURN_COMPLETE || h < hu_floor) return;
    double h_u = c_hu(h, b_c, hu_floor, hu_ceil);
    double t_u = c_interp(tempc, p, h_u, z);
    if (p < p_floor || t_u < tu_floor) return;          // cold-side declared no-burn
    double rho_u = c_interp(rhoc, p, h_u, z);
    double s_l = c_interp(flamec, p, t_u, z);
    double s_t = fmax(s_l * wrinkling, 0.0);
    double d_c, k; c_front(s_t, delta, theta, &d_c, &k);
    *react = rho_u * k * b_c * (1.0 - b_c) * (b_c - IGN_A);
    *rho_dc = rho * d_c;
}

// The combustion source per cell: rate[I_RB] = react + ∇·(ρD_c∇b)/kv. Inputs
// dense [ncell]; face neighbours in-domain gas only (a domain/wall edge is
// zero-flux for b). Writes ONLY the I_RB (=6) slot of the dense NC=7 rate.
#define NC 7
__global__ void k_comb_source(const double* __restrict__ rho, const double* __restrict__ p,
                              const double* __restrict__ z, const double* __restrict__ b,
                              const double* __restrict__ e, const double* __restrict__ gas,
                              int n_r, int n_z, double r_min, double dr, double dz,
                              Col3 tempc, Col3 rhoc, Col3 flamec,
                              double hu_floor, double hu_ceil, double p_floor, double tu_floor,
                              double wrinkling, double theta, double* __restrict__ rate) {
    long c = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long ncell = (long)n_r * n_z;
    if (c >= ncell) return;
    for (int k = 0; k < NC; k++) rate[c*NC+k] = 0.0;
    if (gas[c] == 0.0) return;
    int i_r = (int)(c / n_z), i_z = (int)(c % n_z);
    double delta = sqrt(dr * dz);
    double r_i = r_min + i_r * dr, r_o = r_min + (i_r + 1) * dr;
    double TAU = 6.283185307179586;
    double vol = 0.5 * (r_o*r_o - r_i*r_i) * TAU * dz;
    double kv = vol;                                     // kappa = 1
    double h_c = e[c] + p[c] / rho[c];
    double react, rho_dc_here;
    c_terms(p[c], h_c, z[c], rho[c], b[c], tempc, rhoc, flamec,
            hu_floor, hu_ceil, p_floor, tu_floor, wrinkling, theta, delta, &react, &rho_dc_here);
    double diff = 0.0;
    // 4 faces: r-,r+,z-,z+
    for (int fc = 0; fc < 4; fc++) {
        int nr = i_r, nz = i_z; bool radial = fc < 2;
        if (fc == 0) nr = i_r - 1; else if (fc == 1) nr = i_r + 1;
        else if (fc == 2) nz = i_z - 1; else nz = i_z + 1;
        if (nr < 0 || nr >= n_r || nz < 0 || nz >= n_z) continue;
        long nc = (long)nr * n_z + nz;
        if (gas[nc] == 0.0) continue;
        double area = radial
            ? (r_min + (fc == 1 ? i_r + 1 : i_r) * dr) * TAU * dz     // face_area_r
            : 0.5 * (r_o*r_o - r_i*r_i) * TAU;                        // face_area_z
        double d = radial ? dr : dz;
        double h_n = e[nc] + p[nc] / rho[nc];
        double react_n, rho_dc_n;
        c_terms(p[nc], h_n, z[nc], rho[nc], b[nc], tempc, rhoc, flamec,
                hu_floor, hu_ceil, p_floor, tu_floor, wrinkling, theta, delta, &react_n, &rho_dc_n);
        double rho_dc_face = 0.5 * (rho_dc_here + rho_dc_n);
        diff += area * rho_dc_face * (b[nc] - b[c]) / d;
    }
    double diff_div = kv > 0.0 ? diff / kv : 0.0;
    rate[c*NC + 6] = react + diff_div;                   // I_RB
}

extern "C" void gpu_combustion_source(
    const double* rho, const double* p, const double* z, const double* b, const double* e,
    const double* gas, int n_r, int n_z, double r_min, double dr, double dz,
    // unburnt axes (shared by temp+rho_u) + two data blocks
    const double* up, int unp, int ulp, const double* uh, int unh, int ulh,
    const double* uz, int unz, int ulz, const int* ustr,
    const double* temp_data, int temp_vlog, const double* rhou_data, int rhou_vlog,
    // ignition axes + flame data
    const double* ip, int inp, int ilp, const double* it, int inh, int ilh,
    const double* iz, int inz, int ilz, const int* istr, const double* flame_data, int flame_vlog,
    double hu_floor, double hu_ceil, double p_floor, double tu_floor, double wrinkling, double theta,
    double* rate) {
    long ncell = (long)n_r * n_z; size_t nb = (size_t)ncell * sizeof(double);
    double *d_rho,*d_p,*d_z,*d_b,*d_e,*d_gas,*d_rate;
    double *d_up,*d_uh,*d_uz,*d_temp,*d_rhou,*d_ip,*d_it,*d_iz,*d_flame;
    cudaMalloc(&d_rho,nb);cudaMalloc(&d_p,nb);cudaMalloc(&d_z,nb);cudaMalloc(&d_b,nb);
    cudaMalloc(&d_e,nb);cudaMalloc(&d_gas,nb);cudaMalloc(&d_rate,ncell*NC*sizeof(double));
    long undata = (long)unp*unh*unz, indata = (long)inp*inh*inz;
    cudaMalloc(&d_up,unp*sizeof(double));cudaMalloc(&d_uh,unh*sizeof(double));cudaMalloc(&d_uz,unz*sizeof(double));
    cudaMalloc(&d_temp,undata*sizeof(double));cudaMalloc(&d_rhou,undata*sizeof(double));
    cudaMalloc(&d_ip,inp*sizeof(double));cudaMalloc(&d_it,inh*sizeof(double));cudaMalloc(&d_iz,inz*sizeof(double));
    cudaMalloc(&d_flame,indata*sizeof(double));
    cudaMemcpy(d_rho,rho,nb,cudaMemcpyHostToDevice);cudaMemcpy(d_p,p,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_z,z,nb,cudaMemcpyHostToDevice);cudaMemcpy(d_b,b,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_e,e,nb,cudaMemcpyHostToDevice);cudaMemcpy(d_gas,gas,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_up,up,unp*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_uh,uh,unh*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_uz,uz,unz*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_temp,temp_data,undata*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_rhou,rhou_data,undata*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_ip,ip,inp*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_it,it,inh*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_iz,iz,inz*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(d_flame,flame_data,indata*sizeof(double),cudaMemcpyHostToDevice);
    Col3 tempc{d_up,d_uh,d_uz,unp,unh,unz,ulp,ulh,ulz,ustr[0],ustr[1],ustr[2],d_temp,temp_vlog};
    Col3 rhoc{d_up,d_uh,d_uz,unp,unh,unz,ulp,ulh,ulz,ustr[0],ustr[1],ustr[2],d_rhou,rhou_vlog};
    Col3 flamec{d_ip,d_it,d_iz,inp,inh,inz,ilp,ilh,ilz,istr[0],istr[1],istr[2],d_flame,flame_vlog};
    int tpb = 128; long blk = (ncell + tpb - 1) / tpb;
    k_comb_source<<<blk,tpb>>>(d_rho,d_p,d_z,d_b,d_e,d_gas,n_r,n_z,r_min,dr,dz,
        tempc,rhoc,flamec,hu_floor,hu_ceil,p_floor,tu_floor,wrinkling,theta,d_rate);
    cudaDeviceSynchronize();
    cudaMemcpy(rate,d_rate,ncell*NC*sizeof(double),cudaMemcpyDeviceToHost);
    cudaFree(d_rho);cudaFree(d_p);cudaFree(d_z);cudaFree(d_b);cudaFree(d_e);cudaFree(d_gas);
    cudaFree(d_rate);cudaFree(d_up);cudaFree(d_uh);cudaFree(d_uz);cudaFree(d_temp);cudaFree(d_rhou);
    cudaFree(d_ip);cudaFree(d_it);cudaFree(d_iz);cudaFree(d_flame);
}
