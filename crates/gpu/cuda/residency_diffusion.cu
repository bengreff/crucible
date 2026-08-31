// crucible-gpu — S13b residency: the class-D IMPLICIT DIFFUSION per-component
// symmetric CG on device, bit-for-formula from crucible_solvers::gas_diffusion
// (apply_linear / fill_mass / cg_solve). This is the "single biggest piece" of
// the S13b split (the brief): the fixed-structure Jacobi-preconditioned CG that
// solves (mass - wq*L)·δ = b for one gas component, with the deterministic
// fixed-topology tree reduction that the CG dot products require.
//
// WHY THIS IS THE CRUX. The S13 class-A residency had ZERO reductions (pure
// gather stencils). The class-D CG introduces the first on-device REDUCTIONS
// (the dot products), which is the genuinely new determinism surface: a
// fixed-topology tree reduction (META-1 §2.5 — "fixed-topology tree reductions
// for all grid/ensemble statistics"). Nail it here with a bit-exact CPU
// reference and every downstream residency piece (Robin coupling, stable_dt,
// combustion node solves) reuses the primitive.
//
// SCOPE (S13b, the class-D CG core):
//   N_theta = 1, box world (kappa = aperture = 1), free (zero-flux Neumann)
//   BCs on the domain edge, constant-or-per-cell transport read as the two-cell
//   face average (FaceTr::between). The five components (Ur,Uz,Om,T,C) all use
//   the ONE component-generic cg_solve; only face_coef + mass differ. The
//   Picard-lagged cross terms live in the RHS `b` (assemble_rates), NOT in the
//   CG matrix (module doc), so the b-assembly + Robin-Robin wall coupling +
//   solid-conduction CG + combustion + real HDF5 TableEos + cut/mixed-N_theta
//   geometry + stable_dt are the S13c split (recorded in the SESSION_LOG).
//
// DETERMINISM (META-1 §2.5): gather-only (one writer per cell), no atomics in
// the physics/reduction paths, fixed control flow, fixed-topology reductions.
// Same-build reruns are bit-identical; CPU<->GPU agree to the declared ECT
// (FMA-order + reduction-shape difference — a converged CG solve, not a single
// op). Device symbols prefixed `d_` so they never collide with the other .cu.
#include <cuda_runtime.h>
#include <math.h>   // host-side sqrt in gpu_class_d_cg (NaN via self-inequality)

// Component ids (match GasComp order used by the harness).
#define COMP_UR 0
#define COMP_UZ 1
#define COMP_OM 2
#define COMP_T  3
#define COMP_C  4

#define TAU 6.283185307179586
#define DTPB 256   // reduction/elementwise threads-per-block (power of two)

// ---- geometry (grid/lib.rs, N_theta = 1 so TAU/nt = TAU) ----------------
__device__ __forceinline__ double d_face_radius(double r_min, double dr, int f) {
    return r_min + (double)f * dr;
}
__device__ __forceinline__ double d_r_center(double r_min, double dr, int i_r) {
    return r_min + ((double)i_r + 0.5) * dr;
}
__device__ __forceinline__ double d_cell_volume(double r_min, double dr, double dz, int i_r) {
    double r_i = d_face_radius(r_min, dr, i_r), r_o = d_face_radius(r_min, dr, i_r + 1);
    return 0.5 * (r_o * r_o - r_i * r_i) * TAU * dz;
}
__device__ __forceinline__ double d_area_r(double r_min, double dr, double dz, int i_r, bool outer) {
    double r = d_face_radius(r_min, dr, outer ? i_r + 1 : i_r);
    return r * TAU * dz;
}
__device__ __forceinline__ double d_area_z(double r_min, double dr, int i_r) {
    double r_i = d_face_radius(r_min, dr, i_r), r_o = d_face_radius(r_min, dr, i_r + 1);
    return 0.5 * (r_o * r_o - r_i * r_i) * TAU;
}

// ---- face_coef (gas_diffusion.rs face_coef) -----------------------------
// `radial` = true for an r-face, false for a z-face. `mu_f`,`k_f`,`rhod_f` are
// the two-cell face-averaged transport; `r_face` is the face radius (Om form).
__device__ __forceinline__ double d_face_coef(int comp, bool radial, double r_face,
                                              double mu_f, double k_f, double rhod_f) {
    switch (comp) {
        case COMP_UR: return radial ? (4.0 / 3.0) * mu_f : mu_f;
        case COMP_UZ: return radial ? mu_f : (4.0 / 3.0) * mu_f;
        case COMP_OM: return mu_f * r_face * r_face;
        case COMP_T:  return k_f;
        default:      return rhod_f;   // COMP_C
    }
}

// ---- apply_linear: out = L·x (+ optional diag = diag(L)) ----------------
// Cell = i_r*n_z + i_z. Gas mask `gas` (1.0/0.0). Interior (r,z) faces carry
// the two-point term coef*area*(x_nbr - x_c)/dist; a face to an out-of-domain
// neighbour is a free (zero-flux) BC and contributes nothing (this fixture's
// BCs; cut/solid/no-slip walls are S13c). The Ur negative-definite geometric
// diagonal −(4/3)μ·u_r/r̄·geo·κV is added last (SOLV-1 §3.3 pattern). One writer
// per cell; `diag` is written iff `want_diag`.
__global__ void kd_apply_linear(const double* __restrict__ x, int comp,
                                int n_r, int n_z, double r_min, double dr, double dz,
                                const double* __restrict__ mu, const double* __restrict__ kk,
                                const double* __restrict__ rhod, const double* __restrict__ gas,
                                int want_diag, double* __restrict__ out,
                                double* __restrict__ diag) {
    long c = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long ncell = (long)n_r * n_z;
    if (c >= ncell) return;
    if (gas[c] == 0.0) { out[c] = 0.0; if (want_diag) diag[c] = 0.0; return; }
    int i_r = (int)(c / n_z), i_z = (int)(c % n_z);
    double x_c = x[c];
    double mu_c = mu[c], k_c = kk[c], d_c = rhod[c];
    double acc = 0.0, dg = 0.0;
    // Four (r,z) faces, fixed order r−, r+, z−, z+.
    // r− and r+
    for (int s = 0; s < 2; s++) {
        bool outer = (s == 1);
        int nr = outer ? i_r + 1 : i_r - 1;
        if (nr < 0 || nr >= n_r) continue;             // free BC ⇒ zero-flux
        long nc = (long)nr * n_z + i_z;
        if (gas[nc] == 0.0) continue;                   // (all-gas box; solids = S13c)
        double mu_f = 0.5 * (mu_c + mu[nc]);
        double k_f  = 0.5 * (k_c + kk[nc]);
        double d_f  = 0.5 * (d_c + rhod[nc]);
        double r_face = d_face_radius(r_min, dr, outer ? i_r + 1 : i_r);
        double area = d_area_r(r_min, dr, dz, i_r, outer);
        double coef = d_face_coef(comp, true, r_face, mu_f, k_f, d_f);
        acc += coef * area * (x[nc] - x_c) / dr;
        dg  -= coef * area / dr;
    }
    // z− and z+ (r_face = r_center for the Om form)
    double rbar = d_r_center(r_min, dr, i_r);
    for (int s = 0; s < 2; s++) {
        bool outer = (s == 1);
        int nz = outer ? i_z + 1 : i_z - 1;
        if (nz < 0 || nz >= n_z) continue;
        long nc = (long)i_r * n_z + nz;
        if (gas[nc] == 0.0) continue;
        double mu_f = 0.5 * (mu_c + mu[nc]);
        double k_f  = 0.5 * (k_c + kk[nc]);
        double d_f  = 0.5 * (d_c + rhod[nc]);
        double area = d_area_z(r_min, dr, i_r);
        double coef = d_face_coef(comp, false, rbar, mu_f, k_f, d_f);
        acc += coef * area * (x[nc] - x_c) / dz;
        dg  -= coef * area / dz;
    }
    if (comp == COMP_UR) {
        double vol = d_cell_volume(r_min, dr, dz, i_r);
        double geo = (d_area_r(r_min, dr, dz, i_r, true) - d_area_r(r_min, dr, dz, i_r, false)) / vol;
        double kv = vol;                                // kappa = 1 (box)
        double cc = (4.0 / 3.0) * mu_c * geo / rbar * kv;
        acc -= cc * x_c;
        dg  -= cc;
    }
    out[c] = acc;
    if (want_diag) diag[c] = dg;
}

// ---- fill_mass (gas_diffusion.rs fill_mass), kv = kappa*vol (kappa = 1) --
__global__ void kd_fill_mass(int comp, int n_r, int n_z, double r_min, double dr, double dz,
                             const double* __restrict__ rho, const double* __restrict__ cv,
                             const double* __restrict__ gas, double* __restrict__ mass) {
    long c = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long ncell = (long)n_r * n_z;
    if (c >= ncell) return;
    if (gas[c] == 0.0) { mass[c] = 1.0; return; }
    int i_r = (int)(c / n_z);
    double kv = d_cell_volume(r_min, dr, dz, i_r);      // kappa = 1
    double rho_c = rho[c];
    switch (comp) {
        case COMP_OM: { double r = d_r_center(r_min, dr, i_r); mass[c] = rho_c * r * r * kv; break; }
        case COMP_T:  mass[c] = rho_c * cv[c] * kv; break;
        default:      mass[c] = rho_c * kv;               // Ur, Uz, C
    }
}

// ---- elementwise CG kernels --------------------------------------------
__global__ void kd_copy(double* __restrict__ dst, const double* __restrict__ src, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return; dst[i] = src[i];
}
// diag <- gas ? mass - wq*diagL : 1  (in place: diag currently holds diagL)
__global__ void kd_diag_finish(double* __restrict__ diag, const double* __restrict__ mass,
                               const double* __restrict__ gas, double wq, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return;
    diag[i] = (gas[i] != 0.0) ? (mass[i] - wq * diag[i]) : 1.0;
}
// q <- gas ? mass*p - wq*Lp : 0   (in place: q currently holds Lp)
__global__ void kd_q_finish(double* __restrict__ q, const double* __restrict__ mass,
                            const double* __restrict__ p, const double* __restrict__ gas,
                            double wq, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return;
    double lp = q[i];
    q[i] = (gas[i] != 0.0) ? (mass[i] * p[i] - wq * lp) : 0.0;
}
// z <- r / diag
__global__ void kd_precond(double* __restrict__ z, const double* __restrict__ r,
                           const double* __restrict__ diag, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return; z[i] = r[i] / diag[i];
}
// y <- y + a*x
__global__ void kd_axpy(double* __restrict__ y, double a, const double* __restrict__ x, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return; y[i] += a * x[i];
}
// p <- z + beta*p
__global__ void kd_p_update(double* __restrict__ p, const double* __restrict__ z, double beta, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return; p[i] = z[i] + beta * p[i];
}
// x <- x + delta on GAS cells only (cg_solve applies δ over gas cells; masked
// δ is computed in the iteration but discarded).
__global__ void kd_x_add_delta_gas(double* __restrict__ x, const double* __restrict__ delta,
                                   const double* __restrict__ gas, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x; if (i >= n) return;
    if (gas[i] != 0.0) x[i] += delta[i];
}

// ---- fixed-topology tree reduction (the determinism primitive) ----------
// Stage 1: each block reduces its DTPB masked products a[i]*b[i] via the
// canonical pairwise shared-memory tree (fixed shape for DTPB a power of two,
// every lane padded with 0) into one partial. Stage 2: a single block reduces
// the block partials the same way (each thread first does a fixed-stride serial
// pre-sum, then the pairwise tree). Shape depends only on (ncell, DTPB) — never
// on scheduling ⇒ same-build reruns bit-identical (META-1 §2.5). CPU<->GPU
// differ only at reduction-shape/FMA order (declared ECT).
__global__ void kd_dot_partial(const double* __restrict__ a, const double* __restrict__ b,
                               const double* __restrict__ gas, long ncell,
                               double* __restrict__ partial) {
    __shared__ double sh[DTPB];
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    double v = 0.0;
    if (i < ncell && gas[i] != 0.0) v = a[i] * b[i];
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int s = DTPB / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sh[threadIdx.x] += sh[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) partial[blockIdx.x] = sh[0];
}
__global__ void kd_dot_final(const double* __restrict__ partial, int nblk, double* __restrict__ out) {
    __shared__ double sh[DTPB];
    double v = 0.0;
    for (int j = threadIdx.x; j < nblk; j += DTPB) v += partial[j];   // fixed-stride serial pre-sum
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int s = DTPB / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sh[threadIdx.x] += sh[threadIdx.x + s];
        __syncthreads();
    }
    if (threadIdx.x == 0) *out = sh[0];
}

// Host-side dot: fields resident on device, only the scalar result copies back
// (the CG scalar recurrence — alpha/beta/termination — runs host-side and
// deterministically; the FIELD vectors never leave the device across the CG
// loop, which is the residency contract).
static double dev_dot(const double* a, const double* b, const double* gas, long ncell,
                      double* d_partial, double* d_scalar, int nblk) {
    kd_dot_partial<<<nblk, DTPB>>>(a, b, gas, ncell, d_partial);
    kd_dot_final<<<1, DTPB>>>(d_partial, nblk, d_scalar);
    double h; cudaMemcpy(&h, d_scalar, sizeof(double), cudaMemcpyDeviceToHost);
    return h;
}

// ---- FFI: the resident per-component class-D CG -------------------------
// Solves (mass - wq*L)·δ = b for component `comp`, then x += δ over gas cells,
// bit-for-formula from GasDiffusion::cg_solve. Inputs x0/b/rho/cv/mu/k/rhod/gas
// are dense [ncell]; `x` is inout (x0 in, x0+δ out). If `fixed_iters` > 0 the
// loop runs EXACTLY that many CG iterations with no early break (the clean ECT
// comparison — same work as the CPU's measured count); if <= 0 it runs the
// CPU's data-dependent termination (own reductions) capped at N_CG_ITERS_MAX.
// Writes the iteration count to *iters_out and the terminal relative residual
// to *resid_out. Returns 1 if the CG missed EPS acceptance (never clamps).
#define N_CG_ITERS_MAX 512
#define EPS_CG_RESID 1e-12
extern "C" int gpu_class_d_cg(double* x, const double* b, const double* rho, const double* cv,
                              const double* mu, const double* kk, const double* rhod, const double* gas,
                              int comp, int n_r, int n_z, double r_min, double dr, double dz,
                              double wq, int fixed_iters, int* iters_out, double* resid_out) {
    long ncell = (long)n_r * n_z;
    size_t bytes = ncell * sizeof(double);
    int nblk = (int)((ncell + DTPB - 1) / DTPB);
    double *d_x,*d_b,*d_rho,*d_cv,*d_mu,*d_k,*d_rd,*d_gas;
    double *d_mass,*d_diag,*d_r,*d_z,*d_p,*d_q,*d_delta,*d_partial,*d_scalar;
    cudaMalloc(&d_x,bytes); cudaMalloc(&d_b,bytes); cudaMalloc(&d_rho,bytes); cudaMalloc(&d_cv,bytes);
    cudaMalloc(&d_mu,bytes); cudaMalloc(&d_k,bytes); cudaMalloc(&d_rd,bytes); cudaMalloc(&d_gas,bytes);
    cudaMalloc(&d_mass,bytes); cudaMalloc(&d_diag,bytes); cudaMalloc(&d_r,bytes); cudaMalloc(&d_z,bytes);
    cudaMalloc(&d_p,bytes); cudaMalloc(&d_q,bytes); cudaMalloc(&d_delta,bytes);
    cudaMalloc(&d_partial, nblk * sizeof(double)); cudaMalloc(&d_scalar, sizeof(double));
    cudaMemcpy(d_x,x,bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(d_b,b,bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(d_rho,rho,bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(d_cv,cv,bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(d_mu,mu,bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(d_k,kk,bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(d_rd,rhod,bytes,cudaMemcpyHostToDevice);
    cudaMemcpy(d_gas,gas,bytes,cudaMemcpyHostToDevice);
    long gblk = (ncell + DTPB - 1) / DTPB;

    // mass; diag = mass - wq*diag(L)  (apply_linear x=d_x just for its diag).
    kd_fill_mass<<<gblk,DTPB>>>(comp,n_r,n_z,r_min,dr,dz,d_rho,d_cv,d_gas,d_mass);
    kd_apply_linear<<<gblk,DTPB>>>(d_x,comp,n_r,n_z,r_min,dr,dz,d_mu,d_k,d_rd,d_gas,1,d_q,d_diag);
    kd_diag_finish<<<gblk,DTPB>>>(d_diag,d_mass,d_gas,wq,ncell);

    // delta = 0; r = b; z = r/diag; p = z
    cudaMemset(d_delta,0,bytes);
    kd_copy<<<gblk,DTPB>>>(d_r,d_b,ncell);
    double b_norm2 = dev_dot(d_b,d_b,d_gas,ncell,d_partial,d_scalar,nblk);
    double eps2 = EPS_CG_RESID * EPS_CG_RESID * b_norm2;
    kd_precond<<<gblk,DTPB>>>(d_z,d_r,d_diag,ncell);
    kd_copy<<<gblk,DTPB>>>(d_p,d_z,ncell);
    double rz = dev_dot(d_r,d_z,d_gas,ncell,d_partial,d_scalar,nblk);
    double r_norm2 = dev_dot(d_r,d_r,d_gas,ncell,d_partial,d_scalar,nblk);
    int iters = 0;
    while (true) {
        if (fixed_iters > 0) { if (iters >= fixed_iters) break; }
        else if (!(iters < N_CG_ITERS_MAX && r_norm2 > eps2 && rz > 0.0)) break;
        // q = (mass - wq*L)·p
        kd_apply_linear<<<gblk,DTPB>>>(d_p,comp,n_r,n_z,r_min,dr,dz,d_mu,d_k,d_rd,d_gas,0,d_q,nullptr);
        kd_q_finish<<<gblk,DTPB>>>(d_q,d_mass,d_p,d_gas,wq,ncell);
        double pq = dev_dot(d_p,d_q,d_gas,ncell,d_partial,d_scalar,nblk);
        if (pq != pq || pq <= 0.0) break;               // NaN (self-≠) or SPD breakdown
        double alpha = rz / pq;
        kd_axpy<<<gblk,DTPB>>>(d_delta, alpha, d_p, ncell);
        kd_axpy<<<gblk,DTPB>>>(d_r, -alpha, d_q, ncell);
        kd_precond<<<gblk,DTPB>>>(d_z,d_r,d_diag,ncell);
        double rz_new = dev_dot(d_r,d_z,d_gas,ncell,d_partial,d_scalar,nblk);
        double beta = rz_new / rz;
        rz = rz_new;
        kd_p_update<<<gblk,DTPB>>>(d_p,d_z,beta,ncell);
        r_norm2 = dev_dot(d_r,d_r,d_gas,ncell,d_partial,d_scalar,nblk);
        iters++;
    }
    // resid = sqrt(r_norm2 / max(b_norm2, MIN_POSITIVE)); x += delta on gas.
    double denom = b_norm2 > 2.2250738585072014e-308 ? b_norm2 : 2.2250738585072014e-308;
    double resid = sqrt(r_norm2 / denom);
    kd_x_add_delta_gas<<<gblk,DTPB>>>(d_x, d_delta, d_gas, ncell);   // δ over gas cells (cg_solve)
    cudaDeviceSynchronize();
    cudaMemcpy(x, d_x, bytes, cudaMemcpyDeviceToHost);
    *iters_out = iters;
    *resid_out = resid;
    int fail = ((b_norm2 != b_norm2) || (b_norm2 > 0.0 && ((resid != resid) || resid > EPS_CG_RESID))) ? 1 : 0;

    cudaFree(d_x); cudaFree(d_b); cudaFree(d_rho); cudaFree(d_cv); cudaFree(d_mu); cudaFree(d_k);
    cudaFree(d_rd); cudaFree(d_gas); cudaFree(d_mass); cudaFree(d_diag); cudaFree(d_r); cudaFree(d_z);
    cudaFree(d_p); cudaFree(d_q); cudaFree(d_delta); cudaFree(d_partial); cudaFree(d_scalar);
    return fail;
}
