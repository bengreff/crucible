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
#define NC 7       // NCOMP (rho, m_r, m_theta, m_z, rho_e, rho_c, rho_b) — the
                   // affine class-D rate's dense stride (S13c assemble kernel)

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

// =====================================================================
// S13c — the class-D diffusion FORCING: fill_lag_gradients + assemble_rates
// (the affine viscous-stress physics that BUILDS the CG's RHS b). This is
// the intricate cross-term core of F_visc — compressible viscous stress
// (τ_rr/τ_zz/τ_rz/τ_θθ), Fourier conduction (k∇T), species diffusion
// (ρD∇C) + its enthalpy flux — bit-for-formula from gas_diffusion.rs's
// fill_lag_gradients (~L786) and assemble_rates (~L1054).
//
// SCOPE: N_theta=1 box (kappa=aperture=1), interior cells only (the compared
// set is ≥3 from every edge so the whole lag-gradient + face stencil is real
// interior data), free (FreeSlip/Adiabatic/ZeroFlux) BCs — a domain-edge face
// contributes nothing under those BCs, so "skip a face whose neighbour is
// out of domain" reproduces the CPU exactly for this fixture (real NoSlip/
// Robin walls are the S13c BC leg). Conserved slots: I_RHO0 I_MR1 I_MT2 I_MZ3
// I_EN4 I_RC5 I_RB6. Determinism: gather-only, one writer per cell.
// ---------------------------------------------------------------------

// Cell-centred lag velocity gradients (central where both r/z gas neighbours
// exist, one-sided at a domain edge, 0 if neither — fill_lag_gradients). Box:
// "gas neighbour exists" ⟺ "in domain". At N_theta=1 only these four (of
// lag u_r, u_z) are needed; dom_* / *_dth are the N_theta>1 legs.
__global__ void kd_lag_grads(const double* __restrict__ lag_ur, const double* __restrict__ lag_uz,
                             int n_r, int n_z, double dr, double dz, const double* __restrict__ gas,
                             double* __restrict__ dur_dr, double* __restrict__ dur_dz,
                             double* __restrict__ duz_dr, double* __restrict__ duz_dz) {
    long c = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long ncell = (long)n_r * n_z;
    if (c >= ncell) return;
    dur_dr[c] = dur_dz[c] = duz_dr[c] = duz_dz[c] = 0.0;
    if (gas[c] == 0.0) return;
    int i_r = (int)(c / n_z), i_z = (int)(c % n_z);
    // central/one-sided along a direction given the lo/hi neighbour presence.
    #define DERIV(FLD, LOOK, HIOK, LO, HI, CC, DD) \
        ( (LOOK && HIOK) ? (FLD[HI] - FLD[LO]) / (2.0*(DD)) \
        : (HIOK) ? (FLD[HI] - FLD[CC]) / (DD) \
        : (LOOK) ? (FLD[CC] - FLD[LO]) / (DD) : 0.0 )
    bool rlo = i_r > 0,        rhi = i_r < n_r - 1;
    bool zlo = i_z > 0,        zhi = i_z < n_z - 1;
    long rL = (long)(i_r-1)*n_z + i_z, rH = (long)(i_r+1)*n_z + i_z;
    long zL = (long)i_r*n_z + (i_z-1), zH = (long)i_r*n_z + (i_z+1);
    // a domain neighbour must also be gas (box ⇒ always true, but keep exact)
    rlo = rlo && gas[rL] != 0.0;  rhi = rhi && gas[rH] != 0.0;
    zlo = zlo && gas[zL] != 0.0;  zhi = zhi && gas[zH] != 0.0;
    dur_dr[c] = DERIV(lag_ur, rlo, rhi, rL, rH, c, dr);
    dur_dz[c] = DERIV(lag_ur, zlo, zhi, zL, zH, c, dz);
    duz_dr[c] = DERIV(lag_uz, rlo, rhi, rL, rH, c, dr);
    duz_dz[c] = DERIV(lag_uz, zlo, zhi, zL, zH, c, dz);
    #undef DERIV
}

// The affine class-D rate (assemble_rates), N_theta=1 box interior. sol =
// current iterate (u_r,ω,u_z,T,C); lag_ur feeds e_θθ; the lag gradients feed
// the cross/compressible pieces. Writes the 5 diffusion rates per cell
// (mass slot 0, burn slot 0) in conserved-density-rate units.
__global__ void kd_assemble_rates(const double* __restrict__ s_ur, const double* __restrict__ s_om,
                                  const double* __restrict__ s_uz, const double* __restrict__ s_tt,
                                  const double* __restrict__ s_cc, const double* __restrict__ lag_ur,
                                  const double* __restrict__ dur_dr, const double* __restrict__ dur_dz,
                                  const double* __restrict__ duz_dr, const double* __restrict__ duz_dz,
                                  const double* __restrict__ mu, const double* __restrict__ kk,
                                  const double* __restrict__ rhod, const double* __restrict__ dhdz,
                                  const double* __restrict__ gas, int n_r, int n_z,
                                  double r_min, double dr, double dz, double* __restrict__ rate) {
    long c = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long ncell = (long)n_r * n_z;
    if (c >= ncell) return;
    for (int k = 0; k < NC; k++) rate[c*NC+k] = 0.0;
    if (gas[c] == 0.0) return;
    int i_r = (int)(c / n_z), i_z = (int)(c % n_z);
    double rbar = d_r_center(r_min, dr, i_r);
    double vol  = d_cell_volume(r_min, dr, dz, i_r);
    double kv   = vol;                                   // kappa = 1
    double ur_c = s_ur[c], om_c = s_om[c], uz_c = s_uz[c], tt_c = s_tt[c], cc_c = s_cc[c];
    double tot_mr = 0.0, tot_mz = 0.0, tot_en = 0.0, tot_rc = 0.0, tot_lam = 0.0;
    // 4 faces: 0=r-,1=r+,2=z-,3=z+
    for (int fc = 0; fc < 4; fc++) {
        bool radial = fc < 2;
        double s = (fc & 1) ? 1.0 : -1.0;
        int nr = i_r, nz = i_z;
        if (fc == 0) nr = i_r - 1; else if (fc == 1) nr = i_r + 1;
        else if (fc == 2) nz = i_z - 1; else nz = i_z + 1;
        if (nr < 0 || nr >= n_r || nz < 0 || nz >= n_z) continue;   // free BC ⇒ no term
        long nc = (long)nr * n_z + nz;
        if (gas[nc] == 0.0) continue;
        double d = radial ? dr : dz;
        double r_face = radial ? d_face_radius(r_min, dr, (fc == 1) ? i_r + 1 : i_r) : rbar;
        double area = radial ? d_area_r(r_min, dr, dz, i_r, (fc == 1)) : d_area_z(r_min, dr, i_r);
        double mu_f  = 0.5 * (mu[c]   + mu[nc]);
        double k_f   = 0.5 * (kk[c]   + kk[nc]);
        double d_f   = 0.5 * (rhod[c] + rhod[nc]);
        double dh_f  = 0.5 * (dhdz[c] + dhdz[nc]);
        double two_thirds_mu = (2.0 / 3.0) * mu_f;
        double g_ur = s * (s_ur[nc] - ur_c) / d;
        double g_om = s * (s_om[nc] - om_c) / d;
        double g_uz = s * (s_uz[nc] - uz_c) / d;
        double g_tt = s * (s_tt[nc] - tt_c) / d;
        double g_cc = s * (s_cc[nc] - cc_c) / d;
        // e_θθ face avg: neighbour uses its own ring radius on r-faces, r̄ on z.
        double r_n = radial ? d_r_center(r_min, dr, nr) : rbar;
        double e_thth_f = 0.5 * (lag_ur[c] / rbar + lag_ur[nc] / r_n);
        double e_rr_f   = 0.5 * (dur_dr[c] + dur_dr[nc]);
        double e_zz_f   = 0.5 * (duz_dz[c] + duz_dz[nc]);
        double dur_dz_f = 0.5 * (dur_dz[c] + dur_dz[nc]);
        double duz_dr_f = 0.5 * (duz_dr[c] + duz_dr[nc]);
        double ur_f = 0.5 * (ur_c + s_ur[nc]);
        double uz_f = 0.5 * (uz_c + s_uz[nc]);
        double ut_f = radial ? 0.5 * (om_c * rbar + s_om[nc] * r_n)
                             : 0.5 * (om_c + s_om[nc]) * rbar;
        double f_mr, f_mz, tau_nn, tau_rz, tau_th;
        if (radial) {
            double tau_rr = (4.0/3.0) * mu_f * g_ur - two_thirds_mu * (e_thth_f + e_zz_f);
            double t_rz   = mu_f * (g_uz + dur_dz_f);
            double t_rth  = mu_f * r_face * g_om;
            f_mr = area * tau_rr; f_mz = area * t_rz;
            tau_nn = tau_rr; tau_rz = t_rz; tau_th = t_rth;
        } else {
            double tau_zz = (4.0/3.0) * mu_f * g_uz - two_thirds_mu * (e_rr_f + e_thth_f);
            double t_rz   = mu_f * (g_ur + duz_dr_f);
            double t_thz  = mu_f * rbar * g_om;
            f_mr = area * t_rz; f_mz = area * tau_zz;
            tau_nn = tau_zz; tau_rz = t_rz; tau_th = t_thz;
        }
        tot_mr += s * f_mr;
        tot_mz += s * f_mz;
        double f_lam = area * r_face * tau_th;
        tot_lam += s * f_lam;
        double g_e = (radial ? (ur_f * tau_nn + ut_f * tau_th + uz_f * tau_rz)
                             : (ur_f * tau_rz + ut_f * tau_th + uz_f * tau_nn))
                   + k_f * g_tt + d_f * dh_f * g_cc;
        tot_en += s * area * g_e;
        tot_rc += s * area * d_f * g_cc;
    }
    // −τ_θθ/r volume source of r-momentum (SOLV-1 §3.3), metric-consistent 1/r̄.
    double geo = (d_area_r(r_min, dr, dz, i_r, true) - d_area_r(r_min, dr, dz, i_r, false)) / vol;
    double e_rr_c = dur_dr[c], e_zz_c = duz_dz[c], mu_c = mu[c];
    double tau_thth = (4.0/3.0) * mu_c * (ur_c / rbar) - (2.0/3.0) * mu_c * (e_rr_c + e_zz_c);
    double src_mr = -tau_thth * geo * kv;
    tot_mr += src_mr;
    double inv_kv = 1.0 / kv;
    rate[c*NC + 1] = tot_mr * inv_kv;                    // I_MR
    rate[c*NC + 2] = tot_lam / (rbar * kv);              // I_MT (λ → ρu_θ rate)
    rate[c*NC + 3] = tot_mz * inv_kv;                    // I_MZ
    rate[c*NC + 4] = tot_en * inv_kv;                    // I_EN
    rate[c*NC + 5] = tot_rc * inv_kv;                    // I_RC
}

// FFI: one affine class-D rate evaluation (fill_lag_gradients + assemble_rates)
// on a resident fixture. All inputs dense [ncell]; rate dense [ncell*NC].
extern "C" void gpu_class_d_assemble(const double* s_ur, const double* s_om, const double* s_uz,
                                     const double* s_tt, const double* s_cc, const double* lag_ur,
                                     const double* lag_uz, const double* mu, const double* kk,
                                     const double* rhod, const double* dhdz, const double* gas,
                                     int n_r, int n_z, double r_min, double dr, double dz,
                                     double* rate) {
    long ncell = (long)n_r * n_z;
    size_t nb = ncell * sizeof(double);
    double *d_ur,*d_om,*d_uz,*d_tt,*d_cc,*d_lur,*d_luz,*d_mu,*d_k,*d_rd,*d_dh,*d_gas;
    double *d_grr,*d_grz,*d_gzr,*d_gzz,*d_rate;
    cudaMalloc(&d_ur,nb); cudaMalloc(&d_om,nb); cudaMalloc(&d_uz,nb); cudaMalloc(&d_tt,nb);
    cudaMalloc(&d_cc,nb); cudaMalloc(&d_lur,nb); cudaMalloc(&d_luz,nb); cudaMalloc(&d_mu,nb);
    cudaMalloc(&d_k,nb); cudaMalloc(&d_rd,nb); cudaMalloc(&d_dh,nb); cudaMalloc(&d_gas,nb);
    cudaMalloc(&d_grr,nb); cudaMalloc(&d_grz,nb); cudaMalloc(&d_gzr,nb); cudaMalloc(&d_gzz,nb);
    cudaMalloc(&d_rate, ncell*NC*sizeof(double));
    cudaMemcpy(d_ur,s_ur,nb,cudaMemcpyHostToDevice); cudaMemcpy(d_om,s_om,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_uz,s_uz,nb,cudaMemcpyHostToDevice); cudaMemcpy(d_tt,s_tt,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_cc,s_cc,nb,cudaMemcpyHostToDevice); cudaMemcpy(d_lur,lag_ur,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_luz,lag_uz,nb,cudaMemcpyHostToDevice); cudaMemcpy(d_mu,mu,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_k,kk,nb,cudaMemcpyHostToDevice); cudaMemcpy(d_rd,rhod,nb,cudaMemcpyHostToDevice);
    cudaMemcpy(d_dh,dhdz,nb,cudaMemcpyHostToDevice); cudaMemcpy(d_gas,gas,nb,cudaMemcpyHostToDevice);
    long gblk = (ncell + DTPB - 1) / DTPB;
    kd_lag_grads<<<gblk,DTPB>>>(d_lur,d_luz,n_r,n_z,dr,dz,d_gas,d_grr,d_grz,d_gzr,d_gzz);
    kd_assemble_rates<<<gblk,DTPB>>>(d_ur,d_om,d_uz,d_tt,d_cc,d_lur,d_grr,d_grz,d_gzr,d_gzz,
                                     d_mu,d_k,d_rd,d_dh,d_gas,n_r,n_z,r_min,dr,dz,d_rate);
    cudaDeviceSynchronize();
    cudaMemcpy(rate, d_rate, ncell*NC*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_ur); cudaFree(d_om); cudaFree(d_uz); cudaFree(d_tt); cudaFree(d_cc); cudaFree(d_lur);
    cudaFree(d_luz); cudaFree(d_mu); cudaFree(d_k); cudaFree(d_rd); cudaFree(d_dh); cudaFree(d_gas);
    cudaFree(d_grr); cudaFree(d_grz); cudaFree(d_gzr); cudaFree(d_gzz); cudaFree(d_rate);
}

// =====================================================================
// S13c — the FULL RESIDENT class-D gas Picard iterate: compose the FORCING
// (assemble_rates) + the SOLVER (cg_solve) into one on-device implicit
// diffusion sweep, bit-for-formula from the SDC inner block (sdc.rs ~L1158-
// 1225): fill_lag_gradients → assemble(sol,lag) → {Ur,Uz,Om: fill_rhs, mass,
// cg} → re-assemble(sol',lag) → {T,C: fill_rhs, mass, cg}. The fields stay
// RESIDENT across the whole iterate (operands + dstage + CG work vectors live
// on-device; only the O(1) CG scalars round-trip). This is the resident
// implicit-diffusion solve — the composition the class-D residency exists for.
// ---------------------------------------------------------------------

// fill_gas_rhs (sdc.rs ~L2231): b = wqnew·kv·(dstage−dlag)[k]  (ω: ×r̄;
// T: + ρ·kv·(ke_base − ke_new), the dissipation the momentum solves moved
// into internal energy). k = the comp's conserved slot. ke_new from the
// CURRENT sol velocities (post-velocity-solve for T).
__global__ void kd_fill_gas_rhs(int comp, double wqnew, int n_r, int n_z,
                                double r_min, double dr, double dz,
                                const double* __restrict__ dstage, const double* __restrict__ dlag,
                                const double* __restrict__ s_ur, const double* __restrict__ s_om,
                                const double* __restrict__ s_uz, const double* __restrict__ rho,
                                const double* __restrict__ ke_base, const double* __restrict__ gas,
                                double* __restrict__ b) {
    long c = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long ncell = (long)n_r * n_z;
    if (c >= ncell) return;
    if (gas[c] == 0.0) { b[c] = 0.0; return; }
    int i_r = (int)(c / n_z);
    double kv = d_cell_volume(r_min, dr, dz, i_r);       // kappa = 1
    // comp → conserved slot: Ur=1 Uz=3 Om=2 T=4 C=5
    int k = (comp==COMP_UR)?1 : (comp==COMP_UZ)?3 : (comp==COMP_OM)?2 : (comp==COMP_T)?4 : 5;
    double ddiff = dstage[c*NC+k] - dlag[c*NC+k];
    double val = (comp==COMP_OM) ? wqnew * d_r_center(r_min, dr, i_r) * kv * ddiff
                                 : wqnew * kv * ddiff;
    if (comp == COMP_T) {
        double r = d_r_center(r_min, dr, i_r);
        double ut = s_om[c] * r;
        double ke_new = 0.5 * (s_ur[c]*s_ur[c] + ut*ut + s_uz[c]*s_uz[c]);
        val += rho[c] * kv * (ke_base[c] - ke_new);
    }
    b[c] = val;
}

// ke_base = 0.5(u_r² + u_θ² + u_z²) at the ORIGINAL operands (captured before
// the velocity solves — derive_gas_operands' ke_base).
__global__ void kd_ke_base(int n_r, int n_z, double r_min, double dr,
                           const double* __restrict__ s_ur, const double* __restrict__ s_om,
                           const double* __restrict__ s_uz, const double* __restrict__ gas,
                           double* __restrict__ ke_base) {
    long c = (long)blockIdx.x * blockDim.x + threadIdx.x;
    long ncell = (long)n_r * n_z;
    if (c >= ncell) return;
    if (gas[c] == 0.0) { ke_base[c] = 0.0; return; }
    int i_r = (int)(c / n_z);
    double ut = s_om[c] * d_r_center(r_min, dr, i_r);
    ke_base[c] = 0.5 * (s_ur[c]*s_ur[c] + ut*ut + s_uz[c]*s_uz[c]);
}

// Resident CG for one component on device pointers (no malloc; the SDC-inner
// reuse of cg_solve). `x` (= the sol component buffer) is updated in place:
// x += δ over gas cells. Own data-dependent termination (N_CG_ITERS_MAX cap),
// bit-for-formula from cg_solve. `wnb` structures shared across the iterate.
struct CgWork { double *mass,*diag,*r,*z,*p,*q,*delta,*partial,*scalar; int nblk; };
static void run_cg_resident(double* d_x, const double* d_b, const double* d_rho, const double* d_cv,
                            const double* d_mu, const double* d_k, const double* d_rd, const double* d_gas,
                            int comp, int n_r, int n_z, double r_min, double dr, double dz,
                            double wq, long ncell, long gblk, CgWork w) {
    kd_fill_mass<<<gblk,DTPB>>>(comp,n_r,n_z,r_min,dr,dz,d_rho,d_cv,d_gas,w.mass);
    kd_apply_linear<<<gblk,DTPB>>>(d_x,comp,n_r,n_z,r_min,dr,dz,d_mu,d_k,d_rd,d_gas,1,w.q,w.diag);
    kd_diag_finish<<<gblk,DTPB>>>(w.diag,w.mass,d_gas,wq,ncell);
    cudaMemset(w.delta,0,ncell*sizeof(double));
    kd_copy<<<gblk,DTPB>>>(w.r,d_b,ncell);
    double b_norm2 = dev_dot(d_b,d_b,d_gas,ncell,w.partial,w.scalar,w.nblk);
    double eps2 = EPS_CG_RESID * EPS_CG_RESID * b_norm2;
    kd_precond<<<gblk,DTPB>>>(w.z,w.r,w.diag,ncell);
    kd_copy<<<gblk,DTPB>>>(w.p,w.z,ncell);
    double rz = dev_dot(w.r,w.z,d_gas,ncell,w.partial,w.scalar,w.nblk);
    double r_norm2 = dev_dot(w.r,w.r,d_gas,ncell,w.partial,w.scalar,w.nblk);
    int iters = 0;
    while (iters < N_CG_ITERS_MAX && r_norm2 > eps2 && rz > 0.0) {
        kd_apply_linear<<<gblk,DTPB>>>(w.p,comp,n_r,n_z,r_min,dr,dz,d_mu,d_k,d_rd,d_gas,0,w.q,nullptr);
        kd_q_finish<<<gblk,DTPB>>>(w.q,w.mass,w.p,d_gas,wq,ncell);
        double pq = dev_dot(w.p,w.q,d_gas,ncell,w.partial,w.scalar,w.nblk);
        if (pq != pq || pq <= 0.0) break;
        double alpha = rz / pq;
        kd_axpy<<<gblk,DTPB>>>(w.delta, alpha, w.p, ncell);
        kd_axpy<<<gblk,DTPB>>>(w.r, -alpha, w.q, ncell);
        kd_precond<<<gblk,DTPB>>>(w.z,w.r,w.diag,ncell);
        double rz_new = dev_dot(w.r,w.z,d_gas,ncell,w.partial,w.scalar,w.nblk);
        double beta = rz_new / rz; rz = rz_new;
        kd_p_update<<<gblk,DTPB>>>(w.p,w.z,beta,ncell);
        r_norm2 = dev_dot(w.r,w.r,d_gas,ncell,w.partial,w.scalar,w.nblk);
        iters++;
    }
    kd_x_add_delta_gas<<<gblk,DTPB>>>(d_x, w.delta, d_gas, ncell);
}

// FFI: one full resident class-D gas iterate. sol operands (rho,ur,om,uz,tt,cc)
// are inout (initial in, updated out); lag (ur,uz), transport (mu,k,rhod,dhdz,
// cv), dlag [ncell*NC], gas are in. wqnew = the SDC implicit weight. All fields
// resident across the iterate; returns worst CG iters (diagnostic).
extern "C" int gpu_class_d_iterate(double* sol_rho, double* sol_ur, double* sol_om, double* sol_uz,
                                   double* sol_tt, double* sol_cc, const double* lag_ur,
                                   const double* lag_uz, const double* mu, const double* kk,
                                   const double* rhod, const double* dhdz, const double* cv,
                                   const double* dlag, const double* gas, int n_r, int n_z,
                                   double r_min, double dr, double dz, double wqnew) {
    long ncell = (long)n_r * n_z;
    size_t nb = ncell*sizeof(double), nbc = ncell*NC*sizeof(double);
    long gblk = (ncell + DTPB - 1)/DTPB;
    int nblk = (int)gblk;
    // resident operands + transport + lag + dlag + gas
    double *d_rho,*d_ur,*d_om,*d_uz,*d_tt,*d_cc,*d_lur,*d_luz,*d_mu,*d_k,*d_rd,*d_dh,*d_cv,*d_dlag,*d_gas;
    // lag grads + dstage + ke_base + b + CG work
    double *d_grr,*d_grz,*d_gzr,*d_gzz,*d_dstage,*d_ke,*d_b;
    double *d_mass,*d_diag,*d_r,*d_z,*d_p,*d_q,*d_delta,*d_partial,*d_scalar;
    #define M(p,sz) cudaMalloc(&p,sz)
    M(d_rho,nb);M(d_ur,nb);M(d_om,nb);M(d_uz,nb);M(d_tt,nb);M(d_cc,nb);
    M(d_lur,nb);M(d_luz,nb);M(d_mu,nb);M(d_k,nb);M(d_rd,nb);M(d_dh,nb);M(d_cv,nb);
    M(d_dlag,nbc);M(d_gas,nb);M(d_grr,nb);M(d_grz,nb);M(d_gzr,nb);M(d_gzz,nb);
    M(d_dstage,nbc);M(d_ke,nb);M(d_b,nb);
    M(d_mass,nb);M(d_diag,nb);M(d_r,nb);M(d_z,nb);M(d_p,nb);M(d_q,nb);M(d_delta,nb);
    M(d_partial,nblk*sizeof(double));M(d_scalar,sizeof(double));
    #undef M
    #define CPY(d,s,sz) cudaMemcpy(d,s,sz,cudaMemcpyHostToDevice)
    CPY(d_rho,sol_rho,nb);CPY(d_ur,sol_ur,nb);CPY(d_om,sol_om,nb);CPY(d_uz,sol_uz,nb);
    CPY(d_tt,sol_tt,nb);CPY(d_cc,sol_cc,nb);CPY(d_lur,lag_ur,nb);CPY(d_luz,lag_uz,nb);
    CPY(d_mu,mu,nb);CPY(d_k,kk,nb);CPY(d_rd,rhod,nb);CPY(d_dh,dhdz,nb);CPY(d_cv,cv,nb);
    CPY(d_dlag,dlag,nbc);CPY(d_gas,gas,nb);
    #undef CPY
    CgWork w = {d_mass,d_diag,d_r,d_z,d_p,d_q,d_delta,d_partial,d_scalar,nblk};

    // ke_base from the ORIGINAL velocities; lag gradients from the fixed lag.
    kd_ke_base<<<gblk,DTPB>>>(n_r,n_z,r_min,dr,d_ur,d_om,d_uz,d_gas,d_ke);
    kd_lag_grads<<<gblk,DTPB>>>(d_lur,d_luz,n_r,n_z,dr,dz,d_gas,d_grr,d_grz,d_gzr,d_gzz);

    // assemble(sol,lag) → dstage; solve the velocity block {Ur,Uz,Om}.
    kd_assemble_rates<<<gblk,DTPB>>>(d_ur,d_om,d_uz,d_tt,d_cc,d_lur,d_grr,d_grz,d_gzr,d_gzz,
                                     d_mu,d_k,d_rd,d_dh,d_gas,n_r,n_z,r_min,dr,dz,d_dstage);
    int comps_v[3] = {COMP_UR, COMP_UZ, COMP_OM};
    double* xs_v[3] = {d_ur, d_uz, d_om};
    for (int i=0;i<3;i++) {
        kd_fill_gas_rhs<<<gblk,DTPB>>>(comps_v[i],wqnew,n_r,n_z,r_min,dr,dz,d_dstage,d_dlag,
                                       d_ur,d_om,d_uz,d_rho,d_ke,d_gas,d_b);
        run_cg_resident(xs_v[i],d_b,d_rho,d_cv,d_mu,d_k,d_rd,d_gas,comps_v[i],
                        n_r,n_z,r_min,dr,dz,wqnew,ncell,gblk,w);
    }
    // re-assemble at the NEW velocities (lag grads unchanged); solve {T,C}.
    kd_assemble_rates<<<gblk,DTPB>>>(d_ur,d_om,d_uz,d_tt,d_cc,d_lur,d_grr,d_grz,d_gzr,d_gzz,
                                     d_mu,d_k,d_rd,d_dh,d_gas,n_r,n_z,r_min,dr,dz,d_dstage);
    int comps_t[2] = {COMP_T, COMP_C};
    double* xs_t[2] = {d_tt, d_cc};
    for (int i=0;i<2;i++) {
        kd_fill_gas_rhs<<<gblk,DTPB>>>(comps_t[i],wqnew,n_r,n_z,r_min,dr,dz,d_dstage,d_dlag,
                                       d_ur,d_om,d_uz,d_rho,d_ke,d_gas,d_b);
        run_cg_resident(xs_t[i],d_b,d_rho,d_cv,d_mu,d_k,d_rd,d_gas,comps_t[i],
                        n_r,n_z,r_min,dr,dz,wqnew,ncell,gblk,w);
    }
    cudaDeviceSynchronize();
    // download updated sol (ur,uz,om,tt,cc; rho unchanged).
    cudaMemcpy(sol_ur,d_ur,nb,cudaMemcpyDeviceToHost);
    cudaMemcpy(sol_uz,d_uz,nb,cudaMemcpyDeviceToHost);
    cudaMemcpy(sol_om,d_om,nb,cudaMemcpyDeviceToHost);
    cudaMemcpy(sol_tt,d_tt,nb,cudaMemcpyDeviceToHost);
    cudaMemcpy(sol_cc,d_cc,nb,cudaMemcpyDeviceToHost);
    double* all[] = {d_rho,d_ur,d_om,d_uz,d_tt,d_cc,d_lur,d_luz,d_mu,d_k,d_rd,d_dh,d_cv,d_dlag,
        d_gas,d_grr,d_grz,d_gzr,d_gzz,d_dstage,d_ke,d_b,d_mass,d_diag,d_r,d_z,d_p,d_q,d_delta,
        d_partial,d_scalar};
    for (double* q : all) cudaFree(q);
    return 0;
}
