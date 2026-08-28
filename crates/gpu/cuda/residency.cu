// crucible-gpu — S13 residency: the class-A (explicit hyperbolic) SDC RHS on
// device, generalizing the S12 spike's z-only-uniform sweep to the REAL
// operator: PPM reconstruction + HLLC-Batten flux + the exact cylindrical
// face_radius metric (r-sweep area-weighted, z-sweep metric-ratio) + the
// SOLV-1 §3.3 geometric sources, all bit-for-formula identical to
// crucible_solvers::euler (recon.rs / hllc.rs / mod.rs).
//
// SCOPE of this kernel set (S13, the adiabatic-hyperbolic residency core):
//   N_theta = 1, box world (kappa = aperture = 1), GammaLaw EOS.
// The cross-check compares INTERIOR cells (>= NGHOST=3 from every domain
// edge), whose compact PPM stencil is all real interior data — so the BC /
// mixed-N_theta reflux / axis-parity / SRD machinery (carried to S13b) does
// not enter the compared cells. Cut apertures + real HDF5 TableEos + the
// class-D diffusion CG are S13b (recorded in the SESSION_LOG split).
//
// Determinism (META-1 §2.5): gather-only, one writer per cell, no atomics,
// fixed control flow. Same-build reruns are bit-identical; CPU<->GPU agree to
// the declared ECT (FMA-order). All device functions are prefixed `rz_` so
// they never collide with hllc_kernel.cu's symbols at link time.
#include <cuda_runtime.h>

#define NP 9   // NPRIM  (rho,u_r,u_theta,u_z,p,c,b,e,Gamma1)
#define NC 7   // NCOMP  (rho, m_r, m_theta, m_z, rho_e, rho_c, rho_b)
#define NGH 3  // NGHOST

// TAU = 2*pi as an f64 literal (matches Rust std::f64::consts::TAU to <1 ulp;
// the cross-check is ECT, not bit-exact CPU<->GPU).
__device__ __constant__ double RZ_TAU = 6.283185307179586;

// ---- GammaLaw closures (euler/mod.rs GammaLaw) --------------------------
__device__ __forceinline__ double rz_total_energy(const double* w, double g) {
    return w[4] / (g - 1.0) + 0.5 * w[0] * (w[1]*w[1] + w[2]*w[2] + w[3]*w[3]);
}
__device__ __forceinline__ double rz_sound(const double* w, double g) {
    return sqrt(g * w[4] / w[0]);
}
// Conserved -> primitive (fills p; aux slots e,Gamma1 unused by GammaLaw).
// Returns 0 on non-physical (rho<=0 or p<=0), matching the CPU halt predicate.
__device__ int rz_prim(const double* u, double g, double* w) {
    double rho = u[0];
    if (!isfinite(rho) || rho <= 0.0) return 0;
    double inv = 1.0 / rho;
    double ur = u[1]*inv, ut = u[2]*inv, uz = u[3]*inv;
    double ke = 0.5 * rho * (ur*ur + ut*ut + uz*uz);
    double p = (g - 1.0) * (u[4] - ke);
    if (!isfinite(p) || p <= 0.0) return 0;
    w[0]=rho; w[1]=ur; w[2]=ut; w[3]=uz; w[4]=p;
    w[5]=u[5]*inv; w[6]=u[6]*inv; w[7]=0.0; w[8]=0.0;
    return 1;
}

// ---- HLLC-Batten (euler/hllc.rs), GammaLaw, gamma passed in ------------
__device__ void rz_physical_flux(const double* w, int n, double g, double* f) {
    double rho=w[0], p=w[4], un=w[n], m=rho*un, e=rz_total_energy(w,g);
    f[0]=m; f[1]=m*w[1]; f[2]=m*w[2]; f[3]=m*w[3]; f[n]+=p;
    f[4]=un*(e+p); f[5]=m*w[5]; f[6]=m*w[6];
}
__device__ void rz_hllc(const double* wl, const double* wr, int n, double g, double* f) {
    double rho_l=wl[0], p_l=wl[4], rho_r=wr[0], p_r=wr[4];
    double un_l=wl[n], un_r=wr[n];
    double c_l=rz_sound(wl,g), c_r=rz_sound(wr,g);
    double sql=sqrt(rho_l), sqr=sqrt(rho_r), inv=1.0/(sql+sqr);
    double u1=(sql*wl[1]+sqr*wr[1])*inv;
    double u2=(sql*wl[2]+sqr*wr[2])*inv;
    double u3=(sql*wl[3]+sqr*wr[3])*inv;
    double urn=(n==1?u1:(n==2?u2:u3));
    double h_l=(rz_total_energy(wl,g)+p_l)/rho_l;
    double h_r=(rz_total_energy(wr,g)+p_r)/rho_r;
    double h_roe=(sql*h_l+sqr*h_r)*inv;
    double q2=u1*u1+u2*u2+u3*u3;
    double arg=(g-1.0)*(h_roe-0.5*q2);
    double c_roe=sqrt(arg>0.0?arg:0.0);
    double s_l=fmin(un_l-c_l, urn-c_roe);
    double s_r=fmax(un_r+c_r, urn+c_roe);
    double ml=rho_l*(s_l-un_l), mr=rho_r*(s_r-un_r);
    double s_m=(mr*un_r-ml*un_l+p_l-p_r)/(mr-ml);
    for (int k=0;k<NC;k++) f[k]=0.0;
    if (s_l>=0.0) { rz_physical_flux(wl,n,g,f); return; }
    if (s_r<=0.0) { rz_physical_flux(wr,n,g,f); return; }
    const double* w=(s_m>=0.0)?wl:wr; double s_k=(s_m>=0.0)?s_l:s_r;
    double rho=w[0], p=w[4], un=w[n], e=rz_total_energy(w,g);
    double p_star=rho*(un-s_k)*(un-s_m)+p;
    double u_k[NC]={rho, rho*w[1], rho*w[2], rho*w[3], e, rho*w[5], rho*w[6]};
    double fac=(s_k-un)/(s_k-s_m), rho_s=rho*fac;
    double u_s[NC];
    u_s[0]=rho_s; u_s[1]=rho_s*w[1]; u_s[2]=rho_s*w[2]; u_s[3]=rho_s*w[3];
    u_s[n]=rho_s*s_m;
    u_s[4]=fac*e+(p_star*s_m-p*un)/(s_k-s_m);
    u_s[5]=rho_s*w[5]; u_s[6]=rho_s*w[6];
    double fk[NC]; for(int k=0;k<NC;k++) fk[k]=0.0; rz_physical_flux(w,n,g,fk);
    for (int k=0;k<NC;k++) f[k]=fk[k]+s_k*(u_s[k]-u_k[k]);
}

// ---- PPM (euler/recon.rs) ----------------------------------------------
__device__ __forceinline__ double rz_mc_slope(double wm, double w0, double wp) {
    double dl=w0-wm, dr=wp-w0;
    if (dl*dr <= 0.0) return 0.0;
    double dc=0.5*(dl+dr);
    double a=fabs(dc); a=fmin(a, 2.0*fabs(dl)); a=fmin(a, 2.0*fabs(dr));
    return copysign(a, dc);
}
__device__ void rz_slope(const double* wm, const double* w0, const double* wp, double* s) {
    for (int k=0;k<NP;k++) s[k]=rz_mc_slope(wm[k], w0[k], wp[k]);
}
__device__ void rz_iface(const double* wa, const double* wb,
                         const double* sa, const double* sb, double* out) {
    for (int k=0;k<NP;k++) out[k]=0.5*(wa[k]+wb[k]) - (sb[k]-sa[k])/6.0;
}
// Monotonized left/right edge of a cell (CW84 eq. 1.10).
__device__ void rz_edge(const double* lo_if, const double* hi_if, const double* c,
                        double* elo, double* ehi) {
    for (int k=0;k<NP;k++) {
        double lo=lo_if[k], hi=hi_if[k], cc=c[k];
        if ((hi-cc)*(cc-lo) <= 0.0) { lo=cc; hi=cc; }
        else {
            double d=hi-lo, six=6.0*(cc-0.5*(lo+hi));
            if (d*six > d*d) lo=3.0*cc-2.0*hi;
            else if (d*six < -(d*d)) hi=3.0*cc-2.0*lo;
        }
        elo[k]=lo; ehi[k]=hi;
    }
}
// Faces bounding the centre cell p[3] of a 7-cell pencil p[0..6] (each NP wide).
// Writes the left face's (fl,fr) and the right face's (fl,fr).
__device__ void rz_pencil_faces(const double* p,
                                double* lfl, double* lfr, double* rfl, double* rfr) {
    double s1[NP],s2[NP],s3[NP],s4[NP],s5[NP];
    rz_slope(p+0*NP, p+1*NP, p+2*NP, s1);
    rz_slope(p+1*NP, p+2*NP, p+3*NP, s2);
    rz_slope(p+2*NP, p+3*NP, p+4*NP, s3);
    rz_slope(p+3*NP, p+4*NP, p+5*NP, s4);
    rz_slope(p+4*NP, p+5*NP, p+6*NP, s5);
    double if1[NP],if2[NP],if3[NP],if4[NP];
    rz_iface(p+1*NP, p+2*NP, s1, s2, if1);
    rz_iface(p+2*NP, p+3*NP, s2, s3, if2);
    rz_iface(p+3*NP, p+4*NP, s3, s4, if3);
    rz_iface(p+4*NP, p+5*NP, s4, s5, if4);
    double e2lo[NP],e2hi[NP],e3lo[NP],e3hi[NP],e4lo[NP],e4hi[NP];
    rz_edge(if1, if2, p+2*NP, e2lo, e2hi);   // cell p[2]
    rz_edge(if2, if3, p+3*NP, e3lo, e3hi);   // cell p[3] (centre)
    rz_edge(if3, if4, p+4*NP, e4lo, e4hi);   // cell p[4]
    for (int k=0;k<NP;k++) {
        lfl[k]=e2hi[k]; lfr[k]=e3lo[k];      // left  face: hi(p2) | lo(p3)
        rfl[k]=e3hi[k]; rfr[k]=e4lo[k];      // right face: hi(p3) | lo(p4)
    }
}

// ---- kernels -----------------------------------------------------------
// cell index = i_r * n_z + i_z ; cons/prim/rate strides NC/NP/NC.
__global__ void k_fill_prims(const double* __restrict__ cons, int ncell, double g,
                             double* __restrict__ prim, int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= ncell) return;
    double w[NP];
    if (!rz_prim(cons + c*NC, g, w)) { atomicExch(bad, 1); return; }
    for (int k=0;k<NP;k++) prim[c*NP+k]=w[k];
}

// STAGED rate kernels (S13 register-reduction lever, the brief's STEP 1c):
// one direction per kernel holds ONE 7-cell pencil, not both + the two flux
// sets at once. k_rate_r WRITES the r-sweep divergence + the geometric
// sources; k_rate_z ADDS the z-sweep divergence. One writer per kernel per
// cell, launched in sequence (no atomics, deterministic). Interior cells only
// (>= NGH from every edge); boundary cells are zeroed by k_rate_r (the S13b BC
// wave owns them) and left untouched by k_rate_z.
__global__ void k_rate_r(const double* __restrict__ prim, int n_r, int n_z,
                         double r_min, double dr, double z_min, double dz,
                         double g, double* __restrict__ rate) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    long ncell = (long)n_r*n_z;
    if (c >= ncell) return;
    int i_r = (int)(c / n_z), i_z = (int)(c % n_z);
    if (i_r < NGH || i_r >= n_r-NGH || i_z < NGH || i_z >= n_z-NGH) {
        for (int k=0;k<NC;k++) rate[c*NC+k]=0.0;
        return;
    }
    double p[7*NP];
    for (int m=0;m<7;m++) {
        long cc = (long)(i_r-3+m)*n_z + i_z;
        for (int k=0;k<NP;k++) p[m*NP+k]=prim[cc*NP+k];
    }
    double lfl[NP],lfr[NP],rfl[NP],rfr[NP];
    rz_pencil_faces(p, lfl, lfr, rfl, rfr);
    double fL[NC], fR[NC];
    rz_hllc(lfl, lfr, 1, g, fL);
    rz_hllc(rfl, rfr, 1, g, fR);
    double ri = r_min + (double)i_r*dr;
    double ro = r_min + (double)(i_r+1)*dr;
    double a_in  = ri * RZ_TAU * dz;
    double a_out = ro * RZ_TAU * dz;
    double vol   = 0.5*(ro*ro - ri*ri) * RZ_TAU * dz;
    double out[NC];
    for (int k=0;k<NC;k++) out[k] = (a_in*fL[k] - a_out*fR[k]) / vol;
    // geometric sources (add_sources, box world, no external src)
    double geo = (a_out - a_in) / vol;
    const double* wc = prim + c*NP;
    double rho=wc[0], ur=wc[1], ut=wc[2], pp=wc[4];
    out[1] += (a_out*pp - a_in*pp)/vol + rho*ut*ut*geo;   // s_mr
    out[2] -= rho*ur*ut*geo;                              // s_mt
    for (int k=0;k<NC;k++) rate[c*NC+k]=out[k];
}
__global__ void k_rate_z(const double* __restrict__ prim, int n_r, int n_z,
                         double r_min, double dr, double z_min, double dz,
                         double g, double* __restrict__ rate) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    long ncell = (long)n_r*n_z;
    if (c >= ncell) return;
    int i_r = (int)(c / n_z), i_z = (int)(c % n_z);
    if (i_r < NGH || i_r >= n_r-NGH || i_z < NGH || i_z >= n_z-NGH) return;
    double p[7*NP];
    for (int m=0;m<7;m++) {
        long cc = (long)i_r*n_z + (i_z-3+m);
        for (int k=0;k<NP;k++) p[m*NP+k]=prim[cc*NP+k];
    }
    double lfl[NP],lfr[NP],rfl[NP],rfr[NP];
    rz_pencil_faces(p, lfl, lfr, rfl, rfr);
    double fL[NC], fR[NC];
    rz_hllc(lfl, lfr, 3, g, fL);
    rz_hllc(rfl, rfr, 3, g, fR);
    double inv_dz = 1.0/dz;
    for (int k=0;k<NC;k++) rate[c*NC+k] += (fL[k] - fR[k]) * inv_dz;
}
// Convenience: the full class-A rate = k_rate_r then k_rate_z.
static inline void launch_rate(long cblk, int tpb, const double* prim, int n_r, int n_z,
                               double r_min, double dr, double z_min, double dz,
                               double g, double* rate) {
    k_rate_r<<<cblk,tpb>>>(prim,n_r,n_z,r_min,dr,z_min,dz,g,rate);
    k_rate_z<<<cblk,tpb>>>(prim,n_r,n_z,r_min,dr,z_min,dz,g,rate);
}

// ---- SDC composition kernels (compose_gas, flow-only) ------------------
// per scalar element (ncell*NC). Predictor: U = u0 + dt*A(u0).
__global__ void k_compose_pred(const double* __restrict__ u0, const double* __restrict__ r0,
                               double dt, long n, double* __restrict__ cons) {
    long i=(long)blockIdx.x*blockDim.x+threadIdx.x; if (i>=n) return;
    cons[i]=u0[i]+dt*r0[i];
}
// Correction: U = u0 + (dt/2)*A(u0) + (dt/2)*A(U_prev).
__global__ void k_compose_corr(const double* __restrict__ u0, const double* __restrict__ r0,
                               const double* __restrict__ rl, double half_dt, long n,
                               double* __restrict__ cons) {
    long i=(long)blockIdx.x*blockDim.x+threadIdx.x; if (i>=n) return;
    cons[i]=u0[i]+half_dt*r0[i]+half_dt*rl[i];
}

// ---- FFI: the RESIDENT class-A SDC march --------------------------------
// State lives on-device for the whole `nsteps` loop; the host uploads once at
// entry and downloads once at exit (cons is inout). 2-node Lobatto IMEX-SDC,
// explicit-only (D,R absent): 1 predictor + N_SDC_CORRECTIONS(=2) corrections
// per step, bit-for-formula identical to Sdc::step_flow's compose_gas. Fixed
// dt (the caller feeds CPU and GPU the same dt — stable_dt-on-device is S13b).
// Splitting a march into march(a)+march(b) inserts an FND-6 checkpoint at step
// a (the host round-trip IS the bit-faithful serialize/resume).
extern "C" int gpu_class_a_march(double* cons, int n_r, int n_z,
                                 double r_min, double dr, double z_min, double dz,
                                 double gamma, double dt, int nsteps) {
    long ncell=(long)n_r*n_z, nscal=ncell*NC;
    double *d_cons,*d_u0,*d_prim,*d_r0,*d_rl; int *d_bad;
    cudaMalloc(&d_cons, nscal*sizeof(double));
    cudaMalloc(&d_u0,   nscal*sizeof(double));
    cudaMalloc(&d_prim, ncell*NP*sizeof(double));
    cudaMalloc(&d_r0,   nscal*sizeof(double));
    cudaMalloc(&d_rl,   nscal*sizeof(double));
    cudaMalloc(&d_bad, sizeof(int));
    cudaMemcpy(d_cons, cons, nscal*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_bad, 0, sizeof(int));
    int tpb=256; long cblk=(ncell+tpb-1)/tpb, sblk=(nscal+tpb-1)/tpb;
    double half=0.5*dt;
    for (int step=0; step<nsteps; step++) {
        // node 0: u0 = U^n, rate0 = A(U^n)
        cudaMemcpy(d_u0, d_cons, nscal*sizeof(double), cudaMemcpyDeviceToDevice);
        k_fill_prims<<<cblk,tpb>>>(d_cons,(int)ncell,gamma,d_prim,d_bad);
        launch_rate(cblk,tpb,d_prim,n_r,n_z,r_min,dr,z_min,dz,gamma,d_r0);
        // predictor
        k_compose_pred<<<sblk,tpb>>>(d_u0,d_r0,dt,nscal,d_cons);
        // 2 corrections
        for (int corr=0; corr<2; corr++) {
            k_fill_prims<<<cblk,tpb>>>(d_cons,(int)ncell,gamma,d_prim,d_bad);
            launch_rate(cblk,tpb,d_prim,n_r,n_z,r_min,dr,z_min,dz,gamma,d_rl);
            k_compose_corr<<<sblk,tpb>>>(d_u0,d_r0,d_rl,half,nscal,d_cons);
        }
    }
    cudaDeviceSynchronize();
    int bad=0; cudaMemcpy(&bad,d_bad,sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(cons, d_cons, nscal*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_cons); cudaFree(d_u0); cudaFree(d_prim); cudaFree(d_r0); cudaFree(d_rl); cudaFree(d_bad);
    return bad;
}

// ---- FFI: single-shot class-A RHS (milestone-1 cross-check oracle) ------
// cons: dense [ncell*NC], row-major cell = i_r*n_z+i_z. Writes rate dense.
// Returns 0 on success, 1 if any cell was non-physical (halt, never clamp).
extern "C" int gpu_class_a_rhs(const double* cons, int n_r, int n_z,
                               double r_min, double dr, double z_min, double dz,
                               double gamma, double* rate) {
    long ncell = (long)n_r*n_z;
    double *d_cons,*d_prim,*d_rate; int *d_bad;
    cudaMalloc(&d_cons, ncell*NC*sizeof(double));
    cudaMalloc(&d_prim, ncell*NP*sizeof(double));
    cudaMalloc(&d_rate, ncell*NC*sizeof(double));
    cudaMalloc(&d_bad, sizeof(int));
    cudaMemcpy(d_cons, cons, ncell*NC*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_bad, 0, sizeof(int));
    int tpb=256;
    long blocks=(ncell+tpb-1)/tpb;
    k_fill_prims<<<blocks,tpb>>>(d_cons, (int)ncell, gamma, d_prim, d_bad);
    launch_rate(blocks, tpb, d_prim, n_r, n_z, r_min, dr, z_min, dz, gamma, d_rate);
    cudaDeviceSynchronize();
    int bad=0; cudaMemcpy(&bad, d_bad, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(rate, d_rate, ncell*NC*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_cons); cudaFree(d_prim); cudaFree(d_rate); cudaFree(d_bad);
    return bad;
}

// ---- FFI: throughput of the staged class-A rate (STEP 4 tuning) ----------
// Returns elapsed milliseconds for `iters` full class-A rate evaluations
// (k_rate_r + k_rate_z) on a resident prim buffer; the caller derives
// cell-RHS/s. fill_prims runs once (not timed); one warmup pass excluded.
extern "C" double gpu_class_a_bench(const double* cons, int n_r, int n_z,
                                    double r_min, double dr, double z_min, double dz,
                                    double gamma, int iters) {
    long ncell=(long)n_r*n_z, nscal=ncell*NC;
    double *d_cons,*d_prim,*d_rate; int *d_bad;
    cudaMalloc(&d_cons, nscal*sizeof(double));
    cudaMalloc(&d_prim, ncell*NP*sizeof(double));
    cudaMalloc(&d_rate, nscal*sizeof(double));
    cudaMalloc(&d_bad, sizeof(int));
    cudaMemcpy(d_cons, cons, nscal*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_bad, 0, sizeof(int));
    int tpb=256; long cblk=(ncell+tpb-1)/tpb;
    k_fill_prims<<<cblk,tpb>>>(d_cons,(int)ncell,gamma,d_prim,d_bad);
    launch_rate(cblk,tpb,d_prim,n_r,n_z,r_min,dr,z_min,dz,gamma,d_rate); // warmup
    cudaDeviceSynchronize();
    cudaEvent_t a,b; cudaEventCreate(&a); cudaEventCreate(&b);
    cudaEventRecord(a);
    for (int i=0;i<iters;i++)
        launch_rate(cblk,tpb,d_prim,n_r,n_z,r_min,dr,z_min,dz,gamma,d_rate);
    cudaEventRecord(b); cudaEventSynchronize(b);
    float ms=0; cudaEventElapsedTime(&ms,a,b);
    cudaEventDestroy(a); cudaEventDestroy(b);
    cudaFree(d_cons); cudaFree(d_prim); cudaFree(d_rate); cudaFree(d_bad);
    return (double)ms;
}
