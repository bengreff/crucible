// crucible-gpu — S13c residency: the class-A (explicit hyperbolic) step on the
// REAL 3-D CUT GEOMETRY — the ◆C3 world: uniform N_θ > 1 on a revolved cut
// contour (FND-3 §3.3 per-sector κ + six apertures), the r=0 axis (FND-2 §3.2
// θ↔θ+π parity-pair gather), the periodic θ-sweep, interior wall run
// boundaries (grid-aligned mirror / slip-reflect about the true wall normal),
// the SOLV-1 §3.3 geometric sources + the per-sector embedded-interface
// pressure closure (Grid::wall_closure_cell), and the Berger–Giuliani State
// Redistribution pass (Euler::srd) on the composed state. Bit-for-formula from
// crucible_solvers::euler (mod.rs sweep_r / sweep_theta / sweep_z /
// add_sources / srd / fill_ghosts_* / wall_ghosts_* / recon.rs / hllc.rs).
//
// LAYOUT: θ-plane-major dense cells, c = j·NRZ + i_r·n_z + i_z (NRZ = n_r·n_z),
// mirroring the brick's θ-plane-major storage. The (r,z) activity map and the
// pencil RUN tables (maximal active runs per line + the ghost KIND at each run
// end + the wall normal there) are host-precomputed geometry-time data (they
// are pure functions of the grid — Euler::scratch's `act` + the sweeps' run
// decomposition), uploaded once. Uniform N_θ only (the ◆C3 contour pins it
// structurally — mixed-N_θ reflux is off the ◆C3-parity path).
//
// SCOPE: interior + reflective walls + the two domain BCs the fixture needs
// (Reflecting / Transmissive); the injector mass-flow inflow + pressure
// outflow ghosts ride the composed-step session (Ben, session 28). GammaLaw
// EOS (the general-EOS HLLC aux-slot path = the composed-step session with
// the blend). The COUP-2 ledger reductions are not carried here (the S13b
// fixed-topology reduction primitive; the composed-step session).
//
// Determinism (META-1 §2.5): every kernel is a per-cell gather — a face's flux
// is computed identically by both its cells (same formula, same inputs, same
// ghost rule), so the 2×-recompute is bit-consistent and one writer per cell
// holds. SRD is two gathers over host-fixed member lists in the CPU's own
// order. Same-build reruns bit-identical; CPU↔GPU = FMA-order ECT.
#include <cuda_runtime.h>

#define NP 9
#define NC 7
#define NGH 3
#define KAPPA_SRD 0.5

__device__ __constant__ double G3_TAU = 6.283185307179586;

// ---- ghost-rule kinds at a run end (host-encoded) --------------------------
#define GK_DOMAIN_REFLECT 0
#define GK_DOMAIN_TRANSMISSIVE 1
#define GK_WALL_MIRROR 2
#define GK_WALL_SLIP 3
#define GK_AXIS 4

// ---- GammaLaw closures (euler/mod.rs GammaLaw) ------------------------------
__device__ __forceinline__ double g3_total_energy(const double* w, double g) {
    return w[4] / (g - 1.0) + 0.5 * w[0] * (w[1]*w[1] + w[2]*w[2] + w[3]*w[3]);
}
__device__ __forceinline__ double g3_sound(const double* w, double g) {
    return sqrt(g * w[4] / w[0]);
}
__device__ int g3_prim(const double* u, double g, double* w) {
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

// ---- HLLC-Batten (euler/hllc.rs), GammaLaw ----------------------------------
__device__ void g3_physical_flux(const double* w, int n, double g, double* f) {
    double rho=w[0], p=w[4], un=w[n], m=rho*un, e=g3_total_energy(w,g);
    f[0]=m; f[1]=m*w[1]; f[2]=m*w[2]; f[3]=m*w[3]; f[n]+=p;
    f[4]=un*(e+p); f[5]=m*w[5]; f[6]=m*w[6];
}
__device__ void g3_hllc(const double* wl, const double* wr, int n, double g, double* f) {
    double rho_l=wl[0], p_l=wl[4], rho_r=wr[0], p_r=wr[4];
    double un_l=wl[n], un_r=wr[n];
    double c_l=g3_sound(wl,g), c_r=g3_sound(wr,g);
    double sql=sqrt(rho_l), sqr=sqrt(rho_r), inv=1.0/(sql+sqr);
    double u1=(sql*wl[1]+sqr*wr[1])*inv;
    double u2=(sql*wl[2]+sqr*wr[2])*inv;
    double u3=(sql*wl[3]+sqr*wr[3])*inv;
    double urn=(n==1?u1:(n==2?u2:u3));
    double h_l=(g3_total_energy(wl,g)+p_l)/rho_l;
    double h_r=(g3_total_energy(wr,g)+p_r)/rho_r;
    double h_roe=(sql*h_l+sqr*h_r)*inv;
    double q2=u1*u1+u2*u2+u3*u3;
    double arg=(g-1.0)*(h_roe-0.5*q2);
    double c_roe=sqrt(arg>0.0?arg:0.0);
    double s_l=fmin(un_l-c_l, urn-c_roe);
    double s_r=fmax(un_r+c_r, urn+c_roe);
    double ml=rho_l*(s_l-un_l), mr=rho_r*(s_r-un_r);
    double s_m=(mr*un_r-ml*un_l+p_l-p_r)/(mr-ml);
    for (int k=0;k<NC;k++) f[k]=0.0;
    if (s_l>=0.0) { g3_physical_flux(wl,n,g,f); return; }
    if (s_r<=0.0) { g3_physical_flux(wr,n,g,f); return; }
    const double* w=(s_m>=0.0)?wl:wr; double s_k=(s_m>=0.0)?s_l:s_r;
    double rho=w[0], p=w[4], un=w[n], e=g3_total_energy(w,g);
    double p_star=rho*(un-s_k)*(un-s_m)+p;
    double u_k[NC]={rho, rho*w[1], rho*w[2], rho*w[3], e, rho*w[5], rho*w[6]};
    double fac=(s_k-un)/(s_k-s_m), rho_s=rho*fac;
    double u_s[NC];
    u_s[0]=rho_s; u_s[1]=rho_s*w[1]; u_s[2]=rho_s*w[2]; u_s[3]=rho_s*w[3];
    u_s[n]=rho_s*s_m;
    u_s[4]=fac*e+(p_star*s_m-p*un)/(s_k-s_m);
    u_s[5]=rho_s*w[5]; u_s[6]=rho_s*w[6];
    double fk[NC]; for(int k=0;k<NC;k++) fk[k]=0.0; g3_physical_flux(w,n,g,fk);
    for (int k=0;k<NC;k++) f[k]=fk[k]+s_k*(u_s[k]-u_k[k]);
}

// ---- PPM (euler/recon.rs) on a 7-cell pencil ---------------------------------
__device__ __forceinline__ double g3_mc_slope(double wm, double w0, double wp) {
    double dl=w0-wm, dr=wp-w0;
    if (dl*dr <= 0.0) return 0.0;
    double dc=0.5*(dl+dr);
    double a=fabs(dc); a=fmin(a, 2.0*fabs(dl)); a=fmin(a, 2.0*fabs(dr));
    return copysign(a, dc);
}
__device__ void g3_slope(const double* wm, const double* w0, const double* wp, double* s) {
    for (int k=0;k<NP;k++) s[k]=g3_mc_slope(wm[k], w0[k], wp[k]);
}
__device__ void g3_iface(const double* wa, const double* wb, const double* sa, const double* sb, double* out) {
    for (int k=0;k<NP;k++) out[k]=0.5*(wa[k]+wb[k]) - (sb[k]-sa[k])/6.0;
}
__device__ void g3_edge(const double* lo_if, const double* hi_if, const double* c, double* elo, double* ehi) {
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
// Faces bounding the centre cell p[3] of the 7-cell pencil p[0..6].
__device__ void g3_pencil_faces(const double* p, double* lfl, double* lfr, double* rfl, double* rfr) {
    double s1[NP],s2[NP],s3[NP],s4[NP],s5[NP];
    g3_slope(p+0*NP, p+1*NP, p+2*NP, s1);
    g3_slope(p+1*NP, p+2*NP, p+3*NP, s2);
    g3_slope(p+2*NP, p+3*NP, p+4*NP, s3);
    g3_slope(p+3*NP, p+4*NP, p+5*NP, s4);
    g3_slope(p+4*NP, p+5*NP, p+6*NP, s5);
    double if1[NP],if2[NP],if3[NP],if4[NP];
    g3_iface(p+1*NP, p+2*NP, s1, s2, if1);
    g3_iface(p+2*NP, p+3*NP, s2, s3, if2);
    g3_iface(p+3*NP, p+4*NP, s3, s4, if3);
    g3_iface(p+4*NP, p+5*NP, s4, s5, if4);
    double e2lo[NP],e2hi[NP],e3lo[NP],e3hi[NP],e4lo[NP],e4hi[NP];
    g3_edge(if1, if2, p+2*NP, e2lo, e2hi);
    g3_edge(if2, if3, p+3*NP, e3lo, e3hi);
    g3_edge(if3, if4, p+4*NP, e4lo, e4hi);
    for (int k=0;k<NP;k++) { lfl[k]=e2hi[k]; lfr[k]=e3lo[k]; rfl[k]=e3hi[k]; rfr[k]=e4lo[k]; }
}

// ---- the world: geometry + run tables (host-precomputed) --------------------
struct World {
    int n_r, n_z, nt;            // NRZ = n_r*n_z; N = nt*NRZ
    double r_min, dr, z_min, dz, gamma;
    const int* act;              // [NRZ] (r,z) activity
    const double* kappa;         // [N] per-sector κ
    const double* ap;            // [6*N] apertures: d*N + c, d ∈ {r−,r+,z−,z+,θ−,θ+}
    int has_geom;                // wall-closure sources armed (cut worlds)
    // r-direction runs (per rz cell): run start/len along i_r, ghost kinds +
    // wall normals (n_r, n_z) at the low/high run ends.
    const int* rs_r; const int* rl_r; const int* klo_r; const int* khi_r;
    const double* nlo_r; const double* nhi_r;     // [2*NRZ] each: (n_r, n_z)
    // z-direction runs (per rz cell): along i_z.
    const int* rs_z; const int* rl_z; const int* klo_z; const int* khi_z;
    const double* nlo_z; const double* nhi_z;
};

// Metric (grid/lib.rs — the same expression order as face_area_r /
// cell_volume / face_area_z / face_area_theta).
__device__ __forceinline__ double g3_face_r(const World& W, int f) { return W.r_min + (double)f * W.dr; }
__device__ __forceinline__ double g3_area_r(const World& W, int f) {
    return g3_face_r(W, f) * (G3_TAU / (double)W.nt) * W.dz;
}
__device__ __forceinline__ double g3_vol(const World& W, int i_r) {
    double ri = g3_face_r(W, i_r), ro = g3_face_r(W, i_r + 1);
    return 0.5 * (ro*ro - ri*ri) * (G3_TAU / (double)W.nt) * W.dz;
}
__device__ __forceinline__ double g3_area_z(const World& W, int i_r) {
    double ri = g3_face_r(W, i_r), ro = g3_face_r(W, i_r + 1);
    return 0.5 * (ro*ro - ri*ri) * (G3_TAU / (double)W.nt);
}
__device__ __forceinline__ double g3_area_th(const World& W) { return W.dr * W.dz; }

// slip_reflect: v ← v − 2(v·n̂)n̂ on (u_r, u_z).
__device__ __forceinline__ void g3_slip(double* m, double nr, double nz) {
    double vn = m[1]*nr + m[3]*nz;
    m[1] -= 2.0*vn*nr; m[3] -= 2.0*vn*nz;
}

// Gather one pencil entry of a MERIDIONAL sweep (dir: 1 = r, 3 = z) for the
// cell at (j, i_r, i_z) at pencil offset `m` ∈ [−3, 3] — the cell itself, an
// in-run neighbour, or a ghost by the run end's rule (fill_ghosts_low/high,
// wall_ghosts_low/high, the axis parity pair). `normal` = the velocity slot.
__device__ void g3_gather_rz(const World& W, const double* __restrict__ prim,
                             int j, int i_r, int i_z, int dir, int m, double* out) {
    long NRZ = (long)W.n_r * W.n_z;
    int rz = i_r * W.n_z + i_z;
    int start, len, klo, khi; const double *nlo, *nhi; int normal;
    int pos;   // this cell's coordinate along the pencil
    if (dir == 1) { start = W.rs_r[rz]; len = W.rl_r[rz]; klo = W.klo_r[rz]; khi = W.khi_r[rz];
                    nlo = W.nlo_r + 2*rz; nhi = W.nhi_r + 2*rz; normal = 1; pos = i_r; }
    else          { start = W.rs_z[rz]; len = W.rl_z[rz]; klo = W.klo_z[rz]; khi = W.khi_z[rz];
                    nlo = W.nlo_z + 2*rz; nhi = W.nhi_z + 2*rz; normal = 3; pos = i_z; }
    int q = pos - start + m;   // run-relative coordinate of the wanted entry
    // (source cell coordinate along the line, source θ-plane, flip mask)
    int src; int jsrc = j; int kind = -1; const double* nh = nlo;
    if (q >= 0 && q < len) { src = start + q; }
    else if (q < 0) {
        int k = -q;                                   // ghost k = 1..3 beyond the low end
        kind = klo; nh = nlo;
        src = (kind == GK_DOMAIN_TRANSMISSIVE) ? start : start + min(k - 1, len - 1);
        if (kind == GK_AXIS) jsrc = (j + W.nt / 2) % W.nt;
    } else {
        int k = q - len + 1;                          // ghost k = 1..3 beyond the high end
        kind = khi; nh = nhi;
        src = (kind == GK_DOMAIN_TRANSMISSIVE) ? start + len - 1 : start + max(len - k, 0);
    }
    long cs = (dir == 1) ? ((long)jsrc * NRZ + (long)src * W.n_z + i_z)
                         : ((long)jsrc * NRZ + (long)i_r * W.n_z + src);
    for (int k = 0; k < NP; k++) out[k] = prim[cs*NP + k];
    if (kind == GK_DOMAIN_REFLECT || kind == GK_WALL_MIRROR) out[normal] = -out[normal];
    else if (kind == GK_WALL_SLIP) g3_slip(out, nh[0], nh[1]);
    else if (kind == GK_AXIS) { out[1] = -out[1]; out[2] = -out[2]; }
}

// ---- kernels ------------------------------------------------------------------
__global__ void k3_fill_prims(const double* __restrict__ cons, const World W, long N,
                              double* __restrict__ prim, int* __restrict__ bad) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int rz = (int)(c % NRZ);
    if (!W.act[rz]) { for (int k=0;k<NP;k++) prim[c*NP+k]=0.0; return; }
    double w[NP];
    if (!g3_prim(cons + c*NC, W.gamma, w)) { atomicExch(bad, 1); return; }
    for (int k=0;k<NP;k++) prim[c*NP+k]=w[k];
}

// r-sweep (WRITES the rate): the single difference (af_lo − af_hi)/(κV) with
// af = (A·ap)·F; the face aperture is the CPU's canonical read — the r− of the
// face's high-side cell for an in-run face, the r+ of the last cell at the run's
// high end.
__global__ void k3_rate_r(const double* __restrict__ prim, const World W, long N,
                          double* __restrict__ rate) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z, i_z = rz % W.n_z;
    if (!W.act[rz]) { for (int k=0;k<NC;k++) rate[c*NC+k]=0.0; return; }
    double p[7*NP];
    for (int m=-3;m<=3;m++) g3_gather_rz(W, prim, j, i_r, i_z, 1, m, p + (m+3)*NP);
    double lfl[NP],lfr[NP],rfl[NP],rfr[NP];
    g3_pencil_faces(p, lfl, lfr, rfl, rfr);
    double fL[NC], fR[NC];
    g3_hllc(lfl, lfr, 1, W.gamma, fL);
    g3_hllc(rfl, rfr, 1, W.gamma, fR);
    int start = W.rs_r[rz], len = W.rl_r[rz];
    double ap_lo = W.ap[0*N + c];
    double ap_hi = (i_r + 1 < start + len) ? W.ap[0*N + (c + W.n_z)] : W.ap[1*N + c];
    double a_lo = g3_area_r(W, i_r) * ap_lo;
    double a_hi = g3_area_r(W, i_r + 1) * ap_hi;
    double kv = W.kappa[c] * g3_vol(W, i_r);
    for (int k=0;k<NC;k++) rate[c*NC+k] = (a_lo*fL[k] - a_hi*fR[k]) / kv;
}

// θ-sweep (ADDS): periodic ring pencil; face fi reads the θ+ aperture of ring
// cell (fi − 1) mod n; rate += (af_lo − af_hi)·(A_θ/V)/κ.
__global__ void k3_rate_theta(const double* __restrict__ prim, const World W, long N,
                              double* __restrict__ rate) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z;
    if (!W.act[rz]) return;
    int n = W.nt;
    double p[7*NP];
    for (int m=-3;m<=3;m++) {
        int jj = ((j + m) % n + n) % n;
        long cs = (long)jj * NRZ + rz;
        for (int k=0;k<NP;k++) p[(m+3)*NP+k] = prim[cs*NP+k];
    }
    double lfl[NP],lfr[NP],rfl[NP],rfr[NP];
    g3_pencil_faces(p, lfl, lfr, rfl, rfr);
    double fL[NC], fR[NC];
    g3_hllc(lfl, lfr, 2, W.gamma, fL);
    g3_hllc(rfl, rfr, 2, W.gamma, fR);
    int jm = (j + n - 1) % n;
    double ap_lo = W.ap[5*N + ((long)jm * NRZ + rz)];   // θ+ of the cell before face j
    double ap_hi = W.ap[5*N + c];                        // θ+ of this cell (face j+1)
    double inv = g3_area_th(W) / g3_vol(W, i_r);
    double kap = W.kappa[c];
    for (int k=0;k<NC;k++) rate[c*NC+k] += ((fL[k]*ap_lo) - (fR[k]*ap_hi)) * inv / kap;
}

// z-sweep (ADDS): metric-ratio form (F_lo − F_hi)/dz/κ with F = flux·ap.
__global__ void k3_rate_z(const double* __restrict__ prim, const World W, long N,
                          double* __restrict__ rate) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int j = (int)(c / NRZ), rz = (int)(c % NRZ);
    int i_r = rz / W.n_z, i_z = rz % W.n_z;
    if (!W.act[rz]) return;
    double p[7*NP];
    for (int m=-3;m<=3;m++) g3_gather_rz(W, prim, j, i_r, i_z, 3, m, p + (m+3)*NP);
    double lfl[NP],lfr[NP],rfl[NP],rfr[NP];
    g3_pencil_faces(p, lfl, lfr, rfl, rfr);
    double fL[NC], fR[NC];
    g3_hllc(lfl, lfr, 3, W.gamma, fL);
    g3_hllc(rfl, rfr, 3, W.gamma, fR);
    int start = W.rs_z[rz], len = W.rl_z[rz];
    double ap_lo = W.ap[2*N + c];
    double ap_hi = (i_z + 1 < start + len) ? W.ap[2*N + (c + 1)] : W.ap[3*N + c];
    double inv_dz = 1.0 / W.dz;
    double kap = W.kappa[c];
    for (int k=0;k<NC;k++) rate[c*NC+k] += ((fL[k]*ap_lo) - (fR[k]*ap_hi)) * inv_dz / kap;
}

// Sources (ADDS): the SOLV-1 §3.3 geometric terms + (cut worlds) the per-sector
// wall-closure pressure force p·W/(κV) — add_sources with no external intake.
__global__ void k3_sources(const double* __restrict__ prim, const World W, long N,
                           double* __restrict__ rate) {
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    if (c >= N) return;
    long NRZ = (long)W.n_r * W.n_z;
    int rz = (int)(c % NRZ);
    int i_r = rz / W.n_z;
    if (!W.act[rz]) return;
    double a_in = g3_area_r(W, i_r), a_out = g3_area_r(W, i_r + 1);
    double vol = g3_vol(W, i_r);
    double geo = (a_out - a_in) / vol;
    const double* wc = prim + c*NP;
    double rho=wc[0], ur=wc[1], ut=wc[2], pp=wc[4];
    double s_mr = (a_out*pp - a_in*pp)/vol + rho*ut*ut*geo;
    double s_mt = rho*ur*ut*geo;
    rate[c*NC+1] += s_mr;
    rate[c*NC+2] -= s_mt;
    if (W.has_geom) {
        double kappa = W.kappa[c];
        double w_r = W.ap[1*N+c]*a_out - W.ap[0*N+c]*a_in - kappa*(a_out - a_in);
        double w_th = (W.ap[5*N+c] - W.ap[4*N+c]) * g3_area_th(W);
        double w_z = (W.ap[3*N+c] - W.ap[2*N+c]) * g3_area_z(W, i_r);
        double inv_kv = 1.0 / (kappa * vol);
        double wr_kv = w_r * inv_kv, wz_kv = w_z * inv_kv;
        rate[c*NC+1] += pp * wr_kv;
        rate[c*NC+3] += pp * wz_kv;
        if (w_th != 0.0) { double wt_kv = w_th * inv_kv; rate[c*NC+2] += pp * wt_kv; }
    }
}

static inline void launch_rate3(long cblk, int tpb, const double* prim, const World& W, long N, double* rate) {
    k3_rate_r<<<cblk,tpb>>>(prim, W, N, rate);
    k3_rate_theta<<<cblk,tpb>>>(prim, W, N, rate);
    k3_rate_z<<<cblk,tpb>>>(prim, W, N, rate);
    k3_sources<<<cblk,tpb>>>(prim, W, N, rate);
}

// ---- State Redistribution (Euler::srd) ---------------------------------------
// Host-fixed structure in the CPU's own order: `ns` small neighbourhoods
// (CSR: mem_off[s]..mem_off[s+1] member cells `mem[]` with weights `mem_kv[]`,
// members[0] = owner), overlap counts `cnt[c]`, and the INVERSE map for the
// `na` affected cells: aff[a] = cell, aff_owner[a], and the neighbourhoods
// containing it (CSR: inv_off[a]..inv_off[a+1] → inv[]), in small-list order.
struct SrdTab {
    int ns, na;
    const int* mem_off; const int* mem; const double* mem_kv;
    const int* cnt;
    const int* aff; const int* aff_owner; const int* inv_off; const int* inv;
};
// Q_s = Σ_m (kv_m/n_m)·U_m / Σ_m (kv_m/n_m), member order.
__global__ void k3_srd_q(const double* __restrict__ cons, const SrdTab T, double* __restrict__ q) {
    int s = blockIdx.x*blockDim.x + threadIdx.x;
    if (s >= T.ns) return;
    double num[NC]; for (int k=0;k<NC;k++) num[k]=0.0;
    double den = 0.0;
    for (int m = T.mem_off[s]; m < T.mem_off[s+1]; m++) {
        int c = T.mem[m];
        double w = T.mem_kv[m] / (double)T.cnt[c];
        for (int k=0;k<NC;k++) num[k] += w * cons[(long)c*NC+k];
        den += w;
    }
    for (int k=0;k<NC;k++) q[(long)s*NC+k] = num[k] / den;
}
// U_c = ((owner ? 0 : U_c) + Σ_{s∋c} Q_s) / n_c — one writer per affected cell.
__global__ void k3_srd_apply(const SrdTab T, const double* __restrict__ q, double* __restrict__ cons) {
    int a = blockIdx.x*blockDim.x + threadIdx.x;
    if (a >= T.na) return;
    int c = T.aff[a];
    double acc[NC];
    for (int k=0;k<NC;k++) acc[k] = T.aff_owner[a] ? 0.0 : cons[(long)c*NC+k];
    for (int i = T.inv_off[a]; i < T.inv_off[a+1]; i++) {
        int s = T.inv[i];
        for (int k=0;k<NC;k++) acc[k] += q[(long)s*NC+k];
    }
    double n = (double)T.cnt[c];
    for (int k=0;k<NC;k++) cons[(long)c*NC+k] = acc[k] / n;
}
static inline void launch_srd(const SrdTab& T, double* d_q, double* d_cons) {
    if (T.ns == 0) return;
    int tpb = 128;
    k3_srd_q<<<(T.ns + tpb - 1)/tpb, tpb>>>(d_cons, T, d_q);
    k3_srd_apply<<<(T.na + tpb - 1)/tpb, tpb>>>(T, d_q, d_cons);
}

// ---- SDC composition (compose_gas, flow-only) --------------------------------
__global__ void k3_compose_pred(const double* __restrict__ u0, const double* __restrict__ r0,
                                double dt, long n, double* __restrict__ cons) {
    long i=(long)blockIdx.x*blockDim.x+threadIdx.x; if (i>=n) return;
    cons[i]=u0[i]+dt*r0[i];
}
__global__ void k3_compose_corr(const double* __restrict__ u0, const double* __restrict__ r0,
                                const double* __restrict__ rl, double half_dt, long n,
                                double* __restrict__ cons) {
    long i=(long)blockIdx.x*blockDim.x+threadIdx.x; if (i>=n) return;
    cons[i]=u0[i]+half_dt*r0[i]+half_dt*rl[i];
}

// ---- host-side upload helpers -------------------------------------------------
template <typename T> static T* up(const T* h, long n, void** keep, int* nk) {
    T* d; cudaMalloc(&d, n * sizeof(T));
    cudaMemcpy(d, h, n * sizeof(T), cudaMemcpyHostToDevice);
    keep[(*nk)++] = d; return d;
}
struct HostWorld {
    int n_r, n_z, nt; double r_min, dr, z_min, dz, gamma; int has_geom;
    const int* act; const double* kappa; const double* ap;
    const int* rs_r; const int* rl_r; const int* klo_r; const int* khi_r; const double* nlo_r; const double* nhi_r;
    const int* rs_z; const int* rl_z; const int* klo_z; const int* khi_z; const double* nlo_z; const double* nhi_z;
};
struct HostSrd {
    int ns, na; const int* mem_off; const int* mem; const double* mem_kv; const int* cnt;
    const int* aff; const int* aff_owner; const int* inv_off; const int* inv;
};
static World upload_world(const HostWorld& h, void** keep, int* nk) {
    long NRZ = (long)h.n_r * h.n_z, N = NRZ * h.nt;
    World W;
    W.n_r=h.n_r; W.n_z=h.n_z; W.nt=h.nt; W.r_min=h.r_min; W.dr=h.dr; W.z_min=h.z_min; W.dz=h.dz;
    W.gamma=h.gamma; W.has_geom=h.has_geom;
    W.act = up(h.act, NRZ, keep, nk);
    W.kappa = up(h.kappa, N, keep, nk);
    W.ap = up(h.ap, 6*N, keep, nk);
    W.rs_r = up(h.rs_r, NRZ, keep, nk); W.rl_r = up(h.rl_r, NRZ, keep, nk);
    W.klo_r = up(h.klo_r, NRZ, keep, nk); W.khi_r = up(h.khi_r, NRZ, keep, nk);
    W.nlo_r = up(h.nlo_r, 2*NRZ, keep, nk); W.nhi_r = up(h.nhi_r, 2*NRZ, keep, nk);
    W.rs_z = up(h.rs_z, NRZ, keep, nk); W.rl_z = up(h.rl_z, NRZ, keep, nk);
    W.klo_z = up(h.klo_z, NRZ, keep, nk); W.khi_z = up(h.khi_z, NRZ, keep, nk);
    W.nlo_z = up(h.nlo_z, 2*NRZ, keep, nk); W.nhi_z = up(h.nhi_z, 2*NRZ, keep, nk);
    return W;
}
static SrdTab upload_srd(const HostSrd& h, long N, void** keep, int* nk) {
    SrdTab T; T.ns = h.ns; T.na = h.na;
    if (h.ns == 0) { T.mem_off=0; T.mem=0; T.mem_kv=0; T.cnt=0; T.aff=0; T.aff_owner=0; T.inv_off=0; T.inv=0; return T; }
    long nmem = h.mem_off[h.ns], ninv = h.inv_off[h.na];
    T.mem_off = up(h.mem_off, h.ns + 1, keep, nk);
    T.mem = up(h.mem, nmem, keep, nk);
    T.mem_kv = up(h.mem_kv, nmem, keep, nk);
    T.cnt = up(h.cnt, N, keep, nk);
    T.aff = up(h.aff, h.na, keep, nk);
    T.aff_owner = up(h.aff_owner, h.na, keep, nk);
    T.inv_off = up(h.inv_off, h.na + 1, keep, nk);
    T.inv = up(h.inv, ninv, keep, nk);
    return T;
}
static void free_keep(void** keep, int nk) { for (int i = 0; i < nk; i++) cudaFree(keep[i]); }

// ---- FFI: single-shot class-A RHS on the 3-D cut world ------------------------
extern "C" int gpu_class_a_rhs_3d(const double* cons, const HostWorld* hw, double* rate) {
    void* keep[32]; int nk = 0;
    World W = upload_world(*hw, keep, &nk);
    long N = (long)hw->n_r * hw->n_z * hw->nt, nscal = N * NC;
    double *d_cons, *d_prim, *d_rate; int* d_bad;
    cudaMalloc(&d_cons, nscal*sizeof(double)); cudaMalloc(&d_prim, N*NP*sizeof(double));
    cudaMalloc(&d_rate, nscal*sizeof(double)); cudaMalloc(&d_bad, sizeof(int));
    cudaMemcpy(d_cons, cons, nscal*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_bad, 0, sizeof(int));
    int tpb = 128; long cblk = (N + tpb - 1) / tpb;
    k3_fill_prims<<<cblk,tpb>>>(d_cons, W, N, d_prim, d_bad);
    launch_rate3(cblk, tpb, d_prim, W, N, d_rate);
    cudaDeviceSynchronize();
    int bad = 0; cudaMemcpy(&bad, d_bad, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(rate, d_rate, nscal*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_cons); cudaFree(d_prim); cudaFree(d_rate); cudaFree(d_bad);
    free_keep(keep, nk);
    return bad;
}

// ---- FFI: one SRD pass on a state (cons inout) -----------------------------------
extern "C" void gpu_srd_3d(double* cons, const HostWorld* hw, const HostSrd* hs) {
    void* keep[32]; int nk = 0;
    long N = (long)hw->n_r * hw->n_z * hw->nt, nscal = N * NC;
    SrdTab T = upload_srd(*hs, N, keep, &nk);
    double *d_cons, *d_q;
    cudaMalloc(&d_cons, nscal*sizeof(double)); cudaMalloc(&d_q, (long)(hs->ns > 0 ? hs->ns : 1)*NC*sizeof(double));
    cudaMemcpy(d_cons, cons, nscal*sizeof(double), cudaMemcpyHostToDevice);
    launch_srd(T, d_q, d_cons);
    cudaDeviceSynchronize();
    cudaMemcpy(cons, d_cons, nscal*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_cons); cudaFree(d_q);
    free_keep(keep, nk);
}

// ---- FFI: the RESIDENT class-A SDC march on the 3-D cut world -------------------
// Per step: u0 = U^n; A(u0) → predictor compose → SRD; then N_SDC_CORRECTIONS
// = 2 × { A(U) at the SRD'd iterate → correction compose → SRD } — Sdc::step_flow
// with the stagewise SRD (apply_srd after every composition). cons is inout;
// state resident for the whole loop.
extern "C" int gpu_class_a_march_3d(double* cons, const HostWorld* hw, const HostSrd* hs,
                                    double dt, int nsteps) {
    void* keep[48]; int nk = 0;
    World W = upload_world(*hw, keep, &nk);
    long N = (long)hw->n_r * hw->n_z * hw->nt, nscal = N * NC;
    SrdTab T = upload_srd(*hs, N, keep, &nk);
    double *d_cons, *d_u0, *d_prim, *d_r0, *d_rl, *d_q; int* d_bad;
    cudaMalloc(&d_cons, nscal*sizeof(double)); cudaMalloc(&d_u0, nscal*sizeof(double));
    cudaMalloc(&d_prim, N*NP*sizeof(double));
    cudaMalloc(&d_r0, nscal*sizeof(double)); cudaMalloc(&d_rl, nscal*sizeof(double));
    cudaMalloc(&d_q, (long)(hs->ns > 0 ? hs->ns : 1)*NC*sizeof(double));
    cudaMalloc(&d_bad, sizeof(int));
    cudaMemcpy(d_cons, cons, nscal*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_bad, 0, sizeof(int));
    int tpb = 128; long cblk = (N + tpb - 1)/tpb, sblk = (nscal + tpb - 1)/tpb;
    double half = 0.5 * dt;
    for (int step = 0; step < nsteps; step++) {
        cudaMemcpy(d_u0, d_cons, nscal*sizeof(double), cudaMemcpyDeviceToDevice);
        k3_fill_prims<<<cblk,tpb>>>(d_cons, W, N, d_prim, d_bad);
        launch_rate3(cblk, tpb, d_prim, W, N, d_r0);
        k3_compose_pred<<<sblk,tpb>>>(d_u0, d_r0, dt, nscal, d_cons);
        launch_srd(T, d_q, d_cons);
        for (int corr = 0; corr < 2; corr++) {
            k3_fill_prims<<<cblk,tpb>>>(d_cons, W, N, d_prim, d_bad);
            launch_rate3(cblk, tpb, d_prim, W, N, d_rl);
            k3_compose_corr<<<sblk,tpb>>>(d_u0, d_r0, d_rl, half, nscal, d_cons);
            launch_srd(T, d_q, d_cons);
        }
    }
    cudaDeviceSynchronize();
    int bad = 0; cudaMemcpy(&bad, d_bad, sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(cons, d_cons, nscal*sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_cons); cudaFree(d_u0); cudaFree(d_prim); cudaFree(d_r0); cudaFree(d_rl); cudaFree(d_q); cudaFree(d_bad);
    free_keep(keep, nk);
    return bad;
}

// ---- stable_dt on the 3-D world (Euler::stable_dt_inner, no combustion) -------
// σ = (|u_r|+c)/dr + (|u_z|+c)/dz + [N_θ>1] (|u_θ|+c)/(r̄·Δθ) over active cells;
// Δt = cfl / max σ. MAX is exactly order-independent, so the tree shape is
// immaterial to determinism; CPU↔GPU differ only at the per-cell σ's FMA order.
#define SDTPB 256
__global__ void k3_sig_partial(const double* __restrict__ cons, const World W, long N,
                               double* __restrict__ partial, int* __restrict__ bad) {
    __shared__ double sh[SDTPB];
    long c = (long)blockIdx.x*blockDim.x + threadIdx.x;
    double v = 0.0;
    if (c < N) {
        long NRZ = (long)W.n_r * W.n_z;
        int rz = (int)(c % NRZ);
        if (W.act[rz]) {
            double w[NP];
            if (g3_prim(cons + c*NC, W.gamma, w)) {
                double cs = g3_sound(w, W.gamma);
                v = (fabs(w[1]) + cs)/W.dr + (fabs(w[3]) + cs)/W.dz;
                if (W.nt > 1) {
                    int i_r = rz / W.n_z;
                    double r_c = W.r_min + ((double)i_r + 0.5) * W.dr;
                    double dth = G3_TAU / (double)W.nt;
                    v += (fabs(w[2]) + cs) / (r_c * dth);
                }
            } else {
                atomicExch(bad, 1);   // halt flag (not a physics-path atomic)
            }
        }
    }
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int s = SDTPB/2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sh[threadIdx.x] = fmax(sh[threadIdx.x], sh[threadIdx.x+s]);
        __syncthreads();
    }
    if (threadIdx.x == 0) partial[blockIdx.x] = sh[0];
}
__global__ void k3_max_final(const double* __restrict__ partial, int nblk, double* __restrict__ out) {
    __shared__ double sh[SDTPB];
    double v = 0.0;
    for (int j = threadIdx.x; j < nblk; j += SDTPB) v = fmax(v, partial[j]);
    sh[threadIdx.x] = v;
    __syncthreads();
    for (int s = SDTPB/2; s > 0; s >>= 1) {
        if (threadIdx.x < s) sh[threadIdx.x] = fmax(sh[threadIdx.x], sh[threadIdx.x+s]);
        __syncthreads();
    }
    if (threadIdx.x == 0) *out = sh[0];
}
extern "C" double gpu_stable_dt_3d(const double* cons, const HostWorld* hw, double cfl, int* bad) {
    void* keep[32]; int nk = 0;
    World W = upload_world(*hw, keep, &nk);
    long N = (long)hw->n_r * hw->n_z * hw->nt, nscal = N * NC;
    int nblk = (int)((N + SDTPB - 1)/SDTPB);
    double *d_cons, *d_partial, *d_scalar; int* d_bad;
    cudaMalloc(&d_cons, nscal*sizeof(double)); cudaMalloc(&d_partial, nblk*sizeof(double));
    cudaMalloc(&d_scalar, sizeof(double)); cudaMalloc(&d_bad, sizeof(int));
    cudaMemcpy(d_cons, cons, nscal*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_bad, 0, sizeof(int));
    k3_sig_partial<<<nblk,SDTPB>>>(d_cons, W, N, d_partial, d_bad);
    k3_max_final<<<1,SDTPB>>>(d_partial, nblk, d_scalar);
    cudaDeviceSynchronize();
    double max_sig = 0.0; cudaMemcpy(&max_sig, d_scalar, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(bad, d_bad, sizeof(int), cudaMemcpyDeviceToHost);
    cudaFree(d_cons); cudaFree(d_partial); cudaFree(d_scalar); cudaFree(d_bad);
    free_keep(keep, nk);
    if (*bad || !(max_sig > 0.0)) return -1.0;
    return cfl / max_sig;
}

// ---- FFI: throughput of the 3-D rate on the resident world ----------------------
extern "C" double gpu_class_a_bench_3d(const double* cons, const HostWorld* hw, int iters) {
    void* keep[32]; int nk = 0;
    World W = upload_world(*hw, keep, &nk);
    long N = (long)hw->n_r * hw->n_z * hw->nt, nscal = N * NC;
    double *d_cons, *d_prim, *d_rate; int* d_bad;
    cudaMalloc(&d_cons, nscal*sizeof(double)); cudaMalloc(&d_prim, N*NP*sizeof(double));
    cudaMalloc(&d_rate, nscal*sizeof(double)); cudaMalloc(&d_bad, sizeof(int));
    cudaMemcpy(d_cons, cons, nscal*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemset(d_bad, 0, sizeof(int));
    int tpb = 128; long cblk = (N + tpb - 1)/tpb;
    k3_fill_prims<<<cblk,tpb>>>(d_cons, W, N, d_prim, d_bad);
    launch_rate3(cblk, tpb, d_prim, W, N, d_rate);
    cudaDeviceSynchronize();
    cudaEvent_t a, b; cudaEventCreate(&a); cudaEventCreate(&b);
    cudaEventRecord(a);
    for (int i = 0; i < iters; i++) launch_rate3(cblk, tpb, d_prim, W, N, d_rate);
    cudaEventRecord(b); cudaEventSynchronize(b);
    float ms = 0; cudaEventElapsedTime(&ms, a, b);
    cudaEventDestroy(a); cudaEventDestroy(b);
    cudaFree(d_cons); cudaFree(d_prim); cudaFree(d_rate); cudaFree(d_bad);
    free_keep(keep, nk);
    return (double)ms;
}
