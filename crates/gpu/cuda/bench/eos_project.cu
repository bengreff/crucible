// S12 GPU spike — hot kernel (b): the per-cell TableEos (p,h,Z) projection
// (SOLV-1 §3.4 / FND-5), ported from euler/table_eos.rs illinois_root +
// project_pressure. Each thread roots  rho_tab(p, e_q + p/rho, Z) - rho = 0
// over the p-bracket by Illinois regula-falsi (fixed 48-iter cap) with a
// 64-iter bisection backstop — same fixed-order comparison ladder as the CPU,
// so per-thread bit-deterministic, gather-only, no atomics in the solve.
// The surface here is a SYNTHETIC trilinear (p,h,Z) grid (rho = A(Z)*p/h) that
// exercises the identical control flow + 8-point trilinear gather; binding the
// real HDF5 equilibrium surface is S13 residency. Correctness is self-checked:
// each cell's target is built at a known p_true, and the solver must recover it.
#include <cstdio>
#include <cstring>
#include <cuda_runtime.h>

#define NP 128
#define NH 128
#define NZ 8
#define N_P_ITER_MAX 48
#define N_P_BISECT   64
#define EPS_P_PROJECTION 1e-11
#define PMIN 1.0e3
#define PMAX 1.0e7
#define HMIN 1.0e5
#define HMAX 1.0e7

#define CK(x) do { cudaError_t e=(x); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %d\n",cudaGetErrorString(e),__LINE__); return 1;} } while(0)

// analytic node density (also used to fill the grid host-side). MULTILINEAR in
// (p,h,Z) so trilinear interp reproduces it exactly ⇒ g(p) is smooth + monotone
// increasing (dρ/dp>0, dρ/dh>0) ⇒ unique root, clean straddle — the fast-path
// precondition the real project_pressure's uniqueness guard verifies.
__host__ __device__ __forceinline__ double node_rho(double p,double h,double zn){
    double pn=(p-PMIN)/(PMAX-PMIN), hn=(h-HMIN)/(HMAX-HMIN);
    return 1.0 + 4.0*pn + 0.5*hn + 0.3*zn + 0.2*pn*zn;
}
// trilinear interp of the density grid on uniform axes (fixed-order).
__device__ double interp_rho(const double* __restrict__ g,double p,double h,double z){
    double fp=(p-PMIN)/(PMAX-PMIN)*(NP-1); if(fp<0)fp=0; if(fp>NP-1)fp=NP-1;
    double fh=(h-HMIN)/(HMAX-HMIN)*(NH-1); if(fh<0)fh=0; if(fh>NH-1)fh=NH-1;
    double fz=z*(NZ-1);                     if(fz<0)fz=0; if(fz>NZ-1)fz=NZ-1;
    int ip=(int)fp, ih=(int)fh, iz=(int)fz;
    if(ip>NP-2)ip=NP-2; if(ih>NH-2)ih=NH-2; if(iz>NZ-2)iz=NZ-2;
    double tp=fp-ip, th=fh-ih, tz=fz-iz;
    #define G(a,b,c) g[((a)*NH+(b))*NZ+(c)]
    double c00=G(ip,ih,iz)*(1-tp)+G(ip+1,ih,iz)*tp;
    double c01=G(ip,ih,iz+1)*(1-tp)+G(ip+1,ih,iz+1)*tp;
    double c10=G(ip,ih+1,iz)*(1-tp)+G(ip+1,ih+1,iz)*tp;
    double c11=G(ip,ih+1,iz+1)*(1-tp)+G(ip+1,ih+1,iz+1)*tp;
    #undef G
    double c0=c00*(1-th)+c10*th, c1=c01*(1-th)+c11*th;
    return c0*(1-tz)+c1*tz;
}
// residual g(p) = rho_tab(p, e_q + p/rho, z) - rho
__device__ __forceinline__ double gres(const double* grid,double p,double e_q,double rho,double z){
    return interp_rho(grid,p,e_q+p/rho,z)-rho;
}
// Illinois regula-falsi + bisection backstop (table_eos.rs illinois_root).
__device__ double illinois(const double* grid,double a,double b,double ga,double gb,
                           double e_q,double rho,double z){
    for(int it=0;it<N_P_ITER_MAX;it++){
        if(fabs(b-a)<=EPS_P_PROJECTION*fmax(fabs(a),fabs(b))) return 0.5*(a+b);
        double denom=gb-ga;
        double p=(denom!=0.0)? b-gb*(b-a)/denom : 0.5*(a+b);
        double lo=fmin(a,b),hi=fmax(a,b);
        if(!(p>lo&&p<hi)) p=0.5*(a+b);
        double gp=gres(grid,p,e_q,rho,z);
        if(gp==0.0) return p;
        if(gp*gb<0.0){ a=b; ga=gb; } else { ga*=0.5; }
        b=p; gb=gp;
    }
    if(fabs(b-a)<=EPS_P_PROJECTION*fmax(fabs(a),fabs(b))*10.0) return 0.5*(a+b);
    double aa=(a<b)?a:b, bb=(a<b)?b:a, gaa=(a<b)?ga:gb;
    for(int it=0;it<N_P_BISECT;it++){
        double m=0.5*(aa+bb);
        if(fabs(bb-aa)<=EPS_P_PROJECTION*fmax(fabs(aa),fabs(bb))) return m;
        double gm=gres(grid,m,e_q,rho,z);
        if(gm==0.0) return m;
        if(gaa*gm<0.0) bb=m; else { aa=m; gaa=gm; }
    }
    return 0.5*(aa+bb);
}
__global__ void project(const double* __restrict__ grid,long n,
                        double* __restrict__ p_out,double* __restrict__ p_true_out){
    long i=(long)blockIdx.x*blockDim.x+threadIdx.x;
    if(i>=n) return;
    // deterministic per-cell target built at a known p_true
    double u=(double)((i*2654435761u)&0xffffff)/16777216.0;
    double v=(double)((i*40503u+12345u)&0xffffff)/16777216.0;
    double w=(double)((i*2246822519u)&0xffffff)/16777216.0;
    double p_true=PMIN+(PMAX-PMIN)*(0.05+0.9*u);
    double h_true=HMIN+(HMAX-HMIN)*(0.05+0.9*v);
    double zn=w;
    double rho=interp_rho(grid,p_true,h_true,zn);
    double e_q=h_true-p_true/rho;
    // admissible bracket: h(p)=e_q+p/rho must lie in the h-envelope
    // (project_pressure's closed-form bracket; outside it the surface refuses).
    double a=fmax(PMIN,(HMIN-e_q)*rho);
    double b=fmin(PMAX,(HMAX-e_q)*rho);
    double ga=gres(grid,a,e_q,rho,zn),gb=gres(grid,b,e_q,rho,zn);
    double p=illinois(grid,a,b,ga,gb,e_q,rho,zn);
    p_out[i]=p; p_true_out[i]=p_true;
}

int main(){
    cudaDeviceProp pr; CK(cudaGetDeviceProperties(&pr,0));
    // fill the synthetic surface (host), upload
    long ng=(long)NP*NH*NZ; double* hg=(double*)malloc(ng*sizeof(double));
    for(int ip=0;ip<NP;ip++)for(int ih=0;ih<NH;ih++)for(int iz=0;iz<NZ;iz++){
        double p=PMIN+(PMAX-PMIN)*ip/(NP-1), h=HMIN+(HMAX-HMIN)*ih/(NH-1), zn=(double)iz/(NZ-1);
        hg[(ip*NH+ih)*NZ+iz]=node_rho(p,h,zn);
    }
    double* dg; CK(cudaMalloc(&dg,ng*sizeof(double))); CK(cudaMemcpy(dg,hg,ng*sizeof(double),cudaMemcpyHostToDevice));
    long n=4194304; size_t nb=n*sizeof(double);
    double *dpo,*dpt; CK(cudaMalloc(&dpo,nb)); CK(cudaMalloc(&dpt,nb));
    int blk=128, grd=(int)((n+blk-1)/blk);
    project<<<grd,blk>>>(dg,n,dpo,dpt); CK(cudaDeviceSynchronize());   // warmup
    int iters=20; cudaEvent_t t0,t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    CK(cudaEventRecord(t0));
    for(int it=0;it<iters;it++) project<<<grd,blk>>>(dg,n,dpo,dpt);
    CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
    float ms; cudaEventElapsedTime(&ms,t0,t1); double sec=ms/1e3/iters;
    printf("Device %s\nEOS project f64: %.2f ms/iter, %.3e cell-projections/s (%ld cells)\n",
           pr.name, ms/iters, (double)n/sec, n);
    // correctness + determinism
    double *r1=(double*)malloc(nb),*pt=(double*)malloc(nb);
    CK(cudaMemcpy(r1,dpo,nb,cudaMemcpyDeviceToHost)); CK(cudaMemcpy(pt,dpt,nb,cudaMemcpyDeviceToHost));
    double maxrel=0; long nbad=0; for(long i=0;i<n;i++){ double e=fabs(r1[i]-pt[i])/pt[i]; if(e>maxrel)maxrel=e; if(e>1e-6)nbad++; }
    printf("cells with rel err > 1e-6: %ld / %ld (%.4f%%)\n", nbad, n, 100.0*nbad/n);
    project<<<grd,blk>>>(dg,n,dpo,dpt); CK(cudaDeviceSynchronize());
    double* r2=(double*)malloc(nb); CK(cudaMemcpy(r2,dpo,nb,cudaMemcpyDeviceToHost));
    printf("recovered p vs p_true: max rel err %.2e (root-find correctness)\n", maxrel);
    printf("determinism (run-twice bitwise): %s\n", memcmp(r1,r2,nb)==0?"BIT-IDENTICAL":"*** DIFFERS ***");
    free(hg);free(r1);free(r2);free(pt); cudaFree(dg);cudaFree(dpo);cudaFree(dpt); return 0;
}
