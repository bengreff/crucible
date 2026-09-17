// S12 GPU spike — hot kernel (a), shared-memory z-pencil formulation. Same
// PPM+HLLC math as kernel_flux.cu, but each block loads a z-tile of prim into
// shared memory ONCE (cons->prim done once per cell, not 7x), and threads read
// their reconstruction stencils from shared — no local-memory spills, coalesced
// global traffic. This is the bandwidth-bound production pattern (S13 residency
// generalizes it to r/theta sweeps + the cut metric). Determinism unchanged:
// gather-only shared reads, fixed-order per-thread arithmetic, no atomics.
#include <cstdio>
#include <cstring>
#include <cuda_runtime.h>

#define NC 7
#define GAMMA 1.4
#define TILE 128
#define HALO 3

#define CK(x) do { cudaError_t e=(x); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %d\n",cudaGetErrorString(e),__LINE__); return 1;} } while(0)

__device__ __forceinline__ double total_energy(const double w[NC]){
    return w[4]/(GAMMA-1.0)+0.5*w[0]*(w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
}
__device__ __forceinline__ double sound_speed_w(const double w[NC]){ return sqrt(GAMMA*w[4]/w[0]); }
__device__ void physical_flux(const double w[NC],int n,double f[NC]){
    double rho=w[0],p=w[4],un=w[n],m=rho*un,e=total_energy(w);
    f[0]=m; f[1]=m*w[1]; f[2]=m*w[2]; f[3]=m*w[3]; f[n]+=p; f[4]=un*(e+p); f[5]=m*w[5]; f[6]=m*w[6];
}
__device__ void hllc_flux(const double wl[NC],const double wr[NC],int n,double f[NC]){
    double rho_l=wl[0],p_l=wl[4],rho_r=wr[0],p_r=wr[4],un_l=wl[n],un_r=wr[n];
    double c_l=sound_speed_w(wl),c_r=sound_speed_w(wr);
    double sql=sqrt(rho_l),sqr=sqrt(rho_r),inv=1.0/(sql+sqr);
    double u1=(sql*wl[1]+sqr*wr[1])*inv,u2=(sql*wl[2]+sqr*wr[2])*inv,u3=(sql*wl[3]+sqr*wr[3])*inv;
    double urn=(n==1?u1:(n==2?u2:u3));
    double h_l=(total_energy(wl)+p_l)/rho_l,h_r=(total_energy(wr)+p_r)/rho_r;
    double h_roe=(sql*h_l+sqr*h_r)*inv,q2=u1*u1+u2*u2+u3*u3;
    double arg=(GAMMA-1.0)*(h_roe-0.5*q2),c_roe=sqrt(arg>0.0?arg:0.0);
    double s_l=fmin(un_l-c_l,urn-c_roe),s_r=fmax(un_r+c_r,urn+c_roe);
    double ml=rho_l*(s_l-un_l),mr=rho_r*(s_r-un_r),s_m=(mr*un_r-ml*un_l+p_l-p_r)/(mr-ml);
    if(s_l>=0.0){ physical_flux(wl,n,f); return; }
    if(s_r<=0.0){ physical_flux(wr,n,f); return; }
    const double* w=(s_m>=0.0)?wl:wr; double s_k=(s_m>=0.0)?s_l:s_r;
    double rho=w[0],p=w[4],un=w[n],e=total_energy(w),p_star=rho*(un-s_k)*(un-s_m)+p;
    double u_k[NC]={rho,rho*w[1],rho*w[2],rho*w[3],e,rho*w[5],rho*w[6]};
    double fac=(s_k-un)/(s_k-s_m),rho_s=rho*fac,u_s[NC];
    u_s[0]=rho_s; u_s[1]=rho_s*w[1]; u_s[2]=rho_s*w[2]; u_s[3]=rho_s*w[3]; u_s[n]=rho_s*s_m;
    u_s[4]=fac*e+(p_star*s_m-p*un)/(s_k-s_m); u_s[5]=rho_s*w[5]; u_s[6]=rho_s*w[6];
    double fk[NC]; physical_flux(w,n,fk);
    for(int i=0;i<NC;i++) f[i]=fk[i]+s_k*(u_s[i]-u_k[i]);
}
__device__ __forceinline__ double mc_slope(double wm,double w0,double wp){
    double dl=w0-wm,dr=wp-w0; if(dl*dr<=0.0) return 0.0;
    double dc=0.5*(dl+dr),a=fabs(dc),b=2.0*fabs(dl),c=2.0*fabs(dr);
    double mn=a<b?a:b; mn=mn<c?mn:c; return mn*(dc>0.0?1.0:(dc<0.0?-1.0:0.0));
}
// PPM edges (lo/hi) of the cell whose 5-cell window is p[s-2..s+2] in shared;
// sm is component-major: sm[k*W + idx].
__device__ void edge_hilo_sm(const double* sm,int W,int s,double lo[NC],double hi[NC]){
    for(int k=0;k<NC;k++){
        const double* col=sm+(long)k*W;
        double s_im1=mc_slope(col[s-2],col[s-1],col[s]);
        double s_i  =mc_slope(col[s-1],col[s],col[s+1]);
        double s_ip1=mc_slope(col[s],col[s+1],col[s+2]);
        double if_im1=0.5*(col[s-1]+col[s])-(s_i-s_im1)/6.0;
        double if_i  =0.5*(col[s]+col[s+1])-(s_ip1-s_i)/6.0;
        double c=col[s],lov=if_im1,hiv=if_i;
        if((hiv-c)*(c-lov)<=0.0){lov=c;hiv=c;}
        else{ double d=hiv-lov,six=6.0*(c-0.5*(lov+hiv));
              if(d*six>d*d) lov=3.0*c-2.0*hiv; else if(d*six<-(d*d)) hiv=3.0*c-2.0*lov; }
        lo[k]=lov; hi[k]=hiv;
    }
}
__global__ void sweep_z_sm(const double* __restrict__ cons, double* __restrict__ out,
                           int nr,int nz,double dt_dz){
    const int W=TILE+2*HALO;
    __shared__ double sm[NC*W];           // prim tile, component-major
    int ir=blockIdx.x;                     // one radial index per block-column
    int tile0=blockIdx.y*TILE;             // first interior z of this tile
    long ncell=(long)nr*nz;
    // cooperative load: shared idx t in [0,W) maps to global z = tile0-HALO+t
    for(int t=threadIdx.x;t<W;t+=blockDim.x){
        int gz=tile0-HALO+t;
        if(gz<0||gz>=nz){ for(int k=0;k<NC;k++) sm[k*W+t]=0.0; continue; }
        long cell=(long)gz*nr+ir;
        double u[NC]; for(int k=0;k<NC;k++) u[k]=cons[(long)k*ncell+cell];
        double rho=u[0],iv=1.0/rho,ur=u[1]*iv,ut=u[2]*iv,uz=u[3]*iv;
        double ke=0.5*rho*(ur*ur+ut*ut+uz*uz);
        double w[NC]={rho,ur,ut,uz,(GAMMA-1.0)*(u[4]-ke),u[5]*iv,u[6]*iv};
        for(int k=0;k<NC;k++) sm[k*W+t]=w[k];
    }
    __syncthreads();
    int li=threadIdx.x; int iz=tile0+li; int s=li+HALO;   // shared index of iz
    if(iz<3||iz>=nz-3||li>=TILE) return;
    double loB[NC],hiB[NC],loA[NC],hiA[NC],loC[NC],hiC[NC];
    edge_hilo_sm(sm,W,s,  loB,hiB);   // cell iz
    edge_hilo_sm(sm,W,s-1,loA,hiA);   // cell iz-1
    edge_hilo_sm(sm,W,s+1,loC,hiC);   // cell iz+1
    double fm[NC],fp[NC];
    hllc_flux(hiA,loB,3,fm);          // flux iz-1/2
    hllc_flux(hiB,loC,3,fp);          // flux iz+1/2
    long cell=(long)iz*nr+ir;
    for(int k=0;k<NC;k++) out[(long)k*ncell+cell]=cons[(long)k*ncell+cell]-dt_dz*(fp[k]-fm[k]);
}

int main(){
    cudaDeviceProp p; CK(cudaGetDeviceProperties(&p,0));
    int nr=2048,nz=2048; long ncell=(long)nr*nz,ninner=(long)nr*(nz-6);
    size_t bytes=(size_t)NC*ncell*sizeof(double);
    printf("Device %s | field %dx%d = %ld cells\n",p.name,nr,nz,ncell);
    double* h=(double*)malloc(bytes);
    for(long c=0;c<ncell;c++){ int iz=c/nr,ir=c%nr; double x=(double)iz/nz,y=(double)ir/nr;
        double rho=1.0+0.3*sin(6.2831*x)+0.2*cos(6.2831*y)+(iz>nz/2?0.5:0.0);
        double ur=0.1*sin(3.0*x),ut=0.05*cos(2.0*y),uz=0.2+0.1*sin(4.0*y);
        double pr=1.0+0.4*cos(6.2831*x)+(iz>nz/2?0.3:0.0),cc=0.5+0.5*sin(5.0*y),bb=0.3;
        double E=pr/(GAMMA-1.0)+0.5*rho*(ur*ur+ut*ut+uz*uz);
        double U[NC]={rho,rho*ur,rho*ut,rho*uz,E,rho*cc,rho*bb};
        for(int k=0;k<NC;k++) h[(long)k*ncell+c]=U[k]; }
    double *d_a,*d_b; CK(cudaMalloc(&d_a,bytes)); CK(cudaMalloc(&d_b,bytes));
    dim3 blk(TILE), grd(nr,(nz+TILE-1)/TILE); double dt_dz=0.05;
    CK(cudaMemcpy(d_a,h,bytes,cudaMemcpyHostToDevice));
    sweep_z_sm<<<grd,blk>>>(d_a,d_b,nr,nz,dt_dz); CK(cudaDeviceSynchronize());
    int iters=200; cudaEvent_t t0,t1; cudaEventCreate(&t0); cudaEventCreate(&t1);
    CK(cudaMemcpy(d_a,h,bytes,cudaMemcpyHostToDevice));
    double *A=d_a,*B=d_b; CK(cudaEventRecord(t0));
    for(int it=0;it<iters;it++){ sweep_z_sm<<<grd,blk>>>(A,B,nr,nz,dt_dz); double*t=A;A=B;B=t; }
    CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
    float ms; cudaEventElapsedTime(&ms,t0,t1); double sec=ms/1e3;
    printf("sweep_z_sm f64: %.2f ms/iter, %.3e cell-updates/s\n",ms/iters,(double)ninner*iters/sec);
    double* r1=(double*)malloc(bytes); CK(cudaMemcpy(r1,A,bytes,cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(d_a,h,bytes,cudaMemcpyHostToDevice)); A=d_a;B=d_b;
    for(int it=0;it<iters;it++){ sweep_z_sm<<<grd,blk>>>(A,B,nr,nz,dt_dz); double*t=A;A=B;B=t; }
    CK(cudaDeviceSynchronize());
    double* r2=(double*)malloc(bytes); CK(cudaMemcpy(r2,A,bytes,cudaMemcpyDeviceToHost));
    printf("determinism (run-twice bitwise): %s\n", memcmp(r1,r2,bytes)==0?"BIT-IDENTICAL":"*** DIFFERS ***");
    free(h);free(r1);free(r2);cudaFree(d_a);cudaFree(d_b); return 0;
}
