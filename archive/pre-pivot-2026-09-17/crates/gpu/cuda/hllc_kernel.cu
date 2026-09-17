// crucible-gpu — HLLC-Batten flux kernel (GammaLaw), the GPU port cross-checked
// against crucible_solvers::euler::hllc_flux. Device math is byte-for-formula
// identical to hllc.rs / euler/mod.rs; cross-device bit-identity is NOT promised
// (FMA-contraction differs CPU vs GPU) — the Rust harness measures the tolerance.
// Input prim stride is NPRIM=9 (the CPU Prim), of which slots 0..6 are read
// (rho, u_r, u_theta, u_z, p, c, b); output stride is NCOMP=7.
#include <cuda_runtime.h>

#define NPRIM 9
#define NC 7
#define GAMMA 1.4

__device__ __forceinline__ double total_energy(const double* w){
    return w[4]/(GAMMA-1.0)+0.5*w[0]*(w[1]*w[1]+w[2]*w[2]+w[3]*w[3]);
}
__device__ __forceinline__ double sound_speed_w(const double* w){ return sqrt(GAMMA*w[4]/w[0]); }
__device__ void physical_flux(const double* w,int n,double f[NC]){
    double rho=w[0],p=w[4],un=w[n],m=rho*un,e=total_energy(w);
    f[0]=m; f[1]=m*w[1]; f[2]=m*w[2]; f[3]=m*w[3]; f[n]+=p; f[4]=un*(e+p); f[5]=m*w[5]; f[6]=m*w[6];
}
__device__ void hllc_flux(const double* wl,const double* wr,int n,double f[NC]){
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
__global__ void hllc_batch(const double* __restrict__ wl,const double* __restrict__ wr,
                           int n,int dir,double* __restrict__ out){
    long i=(long)blockIdx.x*blockDim.x+threadIdx.x; if(i>=n) return;
    double f[NC]; hllc_flux(wl+i*NPRIM, wr+i*NPRIM, dir, f);
    for(int k=0;k<NC;k++) out[i*NC+k]=f[k];
}

extern "C" void gpu_hllc(const double* wl,const double* wr,int n,int dir,double* out){
    double *dl,*dr,*dout;
    cudaMalloc(&dl,(size_t)n*NPRIM*sizeof(double));
    cudaMalloc(&dr,(size_t)n*NPRIM*sizeof(double));
    cudaMalloc(&dout,(size_t)n*NC*sizeof(double));
    cudaMemcpy(dl,wl,(size_t)n*NPRIM*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(dr,wr,(size_t)n*NPRIM*sizeof(double),cudaMemcpyHostToDevice);
    hllc_batch<<<(n+255)/256,256>>>(dl,dr,n,dir,dout);
    cudaDeviceSynchronize();
    cudaMemcpy(out,dout,(size_t)n*NC*sizeof(double),cudaMemcpyDeviceToHost);
    cudaFree(dl); cudaFree(dr); cudaFree(dout);
}
