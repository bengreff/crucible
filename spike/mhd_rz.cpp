// Architecture spike, 30 September 2026. Throwaway code: not production, not a framework.
// Question: can the body-fitted RZ layout of core/flow carry magnetic fields?
// Ideal MHD with swirl and toroidal field on the same node/cell/face layout and exact ring
// measures as core/flow. Field unit b = B/sqrt(mu0), so magnetic pressure is b^2/2 and J x b
// is the SI force density. HLLD fluxes (Miyoshi & Kusano 2005). Three div B treatments:
//   none : cell-centred B, normal component averaged at faces;
//   glm  : cell-centred B with mixed GLM cleaning (Dedner et al. 2002);
//   ct   : poloidal B from the flux function psi = r A_theta stored at mesh nodes, advanced by
//          node EMFs averaged from the HLLD face fluxes. Face fluxes 2 pi dpsi are exactly
//          solenoidal on any body-fitted quadrilateral mesh. B_theta stays cell-centred.
//   split: ct for the induced field b1 only; the coil (vacuum) field b0 is steady and curl-free,
//          analytic at faces (Tanaka 1994). Energy slot holds E1 = p/(g-1) + rho u^2/2 + b1^2/2,
//          so pressure never comes from subtracting b0^2/2. The b0 stress T(b0) (divergence-free)
//          is dropped; the cross stress b0 b1 becomes a cell body force (J1 x b0, the force of the
//          coil field on the plasma), whose axial integral is bodyAxialForce.
#include "core/flow.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi=std::numbers::pi, mu0=4e-7*pi;
constexpr int NV=10;
using U=std::array<double,NV>;
// Conserved slots; primitives reuse them as rho,uz,ur,ut,p,bz,br,bt,psi,s. ENT carries rho s with
// s = p / rho^gamma, the auxiliary (entropy) variable for pressure recovery at low beta.
enum {RHO,MZ,MR,MT,EN,BZ,BR,BT,PS,ENT};
double sq(double x) { return x*x; }
double minmod(double a,double b) { return a*b<=0?0:std::copysign(std::min(std::abs(a),std::abs(b)),a); }
enum class Clean { None, GLM, CT, Split };
// Magnetic condition at the side wall for the induced field b1 (split only):
//   transparent: ghost copies b1 (the first spike; sliding footpoints, no exterior field);
//   insulating : ghost b1 is the exterior vacuum field matched to the wall flux psi1 (no wall
//                current, coil current held fixed, so b0 is maintained and b1 decays outside);
//   conducting : perfectly conducting wall, flux through it frozen (psi1 = 0 at wall nodes).
enum class WallB { Transparent, Insulating, Conducting };
const char* name(WallB w) { return w==WallB::Transparent?"transparent":w==WallB::Insulating?"insulating":"conducting"; }
const char* name(Clean c) { return c==Clean::None?"none":c==Clean::GLM?"glm":c==Clean::CT?"ct":"split"; }

// Same construction as core Mesh, plus an optional inner radius (annulus).
struct Grid {
    int nz,nr; double dz,inner;
    std::vector<double> rn,vol,meas;
    Grid(int nz_,int nr_,double length,const std::function<double(double)>& outer,double inner_=0)
        :nz(nz_),nr(nr_),dz(length/nz_),inner(inner_) {
        rn.resize(static_cast<std::size_t>(nz+1)*(nr+1));
        for(int i=0;i<=nz;++i) { double R=outer(i*dz); for(int j=0;j<=nr;++j) rn[node(i,j)]=inner+(R-inner)*j/nr; }
        vol.resize(static_cast<std::size_t>(nz)*nr); meas.resize(vol.size());
        for(int i=0;i<nz;++i) for(int j=0;j<nr;++j) {
            double a0=r(i,j+1),a1=r(i+1,j+1),b0=r(i,j),b1=r(i+1,j);
            vol[cell(i,j)]=pi*dz/3*(a0*a0+a0*a1+a1*a1-b0*b0-b0*b1-b1*b1);
            meas[cell(i,j)]=pi*dz*(a0+a1-b0-b1);
        }
    }
    std::size_t node(int i,int j) const { return static_cast<std::size_t>(i)*(nr+1)+j; }
    std::size_t cell(int i,int j) const { return static_cast<std::size_t>(i)*nr+j; }
    double r(int i,int j) const { return rn[node(i,j)]; }
    double axialArea(int i,int j) const { return pi*(sq(r(i,j+1))-sq(r(i,j))); }
    std::array<double,2> radialArea(int i,int j) const {
        double r0=r(i,j),r1=r(i+1,j); return {-pi*(r0+r1)*(r1-r0),pi*(r0+r1)*dz};
    }
};

// Dual energy (Bryan et al. 1995; FLASH): total energy stays the conserved variable, but where the
// thermal part is below dualEta of it the pressure is recovered from the advected entropy instead.
bool dualEnergy=false; double dualEta=1e-3;
U toCons(const U& w,double g) {
    U u=w; u[ENT]=w[4]/std::pow(w[RHO],g-1); u[MZ]=w[RHO]*w[1]; u[MR]=w[RHO]*w[2]; u[MT]=w[RHO]*w[3];
    u[EN]=w[4]/(g-1)+0.5*w[RHO]*(sq(w[1])+sq(w[2])+sq(w[3]))+0.5*(sq(w[BZ])+sq(w[BR])+sq(w[BT]));
    return u;
}
bool toPrim(const U& u,double g,U& w) {
    for(double x:u) if(!std::isfinite(x)) return false;
    if(!(u[RHO]>0)) return false;
    w=u; w[1]=u[MZ]/u[RHO]; w[2]=u[MR]/u[RHO]; w[3]=u[MT]/u[RHO];
    w[4]=(g-1)*(u[EN]-0.5*(sq(u[MZ])+sq(u[MR])+sq(u[MT]))/u[RHO]-0.5*(sq(u[BZ])+sq(u[BR])+sq(u[BT])));
    w[ENT]=u[ENT]/u[RHO];
    if(dualEnergy && !(w[4]>dualEta*(g-1)*u[EN])) w[4]=u[ENT]*std::pow(u[RHO],g-1);
    return w[4]>0;
}

// ---- 1D HLLD in the face frame: (rho, un, ut1, ut2, p, bt1, bt2) with normal field bx.
struct S1 { double rho,u,v,w,p,by,bz; };
using V7=std::array<double,7>;
double ptot(const S1& s,double bx) { return s.p+0.5*(bx*bx+s.by*s.by+s.bz*s.bz); }
V7 cons1(const S1& s,double bx,double g) {
    return {s.rho,s.rho*s.u,s.rho*s.v,s.rho*s.w,
        s.p/(g-1)+0.5*s.rho*(s.u*s.u+s.v*s.v+s.w*s.w)+0.5*(bx*bx+s.by*s.by+s.bz*s.bz),s.by,s.bz};
}
V7 flux1(const S1& s,double bx,double g) {
    auto c=cons1(s,bx,g); double pt=ptot(s,bx),vb=s.u*bx+s.v*s.by+s.w*s.bz;
    return {s.rho*s.u,s.rho*s.u*s.u+pt-bx*bx,s.rho*s.u*s.v-bx*s.by,s.rho*s.u*s.w-bx*s.bz,
        (c[4]+pt)*s.u-bx*vb,s.by*s.u-bx*s.v,s.bz*s.u-bx*s.w};
}
double fastSpeed(const S1& s,double bx,double g) {
    double a2=g*s.p/s.rho,b2=(bx*bx+s.by*s.by+s.bz*s.bz)/s.rho;
    return std::sqrt(0.5*(a2+b2+std::sqrt(std::max(0.0,sq(a2+b2)-4*a2*bx*bx/s.rho))));
}
long hllFallbacks=0;
bool forceHll=false;
V7 hlld(const S1& L,const S1& R,double bx,double g) {
    double cf=std::max(fastSpeed(L,bx,g),fastSpeed(R,bx,g));
    double SL=std::min(L.u,R.u)-cf,SR=std::max(L.u,R.u)+cf;
    auto UL=cons1(L,bx,g),UR=cons1(R,bx,g),FL=flux1(L,bx,g),FR=flux1(R,bx,g);
    if(SL>=0) return FL;
    if(SR<=0) return FR;
    auto hll=[&] { ++hllFallbacks; V7 f{}; for(int k=0;k<7;++k) f[k]=(SR*FL[k]-SL*FR[k]+SL*SR*(UR[k]-UL[k]))/(SR-SL); return f; };
    if(forceHll) return hll();
    double dl=SL-L.u,dr=SR-R.u,den=dr*R.rho-dl*L.rho;
    double SM=(dr*R.rho*R.u-dl*L.rho*L.u-ptot(R,bx)+ptot(L,bx))/den;
    double pts=(dr*R.rho*ptot(L,bx)-dl*L.rho*ptot(R,bx)+L.rho*R.rho*dr*dl*(R.u-L.u))/den;
    struct St { double rho,v,w,by,bz,e; };
    bool bad=false;
    auto star=[&](const S1& s,const V7& u,double S) {
        double d=S-s.u; St o{s.rho*d/(S-SM),s.v,s.w,s.by,s.bz,0};
        double den2=s.rho*d*(S-SM)-bx*bx;
        // Degenerate (tangential field and switch-on) case: keep tangential state.
        if(std::abs(den2)>1e-10*(std::abs(s.rho*d*(S-SM))+bx*bx)) {
            double f1=bx*(SM-s.u)/den2,f2=(s.rho*d*d-bx*bx)/den2;
            o.v=s.v-s.by*f1; o.w=s.w-s.bz*f1; o.by=s.by*f2; o.bz=s.bz*f2;
        }
        double vb=s.u*bx+s.v*s.by+s.w*s.bz,vbs=SM*bx+o.v*o.by+o.w*o.bz;
        o.e=(d*u[4]-ptot(s,bx)*s.u+pts*SM+bx*(vb-vbs))/(S-SM);
        double p=(g-1)*(o.e-0.5*o.rho*(SM*SM+o.v*o.v+o.w*o.w)-0.5*(bx*bx+o.by*o.by+o.bz*o.bz));
        if(!(o.rho>0) || !(p>0) || !std::isfinite(o.e)) bad=true;
        return o;
    };
    St sL=star(L,UL,SL),sR=star(R,UR,SR);
    if(bad || !(SL<SM && SM<SR)) return hll();
    auto pack=[&](const St& o,double rho) { return V7{rho,rho*SM,rho*o.v,rho*o.w,o.e,o.by,o.bz}; };
    V7 UsL=pack(sL,sL.rho),UsR=pack(sR,sR.rho),FsL{},FsR{};
    for(int k=0;k<7;++k) { FsL[k]=FL[k]+SL*(UsL[k]-UL[k]); FsR[k]=FR[k]+SR*(UsR[k]-UR[k]); }
    double rl=std::sqrt(sL.rho),rr=std::sqrt(sR.rho),sg=bx>=0?1:-1;
    double SLs=SM-std::abs(bx)/rl,SRs=SM+std::abs(bx)/rr;
    if(SLs>=0) return FsL;
    if(SRs<=0) return FsR;
    St ss{0,(rl*sL.v+rr*sR.v+(sR.by-sL.by)*sg)/(rl+rr),(rl*sL.w+rr*sR.w+(sR.bz-sL.bz)*sg)/(rl+rr),
        (rl*sR.by+rr*sL.by+rl*rr*(sR.v-sL.v)*sg)/(rl+rr),(rl*sR.bz+rr*sL.bz+rl*rr*(sR.w-sL.w)*sg)/(rl+rr),0};
    double vbss=SM*bx+ss.v*ss.by+ss.w*ss.bz;
    double vbL=SM*bx+sL.v*sL.by+sL.w*sL.bz,vbR=SM*bx+sR.v*sR.by+sR.w*sR.bz;
    V7 F{};
    if(SM>=0) { St o=ss; o.e=sL.e-rl*(vbL-vbss)*sg; auto Uss=pack(o,sL.rho); for(int k=0;k<7;++k) F[k]=FsL[k]+SLs*(Uss[k]-UsL[k]); }
    else      { St o=ss; o.e=sR.e+rr*(vbR-vbss)*sg; auto Uss=pack(o,sR.rho); for(int k=0;k<7;++k) F[k]=FsR[k]+SRs*(Uss[k]-UsR[k]); }
    return F;
}
struct FaceFlux { U f; double emf; };
// Face with unit normal (nz,nr) in the meridional plane; t1=(-nr,nz), t2=theta.
// emf is (u x b)_theta = flux of b_t1, the axisymmetric CT electric field (sign: dA_theta/dt).
FaceFlux mhdFlux(const U& l,const U& r,double nz,double nr,double ch,double g) {
    auto un=[&](const U& w){return w[1]*nz+w[2]*nr;}; auto ut=[&](const U& w){return -w[1]*nr+w[2]*nz;};
    auto bn=[&](const U& w){return w[BZ]*nz+w[BR]*nr;}; auto bt=[&](const U& w){return -w[BZ]*nr+w[BR]*nz;};
    double b=0.5*(bn(l)+bn(r)),psi=0.5*(l[PS]+r[PS]);
    if(ch>0) { b-=0.5*(r[PS]-l[PS])/ch; psi-=0.5*ch*(bn(r)-bn(l)); }
    S1 L{l[RHO],un(l),ut(l),l[3],l[4],bt(l),l[BT]},R{r[RHO],un(r),ut(r),r[3],r[4],bt(r),r[BT]};
    auto f=hlld(L,R,b,g);
    FaceFlux o{};
    o.f[RHO]=f[0]; o.f[MZ]=f[1]*nz-f[2]*nr; o.f[MR]=f[1]*nr+f[2]*nz; o.f[MT]=f[3]; o.f[EN]=f[4];
    o.f[BZ]=psi*nz-f[5]*nr; o.f[BR]=psi*nr+f[5]*nz; o.f[BT]=f[6]; o.f[PS]=ch*ch*b; o.emf=f[5];
    o.f[ENT]=f[0]*(f[0]>0?l[ENT]:r[ENT]);  // upwind entropy with the mass flux
    return o;
}

// Current loop of radius a at z=zc carrying I: psi = r A_phi in b units (T m^2 / sqrt(mu0)).
double loopPsi(double r,double z,double a,double zc,double I) {
    if(r<=0) return 0;
    double k2=4*a*r/(sq(a+r)+sq(z-zc));
    double A=1,B=std::sqrt(1-k2),sum=k2/2,p2=0.5;
    for(int n=0;n<60;++n) { double c=(A-B)/2,an=(A+B)/2; B=std::sqrt(A*B); A=an; p2*=2; sum+=p2*c*c; if(c<1e-17) break; }
    double K=pi/(2*A),E=K*(1-sum),k=std::sqrt(k2);
    double Aphi=mu0*I/(pi*k)*std::sqrt(a/r)*((1-k2/2)*K-E);
    return r*Aphi/std::sqrt(mu0);
}

// Exact ideal-gas Riemann density (Toro ch. 4), for the field-aligned check.
struct G1 { double r,u,p; };
double exactDensity(G1 L,G1 R,double g,double xi) {
    auto f=[&](double p,const G1& s,double& d) {
        double a=std::sqrt(g*s.p/s.r);
        if(p>s.p) { double A=2/((g+1)*s.r),B=(g-1)/(g+1)*s.p; d=std::sqrt(A/(B+p))*(1-(p-s.p)/(2*(B+p))); return (p-s.p)*std::sqrt(A/(p+B)); }
        d=1/(s.r*a)*std::pow(p/s.p,-(g+1)/(2*g)); return 2*a/(g-1)*(std::pow(p/s.p,(g-1)/(2*g))-1);
    };
    double p=0.5*(L.p+R.p),dL,dR;
    for(int n=0;n<100;++n) { double F=f(p,L,dL)+f(p,R,dR)+R.u-L.u; p=std::max(1e-14,p-F/(dL+dR)); }
    double fl=f(p,L,dL),fr=f(p,R,dR),us=0.5*(L.u+R.u)+0.5*(fr-fl);
    auto side=[&](G1 s,double x,double sign) {  // sign=+1 left, -1 right (mirror)
        double a=std::sqrt(g*s.p/s.r),u=sign*s.u,xx=sign*x,uss=sign*us;
        if(p>s.p) { double S=u-a*std::sqrt((g+1)/(2*g)*p/s.p+(g-1)/(2*g));
            return xx<S?s.r:s.r*((p/s.p+(g-1)/(g+1))/((g-1)/(g+1)*p/s.p+1)); }
        double as=a*std::pow(p/s.p,(g-1)/(2*g));
        if(xx<u-a) return s.r;
        if(xx>uss-as) return s.r*std::pow(p/s.p,1/g);
        double c=2/(g+1)*(a+(g-1)/2*(u-xx)); return s.r*std::pow(c/a,2/(g-1));
    };
    return xi<=us?side(L,xi,1):side(R,xi,-1);
}

// Exterior vacuum field for the insulating wall. Outside the nozzle wall there is no plasma and no
// current except the coil, whose field is b0, so the induced field there is curl-free: in flux form
// div((1/r) grad psi1) = 0 in the meridional plane. Solved with linear triangles on an annulus from
// the wall r = R(z) out to a far radius Rf (psi1 = 0 there, a distant conducting shell) over the
// nozzle length, with b_r1 = 0 (natural condition) on the end planes z = 0 and z = L. Upstream
// plenum and downstream plume are not vacuum and are not modelled here. The problem is linear and
// fixed, so it is reduced once to a matrix from the wall node values psi1_k to the exterior field
// (b_z, b_r) at each wall face (one-sided second-order radial derivative).
struct Vacuum {
    int nz{},nv{}; double Rf{};
    std::vector<double> Dz,Dr;  // nz x (nz+1), row-major
    std::array<double,2> field(int i,const std::vector<double>& wallPsi) const {
        double bz=0,br=0; const double* a=&Dz[static_cast<std::size_t>(i)*(nz+1)]; const double* b=&Dr[static_cast<std::size_t>(i)*(nz+1)];
        for(int k=0;k<=nz;++k) { bz+=a[k]*wallPsi[k]; br+=b[k]*wallPsi[k]; }
        return {bz,br};
    }
};
Vacuum buildVacuum(int nz,double dz,const std::function<double(double)>& wall,double Rf,int nv) {
    Vacuum v; v.nz=nz; v.nv=nv; v.Rf=Rf;
    // Geometric radial spacing with the first layer about dz thick at the mean wall radius.
    double Rm=0; for(int i=0;i<=nz;++i) Rm+=wall(i*dz)/(nz+1);
    double lo=1.0000001,hi=2; for(int n=0;n<200;++n) { double q=(lo+hi)/2; if((q-1)/(std::pow(q,nv)-1)*(Rf-Rm)>dz) lo=q; else hi=q; }
    double q=(lo+hi)/2; std::vector<double> sfrac(nv+1); for(int j=0;j<=nv;++j) sfrac[j]=(std::pow(q,j)-1)/(std::pow(q,nv)-1);
    auto R=[&](int i,int j){ double w=wall(i*dz); return w+(Rf-w)*sfrac[j]; };
    // Unknowns: nodes j = 1..nv-1, ordered i-major; half bandwidth nv.
    const int m=nv-1,N=(nz+1)*m,bw=nv;
    auto id=[&](int i,int j){ return i*m+(j-1); };
    std::vector<double> A(static_cast<std::size_t>(N)*(bw+1),0.0);  // A[n*(bw+1)+(n-c)] for c<=n
    std::vector<std::array<double,3>> KIB;  // (unknown n, wall node k, K entry)
    auto add=[&](int i0,int j0,int i1,int j1,double val) {
        if(j0==0||j0==nv) return;  // rows only for unknowns
        int a=id(i0,j0);
        if(j1==nv) return;
        if(j1==0) { KIB.push_back({double(a),double(i1),val}); return; }
        int b=id(i1,j1); if(b>a) return;  // lower triangle
        A[static_cast<std::size_t>(a)*(bw+1)+(a-b)]+=val;
    };
    auto tri=[&](std::array<std::array<int,2>,3> t) {
        double z[3],r[3]; for(int k=0;k<3;++k) { z[k]=t[k][0]*dz; r[k]=R(t[k][0],t[k][1]); }
        double a2=(z[1]-z[0])*(r[2]-r[0])-(z[2]-z[0])*(r[1]-r[0]),rc=(r[0]+r[1]+r[2])/3;
        double gz[3],gr[3]; for(int k=0;k<3;++k) { gz[k]=(r[(k+1)%3]-r[(k+2)%3])/a2; gr[k]=(z[(k+2)%3]-z[(k+1)%3])/a2; }
        for(int a=0;a<3;++a) for(int b=0;b<3;++b) add(t[a][0],t[a][1],t[b][0],t[b][1],0.5*std::abs(a2)/rc*(gz[a]*gz[b]+gr[a]*gr[b]));
    };
    for(int i=0;i<nz;++i) for(int j=0;j<nv;++j) { tri({{{i,j},{i+1,j},{i+1,j+1}}}); tri({{{i,j},{i+1,j+1},{i,j+1}}}); }
    // Banded Cholesky in place.
    auto L=[&](int r_,int c)->double& { return A[static_cast<std::size_t>(r_)*(bw+1)+(r_-c)]; };
    for(int j=0;j<N;++j) for(int k=std::max(0,j-bw);k<=j;++k) {
        double s=L(j,k); for(int t=std::max(0,j-bw);t<k;++t) s-=L(j,t)*L(k,t);
        if(k==j) { if(!(s>0)) throw std::runtime_error("vacuum matrix not positive definite"); L(j,j)=std::sqrt(s); } else L(j,k)=s/L(k,k);
    }
    v.Dz.assign(static_cast<std::size_t>(nz)*(nz+1),0); v.Dr.assign(v.Dz.size(),0);
    std::vector<double> x(N),full(static_cast<std::size_t>(nz+1)*(nv+1));
    for(int k=0;k<=nz;++k) {
        std::fill(x.begin(),x.end(),0.0);
        for(const auto& e:KIB) if(int(e[1])==k) x[int(e[0])]-=e[2];
        for(int j=0;j<N;++j) { double s=x[j]; for(int t=std::max(0,j-bw);t<j;++t) s-=L(j,t)*x[t]; x[j]=s/L(j,j); }
        for(int j=N-1;j>=0;--j) { double s=x[j]; for(int t=j+1;t<=std::min(N-1,j+bw);++t) s-=L(t,j)*x[t]; x[j]=s/L(j,j); }
        auto P=[&](int i,int j){ if(j==0) return i==k?1.0:0.0; if(j==nv) return 0.0; return x[id(i,j)]; };
        for(int i=0;i<nz;++i) {
            // psi_r at the wall from the quadratic through the first three exterior nodes (second
            // order), averaged over the face's two columns; psi_z from the wall nodes, minus the
            // wall-slope part: d psi/dz along the wall = psi_z + R'(z) psi_r.
            auto dr=[&](int c) { double h1=R(c,1)-R(c,0),h2=R(c,2)-R(c,1);
                return -(2*h1+h2)/(h1*(h1+h2))*P(c,0)+(h1+h2)/(h1*h2)*P(c,1)-h1/(h2*(h1+h2))*P(c,2); };
            double gr=0.5*(dr(i)+dr(i+1)),gz=(P(i+1,0)-P(i,0))/dz-(R(i+1,0)-R(i,0))/dz*gr;
            double rf=0.5*(R(i,0)+R(i+1,0));
            v.Dz[static_cast<std::size_t>(i)*(nz+1)+k]=gr/rf; v.Dr[static_cast<std::size_t>(i)*(nz+1)+k]=-gz/rf;
        }
    }
    return v;
}

struct Rates { double mass{},energy{},inMom{},outMom{},wallMom{},inMass{},outMass{},body{}; };
struct Report {
    double t{},mass{},energy{},Pz{},massErr{},energyErr{},momErr{},deviceThrust{},exitPlaneThrust{};
    double inMom{},outMom{},wallMom{},ambient{},inMdot{},outMdot{};
    double divMax{},divMean{},ctFaceDivMax{},maxMach{},minBeta{},lorentzZ{},coilLorentzZ{},maxSpeed{},body{};
    std::uint64_t steps{},rejected{},robustSteps{};
    double sync{};
};

class Mhd {
public:
    Grid m; double g; Clean clean; bool nozzle{false};
    double p0{},T0{},Rg{287.05},pBack{},cfl{0.4},glmAlpha{0.4};
    bool axisParity{true};
    WallB wallB{WallB::Transparent}; Vacuum vac; double vacuumRadius{1.0}; int vacuumLayers{60};
    bool robust{false};       // dual energy + first-order HLL retry of a failed step before halving dt
    bool firstOrder{false};   // current attempt uses zero slopes
    std::vector<double> wallPsi;
    std::uint64_t robustSteps{}; double intSync{};
    std::vector<U> u,stage,next,du,w,sz,sr;
    std::vector<double> psi,psiStage,psiNext,dpsi,emfA,emfR;
    std::vector<std::array<double,2>> b0;  // applied (coil) field at cells
    std::vector<std::array<double,2>> b0A,b0R,cross;  // split: b0 at axial/radial face centres; cell cross force
    bool ct() const { return clean==Clean::CT || clean==Clean::Split; }
    bool split() const { return clean==Clean::Split; }
    std::array<double,2> total(const U& x,std::size_t q) const { return split()?std::array<double,2>{x[BZ]+b0[q][0],x[BR]+b0[q][1]}:std::array<double,2>{x[BZ],x[BR]}; }
    double t{},ch{};
    std::uint64_t steps{},rejected{};
    Rates last{};
    double M0{},E0{},P0{},intM{},intE{},intP{},intPg{};

    Mhd(Grid grid,double gamma,Clean c):m(std::move(grid)),g(gamma),clean(c) {
        auto n=m.vol.size();
        for(auto* v:{&u,&stage,&next,&du,&w,&sz,&sr}) v->assign(n,U{});
        auto nn=m.rn.size();
        for(auto* v:{&psi,&psiStage,&psiNext,&dpsi}) v->assign(nn,0.0);
        emfA.assign(static_cast<std::size_t>(m.nz+1)*m.nr,0); emfR.assign(static_cast<std::size_t>(m.nz)*(m.nr+1),0);
        b0.assign(n,{0,0}); cross.assign(n,{0,0});
        b0A.assign(static_cast<std::size_t>(m.nz+1)*m.nr,{0,0}); b0R.assign(static_cast<std::size_t>(m.nz)*(m.nr+1),{0,0});
    }
    // Cell-volume-averaged poloidal field from nodal psi (b = curl(psi/r theta)):
    //   V<b_z> = 2 pi int (psi_top - psi_bottom) dz   (psi linear along the sloped edges),
    //   V<b_r> = -2 pi (loop integral of psi dr)     (psi linear in r^2 along the vertical edges,
    //                                                    exact for psi ~ r^2 near the axis).
    // An earlier area-weighted face fit returned b_r at the face, not the cell average, which
    // left an O(1) non-converging force residual in the axis row.
    std::array<double,2> cellB(const std::vector<double>& ps,int i,int j) const {
        auto P=[&](int ii,int jj){return ps[m.node(ii,jj)];};
        double bz=2*pi*m.dz*0.5*((P(i,j+1)-P(i,j))+(P(i+1,j+1)-P(i+1,j)));
        auto colInt=[&](int ii) { double ra=m.r(ii,j),rb=m.r(ii,j+1),pa=P(ii,j),pb=P(ii,j+1),dr=rb-ra;
            return pa*dr+(pb-pa)/((rb+ra)*dr)*((rb*rb*rb-ra*ra*ra)/3-ra*ra*dr); };
        // Green: int int d(psi)/dz dz dr = loop integral of psi dr, including the sloped edges.
        double edges=0.5*(P(i,j)+P(i+1,j))*(m.r(i+1,j)-m.r(i,j))+0.5*(P(i+1,j+1)+P(i,j+1))*(m.r(i,j+1)-m.r(i+1,j+1));
        double br=-2*pi*(colInt(i+1)-colInt(i)+edges);
        double V=m.vol[m.cell(i,j)];
        return {bz/V,br/V};
    }
    void poloidalFromPsi(const std::vector<double>& ps,std::vector<U>& state) const {
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) { auto b=cellB(ps,i,j); state[m.cell(i,j)][BZ]=b[0]; state[m.cell(i,j)][BR]=b[1]; }
    }
    // Initial state from primitives; with a coil the poloidal field is taken from psi nodes.
    void init(const std::function<U(int,int)>& prim,const std::function<double(double,double)>& psiField) {
        if(psiField) for(int i=0;i<=m.nz;++i) for(int j=0;j<=m.nr;++j) psi[m.node(i,j)]=psiField(i*m.dz,m.r(i,j));
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
            U p=prim(i,j);
            if(psiField) { b0[m.cell(i,j)]=cellB(psi,i,j); if(!split()) { p[BZ]=b0[m.cell(i,j)][0]; p[BR]=b0[m.cell(i,j)][1]; } }
            u[m.cell(i,j)]=toCons(p,g);
        }
        if(split() && psiField) {
            // b = curl(psi/r theta): b_z = psi_r / r, b_r = -psi_z / r, by central differences.
            auto field=[&](double z,double r){ double h=1e-7;
                return std::array<double,2>{(psiField(z,r+h)-psiField(z,r-h))/(2*h*r),-(psiField(z+h,r)-psiField(z-h,r))/(2*h*r)}; };
            for(int i=0;i<=m.nz;++i) for(int j=0;j<m.nr;++j) b0A[static_cast<std::size_t>(i)*m.nr+j]=field(i*m.dz,0.5*(m.r(i,j)+m.r(i,j+1)));
            for(int i=0;i<m.nz;++i) for(int j=1;j<=m.nr;++j) b0R[static_cast<std::size_t>(i)*(m.nr+1)+j]=field((i+0.5)*m.dz,0.5*(m.r(i,j)+m.r(i+1,j)));
            std::fill(psi.begin(),psi.end(),0.0);  // evolved flux function is psi1
        }
        if(split() && wallB==WallB::Insulating) { vac=buildVacuum(m.nz,m.dz,[&](double z){return m.r(std::clamp(int(std::lround(z/m.dz)),0,m.nz),m.nr);},vacuumRadius,vacuumLayers); wallPsi.assign(m.nz+1,0); }
        M0=E0=P0=0;
        for(std::size_t q=0;q<u.size();++q) { M0+=u[q][RHO]*m.vol[q]; E0+=u[q][EN]*m.vol[q]; P0+=u[q][MZ]*m.vol[q]; }
    }
    U inletGhost(const U& in) const {
        if(!nozzle) return in;
        const double a0=std::sqrt(g*Rg*T0),inv=in[1]-2*std::sqrt(g*in[4]/in[RHO])/(g-1);
        double lo=0,hi=a0/std::sqrt(1+(g-1)/2);
        for(int n=0;n<40;++n) { double v=(lo+hi)/2,a=std::sqrt(a0*a0-(g-1)*v*v/2); if(v-2*a/(g-1)>inv) hi=v; else lo=v; }
        double v=(lo+hi)/2,T=T0-v*v*(g-1)/(2*g*Rg),p=p0*std::pow(T/T0,g/(g-1));
        U o=in; o[RHO]=p/(Rg*T); o[1]=v; o[2]=0; o[3]=0; o[4]=p; return o;
    }
    U outletGhost(const U& in) const {
        if(!nozzle) return in;
        double a=std::sqrt(g*in[4]/in[RHO]);
        if(in[1]>=a) return in;
        U o=in; o[4]=pBack; o[RHO]=in[RHO]*std::pow(pBack/in[4],1/g);
        o[1]=in[1]+2*(a-std::sqrt(g*pBack/o[RHO]))/(g-1); return o;
    }
    static U wallGhost(U x,double nz,double nr) {
        double un=x[1]*nz+x[2]*nr; x[1]-=2*un*nz; x[2]-=2*un*nr; return x;  // B copied: magnetically transparent wall
    }
    double stableDt() {
        double dt=1e300; ch=0;
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
            auto q=m.cell(i,j); U p; if(!toPrim(u[q],g,p)) throw std::runtime_error("bad state");
            auto bt=total(p,q); double c=std::sqrt(g*p[4]/p[RHO]+(sq(bt[0])+sq(bt[1])+sq(p[BT]))/p[RHO]);
            ch=std::max(ch,std::hypot(p[1],p[2])+c);
            double rate=(m.axialArea(i,j)+m.axialArea(i+1,j))*(std::abs(p[1])+c);
            for(int f:{j,j+1}) { auto a=m.radialArea(i,f); rate+=std::abs(p[1]*a[0]+p[2]*a[1])+c*std::hypot(a[0],a[1]); }
            dt=std::min(dt,cfl*m.vol[q]/rate);
        }
        if(clean!=Clean::GLM) ch=0;
        return dt;
    }
    Rates rhs(const std::vector<U>& s,const std::vector<double>& ps,std::vector<U>& d,std::vector<double>& dps) {
        std::fill(d.begin(),d.end(),U{});
        for(std::size_t q=0;q<s.size();++q) if(!toPrim(s[q],g,w[q])) throw std::runtime_error("bad state in rhs");
        std::fill(sz.begin(),sz.end(),U{}); std::fill(sr.begin(),sr.end(),U{});
        for(int i=1;i<m.nz-1;++i) for(int j=0;j<m.nr;++j) { auto q=m.cell(i,j);
            for(int k=0;k<NV;++k) sz[q][k]=minmod(w[q][k]-w[m.cell(i-1,j)][k],w[m.cell(i+1,j)][k]-w[q][k]); }
        for(int i=0;i<m.nz;++i) for(int j=1;j<m.nr-1;++j) { auto q=m.cell(i,j);
            for(int k=0;k<NV;++k) sr[q][k]=minmod(w[q][k]-w[m.cell(i,j-1)][k],w[m.cell(i,j+1)][k]-w[q][k]); }
        auto oneSided=[&](std::size_t q,std::size_t a,std::size_t b,double sign,U& sl) {
            for(int k=0;k<NV;++k) sl[k]=sign*minmod(w[q][k]-w[a][k],w[a][k]-w[b][k]);
            for(int k:{int(RHO),4}) sl[k]=std::clamp(sl[k],-w[q][k],w[q][k]);
        };
        for(int j=0;j<m.nr;++j) {
            oneSided(m.cell(0,j),m.cell(1,j),m.cell(2,j),-1,sz[m.cell(0,j)]);
            oneSided(m.cell(m.nz-1,j),m.cell(m.nz-2,j),m.cell(m.nz-3,j),1,sz[m.cell(m.nz-1,j)]);
        }
        // Axis row: mirror cell across r=0 with parity (odd: u_r, u_theta, b_r, b_theta). A zero
        // slope here makes odd quantities O(1) wrong at the first radial face (b_r ~ r), which
        // measured as a non-converging force residual of the magnetostatic equilibrium.
        if(m.inner<=0 && m.nr>=2 && axisParity) for(int i=0;i<m.nz;++i) { auto q=m.cell(i,0),n=m.cell(i,1);
            for(int k:{2,3,int(BR),int(BT)}) sr[q][k]=minmod(2*w[q][k],w[n][k]-w[q][k]); }
        if(m.nr>=3) for(int i=0;i<m.nz;++i) {
            oneSided(m.cell(i,m.nr-1),m.cell(i,m.nr-2),m.cell(i,m.nr-3),1,sr[m.cell(i,m.nr-1)]);
            if(m.inner>0) oneSided(m.cell(i,0),m.cell(i,1),m.cell(i,2),-1,sr[m.cell(i,0)]);
        }
        if(firstOrder) { std::fill(sz.begin(),sz.end(),U{}); std::fill(sr.begin(),sr.end(),U{}); }
        auto rec=[&](std::size_t q,bool axial,double dir) { U x=w[q]; const U& sl=axial?sz[q]:sr[q]; for(int k=0;k<NV;++k) x[k]+=dir*0.5*sl[k]; return x; };
        Rates rt{};
        std::fill(cross.begin(),cross.end(),std::array<double,2>{0,0});
        // Split: add the face b0 to both states, then remove T(b0) and the cross stress from the
        // momentum flux and b0 . (induction flux) from the energy flux. Returns the cross stress . n.
        auto faceFlux=[&](U l,U r,double nz,double nr,std::array<double,2> f0,std::array<double,2>& crossN) {
            crossN={0,0};
            if(!split()) return mhdFlux(l,r,nz,nr,ch,g);
            double b1z=0.5*(l[BZ]+r[BZ]),b1r=0.5*(l[BR]+r[BR]);
            for(U* x:{&l,&r}) { (*x)[BZ]+=f0[0]; (*x)[BR]+=f0[1]; }
            auto ff=mhdFlux(l,r,nz,nr,ch,g);
            double b0n=f0[0]*nz+f0[1]*nr,b1n=b1z*nz+b1r*nr,dot=f0[0]*b1z+f0[1]*b1r,b02=sq(f0[0])+sq(f0[1]);
            std::array<double,2> T0{f0[0]*b0n-0.5*b02*nz,f0[1]*b0n-0.5*b02*nr};
            crossN={dot*nz-f0[0]*b1n-b1z*b0n,dot*nr-f0[1]*b1n-b1r*b0n};
            ff.f[MZ]+=T0[0]-crossN[0]; ff.f[MR]+=T0[1]-crossN[1];
            ff.f[EN]-=(-f0[0]*nr+f0[1]*nz)*ff.emf;
            return ff;
        };
        for(int i=0;i<=m.nz;++i) for(int j=0;j<m.nr;++j) {
            auto il=m.cell(std::max(0,i-1),j),ir=m.cell(std::min(m.nz-1,i),j);
            U l=rec(il,true,1),r=rec(ir,true,-1);
            if(i==0) { l=inletGhost(r); if(split() && wallB!=WallB::Transparent) { l[BZ]=0; l[BR]=0; } }  // inflow carries only b0
            if(i==m.nz) r=outletGhost(l);
            std::array<double,2> cn; auto ff=faceFlux(l,r,1,0,b0A[static_cast<std::size_t>(i)*m.nr+j],cn); double A=m.axialArea(i,j);
            for(int k=0;k<NV;++k) { if(i>0) d[il][k]-=A*ff.f[k]; if(i<m.nz) d[ir][k]+=A*ff.f[k]; }
            for(int k=0;k<2;++k) { if(i>0) cross[il][k]-=A*cn[k]; if(i<m.nz) cross[ir][k]+=A*cn[k]; }
            emfA[static_cast<std::size_t>(i)*m.nr+j]=ff.emf;
            if(i==0) { rt.mass+=A*ff.f[RHO]; rt.energy+=A*ff.f[EN]; rt.inMom+=A*ff.f[MZ]; rt.inMass+=A*ff.f[RHO]; }
            if(i==m.nz) { rt.mass-=A*ff.f[RHO]; rt.energy-=A*ff.f[EN]; rt.outMom+=A*ff.f[MZ]; rt.outMass+=A*ff.f[RHO]; }
        }
        if(split() && wallB==WallB::Insulating) for(int i=0;i<=m.nz;++i) wallPsi[i]=ps[m.node(i,m.nr)];
        for(int i=0;i<m.nz;++i) for(int j=0;j<=m.nr;++j) {
            if(j==0 && m.inner<=0) continue;  // axis: zero-area face
            auto a=m.radialArea(i,j); double A=std::hypot(a[0],a[1]),nz=a[0]/A,nr=a[1]/A;
            U l,r;
            if(j==0) { r=rec(m.cell(i,0),false,-1); l=wallGhost(r,nz,nr); }
            else if(j==m.nr) { l=rec(m.cell(i,m.nr-1),false,1); r=wallGhost(l,nz,nr);
                if(split() && wallB==WallB::Conducting) { double bn=r[BZ]*nz+r[BR]*nr; r[BZ]-=2*bn*nz; r[BR]-=2*bn*nr; } }
            else { l=rec(m.cell(i,j-1),false,1); r=rec(m.cell(i,j),false,-1); }
            std::array<double,2> cn; auto ff=faceFlux(l,r,nz,nr,b0R[static_cast<std::size_t>(i)*(m.nr+1)+j],cn);
            if(j==m.nr && split() && wallB==WallB::Insulating) {
                // The Riemann solve above sees no b1 jump (a jump there is dissipated every step as
                // numerical Joule heating). The exterior field enters through the face stress and the
                // Poynting flux instead: a b1_t jump at the wall is a plasma-side current sheet, whose
                // force acts on the wall cell. b1_n is the interior value (continuous).
                auto e=vac.field(i,wallPsi); auto f0=b0R[static_cast<std::size_t>(i)*(m.nr+1)+j];
                double bin=l[BZ]*nz+l[BR]*nr,bit=-l[BZ]*nr+l[BR]*nz,bet=-e[0]*nr+e[1]*nz;
                std::array<double,2> bi{l[BZ],l[BR]},be{bin*nz-bet*nr,bin*nr+bet*nz};
                auto M=[&](std::array<double,2> b) { double bb=0.5*(sq(b[0])+sq(b[1])),bn=b[0]*nz+b[1]*nr; return std::array<double,2>{bb*nz-b[0]*bn,bb*nr-b[1]*bn}; };
                auto C=[&](std::array<double,2> b) { double dot=f0[0]*b[0]+f0[1]*b[1],b0n=f0[0]*nz+f0[1]*nr,bn=b[0]*nz+b[1]*nr;
                    return std::array<double,2>{dot*nz-f0[0]*bn-b[0]*b0n,dot*nr-f0[1]*bn-b[1]*b0n}; };
                auto Mi=M(bi),Me=M(be),Ce=C(be);
                for(int k=0;k<2;++k) { ff.f[MZ+k]+=Me[k]-Mi[k]+cn[k]-Ce[k]; cn[k]=Ce[k]; }
                ff.f[EN]+=ff.emf*(bet-bit);
            }
            for(int k=0;k<2;++k) { if(j>0) cross[m.cell(i,j-1)][k]-=A*cn[k]; if(j<m.nr) cross[m.cell(i,j)][k]+=A*cn[k]; }
            bool wall=j==0||j==m.nr;
            if(wall) { ff.f[RHO]=0; ff.f[ENT]=0; }  // impermeable
            if(j>0) for(int k=0;k<NV;++k) d[m.cell(i,j-1)][k]-=A*ff.f[k];
            if(j<m.nr) for(int k=0;k<NV;++k) d[m.cell(i,j)][k]+=A*ff.f[k];
            emfR[static_cast<std::size_t>(i)*(m.nr+1)+j]=ff.emf;
            if(j==m.nr) { rt.wallMom-=A*ff.f[MZ]; rt.energy-=A*ff.f[EN]; }
            if(j==0)    { rt.wallMom+=A*ff.f[MZ]; rt.energy+=A*ff.f[EN]; }
        }
        // Axisymmetric geometric sources, all with the exact meridional measure (well-balanced for uniform p + b^2/2).
        for(std::size_t q=0;q<s.size();++q) {
            const U& x=w[q]; double M=m.meas[q];
            d[q][MR]+=(x[RHO]*x[3]*x[3]+x[4]+0.5*(sq(x[BZ])+sq(x[BR])+sq(x[BT]))-x[BT]*x[BT])*M;
            d[q][MT]-=(x[RHO]*x[2]*x[3]-x[BR]*x[BT])*M;
            d[q][BT]+=(x[2]*x[BT]-x[3]*x[BR])*M;
            d[q][BR]+=x[PS]*M;
            if(split()) {  // cross body force: -div(C) with C_theta_theta = b0 . b1
                double fz=cross[q][0],fr=cross[q][1]+(b0[q][0]*x[BZ]+b0[q][1]*x[BR])*M;
                d[q][MZ]+=fz; d[q][MR]+=fr; d[q][EN]+=0; rt.body+=fz;
            }
            for(double& v:d[q]) v/=m.vol[q];
        }
        std::fill(dps.begin(),dps.end(),0.0);
        if(ct()) for(int i=0;i<=m.nz;++i) for(int j=0;j<=m.nr;++j) {
            double r0=m.r(i,j); if(r0<=0) continue;  // psi = 0 on the axis
            if(j==m.nr && split() && wallB==WallB::Conducting) continue;  // flux through the wall frozen
            double e=0; int n=0;
            if(j>0)    { e+=emfA[static_cast<std::size_t>(i)*m.nr+j-1]; ++n; }
            if(j<m.nr) { e+=emfA[static_cast<std::size_t>(i)*m.nr+j]; ++n; }
            if(i>0)    { e+=emfR[static_cast<std::size_t>(i-1)*(m.nr+1)+j]; ++n; }
            if(i<m.nz) { e+=emfR[static_cast<std::size_t>(i)*(m.nr+1)+j]; ++n; }
            dps[m.node(i,j)]=r0*e/n;
        }
        return rt;
    }
    bool admissible(const std::vector<U>& s) { U p; for(const auto& x:s) if(!toPrim(x,g,p)) return false; return true; }
    double step(double maxDt) {
        double dt=std::min(stableDt(),maxDt);
        firstOrder=forceHll=false;
        for(int attempt=0;attempt<14;++attempt) {
            try {
                auto r1=rhs(u,psi,du,dpsi);
                for(std::size_t q=0;q<u.size();++q) for(int k=0;k<NV;++k) stage[q][k]=u[q][k]+dt*du[q][k];
                for(std::size_t n=0;n<psi.size();++n) psiStage[n]=psi[n]+dt*dpsi[n];
                if(ct()) poloidalFromPsi(psiStage,stage);
                if(admissible(stage)) {
                    auto r2=rhs(stage,psiStage,du,dpsi);
                    for(std::size_t q=0;q<u.size();++q) for(int k=0;k<NV;++k) next[q][k]=0.5*(u[q][k]+stage[q][k]+dt*du[q][k]);
                    for(std::size_t n=0;n<psi.size();++n) psiNext[n]=0.5*(psi[n]+psiStage[n]+dt*dpsi[n]);
                    if(ct()) poloidalFromPsi(psiNext,next);
                    if(admissible(next)) {
                        if(clean==Clean::GLM) for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
                            double h=std::min(m.dz,(m.r(i,j+1)-m.r(i,j)+m.r(i+1,j+1)-m.r(i+1,j))/2);
                            next[m.cell(i,j)][PS]*=std::exp(-glmAlpha*ch*dt/h);
                        }
                        if(robust) synchronise(next);
                        if(firstOrder) ++robustSteps;
                        firstOrder=forceHll=false;
                        u.swap(next); psi.swap(psiNext); t+=dt; ++steps;
                        last={0.5*(r1.mass+r2.mass),0.5*(r1.energy+r2.energy),0.5*(r1.inMom+r2.inMom),0.5*(r1.outMom+r2.outMom),
                              0.5*(r1.wallMom+r2.wallMom),0.5*(r1.inMass+r2.inMass),0.5*(r1.outMass+r2.outMass),0.5*(r1.body+r2.body)};
                        intM+=dt*last.mass; intE+=dt*last.energy;
                        intP+=dt*(last.inMom-last.outMom+last.wallMom+last.body);
                        intPg+=dt*(std::abs(last.inMom)+std::abs(last.outMom)+std::abs(last.wallMom)+std::abs(last.body));
                        return dt;
                    }
                }
            } catch(const std::runtime_error&) {}
            ++rejected;
            if(robust && !firstOrder) { firstOrder=forceHll=true; continue; }  // same dt, first order + HLL
            dt*=0.5;
        }
        firstOrder=forceHll=false;
        throw std::runtime_error("could not advance an admissible state");
    }
    // Dual-energy synchronisation: where the total-energy pressure is reliable, reset the entropy to
    // it; elsewhere reset the total energy to the entropy pressure and book the change (not hidden).
    void synchronise(std::vector<U>& s) {
        for(std::size_t q=0;q<s.size();++q) { U& x=s[q];
            double ke=0.5*(sq(x[MZ])+sq(x[MR])+sq(x[MT]))/x[RHO],mag=0.5*(sq(x[BZ])+sq(x[BR])+sq(x[BT]));
            double pE=(g-1)*(x[EN]-ke-mag);
            if(pE>dualEta*(g-1)*x[EN]) x[ENT]=pE/std::pow(x[RHO],g-1);
            else { double e=x[ENT]*std::pow(x[RHO],g-1)/(g-1)+ke+mag; intSync+=(e-x[EN])*m.vol[q]; x[EN]=e; }
        }
    }
    void advanceTo(double T) { while(t<T*(1-1e-14)) step(T-t); }
    // Discrete J_theta from the circulation of face-averaged b around the meridional cell loop.
    double jTheta(int i,int j) const {
        auto bAt=[&](int ii,int jj) { auto q=m.cell(std::clamp(ii,0,m.nz-1),std::clamp(jj,0,m.nr-1)); return std::array<double,2>{u[q][BZ],u[q][BR]}; };
        auto c=bAt(i,j);
        auto face=[&](int ii,int jj) { auto o=bAt(ii,jj); bool out=ii<0||ii>=m.nz||jj<0||jj>=m.nr; return out?c:std::array<double,2>{0.5*(c[0]+o[0]),0.5*(c[1]+o[1])}; };
        auto seg=[&](std::array<double,2> b,double z0,double r0,double z1,double r1) { return b[0]*(z1-z0)+b[1]*(r1-r0); };
        double z0=i*m.dz,z1=(i+1)*m.dz;
        double circ=seg(face(i,j-1),z0,m.r(i,j),z1,m.r(i+1,j))+seg(face(i+1,j),z1,m.r(i+1,j),z1,m.r(i+1,j+1))
                   +seg(face(i,j+1),z1,m.r(i+1,j+1),z0,m.r(i,j+1))+seg(face(i-1,j),z0,m.r(i,j+1),z0,m.r(i,j));
        return circ/(m.meas[m.cell(i,j)]/(2*pi));
    }
    Report report() const {
        Report R{}; R.t=t; R.steps=steps; R.rejected=rejected; R.minBeta=1e300;
        double bmax=1e-300,exitArea=0;
        for(const auto& x:u) bmax=std::max(bmax,std::sqrt(sq(x[BZ])+sq(x[BR])+sq(x[BT])));
        double vsum=0;
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
            auto q=m.cell(i,j); U p; toPrim(u[q],g,p);
            R.mass+=u[q][RHO]*m.vol[q]; R.energy+=u[q][EN]*m.vol[q]; R.Pz+=u[q][MZ]*m.vol[q];
            auto bt=total(p,q); double a=std::sqrt(g*p[4]/p[RHO]),b2=sq(bt[0])+sq(bt[1])+sq(p[BT]);
            R.maxMach=std::max(R.maxMach,std::hypot(p[1],p[2])/a); R.maxSpeed=std::max(R.maxSpeed,std::sqrt(sq(p[1])+sq(p[2])+sq(p[3])));
            if(b2>0) R.minBeta=std::min(R.minBeta,2*p[4]/b2);
            // Cell-centred divergence from face-averaged b (the discrete constraint of none/glm).
            double div=0;
            auto nb=[&](int ii,int jj,const std::array<double,2>& a,double sign) {
                bool out=ii<0||ii>=m.nz||jj<0||jj>=m.nr; auto o=out?q:m.cell(ii,jj);
                div+=sign*(0.5*(u[q][BZ]+u[o][BZ])*a[0]+0.5*(u[q][BR]+u[o][BR])*a[1]);
            };
            nb(i-1,j,{m.axialArea(i,j),0},-1); nb(i+1,j,{m.axialArea(i+1,j),0},1);
            nb(i,j-1,m.radialArea(i,j),-1); nb(i,j+1,m.radialArea(i,j+1),1);
            double h=std::min(m.dz,(m.r(i,j+1)-m.r(i,j)+m.r(i+1,j+1)-m.r(i+1,j))/2);
            double e=std::abs(div)/m.vol[q]*h/bmax;
            R.divMax=std::max(R.divMax,e); R.divMean+=e*m.vol[q]; vsum+=m.vol[q];
            if(ct()) {
                double f=-2*pi*(psi[m.node(i,j+1)]-psi[m.node(i,j)])+2*pi*(psi[m.node(i+1,j+1)]-psi[m.node(i+1,j)])
                        -(-2*pi*(psi[m.node(i+1,j)]-psi[m.node(i,j)]))+(-2*pi*(psi[m.node(i+1,j+1)]-psi[m.node(i,j+1)]));
                R.ctFaceDivMax=std::max(R.ctFaceDivMax,std::abs(f)/m.vol[q]*h/bmax);
            }
            double jt=jTheta(i,j);
            R.lorentzZ+=-jt*bt[1]*m.vol[q];            // (J x b)_z with J_r b_theta = 0 here
            R.coilLorentzZ+=-jt*b0[q][1]*m.vol[q];     // part exerted by the applied (coil) field
            if(i==m.nz-1) exitArea+=m.axialArea(m.nz,j);
        }
        R.divMean/=vsum;
        R.inMom=last.inMom; R.outMom=last.outMom; R.wallMom=last.wallMom; R.inMdot=last.inMass; R.outMdot=last.outMass;
        R.ambient=pBack*exitArea;
        R.massErr=(R.mass-M0-intM)/M0; R.energyErr=(R.energy-E0-intE-intSync)/E0; R.sync=intSync; R.robustSteps=robustSteps;
        R.momErr=(R.Pz-P0-intP)/std::max(std::abs(P0)+intPg,1e-300);
        R.body=last.body; R.deviceThrust=R.inMom+R.wallMom+R.body-R.ambient; R.exitPlaneThrust=R.outMom-R.ambient;
        return R;
    }
};

// ---------------------------------------------------------------- cases
double seconds(std::chrono::steady_clock::time_point s) { return std::chrono::duration<double>(std::chrono::steady_clock::now()-s).count(); }

// MHD shock tube (Brio & Wu 1988, gamma=2) in a thin annulus far from the axis (inner radius
// 1e6, so curvature terms are ~1e-6). Axial: shock normal z, transverse field b_theta, run in
// the 4 radial rows. Radial: shock normal r (b_r = 0.75), transverse field b_z, run across the
// annulus in 4 axial columns; this exercises the meridional face rotation and the r-weighted
// radial fluxes and sources. (A b_r-transverse tube along z is not a valid 1D test: its
// transverse velocity u_r is blocked by the annulus walls.)
void shockTube(const std::vector<int>& levels) {
    const double g=2,inner=1e6;
    auto run=[&](int n,bool radial,Clean c,std::vector<U>& out,int across=4) {
        int nz=radial?across:n,nr=radial?n:across;
        Grid grid(nz,nr,radial?0.01:1.0,[&](double){return inner+(radial?1.0:0.01);},inner);
        Mhd s(grid,g,c);
        s.init([&](int i,int j){ U p{}; bool left=(radial?j:i)<n/2; p[RHO]=left?1:0.125; p[4]=left?1:0.1;
            if(radial) { p[BR]=0.75; p[BZ]=left?-1:1; } else { p[BZ]=0.75; p[BT]=left?1:-1; } return p; },nullptr);
        try { s.advanceTo(0.1); } catch(const std::exception&) {
            std::fprintf(stderr,"shock tube %s n=%d failed at t=%.5f after %llu steps\n",radial?"radial":"axial",n,s.t,(unsigned long long)s.steps); throw; }
        // Average across the transverse cells; report (rho, normal u, transverse b) in the tube frame.
        out.assign(n,U{}); int m=radial?nz:nr;
        for(int k=0;k<n;++k) for(int l=0;l<m;++l) { U p; toPrim(s.u[radial?grid.cell(l,k):grid.cell(k,l)],g,p);
            out[k][RHO]+=p[RHO]/m; out[k][1]+=(radial?p[2]:p[1])/m; out[k][BT]+=(radial?-p[BZ]:p[BT])/m; out[k][4]+=p[4]/m; out[k][3]+=(radial?-p[1]:p[3])/m; }
        return s.report();
    };
    int refN=12800; std::vector<U> ref; auto rr=run(refN,false,Clean::None,ref,1);
    std::fflush(stdout);
    std::printf("# Brio-Wu reference: axial tube, nz=%d, steps=%llu, hll fallbacks=%ld\n",refN,(unsigned long long)rr.steps,hllFallbacks);
    std::printf("case,direction,n,L1_rho,L1_un,L1_btransverse,mass_res,energy_res,steps,hll_fallbacks\n");
    for(int n:levels) for(bool radial:{false,true}) {
        hllFallbacks=0; std::vector<U> sol; auto R=run(n,radial,Clean::None,sol);
        double e[3]={0,0,0}; int ratio=refN/n;
        for(int i=0;i<n;++i) { double a[3]={0,0,0};
            for(int k=0;k<ratio;++k){const U& x=ref[i*ratio+k]; a[0]+=x[RHO]/ratio; a[1]+=x[1]/ratio; a[2]+=x[BT]/ratio;}
            e[0]+=std::abs(sol[i][RHO]-a[0])/n; e[1]+=std::abs(sol[i][1]-a[1])/n; e[2]+=std::abs(sol[i][BT]-a[2])/n; }
        std::printf("briowu,%s,%d,%.4e,%.4e,%.4e,%.1e,%.1e,%llu,%ld\n",radial?"radial":"axial",n,e[0],e[1],e[2],R.massErr,R.energyErr,(unsigned long long)R.steps,hllFallbacks);
        std::fflush(stdout);
    }
    for(bool radial:{false,true}) {
        std::vector<U> sol; run(800,radial,Clean::None,sol);
        std::FILE* f=std::fopen(radial?"/tmp/crucible_mhd/briowu_radial_800.csv":"/tmp/crucible_mhd/briowu_axial_800.csv","w");
        if(f) { std::fprintf(f,"x,rho,un,ut,p,bt,ref_rho,ref_un,ref_bt\n");
            for(int i=0;i<800;++i){ const U& a=sol[i]; U r{}; for(int k=0;k<16;++k) for(int v=0;v<NV;++v) r[v]+=ref[i*16+k][v]/16;
                std::fprintf(f,"%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",(i+0.5)/800.,a[RHO],a[1],a[3],a[4],a[BT],r[RHO],r[1],r[BT]); } std::fclose(f); }
    }
}

// Field-aligned check: Sod along z in a full cylinder (with axis) threaded by a strong uniform
// axial field. The field exerts no force on flow along it, so the exact answer is gas Sod.
void aligned(const std::vector<int>& levels) {
    const double g=1.4;
    std::printf("case,b_z,beta_right,nz,nr,L1_rho_vs_exact,max_abs_ur,max_abs_br,max_abs_bt,div_max,mass_res,energy_res,hll_fallbacks\n");
    for(double b:{0.0,4.47213595,14.1421356}) for(int nz:levels) for(Clean c:{Clean::None,Clean::CT}) {
        if(b==0 && c==Clean::CT) continue;
        hllFallbacks=0;
        Grid grid(nz,8,1.0,[](double){return 0.1;});
        Mhd s(grid,g,c);
        s.init([&](int i,int){ U p{}; bool left=i<nz/2; p[RHO]=left?1:0.125; p[4]=left?1:0.1; p[BZ]=b; return p; },
               c==Clean::CT?std::function<double(double,double)>([&](double,double r){return b*r*r/2;}):nullptr);
        s.advanceTo(0.2);
        double err=0,ur=0,br=0,bt=0;
        for(int i=0;i<nz;++i) for(int j=0;j<8;++j) { U p; toPrim(s.u[grid.cell(i,j)],g,p);
            double x=(i+0.5)/nz-0.5; double ex=0; for(int k=0;k<8;++k) ex+=exactDensity({1,0,1},{0.125,0,0.1},g,(x+(k-3.5)/(8.0*nz))/0.2)/8;
            err+=std::abs(p[RHO]-ex)*grid.vol[grid.cell(i,j)]; ur=std::max(ur,std::abs(p[2])); br=std::max(br,std::abs(p[BR])); bt=std::max(bt,std::abs(p[BT])); }
        double V=0; for(double v:grid.vol) V+=v;
        auto R=s.report();
        std::printf("aligned-%s,%.3f,%.3g,%d,8,%.4e,%.1e,%.1e,%.1e,%.1e,%.1e,%.1e,%ld\n",name(c),b,b>0?2*0.1/(b*b):0.0,nz,err/V,ur,br,bt,R.divMax,R.massErr,R.energyErr,hllFallbacks);
    }
}

double nozzleRadius(double z) {
    const double L=0.6,ri=0.035,rt=0.020,re=0.035,tf=0.36; double x=z/L;
    if(x<tf) { double t=x/tf; return rt+(ri-rt)*(1+std::cos(pi*t))/2; }
    double t=(x-tf)/(1-tf); return rt+(re-rt)*(1-std::cos(pi*t))/2;
}
// Coil: one loop, radius 0.05 m, in the throat plane; current set by the field at its centre.
constexpr double coilA=0.05,coilZ=0.216;
double coilCurrent(double tesla) { return 2*coilA*tesla/mu0; }

void dumpField(const Mhd& s,const std::string& path) {
    const Grid& grid=s.m; std::FILE* f=std::fopen(path.c_str(),"w"); if(!f) return;
    std::fprintf(f,"i,j,z0,z1,rb0,rb1,rt0,rt1,rho,uz,ur,ut,p,bz,br,bt,mach,b0z,b0r,jtheta\n");
    for(int i=0;i<grid.nz;++i) for(int j=0;j<grid.nr;++j) { U p; toPrim(s.u[grid.cell(i,j)],s.g,p);
        std::fprintf(f,"%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n",i,j,i*grid.dz,(i+1)*grid.dz,
            grid.r(i,j),grid.r(i+1,j),grid.r(i,j+1),grid.r(i+1,j+1),p[RHO],p[1],p[2],p[3],p[4],p[BZ]*std::sqrt(mu0),p[BR]*std::sqrt(mu0),p[BT]*std::sqrt(mu0),
            std::hypot(p[1],p[2])/std::sqrt(s.g*p[4]/p[RHO]),s.b0[grid.cell(i,j)][0]*std::sqrt(mu0),s.b0[grid.cell(i,j)][1]*std::sqrt(mu0),s.jTheta(i,j)*std::sqrt(mu0)/mu0); }
    std::fclose(f);
    // Nodal flux function (T m^2) for field lines; for none/glm it is the initial (applied) psi.
    if(std::FILE* n=std::fopen((path+".psi").c_str(),"w")) { std::fprintf(n,"i,j,z,r,psi\n");
        for(int i=0;i<=grid.nz;++i) for(int j=0;j<=grid.nr;++j) std::fprintf(n,"%d,%d,%.6e,%.6e,%.8e\n",i,j,i*grid.dz,grid.r(i,j),s.psi[grid.node(i,j)]*std::sqrt(mu0));
        std::fclose(n); }
    // Insulating wall: interior wall-cell b1 next to the matched exterior (ghost) b1, in tesla.
    if(s.split() && s.wallB==WallB::Insulating) if(std::FILE* n=std::fopen((path+".wall").c_str(),"w")) {
        std::vector<double> wp(grid.nz+1); for(int i=0;i<=grid.nz;++i) wp[i]=s.psi[grid.node(i,grid.nr)];
        std::fprintf(n,"i,z,b1z_cell,b1r_cell,b1z_ext,b1r_ext,psi1_wall\n");
        for(int i=0;i<grid.nz;++i) { auto b=s.vac.field(i,wp); const U& x=s.u[grid.cell(i,grid.nr-1)];
            std::fprintf(n,"%d,%.5f,%.5f,%.5f,%.5f,%.5f,%.5e\n",i,(i+0.5)*grid.dz,x[BZ]*std::sqrt(mu0),x[BR]*std::sqrt(mu0),b[0]*std::sqrt(mu0),b[1]*std::sqrt(mu0),wp[i]*std::sqrt(mu0)); }
        std::fclose(n); }
}
// Magnetostatic equilibrium: gas at rest, uniform pressure, vacuum coil field (J = 0).
// Exact answer: nothing moves. Measures the discrete force imbalance of a curl-free field.
void staticCoil(const std::vector<int>& levels,double tesla) {
    const double g=1.4,p=1e5,rho=p/(287.05*300);
    std::printf("case,B_centre_T,nz,nr,clean,t_ms,max_speed_m_s,max_speed_over_vA,div_max,div_mean,ct_face_div,lorentz_z_N,steps\n");
    double I=coilCurrent(tesla),vA=tesla/std::sqrt(mu0*rho);
    for(int nz:levels) for(Clean c:{Clean::None,Clean::GLM,Clean::CT,Clean::Split}) {
        int nr=nz*3/20; Grid grid(nz,nr,0.6,nozzleRadius);
        Mhd s(grid,g,c);
        s.init([&](int,int){ U w{}; w[RHO]=rho; w[4]=p; return w; },[&](double z,double r){return loopPsi(r,z,coilA,coilZ,I);});
        try { s.advanceTo(1e-3); } catch(const std::exception& e) { std::printf("static,%.2f,%d,%d,%s,FAILED at t=%.3e ms after %llu steps\n",tesla,nz,nr,name(c),s.t*1e3,(unsigned long long)s.steps); }
        dumpField(s,"/tmp/crucible_mhd/static_"+std::string(name(c))+"_"+std::to_string(nz)+".csv");
        auto R=s.report();
        std::printf("static,%.2f,%d,%d,%s,1,%.4e,%.3e,%.2e,%.2e,%.1e,%.4e,%llu\n",tesla,nz,nr,name(c),R.maxSpeed,R.maxSpeed/vA,R.divMax,R.divMean,R.ctFaceDivMax,R.lorentzZ,(unsigned long long)R.steps);
    }
}

// Initial force residual of the magnetostatic equilibrium: acceleration du/dt at t=0.
void residual(const std::vector<int>& levels,double tesla) {
    const double g=1.4,p=1e5,rho=p/(287.05*300); double I=coilCurrent(tesla);
    std::printf("case,nz,nr,max_accel_m_s2,at_i,at_j,rms_accel,max_dEdt_W_m3,accel_scale_b2_over_rhoL\n");
    for(bool parity:{false,true}) for(int nz:levels) { int nr=nz*3/20; Grid grid(nz,nr,0.6,nozzleRadius); Mhd s(grid,g,Clean::None); s.axisParity=parity;
        s.init([&](int,int){ U w{}; w[RHO]=rho; w[4]=p; return w; },[&](double z,double r){return loopPsi(r,z,coilA,coilZ,I);});
        s.rhs(s.u,s.psi,s.du,s.dpsi);
        double mx=0,rms=0,de=0; int ai=0,aj=0;
        for(int i=0;i<nz;++i) for(int j=0;j<nr;++j) { auto q=grid.cell(i,j); double a=std::hypot(s.du[q][MZ],s.du[q][MR])/rho;
            rms+=a*a*grid.vol[q]; if(a>mx){mx=a;ai=i;aj=j;} de=std::max(de,std::abs(s.du[q][EN])); }
        double V=0; for(double v:grid.vol) V+=v;
        std::printf("residual-%s,%d,%d,%.4e,%d,%d,%.4e,%.2e,%.3e\n",parity?"parity":"zero-slope",nz,nr,mx,ai,aj,std::sqrt(rms/V),de,sq(tesla)/mu0/rho/0.05);
    }
}
// Diagnostic only: z-momentum budget of one axis cell, split by face, for the static coil state.
void axisProbe(int nz,double zProbe) {
    const double g=1.4,p=1e5,rho=p/(287.05*300); double I=coilCurrent(1.0); int nr=nz*3/20;
    Grid grid(nz,nr,0.6,nozzleRadius); Mhd s(grid,g,Clean::None);
    s.init([&](int,int){ U w{}; w[RHO]=rho; w[4]=p; return w; },[&](double z,double r){return loopPsi(r,z,coilA,coilZ,I);});
    s.rhs(s.u,s.psi,s.du,s.dpsi);
    int i=int(zProbe/grid.dz);
    for(int j=0;j<3;++j) {
        auto q=grid.cell(i,j); const U& c=s.u[q];
        // Face fluxes of z-momentum with exact field values at face centres for comparison.
        double zl=i*grid.dz,zr=(i+1)*grid.dz,rm=[&]{return 0.25*(grid.r(i,j)+grid.r(i+1,j)+grid.r(i,j+1)+grid.r(i+1,j+1));}();
        auto exact=[&](double r,double z){ double h=1e-7; double bz=r>1e-9?(loopPsi(r+h,z,coilA,coilZ,I)-loopPsi(r-h,z,coilA,coilZ,I))/(2*h*r):2*loopPsi(1e-4,z,coilA,coilZ,I)/1e-8;
            double br=r>1e-9?-(loopPsi(r,z+h,coilA,coilZ,I)-loopPsi(r,z-h,coilA,coilZ,I))/(2*h*r):0; return std::array<double,2>{bz,br}; };
        auto e=exact(rm,0.5*(zl+zr));
        std::printf("cell j=%d: cell bz=%.6g br=%.6g | exact at centre bz=%.6g br=%.6g | accel_z=%.4e accel_r=%.4e\n",j,c[BZ],c[BR],e[0],e[1],s.du[q][MZ]/rho,s.du[q][MR]/rho);
        // Split: stress through the axial faces vs the radial faces, using physical flux of cell-face values.
        double top=grid.r(i,j+1),topR=grid.r(i+1,j+1);
        auto et=exact(0.5*(top+topR),0.5*(zl+zr));
        std::printf("   exact top-face (r=%.5f) bz=%.6g br=%.6g ; T_rz=-bz*br=%.6g\n",0.5*(top+topR),et[0],et[1],-et[0]*et[1]);
    }
}
// Default nozzle (core Definition) with prepared 1D flow and a coil at the throat.
void nozzleCoil(const std::vector<int>& levels,const std::vector<double>& teslas,double tEnd,bool dump,std::vector<Clean> cleans={Clean::None,Clean::GLM,Clean::CT,Clean::Split},
                std::vector<WallB> walls={WallB::Transparent},bool robust=false) {
    dualEnergy=robust;
    crucible::Definition d; const double g=d.gas.gamma;
    std::printf("case,B_centre_T,nz,nr,clean,t_ms,device_thrust_N,exit_plane_thrust_N,core_gas_thrust_N,inlet_mom_N,wall_force_on_gas_N,mdot_out,"
                "mass_res,energy_res,mom_res,div_max,div_mean,ct_face_div,lorentz_z_N,coil_lorentz_z_N,body_axial_N,min_beta,steps,rejected,hll_fallbacks,wall_s,cell_updates_per_s,"
                "wall_b,robust,sync_J,robust_steps,wall_flux_mWb_z044\n");
    for(int nz:levels) {
        d.nz=nz; d.nr=nz*3/20;
        crucible::Flow core(d); core.advanceTo(tEnd); double coreThrust=core.measurements().deviceThrust;
        for(double tesla:teslas) for(Clean c:cleans) for(WallB wb:walls) {
            if(tesla==0 && c!=Clean::None) continue;
            if(c!=Clean::Split && wb!=walls.front()) continue;
            hllFallbacks=0;
            Grid grid(d.nz,d.nr,d.length,nozzleRadius);
            Mhd s(grid,g,c); s.wallB=wb; s.robust=robust; s.nozzle=true; s.p0=d.totalPressure; s.T0=d.totalTemperature; s.Rg=d.gas.specificR; s.pBack=d.backPressure;
            double I=coilCurrent(tesla);
            s.init([&](int i,int j){
                double rad=(grid.r(i,grid.nr)+grid.r(i+1,grid.nr))/2;
                double M=crucible::machFromArea(sq(rad/d.throatRadius),(i+0.5)*grid.dz>d.length*d.throatFraction,g);
                double f=1+(g-1)/2*M*M,T=d.totalTemperature/f; U w{};
                w[4]=d.totalPressure/std::pow(f,g/(g-1)); w[RHO]=w[4]/(d.gas.specificR*T); w[1]=M*std::sqrt(g*d.gas.specificR*T);
                double rc=(grid.r(i,j)+grid.r(i,j+1)+grid.r(i+1,j)+grid.r(i+1,j+1))/4;
                w[2]=w[1]*rc/rad*(grid.r(i+1,grid.nr)-grid.r(i,grid.nr))/grid.dz; return w; },
                tesla>0?std::function<double(double,double)>([&](double z,double r){return loopPsi(r,z,coilA,coilZ,I);}):nullptr);
            auto start=std::chrono::steady_clock::now(); std::string failure;
            try { s.advanceTo(tEnd); } catch(const std::exception& e) { failure=e.what(); }
            double wall=seconds(start); auto R=s.report();
            // Induced flux through the wall circle at z = 0.44 m (where the 1 T transparent run failed).
            double wflux=2*pi*s.psi[grid.node(std::min(grid.nz,int(std::lround(0.44/grid.dz))),grid.nr)]*std::sqrt(mu0)*1e3;
            std::printf("coil,%.2f,%d,%d,%s,%.3f,%.4f,%.4f,%.4f,%.3f,%.3f,%.6f,%.1e,%.1e,%.1e,%.2e,%.2e,%.1e,%.3f,%.3f,%.3f,%.3g,%llu,%llu,%ld,%.1f,%.3g,%s,%d,%.4g,%llu,%.4f%s%s\n",
                tesla,d.nz,d.nr,name(c),R.t*1e3,R.deviceThrust,R.exitPlaneThrust,coreThrust,R.inMom,R.wallMom,R.outMdot,R.massErr,R.energyErr,R.momErr,
                R.divMax,R.divMean,R.ctFaceDivMax,R.lorentzZ,R.coilLorentzZ,R.body,R.minBeta,(unsigned long long)R.steps,(unsigned long long)R.rejected,hllFallbacks,wall,
                double(R.steps)*2*d.nz*d.nr/wall,name(wb),int(robust),R.sync,(unsigned long long)R.robustSteps,wflux,failure.empty()?"":",FAILED: ",failure.c_str());
            std::fflush(stdout);
            if(dump) dumpField(s,"/tmp/crucible_mhd/coil_"+std::to_string(int(tesla*100))+"_"+name(c)+(c==Clean::Split?std::string("_")+name(wb)+(robust?"_robust":""):"")+"_"+std::to_string(d.nz)+".csv");
        }
    }
}
// Check of the exterior vacuum matrix: a current loop inside the nozzle (r = 0.01 m, z = 0.3 m)
// has a vacuum field outside the wall. Feed its exact wall psi in, compare the matched exterior
// field at the wall faces with the exact loop field there. Errors come from the triangle gradient
// (first order) and from truncating the exterior (far shell at Rf, b_r = 0 on the end planes).
void vacuumCheck(const std::vector<int>& levels) {
    std::printf("case,nz,Rf_m,layers,max_err_over_max_b,rms_err_over_max_b,err_mid_over_max_b,build_s\n");
    for(int nz:levels) for(double Rf:{0.3,1.0,3.0}) for(int nv:{40,60}) {
        double dz=0.6/nz; auto start=std::chrono::steady_clock::now();
        Vacuum v=buildVacuum(nz,dz,nozzleRadius,Rf,nv); double bs=seconds(start);
        auto P=[&](double z,double r){return loopPsi(r,z,0.01,0.3,1e3);};
        std::vector<double> wp(nz+1); for(int i=0;i<=nz;++i) wp[i]=P(i*dz,nozzleRadius(i*dz));
        double mx=0,rms=0,bmax=0,mid=0;
        std::vector<std::array<double,2>> ex(nz);
        for(int i=0;i<nz;++i) { double z=(i+0.5)*dz,r=0.5*(nozzleRadius(i*dz)+nozzleRadius((i+1)*dz)),h=1e-7;
            ex[i]={(P(z,r+h)-P(z,r-h))/(2*h*r),-(P(z+h,r)-P(z-h,r))/(2*h*r)}; bmax=std::max(bmax,std::hypot(ex[i][0],ex[i][1])); }
        for(int i=0;i<nz;++i) { auto b=v.field(i,wp); double e=std::hypot(b[0]-ex[i][0],b[1]-ex[i][1]); mx=std::max(mx,e); rms+=e*e/nz;
            if(i==nz/2) mid=e; }
        std::printf("vacuum,%d,%.1f,%d,%.4e,%.4e,%.4e,%.2f\n",nz,Rf,nv,mx/bmax,std::sqrt(rms)/bmax,mid/bmax,bs);
    }
}
std::vector<int> ints(const std::string& s) { std::vector<int> v; std::size_t p=0; while(p<s.size()){auto q=s.find(',',p); v.push_back(std::stoi(s.substr(p,q-p))); if(q==std::string::npos)break; p=q+1;} return v; }
std::vector<double> reals(const std::string& s) { std::vector<double> v; for(int x:ints(s)) v.push_back(x/100.0); return v; }
}

int main(int argc,char** argv) {
    std::string mode=argc>1?argv[1]:"";
    try {
        if(mode=="briowu") shockTube(ints(argc>2?argv[2]:"200,400,800,1600"));
        else if(mode=="aligned") aligned(ints(argc>2?argv[2]:"100,200,400"));
        else if(mode=="static") staticCoil(ints(argc>2?argv[2]:"40,80,160"),argc>3?std::stod(argv[3]):1.0);
        else if(mode=="residual") residual(ints(argc>2?argv[2]:"40,80,160,320"),argc>3?std::stod(argv[3]):1.0);
        else if(mode=="axisprobe") axisProbe(std::stoi(argv[2]),std::stod(argv[3]));
        else if(mode=="coil") nozzleCoil(ints(argc>2?argv[2]:"80"),reals(argc>3?argv[3]:"0,30,100"),argc>4?std::stod(argv[4]):0.008,argc>5);
        else if(mode=="split") nozzleCoil(ints(argc>2?argv[2]:"80"),reals(argc>3?argv[3]:"30,100"),argc>4?std::stod(argv[4]):0.008,argc>5,{Clean::CT,Clean::Split});
        else if(mode=="vacuum") vacuumCheck(ints(argc>2?argv[2]:"80,160"));
        else if(mode=="wall") {  // wall levels tesla*100 t_end dump|nodump transparent,insulating,conducting robust(0|1)
            std::vector<WallB> ws; std::string list=argc>6?argv[6]:"transparent,insulating,conducting";
            for(auto [n,w]:{std::pair{"transparent",WallB::Transparent},{"insulating",WallB::Insulating},{"conducting",WallB::Conducting}}) if(list.find(n)!=std::string::npos) ws.push_back(w);
            nozzleCoil(ints(argc>2?argv[2]:"80"),reals(argc>3?argv[3]:"100"),argc>4?std::stod(argv[4]):0.008,argc>5&&std::string(argv[5])=="dump",{Clean::Split},ws,argc>7&&std::string(argv[7])=="1");
        }
        else { std::fprintf(stderr,"crucible_mhd_spike briowu|aligned|static|coil [levels] [tesla*100 list] [t_end] [dump]\n"); return 1; }
    } catch(const std::exception& e) { std::fprintf(stderr,"error: %s\n",e.what()); return 2; }
}
