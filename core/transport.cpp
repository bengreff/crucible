// Molecular transport on the meridional mesh: viscous stress, Fourier conduction and mixture-averaged
// species diffusion with the correction velocity that makes the diffusive mass fluxes sum to zero
// (the flux of Cantera's mixture-averaged model with Soret diffusion off):
//   j_k = -rho (W_k / W) D_km grad X_k + Y_k sum_l rho (W_l / W) D_lm grad X_l,
//   q = -lambda grad T + sum_k h_k j_k,
//   tau = mu (grad u + grad u^T - 2/3 div u I), div u = du_z/dz + du_r/dr + u_r / r (no swirl).
// The face flux per unit area is {0, -(tau.n)_z, -(tau.n)_r, -(tau.n).u + q.n} and j_k.n.
//
// Gradients: weighted least squares (weights 1/d^2) over the four face neighbours of each cell,
// exact for linear fields on the body-fitted mesh. A missing neighbour is replaced by a boundary
// point: the mirror image of the centroid across the axis (even fields keep their value, u_r
// changes sign); the wall face midpoint with the wall value (no-slip: u = 0; slip: the tangential
// part of the cell velocity; isothermal: T_w; otherwise the cell value); and the face midpoint with
// the cell value at open boundaries (supply rings, nozzle inlet, outlet). Face gradients are the
// mean of the two cell gradients with the component along the line of centres replaced by the
// compact difference; a wall face uses the cell gradient and the difference to the wall value.
// Across the first face off the axis the even fields (u_z, T, X_k) are differenced in r^2 with the
// cells' exact second moments, d phi/dr = 2 r (phi_1 - phi_0) / (<r^2>_1 - <r^2>_0), exact for
// phi = a + b r^2 (a difference in r would give the slope at the mean centroid radius, 11% off).
// u_r / r is interpolated to faces as itself (it is smooth and even; the mean of u_r divided by the
// face radius is O(1) wrong next to the axis). Face properties are the mean of the two cells (a wall
// face uses the cell's). The hoop stress enters the radial momentum as -tau_thth integral(dV / r),
// with tau_thth = mu (2 u_r / r - 2/3 div u) at the centroid, as the pressure enters as +p integral(dV / r).
// Open boundaries (nozzle inlet, outlet) take a zero normal gradient (fully developed flow): only
// the stress of the cell's tangential gradients crosses them. Supply rings deliver plug flow with
// their imposed mass, enthalpy and composition and exchange no diffusive flux (with turbulence the
// stream's own Reynolds normal stress 2/3 rho k is in its convective flux, Flow::supplyFace).
//
// SST-2003 (Definition::Turbulence; NASA TMR sst.html, the 2003 constants): with the eddy viscosity
// mu_t = rho a1 k / max(a1 omega, S F2), S = sqrt(2 S_ij S_ij) with the hoop strain, the stress
// becomes (mu + mu_t)(grad u + grad u^T - 2/3 div u I) - 2/3 rho k I, the conductivity
// lambda + cp mu_t / Pr_t, every D_km gains mu_t / (rho Sc_t) (with the correction velocity this
// adds exactly -(mu_t / Sc_t) grad Y_k), and k and omega diffuse with mu + sigma mu_t. The energy
// flux carries the k flux, since k is part of E. Gradients of k and omega are taken like the other
// even fields. A no-slip wall face has mu_t = 0, k = 0 and the wall omega of Definition::Turbulence;
// its k flux (zero in the continuum, where k ~ y^2) stays in the cell as heat, so the energy
// through the wall is the conduction alone. A slip wall passes no k or omega flux, and its normal
// stress includes -2/3 rho k of the cell. Face eddy properties are the mean of the two cells.
// Sources: P = mu_t (S^2 - 2/3 div^2) - 2/3 rho k div, Pt = min(P, 10 beta* rho omega k),
//   d(rho k)/dt = Pt - beta* rho omega k,
//   d(rho omega)/dt = gamma rho Pt / mu_t - beta rho omega^2 + 2 (1 - F1) rho sigma_w2 / omega grad k . grad omega,
// integrated per cell at fixed rho and E by the second-order positive modified Patankar
// Runge-Kutta scheme MPRK22 (Kopecz and Meister, BIT 58, 2018; production explicit, destruction
// weighted by the new over the old value), with the strain, divergence, grad k . grad omega, nu
// and the wall distance frozen and F1, F2 and mu_t evaluated at each stage's k and omega.
#include "core/flow.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace crucible {
namespace {
double sq(double x) { return x*x; }
// SST-2003 constants (sets 1 and 2, blended by F1).
constexpr double sigmaK1=0.85,sigmaK2=1.0,sigmaW1=0.5,sigmaW2=0.856,beta1=0.075,beta2=0.0828;
constexpr double gamma1=5.0/9,gamma2=0.44,betaStar=0.09,a1=0.31;
struct Blending { double f1, f2; };
// F1 and F2 (sst.html); inverseDistance = 0 away from every wall gives F1 = F2 = 0.
Blending blending(double k,double omega,double rho,double nu,double inverseDistance,double cross) {
    const double root=std::sqrt(k),d2=sq(inverseDistance);
    const double cd=std::max(2*rho*sigmaW2/omega*cross,1e-10);
    const double viscous=500*nu*d2/omega;
    const double arg1=std::min(std::max(root/(betaStar*omega)*inverseDistance,viscous),4*rho*sigmaW2*k*d2/cd);
    const double arg2=std::max(2*root/(betaStar*omega)*inverseDistance,viscous);
    return {std::tanh(sq(sq(arg1))),std::tanh(sq(arg2))};
}
// A least-squares neighbour of a cell: an interior cell, the mirror image across the axis, or a
// boundary face midpoint (a wall or an open boundary). (dz, dr) is its displacement from the
// centroid; (nz, nr) is a wall's unit normal (either orientation).
struct Neighbour {
    enum Kind { Cell, Axis, Wall, Open, Supply } kind;
    std::size_t cell;
    double dz, dr, nz, nr;
};
std::array<Neighbour, 4> neighbours(const Mesh& m, bool chamber, const std::vector<int>& faceSupply, int i, int j) {
    const auto& c=m.cells[m.index(i,j)];
    std::array<Neighbour, 4> out{};
    auto cell=[&](int a,int b){ auto q=m.index(a,b);return Neighbour{Neighbour::Cell,q,m.cells[q].z-c.z,m.cells[q].r-c.r,0,0}; };
    // Axial faces are normal to z. A wall (the injector plate) has its value at the face midpoint; an
    // open or supply face has zero normal gradient, which is the cell's own value on its axial line,
    // so that neighbour sits at the face at the centroid's radius (at the midpoint it would also
    // claim a radial derivative of zero and bias the end columns' radial gradients).
    auto axial=[&](int face) {
        auto kind=face==0 && chamber?(faceSupply[j]<0?Neighbour::Wall:Neighbour::Supply):Neighbour::Open;
        return Neighbour{kind,0,face*m.dz-c.z,kind==Neighbour::Wall?m.ringMiddle(face,j)-c.r:0.0,-1,0};
    };
    out[0]=i>0?cell(i-1,j):axial(0);
    out[1]=i<m.nz-1?cell(i+1,j):axial(m.nz);
    out[2]=j>0?cell(i,j-1):Neighbour{Neighbour::Axis,0,0,-2*c.r,0,0};
    if(j<m.nr-1) out[3]=cell(i,j+1);
    else {
        auto ar=m.radialAreaVector(i,m.nr);double area=std::hypot(ar[0],ar[1]);
        out[3]={Neighbour::Wall,0,(i+0.5)*m.dz-c.z,0.5*(m.radius[i]+m.radius[i+1])-c.r,ar[0]/area,ar[1]/area};
    }
    return out;
}
}

void Flow::prepareTransport() {
    const auto count=mesh_.cells.size();
    viscosity_.resize(count);conductivity_.resize(count);
    diffusion_.resize(count*ns_);moles_.resize(count*ns_);faceEnthalpy_.resize(ns_);
    gradients_.resize(count*(3+ns_+nt_)*2);leastSquares_.resize(count);
    if(nt_) {
        eddy_.assign(count,0);eddyConductivity_.assign(count,0);eddyDiffusion_.assign(count*2,0);
        wallOmega_.assign(count,0);wallDiffusion_.resize(ns_);sources_.resize(count);
    }
    const bool chamber=definition_.experiment==Case::Chamber;
    for(int i=0;i<mesh_.nz;++i) for(int j=0;j<mesh_.nr;++j) {
        double a=0,b=0,c=0;
        for(const auto& n:neighbours(mesh_,chamber,faceSupply_,i,j)) {
            double w=1/(n.dz*n.dz+n.dr*n.dr);
            a+=w*n.dz*n.dz;b+=w*n.dz*n.dr;c+=w*n.dr*n.dr;
        }
        double det=a*c-b*b;
        if(!(det>0)) throw std::runtime_error("Degenerate least-squares gradient stencil.");
        leastSquares_[mesh_.index(i,j)]={c/det,-b/det,a/det};
    }
}
void Flow::transportProperties() {
    const auto& sp=medium_.species();const auto& d=definition_;
    for(std::size_t q=0;q<state_.size();++q) {
        const double* y=fractions_.data()+q*nw_;
        auto t=medium_.transport(temperature_[q],primitives_[q].p,y,diffusion_.data()+q*ns_,transportWork_);
        viscosity_[q]=t.viscosity;conductivity_[q]=t.conductivity;
        double moles=0;
        for(std::size_t k=0;k<ns_;++k) moles+=y[k]/sp[k].molarMass;
        for(std::size_t k=0;k<ns_;++k) moles_[q*ns_+k]=y[k]/sp[k].molarMass/moles;
    }
    if(!nt_ || d.wallSlip) return;
    // Wall-face omega = factor 6 nu_w / (beta1 d1^2) of the cells next to a no-slip wall, nu_w at the
    // wall temperature and the cell's pressure and composition (the cell's nu if adiabatic).
    auto wall=[&](std::size_t q) {
        const double* y=fractions_.data()+q*nw_;
        double nu=viscosity_[q]/primitives_[q].rho;
        if(d.wallTemperature>0) {
            auto t=medium_.transport(d.wallTemperature,primitives_[q].p,y,wallDiffusion_.data(),transportWork_);
            nu=t.viscosity*medium_.gasConstant(y)*d.wallTemperature/primitives_[q].p;
        }
        wallOmega_[q]=d.turbulence.wallOmegaFactor*6*nu/(beta1*sq(wallDistance_[q]));
    };
    for(int i=0;i<mesh_.nz;++i) wall(mesh_.index(i,mesh_.nr-1));
    if(d.experiment==Case::Chamber) for(int j=0;j<mesh_.nr;++j) if(faceSupply_[j]<0) wall(mesh_.index(0,j));
}
double Flow::transportValue(std::size_t q,std::size_t f) const {
    return f==0?primitives_[q].uz:f==1?primitives_[q].ur:f==2?temperature_[q]:f<3+ns_?moles_[q*ns_+f-3]:fractions_[q*nw_+f-3];
}
double Flow::wallValue(std::size_t q,std::size_t f,double nz,double nr) const {
    const auto& d=definition_;
    if(f<2) {
        if(!d.wallSlip) return 0.0;
        const auto& w=primitives_[q];double un=w.uz*nz+w.ur*nr;
        return f==0?w.uz-un*nz:w.ur-un*nr;
    }
    if(f==2 && d.wallTemperature>0) return d.wallTemperature;
    if(f>=3+ns_ && !d.wallSlip) return f==3+ns_?0.0:wallOmega_[q];
    return transportValue(q,f);
}
void Flow::transportGradients() {
    const auto& m=mesh_;
    const std::size_t nf=3+ns_+nt_;
    const bool chamber=definition_.experiment==Case::Chamber;
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
        auto q=m.index(i,j);auto stencil=neighbours(m,chamber,faceSupply_,i,j);const auto& inverse=leastSquares_[q];
        for(std::size_t f=0;f<nf;++f) {
            double v=transportValue(q,f),sz=0,sr=0;
            for(const auto& n:stencil) {
                double delta=0;
                switch(n.kind) {
                    case Neighbour::Cell: delta=transportValue(n.cell,f)-v;break;
                    case Neighbour::Axis: delta=f==1?-2*v:0;break;
                    case Neighbour::Wall: delta=wallValue(q,f,n.nz,n.nr)-v;break;
                    case Neighbour::Open: case Neighbour::Supply: break;
                }
                double w=delta/(n.dz*n.dz+n.dr*n.dr);
                sz+=w*n.dz;sr+=w*n.dr;
            }
            gradients_[(q*nf+f)*2]=inverse[0]*sz+inverse[1]*sr;
            gradients_[(q*nf+f)*2+1]=inverse[1]*sz+inverse[2]*sr;
        }
    }
}
void Flow::eddyViscosity(bool withSources) {
    const auto& tu=definition_.turbulence;
    const std::size_t nf=3+ns_+nt_;
    for(std::size_t q=0;q<state_.size();++q) {
        const auto& w=primitives_[q];const double* g=&gradients_[q*nf*2];const double* y=fractions_.data()+q*nw_;
        const double k=y[ns_],omega=y[ns_+1],hoop=w.ur/mesh_.cells[q].r;
        const double div=g[0]+g[3]+hoop;
        const double strain2=2*(sq(g[0])+sq(g[3])+sq(hoop))+sq(g[1]+g[2]),strain=std::sqrt(strain2);
        const double* gk=g+2*(3+ns_);const double* gw=g+2*(4+ns_);
        const double cross=gk[0]*gw[0]+gk[1]*gw[1],nu=viscosity_[q]/w.rho,inverse=1/wallDistance_[q];
        const auto b=blending(k,omega,w.rho,nu,inverse,cross);
        const double mut=w.rho*a1*k/std::max(a1*omega,strain*b.f2);
        auto props=medium_.properties(temperature_[q],y);
        eddy_[q]=mut;eddyConductivity_[q]=(props.cv+props.r)*mut/tu.prandtl;
        eddyDiffusion_[2*q]=(b.f1*sigmaK1+(1-b.f1)*sigmaK2)*mut;
        eddyDiffusion_[2*q+1]=(b.f1*sigmaW1+(1-b.f1)*sigmaW2)*mut;
        if(withSources) sources_[q]={strain2-2.0/3*sq(div),strain,div,cross,nu,inverse};
    }
}
void Flow::mixingInputs(double cmix,const std::vector<std::size_t>& species,std::vector<double>& time,
                        std::vector<double>& segregation) {
    if(!nt_) throw std::logic_error("The PaSR closure needs turbulence.");
    for(std::size_t i:species) if(i>=ns_) throw std::invalid_argument("PaSR species index out of range.");
    refresh(state_,species_,turbulence_);
    transportProperties();transportGradients();eddyViscosity(false);
    const auto& tu=definition_.turbulence;
    const std::size_t nf=3+ns_+nt_,count=state_.size();
    time.resize(count);segregation.resize(count);
    for(std::size_t q=0;q<count;++q) {
        const double rho=primitives_[q].rho,k=fractions_[q*nw_+ns_],omega=fractions_[q*nw_+ns_+1];
        const double epsilon=betaStar*k*omega,nuT=eddy_[q]/rho;
        time[q]=epsilon>0?cmix*std::sqrt((viscosity_[q]+eddy_[q])/rho/epsilon):std::numeric_limits<double>::infinity();
        double s=0;
        for(std::size_t i:species) {
            const double x=moles_[q*ns_+i];
            if(!(x>0 && x<1)) continue;
            const double* g=&gradients_[(q*nf+3+i)*2];
            s=std::max(s,std::min(1.0,nuT/tu.schmidt*(sq(g[0])+sq(g[1]))/(betaStar*omega*x*(1-x))));
        }
        segregation[q]=s;
    }
}
// MPRK22 on specific k and omega (rho is fixed): with production P and destruction D = d y (d >= 0),
//   y1 = (y0 + tau P0) / (1 + tau d0),  y = (y0 + tau/2 (P0 + P1)) / (1 + tau/2 (d0 y0 + d1 y1) / y1).
// Each source term goes to P or D by its sign; positivity of k and omega is unconditional.
void Flow::turbulenceSource(const std::vector<Conserved>& state,std::vector<double>& turbulence,
                            const std::vector<SourceCoefficients>& coefficients,double tau) const {
    for(std::size_t q=0;q<state.size();++q) {
        const auto& c=coefficients[q];const double rho=state[q][0];
        // Production and destruction coefficient of k and omega at (k, omega).
        auto rates=[&](double k,double omega,double& pk,double& dk,double& pw,double& dw) {
            const auto b=blending(k,omega,rho,c.nu,c.inverseDistance,c.crossGradient);
            const double m=std::max(a1*omega,c.strain*b.f2);
            // Pt / (rho k): mu_t / (rho k) = a1 / m.
            const double rate=std::min(a1*c.strain2/m-2.0/3*c.divergence,10*betaStar*omega);
            const double gamma=b.f1*gamma1+(1-b.f1)*gamma2,beta=b.f1*beta1+(1-b.f1)*beta2;
            const double production=gamma*rate*m/a1,cross=2*(1-b.f1)*sigmaW2*c.crossGradient;
            pk=k*std::max(rate,0.0);dk=std::max(-rate,0.0)+betaStar*omega;
            pw=std::max(production,0.0)+std::max(cross,0.0)/omega;
            dw=std::max(-production,0.0)/omega+beta*omega+std::max(-cross,0.0)/sq(omega);
        };
        const double k0=turbulence[q*nt_]/rho,w0=turbulence[q*nt_+1]/rho;
        double pk0,dk0,pw0,dw0,pk1,dk1,pw1,dw1;
        rates(k0,w0,pk0,dk0,pw0,dw0);
        const double k1=(k0+tau*pk0)/(1+tau*dk0),w1=(w0+tau*pw0)/(1+tau*dw0);
        rates(k1,w1,pk1,dk1,pw1,dw1);
        const double k=(k0+0.5*tau*(pk0+pk1))/(1+0.5*tau*(k1>0?(dk0*k0+dk1*k1)/k1:dk1));
        const double w=(w0+0.5*tau*(pw0+pw1))/(1+0.5*tau*(dw0*w0+dw1*w1)/w1);
        turbulence[q*nt_]=rho*k;turbulence[q*nt_+1]=rho*w;
    }
}
void Flow::transportFluxes(std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative,
                           std::vector<double>& turbulenceDerivative,BoundaryRates& rates) {
    const auto& m=mesh_;const auto& d=definition_;const auto& sp=medium_.species();
    const std::size_t nf=3+ns_+nt_;
    const bool chamber=d.experiment==Case::Chamber;
    transportProperties();
    transportGradients();
    if(nt_) eddyViscosity(false);
    auto value=[&](std::size_t q,std::size_t f) { return transportValue(q,f); };
    auto rhoK=[&](std::size_t q) { return primitives_[q].rho*fractions_[q*nw_+ns_]; };
    // Flux per unit area through a face with unit normal (nz, nr) pointing from cell a to cell b, or
    // out of cell a through a wall (wall != nullptr) whose midpoint is at its displacement. rf is the
    // face radius; axis marks the first face off the axis. Species fluxes go to jn, k and omega fluxes to tn.
    std::vector<double> yFace(ns_),jn(ns_);
    std::array<double,2> tn{};
    auto faceFlux=[&](std::size_t a,std::size_t b,const Neighbour* wall,double nz,double nr,double rf,bool axis) {
        const bool boundary=wall!=nullptr;
        double dz=boundary?wall->dz:m.cells[b].z-m.cells[a].z,dr=boundary?wall->dr:m.cells[b].r-m.cells[a].r;
        double len=std::hypot(dz,dr),ez=dz/len,er=dr/len;
        auto other=[&](std::size_t f){ return boundary?wallValue(a,f,wall->nz,wall->nr):value(b,f); };
        auto gradient=[&](std::size_t f,double& gz,double& gr) {
            const double* ga=&gradients_[(a*nf+f)*2];
            gz=ga[0];gr=ga[1];
            if(!boundary) { const double* gb=&gradients_[(b*nf+f)*2];gz=0.5*(gz+gb[0]);gr=0.5*(gr+gb[1]); }
            double along=(other(f)-value(a,f))/len;
            if(axis && f!=1) along=2*rf*(value(b,f)-value(a,f))/(m.cells[b].radialSecondMoment-m.cells[a].radialSecondMoment)*er;
            double c=along-(gz*ez+gr*er);
            gz+=c*ez;gr+=c*er;
        };
        auto mean=[&](double x,double y){ return boundary?x:0.5*(x+y); };
        double uzz,uzr,urz,urr,tz,tr;
        gradient(0,uzz,uzr);gradient(1,urz,urr);gradient(2,tz,tr);
        double uz=boundary?other(0):0.5*(value(a,0)+value(b,0)),ur=boundary?other(1):0.5*(value(a,1)+value(b,1));
        double urOverR=boundary?ur/rf:0.5*(value(a,1)/m.cells[a].r+value(b,1)/m.cells[b].r);
        double mu=mean(viscosity_[a],boundary?0:viscosity_[b]),lambda=mean(conductivity_[a],boundary?0:conductivity_[b]);
        double div=uzz+urr+urOverR;
        double tzz=mu*(2*uzz-2.0/3*div),trr=mu*(2*urr-2.0/3*div),tzr=mu*(uzr+urz);
        double normalGradientT=tz*nz+tr*nr,muT=0;
        tn={0,0};
        if(nt_) {
            // A no-slip wall face takes the wall's mu_t = 0 and k = 0; a slip wall face the cell's.
            const bool noSlip=boundary && !d.wallSlip;
            muT=noSlip?0:mean(eddy_[a],boundary?0:eddy_[b]);
            double k=noSlip?0:mean(rhoK(a),boundary?0:rhoK(b));
            tzz+=muT*(2*uzz-2.0/3*div)-2.0/3*k;trr+=muT*(2*urr-2.0/3*div)-2.0/3*k;tzr+=muT*(uzr+urz);
            lambda+=noSlip?0:mean(eddyConductivity_[a],boundary?0:eddyConductivity_[b]);
            if(!(boundary && d.wallSlip)) for(std::size_t t=0;t<nt_;++t) {
                double gz,gr;gradient(3+ns_+t,gz,gr);
                double diffusivity=mu+(noSlip?0:mean(eddyDiffusion_[2*a+t],boundary?0:eddyDiffusion_[2*b+t]));
                tn[t]=-diffusivity*(gz*nz+gr*nr);
            }
        }
        double fz=tzz*nz+tzr*nr,fr=tzr*nz+trr*nr,heat=-lambda*normalGradientT;
        std::fill(jn.begin(),jn.end(),0.0);
        if(boundary) {
            // A slip wall carries the normal stress only; an adiabatic wall no heat; no wall passes species.
            if(d.wallSlip) { double normal=fz*nz+fr*nr;fz=normal*nz;fr=normal*nr; }
            if(!(d.wallTemperature>0)) heat=0;
        } else {
            if(nt_) heat+=tn[0];  // k is part of E
            if(ns_>1) {
                double rho=0.5*(primitives_[a].rho+primitives_[b].rho),moles=0,sum=0;
                for(std::size_t k=0;k<ns_;++k) { yFace[k]=0.5*(fractions_[a*nw_+k]+fractions_[b*nw_+k]);moles+=yFace[k]/sp[k].molarMass; }
                for(std::size_t k=0;k<ns_;++k) {
                    double gz,gr;gradient(3+k,gz,gr);
                    jn[k]=-rho*sp[k].molarMass*moles*0.5*(diffusion_[a*ns_+k]+diffusion_[b*ns_+k])*(gz*nz+gr*nr);
                    if(nt_) jn[k]-=sp[k].molarMass*moles*muT/d.turbulence.schmidt*(gz*nz+gr*nr);
                    sum+=jn[k];
                }
                medium_.speciesEnthalpies(0.5*(temperature_[a]+temperature_[b]),faceEnthalpy_.data());
                for(std::size_t k=0;k<ns_;++k) { jn[k]-=yFace[k]*sum;heat+=faceEnthalpy_[k]*jn[k]; }
            }
        }
        return Conserved{0,-fz,-fr,-(fz*uz+fr*ur)+heat};
    };
    auto exchange=[&](std::size_t a,std::size_t b,double area,const Conserved& flux) {
        for(int k=1;k<4;++k) { derivative[a][k]-=area*flux[k];derivative[b][k]+=area*flux[k]; }
        for(std::size_t k=0;k<ns_;++k) { speciesDerivative[a*ns_+k]-=area*jn[k];speciesDerivative[b*ns_+k]+=area*jn[k]; }
        for(std::size_t k=0;k<nt_;++k) { turbulenceDerivative[a*nt_+k]-=area*tn[k];turbulenceDerivative[b*nt_+k]+=area*tn[k]; }
    };
    for(int i=0;i<=m.nz;++i) for(int j=0;j<m.nr;++j) {
        double area=m.axialArea(i,j),rf=m.ringMiddle(i,j);
        if(i>0 && i<m.nz) { exchange(m.index(i-1,j),m.index(i,j),area,faceFlux(m.index(i-1,j),m.index(i,j),nullptr,1,0,rf,false));continue; }
        auto q=m.index(i==0?0:m.nz-1,j);
        auto n=neighbours(m,chamber,faceSupply_,i==0?0:m.nz-1,j)[i==0?0:1];
        Conserved flux{};
        if(n.kind==Neighbour::Supply) continue;
        // The flux is along +z: it enters the cell at i = 0 and leaves it at i = nz.
        double sign=i==0?1:-1;
        if(n.kind==Neighbour::Wall) {  // the injector plate
            flux=faceFlux(q,0,&n,1,0,rf,false);
            for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]+=sign*area*tn[k];
        } else {
            // Zero normal gradient: tau_zz = -2/3 mu (du_r/dr + u_r/r), tau_zr = mu du_z/dr, no heat,
            // species, k or omega flux; with turbulence mu + mu_t and the normal stress -2/3 rho k.
            const double* g=&gradients_[q*nf*2];
            double mu=viscosity_[q],ur=primitives_[q].ur;
            if(nt_) mu+=eddy_[q];
            double tzz=-2.0/3*mu*(g[3]+ur/m.cells[q].r),tzr=mu*g[1];
            if(nt_) tzz-=2.0/3*rhoK(q);
            flux={0,-tzz,-tzr,-(tzz*primitives_[q].uz+tzr*ur)};
        }
        for(int k=1;k<4;++k) derivative[q][k]+=sign*area*flux[k];
        rates.energy+=sign*area*flux[3];
        if(n.kind==Neighbour::Wall) { rates.wallAxial+=area*flux[1];rates.wallHeat+=area*flux[3]; }
        else if(i==0) rates.inletMomentum+=area*flux[1];
        else rates.outletMomentum+=area*flux[1];
    }
    for(int i=0;i<m.nz;++i) for(int j=1;j<=m.nr;++j) {
        auto ar=m.radialAreaVector(i,j);double area=std::hypot(ar[0],ar[1]);
        double nz=ar[0]/area,nr=ar[1]/area,rf=m.radialFaceMiddle(i,j);
        auto a=m.index(i,j-1);
        if(j<m.nr) { exchange(a,m.index(i,j),area,faceFlux(a,m.index(i,j),nullptr,nz,nr,rf,j==1));continue; }
        auto n=neighbours(m,chamber,faceSupply_,i,j-1)[3];
        auto flux=faceFlux(a,0,&n,nz,nr,rf,false);
        for(int k=1;k<4;++k) derivative[a][k]-=area*flux[k];
        for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[a*nt_+k]-=area*tn[k];
        rates.wallAxial-=area*flux[1];rates.energy-=area*flux[3];rates.wallHeat-=area*flux[3];
    }
    // Hoop stress: tau_thth = (mu + mu_t)(2 u_r / r - 2/3 div u) - 2/3 rho k at the centroid.
    for(std::size_t q=0;q<state_.size();++q) {
        const auto& c=m.cells[q];
        double ur=primitives_[q].ur,div=gradients_[q*nf*2]+gradients_[(q*nf+1)*2+1]+ur/c.r;
        double hoop=viscosity_[q]*(2*ur/c.r-2.0/3*div);
        if(nt_) hoop+=eddy_[q]*(2*ur/c.r-2.0/3*div)-2.0/3*rhoK(q);
        derivative[q][2]-=hoop*c.radialPressureMeasure;
    }
}
void Flow::transportDerivative(std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative) {
    std::vector<double> turbulenceDerivative;
    transportDerivative(derivative,speciesDerivative,turbulenceDerivative);
}
void Flow::transportDerivative(std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative,
                               std::vector<double>& turbulenceDerivative) {
    if(!medium_.hasTransport()) throw std::logic_error("The medium has no transport data.");
    derivative.assign(state_.size(),Conserved{});speciesDerivative.assign(species_.size(),0.0);
    turbulenceDerivative.assign(turbulence_.size(),0.0);
    refresh(state_,species_,turbulence_);
    BoundaryRates rates{};
    transportFluxes(derivative,speciesDerivative,turbulenceDerivative,rates);
    for(std::size_t q=0;q<state_.size();++q) {
        for(double& v:derivative[q]) v/=mesh_.cells[q].volume;
        for(std::size_t k=0;k<ns_;++k) speciesDerivative[q*ns_+k]/=mesh_.cells[q].volume;
        for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]/=mesh_.cells[q].volume;
    }
}
} // namespace crucible
