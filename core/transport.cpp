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
// their imposed mass, enthalpy and composition and exchange no diffusive flux.
#include "core/flow.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace crucible {
namespace {
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
    auto axial=[&](int face) {
        auto kind=face==0 && chamber?(faceSupply[j]<0?Neighbour::Wall:Neighbour::Supply):Neighbour::Open;
        return Neighbour{kind,0,face*m.dz-c.z,m.radius[face]*(j+0.5)/m.nr-c.r,-1,0};
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
    gradients_.resize(count*(3+ns_)*2);leastSquares_.resize(count);
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
    const auto& sp=medium_.species();
    for(std::size_t q=0;q<state_.size();++q) {
        const double* y=fractions_.data()+q*ns_;
        auto t=medium_.transport(temperature_[q],primitives_[q].p,y,diffusion_.data()+q*ns_,transportWork_);
        viscosity_[q]=t.viscosity;conductivity_[q]=t.conductivity;
        double moles=0;
        for(std::size_t k=0;k<ns_;++k) moles+=y[k]/sp[k].molarMass;
        for(std::size_t k=0;k<ns_;++k) moles_[q*ns_+k]=y[k]/sp[k].molarMass/moles;
    }
}
void Flow::transportFluxes(std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative,BoundaryRates& rates) {
    const auto& m=mesh_;const auto& d=definition_;const auto& sp=medium_.species();
    const std::size_t nf=3+ns_;
    const bool chamber=d.experiment==Case::Chamber;
    transportProperties();
    auto value=[&](std::size_t q,std::size_t f) {
        return f==0?primitives_[q].uz:f==1?primitives_[q].ur:f==2?temperature_[q]:moles_[q*ns_+f-3];
    };
    auto wallValue=[&](std::size_t q,std::size_t f,double nz,double nr) {
        if(f<2) {
            if(!d.wallSlip) return 0.0;
            const auto& w=primitives_[q];double un=w.uz*nz+w.ur*nr;
            return f==0?w.uz-un*nz:w.ur-un*nr;
        }
        if(f==2 && d.wallTemperature>0) return d.wallTemperature;
        return value(q,f);
    };
    for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
        auto q=m.index(i,j);auto stencil=neighbours(m,chamber,faceSupply_,i,j);const auto& inverse=leastSquares_[q];
        for(std::size_t f=0;f<nf;++f) {
            double v=value(q,f),sz=0,sr=0;
            for(const auto& n:stencil) {
                double delta=0;
                switch(n.kind) {
                    case Neighbour::Cell: delta=value(n.cell,f)-v;break;
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
    // Flux per unit area through a face with unit normal (nz, nr) pointing from cell a to cell b, or
    // out of cell a through a wall (wall != nullptr) whose midpoint is at its displacement. rf is the
    // face radius; axis marks the first face off the axis. Species fluxes go to jn.
    std::vector<double> yFace(ns_),jn(ns_);
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
        double fz=tzz*nz+tzr*nr,fr=tzr*nz+trr*nr,heat=-lambda*(tz*nz+tr*nr);
        std::fill(jn.begin(),jn.end(),0.0);
        if(boundary) {
            // A slip wall carries the normal stress only; an adiabatic wall no heat; no wall passes species.
            if(d.wallSlip) { double normal=fz*nz+fr*nr;fz=normal*nz;fr=normal*nr; }
            if(!(d.wallTemperature>0)) heat=0;
        } else if(ns_>1) {
            double rho=0.5*(primitives_[a].rho+primitives_[b].rho),moles=0,sum=0;
            for(std::size_t k=0;k<ns_;++k) { yFace[k]=0.5*(fractions_[a*ns_+k]+fractions_[b*ns_+k]);moles+=yFace[k]/sp[k].molarMass; }
            for(std::size_t k=0;k<ns_;++k) {
                double gz,gr;gradient(3+k,gz,gr);
                jn[k]=-rho*sp[k].molarMass*moles*0.5*(diffusion_[a*ns_+k]+diffusion_[b*ns_+k])*(gz*nz+gr*nr);
                sum+=jn[k];
            }
            medium_.speciesEnthalpies(0.5*(temperature_[a]+temperature_[b]),faceEnthalpy_.data());
            for(std::size_t k=0;k<ns_;++k) { jn[k]-=yFace[k]*sum;heat+=faceEnthalpy_[k]*jn[k]; }
        }
        return Conserved{0,-fz,-fr,-(fz*uz+fr*ur)+heat};
    };
    auto exchange=[&](std::size_t a,std::size_t b,double area,const Conserved& flux) {
        for(int k=1;k<4;++k) { derivative[a][k]-=area*flux[k];derivative[b][k]+=area*flux[k]; }
        for(std::size_t k=0;k<ns_;++k) { speciesDerivative[a*ns_+k]-=area*jn[k];speciesDerivative[b*ns_+k]+=area*jn[k]; }
    };
    for(int i=0;i<=m.nz;++i) for(int j=0;j<m.nr;++j) {
        double area=m.axialArea(i,j),rf=m.radius[i]*(j+0.5)/m.nr;
        if(i>0 && i<m.nz) { exchange(m.index(i-1,j),m.index(i,j),area,faceFlux(m.index(i-1,j),m.index(i,j),nullptr,1,0,rf,false));continue; }
        auto q=m.index(i==0?0:m.nz-1,j);
        auto n=neighbours(m,chamber,faceSupply_,i==0?0:m.nz-1,j)[i==0?0:1];
        Conserved flux{};
        if(n.kind==Neighbour::Supply) continue;
        if(n.kind==Neighbour::Wall) flux=faceFlux(q,0,&n,1,0,rf,false);  // the injector plate
        else {
            // Zero normal gradient: tau_zz = -2/3 mu (du_r/dr + u_r/r), tau_zr = mu du_z/dr, no heat or species.
            const double* g=&gradients_[q*nf*2];
            double mu=viscosity_[q],ur=primitives_[q].ur;
            double tzz=-2.0/3*mu*(g[3]+ur/m.cells[q].r),tzr=mu*g[1];
            flux={0,-tzz,-tzr,-(tzz*primitives_[q].uz+tzr*ur)};
        }
        // The flux is along +z: it enters the cell at i = 0 and leaves it at i = nz.
        double sign=i==0?1:-1;
        for(int k=1;k<4;++k) derivative[q][k]+=sign*area*flux[k];
        rates.energy+=sign*area*flux[3];
        if(n.kind==Neighbour::Wall) { rates.wallAxial+=area*flux[1];rates.wallHeat+=area*flux[3]; }
        else if(i==0) rates.inletMomentum+=area*flux[1];
        else rates.outletMomentum+=area*flux[1];
    }
    for(int i=0;i<m.nz;++i) for(int j=1;j<=m.nr;++j) {
        auto ar=m.radialAreaVector(i,j);double area=std::hypot(ar[0],ar[1]);
        double nz=ar[0]/area,nr=ar[1]/area,rf=0.5*(m.radius[i]+m.radius[i+1])*j/m.nr;
        auto a=m.index(i,j-1);
        if(j<m.nr) { exchange(a,m.index(i,j),area,faceFlux(a,m.index(i,j),nullptr,nz,nr,rf,j==1));continue; }
        auto n=neighbours(m,chamber,faceSupply_,i,j-1)[3];
        auto flux=faceFlux(a,0,&n,nz,nr,rf,false);
        for(int k=1;k<4;++k) derivative[a][k]-=area*flux[k];
        rates.wallAxial-=area*flux[1];rates.energy-=area*flux[3];rates.wallHeat-=area*flux[3];
    }
    for(std::size_t q=0;q<state_.size();++q) {
        const auto& c=m.cells[q];
        double ur=primitives_[q].ur,div=gradients_[q*nf*2]+gradients_[(q*nf+1)*2+1]+ur/c.r;
        derivative[q][2]-=viscosity_[q]*(2*ur/c.r-2.0/3*div)*c.radialPressureMeasure;
    }
}
void Flow::transportDerivative(std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative) {
    if(!medium_.hasTransport()) throw std::logic_error("The medium has no transport data.");
    derivative.assign(state_.size(),Conserved{});speciesDerivative.assign(species_.size(),0.0);
    refresh(state_,species_);
    BoundaryRates rates{};
    transportFluxes(derivative,speciesDerivative,rates);
    for(std::size_t q=0;q<state_.size();++q) {
        for(double& v:derivative[q]) v/=mesh_.cells[q].volume;
        for(std::size_t k=0;k<ns_;++k) speciesDerivative[q*ns_+k]/=mesh_.cells[q].volume;
    }
}
} // namespace crucible
