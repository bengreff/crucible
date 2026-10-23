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
// Sources, in TMR's SSTs form (the -2/3 rho k delta_ij term kept in tau_ij for the momentum and
// energy, the production approximated by P = mu_t S^2 in both equations; adopted 5 October 2026 in
// place of the exact P, whose dilatation part -2/3 rho k div collapsed omega in the expanding nozzle,
// docs/evidence/PASR_C2.md): Pt = min(mu_t S^2, 10 beta* rho omega k),
//   d(rho k)/dt = Pt - beta* rho omega k,
//   d(rho omega)/dt = gamma rho Pt / mu_t - beta rho omega^2 + 2 (1 - F1) rho sigma_w2 / omega grad k . grad omega.
// Both productions are non-negative. The turbulent pressure's expansion work u . grad(2/3 rho k) is
// in the energy flux and does not enter k; since k is part of E, the total energy is unchanged. Sources are
// integrated per cell at fixed rho and E by the second-order positive modified Patankar
// Runge-Kutta scheme MPRK22 (Kopecz and Meister, BIT 58, 2018; production explicit, destruction
// weighted by the new over the old value), with S^2, S, grad k . grad omega, nu
// and the wall distance frozen and F1, F2 and mu_t evaluated at each stage's k and omega.
//
// Wall functions (Definition::Turbulence::wallFunctions; docs/evidence/WALL_FUNCTIONS.md): each no-slip
// wall face takes Nichols and Nelson's law (crucible::wallLaw::solveIsothermal) at the first cell's
// state, and its viscous flux is replaced by the wall's traction tau_w against the cell's tangential
// velocity and the heat q_w, with no species, k or omega flux. rho_w is the cell's pressure at T_w,
// mu_w and k_w are the mixture-averaged transport at T_w with the cell's composition, cp is the layer
// mean (h(T_1) - h(T_w)) / (T_1 - T_w) and r = Pr_w^(1/3). T_w and rho_w are never written to the
// state. The first cell's k and omega are prescribed after every stage (Flow::step): mu_t from
// eq. 10.12 (negative values set to zero), omega = (omega_i^2 + omega_o^2)^(1/2) with
// omega_i = 6 mu_w / (beta1 rho_w y^2) and omega_o = u_tau / (sqrt(beta*) kappa y), and k = omega mu_t / rho
// (eqs. 10.17 to 10.20); its turbulence equations are not solved, and its eddy viscosity is rho k / omega,
// which is that mu_t (the stress limiter is the resolved model's; eq. 10.12 is the constant-stress
// layer's mu_t). For the gradients the wall face carries u = 0 and T_w as before, and the cell's own k
// and omega (no k or omega flux). The law folds back above T_1 / T_w = 11 (WALL_FUNCTIONS.md); there the
// engine throws.
//
// Threads (Flow::pool): the per-cell loops run in blocks of cells. Each face's flux is computed alone
// into a face array, and each cell then sums its faces with the serial loops' arithmetic and order
// (axial face i, i + 1, radial face j, j + 1, the hoop stress); the boundary rates are summed serially
// in the order of those loops. The result is the same bit for bit on any number of threads.
#include "core/flow.hpp"
#include "core/walls.hpp"
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
// D^ic = D(R_t, M_t) / D(R_t, 0), D = [1 - exp(-R_t / (3.5 + 0.39 M_t^0.77))]^2 (Hasan et al. eqs. 4.3 to 4.5).
double intrinsicDamping(double rt,double mt) {
    rt=std::max(rt,1e-12);
    return sq(std::expm1(-rt/(3.5+0.39*std::pow(mt,0.77)))/std::expm1(-rt/3.5));
}
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
// centroid; (nz, nr) is a wall's unit normal (either orientation); plate marks the injector plate.
struct Neighbour {
    enum Kind { Cell, Axis, Wall, Open, Supply } kind;
    std::size_t cell;
    double dz, dr, nz, nr;
    bool plate{false};
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
        return Neighbour{kind,0,face*m.dz-c.z,kind==Neighbour::Wall?m.ringMiddle(face,j)-c.r:0.0,-1,0,kind==Neighbour::Wall};
    };
    out[0]=i>0?cell(i-1,j):axial(0);
    out[1]=i<m.nz-1?cell(i+1,j):axial(m.nz);
    out[2]=j>0?cell(i,j-1):Neighbour{Neighbour::Axis,0,0,-2*c.r,0,0};
    if(j<m.nr-1) out[3]=cell(i,j+1);
    else {
        // The side wall's value sits where the cell's radial line meets it, at the wall radius at the
        // centroid's z (as the open faces' neighbours sit on the cell's axial line, 1e0fd9e).
        auto ar=m.radialAreaVector(i,m.nr);double area=std::hypot(ar[0],ar[1]);
        double wall=m.wallRadiusAt(i,c.z);
        out[3]={Neighbour::Wall,0,0,wall-c.r,ar[0]/area,ar[1]/area};
    }
    return out;
}
}

void Flow::prepareTransport() {
    const auto count=mesh_.cells.size();
    viscosity_.resize(count);conductivity_.resize(count);
    diffusion_.resize(count*ns_);moles_.resize(count*ns_);
    gradients_.resize(count*gradientFields()*2);leastSquares_.resize(count);
    if(corrections_) {
        correctionFields_.assign(count*4,0);wallDensity_.assign(count,0);wallViscosity_.assign(count,0);
        snOverMu_.assign(count,0);blendF1_.assign(count,0);snFloored_.assign(count,0);correctionCell_.assign(count,0);
        axialCorrection_.assign(axialFlux_.size()*4,0);radialCorrection_.assign(radialFlux_.size()*4,0);
    }
    axialViscous_.resize(axialFlux_.size());axialViscousTransported_.resize(axialFlux_.size()*(ns_+nt_));
    radialViscous_.resize(radialFlux_.size());radialViscousTransported_.resize(radialFlux_.size()*(ns_+nt_));
    const bool chamber=definition_.experiment==Case::Chamber;
    if(nt_) {
        eddy_.assign(count,0);eddyConductivity_.assign(count,0);eddyDiffusion_.assign(count*2,0);
        if(definition_.turbulence.consistentFlux) { staticEnthalpy_.assign(count,0);heatCapacity_.assign(count,0); }
        wallOmega_.assign(count,0);sources_.resize(count);
        // The cells next to a no-slip wall: the side wall's, then the injector plate's not already listed.
        wallCells_.clear();
        if(!definition_.wallSlip) {
            for(int i=0;i<mesh_.nz;++i) wallCells_.push_back(mesh_.index(i,mesh_.nr-1));
            if(chamber) for(int j=0;j<mesh_.nr;++j) if(faceSupply_[j]<0 && j!=mesh_.nr-1) wallCells_.push_back(mesh_.index(0,j));
        }
        wallFaces_.clear();wallCellFace_.clear();plateWallFace_.assign(mesh_.nr,-1);prescribedCell_.assign(count,0);
        prescribedChange_.assign(count,0);
        if(definition_.turbulence.wallFunctions) {
            auto face=[&](int i,int j,int side,bool plate) {
                const auto n=neighbours(mesh_,chamber,faceSupply_,i,j)[side];
                const double nz=plate?-1.0:n.nz,nr=plate?0.0:n.nr;
                wallFaces_.push_back({mesh_.index(i,j),plate,nz,nr,n.dz*nz+n.dr*nr});
            };
            for(int i=0;i<mesh_.nz;++i) face(i,mesh_.nr-1,3,false);
            // Open-area supplies: a supply ring with a wall share has its plate face too, its cell not prescribed.
            if(chamber) for(int j=0;j<mesh_.nr;++j)
                if(faceSupply_[j]<0 || (definition_.openAreaSupplies && plateShare_[j]>0)) { plateWallFace_[j]=static_cast<int>(wallFaces_.size());face(0,j,0,true); }
            for(const auto& f:wallFaces_) if(!(f.distance>0)) throw std::runtime_error("A wall face's normal distance is not positive.");
            for(std::size_t q:wallCells_) {
                const int i=static_cast<int>(q/mesh_.nr),j=static_cast<int>(q%mesh_.nr);
                std::size_t f=j==mesh_.nr-1?static_cast<std::size_t>(i):static_cast<std::size_t>(plateWallFace_[j]);
                if(j==mesh_.nr-1 && i==0 && faceSupply_[j]<0 && plateWallFace_[j]>=0 && wallFaces_[plateWallFace_[j]].distance<wallFaces_[f].distance)
                    f=static_cast<std::size_t>(plateWallFace_[j]);
                wallCellFace_.push_back(f);prescribedCell_[q]=1;
            }
            wallSolutions_.assign(wallFaces_.size(),WallSolution{});
        }
    }
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
void Flow::transportProperties(int depth) {
    const auto& sp=medium_.species();const auto& d=definition_;
    eachCell(region(depth),[&](std::size_t q,int worker) {
        const double* y=fractions_.data()+q*nw_;
        auto t=medium_.transport(temperature_[q],primitives_[q].p,y,diffusion_.data()+q*ns_,transportScratch_[worker].work);
        viscosity_[q]=t.viscosity;conductivity_[q]=t.conductivity;
        double moles=0;
        for(std::size_t k=0;k<ns_;++k) moles+=y[k]/sp[k].molarMass;
        for(std::size_t k=0;k<ns_;++k) moles_[q*ns_+k]=y[k]/sp[k].molarMass/moles;
    });
    if(corrections_) eachCell(region(depth),[&](std::size_t q,int) {
        const double rho=primitives_[q].rho,root=std::sqrt(rho),mu=viscosity_[q];
        const double k=fractions_[q*nw_+ns_],omega=fractions_[q*nw_+ns_+1];
        double* h=&correctionFields_[q*4];
        h[0]=root/mu;h[1]=rho*k;h[2]=root*omega;h[3]=mu*omega;
    });
    if(!nt_ || d.wallSlip || d.turbulence.wallFunctions) return;
    // Wall-face omega = factor 6 nu_w / (beta1 d1^2) of the cells next to a no-slip wall, nu_w at the
    // wall temperature and the cell's pressure and composition (the cell's nu if adiabatic).
    pool_->blocks(wallCells_.size(),kFaces,[&](std::size_t begin,std::size_t end,int worker) {
        auto& s=transportScratch_[worker];
        for(std::size_t n=begin;n<end;++n) {
            const std::size_t q=wallCells_[n];
            if(!inRegion(q,depth)) continue;
            const double* y=fractions_.data()+q*nw_;
            double nu=viscosity_[q]/primitives_[q].rho;
            if(corrections_) { wallDensity_[q]=primitives_[q].rho;wallViscosity_[q]=viscosity_[q]; }
            // The side wall's cells (the corner's too) at wallTemperature, the plate's at its own.
            const double tw=static_cast<int>(q%mesh_.nr)==mesh_.nr-1?d.wallTemperature:plateTemperature();
            if(tw>0) {
                auto t=medium_.transport(tw,primitives_[q].p,y,s.diffusion.data(),s.work);
                nu=t.viscosity*medium_.gasConstant(y)*tw/primitives_[q].p;
                if(corrections_) { wallDensity_[q]=primitives_[q].p/(medium_.gasConstant(y)*tw);wallViscosity_[q]=t.viscosity; }
            }
            wallOmega_[q]=d.turbulence.wallOmegaFactor*6*nu/(beta1*sq(wallDistance_[q]));
        }
    });
}
double Flow::transportValue(std::size_t q,std::size_t f) const {
    if(f>=3+ns_+nt_) return correctionFields_[q*4+f-(3+ns_+nt_)];
    return f==0?primitives_[q].uz:f==1?primitives_[q].ur:f==2?temperature_[q]:f<3+ns_?moles_[q*ns_+f-3]:fractions_[q*nw_+f-3];
}
double Flow::wallValue(std::size_t q,std::size_t f,double nz,double nr,bool plate) const {
    const auto& d=definition_;
    if(f<2) {
        if(!d.wallSlip) return 0.0;
        const auto& w=primitives_[q];double un=w.uz*nz+w.ur*nr;
        return f==0?w.uz-un*nz:w.ur-un*nr;
    }
    if(const double tw=plate?plateTemperature():d.wallTemperature;f==2 && tw>0) return tw;
    // The correction fields at a no-slip wall: psi_w, rho k = 0, sqrt(rho_w) omega_w and mu_w omega_w.
    if(f>=3+ns_+nt_ && !d.wallSlip && !d.turbulence.wallFunctions) {
        const double root=std::sqrt(wallDensity_[q]);
        switch(f-(3+ns_+nt_)) {
            case 0: return root/wallViscosity_[q];
            case 1: return 0.0;
            case 2: return root*wallOmega_[q];
            default: return wallViscosity_[q]*wallOmega_[q];
        }
    }
    if(f>=3+ns_ && !d.wallSlip && !d.turbulence.wallFunctions) return f==3+ns_?0.0:wallOmega_[q];
    return transportValue(q,f);
}
// The field values over N_(depth + 1) and the gradients over N_depth (depth kWhole: every cell).
void Flow::transportGradients(int depth) {
    const auto& m=mesh_;
    const std::size_t nf=gradientFields();
    const bool chamber=definition_.experiment==Case::Chamber;
    // The fields are gathered once per cell into one contiguous array, then each neighbour's
    // differences are taken for all fields at once; every field's sums see the same operations in the
    // same order as one field at a time, so the gradients are unchanged bit for bit.
    fieldValues_.resize(state_.size()*nf);
    eachCell(region(depth<0?depth:depth+1),[&](std::size_t q,int) { for(std::size_t f=0;f<nf;++f) fieldValues_[q*nf+f]=transportValue(q,f); });
    eachCell(region(depth),[&](std::size_t q,int worker) {
        auto& scratch=transportScratch_[worker].gradient;
        scratch.resize(3*nf);
        double *sz=scratch.data(),*sr=sz+nf,*delta=sr+nf;
        const int i=static_cast<int>(q/m.nr),j=static_cast<int>(q%m.nr);
        auto stencil=neighbours(m,chamber,faceSupply_,i,j);const auto& inverse=leastSquares_[q];
        const double* v=&fieldValues_[q*nf];
        std::fill(sz,sz+2*nf,0.0);
        for(const auto& n:stencil) {
            switch(n.kind) {
                case Neighbour::Cell: { const double* u=&fieldValues_[n.cell*nf];for(std::size_t f=0;f<nf;++f) delta[f]=u[f]-v[f];break; }
                case Neighbour::Axis: for(std::size_t f=0;f<nf;++f) delta[f]=f==1?-2*v[f]:0;break;
                case Neighbour::Wall: for(std::size_t f=0;f<nf;++f) delta[f]=wallValue(q,f,n.nz,n.nr,n.plate)-v[f];break;
                case Neighbour::Open: case Neighbour::Supply: std::fill(delta,delta+nf,0.0);break;
            }
            const double d2=n.dz*n.dz+n.dr*n.dr;
            for(std::size_t f=0;f<nf;++f) { double w=delta[f]/d2;sz[f]+=w*n.dz;sr[f]+=w*n.dr; }
        }
        for(std::size_t f=0;f<nf;++f) {
            gradients_[(q*nf+f)*2]=inverse[0]*sz[f]+inverse[1]*sr[f];
            gradients_[(q*nf+f)*2+1]=inverse[1]*sz[f]+inverse[2]*sr[f];
        }
    });
}
void Flow::eddyViscosity(bool withSources,int depth) {
    const auto& tu=definition_.turbulence;
    const std::size_t nf=gradientFields();
    eachCell(region(depth),[&](std::size_t q,int) {
        const auto& w=primitives_[q];const double* g=&gradients_[q*nf*2];const double* y=fractions_.data()+q*nw_;
        const double k=y[ns_],omega=y[ns_+1],hoop=w.ur/mesh_.cells[q].r;
        const double strain2=2*(sq(g[0])+sq(g[3])+sq(hoop))+sq(g[1]+g[2]),strain=std::sqrt(strain2);
        const double* gk=g+2*(3+ns_);const double* gw=g+2*(4+ns_);
        const double cross=gk[0]*gw[0]+gk[1]*gw[1],nu=viscosity_[q]/w.rho,inverse=1/wallDistance_[q];
        const auto b=blending(k,omega,w.rho,nu,inverse,cross);
        double mut=prescribedCell_[q]?w.rho*k/omega:w.rho*a1*k/std::max(a1*omega,strain*b.f2);
        double crossSource=cross;
        if(corrections_) {
            // S_n = 1 / (psi + l n . grad psi), l the wall distance and n the unit vector away from the
            // nearest wall, its denominator floored at psi / 10; the corrected cross-diffusion's product;
            // D^ic on mu_t (not on the wall functions' prescribed mu_t).
            const double* h=&correctionFields_[q*4];const double* gp=g+2*(3+ns_+nt_);
            const auto& n=wallDirection_[q];
            const double den=h[0]+wallDistance_[q]*(n[0]*gp[0]+n[1]*gp[1]);
            snFloored_[q]=den<0.1*h[0];
            snOverMu_[q]=1/(std::max(den,0.1*h[0])*viscosity_[q]);
            blendF1_[q]=b.f1;
            crossSource=(gp[2]*gp[4]+gp[3]*gp[5])/(w.rho*std::sqrt(w.rho));
            if(!prescribedCell_[q]) mut*=intrinsicDamping(k/(nu*omega),std::sqrt(2*k)/sound_[q]);
        }
        auto props=medium_.properties(temperature_[q],y);
        eddy_[q]=mut;eddyConductivity_[q]=(props.cv+props.r)*mut/tu.prandtl;
        if(tu.consistentFlux) { staticEnthalpy_[q]=props.e+props.r*temperature_[q];heatCapacity_[q]=props.cv+props.r; }
        eddyDiffusion_[2*q]=(b.f1*sigmaK1+(1-b.f1)*sigmaK2)*mut;
        eddyDiffusion_[2*q+1]=(b.f1*sigmaW1+(1-b.f1)*sigmaW2)*mut;
        if(withSources) sources_[q]={strain2,strain,cross,nu,inverse,crossSource,sound_[q]};
    });
}
void Flow::mixingInputs(double cmix,const std::vector<std::size_t>& species,std::vector<double>& time,
                        std::vector<double>& segregation) {
    if(!nt_) throw std::logic_error("The PaSR closure needs turbulence.");
    for(std::size_t i:species) if(i>=ns_) throw std::invalid_argument("PaSR species index out of range.");
    refresh(state_,species_,turbulence_,kWhole);
    transportProperties(kWhole);transportGradients(kWhole);eddyViscosity(false,kWhole);
    const auto& tu=definition_.turbulence;
    const std::size_t nf=gradientFields(),count=state_.size();
    time.resize(count);segregation.resize(count);
    pool_->blocks(count,kCells,[&](std::size_t begin,std::size_t end,int) {
      for(std::size_t q=begin;q<end;++q) {
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
    });
}
// MPRK22 on specific k and omega (rho is fixed): with production P and destruction D = d y (d >= 0),
//   y1 = (y0 + tau P0) / (1 + tau d0),  y = (y0 + tau/2 (P0 + P1)) / (1 + tau/2 (d0 y0 + d1 y1) / y1).
// Each source term goes to P or D by its sign; positivity of k and omega is unconditional.
void Flow::turbulenceSource(const std::vector<Conserved>& state,std::vector<double>& turbulence,
                            const std::vector<SourceCoefficients>& coefficients,double tau,int depth) const {
    eachCell(region(depth),[&](std::size_t q,int) {
        if(prescribedCell_[q]) return;
        const auto& c=coefficients[q];const double rho=state[q][0];
        // Production and destruction coefficient of k and omega at (k, omega).
        auto rates=[&](double k,double omega,double& pk,double& dk,double& pw,double& dw) {
            const auto b=blending(k,omega,rho,c.nu,c.inverseDistance,c.crossGradient);
            const double m=std::max(a1*omega,c.strain*b.f2);
            // Pt / (rho k), with mu_t / (rho k) = D a1 / m (D = D^ic with the corrections, else 1), and gamma Pt / nu_t.
            const double damping=corrections_?intrinsicDamping(k/(c.nu*omega),std::sqrt(2*k)/c.sound):1.0;
            const double rate=std::min(damping*a1*c.strain2/m,10*betaStar*omega);
            const double gamma=b.f1*gamma1+(1-b.f1)*gamma2,beta=b.f1*beta1+(1-b.f1)*beta2;
            const double production=gamma*std::min(c.strain2,10*betaStar*omega*m/(a1*damping));
            const double cross=2*(1-b.f1)*sigmaW2*(corrections_?c.crossSource:c.crossGradient);
            pk=k*rate;dk=betaStar*omega;
            pw=production+std::max(cross,0.0)/omega;
            dw=beta*omega+std::max(-cross,0.0)/sq(omega);
        };
        const double k0=turbulence[q*nt_]/rho,w0=turbulence[q*nt_+1]/rho;
        double pk0,dk0,pw0,dw0,pk1,dk1,pw1,dw1;
        rates(k0,w0,pk0,dk0,pw0,dw0);
        const double k1=(k0+tau*pk0)/(1+tau*dk0),w1=(w0+tau*pw0)/(1+tau*dw0);
        rates(k1,w1,pk1,dk1,pw1,dw1);
        const double k=(k0+0.5*tau*(pk0+pk1))/(1+0.5*tau*(k1>0?(dk0*k0+dk1*k1)/k1:dk1));
        const double w=(w0+0.5*tau*(pw0+pw1))/(1+0.5*tau*(dw0*w0+dw1*w1)/w1);
        turbulence[q*nt_]=rho*k;turbulence[q*nt_+1]=rho*w;
    });
}
void Flow::transportFluxes(std::vector<Conserved>& derivative,std::vector<double>& speciesDerivative,
                           std::vector<double>& turbulenceDerivative,BoundaryRates& rates) {
    const auto& m=mesh_;const auto& d=definition_;const auto& sp=medium_.species();
    const std::size_t nf=gradientFields(),nv=ns_+nt_;
    const bool chamber=d.experiment==Case::Chamber;
    // AMR step B2: properties over N_2, gradients and mu_t over N_1, the faces with an active cell on a
    // side and the sums over the active cells.
    transportProperties(2);
    const bool law=nt_ && d.turbulence.wallFunctions,open=chamber && d.openAreaSupplies;
    if(law) pool_->blocks(wallFaces_.size(),kFaces,[&](std::size_t begin,std::size_t end,int worker) {
        for(std::size_t n=begin;n<end;++n) {
            const std::size_t q=wallFaces_[n].cell;
            // Under AMR a wall face outside N_2 has no computed state, and no reader takes its solution.
            if(!inRegion(q,2)) continue;
            wallSolutions_[n]=wallSolve(wallFaces_[n],primitives_[q],temperature_[q],fractions_.data()+q*nw_,viscosity_[q],transportScratch_[worker]);
        }
    });
    transportGradients(1);
    if(nt_) eddyViscosity(false,1);
    const bool restricted=!within_.empty();
    auto value=[&](std::size_t q,std::size_t f) { return fieldValues_[q*nf+f]; };  // gathered by transportGradients above
    auto rhoK=[&](std::size_t q) { return primitives_[q].rho*fractions_[q*nw_+ns_]; };
    // Flux per unit area through a face with unit normal (nz, nr) pointing from cell a to cell b, or
    // out of cell a through a wall (wall != nullptr) whose midpoint is at its displacement. rf is the
    // face radius; axis marks the first face off the axis. Species fluxes go to jn, k and omega fluxes
    // to tn, and with the property corrections the inner and outer k and omega fluxes to hn; s is the
    // worker's scratch. With parts (the budget probe), the flux's parts in the order of Flow::FaceTerm from
    // kViscousWork to kTurbulentTrackedDiffusion, per unit area, then 2 ns_ of scratch.
    const bool consistent=nt_ && d.turbulence.consistentFlux;
    auto faceFlux=[&](std::size_t a,std::size_t b,const Neighbour* wall,double nz,double nr,double rf,bool axis,
                      double* jn,double* tn,double* hn,TransportScratch& s,double* parts=nullptr) {
        const bool boundary=wall!=nullptr;
        double dz=boundary?wall->dz:m.cells[b].z-m.cells[a].z,dr=boundary?wall->dr:m.cells[b].r-m.cells[a].r;
        double len=std::hypot(dz,dr),ez=dz/len,er=dr/len;
        auto other=[&](std::size_t f){ return boundary?wallValue(a,f,wall->nz,wall->nr,wall->plate):value(b,f); };
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
        double normalGradientT=tz*nz+tr*nr,muT=0,eddyLambda=0;
        for(std::size_t t=0;t<nt_;++t) tn[t]=0;
        if(nt_) {
            // A no-slip wall face takes the wall's mu_t = 0 and k = 0; a slip wall face the cell's.
            const bool noSlip=boundary && !d.wallSlip;
            muT=noSlip?0:mean(eddy_[a],boundary?0:eddy_[b]);
            double k=noSlip?0:mean(rhoK(a),boundary?0:rhoK(b));
            tzz+=muT*(2*uzz-2.0/3*div)-2.0/3*k;trr+=muT*(2*urr-2.0/3*div)-2.0/3*k;tzr+=muT*(uzr+urz);
            eddyLambda=noSlip || (consistent && !boundary)?0:mean(eddyConductivity_[a],boundary?0:eddyConductivity_[b]);
            lambda+=eddyLambda;
            if(!(boundary && d.wallSlip)) for(std::size_t t=0;t<nt_;++t) {
                double gz,gr;gradient(3+ns_+t,gz,gr);
                double diffusivity=mu+(noSlip?0:mean(eddyDiffusion_[2*a+t],boundary?0:eddyDiffusion_[2*b+t]));
                tn[t]=-diffusivity*(gz*nz+gr*nr);
            }
            // Corrected diffusion (Hasan et al.): the inner form's flux of rho k with (mu + sigma_k mu_t) S_n / mu
            // and of mu omega with (mu + sigma_w mu_t) S_n / mu, the outer form's of rho k and sqrt(rho) omega
            // with (mu + sigma mu_t) / sqrt(rho); face coefficients are two-cell means of the cell products
            // (a no-slip wall face the cell's, with mu_t = 0), as for the conventional terms.
            if(hn) {
                auto coefficient=[&](std::size_t t,bool inner) {
                    auto one=[&](std::size_t c) {
                        const double diffusivity=viscosity_[c]+(noSlip?0:eddyDiffusion_[2*c+t]);
                        return diffusivity*(inner?snOverMu_[c]:1/std::sqrt(primitives_[c].rho));
                    };
                    return boundary?one(a):0.5*(one(a)+one(b));
                };
                const std::size_t base=3+ns_+nt_;
                double gz,gr;
                gradient(base+1,gz,gr);
                hn[0]=-coefficient(0,true)*(gz*nz+gr*nr);hn[1]=-coefficient(0,false)*(gz*nz+gr*nr);
                gradient(base+3,gz,gr);hn[2]=-coefficient(1,true)*(gz*nz+gr*nr);
                gradient(base+2,gz,gr);hn[3]=-coefficient(1,false)*(gz*nz+gr*nr);
            }
        }
        double fz=tzz*nz+tzr*nr,fr=tzr*nz+trr*nr,heat=-lambda*normalGradientT;
        std::fill(jn,jn+ns_,0.0);
        const std::size_t np=Flow::kTurbulentTrackedDiffusion-Flow::kViscousWork+1;
        if(parts) { std::fill(parts,parts+np+2*ns_,0.0);parts[1]=-(lambda-eddyLambda)*normalGradientT;parts[2]=-eddyLambda*normalGradientT; }
        if(boundary) {
            // A slip wall carries the normal stress only; an adiabatic wall no heat; no wall passes species.
            if(d.wallSlip) { double normal=fz*nz+fr*nr;fz=normal*nz;fr=normal*nr; }
            if(!((wall->plate?plateTemperature():d.wallTemperature)>0)) { heat=0;if(parts) parts[1]=parts[2]=0; }
        } else {
            if(nt_) heat+=tn[0];  // k is part of E
            if(parts && nt_) parts[3]=tn[0];
            // The probe's species parts: the molecular and turbulent fluxes, each with its own share of the correction.
            auto bookSpecies=[&]() {
                double molecular=0,turbulent=0;
                for(std::size_t k=0;k<ns_;++k) { molecular+=parts[np+k];turbulent+=parts[np+ns_+k]; }
                for(std::size_t k=0;k<ns_;++k) {
                    const double jm=parts[np+k]-s.face[k]*molecular,jt=parts[np+ns_+k]-s.face[k]*turbulent;
                    parts[4]+=s.enthalpy[k]*jm;parts[5]+=s.enthalpy[k]*jt;
                    parts[6]+=probe_.weights[k]*jm;parts[7]+=probe_.weights[k]*jt;
                }
            };
            if(consistent) {
                // Definition::Turbulence::consistentFlux: Y_k and h differenced between the two cells along
                // the line of centres (in r^2 across the first face off the axis), the rest of each normal
                // derivative from the mean cell gradients along s = n - (e.n) e, as gradient() composes it;
                // the tangential part of grad Y_k from the X_k gradients at the face composition, and of
                // grad h as cp grad T + sum_k h_k grad Y_k.
                const auto& tu=d.turbulence;
                const double en=ez*nz+er*nr,sz=nz-en*ez,sr=nr-en*er;
                auto along=[&](double va,double vb) {
                    return en*(axis?2*rf*(vb-va)/(m.cells[b].radialSecondMoment-m.cells[a].radialSecondMoment)*er:(vb-va)/len);
                };
                auto across=[&](std::size_t f) {
                    const double* ga=&gradients_[(a*nf+f)*2];const double* gb=&gradients_[(b*nf+f)*2];
                    return 0.5*((ga[0]+gb[0])*sz+(ga[1]+gb[1])*sr);
                };
                double dh=along(staticEnthalpy_[a],staticEnthalpy_[b])+0.5*(heatCapacity_[a]+heatCapacity_[b])*across(2);
                double turbulentSpecies=0,molecularSpecies=0;
                // Definition::unityLewis: every molecular diffusivity lambda / (rho cp), lambda here molecular
                // (eddyLambda is zero on interior faces under consistentFlux).
                const double unity=d.unityLewis?lambda/(0.5*(heatCapacity_[a]+heatCapacity_[b])):0;
                if(ns_>1) {
                    auto& yFace=s.face;auto& jt=s.turbulent;
                    double rho=0.5*(primitives_[a].rho+primitives_[b].rho),moles=0,mixed=0,molecularSum=0,turbulentSum=0;
                    for(std::size_t k=0;k<ns_;++k) { yFace[k]=0.5*(fractions_[a*nw_+k]+fractions_[b*nw_+k]);moles+=yFace[k]/sp[k].molarMass; }
                    for(std::size_t k=0;k<ns_;++k) mixed+=sp[k].molarMass*moles*across(3+k);
                    medium_.speciesEnthalpies(0.5*(temperature_[a]+temperature_[b]),s.enthalpy.data());
                    for(std::size_t k=0;k<ns_;++k) {
                        const double tangential=sp[k].molarMass*moles*across(3+k)-yFace[k]*mixed;
                        const double normalY=along(fractions_[a*nw_+k],fractions_[b*nw_+k])+tangential;
                        if(d.unityLewis) jn[k]=-unity*normalY;
                        else {
                            double gz,gr;gradient(3+k,gz,gr);
                            jn[k]=-rho*sp[k].molarMass*moles*0.5*(diffusion_[a*ns_+k]+diffusion_[b*ns_+k])*(gz*nz+gr*nr);
                        }
                        molecularSum+=jn[k];
                        jt[k]=-muT/tu.schmidt*normalY;
                        turbulentSum+=jt[k];
                        dh+=s.enthalpy[k]*tangential;
                        if(parts) { parts[np+k]=jn[k];parts[np+ns_+k]=jt[k]; }
                    }
                    for(std::size_t k=0;k<ns_;++k) {
                        jn[k]-=yFace[k]*molecularSum;jt[k]-=yFace[k]*turbulentSum;
                        heat+=s.enthalpy[k]*jn[k];molecularSpecies+=s.enthalpy[k]*jn[k];
                        turbulentSpecies+=s.enthalpy[k]*jt[k];jn[k]+=jt[k];
                    }
                    if(parts) bookSpecies();
                }
                if(d.unityLewis) {
                    // The molecular energy flux -(lambda / cp) dh/dn in place of -lambda dT/dn + sum_k h_k j_k.
                    const double molecularHeat=-unity*dh-molecularSpecies;
                    heat+=lambda*normalGradientT+molecularHeat;
                    if(parts) parts[1]=molecularHeat;
                }
                // The turbulent conduction: the part of -(mu_t / Pr_t) dh/dn + (1 - Sc_t / Pr_t) sum h_k j_t
                // that is not species enthalpy.
                const double turbulentHeat=-muT/tu.prandtl*dh-tu.schmidt/tu.prandtl*turbulentSpecies;
                heat+=turbulentHeat+turbulentSpecies;
                if(parts) parts[2]=turbulentHeat;
            } else if(ns_>1) {
                auto& yFace=s.face;
                double rho=0.5*(primitives_[a].rho+primitives_[b].rho),moles=0,sum=0;
                for(std::size_t k=0;k<ns_;++k) { yFace[k]=0.5*(fractions_[a*nw_+k]+fractions_[b*nw_+k]);moles+=yFace[k]/sp[k].molarMass; }
                for(std::size_t k=0;k<ns_;++k) {
                    double gz,gr;gradient(3+k,gz,gr);
                    jn[k]=-rho*sp[k].molarMass*moles*0.5*(diffusion_[a*ns_+k]+diffusion_[b*ns_+k])*(gz*nz+gr*nr);
                    if(parts) parts[np+k]=jn[k];
                    if(nt_) jn[k]-=sp[k].molarMass*moles*muT/d.turbulence.schmidt*(gz*nz+gr*nr);
                    if(parts) parts[np+ns_+k]=jn[k]-parts[np+k];
                    sum+=jn[k];
                }
                medium_.speciesEnthalpies(0.5*(temperature_[a]+temperature_[b]),s.enthalpy.data());
                for(std::size_t k=0;k<ns_;++k) { jn[k]-=yFace[k]*sum;heat+=s.enthalpy[k]*jn[k]; }
                if(parts) bookSpecies();
            }
        }
        if(parts) parts[0]=-(fz*uz+fr*ur);
        return Conserved{0,-fz,-fr,-(fz*uz+fr*ur)+heat};
    };
    // A wall face under the wall functions: the traction and q_w of wall face n, along the face's normal
    // (out of the gas at the side wall, into it along +z at the plate).
    auto lawFlux=[&](std::size_t n,double* jn) {
        std::fill(jn,jn+nv,0.0);
        const auto& w=wallSolutions_[n];const double sign=wallFaces_[n].plate?-1:1;
        return Conserved{0,-sign*w.tz,-sign*w.tr,sign*w.heat};
    };
    // The kind of the end face of row j at axial face i = 0 or nz.
    auto endKind=[&](int i,int j) { return neighbours(m,chamber,faceSupply_,i==0?0:m.nz-1,j)[i==0?0:1].kind; };
    // Axial faces (i, j), i = 0..nz: between cells (i - 1, j) and (i, j), or at the ends the injector
    // plate (a wall), a supply ring (no diffusive flux), the nozzle inlet or the outlet.
    pool_->blocks(restricted?axialFaces_:axialViscous_.size(),kFaces,[&](std::size_t begin,std::size_t end,int worker) {
        auto& s=transportScratch_[worker];
        for(std::size_t e=begin;e<end;++e) {
            const std::size_t f=restricted?faces_[e]:e;
            const int i=static_cast<int>(f/m.nr),j=static_cast<int>(f%m.nr);
            double* jn=&axialViscousTransported_[f*nv];double* tn=jn+ns_;
            double* hn=corrections_?&axialCorrection_[f*4]:nullptr;
            if(hn) std::fill(hn,hn+4,0.0);
            double rf=m.ringMiddle(i,j);
            if(i>0 && i<m.nz) { axialViscous_[f]=faceFlux(m.index(i-1,j),m.index(i,j),nullptr,1,0,rf,false,jn,tn,hn,s);continue; }
            auto q=m.index(i==0?0:m.nz-1,j);
            auto n=neighbours(m,chamber,faceSupply_,i==0?0:m.nz-1,j)[i==0?0:1];
            Conserved flux{};
            if(i==0 && open) {
                // Open-area supplies: w_j of the plate's wall flux; a plate ring's whole, a supply ring's
                // momentum and heat only (its k and omega follow the supply share's zero-gradient face).
                std::fill(jn,jn+nv,0.0);
                if(const double w=plateShare_[j];w>0) {
                    if(law) flux=lawFlux(static_cast<std::size_t>(plateWallFace_[j]),jn);
                    else if(n.kind==Neighbour::Wall) flux=faceFlux(q,0,&n,1,0,rf,false,jn,tn,hn,s);
                    else {
                        const Neighbour wall{Neighbour::Wall,0,n.dz,m.ringMiddle(0,j)-m.cells[q].r,-1,0,true};
                        flux=faceFlux(q,0,&wall,1,0,rf,false,jn,tn,nullptr,s);
                        std::fill(jn,jn+nv,0.0);
                    }
                    for(int k=1;k<4;++k) flux[k]*=w;
                    for(std::size_t k=0;k<nv;++k) jn[k]*=w;
                    if(hn) for(int k=0;k<4;++k) hn[k]*=w;
                }
                axialViscous_[f]=flux;
                continue;
            }
            if(n.kind==Neighbour::Wall) flux=law?lawFlux(static_cast<std::size_t>(plateWallFace_[j]),jn):faceFlux(q,0,&n,1,0,rf,false,jn,tn,hn,s);
            else if(n.kind!=Neighbour::Supply) {
                // Zero normal gradient: tau_zz = -2/3 mu (du_r/dr + u_r/r), tau_zr = mu du_z/dr, no heat,
                // species, k or omega flux; with turbulence mu + mu_t and the normal stress -2/3 rho k.
                const double* g=&gradients_[q*nf*2];
                double mu=viscosity_[q],ur=primitives_[q].ur;
                if(nt_) mu+=eddy_[q];
                double tzz=-2.0/3*mu*(g[3]+ur/m.cells[q].r),tzr=mu*g[1];
                if(nt_) tzz-=2.0/3*rhoK(q);
                flux={0,-tzz,-tzr,-(tzz*primitives_[q].uz+tzr*ur)};
            }
            axialViscous_[f]=flux;
        }
    });
    // Radial faces j = 1..nr of column i: between cells (i, j - 1) and (i, j), or the side wall (j = nr).
    const std::size_t na=axialViscous_.size();
    pool_->blocks(restricted?faces_.size()-axialFaces_:radialViscous_.size(),kFaces,[&](std::size_t begin,std::size_t end,int worker) {
        auto& s=transportScratch_[worker];
        for(std::size_t e=begin;e<end;++e) {
            const std::size_t f=restricted?faces_[axialFaces_+e]-na:e;
            const int i=static_cast<int>(f/(m.nr+1)),j=static_cast<int>(f%(m.nr+1));
            if(j==0) continue;
            double* jn=&radialViscousTransported_[f*nv];double* tn=jn+ns_;
            double* hn=corrections_?&radialCorrection_[f*4]:nullptr;
            if(hn) std::fill(hn,hn+4,0.0);
            auto ar=m.radialAreaVector(i,j);double area=std::hypot(ar[0],ar[1]);
            double nz=ar[0]/area,nr=ar[1]/area,rf=m.radialFaceMiddle(i,j);
            auto a=m.index(i,j-1);
            if(j<m.nr) { radialViscous_[f]=faceFlux(a,m.index(i,j),nullptr,nz,nr,rf,j==1,jn,tn,hn,s);continue; }
            auto n=neighbours(m,chamber,faceSupply_,i,j-1)[3];
            radialViscous_[f]=law?lawFlux(static_cast<std::size_t>(i),jn):faceFlux(a,0,&n,nz,nr,rf,false,jn,tn,hn,s);
        }
    });
    // Each cell sums its faces: the + side of an interior face gains area times the flux, the - side
    // loses it; an end face enters with sign +1 at i = 0 and -1 at i = nz.
    eachCell(region(0),[&](std::size_t q,int) {
      const int i=static_cast<int>(q/m.nr),j=static_cast<int>(q%m.nr);
      const std::size_t a0=q,a1=q+m.nr,r0=static_cast<std::size_t>(i)*(m.nr+1)+j,r1=r0+1;
      auto gain=[&](double area,const Conserved& flux,const double* t) {
          for(int k=1;k<4;++k) derivative[q][k]+=area*flux[k];
          for(std::size_t k=0;k<ns_;++k) speciesDerivative[q*ns_+k]+=area*t[k];
          for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]+=area*t[ns_+k];
      };
      auto lose=[&](double area,const Conserved& flux,const double* t) {
          for(int k=1;k<4;++k) derivative[q][k]-=area*flux[k];
          for(std::size_t k=0;k<ns_;++k) speciesDerivative[q*ns_+k]-=area*t[k];
          for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]-=area*t[ns_+k];
      };
      auto endFace=[&](int face,std::size_t f) {
          const auto kind=endKind(face,j);
          if(kind==Neighbour::Supply) return;
          const double sign=face==0?1:-1,area=axialArea_[f];const auto& flux=axialViscous_[f];
          for(int k=1;k<4;++k) derivative[q][k]+=sign*area*flux[k];
          if(kind==Neighbour::Wall) for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]+=sign*area*axialViscousTransported_[f*nv+ns_+k];
      };
      if(i>0) gain(axialArea_[a0],axialViscous_[a0],&axialViscousTransported_[a0*nv]); else endFace(0,a0);
      if(i<m.nz-1) lose(axialArea_[a1],axialViscous_[a1],&axialViscousTransported_[a1*nv]); else endFace(m.nz,a1);
      if(j>0) gain(radialArea_[r0],radialViscous_[r0],&radialViscousTransported_[r0*nv]);
      if(j<m.nr-1) lose(radialArea_[r1],radialViscous_[r1],&radialViscousTransported_[r1*nv]);
      else {  // the side wall
          const double area=radialArea_[r1];const auto& flux=radialViscous_[r1];
          for(int k=1;k<4;++k) derivative[q][k]-=area*flux[k];
          for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]-=area*radialViscousTransported_[r1*nv+ns_+k];
      }
      // Hoop stress: tau_thth = (mu + mu_t)(2 u_r / r - 2/3 div u) - 2/3 rho k at the centroid.
      const auto& c=m.cells[q];
      double ur=primitives_[q].ur,div=gradients_[q*nf*2]+gradients_[(q*nf+1)*2+1]+ur/c.r;
      double hoop=viscosity_[q]*(2*ur/c.r-2.0/3*div);
      if(nt_) hoop+=eddy_[q]*(2*ur/c.r-2.0/3*div)-2.0/3*rhoK(q);
      derivative[q][2]-=hoop*c.radialPressureMeasure;
      // Property corrections: the corrected k and omega diffusion, F1 (S_n / mu) sum(inner k) + (1 - F1)
      // sum(outer k) / sqrt(rho) and F1 (rho / mu)(S_n / mu) sum(inner omega) + (1 - F1) sum(outer omega),
      // replace the conventional face sums (the same faces, signs and areas); the difference for k,
      // Phi_k, also enters E (k is part of E), its integral booked in correctionEnergy. The wall
      // functions' prescribed cells have no k or omega equation.
      if(corrections_) {
          correctionCell_[q]=0;
          if(!prescribedCell_[q]) {
              double conventional[2]={0,0},corrected[4]={0,0,0,0};
              auto add=[&](double signedArea,const double* t,const double* h) {
                  for(int k=0;k<2;++k) conventional[k]+=signedArea*t[ns_+k];
                  for(int k=0;k<4;++k) corrected[k]+=signedArea*h[k];
              };
              auto addEnd=[&](int face,std::size_t f) {
                  if(endKind(face,j)!=Neighbour::Wall) return;
                  add((face==0?1:-1)*axialArea_[f],&axialViscousTransported_[f*nv],&axialCorrection_[f*4]);
              };
              if(i>0) add(axialArea_[a0],&axialViscousTransported_[a0*nv],&axialCorrection_[a0*4]); else addEnd(0,a0);
              if(i<m.nz-1) add(-axialArea_[a1],&axialViscousTransported_[a1*nv],&axialCorrection_[a1*4]); else addEnd(m.nz,a1);
              if(j>0) add(radialArea_[r0],&radialViscousTransported_[r0*nv],&radialCorrection_[r0*4]);
              add(-radialArea_[r1],&radialViscousTransported_[r1*nv],&radialCorrection_[r1*4]);
              const double rho=primitives_[q].rho,f1=blendF1_[q],sn=snOverMu_[q];
              const double k=f1*sn*corrected[0]+(1-f1)*corrected[1]/std::sqrt(rho);
              const double omega=f1*rho/viscosity_[q]*sn*corrected[2]+(1-f1)*corrected[3];
              const double phiK=k-conventional[0];
              turbulenceDerivative[q*nt_]+=phiK;turbulenceDerivative[q*nt_+1]+=omega-conventional[1];
              derivative[q][3]+=phiK;
              correctionCell_[q]=phiK;
          }
      }
    });
    // The budget probe's viscous terms (Flow::setBudgetProbe): each probe cell's faces again, serially,
    // with their parts kept apart; the fluxes the cells used above are read, not changed.
    if(!probe_.cells.empty()) {
        auto& s=transportScratch_[0];
        const std::size_t np=kTurbulentTrackedDiffusion-kViscousWork+1;
        std::vector<double> jn(nv),hn(4),parts(np+2*ns_);
        double* hp=corrections_?hn.data():nullptr;
        auto trackedSum=[&](const double* t) { double x=0;for(std::size_t k=0;k<ns_;++k) x+=probe_.weights[k]*t[k];return x; };
        for(std::size_t c=0;c<probe_.cells.size();++c) {
            const std::size_t q=probe_.cells[c];
            const int i=static_cast<int>(q/m.nr),j=static_cast<int>(q%m.nr);
            const std::size_t a0=q,a1=q+m.nr;
            double* out=probeStage_.data()+c*kProbeStride;
            auto book=[&](int face,double sign,double area,const Conserved& used,const double* moved) {
                double* o=out+face*kFaceTerms;
                for(std::size_t p=0;p<np;++p) o[kViscousWork+p]+=sign*area*parts[p];
                o[kViscousEnergy]+=sign*area*used[3];
                if(moved) o[kTrackedTransport]+=sign*area*trackedSum(moved);
            };
            auto end=[&](int face,std::size_t f) {
                const auto kind=endKind(face,j);
                std::fill(parts.begin(),parts.end(),0.0);
                if(kind==Neighbour::Supply) return;
                const auto& flux=axialViscous_[f];
                if(kind==Neighbour::Wall) {
                    if(law) parts[1]=flux[3];
                    else { auto n=neighbours(m,chamber,faceSupply_,face==0?0:m.nz-1,j)[face==0?0:1];
                           faceFlux(q,0,&n,1,0,m.ringMiddle(face,j),false,jn.data(),jn.data()+ns_,hp,s,parts.data()); }
                } else parts[0]=flux[3];
                book(face==0?0:1,face==0?1:-1,axialArea_[f],flux,nullptr);
            };
            if(i>0) { faceFlux(m.index(i-1,j),q,nullptr,1,0,m.ringMiddle(i,j),false,jn.data(),jn.data()+ns_,hp,s,parts.data());
                      book(0,1,axialArea_[a0],axialViscous_[a0],&axialViscousTransported_[a0*nv]); }
            else end(0,a0);
            if(i<m.nz-1) { faceFlux(q,m.index(i+1,j),nullptr,1,0,m.ringMiddle(i+1,j),false,jn.data(),jn.data()+ns_,hp,s,parts.data());
                           book(1,-1,axialArea_[a1],axialViscous_[a1],&axialViscousTransported_[a1*nv]); }
            else end(m.nz,a1);
            auto radial=[&](int face,int jf,double sign) {
                const std::size_t f=static_cast<std::size_t>(i)*(m.nr+1)+jf;
                auto ar=m.radialAreaVector(i,jf);double area=std::hypot(ar[0],ar[1]);
                double nz=ar[0]/area,nr=ar[1]/area,rf=m.radialFaceMiddle(i,jf);
                const auto lo=m.index(i,jf-1);
                std::fill(parts.begin(),parts.end(),0.0);
                if(jf<m.nr) { faceFlux(lo,m.index(i,jf),nullptr,nz,nr,rf,jf==1,jn.data(),jn.data()+ns_,hp,s,parts.data());
                              book(face,sign,radialArea_[f],radialViscous_[f],&radialViscousTransported_[f*nv]);return; }
                if(law) parts[1]=radialViscous_[f][3];
                else { auto n=neighbours(m,chamber,faceSupply_,i,jf-1)[3];
                       faceFlux(lo,0,&n,nz,nr,rf,false,jn.data(),jn.data()+ns_,hp,s,parts.data()); }
                book(face,sign,radialArea_[f],radialViscous_[f],nullptr);
            };
            if(j>0) radial(2,j,1);
            radial(3,j+1,-1);
            if(corrections_) out[4*kFaceTerms+kCorrection]+=correctionCell_[q];
        }
    }
    if(corrections_) for(std::size_t q=0;q<correctionCell_.size();++q) if(isOwned(q)) rates.correction+=correctionCell_[q];
    // Boundary rates in the order of the serial face loops: the end faces at i = 0, then at i = nz,
    // then the side wall by column; under AMR only beside owned cells.
    for(int i:{0,m.nz}) for(int j=0;j<m.nr;++j) {
        const auto kind=endKind(i,j);
        const bool share=i==0 && open;  // open-area supplies: the face's viscous flux is its wall share's
        if((share?!(plateShare_[j]>0):kind==Neighbour::Supply) || !isOwned(m.index(i==0?0:m.nz-1,j))) continue;
        const auto f=static_cast<std::size_t>(i)*m.nr+j;
        const double sign=i==0?1:-1,area=axialArea_[f];const auto& flux=axialViscous_[f];
        rates.energy+=sign*area*flux[3];
        if(share || kind==Neighbour::Wall) { rates.wallAxial+=area*flux[1];rates.wallHeat+=area*flux[3];if(i==0) rates.plateHeat+=area*flux[3]; }
        else if(i==0) rates.inletMomentum+=area*flux[1];
        else rates.outletMomentum+=area*flux[1];
    }
    for(int i=0;i<m.nz;++i) {
        if(!isOwned(m.index(i,m.nr-1))) continue;
        const auto f=static_cast<std::size_t>(i)*(m.nr+1)+m.nr;
        const double area=radialArea_[f];const auto& flux=radialViscous_[f];
        rates.wallAxial-=area*flux[1];rates.energy-=area*flux[3];rates.wallHeat-=area*flux[3];
    }
}
Flow::WallSolution Flow::wallSolve(const WallFace& f,const Primitive& w,double t1,const double* y,double mu1,
                                   TransportScratch& s) const {
    const auto& d=definition_;const auto& tu=d.turbulence;
    const double fixed=f.plate?plateTemperature():d.wallTemperature;
    const double normal=w.uz*f.nz+w.ur*f.nr,uz=w.uz-normal*f.nz,ur=w.ur-normal*f.nr,u1=std::hypot(uz,ur);
    // The layer at wall temperature tw: mu_w and k_w, cp_w, and the layer's mean cp (below 1e-3 T_w apart
    // the difference of enthalpies would round, and cp at the mean temperature is used).
    Medium::Transport wall{};double cpW=0,rhoW=0,cp=0;
    auto layer=[&](double tw) {
        wall=medium_.transport(tw,w.p,y,s.diffusion.data(),s.work);
        const auto atWall=medium_.properties(tw,y);
        cpW=atWall.cv+atWall.r;rhoW=w.p/(atWall.r*tw);
        if(std::abs(t1-tw)>1e-3*tw) cp=(medium_.enthalpy(t1,y)-medium_.enthalpy(tw,y))/(t1-tw);
        else { const auto m=medium_.properties(0.5*(t1+tw),y);cp=m.cv+m.r; }
    };
    // An adiabatic wall (temperature 0): beta = 0, so Crocco-Busemann at the first cell puts T_w at the
    // recovery temperature T_1 + r u_1^2 / (2 cp), r = Pr_w^(1/3) (WALL_FUNCTIONS.md, eq. 10.8), solved by
    // fixed point with the layer's properties at T_w; the isothermal law at that T_w then has beta = 0
    // and q_w = 0 to rounding.
    double tw=fixed;
    if(!(fixed>0)) {
        tw=t1;
        for(int n=0;;++n) {
            layer(tw);
            const double next=t1+std::cbrt(cpW*wall.viscosity/wall.conductivity)*u1*u1/(2*cp);
            const bool done=std::abs(next-tw)<=1e-13*next;
            tw=next;
            if(done) break;
            if(n==50) throw std::runtime_error("Wall function: the adiabatic wall's recovery temperature did not converge.");
        }
    }
    if(t1>11*tw) throw std::runtime_error("Wall function: T_1 / T_w above 11, where the law folds back.");
    layer(tw);
    const wallLaw::Constants c{tu.wallKappa,tu.wallB};
    const auto law=wallLaw::solveIsothermal({u1,f.distance,t1,tw,rhoW,wall.viscosity,wall.conductivity,cp,
                                             std::cbrt(cpW*wall.viscosity/wall.conductivity)},c);
    WallSolution out;out.temperature=tw;
    if(u1>0) { const double scale=law.shear/u1;out.tz=-scale*uz;out.tr=-scale*ur; }
    out.heat=law.heat;out.yPlus=law.yPlus;
    const double ratio=mu1/wall.viscosity;
    const double eddy=tu.printedDerivative?wallLaw::yPlusDerivativePrinted(law.uPlus,law.gamma,law.beta,c)-ratio
                                          :wallLaw::eddyRatio(law.uPlus,law.gamma,law.beta,ratio,c);
    out.eddy=wall.viscosity*std::max(eddy,0.0);
    const double inner=6*wall.viscosity/(beta1*rhoW*sq(f.distance)),outer=law.uTau/(std::sqrt(betaStar)*c.kappa*f.distance);
    out.omega=std::hypot(inner,outer);
    return out;
}
std::vector<Flow::WallYPlus> Flow::wallYPlus() const {
    std::vector<WallYPlus> out;out.reserve(wallFaces_.size());
    for(std::size_t n=0;n<wallFaces_.size();++n) out.push_back({wallFaces_[n].cell,wallFaces_[n].plate,wallSolutions_[n].yPlus,wallSolutions_[n].heat,wallSolutions_[n].temperature});
    return out;
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
    refresh(state_,species_,turbulence_,kWhole);
    BoundaryRates rates{};
    transportFluxes(derivative,speciesDerivative,turbulenceDerivative,rates);
    for(std::size_t q=0;q<state_.size();++q) {
        for(double& v:derivative[q]) v/=mesh_.cells[q].volume;
        for(std::size_t k=0;k<ns_;++k) speciesDerivative[q*ns_+k]/=mesh_.cells[q].volume;
        for(std::size_t k=0;k<nt_;++k) turbulenceDerivative[q*nt_+k]/=mesh_.cells[q].volume;
    }
}
} // namespace crucible
