// The implicit prototype (docs/evidence/IMPLICIT_PROTOTYPE.md, its criteria written before this file): BDF2 in
// physical time, solved by dual time with LU-SGS sweeps. Kept apart from flow.cpp, so the explicit march's code
// and arithmetic are not touched.
#include "core/flow.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>

namespace crucible {
namespace {
double sq(double x) { return x*x; }
// The change of the convective flux through the area vector (sz, sr) of a cell at (w, a, E, y) when its conserved
// state changes by du (bulk, then ns partial densities, then nt turbulence variables rho k, rho omega): the Jacobian's
// product at fixed composition, the pressure's change from dp/d(rho e) = gamma - 1 and dp/d(rho) = p / rho - (gamma - 1) e,
// gamma = rho a^2 / p the frozen one, rho e = E - rho |u|^2 / 2 - rho k. The sweeps' approximation only; the residual keeps
// the scheme's flux and the full thermodynamics.
void fluxChange(const Primitive& w, double a, double energy, const double* y, std::size_t ns, const double* turbulence,
                std::size_t nt, const double* du, double sz, double sr, double* out) {
    const double un=w.uz*sz+w.ur*sr, gamma=w.rho*a*a/w.p, kinetic=0.5*(sq(w.uz)+sq(w.ur));
    double e=energy/w.rho-kinetic;
    if(nt) e-=turbulence[0]/w.rho;
    const double dm=du[1]*sz+du[2]*sr, dun=(dm-un*du[0])/w.rho;
    double dp=(gamma-1)*(du[3]-w.uz*du[1]-w.ur*du[2]+kinetic*du[0])+(w.p/w.rho-(gamma-1)*e)*du[0];
    if(nt) dp-=(gamma-1)*du[4+ns];
    out[0]=dm;
    out[1]=du[1]*un+w.rho*w.uz*dun+dp*sz;
    out[2]=du[2]*un+w.rho*w.ur*dun+dp*sr;
    out[3]=(du[3]+dp)*un+(energy+w.p)*dun;
    for(std::size_t k=0;k<ns;++k) out[4+k]=du[4+k]*un+w.rho*y[k]*dun;
    for(std::size_t k=0;k<nt;++k) out[4+ns+k]=du[4+ns+k]*un+turbulence[k]*dun;
}
}

// LU-SGS (Yoon and Jameson) on the Rusanov split: a face's flux changes by (J_S(U_q) + l) dU_q / 2 from its own
// cell q and (J_S(U_p) - l) dU_p / 2 from the other, p, with S its area vector out of q and l its spectral radius
// times its area, the larger of the two cells' |u . S| + a |S|. The own-cell Jacobians are dropped (their sum over
// a cell's faces vanishes but for the axisymmetric pressure term), leaving the scalar diagonal
// V / (a dt) + (1/2 + 1/pseudoCfl) sum(l); the off-diagonal products are evaluated at the iterate's cell states.
// With transport each face also carries the viscous flux's thin-layer estimate l_v (U_q - U_p), l_v = nu |S|^2 / V
// with nu the larger cell's largest diffusivity as Flow::stableDt takes it and V the two cells' mean volume (the
// cell's own at a boundary face): l_v on the diagonal, 2 l_v / pseudoCfl more on it for the pseudo-time step (stableDt's
// diffusive rate), and -l_v dU_p off it.
// On entry implicitDelta_ holds the right-hand side -V d / (a dt); on return the update.
void Flow::implicitSweeps(double diagonalTime,double pseudoCfl) {
    const auto& m=mesh_;const int nz=m.nz,nr=m.nr;const std::size_t nv=4+ns_+nt_;
    const bool viscous=medium_.hasTransport();
    auto rate=[&](std::size_t q,double sz,double sr) {
        const auto& w=primitives_[q];
        return std::abs(w.uz*sz+w.ur*sr)+sound_[q]*std::hypot(sz,sr);
    };
    if(viscous) {
        implicitNu_.resize(m.cells.size());
        for(std::size_t q=0;q<m.cells.size();++q) {
            const double* y=fractions_.data()+q*nw_;const double rho=primitives_[q].rho;
            double nu=std::max(4.0/3*viscosity_[q],conductivity_[q]/medium_.cv(temperature_[q],y))/rho;
            for(std::size_t k=0;k<ns_;++k) nu=std::max(nu,diffusion_[q*ns_+k]);
            if(nt_) {
                const auto& tu=definition_.turbulence;
                nu=std::max({nu,4.0/3*(viscosity_[q]+eddy_[q])/rho,(conductivity_[q]+eddyConductivity_[q])/(medium_.cv(temperature_[q],y)*rho)});
                for(std::size_t k=0;k<ns_;++k) nu=std::max(nu,diffusion_[q*ns_+k]+eddy_[q]/(rho*tu.schmidt));
            }
            implicitNu_[q]=nu;
        }
        implicitAxialViscous_.resize(implicitAxialRate_.size());implicitRadialViscous_.resize(implicitRadialRate_.size());
    }
    auto viscousRate=[&](std::ptrdiff_t p,std::ptrdiff_t q,double area2) {
        double nu=0,volume=0;int n=0;
        for(std::ptrdiff_t c:{p,q}) if(c>=0) { nu=std::max(nu,implicitNu_[static_cast<std::size_t>(c)]);volume+=m.cells[static_cast<std::size_t>(c)].volume;++n; }
        return nu*area2*n/volume;
    };
    for(int i=0;i<=nz;++i) for(int j=0;j<nr;++j) {
        const double area=m.axialArea(i,j);double r=0;
        if(i>0) r=rate(m.index(i-1,j),area,0);
        if(i<nz) r=std::max(r,rate(m.index(i,j),area,0));
        implicitAxialRate_[static_cast<std::size_t>(i)*nr+j]=r;
        if(viscous) implicitAxialViscous_[static_cast<std::size_t>(i)*nr+j]=viscousRate(i>0?static_cast<std::ptrdiff_t>(m.index(i-1,j)):-1,
                                                                                     i<nz?static_cast<std::ptrdiff_t>(m.index(i,j)):-1,sq(area));
    }
    for(int i=0;i<nz;++i) for(int j=0;j<=nr;++j) {
        const auto s=m.radialAreaVector(i,j);double r=0;
        if(j>0) r=rate(m.index(i,j-1),s[0],s[1]);
        if(j<nr) r=std::max(r,rate(m.index(i,j),s[0],s[1]));
        implicitRadialRate_[static_cast<std::size_t>(i)*(nr+1)+j]=r;
        if(viscous) implicitRadialViscous_[static_cast<std::size_t>(i)*(nr+1)+j]=viscousRate(j>0?static_cast<std::ptrdiff_t>(m.index(i,j-1)):-1,
                                                                                         j<nr?static_cast<std::ptrdiff_t>(m.index(i,j)):-1,sq(s[0])+sq(s[1]));
    }
    auto axial=[&](int i,int j) { return implicitAxialRate_[static_cast<std::size_t>(i)*nr+j]; };
    auto radial=[&](int i,int j) { return implicitRadialRate_[static_cast<std::size_t>(i)*(nr+1)+j]; };
    auto axialViscous=[&](int i,int j) { return viscous?implicitAxialViscous_[static_cast<std::size_t>(i)*nr+j]:0.0; };
    auto radialViscous=[&](int i,int j) { return viscous?implicitRadialViscous_[static_cast<std::size_t>(i)*(nr+1)+j]:0.0; };
    for(int i=0;i<nz;++i) for(int j=0;j<nr;++j) {
        const std::size_t q=m.index(i,j);
        implicitDiagonal_[q]=m.cells[q].volume/diagonalTime+(0.5+1/pseudoCfl)*(axial(i,j)+axial(i+1,j)+radial(i,j)+radial(i,j+1));
        if(viscous) implicitDiagonal_[q]+=(1+2/pseudoCfl)*(axialViscous(i,j)+axialViscous(i+1,j)+radialViscous(i,j)+radialViscous(i,j+1));
    }
    std::vector<double> change(nv),upper(nv);
    // acc += sign ((J_S(U_p) du - l du) / 2 - l_v du).
    auto offDiagonal=[&](std::size_t p,double sz,double sr,double l,double lv,double sign,double* acc) {
        const double* du=implicitDelta_.data()+p*nv;
        fluxChange(primitives_[p],sound_[p],stage_[p][3],fractions_.data()+p*nw_,ns_,turbulenceStage_.data()+p*nt_,nt_,du,sz,sr,change.data());
        if(!viscous) { for(std::size_t k=0;k<nv;++k) acc[k]+=sign*0.5*(change[k]-l*du[k]); }
        else for(std::size_t k=0;k<nv;++k) acc[k]+=sign*(0.5*(change[k]-l*du[k])-lv*du[k]);
    };
    // Forward: (D + L) dU* = b, the lower neighbours (i - 1, j) and (i, j - 1) already swept.
    for(int i=0;i<nz;++i) for(int j=0;j<nr;++j) {
        const std::size_t q=m.index(i,j);double* d=implicitDelta_.data()+q*nv;
        if(i>0) offDiagonal(m.index(i-1,j),-m.axialArea(i,j),0,axial(i,j),axialViscous(i,j),-1,d);
        if(j>0) { const auto s=m.radialAreaVector(i,j);offDiagonal(m.index(i,j-1),-s[0],-s[1],radial(i,j),radialViscous(i,j),-1,d); }
        for(std::size_t k=0;k<nv;++k) d[k]/=implicitDiagonal_[q];
    }
    // Backward: dU = dU* - D^-1 U dU, the upper neighbours (i + 1, j) and (i, j + 1) already final.
    for(int i=nz-1;i>=0;--i) for(int j=nr-1;j>=0;--j) {
        const std::size_t q=m.index(i,j);double* d=implicitDelta_.data()+q*nv;
        std::fill(upper.begin(),upper.end(),0.0);
        if(i<nz-1) offDiagonal(m.index(i+1,j),m.axialArea(i+1,j),0,axial(i+1,j),axialViscous(i+1,j),1,upper.data());
        if(j<nr-1) { const auto s=m.radialAreaVector(i,j+1);offDiagonal(m.index(i,j+1),s[0],s[1],radial(i,j+1),radialViscous(i,j+1),1,upper.data()); }
        for(std::size_t k=0;k<nv;++k) d[k]-=upper[k]/implicitDiagonal_[q];
    }
}

double Flow::implicitStep(double dt,const ImplicitControl& control) {
    if(!(dt>0)) throw std::invalid_argument("Step interval must be positive.");
    if(!active_.empty() || !within_.empty()) throw std::runtime_error("The implicit step runs one level without regions.");
    const auto& m=mesh_;const std::size_t n=state_.size(),nv=4+ns_+nt_;
    const bool viscous=medium_.hasTransport(),walls=nt_ && definition_.turbulence.wallFunctions;
    const bool history=implicitPreviousDt_>0 && implicitEndTime_==time_ && implicitEndSteps_==steps_;
    // The explicit step's stable interval, for the record, and with turbulence the source coefficients at the start,
    // which its refresh computes, as the explicit step takes them. Without transport or turbulence nothing is
    // refreshed, as before.
    implicitStats_.lastExplicitDt=(nt_ || viscous)?stableDt():0.0;
    beginStep();
    const std::vector<Conserved> start=state_;
    const std::vector<double> startSpecies=species_;
    implicitKnown_.resize(n*nv);implicitDelta_.resize(n*nv);implicitDiagonal_.resize(n);
    implicitAxialRate_.resize(static_cast<std::size_t>(m.nz+1)*m.nr);
    implicitRadialRate_.resize(static_cast<std::size_t>(m.nz)*(m.nr+1));
    double total=0;
    for(const auto& c:m.cells) total+=c.volume;
    std::vector<double> trial(ns_),trialTurbulence(nt_);
    const auto increment=[&](std::size_t q,std::size_t k) { return implicitPreviousIncrement_[q*nv+k]; };
    const auto prescribed=[&](std::size_t q) { return walls && prescribedCell_[q]; };
    // A closure's projection with the state (bulk, species, turbulence) swapped in.
    const auto project=[&](std::vector<Conserved>& u,std::vector<double>& y,std::vector<double>& t) {
        state_.swap(u);species_.swap(y);turbulence_.swap(t);
        try { control.project(); } catch(...) { state_.swap(u);species_.swap(y);turbulence_.swap(t);throw; }
        state_.swap(u);species_.swap(y);turbulence_.swap(t);
    };
    // Why each attempt failed, for the error: the stage and the first cell it failed at.
    std::string firstWhy,lastWhy;
    for(int attempt=0;attempt<14;++attempt) {
        double retry=0.5*dt;  // the next attempt's step: half, or the time error's choice
        // BDF2 on the steps h_(n-1), h = dt, w = h / h_(n-1): U - U^n - c (increment) = a h R(U), a = (1 + w) / (1 + 2 w),
        // c = w^2 / (1 + 2 w) (a = 2/3, c = 1/3 at w = 1); backward Euler (a = 1, c = 0) without a history.
        double a=1,c=0;
        if(history) { const double w=dt/implicitPreviousDt_;a=(1+w)/(1+2*w);c=w*w/(1+2*w); }
        Attempt booked{};
        BoundaryRates rates;
        bool failed=false;
        std::string why;
        const auto fail=[&](const char* stage,std::size_t q) { why=std::string(stage)+(q<n?" at cell "+std::to_string(q):""); };
        // A cell's density, thermal energy over the floor, rho k and rho omega, for the error.
        const auto show=[&](const char* label,const Conserved& u,const double* y,const double* t) {
            double sum=0;
            for(std::size_t k=0;k<ns_;++k) sum+=std::max(y[k],0.0);
            const double floor=sum>0?medium_.energyFloorOfClipped(y)/sum:0.0;
            char text[256];
            std::snprintf(text,sizeof text," %s: rho %.6e, e - floor - k %.6e, rho k %.6e, rho omega %.6e;",label,u[0],
                          (u[3]-(sq(u[1])+sq(u[2]))/(2*u[0]))/u[0]-floor-(nt_?t[0]/u[0]:0.0),nt_?t[0]:0.0,nt_?t[1]:0.0);
            why+=text;
        };
        int it=0;
        double defect=0,first=0,turbulenceDefect=0;
        // The first half source, from the start's rho k and rho omega with the start's coefficients.
        if(nt_) {
            turbulence_=turbulenceStart_;
            turbulenceSource(state_,turbulence_,sourcesStart_,0.5*dt,ghostWidth()-1);
            for(std::size_t q=0;q<n && !failed;++q)
                if((failed=!admissible(state_[q],species_.data()+q*ns_,turbulence_.data()+q*nt_))) {
                    fail("the first half source",q);
                    show("after",state_[q],species_.data()+q*ns_,turbulence_.data()+q*nt_);
                    show("start",state_[q],species_.data()+q*ns_,turbulenceStart_.data()+q*nt_);
                    char text[64];
                    std::snprintf(text,sizeof text," T %.6g K, dt %.6e s",temperature_[q],dt);
                    why+=text;
                }
        }
        if(!failed) for(std::size_t q=0;q<n;++q) {
            double* known=implicitKnown_.data()+q*nv;
            // As U^n + c increment: (1 + c) U^n - c U^(n-1) rounds with a bias (about -6e-17 of U at c = 1/3) that does
            // not vanish when U^n = U^(n-1), and the budgets drift by it every step.
            for(int k=0;k<4;++k) known[k]=start[q][k]+(history?c*increment(q,k):0.0);
            for(std::size_t k=0;k<ns_;++k) known[4+k]=startSpecies[q*ns_+k]+(history?c*increment(q,4+k):0.0);
            for(std::size_t k=0;k<nt_;++k) known[4+ns_+k]=turbulence_[q*nt_+k]+(history?c*increment(q,4+ns_+k):0.0);
            // Where rho omega fell more than (1 + c) / c times over the last step, the history alone asks for rho omega
            // <= 0, the damping then halves it every iteration, and the cross-diffusion's destruction (|CD| / omega,
            // singular at omega = 0) takes the rest: 5.4 to 1.7e-118 in one step, then a NaN in the next first half
            // source at every step size (the run i20u at 8.975 ms, SPEED_PLAN.md 7.4). Half the step shrinks c at least
            // twofold.
            if(nt_ && !failed && !prescribed(q) && !(known[4+ns_+1]>0)) {
                failed=true;fail("the history's rho omega at or below 0",q);
                char text[96];
                std::snprintf(text,sizeof text," (%.6e from %.6e, c %.4g)",known[4+ns_+1],turbulence_[q*nt_+1],c);
                why+=text;
            }
        }
        // The igniter's mean power over [t, t + h], as the explicit step's.
        const auto& ig=definition_.igniter;
        igniterPower_=0;
        if(ig.energy>0) {
            const double overlap=std::min(time_+dt,ig.start+ig.duration)-std::max(time_,ig.start);
            if(overlap>0) igniterPower_=ig.energy/ig.duration*overlap/dt;
        }
        stage_=start;speciesStage_=startSpecies;
        if(nt_) turbulenceStage_=turbulence_;
        if(!failed && history && control.extrapolate) {
            const double w=dt/implicitPreviousDt_;
            for(std::size_t q=0;q<n;++q) {
                Conserved u{};
                for(int k=0;k<4;++k) u[k]=start[q][k]+w*increment(q,k);
                for(std::size_t k=0;k<ns_;++k) trial[k]=startSpecies[q*ns_+k]+w*increment(q,4+k);
                for(std::size_t k=0;k<nt_;++k) trialTurbulence[k]=turbulence_[q*nt_+k]+w*increment(q,4+ns_+k);
                if(!admissible(u,trial.data(),nt_?trialTurbulence.data():nullptr)) continue;
                stage_[q]=u;std::copy(trial.begin(),trial.end(),speciesStage_.begin()+static_cast<std::ptrdiff_t>(q*ns_));
                std::copy(trialTurbulence.begin(),trialTurbulence.end(),turbulenceStage_.begin()+static_cast<std::ptrdiff_t>(q*nt_));
            }
        }
        std::array<double,4> scale{};
        std::array<double,2> turbulenceScale{};
        if(!failed) try {
            for(;;++it) {
                if(control.project) project(stage_,speciesStage_,turbulenceStage_);
                if(walls) { Attempt scratch;if(!prescribe(stage_,speciesStage_,turbulenceStage_,scratch)) { failed=true;fail("the wall functions on an iterate",n);break; } }
                rates=rhs(stage_,speciesStage_,turbulenceStage_,rhs_,speciesRhs_,turbulenceRhs_,time_+dt);
                if(it==0) {
                    double r2=0,ra2=0,e2=0;std::array<double,2> t2{};
                    for(std::size_t q=0;q<n;++q) {
                        const double v=m.cells[q].volume;
                        r2+=v*sq(stage_[q][0]);ra2+=v*sq(stage_[q][0]*sound_[q]);e2+=v*sq(stage_[q][3]);
                        // Over the cells the transport sets: a wall cell's rho omega, the wall functions', is orders above.
                        if(!prescribed(q)) for(std::size_t k=0;k<nt_;++k) t2[k]+=v*sq(turbulence_[q*nt_+k]);
                    }
                    scale={std::sqrt(r2/total),std::sqrt(ra2/total),std::sqrt(ra2/total),std::sqrt(e2/total)};
                    for(std::size_t k=0;k<nt_;++k) turbulenceScale[k]=std::sqrt(t2[k]/total);
                }
                std::array<double,4> d2{};std::array<double,2> dt2{};
                for(std::size_t q=0;q<n;++q) {
                    const double* known=implicitKnown_.data()+q*nv;double* b=implicitDelta_.data()+q*nv;
                    const double v=m.cells[q].volume;
                    for(int k=0;k<4;++k) {
                        const double d=stage_[q][k]-known[k]-a*dt*rhs_[q][k];
                        d2[k]+=v*d*d;b[k]=-v/(a*dt)*d;
                    }
                    for(std::size_t k=0;k<ns_;++k)
                        b[4+k]=-v/(a*dt)*(speciesStage_[q*ns_+k]-known[4+k]-a*dt*speciesRhs_[q*ns_+k]);
                    // A prescribed cell's rho k and rho omega are the wall functions', not the transport's.
                    for(std::size_t k=0;k<nt_;++k) {
                        if(prescribed(q)) { b[4+ns_+k]=0;continue; }
                        const double d=turbulenceStage_[q*nt_+k]-known[4+ns_+k]-a*dt*turbulenceRhs_[q*nt_+k];
                        dt2[k]+=v*d*d;b[4+ns_+k]=-v/(a*dt)*d;
                    }
                }
                defect=0;
                for(int k=0;k<4;++k) defect=std::max(defect,std::sqrt(d2[k]/total)/scale[k]);
                turbulenceDefect=0;
                for(std::size_t k=0;k<nt_;++k) {
                    const double d=std::sqrt(dt2[k]/total);
                    turbulenceDefect=std::max(turbulenceDefect,turbulenceScale[k]>0?d/turbulenceScale[k]:d);
                }
                if(it==0) first=defect;
                if(control.trace) {
                    if(nt_) std::printf("  implicit t %.6e attempt %d iteration %d defect %.4e turbulence %.4e\n",time_+dt,attempt,it,defect,turbulenceDefect);
                    else std::printf("  implicit t %.6e attempt %d iteration %d defect %.4e\n",time_+dt,attempt,it,defect);
                }
                if(!std::isfinite(defect) || !std::isfinite(turbulenceDefect)) { failed=true;fail("a non-finite defect",n);break; }
                if(std::max(defect,turbulenceDefect)<=control.tolerance || it==control.maxIterations) break;
                implicitSweeps(a*dt,control.pseudoCfl);
                // Each cell's update, damped until the cell is admissible (density and rho omega at most halved, rho k
                // floored at 0 in the iterate; the final state's clip is booked).
                for(std::size_t q=0;q<n;++q) {
                    const double* du=implicitDelta_.data()+q*nv;
                    const bool held=prescribed(q);
                    double theta=1;
                    if(du[0]<-0.5*stage_[q][0]) theta=0.5*stage_[q][0]/(-du[0]);
                    if(nt_ && !held) {
                        const double omega=turbulenceStage_[q*nt_+1],dw=du[4+ns_+1];
                        if(dw<-0.5*omega) theta=std::min(theta,0.5*omega/(-dw));
                    }
                    Conserved u{};
                    for(int halve=0;;++halve) {
                        for(int k=0;k<4;++k) u[k]=stage_[q][k]+theta*du[k];
                        for(std::size_t k=0;k<ns_;++k) trial[k]=speciesStage_[q*ns_+k]+theta*du[4+k];
                        for(std::size_t k=0;k<nt_;++k) trialTurbulence[k]=turbulenceStage_[q*nt_+k]+(held?0.0:theta*du[4+ns_+k]);
                        if(nt_) trialTurbulence[0]=std::max(trialTurbulence[0],0.0);
                        if(admissible(u,trial.data(),nt_?trialTurbulence.data():nullptr)) break;
                        if(halve==20) { theta=0;break; }
                        theta*=0.5;
                    }
                    if(theta>0) {
                        stage_[q]=u;std::copy(trial.begin(),trial.end(),speciesStage_.begin()+static_cast<std::ptrdiff_t>(q*ns_));
                        std::copy(trialTurbulence.begin(),trialTurbulence.end(),turbulenceStage_.begin()+static_cast<std::ptrdiff_t>(q*nt_));
                    }
                }
            }
        } catch(const std::runtime_error& e) { failed=true;why=std::string("iteration ")+std::to_string(it)+": "+e.what(); }
        // With the time error on, a step the inner iterations leave over their tolerance is redone at half: the time
        // error does not see rho k or rho omega. One accepted at 400 iterations took a thin nozzle cell's rho k from
        // 5.6 to 0.11 and its rho omega from 31 to 650, where the explicit run holds 5.5 and 29, and omega's collapse
        // followed 9 steps later (the run i20u at 8.89 ms, SPEED_PLAN.md 7.4).
        if(!failed && control.timeTolerance>0 && std::max(defect,turbulenceDefect)>control.tolerance) {
            failed=true;++implicitStats_.unconverged;
            char text[128];
            std::snprintf(text,sizeof text,"no convergence in %d iterations (defect %.3e, turbulence %.3e)",it,defect,turbulenceDefect);
            why=text;
        }
        if(!failed) {
            // The step from the last iterate's residual: U^(n+1) = U_known + a h R(U^(m)), and the increment U^(n+1) - U^n,
            // kept as the next step's history, before anything split after it: rho k and rho omega before the clip, the
            // bulk and species after the projection inside (part of the step's equations; it keeps rho, E and the
            // elements, so it moves the composition only along the equilibrium manifold).
            bool ok=true;
            for(std::size_t q=0;q<n;++q) {
                const double* known=implicitKnown_.data()+q*nv;double* step=implicitDelta_.data()+q*nv;
                for(int k=0;k<4;++k) next_[q][k]=known[k]+a*dt*rhs_[q][k];
                for(std::size_t k=0;k<ns_;++k) speciesNext_[q*ns_+k]=known[4+k]+a*dt*speciesRhs_[q*ns_+k];
                // rho k and rho omega are the last iterate's, admissible by its damping: they carry no budget, and the
                // update form's one further Jacobi step can take rho omega below 0 where the iterate is near tolerance
                // (the agreement run i5, SPEED_PLAN.md 7.3).
                for(std::size_t k=0;k<nt_;++k) {
                    turbulenceNext_[q*nt_+k]=turbulenceStage_[q*nt_+k];
                    step[4+ns_+k]=turbulenceNext_[q*nt_+k]-turbulence_[q*nt_+k];
                }
                clipTurbulence(turbulenceNext_,q);
            }
            sumClipped(booked);
            // With the projection inside, the final update's transported composition (unreacted over the step) is
            // projected before it is judged, as the iterates were: at the reacted gas's energy it can fall below the floor.
            try {
                if(control.project) project(next_,speciesNext_,turbulenceNext_);
                for(std::size_t q=0;q<n;++q) {
                    double* step=implicitDelta_.data()+q*nv;
                    for(int k=0;k<4;++k) step[k]=next_[q][k]-start[q][k];
                    for(std::size_t k=0;k<ns_;++k) step[4+k]=speciesNext_[q*ns_+k]-startSpecies[q*ns_+k];
                }
                for(std::size_t q=0;q<n && ok;++q)
                    if(!admissible(next_[q],speciesNext_.data()+q*ns_,nt_?turbulenceNext_.data()+q*nt_:nullptr)) {
                        ok=false;fail("the final update",q);
                        // What failed, with the start's and the last iterate's alongside.
                        show("final",next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_);
                        show("iterate",stage_[q],speciesStage_.data()+q*ns_,turbulenceStage_.data()+q*nt_);
                        show("start",start[q],startSpecies.data()+q*ns_,turbulence_.data()+q*nt_);
                        why+=" dt "+std::to_string(dt)+", iterations "+std::to_string(it)+", defect "+std::to_string(defect);
                    }
                if(ok && walls && !(ok=prescribe(next_,speciesNext_,turbulenceNext_,booked))) fail("the wall functions on the final state",n);
                // The second half source, with the coefficients at the new state (the explicit second stage's).
                if(ok && nt_) {
                    refresh(next_,speciesNext_,turbulenceNext_,2);
                    transportProperties(2);transportGradients(1);eddyViscosity(true,1);
                    turbulenceSource(next_,turbulenceNext_,sources_,0.5*dt,1);
                    for(std::size_t q=0;q<n && ok;++q)
                        if(!(ok=admissible(next_[q],speciesNext_.data()+q*ns_,turbulenceNext_.data()+q*nt_))) fail("the second half source",q);
                }
            } catch(const std::runtime_error& e) { ok=false;why=std::string("the final state: ")+e.what(); }
            // The time error: the step's own increment against the predictor U^n + w (the last increment), per cell
            // over the start's own scales (rho; |m| with rho times the thermal energy; the internal energy above the
            // floor, k in it), as an RMS over the cells, each cell alike (a front in the thin nozzle gas counts as much
            // as the chamber). The predictor's error is O(h^2), so the next step is h (tolerance / error)^(1/2) times
            // 0.9, within [0.2, 2] h, and a step over the tolerance is redone at that.
            double timeError=0,suggested=dt;
            if(ok && control.timeTolerance>0 && history) {
                const double w=dt/implicitPreviousDt_;
                double sum=0;
                for(std::size_t q=0;q<n;++q) {
                    const Conserved& u0=start[q];
                    const double* step=implicitDelta_.data()+q*nv;
                    double partial=0;
                    for(std::size_t k=0;k<ns_;++k) partial+=std::max(startSpecies[q*ns_+k],0.0);
                    const double floor=partial>0?medium_.energyFloorOfClipped(startSpecies.data()+q*ns_)/partial:0.0;
                    const auto thermal=[&](const Conserved& u) { return u[3]-(sq(u[1])+sq(u[2]))/(2*u[0])-u[0]*floor; };
                    Conserved p{},e{};
                    for(int k=0;k<4;++k) { p[k]=u0[k]+w*increment(q,k);e[k]=u0[k]+step[k]; }
                    const double t0=thermal(u0);
                    const double er=(e[0]-p[0])/u0[0];
                    const double em2=(sq(e[1]-p[1])+sq(e[2]-p[2]))/(sq(u0[1])+sq(u0[2])+u0[0]*t0);
                    const double et=(p[0]>0?thermal(e)-thermal(p):std::numeric_limits<double>::infinity())/t0;
                    sum+=er*er+em2+et*et;
                }
                timeError=std::sqrt(sum/(3.0*static_cast<double>(n)));
                if(!std::isfinite(timeError)) timeError=std::numeric_limits<double>::infinity();
                suggested=dt*std::clamp(0.9*std::sqrt(control.timeTolerance/std::max(timeError,1e-300)),0.2,2.0);
                if(timeError>control.timeTolerance) {
                    ok=false;retry=suggested;++implicitStats_.timeRejected;
                    why="the time error "+std::to_string(timeError)+" over "+std::to_string(control.timeTolerance);
                }
            }
            if(ok) {
                implicitPreviousIncrement_.swap(implicitDelta_);implicitPreviousDt_=dt;
                state_.swap(next_);species_.swap(speciesNext_);
                if(nt_) turbulence_.swap(turbulenceNext_);
                time_+=dt;dt_=dt;++steps_;
                lastRates_=rates;
                integratedRates_.accumulate(rates,dt);
                if(nt_) { clippedTurbulentEnergy_+=booked.clipped;prescribedTurbulentEnergy_+=booked.prescribed; }
                // The budgets follow the state's recurrence: I(n+1) - I(n) = c (the last increment) + a h rate, so
                // mass - I stays where it was (the projection after the step keeps rho, E and the elements).
                const std::array<double,6> r={rates.heat,rates.mass,rates.energy,rates.correction,rates.netMomentum(),rates.grossMomentum()};
                for(std::size_t k=0;k<6;++k) implicitIncrement_[k]=(history?c*implicitIncrement_[k]:0.0)+a*dt*r[k];
                auto add=[](double& plain,double& compensation,double term) {
                    const double sum=plain+term;
                    compensation+=std::abs(plain)>=std::abs(term)?(plain-sum)+term:(term-sum)+plain;
                    plain=sum;
                };
                add(integratedHeat_,heatCompensation_,implicitIncrement_[0]);
                add(integratedMassFlux_,massCompensation_,implicitIncrement_[1]);
                add(integratedEnergyFlux_,energyCompensation_,implicitIncrement_[2]);
                add(integratedCorrection_,correctionCompensation_,implicitIncrement_[3]);
                add(integratedMomentumSource_,momentumCompensation_,implicitIncrement_[4]);
                integratedMomentumGross_+=implicitIncrement_[5];
                implicitEndTime_=time_;implicitEndSteps_=steps_;
                auto& s=implicitStats_;
                ++s.steps;s.iterations+=static_cast<std::uint64_t>(it);s.lastIterations=it;s.mostIterations=std::max(s.mostIterations,it);
                s.lastDefect=defect;s.lastFirstDefect=first;s.lastTurbulenceDefect=turbulenceDefect;s.lastTimeError=timeError;
                implicitSuggestedDt_=control.timeTolerance>0 && history?suggested:0.0;
                if(std::max(defect,turbulenceDefect)>control.tolerance) ++s.unconverged;
                return dt;
            }
        }
        ++rejectedSteps_;++implicitStats_.rejected;implicitStats_.rejectedIterations+=static_cast<std::uint64_t>(it);
        if(control.trace) std::printf("  implicit t %.6e attempt %d failed: %s\n",time_+dt,attempt,why.c_str());
        if(attempt==0) firstWhy=why;
        lastWhy=why;
        dt=retry;
    }
    throw std::runtime_error("Could not advance an admissible gas state after 14 implicit step reductions (the first: "+firstWhy+
                             "; the last: "+lastWhy+").");
}
}
