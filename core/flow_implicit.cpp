// The implicit prototype (docs/evidence/IMPLICIT_PROTOTYPE.md, its criteria written before this file): BDF2 in
// physical time, solved by dual time with LU-SGS sweeps. Kept apart from flow.cpp, so the explicit march's code
// and arithmetic are not touched.
#include "core/flow.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace crucible {
namespace {
double sq(double x) { return x*x; }
// The change of the convective flux through the area vector (sz, sr) of a cell at (w, a, E, y) when its conserved
// state changes by du (bulk, then ns partial densities): the Jacobian's product at fixed composition, the pressure's
// change from dp/d(rho e) = gamma - 1 and dp/d(rho) = p / rho - (gamma - 1) e, gamma = rho a^2 / p the frozen one.
// The sweeps' approximation only; the residual keeps the scheme's flux and the full thermodynamics.
void fluxChange(const Primitive& w, double a, double energy, const double* y, std::size_t ns, const double* du,
                double sz, double sr, double* out) {
    const double un=w.uz*sz+w.ur*sr, gamma=w.rho*a*a/w.p, kinetic=0.5*(sq(w.uz)+sq(w.ur));
    const double e=energy/w.rho-kinetic;
    const double dm=du[1]*sz+du[2]*sr, dun=(dm-un*du[0])/w.rho;
    const double dp=(gamma-1)*(du[3]-w.uz*du[1]-w.ur*du[2]+kinetic*du[0])+(w.p/w.rho-(gamma-1)*e)*du[0];
    out[0]=dm;
    out[1]=du[1]*un+w.rho*w.uz*dun+dp*sz;
    out[2]=du[2]*un+w.rho*w.ur*dun+dp*sr;
    out[3]=(du[3]+dp)*un+(energy+w.p)*dun;
    for(std::size_t k=0;k<ns;++k) out[4+k]=du[4+k]*un+w.rho*y[k]*dun;
}
}

// LU-SGS (Yoon and Jameson) on the Rusanov split: a face's flux changes by (J_S(U_q) + l) dU_q / 2 from its own
// cell q and (J_S(U_p) - l) dU_p / 2 from the other, p, with S its area vector out of q and l its spectral radius
// times its area, the larger of the two cells' |u . S| + a |S|. The own-cell Jacobians are dropped (their sum over
// a cell's faces vanishes but for the axisymmetric pressure term), leaving the scalar diagonal
// V / (a dt) + (1/2 + 1/pseudoCfl) sum(l); the off-diagonal products are evaluated at the iterate's cell states.
// On entry implicitDelta_ holds the right-hand side -V d / (a dt); on return the update.
void Flow::implicitSweeps(double diagonalTime,double pseudoCfl) {
    const auto& m=mesh_;const int nz=m.nz,nr=m.nr;const std::size_t nv=4+ns_;
    auto rate=[&](std::size_t q,double sz,double sr) {
        const auto& w=primitives_[q];
        return std::abs(w.uz*sz+w.ur*sr)+sound_[q]*std::hypot(sz,sr);
    };
    for(int i=0;i<=nz;++i) for(int j=0;j<nr;++j) {
        const double area=m.axialArea(i,j);double r=0;
        if(i>0) r=rate(m.index(i-1,j),area,0);
        if(i<nz) r=std::max(r,rate(m.index(i,j),area,0));
        implicitAxialRate_[static_cast<std::size_t>(i)*nr+j]=r;
    }
    for(int i=0;i<nz;++i) for(int j=0;j<=nr;++j) {
        const auto s=m.radialAreaVector(i,j);double r=0;
        if(j>0) r=rate(m.index(i,j-1),s[0],s[1]);
        if(j<nr) r=std::max(r,rate(m.index(i,j),s[0],s[1]));
        implicitRadialRate_[static_cast<std::size_t>(i)*(nr+1)+j]=r;
    }
    auto axial=[&](int i,int j) { return implicitAxialRate_[static_cast<std::size_t>(i)*nr+j]; };
    auto radial=[&](int i,int j) { return implicitRadialRate_[static_cast<std::size_t>(i)*(nr+1)+j]; };
    for(int i=0;i<nz;++i) for(int j=0;j<nr;++j) {
        const std::size_t q=m.index(i,j);
        implicitDiagonal_[q]=m.cells[q].volume/diagonalTime+(0.5+1/pseudoCfl)*(axial(i,j)+axial(i+1,j)+radial(i,j)+radial(i,j+1));
    }
    std::vector<double> change(nv),upper(nv);
    // acc += sign (J_S(U_p) du - l du) / 2.
    auto offDiagonal=[&](std::size_t p,double sz,double sr,double l,double sign,double* acc) {
        const double* du=implicitDelta_.data()+p*nv;
        fluxChange(primitives_[p],sound_[p],stage_[p][3],fractions_.data()+p*nw_,ns_,du,sz,sr,change.data());
        for(std::size_t k=0;k<nv;++k) acc[k]+=sign*0.5*(change[k]-l*du[k]);
    };
    // Forward: (D + L) dU* = b, the lower neighbours (i - 1, j) and (i, j - 1) already swept.
    for(int i=0;i<nz;++i) for(int j=0;j<nr;++j) {
        const std::size_t q=m.index(i,j);double* d=implicitDelta_.data()+q*nv;
        if(i>0) offDiagonal(m.index(i-1,j),-m.axialArea(i,j),0,axial(i,j),-1,d);
        if(j>0) { const auto s=m.radialAreaVector(i,j);offDiagonal(m.index(i,j-1),-s[0],-s[1],radial(i,j),-1,d); }
        for(std::size_t k=0;k<nv;++k) d[k]/=implicitDiagonal_[q];
    }
    // Backward: dU = dU* - D^-1 U dU, the upper neighbours (i + 1, j) and (i, j + 1) already final.
    for(int i=nz-1;i>=0;--i) for(int j=nr-1;j>=0;--j) {
        const std::size_t q=m.index(i,j);double* d=implicitDelta_.data()+q*nv;
        std::fill(upper.begin(),upper.end(),0.0);
        if(i<nz-1) offDiagonal(m.index(i+1,j),m.axialArea(i+1,j),0,axial(i+1,j),1,upper.data());
        if(j<nr-1) { const auto s=m.radialAreaVector(i,j+1);offDiagonal(m.index(i,j+1),s[0],s[1],radial(i,j+1),1,upper.data()); }
        for(std::size_t k=0;k<nv;++k) d[k]-=upper[k]/implicitDiagonal_[q];
    }
}

double Flow::implicitStep(double dt,const ImplicitControl& control) {
    if(!(dt>0)) throw std::invalid_argument("Step interval must be positive.");
    if(nt_ || medium_.hasTransport() || !active_.empty() || !within_.empty())
        throw std::runtime_error("The implicit prototype runs one inviscid level without turbulence.");
    const auto& m=mesh_;const std::size_t n=state_.size(),nv=4+ns_;
    const bool history=implicitPreviousDt_>0 && implicitEndTime_==time_ && implicitEndSteps_==steps_;
    const std::vector<Conserved> start=state_;
    const std::vector<double> startSpecies=species_;
    implicitKnown_.resize(n*nv);implicitDelta_.resize(n*nv);implicitDiagonal_.resize(n);
    implicitAxialRate_.resize(static_cast<std::size_t>(m.nz+1)*m.nr);
    implicitRadialRate_.resize(static_cast<std::size_t>(m.nz)*(m.nr+1));
    double total=0;
    for(const auto& c:m.cells) total+=c.volume;
    std::vector<double> trial(ns_);
    for(int attempt=0;attempt<14;++attempt,dt*=0.5) {
        // BDF2 on the steps h_(n-1), h = dt, w = h / h_(n-1): U - (1 + c) U^n + c U^(n-1) = a h R(U), a = (1 + w) / (1 + 2 w),
        // c = w^2 / (1 + 2 w) (a = 2/3, c = 1/3 at w = 1); backward Euler (a = 1, c = 0) without a history.
        double a=1,c=0;
        if(history) { const double w=dt/implicitPreviousDt_;a=(1+w)/(1+2*w);c=w*w/(1+2*w); }
        for(std::size_t q=0;q<n;++q) {
            double* known=implicitKnown_.data()+q*nv;
            // As U^n + c (U^n - U^(n-1)): (1 + c) U^n - c U^(n-1) rounds with a bias (about -6e-17 of U at c = 1/3)
            // that does not vanish when U^n = U^(n-1), and the budgets drift by it every step.
            for(int k=0;k<4;++k) known[k]=start[q][k]+(history?c*(start[q][k]-implicitPrevious_[q][k]):0.0);
            for(std::size_t k=0;k<ns_;++k)
                known[4+k]=startSpecies[q*ns_+k]+(history?c*(startSpecies[q*ns_+k]-implicitPreviousSpecies_[q*ns_+k]):0.0);
        }
        // The igniter's mean power over [t, t + h], as the explicit step's.
        const auto& ig=definition_.igniter;
        igniterPower_=0;
        if(ig.energy>0) {
            const double overlap=std::min(time_+dt,ig.start+ig.duration)-std::max(time_,ig.start);
            if(overlap>0) igniterPower_=ig.energy/ig.duration*overlap/dt;
        }
        stage_=start;speciesStage_=startSpecies;
        if(history && control.extrapolate) {
            const double w=dt/implicitPreviousDt_;
            for(std::size_t q=0;q<n;++q) {
                Conserved u{};
                for(int k=0;k<4;++k) u[k]=start[q][k]+w*(start[q][k]-implicitPrevious_[q][k]);
                for(std::size_t k=0;k<ns_;++k) trial[k]=startSpecies[q*ns_+k]+w*(startSpecies[q*ns_+k]-implicitPreviousSpecies_[q*ns_+k]);
                if(!admissible(u,trial.data(),nullptr)) continue;
                stage_[q]=u;std::copy(trial.begin(),trial.end(),speciesStage_.begin()+static_cast<std::ptrdiff_t>(q*ns_));
            }
        }
        BoundaryRates rates;
        bool failed=false;
        int it=0;
        double defect=0,first=0;
        std::array<double,4> scale{};
        try {
            for(;;++it) {
                if(control.project) {
                    state_.swap(stage_);species_.swap(speciesStage_);
                    try { control.project(); } catch(...) { state_.swap(stage_);species_.swap(speciesStage_);throw; }
                    state_.swap(stage_);species_.swap(speciesStage_);
                }
                rates=rhs(stage_,speciesStage_,turbulence_,rhs_,speciesRhs_,turbulenceRhs_,time_+dt);
                if(it==0) {
                    double r2=0,ra2=0,e2=0;
                    for(std::size_t q=0;q<n;++q) {
                        const double v=m.cells[q].volume;
                        r2+=v*sq(stage_[q][0]);ra2+=v*sq(stage_[q][0]*sound_[q]);e2+=v*sq(stage_[q][3]);
                    }
                    scale={std::sqrt(r2/total),std::sqrt(ra2/total),std::sqrt(ra2/total),std::sqrt(e2/total)};
                }
                std::array<double,4> d2{};
                for(std::size_t q=0;q<n;++q) {
                    const double* known=implicitKnown_.data()+q*nv;double* b=implicitDelta_.data()+q*nv;
                    const double v=m.cells[q].volume;
                    for(int k=0;k<4;++k) {
                        const double d=stage_[q][k]-known[k]-a*dt*rhs_[q][k];
                        d2[k]+=v*d*d;b[k]=-v/(a*dt)*d;
                    }
                    for(std::size_t k=0;k<ns_;++k)
                        b[4+k]=-v/(a*dt)*(speciesStage_[q*ns_+k]-known[4+k]-a*dt*speciesRhs_[q*ns_+k]);
                }
                defect=0;
                for(int k=0;k<4;++k) defect=std::max(defect,std::sqrt(d2[k]/total)/scale[k]);
                if(it==0) first=defect;
                if(control.trace) std::printf("  implicit t %.6e attempt %d iteration %d defect %.4e\n",time_+dt,attempt,it,defect);
                if(!std::isfinite(defect)) { failed=true;break; }
                if(defect<=control.tolerance || it==control.maxIterations) break;
                implicitSweeps(a*dt,control.pseudoCfl);
                // Each cell's update, damped until the cell is admissible (density at most halved).
                for(std::size_t q=0;q<n;++q) {
                    const double* du=implicitDelta_.data()+q*nv;
                    double theta=1;
                    if(du[0]<-0.5*stage_[q][0]) theta=0.5*stage_[q][0]/(-du[0]);
                    Conserved u{};
                    for(int halve=0;;++halve) {
                        for(int k=0;k<4;++k) u[k]=stage_[q][k]+theta*du[k];
                        for(std::size_t k=0;k<ns_;++k) trial[k]=speciesStage_[q*ns_+k]+theta*du[4+k];
                        if(admissible(u,trial.data(),nullptr)) break;
                        if(halve==20) { theta=0;break; }
                        theta*=0.5;
                    }
                    if(theta>0) { stage_[q]=u;std::copy(trial.begin(),trial.end(),speciesStage_.begin()+static_cast<std::ptrdiff_t>(q*ns_)); }
                }
            }
        } catch(const std::runtime_error&) { failed=true; }
        if(!failed) {
            // The step from the last iterate's residual: U^(n+1) = U_known + a h R(U^(m)).
            bool ok=true;
            for(std::size_t q=0;q<n;++q) {
                const double* known=implicitKnown_.data()+q*nv;
                for(int k=0;k<4;++k) next_[q][k]=known[k]+a*dt*rhs_[q][k];
                for(std::size_t k=0;k<ns_;++k) speciesNext_[q*ns_+k]=known[4+k]+a*dt*speciesRhs_[q*ns_+k];
            }
            // With the projection inside, the final update's transported composition (unreacted over the step) is
            // projected before it is judged, as the iterates were: at the reacted gas's energy it can fall below the floor.
            if(control.project) {
                state_.swap(next_);species_.swap(speciesNext_);
                try { control.project(); } catch(const std::runtime_error&) { ok=false; }
                state_.swap(next_);species_.swap(speciesNext_);
            }
            for(std::size_t q=0;q<n && ok;++q)
                if(!admissible(next_[q],speciesNext_.data()+q*ns_,nullptr)) ok=false;
            if(ok) {
                implicitPrevious_=start;implicitPreviousSpecies_=startSpecies;implicitPreviousDt_=dt;
                state_.swap(next_);species_.swap(speciesNext_);
                time_+=dt;dt_=dt;++steps_;
                lastRates_=rates;
                // The budgets follow the state's recurrence: I(n+1) - I(n) = c (I(n) - I(n-1)) + a h rate, so
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
                s.lastDefect=defect;s.lastFirstDefect=first;
                if(defect>control.tolerance) ++s.unconverged;
                return dt;
            }
        }
        ++rejectedSteps_;++implicitStats_.rejected;implicitStats_.rejectedIterations+=static_cast<std::uint64_t>(it);
    }
    throw std::runtime_error("Could not advance an admissible gas state after 14 implicit step reductions.");
}
}
