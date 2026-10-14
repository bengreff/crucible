#include "core/amr.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace crucible {

Amr::Amr(const Definition& base,int finer) {
    if(finer<1 || finer>3) throw std::invalid_argument("AMR takes 1 to 3 levels above the base.");
    // The finest mesh from the definition, each coarser one the exact union of its children.
    Definition finest=base;finest.nz=base.nz<<finer;finest.nr=base.nr<<finer;
    std::vector<Mesh> meshes;
    meshes.push_back(Mesh(finest));
    for(int l=finer-1;l>=0;--l) meshes.push_back(Mesh::coarsened(meshes.back()));
    std::reverse(meshes.begin(),meshes.end());
    for(int l=0;l<=finer;++l) {
        Definition d=base;d.nz=base.nz<<l;d.nr=base.nr<<l;
        Level v;v.flow=std::make_unique<Flow>(d,std::move(meshes[static_cast<std::size_t>(l)]));
        levels_.push_back(std::move(v));
    }
    const Flow& b=level(0);
    if(b.medium_.hasTransport() || b.nt_)
        throw std::invalid_argument("AMR step B1 is inviscid: transport and turbulence are not yet refluxed.");
    for(int l=0;l<=finer;++l) {
        Flow& f=level(l);
        f.recordFaces_=true;
        f.faceIntegral_.assign((f.axialFlux_.size()+f.radialFlux_.size())*(4+f.nw_),0.0);
        if(l==0) continue;
        levels_[static_cast<std::size_t>(l)].fineSum.assign(f.faceIntegral_.size(),0.0);
        f.active_.assign(f.mesh_.cells.size(),0);
        // A fine injector ring takes its base parent ring's supply, so the supply areas are the base's;
        // a fine cell over a base igniter cell takes the base's power density.
        if(base.experiment==Case::Chamber) {
            std::fill(f.supplyArea_.begin(),f.supplyArea_.end(),0.0);
            for(int j=0;j<f.mesh_.nr;++j) {
                const int s=f.faceSupply_[static_cast<std::size_t>(j)]=b.faceSupply_[static_cast<std::size_t>(j>>l)];
                if(s>=0) f.supplyArea_[static_cast<std::size_t>(s)]+=f.mesh_.axialArea(0,j);
            }
            if(base.igniter.energy>0) {
                f.igniterCells_.clear();
                for(int i=0;i<f.mesh_.nz;++i) for(int j=0;j<f.mesh_.nr;++j) {
                    const auto q=f.mesh_.index(i,j);
                    f.igniterCell_[q]=b.igniterCell_[b.mesh_.index(i>>l,j>>l)];
                    if(f.igniterCell_[q]) f.igniterCells_.push_back(q);
                }
                f.igniterVolume_=b.igniterVolume_;
            }
        }
    }
    time_=b.time_;
}
void Amr::setThreads(int threads) {
    for(auto& v:levels_) v.flow->setThreads(threads);
}
void Amr::setActive(int l,std::vector<unsigned char> active) {
    if(l<1 || l>=levels()) throw std::invalid_argument("Only the levels above the base take an active region.");
    Flow& f=level(l);
    if(active.size()!=f.mesh_.cells.size()) throw std::invalid_argument("The active region needs one entry per cell of its level.");
    for(auto& a:active) a=a?1:0;
    f.active_=std::move(active);
    checked_=false;
}
bool Amr::covered(int l,std::size_t q) const {
    if(l+1>=levels()) return false;
    const Flow& c=level(l);const Flow& f=level(l+1);
    const int i=static_cast<int>(q/static_cast<std::size_t>(c.mesh_.nr)),j=static_cast<int>(q%static_cast<std::size_t>(c.mesh_.nr));
    return f.isActive(f.mesh_.index(2*i,2*j));
}
// Whole parent cells, and proper nesting: every cell within the ghost width of an active cell of
// level l (clipped to the domain) lies over an active cell of level l - 1.
void Amr::check() {
    top_=0;
    for(int l=1;l<levels();++l) {
        const Flow& f=level(l);const Flow& c=level(l-1);
        const Mesh& m=f.mesh_;
        bool any=false;
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
            const bool a=f.active_[m.index(i,j)]!=0;
            if(a!=(f.active_[m.index(i&~1,j&~1)]!=0))
                throw std::invalid_argument("An active region must consist of whole parent cells.");
            if(!a) continue;
            any=true;
            for(int ii=std::max(0,i-ghost_);ii<=std::min(m.nz-1,i+ghost_);++ii)
                for(int jj=std::max(0,j-ghost_);jj<=std::min(m.nr-1,j+ghost_);++jj)
                    if(!c.isActive(c.mesh_.index(ii>>1,jj>>1)))
                        throw std::invalid_argument("Level "+std::to_string(l)+" is not properly nested in level "+std::to_string(l-1)+".");
        }
        if(any) top_=l;
    }
    checked_=true;
}
void Amr::averageDown() {
    for(int l=levels()-2;l>=0;--l) averageDown(l);
}
void Amr::averageDown(int l) {
    Flow& c=level(l);const Flow& f=level(l+1);
    const Mesh& cm=c.mesh_;const Mesh& fm=f.mesh_;
    const std::size_t ns=c.ns_,nt=c.nt_;
    c.pool_->blocks(cm.cells.size(),Flow::kCells,[&](std::size_t begin,std::size_t end,int) {
        for(std::size_t p=begin;p<end;++p) {
            const int I=static_cast<int>(p/static_cast<std::size_t>(cm.nr)),J=static_cast<int>(p%static_cast<std::size_t>(cm.nr));
            const std::size_t child[4]={fm.index(2*I,2*J),fm.index(2*I,2*J+1),fm.index(2*I+1,2*J),fm.index(2*I+1,2*J+1)};
            if(!f.isActive(child[0])) continue;
            const double v=cm.cells[p].volume;
            for(int k=0;k<4;++k) {
                double s=0;
                for(auto q:child) s+=fm.cells[q].volume*f.state_[q][k];
                c.state_[p][k]=s/v;
            }
            for(std::size_t k=0;k<ns;++k) {
                double s=0;
                for(auto q:child) s+=fm.cells[q].volume*f.species_[q*ns+k];
                c.species_[p*ns+k]=s/v;
            }
            for(std::size_t k=0;k<nt;++k) {
                double s=0;
                for(auto q:child) s+=fm.cells[q].volume*f.turbulence_[q*nt+k];
                c.turbulence_[p*nt+k]=s/v;
            }
        }
    });
}
void Amr::fill(int l,std::vector<Conserved>& state,std::vector<double>& species,std::vector<double>& turbulence,double theta) {
    Level& cv=levels_[static_cast<std::size_t>(l-1)];
    const Flow& c=*cv.flow;const Flow& f=level(l);
    const Mesh& cm=c.mesh_;const Mesh& fm=f.mesh_;
    const std::size_t ns=c.ns_,nt=c.nt_,stride=4+ns+nt,count=cm.cells.size();
    cv.values.resize(count*stride);
    // Level l - 1 at theta, linear in time (its old state at 0 and its new state at 1 exactly). A
    // convex combination of admissible states is admissible (the internal energy is concave in the
    // conserved variables).
    c.pool_->blocks(count,Flow::kCells,[&](std::size_t begin,std::size_t end,int) {
        for(std::size_t q=begin;q<end;++q) {
            double* v=cv.values.data()+q*stride;
            for(int k=0;k<4;++k) v[k]=(1-theta)*cv.oldState[q][k]+theta*c.state_[q][k];
            for(std::size_t k=0;k<ns;++k) v[4+k]=(1-theta)*cv.oldSpecies[q*ns+k]+theta*c.species_[q*ns+k];
            for(std::size_t k=0;k<nt;++k) v[4+ns+k]=(1-theta)*cv.oldTurbulence[q*nt+k]+theta*c.turbulence_[q*nt+k];
        }
    });
    // Each parent of cells that are not active: rho, rho u_z, rho u_r and E linear with the parent's
    // least-squares gradient over its face neighbours (across the axis its mirror image, rho u_r odd),
    // limited by Barth and Jespersen to the neighbours' range, times the child's offset from the
    // children's volume-weighted centroid, so the children's content is the parent's. The partial
    // densities, rho k and rho omega scale with the child's density. If a child is not admissible, the
    // parent's children all take its values.
    c.pool_->blocks(count,Flow::kCells,[&](std::size_t begin,std::size_t end,int) {
        for(std::size_t p=begin;p<end;++p) {
            const int I=static_cast<int>(p/static_cast<std::size_t>(cm.nr)),J=static_cast<int>(p%static_cast<std::size_t>(cm.nr));
            const std::size_t child[4]={fm.index(2*I,2*J),fm.index(2*I,2*J+1),fm.index(2*I+1,2*J),fm.index(2*I+1,2*J+1)};
            if(f.isActive(child[0])) continue;
            const double* vp=cv.values.data()+p*stride;
            const Cell& cp=cm.cells[p];
            double m00=0,m01=0,m11=0,b[4][2]={},lo[4],hi[4];
            for(int k=0;k<4;++k) lo[k]=hi[k]=vp[k];
            auto add=[&](double dz,double dr,const double* vn) {
                m00+=dz*dz;m01+=dz*dr;m11+=dr*dr;
                for(int k=0;k<4;++k) {
                    const double dv=vn[k]-vp[k];
                    b[k][0]+=dz*dv;b[k][1]+=dr*dv;
                    lo[k]=std::min(lo[k],vn[k]);hi[k]=std::max(hi[k],vn[k]);
                }
            };
            auto neighbour=[&](int i,int j) {
                const std::size_t n=cm.index(i,j);
                add(cm.cells[n].z-cp.z,cm.cells[n].r-cp.r,cv.values.data()+n*stride);
            };
            if(I>0) neighbour(I-1,J);
            if(I<cm.nz-1) neighbour(I+1,J);
            if(J>0) neighbour(I,J-1);
            else { const double mirror[4]={vp[0],vp[1],-vp[2],vp[3]};add(0,-2*cp.r,mirror); }
            if(J<cm.nr-1) neighbour(I,J+1);
            const double det=m00*m11-m01*m01;
            double g[4][2]={};
            if(det>0) for(int k=0;k<4;++k) {
                g[k][0]=(m11*b[k][0]-m01*b[k][1])/det;
                g[k][1]=(m00*b[k][1]-m01*b[k][0])/det;
            }
            double volume=0,zbar=0,rbar=0;
            for(auto q:child) { const Cell& cc=fm.cells[q];volume+=cc.volume;zbar+=cc.volume*cc.z;rbar+=cc.volume*cc.r; }
            zbar/=volume;rbar/=volume;
            double change[4][4],phi[4]={1,1,1,1};
            for(int n=0;n<4;++n) {
                const Cell& cc=fm.cells[child[n]];
                for(int k=0;k<4;++k) {
                    const double d=g[k][0]*(cc.z-zbar)+g[k][1]*(cc.r-rbar);
                    change[n][k]=d;
                    if(d>0) phi[k]=std::min(phi[k],(hi[k]-vp[k])/d);
                    else if(d<0) phi[k]=std::min(phi[k],(lo[k]-vp[k])/d);
                }
            }
            bool ok=true;
            for(int n=0;n<4 && ok;++n) {
                const auto q=child[n];
                for(int k=0;k<4;++k) state[q][k]=vp[k]+phi[k]*change[n][k];
                const double scale=state[q][0]/vp[0];
                for(std::size_t k=0;k<ns;++k) species[q*ns+k]=vp[4+k]*scale;
                for(std::size_t k=0;k<nt;++k) turbulence[q*nt+k]=vp[4+ns+k]*scale;
                ok=f.admissible(state[q],species.data()+q*ns,nt?turbulence.data()+q*nt:nullptr);
            }
            if(!ok) for(auto q:child) {
                for(int k=0;k<4;++k) state[q][k]=vp[k];
                for(std::size_t k=0;k<ns;++k) species[q*ns+k]=vp[4+k];
                for(std::size_t k=0;k<nt;++k) turbulence[q*nt+k]=vp[4+ns+k];
            }
        }
    });
}
// Each uncovered level-l cell beside a covered one takes, at the face between them, the sum of
// level l + 1's face integrals over both substeps and the two fine faces less level l's own, divided
// by its volume: gained at its low face, lost at its high face. Serial, in face order.
bool Amr::reflux(int l) {
    Flow& c=level(l);const Flow& f=level(l+1);
    const Mesh& cm=c.mesh_;const Mesh& fm=f.mesh_;
    const std::size_t ns=c.ns_,nt=c.nt_,stride=4+c.nw_;
    const auto& fine=levels_[static_cast<std::size_t>(l+1)].fineSum;const auto& coarse=c.faceIntegral_;
    const std::size_t naC=static_cast<std::size_t>(cm.nz+1)*cm.nr,naF=static_cast<std::size_t>(fm.nz+1)*fm.nr;
    std::vector<std::size_t> touched;
    auto correct=[&](std::size_t u,double sign,std::size_t fc,std::size_t f0,std::size_t f1) {
        const double v=cm.cells[u].volume;
        for(std::size_t k=0;k<stride;++k) {
            const double change=sign*((fine[f0*stride+k]+fine[f1*stride+k])-coarse[fc*stride+k])/v;
            if(k<4) c.state_[u][k]+=change;
            else if(k<4+ns) c.species_[u*ns+k-4]+=change;
            else c.turbulence_[u*nt+k-4-ns]+=change;
        }
        touched.push_back(u);
    };
    auto isCovered=[&](int i,int j) { return f.isActive(fm.index(2*i,2*j)); };
    for(int I=1;I<cm.nz;++I) for(int J=0;J<cm.nr;++J) {
        const bool low=isCovered(I-1,J),high=isCovered(I,J);
        if(low==high) continue;
        const std::size_t fc=static_cast<std::size_t>(I)*cm.nr+J,f0=static_cast<std::size_t>(2*I)*fm.nr+2*J;
        if(low) correct(cm.index(I,J),1,fc,f0,f0+1); else correct(cm.index(I-1,J),-1,fc,f0,f0+1);
    }
    for(int I=0;I<cm.nz;++I) for(int J=1;J<cm.nr;++J) {
        const bool low=isCovered(I,J-1),high=isCovered(I,J);
        if(low==high) continue;
        const std::size_t fc=naC+static_cast<std::size_t>(I)*(cm.nr+1)+J;
        const std::size_t f0=naF+static_cast<std::size_t>(2*I)*(fm.nr+1)+2*J,f1=f0+(fm.nr+1);
        if(low) correct(cm.index(I,J),1,fc,f0,f1); else correct(cm.index(I,J-1),-1,fc,f0,f1);
    }
    for(auto u:touched)
        if(!c.admissible(c.state_[u],c.species_.data()+u*ns,nt?c.turbulence_.data()+u*nt:nullptr)) return false;
    return true;
}
// One step of level l from t over dt, half the substep (0 or 1) of level l - 1's step it is, then
// level l + 1's two substeps, the average down onto level l and the reflux.
bool Amr::advance(int l,double t,double dt,int half) {
    Level& v=levels_[static_cast<std::size_t>(l)];
    Flow& f=*v.flow;
    if(l>0) fill(l,f.state_,f.species_,f.turbulence_,0.5*half);
    v.oldState=f.state_;v.oldSpecies=f.species_;v.oldTurbulence=f.turbulence_;
    f.time_=t;
    f.beginStep();
    Flow::Attempt attempt;
    if(!f.firstStage(dt,attempt)) return false;
    if(l>0) fill(l,f.stage_,f.speciesStage_,f.turbulenceStage_,0.5*(half+1));
    if(!f.secondStage(dt,attempt)) return false;
    f.accept(dt,t+dt,attempt);
    if(l>0) for(std::size_t k=0;k<v.fineSum.size();++k) v.fineSum[k]+=f.faceIntegral_[k];
    if(l<top_) {
        auto& sum=levels_[static_cast<std::size_t>(l+1)].fineSum;
        std::fill(sum.begin(),sum.end(),0.0);
        if(!advance(l+1,t,0.5*dt,0) || !advance(l+1,t+0.5*dt,0.5*dt,1)) return false;
        averageDown(l);
        if(!reflux(l)) return false;
    }
    return true;
}
void Amr::save() {
    for(auto& v:levels_) {
        const Flow& f=*v.flow;auto& s=v.saved;
        s.state=f.state_;s.species=f.species_;s.turbulence=f.turbulence_;s.temperature=f.temperature_;
        s.time=f.time_;s.dt=f.dt_;s.heat=f.integratedHeat_;s.mass=f.integratedMassFlux_;s.energy=f.integratedEnergyFlux_;
        s.correction=f.integratedCorrection_;s.momentum=f.integratedMomentumSource_;s.gross=f.integratedMomentumGross_;
        s.clipped=f.clippedTurbulentEnergy_;s.prescribed=f.prescribedTurbulentEnergy_;s.steps=f.steps_;
    }
}
void Amr::restore() {
    for(auto& v:levels_) {
        Flow& f=*v.flow;const auto& s=v.saved;
        f.state_=s.state;f.species_=s.species;f.turbulence_=s.turbulence;f.temperature_=s.temperature;
        f.time_=s.time;f.dt_=s.dt;f.integratedHeat_=s.heat;f.integratedMassFlux_=s.mass;f.integratedEnergyFlux_=s.energy;
        f.integratedCorrection_=s.correction;f.integratedMomentumSource_=s.momentum;f.integratedMomentumGross_=s.gross;
        f.clippedTurbulentEnergy_=s.clipped;f.prescribedTurbulentEnergy_=s.prescribed;f.steps_=s.steps;
    }
}
double Amr::step(double maxDt) {
    if(!(maxDt>0)) throw std::invalid_argument("Step interval must be positive.");
    if(!checked_) check();
    for(auto& v:levels_)
        if(!v.flow->bodyForce_.empty() || !v.flow->heating_.empty()) throw std::invalid_argument("AMR step B1 takes no body force or heating.");
    double dt=maxDt;
    for(int l=0;l<=top_;++l) dt=std::min(dt,std::ldexp(level(l).stableDt(),l));
    const double event=level(0).nextEvent(time_);
    bool landing=false;
    if(event-time_<=dt) { dt=event-time_;landing=true; }
    save();
    for(int attempt=0;attempt<14;++attempt) {
        if(advance(0,time_,dt,0)) {
            time_=landing?event:time_+dt;++steps_;
            for(auto& v:levels_) v.flow->time_=time_;
            return dt;
        }
        restore();++rejected_;dt*=0.5;landing=false;
    }
    throw std::runtime_error("AMR: could not advance an admissible state after 14 base-step reductions.");
}
void Amr::advanceTo(double target) {
    if(!std::isfinite(target) || target<time_) throw std::invalid_argument("Invalid target time.");
    while(time_<target) step(target-time_);
}

}  // namespace crucible
