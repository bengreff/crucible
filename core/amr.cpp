#include "core/amr.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
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
    ghost_=b.nt_?3:2;
    if(b.nt_ && (base.turbulence.wallFunctions || base.turbulence.propertyCorrections))
        throw std::invalid_argument("AMR does not yet carry the wall functions' prescription or the SST property corrections across levels.");
    for(int l=0;l<=finer;++l) {
        Flow& f=level(l);
        f.recordFaces_=true;
        f.faceIntegral_.assign((f.axialFlux_.size()+f.radialFlux_.size())*(4+f.nw_),0.0);
        if(l==0) continue;
        levels_[static_cast<std::size_t>(l)].fineSum.assign(f.faceIntegral_.size(),0.0);
        f.active_.assign(f.mesh_.cells.size(),0);
        // A fine injector ring takes its base parent ring's supply, so the supply areas are the base's.
        if(base.experiment==Case::Chamber) {
            std::fill(f.supplyArea_.begin(),f.supplyArea_.end(),0.0);
            for(int j=0;j<f.mesh_.nr;++j) {
                const int s=f.faceSupply_[static_cast<std::size_t>(j)]=b.faceSupply_[static_cast<std::size_t>(j>>l)];
                if(s>=0) f.supplyArea_[static_cast<std::size_t>(s)]+=f.mesh_.axialArea(0,j);
            }
        }
    }
    // The igniter's region is the finest level's cells by centroid, as a single level on that mesh
    // chooses them; a coarser cell's igniter volume is the sum of its children's, so the owned cells'
    // shares of the power (power times the cell's igniter volume over the total) sum to the whole.
    if(base.experiment==Case::Chamber && base.igniter.energy>0) {
        const auto& g=base.igniter;
        Flow& top=level(finer);
        const Mesh& tm=top.mesh_;
        top.igniterCells_.clear();top.igniterVolume_=0;
        for(std::size_t q=0;q<tm.cells.size();++q) {
            const Cell& c=tm.cells[q];
            const bool in=c.z>=g.zMin && c.z<=g.zMax && c.r<=g.rMax;
            top.igniterCell_[q]=in?1:0;top.igniterCellVolume_[q]=in?c.volume:0.0;
            if(in) { top.igniterCells_.push_back(q);top.igniterVolume_+=c.volume; }
        }
        if(top.igniterCells_.empty()) throw std::invalid_argument("The igniter region contains no centroid of the finest level.");
        for(int l=finer-1;l>=0;--l) {
            Flow& c=level(l);const Flow& f=level(l+1);
            const Mesh& cm=c.mesh_;const Mesh& fm=f.mesh_;
            c.igniterCells_.clear();c.igniterVolume_=top.igniterVolume_;
            for(int I=0;I<cm.nz;++I) for(int J=0;J<cm.nr;++J) {
                const auto p=cm.index(I,J);
                double v=0;
                for(auto q:{fm.index(2*I,2*J),fm.index(2*I,2*J+1),fm.index(2*I+1,2*J),fm.index(2*I+1,2*J+1)}) v+=f.igniterCellVolume_[q];
                c.igniterCell_[p]=v>0?1:0;c.igniterCellVolume_[p]=v;
                if(v>0) c.igniterCells_.push_back(p);
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
// level l (clipped to the domain) lies over an active cell of level l - 1. Then each level's owned
// cells (active and not covered), where it books its boundary exchange and sources, and in B2 the
// regions of the levels above the base.
void Amr::check() {
    if(ghost_<2 || (restricted_ && level(0).nt_ && ghost_<3))
        throw std::invalid_argument("The ghost width must be at least 2, and 3 with turbulence in B2.");
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
    for(int l=0;l<levels();++l) {
        Flow& f=level(l);
        f.owned_.assign(f.mesh_.cells.size(),0);
        for(std::size_t q=0;q<f.owned_.size();++q) f.owned_[q]=f.isActive(q) && !covered(l,q)?1:0;
    }
    // B2: the parents of the ghosts (N_g less the active cells), all active by the nesting; their
    // children are filled, and the other non-active cells (off cells) keep their state through the
    // stages, so their stage and new values are set to it here once.
    for(int l=1;l<levels();++l) {
        Level& v=levels_[static_cast<std::size_t>(l)];
        Flow& f=*v.flow;const Flow& c=level(l-1);
        const Mesh& m=f.mesh_;const Mesh& cm=c.mesh_;
        v.parents.clear();v.around.clear();
        if(!restricted_) { f.clearRegions();continue; }
        f.setRegions(ghost_);
        const auto nr=static_cast<std::size_t>(m.nr);
        auto parent=[&](std::size_t q) { return cm.index(static_cast<int>(q/nr)>>1,static_cast<int>(q%nr)>>1); };
        std::vector<unsigned char> mark(cm.cells.size(),0),near(cm.cells.size(),0);
        for(std::size_t q:f.within_[static_cast<std::size_t>(ghost_)]) if(!f.isActive(q)) mark[parent(q)]=1;
        for(std::size_t p=0;p<mark.size();++p) {
            if(!mark[p]) continue;
            v.parents.push_back(p);
            const int I=static_cast<int>(p/static_cast<std::size_t>(cm.nr)),J=static_cast<int>(p%static_cast<std::size_t>(cm.nr));
            near[p]=1;
            if(I>0) near[cm.index(I-1,J)]=1;
            if(I<cm.nz-1) near[cm.index(I+1,J)]=1;
            if(J>0) near[cm.index(I,J-1)]=1;
            if(J<cm.nr-1) near[cm.index(I,J+1)]=1;
        }
        for(std::size_t p=0;p<near.size();++p) if(near[p]) v.around.push_back(p);
        f.filled_.clear();
        const std::size_t ns=f.ns_,nt=f.nt_;
        for(std::size_t q=0;q<m.cells.size();++q) {
            if(f.isActive(q)) continue;
            if(mark[parent(q)]) f.filled_.push_back(q);
            f.stage_[q]=f.next_[q]=f.state_[q];
            for(std::size_t k=q*ns;k<(q+1)*ns;++k) f.speciesStage_[k]=f.speciesNext_[k]=f.species_[k];
            for(std::size_t k=q*nt;k<(q+1)*nt;++k) f.turbulenceStage_[k]=f.turbulenceNext_[k]=f.turbulence_[k];
        }
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
    Level& cv=levels_[static_cast<std::size_t>(l-1)];const Level& fv=levels_[static_cast<std::size_t>(l)];
    const Flow& c=*cv.flow;Flow& f=*fv.flow;
    const Mesh& cm=c.mesh_;const Mesh& fm=f.mesh_;
    const std::size_t ns=c.ns_,nt=c.nt_,stride=4+ns+nt,count=cm.cells.size();
    cv.values.resize(count*stride);
    // Level l - 1 at theta, linear in time (its old state at 0 and its new state at 1 exactly). A
    // convex combination of admissible states is admissible (the internal energy is concave in the
    // conserved variables).
    const auto* around=restricted_?&fv.around:nullptr;
    c.pool_->blocks(around?around->size():count,Flow::kCells,[&](std::size_t begin,std::size_t end,int) {
        for(std::size_t n=begin;n<end;++n) {
            const std::size_t q=around?(*around)[n]:n;
            double* v=cv.values.data()+q*stride;
            for(int k=0;k<4;++k) v[k]=(1-theta)*cv.oldState[q][k]+theta*c.state_[q][k];
            for(std::size_t k=0;k<ns;++k) v[4+k]=(1-theta)*cv.oldSpecies[q*ns+k]+theta*c.species_[q*ns+k];
            for(std::size_t k=0;k<nt;++k) v[4+ns+k]=(1-theta)*cv.oldTurbulence[q*nt+k]+theta*c.turbulence_[q*nt+k];
        }
    });
    // Each parent of cells that are not active, by spread.
    const auto* parents=restricted_?&fv.parents:nullptr;
    c.pool_->blocks(parents?parents->size():count,Flow::kCells,[&](std::size_t begin,std::size_t end,int) {
        for(std::size_t n=begin;n<end;++n) {
            const std::size_t p=parents?(*parents)[n]:n;
            const int I=static_cast<int>(p/static_cast<std::size_t>(cm.nr)),J=static_cast<int>(p%static_cast<std::size_t>(cm.nr));
            if(f.isActive(fm.index(2*I,2*J))) continue;
            spread(l,p,state,species,turbulence);
        }
    });
}
// rho, rho u_z, rho u_r and E linear with the parent's least-squares gradient over its face neighbours
// (across the axis its mirror image, rho u_r odd), limited by Barth and Jespersen to the neighbours'
// range, times the child's offset from the children's volume-weighted centroid, so the children's
// content is the parent's. The partial densities, rho k and rho omega scale with the child's density.
// If a child is not admissible, the parent's children all take its values.
void Amr::spread(int l,std::size_t p,std::vector<Conserved>& state,std::vector<double>& species,std::vector<double>& turbulence) {
    const Level& cv=levels_[static_cast<std::size_t>(l-1)];
    const Flow& c=*cv.flow;Flow& f=level(l);
    const Mesh& cm=c.mesh_;const Mesh& fm=f.mesh_;
    const std::size_t ns=c.ns_,nt=c.nt_,stride=4+ns+nt;
    const int I=static_cast<int>(p/static_cast<std::size_t>(cm.nr)),J=static_cast<int>(p%static_cast<std::size_t>(cm.nr));
    const std::size_t child[4]={fm.index(2*I,2*J),fm.index(2*I,2*J+1),fm.index(2*I+1,2*J),fm.index(2*I+1,2*J+1)};
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
    for(auto q:child) f.temperature_[q]=c.temperature_[p];
}
void Amr::setTags(Tags tags) {
    if(tags.interval<1 || tags.buffer<0) throw std::invalid_argument("The regrid interval must be at least 1 and the buffer at least 0.");
    if(!tags.hydrogen.empty() && tags.hydrogen.size()!=level(0).ns_) throw std::invalid_argument("The f_H weights need one entry per species.");
    tags_=std::move(tags);tagged_=true;
}
// Tags on level l over its active cells (AMR.md, Step C: setup), from its refreshed primitives:
// across each face to an active face neighbour (f_H, T and p, both cells), vorticity times the side of
// the cell's meridional square against 0.05 of its sound speed, the throat and the wall strip; then
// buffered. The level-l set is those together with the parents of the level-(l + 1) cells within the
// ghost width of level l + 2's new active cells, so each new level nests in the one below by
// construction, and level l + 1's new active cells are its children. New cells are filled from the
// bottom up by spread from the level below at the regrid time; cells leaving a level keep nothing.
void Amr::regrid() {
    if(!checked_) check();
    const int top=levels()-1;
    const Tags& tg=tags_;
    const bool first=!regridded_;
    std::vector<std::vector<unsigned char>> set(static_cast<std::size_t>(top));
    for(int l=top-1;l>=0;--l) {
        Flow& f=level(l);const Mesh& m=f.mesh_;
        const std::size_t count=m.cells.size(),ns=f.ns_,nw=f.nw_;
        std::vector<unsigned char> tag(count,0);
        f.refresh(f.state_,f.species_,f.turbulence_,0);
        auto on=[&](int i,int j) { return i>=0 && i<m.nz && j>=0 && j<m.nr && f.isActive(m.index(i,j)); };
        auto mixture=[&](std::size_t q) {
            const double* y=f.fractions_.data()+q*nw;
            double s=0;
            for(std::size_t k=0;k<ns;++k) s+=tg.hydrogen[k]*y[k];
            return s;
        };
        auto across=[&](std::size_t a,std::size_t b) {
            const double ta=f.temperature_[a],tb=f.temperature_[b],pa=f.primitives_[a].p,pb=f.primitives_[b].p;
            if((!tg.hydrogen.empty() && tg.mixture>0 && std::abs(mixture(a)-mixture(b))>tg.mixture) ||
               (tg.temperature>0 && std::abs(ta-tb)>tg.temperature*std::min(ta,tb)) ||
               (tg.pressure>0 && std::abs(pa-pb)>tg.pressure*std::min(pa,pb))) tag[a]=tag[b]=1;
        };
        // The throat: the wall's narrowest station strictly inside the domain (none if the narrowest is an end).
        int throat=-1;double rt=0;
        if(tg.throatHalfWidth>0 && l+1<=tg.throatLevel) {
            int it=0;
            for(int i=1;i<=m.nz;++i) if(m.radius[static_cast<std::size_t>(i)]<m.radius[static_cast<std::size_t>(it)]) it=i;
            rt=m.radius[static_cast<std::size_t>(it)];
            if(it>0 && it<m.nz && rt<m.radius.front() && rt<m.radius.back()) throat=it;
        }
        const bool wall=tg.wallYPlus>0 && l+1<=tg.throatLevel;
        if(wall && !first && f.viscosity_.size()<count) throw std::invalid_argument("The wall strip needs molecular transport.");
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
            if(!on(i,j)) continue;
            const std::size_t q=m.index(i,j);
            const Cell& c=m.cells[q];
            if(on(i+1,j)) across(q,m.index(i+1,j));
            if(on(i,j+1)) across(q,m.index(i,j+1));
            if(tg.vorticity>0) {
                // The fill's unlimited least-squares fit of u_z and u_r over the active face neighbours.
                const Primitive& w=f.primitives_[q];
                double m00=0,m01=0,m11=0,bz[2]={},br[2]={};
                auto add=[&](double dz,double dr,double uz,double ur) {
                    m00+=dz*dz;m01+=dz*dr;m11+=dr*dr;
                    bz[0]+=dz*(uz-w.uz);bz[1]+=dr*(uz-w.uz);br[0]+=dz*(ur-w.ur);br[1]+=dr*(ur-w.ur);
                };
                auto neighbour=[&](int ii,int jj) {
                    if(!on(ii,jj)) return;
                    const std::size_t n=m.index(ii,jj);
                    add(m.cells[n].z-c.z,m.cells[n].r-c.r,f.primitives_[n].uz,f.primitives_[n].ur);
                };
                neighbour(i-1,j);neighbour(i+1,j);neighbour(i,j+1);
                if(j>0) neighbour(i,j-1); else add(0,-2*c.r,w.uz,-w.ur);
                const double det=m00*m11-m01*m01;
                if(det>0) {
                    const double durdz=(m11*br[0]-m01*br[1])/det,duzdr=(m00*bz[1]-m01*bz[0])/det;
                    const double side=std::sqrt(c.volume/(2*std::numbers::pi*c.r));
                    if(std::abs(durdz-duzdr)*side>tg.vorticity*f.sound_[q]) tag[q]=1;
                }
            }
            if(throat>=0 && std::abs(c.z-throat*m.dz)<=tg.throatHalfWidth*rt) tag[q]=1;
            if(wall && j==m.nr-1) {
                if(first) tag[q]=1;
                else {
                    // Laminar y+: u_tau from mu |u_t| / y, y the centroid's distance to the wall segment.
                    const double r0=m.radius[static_cast<std::size_t>(i)],r1=m.radius[static_cast<std::size_t>(i+1)];
                    const double len=std::hypot(m.dz,r1-r0);
                    const Primitive& w=f.primitives_[q];
                    const double ut=std::abs(w.uz*m.dz+w.ur*(r1-r0))/len;
                    const double y=(r0+(r1-r0)*(c.z-i*m.dz)/m.dz-c.r)*m.dz/len;
                    if(std::sqrt(w.rho*ut*y/f.viscosity_[q])>tg.wallYPlus) tag[q]=1;
                }
            }
        }
        auto& s=set[static_cast<std::size_t>(l)];
        s.assign(count,0);
        const int b=tg.buffer;
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) {
            if(!tag[m.index(i,j)]) continue;
            for(int ii=std::max(0,i-b);ii<=std::min(m.nz-1,i+b);++ii)
                for(int jj=std::max(0,j-b);jj<=std::min(m.nr-1,j+b);++jj) s[m.index(ii,jj)]=1;
        }
        // Nesting: level l + 2's new active cells are the children of set[l + 1].
        if(l+1<top) {
            const auto& above=set[static_cast<std::size_t>(l+1)];
            const Mesh& am=level(l+1).mesh_;const Mesh& tm=level(l+2).mesh_;
            for(int I=0;I<am.nz;++I) for(int J=0;J<am.nr;++J) {
                if(!above[am.index(I,J)]) continue;
                for(int ii=std::max(0,2*I-ghost_);ii<=std::min(tm.nz-1,2*I+1+ghost_);++ii)
                    for(int jj=std::max(0,2*J-ghost_);jj<=std::min(tm.nr-1,2*J+1+ghost_);++jj) s[m.index(ii>>2,jj>>2)]=1;
            }
        }
    }
    // The new active sets; nothing is touched if none changes.
    std::vector<std::vector<unsigned char>> next(static_cast<std::size_t>(levels()));
    bool changed=false;
    for(int l=1;l<=top;++l) {
        const Mesh& m=level(l).mesh_;const Mesh& cm=level(l-1).mesh_;
        auto& a=next[static_cast<std::size_t>(l)];
        a.assign(m.cells.size(),0);
        for(int i=0;i<m.nz;++i) for(int j=0;j<m.nr;++j) a[m.index(i,j)]=set[static_cast<std::size_t>(l-1)][cm.index(i>>1,j>>1)];
        changed=changed || a!=level(l).active_;
    }
    regridded_=true;
    if(!changed) return;
    fillError_=fillGross_=0;
    for(int l=1;l<=top;++l) {
        Level& cv=levels_[static_cast<std::size_t>(l-1)];
        Flow& f=level(l);const Flow& c=*cv.flow;
        const Mesh& m=f.mesh_;const Mesh& cm=c.mesh_;
        const std::size_t ns=c.ns_,nt=c.nt_,stride=4+ns+nt,count=cm.cells.size();
        cv.values.resize(count*stride);
        for(std::size_t q=0;q<count;++q) {
            double* v=cv.values.data()+q*stride;
            for(int k=0;k<4;++k) v[k]=c.state_[q][k];
            for(std::size_t k=0;k<ns;++k) v[4+k]=c.species_[q*ns+k];
            for(std::size_t k=0;k<nt;++k) v[4+ns+k]=c.turbulence_[q*nt+k];
        }
        const auto& a=next[static_cast<std::size_t>(l)];
        for(int I=0;I<cm.nz;++I) for(int J=0;J<cm.nr;++J) {
            const std::size_t p=cm.index(I,J),q0=m.index(2*I,2*J);
            if(!a[q0] || f.isActive(q0)) continue;
            spread(l,p,f.state_,f.species_,f.turbulence_);
            // Check 2: the children's content against the parent's.
            const double* vp=cv.values.data()+p*stride;
            const double vol=cm.cells[p].volume;
            const std::size_t child[4]={q0,m.index(2*I,2*J+1),m.index(2*I+1,2*J),m.index(2*I+1,2*J+1)};
            // As stated, over the parent's own content (momentum's magnitude for both components); and over the
            // children's gross content, the sum of their |content|, which stays well conditioned where a content
            // nears zero (the total energy, which carries the formation enthalpies, changes sign in this flow).
            auto error=[&](double sum,double parent,double own,double gross) {
                if(own>0) fillError_=std::max(fillError_,std::abs(sum-parent)/own);
                if(gross>0) fillGross_=std::max(fillGross_,std::abs(sum-parent)/gross);
            };
            const double momentum=std::hypot(vp[1],vp[2])*vol;
            for(int k=0;k<4;++k) {
                double sum=0,gross=0;
                for(auto q:child) { const double v=m.cells[q].volume*f.state_[q][k];sum+=v;gross+=std::abs(v); }
                error(sum,vol*vp[k],k==1 || k==2?momentum:std::abs(vol*vp[k]),gross);
            }
            for(std::size_t k=0;k<ns;++k) {
                double sum=0,gross=0;
                for(auto q:child) { const double v=m.cells[q].volume*f.species_[q*ns+k];sum+=v;gross+=std::abs(v); }
                error(sum,vol*vp[4+k],std::abs(vol*vp[4+k]),gross);
            }
        }
        f.active_=a;
    }
    ++regrids_;checked_=false;
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
// level l + 1's two substeps, the average down onto level l and the reflux. With turbulence the SST
// source coefficients are taken at the step's start, filled cells included (and again after a redone
// base step), by stableDt as Flow::step takes them.
bool Amr::advance(int l,double t,double dt,int half) {
    Level& v=levels_[static_cast<std::size_t>(l)];
    Flow& f=*v.flow;
    if(l>0) fill(l,f.state_,f.species_,f.turbulence_,0.5*half);
    v.oldState=f.state_;v.oldSpecies=f.species_;v.oldTurbulence=f.turbulence_;
    f.time_=t;
    if(f.nt_) f.stableDt();
    f.beginStep();
    Flow::Attempt attempt;
    if(!f.firstStage(dt,attempt)) return false;
    if(l>0) fill(l,f.stage_,f.speciesStage_,f.turbulenceStage_,0.5*(half+1));
    if(!f.secondStage(dt,attempt)) return false;
    f.accept(dt,t+dt,attempt);
    if(afterStep_) afterStep_(l,dt);
    // In B2 the registers hold only the faces a level computes (Flow::faces_), which contain every
    // coarse-fine face.
    const std::size_t stride=4+f.nw_;
    if(l>0) {
        if(!restricted_) for(std::size_t k=0;k<v.fineSum.size();++k) v.fineSum[k]+=f.faceIntegral_[k];
        else for(std::size_t e:f.faces_) for(std::size_t k=e*stride;k<(e+1)*stride;++k) v.fineSum[k]+=f.faceIntegral_[k];
    }
    if(l<top_) {
        Level& u=levels_[static_cast<std::size_t>(l+1)];
        auto& sum=u.fineSum;
        if(!restricted_) std::fill(sum.begin(),sum.end(),0.0);
        else for(std::size_t e:u.flow->faces_) std::fill_n(sum.begin()+static_cast<std::ptrdiff_t>(e*stride),stride,0.0);
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
        s.time=f.time_;s.dt=f.dt_;s.heat=f.integratedHeat_;s.heatCompensation=f.heatCompensation_;s.mass=f.integratedMassFlux_;s.energy=f.integratedEnergyFlux_;
        s.correction=f.integratedCorrection_;s.momentum=f.integratedMomentumSource_;s.gross=f.integratedMomentumGross_;
        s.clipped=f.clippedTurbulentEnergy_;s.prescribed=f.prescribedTurbulentEnergy_;s.steps=f.steps_;
    }
}
void Amr::restore() {
    for(auto& v:levels_) {
        Flow& f=*v.flow;const auto& s=v.saved;
        f.state_=s.state;f.species_=s.species;f.turbulence_=s.turbulence;f.temperature_=s.temperature;
        f.time_=s.time;f.dt_=s.dt;f.integratedHeat_=s.heat;f.heatCompensation_=s.heatCompensation;f.integratedMassFlux_=s.mass;f.integratedEnergyFlux_=s.energy;
        f.integratedCorrection_=s.correction;f.integratedMomentumSource_=s.momentum;f.integratedMomentumGross_=s.gross;
        f.clippedTurbulentEnergy_=s.clipped;f.prescribedTurbulentEnergy_=s.prescribed;f.steps_=s.steps;
    }
}
Amr::Budget Amr::content() const {
    const Flow& b=level(0);const Mesh& m=b.mesh_;
    Budget c;
    for(std::size_t q=0;q<m.cells.size();++q) {
        const double v=m.cells[q].volume;
        c.mass+=b.state_[q][0]*v;c.energy+=b.state_[q][3]*v;c.momentum+=b.state_[q][1]*v;
    }
    return c;
}
Amr::Budget Amr::budget() const {
    Budget c=content();
    c.massResidual=c.mass-initial_.mass;c.energyResidual=c.energy-initial_.energy;c.momentumResidual=c.momentum-initial_.momentum;
    for(const auto& v:levels_) {
        const Flow& f=*v.flow;
        c.massResidual-=f.integratedMassFlux_+f.exchangedMass_;
        c.energyResidual-=f.integratedEnergyFlux_+f.integratedCorrection_+f.exchangedEnergy_;
        c.momentumResidual-=f.integratedMomentumSource_+f.exchangedMomentum_;
        c.momentumGross+=f.integratedMomentumGross_;c.igniterEnergy+=f.deliveredHeat();
    }
    return c;
}
double Amr::step(double maxDt) {
    if(!(maxDt>0)) throw std::invalid_argument("Step interval must be positive.");
    if(tagged_ && !regridded_) regrid();
    if(!checked_) check();
    if(!started_) { initial_=content();started_=true; }
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
            if(tagged_ && steps_%static_cast<std::uint64_t>(tags_.interval)==0) regrid();
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
