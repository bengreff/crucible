#include "core/equilibrium_table.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>

namespace crucible {
namespace {
constexpr char kMagic[]="CRUCIBLE equilibrium table 3\n";
// Index i of the interval [x_i, x_i+1] holding v (clamped to the axis) and the weight of x_i+1.
std::size_t interval(const std::vector<double>& x, double v, double& w) {
    v=std::clamp(v,x.front(),x.back());
    std::size_t i=static_cast<std::size_t>(std::upper_bound(x.begin(),x.end(),v)-x.begin());
    i=std::clamp<std::size_t>(i,1,x.size()-1)-1;
    w=(v-x[i])/(x[i+1]-x[i]);
    return i;
}
template<class T> void put(std::ofstream& o,const T& v) { o.write(reinterpret_cast<const char*>(&v),sizeof v); }
template<class T> void putVector(std::ofstream& o,const std::vector<T>& v) {
    put(o,static_cast<std::uint64_t>(v.size()));
    o.write(reinterpret_cast<const char*>(v.data()),static_cast<std::streamsize>(v.size()*sizeof(T)));
}
void putString(std::ofstream& o,const std::string& s) { put(o,static_cast<std::uint64_t>(s.size()));o.write(s.data(),static_cast<std::streamsize>(s.size())); }
template<class T> T get(std::ifstream& in) { T v{};in.read(reinterpret_cast<char*>(&v),sizeof v);return v; }
template<class T> std::vector<T> getVector(std::ifstream& in) {
    std::vector<T> v(get<std::uint64_t>(in));
    in.read(reinterpret_cast<char*>(v.data()),static_cast<std::streamsize>(v.size()*sizeof(T)));
    return v;
}
std::string getString(std::ifstream& in) { std::string s(get<std::uint64_t>(in),'\0');in.read(s.data(),static_cast<std::streamsize>(s.size()));return s; }
}

EquilibriumTable::Result EquilibriumTable::lookup(const double* y, double e, double rho, double* out, double& t) const {
    const std::size_t ns=species.size(),nv=width(),stride=static_cast<std::size_t>(nRho);
    double zh=0,zo=0;
    for(std::size_t k=0;k<ns;++k) { zh+=y[k]*elements[k][0];zo+=y[k]*elements[k][1]; }
    const double reacting=zh+zo;
    double wf;
    const std::size_t i=interval(f,reacting>0 ? zh/reacting : 0.5,wf);
    double lr=std::log(rho*reacting);
    if(lr>lnRhoMax || !(rho>0)) return Result::Outside;
    bool clamped=false;
    if(!(lr>=lnRhoMin)) { lr=lnRhoMin;clamped=true; }
    const double sr=(lr-lnRhoMin)/(lnRhoMax-lnRhoMin)*(nRho-1);
    const std::size_t l=std::min(static_cast<std::size_t>(sr),static_cast<std::size_t>(nRho-2));
    const double wr=sr-static_cast<double>(l);
    // The four (f, rho_r) corners: offsets of their T = tMin nodes, and weights per unit reacting mass.
    std::size_t base[4];
    double weight[4];
    for(int c=0;c<4;++c) {
        base[c]=node(i+(c&1),0,l+(c>>1));
        weight[c]=reacting*((c&1) ? wf : 1-wf)*((c&2) ? wr : 1-wr);
    }
    auto energy=[&](std::size_t k) {
        double v=0;
        for(int c=0;c<4;++c) v+=weight[c]*values[(base[c]+k*stride)*nv];
        for(std::size_t q=0;q<inert.size();++q) v+=y[inert[q]]*inertEnergy[q*static_cast<std::size_t>(nT)+k];
        return v;
    };
    std::size_t lo=0,hi=static_cast<std::size_t>(nT-1);
    const double eLo=energy(lo),eHi=energy(hi);
    if(!(e<=eHi)) return Result::Outside;
    double wt=0;
    if(e<eLo) clamped=true;
    else {
        double a=eLo,b=eHi;
        while(hi-lo>1) {
            const std::size_t mid=(lo+hi)/2;
            const double m=energy(mid);
            if(m<=e) { lo=mid;a=m; } else { hi=mid;b=m; }
        }
        wt=(e-a)/(b-a);
    }
    t=temperature(lo)+wt*(temperature(lo+1)-temperature(lo));
    std::fill(out,out+ns,0.0);
    for(int c=0;c<4;++c)
        for(int d=0;d<2;++d) {
            const double w=weight[c]*(d ? wt : 1-wt);
            if(w==0) continue;
            const double* p=values.data()+(base[c]+(lo+static_cast<std::size_t>(d))*stride)*nv+1;
            for(std::size_t s=0;s<ns;++s) out[s]+=w*p[s];
        }
    for(std::size_t q:inert) out[q]=y[q];
    return clamped ? Result::Clamped : Result::Table;
}

void EquilibriumTable::save(const std::string& path) const {
    std::ofstream o(path,std::ios::binary);
    if(!o) throw std::runtime_error("Cannot write "+path);
    o.write(kMagic,sizeof kMagic-1);
    putString(o,provenance);
    put(o,static_cast<std::uint64_t>(species.size()));
    for(const auto& s:species) putString(o,s);
    putVector(o,elements);putVector(o,inert);putVector(o,f);
    put(o,nT);put(o,nRho);put(o,tMin);put(o,tMax);put(o,lnRhoMin);put(o,lnRhoMax);
    putVector(o,inertEnergy);putVector(o,values);
    if(!o) throw std::runtime_error("Write failed: "+path);
}

EquilibriumTable EquilibriumTable::load(const std::string& path) {
    std::ifstream in(path,std::ios::binary);
    if(!in) throw std::runtime_error("Cannot read equilibrium table "+path);
    std::string magic(sizeof kMagic-1,'\0');
    in.read(magic.data(),static_cast<std::streamsize>(magic.size()));
    if(magic!=kMagic) throw std::runtime_error(path+" is not an equilibrium table.");
    EquilibriumTable t;
    t.provenance=getString(in);
    t.species.resize(get<std::uint64_t>(in));
    for(auto& s:t.species) s=getString(in);
    t.elements=getVector<std::array<double,2>>(in);t.inert=getVector<std::size_t>(in);t.f=getVector<double>(in);
    t.nT=get<int>(in);t.nRho=get<int>(in);t.tMin=get<double>(in);t.tMax=get<double>(in);
    t.lnRhoMin=get<double>(in);t.lnRhoMax=get<double>(in);
    t.inertEnergy=getVector<double>(in);t.values=getVector<double>(in);
    bool inertOk=true;
    for(std::size_t q:t.inert) inertOk=inertOk && q<t.species.size();
    if(!in || !inertOk || t.elements.size()!=t.species.size() || t.f.size()<2 || t.nT<2 || t.nRho<2 ||
       t.inertEnergy.size()!=t.inert.size()*static_cast<std::size_t>(t.nT) ||
       t.values.size()!=t.f.size()*static_cast<std::size_t>(t.nT)*static_cast<std::size_t>(t.nRho)*t.width())
        throw std::runtime_error(path+" is truncated or inconsistent.");
    return t;
}
} // namespace crucible
