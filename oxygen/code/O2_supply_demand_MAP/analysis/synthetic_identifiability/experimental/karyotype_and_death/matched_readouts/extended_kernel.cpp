#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

namespace {
constexpr int N = 133, M = 14, NC = 46;
constexpr double dt = 0.05, scale = 1e12, decay = 1.-dt*1e-4;
using V = std::array<double,N>;
using D = std::array<double,M>;
struct Dual {
  double v;
  D d{};
  Dual(double v_ = 0):v(v_) {}
  static Dual par(double x,int j){Dual y(x); y.d[j]=1; return y;}
};
Dual operator+(Dual x, Dual y){Dual r(x.v+y.v); for(int k=0;k<M;++k)r.d[k]=x.d[k]+y.d[k];return r;}
Dual operator-(Dual x, Dual y){Dual r(x.v-y.v); for(int k=0;k<M;++k)r.d[k]=x.d[k]-y.d[k];return r;}
Dual operator*(Dual x, Dual y){Dual r(x.v*y.v); for(int k=0;k<M;++k)r.d[k]=x.d[k]*y.v+y.d[k]*x.v;return r;}
Dual operator/(Dual x, Dual y){Dual r(x.v/y.v); for(int k=0;k<M;++k)r.d[k]=(x.d[k]-r.v*y.d[k])/y.v;return r;}
Dual expd(Dual x){Dual r(std::exp(x.v));for(int k=0;k<M;++k)r.d[k]=r.v*x.d[k];return r;}
Dual powd(Dual x,int m){if(m==0)return 1;double a=std::pow(x.v,m-1);Dual r(a*x.v);for(int k=0;k<M;++k)r.d[k]=m*a*x.d[k];return r;}

struct Entry {int to; double v;};
struct Rates {
  std::array<std::vector<Entry>,N> g;
  std::array<std::vector<D>,N> dg;
  V mu,b;
  std::array<D,N> dm,db;
  unsigned long long mask=1469598103934665603ULL;
};
struct Slot {int tick,row;double w;};
struct Segment {
  int oxygen,cohort,parent,donor,last,tilt;
  double seed;
  std::vector<Slot> slots;
  std::vector<double> demand;
  std::vector<int> demand_index;
  std::vector<V> x, external;
  int adjoint_last=-1;
  V before{}, base_pi{}, init_pi{};
};
struct Design {
  std::vector<Segment> passages;
  std::vector<double> oxygens;
  std::array<V,2> histogram{};
  int nslope=0;
  int nobs=0;
  std::vector<int> group;
  std::vector<double> weights;
  std::vector<double> logc, cells, margins, live, visible_h, visible_m;
  std::vector<Rates> cached_rates;
  std::array<double,NC> cached_p{};
  bool cache_valid=false;
  bool cached_gaussian=false;
};

double sum(const V& v){ double r=0;for(double x:v) r+=x;return r;}
double dot(const V& a,const V& b){double r=0;for(int i=0;i<N;++i)r+=a[i]*b[i];return r;}

Rates build(double oxygen,const double* p) {
  auto P=[&](int k){return Dual::par(p[k],k);};
  Dual w=P(4);
  Dual h=1.;
  if(oxygen>0.) {
    // Differentiate h explicitly to keep the positive Hill branch stable.
    h.v=1./(1.+std::exp(p[13]*std::log(oxygen/p[12])));
    h.d.fill(0.);
    h.d[12]=p[13]*h.v*(1-h.v)/p[12];
    h.d[13]=h.v*(1-h.v)*std::log(p[12]/oxygen);
  }
  Rates r;
  for(int j=0;j<N;++j){
    const int n=j+22;
    const double lu=std::log(n/44.);
    Dual lam=P(0)/(1.+P(8)*h*expd(P(9)*lu));
    Dual mu=P(10)*h*expd(P(11)*lu);
    Dual pm=P(1)+P(2)*mu/(P(3)+mu);
    Dual sv=P(5)*expd((0.-P(6))*expd((0.-P(7))*lu));
    Dual ord=lam*(1.-w);
    Dual bio=2.*ord*(1.-powd(1.-pm*(1.-sv),n));
    r.mu[j]=mu.v; r.dm[j]=mu.d; r.b[j]=bio.v; r.db[j]=bio.d;
    r.g[j].reserve(30);
    r.dg[j].reserve(30);
    auto add=[&](int i,Dual x){r.g[j].push_back({i,x.v});r.dg[j].push_back(x.d);};
    add(j,0.-lam);
    if(2*n<=154)add(2*n-22,lam*w);
    Dual q=powd(1.-pm,n), survival=1;
    Dual odds=pm/(1.-pm);
    for(int m=0;m<=n;++m){
      if(q.v>=1e-8) {
        Dual flow=ord*q*survival;
        if(j+m<N)add(j+m,flow);
        if(j-m>=0)add(j-m,flow);
        r.mask=(r.mask^static_cast<unsigned long long>(1+m+j*154))*1099511628211ULL;
      }
      // Once below cutoff after the mode, the remaining tail is also below.
      if(q.v<1e-8 && m+1>(n+1)*pm.v)break;
      q=q*(double(n-m)/(m+1))*odds;
      survival=survival*sv;
    }
  }
  return r;
}

void mult(const Rates &r, const V& x, V& gx) {
  gx.fill(0.);
  for(int i=0;i<N;++i) {
    double z=x[i];
    for(auto& a:r.g[i])gx[a.to]+=a.v*z;
  }
}

void init(Design& d,int index,const double* p,bool gaussian) {
  auto& s=d.passages[index];
  auto& v=s.x[0];
  if(s.parent<0){
    if(gaussian){
      double m=p[18+2*s.cohort], sd=p[19+2*s.cohort];
      for(int i=0;i<N;++i){double a=(i+22-m)/sd;v[i]=std::exp(-0.5*a*a);}
    } else v=d.histogram[s.cohort];
    const double z=sum(v);
    for(double& x:v)x*=s.seed/z;
    return;
  }
  const auto& x=d.passages[s.parent].x[s.donor];
  const double z=sum(x);
  if(!(z>0))throw std::runtime_error("zero donor");
  for(int i=0;i<N;++i){
    s.base_pi[i]=x[i]/z;
    const double n=(i+22-88)/44.;
    v[i]=s.base_pi[i]*std::exp(s.tilt<0 ? 0 : p[22+2*s.tilt]*n+p[23+2*s.tilt]*n*n);
  }
  const double zz=sum(v);
  for(int i=0;i<N;++i){
    s.init_pi[i]=v[i]/zz;
    v[i]=s.seed*s.init_pi[i];
  }
}

int forward(Design& d,const double* p,bool gaussian,std::vector<Rates>& r,double* slopes,
            double* counts,double& supply_pen,double& min_supply) {
  std::fill(slopes,slopes+d.nslope,0.);
  supply_pen=0.; min_supply=1e30;
  for(size_t i=0;i<d.passages.size();++i){
    auto& s=d.passages[i];
    init(d,i,p,gaussian);
    double bh=0., bm=0.;
    size_t k=0;
    const Rates& o=r[s.oxygen];
    const double eh=p[14+2*s.cohort], em=p[15+2*s.cohort];
    for(int t=0;t<=s.last;++t) {
      const V& v=s.x[t];
      const double L=sum(v);
      if(s.demand[t]>0.){
        double gap=(L-s.demand[t])/s.demand[t];
        d.margins[s.demand_index[t]]=gap;
        min_supply=std::min(min_supply,gap);
        if(gap<0)supply_pen+=gap*gap;
      }
      if(k<s.slots.size() && s.slots[k].tick==t){
        const int row=s.slots[k].row;
        d.live[row]=L; d.visible_h[row]=bh; d.visible_m[row]=bm;
        double C=(L+(bh+bm))*scale;
        if(!std::isfinite(C))return 0;
        if(counts)counts[s.slots[k].row]=C;
        // All slopes use the original calendar-day OLS weights.
        slopes[s.slots[k].row]=std::log(std::max(1.,C));
        ++k;
      }
      if(t==s.last)break;
      const double dtc=dt*std::max(0.,1-L);
      V gx;
      mult(o,v,gx);
      for(int j=0;j<N;++j)s.x[t+1][j]=std::max(0.,v[j]*(1-dt*o.mu[j])+dtc*gx[j]);
      bh=decay*bh+dt*eh*dot(o.mu,v);
      bm=decay*bm+dtc*em*dot(o.b,v);
    }
  }
  return 1;
}

struct RateAdj {
  std::array<std::vector<double>,N> g;
  V mu{},b{};
  explicit RateAdj(const Rates& r){for(int i=0;i<N;++i)g[i].resize(r.g[i].size());}
};

void reverse(Design& d,const double* p,bool gaussian,const std::vector<Rates>& r,
             const double* source,double* grad,bool feasibility,const double* constraint_source=nullptr, const double* live_source=nullptr,
             const double* h_source=nullptr, const double* m_source=nullptr,
             const double* live_states_source=nullptr) {
  std::fill(grad,grad+NC,0.);
  std::vector<RateAdj> accum;
  accum.reserve(r.size());
  for(auto& op:r)accum.emplace_back(op);
  for(auto& s:d.passages){
    s.adjoint_last=-1;
    for(int t=0;t<=s.last;++t)if(s.demand_index[t]>=0)s.external[t].fill(0.);
  }
  for(int index=int(d.passages.size())-1;index>=0;--index){
    auto& s=d.passages[index];
    int end=s.last;
    if(!source && !feasibility) {
      end=s.adjoint_last;
      for(int t=0;t<=s.last;++t)if(s.demand_index[t]>=0 && constraint_source[s.demand_index[t]]!=0.)
        end=std::max(end,t);
      if(end<0)continue;
    }
    const Rates& o=r[s.oxygen];
    auto& q=accum[s.oxygen];
    const double eh=p[14+2*s.cohort],em=p[15+2*s.cohort];
    V xj{};
    double bjh=0., bjm=0.;
    int slot=int(s.slots.size())-1;
    while(slot>=0 && s.slots[slot].tick>end)--slot;
    for(int t=end;t>=0;--t){
      const V& x=s.x[t];
      if(slot>=0 && s.slots[slot].tick==t) {
        double g=source ? source[s.slots[slot].row] : 0.;
        const int row=s.slots[slot].row;
        bjh+=g+(h_source?h_source[row]:0.);
        bjm+=g+(m_source?m_source[row]:0.);
        const double lj=g+(live_source?live_source[row]:0.);
        for(int i=0;i<N;++i)
          xj[i]+=lj+(live_states_source?live_states_source[row*N+i]:0.);
        --slot;
      }
      if(feasibility && s.demand[t]>0.) {
        const double g=(sum(x)-s.demand[t])/s.demand[t];
        if(g<0)for(double& a:xj)a+=2*g/s.demand[t];
      }
      if(constraint_source && s.demand[t]>0.){
        const double g=constraint_source[s.demand_index[t]]/s.demand[t];
        for(double& a:xj)a+=g;
      }
      for(int i=0;i<N;++i)xj[i]+=s.external[t][i];
      if(t==0)break;
      const V& prev=s.x[t-1];
      for(int j=0;j<N;++j)if(x[j]==0.)xj[j]=0.;
      const double L=sum(prev), c=std::max(0.,1-L), dc=dt*c;
      double crowd_adj=bjm*em*dot(o.b,prev);
      grad[14+2*s.cohort] += dt*bjh*dot(o.mu,prev);
      grad[15+2*s.cohort] += dc*bjm*dot(o.b,prev);
      V tmp;
      for(int i=0;i<N;++i){
        const double a=prev[i], sw=bjm*em;
        q.mu[i]+=dt*a*(bjh*eh-xj[i]);
        q.b[i]+=dc*sw*a;
        double back=0.;
        for(size_t j=0;j<o.g[i].size();++j) {
          const auto& e=o.g[i][j];
          double ad=xj[e.to];
          back+=e.v*ad;
          q.g[i][j]+=dc*ad*a;
        }
        crowd_adj+=a*back;
        tmp[i]=xj[i]*(1-dt*o.mu[i])+dc*back+bjh*dt*eh*o.mu[i]+bjm*dc*em*o.b[i];
      }
      for(int i=0;i<N;++i)xj[i]=tmp[i]-(c>0 ? dt*crowd_adj : 0.);
      bjh*=decay; bjm*=decay;
    }
    if(s.parent<0){
      if(gaussian){
        double m=p[18+2*s.cohort], sd=p[19+2*s.cohort], Lj=dot(xj,s.x[0])/s.seed;
        for(int j=0;j<N;++j){
          double z=(j+22-m)/sd;
          double qj=s.x[0][j]*(xj[j]-Lj)/sd;
          grad[18+2*s.cohort]+=qj*z;
          grad[19+2*s.cohort]+=qj*z*z;
        }
      }
    }else {
      const auto& parent_x=d.passages[s.parent].x[s.donor];
      d.passages[s.parent].adjoint_last=std::max(d.passages[s.parent].adjoint_last,s.donor);
      V& parent_j=d.passages[s.parent].external[s.donor];
      const double L=sum(parent_x);
      const double jmean=dot(xj,s.init_pi);
      double factor=0.,mean_grad=0.,variance_grad=0.;
      if(s.tilt>=0){
        for(int j=0;j<N;++j){double z=(j+22-88)/44.;factor+=s.base_pi[j]*std::exp(p[22+2*s.tilt]*z+p[23+2*s.tilt]*z*z);}
      }
      for(int j=0;j<N;++j){
        const double z=(j+22-88)/44.;
        const double j0=s.seed*(xj[j]-jmean);
        const double f=s.tilt<0 ? 1. : std::exp(p[22+2*s.tilt]*z+p[23+2*s.tilt]*z*z)/factor;
        parent_j[j]+=j0*f/L;
        if(s.tilt>=0){mean_grad+=j0*s.init_pi[j]*z;variance_grad+=j0*s.init_pi[j]*z*z;}
      }
      if(s.tilt>=0){grad[22+2*s.tilt]+=mean_grad;grad[23+2*s.tilt]+=variance_grad;}
    }
  }
  for(size_t op=0;op<r.size();++op)for(int i=0;i<N;++i){
    const auto& v=r[op];
    const auto& a=accum[op];
    for(int p=0;p<M;++p)grad[p]+=v.dm[i][p]*a.mu[i]+v.db[i][p]*a.b[i];
    for(size_t j=0;j<v.g[i].size();++j)for(int p=0;p<M;++p)grad[p]+=v.dg[i][j][p]*a.g[i][j];
  }
}
int prepare(Design& d,const double* p,bool gaussian){
  if(d.cache_valid && d.cached_gaussian==gaussian && std::equal(d.cached_p.begin(),d.cached_p.end(),p))
    return 1;
  d.cache_valid=false;
  d.cached_rates.clear();
  for(double x:d.oxygens)d.cached_rates.emplace_back(build(x,p));
  double penalty=0., margin=0.;
  int ok=forward(d,p,gaussian,d.cached_rates,d.logc.data(),d.cells.data(),penalty,margin);
  if(ok){
    std::copy(p,p+NC,d.cached_p.begin());
    d.cached_gaussian=gaussian;
    d.cache_valid=true;
  }
  return ok;
}
}

extern "C" {
void* ols_design(int ns,int no,int ng,const double* oxygens,const int* passage_meta,
                 const double* amounts,int nk,const int* sample_meta,const double* weights,
                 const double* hist) {
  try {
    auto d=std::make_unique<Design>();
    d->oxygens.assign(oxygens,oxygens+no);
    d->passages.resize(ns);
    d->nslope=nk;
    d->nobs=ng;
    d->group.resize(nk);
    d->weights.assign(weights,weights+nk);
    for(int i=0;i<ns;++i){
      auto& p=d->passages[i];
      const int* m=passage_meta+i*6;
      p.oxygen=m[0];p.cohort=m[1];p.parent=m[2];p.donor=m[3];p.last=m[4];p.tilt=m[5];
      p.seed=amounts[i]/scale;
      p.demand.resize(p.last+1);
      p.demand_index.assign(p.last+1,-1);
      p.x.resize(p.last+1);
      p.external.resize(p.last+1);
      if(p.parent>=0)d->passages[p.parent].demand[p.donor]+=p.seed;
    }
    int nc=0;
    for(auto& p:d->passages)for(int tick=0;tick<=p.last;++tick)
      if(p.demand[tick]>0.)p.demand_index[tick]=nc++;
    d->margins.resize(nc);
    d->logc.resize(nk);
    d->cells.resize(nk);
    d->live.resize(nk); d->visible_h.resize(nk); d->visible_m.resize(nk);
    for(int k=0;k<nk;++k){
      const int* m=sample_meta+k*3;
      d->passages[m[0]].slots.push_back({m[1],k,weights[k]});
      d->group[k]=m[2];
    }
    for(int i=0;i<2;++i)std::copy(hist+i*N,hist+(i+1)*N,d->histogram[i].begin());
    return d.release();
  } catch (...) {return nullptr;}
}

void ols_destroy(void* handle){delete static_cast<Design*>(handle);}

int ols_eval(void* handle,const double* p,int gaussian,const double* source_raw,
             double* logc,double* cells,double* grad,double* diagnostics,unsigned long long* mask) {
  try {
    auto& d=*static_cast<Design*>(handle);
    d.cache_valid=false;
    std::vector<Rates> r;
    r.reserve(d.oxygens.size());
    for(double x:d.oxygens)r.emplace_back(build(x,p));
    if(mask)for(size_t j=0;j<r.size();++j)mask[j]=r[j].mask;
    double pen=0,min_supply=0;
    if(!forward(d,p,gaussian!=0,r,logc,cells,pen,min_supply))return 0;
    diagnostics[0]=pen;diagnostics[1]=min_supply;
    if(grad){
      if(pen>0){
        std::vector<double> zero(d.nslope,0.);
        reverse(d,p,gaussian!=0,r,zero.data(),grad,true);
      }else if(source_raw){
        std::vector<double> source(d.nslope);
        for(int i=0;i<d.nslope;++i)source[i]=(cells[i]>1.?source_raw[i]*scale/cells[i]:0.);
        reverse(d,p,gaussian!=0,r,source.data(),grad,false);
      }
    }
    return 1;
  }catch(...){return 0;}
}

int ols_loss(void* handle,const double* p,int gaussian,const double* target,const double* precision,
             int dense,double* output,double* grad) {
  try {
    auto& d=*static_cast<Design*>(handle);
    d.cache_valid=false;
    std::vector<Rates> r;
    r.reserve(d.oxygens.size());
    for(double x:d.oxygens)r.emplace_back(build(x,p));
    std::vector<double> logc(d.nslope), cells(d.nslope), source(d.nslope), f(d.nobs), a(d.nobs);
    if(!forward(d,p,gaussian!=0,r,logc.data(),cells.data(),output[1],output[2]))return 0;
    if(output[1]>0){
      reverse(d,p,gaussian!=0,r,source.data(),grad,true);
      output[0]=100. + 100.*output[1];
      for(int j=0;j<NC;++j)grad[j]*=100.;
      return 1;
    }
    for(int i=0;i<d.nslope;++i)f[d.group[i]]+=logc[i]*d.weights[i];
    for(int i=0;i<d.nobs;++i)f[i]-=target[i];
    output[0]=0.;
    for(int i=0;i<d.nobs;++i){
      if(dense){
        for(int j=0;j<d.nobs;++j)a[i]+=precision[i*d.nobs+j]*f[j];
      }else a[i]=precision[i]*f[i];
      output[0]+=0.5*a[i]*f[i];
    }
    for(int i=0;i<d.nslope;++i)if(cells[i]>1.)source[i]=scale/cells[i]*d.weights[i]*a[d.group[i]];
    reverse(d,p,gaussian!=0,r,source.data(),grad,false);
    return 1;
  }catch(...){return 0;}
}
int ols_nconstraints(void* handle){return static_cast<Design*>(handle)->margins.size();}

int ols_constrained_loss(void* handle,const double* p,int gaussian,const double* target,
                         const double* precision,int dense,double* f,double* grad){
  try{
    auto& d=*static_cast<Design*>(handle);
    if(!prepare(d,p,gaussian!=0))return 0;
    std::vector<double> source(d.nslope), residual(d.nobs), a(d.nobs);
    for(int i=0;i<d.nslope;++i)residual[d.group[i]]+=d.logc[i]*d.weights[i];
    for(int i=0;i<d.nobs;++i)residual[i]-=target[i];
    f[0]=0.;
    for(int i=0;i<d.nobs;++i) {
      if(dense){
        for(int j=0;j<d.nobs;++j)a[i]+=precision[i*d.nobs+j]*residual[j];
      }else a[i]=precision[i]*residual[i];
      f[0]+=0.5*a[i]*residual[i];
    }
    f[1]=*std::min_element(d.margins.begin(),d.margins.end());
    if(grad){
      for(int i=0;i<d.nslope;++i)if(d.cells[i]>1.)source[i]=scale/d.cells[i]*d.weights[i]*a[d.group[i]];
      reverse(d,p,gaussian!=0,d.cached_rates,source.data(),grad,false);
    }
    return 1;
  }catch(...){return 0;}
}

int ols_supply(void* handle,const double* p,int gaussian,double* margins,int ng,const int* pick,double* grad){
  try{
    auto& d=*static_cast<Design*>(handle);
    if(!prepare(d,p,gaussian!=0))return 0;
    std::copy(d.margins.begin(),d.margins.end(),margins);
    if(grad){
      std::vector<double> source(d.margins.size());
      for(int i=0;i<ng;++i){
        source[pick[i]]=1.;
        reverse(d,p,gaussian!=0,d.cached_rates,nullptr,grad+i*NC,false,source.data());
        source[pick[i]]=0.;
      }
    }
    return 1;
  }catch(...){return 0;}
}
}


extern "C" int count_loss(void* handle, const double* p, int gaussian, int n,
                          const int* slots, const double* target, const double* precision,
                          int dense, double* out, double* grad) {
  try {
    auto& d = *static_cast<Design*>(handle);
    if (!prepare(d, p, gaussian != 0)) return 0;
    std::vector<double> residual(n), weighted(n), source(d.nslope);
    for (int j = 0; j < n; ++j) residual[j] = d.logc[slots[j]] - target[j];
    out[0] = 0.;
    for (int j = 0; j < n; ++j) {
      if (dense) {
        for (int k = 0; k < n; ++k) weighted[j] += precision[j*n+k]*residual[k];
      } else {
        weighted[j] = precision[j]*residual[j];
      }
      out[0] += .5*weighted[j]*residual[j];
      const int slot = slots[j];
      if (grad && d.cells[slot] > 1.) source[slot] = scale/d.cells[slot]*weighted[j];
    }
    out[1] = *std::min_element(d.margins.begin(), d.margins.end());
    if (grad) reverse(d, p, gaussian != 0, d.cached_rates, source.data(), grad, false);
    return 1;
  } catch (...) {
    return 0;
  }
}


extern "C" int ar1_count_loss(void* handle, const double* p, int n, const int* slots,
                              const int* previous, const int* gap, const double* target,
                              double* out, double* grad) {
  try {
    auto& d = *static_cast<Design*>(handle);
    if (!prepare(d, p, false)) return 0;
    const double rho = p[NC];
    std::vector<double> residual(n), weighted(n), source(d.nslope);
    double q = 0., logdet = 0., dq = 0., dlogdet = 0.;
    for (int i = 0; i < n; ++i) residual[i] = d.logc[slots[i]] - target[i];
    for (int i = 0; i < n; ++i) {
      const int j = previous[i];
      if (j < 0) {
        q += residual[i] * residual[i];
        weighted[i] += residual[i];
      } else {
        const double a = std::pow(rho, gap[i]);
        const double da = gap[i] == 1 ? 1. : gap[i] * std::pow(rho, gap[i]-1);
        const double v = 1.-a*a;
        const double e = residual[i] - a*residual[j];
        const double z = e/v;
        q += e*z;
        logdet += std::log1p(-a*a);
        dq += 2.*z*da*(-residual[j]+a*z);
        dlogdet -= 2.*a*da/v;
        weighted[i] += z;
        weighted[j] -= a*z;
      }
    }
    const double factor = std::exp(logdet/n);
    out[0] = .5*q*factor;
    out[1] = *std::min_element(d.margins.begin(),d.margins.end());
    out[2] = q;
    out[3] = logdet;
    if (grad) {
      for (int i = 0; i < n; ++i) {
        const int j = slots[i];
        if (d.cells[j]>1.) source[j] = factor * weighted[i]*scale/d.cells[j];
      }
      reverse(d,p,false,d.cached_rates,source.data(),grad,false);
      grad[NC] = .5*factor*(dq+q*dlogdet/n);
    }
    return std::isfinite(out[0]);
  } catch (...) {
    return 0;
  }
}

extern "C" int extension_keep_original_supply(void* handle, int last) {
  auto& d=*static_cast<Design*>(handle);
  for (int i=int(d.passages.size())-last;i<int(d.passages.size());++i) {
    const auto& branch=d.passages.at(i);
    auto& demand=d.passages.at(branch.parent).demand.at(branch.donor);
    demand-=branch.seed;
    if (!(demand>0.)) return 0;
  }
  d.cache_valid=false;
  return 1;
}
extern "C" int extension_readouts(void* handle, const double* p, int nk, const int* kary_rows,
                                  double* channels, double* kary_states) {
  try {
    auto& d=*static_cast<Design*>(handle);
    if (!prepare(d,p,false)) return 0;
    std::copy(d.live.begin(),d.live.end(),channels);
    std::copy(d.visible_h.begin(),d.visible_h.end(),channels+d.nslope);
    std::copy(d.visible_m.begin(),d.visible_m.end(),channels+2*d.nslope);
    for (auto const& seg:d.passages) for (auto const& s:seg.slots) for(int j=0;j<nk;++j)
      if (s.row==kary_rows[j]) std::copy(seg.x[s.tick].begin(),seg.x[s.tick].end(),kary_states+j*N);
    return 1;
  } catch (...) { return 0; }
}
extern "C" int extension_deviance_adjoint(
    void* handle, const double* p, int nk, const int* kary_rows,
    const double* channels, const double* kary_source,
    int n, const int* slots, const int* previous, const int* gaps, const double* target,
    double quadratic, double* grad) {
  try {
    auto& d=*static_cast<Design*>(handle);
    if (!prepare(d,p,false)) return 0;
    std::vector<double> residual(n), weight(n), count_source(d.nslope);
    const double rho=p[NC];
    double dq=0., dlogdet=0.;
    for(int i=0;i<n;++i) residual[i]=d.logc[slots[i]]-target[i];
    for(int i=0;i<n;++i) {
      const int j=previous[i];
      if(j<0) { weight[i]+=residual[i]; continue; }
      const double a=std::pow(rho,gaps[i]);
      const double da=gaps[i]==1?1.:gaps[i]*std::pow(rho,gaps[i]-1);
      const double v=1.-a*a;
      const double r=(residual[i]-a*residual[j])/v;
      weight[i]+=r;
      weight[j]-=a*r;
      dq+=2.*r*da*(-residual[j]+a*r);
      dlogdet-=2.*a*da/v;
    }
    for(int i=0;i<n;++i)
      if(d.cells[slots[i]]>1.) count_source[slots[i]]=2.*n/quadratic*weight[i]*scale/d.cells[slots[i]];
    std::vector<double> states(d.nslope*N);
    for(int j=0;j<nk;++j)
      std::copy_n(kary_source+j*N,N,states.begin()+kary_rows[j]*N);
    reverse(d,p,false,d.cached_rates,count_source.data(),grad,false,nullptr,
            channels,channels+d.nslope,channels+2*d.nslope,states.data());
    grad[NC]=n*dq/quadratic+dlogdet;
    return 1;
  } catch (...) { return 0; }
}
