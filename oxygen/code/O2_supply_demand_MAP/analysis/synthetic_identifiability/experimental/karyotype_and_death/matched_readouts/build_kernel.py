"""Extend the frozen propagation and adjoint with two visible-object channels."""
from pathlib import Path
import os
import platform
import subprocess

HERE = Path(__file__).resolve().parent
parent = HERE.parent
base = parent/"20261118_experimental_OLS_practical_01/ols_kernel.cpp"
text = base.read_text()


def replace(old, new):
    global text
    assert text.count(old)==1, old
    text = text.replace(old,new)


replace("std::vector<double> logc, cells, margins;",
        "std::vector<double> logc, cells, margins, live, visible_h, visible_m;")
replace("    double b=0.;", "    double bh=0., bm=0.;")
replace("        double C=(L+b)*scale;",
        """        const int row=s.slots[k].row;
        d.live[row]=L; d.visible_h[row]=bh; d.visible_m[row]=bm;
        double C=(L+(bh+bm))*scale;""")
replace("      b=decay*b+dt*eh*dot(o.mu,v)+dtc*em*dot(o.b,v);",
        "      bh=decay*bh+dt*eh*dot(o.mu,v);\n      bm=decay*bm+dtc*em*dot(o.b,v);")
replace("const double* constraint_source=nullptr) {",
        """const double* constraint_source=nullptr, const double* live_source=nullptr,
             const double* h_source=nullptr, const double* m_source=nullptr,
             const double* live_states_source=nullptr) {""")
replace("    double bj=0.;", "    double bjh=0., bjm=0.;")
replace("        bj+=g;\n        for(double& a:xj)a+=g;",
        """        const int row=s.slots[slot].row;
        bjh+=g+(h_source?h_source[row]:0.);
        bjm+=g+(m_source?m_source[row]:0.);
        const double lj=g+(live_source?live_source[row]:0.);
        for(int i=0;i<N;++i)
          xj[i]+=lj+(live_states_source?live_states_source[row*N+i]:0.);""")
replace("      double crowd_adj=bj*em*dot(o.b,prev);","      double crowd_adj=bjm*em*dot(o.b,prev);")
replace(" += dt*bj*dot"," += dt*bjh*dot")
replace(" += dc*bj*dot"," += dc*bjm*dot")
replace("const double a=prev[i], sw=bj*em;", "const double a=prev[i], sw=bjm*em;")
replace("(bj*eh-xj[i])","(bjh*eh-xj[i])")
replace("+bj*(dt*eh*o.mu[i]+dc*em*o.b[i]);","+bjh*dt*eh*o.mu[i]+bjm*dc*em*o.b[i];")
replace("      bj*=decay;","      bjh*=decay; bjm*=decay;")
replace("    d->cells.resize(nk);","    d->cells.resize(nk);\n    d->live.resize(nk); d->visible_h.resize(nk); d->visible_m.resize(nk);")

sources=[(parent/"20261118_experimental_counts_practical_01/count_kernel.cpp").read_text(),
         (parent/"20261118_experimental_counts_ar1_practical_01/ar1_kernel.cpp").read_text()]
for src in sources:
    text += "\n" + src.split("\n",1)[1]
text += """
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
"""
path=HERE/"extended_kernel.cpp"
path.write_text(text)
library=Path(os.environ.get("CULTURE_DEATH_A7K_LIB","/private/tmp/culture_death_A7K_20261119.dylib"))
temporary=library.with_suffix(".next"+library.suffix)
flags=["-dynamiclib"] if platform.system()=="Darwin" else ["-shared","-fPIC"]
subprocess.run(["c++","-std=c++17","-O3",*flags,str(path),"-o",
                str(temporary)], check=True)
temporary.replace(library)
