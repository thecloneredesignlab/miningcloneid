"""Extend the frozen A7K kernel with four reset tilts and exact TV gradients."""
import os
from pathlib import Path
import platform
import subprocess

HERE=Path(__file__).resolve().parent
ORIGINAL=HERE.parent/"20261119_experimental_death_A7K_comparison_01/extended_kernel.cpp"
source=ORIGINAL.read_text()


def replace(old,new):
    global source
    assert source.count(old)==1,old
    source=source.replace(old,new)


replace("NC = 46;", "NC = 54;")
replace("if(s.demand_index[t]>=0 && constraint_source[s.demand_index[t]]!=0.)",
        "if(constraint_source && s.demand_index[t]>=0 && constraint_source[s.demand_index[t]]!=0.)")
replace("const double* live_states_source=nullptr) {",
        """const double* live_states_source=nullptr,
             const int* direct_pairs=nullptr,const double* direct_weights=nullptr,
             int ndirect=0) {""")
replace("  for(int index=int(d.passages.size())-1;",
        """  for(int j=0;j<ndirect;++j) {
    auto& seg=d.passages[direct_pairs[2*j]];
    const int t=direct_pairs[2*j+1];
    seg.adjoint_last=std::max(seg.adjoint_last,t);
    for(int i=0;i<N;++i) seg.external[t][i]+=direct_weights[j*N+i];
  }
  for(int index=int(d.passages.size())-1;""")
replace("  for(size_t op=0;op<r.size();++op)for(int i=0;i<N;++i){",
        """  for(int j=0;j<ndirect;++j)
    d.passages[direct_pairs[2*j]].external[direct_pairs[2*j+1]].fill(0.);
  for(size_t op=0;op<r.size();++op)for(int i=0;i<N;++i){""")

source+="""
extern "C" int extension_enable_a7k_tilts(void* handle) {
  auto& d=*static_cast<Design*>(handle);
  if (d.passages.size()!=120) return 0;
  for(int i=0;i<4;++i) {
    auto& branch=d.passages[116+i];
    if(branch.last!=80 || branch.donor!=100 || branch.parent<0) return 0;
    branch.tilt=12+i;
  }
  d.cache_valid=false;
  return 1;
}

extern "C" int extension_a7k_composition(void* handle,const double* p,double* tv,
                                         double* base,double* after,double* seeds) {
  try {
    auto& d=*static_cast<Design*>(handle);
    if(!prepare(d,p,false)) return 0;
    for(int j=0;j<4;++j) {
      const auto& branch=d.passages[116+j];
      tv[j]=0.;
      for(int i=0;i<N;++i) tv[j]+=.5*std::abs(branch.init_pi[i]-branch.base_pi[i]);
      if(base) std::copy(branch.base_pi.begin(),branch.base_pi.end(),base+j*N);
      if(after) std::copy(branch.init_pi.begin(),branch.init_pi.end(),after+j*N);
      if(seeds) seeds[j]=sum(branch.x[0])*scale;
    }
    return 1;
  } catch(...) { return 0; }
}

extern "C" int extension_tv_adjoint(void* handle,const double* p,int nt,
                                    const int* picks,double* grad) {
  try {
    auto& d=*static_cast<Design*>(handle);
    if(!prepare(d,p,false)) return 0;
    for(int j=0;j<nt;++j) {
      const int branch_id=116+picks[j];
      const auto& branch=d.passages.at(branch_id);
      const auto& parent=d.passages.at(branch.parent);
      const double parent_mass=sum(parent.x.at(branch.donor));
      std::array<double,2*N> force{};
      double mean=0.;
      for(int i=0;i<N;++i) {
        const double diff=branch.init_pi[i]-branch.base_pi[i];
        const double g=.5*((diff>1e-16)-(diff< -1e-16));
        force[i]=g/branch.seed;
        force[N+i]=-g;
        mean+=g*branch.base_pi[i];
      }
      for(int i=0;i<N;++i) force[N+i]=(force[N+i]+mean)/parent_mass;
      const int pairs[4]={branch_id,0,branch.parent,branch.donor};
      reverse(d,p,false,d.cached_rates,nullptr,grad+j*NC,false,nullptr,
              nullptr,nullptr,nullptr,nullptr,pairs,force.data(),2);
    }
    return 1;
  } catch(...) { return 0; }
}
"""
cpp=HERE/"tilt_kernel.cpp"
cpp.write_text(source)
fallback="/private/tmp/culture_a7k_tv_20261119.dylib"
lib=Path(os.environ.get("CULTURE_A7K_TV_LIB",fallback))
tmp=lib.with_suffix(".next"+lib.suffix)
flags=["-dynamiclib"] if platform.system()=="Darwin" else ["-shared","-fPIC"]
subprocess.run(["c++","-std=c++17","-O3",*flags,str(cpp),"-o",str(tmp)],check=True)
tmp.replace(lib)
