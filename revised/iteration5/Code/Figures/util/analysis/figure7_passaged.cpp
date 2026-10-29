#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Figure7-only Euler instrumentation. All biological coefficients come from
// the locked original model. R chooses native passage times and next inocula.
// X is normalized only for numerical scaling; log_cells retains absolute size.
// D: surviving first-generation H->L daughters, label NOT inherited.
// U: cells whose lineage has never crossed H->L. Passive labels do not alter X.
// [[Rcpp::export(rng = false)]]
List f7p_segment_cpp(S4 G,S4 KH,S4 GU,NumericVector lambda,NumericVector mu,
    NumericVector hl,NumericVector ll,NumericVector upper,NumericVector lower,
    NumericVector ploidy,LogicalVector low,NumericVector initial,
    NumericVector direct,NumericVector never,int days,double dt,
    double initial_cells,bool crowding,double K) {
  int n=initial.size(),steps=std::lround(1/dt);
  if(n<1||days<0||steps<1||std::abs(steps*dt-1)>1e-12||!(initial_cells>0))
    stop("Invalid Figure7 segment inputs");
  IntegerVector gi=G.slot("i"),gp=G.slot("p"),hi=KH.slot("i"),hp=KH.slot("p"),
    ui=GU.slot("i"),up=GU.slot("p");
  NumericVector gv=G.slot("x"),hv=KH.slot("x"),uv=GU.slot("x");
  std::vector<double>x(initial.begin(),initial.end()),d(direct.begin(),direct.end()),
    u(never.begin(),never.end()),gx(n),gu(n),kx(n);
  NumericMatrix X(n,days+1),D(n,days+1),U(n,days+1),metrics(days+1,32);
  colnames(metrics)=CharacterVector::create("log_cells","mean_ploidy","low_fraction",
    "A","B","HL_per_start_cell","LL_per_start_cell","upper_loss_fraction",
    "direct_share","never_share","all_ploidy","direct_ploidy","other_ploidy",
    "historical_ploidy","never_ploidy","all_division_rate","direct_division_rate",
    "other_division_rate","historical_division_rate","never_division_rate",
    "all_low_birth_rate","direct_low_birth_rate","other_low_birth_rate",
    "historical_low_birth_rate","never_low_birth_rate","direct_birth_share",
    "historical_birth_share","crowding","net_rate","upper_lost_births_per_cell",
    "retained_births_per_cell","lower_lost_births_per_cell");
  double logn=std::log(initial_cells),max_subset=0;int corrections=0;
  auto cfactor=[&](){return crowding?std::max(0.0,1-std::exp(logn-std::log(K))):1.0;};
  auto record=[&](int day,double H,double L) {
    double cf=cfactor(),mp=0,lm=0,ul=0,lo=0,birth=0,net=0;
    double mass[5]={0},pm[5]={0},div[5]={0},br[5]={0};
    for(int j=0;j<n;j++) {
      X(j,day)=x[j];D(j,day)=d[j];U(j,day)=u[j];mp+=ploidy[j]*x[j];
      ul+=upper[j]*x[j];lo+=lower[j]*x[j];
      double gs=0;for(int q=gp[j];q<gp[j+1];q++)gs+=gv[q];
      birth+=(gs+lambda[j])*x[j];net+=(cf*gs-mu[j])*x[j];
      if(low[j]) {
        lm+=x[j];double group[5]={x[j],d[j],x[j]-d[j],x[j]-u[j],u[j]};
        for(int g=0;g<5;g++){mass[g]+=group[g];pm[g]+=ploidy[j]*group[g];
          div[g]+=cf*lambda[j]*group[g];br[g]+=cf*ll[j]*group[g];}
      }
    }
    metrics(day,0)=logn;metrics(day,1)=mp;metrics(day,2)=lm;
    metrics(day,3)=day&&H+L>0?100*H/(H+L):NA_REAL;
    metrics(day,4)=day&&H+L>0?100*L/(H+L):NA_REAL;
    metrics(day,5)=day?H:NA_REAL;metrics(day,6)=day?L:NA_REAL;
    metrics(day,7)=birth+ul+lo>0?ul/(birth+ul+lo):NA_REAL;
    metrics(day,8)=lm>0?100*mass[1]/lm:NA_REAL;
    metrics(day,9)=lm>0?100*mass[4]/lm:NA_REAL;
    for(int g=0;g<5;g++) {
      metrics(day,10+g)=mass[g]>0?pm[g]/mass[g]:NA_REAL;
      metrics(day,15+g)=mass[g]>0?div[g]/mass[g]:NA_REAL;
      metrics(day,20+g)=mass[g]>0?br[g]/mass[g]:NA_REAL;
    }
    metrics(day,25)=br[0]>0?100*br[1]/br[0]:NA_REAL;
    metrics(day,26)=br[0]>0?100*br[3]/br[0]:NA_REAL;
    metrics(day,27)=cf;metrics(day,28)=net;metrics(day,29)=cf*ul;
    metrics(day,30)=cf*birth;metrics(day,31)=cf*lo;
  };
  record(0,0,0);
  for(int day=1;day<=days;day++) {
    double start=logn,H=0,L=0;
    for(int step=0;step<steps;step++) {
      std::fill(gx.begin(),gx.end(),0);std::fill(gu.begin(),gu.end(),0);
      std::fill(kx.begin(),kx.end(),0);
      for(int j=0;j<n;j++) {
        for(int q=gp[j];q<gp[j+1];q++)gx[gi[q]]+=gv[q]*x[j];
        for(int q=up[j];q<up[j+1];q++)gu[ui[q]]+=uv[q]*u[j];
        for(int q=hp[j];q<hp[j+1];q++)kx[hi[q]]+=hv[q]*x[j];
      }
      double cf=cfactor(),scale=std::exp(logn-start),total=0;
      for(int j=0;j<n;j++) {
        H+=dt*cf*hl[j]*x[j]*scale;L+=dt*cf*ll[j]*x[j]*scale;
        double nx=x[j]+dt*(cf*gx[j]-mu[j]*x[j]);
        double nd=d[j]+dt*(cf*kx[j]-(cf*lambda[j]+mu[j])*d[j]);
        // Form G-KH before multiplication: avoid cancellation of two tiny
        // products for a vanishing never-HL group. This is the same equation.
        double nu=u[j]+dt*(cf*gu[j]-mu[j]*u[j]);
        if(!R_finite(nx)||!R_finite(nd)||!R_finite(nu)||nx< -1e-12||nd< -1e-12||nu< -1e-12)
          stop("Invalid Figure7 numerical state");
        max_subset=std::max(max_subset,std::max(nd-nx,nu-nx));
        if(nd>nx+1e-10||nu>nx+1e-10||(!low[j]&&nd>1e-12))stop("Invalid passive source subset");
        if(nx<0||nd<0||nu<0||nd>nx||nu>nx)corrections++;
        x[j]=std::max(0.0,nx);d[j]=std::max(0.0,std::min(nd,x[j]));
        u[j]=std::max(0.0,std::min(nu,x[j]));total+=x[j];
      }
      if(!(total>0)||!R_finite(total))stop("Numerical failure: no positive finite mass");
      logn+=std::log(total);
      for(int j=0;j<n;j++){x[j]/=total;d[j]/=total;u[j]/=total;}
    }
    record(day,H,L);checkUserInterrupt();
  }
  return List::create(_["X"]=X,_["D"]=D,_["U"]=U,_["metrics"]=metrics,
    _["subset_error"]=max_subset,_["roundoff_corrections"]=corrections);
}
