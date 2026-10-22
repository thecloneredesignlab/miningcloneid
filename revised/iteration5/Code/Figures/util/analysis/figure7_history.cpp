#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

namespace {
void product(const NumericMatrix& A,const NumericMatrix& x,NumericMatrix& y) {
  int n=A.nrow(),k=x.ncol();const char tr='N';const double one=1,zero=0;
  F77_CALL(dgemm)(&tr,&tr,&n,&k,&n,&one,A.begin(),&n,x.begin(),&n,&zero,y.begin(),&n FCONE FCONE);
}
double crowd(double logn,bool on,int mode,double K) {
  if(!on)return 1;
  double ratio=std::exp(logn-std::log(K));
  return mode==0?std::max(0.0,1-ratio):std::exp(-ratio);
}
double error(double a,double b) {
  if(!R_finite(a)||!R_finite(b))return R_finite(a)==R_finite(b)?0:R_PosInf;
  return std::abs(a-b)/(1+std::abs(b));
}
}

// Telemetry-only replay of authenticated operators. No RNG or new passage decisions.
// Culture follows the saved fitted-DT Euler update and saved X/U sampling anchors.
// Tumor follows the original one-day exp(M) and exp(M-K_HL) operators.
// [[Rcpp::export(rng = false)]]
List f7_history_cpp(List op,List reference,NumericVector initial,LogicalVector low,
    NumericVector ploidy,NumericVector hl,NumericVector ll,bool culture,double dt,
    bool crowding,int mode,double K,double seed_cells,int days) {
  NumericMatrix G=op["G"],KH=op["KH"],P=op["P"],PU=op["PU"];
  NumericVector mu=op["mu"];
  NumericMatrix rf=reference["flux_percent"],rm=reference["mean_ploidy"],rt=reference["tag_fraction"];
  NumericMatrix endx=reference["final_state"],endu=reference["final_untagged_state"];
  DataFrame events=as<DataFrame>(reference["passage_events"]);
  IntegerVector eday=events["day"],erep=events["repeat_id"];
  NumericMatrix prex=reference["pre_passage_state"],preu=reference["pre_passage_untagged_state"],
    postx=reference["post_passage_state"],postu=reference["post_passage_untagged_state"];
  int n=initial.size(),k=rf.nrow(),substeps=culture?std::lround(1/dt):1,ne=eday.size(),event=0;
  if(n<1||G.nrow()!=n||G.ncol()!=n||KH.nrow()!=n||KH.ncol()!=n||mu.size()!=n||low.size()!=n||
    hl.size()!=n||ll.size()!=n||ploidy.size()!=n||rf.ncol()!=days+1||rm.nrow()!=k||rm.ncol()!=days+1||
    rt.nrow()!=k||rt.ncol()!=days+1||endx.nrow()!=n||endx.ncol()!=k||endu.nrow()!=n||endu.ncol()!=k||
    substeps<1||(culture&&std::abs(substeps*dt-1)>1e-12)||erep.size()!=ne||
    prex.ncol()!=ne||preu.ncol()!=ne||postx.ncol()!=ne||postu.ncol()!=ne)
    stop("Invalid authenticated replay dimensions");
  if(!culture&&(ne!=0||P.nrow()!=n||P.ncol()!=n||PU.nrow()!=n||PU.ncol()!=n))
    stop("Continuous tumor replay cannot include culture passages");
  NumericMatrix x(n,k),u(n,k),gx(n,k),gu(n,k),AU=clone(G);
  for(int r=0;r<n*n;++r)AU[r]-=KH[r];
  for(int c=0;c<k;++c)for(int r=0;r<n;++r)x(r,c)=u(r,c)=initial[r];
  NumericVector logn(k,std::log(seed_cells));
  NumericMatrix shares(k,days+1),growth(k,days+1),lowmass(k,days+1),umass(k,days+1);
  double metric_error=0,anchor_error=0,balance_error=0,final_error=0;int clipping=0;
  auto record=[&](int c,int day) {
    double cf=culture?crowd(logn[c],crowding,mode,K):1,m=0,ul=0,H=0,L=0,TL=0,mp=0,tag=0,total=0;
    for(int r=0;r<n;++r) {
      double t=x(r,c)-u(r,c);total+=x(r,c);mp+=ploidy[r]*x(r,c);tag+=t;
      H+=cf*hl[r]*x(r,c);L+=cf*ll[r]*x(r,c);TL+=cf*ll[r]*t;
      if(low[r]){m+=x(r,c);ul+=u(r,c);}
    }
    shares(c,day)=m>1e-14?100*(m-ul)/m:NA_REAL;
    growth(c,day)=culture?(L>1e-14?100*TL/L:NA_REAL):
      (m>1e-14&&L/m>1e-8?100*TL/L:NA_REAL);
    lowmass(c,day)=m;umass(c,day)=ul;
    metric_error=std::max(metric_error,error(H+L>1e-14?100*H/(H+L):NA_REAL,rf(c,day)));
    metric_error=std::max(metric_error,error(mp,rm(c,day)));
    metric_error=std::max(metric_error,error(tag,rt(c,day)));
    balance_error=std::max(balance_error,std::abs(total-1));
    if(R_finite(shares(c,day))&&(shares(c,day)<-1e-8||shares(c,day)>100+1e-8))stop("Lineage stock share outside [0,100]");
    if(R_finite(growth(c,day))&&(growth(c,day)<-1e-8||growth(c,day)>100+1e-8))stop("Lineage proliferation share outside [0,100]");
  };
  for(int c=0;c<k;++c)record(c,0);
  for(int day=1;day<=days;++day) {
    for(int s=0;s<substeps;++s) {
      product(culture?G:P,x,gx);product(culture?AU:PU,u,gu);
      for(int c=0;c<k;++c) {
        double cf=culture?crowd(logn[c],crowding,mode,K):1,total=0;
        for(int r=0;r<n;++r) {
          double v=culture?x(r,c)+dt*(cf*gx(r,c)-mu[r]*x(r,c)):gx(r,c);
          double w=culture?u(r,c)+dt*(cf*gu(r,c)-mu[r]*u(r,c)):gu(r,c);
          if(v<0){if(!culture&&v<-1e-10)stop("Negative tumor state");v=0;clipping++;}
          if(w<0){if(!culture&&w<-1e-10)stop("Negative tumor untagged state");w=0;clipping++;}
          if(!R_finite(v)||!R_finite(w))stop("Nonfinite replay state");
          balance_error=std::max(balance_error,std::max(0.0,w-v));
          if(w>v+1e-9)stop("Replay U exceeds total X");
          x(r,c)=v;u(r,c)=std::min(v,w);total+=v;
        }
        if(!R_finite(total)||total<=0)stop("No live replay estimand");
        for(int r=0;r<n;++r){x(r,c)/=total;u(r,c)/=total;}
        logn[c]+=std::log(total);
      }
    }
    while(event<ne&&eday[event]==day) {
      int c=erep[event]-1;if(c<0||c>=k)stop("Invalid saved passage repeat");
      for(int r=0;r<n;++r) {
        anchor_error=std::max(anchor_error,std::abs(x(r,c)-prex(r,event)));
        anchor_error=std::max(anchor_error,std::abs(u(r,c)-preu(r,event)));
        x(r,c)=postx(r,event);u(r,c)=postu(r,event);
      }
      logn[c]=std::log(seed_cells);event++;
    }
    for(int c=0;c<k;++c)record(c,day);
    if(day%50==0)checkUserInterrupt();
  }
  if(event!=ne)stop("Saved passage anchors not consumed");
  for(int r=0;r<n*k;++r) {
    final_error=std::max(final_error,std::abs(x[r]-endx[r]));
    final_error=std::max(final_error,std::abs(u[r]-endu[r]));
  }
  return List::create(_["lineage_cell_share"]=shares,_["lineage_growth_share"]=growth,
    _["low_fraction"]=lowmass,_["untagged_low_fraction"]=umass,
    _["maximum_reference_metric_error"]=metric_error,_["maximum_pre_passage_error"]=anchor_error,
    _["maximum_balance_error"]=balance_error,_["maximum_final_state_error"]=final_error,
    _["negative_state_clipping_count"]=clipping,_["n_saved_passages"]=ne);
}
