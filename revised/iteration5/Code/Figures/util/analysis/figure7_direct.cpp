#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
namespace {
void multiply(const NumericMatrix& A,const NumericMatrix& x,NumericMatrix& y) {
  int n=A.nrow(),k=x.ncol();const char tr='N';const double one=1,zero=0;
  F77_CALL(dgemm)(&tr,&tr,&n,&k,&n,&one,A.begin(),&n,x.begin(),&n,&zero,y.begin(),&n FCONE FCONE);
}
double crowd_factor(double logn,bool on,int mode,double K) {
  if(!on)return 1;
  double r=std::exp(logn-std::log(K));return mode==0?std::max(0.0,1-r):std::exp(-r);
}
double discrepancy(double a,double b) {
  if(!R_finite(a)||!R_finite(b))return R_finite(a)==R_finite(b)?0:R_PosInf;
  return std::abs(a-b)/(1+std::abs(b));
}
}

// A passive, NON-INHERITED birth-origin label. D has no J_LL * D birth term.
// The X dynamics, density feedback and saved passage samples remain unchanged.
// [[Rcpp::export(rng = false)]]
List f7_direct_cpp(List op,List reference,NumericVector initial,LogicalVector low,
    NumericVector ploidy,NumericVector hl,NumericVector ll,bool culture,double dt,
    bool crowding,int mode,double K,double seed_cells,int days) {
  NumericMatrix G=op["G"],KH=op["KH"],P=op["P"],Q=op["Q"];
  NumericVector lambda=op["lambda"],mu=op["mu"],survive=op["survive"];
  NumericMatrix rf=reference["flux_percent"],rm=reference["mean_ploidy"],endx=reference["final_state"];
  DataFrame events=as<DataFrame>(reference["passage_events"]);
  IntegerVector eday=events["day"],erep=events["repeat_id"];
  NumericMatrix pre=reference["pre_passage_state"],post=reference["post_passage_state"];
  int n=initial.size(),k=rf.nrow(),steps=culture?std::lround(1/dt):1,ne=eday.size(),event=0;
  if(n<1||G.nrow()!=n||G.ncol()!=n||KH.nrow()!=n||KH.ncol()!=n||lambda.size()!=n||mu.size()!=n||
     low.size()!=n||ploidy.size()!=n||hl.size()!=n||ll.size()!=n||rf.ncol()!=days+1||
     rm.nrow()!=k||rm.ncol()!=days+1||endx.nrow()!=n||endx.ncol()!=k||pre.ncol()!=ne||post.ncol()!=ne||
     steps<1||(culture&&std::abs(steps*dt-1)>1e-12))stop("Invalid direct-source replay dimensions");
  if(!culture&&(ne!=0||P.nrow()!=n||P.ncol()!=n||Q.nrow()!=n||Q.ncol()!=n||survive.size()!=n))
    stop("Invalid exact continuous direct-source propagator");
  NumericMatrix x(n,k),d(n,k),gx(n,k),inflow(n,k);
  for(int c=0;c<k;++c)for(int i=0;i<n;++i)x(i,c)=initial[i];
  NumericVector logn(k,std::log(seed_cells));
  NumericMatrix growth(k,days+1),cell_share(k,days+1),all_births(k,days+1),direct_births(k,days+1);
  double flux_error=0,ploidy_error=0,anchor_error=0,final_error=0,subset_error=0,balance_error=0;
  int corrections=0;
  auto record=[&](int c,int day) {
    double cf=culture?crowd_factor(logn[c],crowding,mode,K):1,H=0,L=0,D=0,m=0,dm=0,mp=0;
    for(int i=0;i<n;++i) {
      H+=cf*hl[i]*x(i,c);L+=cf*ll[i]*x(i,c);D+=cf*ll[i]*d(i,c);mp+=ploidy[i]*x(i,c);
      if(low[i]){m+=x(i,c);dm+=d(i,c);}
      subset_error=std::max(subset_error,std::max(0.0,d(i,c)-x(i,c)));
      if(!low[i]&&std::abs(d(i,c))>1e-12)stop("Direct-origin cells cannot occupy H");
    }
    growth(c,day)=culture?(L>1e-14?100*D/L:NA_REAL):
      (m>1e-14&&L/m>1e-8?100*D/L:NA_REAL);
    cell_share(c,day)=m>1e-14?100*dm/m:NA_REAL;
    all_births(c,day)=L;direct_births(c,day)=D;
    balance_error=std::max(balance_error,std::abs(D+(L-D)-L));
    flux_error=std::max(flux_error,discrepancy(H+L>1e-14?100*H/(H+L):NA_REAL,rf(c,day)));
    ploidy_error=std::max(ploidy_error,discrepancy(mp,rm(c,day)));
    if(R_finite(growth(c,day))&&(growth(c,day)<-1e-8||growth(c,day)>100+1e-8))
      stop("Direct-source growth percentage outside its natural range");
  };
  for(int c=0;c<k;++c)record(c,0);
  for(int day=1;day<=days;++day) {
    for(int s=0;s<steps;++s) {
      multiply(culture?G:P,x,gx);multiply(culture?KH:Q,x,inflow);
      for(int c=0;c<k;++c) {
        double cf=culture?crowd_factor(logn[c],crowding,mode,K):1,total=0;
        for(int i=0;i<n;++i) {
          double v=culture?x(i,c)+dt*(cf*gx(i,c)-mu[i]*x(i,c)):gx(i,c);
          double w=culture?d(i,c)+dt*(cf*inflow(i,c)-(cf*lambda[i]+mu[i])*d(i,c)):
            inflow(i,c)+survive[i]*d(i,c);
          if(!R_finite(v)||!R_finite(w))stop("Nonfinite direct-source state");
          if(v<-1e-10||w<-1e-10||w>v+1e-9)stop("Invalid direct-source population subset");
          if(v<0){v=0;corrections++;}if(w<0){w=0;corrections++;}
          if(w>v){w=v;corrections++;}
          subset_error=std::max(subset_error,std::max(0.0,w-v));
          x(i,c)=v;d(i,c)=w;total+=v;
        }
        if(total<=0||!R_finite(total))stop("No live replay population");
        for(int i=0;i<n;++i){x(i,c)/=total;d(i,c)/=total;}
        logn[c]+=std::log(total);
      }
    }
    while(event<ne&&eday[event]==day) {
      int c=erep[event]-1;if(c<0||c>=k)stop("Invalid saved passage repeat");
      for(int i=0;i<n;++i) {
        anchor_error=std::max(anchor_error,std::abs(x(i,c)-pre(i,event)));
        // Exact conditional expectation given the saved within-state sample.
        double direct_probability=x(i,c)>0?d(i,c)/x(i,c):0;
        if(direct_probability<0||direct_probability>1+1e-12)stop("Invalid direct passage probability");
        d(i,c)=direct_probability*post(i,event);x(i,c)=post(i,event);
      }
      logn[c]=std::log(seed_cells);event++;
    }
    for(int c=0;c<k;++c)record(c,day);
    if(day%25==0)checkUserInterrupt();
  }
  if(event!=ne)stop("Unconsumed saved passages");
  for(int i=0;i<n*k;++i)final_error=std::max(final_error,std::abs(x[i]-endx[i]));
  return List::create(_["growth_percent"]=growth,_["direct_cell_share"]=cell_share,
    _["all_low_births"]=all_births,_["direct_low_births"]=direct_births,
    _["final_direct_state"]=d,_["reference_flux_error"]=flux_error,
    _["reference_ploidy_error"]=ploidy_error,_["saved_passage_error"]=anchor_error,
    _["final_state_error"]=final_error,_["subset_error"]=subset_error,
    _["growth_balance_error"]=balance_error,_["roundoff_corrections"]=corrections,
    _["n_saved_passages"]=ne);
}

// Minimal unit-test hook for the same non-inherited cohort Euler equation.
// [[Rcpp::export(rng = false)]]
NumericVector f7_direct_cohort_step_cpp(NumericVector d,NumericVector incoming,
    NumericVector lambda,NumericVector mu,double cf,double dt) {
  int n=d.size();NumericVector answer(n);
  for(int i=0;i<n;++i)answer[i]=d[i]+dt*(cf*incoming[i]-(cf*lambda[i]+mu[i])*d[i]);
  return answer;
}
