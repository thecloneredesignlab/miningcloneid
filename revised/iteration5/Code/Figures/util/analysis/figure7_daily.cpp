#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
namespace {
void mult(const NumericMatrix& A,const NumericMatrix& x,NumericMatrix& y) {
  int n=A.nrow(),k=x.ncol();const char tr='N';const double one=1,zero=0;
  F77_CALL(dgemm)(&tr,&tr,&n,&k,&n,&one,A.begin(),&n,x.begin(),&n,&zero,y.begin(),&n FCONE FCONE);
}
double crowd(double logn,bool on,int mode,double K) {
  if(!on)return 1;
  double r=std::exp(logn-std::log(K));return mode==0?std::max(0.0,1-r):std::exp(-r);
}
}

// Passive daily daughter counters and current-cell composition. No new RNG,
// passage decision, source-specific fitness, or alteration of the total model.
// Counts are per live cell at the START of each daily interval. The saved
// log population supplies a common absolute scale; normalization never weights
// different substeps equally when their actual population sizes differ.
// [[Rcpp::export(rng = false)]]
List f7_daily_cpp(List op,List reference,NumericVector initial,LogicalVector low,
    NumericVector ploidy,NumericVector hl,NumericVector ll,NumericMatrix count_rows,
    bool culture,double dt,bool crowding,int mode,double K,double seed_cells,
    int days,bool profiles,IntegerVector snapshot_days) {
  NumericMatrix G=op["G"],KH=op["KH"],P=op["P"],PU=op["PU"],Q=op["Q"];
  NumericVector lambda=op["lambda"],mu=op["mu"],survive=op["survive"];
  NumericMatrix endx=reference["final_state"],endu=reference["final_untagged_state"];
  DataFrame events=as<DataFrame>(reference["passage_events"]);
  IntegerVector eday=events["day"],erep=events["repeat_id"];
  NumericMatrix pre=reference["pre_passage_state"],post=reference["post_passage_state"],
    preu=reference["pre_passage_untagged_state"],postu=reference["post_passage_untagged_state"];
  int n=initial.size(),k=endx.ncol(),steps=culture?std::lround(1/dt):1,ne=eday.size(),event=0;
  if(n<1||G.nrow()!=n||G.ncol()!=n||KH.nrow()!=n||KH.ncol()!=n||
    low.size()!=n||ploidy.size()!=n||hl.size()!=n||ll.size()!=n||lambda.size()!=n||mu.size()!=n||
    endx.nrow()!=n||endu.nrow()!=n||endu.ncol()!=k||steps<1||
    (culture&&std::abs(steps*dt-1)>1e-12)||(!culture&&(count_rows.nrow()!=2||count_rows.ncol()!=n)))
    stop("Invalid daily source accounting dimensions");
  NumericMatrix x(n,k),u(n,k),d(n,k),gx(n,k),gu(n,k),incoming(n,k),AU=clone(G);
  for(int i=0;i<n*n;++i)AU[i]-=KH[i];
  for(int c=0;c<k;++c)for(int i=0;i<n;++i)x(i,c)=u(i,c)=initial[i];
  NumericVector logn(k,std::log(seed_cells));
  NumericMatrix H(k,days+1),L(k,days+1),log_start(k,days+1),A(k,days+1),B(k,days+1);
  std::fill(H.begin(),H.end(),NA_REAL);std::fill(L.begin(),L.end(),NA_REAL);
  std::fill(A.begin(),A.end(),NA_REAL);std::fill(B.begin(),B.end(),NA_REAL);
  std::fill(log_start.begin(),log_start.end(),NA_REAL);
  // Five disjoint/overlapping explicitly named groups, not an implicit ancestry
  // interpretation of the complement of the first-generation compartment.
  CharacterVector group_names=CharacterVector::create("all_low","direct_first_generation","other_low",
    "historical_lineage","never_HL_lineage");
  CharacterVector metric_names=CharacterVector::create("cell_share_percent","low_birth_share_percent",
    "mean_low_birth_rate","mean_division_rate","mean_ploidy");
  List metrics(profiles?25:0);
  if(profiles)for(int i=0;i<25;++i)metrics[i]=NumericMatrix(k,days+1);
  int ns=profiles?snapshot_days.size():0;
  NumericMatrix sx(n,k*ns),su(n,k*ns),sd(n,k*ns);NumericVector scf(k*ns);
  double anchor_error=0,final_error=0,subset_error=0,complement_error=0;
  auto record=[&](int day) {
    if(!profiles)return;
    for(int c=0;c<k;++c) {
      double cf=culture?crowd(logn[c],crowding,mode,K):1;
      double mass[5]={0},births[5]={0},divs[5]={0},pl[5]={0};
      for(int i=0;i<n;++i)if(low[i]) {
        double w[5]={x(i,c),d(i,c),x(i,c)-d(i,c),x(i,c)-u(i,c),u(i,c)};
        for(int g=0;g<5;++g) {
          if(w[g]<-1e-10)stop("Negative source-profile subset");
          double v=std::max(0.0,w[g]);mass[g]+=v;births[g]+=cf*ll[i]*v;
          divs[g]+=cf*lambda[i]*v;pl[g]+=ploidy[i]*v;
        }
      }
      for(int g=0;g<5;++g) {
        double val[5]={mass[0]>0?100*mass[g]/mass[0]:NA_REAL,
          births[0]>0?100*births[g]/births[0]:NA_REAL,
          mass[g]>0?births[g]/mass[g]:NA_REAL,
          mass[g]>0?divs[g]/mass[g]:NA_REAL,mass[g]>0?pl[g]/mass[g]:NA_REAL};
        for(int m=0;m<5;++m){NumericMatrix z=metrics[g*5+m];z(c,day)=val[m];}
      }
      for(int s=0;s<ns;++s)if(snapshot_days[s]==day) {
        int col=s*k+c;scf[col]=cf;
        for(int i=0;i<n;++i){sx(i,col)=x(i,c);su(i,col)=u(i,c);sd(i,col)=d(i,c);}
      }
    }
  };
  record(0);
  for(int day=1;day<=days;++day) {
    NumericVector start=clone(logn),ch(k),cl(k);
    if(!culture)for(int c=0;c<k;++c)for(int i=0;i<n;++i) {
      ch[c]+=count_rows(0,i)*x(i,c);cl[c]+=count_rows(1,i)*x(i,c);
    }
    for(int s=0;s<steps;++s) {
      mult(culture?G:P,x,gx);mult(culture?AU:PU,u,gu);mult(culture?KH:Q,x,incoming);
      for(int c=0;c<k;++c) {
        double cf=culture?crowd(logn[c],crowding,mode,K):1,total=0;
        if(culture) {
          double scale=std::exp(logn[c]-start[c]);
          for(int i=0;i<n;++i){ch[c]+=dt*cf*hl[i]*x(i,c)*scale;cl[c]+=dt*cf*ll[i]*x(i,c)*scale;}
        }
        for(int i=0;i<n;++i) {
          double v=culture?x(i,c)+dt*(cf*gx(i,c)-mu[i]*x(i,c)):gx(i,c);
          double w=culture?u(i,c)+dt*(cf*gu(i,c)-mu[i]*u(i,c)):gu(i,c);
          double z=culture?d(i,c)+dt*(cf*incoming(i,c)-(cf*lambda[i]+mu[i])*d(i,c)):
            incoming(i,c)+survive[i]*d(i,c);
          if(!R_finite(v)||!R_finite(w)||!R_finite(z)||v< -1e-10||w< -1e-10||z< -1e-10||
            w>v+1e-9||z>v+1e-9)stop("Invalid daily replay population subset");
          subset_error=std::max(subset_error,std::max(0.0,std::max(w-v,z-v)));
          x(i,c)=std::max(0.0,v);u(i,c)=std::min(x(i,c),std::max(0.0,w));
          d(i,c)=std::min(x(i,c),std::max(0.0,z));total+=x(i,c);
        }
        if(total<=0||!R_finite(total))stop("No live daily replay population");
        for(int i=0;i<n;++i){x(i,c)/=total;u(i,c)/=total;d(i,c)/=total;}
        logn[c]+=std::log(total);
      }
    }
    for(int c=0;c<k;++c) {
      H(c,day)=ch[c];L(c,day)=cl[c];log_start(c,day)=start[c];
      double total=ch[c]+cl[c];
      if(!R_finite(total)||ch[c]<-1e-12||cl[c]<-1e-12)stop("Invalid daily production count");
      if(total>0){A(c,day)=100*ch[c]/total;B(c,day)=100*cl[c]/total;
        complement_error=std::max(complement_error,std::abs(A(c,day)+B(c,day)-100));}
    }
    // Culture births above end immediately before the integer-day passage.
    // Sampling is not daughter production; the next interval starts post-passaging.
    while(event<ne&&eday[event]==day) {
      int c=erep[event]-1;if(c<0||c>=k)stop("Invalid daily replay passage index");
      for(int i=0;i<n;++i) {
        anchor_error=std::max(anchor_error,std::abs(x(i,c)-pre(i,event)));
        anchor_error=std::max(anchor_error,std::abs(u(i,c)-preu(i,event)));
        double prob=x(i,c)>0?d(i,c)/x(i,c):0;
        d(i,c)=prob*post(i,event);x(i,c)=post(i,event);u(i,c)=postu(i,event);
      }
      logn[c]=std::log(seed_cells);event++;
    }
    record(day);if(day%25==0)checkUserInterrupt();
  }
  if(event!=ne)stop("Daily replay did not consume all passages");
  for(int i=0;i<n*k;++i)final_error=std::max(final_error,
    std::max(std::abs(x[i]-endx[i]),std::abs(u[i]-endu[i])));
  return List::create(_["high_parent_percent"]=A,_["low_parent_percent"]=B,
    _["HL_per_start_live_cell"]=H,_["LL_per_start_live_cell"]=L,_["log_start_live_cells"]=log_start,
    _["metrics"]=metrics,_["group_names"]=group_names,_["metric_names"]=metric_names,
    _["snapshot_days"]=snapshot_days,_["snapshot_X"]=sx,_["snapshot_U"]=su,_["snapshot_D"]=sd,
    _["snapshot_crowding"]=scf,_["final_direct_state"]=d,
    _["passage_error"]=anchor_error,_["final_state_error"]=final_error,
    _["subset_error"]=subset_error,_["complement_error"]=complement_error);
}
