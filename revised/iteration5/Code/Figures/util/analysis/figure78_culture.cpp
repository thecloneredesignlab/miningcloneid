#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace Rcpp;

namespace {
void mult(const NumericMatrix& A, const NumericMatrix& x, NumericMatrix& y) {
  int n=A.nrow(), k=x.ncol(); const char tr='N'; const double one=1,zero=0;
  F77_CALL(dgemm)(&tr,&tr,&n,&k,&n,&one,A.begin(),&n,x.begin(),&n,&zero,y.begin(),&n FCONE FCONE);
}
NumericVector counts_from_fraction(const NumericVector& p,double population) {
  if(!R_finite(population)||population<1||population>9007199254740991.0)
    stop("Population cannot be integerized exactly");
  double target=std::floor(population+0.5),offset=R::unif_rand(),previous=0;
  long double cumulative=0;NumericVector counts(p.size());
  for(int i=0;i<p.size();++i) {
    cumulative+=static_cast<long double>(target)*p[i];
    double upper=i+1==p.size()?target:std::min(target,static_cast<double>(std::floor(cumulative+offset)));
    counts[i]=upper-previous;previous=upper;if(counts[i]<0)stop("Negative state count");
  }
  return counts;
}
NumericVector sample_without_replacement(const NumericVector& counts,double inoculum) {
  double remaining=sum(counts),needed=inoculum;NumericVector sampled(counts.size());
  if(inoculum!=std::floor(inoculum)||inoculum<1||inoculum>remaining)stop("Invalid inoculum");
  for(int i=0;i<counts.size();++i) {
    double good=counts[i],bad=remaining-good,value=0;
    if(needed>0&&good>0)value=bad<=0?needed:(needed==remaining?good:R::rhyper(good,bad,needed));
    if(!R_finite(value)||value<0||value>good||value>needed)stop("Invalid passage draw");
    sampled[i]=value;remaining-=good;needed-=value;
  }
  if(needed!=0||sum(sampled)!=inoculum)stop("Passage count not conserved");
  return sampled;
}
double crowd_factor(double logpop,bool enabled,int mode,double K) {
  if(!enabled)return 1;
  double ratio=std::exp(logpop-std::log(K));
  return mode==0?std::max(0.0,1-ratio):std::exp(-ratio);
}
NumericMatrix as_matrix(const std::vector<double>& values,int n) {
  NumericMatrix out(n,values.size()/n);std::copy(values.begin(),values.end(),out.begin());return out;
}
}

// Same live-cell Euler update as cpp_o2simps_simulate_one. Only the externally
// imposed fixed oxygen and the approved threshold-random-passage design differ.
// G contains ordinary/WGD births and maternal removal; mu is direct death.
// K_HL contains ordinary H->L births; U'=(cG-mu-cK_HL)U. T=x-U.
// [[Rcpp::export(rng = false)]]
List f78_culture_cpp(NumericMatrix G,NumericVector mu,NumericMatrix K_HL,
    NumericVector initial,NumericVector ploidy,NumericVector hl,NumericVector ll,
    LogicalVector low,int n_day,double dt,bool crowding,int crowding_mode,double K,
    double seed_cells,double target_cells,IntegerMatrix rng_states,
    IntegerMatrix label_rng_states,bool passage=true,bool tagged=false,
    bool save_daily=false,bool save_events=false) {
  int n=G.nrow(),k=rng_states.ncol(),substeps=std::lround(1/dt);
  if(n<1||G.ncol()!=n||mu.size()!=n||initial.size()!=n||ploidy.size()!=n||
     hl.size()!=n||ll.size()!=n||low.size()!=n||n_day<0||dt<=0||substeps<1||
     std::abs(substeps*dt-1)>1e-12||rng_states.nrow()!=7||k<1||
     label_rng_states.nrow()!=7||label_rng_states.ncol()!=k||K<=0||!R_finite(K)||
     (crowding_mode!=0&&crowding_mode!=1)||seed_cells<1||seed_cells!=std::floor(seed_cells)||
     target_cells<=seed_cells||target_cells!=std::floor(target_cells))stop("Invalid culture contract");
  if(tagged&&(K_HL.nrow()!=n||K_HL.ncol()!=n))stop("Invalid lineage kernel");
  if(std::abs(sum(initial)-1)>1e-10)stop("Initial fractions do not sum to one");
  for(int i=0;i<n;++i)if(!R_finite(initial[i])||initial[i]<0||!R_finite(mu[i])||mu[i]<0)
    stop("Invalid state or direct-death rate");
  NumericMatrix x(n,k),u(n,k),gx(n,k),gu(n,k),A_U;
  if(tagged) {A_U=clone(G);for(int i=0;i<n*n;++i)A_U[i]-=K_HL[i];}
  for(int c=0;c<k;++c)for(int r=0;r<n;++r)x(r,c)=u(r,c)=initial[r];
  IntegerMatrix streams=clone(rng_states),labelstreams=clone(label_rng_states);
  NumericVector logpop(k,std::log(seed_cells)),mean(k),net(k),last_crowd(k);
  NumericVector lowfrac(k),flux(k),growth(k),tags(k),lowbirth(k),tagbirth(k);
  IntegerVector passes(k),first(k,NA_INTEGER),last(k,NA_INTEGER);
  NumericMatrix means(k,save_daily?n_day+1:1),nets(k,save_daily?n_day+1:1),
    fluxes(k,save_daily?n_day+1:1),growths(k,save_daily?n_day+1:1),
    tags_daily(k,save_daily?n_day+1:1),low_daily(k,save_daily?n_day+1:1);
  NumericVector colG(n);for(int c=0;c<n;++c)for(int r=0;r<n;++r)colG[c]+=G(r,c);
  std::vector<int> event_rep,event_day;
  std::vector<double> event_before,event_after,event_pre_mean,event_post_mean,preX,postX,preU,postU;
  double max_balance=0,min_raw=0;int clipping=0;
  auto stats=[&](int c,int day) {
    double cf=crowd_factor(logpop[c],crowding,crowding_mode,K),H=0,L=0,birth=0,tbirth=0,lf=0,m=0,rate=0,tag=0;
    for(int r=0;r<n;++r) {
      double t=tagged?x(r,c)-u(r,c):0;
      m+=ploidy[r]*x(r,c);rate+=(cf*colG[r]-mu[r])*x(r,c);
      H+=cf*hl[r]*x(r,c);L+=cf*ll[r]*x(r,c);
      // ll only includes low daughters from currently low mothers.
      birth+=cf*ll[r]*x(r,c);tbirth+=cf*ll[r]*t;
      if(low[r])lf+=x(r,c);tag+=t;
    }
    mean[c]=m;net[c]=rate;last_crowd[c]=cf;lowfrac[c]=lf;tags[c]=tag;
    flux[c]=H+L>1e-14?100*H/(H+L):NA_REAL;
    growth[c]=tagged&&birth>1e-14?100*tbirth/birth:NA_REAL;
    if(R_finite(growth[c])&&(growth[c]<-1e-8||growth[c]>100+1e-8))stop("Low-cell birth contribution outside [0,100]");
    lowbirth[c]=birth;tagbirth[c]=tbirth;
    int d=save_daily?day:0;
    means(c,d)=m;nets(c,d)=rate;fluxes(c,d)=flux[c];growths(c,d)=growth[c];tags_daily(c,d)=tag;low_daily(c,d)=lf;
  };
  for(int c=0;c<k;++c)stats(c,0);
  Environment global=Environment::global_env();bool had=global.exists(".Random.seed");
  RObject caller=R_NilValue;if(had)caller=clone(as<IntegerVector>(global[".Random.seed"]));
  for(int day=1;day<=n_day;++day) {
    for(int s=0;s<substeps;++s) {
      mult(G,x,gx);if(tagged)mult(A_U,u,gu);
      for(int c=0;c<k;++c) {
        double cf=crowd_factor(logpop[c],crowding,crowding_mode,K),total=0;
        for(int r=0;r<n;++r) {
          double value=x(r,c)+dt*(cf*gx(r,c)-mu[r]*x(r,c));
          if(!R_finite(value))stop("Nonfinite Euler live state");
          min_raw=std::min(min_raw,value);
          if(value<0){clipping++;value=0;} // Canonical core's positivity guard; counted explicitly.
          x(r,c)=value;total+=value;
          if(tagged) {
            double v=u(r,c)+dt*(cf*gu(r,c)-mu[r]*u(r,c));
            if(!R_finite(v))stop("Nonfinite Euler untagged state");
            if(v<0){clipping++;v=0;}
            max_balance=std::max(max_balance,std::max(0.0,v-value));
            if(v>value+1e-9)stop("Lineage conservation failed");
            u(r,c)=std::min(v,value);
          }
        }
        if(!R_finite(total)||total<=0)stop("Extinct live state has no ploidy estimand");
        for(int r=0;r<n;++r){x(r,c)/=total;if(tagged)u(r,c)/=total;}
        logpop[c]+=std::log(total);
      }
    }
    for(int c=0;c<k;++c) {
      if(passage&&logpop[c]>=std::log(target_cells)) {
        double population=std::exp(logpop[c]),pre_mean=0;
        for(int r=0;r<n;++r)pre_mean+=ploidy[r]*x(r,c);
        if(save_events)for(int r=0;r<n;++r){preX.push_back(x(r,c));if(tagged)preU.push_back(u(r,c));}
        global[".Random.seed"]=IntegerVector(streams(_,c));GetRNGstate();
        NumericVector p=x(_,c),counts=counts_from_fraction(p,population),sample=sample_without_replacement(counts,seed_cells);
        PutRNGstate();streams(_,c)=as<IntegerVector>(global[".Random.seed"]);
        if(tagged) {
          global[".Random.seed"]=IntegerVector(labelstreams(_,c));GetRNGstate();
          for(int r=0;r<n;++r) {
            double ratio=x(r,c)>0?std::max(0.0,std::min(1.0,u(r,c)/x(r,c))):0;
            double expected=counts[r]*ratio,good=std::floor(expected);
            if(expected>good&&R::unif_rand()<expected-good)good++;
            good=std::min(counts[r],good);double chosen=0;
            if(sample[r]>0&&good>0)chosen=good==counts[r]?sample[r]:(sample[r]==counts[r]?good:R::rhyper(good,counts[r]-good,sample[r]));
            u(r,c)=chosen/seed_cells;
          }
          PutRNGstate();labelstreams(_,c)=as<IntegerVector>(global[".Random.seed"]);
        }
        for(int r=0;r<n;++r)x(r,c)=sample[r]/seed_cells;
        passes[c]++;if(first[c]==NA_INTEGER)first[c]=day;last[c]=day;logpop[c]=std::log(seed_cells);
        if(save_events) {
          double post_mean=0;for(int r=0;r<n;++r){post_mean+=ploidy[r]*x(r,c);postX.push_back(x(r,c));if(tagged)postU.push_back(u(r,c));}
          event_rep.push_back(c+1);event_day.push_back(day);event_before.push_back(population);event_after.push_back(seed_cells);
          event_pre_mean.push_back(pre_mean);event_post_mean.push_back(post_mean);
        }
      }
      stats(c,day);
    }
    if(day%50==0)checkUserInterrupt();
  }
  if(had)global[".Random.seed"]=caller;else global.remove(".Random.seed");
  DataFrame events=DataFrame::create(_["repeat_id"]=event_rep,_["day"]=event_day,
    _["cells_before"]=event_before,_["cells_after"]=event_after,
    _["ploidy_before"]=event_pre_mean,_["ploidy_after"]=event_post_mean);
  return List::create(_["mean_ploidy"]=means,_["net_growth_rate"]=nets,
    _["flux_percent"]=fluxes,_["historical_growth_percent"]=growths,_["tag_fraction"]=tags_daily,
    _["low_fraction"]=low_daily,_["final_state"]=x,_["final_untagged_state"]=u,
    _["final_log_population"]=logpop,_["final_crowding_factor"]=last_crowd,
    _["low_birth_rate"]=lowbirth,_["tagged_low_birth_rate"]=tagbirth,
    _["passage_count"]=passes,_["first_passage_day"]=first,_["last_passage_day"]=last,
    _["passage_events"]=events,_["pre_passage_state"]=as_matrix(preX,n),_["post_passage_state"]=as_matrix(postX,n),
    _["pre_passage_untagged_state"]=as_matrix(preU,n),_["post_passage_untagged_state"]=as_matrix(postU,n),
    _["maximum_lineage_balance_error"]=max_balance,_["negative_state_clipping_count"]=clipping,
    _["minimum_raw_state"]=min_raw,_["rng_final"]=streams,_["label_rng_final"]=labelstreams);
}
