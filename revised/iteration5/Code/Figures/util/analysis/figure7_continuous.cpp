#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace {
// Same systematic integerization and without-replacement sampling as Figure6.
NumericVector integerize(const NumericVector& p, double total) {
  if (!R_finite(total)||total<1||total>9007199254740991.0)stop("Invalid passage population");
  const double target=std::floor(total+0.5),offset=R::unif_rand();
  NumericVector count(p.size());long double cumulative=0;double previous=0;
  for(int i=0;i<p.size();++i) {
    cumulative+=static_cast<long double>(target)*p[i];
    double upper=i+1==p.size()?target:std::min(target,static_cast<double>(std::floor(cumulative+offset)));
    count[i]=upper-previous;previous=upper;if(count[i]<0)stop("Negative integerized state");
  }
  return count;
}
NumericVector sample_counts(const NumericVector& count,double size) {
  double remaining=sum(count),needed=size;NumericVector sample(count.size());
  if(size!=std::floor(size)||size<1||size>remaining)stop("Invalid inoculum");
  for(int i=0;i<count.size();++i) {
    double good=count[i],bad=remaining-good,value=0;
    if(needed>0&&good>0)value=bad<=0?needed:(needed==remaining?good:R::rhyper(good,bad,needed));
    if(!R_finite(value)||value<0||value>good||value>needed)stop("Invalid hypergeometric draw");
    sample[i]=value;needed-=value;remaining-=good;
  }
  if(needed!=0||sum(sample)!=size)stop("Inoculum not conserved");return sample;
}
void multiply(const NumericMatrix& P,const NumericMatrix& x,NumericMatrix& y) {
  int n=P.nrow(),k=x.ncol();const char tr='N';const double one=1,zero=0;
  F77_CALL(dgemm)(&tr,&tr,&n,&k,&n,&one,P.begin(),&n,x.begin(),&n,&zero,y.begin(),&n FCONE FCONE);
}
struct Stats {double flux,net,tag,ploidy,hl,ll;};
Stats lineage_stats(const NumericMatrix& x,const NumericMatrix& u,int c,
    const NumericVector& pg,const NumericVector& rate,const NumericVector& hl,
    const NumericVector& ll) {
  double H=0,L=0,net=0,tag=0,ploidy=0;
  for(int r=0;r<x.nrow();++r) {
    H+=hl[r]*x(r,c);L+=ll[r]*x(r,c);net+=rate[r]*x(r,c);
    tag+=x(r,c)-u(r,c);ploidy+=pg[r]*x(r,c);
  }
  return Stats{H+L>1e-14?100*H/(H+L):NA_REAL,
    net,tag,ploidy,H,L};
}
NumericMatrix state_matrix(const std::vector<double>& x,int n) {
  int k=x.size()/n;NumericMatrix out(n,k);std::copy(x.begin(),x.end(),out.begin());return out;
}
}

// U'=(A-K)U, T'=AT+KU, x=U+T and x'=Ax. K contains ordinary H->L production.
// Separate label RNG leaves chromosome-state passage draws and their streams unchanged.
// [[Rcpp::export(rng = false)]]
List f7_continuous_cpp(NumericMatrix step,NumericMatrix untagged_step,
    NumericMatrix initial,NumericVector ploidy_grid,NumericVector net_rate,
    NumericVector hl,NumericVector ll,int n_day,double seed_cells,double target_cells,
    IntegerMatrix rng_states,IntegerMatrix label_rng_states,bool passage,
    bool save_passage_states=true) {
  int n=step.nrow(),k=initial.ncol();
  if(n<1||step.ncol()!=n||untagged_step.nrow()!=n||untagged_step.ncol()!=n||
      initial.nrow()!=n||ploidy_grid.size()!=n||net_rate.size()!=n||hl.size()!=n||ll.size()!=n||
      rng_states.nrow()!=7||rng_states.ncol()!=k||label_rng_states.nrow()!=7||label_rng_states.ncol()!=k||
      n_day<0||seed_cells<1||seed_cells!=std::floor(seed_cells)||target_cells<=seed_cells||
      target_cells!=std::floor(target_cells))stop("Invalid tagged propagation contract");
  NumericMatrix x=clone(initial),u=clone(initial),xn(n,k),un(n,k);
  NumericMatrix flux(k,n_day+1),net(k,n_day+1),tag(k,n_day+1),ploidy(k,n_day+1);
  IntegerMatrix streams=clone(rng_states),labelStreams=clone(label_rng_states);
  NumericVector logpop(k,std::log(seed_cells));IntegerVector passage_count(k);
  std::vector<int> evRep,evDay;
  std::vector<double> evBefore,evAfter,evPreP,evPostP,evPreTag,evPostTag,evPreFlux,evPostFlux;
  std::vector<double> preX,postX,preU,postU;
  double maximumError=0;
  auto record=[&](int c,int day) {
    Stats s=lineage_stats(x,u,c,ploidy_grid,net_rate,hl,ll);
    flux(c,day)=s.flux;net(c,day)=s.net;tag(c,day)=s.tag;ploidy(c,day)=s.ploidy;
  };
  for(int c=0;c<k;++c) {
    double total=0;for(int r=0;r<n;++r){if(x(r,c)<0||!R_finite(x(r,c)))stop("Invalid initial composition");total+=x(r,c);}
    if(std::abs(total-1)>1e-10)stop("Initial composition not normalized");record(c,0);
  }
  Environment global=Environment::global_env();bool had=global.exists(".Random.seed");
  RObject caller=R_NilValue;if(had)caller=clone(as<IntegerVector>(global[".Random.seed"]));
  for(int day=1;day<=n_day;++day) {
    multiply(step,x,xn);multiply(untagged_step,u,un);
    for(int c=0;c<k;++c) {
      double total=0;
      for(int r=0;r<n;++r) {
        if(xn(r,c)<0&&xn(r,c)>-1e-10)xn(r,c)=0;
        if(un(r,c)<0&&un(r,c)>-1e-10)un(r,c)=0;
        if(!R_finite(xn(r,c))||xn(r,c)<0||!R_finite(un(r,c))||un(r,c)<0)stop("Invalid propagated lineage state");
        total+=xn(r,c);
      }
      if(!R_finite(total)||total<=0)stop("Non-positive growth factor");
      for(int r=0;r<n;++r) {
        xn(r,c)/=total;un(r,c)/=total;
        double error=std::max(0.0,un(r,c)-xn(r,c));maximumError=std::max(maximumError,error);
        if(error>1e-9)stop("U+T lineage balance failed");
        if(un(r,c)>xn(r,c))un(r,c)=xn(r,c);
      }
      logpop[c]+=std::log(total);
      if(passage&&logpop[c]>=std::log(target_cells)) {
        double population=std::exp(logpop[c]);
        Stats before=lineage_stats(xn,un,c,ploidy_grid,net_rate,hl,ll);
        if(save_passage_states)for(int r=0;r<n;++r){preX.push_back(xn(r,c));preU.push_back(un(r,c));}
        global[".Random.seed"]=IntegerVector(streams(_,c));GetRNGstate();
        NumericVector p=xn(_,c),counts=integerize(p,population),sample=sample_counts(counts,seed_cells);
        PutRNGstate();streams(_,c)=as<IntegerVector>(global[".Random.seed"]);
        // Conditional source-label sampling inside each chromosome state.
        global[".Random.seed"]=IntegerVector(labelStreams(_,c));GetRNGstate();
        for(int r=0;r<n;++r) {
          double ratio=xn(r,c)>0?std::max(0.0,std::min(1.0,un(r,c)/xn(r,c))):0;
          double expected=counts[r]*ratio,good=std::floor(expected);
          if(expected>good&&R::unif_rand()<expected-good)good+=1;
          good=std::min(counts[r],good);double chosen=0;
          if(sample[r]>0&&good>0)chosen=good==counts[r]?sample[r]:
            (sample[r]==counts[r]?good:R::rhyper(good,counts[r]-good,sample[r]));
          xn(r,c)=sample[r]/seed_cells;un(r,c)=chosen/seed_cells;
        }
        PutRNGstate();labelStreams(_,c)=as<IntegerVector>(global[".Random.seed"]);
        Stats after=lineage_stats(xn,un,c,ploidy_grid,net_rate,hl,ll);
        evRep.push_back(c+1);evDay.push_back(day);evBefore.push_back(population);evAfter.push_back(seed_cells);
        evPreP.push_back(before.ploidy);evPostP.push_back(after.ploidy);
        evPreTag.push_back(before.tag);evPostTag.push_back(after.tag);
        evPreFlux.push_back(before.flux);evPostFlux.push_back(after.flux);
        if(save_passage_states)for(int r=0;r<n;++r){postX.push_back(xn(r,c));postU.push_back(un(r,c));}
        passage_count[c]++;logpop[c]=std::log(seed_cells);
      }
    }
    std::copy(xn.begin(),xn.end(),x.begin());std::copy(un.begin(),un.end(),u.begin());
    for(int c=0;c<k;++c)record(c,day);
    if(day%50==0)checkUserInterrupt();
  }
  if(had)global[".Random.seed"]=caller;else global.remove(".Random.seed");
  DataFrame events=DataFrame::create(_["repeat_id"]=evRep,_["day"]=evDay,
    _["cells_before"]=evBefore,_["cells_after"]=evAfter,_["ploidy_before"]=evPreP,
    _["ploidy_after"]=evPostP,_["tag_fraction_before"]=evPreTag,_["tag_fraction_after"]=evPostTag,
    _["flux_percent_before"]=evPreFlux,_["flux_percent_after"]=evPostFlux);
  return List::create(_["flux_percent"]=flux,_["net_growth"]=net,
    _["tag_fraction"]=tag,_["mean_ploidy"]=ploidy,_["final_state"]=x,_["final_untagged_state"]=u,
    _["final_log_population"]=logpop,_["rng_final"]=streams,_["label_rng_final"]=labelStreams,
    _["passage_count"]=passage_count,_["passage_events"]=events,
    _["pre_passage_state"]=state_matrix(preX,n),_["post_passage_state"]=state_matrix(postX,n),
    _["pre_passage_untagged_state"]=state_matrix(preU,n),_["post_passage_untagged_state"]=state_matrix(postU,n),
    _["maximum_lineage_balance_error"]=maximumError);
}
