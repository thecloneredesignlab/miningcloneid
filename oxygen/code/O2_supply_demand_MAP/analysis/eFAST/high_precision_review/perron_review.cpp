// Same fixed double matrix; multiprecision linear algebra, decimal output.
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <vector>
#include <sstream>
#include <iomanip>

template<unsigned D> Rcpp::List solve(const Rcpp::NumericMatrix& input,
 const Rcpp::NumericVector& states, double unit, int limit) {
 using R=boost::multiprecision::number<boost::multiprecision::cpp_dec_float<D>>;
 const int n=input.nrow();
 if(n<2||input.ncol()!=n||states.size()!=n||unit<=0) Rcpp::stop("Invalid dimensions/unit");
 std::vector<R> m(n*n),x(n,R(1)/n),mx(n),y(n),lu(n*n),rhs(n);
 R norm=0;
 for(int i=0;i<n;++i){R sum=0;
  for(int j=0;j<n;++j){double a=input(i,j);
   if(!R_finite(a)||(i!=j&&a<0))Rcpp::stop("Not a finite Metzler matrix");
   m[i*n+j]=R(a);sum+=abs(m[i*n+j]);}
  if(sum>norm)norm=sum;}
 if(norm==0)Rcpp::stop("Zero matrix: nonunique eigenvector");
 const R shift_floor=norm*pow(R(10),-int(D*8/10));
 const R tolerance=pow(R(10),-int(D*65/100));
 const R vector_tolerance=pow(R(10),-int(D*6/10));
 const R component_tolerance=pow(R(10),-int(D*4/10));
 R upper=0,lower=0,lambda=0,residual=0,relative=0,change=0,component_change=0,last_sum=1;
 auto measure=[&](){lambda=0;residual=0;R xmax=0;
  for(int i=0;i<n;++i){mx[i]=0;for(int j=0;j<n;++j)mx[i]+=m[i*n+j]*x[j];
   if(!(x[i]>0))Rcpp::stop("Lost strict positivity");
   R ratio=mx[i]/x[i];if(i==0||ratio>upper)upper=ratio;if(i==0||ratio<lower)lower=ratio;
   lambda+=mx[i];if(x[i]>xmax)xmax=x[i];}
  for(int i=0;i<n;++i){R v=abs(mx[i]-lambda*x[i]);if(v>residual)residual=v;}
  relative=residual/((norm+abs(lambda))*xmax);};
 measure();
 int it=0;bool converged=false;
 for(it=1;it<=limit;++it){R shift=upper+shift_floor;
  for(int i=0;i<n;++i){rhs[i]=x[i];for(int j=0;j<n;++j)lu[i*n+j]=-m[i*n+j];lu[i*n+i]+=shift;}
  for(int k=0;k<n;++k){if(!(lu[k*n+k]>0))Rcpp::stop("Nonpositive safeguarded M-matrix pivot");
   for(int i=k+1;i<n;++i){if(lu[i*n+k]==0)continue;R factor=lu[i*n+k]/lu[k*n+k];lu[i*n+k]=0;
    for(int j=k+1;j<n;++j)lu[i*n+j]-=factor*lu[k*n+j];rhs[i]-=factor*rhs[k];}}
  last_sum=0;
  for(int i=n-1;i>=0;--i){R value=rhs[i];for(int j=i+1;j<n;++j)value-=lu[i*n+j]*y[j];
   y[i]=value/lu[i*n+i];if(!(y[i]>0))Rcpp::stop("Inverse failed to preserve positivity");last_sum+=y[i];}
  change=0;component_change=0;for(int i=0;i<n;++i){y[i]/=last_sum;R delta=abs(y[i]-x[i]);change+=delta;
   R relative_component=delta/(y[i]>x[i]?y[i]:x[i]);if(relative_component>component_change)component_change=relative_component;x[i]=y[i];}
  measure();
  if(relative<tolerance&&abs(upper-lambda)/norm<tolerance&&change<vector_tolerance&&component_change<component_tolerance){converged=true;break;}
  if(it%5==0)Rcpp::checkUserInterrupt();}
 auto str=[](const R& z){return z.str(D,std::ios_base::scientific);};
 R sum=0,ploidy=0;Rcpp::CharacterVector vector(n);
 for(int i=0;i<n;++i){sum+=x[i];ploidy+=R(states[i])*x[i]/R(unit);vector[i]=str(x[i]);}
 return Rcpp::List::create(Rcpp::Named("decimal_digits")=D,Rcpp::Named("converged")=converged,
  Rcpp::Named("iterations")=std::min(it,limit),Rcpp::Named("vector_decimal")=vector,
  Rcpp::Named("lambda")=str(lambda),Rcpp::Named("ploidy")=str(ploidy),Rcpp::Named("sum")=str(sum),
  Rcpp::Named("residual_inf")=str(residual),Rcpp::Named("relative_residual_inf")=str(relative),
  Rcpp::Named("matrix_norm_inf")=str(norm),Rcpp::Named("lower_bound")=str(lower),Rcpp::Named("upper_bound")=str(upper),
  Rcpp::Named("relative_upper_minus_lambda")=str(abs(upper-lambda)/norm),
  Rcpp::Named("vector_change_l1")=str(change),Rcpp::Named("max_relative_component_change")=str(component_change),
  Rcpp::Named("component_tolerance")=str(component_tolerance),Rcpp::Named("inverse_sum_before_normalization")=str(last_sum),
  Rcpp::Named("relative_tolerance")=str(tolerance),Rcpp::Named("vector_tolerance")=str(vector_tolerance));
}

// [[Rcpp::export]]
Rcpp::List review_perron(Rcpp::NumericMatrix input,Rcpp::NumericVector states,
 double unit=22,int digits=50,int max_iterations=250){
 if(digits==50)return solve<50>(input,states,unit,max_iterations);
 if(digits==100)return solve<100>(input,states,unit,max_iterations);
 Rcpp::stop("Only 50 and 100 decimal digits are supported");return Rcpp::List();
}
