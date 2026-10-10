// Positivity preserving Noda inverse iteration on the original double matrix.
// No rates, graph edges, parameter ranges, or matrix entries are regularized.
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <vector>
#include <algorithm>

template <unsigned Digits>
Rcpp::List perron_impl(const Rcpp::NumericMatrix& input, int max_iterations) {
  using Real = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<Digits>>;
  const int n = input.nrow();
  if (n < 2 || input.ncol() != n) Rcpp::stop("Require a square Metzler matrix");
  std::vector<Real> m(n*n), x(n, Real(1)/n), mx(n), y(n), lu(n*n), rhs(n);
  Real scale = 1;
  for (int i=0; i<n; ++i) for (int j=0; j<n; ++j) {
    double a=input(i,j);
    if (!R_finite(a) || (i!=j && a<0)) Rcpp::stop("Not a finite Metzler matrix");
    m[i*n+j]=Real(a);
    if (abs(m[i*n+j])>scale) scale=abs(m[i*n+j]);
  }
  const Real shift_floor = scale*Real("1e-25");
  const Real residual_tolerance = scale*Real("1e-24");
  Real upper=0, lambda=0, residual=0, change=0;
  auto measure = [&]() {
    lambda=0; residual=0;
    for (int i=0; i<n; ++i) {
      mx[i]=0;
      for (int j=0; j<n; ++j) mx[i]+=m[i*n+j]*x[j];
      if (!(x[i]>0)) Rcpp::stop("Lost positivity in multiprecision iterate");
      Real ratio=mx[i]/x[i];
      if (i==0 || ratio>upper) upper=ratio;
      lambda+=mx[i];
    }
    for (int i=0; i<n; ++i) {
      Real r=abs(mx[i]-lambda*x[i]);
      if (r>residual) residual=r;
    }
  };
  measure();
  for (int iteration=1; iteration<=max_iterations; ++iteration) {
    Real shift=upper+shift_floor;
    for (int i=0; i<n; ++i) {
      rhs[i]=x[i];
      for (int j=0; j<n; ++j) lu[i*n+j]=-m[i*n+j];
      lu[i*n+i]+=shift;
    }
    // Dense M-matrix LU without row exchanges in decimal arithmetic; do not perturb pivots.
    for (int k=0; k<n; ++k) {
      if (!(lu[k*n+k]>0)) Rcpp::stop("Nonpositive M-matrix pivot at the safeguarded shift");
      for (int i=k+1; i<n; ++i) {
        if (lu[i*n+k]==0) continue;
        Real factor=lu[i*n+k]/lu[k*n+k];
        lu[i*n+k]=0;
        for (int j=k+1; j<n; ++j) lu[i*n+j]-=factor*lu[k*n+j];
        rhs[i]-=factor*rhs[k];
      }
    }
    Real sum=0;
    for (int i=n-1; i>=0; --i) {
      Real value=rhs[i];
      for (int j=i+1; j<n; ++j) value-=lu[i*n+j]*y[j];
      y[i]=value/lu[i*n+i];
      if (!(y[i]>0)) Rcpp::stop("Shifted inverse did not preserve positivity");
      sum+=y[i];
    }
    change=0;
    for (int i=0; i<n; ++i) {
      y[i]/=sum; change+=abs(y[i]-x[i]); x[i]=y[i];
    }
    measure();
    if (residual<residual_tolerance && abs(upper-lambda)<residual_tolerance && change<Real("1e-18")) {
      Rcpp::NumericVector vector(n);
      for (int i=0; i<n; ++i) vector[i]=x[i].template convert_to<double>();
      return Rcpp::List::create(Rcpp::Named("vector")=vector,
        Rcpp::Named("lambda")=lambda.template convert_to<double>(),
        Rcpp::Named("upper_bound")=upper.template convert_to<double>(),
        Rcpp::Named("upper_minus_lambda")=Real(abs(upper-lambda)).template convert_to<double>(),
        Rcpp::Named("residual")=residual.template convert_to<double>(),
        Rcpp::Named("vector_change_l1")=change.template convert_to<double>(),
        Rcpp::Named("iterations")=iteration, Rcpp::Named("decimal_digits")=Digits);
    }
    if (iteration % 10 == 0) Rcpp::checkUserInterrupt();
  }
  Rcpp::stop("Noda iteration did not meet the residual, leading-root bound and vector stability criteria");
  return Rcpp::List();
}

// [[Rcpp::export]]
Rcpp::List perron_high_precision(Rcpp::NumericMatrix input, int digits=50, int max_iterations=200) {
  if (digits==50) return perron_impl<50>(input,max_iterations);
  if (digits==100) return perron_impl<100>(input,max_iterations);
  Rcpp::stop("Supported decimal precision: 50 or 100");
  return Rcpp::List();
}
