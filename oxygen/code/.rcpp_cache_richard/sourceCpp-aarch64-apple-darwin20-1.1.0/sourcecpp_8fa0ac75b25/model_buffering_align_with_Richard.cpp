#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

namespace {

inline int boundary_mode(const std::string& boundary) {
  if (boundary == "drop") return 0;
  if (boundary == "absorb_minmax") return 1;
  stop("boundary must be one of: drop, absorb_minmax");
}

inline void append_with_boundary(
    int Np,
    int row_min,
    int row_max,
    int col_1based,
    double value,
    int bmode,
    std::vector<int>& ii,
    std::vector<int>& jj,
    std::vector<double>& xx
) {
  if (Np < row_min || Np > row_max) {
    if (bmode == 1) {
      const int Np2 = std::max(std::min(Np, row_max), row_min);
      ii.push_back(Np2 - row_min + 1);
      jj.push_back(col_1based);
      xx.push_back(value);
    }
    return;
  }
  ii.push_back(Np - row_min + 1);
  jj.push_back(col_1based);
  xx.push_back(value);
}

void richard_pr_delta_internal(
    int N,
    double p,
    double eps_tail,
    double beta_buffer,
    double n_exp,
    double smax,
    int N_unit,
    std::vector<int>& ts_out,
    std::vector<double>& prob_out,
    double& mass_dropped
) {
  (void) eps_tail; // kept for API compatibility with the R implementation

  if (p <= 0.0 || N <= 0) {
    ts_out.assign(1, 0);
    prob_out.assign(1, 1.0);
    mass_dropped = 0.0;
    return;
  }

  const double sd = std::sqrt(static_cast<double>(N) * p);
  if (sd == 0.0) {
    ts_out.assign(1, 0);
    prob_out.assign(1, 1.0);
    mass_dropped = 0.0;
    return;
  }

  const double n_d = static_cast<double>(N);
  const double n_unit_d = static_cast<double>(N_unit);
  const double sN = smax * std::exp(-beta_buffer * std::pow((2.0 * n_unit_d) / n_d, n_exp));

  const double z = 9.0;
  const int T = std::min(N, std::max(0, static_cast<int>(std::ceil(z * sd))));
  const int len = 2 * T + 1;

  ts_out.resize(len);
  prob_out.assign(len, 0.0);

  for (int idx = 0; idx < len; ++idx) {
    const int t = idx - T;
    ts_out[idx] = t;
    const int k_start = std::abs(t);
    double acc = 0.0;

    for (int ks = k_start; ks <= N; ks += 2) {
      const double pk = R::dbinom(ks, N, p, false);
      const double m = (static_cast<double>(ks) + static_cast<double>(t)) / 2.0;
      const double qm = R::dbinom(m, ks, 0.5, false);
      const double s_pow = std::pow(sN, static_cast<double>(ks));
      acc += pk * qm * s_pow;
    }
    prob_out[idx] = acc;
  }

  const double total = std::accumulate(prob_out.begin(), prob_out.end(), 0.0);
  mass_dropped = std::max(0.0, 1.0 - total);
}

} // namespace

// [[Rcpp::export]]
List cpp_richard_pr_delta_vec(
    int N,
    double p,
    double eps_tail = 1e-8,
    double beta_buffer = 0.0,
    double n_exp = 1.0,
    double smax = 1.0,
    int N_unit = 22
) {
  std::vector<int> ts;
  std::vector<double> prob;
  double mass_dropped = 0.0;

  richard_pr_delta_internal(
    N,
    p,
    eps_tail,
    beta_buffer,
    n_exp,
    smax,
    N_unit,
    ts,
    prob,
    mass_dropped
  );

  return List::create(
    _["ts"] = IntegerVector(ts.begin(), ts.end()),
    _["prob"] = NumericVector(prob.begin(), prob.end()),
    _["mass_dropped"] = mass_dropped
  );
}

// [[Rcpp::export]]
List cpp_richard_build_B_total_triplet(
    int Nmin,
    int Nmax,
    NumericVector p_vec,
    std::string boundary = "drop",
    double eps_tail = 1e-8,
    double beta_buffer = 0.0,
    double n_exp = 1.0,
    double smax = 1.0,
    int N_unit = 22
) {
  const int R = Nmax - Nmin + 1;
  if (R <= 0) stop("Nmax must be >= Nmin");

  const int p_len = p_vec.size();
  if (!(p_len == 1 || p_len == R)) stop("p_vec length must be 1 or R");

  const int bmode = boundary_mode(boundary);

  std::vector<int> ii;
  std::vector<int> jj;
  std::vector<double> xx;
  ii.reserve(static_cast<size_t>(R) * 12);
  jj.reserve(static_cast<size_t>(R) * 12);
  xx.reserve(static_cast<size_t>(R) * 12);

  for (int col = 0; col < R; ++col) {
    const int N = Nmin + col;
    double pN = (p_len == 1) ? p_vec[0] : p_vec[col];
    pN = std::max(0.0, std::min(1.0, pN));

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    richard_pr_delta_internal(
      N,
      pN,
      eps_tail,
      beta_buffer,
      n_exp,
      smax,
      N_unit,
      ts,
      pr,
      mass_dropped
    );

    const int col_1based = col + 1;
    const int K = static_cast<int>(ts.size());
    for (int k = 0; k < K; ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;

      if (t == 0) {
        append_with_boundary(
          N,
          Nmin,
          Nmax,
          col_1based,
          2.0 * w,
          bmode,
          ii,
          jj,
          xx
        );
      } else {
        append_with_boundary(
          N + t,
          Nmin,
          Nmax,
          col_1based,
          w,
          bmode,
          ii,
          jj,
          xx
        );
        append_with_boundary(
          N - t,
          Nmin,
          Nmax,
          col_1based,
          w,
          bmode,
          ii,
          jj,
          xx
        );
      }
    }
  }

  return List::create(
    _["i"] = IntegerVector(ii.begin(), ii.end()),
    _["j"] = IntegerVector(jj.begin(), jj.end()),
    _["x"] = NumericVector(xx.begin(), xx.end()),
    _["nrow"] = R,
    _["ncol"] = R
  );
}

// [[Rcpp::export]]
List cpp_richard_build_B_WGD_triplet(
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    std::string boundary = "drop",
    double wgd_value = 1.0
) {
  const int R0 = N0max - N0min + 1;
  const int R1 = N1max - N1min + 1;
  if (R0 <= 0 || R1 <= 0) stop("Nmax must be >= Nmin for both layers");

  const int bmode = boundary_mode(boundary);
  const double val = wgd_value;

  std::vector<int> ii;
  std::vector<int> jj;
  std::vector<double> xx;
  ii.reserve(static_cast<size_t>(R0));
  jj.reserve(static_cast<size_t>(R0));
  xx.reserve(static_cast<size_t>(R0));

  for (int col = 0; col < R0; ++col) {
    const int N0 = N0min + col;
    const int Np = 2 * N0;
    const int col_1based = col + 1;
    append_with_boundary(
      Np,
      N1min,
      N1max,
      col_1based,
      val,
      bmode,
      ii,
      jj,
      xx
    );
  }

  return List::create(
    _["i"] = IntegerVector(ii.begin(), ii.end()),
    _["j"] = IntegerVector(jj.begin(), jj.end()),
    _["x"] = NumericVector(xx.begin(), xx.end()),
    _["nrow"] = R1,
    _["ncol"] = R0
  );
}


#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_richard_pr_delta_vec
List cpp_richard_pr_delta_vec(int N, double p, double eps_tail, double beta_buffer, double n_exp, double smax, int N_unit);
RcppExport SEXP sourceCpp_3_cpp_richard_pr_delta_vec(SEXP NSEXP, SEXP pSEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type eps_tail(eps_tailSEXP);
    Rcpp::traits::input_parameter< double >::type beta_buffer(beta_bufferSEXP);
    Rcpp::traits::input_parameter< double >::type n_exp(n_expSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< int >::type N_unit(N_unitSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_richard_pr_delta_vec(N, p, eps_tail, beta_buffer, n_exp, smax, N_unit));
    return rcpp_result_gen;
END_RCPP
}
// cpp_richard_build_B_total_triplet
List cpp_richard_build_B_total_triplet(int Nmin, int Nmax, NumericVector p_vec, std::string boundary, double eps_tail, double beta_buffer, double n_exp, double smax, int N_unit);
RcppExport SEXP sourceCpp_3_cpp_richard_build_B_total_triplet(SEXP NminSEXP, SEXP NmaxSEXP, SEXP p_vecSEXP, SEXP boundarySEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Nmin(NminSEXP);
    Rcpp::traits::input_parameter< int >::type Nmax(NmaxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_vec(p_vecSEXP);
    Rcpp::traits::input_parameter< std::string >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< double >::type eps_tail(eps_tailSEXP);
    Rcpp::traits::input_parameter< double >::type beta_buffer(beta_bufferSEXP);
    Rcpp::traits::input_parameter< double >::type n_exp(n_expSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< int >::type N_unit(N_unitSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_richard_build_B_total_triplet(Nmin, Nmax, p_vec, boundary, eps_tail, beta_buffer, n_exp, smax, N_unit));
    return rcpp_result_gen;
END_RCPP
}
// cpp_richard_build_B_WGD_triplet
List cpp_richard_build_B_WGD_triplet(int N0min, int N0max, int N1min, int N1max, std::string boundary, double wgd_value);
RcppExport SEXP sourceCpp_3_cpp_richard_build_B_WGD_triplet(SEXP N0minSEXP, SEXP N0maxSEXP, SEXP N1minSEXP, SEXP N1maxSEXP, SEXP boundarySEXP, SEXP wgd_valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N0min(N0minSEXP);
    Rcpp::traits::input_parameter< int >::type N0max(N0maxSEXP);
    Rcpp::traits::input_parameter< int >::type N1min(N1minSEXP);
    Rcpp::traits::input_parameter< int >::type N1max(N1maxSEXP);
    Rcpp::traits::input_parameter< std::string >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< double >::type wgd_value(wgd_valueSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_richard_build_B_WGD_triplet(N0min, N0max, N1min, N1max, boundary, wgd_value));
    return rcpp_result_gen;
END_RCPP
}
