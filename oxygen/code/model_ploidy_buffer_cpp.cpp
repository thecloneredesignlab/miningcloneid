#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

namespace {

void pr_delta_internal(
    int N,
    double p,
    double eps_tail,
    double mr_lethality,
    std::vector<int>& ts_out,
    std::vector<double>& prob_out,
    double& mass_dropped
) {
  if (p <= 0.0 || N == 0) {
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

  const double z = R::qnorm5(1.0 - eps_tail / 2.0, 0.0, 1.0, 1, 0);
  const int T = std::min(N, std::max(0, static_cast<int>(std::ceil(z * sd))));
  const int len = 2 * T + 1;

  ts_out.resize(len);
  prob_out.assign(len, 0.0);
  std::vector<double> out_raw(len, 0.0);

  for (int idx = 0; idx < len; ++idx) {
    const int t = idx - T;
    ts_out[idx] = t;
    const int k_start = std::abs(t);

    for (int ks = k_start; ks <= N; ks += 2) {
      const double pk = R::dbinom(ks, N, p, false);
      const double m = (static_cast<double>(ks) + static_cast<double>(t)) / 2.0;
      const double qm = R::dbinom(m, ks, 0.5, false);
      const double base = pk * qm;
      out_raw[idx] += base;
      if (base > 0.0) {
        const double surv = std::pow(1.0 - mr_lethality, static_cast<double>(ks));
        prob_out[idx] += base * surv;
      }
    }
  }

  const double total_raw = std::accumulate(out_raw.begin(), out_raw.end(), 0.0);
  mass_dropped = std::max(0.0, 1.0 - total_raw);

  if (total_raw > 0.0) {
    for (int idx = 0; idx < len; ++idx) {
      const double out_raw_norm = out_raw[idx] / total_raw;
      if (out_raw_norm > 0.0) {
        const double surv_ratio = prob_out[idx] / (out_raw_norm * total_raw);
        prob_out[idx] = out_raw_norm * surv_ratio;
      } else {
        prob_out[idx] = 0.0;
      }
    }
  }
}

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

} // namespace

// [[Rcpp::export]]
List cpp_pr_delta_vec(
    int N,
    double p,
    double eps_tail = 1e-8,
    double mr_lethality = 0.9
) {
  std::vector<int> ts;
  std::vector<double> prob;
  double mass_dropped = 0.0;
  pr_delta_internal(N, p, eps_tail, mr_lethality, ts, prob, mass_dropped);

  return List::create(
    _["ts"] = IntegerVector(ts.begin(), ts.end()),
    _["prob"] = NumericVector(prob.begin(), prob.end()),
    _["mass_dropped"] = mass_dropped
  );
}

// [[Rcpp::export]]
List cpp_build_B_total_triplet(
    int Nmin,
    int Nmax,
    NumericVector p_vec,
    NumericVector mr_lethality,
    std::string boundary = "drop",
    double eps_tail = 1e-8
) {
  const int R = Nmax - Nmin + 1;
  if (R <= 0) stop("Nmax must be >= Nmin");

  const int p_len = p_vec.size();
  const int mr_len = mr_lethality.size();
  if (!(p_len == 1 || p_len == R)) stop("p_vec length must be 1 or R");
  if (!(mr_len == 1 || mr_len == R)) stop("mr_lethality length must be 1 or R");

  const int bmode = boundary_mode(boundary);

  std::vector<int> ii;
  std::vector<int> jj;
  std::vector<double> xx;
  ii.reserve(static_cast<size_t>(R) * 12);
  jj.reserve(static_cast<size_t>(R) * 12);
  xx.reserve(static_cast<size_t>(R) * 12);

  for (int col = 0; col < R; ++col) {
    const int N = Nmin + col;
    const double pN = (p_len == 1) ? p_vec[0] : p_vec[col];
    const double mr = (mr_len == 1) ? mr_lethality[0] : mr_lethality[col];

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    pr_delta_internal(N, pN, eps_tail, mr, ts, pr, mass_dropped);

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
List cpp_build_B_WGD_triplet(
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    std::string boundary = "drop"
) {
  const int R0 = N0max - N0min + 1;
  const int R1 = N1max - N1min + 1;
  if (R0 <= 0 || R1 <= 0) stop("Nmax must be >= Nmin for both layers");

  const int bmode = boundary_mode(boundary);

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
      2.0,
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
