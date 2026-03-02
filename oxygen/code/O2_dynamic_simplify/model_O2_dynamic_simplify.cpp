#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

namespace {

inline double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

inline double clamp_o2_pct(double x) {
  if (x < 0.0) return 0.0;
  if (x > 100.0) return 100.0;
  return x;
}

inline double sigmoid01(double z) {
  return 1.0 / (1.0 + std::exp(-z));
}

inline int quantize_o2_key_1e3(double o2_pct) {
  const double o2_use = clamp_o2_pct(o2_pct);
  int k = static_cast<int>(std::llround(o2_use * 1000.0));
  if (k < 0) k = 0;
  if (k > 100000) k = 100000;
  return k;
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

void o2simps_pr_delta_internal(
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
List cpp_o2simps_pr_delta_vec(
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

  o2simps_pr_delta_internal(
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
NumericVector cpp_o2simps_o2_window_supply(
    NumericVector Ntot,
    std::string curve_type = "gompertz",
    double O2_cap = 5.0,
    double o2_init = 0.5,
    double o2_rate = 1.0,
    double o2_shape_v = 1.0,
    double o2_anchor_N = 1e6,
    double o2_logN_eps = 1.0
) {
  const int n = Ntot.size();
  NumericVector out(n);

  const bool use_glogistic = (curve_type == "glogistic");
  if (!(curve_type == "gompertz" || use_glogistic)) {
    stop("curve_type must be one of: gompertz, glogistic");
  }

  const double O2_cap_use = clamp_o2_pct(O2_cap);
  const double eps_o2 = 1e-9;
  const double o2_init_use = std::max(eps_o2, std::min(O2_cap_use - eps_o2, o2_init));
  const double rate_use = (std::isfinite(o2_rate) && o2_rate > 0.0) ? o2_rate : 1.0;
  const double v_use = (std::isfinite(o2_shape_v) && o2_shape_v > 0.0) ? o2_shape_v : 1.0;
  const double anchor_use = (std::isfinite(o2_anchor_N) && o2_anchor_N >= 0.0) ? o2_anchor_N : 1e6;
  const double eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;
  const double x0 = std::log10(anchor_use + eps_use);

  const double bg = std::log(-std::log(o2_init_use / std::max(O2_cap_use, eps_o2)));
  const double bl = std::log(std::pow(std::max(O2_cap_use, eps_o2) / o2_init_use, v_use) - 1.0);
  for (int i = 0; i < n; ++i) {
    const double n_raw = Ntot[i];
    const double n_use = (std::isfinite(n_raw) && n_raw > 0.0) ? n_raw : 0.0;
    const double x = std::log10(n_use + eps_use) - x0;
    double o2 = O2_cap_use;
    if (use_glogistic) {
      o2 = O2_cap_use / std::pow(1.0 + std::exp(-rate_use * x + bl), 1.0 / v_use);
    } else {
      o2 = O2_cap_use * std::exp(-std::exp(-rate_use * x + bg));
    }
    out[i] = clamp_o2_pct(o2);
  }

  return out;
}

// [[Rcpp::export]]
List cpp_o2simps_build_B_total_triplet(
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
    o2simps_pr_delta_internal(
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
List cpp_o2simps_build_B_WGD_triplet(
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

namespace {

inline void append_block_with_boundary(
    int Np,
    int row_min,
    int row_max,
    int row_offset_1based,
    int col_1based,
    double value,
    int bmode,
    std::vector<int>& ii,
    std::vector<int>& jj,
    std::vector<double>& xx
) {
  if (value == 0.0) return;

  int Np_use = Np;
  if (Np_use < row_min || Np_use > row_max) {
    if (bmode == 0) return;
    Np_use = std::max(std::min(Np_use, row_max), row_min);
  }

  const int row_local_1based = Np_use - row_min + 1;
  ii.push_back(row_offset_1based + row_local_1based - 1);
  jj.push_back(col_1based);
  xx.push_back(value);
}

inline double resolve_pmis_for_o2(
    double O2_pct,
    bool has_p_misseg,
    double p_misseg,
    double k_o_mis,
    bool has_pmis_endpoints,
    double pmis_O2_0,
    double pmis_O2_1,
    double p_const
) {
  if (has_p_misseg) {
    const double p0 = std::isfinite(p_misseg) ? p_misseg : 0.0;
    const double k_use = (std::isfinite(k_o_mis) && k_o_mis > 0.0) ? k_o_mis : 1e-12;
    const double frac = O2_pct / (O2_pct + k_use);
    return clamp01(p0 * (1.0 - frac));
  }

  if (has_pmis_endpoints) {
    if (!(pmis_O2_0 > 0.0 && pmis_O2_1 > 0.0)) return 0.0;
    const double frac = clamp01(O2_pct / 100.0);
    const double logp = (1.0 - frac) * std::log10(pmis_O2_0) + frac * std::log10(pmis_O2_1);
    return clamp01(std::pow(10.0, logp));
  }

  return clamp01(p_const);
}

} // namespace

// [[Rcpp::export]]
List cpp_o2simps_build_G_for_o2_triplet(
    double O2,
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    double lam_min,
    double lam_max,
    double k_o,
    bool has_p_misseg,
    double p_misseg,
    double k_o_mis,
    bool has_pmis_endpoints,
    double pmis_O2_0,
    double pmis_O2_1,
    double p_const,
    double p_wgd,
    std::string boundary = "drop",
    double eps_tail = 1e-8,
    double beta_buffer = 0.0,
    double n_exp = 1.0,
    double smax = 1.0,
    int N_unit = 22
) {
  const int R0 = N0max - N0min + 1;
  const int R1 = N1max - N1min + 1;
  if (R0 <= 0 || R1 <= 0) stop("Nmax must be >= Nmin for both layers");

  const int bmode = boundary_mode(boundary);

  const double O2_use = clamp_o2_pct(O2);
  const double lam_min_use = std::isfinite(lam_min) ? lam_min : 1.0;
  const double lam_max_use = std::isfinite(lam_max) ? lam_max : lam_min_use;
  const double k_o_use = (std::isfinite(k_o) && k_o > 0.0) ? k_o : 1e-12;
  const double frac = O2_use / (O2_use + k_o_use);
  double lam = lam_min_use + (lam_max_use - lam_min_use) * frac;
  if (!std::isfinite(lam) || lam < 0.0) lam = 0.0;

  const double p_mis = resolve_pmis_for_o2(
    O2_use,
    has_p_misseg,
    p_misseg,
    k_o_mis,
    has_pmis_endpoints,
    pmis_O2_0,
    pmis_O2_1,
    p_const
  );
  const double pw = clamp01(p_wgd);

  std::vector<int> ii;
  std::vector<int> jj;
  std::vector<double> xx;
  ii.reserve(static_cast<size_t>(R0 + R1) * 20);
  jj.reserve(static_cast<size_t>(R0 + R1) * 20);
  xx.reserve(static_cast<size_t>(R0 + R1) * 20);

  for (int c0 = 0; c0 < R0; ++c0) {
    const int N = N0min + c0;
    const int col_1based = c0 + 1;

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    o2simps_pr_delta_internal(
      N,
      p_mis,
      eps_tail,
      beta_buffer,
      n_exp,
      smax,
      N_unit,
      ts,
      pr,
      mass_dropped
    );
    (void)mass_dropped;

    const double scale_pre = lam * (1.0 - pw);
    for (size_t k = 0; k < ts.size(); ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;
      if (t == 0) {
        append_block_with_boundary(
          N,
          N0min,
          N0max,
          1,
          col_1based,
          scale_pre * (2.0 * w),
          bmode,
          ii,
          jj,
          xx
        );
      } else {
        append_block_with_boundary(
          N + t,
          N0min,
          N0max,
          1,
          col_1based,
          scale_pre * w,
          bmode,
          ii,
          jj,
          xx
        );
        append_block_with_boundary(
          N - t,
          N0min,
          N0max,
          1,
          col_1based,
          scale_pre * w,
          bmode,
          ii,
          jj,
          xx
        );
      }
    }

    ii.push_back(col_1based);
    jj.push_back(col_1based);
    xx.push_back(-lam);

    append_block_with_boundary(
      2 * N,
      N1min,
      N1max,
      R0 + 1,
      col_1based,
      lam * pw,
      bmode,
      ii,
      jj,
      xx
    );
  }

  for (int c1 = 0; c1 < R1; ++c1) {
    const int N = N1min + c1;
    const int col_1based = R0 + c1 + 1;

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    o2simps_pr_delta_internal(
      N,
      p_mis,
      eps_tail,
      beta_buffer,
      n_exp,
      smax,
      N_unit,
      ts,
      pr,
      mass_dropped
    );
    (void)mass_dropped;

    for (size_t k = 0; k < ts.size(); ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;
      if (t == 0) {
        append_block_with_boundary(
          N,
          N1min,
          N1max,
          R0 + 1,
          col_1based,
          lam * (2.0 * w),
          bmode,
          ii,
          jj,
          xx
        );
      } else {
        append_block_with_boundary(
          N + t,
          N1min,
          N1max,
          R0 + 1,
          col_1based,
          lam * w,
          bmode,
          ii,
          jj,
          xx
        );
        append_block_with_boundary(
          N - t,
          N1min,
          N1max,
          R0 + 1,
          col_1based,
          lam * w,
          bmode,
          ii,
          jj,
          xx
        );
      }
    }

    ii.push_back(col_1based);
    jj.push_back(col_1based);
    xx.push_back(-lam);
  }

  return List::create(
    _["i"] = IntegerVector(ii.begin(), ii.end()),
    _["j"] = IntegerVector(jj.begin(), jj.end()),
    _["x"] = NumericVector(xx.begin(), xx.end()),
    _["nrow"] = R0 + R1,
    _["ncol"] = R0 + R1
  );
}

namespace {

struct SparseCacheEntry {
  std::vector<int> row0;
  std::vector<int> col0;
  std::vector<double> val;
};

template <typename T>
inline void hash_combine_cpp(std::size_t& seed, const T& value) {
  seed ^= std::hash<T>{}(value) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

inline std::uint64_t bits_of_double_cpp(double x) {
  std::uint64_t out = 0ULL;
  std::memcpy(&out, &x, sizeof(double));
  return out;
}

inline std::size_t g_cache_signature_cpp(
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    double lam_min,
    double lam_max,
    double k_o,
    bool has_p_misseg,
    double p_misseg,
    double k_o_mis,
    bool has_pmis_endpoints,
    double pmis_O2_0,
    double pmis_O2_1,
    double p_const,
    double p_wgd,
    const std::string& boundary,
    double eps_tail,
    double beta_buffer,
    double n_exp,
    double smax,
    int N_unit
) {
  std::size_t seed = 0ULL;
  hash_combine_cpp(seed, N0min);
  hash_combine_cpp(seed, N0max);
  hash_combine_cpp(seed, N1min);
  hash_combine_cpp(seed, N1max);
  hash_combine_cpp(seed, bits_of_double_cpp(lam_min));
  hash_combine_cpp(seed, bits_of_double_cpp(lam_max));
  hash_combine_cpp(seed, bits_of_double_cpp(k_o));
  hash_combine_cpp(seed, has_p_misseg ? 1 : 0);
  hash_combine_cpp(seed, bits_of_double_cpp(p_misseg));
  hash_combine_cpp(seed, bits_of_double_cpp(k_o_mis));
  hash_combine_cpp(seed, has_pmis_endpoints ? 1 : 0);
  hash_combine_cpp(seed, bits_of_double_cpp(pmis_O2_0));
  hash_combine_cpp(seed, bits_of_double_cpp(pmis_O2_1));
  hash_combine_cpp(seed, bits_of_double_cpp(p_const));
  hash_combine_cpp(seed, bits_of_double_cpp(p_wgd));
  hash_combine_cpp(seed, boundary);
  hash_combine_cpp(seed, bits_of_double_cpp(eps_tail));
  hash_combine_cpp(seed, bits_of_double_cpp(beta_buffer));
  hash_combine_cpp(seed, bits_of_double_cpp(n_exp));
  hash_combine_cpp(seed, bits_of_double_cpp(smax));
  hash_combine_cpp(seed, N_unit);
  return seed;
}

inline double vector_sum_cpp(const std::vector<double>& x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}

inline void sparse_mv_cpp(
    const SparseCacheEntry& G,
    const std::vector<double>& x,
    std::vector<double>& y
) {
  std::fill(y.begin(), y.end(), 0.0);
  const size_t nnz = G.val.size();
  for (size_t e = 0; e < nnz; ++e) {
    y[static_cast<size_t>(G.row0[e])] += G.val[e] * x[static_cast<size_t>(G.col0[e])];
  }
}

inline double o2_window_supply_scalar_cpp(
    double Ntot,
    const std::string& curve_type,
    double O2_cap,
    double o2_init,
    double o2_rate,
    double o2_shape_v,
    double o2_anchor_N,
    double o2_logN_eps
) {
  const double n_use = (std::isfinite(Ntot) && Ntot > 0.0) ? Ntot : 0.0;
  const bool use_glogistic = (curve_type == "glogistic");
  if (!(curve_type == "gompertz" || use_glogistic)) {
    stop("curve_type must be one of: gompertz, glogistic");
  }
  const double O2_cap_use = clamp_o2_pct(O2_cap);
  const double eps_o2 = 1e-9;
  const double o2_init_use = std::max(eps_o2, std::min(O2_cap_use - eps_o2, o2_init));
  const double rate_use = (std::isfinite(o2_rate) && o2_rate > 0.0) ? o2_rate : 1.0;
  const double v_use = (std::isfinite(o2_shape_v) && o2_shape_v > 0.0) ? o2_shape_v : 1.0;
  const double anchor_use = (std::isfinite(o2_anchor_N) && o2_anchor_N >= 0.0) ? o2_anchor_N : 1e6;
  const double eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;
  const double x0 = std::log10(anchor_use + eps_use);
  const double x = std::log10(n_use + eps_use) - x0;
  if (use_glogistic) {
    const double bl = std::log(std::pow(std::max(O2_cap_use, eps_o2) / o2_init_use, v_use) - 1.0);
    const double o2 = O2_cap_use / std::pow(1.0 + std::exp(-rate_use * x + bl), 1.0 / v_use);
    return clamp_o2_pct(o2);
  }
  const double bg = std::log(-std::log(o2_init_use / std::max(O2_cap_use, eps_o2)));
  const double o2 = O2_cap_use * std::exp(-std::exp(-rate_use * x + bg));
  return clamp_o2_pct(o2);
}

inline SparseCacheEntry build_sparse_cache_entry_from_triplet(const List& tri) {
  IntegerVector ii = tri["i"];
  IntegerVector jj = tri["j"];
  NumericVector xx = tri["x"];
  const int n = xx.size();
  if (ii.size() != n || jj.size() != n) {
    stop("Triplet i/j/x length mismatch.");
  }
  SparseCacheEntry out;
  out.row0.resize(static_cast<size_t>(n));
  out.col0.resize(static_cast<size_t>(n));
  out.val.resize(static_cast<size_t>(n));
  for (int e = 0; e < n; ++e) {
    out.row0[static_cast<size_t>(e)] = ii[e] - 1;
    out.col0[static_cast<size_t>(e)] = jj[e] - 1;
    out.val[static_cast<size_t>(e)] = xx[e];
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
List cpp_o2simps_simulate_one(
    NumericVector init_state,
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    IntegerVector obs_steps,
    int sim_end_step,
    double DT,
    double dose,
    double dose_ref,
    double treat_day,
    bool fit_treatment,
    double alpha,
    double gamma,
    double tx_mult_min,
    std::string crowding,
    double K,
    double min_pop,
    double O2_cap,
    bool o2_feedback,
    std::string o2_curve_type,
    double o2_init,
    double o2_rate,
    double o2_shape_v,
    double o2_anchor_N,
    double o2_logN_eps,
    double lam_min,
    double lam_max,
    double k_o,
    bool has_p_misseg,
    double p_misseg,
    double k_o_mis,
    bool has_pmis_endpoints,
    double pmis_O2_0,
    double pmis_O2_1,
    double p_const,
    double p_wgd,
    std::string boundary,
    double eps_tail,
    double beta_buffer,
    double n_exp,
    double smax,
    int N_unit,
    NumericVector vol_by_N,
    double burden_floor
) {
  const int R0 = N0max - N0min + 1;
  const int R1 = N1max - N1min + 1;
  if (R0 <= 0 || R1 <= 0) stop("Nmax must be >= Nmin for both layers.");
  if (init_state.size() != (R0 + R1)) stop("init_state length mismatch.");
  if (vol_by_N.size() != R0) stop("vol_by_N length mismatch.");

  const bool crowd_logistic = (crowding == "logistic");
  const bool crowd_gompertz = (crowding == "gompertz");
  if (!crowd_logistic && !crowd_gompertz) stop("crowding must be logistic or gompertz.");

  const double DT_use = (std::isfinite(DT) && DT > 0.0) ? DT : 0.5;
  const double K_use = (std::isfinite(K) && K > 0.0) ? K : 1e12;
  const double min_pop_use = (std::isfinite(min_pop) && min_pop > 0.0) ? min_pop : 1e-12;
  const double dose_ref_use = (std::isfinite(dose_ref) && dose_ref > 0.0) ? dose_ref : 30.0;
  const double dose_use = (std::isfinite(dose) && dose > 0.0) ? dose : 0.0;
  const double tx_min_use = clamp01(tx_mult_min);
  const double burden_floor_use = (std::isfinite(burden_floor) && burden_floor >= 0.0) ? burden_floor : 0.0;

  std::vector<int> obs_v(obs_steps.size());
  for (int i = 0; i < obs_steps.size(); ++i) obs_v[static_cast<size_t>(i)] = obs_steps[i];
  std::vector<int> step_unique = obs_v;
  std::sort(step_unique.begin(), step_unique.end());
  step_unique.erase(std::unique(step_unique.begin(), step_unique.end()), step_unique.end());
  const int max_obs_step = step_unique.empty() ? 0 : step_unique.back();
  const int final_step = std::max(sim_end_step, max_obs_step);

  std::unordered_map<int, int> step_to_idx;
  step_to_idx.reserve(step_unique.size() * 2 + 1);
  for (size_t i = 0; i < step_unique.size(); ++i) {
    step_to_idx[step_unique[i]] = static_cast<int>(i);
  }

  std::vector<double> Ntot_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_at_step(step_unique.size(), NA_REAL);

  std::vector<double> v(init_state.begin(), init_state.end());
  std::vector<double> growth(static_cast<size_t>(R0 + R1), 0.0);

  // Shared across scenario calls in the same worker process.
  // We keep one active parameter signature at a time so cache is reused
  // within one objective (same params), then reset when params change.
  static std::size_t active_sig = std::numeric_limits<std::size_t>::max();
  static std::unordered_map<int, SparseCacheEntry> shared_G_cache;

  const std::size_t cur_sig = g_cache_signature_cpp(
    N0min,
    N0max,
    N1min,
    N1max,
    lam_min,
    lam_max,
    k_o,
    has_p_misseg,
    p_misseg,
    k_o_mis,
    has_pmis_endpoints,
    pmis_O2_0,
    pmis_O2_1,
    p_const,
    p_wgd,
    boundary,
    eps_tail,
    beta_buffer,
    n_exp,
    smax,
    N_unit
  );
  if (cur_sig != active_sig) {
    shared_G_cache.clear();
    shared_G_cache.reserve(256);
    active_sig = cur_sig;
  }

  const double O2_cap_use = clamp_o2_pct(O2_cap);
  const bool o2_glogistic = (o2_curve_type == "glogistic");
  if (!(o2_curve_type == "gompertz" || o2_glogistic)) {
    stop("o2_curve_type must be one of: gompertz, glogistic");
  }
  const double o2_init_use = (std::isfinite(o2_init) ? o2_init : 0.5);
  const double o2_rate_use = (std::isfinite(o2_rate) ? o2_rate : 1.0);
  const double o2_shape_v_use = (std::isfinite(o2_shape_v) ? o2_shape_v : 1.0);
  const double o2_anchor_use = (std::isfinite(o2_anchor_N) && o2_anchor_N >= 0.0) ? o2_anchor_N : 1e6;
  const double o2_eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;

  for (int step = 0; step <= final_step; ++step) {
    auto it_obs = step_to_idx.find(step);
    if (it_obs != step_to_idx.end()) {
      const int idx = it_obs->second;
      const double Ntot_now = vector_sum_cpp(v);
      Ntot_at_step[static_cast<size_t>(idx)] = Ntot_now;
      double burden_now = 0.0;
      for (int i = 0; i < R0; ++i) {
        const double n_i = v[static_cast<size_t>(i)] + v[static_cast<size_t>(R0 + i)];
        burden_now += n_i * vol_by_N[i];
      }
      Vmm3_at_step[static_cast<size_t>(idx)] = burden_now;
    }
    if (step >= final_step) break;

    const double t = static_cast<double>(step) * DT_use;
    double tx_mult = 1.0;
    if (fit_treatment) {
      if (!(t < treat_day) && dose_use > 0.0) {
        double dose_scaled = dose_use / dose_ref_use;
        if (!std::isfinite(dose_scaled) || dose_scaled < 0.0) dose_scaled = 0.0;
        tx_mult = std::exp(-alpha * std::pow(dose_scaled, gamma));
      } else {
        tx_mult = 1.0;
      }
      if (!std::isfinite(tx_mult)) tx_mult = tx_min_use;
      if (tx_mult < tx_min_use) tx_mult = tx_min_use;
      if (tx_mult > 1.0) tx_mult = 1.0;
    }

    const double Ntot = vector_sum_cpp(v);
    double O2_eff = O2_cap_use;
    if (o2_feedback) {
      O2_eff = o2_window_supply_scalar_cpp(
        Ntot,
        o2_curve_type,
        O2_cap_use,
        o2_init_use,
        o2_rate_use,
        o2_shape_v_use,
        o2_anchor_use,
        o2_eps_use
      );
    }
    O2_eff = clamp_o2_pct(O2_eff);

    const int gkey = quantize_o2_key_1e3(O2_eff);
    auto itG = shared_G_cache.find(gkey);
    if (itG == shared_G_cache.end()) {
      const List tri = cpp_o2simps_build_G_for_o2_triplet(
        O2_eff,
        N0min,
        N0max,
        N1min,
        N1max,
        lam_min,
        lam_max,
        k_o,
        has_p_misseg,
        p_misseg,
        k_o_mis,
        has_pmis_endpoints,
        pmis_O2_0,
        pmis_O2_1,
        p_const,
        p_wgd,
        boundary,
        eps_tail,
        beta_buffer,
        n_exp,
        smax,
        N_unit
      );
      SparseCacheEntry entry = build_sparse_cache_entry_from_triplet(tri);
      auto insert_res = shared_G_cache.emplace(gkey, std::move(entry));
      itG = insert_res.first;
    }

    sparse_mv_cpp(itG->second, v, growth);
    const double crowd = crowd_logistic ? std::max(0.0, 1.0 - Ntot / K_use) : std::exp(-Ntot / K_use);
    const double scalar = DT_use * crowd * tx_mult;
    for (size_t i = 0; i < v.size(); ++i) {
      const double next_v = v[i] + scalar * growth[i];
      if (!std::isfinite(next_v) || next_v < 0.0) {
        v[i] = 0.0;
      } else {
        v[i] = next_v;
      }
    }
    if (vector_sum_cpp(v) <= min_pop_use) break;
  }

  NumericVector Ntot_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_obs(obs_v.size(), NA_REAL);
  for (int i = 0; i < static_cast<int>(obs_v.size()); ++i) {
    auto it = step_to_idx.find(obs_v[static_cast<size_t>(i)]);
    if (it == step_to_idx.end()) {
      Ntot_obs[i] = min_pop_use;
      Vmm3_obs[i] = burden_floor_use;
      continue;
    }
    const int idx = it->second;
    double nv = Ntot_at_step[static_cast<size_t>(idx)];
    double bv = Vmm3_at_step[static_cast<size_t>(idx)];
    if (!std::isfinite(nv)) nv = min_pop_use;
    if (!std::isfinite(bv)) bv = burden_floor_use;
    Ntot_obs[i] = nv;
    Vmm3_obs[i] = bv;
  }

  NumericVector frac_N(R0, 0.0);
  double total_frac = 0.0;
  for (int i = 0; i < R0; ++i) {
    const double val = v[static_cast<size_t>(i)] + v[static_cast<size_t>(R0 + i)];
    frac_N[i] = val;
    total_frac += val;
  }
  if (total_frac > 0.0 && std::isfinite(total_frac)) {
    for (int i = 0; i < R0; ++i) frac_N[i] = frac_N[i] / total_frac;
  } else {
    const double u = 1.0 / static_cast<double>(R0);
    for (int i = 0; i < R0; ++i) frac_N[i] = u;
  }

  return List::create(
    _["Ntot_obs"] = Ntot_obs,
    _["Vmm3_obs"] = Vmm3_obs,
    _["frac_N"] = frac_N
  );
}
