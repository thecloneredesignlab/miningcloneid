#include <RcppEigen.h>
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
// [[Rcpp::depends(RcppEigen)]]

namespace {

// -----------------------------------------------------------------------------
// Function: clamp01
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - x: Input value or vector to process.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double clamp01(double x) {
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

// -----------------------------------------------------------------------------
// Function: clamp_o2_pct
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - x: Input value or vector to process.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double clamp_o2_pct(double x) {
  if (x < 0.0) return 0.0;
  if (x > 100.0) return 100.0;
  return x;
}

// -----------------------------------------------------------------------------
// Function: clamp_o2_pct_cap
// Purpose: Clamp oxygen percentage to [0, o2_cap] where o2_cap is in [0,100].
// Parameters:
//   - x: Input oxygen value in percent.
//   - o2_cap: Maximum oxygen cap in percent.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double clamp_o2_pct_cap(double x, double o2_cap) {
  const double cap = clamp_o2_pct(o2_cap);
  if (x < 0.0) return 0.0;
  if (x > cap) return cap;
  return x;
}

// -----------------------------------------------------------------------------
// Function: sigmoid01
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - z: Function-specific input argument.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double sigmoid01(double z) {
  return 1.0 / (1.0 + std::exp(-z));
}

// -----------------------------------------------------------------------------
// Function: quantize_o2_key
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - o2_pct: Oxygen level in percent scale (0-100).
//   - bin_pct: Function-specific input argument.
// Returns:
//   int return value containing the computed result.
// -----------------------------------------------------------------------------
inline int quantize_o2_key(double o2_pct, double bin_pct) {
  const double o2_use = clamp_o2_pct(o2_pct);
  const double bin_use = (std::isfinite(bin_pct) && bin_pct > 0.0) ? bin_pct : 1e-3;
  const double raw = o2_use / bin_use;
  const double cap = static_cast<double>(std::numeric_limits<int>::max() / 4);
  const double clamped = std::min(std::max(raw, 0.0), cap);
  return static_cast<int>(std::llround(clamped));
}

// -----------------------------------------------------------------------------
// Function: boundary_mode
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
// Returns:
//   int return value containing the computed result.
// -----------------------------------------------------------------------------
inline int boundary_mode(const std::string& boundary) {
  if (boundary == "drop") return 0;
  if (boundary == "absorb_minmax") return 1;
  stop("boundary must be one of: drop, absorb_minmax");
}

// -----------------------------------------------------------------------------
// Function: append_with_boundary
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - Np: Function-specific input argument.
//   - row_min: Minimum valid row state index for sparse insertion.
//   - row_max: Maximum valid row state index for sparse insertion.
//   - col_1based: 1-based sparse column index for sparse insertion.
//   - value: Numeric transition value to append.
//   - bmode: Encoded boundary mode used internally by C++ helpers.
//   - ii: Sparse-triplet row index accumulator.
//   - jj: Sparse-triplet column index accumulator.
//   - xx: Sparse-triplet value accumulator.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
// Function: o2invivo_pr_delta_internal
// Purpose: Compute missegregation delta-kernel probabilities over ploidy shifts.
// Parameters:
//   - N: Ploidy state value or chromosome-copy count.
//   - p: Missegregation probability parameter.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - beta_buffer: Buffer exponent controlling ploidy-dependence of missegregation survival.
//   - n_exp: Exponent controlling ploidy scaling in buffering term.
//   - smax: Maximum survival factor for missegregation events.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
//   - ts_out: Function-specific input argument.
//   - prob_out: Function-specific input argument.
//   - mass_dropped: Function-specific input argument.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
void o2invivo_pr_delta_internal(
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

// -----------------------------------------------------------------------------
// Function: cpp_o2invivo_pr_delta_vec
// Purpose: Compute missegregation delta-kernel probabilities over ploidy shifts.
// Parameters:
//   - N: Ploidy state value or chromosome-copy count.
//   - p: Missegregation probability parameter.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - beta_buffer: Buffer exponent controlling ploidy-dependence of missegregation survival.
//   - n_exp: Exponent controlling ploidy scaling in buffering term.
//   - smax: Maximum survival factor for missegregation events.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2invivo_pr_delta_vec(
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

  o2invivo_pr_delta_internal(
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

// -----------------------------------------------------------------------------
// Function: cpp_o2invivo_o2_window_supply
// Purpose: Compute oxygen supply fraction/level from burden under the selected O2 model.
// Parameters:
//   - Ntot: Total predicted cell count (or burden proxy) at current time.
//   - O2_base: Function-specific input argument.
//   - o2_min: Minimum oxygen floor in percent.
//   - h_down: Shape exponent controlling steepness of O2 down-regulation.
//   - K_down: Burden scale controlling O2 down-regulation window.
//   - A_ang: Angiogenesis amplitude parameter in dynamic O2 window model.
//   - m_on: Center of angiogenesis activation window in log-burden space.
//   - m_off: Function-specific input argument.
//   - s_on: Slope of angiogenesis activation edge.
//   - s_off: Slope of angiogenesis deactivation edge.
//   - o2_logN_eps: Function-specific input argument.
// Returns:
//   NumericVector return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cpp_o2invivo_o2_window_supply(
    NumericVector Ntot,
    double O2_base = 5.0,
    double o2_cap = 5.0,
    double o2_min = 0.0,
    double h_down = 1.0,
    double K_down = 1e12,
    double A_ang = 0.0,
    double m_on = 9.0,
    double m_off = 10.0,
    double s_on = 0.3,
    double s_off = 0.3,
    double o2_logN_eps = 1.0
) {
  const int n = Ntot.size();
  NumericVector out(n);

  const double o2_cap_use = clamp_o2_pct(o2_cap);
  const double O2_base_use = clamp_o2_pct_cap(O2_base, o2_cap_use);
  const double o2_floor = clamp_o2_pct(o2_min);
  const double o2_floor_use = std::min(o2_floor, o2_cap_use);
  const double h_use = (std::isfinite(h_down) && h_down > 0.0) ? h_down : 1.0;
  const double K_down_use = (std::isfinite(K_down) && K_down > 0.0) ? K_down : 1e12;
  const double A_ang_use = clamp_o2_pct(A_ang);
  const double m_on_use = std::isfinite(m_on) ? m_on : 9.0;
  const double m_off_use = (std::isfinite(m_off) && m_off > m_on_use) ? m_off : (m_on_use + 1.0);
  const double s_on_use = (std::isfinite(s_on) && s_on > 0.0) ? s_on : 0.3;
  const double s_off_use = (std::isfinite(s_off) && s_off > 0.0) ? s_off : 0.3;
  const double eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;

  for (int i = 0; i < n; ++i) {
    const double n_raw = Ntot[i];
    const double n_use = (std::isfinite(n_raw) && n_raw > 0.0) ? n_raw : 0.0;
    const double x = std::log10(n_use + eps_use);
    const double down = o2_floor_use + (O2_base_use - o2_floor_use) / (1.0 + std::pow(n_use / K_down_use, h_use));
    const double win = sigmoid01((x - m_on_use) / s_on_use) * (1.0 - sigmoid01((x - m_off_use) / s_off_use));
    out[i] = clamp_o2_pct_cap(down + A_ang_use * win, o2_cap_use);
  }

  return out;
}

// -----------------------------------------------------------------------------
// Function: cpp_o2invivo_build_B_total_triplet
// Purpose: Build total missegregation transition operator on ploidy grid.
// Parameters:
//   - Nmin: Minimum ploidy state on source grid.
//   - Nmax: Maximum ploidy state on source grid.
//   - p_vec: State-specific missegregation probability vector.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - beta_buffer: Buffer exponent controlling ploidy-dependence of missegregation survival.
//   - n_exp: Exponent controlling ploidy scaling in buffering term.
//   - smax: Maximum survival factor for missegregation events.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2invivo_build_B_total_triplet(
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
    o2invivo_pr_delta_internal(
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

// -----------------------------------------------------------------------------
// Function: cpp_o2invivo_build_B_WGD_triplet
// Purpose: Build WGD transition operator between source and doubled-ploidy grids.
// Parameters:
//   - N0min: Minimum ploidy state on pre-WGD source grid.
//   - N0max: Maximum ploidy state on pre-WGD source grid.
//   - N1min: Minimum ploidy state on WGD target grid.
//   - N1max: Maximum ploidy state on WGD target grid.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - wgd_value: Function-specific input argument.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2invivo_build_B_WGD_triplet(
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

// -----------------------------------------------------------------------------
// Function: append_block_with_boundary
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - Np: Function-specific input argument.
//   - row_min: Minimum valid row state index for sparse insertion.
//   - row_max: Maximum valid row state index for sparse insertion.
//   - row_offset_1based: Function-specific input argument.
//   - col_1based: 1-based sparse column index for sparse insertion.
//   - value: Numeric transition value to append.
//   - bmode: Encoded boundary mode used internally by C++ helpers.
//   - ii: Sparse-triplet row index accumulator.
//   - jj: Sparse-triplet column index accumulator.
//   - xx: Sparse-triplet value accumulator.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
// Function: resolve_pmis_for_o2
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - O2_pct: Function-specific input argument.
//   - has_p_misseg: Function-specific input argument.
//   - p_misseg: Function-specific input argument.
//   - k_o_mis: Oxygen-sensitivity parameter for missegregation rate.
//   - has_pmis_endpoints: Function-specific input argument.
//   - pmis_O2_0: Function-specific input argument.
//   - pmis_O2_1: Function-specific input argument.
//   - p_const: Function-specific input argument.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
// Function: cpp_o2invivo_build_G_for_o2_triplet
// Purpose: Build generator matrix at the current oxygen/burden condition.
// Parameters:
//   - O2: Oxygen level used by model rate functions.
//   - N0min: Minimum ploidy state on pre-WGD source grid.
//   - N0max: Maximum ploidy state on pre-WGD source grid.
//   - N1min: Minimum ploidy state on WGD target grid.
//   - N1max: Maximum ploidy state on WGD target grid.
//   - lam_min: Lower asymptote of proliferation rate.
//   - lam_max: Upper asymptote of proliferation rate.
//   - k_o: Oxygen-sensitivity parameter for proliferation rate.
//   - has_p_misseg: Function-specific input argument.
//   - p_misseg: Function-specific input argument.
//   - k_o_mis: Oxygen-sensitivity parameter for missegregation rate.
//   - has_pmis_endpoints: Function-specific input argument.
//   - pmis_O2_0: Function-specific input argument.
//   - pmis_O2_1: Function-specific input argument.
//   - p_const: Function-specific input argument.
//   - p_wgd: Function-specific input argument.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - beta_buffer: Buffer exponent controlling ploidy-dependence of missegregation survival.
//   - n_exp: Exponent controlling ploidy scaling in buffering term.
//   - smax: Maximum survival factor for missegregation events.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2invivo_build_G_for_o2_triplet(
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
    int N_unit = 22,
    double beta_size = 0.0,
    double alpha_o2 = 0.0,
    double o2_ref_pct = 0.0,
    double gamma_growth = 1.0,
    double mu_hp = 0.0
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
  double lam_base = lam_min_use + (lam_max_use - lam_min_use) * frac;
  if (!std::isfinite(lam_base) || lam_base < 0.0) lam_base = 0.0;
  const double beta_size_use = std::isfinite(beta_size) ? beta_size : 0.0;
  const double alpha_o2_use = (std::isfinite(alpha_o2) && alpha_o2 >= 0.0) ? alpha_o2 : 0.0;
  const double o2_ref_use = clamp_o2_pct(o2_ref_pct);
  const double gamma_growth_use = (std::isfinite(gamma_growth) && gamma_growth > 0.0) ? gamma_growth : 1.0;
  const double mu_hp_use = (std::isfinite(mu_hp) && mu_hp > 0.0) ? mu_hp : 0.0;
  const int N_unit_safe = (N_unit > 0) ? N_unit : 1;
  auto lam_for_N = [&](int N_state) -> double {
    if (lam_base <= 0.0) return 0.0;
    const double d = std::max(0.0, static_cast<double>(N_state) / static_cast<double>(N_unit_safe) - 2.0);
    const double size_penalty = std::exp(-beta_size_use * std::pow(d, gamma_growth_use));
    const double hypoxia_deficit = std::max(0.0, o2_ref_use - O2_use);
    const double hypoxia_penalty = 1.0 / (1.0 + alpha_o2_use * d * hypoxia_deficit);
    double lam_eff = lam_base * size_penalty * hypoxia_penalty;
    if (!std::isfinite(lam_eff) || lam_eff < 0.0) lam_eff = 0.0;
    return lam_eff;
  };
  auto mu_for_N = [&](int N_state) -> double {
    if (mu_hp_use <= 0.0) return 0.0;
    const double d = std::max(0.0, static_cast<double>(N_state) / static_cast<double>(N_unit_safe) - 2.0);
    if (d <= 0.0) return 0.0;
    const double hypoxia_deficit = std::max(0.0, o2_ref_use - O2_use);
    if (hypoxia_deficit <= 0.0) return 0.0;
    double mu_eff = mu_hp_use * d * hypoxia_deficit;
    if (!std::isfinite(mu_eff) || mu_eff < 0.0) mu_eff = 0.0;
    return mu_eff;
  };

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
    const double lam_N = lam_for_N(N);
    const double mu_N = mu_for_N(N);
    const int col_1based = c0 + 1;

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    o2invivo_pr_delta_internal(
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

    const double scale_pre = lam_N * (1.0 - pw);
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
    xx.push_back(-(lam_N + mu_N));

    append_block_with_boundary(
      2 * N,
      N1min,
      N1max,
      R0 + 1,
      col_1based,
      lam_N * pw,
      bmode,
      ii,
      jj,
      xx
    );
  }

  for (int c1 = 0; c1 < R1; ++c1) {
    const int N = N1min + c1;
    const double lam_N = lam_for_N(N);
    const double mu_N = mu_for_N(N);
    const int col_1based = R0 + c1 + 1;

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    o2invivo_pr_delta_internal(
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
          lam_N * (2.0 * w),
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
          lam_N * w,
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
          lam_N * w,
          bmode,
          ii,
          jj,
          xx
        );
      }
    }

    ii.push_back(col_1based);
    jj.push_back(col_1based);
    xx.push_back(-(lam_N + mu_N));
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

using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

struct SparseCacheEntry {
  SpMat mat;
};

template <typename T>
// -----------------------------------------------------------------------------
// Function: hash_combine_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - seed: Random-seed value used for reproducible initialization.
//   - value: Numeric transition value to append.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
inline void hash_combine_cpp(std::size_t& seed, const T& value) {
  seed ^= std::hash<T>{}(value) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
}

// -----------------------------------------------------------------------------
// Function: bits_of_double_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - x: Input value or vector to process.
// Returns:
//   std::uint64_t return value containing the computed result.
// -----------------------------------------------------------------------------
inline std::uint64_t bits_of_double_cpp(double x) {
  std::uint64_t out = 0ULL;
  std::memcpy(&out, &x, sizeof(double));
  return out;
}

// -----------------------------------------------------------------------------
// Function: g_cache_signature_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - N0min: Minimum ploidy state on pre-WGD source grid.
//   - N0max: Maximum ploidy state on pre-WGD source grid.
//   - N1min: Minimum ploidy state on WGD target grid.
//   - N1max: Maximum ploidy state on WGD target grid.
//   - lam_min: Lower asymptote of proliferation rate.
//   - lam_max: Upper asymptote of proliferation rate.
//   - k_o: Oxygen-sensitivity parameter for proliferation rate.
//   - has_p_misseg: Function-specific input argument.
//   - p_misseg: Function-specific input argument.
//   - k_o_mis: Oxygen-sensitivity parameter for missegregation rate.
//   - has_pmis_endpoints: Function-specific input argument.
//   - pmis_O2_0: Function-specific input argument.
//   - pmis_O2_1: Function-specific input argument.
//   - p_const: Function-specific input argument.
//   - p_wgd: Function-specific input argument.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - beta_buffer: Buffer exponent controlling ploidy-dependence of missegregation survival.
//   - n_exp: Exponent controlling ploidy scaling in buffering term.
//   - smax: Maximum survival factor for missegregation events.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
// Returns:
//   std::size_t return value containing the computed result.
// -----------------------------------------------------------------------------
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
    double beta_size,
    double alpha_o2,
    double o2_ref_pct,
    double gamma_growth,
    double mu_hp,
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
  hash_combine_cpp(seed, bits_of_double_cpp(beta_size));
  hash_combine_cpp(seed, bits_of_double_cpp(alpha_o2));
  hash_combine_cpp(seed, bits_of_double_cpp(o2_ref_pct));
  hash_combine_cpp(seed, bits_of_double_cpp(gamma_growth));
  hash_combine_cpp(seed, bits_of_double_cpp(mu_hp));
  hash_combine_cpp(seed, N_unit);
  return seed;
}

// -----------------------------------------------------------------------------
// Function: vector_sum_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - x: Input value or vector to process.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double vector_sum_cpp(const std::vector<double>& x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}

// -----------------------------------------------------------------------------
// Function: sparse_mv_cpp
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - G: Generator or transition matrix used for state propagation.
//   - x: Input value or vector to process.
//   - y: Function-specific input argument.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
inline void sparse_mv_cpp(
    const SparseCacheEntry& G,
    const std::vector<double>& x,
    std::vector<double>& y
) {
  if (static_cast<int>(x.size()) != G.mat.cols() || static_cast<int>(y.size()) != G.mat.rows()) {
    stop("Sparse matvec dimension mismatch.");
  }
  Eigen::Map<const Eigen::VectorXd> xmap(x.data(), static_cast<int>(x.size()));
  Eigen::Map<Eigen::VectorXd> ymap(y.data(), static_cast<int>(y.size()));
  ymap.noalias() = G.mat * xmap;
}

// -----------------------------------------------------------------------------
// Function: o2_window_supply_scalar_cpp
// Purpose: Compute oxygen supply fraction/level from burden under the selected O2 model.
// Parameters:
//   - Ntot: Total predicted cell count (or burden proxy) at current time.
//   - O2_base: Function-specific input argument.
//   - o2_min: Minimum oxygen floor in percent.
//   - h_down: Shape exponent controlling steepness of O2 down-regulation.
//   - K_down: Burden scale controlling O2 down-regulation window.
//   - A_ang: Angiogenesis amplitude parameter in dynamic O2 window model.
//   - m_on: Center of angiogenesis activation window in log-burden space.
//   - m_off: Function-specific input argument.
//   - s_on: Slope of angiogenesis activation edge.
//   - s_off: Slope of angiogenesis deactivation edge.
//   - o2_logN_eps: Function-specific input argument.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double o2_window_supply_scalar_cpp(
    double Ntot,
    double O2_base,
    double o2_cap,
    double o2_min,
    double h_down,
    double K_down,
    double A_ang,
    double m_on,
    double m_off,
    double s_on,
    double s_off,
    double o2_logN_eps
) {
  const double n_use = (std::isfinite(Ntot) && Ntot > 0.0) ? Ntot : 0.0;
  const double o2_cap_use = clamp_o2_pct(o2_cap);
  const double O2_base_use = clamp_o2_pct_cap(O2_base, o2_cap_use);
  const double o2_floor = clamp_o2_pct(o2_min);
  const double o2_floor_use = std::min(o2_floor, o2_cap_use);
  const double h_use = (std::isfinite(h_down) && h_down > 0.0) ? h_down : 1.0;
  const double K_down_use = (std::isfinite(K_down) && K_down > 0.0) ? K_down : 1e12;
  const double A_ang_use = clamp_o2_pct(A_ang);
  const double m_on_use = std::isfinite(m_on) ? m_on : 9.0;
  const double m_off_use = (std::isfinite(m_off) && m_off > m_on_use) ? m_off : (m_on_use + 1.0);
  const double s_on_use = (std::isfinite(s_on) && s_on > 0.0) ? s_on : 0.3;
  const double s_off_use = (std::isfinite(s_off) && s_off > 0.0) ? s_off : 0.3;
  const double eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;

  const double x = std::log10(n_use + eps_use);
  const double down = o2_floor_use + (O2_base_use - o2_floor_use) / (1.0 + std::pow(n_use / K_down_use, h_use));
  const double win = sigmoid01((x - m_on_use) / s_on_use) * (1.0 - sigmoid01((x - m_off_use) / s_off_use));
  return clamp_o2_pct_cap(down + A_ang_use * win, o2_cap_use);
}

// -----------------------------------------------------------------------------
// Function: build_sparse_cache_entry_from_triplet
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - tri: Function-specific input argument.
// Returns:
//   SparseCacheEntry return value containing the computed result.
// -----------------------------------------------------------------------------
inline SparseCacheEntry build_sparse_cache_entry_from_triplet(const List& tri) {
  IntegerVector ii = tri["i"];
  IntegerVector jj = tri["j"];
  NumericVector xx = tri["x"];
  const int nrow = as<int>(tri["nrow"]);
  const int ncol = as<int>(tri["ncol"]);
  const int n = xx.size();
  if (ii.size() != n || jj.size() != n) {
    stop("Triplet i/j/x length mismatch.");
  }
  SparseCacheEntry out;
  out.mat.resize(nrow, ncol);
  std::vector<Eigen::Triplet<double, int>> triplets;
  triplets.reserve(static_cast<size_t>(n));
  for (int e = 0; e < n; ++e) {
    const int r = ii[e] - 1;
    const int c = jj[e] - 1;
    if (r < 0 || r >= nrow || c < 0 || c >= ncol) {
      stop("Triplet index out of bounds.");
    }
    triplets.emplace_back(r, c, xx[e]);
  }
  out.mat.setFromTriplets(triplets.begin(), triplets.end());
  out.mat.makeCompressed();
  return out;
}

} // namespace

// -----------------------------------------------------------------------------
// Function: cpp_o2invivo_simulate_one
// Purpose: Run one forward simulation trajectory for a single scenario.
// Parameters:
//   - init_state: Function-specific input argument.
//   - N0min: Minimum ploidy state on pre-WGD source grid.
//   - N0max: Maximum ploidy state on pre-WGD source grid.
//   - N1min: Minimum ploidy state on WGD target grid.
//   - N1max: Maximum ploidy state on WGD target grid.
//   - obs_steps: Function-specific input argument.
//   - sim_end_step: Function-specific input argument.
//   - DT: Function-specific input argument.
//   - dose: Function-specific input argument.
//   - dose_ref: Function-specific input argument.
//   - treat_day: Function-specific input argument.
//   - fit_treatment: Logical flag indicating whether treatment-effect parameters are optimized.
//   - alpha: Function-specific input argument.
//   - gamma: Function-specific input argument.
//   - tx_mult_min: Function-specific input argument.
//   - crowding: Function-specific input argument.
//   - K: Function-specific input argument.
//   - min_pop: Function-specific input argument.
//   - O2_base: Function-specific input argument.
//   - o2_feedback: Function-specific input argument.
//   - o2_min: Minimum oxygen floor in percent.
//   - h_O2: Function-specific input argument.
//   - K_down: Burden scale controlling O2 down-regulation window.
//   - A_ang: Angiogenesis amplitude parameter in dynamic O2 window model.
//   - m_on: Center of angiogenesis activation window in log-burden space.
//   - m_off: Function-specific input argument.
//   - s_on: Slope of angiogenesis activation edge.
//   - s_off: Slope of angiogenesis deactivation edge.
//   - tau_O2: Relaxation time constant controlling lag from O2 target to O2 effective.
//   - o2_logN_eps: Function-specific input argument.
//   - o2_cache_bin_pct: Function-specific input argument.
//   - o2_cache_hysteresis_pct: Function-specific input argument.
//   - o2_cache_profile: Function-specific input argument.
//   - lam_min: Lower asymptote of proliferation rate.
//   - lam_max: Upper asymptote of proliferation rate.
//   - k_o: Oxygen-sensitivity parameter for proliferation rate.
//   - has_p_misseg: Function-specific input argument.
//   - p_misseg: Function-specific input argument.
//   - k_o_mis: Oxygen-sensitivity parameter for missegregation rate.
//   - has_pmis_endpoints: Function-specific input argument.
//   - pmis_O2_0: Function-specific input argument.
//   - pmis_O2_1: Function-specific input argument.
//   - p_const: Function-specific input argument.
//   - p_wgd: Function-specific input argument.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - beta_buffer: Buffer exponent controlling ploidy-dependence of missegregation survival.
//   - n_exp: Exponent controlling ploidy scaling in buffering term.
//   - smax: Maximum survival factor for missegregation events.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
//   - vol_by_N: Optional precomputed per-state cell volume lookup.
//   - burden_floor: Function-specific input argument.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2invivo_simulate_one(
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
    double O2_base,
    double o2_cap,
    bool o2_feedback,
    double o2_min,
    double h_O2,
    double K_down,
    double A_ang,
    double m_on,
    double m_off,
    double s_on,
    double s_off,
    double tau_O2,
    double o2_logN_eps,
    double o2_cache_bin_pct,
    double o2_cache_hysteresis_pct,
    bool o2_cache_profile,
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
    double beta_size,
    double alpha_o2,
    double o2_ref_pct,
    double gamma_growth,
    double mu_hp,
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
    beta_size,
    alpha_o2,
    o2_ref_pct,
    gamma_growth,
    mu_hp,
    N_unit
  );
  if (cur_sig != active_sig) {
    shared_G_cache.clear();
    shared_G_cache.reserve(256);
    active_sig = cur_sig;
  }

  const double o2_cap_use = clamp_o2_pct(o2_cap);
  const double O2_base_use = clamp_o2_pct_cap(O2_base, o2_cap_use);
  const double o2_min_use = std::min(clamp_o2_pct(o2_min), o2_cap_use);
  const double h_use = (std::isfinite(h_O2) && h_O2 > 0.0) ? h_O2 : 1.0;
  const double K_down_use = (std::isfinite(K_down) && K_down > 0.0) ? K_down : 1e12;
  const double A_ang_use = clamp_o2_pct(A_ang);
  const double m_on_use = std::isfinite(m_on) ? m_on : 9.0;
  const double m_off_use = (std::isfinite(m_off) && m_off > m_on_use) ? m_off : (m_on_use + 1.0);
  const double s_on_use = (std::isfinite(s_on) && s_on > 0.0) ? s_on : 0.3;
  const double s_off_use = (std::isfinite(s_off) && s_off > 0.0) ? s_off : 0.3;
  const double tau_use = (std::isfinite(tau_O2) && tau_O2 > 0.0) ? tau_O2 : 2.0;
  const double alpha_tau = 1.0 - std::exp(-DT_use / tau_use);
  const double o2_eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;
  const double o2_bin_use = (std::isfinite(o2_cache_bin_pct) && o2_cache_bin_pct > 0.0) ? o2_cache_bin_pct : 1e-3;
  const double o2_hyst_use = (std::isfinite(o2_cache_hysteresis_pct) && o2_cache_hysteresis_pct >= 0.0) ? o2_cache_hysteresis_pct : 0.0;
  (void) o2_cache_profile;
  int cache_g_build = 0;
  int cache_g_hit = 0;
  int cache_g_hysteresis = 0;
  bool has_last_key = false;
  int last_key = 0;
  double last_o2_eff = 0.0;
  double O2_state = O2_base_use;
  if (o2_feedback) {
    O2_state = o2_window_supply_scalar_cpp(
      vector_sum_cpp(v),
      O2_base_use,
      o2_cap_use,
      o2_min_use,
      h_use,
      K_down_use,
      A_ang_use,
      m_on_use,
      m_off_use,
      s_on_use,
      s_off_use,
      o2_eps_use
    );
    O2_state = clamp_o2_pct_cap(O2_state, o2_cap_use);
  }

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
    double O2_target = O2_base_use;
    if (o2_feedback) {
      O2_target = o2_window_supply_scalar_cpp(
        Ntot,
        O2_base_use,
        o2_cap_use,
        o2_min_use,
        h_use,
        K_down_use,
        A_ang_use,
        m_on_use,
        m_off_use,
        s_on_use,
        s_off_use,
        o2_eps_use
      );
    }
    O2_target = clamp_o2_pct_cap(O2_target, o2_cap_use);
    O2_state = O2_state + alpha_tau * (O2_target - O2_state);
    double O2_eff = clamp_o2_pct_cap(O2_state, o2_cap_use);

    int gkey = quantize_o2_key(O2_eff, o2_bin_use);
    if (o2_hyst_use > 0.0 && has_last_key && std::abs(O2_eff - last_o2_eff) <= o2_hyst_use) {
      gkey = last_key;
      ++cache_g_hysteresis;
    }
    auto itG = shared_G_cache.find(gkey);
    if (itG == shared_G_cache.end()) {
      const List tri = cpp_o2invivo_build_G_for_o2_triplet(
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
        N_unit,
        beta_size,
        alpha_o2,
        o2_ref_pct,
        gamma_growth,
        mu_hp
      );
      SparseCacheEntry entry = build_sparse_cache_entry_from_triplet(tri);
      auto insert_res = shared_G_cache.emplace(gkey, std::move(entry));
      itG = insert_res.first;
      ++cache_g_build;
    } else {
      ++cache_g_hit;
    }
    has_last_key = true;
    last_key = gkey;
    last_o2_eff = O2_eff;

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
    _["frac_N"] = frac_N,
    _["cache_g_build"] = cache_g_build,
    _["cache_g_hit"] = cache_g_hit,
    _["cache_g_hysteresis"] = cache_g_hysteresis,
    _["cache_bin_pct"] = o2_bin_use,
    _["cache_hysteresis_pct"] = o2_hyst_use
  );
}

// -----------------------------------------------------------------------------
// Function: cpp_o2invivo_objective_components_map
// Purpose: Compute MAP objective components using log-normal burden likelihood
//   and continuous single-cell ploidy mixture likelihood.
// Parameters:
//   - ploidy_z_list: Per-tumor continuous single-cell ploidy observations.
//   - mu_by_N: Representative ploidy value for each modeled N state.
//   - sigma_burden: Log-normal observation SD for burden.
//   - sigma_ploidy: Gaussian observation SD for single-cell ploidy.
// Returns:
//   List return value containing per-modality mean NLL components.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2invivo_objective_components_map(
    IntegerVector cohort_code,
    NumericVector dose_vec,
    NumericVector treat_day_vec,
    List obs_steps_list,
    IntegerVector sim_end_step_vec,
    List obs_burden_list,
    List keep_burden_list,
    List ploidy_z_list,
    NumericVector mu_by_N,
    double sigma_burden,
    double sigma_ploidy,
    NumericVector init_state_2N,
    NumericVector init_state_4N,
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    double DT,
    double dose_ref,
    bool fit_treatment,
    double alpha,
    double gamma,
    double tx_mult_min,
    std::string crowding,
    double K,
    double min_pop,
    double O2_base,
    double o2_cap,
    bool o2_feedback,
    double o2_min,
    double h_O2,
    double K_down,
    double A_ang,
    double m_on,
    double m_off,
    double s_on,
    double s_off,
    double tau_O2,
    double o2_logN_eps,
    double o2_cache_bin_pct,
    double o2_cache_hysteresis_pct,
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
    NumericVector growth_oxy_par,
    NumericVector vol_by_N,
    double burden_floor,
    double burden_log_eps
) {
  const int n_sc = cohort_code.size();
  if (dose_vec.size() != n_sc || treat_day_vec.size() != n_sc ||
      obs_steps_list.size() != n_sc || sim_end_step_vec.size() != n_sc ||
      obs_burden_list.size() != n_sc || keep_burden_list.size() != n_sc ||
      ploidy_z_list.size() != n_sc) {
    stop("Scenario containers must have consistent length.");
  }

  const double log_eps_use =
    (std::isfinite(burden_log_eps) && burden_log_eps > 0.0) ? burden_log_eps : 1e-15;
  const double sigma_b_use =
    (std::isfinite(sigma_burden) && sigma_burden > 0.0) ? sigma_burden : 0.35;
  const double sigma_p_use =
    (std::isfinite(sigma_ploidy) && sigma_ploidy > 0.0) ? sigma_ploidy : 0.08;
  const double prob_eps = 1e-300;
  const double beta_size_use =
    (growth_oxy_par.size() > 0 && std::isfinite(growth_oxy_par[0])) ? static_cast<double>(growth_oxy_par[0]) : 0.0;
  const double alpha_o2_use =
    (growth_oxy_par.size() > 1 && std::isfinite(growth_oxy_par[1])) ? static_cast<double>(growth_oxy_par[1]) : 0.0;
  const double o2_ref_pct_use =
    (growth_oxy_par.size() > 2 && std::isfinite(growth_oxy_par[2])) ? static_cast<double>(growth_oxy_par[2]) : 0.0;
  const double gamma_growth_use =
    (growth_oxy_par.size() > 3 && std::isfinite(growth_oxy_par[3])) ? static_cast<double>(growth_oxy_par[3]) : 1.0;
  const double mu_hp_use =
    (growth_oxy_par.size() > 4 && std::isfinite(growth_oxy_par[4])) ? static_cast<double>(growth_oxy_par[4]) : 0.0;

  std::vector<double> burden_losses;
  std::vector<double> ploidy_losses;
  burden_losses.reserve(static_cast<size_t>(n_sc));
  ploidy_losses.reserve(static_cast<size_t>(n_sc));

  int cache_g_build = 0;
  int cache_g_hit = 0;
  int cache_g_hysteresis = 0;

  for (int i = 0; i < n_sc; ++i) {
    const int cohort = cohort_code[i];
    NumericVector init_state = (cohort == 0) ? init_state_2N : init_state_4N;
    IntegerVector obs_steps = as<IntegerVector>(obs_steps_list[i]);
    NumericVector obs_burden = as<NumericVector>(obs_burden_list[i]);
    LogicalVector keep_day = as<LogicalVector>(keep_burden_list[i]);
    NumericVector ploidy_z = as<NumericVector>(ploidy_z_list[i]);

    List sim = cpp_o2invivo_simulate_one(
      init_state,
      N0min,
      N0max,
      N1min,
      N1max,
      obs_steps,
      sim_end_step_vec[i],
      DT,
      dose_vec[i],
      dose_ref,
      treat_day_vec[i],
      fit_treatment,
      alpha,
      gamma,
      tx_mult_min,
      crowding,
      K,
      min_pop,
      O2_base,
      o2_cap,
      o2_feedback,
      o2_min,
      h_O2,
      K_down,
      A_ang,
      m_on,
      m_off,
      s_on,
      s_off,
      tau_O2,
      o2_logN_eps,
      o2_cache_bin_pct,
      o2_cache_hysteresis_pct,
      false,
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
      N_unit,
      beta_size_use,
      alpha_o2_use,
      o2_ref_pct_use,
      gamma_growth_use,
      mu_hp_use,
      vol_by_N,
      burden_floor
    );

    NumericVector pred_burden = sim["Vmm3_obs"];
    NumericVector frac_N = sim["frac_N"];
    cache_g_build += as<int>(sim["cache_g_build"]);
    cache_g_hit += as<int>(sim["cache_g_hit"]);
    cache_g_hysteresis += as<int>(sim["cache_g_hysteresis"]);

    if (mu_by_N.size() != frac_N.size()) {
      stop("mu_by_N length must equal simulated terminal state vector length.");
    }

    // Burden log-normal NLL per tumor (mean across available time points).
    const int nb = std::min(obs_burden.size(), pred_burden.size());
    double burden_nll_sum = 0.0;
    int burden_n = 0;
    for (int j = 0; j < nb; ++j) {
      const bool keepj = (keep_day.size() == nb) ? static_cast<bool>(keep_day[j]) : true;
      if (!keepj) continue;
      const double obs = obs_burden[j];
      const double pred = pred_burden[j];
      if (!std::isfinite(obs) || !std::isfinite(pred) || obs <= 0.0 || pred <= 0.0) continue;
      const double resid = std::log(std::max(obs, log_eps_use)) - std::log(std::max(pred, log_eps_use));
      const double z = resid / sigma_b_use;
      burden_nll_sum += std::log(sigma_b_use) + 0.5 * std::log(2.0 * 3.14159265358979323846) + 0.5 * z * z;
      ++burden_n;
    }
    if (burden_n > 0) {
      burden_losses.push_back(burden_nll_sum / static_cast<double>(burden_n));
    }

    // Continuous single-cell ploidy NLL per tumor:
    // p(z) = sum_j pi_j * Normal(z; mu_j, sigma_ploidy^2), then average -log p(z).
    double p_sum = 0.0;
    for (int j = 0; j < frac_N.size(); ++j) {
      const double pv = frac_N[j];
      if (std::isfinite(pv) && pv > 0.0) p_sum += pv;
    }
    if (ploidy_z.size() > 0 && p_sum > 0.0) {
      double ploidy_nll_sum = 0.0;
      int ploidy_n = 0;
      for (int c = 0; c < ploidy_z.size(); ++c) {
        const double z_obs = ploidy_z[c];
        if (!std::isfinite(z_obs)) continue;
        double max_log = -std::numeric_limits<double>::infinity();
        std::vector<double> comp_log;
        comp_log.reserve(static_cast<size_t>(frac_N.size()));
        for (int j = 0; j < frac_N.size(); ++j) {
          const double pv = frac_N[j];
          if (!std::isfinite(pv) || pv <= 0.0) continue;
          const double pj = pv / p_sum;
          const double mu_j = mu_by_N[j];
          if (!std::isfinite(mu_j)) continue;
          const double log_comp =
            std::log(std::max(pj, prob_eps)) + R::dnorm4(z_obs, mu_j, sigma_p_use, /*give_log=*/1);
          comp_log.push_back(log_comp);
          if (log_comp > max_log) max_log = log_comp;
        }
        if (comp_log.empty() || !std::isfinite(max_log)) continue;
        double sum_exp = 0.0;
        for (size_t t = 0; t < comp_log.size(); ++t) {
          sum_exp += std::exp(comp_log[t] - max_log);
        }
        if (!std::isfinite(sum_exp) || sum_exp <= 0.0) continue;
        const double log_mix = max_log + std::log(sum_exp);
        ploidy_nll_sum += -log_mix;
        ++ploidy_n;
      }
      if (ploidy_n > 0) {
        ploidy_losses.push_back(ploidy_nll_sum / static_cast<double>(ploidy_n));
      }
    }
  }

  const double L_b = burden_losses.empty()
    ? 0.0
    : std::accumulate(burden_losses.begin(), burden_losses.end(), 0.0) /
        static_cast<double>(burden_losses.size());
  const double L_p = ploidy_losses.empty()
    ? 0.0
    : std::accumulate(ploidy_losses.begin(), ploidy_losses.end(), 0.0) /
        static_cast<double>(ploidy_losses.size());

  return List::create(
    _["L_b"] = L_b,
    _["L_p"] = L_p,
    _["n_burden"] = static_cast<int>(burden_losses.size()),
    _["n_ploidy"] = static_cast<int>(ploidy_losses.size()),
    _["cache_g_build"] = cache_g_build,
    _["cache_g_hit"] = cache_g_hit,
    _["cache_g_hysteresis"] = cache_g_hysteresis
  );
}


#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_o2invivo_pr_delta_vec
List cpp_o2invivo_pr_delta_vec(int N, double p, double eps_tail, double beta_buffer, double n_exp, double smax, int N_unit);
RcppExport SEXP sourceCpp_1_cpp_o2invivo_pr_delta_vec(SEXP NSEXP, SEXP pSEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(cpp_o2invivo_pr_delta_vec(N, p, eps_tail, beta_buffer, n_exp, smax, N_unit));
    return rcpp_result_gen;
END_RCPP
}
// cpp_o2invivo_o2_window_supply
NumericVector cpp_o2invivo_o2_window_supply(NumericVector Ntot, double O2_base, double o2_cap, double o2_min, double h_down, double K_down, double A_ang, double m_on, double m_off, double s_on, double s_off, double o2_logN_eps);
RcppExport SEXP sourceCpp_1_cpp_o2invivo_o2_window_supply(SEXP NtotSEXP, SEXP O2_baseSEXP, SEXP o2_capSEXP, SEXP o2_minSEXP, SEXP h_downSEXP, SEXP K_downSEXP, SEXP A_angSEXP, SEXP m_onSEXP, SEXP m_offSEXP, SEXP s_onSEXP, SEXP s_offSEXP, SEXP o2_logN_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Ntot(NtotSEXP);
    Rcpp::traits::input_parameter< double >::type O2_base(O2_baseSEXP);
    Rcpp::traits::input_parameter< double >::type o2_cap(o2_capSEXP);
    Rcpp::traits::input_parameter< double >::type o2_min(o2_minSEXP);
    Rcpp::traits::input_parameter< double >::type h_down(h_downSEXP);
    Rcpp::traits::input_parameter< double >::type K_down(K_downSEXP);
    Rcpp::traits::input_parameter< double >::type A_ang(A_angSEXP);
    Rcpp::traits::input_parameter< double >::type m_on(m_onSEXP);
    Rcpp::traits::input_parameter< double >::type m_off(m_offSEXP);
    Rcpp::traits::input_parameter< double >::type s_on(s_onSEXP);
    Rcpp::traits::input_parameter< double >::type s_off(s_offSEXP);
    Rcpp::traits::input_parameter< double >::type o2_logN_eps(o2_logN_epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_o2invivo_o2_window_supply(Ntot, O2_base, o2_cap, o2_min, h_down, K_down, A_ang, m_on, m_off, s_on, s_off, o2_logN_eps));
    return rcpp_result_gen;
END_RCPP
}
// cpp_o2invivo_build_B_total_triplet
List cpp_o2invivo_build_B_total_triplet(int Nmin, int Nmax, NumericVector p_vec, std::string boundary, double eps_tail, double beta_buffer, double n_exp, double smax, int N_unit);
RcppExport SEXP sourceCpp_1_cpp_o2invivo_build_B_total_triplet(SEXP NminSEXP, SEXP NmaxSEXP, SEXP p_vecSEXP, SEXP boundarySEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP) {
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
    rcpp_result_gen = Rcpp::wrap(cpp_o2invivo_build_B_total_triplet(Nmin, Nmax, p_vec, boundary, eps_tail, beta_buffer, n_exp, smax, N_unit));
    return rcpp_result_gen;
END_RCPP
}
// cpp_o2invivo_build_B_WGD_triplet
List cpp_o2invivo_build_B_WGD_triplet(int N0min, int N0max, int N1min, int N1max, std::string boundary, double wgd_value);
RcppExport SEXP sourceCpp_1_cpp_o2invivo_build_B_WGD_triplet(SEXP N0minSEXP, SEXP N0maxSEXP, SEXP N1minSEXP, SEXP N1maxSEXP, SEXP boundarySEXP, SEXP wgd_valueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N0min(N0minSEXP);
    Rcpp::traits::input_parameter< int >::type N0max(N0maxSEXP);
    Rcpp::traits::input_parameter< int >::type N1min(N1minSEXP);
    Rcpp::traits::input_parameter< int >::type N1max(N1maxSEXP);
    Rcpp::traits::input_parameter< std::string >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< double >::type wgd_value(wgd_valueSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_o2invivo_build_B_WGD_triplet(N0min, N0max, N1min, N1max, boundary, wgd_value));
    return rcpp_result_gen;
END_RCPP
}
// cpp_o2invivo_build_G_for_o2_triplet
List cpp_o2invivo_build_G_for_o2_triplet(double O2, int N0min, int N0max, int N1min, int N1max, double lam_min, double lam_max, double k_o, bool has_p_misseg, double p_misseg, double k_o_mis, bool has_pmis_endpoints, double pmis_O2_0, double pmis_O2_1, double p_const, double p_wgd, std::string boundary, double eps_tail, double beta_buffer, double n_exp, double smax, int N_unit, double beta_size, double alpha_o2, double o2_ref_pct, double gamma_growth, double mu_hp);
RcppExport SEXP sourceCpp_1_cpp_o2invivo_build_G_for_o2_triplet(SEXP O2SEXP, SEXP N0minSEXP, SEXP N0maxSEXP, SEXP N1minSEXP, SEXP N1maxSEXP, SEXP lam_minSEXP, SEXP lam_maxSEXP, SEXP k_oSEXP, SEXP has_p_missegSEXP, SEXP p_missegSEXP, SEXP k_o_misSEXP, SEXP has_pmis_endpointsSEXP, SEXP pmis_O2_0SEXP, SEXP pmis_O2_1SEXP, SEXP p_constSEXP, SEXP p_wgdSEXP, SEXP boundarySEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP, SEXP beta_sizeSEXP, SEXP alpha_o2SEXP, SEXP o2_ref_pctSEXP, SEXP gamma_growthSEXP, SEXP mu_hpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type O2(O2SEXP);
    Rcpp::traits::input_parameter< int >::type N0min(N0minSEXP);
    Rcpp::traits::input_parameter< int >::type N0max(N0maxSEXP);
    Rcpp::traits::input_parameter< int >::type N1min(N1minSEXP);
    Rcpp::traits::input_parameter< int >::type N1max(N1maxSEXP);
    Rcpp::traits::input_parameter< double >::type lam_min(lam_minSEXP);
    Rcpp::traits::input_parameter< double >::type lam_max(lam_maxSEXP);
    Rcpp::traits::input_parameter< double >::type k_o(k_oSEXP);
    Rcpp::traits::input_parameter< bool >::type has_p_misseg(has_p_missegSEXP);
    Rcpp::traits::input_parameter< double >::type p_misseg(p_missegSEXP);
    Rcpp::traits::input_parameter< double >::type k_o_mis(k_o_misSEXP);
    Rcpp::traits::input_parameter< bool >::type has_pmis_endpoints(has_pmis_endpointsSEXP);
    Rcpp::traits::input_parameter< double >::type pmis_O2_0(pmis_O2_0SEXP);
    Rcpp::traits::input_parameter< double >::type pmis_O2_1(pmis_O2_1SEXP);
    Rcpp::traits::input_parameter< double >::type p_const(p_constSEXP);
    Rcpp::traits::input_parameter< double >::type p_wgd(p_wgdSEXP);
    Rcpp::traits::input_parameter< std::string >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< double >::type eps_tail(eps_tailSEXP);
    Rcpp::traits::input_parameter< double >::type beta_buffer(beta_bufferSEXP);
    Rcpp::traits::input_parameter< double >::type n_exp(n_expSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< int >::type N_unit(N_unitSEXP);
    Rcpp::traits::input_parameter< double >::type beta_size(beta_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_o2(alpha_o2SEXP);
    Rcpp::traits::input_parameter< double >::type o2_ref_pct(o2_ref_pctSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_growth(gamma_growthSEXP);
    Rcpp::traits::input_parameter< double >::type mu_hp(mu_hpSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_o2invivo_build_G_for_o2_triplet(O2, N0min, N0max, N1min, N1max, lam_min, lam_max, k_o, has_p_misseg, p_misseg, k_o_mis, has_pmis_endpoints, pmis_O2_0, pmis_O2_1, p_const, p_wgd, boundary, eps_tail, beta_buffer, n_exp, smax, N_unit, beta_size, alpha_o2, o2_ref_pct, gamma_growth, mu_hp));
    return rcpp_result_gen;
END_RCPP
}
// cpp_o2invivo_simulate_one
List cpp_o2invivo_simulate_one(NumericVector init_state, int N0min, int N0max, int N1min, int N1max, IntegerVector obs_steps, int sim_end_step, double DT, double dose, double dose_ref, double treat_day, bool fit_treatment, double alpha, double gamma, double tx_mult_min, std::string crowding, double K, double min_pop, double O2_base, double o2_cap, bool o2_feedback, double o2_min, double h_O2, double K_down, double A_ang, double m_on, double m_off, double s_on, double s_off, double tau_O2, double o2_logN_eps, double o2_cache_bin_pct, double o2_cache_hysteresis_pct, bool o2_cache_profile, double lam_min, double lam_max, double k_o, bool has_p_misseg, double p_misseg, double k_o_mis, bool has_pmis_endpoints, double pmis_O2_0, double pmis_O2_1, double p_const, double p_wgd, std::string boundary, double eps_tail, double beta_buffer, double n_exp, double smax, int N_unit, double beta_size, double alpha_o2, double o2_ref_pct, double gamma_growth, double mu_hp, NumericVector vol_by_N, double burden_floor);
RcppExport SEXP sourceCpp_1_cpp_o2invivo_simulate_one(SEXP init_stateSEXP, SEXP N0minSEXP, SEXP N0maxSEXP, SEXP N1minSEXP, SEXP N1maxSEXP, SEXP obs_stepsSEXP, SEXP sim_end_stepSEXP, SEXP DTSEXP, SEXP doseSEXP, SEXP dose_refSEXP, SEXP treat_daySEXP, SEXP fit_treatmentSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP tx_mult_minSEXP, SEXP crowdingSEXP, SEXP KSEXP, SEXP min_popSEXP, SEXP O2_baseSEXP, SEXP o2_capSEXP, SEXP o2_feedbackSEXP, SEXP o2_minSEXP, SEXP h_O2SEXP, SEXP K_downSEXP, SEXP A_angSEXP, SEXP m_onSEXP, SEXP m_offSEXP, SEXP s_onSEXP, SEXP s_offSEXP, SEXP tau_O2SEXP, SEXP o2_logN_epsSEXP, SEXP o2_cache_bin_pctSEXP, SEXP o2_cache_hysteresis_pctSEXP, SEXP o2_cache_profileSEXP, SEXP lam_minSEXP, SEXP lam_maxSEXP, SEXP k_oSEXP, SEXP has_p_missegSEXP, SEXP p_missegSEXP, SEXP k_o_misSEXP, SEXP has_pmis_endpointsSEXP, SEXP pmis_O2_0SEXP, SEXP pmis_O2_1SEXP, SEXP p_constSEXP, SEXP p_wgdSEXP, SEXP boundarySEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP, SEXP beta_sizeSEXP, SEXP alpha_o2SEXP, SEXP o2_ref_pctSEXP, SEXP gamma_growthSEXP, SEXP mu_hpSEXP, SEXP vol_by_NSEXP, SEXP burden_floorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type init_state(init_stateSEXP);
    Rcpp::traits::input_parameter< int >::type N0min(N0minSEXP);
    Rcpp::traits::input_parameter< int >::type N0max(N0maxSEXP);
    Rcpp::traits::input_parameter< int >::type N1min(N1minSEXP);
    Rcpp::traits::input_parameter< int >::type N1max(N1maxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type obs_steps(obs_stepsSEXP);
    Rcpp::traits::input_parameter< int >::type sim_end_step(sim_end_stepSEXP);
    Rcpp::traits::input_parameter< double >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< double >::type dose(doseSEXP);
    Rcpp::traits::input_parameter< double >::type dose_ref(dose_refSEXP);
    Rcpp::traits::input_parameter< double >::type treat_day(treat_daySEXP);
    Rcpp::traits::input_parameter< bool >::type fit_treatment(fit_treatmentSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type tx_mult_min(tx_mult_minSEXP);
    Rcpp::traits::input_parameter< std::string >::type crowding(crowdingSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type min_pop(min_popSEXP);
    Rcpp::traits::input_parameter< double >::type O2_base(O2_baseSEXP);
    Rcpp::traits::input_parameter< double >::type o2_cap(o2_capSEXP);
    Rcpp::traits::input_parameter< bool >::type o2_feedback(o2_feedbackSEXP);
    Rcpp::traits::input_parameter< double >::type o2_min(o2_minSEXP);
    Rcpp::traits::input_parameter< double >::type h_O2(h_O2SEXP);
    Rcpp::traits::input_parameter< double >::type K_down(K_downSEXP);
    Rcpp::traits::input_parameter< double >::type A_ang(A_angSEXP);
    Rcpp::traits::input_parameter< double >::type m_on(m_onSEXP);
    Rcpp::traits::input_parameter< double >::type m_off(m_offSEXP);
    Rcpp::traits::input_parameter< double >::type s_on(s_onSEXP);
    Rcpp::traits::input_parameter< double >::type s_off(s_offSEXP);
    Rcpp::traits::input_parameter< double >::type tau_O2(tau_O2SEXP);
    Rcpp::traits::input_parameter< double >::type o2_logN_eps(o2_logN_epsSEXP);
    Rcpp::traits::input_parameter< double >::type o2_cache_bin_pct(o2_cache_bin_pctSEXP);
    Rcpp::traits::input_parameter< double >::type o2_cache_hysteresis_pct(o2_cache_hysteresis_pctSEXP);
    Rcpp::traits::input_parameter< bool >::type o2_cache_profile(o2_cache_profileSEXP);
    Rcpp::traits::input_parameter< double >::type lam_min(lam_minSEXP);
    Rcpp::traits::input_parameter< double >::type lam_max(lam_maxSEXP);
    Rcpp::traits::input_parameter< double >::type k_o(k_oSEXP);
    Rcpp::traits::input_parameter< bool >::type has_p_misseg(has_p_missegSEXP);
    Rcpp::traits::input_parameter< double >::type p_misseg(p_missegSEXP);
    Rcpp::traits::input_parameter< double >::type k_o_mis(k_o_misSEXP);
    Rcpp::traits::input_parameter< bool >::type has_pmis_endpoints(has_pmis_endpointsSEXP);
    Rcpp::traits::input_parameter< double >::type pmis_O2_0(pmis_O2_0SEXP);
    Rcpp::traits::input_parameter< double >::type pmis_O2_1(pmis_O2_1SEXP);
    Rcpp::traits::input_parameter< double >::type p_const(p_constSEXP);
    Rcpp::traits::input_parameter< double >::type p_wgd(p_wgdSEXP);
    Rcpp::traits::input_parameter< std::string >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< double >::type eps_tail(eps_tailSEXP);
    Rcpp::traits::input_parameter< double >::type beta_buffer(beta_bufferSEXP);
    Rcpp::traits::input_parameter< double >::type n_exp(n_expSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< int >::type N_unit(N_unitSEXP);
    Rcpp::traits::input_parameter< double >::type beta_size(beta_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_o2(alpha_o2SEXP);
    Rcpp::traits::input_parameter< double >::type o2_ref_pct(o2_ref_pctSEXP);
    Rcpp::traits::input_parameter< double >::type gamma_growth(gamma_growthSEXP);
    Rcpp::traits::input_parameter< double >::type mu_hp(mu_hpSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vol_by_N(vol_by_NSEXP);
    Rcpp::traits::input_parameter< double >::type burden_floor(burden_floorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_o2invivo_simulate_one(init_state, N0min, N0max, N1min, N1max, obs_steps, sim_end_step, DT, dose, dose_ref, treat_day, fit_treatment, alpha, gamma, tx_mult_min, crowding, K, min_pop, O2_base, o2_cap, o2_feedback, o2_min, h_O2, K_down, A_ang, m_on, m_off, s_on, s_off, tau_O2, o2_logN_eps, o2_cache_bin_pct, o2_cache_hysteresis_pct, o2_cache_profile, lam_min, lam_max, k_o, has_p_misseg, p_misseg, k_o_mis, has_pmis_endpoints, pmis_O2_0, pmis_O2_1, p_const, p_wgd, boundary, eps_tail, beta_buffer, n_exp, smax, N_unit, beta_size, alpha_o2, o2_ref_pct, gamma_growth, mu_hp, vol_by_N, burden_floor));
    return rcpp_result_gen;
END_RCPP
}
// cpp_o2invivo_objective_components_map
List cpp_o2invivo_objective_components_map(IntegerVector cohort_code, NumericVector dose_vec, NumericVector treat_day_vec, List obs_steps_list, IntegerVector sim_end_step_vec, List obs_burden_list, List keep_burden_list, List ploidy_z_list, NumericVector mu_by_N, double sigma_burden, double sigma_ploidy, NumericVector init_state_2N, NumericVector init_state_4N, int N0min, int N0max, int N1min, int N1max, double DT, double dose_ref, bool fit_treatment, double alpha, double gamma, double tx_mult_min, std::string crowding, double K, double min_pop, double O2_base, double o2_cap, bool o2_feedback, double o2_min, double h_O2, double K_down, double A_ang, double m_on, double m_off, double s_on, double s_off, double tau_O2, double o2_logN_eps, double o2_cache_bin_pct, double o2_cache_hysteresis_pct, double lam_min, double lam_max, double k_o, bool has_p_misseg, double p_misseg, double k_o_mis, bool has_pmis_endpoints, double pmis_O2_0, double pmis_O2_1, double p_const, double p_wgd, std::string boundary, double eps_tail, double beta_buffer, double n_exp, double smax, int N_unit, NumericVector growth_oxy_par, NumericVector vol_by_N, double burden_floor, double burden_log_eps);
RcppExport SEXP sourceCpp_1_cpp_o2invivo_objective_components_map(SEXP cohort_codeSEXP, SEXP dose_vecSEXP, SEXP treat_day_vecSEXP, SEXP obs_steps_listSEXP, SEXP sim_end_step_vecSEXP, SEXP obs_burden_listSEXP, SEXP keep_burden_listSEXP, SEXP ploidy_z_listSEXP, SEXP mu_by_NSEXP, SEXP sigma_burdenSEXP, SEXP sigma_ploidySEXP, SEXP init_state_2NSEXP, SEXP init_state_4NSEXP, SEXP N0minSEXP, SEXP N0maxSEXP, SEXP N1minSEXP, SEXP N1maxSEXP, SEXP DTSEXP, SEXP dose_refSEXP, SEXP fit_treatmentSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP tx_mult_minSEXP, SEXP crowdingSEXP, SEXP KSEXP, SEXP min_popSEXP, SEXP O2_baseSEXP, SEXP o2_capSEXP, SEXP o2_feedbackSEXP, SEXP o2_minSEXP, SEXP h_O2SEXP, SEXP K_downSEXP, SEXP A_angSEXP, SEXP m_onSEXP, SEXP m_offSEXP, SEXP s_onSEXP, SEXP s_offSEXP, SEXP tau_O2SEXP, SEXP o2_logN_epsSEXP, SEXP o2_cache_bin_pctSEXP, SEXP o2_cache_hysteresis_pctSEXP, SEXP lam_minSEXP, SEXP lam_maxSEXP, SEXP k_oSEXP, SEXP has_p_missegSEXP, SEXP p_missegSEXP, SEXP k_o_misSEXP, SEXP has_pmis_endpointsSEXP, SEXP pmis_O2_0SEXP, SEXP pmis_O2_1SEXP, SEXP p_constSEXP, SEXP p_wgdSEXP, SEXP boundarySEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP, SEXP growth_oxy_parSEXP, SEXP vol_by_NSEXP, SEXP burden_floorSEXP, SEXP burden_log_epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type cohort_code(cohort_codeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dose_vec(dose_vecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type treat_day_vec(treat_day_vecSEXP);
    Rcpp::traits::input_parameter< List >::type obs_steps_list(obs_steps_listSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sim_end_step_vec(sim_end_step_vecSEXP);
    Rcpp::traits::input_parameter< List >::type obs_burden_list(obs_burden_listSEXP);
    Rcpp::traits::input_parameter< List >::type keep_burden_list(keep_burden_listSEXP);
    Rcpp::traits::input_parameter< List >::type ploidy_z_list(ploidy_z_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu_by_N(mu_by_NSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_burden(sigma_burdenSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_ploidy(sigma_ploidySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init_state_2N(init_state_2NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init_state_4N(init_state_4NSEXP);
    Rcpp::traits::input_parameter< int >::type N0min(N0minSEXP);
    Rcpp::traits::input_parameter< int >::type N0max(N0maxSEXP);
    Rcpp::traits::input_parameter< int >::type N1min(N1minSEXP);
    Rcpp::traits::input_parameter< int >::type N1max(N1maxSEXP);
    Rcpp::traits::input_parameter< double >::type DT(DTSEXP);
    Rcpp::traits::input_parameter< double >::type dose_ref(dose_refSEXP);
    Rcpp::traits::input_parameter< bool >::type fit_treatment(fit_treatmentSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type tx_mult_min(tx_mult_minSEXP);
    Rcpp::traits::input_parameter< std::string >::type crowding(crowdingSEXP);
    Rcpp::traits::input_parameter< double >::type K(KSEXP);
    Rcpp::traits::input_parameter< double >::type min_pop(min_popSEXP);
    Rcpp::traits::input_parameter< double >::type O2_base(O2_baseSEXP);
    Rcpp::traits::input_parameter< double >::type o2_cap(o2_capSEXP);
    Rcpp::traits::input_parameter< bool >::type o2_feedback(o2_feedbackSEXP);
    Rcpp::traits::input_parameter< double >::type o2_min(o2_minSEXP);
    Rcpp::traits::input_parameter< double >::type h_O2(h_O2SEXP);
    Rcpp::traits::input_parameter< double >::type K_down(K_downSEXP);
    Rcpp::traits::input_parameter< double >::type A_ang(A_angSEXP);
    Rcpp::traits::input_parameter< double >::type m_on(m_onSEXP);
    Rcpp::traits::input_parameter< double >::type m_off(m_offSEXP);
    Rcpp::traits::input_parameter< double >::type s_on(s_onSEXP);
    Rcpp::traits::input_parameter< double >::type s_off(s_offSEXP);
    Rcpp::traits::input_parameter< double >::type tau_O2(tau_O2SEXP);
    Rcpp::traits::input_parameter< double >::type o2_logN_eps(o2_logN_epsSEXP);
    Rcpp::traits::input_parameter< double >::type o2_cache_bin_pct(o2_cache_bin_pctSEXP);
    Rcpp::traits::input_parameter< double >::type o2_cache_hysteresis_pct(o2_cache_hysteresis_pctSEXP);
    Rcpp::traits::input_parameter< double >::type lam_min(lam_minSEXP);
    Rcpp::traits::input_parameter< double >::type lam_max(lam_maxSEXP);
    Rcpp::traits::input_parameter< double >::type k_o(k_oSEXP);
    Rcpp::traits::input_parameter< bool >::type has_p_misseg(has_p_missegSEXP);
    Rcpp::traits::input_parameter< double >::type p_misseg(p_missegSEXP);
    Rcpp::traits::input_parameter< double >::type k_o_mis(k_o_misSEXP);
    Rcpp::traits::input_parameter< bool >::type has_pmis_endpoints(has_pmis_endpointsSEXP);
    Rcpp::traits::input_parameter< double >::type pmis_O2_0(pmis_O2_0SEXP);
    Rcpp::traits::input_parameter< double >::type pmis_O2_1(pmis_O2_1SEXP);
    Rcpp::traits::input_parameter< double >::type p_const(p_constSEXP);
    Rcpp::traits::input_parameter< double >::type p_wgd(p_wgdSEXP);
    Rcpp::traits::input_parameter< std::string >::type boundary(boundarySEXP);
    Rcpp::traits::input_parameter< double >::type eps_tail(eps_tailSEXP);
    Rcpp::traits::input_parameter< double >::type beta_buffer(beta_bufferSEXP);
    Rcpp::traits::input_parameter< double >::type n_exp(n_expSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< int >::type N_unit(N_unitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type growth_oxy_par(growth_oxy_parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vol_by_N(vol_by_NSEXP);
    Rcpp::traits::input_parameter< double >::type burden_floor(burden_floorSEXP);
    Rcpp::traits::input_parameter< double >::type burden_log_eps(burden_log_epsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_o2invivo_objective_components_map(cohort_code, dose_vec, treat_day_vec, obs_steps_list, sim_end_step_vec, obs_burden_list, keep_burden_list, ploidy_z_list, mu_by_N, sigma_burden, sigma_ploidy, init_state_2N, init_state_4N, N0min, N0max, N1min, N1max, DT, dose_ref, fit_treatment, alpha, gamma, tx_mult_min, crowding, K, min_pop, O2_base, o2_cap, o2_feedback, o2_min, h_O2, K_down, A_ang, m_on, m_off, s_on, s_off, tau_O2, o2_logN_eps, o2_cache_bin_pct, o2_cache_hysteresis_pct, lam_min, lam_max, k_o, has_p_misseg, p_misseg, k_o_mis, has_pmis_endpoints, pmis_O2_0, pmis_O2_1, p_const, p_wgd, boundary, eps_tail, beta_buffer, n_exp, smax, N_unit, growth_oxy_par, vol_by_N, burden_floor, burden_log_eps));
    return rcpp_result_gen;
END_RCPP
}
