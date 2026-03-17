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
#include <utility>
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
// Function: representative_chromosome_vector
// Purpose: Build deterministic temporary chromosome-vector representative for a
//   coarse total-chromosome state N. Entries are as evenly distributed as
//   possible and sum exactly to N.
// Parameters:
//   - N: Coarse total-chromosome state.
//   - n_chr: Number of modeled chromosome classes (22 autosomes).
// Returns:
//   std::vector<int> return value containing the representative copy vector.
// -----------------------------------------------------------------------------
inline std::vector<int> representative_chromosome_vector(int N, int n_chr) {
  const int n_chr_use = std::max(1, n_chr);
  const int N_use = std::max(0, N);
  const int base = N_use / n_chr_use;
  const int rem = N_use % n_chr_use;
  std::vector<int> out(static_cast<size_t>(n_chr_use), base);
  for (int i = 0; i < rem; ++i) out[static_cast<size_t>(i)] += 1;
  return out;
}

// -----------------------------------------------------------------------------
// Function: select_affected_chromosomes
// Purpose: Deterministic expected-kernel approximation for choosing affected
//   chromosomes. Lower-copy chromosomes are selected first to preserve strong
//   loss-side lethality in the coarse projection.
// Parameters:
//   - n_rep: Representative chromosome vector for the mother state.
//   - k: Number of affected chromosomes requested.
// Returns:
//   std::vector<int> return value containing selected chromosome indices.
// -----------------------------------------------------------------------------
inline std::vector<int> select_affected_chromosomes(
    const std::vector<int>& n_rep,
    int k
) {
  const int n_chr = static_cast<int>(n_rep.size());
  std::vector<int> idx;
  idx.reserve(static_cast<size_t>(n_chr));
  for (int c = 0; c < n_chr; ++c) {
    if (n_rep[static_cast<size_t>(c)] > 0) idx.push_back(c);
  }
  std::stable_sort(
    idx.begin(),
    idx.end(),
    [&](int a, int b) {
      const int na = n_rep[static_cast<size_t>(a)];
      const int nb = n_rep[static_cast<size_t>(b)];
      if (na != nb) return na < nb;
      return a < b;
    }
  );
  if (k < static_cast<int>(idx.size())) idx.resize(static_cast<size_t>(std::max(0, k)));
  return idx;
}

// -----------------------------------------------------------------------------
// Function: daughter_viability
// Purpose: Evaluate daughter-specific viability with hard null-chromosome
//   lethality and asymmetric loss-vs-gain penalties relative to diploid copy=2.
// Parameters:
//   - n_vec: Daughter chromosome copy vector.
//   - beta_loss: Loss-burden penalty coefficient (dominant penalty).
//   - beta_gain: Gain-burden penalty coefficient.
// Returns:
//   double return value containing viability in [0, 1].
// -----------------------------------------------------------------------------
inline double daughter_viability(
    const std::vector<int>& n_vec,
    double beta_loss,
    double beta_gain
) {
  int z0 = 0;
  double L = 0.0;
  double G = 0.0;
  for (size_t i = 0; i < n_vec.size(); ++i) {
    const int nc = n_vec[i];
    if (nc == 0) ++z0;
    const int loss = std::max(0, 2 - nc);
    const int gain = std::max(0, nc - 2);
    L += static_cast<double>(loss * loss);
    G += static_cast<double>(gain);
  }
  if (z0 > 0) return 0.0;
  const double beta_loss_use = (std::isfinite(beta_loss) && beta_loss >= 0.0) ? beta_loss : 0.0;
  const double beta_gain_use = (std::isfinite(beta_gain) && beta_gain >= 0.0) ? beta_gain : 0.0;
  const double expo = -beta_loss_use * L - beta_gain_use * G;
  const double expo_clamped = std::max(-700.0, std::min(0.0, expo));
  return std::exp(expo_clamped);
}

// -----------------------------------------------------------------------------
// Function: o2simps_pr_delta_internal
// Purpose: Build coarse ploidy transition weights using a temporary chromosome-
//   vector representation with asymmetric daughter-specific post-missegregation
//   viability and explicit nonviable-daughter mass accounting.
// Parameters:
//   - N: Ploidy state value or chromosome-copy count.
//   - p: Missegregation probability parameter.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - beta_buffer: Buffer exponent controlling ploidy-dependence of missegregation survival.
//   - n_exp: Exponent controlling ploidy scaling in buffering term.
//   - smax: Maximum survival factor for missegregation events.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
//   - beta_loss: Loss-burden penalty coefficient.
//   - beta_gain: Gain-burden penalty coefficient.
//   - ts_out: Function-specific input argument.
//   - prob_out: Function-specific input argument.
//   - mass_dropped: Function-specific input argument.
// Returns:
//   void return value containing the computed result.
// -----------------------------------------------------------------------------
void o2simps_pr_delta_internal(
    int N,
    double p,
    double eps_tail,
    double beta_buffer,
    double n_exp,
    double smax,
    int N_unit,
    double beta_loss,
    double beta_gain,
    std::vector<int>& ts_out,
    std::vector<double>& prob_out,
    double& mass_dropped
) {
  const int n_chr = std::max(1, N_unit);
  const std::vector<int> n_rep = representative_chromosome_vector(N, n_chr);
  if (p <= 0.0 || N <= 0) {
    ts_out.assign(1, 0);
    prob_out.assign(1, 2.0);
    mass_dropped = 0.0;
    return;
  }

  const double p_use = clamp01(p);
  const double eps_use = (std::isfinite(eps_tail) && eps_tail > 0.0) ? eps_tail : 0.0;
  const double n_d = static_cast<double>(std::max(1, N));
  const double n_unit_d = static_cast<double>(std::max(1, N_unit));
  double sN = smax * std::exp(-beta_buffer * std::pow((2.0 * n_unit_d) / n_d, n_exp));
  if (!std::isfinite(sN)) sN = 1.0;
  sN = std::max(0.0, std::min(1.0, sN));

  const int shift_offset = n_chr;
  std::vector<double> shift_mass(static_cast<size_t>(2 * n_chr + 1), 0.0);
  double survivors_total = 0.0;

  for (int ks = 0; ks <= N; ++ks) {
    const double pk = R::dbinom(ks, N, p_use, false);
    if (!std::isfinite(pk) || pk <= 0.0) continue;
    if (eps_use > 0.0 && pk < eps_use) continue;

    if (ks == 0) {
      shift_mass[static_cast<size_t>(shift_offset)] += 2.0 * pk;
      survivors_total += 2.0 * pk;
      continue;
    }

    const std::vector<int> affected = select_affected_chromosomes(n_rep, ks);
    const int k_eff = static_cast<int>(affected.size());
    if (k_eff <= 0) {
      shift_mass[static_cast<size_t>(shift_offset)] += 2.0 * pk;
      survivors_total += 2.0 * pk;
      continue;
    }

    std::vector<int> daughter_gain = n_rep;
    std::vector<int> daughter_loss = n_rep;
    for (int i = 0; i < k_eff; ++i) {
      const int c = affected[static_cast<size_t>(i)];
      daughter_gain[static_cast<size_t>(c)] += 1;
      daughter_loss[static_cast<size_t>(c)] = std::max(0, daughter_loss[static_cast<size_t>(c)] - 1);
    }

    const double V_gain = daughter_viability(daughter_gain, beta_loss, beta_gain);
    const double V_loss = daughter_viability(daughter_loss, beta_loss, beta_gain);
    const double bg = std::pow(sN, static_cast<double>(k_eff));
    const double bg_use = (std::isfinite(bg) && bg >= 0.0) ? std::min(bg, 1.0) : 0.0;
    const double w_gain = pk * bg_use * V_gain;
    const double w_loss = pk * bg_use * V_loss;

    const int shift_gain = std::min(n_chr, k_eff);
    const int shift_loss = -std::min(n_chr, k_eff);
    shift_mass[static_cast<size_t>(shift_gain + shift_offset)] += w_gain;
    shift_mass[static_cast<size_t>(shift_loss + shift_offset)] += w_loss;
    survivors_total += (w_gain + w_loss);
  }

  if (!std::isfinite(survivors_total) || survivors_total < 0.0) survivors_total = 0.0;
  if (survivors_total > 2.0) {
    const double scale = 2.0 / survivors_total;
    for (size_t i = 0; i < shift_mass.size(); ++i) shift_mass[i] *= scale;
    survivors_total = 2.0;
  }

  const double dead_daughters = std::max(0.0, 2.0 - survivors_total);
  mass_dropped = std::max(0.0, std::min(1.0, dead_daughters / 2.0));

  ts_out.clear();
  prob_out.clear();
  ts_out.reserve(shift_mass.size());
  prob_out.reserve(shift_mass.size());
  for (int t = -n_chr; t <= n_chr; ++t) {
    const double w = shift_mass[static_cast<size_t>(t + shift_offset)];
    if (!std::isfinite(w) || w <= 0.0) continue;
    ts_out.push_back(t);
    prob_out.push_back(w);
  }

  if (ts_out.empty()) {
    ts_out.assign(1, 0);
    prob_out.assign(1, 0.0);
  }
}

} // namespace

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_representative_vector
// Purpose: Expose representative chromosome-vector construction used by the
//   asymmetric kernel for validation/debugging.
// Parameters:
//   - N: Coarse total-chromosome state.
//   - N_unit: Number of chromosome classes.
// Returns:
//   IntegerVector return value containing n_rep(N).
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
IntegerVector cpp_o2simps_representative_vector(
    int N,
    int N_unit = 22
) {
  const std::vector<int> rep = representative_chromosome_vector(N, std::max(1, N_unit));
  return IntegerVector(rep.begin(), rep.end());
}

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_asymmetric_daughters_example
// Purpose: Build one deterministic asymmetric daughter pair for validation of
//   chromosome-aware missegregation handling.
// Parameters:
//   - N: Coarse mother state.
//   - k_s: Number of affected chromosomes.
//   - N_unit: Number of chromosome classes.
//   - beta_loss: Loss-burden penalty coefficient.
//   - beta_gain: Gain-burden penalty coefficient.
// Returns:
//   List return value containing mother/daughter vectors and viability metrics.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_asymmetric_daughters_example(
    int N,
    int k_s,
    int N_unit = 22,
    double beta_loss = 0.25,
    double beta_gain = 0.05
) {
  const std::vector<int> mother = representative_chromosome_vector(N, std::max(1, N_unit));
  const std::vector<int> affected = select_affected_chromosomes(mother, std::max(0, k_s));
  std::vector<int> daughter_gain = mother;
  std::vector<int> daughter_loss = mother;
  for (size_t i = 0; i < affected.size(); ++i) {
    const int c = affected[i];
    daughter_gain[static_cast<size_t>(c)] += 1;
    daughter_loss[static_cast<size_t>(c)] = std::max(0, daughter_loss[static_cast<size_t>(c)] - 1);
  }
  int z0_gain = 0;
  int z0_loss = 0;
  for (size_t i = 0; i < daughter_gain.size(); ++i) {
    if (daughter_gain[i] == 0) ++z0_gain;
    if (daughter_loss[i] == 0) ++z0_loss;
  }
  return List::create(
    _["mother"] = IntegerVector(mother.begin(), mother.end()),
    _["daughter_gain"] = IntegerVector(daughter_gain.begin(), daughter_gain.end()),
    _["daughter_loss"] = IntegerVector(daughter_loss.begin(), daughter_loss.end()),
    _["affected_idx_1based"] = IntegerVector(affected.begin(), affected.end()) + 1,
    _["N_mother"] = std::accumulate(mother.begin(), mother.end(), 0),
    _["N_gain"] = std::accumulate(daughter_gain.begin(), daughter_gain.end(), 0),
    _["N_loss"] = std::accumulate(daughter_loss.begin(), daughter_loss.end(), 0),
    _["z0_gain"] = z0_gain,
    _["z0_loss"] = z0_loss,
    _["V_gain"] = daughter_viability(daughter_gain, beta_loss, beta_gain),
    _["V_loss"] = daughter_viability(daughter_loss, beta_loss, beta_gain)
  );
}

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_pr_delta_vec
// Purpose: Compute missegregation delta-kernel probabilities over ploidy shifts.
// Parameters:
//   - N: Ploidy state value or chromosome-copy count.
//   - p: Missegregation probability parameter.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - beta_buffer: Buffer exponent controlling ploidy-dependence of missegregation survival.
//   - n_exp: Exponent controlling ploidy scaling in buffering term.
//   - smax: Maximum survival factor for missegregation events.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
//   - beta_loss: Loss-burden penalty coefficient.
//   - beta_gain: Gain-burden penalty coefficient.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_pr_delta_vec(
    int N,
    double p,
    double eps_tail = 1e-8,
    double beta_buffer = 0.0,
    double n_exp = 1.0,
    double smax = 1.0,
    int N_unit = 22,
    double beta_loss = 0.25,
    double beta_gain = 0.05
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
    beta_loss,
    beta_gain,
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
// Function: cpp_o2simps_o2_window_supply
// Purpose: Compute oxygen supply fraction/level from burden under the selected O2 model.
// Parameters:
//   - Ntot: Total predicted cell count (or burden proxy) at current time.
//   - curve_type: Function-specific input argument.
//   - O2_cap: Function-specific input argument.
//   - o2_init: Function-specific input argument.
//   - o2_rate: Function-specific input argument.
//   - o2_shape_v: Function-specific input argument.
//   - o2_anchor_N: Function-specific input argument.
//   - o2_logN_eps: Function-specific input argument.
// Returns:
//   NumericVector return value containing the computed result.
// -----------------------------------------------------------------------------
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
  const double O2_min_use = 0.0;
  const double O2_span_use = std::max(0.0, O2_cap_use - O2_min_use);
  const double eps_o2 = 1e-9;
  const double o2_init_use = std::max(O2_min_use + eps_o2, std::min(O2_cap_use - eps_o2, o2_init));
  const double rate_use = (std::isfinite(o2_rate) && o2_rate > 0.0) ? o2_rate : 1.0;
  const double v_use = (std::isfinite(o2_shape_v) && o2_shape_v > 0.0) ? o2_shape_v : 1.0;
  const double anchor_use = (std::isfinite(o2_anchor_N) && o2_anchor_N >= 0.0) ? o2_anchor_N : 1e6;
  const double eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;
  const double log_anchor = std::log10(anchor_use + eps_use);
  const double y_anchor = std::max(
    eps_o2,
    std::min(1.0 - eps_o2, (o2_init_use - O2_min_use) / std::max(O2_span_use, eps_o2))
  );
  const double bg = std::log(std::max(eps_o2, -std::log(y_anchor)));
  const double bl = std::log(std::max(eps_o2, std::pow(1.0 / y_anchor, v_use) - 1.0));
  for (int i = 0; i < n; ++i) {
    const double n_raw = Ntot[i];
    const double n_use = (std::isfinite(n_raw) && n_raw > 0.0) ? n_raw : 0.0;
    const double xi_size = (std::log10(n_use + eps_use) - log_anchor) / 3.0;
    double o2 = O2_min_use;
    if (use_glogistic) {
      const double arg = rate_use * xi_size + bl;
      const double exp_arg = (arg > 700.0) ? std::exp(700.0) : ((arg < -700.0) ? std::exp(-700.0) : std::exp(arg));
      const double denom = std::pow(1.0 + exp_arg, 1.0 / v_use);
      o2 = O2_min_use + O2_span_use / std::max(denom, eps_o2);
    } else {
      const double arg = rate_use * xi_size + bg;
      const double exp_arg = (arg > 700.0) ? std::exp(700.0) : ((arg < -700.0) ? std::exp(-700.0) : std::exp(arg));
      o2 = O2_min_use + O2_span_use * std::exp(-exp_arg);
    }
    out[i] = clamp_o2_pct(o2);
  }

  return out;
}

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_build_B_total_triplet
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
//   - beta_loss: Loss-burden penalty coefficient.
//   - beta_gain: Gain-burden penalty coefficient.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
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
    int N_unit = 22,
    double beta_loss = 0.25,
    double beta_gain = 0.05
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
      beta_loss,
      beta_gain,
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
      // Asymmetric kernel returns daughter-projected signed shifts directly.
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
// Function: cpp_o2simps_build_B_WGD_triplet
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
// Function: cpp_o2simps_build_G_for_o2_triplet
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
//   - beta_loss: Loss-burden penalty coefficient.
//   - beta_gain: Gain-burden penalty coefficient.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_build_G_for_o2_triplet(
    double O2,
    double O2_cap,
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
    double beta_loss = 0.25,
    double beta_gain = 0.05,
    double beta_size = 0.0,
    double alpha_o2 = 0.0,
    double gamma_growth = 1.0,
    bool growth_penalty_ploidy = false,
    bool growth_penalty_hypoxia = false,
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
  const double gamma_growth_use = (std::isfinite(gamma_growth) && gamma_growth > 0.0) ? gamma_growth : 1.0;
  const double o2_cap_use = clamp_o2_pct(O2_cap);
  const double oxygen_norm_eps = 1e-12;
  const double h_o2 = std::max(0.0, 1.0 - O2_use / (o2_cap_use + oxygen_norm_eps));
  const double mu_hp_use = (std::isfinite(mu_hp) && mu_hp > 0.0) ? mu_hp : 0.0;
  auto lam_for_N = [&](int N_state) -> double {
    if (lam_base <= 0.0) return 0.0;
    const double N_state_use = static_cast<double>(N_state);
    const double size_penalty =
      growth_penalty_ploidy ? std::exp(-beta_size_use * std::pow(N_state_use, gamma_growth_use)) : 1.0;
    // Unified oxygen weighting: continuous and zero at oxygen cap.
    const double hypoxia_penalty = growth_penalty_hypoxia ?
      std::exp(-alpha_o2_use * std::pow(N_state_use, gamma_growth_use) * h_o2) : 1.0;
    double lam_eff = lam_base * size_penalty * hypoxia_penalty;
    if (!std::isfinite(lam_eff) || lam_eff < 0.0) lam_eff = 0.0;
    return lam_eff;
  };
  auto mu_for_N = [&](int N_state) -> double {
    if (mu_hp_use <= 0.0) return 0.0;
    if (h_o2 <= 0.0) return 0.0;
    const double N_state_use = static_cast<double>(N_state);
    // This channel represents hypoxia-linked continuous death only.
    // Nonviable post-missegregation daughters are routed separately via dead_buffer_rate.
    double mu_eff = mu_hp_use * N_state_use * h_o2;
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
  std::vector<double> dead_buffer_rate(static_cast<size_t>(R0 + R1), 0.0);
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
    o2simps_pr_delta_internal(
      N,
      p_mis,
      eps_tail,
      beta_buffer,
      n_exp,
      smax,
      N_unit,
      beta_loss,
      beta_gain,
      ts,
      pr,
      mass_dropped
    );

    const double scale_pre = lam_N * (1.0 - pw);
    {
      // Nonviable daughters from asymmetric post-missegregation fate are
      // routed into the dead compartment (not dropped silently).
      double nonviable_rate = 2.0 * scale_pre * mass_dropped;
      if (!std::isfinite(nonviable_rate) || nonviable_rate < 0.0) nonviable_rate = 0.0;
      dead_buffer_rate[static_cast<size_t>(col_1based - 1)] = nonviable_rate;
    }
    for (size_t k = 0; k < ts.size(); ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;
      // Asymmetric daughter projection already encoded as signed shift t.
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
    o2simps_pr_delta_internal(
      N,
      p_mis,
      eps_tail,
      beta_buffer,
      n_exp,
      smax,
      N_unit,
      beta_loss,
      beta_gain,
      ts,
      pr,
      mass_dropped
    );
    {
      // Nonviable daughters from asymmetric post-missegregation fate are
      // routed into the dead compartment (not dropped silently).
      double nonviable_rate = 2.0 * lam_N * mass_dropped;
      if (!std::isfinite(nonviable_rate) || nonviable_rate < 0.0) nonviable_rate = 0.0;
      dead_buffer_rate[static_cast<size_t>(col_1based - 1)] = nonviable_rate;
    }

    for (size_t k = 0; k < ts.size(); ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;
      // Asymmetric daughter projection already encoded as signed shift t.
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
    _["ncol"] = R0 + R1,
    _["dead_buffer_rate"] = NumericVector(dead_buffer_rate.begin(), dead_buffer_rate.end())
  );
}

namespace {

using SpMat = Eigen::SparseMatrix<double, Eigen::RowMajor, int>;

struct SparseCacheEntry {
  SpMat mat;
  std::vector<double> dead_buffer_rate;
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
//   - beta_loss: Loss-burden penalty coefficient.
//   - beta_gain: Gain-burden penalty coefficient.
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
    double beta_loss,
    double beta_gain,
    double beta_size,
    double alpha_o2,
    double O2_cap,
    double gamma_growth,
    bool growth_penalty_ploidy,
    bool growth_penalty_hypoxia,
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
  hash_combine_cpp(seed, bits_of_double_cpp(beta_loss));
  hash_combine_cpp(seed, bits_of_double_cpp(beta_gain));
  hash_combine_cpp(seed, bits_of_double_cpp(beta_size));
  hash_combine_cpp(seed, bits_of_double_cpp(alpha_o2));
  hash_combine_cpp(seed, bits_of_double_cpp(O2_cap));
  hash_combine_cpp(seed, bits_of_double_cpp(gamma_growth));
  hash_combine_cpp(seed, growth_penalty_ploidy ? 1 : 0);
  hash_combine_cpp(seed, growth_penalty_hypoxia ? 1 : 0);
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
//   - curve_type: Function-specific input argument.
//   - O2_cap: Function-specific input argument.
//   - o2_init: Function-specific input argument.
//   - o2_rate: Function-specific input argument.
//   - o2_shape_v: Function-specific input argument.
//   - o2_anchor_N: Function-specific input argument.
//   - o2_logN_eps: Function-specific input argument.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
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
  const double O2_min_use = 0.0;
  const double O2_span_use = std::max(0.0, O2_cap_use - O2_min_use);
  const double eps_o2 = 1e-9;
  const double o2_init_use = std::max(O2_min_use + eps_o2, std::min(O2_cap_use - eps_o2, o2_init));
  const double rate_use = (std::isfinite(o2_rate) && o2_rate > 0.0) ? o2_rate : 1.0;
  const double v_use = (std::isfinite(o2_shape_v) && o2_shape_v > 0.0) ? o2_shape_v : 1.0;
  const double anchor_use = (std::isfinite(o2_anchor_N) && o2_anchor_N >= 0.0) ? o2_anchor_N : 1e6;
  const double eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;
  const double log_anchor = std::log10(anchor_use + eps_use);
  const double xi_size = (std::log10(n_use + eps_use) - log_anchor) / 3.0;
  const double y_anchor = std::max(
    eps_o2,
    std::min(1.0 - eps_o2, (o2_init_use - O2_min_use) / std::max(O2_span_use, eps_o2))
  );
  if (use_glogistic) {
    const double bl = std::log(std::max(eps_o2, std::pow(1.0 / y_anchor, v_use) - 1.0));
    const double arg = rate_use * xi_size + bl;
    const double exp_arg = (arg > 700.0) ? std::exp(700.0) : ((arg < -700.0) ? std::exp(-700.0) : std::exp(arg));
    const double denom = std::pow(1.0 + exp_arg, 1.0 / v_use);
    const double o2 = O2_min_use + O2_span_use / std::max(denom, eps_o2);
    return clamp_o2_pct(o2);
  }
  const double bg = std::log(std::max(eps_o2, -std::log(y_anchor)));
  const double arg = rate_use * xi_size + bg;
  const double exp_arg = (arg > 700.0) ? std::exp(700.0) : ((arg < -700.0) ? std::exp(-700.0) : std::exp(arg));
  const double o2 = O2_min_use + O2_span_use * std::exp(-exp_arg);
  return clamp_o2_pct(o2);
}

// -----------------------------------------------------------------------------
// Function: death_rate_for_N_cpp
// Purpose: Compute live->dead transfer death rate using unified oxygen weighting.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - N_unit: Ploidy scaling unit used to map integer states to N values.
//   - O2_use: Oxygen level used by model rate functions.
//   - O2_cap_use: Oxygen cap used to normalize continuous hypoxia weighting.
//   - mu_hp_use: Hypoxia-linked high-ploidy death strength.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double death_rate_for_N_cpp(
    int N_state,
    int N_unit,
    double O2_use,
    double O2_cap_use,
    double mu_hp_use
) {
  if (!(std::isfinite(mu_hp_use) && mu_hp_use > 0.0)) return 0.0;
  (void)N_unit;
  const double o2_cap_clamped = clamp_o2_pct(O2_cap_use);
  const double oxygen_norm_eps = 1e-12;
  const double h_o2 = std::max(0.0, 1.0 - clamp_o2_pct(O2_use) / (o2_cap_clamped + oxygen_norm_eps));
  if (h_o2 <= 0.0) return 0.0;
  const double N_state_use = static_cast<double>(N_state);
  double mu_eff = mu_hp_use * N_state_use * h_o2;
  if (!std::isfinite(mu_eff) || mu_eff < 0.0) mu_eff = 0.0;
  return mu_eff;
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
  out.dead_buffer_rate.assign(static_cast<size_t>(ncol), 0.0);
  if (tri.containsElementNamed("dead_buffer_rate")) {
    NumericVector db = tri["dead_buffer_rate"];
    if (db.size() != ncol) {
      stop("dead_buffer_rate length mismatch.");
    }
    for (int i = 0; i < ncol; ++i) {
      double v = db[i];
      if (!std::isfinite(v) || v < 0.0) v = 0.0;
      out.dead_buffer_rate[static_cast<size_t>(i)] = v;
    }
  }
  return out;
}

} // namespace

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_simulate_one
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
//   - O2_cap: Function-specific input argument.
//   - o2_feedback: Function-specific input argument.
//   - o2_curve_type: Function-specific input argument.
//   - o2_init: Function-specific input argument.
//   - o2_rate: Function-specific input argument.
//   - o2_shape_v: Function-specific input argument.
//   - tau_O2: Relaxation time constant controlling lag from O2 target to O2 effective.
//   - o2_anchor_N: Function-specific input argument.
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
//   - beta_loss: Loss-burden penalty coefficient.
//   - beta_gain: Gain-burden penalty coefficient.
//   - vol_by_N: Optional precomputed per-state cell volume lookup.
//   - burden_floor: Function-specific input argument.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
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
    double tau_O2,
    double o2_anchor_N,
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
    double beta_loss,
    double beta_gain,
    double beta_size,
    double alpha_o2,
    double gamma_growth,
    bool growth_penalty_ploidy,
    bool growth_penalty_hypoxia,
    double mu_hp,
    double k_clear,
    NumericVector vol_by_N,
    double burden_floor
) {
  const int R0 = N0max - N0min + 1;
  const int R1 = N1max - N1min + 1;
  if (R0 <= 0 || R1 <= 0) stop("Nmax must be >= Nmin for both layers.");
  const int D = R0 + R1;
  if (!(init_state.size() == D || init_state.size() == 2 * D)) stop("init_state length mismatch.");
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

  std::vector<double> Ntot_live_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Ntot_dead_hypoxia_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Ntot_dead_buffer_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Ntot_dead_total_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Ntot_total_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_live_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_dead_hypoxia_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_dead_buffer_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_dead_total_at_step(step_unique.size(), NA_REAL);
  std::vector<double> Vmm3_total_at_step(step_unique.size(), NA_REAL);

  std::vector<double> v_live(static_cast<size_t>(D), 0.0);
  std::vector<double> v_dead_hypoxia(static_cast<size_t>(D), 0.0);
  std::vector<double> v_dead_buffer(static_cast<size_t>(D), 0.0);
  if (init_state.size() == D) {
    std::copy(init_state.begin(), init_state.end(), v_live.begin());
  } else {
    for (int i = 0; i < D; ++i) {
      v_live[static_cast<size_t>(i)] = init_state[i];
      // Preserve backward compatibility for 2*D init_state by placing legacy
      // dead initialization into hypoxia-origin dead pool.
      v_dead_hypoxia[static_cast<size_t>(i)] = init_state[D + i];
    }
  }
  std::vector<double> growth(static_cast<size_t>(D), 0.0);
  std::vector<double> death_flow_hypoxia(static_cast<size_t>(D), 0.0);
  std::vector<double> death_flow_buffer(static_cast<size_t>(D), 0.0);

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
    beta_loss,
    beta_gain,
    beta_size,
    alpha_o2,
    O2_cap,
    gamma_growth,
    growth_penalty_ploidy,
    growth_penalty_hypoxia,
    mu_hp,
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
  const double tau_use = (std::isfinite(tau_O2) && tau_O2 > 0.0) ? tau_O2 : 2.0;
  const double alpha_tau = 1.0 - std::exp(-DT_use / tau_use);
  const double o2_anchor_use = (std::isfinite(o2_anchor_N) && o2_anchor_N >= 0.0) ? o2_anchor_N : 1e6;
  const double o2_eps_use = (std::isfinite(o2_logN_eps) && o2_logN_eps > 0.0) ? o2_logN_eps : 1.0;
  const double o2_bin_use = (std::isfinite(o2_cache_bin_pct) && o2_cache_bin_pct > 0.0) ? o2_cache_bin_pct : 1e-3;
  const double o2_hyst_use = (std::isfinite(o2_cache_hysteresis_pct) && o2_cache_hysteresis_pct >= 0.0) ? o2_cache_hysteresis_pct : 0.0;
  const double k_clear_use = (std::isfinite(k_clear) && k_clear >= 0.0) ? k_clear : 0.0;
  (void) o2_cache_profile;
  int cache_g_build = 0;
  int cache_g_hit = 0;
  int cache_g_hysteresis = 0;
  bool has_last_key = false;
  int last_key = 0;
  double last_o2_eff = 0.0;
  double O2_state = O2_cap_use;
  if (o2_feedback) {
    O2_state = o2_window_supply_scalar_cpp(
      vector_sum_cpp(v_live),
      o2_curve_type,
      O2_cap_use,
      o2_init_use,
      o2_rate_use,
      o2_shape_v_use,
      o2_anchor_use,
      o2_eps_use
    );
    O2_state = clamp_o2_pct(O2_state);
  }

  for (int step = 0; step <= final_step; ++step) {
    auto it_obs = step_to_idx.find(step);
    if (it_obs != step_to_idx.end()) {
      const int idx = it_obs->second;
      const double Ntot_live_now = vector_sum_cpp(v_live);
      const double Ntot_dead_h_now = vector_sum_cpp(v_dead_hypoxia);
      const double Ntot_dead_b_now = vector_sum_cpp(v_dead_buffer);
      const double Ntot_dead_now = Ntot_dead_h_now + Ntot_dead_b_now;
      const double Ntot_total_now = Ntot_live_now + Ntot_dead_now;
      Ntot_live_at_step[static_cast<size_t>(idx)] = Ntot_live_now;
      Ntot_dead_hypoxia_at_step[static_cast<size_t>(idx)] = Ntot_dead_h_now;
      Ntot_dead_buffer_at_step[static_cast<size_t>(idx)] = Ntot_dead_b_now;
      Ntot_dead_total_at_step[static_cast<size_t>(idx)] = Ntot_dead_now;
      Ntot_total_at_step[static_cast<size_t>(idx)] = Ntot_total_now;
      double burden_live_now = 0.0;
      double burden_dead_h_now = 0.0;
      double burden_dead_b_now = 0.0;
      double burden_dead_now = 0.0;
      double burden_total_now = 0.0;
      for (int i = 0; i < R0; ++i) {
        const size_t pre_idx = static_cast<size_t>(i);
        const size_t post_idx = static_cast<size_t>(R0 + i);
        const double n_live_i = v_live[pre_idx] + v_live[post_idx];
        const double n_dead_h_i = v_dead_hypoxia[pre_idx] + v_dead_hypoxia[post_idx];
        const double n_dead_b_i = v_dead_buffer[pre_idx] + v_dead_buffer[post_idx];
        const double n_dead_i = n_dead_h_i + n_dead_b_i;
        const double n_total_i = n_live_i + n_dead_i;
        burden_live_now += n_live_i * vol_by_N[i];
        burden_dead_h_now += n_dead_h_i * vol_by_N[i];
        burden_dead_b_now += n_dead_b_i * vol_by_N[i];
        burden_dead_now += n_dead_i * vol_by_N[i];
        burden_total_now += n_total_i * vol_by_N[i];
      }
      Vmm3_live_at_step[static_cast<size_t>(idx)] = burden_live_now;
      Vmm3_dead_hypoxia_at_step[static_cast<size_t>(idx)] = burden_dead_h_now;
      Vmm3_dead_buffer_at_step[static_cast<size_t>(idx)] = burden_dead_b_now;
      Vmm3_dead_total_at_step[static_cast<size_t>(idx)] = burden_dead_now;
      Vmm3_total_at_step[static_cast<size_t>(idx)] = burden_total_now;
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

    const double Ntot_live = vector_sum_cpp(v_live);
    double O2_target = O2_cap_use;
    if (o2_feedback) {
      O2_target = o2_window_supply_scalar_cpp(
        Ntot_live,
        o2_curve_type,
        O2_cap_use,
        o2_init_use,
        o2_rate_use,
        o2_shape_v_use,
        o2_anchor_use,
        o2_eps_use
      );
    }
    O2_target = clamp_o2_pct(O2_target);
    O2_state = O2_state + alpha_tau * (O2_target - O2_state);
    double O2_eff = clamp_o2_pct(O2_state);

    int gkey = quantize_o2_key(O2_eff, o2_bin_use);
    if (o2_hyst_use > 0.0 && has_last_key && std::abs(O2_eff - last_o2_eff) <= o2_hyst_use) {
      gkey = last_key;
      ++cache_g_hysteresis;
    }
    auto itG = shared_G_cache.find(gkey);
    if (itG == shared_G_cache.end()) {
      const List tri = cpp_o2simps_build_G_for_o2_triplet(
        O2_eff,
        O2_cap,
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
        beta_loss,
        beta_gain,
        beta_size,
        alpha_o2,
        gamma_growth,
        growth_penalty_ploidy,
        growth_penalty_hypoxia,
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

    sparse_mv_cpp(itG->second, v_live, growth);
    const double crowd = crowd_logistic ? std::max(0.0, 1.0 - Ntot_live / K_use) : std::exp(-Ntot_live / K_use);
    const double scalar = DT_use * crowd * tx_mult;
    for (int i = 0; i < D; ++i) {
      const int N_state = (i < R0) ? (N0min + i) : (N1min + (i - R0));
      const double mu_i = death_rate_for_N_cpp(N_state, N_unit, O2_eff, O2_cap, mu_hp);
      const double src_live = v_live[static_cast<size_t>(i)];
      double flow_h_i = scalar * mu_i * src_live;
      if (!std::isfinite(flow_h_i) || flow_h_i < 0.0) flow_h_i = 0.0;
      death_flow_hypoxia[static_cast<size_t>(i)] = flow_h_i;
      double db_rate_i = 0.0;
      if (static_cast<size_t>(i) < itG->second.dead_buffer_rate.size()) {
        db_rate_i = itG->second.dead_buffer_rate[static_cast<size_t>(i)];
      }
      if (!std::isfinite(db_rate_i) || db_rate_i < 0.0) db_rate_i = 0.0;
      // Buffer-derived mitotic catastrophe inflow (not a continuous death hazard).
      double flow_b_i = scalar * db_rate_i * src_live;
      if (!std::isfinite(flow_b_i) || flow_b_i < 0.0) flow_b_i = 0.0;
      death_flow_buffer[static_cast<size_t>(i)] = flow_b_i;
    }
    for (size_t i = 0; i < v_live.size(); ++i) {
      const double next_v = v_live[i] + scalar * growth[i];
      if (!std::isfinite(next_v) || next_v < 0.0) {
        v_live[i] = 0.0;
      } else {
        v_live[i] = next_v;
      }
    }
    for (size_t i = 0; i < v_dead_hypoxia.size(); ++i) {
      const double dead_h_prev = v_dead_hypoxia[i];
      const double dead_h_next = dead_h_prev + death_flow_hypoxia[i] - DT_use * k_clear_use * dead_h_prev;
      if (!std::isfinite(dead_h_next) || dead_h_next < 0.0) {
        v_dead_hypoxia[i] = 0.0;
      } else {
        v_dead_hypoxia[i] = dead_h_next;
      }
      const double dead_b_prev = v_dead_buffer[i];
      const double dead_b_next = dead_b_prev + death_flow_buffer[i] - DT_use * k_clear_use * dead_b_prev;
      if (!std::isfinite(dead_b_next) || dead_b_next < 0.0) {
        v_dead_buffer[i] = 0.0;
      } else {
        v_dead_buffer[i] = dead_b_next;
      }
    }
    if (vector_sum_cpp(v_live) <= min_pop_use &&
        (vector_sum_cpp(v_dead_hypoxia) + vector_sum_cpp(v_dead_buffer)) <= min_pop_use) break;
  }

  NumericVector Ntot_live_obs(obs_v.size(), NA_REAL);
  NumericVector Ntot_dead_hypoxia_obs(obs_v.size(), NA_REAL);
  NumericVector Ntot_dead_buffer_obs(obs_v.size(), NA_REAL);
  NumericVector Ntot_dead_total_obs(obs_v.size(), NA_REAL);
  NumericVector Ntot_total_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_live_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_dead_hypoxia_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_dead_buffer_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_dead_total_obs(obs_v.size(), NA_REAL);
  NumericVector Vmm3_total_obs(obs_v.size(), NA_REAL);
  for (int i = 0; i < static_cast<int>(obs_v.size()); ++i) {
    auto it = step_to_idx.find(obs_v[static_cast<size_t>(i)]);
    if (it == step_to_idx.end()) {
      Ntot_live_obs[i] = min_pop_use;
      Ntot_dead_hypoxia_obs[i] = 0.0;
      Ntot_dead_buffer_obs[i] = 0.0;
      Ntot_dead_total_obs[i] = 0.0;
      Ntot_total_obs[i] = min_pop_use;
      Vmm3_live_obs[i] = burden_floor_use;
      Vmm3_dead_hypoxia_obs[i] = 0.0;
      Vmm3_dead_buffer_obs[i] = 0.0;
      Vmm3_dead_total_obs[i] = 0.0;
      Vmm3_total_obs[i] = burden_floor_use;
      continue;
    }
    const int idx = it->second;
    double nv_live = Ntot_live_at_step[static_cast<size_t>(idx)];
    double nv_dead_h = Ntot_dead_hypoxia_at_step[static_cast<size_t>(idx)];
    double nv_dead_b = Ntot_dead_buffer_at_step[static_cast<size_t>(idx)];
    double nv_dead = Ntot_dead_total_at_step[static_cast<size_t>(idx)];
    double nv_total = Ntot_total_at_step[static_cast<size_t>(idx)];
    double bv_live = Vmm3_live_at_step[static_cast<size_t>(idx)];
    double bv_dead_h = Vmm3_dead_hypoxia_at_step[static_cast<size_t>(idx)];
    double bv_dead_b = Vmm3_dead_buffer_at_step[static_cast<size_t>(idx)];
    double bv_dead = Vmm3_dead_total_at_step[static_cast<size_t>(idx)];
    double bv_total = Vmm3_total_at_step[static_cast<size_t>(idx)];
    if (!std::isfinite(nv_live)) nv_live = min_pop_use;
    if (!std::isfinite(nv_dead_h) || nv_dead_h < 0.0) nv_dead_h = 0.0;
    if (!std::isfinite(nv_dead_b) || nv_dead_b < 0.0) nv_dead_b = 0.0;
    if (!std::isfinite(nv_dead) || nv_dead < 0.0) nv_dead = (nv_dead_h + nv_dead_b);
    if (!std::isfinite(nv_total)) nv_total = min_pop_use;
    if (!std::isfinite(bv_live)) bv_live = burden_floor_use;
    if (!std::isfinite(bv_dead_h) || bv_dead_h < 0.0) bv_dead_h = 0.0;
    if (!std::isfinite(bv_dead_b) || bv_dead_b < 0.0) bv_dead_b = 0.0;
    if (!std::isfinite(bv_dead) || bv_dead < 0.0) bv_dead = (bv_dead_h + bv_dead_b);
    if (!std::isfinite(bv_total)) bv_total = burden_floor_use;
    Ntot_live_obs[i] = nv_live;
    Ntot_dead_hypoxia_obs[i] = nv_dead_h;
    Ntot_dead_buffer_obs[i] = nv_dead_b;
    Ntot_dead_total_obs[i] = nv_dead;
    Ntot_total_obs[i] = nv_total;
    Vmm3_live_obs[i] = bv_live;
    Vmm3_dead_hypoxia_obs[i] = bv_dead_h;
    Vmm3_dead_buffer_obs[i] = bv_dead_b;
    Vmm3_dead_total_obs[i] = bv_dead;
    Vmm3_total_obs[i] = bv_total;
  }

  NumericVector frac_N_live(R0, 0.0);
  double total_frac = 0.0;
  for (int i = 0; i < R0; ++i) {
    const double val = v_live[static_cast<size_t>(i)] + v_live[static_cast<size_t>(R0 + i)];
    frac_N_live[i] = val;
    total_frac += val;
  }
  if (total_frac > 0.0 && std::isfinite(total_frac)) {
    for (int i = 0; i < R0; ++i) frac_N_live[i] = frac_N_live[i] / total_frac;
  } else {
    const double u = 1.0 / static_cast<double>(R0);
    for (int i = 0; i < R0; ++i) frac_N_live[i] = u;
  }

  return List::create(
    // Backward-compatible aliases:
    _["Ntot_obs"] = Ntot_total_obs,
    _["Vmm3_obs"] = Vmm3_total_obs,
    _["frac_N"] = frac_N_live,
    // Explicit live/dead-consistent outputs:
    _["Ntot_live_obs"] = Ntot_live_obs,
    _["Ntot_dead_hypoxia_obs"] = Ntot_dead_hypoxia_obs,
    _["Ntot_dead_buffer_obs"] = Ntot_dead_buffer_obs,
    _["Ntot_dead_total_obs"] = Ntot_dead_total_obs,
    _["Ntot_total_obs"] = Ntot_total_obs,
    _["Vmm3_live_obs"] = Vmm3_live_obs,
    _["Vmm3_dead_hypoxia_obs"] = Vmm3_dead_hypoxia_obs,
    _["Vmm3_dead_buffer_obs"] = Vmm3_dead_buffer_obs,
    _["Vmm3_dead_total_obs"] = Vmm3_dead_total_obs,
    _["Vmm3_total_obs"] = Vmm3_total_obs,
    _["frac_N_live"] = frac_N_live,
    _["cache_g_build"] = cache_g_build,
    _["cache_g_hit"] = cache_g_hit,
    _["cache_g_hysteresis"] = cache_g_hysteresis,
    _["cache_bin_pct"] = o2_bin_use,
    _["cache_hysteresis_pct"] = o2_hyst_use
  );
}

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_objective_components_map
// Purpose: Compute MAP objective components using log-normal burden likelihood
//   and continuous single-cell ploidy mixture likelihood with balanced
//   2N/4N tumor-group aggregation for ploidy loss.
// Parameters:
//   - ploidy_z_list: Per-tumor continuous single-cell ploidy observations.
//   - mu_by_N: Representative ploidy value for each modeled N state.
//   - sigma_burden: Log-normal observation SD for burden.
//   - sigma_ploidy: Gaussian observation SD for single-cell ploidy.
// Returns:
//   List return value containing per-modality mean NLL components.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_objective_components_map(
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
    double O2_cap,
    bool o2_feedback,
    std::string o2_curve_type,
    double o2_init,
    double o2_rate,
    double o2_shape_v,
    double tau_O2,
    double o2_anchor_N,
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
    double p_const,
    double p_wgd,
    std::string boundary,
    double eps_tail,
    double beta_buffer,
    double n_exp,
    double smax,
    int N_unit,
    double beta_loss,
    double beta_gain,
    double beta_size,
    double alpha_o2,
    double gamma_growth,
    int growth_penalty_mode,
    double mu_hp,
    double k_clear,
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
  const int growth_mode_use = growth_penalty_mode;
  const bool growth_penalty_ploidy_use = (growth_mode_use & 1) != 0;
  const bool growth_penalty_hypoxia_use = (growth_mode_use & 2) != 0;

  std::vector<double> burden_losses;
  std::vector<double> ploidy_losses_2N;
  std::vector<double> ploidy_losses_4N;
  burden_losses.reserve(static_cast<size_t>(n_sc));
  ploidy_losses_2N.reserve(static_cast<size_t>(n_sc));
  ploidy_losses_4N.reserve(static_cast<size_t>(n_sc));

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

    List sim = cpp_o2simps_simulate_one(
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
      O2_cap,
      o2_feedback,
      o2_curve_type,
      o2_init,
      o2_rate,
      o2_shape_v,
      tau_O2,
      o2_anchor_N,
      o2_logN_eps,
      o2_cache_bin_pct,
      o2_cache_hysteresis_pct,
      o2_cache_profile,
      lam_min,
      lam_max,
      k_o,
      has_p_misseg,
      p_misseg,
      k_o_mis,
      false,
      0.0,
      0.0,
      p_const,
      p_wgd,
      boundary,
      eps_tail,
      beta_buffer,
      n_exp,
      smax,
      N_unit,
      beta_loss,
      beta_gain,
      beta_size,
      alpha_o2,
      gamma_growth,
      growth_penalty_ploidy_use,
      growth_penalty_hypoxia_use,
      mu_hp,
      k_clear,
      vol_by_N,
      burden_floor
    );

    NumericVector pred_burden = sim["Vmm3_total_obs"];
    NumericVector frac_N = sim["frac_N_live"];
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
        const double tumor_ploidy_loss = ploidy_nll_sum / static_cast<double>(ploidy_n);
        if (cohort == 0) {
          ploidy_losses_2N.push_back(tumor_ploidy_loss);
        } else {
          ploidy_losses_4N.push_back(tumor_ploidy_loss);
        }
      }
    }
  }

  const double L_b = burden_losses.empty()
    ? 0.0
    : std::accumulate(burden_losses.begin(), burden_losses.end(), 0.0) /
        static_cast<double>(burden_losses.size());
  const bool has_2N = !ploidy_losses_2N.empty();
  const bool has_4N = !ploidy_losses_4N.empty();
  const double L_p_2N = has_2N
    ? std::accumulate(ploidy_losses_2N.begin(), ploidy_losses_2N.end(), 0.0) /
        static_cast<double>(ploidy_losses_2N.size())
    : 0.0;
  const double L_p_4N = has_4N
    ? std::accumulate(ploidy_losses_4N.begin(), ploidy_losses_4N.end(), 0.0) /
        static_cast<double>(ploidy_losses_4N.size())
    : 0.0;
  const double L_p = (has_2N && has_4N)
    ? (0.5 * L_p_2N + 0.5 * L_p_4N)
    : (has_2N ? L_p_2N : (has_4N ? L_p_4N : 0.0));
  const int n_ploidy_total = static_cast<int>(ploidy_losses_2N.size() + ploidy_losses_4N.size());

  return List::create(
    _["L_b"] = L_b,
    _["L_p"] = L_p,
    _["n_burden"] = static_cast<int>(burden_losses.size()),
    _["n_ploidy"] = n_ploidy_total,
    _["n_ploidy_2N"] = static_cast<int>(ploidy_losses_2N.size()),
    _["n_ploidy_4N"] = static_cast<int>(ploidy_losses_4N.size()),
    _["L_p_2N"] = L_p_2N,
    _["L_p_4N"] = L_p_4N,
    _["cache_g_build"] = cache_g_build,
    _["cache_g_hit"] = cache_g_hit,
    _["cache_g_hysteresis"] = cache_g_hysteresis
  );
}
