#include <RcppEigen.h>
#include <algorithm>
#include <cctype>
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

constexpr double kNDip = 44.0;
constexpr int kPloidyDeathUniform = 0;
constexpr int kPloidyDeathDiploidNull = 1;
constexpr int kPloidyDeathPloidyRelated = 2;
constexpr int kStartWithPloidy = 0;
constexpr int kStartWithChrNumber = 1;

// -----------------------------------------------------------------------------
// Function: trim_lower_ascii_cpp
// Purpose: Normalize mode strings for robust parsing.
// Parameters:
//   - x: Raw mode string.
// Returns:
//   std::string return value containing lowercase-trimmed ASCII text.
// -----------------------------------------------------------------------------
inline std::string trim_lower_ascii_cpp(const std::string& x) {
  size_t b = 0;
  while (b < x.size() && std::isspace(static_cast<unsigned char>(x[b]))) ++b;
  size_t e = x.size();
  while (e > b && std::isspace(static_cast<unsigned char>(x[e - 1]))) --e;
  std::string out = x.substr(b, e - b);
  std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return out;
}

// -----------------------------------------------------------------------------
// Function: canonical_ploidy_o2_death_mode_cpp
// Purpose: Parse canonical ploidy_O2_death mode.
// Parameters:
//   - mode_raw: Requested mode string.
// Returns:
//   int return value containing one of:
//     0=uniform, 1=diploid_NULL, 2=ploidy_related.
// -----------------------------------------------------------------------------
inline int canonical_ploidy_o2_death_mode_cpp(const std::string& mode_raw) {
  const std::string s = trim_lower_ascii_cpp(mode_raw);
  if (s.empty()) {
    stop(
      "ploidy_O2_death must be supplied as a canonical mode string. ",
      "Allowed values are: uniform, diploid_NULL, ploidy_related."
    );
  }
  if (s == "ploidy_related") {
    return kPloidyDeathPloidyRelated;
  }
  if (s == "uniform") {
    return kPloidyDeathUniform;
  }
  if (s == "diploid_null") {
    return kPloidyDeathDiploidNull;
  }
  stop(
    "Invalid ploidy_O2_death mode: '", mode_raw,
    "'. Allowed canonical values are: uniform, diploid_NULL, ploidy_related."
  );
}

// -----------------------------------------------------------------------------
// Function: ploidy_o2_death_mode_name_cpp
// Purpose: Return canonical mode name for logging/cache consistency.
// Parameters:
//   - mode_code: Integer mode code.
// Returns:
//   std::string return value containing canonical mode name.
// -----------------------------------------------------------------------------
inline std::string ploidy_o2_death_mode_name_cpp(int mode_code) {
  if (mode_code == kPloidyDeathUniform) return "uniform";
  if (mode_code == kPloidyDeathDiploidNull) return "diploid_NULL";
  return "ploidy_related";
}

// -----------------------------------------------------------------------------
// Function: canonical_start_with_mode_cpp
// Purpose: Parse canonical start_with mode.
// Parameters:
//   - mode_raw: Requested mode string.
// Returns:
//   int return value containing one of:
//     0=ploidy, 1=chr_number.
// -----------------------------------------------------------------------------
inline int canonical_start_with_mode_cpp(const std::string& mode_raw) {
  const std::string s = trim_lower_ascii_cpp(mode_raw);
  if (s.empty()) {
    stop(
      "start_with must be supplied as a canonical mode string. ",
      "Allowed values are: ploidy, chr_number."
    );
  }
  if (s == "ploidy") return kStartWithPloidy;
  if (s == "chr_number") return kStartWithChrNumber;
  stop(
    "Invalid start_with mode: '", mode_raw,
    "'. Allowed canonical values are: ploidy, chr_number."
  );
}

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
// Function: hypoxia_weight_cpp
// Purpose: Compute Hill-type hypoxia weight used by growth/death modules.
// Parameters:
//   - O2_use: Oxygen level used by model rate functions.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double hypoxia_weight_cpp(double O2_use, double O2_crit_use, double n_O_use) {
  double o2_crit = (std::isfinite(O2_crit_use) && O2_crit_use >= 0.0) ? O2_crit_use : 1.0;
  o2_crit = std::max(o2_crit, 1e-12);
  const double o2 = clamp_o2_pct(O2_use);
  const double n_O = (std::isfinite(n_O_use) && n_O_use >= 0.0) ? n_O_use : 1.0;
  const double num = std::pow(o2_crit, n_O);
  const double den = num + std::pow(o2, n_O);
  if (!std::isfinite(den) || den <= 0.0) return 0.0;
  const double h = num / den;
  if (!std::isfinite(h)) return 0.0;
  return clamp01(h);
}

// -----------------------------------------------------------------------------
// Function: combined_resource_stress_cpp
// Purpose: Combine oxygen and glucose stress into one bounded resource stress.
// Parameters:
//   - O2_use: Oxygen level used by model rate functions.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double combined_resource_stress_cpp(
    double O2_use,
    double O2_crit_use,
    double n_O_use,
    bool glucose_enabled
) {
  const double h_o2 = hypoxia_weight_cpp(O2_use, O2_crit_use, n_O_use);
  if (!glucose_enabled) {
    return h_o2;
  }
  const double h_g = h_o2;
  const double combined = 1.0 - (1.0 - h_o2) * (1.0 - h_g);
  if (!std::isfinite(combined)) return 0.0;
  return clamp01(combined);
}

// -----------------------------------------------------------------------------
// Function: constant_p_wgd_cpp
// Purpose: Compute the constant per-division WGD probability.
// Parameters:
//   - p_wgd: Per-division whole-genome doubling probability.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double constant_p_wgd_cpp(double p_wgd) {
  if (!std::isfinite(p_wgd) || p_wgd <= 0.0) return 0.0;
  return clamp01(p_wgd);
}

// -----------------------------------------------------------------------------
// Function: lambda_base_from_o2_cpp
// Purpose: Compute baseline proliferation as the maximal growth rate.
// Parameters:
//   - lam_max: Maximal proliferation rate.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double lambda_base_from_o2_cpp(double lam_max) {
  double lam_base = std::isfinite(lam_max) ? lam_max : 0.0;
  if (!std::isfinite(lam_base) || lam_base < 0.0) lam_base = 0.0;
  return lam_base;
}

// -----------------------------------------------------------------------------
// Function: lambda_base_from_resource_cpp
// Purpose: Compute baseline proliferation as the maximal growth rate.
// Parameters:
//   - lam_max: Maximal proliferation rate.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double lambda_base_from_resource_cpp(double lam_max) {
  double lam_base = std::isfinite(lam_max) ? lam_max : 0.0;
  if (!std::isfinite(lam_base) || lam_base < 0.0) lam_base = 0.0;
  return lam_base;
}

// -----------------------------------------------------------------------------
// Function: lambda_eff_soft_o2_only_cpp
// Purpose: Compute effective growth rate with oxygen-stress growth damping.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - O2_use: Oxygen level used by model rate functions.
//   - lam_max: Maximal proliferation rate.
//   - alpha_o2: Oxygen-mediated growth-penalty strength.
//   - gamma_growth: Exponent for oxygen-mediated ploidy growth penalty.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double lambda_eff_soft_o2_only_cpp(
    int N_state,
    double O2_use,
    double lam_max,
    double alpha_o2,
    double gamma_growth,
    double O2_crit_use,
    double n_O_use
) {
  const double lam_base = lambda_base_from_o2_cpp(lam_max);
  if (lam_base <= 0.0) return 0.0;
  const double alpha_use = (std::isfinite(alpha_o2) && alpha_o2 > 0.0) ? alpha_o2 : 0.0;
  const double gamma_use = (std::isfinite(gamma_growth) && gamma_growth > 0.0) ? gamma_growth : 1.0;
  const double N_ratio = std::max(static_cast<double>(N_state) / kNDip, 0.0);
  const double h_o2 = hypoxia_weight_cpp(O2_use, O2_crit_use, n_O_use);
  const double denom = 1.0 + alpha_use * h_o2 * std::pow(N_ratio, gamma_use);
  if (!std::isfinite(denom) || denom <= 0.0) return 0.0;
  const double lam_eff = lam_base / denom;
  if (!std::isfinite(lam_eff) || lam_eff < 0.0) return 0.0;
  return lam_eff;
}

// -----------------------------------------------------------------------------
// Function: lambda_eff_soft_cpp
// Purpose: Compute effective growth rate with soft resource/ploidy penalty.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - O2_use: Oxygen level used by model rate functions.
//   - lam_max: Maximal proliferation rate.
//   - alpha_o2: Resource-mediated growth-penalty strength.
//   - gamma_growth: Exponent for resource-mediated ploidy growth penalty.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double lambda_eff_soft_cpp(
    int N_state,
    double O2_use,
    double lam_max,
    double alpha_o2,
    double gamma_growth,
    double O2_crit_use,
    double n_O_use,
    bool glucose_enabled
) {
  const double lam_base = lambda_base_from_resource_cpp(lam_max);
  if (lam_base <= 0.0) return 0.0;
  const double alpha_use = (std::isfinite(alpha_o2) && alpha_o2 > 0.0) ? alpha_o2 : 0.0;
  const double gamma_use = (std::isfinite(gamma_growth) && gamma_growth > 0.0) ? gamma_growth : 1.0;
  const double N_ratio = std::max(static_cast<double>(N_state) / kNDip, 0.0);
  const double h_resource = combined_resource_stress_cpp(
    O2_use,
    O2_crit_use,
    n_O_use,
    glucose_enabled
  );
  const double denom = 1.0 + alpha_use * h_resource * std::pow(N_ratio, gamma_use);
  if (!std::isfinite(denom) || denom <= 0.0) return 0.0;
  const double lam_eff = lam_base / denom;
  if (!std::isfinite(lam_eff) || lam_eff < 0.0) return 0.0;
  return lam_eff;
}

// -----------------------------------------------------------------------------
// Function: lambda_eff_runtime_cpp
// Purpose: Dispatch proliferation rate according to the O2_growth runtime switch.
// -----------------------------------------------------------------------------
inline double lambda_eff_runtime_cpp(
    int N_state,
    double O2_use,
    double lam_max,
    bool o2_growth,
    double alpha_o2,
    double gamma_growth,
    double O2_crit_use,
    double n_O_use,
    bool glucose_enabled
) {
  if (!glucose_enabled) {
    if (!o2_growth) {
      return lambda_base_from_o2_cpp(lam_max);
    }
    return lambda_eff_soft_o2_only_cpp(
      N_state,
      O2_use,
      lam_max,
      alpha_o2,
      gamma_growth,
      O2_crit_use,
      n_O_use
    );
  }
  if (!o2_growth) {
    return lambda_base_from_resource_cpp(lam_max);
  }
  return lambda_eff_soft_cpp(
    N_state,
    O2_use,
    lam_max,
    alpha_o2,
    gamma_growth,
    O2_crit_use,
    n_O_use,
    glucose_enabled
  );
}

// -----------------------------------------------------------------------------
// Function: mu_eff_soft_cpp
// Purpose: Compute effective hypoxia-linked death rate with optional ploidy modulation.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - O2_use: Oxygen level used by model rate functions.
//   - mu_hp: Hypoxia-linked high-ploidy death strength.
//   - gamma_mu: Exponent for high-ploidy hypoxia death above diploid reference.
//   - O2_crit_use: Hill critical oxygen scale.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
//   - ploidy_O2_death_mode: Mode code parsed from ploidy_O2_death.
//     Allowed values:
//       uniform       -> mu_eff = mu_hp * h(O2)
//       diploid_NULL  -> mu_eff = mu_hp * h(O2) * (1 + max(N/N_dip - 1, 0)^gamma_mu)
//       ploidy_related-> mu_eff = mu_hp * h(O2) * (N/N_dip)^gamma_mu
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double mu_eff_soft_cpp(
    int N_state,
    double O2_use,
    double mu_hp,
    double gamma_mu,
    double O2_crit_use,
    double n_O_use,
    int ploidy_O2_death_mode,
    bool glucose_enabled
) {
  const double mu_hp_use = (std::isfinite(mu_hp) && mu_hp > 0.0) ? mu_hp : 0.0;
  if (mu_hp_use <= 0.0) return 0.0;
  (void)glucose_enabled;
  const double h_o2 = hypoxia_weight_cpp(O2_use, O2_crit_use, n_O_use);
  if (h_o2 <= 0.0) return 0.0;
  if (ploidy_O2_death_mode == kPloidyDeathUniform) {
    const double mu_eff = mu_hp_use * h_o2;
    if (!std::isfinite(mu_eff) || mu_eff < 0.0) return 0.0;
    return mu_eff;
  }
  const double gamma_mu_use = (std::isfinite(gamma_mu) && gamma_mu > 0.0) ? gamma_mu : 1.0;
  const double n_ratio = std::max(static_cast<double>(N_state) / kNDip, 0.0);
  if (ploidy_O2_death_mode == kPloidyDeathDiploidNull) {
    const double above_dip = std::max(n_ratio - 1.0, 0.0);
    const double mu_eff = mu_hp_use * h_o2 * (1.0 + std::pow(above_dip, gamma_mu_use));
    if (!std::isfinite(mu_eff) || mu_eff < 0.0) return 0.0;
    return mu_eff;
  }
  const double mu_eff = mu_hp_use * h_o2 * std::pow(n_ratio, gamma_mu_use);
  if (!std::isfinite(mu_eff) || mu_eff < 0.0) return 0.0;
  return mu_eff;
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
//   - dropped_value: Optional accumulator for out-of-grid dropped transition mass.
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
    std::vector<double>& xx,
    double* dropped_value = nullptr
) {
  if (Np < row_min || Np > row_max) {
    if (bmode == 1) {
      const int Np2 = std::max(std::min(Np, row_max), row_min);
      ii.push_back(Np2 - row_min + 1);
      jj.push_back(col_1based);
      xx.push_back(value);
    } else if (dropped_value != nullptr) {
      // Under boundary=drop, out-of-grid offspring are not written to live
      // states; caller is responsible for routing this mass into dead buffer.
      *dropped_value += value;
    }
    return;
  }
  ii.push_back(Np - row_min + 1);
  jj.push_back(col_1based);
  xx.push_back(value);
}

// -----------------------------------------------------------------------------
// Function: binom_prob_int
// Purpose: Numerically robust binomial PMF evaluator for integer n in [0, N].
// Parameters:
//   - n: Number of successes.
//   - N: Number of Bernoulli trials.
//   - p: Per-trial success probability.
// Returns:
//   double return value containing PMF value in [0,1].
// -----------------------------------------------------------------------------
inline double binom_prob_int(int n, int N, double p) {
  if (N < 0) return 0.0;
  if (n < 0 || n > N) return 0.0;
  const double p_use = clamp01(p);
  const double v = R::dbinom(n, N, p_use, false);
  if (!std::isfinite(v) || v < 0.0) return 0.0;
  return v;
}

// -----------------------------------------------------------------------------
// Function: buffering_survival_modifier
// Purpose: Compute buffering-model survival for any missegregated
//   daughter with m affected chromosome copies.
// Parameters:
//   - q: Mother chromosome count state.
//   - m_misseg: Number of missegregated chromosome copies.
//   - buffer_smax: Maximum per-copy survival factor, constrained to [0,1].
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - n_chr: Number of modeled chromosome classes.
// Returns:
//   double return value containing survival modifier in [0,1].
// -----------------------------------------------------------------------------
inline double buffering_survival_modifier(
    int q,
    int m_misseg,
    double buffer_smax,
    double buffer_beta,
    double buffer_n_exp,
    int n_chr
) {
  if (m_misseg <= 0) return 1.0;
  if (q <= 0) return 1.0;
  const double smax_use = clamp01(std::isfinite(buffer_smax) ? buffer_smax : 1.0);
  const double beta_use = (std::isfinite(buffer_beta) && buffer_beta >= 0.0) ? buffer_beta : 0.0;
  const double n_exp_use = (std::isfinite(buffer_n_exp) && buffer_n_exp >= 0.0) ? buffer_n_exp : 1.0;
  const double n_chr_use = (n_chr > 0) ? static_cast<double>(n_chr) : 22.0;
  const double ratio = (2.0 * n_chr_use) / static_cast<double>(q);
  double sN = smax_use * std::exp(-beta_use * std::pow(std::max(ratio, 0.0), n_exp_use));
  sN = clamp01(std::isfinite(sN) ? sN : 0.0);
  if (sN <= 0.0) return 0.0;
  const double log_survival = static_cast<double>(m_misseg) * std::log(sN);
  if (!std::isfinite(log_survival)) return 0.0;
  const double survival = std::exp(log_survival);
  if (!std::isfinite(survival)) return 0.0;
  return clamp01(survival);
}

// -----------------------------------------------------------------------------
// Function: finalize_pr_delta_internal
// Purpose: Finalize coarse ploidy transition weights and dropped daughter mass.
// -----------------------------------------------------------------------------
inline void finalize_pr_delta_internal(
    int shift_offset,
    const std::vector<double>& shift_mass,
    double survivors_total,
    std::vector<int>& ts_out,
    std::vector<double>& prob_out,
    double& mass_dropped
) {
  if (!std::isfinite(survivors_total) || survivors_total < 0.0) survivors_total = 0.0;

  const double dead_daughters = std::max(0.0, 2.0 - survivors_total);
  mass_dropped = std::max(0.0, std::min(1.0, dead_daughters / 2.0));

  ts_out.clear();
  prob_out.clear();
  ts_out.reserve(shift_mass.size());
  prob_out.reserve(shift_mass.size());
  for (int t = -shift_offset; t <= shift_offset; ++t) {
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

// -----------------------------------------------------------------------------
// Function: o2simps_pr_delta_internal
// Purpose: Build symmetric buffering transition weights, matching the legacy O2
//   buffering model: gain and loss daughters share the same buffering survival.
// -----------------------------------------------------------------------------
void o2simps_pr_delta_internal(
    int N,
    double p,
    double eps_tail,
    double buffer_smax,
    double buffer_beta,
    double buffer_n_exp,
    int N_unit,
    std::vector<int>& ts_out,
    std::vector<double>& prob_out,
    double& mass_dropped
) {
  const int N_use = std::max(0, N);
  const double p_use = clamp01(p);
  const double eps_use = (std::isfinite(eps_tail) && eps_tail > 0.0) ? eps_tail : 0.0;
  const int shift_offset = N_use;
  std::vector<double> shift_mass(static_cast<size_t>(2 * shift_offset + 1), 0.0);

  double survivors_total = 0.0;
  for (int n = 0; n <= N_use; ++n) {
    const double pn = binom_prob_int(n, N_use, p_use);
    if (pn <= 0.0) continue;
    if (eps_use > 0.0 && pn < eps_use) continue;
    const int delta_gain = n;
    const int delta_loss = -n;
    const double s_buf = buffering_survival_modifier(
      N_use,
      n,
      buffer_smax,
      buffer_beta,
      buffer_n_exp,
      N_unit
    );
    const double w_gain = pn * s_buf;
    const double w_loss = pn * s_buf;
    if (w_gain > 0.0) shift_mass[static_cast<size_t>(shift_offset + delta_gain)] += w_gain;
    if (w_loss > 0.0) shift_mass[static_cast<size_t>(shift_offset + delta_loss)] += w_loss;
    survivors_total += (w_gain + w_loss);
  }

  finalize_pr_delta_internal(shift_offset, shift_mass, survivors_total, ts_out, prob_out, mass_dropped);
}

} // namespace

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_pr_delta_vec
// Purpose: Compute missegregation delta-kernel probabilities over ploidy shifts.
// Parameters:
//   - N: Ploidy state value or chromosome-copy count.
//   - p: Missegregation probability parameter.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_pr_delta_vec(
    int N,
    double p,
    double eps_tail = 1e-8,
    double buffer_smax = 1.0,
    double buffer_beta = 0.0,
    double buffer_n_exp = 1.0,
    int N_unit = 22
) {
  std::vector<int> ts;
  std::vector<double> prob;
  double mass_dropped = 0.0;

  o2simps_pr_delta_internal(
    N,
    p,
    eps_tail,
    buffer_smax,
    buffer_beta,
    buffer_n_exp,
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
// Function: cpp_o2simps_o2_window_supply
// Purpose: Compute oxygen target from viable burden using a logarithmic
//   supply-demand form with lower oxygen floor.
// Parameters:
//   - Ntot: Total predicted cell count (or burden proxy) at current time.
//   - o2_S0: Baseline oxygen supply level at near-zero burden (%).
//   - kappa_O: Function-specific input argument.
//   - o2_Nref: Fixed viable-cell scaling constant for demand normalization.
//   - o2_min: Lower floor for oxygen target in logarithmic supply-demand model (%).
// Returns:
//   NumericVector return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cpp_o2simps_o2_window_supply(
    NumericVector Ntot,
    double o2_S0 = 0.5,
    double kappa_O = 1.0,
    double o2_Nref = 1e6,
    double o2_min = 0.0
) {
  const int n = Ntot.size();
  NumericVector out(n);
  const double o2_S0_use = clamp_o2_pct((std::isfinite(o2_S0) && o2_S0 >= 0.0) ? o2_S0 : 0.5);
  const double kappa_use = (std::isfinite(kappa_O) && kappa_O > 0.0) ? kappa_O : 1.0;
  const double Nref_use = (std::isfinite(o2_Nref) && o2_Nref > 0.0) ? o2_Nref : 1e6;
  const double O2_min_use = clamp_o2_pct((std::isfinite(o2_min) && o2_min >= 0.0) ? o2_min : 0.0);

  for (int i = 0; i < n; ++i) {
    const double Nlive = (std::isfinite(Ntot[i]) && Ntot[i] > 0.0) ? Ntot[i] : 0.0;
    const double burden_ratio = Nlive / Nref_use;
    double o2_target = o2_S0_use - kappa_use * std::log1p(burden_ratio);
    if (!std::isfinite(o2_target)) o2_target = O2_min_use;
    o2_target = std::max(O2_min_use, o2_target);
    out[i] = clamp_o2_pct(o2_target);
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
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
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
    double buffer_smax = 1.0,
    double buffer_beta = 0.0,
    double buffer_n_exp = 1.0,
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
      buffer_smax,
      buffer_beta,
      buffer_n_exp,
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
      // Signed-shift contract: each (t, w) already encodes final daughter
      // displacement and mass, so write exactly once to N + t.
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
//   - N0min: Minimum ploidy state on source grid.
//   - N0max: Maximum ploidy state on source grid.
//   - N1min: Minimum ploidy state on doubled-state target grid.
//   - N1max: Maximum ploidy state on doubled-state target grid.
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
  if (R0 <= 0 || R1 <= 0) stop("Nmax must be >= Nmin for source and target grids");

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
//   - dropped_value: Optional accumulator for out-of-grid dropped transition mass.
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
    std::vector<double>& xx,
    double* dropped_value = nullptr
) {
  if (value == 0.0) return;

  int Np_use = Np;
  if (Np_use < row_min || Np_use > row_max) {
    if (bmode == 0) {
      if (dropped_value != nullptr) {
        // Under boundary=drop, out-of-grid offspring mass is accumulated here
        // so caller can add it to dead_buffer_rate.
        *dropped_value += value;
      }
      return;
    }
    Np_use = std::max(std::min(Np_use, row_max), row_min);
  }

  const int row_local_1based = Np_use - row_min + 1;
  ii.push_back(row_offset_1based + row_local_1based - 1);
  jj.push_back(col_1based);
  xx.push_back(value);
}

// -----------------------------------------------------------------------------
// Function: resolve_pmis_for_death
// Purpose: Internal helper used by the model fitting and simulation pipeline.
// Parameters:
//   - mu_eff: Effective death rate for the current ploidy state.
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale for death-linked missegregation.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double resolve_pmis_for_death(
    double mu_eff,
    double p_mis_base,
    double p_misseg,
    double k_o_mis
) {
  const double p_base = clamp01(std::isfinite(p_mis_base) ? p_mis_base : 1e-5);
  const double p_amp = (std::isfinite(p_misseg) && p_misseg > 0.0) ? p_misseg : 0.0;
  const double k_use = (std::isfinite(k_o_mis) && k_o_mis > 0.0) ? k_o_mis : 1e-12;
  const double mu_use = std::max(mu_eff, 0.0);
  const double frac = mu_use / (mu_use + k_use);
  const double delta_p = p_amp * frac;
  return clamp01(p_base + delta_p);
}

} // namespace

// -----------------------------------------------------------------------------
// Function: cpp_o2simps_build_G_for_o2_triplet
// Purpose: Build division-related live-state generator at the current oxygen/burden condition.
// Parameters:
//   - O2: Oxygen level used by model rate functions.
//   - N0min: Minimum ploidy state on the single chromosome-count grid.
//   - N0max: Maximum ploidy state on the single chromosome-count grid.
//   - N1min: Legacy argument kept for interface stability (unused).
//   - N1max: Legacy argument kept for interface stability (unused).
//   - lam_max: Maximal proliferation rate.
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale for death-linked missegregation (mu_eff scale).
//   - p_wgd: Constant per-division WGD probability.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_build_G_for_o2_triplet(
    double O2,
    double O2_crit,
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    double lam_max,
    double p_mis_base,
    double p_misseg,
    double k_o_mis,
    bool glucose = true,
    double p_wgd = 0.0,
    std::string boundary = "drop",
    double eps_tail = 1e-8,
    double buffer_smax = 1.0,
    double buffer_beta = 0.0,
    double buffer_n_exp = 1.0,
    int N_unit = 22,
    double beta_size = 0.0,
    bool O2_growth = true,
    double alpha_o2 = 0.0,
    double gamma_growth = 1.0,
    double mu_hp = 0.0,
    double gamma_mu = 1.0,
    double n_O = 1.0,
    std::string ploidy_O2_death = "diploid_NULL"
) {
  const int R = N0max - N0min + 1;
  if (R <= 0) stop("Nmax must be >= Nmin");
  (void)N1min;
  (void)N1max;

  const int bmode = boundary_mode(boundary);

  const double O2_use = clamp_o2_pct(O2);
  const double o2_crit_use = (std::isfinite(O2_crit) && O2_crit >= 0.0) ? O2_crit : 1.0;
  const bool o2_growth_use = O2_growth;
  const double alpha_o2_use = (std::isfinite(alpha_o2) && alpha_o2 > 0.0) ? alpha_o2 : 0.0;
  const double gamma_growth_use = (std::isfinite(gamma_growth) && gamma_growth > 0.0) ? gamma_growth : 1.0;
  const double mu_hp_use = (std::isfinite(mu_hp) && mu_hp > 0.0) ? mu_hp : 0.0;
  const double gamma_mu_use = (std::isfinite(gamma_mu) && gamma_mu > 0.0) ? gamma_mu : 1.0;
  const double n_O_use = (std::isfinite(n_O) && n_O >= 0.0) ? n_O : 1.0;
  const int ploidy_O2_death_mode_use = canonical_ploidy_o2_death_mode_cpp(ploidy_O2_death);
  const bool glucose_use = glucose;
  (void)beta_size;
  auto lam_for_N = [&](int N_state) -> double {
    return lambda_eff_runtime_cpp(
      N_state,
      O2_use,
      lam_max,
      o2_growth_use,
      alpha_o2_use,
      gamma_growth_use,
      o2_crit_use,
      n_O_use,
      glucose_use
    );
  };
  auto mu_for_N = [&](int N_state) -> double {
    return mu_eff_soft_cpp(
      N_state,
      O2_use,
      mu_hp_use,
      gamma_mu_use,
      o2_crit_use,
      n_O_use,
      ploidy_O2_death_mode_use,
      glucose_use
    );
  };

  const double pw = constant_p_wgd_cpp(p_wgd);

  std::vector<int> ii;
  std::vector<int> jj;
  std::vector<double> xx;
  std::vector<double> dead_buffer_rate(static_cast<size_t>(R), 0.0);
  std::vector<double> misseg_nonviable_rate(static_cast<size_t>(R), 0.0);
  std::vector<double> boundary_dropped_rate_vec(static_cast<size_t>(R), 0.0);
  std::vector<double> misseg_nonviable_division_prob(static_cast<size_t>(R), 0.0);
  std::vector<double> misseg_nonviable_daughters_per_division(static_cast<size_t>(R), 0.0);
  ii.reserve(static_cast<size_t>(R) * 20);
  jj.reserve(static_cast<size_t>(R) * 20);
  xx.reserve(static_cast<size_t>(R) * 20);

  for (int c = 0; c < R; ++c) {
    const int N = N0min + c;
    const double lam_N = lam_for_N(N);
    const double mu_N = mu_for_N(N);
    const double p_mis_N = resolve_pmis_for_death(
      mu_N,
      p_mis_base,
      p_misseg,
      k_o_mis
    );
    const int col_1based = c + 1;

    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    o2simps_pr_delta_internal(
      N,
      p_mis_N,
      eps_tail,
      buffer_smax,
      buffer_beta,
      buffer_n_exp,
      N_unit,
      ts,
      pr,
      mass_dropped
    );

    const double scale_mis = lam_N * (1.0 - pw);
    const double scale_wgd = lam_N * pw;
    double boundary_dropped_rate = 0.0;
    const double dead_daughters_per_division = std::max(0.0, std::min(2.0, 2.0 * mass_dropped));
    {
      // Expected nonviable offspring inflow from missegregation-linked loss.
      double nonviable_daughters_per_division = (1.0 - pw) * dead_daughters_per_division;
      if (!std::isfinite(nonviable_daughters_per_division) || nonviable_daughters_per_division < 0.0) {
        nonviable_daughters_per_division = 0.0;
      }
      if (nonviable_daughters_per_division > 2.0) nonviable_daughters_per_division = 2.0;
      double nonviable_rate = lam_N * nonviable_daughters_per_division;
      if (!std::isfinite(nonviable_rate) || nonviable_rate < 0.0) nonviable_rate = 0.0;
      misseg_nonviable_division_prob[static_cast<size_t>(col_1based - 1)] = std::min(nonviable_daughters_per_division, 1.0);
      misseg_nonviable_daughters_per_division[static_cast<size_t>(col_1based - 1)] = nonviable_daughters_per_division;
      misseg_nonviable_rate[static_cast<size_t>(col_1based - 1)] = nonviable_rate;
      dead_buffer_rate[static_cast<size_t>(col_1based - 1)] = nonviable_rate;
    }
    for (size_t k = 0; k < ts.size(); ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;
      // Signed-shift contract: each (t, w) already encodes final daughter
      // displacement and mass, so write exactly once to N + t.
      append_block_with_boundary(
        N + t,
        N0min,
        N0max,
        1,
        col_1based,
        scale_mis * w,
        bmode,
        ii,
        jj,
        xx,
        &boundary_dropped_rate
      );
    }

    ii.push_back(col_1based);
    jj.push_back(col_1based);
    xx.push_back(-lam_N);

    append_block_with_boundary(
      2 * N,
      N0min,
      N0max,
      1,
      col_1based,
      scale_wgd,
      bmode,
      ii,
      jj,
      xx,
      &boundary_dropped_rate
    );
    if (!std::isfinite(boundary_dropped_rate) || boundary_dropped_rate < 0.0) {
      boundary_dropped_rate = 0.0;
    }
    boundary_dropped_rate_vec[static_cast<size_t>(col_1based - 1)] = boundary_dropped_rate;
    // Boundary-drop losses are also routed to the dead buffer so mass does not
    // disappear from the represented system.
    dead_buffer_rate[static_cast<size_t>(col_1based - 1)] += boundary_dropped_rate;
  }

  return List::create(
    _["i"] = IntegerVector(ii.begin(), ii.end()),
    _["j"] = IntegerVector(jj.begin(), jj.end()),
    _["x"] = NumericVector(xx.begin(), xx.end()),
    _["nrow"] = R,
    _["ncol"] = R,
    _["dead_buffer_rate"] = NumericVector(dead_buffer_rate.begin(), dead_buffer_rate.end()),
    _["misseg_nonviable_rate"] = NumericVector(misseg_nonviable_rate.begin(), misseg_nonviable_rate.end()),
    _["boundary_dropped_rate"] = NumericVector(boundary_dropped_rate_vec.begin(), boundary_dropped_rate_vec.end()),
    _["misseg_nonviable_division_prob"] = NumericVector(
      misseg_nonviable_division_prob.begin(),
      misseg_nonviable_division_prob.end()
    ),
    _["misseg_nonviable_daughters_per_division"] = NumericVector(
      misseg_nonviable_daughters_per_division.begin(),
      misseg_nonviable_daughters_per_division.end()
    )
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
//   - N0min: Minimum ploidy state on source grid.
//   - N0max: Maximum ploidy state on source grid.
//   - N1min: Legacy interface argument (unused in single-layer dynamics).
//   - N1max: Legacy interface argument (unused in single-layer dynamics).
//   - lam_max: Maximal proliferation rate.
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale for death-linked missegregation (mu_eff scale).
//   - p_wgd: Constant per-division WGD probability.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
// Returns:
//   std::size_t return value containing the computed result.
// -----------------------------------------------------------------------------
inline std::size_t g_cache_signature_cpp(
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    double lam_max,
    double p_mis_base,
    double p_misseg,
    double k_o_mis,
    bool glucose,
    double p_wgd,
    const std::string& boundary,
    double eps_tail,
    double buffer_smax,
    double buffer_beta,
    double buffer_n_exp,
    double beta_size,
    bool o2_growth,
    double alpha_o2,
    double O2_crit,
    double gamma_growth,
    double mu_hp,
    double gamma_mu,
    double n_O,
    int ploidy_O2_death_mode,
    int N_unit
) {
  std::size_t seed = 0ULL;
  hash_combine_cpp(seed, N0min);
  hash_combine_cpp(seed, N0max);
  hash_combine_cpp(seed, N1min);
  hash_combine_cpp(seed, N1max);
  hash_combine_cpp(seed, bits_of_double_cpp(lam_max));
  hash_combine_cpp(seed, bits_of_double_cpp(p_mis_base));
  hash_combine_cpp(seed, bits_of_double_cpp(p_misseg));
  hash_combine_cpp(seed, bits_of_double_cpp(k_o_mis));
  hash_combine_cpp(seed, glucose ? 1 : 0);
  hash_combine_cpp(seed, bits_of_double_cpp(p_wgd));
  hash_combine_cpp(seed, boundary);
  hash_combine_cpp(seed, bits_of_double_cpp(eps_tail));
  hash_combine_cpp(seed, bits_of_double_cpp(buffer_smax));
  hash_combine_cpp(seed, bits_of_double_cpp(buffer_beta));
  hash_combine_cpp(seed, bits_of_double_cpp(buffer_n_exp));
  hash_combine_cpp(seed, bits_of_double_cpp(beta_size));
  hash_combine_cpp(seed, o2_growth ? 1 : 0);
  hash_combine_cpp(seed, bits_of_double_cpp(alpha_o2));
  hash_combine_cpp(seed, bits_of_double_cpp(O2_crit));
  hash_combine_cpp(seed, bits_of_double_cpp(gamma_growth));
  hash_combine_cpp(seed, bits_of_double_cpp(mu_hp));
  hash_combine_cpp(seed, bits_of_double_cpp(gamma_mu));
  hash_combine_cpp(seed, bits_of_double_cpp(n_O));
  hash_combine_cpp(seed, ploidy_O2_death_mode);
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
// Purpose: Compute oxygen target from viable burden using a logarithmic
//   supply-demand form with lower oxygen floor.
// Parameters:
//   - Ntot: Total predicted cell count (or burden proxy) at current time.
//   - o2_S0: Baseline oxygen supply level at near-zero burden (%).
//   - kappa_O: Function-specific input argument.
//   - o2_Nref: Fixed viable-cell scaling constant for demand normalization.
//   - o2_min: Lower floor for oxygen target in logarithmic supply-demand model (%).
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double o2_window_supply_scalar_cpp(
    double Ntot,
    double o2_S0,
    double kappa_O,
    double o2_Nref,
    double o2_min
) {
  const double o2_S0_use = clamp_o2_pct((std::isfinite(o2_S0) && o2_S0 >= 0.0) ? o2_S0 : 0.5);
  const double kappa_use = (std::isfinite(kappa_O) && kappa_O > 0.0) ? kappa_O : 1.0;
  const double Nref_use = (std::isfinite(o2_Nref) && o2_Nref > 0.0) ? o2_Nref : 1e6;
  const double O2_min_use = clamp_o2_pct((std::isfinite(o2_min) && o2_min >= 0.0) ? o2_min : 0.0);
  const double Nlive = (std::isfinite(Ntot) && Ntot > 0.0) ? Ntot : 0.0;
  const double burden_ratio = Nlive / Nref_use;
  double o2_target = o2_S0_use - kappa_use * std::log1p(burden_ratio);
  if (!std::isfinite(o2_target)) o2_target = O2_min_use;
  o2_target = std::max(O2_min_use, o2_target);
  return clamp_o2_pct(o2_target);
}

// -----------------------------------------------------------------------------
// Function: death_rate_for_N_cpp
// Purpose: Compute live->dead transfer death rate with optional ploidy modulation.
// Parameters:
//   - N_state: Ploidy state value or chromosome-copy count.
//   - O2_use: Oxygen level used by model rate functions.
//   - O2_crit_use: Hill critical oxygen scale.
//   - mu_hp_use: Hypoxia-linked high-ploidy death strength.
//   - gamma_mu_use: Exponent for high-ploidy hypoxia death above diploid reference.
//   - n_O_use: Hill exponent controlling oxygen-response steepness.
//   - ploidy_O2_death_mode: Parsed mode code for hypoxia-death ploidy dependence.
// Returns:
//   double return value containing the computed result.
// -----------------------------------------------------------------------------
inline double death_rate_for_N_cpp(
    int N_state,
    double O2_use,
    double O2_crit_use,
    double mu_hp_use,
    double gamma_mu_use,
    double n_O_use,
    int ploidy_O2_death_mode,
    bool glucose_enabled
) {
  return mu_eff_soft_cpp(
    N_state,
    O2_use,
    mu_hp_use,
    gamma_mu_use,
    O2_crit_use,
    n_O_use,
    ploidy_O2_death_mode,
    glucose_enabled
  );
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
//   - N0min: Minimum ploidy state on the single chromosome-count grid.
//   - N0max: Maximum ploidy state on the single chromosome-count grid.
//   - N1min: Legacy argument kept for interface stability (unused).
//   - N1max: Legacy argument kept for interface stability (unused).
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
//   - crowding_enabled: When false, disable crowding and force c(N)=1.
//   - crowding: Function-specific input argument.
//   - K: Function-specific input argument.
//   - min_pop: Function-specific input argument.
//   - O2_crit: Hill critical oxygen scale.
//   - o2_feedback: Function-specific input argument.
//   - o2_S0: Baseline oxygen supply level at near-zero burden (%).
//   - kappa_O: Function-specific input argument.
//   - tau_O2: Relaxation time constant controlling lag from O2 target to O2 effective.
//   - o2_Nref: Fixed viable-cell scaling constant for demand normalization.
//   - o2_min: Lower floor for oxygen target in logarithmic supply-demand model (%).
//   - eta_o2: Exponent for ploidy-weighted oxygen demand term (P/2)^eta_o2.
//   - o2_cache_bin_pct: Function-specific input argument.
//   - o2_cache_hysteresis_pct: Function-specific input argument.
//   - o2_cache_profile: Function-specific input argument.
//   - lam_max: Maximal proliferation rate.
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale for death-linked missegregation (mu_eff scale).
//   - p_wgd: Constant per-division WGD probability.
//   - boundary: Boundary handling mode when transitions leave the ploidy grid.
//   - eps_tail: Small truncation threshold for tail probabilities.
//   - buffer_smax: Maximum per-copy survival factor.
//   - buffer_beta: Ploidy-buffering strength.
//   - buffer_n_exp: Ploidy-buffering exponent.
//   - N_unit: Number of modeled chromosome classes for buffering scale.
//   - vol_by_N: Optional precomputed per-state cell volume lookup.
//   - burden_floor: Function-specific input argument.
//   - return_full_trajectory: When true, return per-observation live-state and O2
//     trajectories and do not short-circuit on extinction.
// Returns:
//   List return value containing the computed result.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_simulate_one(List sim_args) {
  NumericVector init_state = as<NumericVector>(sim_args["init_state"]);
  int N0min = as<int>(sim_args["N0min"]);
  int N0max = as<int>(sim_args["N0max"]);
  int N1min = as<int>(sim_args["N1min"]);
  int N1max = as<int>(sim_args["N1max"]);
  IntegerVector obs_steps = as<IntegerVector>(sim_args["obs_steps"]);
  int sim_end_step = as<int>(sim_args["sim_end_step"]);
  double DT = as<double>(sim_args["DT"]);
  double dose = as<double>(sim_args["dose"]);
  double dose_ref = as<double>(sim_args["dose_ref"]);
  double treat_day = as<double>(sim_args["treat_day"]);
  bool fit_treatment = as<bool>(sim_args["fit_treatment"]);
  double alpha = as<double>(sim_args["alpha"]);
  double gamma = as<double>(sim_args["gamma"]);
  double tx_mult_min = as<double>(sim_args["tx_mult_min"]);
  bool crowding_enabled = as<bool>(sim_args["crowding_enabled"]);
  std::string crowding = as<std::string>(sim_args["crowding"]);
  double K = as<double>(sim_args["K"]);
  double min_pop = as<double>(sim_args["min_pop"]);
  double O2_crit = as<double>(sim_args["O2_crit"]);
  bool o2_feedback = as<bool>(sim_args["o2_feedback"]);
  double o2_S0 = as<double>(sim_args["o2_S0"]);
  double kappa_O = as<double>(sim_args["kappa_O"]);
  double tau_O2 = as<double>(sim_args["tau_O2"]);
  double o2_Nref = as<double>(sim_args["o2_Nref"]);
  double o2_min = as<double>(sim_args["o2_min"]);
  double eta_o2 = as<double>(sim_args["eta_o2"]);
  double o2_cache_bin_pct = as<double>(sim_args["o2_cache_bin_pct"]);
  double o2_cache_hysteresis_pct = as<double>(sim_args["o2_cache_hysteresis_pct"]);
  bool o2_cache_profile = as<bool>(sim_args["o2_cache_profile"]);
  bool glucose = as<bool>(sim_args["glucose"]);
  double lam_max = as<double>(sim_args["lam_max"]);
  double p_mis_base = as<double>(sim_args["p_mis_base"]);
  double p_misseg = as<double>(sim_args["p_misseg"]);
  double k_o_mis = as<double>(sim_args["k_o_mis"]);
  double p_wgd = as<double>(sim_args["p_wgd"]);
  std::string boundary = as<std::string>(sim_args["boundary"]);
  double eps_tail = as<double>(sim_args["eps_tail"]);
  double buffer_smax = as<double>(sim_args["buffer_smax"]);
  double buffer_beta = as<double>(sim_args["buffer_beta"]);
  double buffer_n_exp = as<double>(sim_args["buffer_n_exp"]);
  int N_unit = as<int>(sim_args["N_unit"]);
  double beta_size = as<double>(sim_args["beta_size"]);
  bool O2_growth = as<bool>(sim_args["O2_growth"]);
  double alpha_o2 = as<double>(sim_args["alpha_o2"]);
  double gamma_growth = as<double>(sim_args["gamma_growth"]);
  double mu_hp = as<double>(sim_args["mu_hp"]);
  double gamma_mu = as<double>(sim_args["gamma_mu"]);
  double n_O = as<double>(sim_args["n_O"]);
  std::string ploidy_O2_death = as<std::string>(sim_args["ploidy_O2_death"]);
  std::string start_with = as<std::string>(sim_args["start_with"]);
  double k_clear = as<double>(sim_args["k_clear"]);
  NumericVector vol_by_N = as<NumericVector>(sim_args["vol_by_N"]);
  double burden_floor = as<double>(sim_args["burden_floor"]);
  bool return_full_trajectory = as<bool>(sim_args["return_full_trajectory"]);

  const int R = N0max - N0min + 1;
  if (R <= 0) stop("Nmax must be >= Nmin.");
  (void)N1min;
  (void)N1max;
  const int D = R;
  if (init_state.size() != D) {
    stop("init_state length mismatch: expected N0max - N0min + 1 live-state entries.");
  }
  if (vol_by_N.size() != R) stop("vol_by_N length mismatch.");

  if (!(crowding == "logistic" || crowding == "gompertz")) {
    stop("crowding must be logistic or gompertz.");
  }

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
  std::vector<double> O2_target_at_step(step_unique.size(), NA_REAL);
  std::vector<double> O2_eff_at_step(step_unique.size(), NA_REAL);
  std::vector<double> G_target_at_step(step_unique.size(), NA_REAL);
  std::vector<double> G_eff_at_step(step_unique.size(), NA_REAL);
  NumericMatrix live_state_at_step(
    return_full_trajectory ? static_cast<int>(step_unique.size()) : 0,
    return_full_trajectory ? R : 0
  );
  NumericMatrix dead_hypoxia_state_at_step(
    return_full_trajectory ? static_cast<int>(step_unique.size()) : 0,
    return_full_trajectory ? R : 0
  );
  NumericMatrix dead_buffer_state_at_step(
    return_full_trajectory ? static_cast<int>(step_unique.size()) : 0,
    return_full_trajectory ? R : 0
  );

  std::vector<double> v_live(static_cast<size_t>(D), 0.0);
  std::vector<double> v_dead_hypoxia(static_cast<size_t>(D), 0.0);
  std::vector<double> v_dead_buffer(static_cast<size_t>(D), 0.0);
  std::copy(init_state.begin(), init_state.end(), v_live.begin());
  std::vector<double> growth(static_cast<size_t>(D), 0.0);
  std::vector<double> death_flow_hypoxia(static_cast<size_t>(D), 0.0);
  std::vector<double> death_flow_buffer(static_cast<size_t>(D), 0.0);

  // Shared across scenario calls in the same worker process.
  // We keep one active parameter signature at a time so cache is reused
  // within one objective (same params), then reset when params change.
  static std::size_t active_sig = std::numeric_limits<std::size_t>::max();
  static std::unordered_map<std::size_t, SparseCacheEntry> shared_G_cache;

  const int ploidy_O2_death_mode_use = canonical_ploidy_o2_death_mode_cpp(ploidy_O2_death);
  const bool glucose_use = glucose;
  const int start_with_mode_use = canonical_start_with_mode_cpp(start_with);
  const double n_O_use = (std::isfinite(n_O) && n_O >= 0.0) ? n_O : 1.0;
  const std::size_t cur_sig = g_cache_signature_cpp(
    N0min,
    N0max,
    N1min,
    N1max,
    lam_max,
    p_mis_base,
    p_misseg,
    k_o_mis,
    glucose_use,
    p_wgd,
    boundary,
    eps_tail,
    buffer_smax,
    buffer_beta,
    buffer_n_exp,
    beta_size,
    O2_growth,
    alpha_o2,
    O2_crit,
    gamma_growth,
    mu_hp,
    gamma_mu,
    n_O_use,
    ploidy_O2_death_mode_use,
    N_unit
  );
  if (cur_sig != active_sig) {
    shared_G_cache.clear();
    shared_G_cache.reserve(256);
    active_sig = cur_sig;
  }

  const double O2_crit_use = (std::isfinite(O2_crit) && O2_crit >= 0.0) ? O2_crit : 1.0;
  const double o2_S0_use = (std::isfinite(o2_S0) ? o2_S0 : 0.5);
  const double kappa_O_use = (std::isfinite(kappa_O) ? kappa_O : 1.0);
  const double tau_use = (std::isfinite(tau_O2) && tau_O2 > 0.0) ? tau_O2 : 2.0;
  const double alpha_tau = 1.0 - std::exp(-DT_use / tau_use);
  const double o2_Nref_use = (std::isfinite(o2_Nref) && o2_Nref > 0.0) ? o2_Nref : 1e6;
  const double o2_min_use = clamp_o2_pct((std::isfinite(o2_min) && o2_min >= 0.0) ? o2_min : 0.0);
  const double eta_o2_use = (std::isfinite(eta_o2) && eta_o2 >= 0.0) ? eta_o2 : 1.0;
  const double N_unit_use = (N_unit > 0) ? static_cast<double>(N_unit) : 22.0;
  const double o2_bin_use = (std::isfinite(o2_cache_bin_pct) && o2_cache_bin_pct > 0.0) ? o2_cache_bin_pct : 1e-3;
  const double o2_hyst_use = (std::isfinite(o2_cache_hysteresis_pct) && o2_cache_hysteresis_pct >= 0.0) ? o2_cache_hysteresis_pct : 0.0;
  const double k_clear_use = (std::isfinite(k_clear) && k_clear >= 0.0) ? k_clear : 0.0;
  const std::string ploidy_O2_death_mode_name = ploidy_o2_death_mode_name_cpp(ploidy_O2_death_mode_use);
  (void) o2_cache_profile;
  int cache_g_build = 0;
  int cache_g_hit = 0;
  int cache_g_hysteresis = 0;
  bool has_last_key = false;
  std::size_t last_key = 0ULL;
  double last_o2_eff = 0.0;
  double last_g_eff = 0.0;
  std::vector<double> o2_demand_weight(static_cast<size_t>(D), 1.0);
  for (int i = 0; i < D; ++i) {
    const double N_state = static_cast<double>(N0min + i);
    const double ratio = (start_with_mode_use == kStartWithChrNumber)
      ? std::max(0.0, N_state / kNDip)
      : std::max(0.0, (N_state / N_unit_use) / 2.0);
    double w = std::pow(ratio, eta_o2_use);
    if (!std::isfinite(w) || w < 0.0) w = 1.0;
    o2_demand_weight[static_cast<size_t>(i)] = w;
  }
  auto compute_o2_demand_eff = [&](const std::vector<double>& live_state) -> double {
    double s = 0.0;
    for (int i = 0; i < D; ++i) {
      s += live_state[static_cast<size_t>(i)] * o2_demand_weight[static_cast<size_t>(i)];
    }
    if (!std::isfinite(s) || s < 0.0) s = 0.0;
    return s;
  };
  double O2_state = clamp_o2_pct(o2_S0_use);
  if (o2_feedback) {
    O2_state = o2_window_supply_scalar_cpp(
      compute_o2_demand_eff(v_live),
      o2_S0_use,
      kappa_O_use,
      o2_Nref_use,
      o2_min_use
    );
    O2_state = clamp_o2_pct(O2_state);
  }

  for (int step = 0; step <= final_step; ++step) {
    const double Ntot_live_now = vector_sum_cpp(v_live);
    const double Ntot_live_eff_for_o2_now = compute_o2_demand_eff(v_live);
    double O2_target_now = clamp_o2_pct(o2_S0_use);
    if (o2_feedback) {
      O2_target_now = o2_window_supply_scalar_cpp(
        Ntot_live_eff_for_o2_now,
        o2_S0_use,
        kappa_O_use,
        o2_Nref_use,
        o2_min_use
      );
    }
    O2_target_now = clamp_o2_pct(O2_target_now);
    const double O2_eff_now = clamp_o2_pct(O2_state);
    double G_target_now = glucose_use ? O2_target_now : 0.0;
    G_target_now = clamp_o2_pct(G_target_now);
    const double G_eff_now = glucose_use ? O2_eff_now : 0.0;

    auto it_obs = step_to_idx.find(step);
    if (it_obs != step_to_idx.end()) {
      const int idx = it_obs->second;
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
      for (int i = 0; i < R; ++i) {
        const size_t idx_i = static_cast<size_t>(i);
        const double n_live_i = v_live[idx_i];
        const double n_dead_h_i = v_dead_hypoxia[idx_i];
        const double n_dead_b_i = v_dead_buffer[idx_i];
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
      O2_target_at_step[static_cast<size_t>(idx)] = O2_target_now;
      O2_eff_at_step[static_cast<size_t>(idx)] = O2_eff_now;
      G_target_at_step[static_cast<size_t>(idx)] = G_target_now;
      G_eff_at_step[static_cast<size_t>(idx)] = G_eff_now;
      if (return_full_trajectory) {
        for (int i = 0; i < R; ++i) {
          live_state_at_step(idx, i) = v_live[static_cast<size_t>(i)];
          dead_hypoxia_state_at_step(idx, i) = v_dead_hypoxia[static_cast<size_t>(i)];
          dead_buffer_state_at_step(idx, i) = v_dead_buffer[static_cast<size_t>(i)];
        }
      }
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

    O2_state = O2_state + alpha_tau * (O2_target_now - O2_state);
    double O2_eff = clamp_o2_pct(O2_state);
    double G_eff = glucose_use ? O2_eff : 0.0;

    const int o2_key = quantize_o2_key(O2_eff, o2_bin_use);
    const int g_key = o2_key;
    std::size_t gkey = 0ULL;
    hash_combine_cpp(gkey, o2_key);
    hash_combine_cpp(gkey, g_key);
    if (o2_hyst_use > 0.0 && has_last_key &&
        std::abs(O2_eff - last_o2_eff) <= o2_hyst_use &&
        std::abs(G_eff - last_g_eff) <= o2_hyst_use) {
      gkey = last_key;
      ++cache_g_hysteresis;
    }
    auto itG = shared_G_cache.find(gkey);
    if (itG == shared_G_cache.end()) {
      const List tri = cpp_o2simps_build_G_for_o2_triplet(
        O2_eff,
        O2_crit,
        N0min,
        N0max,
        N1min,
        N1max,
        lam_max,
        p_mis_base,
        p_misseg,
        k_o_mis,
        glucose_use,
        p_wgd,
        boundary,
        eps_tail,
        buffer_smax,
        buffer_beta,
        buffer_n_exp,
        N_unit,
        beta_size,
        O2_growth,
        alpha_o2,
        gamma_growth,
        mu_hp,
        gamma_mu,
        n_O_use,
        ploidy_O2_death_mode_name
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
    last_g_eff = G_eff;

    sparse_mv_cpp(itG->second, v_live, growth);
    // Crowding scales division-linked activity when the config switch is enabled.
    double crowd_mult = 1.0;
    if (crowding_enabled) {
      if (crowding == "logistic") {
        crowd_mult = std::max(0.0, 1.0 - (Ntot_live_now / K_use));
      } else {
        crowd_mult = std::exp(-(Ntot_live_now / K_use));
      }
    }
    if (!std::isfinite(crowd_mult) || crowd_mult < 0.0) crowd_mult = 0.0;
    const double scalar = DT_use * crowd_mult * tx_mult;
    for (int i = 0; i < D; ++i) {
      const int N_state = N0min + i;
      const double mu_i = death_rate_for_N_cpp(
        N_state,
        O2_eff,
        O2_crit_use,
        mu_hp,
        gamma_mu,
        n_O_use,
        ploidy_O2_death_mode_use,
        glucose_use
      );
      const double src_live = v_live[static_cast<size_t>(i)];
      // Hypoxia death flow is independent of crowding/treatment scaling.
      double flow_h_i = DT_use * mu_i * src_live;
      if (!std::isfinite(flow_h_i) || flow_h_i < 0.0) flow_h_i = 0.0;
      death_flow_hypoxia[static_cast<size_t>(i)] = flow_h_i;
      double db_rate_i = 0.0;
      if (static_cast<size_t>(i) < itG->second.dead_buffer_rate.size()) {
        db_rate_i = itG->second.dead_buffer_rate[static_cast<size_t>(i)];
      }
      if (!std::isfinite(db_rate_i) || db_rate_i < 0.0) db_rate_i = 0.0;
      // Dead-buffer inflow from mitotic nonviability + boundary-dropped offspring
      // (not a continuous death hazard).
      double flow_b_i = scalar * db_rate_i * src_live;
      if (!std::isfinite(flow_b_i) || flow_b_i < 0.0) flow_b_i = 0.0;
      death_flow_buffer[static_cast<size_t>(i)] = flow_b_i;
    }
    for (size_t i = 0; i < v_live.size(); ++i) {
      const double next_v = v_live[i] + scalar * growth[i] - death_flow_hypoxia[i];
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
    if (!return_full_trajectory &&
        vector_sum_cpp(v_live) <= min_pop_use &&
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
  NumericVector O2_target_obs(obs_v.size(), NA_REAL);
  NumericVector O2_eff_obs(obs_v.size(), NA_REAL);
  NumericVector G_target_obs(obs_v.size(), NA_REAL);
  NumericVector G_eff_obs(obs_v.size(), NA_REAL);
  NumericMatrix live_state_obs(
    return_full_trajectory ? static_cast<int>(obs_v.size()) : 0,
    return_full_trajectory ? R : 0
  );
  NumericMatrix dead_hypoxia_state_obs(
    return_full_trajectory ? static_cast<int>(obs_v.size()) : 0,
    return_full_trajectory ? R : 0
  );
  NumericMatrix dead_buffer_state_obs(
    return_full_trajectory ? static_cast<int>(obs_v.size()) : 0,
    return_full_trajectory ? R : 0
  );
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
      O2_target_obs[i] = o2_window_supply_scalar_cpp(
        0.0,
        o2_S0_use,
        kappa_O_use,
        o2_Nref_use,
        o2_min_use
      );
      O2_eff_obs[i] = O2_target_obs[i];
      G_target_obs[i] = glucose_use ? O2_target_obs[i] : 0.0;
      G_eff_obs[i] = G_target_obs[i];
      if (return_full_trajectory) {
        for (int j = 0; j < R; ++j) {
          live_state_obs(i, j) = 0.0;
          dead_hypoxia_state_obs(i, j) = 0.0;
          dead_buffer_state_obs(i, j) = 0.0;
        }
      }
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
    double o2_target_val = O2_target_at_step[static_cast<size_t>(idx)];
    double o2_eff_val = O2_eff_at_step[static_cast<size_t>(idx)];
    double g_target_val = G_target_at_step[static_cast<size_t>(idx)];
    double g_eff_val = G_eff_at_step[static_cast<size_t>(idx)];
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
    if (!std::isfinite(o2_target_val)) {
      o2_target_val = o2_window_supply_scalar_cpp(
        nv_live,
        o2_S0_use,
        kappa_O_use,
        o2_Nref_use,
        o2_min_use
      );
    }
    if (!std::isfinite(o2_eff_val)) o2_eff_val = o2_target_val;
    if (!std::isfinite(g_target_val)) {
      g_target_val = glucose_use ? o2_target_val : 0.0;
    }
    if (!std::isfinite(g_eff_val)) g_eff_val = glucose_use ? o2_eff_val : 0.0;
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
    O2_target_obs[i] = o2_target_val;
    O2_eff_obs[i] = o2_eff_val;
    G_target_obs[i] = g_target_val;
    G_eff_obs[i] = g_eff_val;
    if (return_full_trajectory) {
      for (int j = 0; j < R; ++j) {
        live_state_obs(i, j) = live_state_at_step(idx, j);
        dead_hypoxia_state_obs(i, j) = dead_hypoxia_state_at_step(idx, j);
        dead_buffer_state_obs(i, j) = dead_buffer_state_at_step(idx, j);
      }
    }
  }

  NumericVector frac_N_live(R, 0.0);
  double total_frac = 0.0;
  for (int i = 0; i < R; ++i) {
    const double val = v_live[static_cast<size_t>(i)];
    frac_N_live[i] = val;
    total_frac += val;
  }
  if (total_frac > 0.0 && std::isfinite(total_frac)) {
    for (int i = 0; i < R; ++i) frac_N_live[i] = frac_N_live[i] / total_frac;
  } else {
    const double u = 1.0 / static_cast<double>(R);
    for (int i = 0; i < R; ++i) frac_N_live[i] = u;
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
    _["O2_target_obs"] = O2_target_obs,
    _["O2_eff_obs"] = O2_eff_obs,
    _["G_target_obs"] = G_target_obs,
    _["G_eff_obs"] = G_eff_obs,
    _["frac_N_live"] = frac_N_live,
    _["live_state_obs"] = live_state_obs,
    _["dead_hypoxia_state_obs"] = dead_hypoxia_state_obs,
    _["dead_buffer_state_obs"] = dead_buffer_state_obs,
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
//   - p_mis_base: Baseline per-chromosome missegregation probability.
//   - p_misseg: Death-linked missegregation amplification scale.
//   - k_o_mis: Half-saturation scale on mu_eff for missegregation amplification.
//   - o2_min: Lower floor for oxygen target in logarithmic supply-demand model (%).
// Returns:
//   List return value containing per-modality mean NLL components.
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
List cpp_o2simps_objective_components_map(
    List scenario_data,
    List objective_data,
    List state_data,
    List sim_args
) {
  IntegerVector cohort_code = as<IntegerVector>(scenario_data["cohort_code"]);
  NumericVector dose_vec = as<NumericVector>(scenario_data["dose_vec"]);
  NumericVector treat_day_vec = as<NumericVector>(scenario_data["treat_day_vec"]);
  List obs_steps_list = as<List>(scenario_data["obs_steps_list"]);
  IntegerVector sim_end_step_vec = as<IntegerVector>(scenario_data["sim_end_step_vec"]);
  List obs_burden_list = as<List>(scenario_data["obs_burden_list"]);
  List keep_burden_list = as<List>(scenario_data["keep_burden_list"]);
  List ploidy_z_list = as<List>(scenario_data["ploidy_z_list"]);
  NumericVector init_mult_vec = scenario_data.containsElementNamed("init_mult_vec")
    ? as<NumericVector>(scenario_data["init_mult_vec"])
    : NumericVector(cohort_code.size(), 1.0);

  NumericVector mu_by_N = as<NumericVector>(objective_data["mu_by_N"]);
  double sigma_burden = as<double>(objective_data["sigma_burden"]);
  double sigma_ploidy = as<double>(objective_data["sigma_ploidy"]);

  NumericVector init_state_2N = as<NumericVector>(state_data["init_state_2N"]);
  NumericVector init_state_4N = as<NumericVector>(state_data["init_state_4N"]);
  int N0min = as<int>(state_data["N0min"]);
  int N0max = as<int>(state_data["N0max"]);
  int N1min = as<int>(state_data["N1min"]);
  int N1max = as<int>(state_data["N1max"]);
  int N_unit = as<int>(state_data["N_unit"]);
  NumericVector vol_by_N = as<NumericVector>(state_data["vol_by_N"]);

  double DT = as<double>(sim_args["DT"]);
  double dose_ref = as<double>(sim_args["dose_ref"]);
  bool fit_treatment = as<bool>(sim_args["fit_treatment"]);
  double alpha = as<double>(sim_args["alpha"]);
  double gamma = as<double>(sim_args["gamma"]);
  double tx_mult_min = as<double>(sim_args["tx_mult_min"]);
  bool crowding_enabled = as<bool>(sim_args["crowding_enabled"]);
  std::string crowding = as<std::string>(sim_args["crowding"]);
  double K = as<double>(sim_args["K"]);
  double min_pop = as<double>(sim_args["min_pop"]);
  double O2_crit = as<double>(sim_args["O2_crit"]);
  bool o2_feedback = as<bool>(sim_args["o2_feedback"]);
  double o2_S0 = as<double>(sim_args["o2_S0"]);
  double kappa_O = as<double>(sim_args["kappa_O"]);
  double tau_O2 = as<double>(sim_args["tau_O2"]);
  double o2_Nref = as<double>(sim_args["o2_Nref"]);
  double o2_min = as<double>(sim_args["o2_min"]);
  double eta_o2 = as<double>(sim_args["eta_o2"]);
  double o2_cache_bin_pct = as<double>(sim_args["o2_cache_bin_pct"]);
  double o2_cache_hysteresis_pct = as<double>(sim_args["o2_cache_hysteresis_pct"]);
  bool o2_cache_profile = as<bool>(sim_args["o2_cache_profile"]);
  bool glucose = as<bool>(sim_args["glucose"]);
  double lam_max = as<double>(sim_args["lam_max"]);
  double p_mis_base = as<double>(sim_args["p_mis_base"]);
  double p_misseg = as<double>(sim_args["p_misseg"]);
  double k_o_mis = as<double>(sim_args["k_o_mis"]);
  double p_wgd = as<double>(sim_args["p_wgd"]);
  std::string boundary = as<std::string>(sim_args["boundary"]);
  double eps_tail = as<double>(sim_args["eps_tail"]);
  double buffer_smax = as<double>(sim_args["buffer_smax"]);
  double buffer_beta = as<double>(sim_args["buffer_beta"]);
  double buffer_n_exp = as<double>(sim_args["buffer_n_exp"]);
  double beta_size = as<double>(sim_args["beta_size"]);
  double alpha_o2 = as<double>(sim_args["alpha_o2"]);
  double gamma_growth = as<double>(sim_args["gamma_growth"]);
  double mu_hp = as<double>(sim_args["mu_hp"]);
  double gamma_mu = as<double>(sim_args["gamma_mu"]);
  double n_O = as<double>(sim_args["n_O"]);
  std::string ploidy_O2_death = as<std::string>(sim_args["ploidy_O2_death"]);
  std::string start_with = as<std::string>(sim_args["start_with"]);
  double k_clear = as<double>(sim_args["k_clear"]);
  double burden_log_eps = as<double>(sim_args["burden_log_eps"]);

  const int n_sc = cohort_code.size();
  if (dose_vec.size() != n_sc || treat_day_vec.size() != n_sc ||
      obs_steps_list.size() != n_sc || sim_end_step_vec.size() != n_sc ||
      obs_burden_list.size() != n_sc || keep_burden_list.size() != n_sc ||
      ploidy_z_list.size() != n_sc || init_mult_vec.size() != n_sc) {
    stop("Scenario containers must have consistent length.");
  }

  const double log_eps_use =
    (std::isfinite(burden_log_eps) && burden_log_eps > 0.0) ? burden_log_eps : 1e-15;
  const double sigma_b_use =
    (std::isfinite(sigma_burden) && sigma_burden > 0.0) ? sigma_burden : 0.35;
  const double sigma_p_use =
    (std::isfinite(sigma_ploidy) && sigma_ploidy > 0.0) ? sigma_ploidy : 0.08;
  const double prob_eps = 1e-300;
  const bool o2_growth_use = !(std::isfinite(alpha_o2) && alpha_o2 < 0.0);
  const double alpha_o2_use = std::fabs(alpha_o2);
  std::vector<double> burden_losses;
  std::vector<double> ploidy_losses_2N;
  std::vector<double> ploidy_losses_4N;
  burden_losses.reserve(static_cast<size_t>(n_sc));
  ploidy_losses_2N.reserve(static_cast<size_t>(n_sc));
  ploidy_losses_4N.reserve(static_cast<size_t>(n_sc));
  double burden_nll_total = 0.0;
  double ploidy_nll_total = 0.0;
  int burden_obs_total = 0;
  int ploidy_obs_total = 0;

  int cache_g_build = 0;
  int cache_g_hit = 0;
  int cache_g_hysteresis = 0;

  for (int i = 0; i < n_sc; ++i) {
    const int cohort = cohort_code[i];
    NumericVector init_state = (cohort == 0) ? init_state_2N : init_state_4N;
    const double init_mult = (std::isfinite(init_mult_vec[i]) && init_mult_vec[i] > 0.0) ? init_mult_vec[i] : 1.0;
    if (init_mult != 1.0) {
      init_state = clone(init_state);
      for (int j = 0; j < init_state.size(); ++j) {
        init_state[j] = init_state[j] * init_mult;
      }
    }
    IntegerVector obs_steps = as<IntegerVector>(obs_steps_list[i]);
    NumericVector obs_burden = as<NumericVector>(obs_burden_list[i]);
    LogicalVector keep_day = as<LogicalVector>(keep_burden_list[i]);
    NumericVector ploidy_z = as<NumericVector>(ploidy_z_list[i]);

    List sim_one_args = clone(sim_args);
    sim_one_args["init_state"] = init_state;
    sim_one_args["N0min"] = N0min;
    sim_one_args["N0max"] = N0max;
    sim_one_args["N1min"] = N1min;
    sim_one_args["N1max"] = N1max;
    sim_one_args["obs_steps"] = obs_steps;
    sim_one_args["sim_end_step"] = sim_end_step_vec[i];
    sim_one_args["dose"] = dose_vec[i];
    sim_one_args["treat_day"] = treat_day_vec[i];
    sim_one_args["fit_treatment"] = fit_treatment;
    sim_one_args["alpha"] = alpha;
    sim_one_args["gamma"] = gamma;
    sim_one_args["tx_mult_min"] = tx_mult_min;
    sim_one_args["crowding_enabled"] = crowding_enabled;
    sim_one_args["crowding"] = crowding;
    sim_one_args["K"] = K;
    sim_one_args["min_pop"] = min_pop;
    sim_one_args["O2_crit"] = O2_crit;
    sim_one_args["o2_feedback"] = o2_feedback;
    sim_one_args["o2_S0"] = o2_S0;
    sim_one_args["kappa_O"] = kappa_O;
    sim_one_args["tau_O2"] = tau_O2;
    sim_one_args["o2_Nref"] = o2_Nref;
    sim_one_args["o2_min"] = o2_min;
    sim_one_args["eta_o2"] = eta_o2;
    sim_one_args["o2_cache_bin_pct"] = o2_cache_bin_pct;
    sim_one_args["o2_cache_hysteresis_pct"] = o2_cache_hysteresis_pct;
    sim_one_args["o2_cache_profile"] = o2_cache_profile;
    sim_one_args["glucose"] = glucose;
    sim_one_args["lam_max"] = lam_max;
    sim_one_args["p_mis_base"] = p_mis_base;
    sim_one_args["p_misseg"] = p_misseg;
    sim_one_args["k_o_mis"] = k_o_mis;
    sim_one_args["p_wgd"] = p_wgd;
    sim_one_args["boundary"] = boundary;
    sim_one_args["eps_tail"] = eps_tail;
    sim_one_args["buffer_smax"] = buffer_smax;
    sim_one_args["buffer_beta"] = buffer_beta;
    sim_one_args["buffer_n_exp"] = buffer_n_exp;
    sim_one_args["N_unit"] = N_unit;
    sim_one_args["beta_size"] = beta_size;
    sim_one_args["O2_growth"] = o2_growth_use;
    sim_one_args["alpha_o2"] = alpha_o2_use;
    sim_one_args["gamma_growth"] = gamma_growth;
    sim_one_args["mu_hp"] = mu_hp;
    sim_one_args["gamma_mu"] = gamma_mu;
    sim_one_args["n_O"] = n_O;
    sim_one_args["ploidy_O2_death"] = ploidy_O2_death;
    sim_one_args["start_with"] = start_with;
    sim_one_args["k_clear"] = k_clear;
    sim_one_args["vol_by_N"] = vol_by_N;
    sim_one_args["burden_floor"] = log_eps_use;
    sim_one_args["return_full_trajectory"] = false;

    List sim = cpp_o2simps_simulate_one(sim_one_args);

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
      burden_nll_total += burden_nll_sum;
      burden_obs_total += burden_n;
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
        ploidy_nll_total += ploidy_nll_sum;
        ploidy_obs_total += ploidy_n;
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
  const double objective_burden_neg2loglik_raw = 2.0 * burden_nll_total;
  const double objective_ploidy_neg2loglik_raw = 2.0 * ploidy_nll_total;

  return List::create(
    _["L_b"] = L_b,
    _["L_p"] = L_p,
    _["burden_nll_total"] = burden_nll_total,
    _["ploidy_nll_total"] = ploidy_nll_total,
    _["objective_burden_neg2loglik_raw"] = objective_burden_neg2loglik_raw,
    _["objective_ploidy_neg2loglik_raw"] = objective_ploidy_neg2loglik_raw,
    _["n_burden"] = static_cast<int>(burden_losses.size()),
    _["n_burden_obs_total"] = burden_obs_total,
    _["n_ploidy"] = n_ploidy_total,
    _["n_ploidy_obs_total"] = ploidy_obs_total,
    _["n_ploidy_2N"] = static_cast<int>(ploidy_losses_2N.size()),
    _["n_ploidy_4N"] = static_cast<int>(ploidy_losses_4N.size()),
    _["L_p_2N"] = L_p_2N,
    _["L_p_4N"] = L_p_4N,
    _["cache_g_build"] = cache_g_build,
    _["cache_g_hit"] = cache_g_hit,
    _["cache_g_hysteresis"] = cache_g_hysteresis
  );
}
