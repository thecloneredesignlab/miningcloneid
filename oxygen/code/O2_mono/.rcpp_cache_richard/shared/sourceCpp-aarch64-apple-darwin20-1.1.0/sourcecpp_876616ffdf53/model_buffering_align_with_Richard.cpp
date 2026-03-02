#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

namespace {

inline int boundary_mode(const std::string& boundary) {
  if (boundary == "drop") return 0;
  if (boundary == "absorb_minmax") return 1;
  stop("boundary must be one of: drop, absorb_minmax");
}

inline double clip01_cpp(double x) {
  if (!std::isfinite(x)) return 0.0;
  if (x < 0.0) return 0.0;
  if (x > 1.0) return 1.0;
  return x;
}

inline double clip_o2_pct_cpp(double x) {
  if (!std::isfinite(x)) return 0.0;
  if (x < 0.0) return 0.0;
  if (x > 100.0) return 100.0;
  return x;
}

inline int quantize_o2_key(double o2) {
  // O2 is represented in percentage [0,100]; quantize at 0.1% bins.
  int k = static_cast<int>(std::llround(o2 * 10.0));
  if (k < 0) k = 0;
  if (k > 1000) k = 1000;
  return k;
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

namespace {

struct SparseEntries {
  std::vector<int> row;
  std::vector<int> col;
  std::vector<double> val;
};

SparseEntries build_G_sparse_for_o2(
    double O2,
    int N0min,
    int N0max,
    int N1min,
    int N1max,
    double lam_min,
    double lam_max,
    double k_o,
    double p_misseg,
    double k_o_mis,
    double beta_buffer,
    double n_exp,
    double smax,
    double p_wgd,
    int N_unit
) {
  const int R0 = N0max - N0min + 1;
  const int R1 = N1max - N1min + 1;

  const double k_o_use = std::max(k_o, 1e-12);
  const double frac = O2 / (O2 + k_o_use);
  const double lam = std::max(0.0, lam_min + (lam_max - lam_min) * frac);
  const double lam0 = lam;
  const double lam1 = lam;

  const double k_o_mis_use = std::max(k_o_mis, 1e-12);
  const double p_mis = clip01_cpp(p_misseg * (1.0 - (O2 / (O2 + k_o_mis_use))));
  const double pw = clip01_cpp(p_wgd);

  SparseEntries G;
  G.row.reserve(static_cast<size_t>(R0 + R1) * 64);
  G.col.reserve(static_cast<size_t>(R0 + R1) * 64);
  G.val.reserve(static_cast<size_t>(R0 + R1) * 64);

  // Pre layer columns (UL + LL)
  for (int c0 = 0; c0 < R0; ++c0) {
    const int N = N0min + c0;
    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    richard_pr_delta_internal(
      N,
      p_mis,
      1e-8,
      beta_buffer,
      n_exp,
      smax,
      N_unit,
      ts,
      pr,
      mass_dropped
    );
    (void)mass_dropped;

    const double scale_pre = lam0 * (1.0 - pw);
    const int col_state = c0;
    for (size_t k = 0; k < ts.size(); ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;

      if (t == 0) {
        const int Np = N;
        if (Np >= N0min && Np <= N0max) {
          G.row.push_back(Np - N0min);
          G.col.push_back(col_state);
          G.val.push_back(scale_pre * (2.0 * w));
        }
      } else {
        const int Np1 = N + t;
        const int Np2 = N - t;
        if (Np1 >= N0min && Np1 <= N0max) {
          G.row.push_back(Np1 - N0min);
          G.col.push_back(col_state);
          G.val.push_back(scale_pre * w);
        }
        if (Np2 >= N0min && Np2 <= N0max) {
          G.row.push_back(Np2 - N0min);
          G.col.push_back(col_state);
          G.val.push_back(scale_pre * w);
        }
      }
    }

    // -L0 diagonal
    G.row.push_back(col_state);
    G.col.push_back(col_state);
    G.val.push_back(-lam0);

    // LL WGD contribution (BW * SW * L0), boundary="drop"
    const int Nw = 2 * N;
    if (Nw >= N1min && Nw <= N1max) {
      G.row.push_back(R0 + (Nw - N1min));
      G.col.push_back(col_state);
      G.val.push_back(lam0 * pw);
    }
  }

  // Post layer columns (LR)
  for (int c1 = 0; c1 < R1; ++c1) {
    const int N = N1min + c1;
    std::vector<int> ts;
    std::vector<double> pr;
    double mass_dropped = 0.0;
    richard_pr_delta_internal(
      N,
      p_mis,
      1e-8,
      beta_buffer,
      n_exp,
      smax,
      N_unit,
      ts,
      pr,
      mass_dropped
    );
    (void)mass_dropped;

    const int col_state = R0 + c1;
    const double scale_post = lam1;
    for (size_t k = 0; k < ts.size(); ++k) {
      const int t = ts[k];
      const double w = pr[k];
      if (w == 0.0) continue;

      if (t == 0) {
        const int Np = N;
        if (Np >= N1min && Np <= N1max) {
          G.row.push_back(R0 + (Np - N1min));
          G.col.push_back(col_state);
          G.val.push_back(scale_post * (2.0 * w));
        }
      } else {
        const int Np1 = N + t;
        const int Np2 = N - t;
        if (Np1 >= N1min && Np1 <= N1max) {
          G.row.push_back(R0 + (Np1 - N1min));
          G.col.push_back(col_state);
          G.val.push_back(scale_post * w);
        }
        if (Np2 >= N1min && Np2 <= N1max) {
          G.row.push_back(R0 + (Np2 - N1min));
          G.col.push_back(col_state);
          G.val.push_back(scale_post * w);
        }
      }
    }

    // -L1 diagonal
    G.row.push_back(col_state);
    G.col.push_back(col_state);
    G.val.push_back(-lam1);
  }

  return G;
}

inline void sparse_mv_accumulate(
    const SparseEntries& G,
    const std::vector<double>& x,
    std::vector<double>& y
) {
  std::fill(y.begin(), y.end(), 0.0);
  const size_t nnz = G.val.size();
  for (size_t e = 0; e < nnz; ++e) {
    y[static_cast<size_t>(G.row[e])] += G.val[e] * x[static_cast<size_t>(G.col[e])];
  }
}

inline double vector_sum(const std::vector<double>& x) {
  return std::accumulate(x.begin(), x.end(), 0.0);
}

} // namespace

// [[Rcpp::export]]
List cpp_richard_simulate_one(
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
    bool o2_feedback,
    double o2_min,
    double h_O2,
    double K_O2,
    double lam_min,
    double lam_max,
    double k_o,
    double p_misseg,
    double k_o_mis,
    double beta_buffer,
    double n_exp,
    double smax,
    double p_wgd,
    int N_unit,
    NumericVector vol_by_N,
    double burden_floor
) {
  const int R0 = N0max - N0min + 1;
  const int R1 = N1max - N1min + 1;
  if (R0 <= 0 || R1 <= 0) stop("Invalid grid bounds.");
  if (init_state.size() != (R0 + R1)) stop("init_state length mismatch.");
  if (vol_by_N.size() != R0) stop("vol_by_N length mismatch.");

  const bool crowd_logistic = (crowding == "logistic");
  const bool crowd_gompertz = (crowding == "gompertz");
  if (!crowd_logistic && !crowd_gompertz) stop("crowding must be logistic or gompertz");

  const double DT_use = (std::isfinite(DT) && DT > 0.0) ? DT : 0.5;
  const double K_use = (std::isfinite(K) && K > 0.0) ? K : 1e12;
  const double min_pop_use = (std::isfinite(min_pop) && min_pop > 0.0) ? min_pop : 1e-12;
  const double dose_ref_use = (std::isfinite(dose_ref) && dose_ref > 0.0) ? dose_ref : 30.0;
  double dose_scaled = dose / dose_ref_use;
  if (!std::isfinite(dose_scaled) || dose_scaled < 0.0) dose_scaled = 0.0;
  const double tx_min_use = clip01_cpp(tx_mult_min);

  const double O2_base_use = clip_o2_pct_cpp(O2_base);
  const double O2_floor = clip_o2_pct_cpp(o2_min);
  const double h_use = (std::isfinite(h_O2) && h_O2 > 0.0) ? h_O2 : 1.0;
  const double K_O2_use = (std::isfinite(K_O2) && K_O2 > 0.0) ? K_O2 : 1e12;
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

  std::unordered_map<int, SparseEntries> G_cache;
  G_cache.reserve(64);

  for (int step = 0; step <= final_step; ++step) {
    auto it_obs = step_to_idx.find(step);
    if (it_obs != step_to_idx.end()) {
      const int idx = it_obs->second;
      const double Ntot_now = vector_sum(v);
      Ntot_at_step[static_cast<size_t>(idx)] = Ntot_now;
      double burden = 0.0;
      for (int i = 0; i < R0; ++i) {
        const double n_i = v[static_cast<size_t>(i)] + v[static_cast<size_t>(R0 + i)];
        burden += n_i * vol_by_N[i];
      }
      Vmm3_at_step[static_cast<size_t>(idx)] = burden;
    }

    if (step >= final_step) break;

    const double t = static_cast<double>(step) * DT_use;
    double tx_mult = 1.0;
    if (fit_treatment) {
      if (!(t < treat_day) && dose > 0.0) {
        const double dose_term = std::pow(dose_scaled, gamma);
        tx_mult = std::exp(-alpha * dose_term);
      } else {
        tx_mult = 1.0;
      }
      if (!std::isfinite(tx_mult)) tx_mult = tx_min_use;
      if (tx_mult < tx_min_use) tx_mult = tx_min_use;
      if (tx_mult > 1.0) tx_mult = 1.0;
    }

    const double Ntot = vector_sum(v);
    double O2_eff = O2_base_use;
    if (o2_feedback) {
      O2_eff = O2_floor + (O2_base_use - O2_floor) / (1.0 + std::pow(Ntot / K_O2_use, h_use));
      O2_eff = clip_o2_pct_cpp(O2_eff);
    }

    const int key = quantize_o2_key(O2_eff);
    auto itG = G_cache.find(key);
    if (itG == G_cache.end()) {
      SparseEntries Gnew = build_G_sparse_for_o2(
        O2_eff,
        N0min,
        N0max,
        N1min,
        N1max,
        lam_min,
        lam_max,
        k_o,
        p_misseg,
        k_o_mis,
        beta_buffer,
        n_exp,
        smax,
        p_wgd,
        N_unit
      );
      itG = G_cache.emplace(key, std::move(Gnew)).first;
    }

    double crowd = 1.0;
    if (crowd_logistic) {
      crowd = std::max(0.0, 1.0 - Ntot / K_use);
    } else {
      crowd = std::exp(-Ntot / K_use);
    }
    double scalar = DT_use * crowd * tx_mult;
    if (!std::isfinite(scalar)) scalar = 0.0;

    sparse_mv_accumulate(itG->second, v, growth);
    for (size_t i = 0; i < v.size(); ++i) {
      v[i] += scalar * growth[i];
      if (!std::isfinite(v[i]) || v[i] < 0.0) v[i] = 0.0;
    }
    if (vector_sum(v) <= min_pop_use) break;
  }

  NumericVector Ntot_obs(obs_v.size());
  NumericVector Vmm3_obs(obs_v.size());
  for (size_t i = 0; i < obs_v.size(); ++i) {
    const int st = obs_v[i];
    auto it = step_to_idx.find(st);
    double nval = NA_REAL;
    double bval = NA_REAL;
    if (it != step_to_idx.end()) {
      const int idx = it->second;
      nval = Ntot_at_step[static_cast<size_t>(idx)];
      bval = Vmm3_at_step[static_cast<size_t>(idx)];
    }
    if (!std::isfinite(nval)) nval = min_pop_use;
    if (!std::isfinite(bval)) bval = burden_floor_use;
    Ntot_obs[static_cast<int>(i)] = nval;
    Vmm3_obs[static_cast<int>(i)] = bval;
  }

  NumericVector frac_N(R0);
  double sum_frac = 0.0;
  for (int i = 0; i < R0; ++i) {
    const double f = v[static_cast<size_t>(i)] + v[static_cast<size_t>(R0 + i)];
    frac_N[i] = f;
    sum_frac += f;
  }
  if (sum_frac > 0.0 && std::isfinite(sum_frac)) {
    for (int i = 0; i < R0; ++i) frac_N[i] = frac_N[i] / sum_frac;
  } else {
    for (int i = 0; i < R0; ++i) frac_N[i] = 1.0 / static_cast<double>(R0);
  }

  return List::create(
    _["Ntot_obs"] = Ntot_obs,
    _["Vmm3_obs"] = Vmm3_obs,
    _["frac_N"] = frac_N
  );
}


#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cpp_richard_pr_delta_vec
List cpp_richard_pr_delta_vec(int N, double p, double eps_tail, double beta_buffer, double n_exp, double smax, int N_unit);
RcppExport SEXP sourceCpp_1_cpp_richard_pr_delta_vec(SEXP NSEXP, SEXP pSEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP) {
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
RcppExport SEXP sourceCpp_1_cpp_richard_build_B_total_triplet(SEXP NminSEXP, SEXP NmaxSEXP, SEXP p_vecSEXP, SEXP boundarySEXP, SEXP eps_tailSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP N_unitSEXP) {
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
RcppExport SEXP sourceCpp_1_cpp_richard_build_B_WGD_triplet(SEXP N0minSEXP, SEXP N0maxSEXP, SEXP N1minSEXP, SEXP N1maxSEXP, SEXP boundarySEXP, SEXP wgd_valueSEXP) {
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
// cpp_richard_simulate_one
List cpp_richard_simulate_one(NumericVector init_state, int N0min, int N0max, int N1min, int N1max, IntegerVector obs_steps, int sim_end_step, double DT, double dose, double dose_ref, double treat_day, bool fit_treatment, double alpha, double gamma, double tx_mult_min, std::string crowding, double K, double min_pop, double O2_base, bool o2_feedback, double o2_min, double h_O2, double K_O2, double lam_min, double lam_max, double k_o, double p_misseg, double k_o_mis, double beta_buffer, double n_exp, double smax, double p_wgd, int N_unit, NumericVector vol_by_N, double burden_floor);
RcppExport SEXP sourceCpp_1_cpp_richard_simulate_one(SEXP init_stateSEXP, SEXP N0minSEXP, SEXP N0maxSEXP, SEXP N1minSEXP, SEXP N1maxSEXP, SEXP obs_stepsSEXP, SEXP sim_end_stepSEXP, SEXP DTSEXP, SEXP doseSEXP, SEXP dose_refSEXP, SEXP treat_daySEXP, SEXP fit_treatmentSEXP, SEXP alphaSEXP, SEXP gammaSEXP, SEXP tx_mult_minSEXP, SEXP crowdingSEXP, SEXP KSEXP, SEXP min_popSEXP, SEXP O2_baseSEXP, SEXP o2_feedbackSEXP, SEXP o2_minSEXP, SEXP h_O2SEXP, SEXP K_O2SEXP, SEXP lam_minSEXP, SEXP lam_maxSEXP, SEXP k_oSEXP, SEXP p_missegSEXP, SEXP k_o_misSEXP, SEXP beta_bufferSEXP, SEXP n_expSEXP, SEXP smaxSEXP, SEXP p_wgdSEXP, SEXP N_unitSEXP, SEXP vol_by_NSEXP, SEXP burden_floorSEXP) {
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
    Rcpp::traits::input_parameter< bool >::type o2_feedback(o2_feedbackSEXP);
    Rcpp::traits::input_parameter< double >::type o2_min(o2_minSEXP);
    Rcpp::traits::input_parameter< double >::type h_O2(h_O2SEXP);
    Rcpp::traits::input_parameter< double >::type K_O2(K_O2SEXP);
    Rcpp::traits::input_parameter< double >::type lam_min(lam_minSEXP);
    Rcpp::traits::input_parameter< double >::type lam_max(lam_maxSEXP);
    Rcpp::traits::input_parameter< double >::type k_o(k_oSEXP);
    Rcpp::traits::input_parameter< double >::type p_misseg(p_missegSEXP);
    Rcpp::traits::input_parameter< double >::type k_o_mis(k_o_misSEXP);
    Rcpp::traits::input_parameter< double >::type beta_buffer(beta_bufferSEXP);
    Rcpp::traits::input_parameter< double >::type n_exp(n_expSEXP);
    Rcpp::traits::input_parameter< double >::type smax(smaxSEXP);
    Rcpp::traits::input_parameter< double >::type p_wgd(p_wgdSEXP);
    Rcpp::traits::input_parameter< int >::type N_unit(N_unitSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vol_by_N(vol_by_NSEXP);
    Rcpp::traits::input_parameter< double >::type burden_floor(burden_floorSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_richard_simulate_one(init_state, N0min, N0max, N1min, N1max, obs_steps, sim_end_step, DT, dose, dose_ref, treat_day, fit_treatment, alpha, gamma, tx_mult_min, crowding, K, min_pop, O2_base, o2_feedback, o2_min, h_O2, K_O2, lam_min, lam_max, k_o, p_misseg, k_o_mis, beta_buffer, n_exp, smax, p_wgd, N_unit, vol_by_N, burden_floor));
    return rcpp_result_gen;
END_RCPP
}
