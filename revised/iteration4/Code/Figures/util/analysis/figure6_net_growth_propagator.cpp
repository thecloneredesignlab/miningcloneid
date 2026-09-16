#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

namespace {

NumericVector f6ng_integerize(const NumericVector& probability, double total) {
  if (!R_finite(total) || total < 1 || total > 9007199254740991.0)
    stop("Passage population outside exact-double integer range.");
  const double target = std::floor(total + 0.5);
  const double offset = R::unif_rand();
  NumericVector count(probability.size());
  long double cumulative = 0;
  double previous = 0;
  for (int i = 0; i < probability.size(); ++i) {
    cumulative += static_cast<long double>(target) * probability[i];
    const double upper = i + 1 == probability.size() ? target :
      std::min(target, static_cast<double>(std::floor(cumulative + offset)));
    count[i] = upper - previous;
    if (count[i] < 0) stop("Negative integerized state.");
    previous = upper;
  }
  return count;
}

NumericVector f6ng_sample_counts(const NumericVector& counts, double size) {
  double remaining = sum(counts), needed = size;
  if (size != std::floor(size) || size < 1 || size > remaining)
    stop("Cannot sample requested inoculum without replacement.");
  NumericVector sample(counts.size());
  for (int i = 0; i < counts.size(); ++i) {
    const double good = counts[i], bad = remaining - good;
    double value = 0;
    if (needed > 0 && good > 0) {
      value = bad <= 0 ? needed :
        (needed == remaining ? good : R::rhyper(good, bad, needed));
    }
    if (!R_finite(value) || value < 0 || value > good || value > needed)
      stop("Invalid hypergeometric sample.");
    sample[i] = value;
    needed -= value;
    remaining -= good;
  }
  if (needed != 0 || sum(sample) != size)
    stop("Sampling total is not conserved.");
  return sample;
}

void f6ng_validate_common(const NumericMatrix& step,
    const NumericMatrix& initial_state, const NumericVector& ploidy_grid,
    const NumericVector& net_live_rate, int n_day) {
  const int n = step.nrow();
  if (step.ncol() != n || initial_state.nrow() != n ||
      ploidy_grid.size() != n || net_live_rate.size() != n || n_day < 0)
    stop("Invalid net-growth propagation dimensions.");
  for (int r = 0; r < n; ++r) {
    if (!R_finite(ploidy_grid[r]) || !R_finite(net_live_rate[r]))
      stop("Non-finite ploidy or net-live-rate vector.");
  }
  for (int c = 0; c < initial_state.ncol(); ++c) {
    double total = 0;
    for (int r = 0; r < n; ++r) {
      if (!R_finite(initial_state(r, c)) || initial_state(r, c) < 0)
        stop("Invalid initial state.");
      total += initial_state(r, c);
    }
    if (!R_finite(total) || std::abs(total - 1) > 1e-10)
      stop("Initial composition must sum to one.");
  }
}

void f6ng_state_summaries(const NumericMatrix& state, int column,
    const NumericVector& ploidy_grid, const NumericVector& net_live_rate,
    double& mean_ploidy, double& net_growth_rate) {
  mean_ploidy = 0;
  net_growth_rate = 0;
  for (int r = 0; r < state.nrow(); ++r) {
    mean_ploidy += state(r, column) * ploidy_grid[r];
    net_growth_rate += state(r, column) * net_live_rate[r];
  }
}

}  // namespace

// The instantaneous population net-live growth rate is 1' M f(t), where
// f(t) is the normalized live-state composition and colSums(M) is supplied as
// net_live_rate. Units therefore match the model generator: day^-1.
// [[Rcpp::export(rng = false)]]
List f6ng_propagate_continuous_cpp(NumericMatrix step,
    NumericMatrix initial_state, NumericVector ploidy_grid,
    NumericVector net_live_rate, int n_day) {
  f6ng_validate_common(step, initial_state, ploidy_grid, net_live_rate, n_day);
  const int n = step.nrow(), k = initial_state.ncol();
  NumericMatrix state = clone(initial_state), next(n, k);
  NumericMatrix mean_ploidy(k, n_day + 1), net_growth_rate(k, n_day + 1);
  for (int c = 0; c < k; ++c) {
    double mean = 0, rate = 0;
    f6ng_state_summaries(state, c, ploidy_grid, net_live_rate, mean, rate);
    mean_ploidy(c, 0) = mean;
    net_growth_rate(c, 0) = rate;
  }
  for (int day = 1; day <= n_day; ++day) {
    const char trans = 'N';
    const double alpha = 1, beta = 0;
    F77_CALL(dgemm)(&trans, &trans, &n, &k, &n, &alpha, step.begin(), &n,
      state.begin(), &n, &beta, next.begin(), &n FCONE FCONE);
    for (int c = 0; c < k; ++c) {
      double total = 0;
      for (int r = 0; r < n; ++r) {
        double value = next(r, c);
        if (value < 0 && value > -1e-10) value = 0;
        if (!R_finite(value) || value < 0) stop("Invalid propagated state.");
        next(r, c) = value;
        total += value;
      }
      if (!R_finite(total) || total <= 0) stop("Non-positive growth factor.");
      for (int r = 0; r < n; ++r) next(r, c) /= total;
      double mean = 0, rate = 0;
      f6ng_state_summaries(next, c, ploidy_grid, net_live_rate, mean, rate);
      mean_ploidy(c, day) = mean;
      net_growth_rate(c, day) = rate;
    }
    std::copy(next.begin(), next.end(), state.begin());
    if (day % 100 == 0) checkUserInterrupt();
  }
  return List::create(_["mean_ploidy"] = mean_ploidy,
    _["net_growth_rate"] = net_growth_rate, _["final_state"] = state);
}

// This reproduces the Figure 6 threshold-triggered stochastic passage rule.
// Rates on a passage day use the sampled post-passage composition. The
// exogenous dilution itself is not counted as negative biological growth.
// [[Rcpp::export(rng = false)]]
List f6ng_propagate_stochastic_cpp(NumericMatrix step,
    NumericMatrix initial_state, NumericVector ploidy_grid,
    NumericVector net_live_rate, int n_day, double seed_cells,
    double target_cells, IntegerMatrix rng_states,
    NumericVector initial_log_population, int start_day = 0) {
  f6ng_validate_common(step, initial_state, ploidy_grid, net_live_rate, n_day);
  const int n = step.nrow(), k = initial_state.ncol();
  if (rng_states.nrow() != 7 || rng_states.ncol() != k ||
      initial_log_population.size() != k || start_day < 0 ||
      seed_cells < 1 || seed_cells != std::floor(seed_cells) ||
      target_cells <= seed_cells || target_cells != std::floor(target_cells))
    stop("Invalid stochastic net-growth propagation contract.");
  NumericMatrix state = clone(initial_state), next(n, k);
  NumericMatrix mean_ploidy(k, n_day + 1), net_growth_rate(k, n_day + 1);
  IntegerMatrix streams = clone(rng_states);
  NumericVector log_population = clone(initial_log_population);
  IntegerVector passage_count(k), first_passage_day(k, NA_INTEGER),
    last_passage_day(k, NA_INTEGER);
  const double log_target = std::log(target_cells), log_seed = std::log(seed_cells);
  for (int c = 0; c < k; ++c) {
    if (!R_finite(log_population[c])) stop("Invalid initial log population.");
    double mean = 0, rate = 0;
    f6ng_state_summaries(state, c, ploidy_grid, net_live_rate, mean, rate);
    mean_ploidy(c, 0) = mean;
    net_growth_rate(c, 0) = rate;
  }
  Environment global = Environment::global_env();
  const bool had_seed = global.exists(".Random.seed");
  RObject caller_seed = R_NilValue;
  if (had_seed) caller_seed = clone(as<IntegerVector>(global[".Random.seed"]));
  for (int day = 1; day <= n_day; ++day) {
    const char trans = 'N';
    const double alpha = 1, beta = 0;
    F77_CALL(dgemm)(&trans, &trans, &n, &k, &n, &alpha, step.begin(), &n,
      state.begin(), &n, &beta, next.begin(), &n FCONE FCONE);
    for (int c = 0; c < k; ++c) {
      double total = 0;
      for (int r = 0; r < n; ++r) {
        double value = next(r, c);
        if (value < 0 && value > -1e-10) value = 0;
        if (!R_finite(value) || value < 0) stop("Invalid propagated state.");
        next(r, c) = value;
        total += value;
      }
      if (!R_finite(total) || total <= 0) stop("Non-positive growth factor.");
      for (int r = 0; r < n; ++r) next(r, c) /= total;
      log_population[c] += std::log(total);
      if (log_population[c] >= log_target) {
        const double actual_population = std::exp(log_population[c]);
        IntegerVector stream = streams(_, c);
        global[".Random.seed"] = stream;
        GetRNGstate();
        NumericVector probability = next(_, c);
        NumericVector rounded = f6ng_integerize(probability, actual_population);
        NumericVector sampled = f6ng_sample_counts(rounded, seed_cells);
        PutRNGstate();
        IntegerVector advanced = global[".Random.seed"];
        streams(_, c) = advanced;
        for (int r = 0; r < n; ++r) next(r, c) = sampled[r] / seed_cells;
        passage_count[c] += 1;
        if (first_passage_day[c] == NA_INTEGER)
          first_passage_day[c] = start_day + day;
        last_passage_day[c] = start_day + day;
        log_population[c] = log_seed;
      }
      double mean = 0, rate = 0;
      f6ng_state_summaries(next, c, ploidy_grid, net_live_rate, mean, rate);
      mean_ploidy(c, day) = mean;
      net_growth_rate(c, day) = rate;
    }
    std::copy(next.begin(), next.end(), state.begin());
    if (day % 100 == 0) checkUserInterrupt();
  }
  if (had_seed) global[".Random.seed"] = caller_seed;
  else global.remove(".Random.seed");
  return List::create(_["mean_ploidy"] = mean_ploidy,
    _["net_growth_rate"] = net_growth_rate, _["final_state"] = state,
    _["final_log_population"] = log_population, _["rng_final"] = streams,
    _["passage_count"] = passage_count,
    _["first_passage_day"] = first_passage_day,
    _["last_passage_day"] = last_passage_day);
}
