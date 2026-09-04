#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Systematic randomized rounding has integer total T and E[count_i]=T*p_i.
// One common uniform preserves rare states in expectation, unlike independent
// deterministic rounding. Only passage-day populations are integerized.
NumericVector f7s_integerize(const NumericVector& probability, double total) {
  if (!R_finite(total) || total < 1 || total > 9007199254740991.0)
    stop("Passage population outside exact-double integer range.");
  const double T = std::floor(total + 0.5), u = R::unif_rand();
  NumericVector count(probability.size());
  long double cumulative = 0;
  double previous = 0;
  for (int i = 0; i < probability.size(); ++i) {
    cumulative += static_cast<long double>(T) * probability[i];
    const double upper = i + 1 == probability.size() ? T :
      std::min(T, static_cast<double>(std::floor(cumulative + u)));
    count[i] = upper - previous;
    if (count[i] < 0) stop("Negative integerized state.");
    previous = upper;
  }
  return count;
}

NumericVector f7s_sample_counts(const NumericVector& counts, double size) {
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
  if (needed != 0 || sum(sample) != size) stop("Sampling total is not conserved.");
  return sample;
}

// [[Rcpp::export]]
List f7s_draw_cpp(NumericVector probability, double total, double size) {
  if (is_true(any(probability < 0)) || !R_finite(sum(probability)) || sum(probability) <= 0)
    stop("Invalid sampling probabilities.");
  probability = probability / sum(probability);
  NumericVector counts = f7s_integerize(probability, total);
  return List::create(_["population"] = counts,
    _["sample"] = f7s_sample_counts(counts, size));
}

// RNG is owned by trajectory, not worker. Explicit Get/PutRNGstate below lets
// batched BLAS share the deterministic operator without sharing random streams.
// A call can cover any integer-day block. Returned state/log size/RNG streams
// are sufficient to resume with exactly the same trajectory and random draws.
// [[Rcpp::export(rng = false)]]
List f7s_propagate_cpp(NumericMatrix step, NumericMatrix initial_state,
    NumericVector ploidy_grid, int n_day, double seed_cells,
    double target_cells, IntegerMatrix rng_states,
    NumericVector initial_log_population, int start_day = 0,
    bool stochastic = true, int trace_columns = 0) {
  const int n = step.nrow(), k = initial_state.ncol();
  if (step.ncol() != n || initial_state.nrow() != n || ploidy_grid.size() != n ||
      rng_states.nrow() != 7 || rng_states.ncol() != k ||
      initial_log_population.size() != k || n_day < 0 || start_day < 0 ||
      seed_cells < 1 || seed_cells != std::floor(seed_cells) ||
      target_cells <= seed_cells || target_cells != std::floor(target_cells))
    stop("Invalid stochastic propagation contract.");
  NumericMatrix state = clone(initial_state), next(n, k), means(k, n_day + 1);
  IntegerMatrix streams = clone(rng_states);
  NumericVector logpop = clone(initial_log_population), max_jump(k);
  IntegerVector count(k), first(k, NA_INTEGER), last(k, NA_INTEGER);
  std::vector<double> event_day, event_column, event_before, event_after;
  std::vector<double> event_mean_before, event_mean_after;
  const double logtarget = std::log(target_cells), logseed = std::log(seed_cells);
  for (int c = 0; c < k; ++c) {
    double total = 0;
    for (int r = 0; r < n; ++r) {
      if (!R_finite(state(r,c)) || state(r,c) < 0) stop("Invalid initial state.");
      total += state(r,c);
    }
    if (std::abs(total-1) > 1e-10 || !R_finite(logpop[c])) stop("Initial composition must sum to one.");
    for (int r = 0; r < n; ++r) {
      means(c,0) += state(r,c) * ploidy_grid[r];
    }
  }
  Environment global = Environment::global_env();
  const bool had_seed = global.exists(".Random.seed");
  RObject caller_seed = R_NilValue;
  if (had_seed) caller_seed = clone(as<IntegerVector>(global[".Random.seed"]));
  for (int d = 1; d <= n_day; ++d) {
    const char trans = 'N';
    const double alpha = 1, beta = 0;
    F77_CALL(dgemm)(&trans, &trans, &n, &k, &n, &alpha, step.begin(), &n,
      state.begin(), &n, &beta, next.begin(), &n FCONE FCONE);
    for (int c = 0; c < k; ++c) {
      double growth = 0, mean = 0;
      for (int r = 0; r < n; ++r) {
        double value = next(r,c);
        if (value < 0 && value > -1e-10) value = 0;
        if (!R_finite(value) || value < 0) stop("Invalid propagated state.");
        next(r,c) = value;
        growth += value;
      }
      if (!R_finite(growth) || growth <= 0) stop("Non-positive growth factor.");
      for (int r = 0; r < n; ++r) {
        next(r,c) /= growth;
        mean += next(r,c) * ploidy_grid[r];
      }
      logpop[c] += std::log(growth);
      if (logpop[c] >= logtarget) {
        const double actual_population = std::exp(logpop[c]);
        const double before = mean;
        if (stochastic) {
          IntegerVector stream = streams(_,c);
          global[".Random.seed"] = stream;
          GetRNGstate();
          NumericVector probability = next(_,c);
          NumericVector rounded = f7s_integerize(probability, actual_population);
          NumericVector sampled = f7s_sample_counts(rounded, seed_cells);
          PutRNGstate();
          IntegerVector advanced = global[".Random.seed"];
          streams(_,c) = advanced;
          mean = 0;
          for (int r = 0; r < n; ++r) {
            next(r,c) = sampled[r] / seed_cells;
            mean += next(r,c) * ploidy_grid[r];
          }
        }
        max_jump[c] = std::max(max_jump[c], std::abs(mean - before));
        count[c]++;
        if (first[c] == NA_INTEGER) first[c] = start_day + d;
        last[c] = start_day + d;
        logpop[c] = logseed;
        if (c < trace_columns) {
          event_day.push_back(start_day+d); event_column.push_back(c+1);
          event_before.push_back(actual_population); event_after.push_back(seed_cells);
          event_mean_before.push_back(before); event_mean_after.push_back(mean);
        }
      }
      means(c,d) = mean; // right-continuous, actual integer-day sampled population
    }
    std::copy(next.begin(), next.end(), state.begin());
    if (d % 100 == 0) checkUserInterrupt();
  }
  if (had_seed) global[".Random.seed"] = caller_seed;
  else global.remove(".Random.seed");
  return List::create(_["mean_ploidy"] = means, _["final_state"] = state,
    _["final_log_population"] = logpop, _["rng_final"] = streams,
    _["passage_count"] = count, _["first_passage_day"] = first,
    _["last_passage_day"] = last, _["maximum_boundary_mean_jump"] = max_jump,
    _["events"] = DataFrame::create(_["day"] = event_day, _["column"] = event_column,
      _["cells_before"] = event_before, _["cells_after"] = event_after,
      _["ploidy_before"] = event_mean_before, _["ploidy_after"] = event_mean_after));
}
