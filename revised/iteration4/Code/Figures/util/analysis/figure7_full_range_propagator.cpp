#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

// Propagate normalized live-cell composition with a one-day expm operator.
// Population size is tracked in log space solely to identify target-triggered
// passage days.  A passage is composition neutral, so reseeding changes the
// population-size clock but not the normalized mean-ploidy trajectory.

// [[Rcpp::export]]
Rcpp::List f7g_propagate_operator_cpp(
    Rcpp::NumericMatrix step,
    Rcpp::NumericMatrix initial_state,
    Rcpp::NumericVector ploidy_grid,
    int n_day,
    double log_seed,
    double log_target,
    bool track_passage) {
  const int n_state = step.nrow();
  const int n_initial = initial_state.ncol();
  if (step.ncol() != n_state || initial_state.nrow() != n_state ||
      ploidy_grid.size() != n_state || n_day < 0) {
    Rcpp::stop("Invalid dimensions for Figure 7 full-range propagation.");
  }
  if (track_passage &&
      (!R_finite(log_seed) || !R_finite(log_target) ||
       !(log_target > log_seed))) {
    Rcpp::stop("Canonical passage thresholds must be finite and increasing.");
  }

  Rcpp::NumericMatrix state = Rcpp::clone(initial_state);
  Rcpp::NumericMatrix next(n_state, n_initial);
  Rcpp::NumericMatrix mean_ploidy(n_initial, n_day + 1);
  Rcpp::IntegerVector passage_count(n_initial, 0);
  Rcpp::IntegerVector first_passage_day(n_initial, NA_INTEGER);
  Rcpp::IntegerVector last_passage_day(n_initial, NA_INTEGER);
  Rcpp::NumericVector final_log_population(n_initial, log_seed);
  Rcpp::NumericVector maximum_pre_post_mean_error(n_initial, 0.0);

  for (int column = 0; column < n_initial; ++column) {
    double total = 0.0;
    double weighted = 0.0;
    for (int row = 0; row < n_state; ++row) {
      total += state(row, column);
      weighted += state(row, column) * ploidy_grid[row];
    }
    if (!R_finite(total) || total <= 0.0) {
      Rcpp::stop("Empty initial state in Figure 7 propagation.");
    }
    mean_ploidy(column, 0) = weighted / total;
  }

  for (int day = 1; day <= n_day; ++day) {
    const char no_transpose = 'N';
    const double alpha = 1.0;
    const double beta = 0.0;
    F77_CALL(dgemm)(
      &no_transpose, &no_transpose,
      &n_state, &n_initial, &n_state,
      &alpha, step.begin(), &n_state,
      state.begin(), &n_state,
      &beta, next.begin(), &n_state FCONE FCONE
    );
    for (int column = 0; column < n_initial; ++column) {
      double growth = 0.0;
      for (int row = 0; row < n_state; ++row) {
        double value = next(row, column);
        if (value < 0.0 && value > -1e-10) value = 0.0;
        if (!R_finite(value) || value < 0.0) {
          Rcpp::stop("Non-finite or negative propagated state.");
        }
        next(row, column) = value;
        growth += value;
      }
      if (!R_finite(growth) || growth <= 0.0) {
        Rcpp::stop("Non-positive one-day population growth factor.");
      }

      double weighted = 0.0;
      for (int row = 0; row < n_state; ++row) {
        next(row, column) /= growth;
        weighted += next(row, column) * ploidy_grid[row];
      }
      mean_ploidy(column, day) = weighted;
      final_log_population[column] += std::log(growth);

      if (track_passage && final_log_population[column] >= log_target) {
        const double before = weighted;
        const double after = weighted;
        maximum_pre_post_mean_error[column] = std::max(
          maximum_pre_post_mean_error[column], std::abs(before - after)
        );
        passage_count[column] += 1;
        if (first_passage_day[column] == NA_INTEGER) {
          first_passage_day[column] = day;
        }
        last_passage_day[column] = day;
        final_log_population[column] = log_seed;
      }
    }
    std::copy(next.begin(), next.end(), state.begin());
  }

  Rcpp::LogicalVector no_crossing(n_initial);
  for (int column = 0; column < n_initial; ++column) {
    no_crossing[column] = passage_count[column] == 0;
  }

  return Rcpp::List::create(
    Rcpp::Named("mean_ploidy") = mean_ploidy,
    Rcpp::Named("passage_count") = passage_count,
    Rcpp::Named("first_passage_day") = first_passage_day,
    Rcpp::Named("last_passage_day") = last_passage_day,
    Rcpp::Named("no_crossing") = no_crossing,
    Rcpp::Named("final_log_population") = final_log_population,
    Rcpp::Named("maximum_pre_post_mean_error") =
      maximum_pre_post_mean_error
  );
}
