#include <Rcpp.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>
#include <algorithm>
#include <cmath>
#include <limits>

using namespace Rcpp;

// Fitting-compatible segment selection. Each segment is simulated in full;
// positive local days compete by absolute live-cell-count distance to target.
// The selected composition seeds the next segment at its EXPERIMENTAL boundary,
// which advances by segment_days, not by the selected model-local day.
// [[Rcpp::export]]
Rcpp::List f7g_propagate_segments_cpp(
    Rcpp::NumericMatrix step, Rcpp::NumericMatrix initial_state,
    Rcpp::NumericVector ploidy_grid, int n_day,
    double log_seed, double log_target, int segment_days,
    bool keep_trace = false) {
  const int n_state = step.nrow(), n_initial = initial_state.ncol();
  if (step.ncol() != n_state || initial_state.nrow() != n_state ||
      ploidy_grid.size() != n_state || n_day < 0 || segment_days < 1 ||
      !R_finite(log_seed) || !R_finite(log_target)) {
    stop("Invalid segment propagation inputs.");
  }
  NumericMatrix state = clone(initial_state), next(n_state, n_initial);
  NumericMatrix selected(n_state, n_initial), means(n_initial, n_day + 1);
  NumericVector log_population(n_initial, log_seed), best_distance(n_initial);
  NumericVector selected_mean(n_initial), error(n_initial), boundary_jump(n_initial);
  NumericVector selected_day_sum(n_initial), selected_distance_sum(n_initial);
  IntegerVector best_day(n_initial), count(n_initial), early_count(n_initial);
  IntegerVector first_day(n_initial, NA_INTEGER), last_day(n_initial, NA_INTEGER);
  IntegerMatrix histogram(n_initial, segment_days);
  IntegerMatrix trace(n_initial, keep_trace ? n_day / segment_days : 0);
  for (int col = 0; col < n_initial; ++col) {
    double total = 0.0;
    for (int row = 0; row < n_state; ++row) total += state(row, col);
    if (!R_finite(total) || total <= 0) stop("Invalid initial population.");
    for (int row = 0; row < n_state; ++row) {
      state(row, col) /= total;
      means(col, 0) += state(row, col) * ploidy_grid[row];
    }
  }
  for (int day = 1; day <= n_day; ++day) {
    const int local_day = (day - 1) % segment_days + 1;
    const char trans = 'N';
    const double alpha = 1.0, beta = 0.0;
    F77_CALL(dgemm)(&trans, &trans, &n_state, &n_initial, &n_state,
      &alpha, step.begin(), &n_state, state.begin(), &n_state,
      &beta, next.begin(), &n_state FCONE FCONE);
    for (int col = 0; col < n_initial; ++col) {
      double growth = 0.0, mean = 0.0;
      for (int row = 0; row < n_state; ++row) {
        double value = next(row, col);
        if (value < 0 && value > -1e-10) value = 0;
        if (!R_finite(value) || value < 0) stop("Invalid propagated state.");
        next(row, col) = value;
        growth += value;
      }
      if (!R_finite(growth) || growth <= 0) stop("Non-positive growth.");
      for (int row = 0; row < n_state; ++row) {
        next(row, col) /= growth;
        mean += next(row, col) * ploidy_grid[row];
      }
      log_population[col] += std::log(growth);
      // Dividing all distances by the same positive target preserves argmin;
      // this is NOT nearest log(cell count), which would change the selector.
      const double distance = std::abs(std::expm1(log_population[col] - log_target));
      if (local_day == 1 || distance < best_distance[col]) {
        best_distance[col] = distance;
        best_day[col] = local_day;
        selected_mean[col] = mean;
        for (int row = 0; row < n_state; ++row) selected(row, col) = next(row, col);
      }
      means(col, day) = mean;
      if (local_day == segment_days) {
        double after = 0.0;
        for (int row = 0; row < n_state; ++row) {
          next(row, col) = selected(row, col);
          after += next(row, col) * ploidy_grid[row];
        }
        error[col] = std::max(error[col], std::abs(after - selected_mean[col]));
        boundary_jump[col] = std::max(boundary_jump[col], std::abs(after - mean));
        means(col, day) = after; // right-continuous reseeded state at the boundary
        log_population[col] = log_seed;
        count[col] += 1;
        if (first_day[col] == NA_INTEGER) first_day[col] = day;
        last_day[col] = day;
        selected_day_sum[col] += best_day[col];
        selected_distance_sum[col] += best_distance[col];
        early_count[col] += best_day[col] < segment_days;
        histogram(col, best_day[col] - 1) += 1;
        if (keep_trace) trace(col, count[col] - 1) = best_day[col];
      }
    }
    std::copy(next.begin(), next.end(), state.begin());
  }
  return List::create(
    Named("mean_ploidy") = means, Named("final_state") = state,
    Named("passage_count") = count, Named("first_passage_day") = first_day,
    Named("last_passage_day") = last_day, Named("no_crossing") = count == 0,
    Named("final_log_population") = log_population,
    Named("maximum_pre_post_mean_error") = error,
    Named("maximum_boundary_mean_jump") = boundary_jump,
    Named("selected_model_day_sum") = selected_day_sum,
    Named("selected_relative_target_distance_sum") = selected_distance_sum,
    Named("earlier_than_segment_end_count") = early_count,
    Named("selected_day_histogram") = histogram,
    Named("selected_day_trace") = trace);
}

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
