# Population-count outputs derived from fitted in-vitro simulations.

ivt_sim_collect_population_tables <- function(components) {
  run_2n <- components$run_2N
  run_4n <- components$run_4N
  if (!is.list(run_2n) || !is.list(run_4n)) {
    stop("In-vitro best_components must contain run_2N and run_4N.")
  }
  list(
    invitro_lineage_summary = as.data.frame(components$summary),
    invitro_daily_counts = dplyr::bind_rows(
      ivt_collect_daily_counts(run_2n),
      ivt_collect_daily_counts(run_4n)
    )
  )
}
