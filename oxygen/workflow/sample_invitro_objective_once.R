args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}
arg_flag <- function(name, default = FALSE) {
  raw <- arg_value(name, NULL)
  if (is.null(raw)) return(isTRUE(default))
  tolower(trimws(raw)) %in% c("1", "true", "t", "yes", "y")
}

seed <- suppressWarnings(as.integer(arg_value("seed", NA)))
if (!is.na(seed)) set.seed(seed)

maxit <- suppressWarnings(as.integer(arg_value("maxit", arg_value("maxiter", 50L))))
if (!is.finite(maxit) || maxit < 0L) stop("--maxit must be a non-negative integer.")

out_dir_arg <- arg_value("out_dir", NULL)
out_tsv_arg <- arg_value("out_tsv", NULL)
parameter_table_arg <- arg_value("parameter_table", NULL)
misseg_loss_survival_arg <- arg_value("misseg_loss_survival", "nullisomy")
write_rds <- arg_flag("write_rds", FALSE)
quiet <- arg_flag("quiet", TRUE)

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1]]) else "workflow/sample_invitro_objective_once.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(repo_root, "code", "in-vitro-utils", "io.R"))
source(file.path(repo_root, "code", "in-vitro-utils", "lineage_adapter.R"))
source(file.path(repo_root, "code", "in-vitro-utils", "runner.R"))
source(file.path(repo_root, "code", "in-vitro-utils", "summaries.R"))
source(file.path(repo_root, "code", "in-vitro-utils", "objective.R"))

ivt_source_map_model(repo_root)

cfg <- ivt_build_default_cfg(
  repo_root = repo_root,
  dt = 0.1,
  init_total_size = 1e6,
  o2_upper_bound = 21,
  fixed_oxygen = TRUE
)
cfg$misseg_loss_survival <- canonical_misseg_loss_survival_mode(
  misseg_loss_survival_arg,
  "nullisomy"
)
if (!is.null(parameter_table_arg) && nzchar(trimws(parameter_table_arg))) {
  cfg$parameter_table <- if (grepl("^/", parameter_table_arg)) {
    parameter_table_arg
  } else {
    file.path(repo_root, parameter_table_arg)
  }
} else {
  cfg$parameter_table <- ivt_parameter_table_for_loss_mode(
    repo_root = repo_root,
    loss_mode = cfg$misseg_loss_survival
  )
}
fit_objects <- ivt_load_fit_objects(repo_root)
optim_spec <- ivt_optimizer_spec(cfg)

sampled_par_t <- setNames(numeric(nrow(optim_spec)), optim_spec$param_name)
for (i in seq_len(nrow(optim_spec))) {
  sampled_par_t[[optim_spec$param_name[[i]]]] <- stats::runif(
    1,
    min = as.numeric(optim_spec$lower[[i]]),
    max = as.numeric(optim_spec$upper[[i]])
  )
}

objective_from_par <- function(par_t) {
  run_params_try <- ivt_optim_par_to_run_params(par_t = par_t, cfg = cfg)
  ivt_objective_components(
    run_params = run_params_try,
    fit_objects = fit_objects,
    cfg = cfg,
    fallback_max_passage_days = 14,
    growth_weight = 1,
    ploidy_weight = 1,
    flow_weight = 1
  )
}

initial_comp <- objective_from_par(sampled_par_t)

optim_fit <- stats::optim(
  par = sampled_par_t,
  fn = function(par_t) objective_from_par(par_t)$objective,
  method = "L-BFGS-B",
  lower = optim_spec$lower,
  upper = optim_spec$upper,
  control = list(maxit = maxit)
)

optimized_par_t <- setNames(as.numeric(optim_fit$par), optim_spec$param_name)
final_comp <- objective_from_par(optimized_par_t)
start_run_params <- ivt_optim_par_to_run_params(sampled_par_t, cfg = cfg)
optimized_run_params <- ivt_optim_par_to_run_params(optimized_par_t, cfg = cfg)

out_dir <- if (is.null(out_dir_arg) || !nzchar(out_dir_arg)) {
  file.path(repo_root, "workflow", "sampled_objective_draws")
} else {
  if (grepl("^/", out_dir_arg)) out_dir_arg else file.path(repo_root, out_dir_arg)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
draw_id <- paste0("draw_", stamp, if (!is.na(seed)) paste0("_seed", seed) else "")

summary_row <- data.frame(
  draw_id = draw_id,
  seed = if (is.na(seed)) NA_integer_ else seed,
  misseg_loss_survival = as.character(cfg$misseg_loss_survival),
  parameter_table = normalizePath(cfg$parameter_table, mustWork = FALSE),
  optimizer = "L-BFGS-B",
  maxit = as.integer(maxit),
  convergence = as.integer(optim_fit$convergence),
  message = if (is.null(optim_fit$message)) NA_character_ else as.character(optim_fit$message),
  initial_objective = as.numeric(initial_comp$objective),
  objective = as.numeric(final_comp$objective),
  total_loglik = as.numeric(final_comp$total_loglik),
  growth_loglik = as.numeric(final_comp$growth_loglik),
  ploidy_loglik = as.numeric(final_comp$ploidy_loglik),
  flow_loglik = as.numeric(final_comp$flow_loglik),
  growth_loglik_sum = as.numeric(final_comp$growth_loglik_sum),
  ploidy_loglik_sum = as.numeric(final_comp$ploidy_loglik_sum),
  flow_loglik_sum = as.numeric(final_comp$flow_loglik_sum),
  sigma_growth = as.numeric(final_comp$sigma_growth),
  sigma_kary = as.numeric(final_comp$sigma_kary),
  sigma_flow_ploidy = as.numeric(final_comp$sigma_flow_ploidy),
  n_growth = as.integer(final_comp$n_growth),
  n_growth_observed = as.integer(final_comp$n_growth_observed),
  n_growth_missing_pred = as.integer(final_comp$n_growth_missing_pred),
  n_growth_negative_pred = as.integer(final_comp$n_growth_negative_pred),
  n_ploidy_passages = as.integer(final_comp$n_ploidy_passages),
  n_kary_cells = as.integer(final_comp$n_kary_cells),
  n_flow_passages = as.integer(final_comp$n_flow_passages),
  n_flow_samples = as.integer(final_comp$n_flow_samples),
  fn_evals = as.integer(if ("function" %in% names(optim_fit$counts)) optim_fit$counts[["function"]] else NA_integer_),
  gr_evals = as.integer(if ("gradient" %in% names(optim_fit$counts)) optim_fit$counts[["gradient"]] else NA_integer_),
  stringsAsFactors = FALSE
)

for (nm in names(sampled_par_t)) summary_row[[paste0("start__", nm)]] <- as.numeric(sampled_par_t[[nm]])
for (nm in names(optimized_par_t)) summary_row[[nm]] <- as.numeric(optimized_par_t[[nm]])

natural_keep <- c(
  "lam_min", "lam_max", "k_o", "p_mis_base", "p_misseg", "k_o_mis",
  "gamma_loss", "buffer_smax", "buffer_beta", "buffer_n_exp",
  "p_wgd", "alpha_o2", "gamma_growth", "mu_hp",
  "gamma_mu", "O2_crit", "n_O", "sigma_growth", "sigma_kary",
  "init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N"
)
for (nm in natural_keep) summary_row[[paste0("start_nat__", nm)]] <- as.numeric(.first_non_null_local(start_run_params[[nm]], NA_real_))
for (nm in natural_keep) summary_row[[paste0("nat__", nm)]] <- as.numeric(.first_non_null_local(optimized_run_params[[nm]], NA_real_))

tsv_path <- file.path(out_dir, paste0(draw_id, ".tsv"))
rds_path <- file.path(out_dir, paste0(draw_id, ".rds"))
append_row_with_lock <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lock_dir <- paste0(path, ".lock")
  lock_ok <- FALSE
  start <- Sys.time()
  while (!lock_ok) {
    lock_ok <- dir.create(lock_dir, showWarnings = FALSE)
    if (!lock_ok) {
      elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
      if (is.finite(elapsed) && elapsed > 60) {
        stop("Timed out waiting for lock: ", lock_dir)
      }
      Sys.sleep(stats::runif(1L, min = 0.02, max = 0.10))
    }
  }
  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)

  write_header <- !file.exists(path) || isTRUE(file.info(path)$size == 0)
  write.table(
    df,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = write_header,
    append = !write_header
  )
}

if (!is.null(out_tsv_arg) && nzchar(out_tsv_arg)) {
  shared_tsv <- if (grepl("^/", out_tsv_arg)) out_tsv_arg else file.path(repo_root, out_tsv_arg)
  append_row_with_lock(summary_row, shared_tsv)
  if (!isTRUE(quiet)) cat("Appended summary:", shared_tsv, "\n")
} else {
  write.table(summary_row, file = tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
  if (!isTRUE(quiet)) cat("Wrote summary:", tsv_path, "\n")
}

if (isTRUE(write_rds)) {
  saveRDS(
    list(
      summary = summary_row,
      start_par_t = sampled_par_t,
      optimized_par_t = optimized_par_t,
      start_run_params = start_run_params,
      optimized_run_params = optimized_run_params,
      initial_objective_components = initial_comp[c("objective", "total_loglik", "growth_loglik", "ploidy_loglik", "flow_loglik", "growth_loglik_sum", "ploidy_loglik_sum", "flow_loglik_sum", "n_growth", "n_growth_observed", "n_growth_missing_pred", "n_growth_negative_pred", "n_ploidy_passages", "n_kary_cells", "n_flow_passages", "n_flow_samples", "sigma_growth", "sigma_kary", "sigma_flow_ploidy")],
      final_objective_components = final_comp[c("objective", "total_loglik", "growth_loglik", "ploidy_loglik", "flow_loglik", "growth_loglik_sum", "ploidy_loglik_sum", "flow_loglik_sum", "n_growth", "n_growth_observed", "n_growth_missing_pred", "n_growth_negative_pred", "n_ploidy_passages", "n_kary_cells", "n_flow_passages", "n_flow_samples", "sigma_growth", "sigma_kary", "sigma_flow_ploidy")]
    ),
    file = rds_path
  )
  if (!isTRUE(quiet)) cat("Wrote payload:", rds_path, "\n")
}

if (!isTRUE(quiet)) print(summary_row)
