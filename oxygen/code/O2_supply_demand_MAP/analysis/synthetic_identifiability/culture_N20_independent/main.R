main <- function() {
  args <- commandArgs(TRUE)
  study <- normalizePath(args[1L])
  task <- args[2L]
  out <- file.path(study, 'verification', task)
  for (name in c('production_bindings.R','interface.R','feasibility.R','search_helpers.R',
                 'n20_model.R','production_score.R','core.R','engineering.R','fit.R')) {
    source(file.path(study, 'release', name))
  }
  on.exit(sf_inventory(out), add = TRUE)
  rows <- list()
  check <- function(name, ok, detail = '') {
    rows[[length(rows) + 1L]] <<- data.frame(check = name, passed = isTRUE(ok), detail = detail)
    sf_emit(do.call(rbind, rows), file.path(out, 'checks.tsv'))
    if (!isTRUE(ok)) stop('Check failed: ', name)
  }
  check('native_slurm_identity',grepl('^[0-9]+$',Sys.getenv('SLURM_JOB_ID')))
  if (task == 'T00') {
    n_setup(study, out, check)
  } else if (task == 'G') {
    n_canaries(study, out, check)
  } else {
    env <- sf_backend(study, out, native = task != 'N_A_import' && !endsWith(task, '_import'))
    if (task == 'E') {
      check('production_score_binding',is.function(cg_score) && is.function(cg_entries) &&
        is.function(env$ivt_collect_lineage_summary) && is.function(env$ivt_growth_loglik_df) &&
        is.function(env$ivt_ploidy_loglik_df))
      n_tests(study, out, env, check)
    }
    else if (startsWith(task, 'N_')) n_fit(study, task, out, env, check)
    else stop('Unauthorized phase ', task)
  }
  sf_emit(data.frame(status = 'PASS', task = task, job=Sys.getenv('SLURM_JOB_ID'), profiles = 0, tumor = 0),
          file.path(out, 'application_result.tsv'))
  sf_emit(data.frame(status = 'PASS', task = task), file.path(out, 'gate.tsv'))
}
tryCatch(main(), error = function(e) {
  cat('N20_INDEPENDENT_PHASE_FAILED: ', conditionMessage(e), '\n', sep = '')
  quit(status = 1, save = 'no')
})
