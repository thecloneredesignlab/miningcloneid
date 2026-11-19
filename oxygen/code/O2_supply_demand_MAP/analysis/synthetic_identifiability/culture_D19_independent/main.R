main <- function() {
  args <- commandArgs(TRUE)
  study <- normalizePath(args[1L])
  task <- args[2L]
  out <- file.path(study, 'verification', task)
  for (name in c('production_bindings.R','interface.R','feasibility.R','search_helpers.R',
                 'core.R','joint_numerics.R','constraints.R','refine.R','engineering.R','fit.R')) {
    source(file.path(study, 'release', name))
  }
  on.exit(sf_inventory(out), add = TRUE)
  rows <- list()
  check <- function(name, ok, detail = '') {
    rows[[length(rows) + 1L]] <<- data.frame(check = name, passed = isTRUE(ok), detail = detail)
    sf_emit(do.call(rbind, rows), file.path(out, 'checks.tsv'))
    if (!isTRUE(ok)) stop('Check failed: ', name)
  }
  if (task == 'T00') {
    blind_setup(study, out, check)
  } else if (task == 'G') {
    blind_canaries(study, out, check)
  } else {
    env <- sf_backend(study, out, native = task != 'D_A_import' && !endsWith(task, '_import'))
    if (task == 'E') d_tests(study, out, env, check)
    else if (startsWith(task, 'D_')) d_fit(study, task, out, env, check)
    else stop('Unauthorized phase ', task)
  }
  sf_emit(data.frame(status = 'PASS', task = task, N20 = 0, profiles = 0, tumor = 0),
          file.path(out, 'application_result.tsv'))
  sf_emit(data.frame(status = 'PASS', task = task), file.path(out, 'gate.tsv'))
}
tryCatch(main(), error = function(e) {
  cat('D19_INDEPENDENT_PHASE_FAILED: ', conditionMessage(e), '\n', sep = '')
  quit(status = 1, save = 'no')
})
