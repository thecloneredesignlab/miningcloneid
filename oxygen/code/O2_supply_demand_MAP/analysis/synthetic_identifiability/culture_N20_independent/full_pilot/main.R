main <- function() {
  args <- commandArgs(TRUE)
  study <- normalizePath(args[1L])
  task <- args[2L]
  stopifnot(grepl('^N_[A-F]_(import|bench|de1b|de[2-6]|local|pattern|verify[12]|replay[12][ab])$',task))
  out <- file.path(study,'verification',task)
  for (f in c('production_bindings.R','interface.R','feasibility.R','search_helpers.R',
              'n20_model.R','production_score.R','core.R','fit.R')) source(file.path(study,'release',f))
  on.exit(sf_inventory(out),add=TRUE)
  rows <- list()
  check <- function(label,condition,detail='') {
    rows[[length(rows)+1L]] <<- data.frame(check=label,passed=isTRUE(condition),detail=detail)
    sf_emit(do.call(rbind,rows),file.path(out,'checks.tsv'))
    if (!isTRUE(condition)) stop('Check failed: ',label,call.=FALSE)
  }
  job <- Sys.getenv('SLURM_JOB_ID')
  check('native_slurm_identity',grepl('^[0-9]+$',job))
  env <- sf_backend(study,out,native=!endsWith(task,'_import'))
  n_fit(study,task,out,env,check)
  sf_emit(data.frame(status='PASS',job=job,task=task,profiles=0,tumor=0,repeated_noisy=0),
    file.path(out,'application_result.tsv'))
  sf_emit(data.frame(status='PASS',job=job,task=task),file.path(out,'gate.tsv'))
}
tryCatch(main(), error=function(e) {
  cat('N20_FULL_PILOT_FAILED: ',conditionMessage(e),'\n',sep='')
  quit(status=1L,save='no')
})
