# Recompute the original production score on archived trajectories with zero propagations.
args <- commandArgs(TRUE)
stopifnot(length(args)==4L)
code <- normalizePath(args[1L])
archive <- normalizePath(args[2L])
case_file <- normalizePath(args[3L])
out <- normalizePath(args[4L])
invisible(lapply(c('production_bindings.R','interface.R'),function(f)source(file.path(archive,'release',f))))
source(file.path(code,'n20_model.R'))
source(file.path(code,'production_score.R'))
legacy <- parse(file.path(dirname(code),'culture_growth_karyotype_recovery','interface.R'))
native <- Filter(function(x) is.call(x) && as.character(x[[1]])=='<-' &&
  is.symbol(x[[2]]) && as.character(x[[2]])=='cg_score',as.list(legacy))
stopifnot(length(native)==1L,identical(body(cg_score),body(eval(native[[1L]][[3L]]))))
env <- sf_backend(archive,out,native=FALSE)
case <- readRDS(case_file)
t <- readRDS(file.path(archive,'inputs/private/truth.rds'))
r <- readRDS(file.path(archive,'inputs/private/reference.rds'))
pred <- n_align_predictions(case,r$pred)
counts <- lapply(case$data$kary,function(x)tabulate(match(x,22:154),133))
context <- list(case=case,hist=counts)
out_rows <- list()
production <- numeric(2)
kary <- numeric(2)
for (i in seq_len(2L)) {
  q <- t
  q$sigma_growth <- q$sigma_growth*c(1,1.5)[i]
  a <- cg_score(env,case$design,r$runs,q,pred,case$data)
  b <- n_loss(context,pred,q)
  independent <- n_independent_loss(context,pred,q)
  err <- abs(c(a$growth,a$kary,a$objective)-unlist(b[c('growth','kary','value')]))
  check <- max(err) <= 1e-8 && max(abs(unlist(b[c('growth','kary','value')])-independent))<=1e-8 &&
    a$growth_rows==114L&&a$kary_rows==12L&&sum(lengths(case$data$kary))==220L
  stopifnot(check)
  production[i] <- a$growth
  kary[i] <- a$kary
  out_rows[[i]] <- data.frame(fixture=c('archived_reference','archived_reference_scale_perturbation')[i],
    growth_rows=a$growth_rows,karyotype_samples=a$kary_rows,karyotype_cells=sum(lengths(case$data$kary)),
    max_component_error=max(err),forward_evaluations=0L,result='PASS')
}
stopifnot(abs(diff(production))>1e-8,abs(diff(kary))<=1e-10,
          is.function(env$ivt_collect_lineage_summary),is.function(env$ivt_growth_loglik_df),
          is.function(env$ivt_ploidy_loglik_df))
sf_emit(do.call(rbind,out_rows),file.path(out,'local_production_score_regression.tsv'))
cat('N20_PRODUCTION_SCORE_REGRESSION_PASS on archived runs with zero new forwards\n')
