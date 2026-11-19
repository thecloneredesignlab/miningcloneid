# Descriptive post-lock assessment with no simulator forwards.
main <- function() {
  args <- commandArgs(TRUE)
  study <- normalizePath(args[1L])
  out <- normalizePath(args[2L])
  release <- file.path(study,'release')
  for(f in c('production_bindings.R','interface.R','feasibility.R','search_helpers.R','n20_model.R'))
    source(file.path(release,f))
  on.exit(sf_inventory(out),add=TRUE)
  lock_path <- file.path(study,'reports/selection_lock.tsv')
  lock_sha <- digest::digest(file=lock_path,algo='sha256')
  lock <- sf_tab(lock_path)
  check <- function(x) stopifnot(isTRUE(x))
  check(nrow(lock)==1L&&file.exists(file.path(study,'controller/SELECTION_LOCKED')))
  quality <- sf_tab(file.path(study,'reports/optimizer_quality.tsv'))
  origins <- sf_tab(file.path(study,'reports/origin_verdicts.tsv'))
  case <- readRDS(file.path(study,'inputs/public/N20.rds'))
  truth <- readRDS(file.path(study,'inputs/private/truth.rds'))
  reference <- readRDS(file.path(study,'inputs/private/reference.rds'))$pred
  reference <- n_align_predictions(case,reference)
  env <- sf_backend(study,out,native=FALSE)
  grid <- 22:154
  kernel <- function(s) {
    k <- outer(grid,grid,function(i,j)dnorm(i,j,s))
    sweep(k,2,colSums(k),'/')
  }
  k0 <- kernel(truth$sigma_kary)
  s0 <- unlist(lapply(refs <- names(reference$states),function(id) {
    state <- reference$states[[id]]
    if (!is.matrix(state)||ncol(state)!=133L) return(integer())
    which(apply(state/pmax(rowSums(state),.Machine$double.xmin),2,max)>=1e-8)
  }),use.names=FALSE)
  check(length(s0)>0L)
  oxygen_values <- sort(unique(vapply(unlist(lapply(case$design$adapters,`[[`,'segments'),recursive=FALSE),
    function(segment)segment$oxygen_pct,numeric(1))))
  diagnostic_grid <- expand.grid(N=grid[sort(unique(s0))],oxygen_pct=oxygen_values)
  # The grid combines recorded oxygen levels and states encountered in the sealed reference.
  answer <- noise <- outputs <- responses <- roots <- list()
  for (i in seq_len(nrow(origins))) {
    row <- origins[i,]
    q <- quality[quality$origin==row$origin & quality$key==row$best_verified_key,,drop=FALSE]
    check(nrow(q)==1L)
    saved <- readRDS(file.path(study,q$state_file))
    pred <- saved$pred
    p <- c(case$design$fixed,saved$candidate$params[cg_roles])
    key <- saved$candidate$key
    check(key==q$key)
    for (j in seq_len(nrow(case$design$bounds))) {
      bound <- case$design$bounds[j,]
      name <- bound$parameter
      x <- p[[name]]
      t <- truth[[name]]
      scale <- function(v) env$o2prof_transform(v,bound)
      answer[[length(answer)+1L]] <- data.frame(origin=row$origin,key=key,parameter=name,
        estimate=x,generating=t,estimated_minus_generating=x-t,
        transformed_fraction_error=(scale(x)-scale(t))/(scale(bound$upper)-scale(bound$lower)),
        lower_hit=x==bound$lower,upper_hit=x==bound$upper,inactive=name%in%r_activity(p))
    }
    qg <- pred$growth-case$data$growth
    component <- function(a,b) sum(abs(a-b))/2
    tv <- vapply(names(pred$kary),function(id)component(pred$kary[[id]],reference$kary[[id]]),numeric(1))
    outputs[[i]] <- data.frame(origin=row$origin,key=key,Q=q$Q,
      growth_RMSE_to_noisy=sqrt(mean(qg^2)),growth_RMSE_to_generating=sqrt(mean((pred$growth-reference$growth)^2)),
      mean_prediction_TV_to_generating=mean(tv),max_prediction_TV_to_generating=max(tv),
      changed_scored_days=sum(pred$selection$day!=reference$selection$day,na.rm=TRUE),
      changed_parent_days=sum(pred$selection$parent_day!=reference$selection$parent_day,na.rm=TRUE))
    k <- kernel(p$sigma_kary)
    noise[[i]] <- data.frame(origin=row$origin,key=key,sigma_kary=p$sigma_kary,generating_sigma_kary=truth$sigma_kary,
      max_kernel_TV_to_generating=max(colSums(abs(k-k0))/2),sigma_growth=p$sigma_growth,
      generating_sigma_growth=truth$sigma_growth,init_sd_2N=p$init_sd_2N,init_sd_4N=p$init_sd_4N)
    for (cohort in c('2N','4N')) {
      x <- env$ivt_gaussian_initial_state(grid,p[[paste0('init_mean_',cohort)]],p[[paste0('init_sd_',cohort)]])
      t <- env$ivt_gaussian_initial_state(grid,truth[[paste0('init_mean_',cohort)]],truth[[paste0('init_sd_',cohort)]])
      roots[[length(roots)+1L]] <- data.frame(origin=row$origin,key=key,cohort=cohort,chromosomes=grid,prediction=x,generating=t)
    }
    d <- diagnostic_grid
    law <- function(b) list(missegregation=env$.pmisseg_of_O2(d$oxygen_pct,b,d$N),
                              death=env$.mu_eff_of_O2(d$oxygen_pct,b,d$N),
                              growth=env$.lambda_eff_of_O2(d$oxygen_pct,b,d$N,O2_growth=case$design$cfg$O2_growth),
                              one_error_survival=b$buffer_smax*exp(-b$buffer_beta*(44/d$N)^b$buffer_n_exp))
    fit <- law(p)
    target <- law(truth)
    responses[[i]] <- do.call(rbind,lapply(names(fit),function(name)data.frame(
      origin=row$origin,key=key,response=name,grid_rows=nrow(d),RMSE=sqrt(mean((fit[[name]]-target[[name]])^2)),
      maximum_error=max(abs(fit[[name]]-target[[name]])))))
  }
  sf_emit(do.call(rbind,answer),file.path(out,'all20_parameter_recovery.tsv'))
  sf_emit(do.call(rbind,noise),file.path(out,'sigma_kary_and_growth_diagnostics.tsv'))
  sf_emit(do.call(rbind,outputs),file.path(out,'prediction_recovery.tsv'))
  sf_emit(do.call(rbind,responses),file.path(out,'response_recovery.tsv'))
  sf_emit(do.call(rbind,roots),file.path(out,'root_distribution_comparison.tsv'))
  check(identical(lock_sha,digest::digest(file=lock_path,algo='sha256')))
  job <- Sys.getenv('SLURM_JOB_ID')
  check(grepl('^[0-9]+$',job))
  sf_emit(data.frame(status='PASS',job=job,selection_sha=lock_sha,full_forward_calls=0,
     profiles=0,tumor=0,repeated_noisy=0),file.path(out,'application_result.tsv'))
  sf_emit(data.frame(status='PASS',job=job),file.path(out,'gate.tsv'))
}
tryCatch(main(),error=function(e) {
  cat('N20_ASSESSMENT_FAILED: ',conditionMessage(e),'\n',sep='')
  quit(status=1L,save='no')
})
