source(file.path(dirname(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])),"common.R"))
stopifnot(length(args)==4L)
case <- args[3];replicate_id <- as.integer(args[4]);stopifnot(replicate_id%in%1:2)
out <- file.path(study,"results/replays",paste0(case,"_rep",replicate_id));dir.create(out,recursive=TRUE,showWarnings=FALSE)
load_defs();guard(out)
row <- tab(file.path(study,"cases.tsv"));row <- row[row$case==case,,drop=FALSE];stopifnot(nrow(row)==1L)
stopifnot(o2prof_hash(file.path(project,row$path))==row$sha256)
ctx <- o2prof_load_culture(wp0,file.path(archive,"implementation",c("helpers.R","runtime.R")))
# Allow one full forward call per fresh process. No fit entry is invoked.
calls <- new.env(parent=emptyenv());calls$n <- 0L
raw <- ctx$raw_objective
ctx$raw_objective <- function(...) {
  calls$n <- calls$n+1L
  if(calls$n>1L)stop("Full-forward call allowance exceeded")
  parent <- file.path(study,"results/full_forward_slots")
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  marker <- file.path(parent,paste0(case,"_rep",replicate_id))
  if(!dir.create(marker,showWarnings=FALSE))stop("Full-forward slot already consumed; never repeat an ambiguous call")
  kv(list(case=case,replicate=replicate_id,job_id=Sys.getenv("SLURM_JOB_ID"),status="consumed_before_forward_call"),file.path(marker,"receipt.tsv"))
  raw(...)
}
fixed <- if(row$target=="none")list() else setNames(list(as.numeric(row$value)),row$target)
map <- original_map(ctx,fixed)
fit <- readRDS(file.path(project,row$path))
projection <- o2prof_project_start(map,fit$best_params)
parameters <- projection$params
stopifnot(!any(projection$projection$reason=="projected_to_bound"))
active_rows <- map$table$status=="active"
stopifnot(all(projection$projection$reason[active_rows]%in%c("unchanged","roundoff_to_bound")))
emit(projection$projection,file.path(out,"requested_actual_parameters.tsv"))
comp <- replay_fit(ctx,map,fit,file.path(out,"comparison"))
stopifnot(calls$n==1L,comp$n_death_observations==0L,comp$death_weight==0)
saveRDS(comp,file.path(out,"replayed_components.rds"),compress=TRUE)
independent <- o2prof_independent_score(ctx,comp)
kv(list(pass=TRUE,case=case,replicate=replicate_id,source_sha256=row$sha256,objective=comp$objective,archived_objective=fit$best_components$objective,independent_objective=independent$objective,score_difference=comp$objective-fit$best_components$objective,model_propagations=calls$n,optimizer_calls=0L,synthetic_datasets=0L,random_draws=0L),file.path(out,"validation.tsv"))
# Re-evaluate actual saved selection and boundary interfaces, without propagation.
passage <- list()
for(cohort in c("2N","4N"))for(i in seq_along(comp[[paste0("run_",cohort)]]$segment_results)) {
  r <- comp[[paste0("run_",cohort)]]$segment_results[[i]];g <- comp[[paste0("run_",cohort)]]$grid_pre
  selected <- ctx$backend$ivt_select_segment_observation(r$sim,r$selection$required_cells,g,r$selection$target_live_cells,r$segment$obs_days_local)
  stopifnot(selected$selected_index==r$selection$selected_index,selected$selected_day==r$selection$selected_day)
  boundary <- r$parent_boundary_selection
  if(!is.null(boundary)) {
    parent <- comp[[paste0("run_",cohort)]]$segment_results[[r$segment$parent_segment_id]]
    checked <- ctx$backend$ivt_extract_passage_end_state(parent$sim,boundary$required_cells,g,parent$selection$target_live_cells,parent$segment$obs_days_local)
    stopifnot(isTRUE(all.equal(checked,boundary,tolerance=1e-12,check.attributes=TRUE)))
  }
  passage[[length(passage)+1L]]<-data.frame(cohort=cohort,segment_id=r$segment$segment_id,selected_day=selected$selected_day,nominal_day=r$segment$duration_days,boundary_present=!is.null(boundary),selection_pass=TRUE,boundary_pass=TRUE)
}
emit(do.call(rbind,passage),file.path(out,"passage_interface_checks.tsv"))
differences <- list()
for(cohort in c("2N","4N")) {
  A<-comp[[paste0("run_",cohort)]]$segment_results;B<-fit$best_components[[paste0("run_",cohort)]]$segment_results
  for(i in seq_along(A))for(n in intersect(names(A[[i]]$sim),names(B[[i]]$sim))) {
    a<-A[[i]]$sim[[n]];b<-B[[i]]$sim[[n]]
    if(!is.numeric(a)||!is.numeric(b)||length(a)!=length(b)||!length(a)||n%in%c("cache_g_build","cache_g_hit","cache_g_hysteresis"))next
    finite<-is.finite(a)&is.finite(b)
    if(!any(finite))next
    differences[[length(differences)+1L]]<-data.frame(cohort=cohort,segment_id=A[[i]]$segment$segment_id,quantity=n,elements=length(a),maximum_absolute_difference=max(abs(a[finite]-b[finite])),maximum_scaled_difference=max(abs(a[finite]-b[finite])/pmax(1,abs(b[finite]))))
  }
}
emit(do.call(rbind,differences),file.path(out,"prediction_state_differences.tsv"))
if(case=="production_seed144" && replicate_id==1L) {
  # Parameter wiring tests evaluate functions/kernels, never full trajectories.
  source(file.path(archive,"implementation/diagnostics.R"))
  p <- parameters;full <- original_map(ctx);spec <- o2prof_spec(full,p);free <- o2prof_to_free(full,p)
  signature <- function(pp) {
    v <- numeric()
    for(oxygen in c(0,.5)) {
      k <- pilotdiag_rates(ctx,pp,ctx$cfg,oxygen)
      v <- c(v,as.numeric(as.matrix(k$generator)),k$rates$lambda,k$rates$mu,k$rates$p_per_chromosome)
    }
    v
  }
  before <- signature(p);roles <- tab(file.path(contract,"parameter_roles.tsv"));wiring <- list()
  for(i in seq_along(free)) {
    a <- free;unit <- (free[i]-spec$lower[i])/(spec$upper[i]-spec$lower[i]);a[i]<-spec$lower[i]+(if(unit>.99)unit-.01 else unit+.01)*(spec$upper[i]-spec$lower[i])
    pp <- o2prof_from_free(full,a);normalized <- ctx$backend$normalize_run_params_common(pp,cfg=ctx$cfg)
    other <- setdiff(roles$parameter,roles$parameter[i]);tol <- 64*.Machine$double.eps*pmax(1,abs(vapply(other,function(n)pp[[n]],numeric(1))))
    stopifnot(all(abs(vapply(other,function(n)normalized[[n]]-pp[[n]],numeric(1)))<=tol),abs(normalized[[roles$parameter[i]]]-pp[[roles$parameter[i]]])<=64*.Machine$double.eps*max(1,abs(pp[[roles$parameter[i]]])))
    n <- roles$parameter[i];effect <- NA_real_;operator <- ""
    if(roles$scientific_role[i]=="mechanism") {effect <- max(abs(signature(normalized)-before));operator<-"R response and native C++ daughter/generator comparison"}
    if(n%in%c("init_mean_2N","init_sd_2N","init_mean_4N","init_sd_4N")) {
      cohort <- if(grepl("2N$",n))"2N" else "4N";grid<-22:154
      A <- ctx$backend$ivt_gaussian_initial_state(grid,p[[paste0("init_mean_",cohort)]],p[[paste0("init_sd_",cohort)]])
      B <- ctx$backend$ivt_gaussian_initial_state(grid,pp[[paste0("init_mean_",cohort)]],pp[[paste0("init_sd_",cohort)]])
      stopifnot(abs(sum(A)-1)<1e-12,abs(sum(B)-1)<1e-12);effect<-max(abs(A-B));operator<-"Normalized initial chromosome composition"
    }
    if(n=="sigma_growth") {
      g<-fit$best_components$growth_df;A<-ctx$backend$ivt_growth_loglik_df(g,p$sigma_growth);B<-ctx$backend$ivt_growth_loglik_df(g,pp$sigma_growth)
      effect<-max(abs(A$loglik-B$loglik));operator<-"Observed Gaussian growth score; no new dynamics"
    }
    if(n=="sigma_kary") {
      g<-22:154;probs<-fit$best_components$run_2N$segment_results[[1]]$selection$selected_frac
      A<-ctx$backend$ivt_smooth_kary_distribution(g,probs,p$sigma_kary);B<-ctx$backend$ivt_smooth_kary_distribution(g,probs,pp$sigma_kary)
      effect<-max(abs(A-B));operator<-"Chromosome observation kernel; near-identity plateau is not inactivity"
    }
    stopifnot(is.finite(effect))
    wiring[[i]]<-data.frame(parameter=n,scientific_role=roles$scientific_role[i],requested_value=pp[[n]],normalized_value=normalized[[n]],other_coordinates_preserved=TRUE,operator=operator,maximum_function_change=effect,interpretation=if(effect==0)"Numerical plateau; no identifiability conclusion" else "Value reaches its declared function; no recovery conclusion",pass=TRUE)
  }
  stopifnot(calls$n==1L);emit(do.call(rbind,wiring),file.path(out,"parameter_wiring_audit.tsv"))
}
record_environment(out)
writeLines("Fixed existing-fit replay and interface checks complete",file.path(out,"COMPLETE"))
