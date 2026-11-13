# Formal execution helpers. The simulator and objective remain frozen.
verify_payload_R <- function(root) {
 stopifnot(file.exists(file.path(root,"COMPLETE")),
  identical(readLines(file.path(root,"COMPLETE"),warn=FALSE),hash_file(file.path(root,"payload_manifest.tsv"))))
 m <- read_tsv(file.path(root,"payload_manifest.tsv"))
 stopifnot(all(vapply(file.path(root,m$path),hash_file,character(1))==m$sha256))
 invisible(TRUE)
}
formal_context <- function(study,spec) {
 ctx <- load_contract(study,TRUE,spec$policy)
 stopifnot(identical(ctx$evaluation_mode,"fresh_fork"))
 ctx$provenance <- getOption("ivi.provenance",list())
 ctx$provenance$task <- spec$task;ctx$provenance$role <- spec$role
 ctx$provenance$family <- spec$family;ctx$provenance$initialization <- spec$init_mode
 ctx$provenance$seed <- spec$seed;ctx$provenance$density_index <- spec$density_index
 ctx$provenance$root_id <- spec$root_id;ctx$provenance$root_family <- spec$root_family
 ctx
}
resolve_start <- function(ctx,spec) {
 record<-NULL;expected<-NULL;origin<-'';source_hash<-''
 if(spec$selector=='reference') {p<-ctx$endpoint;origin<-'canonical_reference'}
 else if(spec$selector=='import') {
  path<-file.path(ctx$study,'reference/imports',paste0(spec$import_id,'.rds'))
  z<-readRDS(path);stopifnot(z$target_policy==spec$policy,identical(z$root_id,spec$root_id))
  original<-z$original_record;p<-original$parameters;record<-original;expected<-z$expected_objective
  source_hash<-hash_file(path);origin<-paste0('import:',spec$import_id)
 } else if(spec$selector=='task') {
  parent<-readRDS(file.path(ctx$study,'manifests/task_specs.rds'))[[spec$source_task]]
  root<-file.path(ctx$study,if(parent$mode%in%c('global','local')) 'fits' else 'tests',spec$source_task);verify_payload_R(root)
  a<-readRDS(file.path(root,'accepted_candidate.rds'));record<-a$record;p<-record$parameters;expected<-a$checked$L
  stopifnot(a$policy==spec$policy,identical(record$source_manifest_sha256,hash_file(file.path(ctx$study,'manifests/input_manifest.tsv'))))
  if(!startsWith(spec$role,'density_') && !startsWith(spec$role,'reserve_density_'))stopifnot(readRDS(file.path(root,'dispatch_spec.rds'))$root_id==spec$root_id)
  source_hash<-hash_file(file.path(root,'payload_manifest.tsv'));origin<-spec$source_task
 } else {
  path<-file.path(ctx$study,'dispatch',getOption('ivi.root.task'),'selection.rds');z<-readRDS(path)
  stopifnot(z$policy==spec$policy,identical(z$record$source_manifest_sha256,hash_file(file.path(ctx$study,'manifests/input_manifest.tsv'))))
  if(!startsWith(spec$role,'density_') && !startsWith(spec$role,'reserve_density_'))stopifnot(z$root_id==spec$root_id)
  record<-z$record;p<-record$parameters;expected<-record$components$L;source_hash<-hash_file(path);origin<-z$source_task
 }
 stopifnot(identical(names(p),names(ctx$endpoint)),all(is.finite(p)),all(p>=ctx$lower),all(p<=ctx$upper))
 list(parameters=p,original_parameters=p,projection=p-p,source=origin,source_sha256=source_hash,record=record,expected=expected,root_id=spec$root_id)
}
validate_point <- function(ctx,p,expected=NULL,record=NULL) {
 z <- score(ctx,p);pred <- predictions(ctx,p,TRUE);ind <- independent_components(ctx,p,pred)
 stopifnot(max(abs(ind-unlist(z[names(ind)])))<1e-8)
 if(!is.null(expected)) stopifnot(abs(z$L-expected)<1e-8)
 if(!is.null(record)) stopifnot(abs(z$L-record$components$L)<1e-8,
  isTRUE(all.equal(pred$tables,record$predictions$tables,tolerance=1e-8)))
 list(components=z,predictions=pred,independent=ind)
}
make_direction_probes <- function(p,lower,upper,fixed,steps) {
 out <- list();free <- setdiff(names(p),names(fixed));span <- upper-lower
 for(step in steps) for(name in free) for(direction in c(-1L,1L)) {
  q <- p;q[name] <- p[name]+direction*step*span[name]
  if(q[name]<lower[name] || q[name]>upper[name]) next
  out[[length(out)+1L]] <- list(parameter=name,step=step,direction=direction,parameters=q)
 }
 out
}
run_acceptance <- function(study,task,spec,do_directions=FALSE) {
 out <- file.path(study,"tests",task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 ctx <- formal_context(study,spec);start <- resolve_start(ctx,spec);p <- start$parameters
 if(length(spec$fixed)) p[names(spec$fixed)] <- spec$fixed
 atomic_rds(start,file.path(out,"actual_start.rds"))
 ctx$provenance$start_sha256 <- hash_file(file.path(out,"actual_start.rds"))
 # A projected neighbor is a new local start; an acceptance point must already match.
 point <- validate_point(ctx,p,expected=start$expected,record=if(spec$selector=="import") NULL else start$record)
 if(spec$selector=="import") stopifnot(isTRUE(all.equal(point$predictions$tables,start$record$predictions$tables,tolerance=1e-8)))
 stopifnot(all(p[names(spec$fixed)]==spec$fixed))
 current_record <- if(spec$selector=="import") list(parameters=p,components=point$components,predictions=point$predictions,policy=spec$policy,free_names=setdiff(names(p),names(spec$fixed)),fixed=spec$fixed,contract_sha256=hash_file(file.path(study,"manifests/contract.rds")),source_manifest_sha256=hash_file(file.path(study,"manifests/input_manifest.tsv")),evaluation_mode="fresh_fork",provenance=ctx$provenance) else start$record
 accepted <- list(record=current_record,checked=point$components,fixed=spec$fixed,policy=spec$policy,
  predictions_full=point$predictions,source_task=start$source,provenance=ctx$provenance)
 if(do_directions) {
  free <- setdiff(names(p),names(spec$fixed))
  evaluator <- new_evaluator(ctx,out,spec$fixed,"coordinate_probe",getOption("ivi.job.deadline"))
  invisible(evaluator(p[free]));probes <- list();n <- 0L
  protocol <- readRDS(file.path(study,"manifests/protocol.rds"))
  for(probe in make_direction_probes(p,ctx$lower,ctx$upper,spec$fixed,protocol$steps)) {
   value <- evaluator(probe$parameters[free]);n <- n+1L
   probes[[n]] <- c(probe,list(objective=value,components=.ivi_state$last_components))
  }
  atomic_rds(probes,file.path(out,"directional_probes.rds"))
  best <- incumbents(out,ctx);verified <- validate_point(ctx,best$parameters,record=best)
  accepted$record <- best;accepted$checked <- verified$components;accepted$predictions_full <- verified$predictions
  improvement <- point$components$L-verified$components$L
  write_tsv(do.call(rbind,lapply(probes,function(z) data.frame(parameter=z$parameter,step=z$step,direction=z$direction,objective=z$objective,
   data=z$components$L_data,prior=z$components$L_prior,burden=z$components$L_b,ploidy=z$components$L_p,necrosis=z$components$L_n,
   improvement=point$components$L-z$objective))),file.path(out,"directional_scores.tsv"))
 } else { n <- 0L;improvement <- 0 }
 atomic_rds(accepted,file.path(out,"accepted_candidate.rds"))
 atomic_rds(list(parameters=p,components=point$components,predictions=point$predictions,source_task=start$source,policy=spec$policy),file.path(out,"fresh_acceptance.rds"))
 write_kv(c(list(pass=TRUE,policy=spec$policy,source_task=start$source,objective=point$components$L,
  difference=if(is.null(start$expected)) NA_real_ else point$components$L-start$expected,
  fixed_residual=if(length(spec$fixed)) max(abs(p[names(spec$fixed)]-spec$fixed)) else 0,
  probes=n,directional_improvement=improvement),component_values(point$components)),file.path(out,"result.tsv"))
 finalize_payload(out)
}
run_reference <- function(study,task,spec) {
 out <- file.path(study,"tests",task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 ctx <- formal_context(study,spec);old <- readRDS(file.path(study,"reference/production_replay.rds"))
 rows <- list()
 for(seed in c(25L,74L,366L,292L,392L)) {
  s <- spec;s$start_seed <- seed;s$selector <- "seed"
  p <- resolve_start(ctx,s)$parameters
  v <- validate_point(ctx,p)
  if(seed==25L) {
   expected <- if(spec$policy=="production") old$components$L else old$components$L_data
   stopifnot(abs(v$components$L-expected)<1e-8,
    isTRUE(all.equal(v$predictions$tables,old$predictions$tables,tolerance=1e-8)))
  }
  atomic_rds(list(parameters=p,components=v$components,predictions=v$predictions,policy=spec$policy),file.path(out,paste0("seed",seed,"_reference.rds")))
  rows[[length(rows)+1L]] <- data.frame(seed=seed,objective=v$components$L,data=v$components$L_data,prior=v$components$L_prior)
 }
 for(name in c("log10_o2_S0","log10_kappa_O","log10_eta_o2")) {
  before <- score(ctx,ctx$endpoint);q <- ctx$endpoint;q[name] <- min(ctx$upper[name],q[name]+.002)
  invisible(score(ctx,q));stopifnot(identical(component_values(before),component_values(score(ctx,ctx$endpoint))))
 }
 for(fixed in c(FALSE,TRUE)) {
  s <- spec;s$selector <- "seed";s$start_seed <- 25L
  s$fixed <- if(fixed) setNames(ctx$rho_grid[4],"log10_rho_2N") else NULL
  s$task <- paste0(task,"_direction_fixture_",if(fixed) "fixed" else "free")
  run_acceptance(study,s$task,s,TRUE)
 }
 write_tsv(do.call(rbind,rows),file.path(out,"reference_scores.tsv"))
 write_kv(list(pass=TRUE,policy=spec$policy,references=5L,resource_history_invariance=TRUE),file.path(out,"result.tsv"));finalize_payload(out)
}
