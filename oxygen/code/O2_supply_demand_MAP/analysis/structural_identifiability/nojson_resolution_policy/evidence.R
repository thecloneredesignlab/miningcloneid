#!/usr/bin/env Rscript
# Fresh-process mapping and endpoint evidence; no optimizer is invoked.
args <- commandArgs(TRUE); stopifnot(length(args)==4L)
project <- normalizePath(args[1]); study <- normalizePath(args[2])
release <- normalizePath(args[3]); batch <- args[4]
modules <- file.path(study,"implementation",c("helpers.R","runtime.R"))
for(path in modules) source(path)
verify_snapshot(study)
manifest <- read_table(file.path(release,"source_manifest.tsv"))
stopifnot(all(vapply(file.path(release,manifest$path),o2prof_hash,character(1))==manifest$sha256))
settings <- read_kv(file.path(study,"settings.tsv"))
stopifnot(identical(Sys.getenv("O2PROF_IMAGE_SHA256"),settings$image_sha256),
          as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))==1L)
ctx <- o2prof_load_culture(file.path(project,settings$wp0_relative_path),modules)
ctx$backend$main <- function(...) stop("Optimization is forbidden in policy evidence")
request_path <- file.path(study,"policy/requests",paste0(batch,".tsv"))
requests <- read_table(request_path)
generation <- trimws(readLines(file.path(study,"state/CURRENT")))
stopifnot(length(generation)==1L,grepl("^g_[0-9a-f]+$",generation))
state_manifest <- read_table(file.path(study,"state",generation,"manifest.tsv"))
stopifnot(all(vapply(file.path(study,"state",generation,state_manifest$path),o2prof_hash,character(1))==state_manifest$sha256))
events <- read_table(file.path(study,"state",generation,"events.tsv"))
expected_request <- paste(batch,o2prof_hash(request_path),o2prof_hash(file.path(release,"source_manifest.tsv")))
stopifnot(expected_request %in% events$detail[events$event=="policy_evidence_request"])
out <- file.path(study,"policy/evidence",batch)
stopifnot(!dir.exists(out));dir.create(out,recursive=TRUE)
objective_calls <- 0L
raw_objective <- ctx$raw_objective
ctx$raw_objective <- function(...) {objective_calls <<- objective_calls+1L;raw_objective(...)}
results <- list()
for(i in seq_len(nrow(requests))) {
  r <- requests[i,];src <- file.path(project,r$source_fit)
  original <- file.path(study,"fits",r$task_id)
  stopifnot(identical(normalizePath(src),normalizePath(file.path(original,"accepted_fit.rds"))),
    o2prof_hash(src)==r$source_sha256,
    o2prof_hash(file.path(original,"result.tsv"))==r$source_result_sha256)
  receipt <- read_kv(file.path(original,"artifact_validation.tsv"))
  stopifnot(isTRUE(receipt$pass),receipt$accepted_fit_sha256==r$source_sha256,
    receipt$result_sha256==r$source_result_sha256,
    receipt$payload_manifest_sha256==o2prof_hash(file.path(original,"payload_manifest.tsv")))
  source_result <- read_kv(file.path(original,"result.tsv"))
  stopifnot(abs(source_result$objective-r$source_objective)<=1e-8)
  fit <- readRDS(src)
  task <- read_kv(file.path(study,"tasks",paste0(r$task_id,".tsv")))
  registry <- read_table(file.path(study,"cases.tsv"))
  fixed <- if(task$point=="base") list() else {
    z <- registry[registry$point==task$point,,drop=FALSE];stopifnot(nrow(z)==1L)
    stopifnot(identical(as.numeric(task$value),z$value[[1]]),identical(task$target,z$target[[1]]))
    setNames(list(task$value),task$target)
  }
  folder <- file.path(out,paste0(r$task_id,"_",r$action));dir.create(folder)
  local_pass <- FALSE
  if(r$action=="free_map") {
    map <- original_map(ctx)
    p <- o2prof_project_start(map,fit$best_params)$params
    stopifnot(sum(map$table$status=="active")==20L)
    comp <- replay_fit(ctx,map,fit,file.path(folder,"free_replay"))
    # Projection may reactivate formerly inactive coordinates, never alter a zero.
    for(n in names(fixed)) stopifnot(identical(p[[n]],fixed[[n]]))
    o2prof_to_free(map,p)
  } else if(r$action=="endpoint") {
    stopifnot(task$mode=="local",!isTRUE(source_result$local_stable))
    map <- original_map(ctx,fixed)
    local <- readRDS(file.path(original,"local_search.rds"))
    stopifnot(!isTRUE(local$budget_limited),local$evaluations<=20000L,
              identical(as.numeric(local$par),as.numeric(tail(local$history,1)[[1]])))
    for(h in c(.001,.0001)) {
      z <- local$sweeps[local$sweeps$mesh==h,,drop=FALSE]
      stopifnot(nrow(z)>0L);z <- tail(z,1)
      stopifnot(isTRUE(z$complete),z$maximum_tested_improvement<=1e-4,z$improvement<1e-5)
    }
    stopifnot(all(local$sweeps$complete),local$value==tail(local$sweeps$after,1))
    spec <- o2prof_spec(map,fit$best_params)
    stopifnot(length(local$par)==length(spec$lower),all(local$par>=0 & local$par<=1))
    p <- o2prof_from_free(map,setNames(spec$lower+local$par*(spec$upper-spec$lower),spec$param_name))
    near <- spec$param_name[local$par<=.001 | local$par>=.999]
    boundary_path <- file.path(original,"boundary_probes.tsv")
    if(length(near)) {
      b <- read_table(boundary_path)
      stopifnot(nrow(b)==3L*length(near),setequal(unique(b$coordinate),near),
                all(abs(b$endpoint-local$value)<=1e-8),all(local$value-b$probe<=1e-5))
      for(n in near) stopifnot(identical(sort(b$inward_fraction[b$coordinate==n]),c(.001,.01,.05)))
    } else stopifnot(!file.exists(boundary_path) || nrow(read_table(boundary_path))==0L)
    comp <- ctx$raw_objective(p,ctx$objects,ctx$cfg)
    o2prof_check_samples(comp,ctx$manifest);assert_predictions(comp)
    independent <- o2prof_independent_score(ctx,comp)
    stopifnot(abs(comp$objective-local$value)<=1e-8,abs(independent$objective-comp$objective)<=1e-8,
              comp$objective-source_result$objective<=1e-4,
              comp$objective-source_result$objective>=-1e-8)
    local_pass <- TRUE
    saveRDS(local$par,file.path(folder,"audited_endpoint_unit.rds"))
    saveRDS(local$sweeps,file.path(folder,"audited_sweeps.rds"))
  } else stop("Unknown policy evidence action")
  verified <- list(best_params=p,best_components=comp,cfg=ctx$cfg)
  path <- file.path(folder,"verified_fit.rds");saveRDS(verified,path)
  saveRDS(p,file.path(folder,"parameters.rds"));saveRDS(map,file.path(folder,"parameter_map.rds"))
  audit_native_strings(verified);audit_native_strings(map)
  # Re-evaluation enforces the original score/state/sample/day comparison policy.
  replay_fit(ctx,map,verified,file.path(folder,"verified_vector_replay"))
  if(i==1L && startsWith(batch,"policy_validation_")) {
    bad <- verified;bad$best_components$objective <- bad$best_components$objective+1e-6
    message <- tryCatch({replay_fit(ctx,map,bad,file.path(folder,"expected_score_rejection"));"NO_REJECTION"},
      error=function(e) conditionMessage(e))
    stopifnot(grepl("score_components",message,fixed=TRUE))
    o2prof_write_tsv(data.frame(pass=TRUE,injected_score_error=1e-6,
      retained_integrity_tolerance=1e-8,expected_error=gsub("[\t\r\n]"," ",message)),
      file.path(out,"integrity_negative_control.tsv"))
  }
  results[[i]] <- data.frame(policy_id="culture_resolution_aware_v1",action=r$action,
    task_id=r$task_id,source_fit=r$source_fit,source_sha256=r$source_sha256,
    source_result_sha256=r$source_result_sha256,source_objective=source_result$objective,
    objective=comp$objective,endpoint_gap=comp$objective-source_result$objective,
    local_evidence=local_pass,verified_fit=substring(path,nchar(project)+2L),
    verified_fit_sha256=o2prof_hash(path),parameters_sha256=o2prof_hash(file.path(folder,"parameters.rds")))
}
o2prof_write_tsv(do.call(rbind,results),file.path(out,"results.tsv"))
assert_runtime(out)
paths <- list.files(out,recursive=TRUE,full.names=TRUE);paths <- paths[!dir.exists(paths)]
o2prof_write_tsv(data.frame(path=substring(paths,nchar(out)+2L),sha256=vapply(paths,o2prof_hash,character(1))),file.path(out,"payload_manifest.tsv"))
write_kv(list(pass=TRUE,policy_sha256=o2prof_hash(file.path(release,"source_manifest.tsv")),
  requests_sha256=o2prof_hash(request_path),rows=nrow(requests),model_objective_calls=objective_calls,
  optimizer_calls=0L,source_fresh_process=TRUE,
  payload_manifest_sha256=o2prof_hash(file.path(out,"payload_manifest.tsv"))),file.path(out,"validation.tsv"))
writeLines("Verified policy mapping and vector-specific endpoint evidence",file.path(out,"COMPLETE"))
cat("Policy evidence passed:",nrow(requests),"rows;",objective_calls,"objective calls; zero optimizer calls.\n")
