#!/usr/bin/env Rscript
# Deliberately separate acceptance entry: never launches a full campaign/draw.
.script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
.code <- dirname(.script);.ws <- normalizePath(file.path(.code,"../.."),mustWork=TRUE)
source(file.path(.code,"util/analysis/figure7_passaged.R"),local=globalenv())
source(file.path(.code,"test_Figure7_passaged.R"),local=globalenv())
.args <- commandArgs(TRUE)
.root_arg <- .args[grepl("^--audit-root=",.args)]
.resume <- "--resume=yes"%in%.args
if(length(.root_arg)!=1L||length(setdiff(.args,c(.root_arg,"--resume=yes"))))
  stop("Required: --audit-root=/absolute/audit/path [--resume=yes]")
.run <- normalizePath(sub("^--audit-root=","",.root_arg),mustWork=TRUE)
.audit_parent <- normalizePath(file.path(.ws,"audit"),mustWork=TRUE)
if(!startsWith(.run,paste0(.audit_parent,"/")))stop("Acceptance output must be an isolated workspace audit directory")
if(!file.exists(file.path(.run,"backup_receipt.json")))stop("Verified backup receipt is required first")
if(file.exists(file.path(.run,"execution_started.json"))&&!.resume)
  stop("Do not overwrite an existing execution; use validated --resume=yes for this same small suite")
.cached_index <- NULL;.previous_inputs <- .previous_locked <- NULL
if(.resume) {
  first <- jsonlite::fromJSON(file.path(.run,"execution_receipt.json"))
  stopifnot(first$timed_out,all(unlist(first$sources_unchanged)))
  # Entry/controller-only updates may add cache recovery or bounded probes. No
  # scientific implementation/config/test may drift while reusing results.
  for(path in c("Code/config/figure7_passaged.json","Code/Figures/util/analysis/figure7_passaged.cpp",
      "Code/Figures/util/analysis/figure7_passaged.R","Code/Figures/test_Figure7_passaged.R"))
    stopifnot(identical(digest::digest(file.path(.ws,path),file=TRUE,algo="sha256"),first$source_sha256[[path]]))
  .cached_index <- data.table::fread(file.path(.run,"results/case_index.tsv"))
  .previous_inputs <- data.table::fread(file.path(.run,"results/authoritative_inputs.tsv"))
  .previous_locked <- data.table::fread(file.path(.run,"results/locked_model_sources.tsv"))
}
cfg <- jsonlite::fromJSON(file.path(.ws,"Code/config/figure7_passaged.json"))
stopifnot(cfg$days==1000L,!cfg$full_campaign_authorized,!cfg$redraw_authorized,
  !cfg$historical_product_writes,!cfg$figure8_changes,!cfg$inoculum_feasibility_stop)
dir.create(file.path(.run,"results/cases"),recursive=TRUE,showWarnings=FALSE)
jsonlite::write_json(list(utc=format(Sys.time(),tz="UTC"),scope="Figure7-only small 1000-day acceptance",
  full_campaign=FALSE,redraw=FALSE,HPC=FALSE,resumed=.resume),
  file.path(.run,if(.resume)"execution_resume_started.json" else "execution_started.json"),auto_unbox=TRUE,pretty=TRUE)
runtime <- f7p_runtime(.ws,.run)
inputs <- lapply(setNames(c("C01","C02"),c("C01","C02")),function(f)f7p_load_fit(.ws,f))
if(.resume) {
  current <- data.table::rbindlist(lapply(inputs,`[[`,"manifest"))
  data.table::setorder(current,path);data.table::setorder(.previous_inputs,path)
  stopifnot(identical(current$path,.previous_inputs$path),identical(current$sha256,.previous_inputs$sha256),
    identical(runtime$locked$relative_path,.previous_locked$relative_path),
    identical(runtime$locked$sha256,.previous_locked$sha256))
}
f7p_write(data.table::rbindlist(lapply(inputs,`[[`,"parity")),file.path(.run,"results/complete_native_parameter_parity.tsv"))
f7p_write(data.table::rbindlist(lapply(inputs,`[[`,"manifest")),file.path(.run,"results/authoritative_inputs.tsv"))
f7p_unit_checks(runtime$native,inputs$C01,.run)
# Four factorial-diagonal representatives at all three candidate grids. Other
# combinations at 141N are coverage probes, NOT extra convergence certificates.
spec <- data.table::data.table(family=c("C01","C01","C02","C02"),
  initial=c(2,4,2,4),O2=c(.5,5,5,.5),pm=c(.1,.3,.3,.1))
spec[,`:=`(beta=0,threshold=66L,block=5L,role="grid_representative")]
grid_rows <- qc_rows <- native_rows <- prefix_rows <- route_rows <- list()
case_index <- summaries <- list();previous <- list();last_result <- NULL
run_one <- function(s,top) {
  id <- paste(s$family,s$initial,s$O2,s$pm,top,s$beta,s$threshold,s$block,sep="_")
  tic <- proc.time()[3];input <- inputs[[s$family]]
  op <- f7p_operator(input,s$O2,s$pm,top,s$beta,s$threshold)
  protocol <- f7p_protocol(input$fit$best_components$invitro[[paste0("run_",s$initial,"N")]]$adapter,s$block)
  prefix <- f7p_prefix(op,runtime$native,protocol,s$initial)
  files <- file.path(.run,"results/cases",paste0(id,"_",cfg$routes,".rds"))
  registered <- if(is.null(.cached_index))NULL else .cached_index[path%in%files]
  reused <- !is.null(registered)&&nrow(registered)>0L
  if(reused) {
    stopifnot(nrow(registered)==3L,setequal(registered$path,files))
    for(file in files)stopifnot(identical(digest::digest(file,file=TRUE,algo="sha256"),
      registered$sha256[match(file,registered$path)]))
    routes <- lapply(setNames(files,cfg$routes),readRDS)
    for(route in names(routes)) {
      r <- routes[[route]]
      stopifnot(identical(as.numeric(r$metrics$day),as.numeric(0:1000)),
        r$meta$family==s$family,r$meta$initial_ploidy==s$initial,r$meta$O2_pct==s$O2,
        r$meta$p_misseg==s$pm,r$meta$maximum_ploidy==top,r$meta$beta==s$beta,
        r$meta$threshold_N==s$threshold,r$meta$block==s$block,r$meta$route==route)
    }
  } else {
    # Preserve any unregistered interrupted writes; never mistake them for a
    # completed case or overwrite them without a recoverable local copy.
    for(file in files[file.exists(files)]) {
      dest <- file.path(.run,"interrupted_unregistered",basename(file))
      dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
      stopifnot(!file.exists(dest),file.rename(file,dest))
    }
    routes <- lapply(setNames(cfg$routes,cfg$routes),function(r)
      f7p_route(op,runtime$native,protocol,prefix,r,cfg$days,s$initial))
  }
  summary <- f7p_summary(routes);A <- summary[metric=="flux_percent"];B <- summary[metric=="growth_percent"]
  ok <- is.finite(A$median)
  complement <- f7p_max(abs(A$median[ok]+B$median[ok]-100))
  iqr_error <- max(f7p_max(abs(A$q25[ok]+B$q75[ok]-100)),f7p_max(abs(A$q75[ok]+B$q25[ok]-100)))
  ev <- data.table::rbindlist(lapply(routes,`[[`,"events"))
  pure <- fp_pure_initial(op$l,s$initial)
  pure_error <- max(vapply(routes,function(r)max(abs(r$X[,1]-pure)),0.0))
  row <- cbind(op$meta,initial_ploidy=s$initial,block=s$block,role=s$role,
    max_A_plus_B_error_pp=complement,max_IQR_complement_error_pp=iqr_error,
    pure_initial_error=pure_error,all_routes_reach_1000=all(vapply(routes,function(r)max(r$metrics$day)==1000,TRUE)),
    uniform_fallbacks=sum(ev$uniform_fallback),max_subset_error=max(ev$subset_error),
    roundoff_corrections=sum(ev$roundoff_corrections),
    stock_expansions_recorded_not_stopped=sum(ev$insufficient_stock_recorded_only),
    max_reseeding_scale=f7p_max(ev$reseeding_scale),
    minimum_expected_cells=min(vapply(routes,function(r)min(exp(r$metrics$log_cells)),0.0)),
    maximum_upper_loss_fraction=max(vapply(routes,function(r)f7p_max(r$metrics$upper_loss_fraction),0.0)),
    passed=complement<1e-10&&iqr_error<1e-10&&pure_error==0&&
      !any(ev$uniform_fallback)&&max(ev$subset_error)<1e-10)
  stopifnot(row$passed)
  qc_rows[[length(qc_rows)+1L]] <<- row
  for(r in names(routes)) {
    file <- file.path(.run,"results/cases",paste0(id,"_",r,".rds"))
    if(!reused)saveRDS(routes[[r]],file,compress=FALSE)
    case_index[[length(case_index)+1L]] <<- cbind(routes[[r]]$meta,role=s$role,
      path=file,sha256=digest::digest(file,file=TRUE,algo="sha256"))
    route_rows[[length(route_rows)+1L]] <<- cbind(routes[[r]]$meta,
      route_end=max(routes[[r]]$metrics$day),prefix_end=routes[[r]]$prefix_end,
      passage_segments=nrow(routes[[r]]$events),extension_segments=sum(routes[[r]]$events$phase=="extrapolation"))
    key <- paste(s$family,s$initial,s$O2,s$pm,s$beta,s$threshold,s$block,r,sep="_")
    if(s$role=="grid_representative"&&!is.null(previous[[key]])) {
      old <- readRDS(previous[[key]])
      grid_rows[[length(grid_rows)+1L]] <<- f7p_compare(old,routes[[r]],cfg$grid_limits)
      rm(old)
    }
    if(s$role=="grid_representative")previous[[key]] <<- file
  }
  # Native kernels only support unblocked dynamics. Check actual late states,
  # not a second run from an unrelated default/root distribution.
  if(top==141L&&s$beta==0&&s$role%in%c("grid_representative","oxygen_p_crossing_probe"))
    native_rows[[length(native_rows)+1L]] <<- f7p_native_checkpoint_checks(op,runtime$native,routes)
  prefix_rows[[length(prefix_rows)+1L]] <<- cbind(op$meta,initial_ploidy=s$initial,block=s$block,
    data.table::rbindlist(lapply(prefix,`[[`,"row")))
  summaries[[length(summaries)+1L]] <<- summary
  f7p_write(data.table::rbindlist(qc_rows),file.path(.run,"results/case_acceptance.tsv"))
  f7p_write(data.table::rbindlist(case_index),file.path(.run,"results/case_index.tsv"))
  if(length(grid_rows))f7p_write(data.table::rbindlist(grid_rows),file.path(.run,"results/grid_comparisons.tsv"))
  cat(if(reused)"F7_1000_CASE_REUSED" else "F7_1000_CASE",id,s$role,"routes",length(routes),
    "seconds",round(proc.time()[3]-tic,2),"\n");flush.console()
  invisible(NULL)
}
for(top in cfg$candidate_grid_upper_N)for(i in seq_len(nrow(spec))) {
  run_one(spec[i],as.integer(top));gc(verbose=FALSE)
}
probes <- data.table::data.table(family=c("C01","C01","C02","C02","C01","C01","C02","C01","C01","C01","C02","C02"),
  initial=c(2,4,2,4,4,4,2,4,4,4,2,2),
  O2=c(5,.5,.5,5,.5,.5,0,20,5,5,.5,.5),
  pm=c(.3,.1,.1,.3,.3,.3,.1,.1,.3,.3,.1,.1),
  beta=c(0,0,0,0,.5,1,0,0,0,0,0,0),
  threshold=c(66L,66L,66L,66L,66L,66L,44L,88L,66L,66L,66L,66L),
  block=c(5L,5L,5L,5L,5L,5L,5L,5L,3L,10L,3L,10L),
  role=c(rep("factorial_coverage",4),rep("blocking_probe",2),rep("oxygen_threshold_endpoint",2),rep("template_sensitivity",4)))
probes <- data.table::rbindlist(list(probes,data.table::data.table(
  family=c("C01","C02"),initial=c(4,2),O2=c(.5,5),pm=c(.3,.1),beta=0,
  threshold=66L,block=5L,role="oxygen_p_crossing_probe")))
for(i in seq_len(nrow(probes))){run_one(probes[i],141L);gc(verbose=FALSE)}
f7p_write(data.table::rbindlist(summaries),file.path(.run,"results/instant_time_summary.tsv"))
f7p_write(data.table::rbindlist(route_rows),file.path(.run,"results/route_inventory.tsv"))
f7p_write(data.table::rbindlist(prefix_rows),file.path(.run,"results/experimental_prefix_segments.tsv"))
f7p_write(data.table::rbindlist(native_rows),file.path(.run,"results/native_checkpoint_validation.tsv"))
if(!all(data.table::rbindlist(native_rows)$passed))stop("Native checkpoint parity failed")
g <- data.table::rbindlist(grid_rows)
status <- list(status=if(all(g$passed))"representative_1000day_acceptance_passed" else
    "implementation_checks_passed_but_1000day_grid_needs_revision",
  candidate_grids=cfg$candidate_grid_upper_N,grid_comparisons=nrow(g),
  grid_comparisons_passed=sum(g$passed),all_tested_grid_comparisons_passed=all(g$passed),
  numerical_cases=length(qc_rows),protocol_routes=length(case_index),
  full_campaign=FALSE,redraw=FALSE,HPC=FALSE,formal_data_promoted=FALSE,
  Figure8_modified=FALSE,historical_products_modified=FALSE)
jsonlite::write_json(status,file.path(.run,"acceptance_summary.json"),auto_unbox=TRUE,pretty=TRUE)
cat("FIGURE7_SMALL_ACCEPTANCE_COMPLETE",status$status,"\n")
