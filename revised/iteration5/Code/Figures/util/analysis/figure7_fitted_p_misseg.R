# Additional fitted-p condition. All biological kernels and native passage
# functions are reused unchanged; fixed-p receipts and their sources stay intact.
.f7fit_dir <- dirname(normalizePath(sys.frame(1)$ofile,mustWork=TRUE))
for(.f7fit_module in c("figure78_pipeline.R","figure78_model.R","figure7_simulation.R",
    "figure7_current_pipeline.R","figure7_population_validity.R"))
  source(file.path(.f7fit_dir,.f7fit_module),local=globalenv())

f7fit_args <- function(x=commandArgs(TRUE)) {
  z <- list(stage="plan",phase="full")
  for(v in x) {
    if(!grepl("^--[^=]+=.+$",v))stop("Use --option=value")
    at <- regexpr("=",v,fixed=TRUE)
    z[[gsub("-","_",substring(v,3,at-1))]] <- substring(v,at+1)
  }
  known <- c("stage","phase","workspace","campaign_root","fit_root","task_ids","task_index",
    "output_root","products","condition_mode")
  if(length(setdiff(names(z),known)))stop("Unknown fitted-p option")
  if(!is.null(z$condition_mode)&&z$condition_mode!="fitted_p_misseg")stop("Wrong condition mode")
  z
}
f7fit_init <- function(script,args) {
  suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite);library(Matrix);library(Rcpp)})
  data.table::setDTthreads(2L)
  ws <- normalizePath(args$workspace %||% file.path(dirname(script),"../.."),mustWork=TRUE)
  cfg_path <- file.path(ws,"Code/config/figure7_fitted_p_misseg.json")
  policy <- jsonlite::fromJSON(cfg_path)
  stopifnot(policy$condition_mode=="fitted_p_misseg",policy$days==1000,
    !policy$inoculum_feasibility_stop,!policy$existing_product_writes,!policy$absolute_count_filter,
    policy$relative_mass_min==1e-8,policy$template_segments==5)
  shared <- jsonlite::fromJSON(file.path(ws,"Code/config/figure78_inputs.json"))
  run <- f78_current_path(args$campaign_root %||% "data/Figures/Figure7/fitted_p_misseg",ws)
  protected <- file.path(ws,"data/Figures",c("Figure7","Figure8"))
  if(run%in%protected)stop("Refusing an existing Figure7/8 data root")
  a <- list(figure=7L,workspace=ws,code=file.path(ws,"Code/Figures"),shared=shared,
    config_path=cfg_path,policy=policy,cfg=jsonlite::fromJSON(file.path(ws,"Code/config/figure7.json")),
    model_root=f78_model_root(ws),fit_root=args$fit_root %||% file.path(ws,shared$fit_root),
    run=run,data_root=dirname(run),profile="full",args=args,
    output_root=args$output_root %||% file.path(run,"rendered"))
  a$cfg$day <- policy$days;a$cfg$main_oxygen <- policy$oxygen
  a$cfg$intervention_q <- policy$beta;a$cfg$snapshot_days <- policy$snapshot_days
  f78_locked_model_manifest(ws)
  a
}
f7fit_sources <- function(a) {
  unique(c(fp_files(a),file.path(a$code,c("run_Figure7_fitted_p_misseg.R",
    "util/analysis/figure7_fitted_p_misseg.R","util/analysis/figure7_passaged.R",
    "util/analysis/figure7_passaged.cpp","util/analysis/figure7_current_pipeline.R",
    "util/analysis/figure7_population_validity.R","test_Figure7_passaged.R")),
    a$config_path,file.path(a$workspace,"Code/config",c("figure78_inputs.json",
      "manifests/figure78_model_sources.tsv"))))
}
f7fit_prepared_sources <- function(a,b) {
  slots <- sub("^.*(/Code/)","Code/",b$code_manifest$path)
  paths <- file.path(a$workspace,slots)
  if(!identical(fp_sha(paths),b$code_manifest$sha256))stop("Fitted-p numerical sources changed")
  stopifnot(identical(b$policy,a$policy),b$condition_mode=="fitted_p_misseg")
  invisible(TRUE)
}
f7fit_bundle <- function(a) {
  b <- readRDS(file.path(a$run,"prepared/analysis_bundle.rds"))
  f7fit_prepared_sources(a,b);b
}
f7fit_prepare <- function(a) {
  dest <- file.path(a$run,"prepared/analysis_bundle.rds")
  if(file.exists(dest))return(invisible(f7fit_bundle(a)))
  if(dir.exists(a$run)&&length(setdiff(list.files(a$run,all.files=TRUE,no..=TRUE),"runtime")))
    stop("Nonempty unprepared campaign; only verified compiler caches may pre-exist")
  fp_preflight(a)
  dir.create(a$run,recursive=TRUE)
  runtime <- f7p_runtime(a$workspace,file.path(a$run,"runtime/preparation"))
  pars <- fp_load_parameters(a,fp_endpoint_spec(a))
  fits <- lapply(setNames(c("C01","C02"),c("C01","C02")),function(f)f7p_load_fit(a$workspace,f,a$fit_root))
  protocols <- list()
  for(f in names(fits))for(start in a$policy$initial_ploidy)
    protocols[[paste(f,start,sep="_")]] <- f7p_protocol(
      fits[[f]]$fit$best_components$invitro[[paste0("run_",start,"N")]]$adapter,5L)
  tasks <- list();id <- 0L
  for(i in seq_len(nrow(pars$endpoints))) {
    e <- pars$endpoints[i];pr <- pars$prepared[[e$prepared_index]]
    pm <- as.numeric(pr$native_parameters["p_misseg"])
    stopifnot(identical(pm,as.numeric(pr$run_params$p_misseg)))
    grids <- if(e$model_context=="in vitro")a$policy$in_vitro_grid_N else a$policy$in_vivo_grid_N
    for(o2 in a$policy$oxygen)for(start in a$policy$initial_ploidy)for(top in grids) {
      id <- id+1L
      # Four diagonal representatives cover both families, contexts, starts and oxygens.
      small <- (e$pair_label=="C01"&&o2==.5&&start==2)||
        (e$pair_label=="C02"&&o2==5&&start==4)
      tasks[[id]] <- data.table::data.table(task_id=id,pair_label=e$pair_label,
        model_context=e$model_context,prepared_index=e$prepared_index,seed_number=e$seed_number,
        O2_pct=o2,p_misseg=pm,initial_ploidy=start,threshold_N=66L,
        maximum_ploidy=as.integer(top),condition_mode="fitted_p_misseg",small_acceptance=small)
    }
  }
  manifest <- data.table::data.table(path=f7fit_sources(a));manifest[,sha256:=fp_sha(path)]
  input <- unique(data.table::rbindlist(c(list(pars$manifest),lapply(fits,`[[`,"manifest")),fill=TRUE),by="path")
  legacy <- fp_passage_rule(a) # exclusively the unchanged in-vivo count scale
  task <- data.table::rbindlist(tasks)
  b <- list(figure=7L,condition_mode="fitted_p_misseg",engine="Figure7_native_fitted_p_misseg",
    profile="full",policy=a$policy,config=a$cfg,shared=a$shared,prepared=pars$prepared,
    endpoints=pars$endpoints,protocols=protocols,rule=legacy$rule,replicates=3L,
    code_manifest=manifest,input_manifest=input,tasks=task)
  b$fingerprint <- digest::digest(list(policy=b$policy,config=b$config,
    code=manifest$sha256,inputs=input$sha256,tasks=task,protocols=protocols,
    parameters=lapply(b$prepared,`[[`,"native_parameters")),algo="sha256")
  fp_save(b,dest)
  fp_write(task,file.path(a$run,"prepared/task_manifest.tsv"))
  fp_write(pars$audit,file.path(a$run,"provenance/full_native_parameter_passthrough.tsv"))
  fp_write(pars$endpoints,file.path(a$run,"prepared/seed_input_manifest.tsv"))
  fp_write(input,file.path(a$run,"provenance/input_manifest.tsv"))
  fp_write(manifest,file.path(a$run,"provenance/current_code_manifest.tsv"))
  fp_write(task[,.(pair_label,model_context,seed_number,p_misseg)][!duplicated(paste(pair_label,model_context))],
    file.path(a$run,"provenance/fitted_p_misseg.tsv"))
  jsonlite::write_json(list(status="prepared",fingerprint=b$fingerprint,tasks=nrow(task),
    small_tasks=sum(task$small_acceptance),condition_mode=b$condition_mode),
    file.path(a$run,"prepared/preparation.json"),auto_unbox=TRUE,pretty=TRUE)
  cat("PREPARED_FITTED_PM",nrow(task),"workers; small",sum(task$small_acceptance),"\n")
  invisible(b)
}
f7fit_operator <- function(pr,family,seed,o2,top,beta,threshold=66L) {
  native <- pr$native_parameters;pm <- as.numeric(native["p_misseg"])
  stopifnot(is.finite(pm),pm>=.005,pm<=.5,beta%in%c(0,.5,1),top>=4,
    identical(as.numeric(pr$run_params$p_misseg),pm))
  pr$config$N_MAX <- as.integer(top*pr$config$N_UNIT)
  pr$run_params$o2_S0 <- o2
  # NA explicitly requests the original parameter rather than an override.
  l <- fp_redirect(fp_ledger(pr,o2,NA_real_),threshold,beta)
  stopifnot(identical(as.numeric(l$p_misseg_used),pm))
  c <- fp_birth_coefficients(l,threshold);loss <- f7p_loss(l,pr,o2,pm)
  sparse <- function(m)as(Matrix::Matrix(m,sparse=TRUE),"dgCMatrix")
  list(pr=pr,l=l,c=c,G=sparse(l$G),KH=sparse(c$KH),GU=sparse(l$G-c$KH),loss=loss,
    meta=data.table::data.table(family=family,seed_number=seed,O2_pct=o2,p_misseg=pm,
      maximum_ploidy=top,beta=beta,threshold_N=threshold,condition_mode="fitted_p_misseg"))
}
f7fit_runtime <- function(a,id) {
  root <- file.path(a$run,"runtime",sprintf("task_%05d",id))
  z <- f7p_runtime(a$workspace,root)
  for(name in c("continuous","daily"))Rcpp::sourceCpp(file.path(.f7fit_dir,
    paste0("figure7_",name,".cpp")),cacheDir=file.path(root,"instrumentation_cache"),showOutput=FALSE)
  z
}
f7fit_vivo <- function(op,b,initial,days) {
  l <- op$l;c <- op$c;x <- fp_pure_initial(l,initial);p <- f7_direct_operator(l,c,FALSE)
  streams <- fp_streams(op$meta$family,op$meta$seed_number,op$meta$O2_pct,
    op$meta$p_misseg,initial,1L)
  r <- f7_continuous_cpp(p$op$P,p$op$PU,matrix(x,ncol=1L),l$ngrid/l$N_unit,
    colSums(l$M),c$hl,c$ll,days,b$rule$seed_cells[1],b$rule$target_live_cells[1],
    streams$state,streams$label,FALSE,TRUE)
  d <- f7_daily_cpp(p$op,r,x,c$low,l$ngrid/l$N_unit,c$hl,c$ll,p$count_rows,
    FALSE,1,FALSE,0L,1e12,b$rule$seed_cells[1],days,TRUE,as.integer(0:days))
  stopifnot(d$passage_error<1e-10,d$final_state_error<1e-9,d$subset_error<1e-9,
    d$complement_error<1e-10,p$error<1e-9)
  X <- d$snapshot_X;D <- d$snapshot_D[c$low,,drop=FALSE];U <- d$snapshot_U[c$low,,drop=FALSE]
  factors <- colSums(p$op$P%*%X[,seq_len(days),drop=FALSE])
  logs <- log(b$rule$seed_cells[1])+c(0,cumsum(log(factors)))
  stopifnot(f7p_distance(as.numeric(d$log_start_live_cells)[-1],head(logs,-1))<1e-8,
    abs(tail(logs,1)-r$final_log_population[1])<1e-8)
  retained <- colSums(op$c$J);upper <- as.numeric(crossprod(op$loss$upper,X))
  lower <- as.numeric(crossprod(op$loss$lower,X));birth <- as.numeric(crossprod(retained,X))
  m <- data.table::data.table(day=0:days,log_cells=logs,
    log_start_live_cells=as.numeric(d$log_start_live_cells),
    mean_ploidy=as.numeric(r$mean_ploidy),low_fraction=colSums(X[c$low,,drop=FALSE]),
    A=as.numeric(d$high_parent_percent),B=as.numeric(d$low_parent_percent),
    HL_per_start_cell=as.numeric(d$HL_per_start_live_cell),LL_per_start_cell=as.numeric(d$LL_per_start_live_cell),
    crowding=1,upper_loss_fraction=ifelse(birth+upper+lower>0,upper/(birth+upper+lower),0))
  prefix <- c("all","direct","other","historical","never")
  for(g in seq_along(prefix))for(v in c("ploidy","division_rate","low_birth_rate")) {
    ix <- (g-1L)*5+switch(v,low_birth_rate=3L,division_rate=4L,ploidy=5L)
    m[,(paste(prefix[g],v,sep="_")):=as.numeric(d$metrics[[ix]])]
  }
  list(meta=cbind(op$meta,initial_ploidy=initial,route="continuous",block=0L),metrics=m,
    events=data.table::data.table(),X=X,D=D,U=U,ngrid=l$ngrid,low=c$low,
    coefficients=list(lambda=l$lambda[c$low],ll=c$ll[c$low]),
    numerical_reference_error=d$final_state_error,propagator_error=p$error)
}
f7fit_qc <- function(op,r,days,initial) {
  m <- r$metrics;ev <- r$events;low <- r$X[op$c$low,,drop=FALSE]
  expected <- ifelse(m$HL_per_start_cell+m$LL_per_start_cell>0,
    100*m$HL_per_start_cell/(m$HL_per_start_cell+m$LL_per_start_cell),NA_real_)
  errs <- c(complement=f7p_max(abs(m$A+m$B-100)),ratio=f7p_distance(expected,m$A),
    norm=max(abs(colSums(r$X)-1)),initial=max(abs(r$X[,1]-fp_pure_initial(op$l,initial))),
    subset=max(f7p_max(r$D-low),f7p_max(r$U-low)))
  culture <- op$pr$model_context=="in vitro"
  protocol_ok <- !culture||(sum(ev$retained_days)==days&&all(ev$selected_day>0)&
    all(ev$retained_days<=ev$selected_day&ev$selected_day<=ev$duration_days)&&!any(ev$uniform_fallback))
  empty <- if(op$meta$beta==1)all(r$D==0) else TRUE
  if(op$meta$beta==1&&initial==4)empty <- empty&&all(low==0)&&all(is.na(m$A))
  passed <- all(errs<1e-9)&&protocol_ok&&empty&&all(is.finite(m$log_cells))&&
    identical(as.numeric(m$day),as.numeric(0:days))&&identical(is.na(m$A),is.na(m$B))&&
    is.na(m$A[1])&&is.na(m$B[1])&&min(r$D)>=0&&min(r$U)>=0
  data.table::data.table(route=r$meta$route,days=days,daily_complement_error=errs[1],
    daily_ratio_error=errs[2],norm_error=errs[3],pure_initial_error=errs[4],subset_error=errs[5],
    maximum_upper_loss_fraction=f7p_max(m$upper_loss_fraction),
    passages=if(culture)sum(ev$passage_after_segment_applied) else 0L,
    native_fallbacks=if(culture)sum(ev$uniform_fallback) else 0L,passed=passed)
}
f7fit_run <- function(a,b,t,runtime,days=1000L,smoke=FALSE) {
  e <- b$endpoints[prepared_index==t$prepared_index];pr <- b$prepared[[t$prepared_index]]
  outputs <- list();source_rows <- list()
  for(beta in a$policy$beta) {
    message("FITTED_PM_BEGIN task=",t$task_id," context=",t$model_context,
      " grid=",t$maximum_ploidy," beta=",beta)
    op <- f7fit_operator(pr,t$pair_label,t$seed_number,t$O2_pct,t$maximum_ploidy,beta)
    stopifnot(identical(as.numeric(t$p_misseg),as.numeric(op$meta$p_misseg)))
    meta <- cbind(t,q=beta)
    if(t$model_context=="in vitro") {
      protocol <- b$protocols[[paste(t$pair_label,t$initial_ploidy,sep="_")]]
      prefix <- f7p_prefix(op,runtime$native,protocol,t$initial_ploidy)
      routes <- lapply(setNames(a$policy$routes,a$policy$routes),function(route)
        f7p_route(op,runtime$native,protocol,prefix,route,days,t$initial_ploidy))
    } else routes <- list(continuous=f7fit_vivo(op,b,t$initial_ploidy,days))
    qc <- data.table::rbindlist(lapply(routes,function(r)cbind(meta,f7fit_qc(op,r,days,t$initial_ploidy))))
    stopifnot(all(qc$passed))
    native <- NULL
    if(t$model_context=="in vitro"&&beta==0&&t$small_acceptance) {
      source(file.path(a$code,"test_Figure7_passaged.R"),local=globalenv())
      native <- f7p_native_checkpoint_checks(op,runtime$native,routes)
      stopifnot(nrow(native)>0,all(native$passed))
    }
    snapshots <- if(smoke)unique(as.integer(c(0,min(10,days),days))) else a$policy$snapshot_days
    d <- f7c_profiles(routes,op$l,op$c,snapshots)
    nominal <- length(routes)
    series <- data.table::rbindlist(lapply(c("flux_percent","growth_percent"),function(metric)
      cbind(meta,analysis="instantaneous",metric=metric,source_accounting="daily_viable_low_daughters",
        fp_summaries(if(metric=="flux_percent")d$high_parent_percent else d$low_parent_percent,nominal))))
    pop <- f7_population_summary(d,meta,nominal)
    dist <- f7_distribution_summary(d,meta,op$l,op$c)
    route_files <- character()
    for(route in names(routes)) {
      r <- routes[[route]]
      path <- file.path(a$run,"routes",sprintf("task_%05d_beta_%s_%s.rds",t$task_id,beta,route))
      r$checkpoints <- NULL;r$fingerprint <- b$fingerprint
      fp_save(r,path);route_files <- c(route_files,path)
    }
    outputs[[length(outputs)+1L]] <- list(metadata=meta,summary=series,population_summary=pop,
      distribution_summary=dist,source_counts_summary=f7_count_summary(d,meta,nominal),qc=qc,
      native_checkpoint_checks=native,route_files=route_files,route_sha256=fp_sha(route_files))
    message("FITTED_PM_BETA_COMPLETE task=",t$task_id," beta=",beta)
    rm(op,routes,d);gc(verbose=FALSE)
  }
  list(status="task_validated",fingerprint=b$fingerprint,task=t,outputs=outputs)
}
f7fit_task_path <- function(a,id)file.path(a$run,"tasks",sprintf("task_%05d.rds",id))
f7fit_simulate <- function(a) {
  b <- f7fit_bundle(a)
  ids <- if(!is.null(a$args$task_index))as.integer(a$args$task_index) else
    if(!is.null(a$args$task_ids))as.integer(strsplit(a$args$task_ids,",",fixed=TRUE)[[1]]) else
      b$tasks$task_id[if(a$args$phase=="small")b$tasks$small_acceptance else TRUE]
  stopifnot(length(ids)>0,!anyNA(ids),all(ids%in%b$tasks$task_id),!anyDuplicated(ids))
  if(a$args$phase!="small") {
    gate <- jsonlite::fromJSON(file.path(a$run,"provenance/small_acceptance.json"))
    stopifnot(gate$status=="passed",gate$fingerprint==b$fingerprint)
  } else stopifnot(all(b$tasks$small_acceptance[match(ids,b$tasks$task_id)]))
  for(id in ids) {
    dest <- f7fit_task_path(a,id)
    if(file.exists(dest)) {
      old <- readRDS(dest);stopifnot(old$status=="task_validated",old$fingerprint==b$fingerprint,
        all(vapply(old$outputs,function(x)all(x$qc$passed),TRUE)))
      receipt <- file.path(a$run,"receipts",sprintf("task_%05d.tsv",id))
      if(!file.exists(receipt))fp_write(data.table::data.table(task_id=id,status=old$status,
        fingerprint=b$fingerprint,sha256=fp_sha(dest),elapsed_seconds=old$elapsed_seconds),receipt)
      next
    }
    t <- b$tasks[task_id==id];tic <- proc.time()[3]
    runtime <- f7fit_runtime(a,id)
    value <- f7fit_run(a,b,t,runtime);value$elapsed_seconds <- proc.time()[3]-tic
    fp_save(value,dest)
    fp_write(data.table::data.table(task_id=id,status=value$status,fingerprint=b$fingerprint,
      sha256=fp_sha(dest),elapsed_seconds=value$elapsed_seconds),
      file.path(a$run,"receipts",sprintf("task_%05d.tsv",id)))
    cat("FITTED_PM_TASK_COMPLETE",id,"seconds",round(value$elapsed_seconds,2),"\n");flush.console()
  }
}
f7fit_read_tasks <- function(a,b,small=FALSE) {
  ids <- b$tasks$task_id[if(small)b$tasks$small_acceptance else TRUE]
  lapply(ids,function(id) {
    path <- f7fit_task_path(a,id);receipt <- data.table::fread(file.path(a$run,"receipts",sprintf("task_%05d.tsv",id)))
    stopifnot(receipt$sha256==fp_sha(path),receipt$fingerprint==b$fingerprint)
    z <- readRDS(path);stopifnot(z$status=="task_validated",z$fingerprint==b$fingerprint);z
  })
}
f7fit_grid <- function(a,b,tasks) {
  outputs <- unlist(lapply(tasks,`[[`,"outputs"),recursive=FALSE)
  specs <- data.table::rbindlist(lapply(outputs,function(x)cbind(x$metadata,
    route_files=list(x$route_files),route_sha256=list(x$route_sha256))))
  rows <- list()
  keys <- c("pair_label","model_context","O2_pct","initial_ploidy","q")
  groups <- unique(specs[,..keys])
  for(i in seq_len(nrow(groups))) {
    g <- groups[i];z <- merge(specs,g,by=keys);data.table::setorder(z,maximum_ploidy)
    stopifnot(nrow(z)==3L)
    for(j in 2:3) {
      old_paths <- z$route_files[[j-1]];new_paths <- z$route_files[[j]]
      stopifnot(identical(fp_sha(old_paths),z$route_sha256[[j-1]]),
        identical(fp_sha(new_paths),z$route_sha256[[j]]))
      for(k in seq_along(old_paths)) {
        old <- readRDS(old_paths[k]);now <- readRDS(new_paths[k])
        stopifnot(old$fingerprint==b$fingerprint,now$fingerprint==b$fingerprint)
        row <- f7p_compare(old,now,a$policy$grid_limits)
        row[,model_context:=g$model_context];rows[[length(rows)+1L]] <- row
        rm(old,now);gc(verbose=FALSE)
      }
    }
  }
  data.table::rbindlist(rows,fill=TRUE)
}
f7fit_acceptance <- function(a) {
  b <- f7fit_bundle(a);tasks <- f7fit_read_tasks(a,b,TRUE)
  grid <- f7fit_grid(a,b,tasks)
  qc <- data.table::rbindlist(lapply(tasks,function(t)data.table::rbindlist(lapply(t$outputs,`[[`,"qc"))))
  native <- data.table::rbindlist(unlist(lapply(tasks,function(t)lapply(t$outputs,`[[`,"native_checkpoint_checks")),recursive=FALSE),fill=TRUE)
  fp_write(grid,file.path(a$run,"provenance/small_grid_comparisons.tsv"))
  fp_write(qc,file.path(a$run,"provenance/small_numerical_checks.tsv"))
  fp_write(native,file.path(a$run,"provenance/small_native_checkpoint_checks.tsv"))
  passed <- nrow(grid)>0&&all(grid$passed)&&all(qc$passed)&&nrow(native)>0&&all(native$passed)
  jsonlite::write_json(list(status=if(passed)"passed" else "failed",fingerprint=b$fingerprint,
    representative_workers=length(tasks),grid_comparisons=nrow(grid),numerical_checks=nrow(qc),
    native_checkpoints=nrow(native),full_calculation_released=passed),
    file.path(a$run,"provenance/small_acceptance.json"),pretty=TRUE,auto_unbox=TRUE)
  if(!passed)stop("Small acceptance/grid gate failed; full calculation remains blocked")
  cat("FITTED_PM_SMALL_ACCEPTANCE_PASSED\n")
}
f7fit_telemetry <- function(r,meta) {
  low <- r$X[r$low,,drop=FALSE];groups <- list(all_low=low,direct_first_generation=r$D,
    other_low=pmax(low-r$D,0),historical_lineage=pmax(low-r$U,0),never_HL_lineage=r$U)
  prefix <- c("all","direct","other","historical","never")
  data.table::rbindlist(lapply(seq_along(groups),function(i) {
    g <- names(groups)[i];m <- r$metrics
    z <- cbind(meta,group=g,route=r$meta$route,day=m$day,mass_fraction=colSums(groups[[i]]),
      log_total_cells=m$log_cells,ancestry_identifiable=if("ancestry_identifiable"%in%names(m))m$ancestry_identifiable else TRUE,
      mean_low_birth_rate=m[[paste0(prefix[i],"_low_birth_rate")]],
      mean_division_rate=m[[paste0(prefix[i],"_division_rate")]],mean_ploidy=m[[paste0(prefix[i],"_ploidy")]])
    cbind(z,f7pv_flags(z,1e-8))
  }))
}
f7fit_aggregate <- function(a) {
  b <- f7fit_bundle(a);gate <- jsonlite::fromJSON(file.path(a$run,"provenance/small_acceptance.json"))
  stopifnot(gate$status=="passed",gate$fingerprint==b$fingerprint)
  tasks <- f7fit_read_tasks(a,b);grid <- f7fit_grid(a,b,tasks)
  fp_write(grid,file.path(a$run,"provenance/grid_comparisons.tsv"))
  if(!all(grid$passed))stop("Full grid validation failed; no publication")
  outputs <- unlist(lapply(tasks,`[[`,"outputs"),recursive=FALSE)
  primary <- Filter(function(x)x$metadata$maximum_ploidy==a$policy$primary_grid_N[[x$metadata$model_context]],outputs)
  tables <- list(instant_time_summary="summary",low_population_time_summary="population_summary",
    low_population_distribution_summary="distribution_summary",daily_source_counts_summary="source_counts_summary")
  for(name in names(tables))fp_write(data.table::rbindlist(lapply(primary,`[[`,tables[[name]]),fill=TRUE),
    file.path(a$run,"tables",paste0(name,".tsv")))
  qc <- data.table::rbindlist(lapply(outputs,`[[`,"qc"));fp_write(qc,file.path(a$run,"provenance/task_numerical_validation.tsv"))
  telemetry <- events <- list()
  for(x in primary)for(path in x$route_files) {
    r <- readRDS(path);telemetry[[length(telemetry)+1L]] <- f7fit_telemetry(r,x$metadata)
    if(nrow(r$events))events[[length(events)+1L]] <- cbind(x$metadata,r$events)
  }
  z <- data.table::rbindlist(telemetry,fill=TRUE)
  fp_save(z,file.path(a$run,"postprocess/route_metrics.rds"))
  fp_write(z[,.(pair_label,model_context,O2_pct,p_misseg,initial_ploidy,q,day,group,route,
    mass_fraction,ancestry_identifiable,mean_low_birth_rate,mean_division_rate,mean_ploidy,
    log_total_cells,valid,reason,structurally_empty,log_expected_cells,expected_cells_lt1,expected_cells_status)],
    file.path(a$run,"postprocess/route_validity.tsv.gz"))
  key <- c("pair_label","model_context","O2_pct","p_misseg","initial_ploidy","q","day","group")
  mask <- z[,.(all_valid=all(valid),n_valid=sum(valid),n_routes=.N),by=key]
  pop <- data.table::fread(file.path(a$run,"tables/low_population_time_summary.tsv"))
  pop <- merge(pop,mask,by=key,all.x=TRUE)
  means <- pop$metric%in%c("mean_low_birth_rate","mean_division_rate","mean_ploidy")
  stopifnot(!anyNA(pop$all_valid))
  pop[means&!all_valid,c("median","q25","q75"):=list(NA_real_,NA_real_,NA_real_)]
  pop[means,n_valid.x:=n_valid.y]
  pop[,c("all_valid","n_valid.y","n_routes"):=NULL]
  data.table::setnames(pop,"n_valid.x","n_valid")
  fp_write(pop,file.path(a$run,"postprocess/low_population_time_summary.tsv"))
  fp_write(data.table::rbindlist(events,fill=TRUE),file.path(a$run,"provenance/passage_events.tsv"))
  fp_write(z[,.(points=.N,valid=sum(valid),valid_below_one=sum(valid&expected_cells_lt1,na.rm=TRUE)),
    by=.(model_context,group)],file.path(a$run,"postprocess/coverage.tsv"))
  sensitivity <- lapply(c(1e-9,1e-8,1e-7),function(rel) {
    f <- f7pv_flags(z,rel);zz <- data.table::copy(z[,..key]);zz[,valid:=f$valid]
    ss <- zz[,.(valid=all(valid)),by=key]
    ss[,.(summary_points=.N,valid_points=sum(valid)),by=.(model_context,group)][,
      `:=`(relative_mass_min=rel,absolute_count_filter=FALSE)]
  })
  fp_write(data.table::rbindlist(sensitivity),file.path(a$run,"postprocess/sensitivity_summary.tsv"))
  jsonlite::write_json(list(status="aggregated",fingerprint=b$fingerprint,primary_conditions=length(primary),
    relative_mass_min=1e-8,absolute_count_filter=FALSE),file.path(a$run,"provenance/aggregation.json"),auto_unbox=TRUE,pretty=TRUE)
  cat("FITTED_PM_AGGREGATED",length(primary),"conditions\n")
}
f7fit_validate <- function(a) {
  b <- f7fit_bundle(a);checks <- list()
  add <- function(n,p)checks[[length(checks)+1L]] <<- data.table::data.table(test=n,passed=isTRUE(p))
  grid <- data.table::fread(file.path(a$run,"provenance/grid_comparisons.tsv"))
  qc <- data.table::fread(file.path(a$run,"provenance/task_numerical_validation.tsv"))
  add("every_grid_and_numerical_check_passed",nrow(grid)>0&&all(grid$passed)&&all(qc$passed))
  for(i in seq_len(nrow(b$endpoints))) {
    e <- b$endpoints[i];pr <- b$prepared[[e$prepared_index]]
    add(paste(e$pair_label,e$model_context,"full_parameter_count"),length(pr$native_parameters)==if(e$model_context=="in vitro")33L else 26L)
    add(paste(e$pair_label,e$model_context,"native_p_misseg_preserved"),
      identical(as.numeric(pr$run_params$p_misseg),as.numeric(pr$native_parameters["p_misseg"])))
  }
  z <- data.table::fread(file.path(a$run,"tables/instant_time_summary.tsv"))
  keys <- setdiff(names(z),c("metric","median","q25","q75","n_valid"))
  A <- z[metric=="flux_percent"];B <- z[metric=="growth_percent"]
  data.table::setorderv(A,c("task_id","q","day"));data.table::setorderv(B,c("task_id","q","day"))
  add("A_B_same_defined_windows",identical(is.na(A$median),is.na(B$median)))
  add("A_plus_B_100",f7p_max(abs(A$median+B$median-100))<1e-10)
  add("IQR_complements",f7p_max(abs(A$q25+B$q75-100))<1e-10&&f7p_max(abs(A$q75+B$q25-100))<1e-10)
  add("day_zero_undefined",all(is.na(z[day==0]$median)))
  add("new_mode_only",all(z$condition_mode=="fitted_p_misseg"))
  add("all_48_primary_blocked_conditions",data.table::uniqueN(z[,.(pair_label,model_context,O2_pct,initial_ploidy,q)])==48L)
  rt <- readRDS(file.path(a$run,"postprocess/route_metrics.rds"))
  fresh <- f7pv_flags(rt[,setdiff(names(rt),names(f7pv_flags(rt))),with=FALSE])
  add("reporting_mask_recomputed",identical(rt$valid,fresh$valid))
  add("below_one_is_not_a_filter",any(rt$valid&rt$expected_cells_lt1,na.rm=TRUE)||!any(rt$expected_cells_lt1,na.rm=TRUE))
  add("theoretical_empty_groups_invalid",!any(rt$valid&rt$structurally_empty))
  ev <- data.table::fread(file.path(a$run,"provenance/passage_events.tsv"))
  add("native_no_fallback",nrow(ev)>0&&!any(ev$uniform_fallback))
  add("all_culture_routes_1000_days",all(ev[,.(total=sum(retained_days)),by=.(task_id,q,route)]$total==1000))
  dist <- data.table::fread(file.path(a$run,"tables/low_population_distribution_summary.tsv"))
  add("CDF_probability_range",all(dist[is.finite(cdf)]$cdf>=-1e-9&dist[is.finite(cdf)]$cdf<=1+1e-9))
  add("empty_CDF_NA",all(is.na(dist[n_valid_repeats==0]$cdf)))
  add("no_historical_or_Figure8_products",!any(grepl("history_",list.files(file.path(a$run,"tables")))))
  checks <- data.table::rbindlist(checks);fp_write(checks,file.path(a$run,"provenance/numerical_validation.tsv"))
  if(!all(checks$passed))stop("Fitted-p validation failed: ",paste(checks[passed==FALSE]$test,collapse=", "))
  files <- list.files(a$run,recursive=TRUE,full.names=TRUE)
  files <- files[file.info(files)$isdir==FALSE&grepl("/(tables|postprocess|provenance)/",files)&
    !grepl("completion.json$|validated_files.tsv$|render|publication",files)]
  manifest <- data.table::data.table(relative_path=substring(files,nchar(a$run)+2L),sha256=fp_sha(files))
  fp_write(manifest,file.path(a$run,"provenance/validated_files.tsv"))
  jsonlite::write_json(list(status="validated_complete",condition_mode="fitted_p_misseg",
    fingerprint=b$fingerprint,bundle_sha256=fp_sha(file.path(a$run,"prepared/analysis_bundle.rds")),
    manifest_sha256=fp_sha(file.path(a$run,"provenance/validated_files.tsv")),checks=nrow(checks),
    existing_product_writes=FALSE),file.path(a$run,"provenance/completion.json"),auto_unbox=TRUE,pretty=TRUE)
  cat("FITTED_PM_VALIDATED",nrow(checks),"checks\n")
}
f7fit_smoke <- function(a) {
  # An explicitly short, isolated test is never eligible for acceptance/publication.
  stopifnot(grepl("/audit/",a$run))
  b <- f7fit_prepare(a);runtime <- f7fit_runtime(a,0L)
  t <- b$tasks[pair_label=="C01"&model_context=="in vitro"&O2_pct==.5&initial_ploidy==2][1]
  t$maximum_ploidy <- 7L
  z <- f7fit_run(a,b,t,runtime,30L,TRUE)
  t2 <- b$tasks[pair_label=="C02"&model_context=="in vivo"&O2_pct==5&initial_ploidy==4][1]
  z2 <- f7fit_run(a,b,t2,runtime,30L,TRUE)
  qc <- data.table::rbindlist(lapply(c(z$outputs,z2$outputs),`[[`,"qc"))
  # Compare the new continuous adapter with the complete established Figure7
  # implementation at exactly the same fitted vector, start, oxygen and horizon.
  for(name in c("history","direct"))Rcpp::sourceCpp(file.path(.f7fit_dir,
    paste0("figure7_",name,".cpp")),cacheDir=file.path(a$run,"runtime/task_00000/instrumentation_cache"),showOutput=FALSE)
  rb <- b;rb$config$day <- 30L;rb$config$snapshot_days <- c(0L,10L,30L)
  ref <- f7_simulate_task(a,rb,t2)
  for(j in seq_along(ref$outputs)) {
    old <- ref$outputs[[j]];new <- z2$outputs[[j]]
    aa <- old$summary[analysis=="instantaneous"];bb <- new$summary
    data.table::setorderv(aa,c("metric","day"));data.table::setorderv(bb,c("metric","day"))
    stopifnot(f7p_distance(aa$median,bb$median)<1e-9)
    rr <- readRDS(new$route_files[1])
    stopifnot(max(abs(old$trajectories$final_state[,1]-rr$X[,31]))<1e-9)
    pa <- old$population_summary;pb <- new$population_summary
    data.table::setorderv(pa,c("group","metric","day"));data.table::setorderv(pb,c("group","metric","day"))
    stopifnot(f7p_distance(pa$median,pb$median)<1e-8)
  }
  fp_write(qc,file.path(a$run,"smoke_checks.tsv"))
  stopifnot(all(qc$passed))
  jsonlite::write_json(list(status="passed",days=30,publication_eligible=FALSE,
    native_parameters=TRUE,checks=nrow(qc)),file.path(a$run,"smoke.json"),auto_unbox=TRUE,pretty=TRUE)
  cat("FITTED_PM_SMOKE_PASSED\n")
}
f7fit_cli <- function(script) {
  args <- f7fit_args();a <- f7fit_init(script,args)
  switch(args$stage,plan={cat("Figure7 fitted-p mode; distinct native parameter for every family/environment.\n")},
    prepare=f7fit_prepare(a),smoke=f7fit_smoke(a),simulate=f7fit_simulate(a),
    acceptance=f7fit_acceptance(a),aggregate=f7fit_aggregate(a),validate=f7fit_validate(a),
    plot={source(file.path(.f7fit_dir,"figure7_fitted_plotting.R"),local=globalenv());f7fit_plot(a)},
    stop("Unknown fitted-p stage"))
  invisible(a)
}
