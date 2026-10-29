# Figure7 current products only. Figure8 and historical products are read-only.
# This adapter uses the independently accepted native-fit passage implementation;
# it never calls fp_culture for an in-vitro task or consumes an old prediction.
.f7c_dir <- dirname(normalizePath(sys.frame(1)$ofile,mustWork=TRUE))
source(file.path(.f7c_dir,"figure7_passaged.R"),local=globalenv())

f7c_policy <- function(a) {
  p <- jsonlite::fromJSON(file.path(a$workspace,"Code/config/figure7_passaged.json"))
  stopifnot(p$days==1000L,p$template_segments==5L,!p$inoculum_feasibility_stop,
    !p$random_sampling,identical(as.numeric(p$p_misseg),c(.1,.3)),
    identical(as.numeric(p$initial_ploidy),c(2,4)),
    identical(as.numeric(p$candidate_grid_upper_N),c(137,139,141)))
  p
}
f7c_tasks <- function(a) {
  z <- f7_tasks(a);z[,base_task_id:=task_id]
  p <- f7c_policy(a)
  out <- data.table::rbindlist(lapply(seq_len(nrow(z)),function(i) {
    t <- z[i]
    t <- t[rep(1L,if(t$model_context=="in vitro")3L else 1L)]
    t[,maximum_ploidy:=if(model_context[1]=="in vitro")as.integer(p$candidate_grid_upper_N) else
      as.integer(a$cfg$model_grid$N_MAX/a$cfg$model_grid$N_UNIT)]
    t
  }))
  out[,task_id:=seq_len(.N)]
  out[,engine:=ifelse(model_context=="in vitro","native_fit_last5","unchanged_continuous")]
  out
}
f7c_sources <- function(a) {
  unique(c(fp_files(a),file.path(a$code,"util/analysis",c("figure7_current_pipeline.R",
    "figure7_passaged.R","figure7_passaged.cpp")),a$config_path,
    file.path(a$workspace,"Code/config",c("figure7_passaged.json","figure78_inputs.json",
      "manifests/figure78_model_sources.tsv"))))
}
f7c_not_formal <- function(a) {
  if(identical(normalizePath(a$run,mustWork=FALSE),
    normalizePath(file.path(a$workspace,"data/Figures/Figure7"),mustWork=TRUE)))
    stop("Generate into a separate campaign root; promotion is a separately gated file allowlist operation")
}
f7c_prepare <- function(a) {
  f7c_not_formal(a);dest <- file.path(a$run,"prepared/analysis_bundle.rds")
  if(file.exists(dest))return(invisible(f7c_bundle(a)))
  if(dir.exists(a$run)&&length(list.files(a$run,all.files=TRUE,no..=TRUE)))
    stop("Refusing to replace a nonempty unprepared run")
  fp_preflight(a);p <- f7c_policy(a)
  dir.create(a$run,recursive=TRUE,showWarnings=FALSE)
  # Retain the normal prepared-model layout for existing container launchers.
  locked <- f78_locked_model_manifest(a$workspace,a$model_root)
  for(i in seq_len(nrow(locked))) {
    target <- file.path(a$run,"prepared/model_source",locked$relative_path[i])
    dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
    if(!file.copy(locked$path[i],target,overwrite=FALSE))stop("Could not stage locked model")
    stopifnot(identical(fp_sha(target),locked$sha256[i]))
  }
  runtime <- f7p_runtime(a$workspace,file.path(a$run,"runtime/preparation"))
  pars <- fp_load_parameters(a,fp_endpoint_spec(a))
  # This legacy scale is used ONLY by the unchanged continuous in-vivo adapter.
  # No in-vitro initial count, target, selection or reseeding is taken from it.
  legacy_vivo <- fp_passage_rule(a)
  fits <- lapply(setNames(c("C01","C02"),c("C01","C02")),function(f)
    f7p_load_fit(a$workspace,f,a$fit_root))
  protocols <- list();inventory <- templates <- list()
  for(f in names(fits))for(initial in c(2,4)) {
    pp <- f7p_protocol(fits[[f]]$fit$best_components$invitro[[paste0("run_",initial,"N")]]$adapter,5L)
    protocols[[paste(f,initial,sep="_")]] <- pp
    inventory[[length(inventory)+1L]] <- cbind(pair_label=f,initial_ploidy=initial,pp$inventory)
    for(route in names(pp$templates))templates[[length(templates)+1L]] <- cbind(
      pair_label=f,initial_ploidy=initial,route=route,template_order=seq_len(5),
      pp$inventory[match(pp$templates[[route]],pp$inventory$segment_id)])
  }
  inputs <- unique(data.table::rbindlist(c(list(pars$manifest,fp_model_manifest(a)),
    lapply(fits,function(x)data.table::copy(x$manifest)[,role:="complete native winning fit"])),fill=TRUE))
  inputs <- inputs[!duplicated(path)]
  extra <- data.table::data.table(path=legacy_vivo$path,sha256=fp_sha(legacy_vivo$path),role="legacy in-vivo count scale only")
  inputs <- unique(data.table::rbindlist(list(inputs,extra)),by="path")
  codes <- data.table::data.table(path=f7c_sources(a));codes[,sha256:=fp_sha(path)]
  tasks <- f7c_tasks(a)
  fingerprint <- digest::digest(list(engine="Figure7_native_fit_last5_current",config=a$cfg,
    policy=p,shared=a$shared,profile=a$profile,codes=codes$sha256,inputs=inputs$sha256,
    endpoints=pars$endpoints,protocols=lapply(protocols,function(x)x[c("inventory","routes","templates")]),
    tasks=tasks),algo="sha256")
  b <- list(figure=7L,engine="Figure7_native_fit_last5_current",config=a$cfg,policy=p,
    shared=a$shared,profile=a$profile,replicates=3L,prepared=pars$prepared,endpoints=pars$endpoints,
    protocols=protocols,tasks=tasks,rule=legacy_vivo$rule,fingerprint=fingerprint,
    code_manifest=codes,input_manifest=inputs,model_manifest=fp_model_manifest(a))
  fp_write(pars$endpoints,file.path(a$run,"prepared/seed_input_manifest.tsv"))
  fp_write(pars$audit,file.path(a$run,"provenance/full_native_parameter_passthrough.tsv"))
  fp_write(tasks,file.path(a$run,"prepared/task_manifest.tsv"))
  fp_write(inputs,file.path(a$run,"provenance/input_manifest.tsv"))
  fp_write(codes,file.path(a$run,"provenance/current_code_manifest.tsv"))
  fp_write(data.table::rbindlist(inventory),file.path(a$run,"provenance/native_passage_inventory.tsv"))
  fp_write(data.table::rbindlist(templates),file.path(a$run,"provenance/native_last5_templates.tsv"))
  fp_write(legacy_vivo$rule,file.path(a$run,"provenance/invivo_legacy_count_scale.tsv"))
  fp_write(data.frame(session=capture.output(sessionInfo())),file.path(a$run,"provenance/runtime_session.tsv"))
  fp_save(b,dest)
  fp_write(data.frame(figure=7L,engine=b$engine,profile=a$profile,fingerprint=fingerprint,status="prepared",
    n_tasks=nrow(tasks),n_base_conditions=data.table::uniqueN(tasks$base_task_id),
    n_context_seed_inputs=nrow(pars$endpoints),historical_product_writes=FALSE,
    figure8_changes=FALSE),file.path(a$run,"provenance/preparation.tsv"))
  message("Prepared ",nrow(tasks)," workers for ",data.table::uniqueN(tasks$base_task_id)," conditions")
  invisible(b)
}
f7c_bundle <- function(a) {
  b <- fp_bundle(a)
  if(!identical(b$engine,"Figure7_native_fit_last5_current")||!identical(b$policy,f7c_policy(a)))
    stop("Not a current Figure7 native-fit passage campaign; old caches cannot be reused")
  b
}
f7c_route_path <- function(a,id,beta,route)file.path(a$run,"routes",
  sprintf("task_%05d_beta_%s_%s.rds",id,format(beta,trim=TRUE),route))
f7c_save_route <- function(r,path,fingerprint) {
  # Keep the full daily state for independent convergence; do not keep repeated
  # full candidate windows once checkpoint verification has been registered.
  r$checkpoints <- NULL;r$fingerprint <- fingerprint
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(".route_",tmpdir=dirname(path));saveRDS(r,tmp,compress="gzip")
  if(!file.rename(tmp,path))stop("Atomic route save failed")
  fp_write(data.frame(sha256=fp_sha(path),fingerprint=fingerprint),paste0(path,".sha.tsv"))
  invisible(r)
}
f7c_profiles <- function(routes,l,c,snapshots) {
  k <- length(routes);nd <- nrow(routes[[1]]$metrics)
  d <- list(high_parent_percent=do.call(rbind,lapply(routes,function(r)r$metrics$A)),
    low_parent_percent=do.call(rbind,lapply(routes,function(r)r$metrics$B)),
    HL_per_start_live_cell=do.call(rbind,lapply(routes,function(r)r$metrics$HL_per_start_cell)),
    LL_per_start_live_cell=do.call(rbind,lapply(routes,function(r)r$metrics$LL_per_start_cell)),
    log_start_live_cells=do.call(rbind,lapply(routes,function(r)r$metrics$log_start_live_cells)),
    group_names=c("all_low","direct_first_generation","other_low","historical_lineage","never_HL_lineage"),
    metric_names=c("cell_share_percent","low_birth_share_percent","mean_low_birth_rate","mean_division_rate","mean_ploidy"),
    metrics=list(),snapshot_days=as.integer(snapshots),snapshot_crowding=numeric())
  prefix <- c("all","direct","other","historical","never")
  for(g in seq_along(prefix))for(metric in d$metric_names) {
    values <- lapply(routes,function(r) {
      X <- r$X[c$low,,drop=FALSE];groups <- list(X,r$D,pmax(X-r$D,0),pmax(X-r$U,0),r$U)
      mass <- colSums(groups[[g]]);total <- colSums(X)
      if(metric=="cell_share_percent")return(ifelse(total>0,100*mass/total,NA_real_))
      if(metric=="low_birth_share_percent") {
        numerator <- colSums(sweep(groups[[g]],1,c$ll[c$low],"*"))*r$metrics$crowding
        denominator <- colSums(sweep(X,1,c$ll[c$low],"*"))*r$metrics$crowding
        return(ifelse(denominator>0,100*numerator/denominator,NA_real_))
      }
      col <- paste0(prefix[g],switch(metric,mean_low_birth_rate="_low_birth_rate",
        mean_division_rate="_division_rate",mean_ploidy="_ploidy"))
      r$metrics[[col]]
    })
    d$metrics[[length(d$metrics)+1L]] <- do.call(rbind,values)
  }
  sx <- sd <- su <- list()
  for(day in snapshots)for(r in routes) {
    j <- day+1L;dx <- ux <- numeric(length(l$ngrid));dx[c$low] <- r$D[,j];ux[c$low] <- r$U[,j]
    sx[[length(sx)+1L]] <- r$X[,j];sd[[length(sd)+1L]] <- dx;su[[length(su)+1L]] <- ux
    d$snapshot_crowding <- c(d$snapshot_crowding,r$metrics$crowding[j])
  }
  d$snapshot_X <- do.call(cbind,sx);d$snapshot_D <- do.call(cbind,sd);d$snapshot_U <- do.call(cbind,su)
  d
}
f7c_diagnostics <- function(op,routes,meta) {
  l <- op$l;low <- op$c$low
  X <- do.call(cbind,lapply(routes,function(r)r$X[,ncol(r$X)]))
  cf <- vapply(routes,function(r)tail(r$metrics$crowding,1),0.0)
  weighted <- rowMeans(sweep(X,2,cf,"*"))
  bins <- floor(l$ngrid/l$N_unit/.25)*.25
  heat <- lapply(c("Ordinary division","WGD"),function(branch) {
    J <- if(branch=="Ordinary division")sweep(l$B,2,l$lambda*(1-l$p_wgd),"*") else
      sweep(l$W,2,l$lambda*l$p_wgd,"*")
    trip <- Matrix::summary(Matrix::Matrix(J,sparse=TRUE))
    h <- data.table::data.table(source=bins[trip$j],daughter=bins[trip$i],
      flux=1000*trip$x*weighted[trip$j])[,.(flux=sum(flux)),by=.(source,daughter)]
    grid <- data.table::CJ(source=sort(unique(bins)),daughter=sort(unique(bins)))
    h <- merge(grid,h,by=c("source","daughter"),all.x=TRUE,sort=FALSE);h[is.na(flux),flux:=0]
    cbind(meta,branch=branch,h)
  })
  H <- sum(op$c$hl*weighted);L <- sum(op$c$ll*weighted)
  removal <- sum(l$lambda[low]*weighted[low])
  death <- mean(colSums(sweep(X,1,l$mu,"*")[low,,drop=FALSE]))
  budget <- cbind(meta,external_HL=H,internal_LL=L,maternal_removal=removal,direct_death=death,
    low_net_change=H+L-removal-death,low_fraction=mean(colSums(X[low,,drop=FALSE])))
  list(heat=data.table::rbindlist(heat),budget=budget)
}
f7c_route_qc <- function(op,routes,meta) {
  checks <- lapply(routes,function(r) {
    m <- r$metrics;low <- r$X[op$c$low,,drop=FALSE];pure <- fp_pure_initial(op$l,meta$initial_ploidy)
    A <- m$A;B <- m$B;valid <- is.finite(A)
    complement <- f7p_max(abs(A[valid]+B[valid]-100))
    expected <- ifelse(m$HL_per_start_cell+m$LL_per_start_cell>0,
      100*m$HL_per_start_cell/(m$HL_per_start_cell+m$LL_per_start_cell),NA_real_)
    subset <- max(f7p_max(r$D-low),f7p_max(r$U-low))
    ev <- r$events
    used <- ev$retained_days;selected <- ev$selected_day
    ok <- identical(as.numeric(m$day),as.numeric(0:1000))&&all(is.finite(m$log_cells))&&
      max(abs(colSums(r$X)-1))<1e-10&&max(abs(r$X[,1]-pure))<1e-12&&
      complement<1e-10&&f7p_distance(expected,A)<1e-10&&
      identical(is.na(A),is.na(B))&&is.na(A[1])&&is.na(B[1])&&
      all(used>0&used<=selected&selected<=ev$duration_days)&&sum(used)==1000&&
      !any(ev$uniform_fallback)&&subset<1e-10&&min(r$D)>=0&&min(r$U)>=0
    if(meta$q==1&&meta$initial_ploidy*op$l$N_unit>=meta$threshold_N)
      ok <- ok&&all(low==0)&&all(is.na(A))&&all(is.na(B))
    cbind(meta,route=r$meta$route,complete_days=nrow(m)-1L,
      daily_complement_error=complement,daily_ratio_error=f7p_distance(expected,A),
      maximum_source_subset_error=subset,max_upper_loss_fraction=f7p_max(m$upper_loss_fraction),
      n_passages=sum(ev$passage_after_segment_applied),
      stock_shortages_recorded_not_stopped=sum(ev$insufficient_stock_recorded_only),
      roundoff_corrections=sum(ev$roundoff_corrections),maximum_subset_roundoff=f7p_max(ev$subset_error),
      native_uniform_fallback=sum(ev$uniform_fallback),passed=ok)
  })
  data.table::rbindlist(checks)
}
f7c_vitro_task <- function(a,b,t,runtime) {
  input <- f7p_load_fit(a$workspace,t$pair_label,a$fit_root)
  protocol <- b$protocols[[paste(t$pair_label,t$initial_ploidy,sep="_")]]
  outputs <- list()
  for(beta in b$config$intervention_q) {
    meta <- cbind(t,seed_number=input$seed,q=beta)
    op <- f7p_operator(input,t$O2_pct,t$p_misseg,t$maximum_ploidy,beta,t$threshold_N)
    files <- setNames(vapply(b$policy$routes,function(route)f7c_route_path(a,t$task_id,beta,route),""),b$policy$routes)
    prefix <- NULL
    routes <- lapply(names(files),function(route) {
      path <- files[[route]]
      if(file.exists(path)) {
        sidecar <- paste0(path,".sha.tsv")
        if(!file.exists(sidecar))stop("Unregistered interrupted route retained at ",path,
          "; move it to a recoverable audit location before retrying")
        registered <- data.table::fread(sidecar)
        stopifnot(nrow(registered)==1L,registered$sha256==fp_sha(path),registered$fingerprint==b$fingerprint)
        r <- readRDS(path)
        stopifnot(identical(r$fingerprint,b$fingerprint),r$meta$route==route,
          r$meta$family==t$pair_label,r$meta$maximum_ploidy==t$maximum_ploidy,
          r$meta$beta==beta,r$meta$threshold_N==t$threshold_N,r$meta$initial_ploidy==t$initial_ploidy,
          r$meta$p_misseg==t$p_misseg,r$meta$O2_pct==t$O2_pct)
        return(r)
      }
      if(is.null(prefix))prefix <<- f7p_prefix(op,runtime$native,protocol,t$initial_ploidy)
      r <- f7p_route(op,runtime$native,protocol,prefix,route,1000L,t$initial_ploidy)
      f7c_save_route(r,path,b$fingerprint)
    });names(routes) <- names(files)
    qc <- f7c_route_qc(op,routes,meta)
    if(!all(qc$passed))stop("Figure7 native passage route validation failed")
    daily <- f7c_profiles(routes,op$l,op$c,b$config$snapshot_days)
    series <- data.table::rbindlist(lapply(c(A="flux_percent",B="growth_percent"),function(metric)
      cbind(meta,analysis="instantaneous",metric=metric,source_accounting="daily_viable_low_daughters",
        interval_type="equal-route protocol IQR; not CI",n_protocol_routes=3L,
        fp_summaries(if(metric=="flux_percent")daily$high_parent_percent else daily$low_parent_percent,3L))))
    profiles <- t$threshold_N==b$config$threshold_N&&t$O2_pct%in%b$config$main_oxygen
    diag <- if(beta==0&&t$threshold_N==66L&&t$initial_ploidy==4&&t$p_misseg==.3&&t$O2_pct%in%c(0,5))
      f7c_diagnostics(op,routes,meta) else NULL
    outputs[[length(outputs)+1L]] <- list(metadata=meta,summary=series,qc=qc,diagnostics=diag,
      source_counts_summary=f7_count_summary(daily,meta,3L),
      population_summary=if(profiles)f7_population_summary(daily,meta,3L) else NULL,
      distribution_summary=if(profiles)f7_distribution_summary(daily,meta,op$l,op$c) else NULL,
      passage_events=cbind(meta,data.table::rbindlist(lapply(routes,`[[`,"events"))),
      route_paths=unname(files),route_sha256=fp_sha(unname(files)))
    rm(op,routes,prefix,daily);gc(verbose=FALSE)
  }
  list(outputs=outputs)
}
f7c_simulate <- function(a) {
  f7c_not_formal(a);b <- f7c_bundle(a)
  ids <- if(!is.null(a$args$task_index))as.integer(a$args$task_index) else
    if(!is.null(a$args$task_ids))as.integer(strsplit(a$args$task_ids,",",fixed=TRUE)[[1]]) else b$tasks$task_id
  stopifnot(length(ids)>0,!anyNA(ids),!anyDuplicated(ids),all(ids%in%b$tasks$task_id))
  for(id in ids) {
    dest <- fp_task_path(a,id)
    if(file.exists(dest)) {
      o <- fp_read(dest);stopifnot(identical(o$fingerprint,b$fingerprint),identical(o$status,"task_validated"))
      receipt <- file.path(a$run,"receipts",sprintf("task_%05d.tsv",id))
      if(!file.exists(receipt))fp_write(data.frame(task_id=id,status=o$status,fingerprint=b$fingerprint,
        sha256=fp_sha(dest),elapsed_seconds=o$elapsed_seconds),receipt)
      stopifnot(data.table::fread(receipt)$sha256==fp_sha(dest))
      message("Reusing validated task ",id);next
    }
    t <- b$tasks[task_id==id];wa <- a
    wa$run <- file.path(a$run,"runtime",sprintf("task_%05d",id))
    tic <- proc.time()[3]
    if(t$model_context=="in vitro") {
      runtime <- f7p_runtime(a$workspace,wa$run)
      o <- f7c_vitro_task(a,b,t,runtime)
    } else {
      fp_engine(wa);o <- f7_simulate_task(wa,b,t)
      o$outputs <- lapply(o$outputs,function(v) {
        v$summary <- v$summary[analysis=="instantaneous"]
        v$passage_events <- NULL
        if(!is.null(v$diagnostics))v$diagnostics <- v$diagnostics[c("heat","budget")]
        v[c("metadata","summary","qc","diagnostics","source_counts_summary","population_summary","distribution_summary")]
      })
    }
    o$fingerprint <- b$fingerprint;o$task <- t;o$status <- "task_validated"
    o$elapsed_seconds <- unname(proc.time()[3]-tic)
    if(!all(vapply(o$outputs,function(v)all(v$qc$passed),TRUE)))stop("Task failed internal validation")
    fp_save(o,dest)
    fp_write(data.frame(task_id=id,status=o$status,fingerprint=b$fingerprint,
      sha256=fp_sha(dest),elapsed_seconds=o$elapsed_seconds),
      file.path(a$run,"receipts",sprintf("task_%05d.tsv",id)))
    message("Completed Figure7 task ",id,"/",nrow(b$tasks)," in ",round(o$elapsed_seconds,1)," s")
  }
  invisible(ids)
}
f7c_collect_tasks <- function(a,b) {
  lapply(b$tasks$task_id,function(id) {
    path <- fp_task_path(a,id);receipt <- file.path(a$run,"receipts",sprintf("task_%05d.tsv",id))
    if(!file.exists(path)||!file.exists(receipt))stop("Missing task/receipt ",id)
    row <- data.table::fread(receipt);o <- fp_read(path)
    stopifnot(nrow(row)==1,row$sha256==fp_sha(path),row$status=="task_validated",
      identical(o$fingerprint,b$fingerprint),identical(o$status,"task_validated"),
      identical(as.list(o$task),as.list(b$tasks[task_id==id])))
    o
  })
}
f7c_grid <- function(a,b,all_tasks) {
  pairs <- list();certificates <- list()
  for(id in unique(b$tasks[model_context=="in vitro"]$base_task_id)) {
    tasks <- b$tasks[base_task_id==id][order(maximum_ploidy)]
    stopifnot(identical(tasks$maximum_ploidy,as.integer(c(137,139,141))))
    for(beta in b$config$intervention_q)for(route in b$policy$routes) {
      rows <- list()
      previous <- NULL
      for(j in seq_len(nrow(tasks))) {
        t <- tasks[j];o <- all_tasks[[t$task_id]]
        v <- o$outputs[[match(beta,vapply(o$outputs,function(v)v$metadata$q,0.0))]]
        path <- f7c_route_path(a,t$task_id,beta,route)
        stopifnot(fp_sha(path)==v$route_sha256[match(path,v$route_paths)])
        r <- fp_read(path);stopifnot(identical(r$fingerprint,b$fingerprint))
        if(!is.null(previous))rows[[length(rows)+1L]] <- cbind(base_task_id=id,
          f7p_compare(previous,r,b$policy$grid_limits))
        previous <- r
      }
      zz <- data.table::rbindlist(rows);pairs[[length(pairs)+1L]] <- zz
      certificates[[length(certificates)+1L]] <- data.table::data.table(base_task_id=id,q=beta,route=route,
        accepted_maximum_ploidy=141L,adjacent_comparisons=nrow(zz),passed=nrow(zz)==2L&&all(zz$passed))
      rm(previous,r);gc(verbose=FALSE)
    }
  }
  comparisons <- data.table::rbindlist(pairs);cert <- data.table::rbindlist(certificates)
  fp_write(comparisons,file.path(a$run,"provenance/grid_comparisons.tsv"))
  fp_write(cert,file.path(a$run,"provenance/grid_certificates.tsv"))
  if(!all(cert$passed))stop("Grid acceptance failed; keep all results isolated; do not publish or expand beyond 141N automatically")
  invisible(cert)
}
f7c_aggregate <- function(a) {
  f7c_not_formal(a);b <- f7c_bundle(a);all_tasks <- f7c_collect_tasks(a,b)
  f7c_grid(a,b,all_tasks)
  accepted <- b$tasks[model_context=="in vivo"|maximum_ploidy==141L]
  outputs <- unlist(lapply(all_tasks[accepted$task_id],`[[`,"outputs"),recursive=FALSE)
  for(spec in list(c("summary","instant_time_summary.tsv"),
    c("source_counts_summary","daily_source_counts_summary.tsv"),
    c("population_summary","low_population_time_summary.tsv"),
    c("distribution_summary","low_population_distribution_summary.tsv")))
    fp_write(data.table::rbindlist(lapply(outputs,`[[`,spec[1]),fill=TRUE),file.path(a$run,"tables",spec[2]))
  fp_write(data.table::rbindlist(lapply(outputs,`[[`,"qc"),fill=TRUE),
    file.path(a$run,"provenance/task_numerical_validation.tsv"))
  diagnostics <- Filter(Negate(is.null),lapply(outputs,`[[`,"diagnostics"))
  for(spec in list(c("heat","source_daughter_flux_binned.tsv"),c("budget","winner_low_population_budget.tsv")))
    fp_write(data.table::rbindlist(lapply(diagnostics,`[[`,spec[1]),fill=TRUE),file.path(a$run,"tables",spec[2]))
  events <- data.table::rbindlist(lapply(outputs,`[[`,"passage_events"),fill=TRUE)
  fp_write(events,file.path(a$run,"provenance/passage_events.tsv"))
  fp_write(accepted,file.path(a$run,"provenance/accepted_task_manifest.tsv"))
  fp_write(data.frame(status="aggregated_not_yet_validated",profile=a$profile,fingerprint=b$fingerprint),
    file.path(a$run,"provenance/completion.tsv"))
  invisible(b)
}
f7c_validate <- function(a) {
  f7c_not_formal(a);b <- f7c_bundle(a);checks <- list()
  add <- function(name,pass)checks[[length(checks)+1L]] <<- data.table::data.table(test=name,passed=isTRUE(pass))
  parameter <- data.table::fread(file.path(a$run,"provenance/full_native_parameter_passthrough.tsv"))
  add("complete_native_parameter_passthrough",all(parameter$passed))
  add("context_parameter_counts",all(b$endpoints$n_parameters==ifelse(b$endpoints$model_context=="in vitro",33L,26L)))
  add("authoritative_culture_config",all(b$endpoints[model_context=="in vitro"]$config_source=="fit_result.rds$ctx$invitro_cfg"))
  templates <- data.table::fread(file.path(a$run,"provenance/native_last5_templates.tsv"))
  add("three_routes_last_five_complete_tuples",all(templates[,.(n=.N),by=.(pair_label,initial_ploidy,route)]$n==5L)&&all(templates$complete))
  tasks <- data.table::fread(file.path(a$run,"provenance/accepted_task_manifest.tsv"))
  grid <- data.table::fread(file.path(a$run,"provenance/grid_certificates.tsv"))
  add("all_invitro_beta_route_grid_certificates",nrow(grid)==sum(tasks$model_context=="in vitro")*9L&&all(grid$passed))
  z <- data.table::fread(file.path(a$run,"tables/instant_time_summary.tsv"))
  add("complete_condition_beta_metric_day_grid",nrow(z)==nrow(tasks)*3L*2L*1001L)
  add("unique_condition_beta_metric_day",!anyDuplicated(z[,.(base_task_id,q,metric,day)]))
  add("daily_viable_mother_source_definition",all(z$source_accounting=="daily_viable_low_daughters"))
  A <- z[metric=="flux_percent"];B <- z[metric=="growth_percent"]
  data.table::setorder(A,base_task_id,q,day);data.table::setorder(B,base_task_id,q,day)
  ok <- is.finite(A$median)
  add("paired_validity",identical(A$n_valid,B$n_valid)&&identical(is.na(A$median),is.na(B$median)))
  add("daily_source_complement",all(abs(A$median[ok]+B$median[ok]-100)<1e-9))
  add("reversed_IQR_complement",all(abs(A$q25[ok]+B$q75[ok]-100)<1e-9)&&all(abs(A$q75[ok]+B$q25[ok]-100)<1e-9))
  add("day_zero_undefined",all(is.na(z[day==0]$median)))
  add("shares_range",all(z[is.finite(median)]$median>=-1e-9&z[is.finite(median)]$median<=100+1e-9))
  add("quantile_order",all(z[is.finite(median)]$q25<=z[is.finite(median)]$median&
    z[is.finite(median)]$median<=z[is.finite(median)]$q75))
  add("three_real_protocol_routes_not_200_repeats",all(z[model_context=="in vitro"]$n_repeats==3L)&&
    all(z[model_context=="in vitro"]$n_protocol_routes==3L))
  qc <- data.table::fread(file.path(a$run,"provenance/task_numerical_validation.tsv"))
  add("every_accepted_task_internal_validation",all(qc$passed))
  ev <- data.table::fread(file.path(a$run,"provenance/passage_events.tsv"))
  add("no_native_zero_mass_fallback",!any(ev$uniform_fallback))
  add("route_elapsed_1000_days",all(ev[,.(days=sum(retained_days)),by=.(task_id,q,route)]$days==1000L))
  add("full_window_before_final_truncation",all(ev$retained_days<=ev$selected_day&ev$selected_day<=ev$duration_days))
  pop <- data.table::fread(file.path(a$run,"tables/low_population_time_summary.tsv"))
  add("five_explicit_origin_groups",setequal(unique(pop$group),c("all_low","direct_first_generation","other_low","historical_lineage","never_HL_lineage")))
  add("empty_means_undefined",all(is.na(pop[n_valid==0]$median)))
  add("source_ploidy_below_threshold",all(pop[metric=="mean_ploidy"&is.finite(median)]$median<3+1e-10))
  dist <- data.table::fread(file.path(a$run,"tables/low_population_distribution_summary.tsv"))
  add("cell_weighted_CDF_range",all(dist[is.finite(cdf)]$cdf>=-1e-9&dist[is.finite(cdf)]$cdf<=1+1e-9))
  add("empty_CDF_undefined",all(is.na(dist[n_valid_repeats==0]$cdf)))
  for(name in c("source_daughter_flux_binned.tsv","winner_low_population_budget.tsv")) {
    d <- data.table::fread(file.path(a$run,"tables",name))
    add(paste0(name,"_all_context_family_oxygen_anchors"),nrow(unique(d[,.(pair_label,model_context,O2_pct)]))==8L)
  }
  out <- data.table::rbindlist(checks)
  fp_write(out,file.path(a$run,"provenance/numerical_validation.tsv"))
  if(!all(out$passed))stop("Figure7 validation failed: ",paste(out[passed==FALSE]$test,collapse=", "))
  files <- list.files(file.path(a$run,"tables"),pattern="\\.tsv$",full.names=TRUE)
  stopifnot(!any(grepl("history",basename(files))))
  fp_write(data.table::data.table(relative_path=substring(files,nchar(a$data_root)+2L),
    logical_name=basename(files),sha256=fp_sha(files)),file.path(a$run,"provenance/plot_inputs.tsv"))
  fp_write(data.frame(status="validated_complete",profile=a$profile,fingerprint=b$fingerprint,
    numerical_checks=nrow(out),production_ensemble=a$profile=="full",engine=b$engine,
    historical_product_writes=FALSE,figure8_changes=FALSE,
    bundle_sha256=fp_sha(file.path(a$run,"prepared/analysis_bundle.rds")),
    validation_sha256=fp_sha(file.path(a$run,"provenance/numerical_validation.tsv")),
    plot_manifest_sha256=fp_sha(file.path(a$run,"provenance/plot_inputs.tsv"))),file.path(a$run,"provenance/completion.tsv"))
  message("PASS Figure7: ",nrow(out)," validation gates")
  invisible(out)
}
f7c_plot <- function(a) {
  if(a$profile!="full")stop("Representative outputs cannot be published")
  args <- c(shQuote(file.path(a$code,"draw_Figure7.R")),"--analysis=instantaneous",
    paste0("--workspace=",shQuote(a$workspace)),paste0("--data-root=",shQuote(a$data_root)),
    paste0("--output-root=",shQuote(a$output_root)),paste0("--formats=",a$args$formats %||% "pdf,png"))
  if(system2(file.path(R.home("bin"),"Rscript"),args)!=0L)stop("Figure7 current-only rendering failed")
}
f7c_cli <- function(script) {
  args <- fp_args()
  if(identical(args$help,"yes")) {
    cat("Figure7 current native-fit passage pipeline; historical/Figure8 products never written.\n",
      "--stage=plan|prepare|simulate|aggregate|validate|plot|all --profile=full|representative\n",
      "--data-root=ISOLATED_CAMPAIGN_DATA --output-root=ISOLATED_FIGURES --task-index=N --task-ids=N,M\n",
      "In vitro: 137/139/141N grid certificates, three last-five protocol routes, no stock gate.\n")
    return(invisible(NULL))
  }
  a <- fp_init(7L,script,args)
  if(args$stage=="plan") {
    fp_preflight(a);tasks <- f7c_tasks(a)
    print(tasks[,.(workers=.N,conditions=data.table::uniqueN(base_task_id)),by=.(model_context,engine)])
    return(invisible(tasks))
  }
  stages <- if(args$stage=="all")c("prepare","simulate","aggregate","validate",if(a$profile=="full")"plot") else args$stage
  if(any(!stages%in%c("prepare","simulate","aggregate","validate","plot")))stop("Unknown stage")
  for(stage in stages)switch(stage,prepare=f7c_prepare(a),simulate=f7c_simulate(a),
    aggregate=f7c_aggregate(a),validate=f7c_validate(a),plot=f7c_plot(a))
  invisible(a)
}
