# Current Figure7/8 forward pipeline. No fit, prediction or code is read from an
# archived Figure7/8 run. Joint fits and the canonical model are upstream inputs.
.fp_module_directory <- dirname(normalizePath(sys.frame(1)$ofile,mustWork=TRUE))
source(file.path(.fp_module_directory,"figure78_paths.R"),local=globalenv())
`%||%` <- function(x,y) if(is.null(x)) y else x

fp_args <- function(x=commandArgs(TRUE)) {
  out <- list(stage="all",profile="full")
  for(s in x) {
    if(!grepl("^--[^=]+=.+$",s)) stop("Use --option=value: ",s)
    p <- regexpr("=",s,fixed=TRUE)
    out[[gsub("-","_",substring(s,3L,p-1L))]] <- substring(s,p+1L)
  }
  known <- c("stage","profile","workspace","data_root","output_root","fit_root",
    "model_root","task_index","task_ids","rebuild","formats","help")
  if(length(setdiff(names(out),known)))stop("Unknown option: ",paste(setdiff(names(out),known),collapse=", "))
  out
}
fp_safe <- function(path) {
  path <- path.expand(path)
  if(any(grepl("(^|/)(arxives|archives|frozen_sources)(/|$)",path)))
    stop("Archived/retired sources are not pipeline inputs: ",paste(path,collapse=", "))
  existing <- path[file.exists(path)]
  if(length(existing)&&any(grepl("(^|/)(arxives|archives|frozen_sources)(/|$)",normalizePath(existing))))
    stop("An input resolves into an archive")
  invisible(path)
}
fp_sha <- function(paths) {
  fp_safe(paths)
  if(any(!file.exists(paths)))stop("Missing required input: ",paste(paths[!file.exists(paths)],collapse=", "))
  unname(vapply(paths,digest::digest,"",file=TRUE,algo="sha256"))
}
fp_write <- function(x,path) {
  fp_safe(path);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(".table_",tmpdir=dirname(path),fileext=if(grepl("gz$",path))".tsv.gz" else ".tsv")
  data.table::fwrite(x,tmp,sep="\t",na="NA",compress=if(grepl("gz$",path))"gzip" else "none")
  if(!file.rename(tmp,path))stop("Atomic table write failed: ",path)
  invisible(path)
}
fp_save <- function(x,path) {
  fp_safe(path);dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- tempfile(".bundle_",tmpdir=dirname(path));saveRDS(x,tmp,compress=FALSE)
  if(!file.rename(tmp,path))stop("Atomic RDS write failed: ",path)
  invisible(path)
}
fp_read <- function(path) {fp_safe(path);readRDS(path)}
fp_files <- function(a) file.path(a$code,c(paste0("run_Figure",a$figure,".R"),
  paste0("data_Figure",a$figure,".R"),"util/analysis/figure78_pipeline.R",
  "util/analysis/figure78_model.R","util/analysis/figure78_culture.cpp",
  "util/analysis/figure78_paths.R",
  paste0("util/analysis/figure",a$figure,"_simulation.R"),
  if(a$figure==7L)c("util/analysis/figure7_source_accounting.R",
    paste0("util/analysis/figure7_",c("continuous","history","direct","daily"),".cpp")) else
    "util/analysis/figure8_empirical.R"))

fp_init <- function(figure,script,args) {
  suppressPackageStartupMessages({library(data.table);library(digest);library(jsonlite);library(Matrix);library(Rcpp)})
  data.table::setDTthreads(2L)
  workspace <- normalizePath(args$workspace %||% file.path(dirname(script),"../.."),mustWork=TRUE)
  cfg_path <- file.path(workspace,"Code/config",paste0("figure",figure,".json"))
  cfg <- jsonlite::fromJSON(cfg_path,simplifyVector=TRUE)
  shared <- jsonlite::fromJSON(file.path(workspace,"Code/config/figure78_inputs.json"))
  profile <- args$profile
  if(!profile%in%c("full","representative"))stop("Profile must be full or representative")
  data_root <- f78_current_path(args$data_root %||% "data/Figures",workspace)
  if(profile!="full"&&identical(normalizePath(data_root,mustWork=FALSE),
    normalizePath(file.path(workspace,"data/Figures"),mustWork=TRUE)))
    stop("Representative tests require a separate --data-root; formal data are protected")
  run <- file.path(data_root,paste0("Figure",figure));fp_safe(c(workspace,data_root,run))
  fit_root <- f78_fit_root(workspace,shared,args$fit_root)
  model_root <- f78_model_root(workspace,args$model_root)
  f78_locked_model_manifest(workspace,model_root)
  fp_safe(c(fit_root,model_root))
  a <- list(figure=figure,workspace=workspace,code=file.path(workspace,"Code/Figures"),
    config_path=cfg_path,cfg=cfg,shared=shared,profile=profile,run=run,data_root=data_root,
    fit_root=normalizePath(fit_root,mustWork=TRUE),model_root=normalizePath(model_root,mustWork=TRUE),
    output_root=args$output_root %||% if(identical(normalizePath(data_root,mustWork=FALSE),
      normalizePath(file.path(workspace,"data/Figures"),mustWork=TRUE)))workspace else file.path(dirname(data_root),"rendered"),args=args)
  a$replicates <- if(profile=="full")cfg$replicates else cfg$representative$replicates
  a
}

fp_model_manifest <- function(a) {
  m <- f78_locked_model_manifest(a$workspace,a$model_root)
  data.table::data.table(path=m$path,sha256=m$sha256,role="canonical model source")
}
fp_endpoint_spec <- function(a) {
  out <- lapply(names(a$shared$families),function(fam) {
    f <- a$shared$families[[fam]]
    seeds <- if(a$figure==7L)f$winner_seed else if(a$profile=="full")seq_len(a$cfg$seeds_per_family) else
      unique(c(1L,f$winner_seed))
    data.table::CJ(pair_label=fam,seed_number=as.integer(seeds),model_context=c("in vivo","in vitro"))[,pair_id:=f$pair_id]
  })
  z <- data.table::rbindlist(out);data.table::setorder(z,pair_label,seed_number,model_context)
  z[,prepared_index:=seq_len(.N)];z
}
fp_required_inputs <- function(a,ep=fp_endpoint_spec(a)) {
  dirs <- unique(file.path(a$fit_root,ep$pair_id,paste0("seed",ep$seed_number)))
  paths <- c(as.vector(outer(dirs,c("joint_best_params_long.tsv","fit_config.rds","fit_result.rds"),file.path)),
    file.path(a$fit_root,a$shared$families$C01$pair_id,"seed1/fit_result.rds"))
  if(a$figure==8L)paths <- c(paths,file.path(a$workspace,"data/share_bundle_figures_3_5",
    c("data/Figure3_PCAWG_included_samples.csv","data/Figure3_DepMap_included_lines.csv","originals/organ_lut.csv")))
  unique(paths)
}
fp_preflight <- function(a) {
  paths <- fp_required_inputs(a);fp_safe(paths)
  missing <- paths[!file.exists(paths)]
  if(length(missing))stop("Raw inputs missing (no fallback):\n",paste(missing,collapse="\n"))
  fp_model_manifest(a)
  invisible(paths)
}

fp_trace_schedule <- function(ids,parents,terminal) {
  path <- character();current <- terminal
  repeat {
    if(current%in%path)stop("Cycle in experimental passage tree")
    path <- c(current,path);i <- match(current,ids)
    if(is.na(i))stop("Missing passage segment: ",current)
    p <- parents[i];if(is.na(p)||!nzchar(p))break;current <- p
  }
  path
}
fp_passage_rule <- function(a) {
  # Experimental adapters, not an old Figure6/7/8 prediction or schedule file.
  path <- file.path(a$fit_root,a$shared$families$C01$pair_id,"seed1/fit_result.rds")
  fit <- fp_read(path);vitro <- fit$best_components$invitro;rows <- list()
  for(cohort in c("2N","4N")) {
    run <- vitro[[paste0("run_",cohort)]];segments <- run$adapter$segments
    if(!length(segments))stop("No experimental passage adapters in ",path)
    ids <- vapply(segments,`[[`,"","segment_id")
    parents <- vapply(segments,function(x)as.character(x$parent_segment_id %||% NA_character_)[1],"")
    terminals <- setdiff(ids,stats::na.omit(parents))
    group <- function(x)if(!grepl("::",x,fixed=TRUE))"control" else if(grepl("::O1$",x))"O1" else if(grepl("::O2$",x))"O2" else NA_character_
    for(g in c("control","O1","O2")) {
      candidate <- terminals[vapply(terminals,group,"")==g]
      candidate <- candidate[!is.na(candidate)]
      if(!length(candidate))stop("Missing experimental lineage ",cohort,"/",g)
      paths <- lapply(candidate,fp_trace_schedule,ids=ids,parents=parents)
      duration <- vapply(paths,function(p)sum(vapply(segments[match(p,ids)],`[[`,0.0,"duration_days")),0.0)
      k <- order(-duration,candidate)[1];p <- paths[[k]];elapsed <- 0
      for(j in seq_along(p)) {
        i <- match(p[j],ids);s <- segments[[i]];selection <- run$segment_results[[i]]$selection
        d <- as.numeric(s$duration_days)
        rows[[length(rows)+1L]] <- data.table::data.table(cohort=cohort,lineage=paste(cohort,g,sep="_"),
          passage_index=j,segment_id=s$segment_id,experimental_day_start=elapsed,
          experimental_day_end=elapsed+d,passage_duration=d,initial_cells=as.numeric(s$initial_cells),
          target_live_cells=as.numeric(selection$target_live_cells))
        elapsed <- elapsed+d
      }
    }
  }
  raw <- data.table::rbindlist(rows)
  valid <- raw[is.finite(initial_cells)&initial_cells>0&is.finite(target_live_cells)&target_live_cells>initial_cells]
  gm <- function(x)exp(mean(log(x)))
  lineages <- valid[,.(n_valid_passage=.N,seed_geomean=gm(initial_cells),
    ratio_geomean=gm(target_live_cells/initial_cells),duration_median=median(passage_duration)),by=.(cohort,lineage)]
  stopifnot(nrow(lineages)==6L,nrow(valid)==114L)
  seed <- gm(lineages$seed_geomean);ratio <- gm(lineages$ratio_geomean)
  rule <- data.table::data.table(seed_cells=round(seed),target_live_cells=ceiling(seed*ratio),
    seed_cells_before_integerization=seed,target_live_cells_before_integerization=seed*ratio,
    expansion_ratio=ratio,n_source_lineage=6L,n_source_maintenance_passage=nrow(valid),
    representative_passage_duration_day=median(lineages$duration_median),
    trigger="first integer day after growth at/above threshold",sampling="systematic integerization; hypergeometric without replacement")
  if(rule$seed_cells!=569799||rule$target_live_cells!=6180783)
    stop("Raw experimental passage rule differs from the approved protocol; do not silently substitute")
  list(rule=rule,raw=raw,lineages=lineages,path=path)
}

fp_load_parameters <- function(a,ep) {
  prepared <- audit <- manifest <- list()
  for(i in seq_len(nrow(ep))) {
    e <- ep[i];dir <- file.path(a$fit_root,e$pair_id,paste0("seed",e$seed_number))
    parfile <- file.path(dir,"joint_best_params_long.tsv")
    culture <- e$model_context=="in vitro"
    configfile <- file.path(dir,if(culture)"fit_result.rds" else "fit_config.rds")
    if(culture)cfg <- fp_read(configfile)$ctx$invitro_cfg else cfg <- fp_read(configfile)
    if(!is.list(cfg))stop("Missing context-specific fit configuration: ",configfile)
    scope <- if(culture)"invitro_effective" else "shared_invivo"
    native <- read_joint_scope_params(parfile,scope)
    if(length(native)!=(if(culture)33L else 26L) || anyDuplicated(names(native)) || any(!is.finite(native)))
      stop("Invalid complete native parameter vector: ",parfile," / ",scope)
    rp <- prepare_run_params(native,if(culture)"invitro" else "joint",cfg,5)
    simcfg <- prepare_sim_cfg(cfg,list(),5,rp)
    after <- vapply(names(native),function(n)as.numeric(rp[[n]]),0.0)
    allow <- names(native)=="o2_S0"
    if(any(abs(after[!allow]-native[!allow])>1e-12))stop("Native fitted parameter changed during preparation")
    if(culture&&!all(vapply(c("DT","Crowding","K","crowding"),function(n)identical(simcfg[[n]],cfg[[n]]),TRUE)))
      stop("Culture configuration was overwritten")
    cfg_origin <- if(culture)"fit_result.rds$ctx$invitro_cfg" else "fit_config.rds"
    expected_grid <- unlist(a$cfg$model_grid,use.names=TRUE)
    if(any(vapply(names(expected_grid),function(n)as.numeric(simcfg[[n]]),0.0)!=expected_grid))
      stop("Upstream chromosome grid differs from the approved prediction grid")
    prepared[[i]] <- list(run_params=rp,config=simcfg,native_parameters=native,
      model_context=e$model_context,original_fit_config=cfg,config_source=cfg_origin)
    signature <- digest::digest(list(native=native,config=cfg),algo="sha256")
    ep[i,`:=`(parameter_scope=scope,n_parameters=length(native),complete_input_signature=signature,
      source_parameter_path=parfile,source_config_path=configfile,config_source=cfg_origin)]
    audit[[i]] <- data.table::data.table(pair_label=e$pair_label,seed_number=e$seed_number,model_context=e$model_context,
      parameter=names(native),fitted_value=as.numeric(native),prepared_value=after,
      allowed_override=allow,passed=allow|abs(after-native)<1e-12)
    manifest[[i]] <- data.table::data.table(path=c(parfile,configfile),sha256=fp_sha(c(parfile,configfile)),role=c(scope,cfg_origin))
    if(i%%50L==0L)message("Prepared ",i," complete context-specific seed inputs")
  }
  list(prepared=prepared,endpoints=ep,audit=data.table::rbindlist(audit),
    manifest=unique(data.table::rbindlist(manifest)))
}

fp_prepare <- function(a) {
  dest <- file.path(a$run,"prepared/analysis_bundle.rds")
  if(file.exists(dest)&&!identical(a$args$rebuild,"yes"))return(invisible(fp_bundle(a)))
  # Never overwrite preserved published results as a side effect of preparation.
  entries <- if(dir.exists(a$run))list.files(a$run,all.files=TRUE,no..=TRUE) else character()
  cache_only <- all(entries%in%c(".cache","prepared"))&&
    all(list.files(file.path(a$run,"prepared"),all.files=TRUE,no..=TRUE)%in%"model_source")
  if(length(entries)&&!cache_only) {
    if(!identical(a$args$rebuild,"yes"))stop("Existing retained data are protected. Use a separate --data-root or explicit --rebuild=yes")
    backup <- file.path(a$workspace,"audit",paste0("figure",a$figure,"_rebuild_",format(Sys.time(),"%Y%m%d_%H%M%S")),"before_data",paste0("Figure",a$figure))
    dir.create(dirname(backup),recursive=TRUE,showWarnings=FALSE)
    if(!file.rename(a$run,backup))stop("Could not move existing data to recoverable audit backup: ",backup)
    message("Previous data moved, not deleted: ",backup)
  }
  fp_preflight(a)
  dir.create(a$run,recursive=TRUE,showWarnings=FALSE)
  fp_engine(a);p <- fp_load_parameters(a,fp_endpoint_spec(a));r <- fp_passage_rule(a)
  tasks <- if(a$figure==7L)f7_tasks(a) else f8_tasks(a)
  inputs <- unique(data.table::rbindlist(list(p$manifest,fp_model_manifest(a),
    data.table::data.table(path=r$path,sha256=fp_sha(r$path),role="raw passage adapters"))))
  if(a$figure==8L) {
    emp <- f8_empirical(a)
    inputs <- unique(data.table::rbindlist(list(inputs,emp$manifest)))
  }
  code <- c(fp_files(a),a$config_path,file.path(a$workspace,"Code/config/figure78_inputs.json"),
    file.path(a$workspace,"Code/config/manifests/figure78_model_sources.tsv"))
  codes <- data.table::data.table(path=code,sha256=fp_sha(code))
  fingerprint <- digest::digest(list(config=a$cfg,shared=a$shared,profile=a$profile,replicates=a$replicates,
    codes=codes$sha256,inputs=inputs,endpoints=p$endpoints,tasks=tasks,rule=r$rule),algo="sha256")
  b <- list(figure=a$figure,config=a$cfg,shared=a$shared,profile=a$profile,replicates=a$replicates,
    prepared=p$prepared,endpoints=p$endpoints,tasks=tasks,rule=r$rule,fingerprint=fingerprint,
    code_manifest=codes,input_manifest=inputs,model_manifest=fp_model_manifest(a))
  fp_write(p$endpoints,file.path(a$run,"prepared/seed_input_manifest.tsv"))
  fp_write(p$audit,file.path(a$run,"provenance/full_native_parameter_passthrough.tsv"))
  fp_write(tasks,file.path(a$run,"prepared/task_manifest.tsv"))
  fp_write(inputs,file.path(a$run,"provenance/input_manifest.tsv"))
  fp_write(codes,file.path(a$run,"provenance/current_code_manifest.tsv"))
  fp_write(r$raw,file.path(a$run,"provenance/source_six_lineage_passage_schedule.tsv"))
  fp_write(r$lineages,file.path(a$run,"provenance/canonical_passage_lineage_summary.tsv"))
  fp_write(r$rule,file.path(a$run,"provenance/canonical_passage_rule.tsv"))
  fp_write(data.frame(session=capture.output(sessionInfo())),file.path(a$run,"provenance/runtime_session.tsv"))
  fp_save(b,dest)
  fp_write(data.frame(figure=a$figure,profile=a$profile,fingerprint=fingerprint,status="prepared",
    n_tasks=nrow(tasks),n_context_seed_inputs=nrow(p$endpoints)),file.path(a$run,"provenance/preparation.tsv"))
  message("Prepared Figure",a$figure,": ",nrow(tasks)," tasks; ",nrow(p$endpoints)," full seed/context inputs")
  invisible(b)
}
fp_bundle <- function(a,verify_inputs=TRUE) {
  b <- fp_read(file.path(a$run,"prepared/analysis_bundle.rds"))
  if(!identical(b$config,a$cfg)||!identical(b$shared,a$shared)||!identical(b$profile,a$profile))
    stop("Prepared configuration differs; use a separate data root for a different computation")
  if(!identical(b$code_manifest$sha256,fp_sha(b$code_manifest$path)))stop("Current numerical code changed since prepare")
  if(verify_inputs&&!identical(b$input_manifest$sha256,fp_sha(b$input_manifest$path)))stop("An upstream scientific input changed")
  b
}
fp_task_path <- function(a,id)file.path(a$run,"tasks",sprintf("task_%05d.rds",as.integer(id)))
fp_simulate <- function(a) {
  b <- fp_bundle(a);fp_engine(a)
  ids <- if(!is.null(a$args$task_index))as.integer(a$args$task_index) else if(!is.null(a$args$task_ids))
    as.integer(strsplit(a$args$task_ids,",",fixed=TRUE)[[1]]) else b$tasks$task_id
  if(!length(ids)||anyNA(ids)||anyDuplicated(ids)||any(!ids%in%b$tasks$task_id))stop("Invalid task IDs")
  for(id in ids) {
    dest <- fp_task_path(a,id)
    if(file.exists(dest)) {
      old <- fp_read(dest)
      if(!identical(old$fingerprint,b$fingerprint))stop("Incompatible task cache: ",dest)
      next
    }
    t <- b$tasks[task_id==id]
    o <- if(a$figure==7L)f7_simulate_task(a,b,t) else f8_simulate_task(a,b,t,dest)
    o$fingerprint <- b$fingerprint;o$task <- t
    fp_save(o,dest)
    checkpoint <- sub(".rds$","_checkpoint.rds",dest)
    if(file.exists(checkpoint))unlink(checkpoint)
    message("Completed Figure",a$figure," task ",id,"/",nrow(b$tasks))
  }
  invisible(ids)
}
fp_summaries <- function(m,nominal) {
  get <- function(p)vapply(seq_len(ncol(m)),function(i) {
    y <- m[,i];if(!any(is.finite(y)))NA_real_ else unname(quantile(y[is.finite(y)],p,type=7))},0.0)
  data.table::data.table(day=seq_len(ncol(m))-1L,median=get(.5),q25=get(.25),q75=get(.75),
    n_valid=colSums(is.finite(m)),n_repeats=nrow(m),nominal_repeats=nominal)
}
fp_quantiles <- function(z,keys) {
  z[,.(median=median(mean_ploidy),q25=unname(quantile(mean_ploidy,.25,type=7)),
    q75=unname(quantile(mean_ploidy,.75,type=7)),n_seeds=.N,
    n_c01=sum(pair_label=="C01"),n_c02=sum(pair_label=="C02"),
    unique_seed_ids=data.table::uniqueN(paste(pair_label,seed_number))),by=keys]
}
fp_aggregate <- function(a) {
  b <- fp_bundle(a);paths <- vapply(b$tasks$task_id,function(id)fp_task_path(a,id),"")
  if(any(!file.exists(paths)))stop("Missing ",sum(!file.exists(paths))," tasks; aggregation is not complete")
  tasks <- lapply(paths,function(path) {
    o <- fp_read(path);if(!identical(o$fingerprint,b$fingerprint))stop("Task provenance mismatch: ",path)
    # Discard large repeat/trajectory arrays before reading the next task.
    if(a$figure==7L)list(outputs=lapply(o$outputs,function(v)c(v[c("metadata","summary","qc","diagnostics",
        "source_counts_summary","population_summary","distribution_summary")],
      list(passage_events=v$trajectories$passage_events)))) else
      list(data=o$data)
  })
  if(a$figure==7L)f7_aggregate(a,b,tasks) else f8_aggregate(a,b,tasks)
  fp_write(data.frame(status="aggregated_not_yet_validated",profile=a$profile,fingerprint=b$fingerprint),
    file.path(a$run,"provenance/completion.tsv"))
  invisible(b)
}
fp_validate <- function(a) {
  b <- fp_bundle(a);checks <- list()
  check <- function(name,pass,value=NA_real_)checks[[length(checks)+1L]] <<-
    data.table::data.table(test=name,value=value,passed=isTRUE(pass))
  pass <- data.table::fread(file.path(a$run,"provenance/full_native_parameter_passthrough.tsv"))
  check("every_complete_native_parameter_passed_through",all(pass$passed))
  ep <- b$endpoints
  check("context_parameter_counts",all(ep$n_parameters==ifelse(ep$model_context=="in vivo",26L,33L)))
  check("authoritative_culture_configuration",all(ep[model_context=="in vitro"]$config_source=="fit_result.rds$ctx$invitro_cfg"))
  for(i in seq_along(b$prepared))if(ep$model_context[i]=="in vitro") {
    pr <- b$prepared[[i]]
    check(paste0("culture_config_passthrough_",i),all(vapply(c("DT","K","Crowding","crowding"),
      function(n)identical(pr$config[[n]],pr$original_fit_config[[n]]),TRUE)))
  }
  if(a$figure==7L)f7_validate(a,b,check) else f8_validate(a,b,check)
  out <- data.table::rbindlist(checks);fp_write(out,file.path(a$run,"provenance/numerical_validation.tsv"))
  if(!all(out$passed))stop("Numerical validation failed: ",paste(out[passed==FALSE]$test,collapse=", "))
  # Register newly computed output bytes; no old-code hash is forged or disabled.
  files <- list.files(file.path(a$run,"tables"),pattern="\\.tsv(\\.gz)?$",full.names=TRUE)
  files <- c(files,list.files(file.path(a$run,"seed_data"),pattern="\\.tsv.gz$",full.names=TRUE))
  fp_write(data.table::data.table(relative_path=substring(files,nchar(a$data_root)+2L),
    logical_name=basename(files),sha256=fp_sha(files)),file.path(a$run,"provenance/plot_inputs.tsv"))
  fp_write(data.frame(status="validated_complete",profile=a$profile,fingerprint=b$fingerprint,
    numerical_checks=nrow(out),production_ensemble=a$profile=="full",
    bundle_sha256=fp_sha(file.path(a$run,"prepared/analysis_bundle.rds")),
    validation_sha256=fp_sha(file.path(a$run,"provenance/numerical_validation.tsv")),
    plot_manifest_sha256=fp_sha(file.path(a$run,"provenance/plot_inputs.tsv"))),
    file.path(a$run,"provenance/completion.tsv"))
  message("PASS Figure",a$figure,": ",nrow(out)," numerical checks; profile=",a$profile)
  invisible(out)
}
fp_plot <- function(a) {
  if(a$profile!="full")stop("Representative data cannot masquerade as 500+500 publication results; inspect the acceptance report")
  script <- file.path(a$code,paste0("draw_Figure",a$figure,".R"))
  args <- c(shQuote(script),paste0("--workspace=",shQuote(a$workspace)),paste0("--data-root=",shQuote(a$data_root)),
    paste0("--output-root=",shQuote(a$output_root)),paste0("--formats=",a$args$formats %||% "pdf,png"))
  if(a$figure==7L)args <- c(args,"--analysis=all")
  status <- system2(file.path(R.home("bin"),"Rscript"),args)
  if(status!=0L)stop("Plot stage failed")
}
fp_cli <- function(figure,script) {
  args <- fp_args()
  if(identical(args$help,"yes")) {
    cat("Figure",figure,"complete forward pipeline\n",
      "--stage=plan|prepare|simulate|aggregate|validate|plot|all\n",
      "--profile=full|representative --data-root=DIR --output-root=DIR\n",
      "--task-index=N (one manifest task) --task-ids=N,M --fit-root=DIR --model-root=DIR\n",
      "No HPC submission. Existing retained data require separate data root or explicit --rebuild=yes.\n")
    return(invisible(NULL))
  }
  a <- fp_init(figure,script,args)
  if(args$stage=="plan") {
    paths <- fp_preflight(a);tasks <- if(figure==7L)f7_tasks(a) else f8_tasks(a)
    print(data.frame(figure=figure,profile=a$profile,tasks=nrow(tasks),
      full_context_seed_inputs=nrow(fp_endpoint_spec(a)),replicates=a$replicates,raw_files=length(paths),
      missing_files=sum(!file.exists(paths)),hpc_jobs_submitted=0L))
    return(invisible(tasks))
  }
  if(!args$stage%in%c("prepare","simulate","aggregate","validate","plot","all"))stop("Unknown stage")
  stages <- if(args$stage=="all")c("prepare","simulate","aggregate","validate",if(a$profile=="full")"plot") else args$stage
  for(stage in stages)switch(stage,prepare=fp_prepare(a),simulate=fp_simulate(a),
    aggregate=fp_aggregate(a),validate=fp_validate(a),plot=fp_plot(a))
  invisible(a)
}
