# Formal native DEoptim/L-BFGS-B control against the frozen tumor objective.
new_evaluator <- function(ctx,out,fixed=NULL,method="evaluation",deadline=Inf) {
  e <- new.env(parent=emptyenv()); e$n <- 0L; e$best <- Inf; e$seconds <- 0; e$improvements <- 0L; e$received <- list()
  free <- setdiff(names(ctx$endpoint),names(fixed))
  pidout <- file.path(out,"evaluators",paste0("pid_",Sys.getpid()))
  dir.create(pidout,recursive=TRUE,showWarnings=FALSE)
  .ivi_state <<- e
  fn <- function(x) {
    if(as.numeric(Sys.time())>=deadline) stop("IVI_SOFT_TIME_LIMIT")
    p <- ctx$endpoint; p[free] <- as.numeric(x); if(length(fixed)) p[names(fixed)] <- fixed
    stopifnot(all(p[names(fixed)]==fixed))
    if(nzchar(Sys.getenv("CONT_TOKEN_PATH"))) take_token(Sys.getenv("CONT_TOKEN_PATH"))
    e$n <- e$n+1L; t <- proc.time()[[3]]
    if(e$n<=256L) e$received[[e$n]] <- list(started=as.numeric(Sys.time()),free=as.numeric(x))
    z <- score(ctx,p);e$last_components <- z; e$seconds <- e$seconds+proc.time()[[3]]-t
    model_pid <- if(exists(".ivi_last_model_pid",envir=.GlobalEnv,inherits=FALSE)) get(".ivi_last_model_pid",envir=.GlobalEnv) else Sys.getpid()
    if(z$L<e$best) {
      pred <- predictions(ctx,p,TRUE)
      e$improvements <- e$improvements+1L
      record <- list(parameters=p,components=z,predictions=pred,method=method,evaluation=e$n,
        policy=ctx$policy,provenance=ctx$provenance,evaluation_mode=ctx$evaluation_mode,model_pid=model_pid,created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
        contract_sha256=hash_file(file.path(ctx$study,"manifests/contract.rds")),fixed=fixed,
        parameter_names=names(ctx$endpoint),free_names=free,source_manifest_sha256=hash_file(file.path(ctx$study,"manifests/input_manifest.tsv")))
      path <- file.path(pidout,sprintf("improvement_%08d.rds",e$improvements))
      if(file.exists(path)) stop("Immutable improvement already exists")
      atomic_rds(record,path)
      e$best_record <- record
      write_kv(list(path=basename(path),sha256=hash_file(path),objective=z$L),file.path(pidout,"BEST.tsv"))
      e$best <- z$L
    }
    if(e$n%%16L==0L || e$n==1L) write_kv(list(evaluations=e$n,best=e$best,objective_seconds=e$seconds,improvements=e$improvements),file.path(pidout,"evaluations.tsv"))
    z$L
  }
  attr(fn,"free") <- free
  fn
}
incumbents <- function(out,ctx) {
  paths <- list.files(out,pattern="^improvement_.*\\.rds$",recursive=TRUE,full.names=TRUE)
  if(!length(paths)) stop("No complete feasible incumbent")
  records <- lapply(paths,readRDS)
  stopifnot(all(vapply(records,function(z) identical(z$source_manifest_sha256,hash_file(file.path(ctx$study,"manifests/input_manifest.tsv"))) && identical(z$contract_sha256,hash_file(file.path(ctx$study,"manifests/contract.rds"))),logical(1))))
  records[[which.min(vapply(records,function(z) z$components$L,numeric(1)))]]
}
run_global <- function(study,entry,task,spec,fixture=FALSE) {
  out <- file.path(study,if(fixture) "tests" else "fits",task); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  policy <- spec$policy;fixed <- spec$fixed
  ctx <- formal_context(study,spec)
  atomic_rds(spec,file.path(out,"actual_spec.rds"))
  stopifnot(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))==8L)
  free <- setdiff(names(ctx$endpoint),names(fixed)); seed <- spec$seed
  set.seed(seed); if(!is.null(spec$rng_state)) .Random.seed <- spec$rng_state; initial_rng <- .Random.seed
  start <- resolve_start(ctx,spec);atomic_rds(start,file.path(out,"actual_start.rds"))
  ctx$provenance$start_sha256 <- hash_file(file.path(out,"actual_start.rds"))
  center <- start$parameters;if(length(fixed)) center[names(fixed)] <- fixed
  plan <- ctx$backend$.build_de_initialpop(if(fixture) 40L else 256L,ctx$lower[free],ctx$upper[free],if(spec$init_mode=="uniform") NULL else center[free],spec$init_mode,.35,if(is.null(spec$sigma_frac)) .08 else spec$sigma_frac)
  if(!is.null(spec$initial_population)) plan<-list(pop=spec$initial_population,mode="population_warm_restart")
  pop <- plan$pop
  if(is.null(pop)) stop("Native population builder contract changed")
  stopifnot(ncol(pop)==length(free),nrow(pop)==if(fixture) 40L else 256L,all(is.finite(pop)))
  atomic_rds(list(free=pop,fixed=fixed,free_names=free,plan=plan,initial_rng=initial_rng,seed=seed),file.path(out,"actual_initial_population.rds"))
  deadline <- getOption("ivi.job.deadline");stopifnot(is.finite(deadline))
  .ivi_evaluator <<- new_evaluator(ctx,out,fixed,"DEoptim",deadline)
  ivi_objective <<- function(x) .ivi_evaluator(x)
  cl <- parallel::makePSOCKcluster(8L,outfile=file.path(out,"workers.log"))
  on.exit(if(!is.null(cl)) try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  initialized <- parallel::clusterCall(cl,function(study_arg,entry_arg,task_arg,spec_arg,provenance_arg,fixture_arg,deadline_arg) {
    study <<- study_arg
    source(file.path(entry_arg,"adapter.R"),local=.GlobalEnv)
    source(file.path(entry_arg,"optimize.R"),local=.GlobalEnv)
    source(file.path(entry_arg,"helpers.R"),local=.GlobalEnv)
    source(file.path(entry_arg,"control.R"),local=.GlobalEnv)
    ctx <- formal_context(study_arg,spec_arg);ctx$provenance <- provenance_arg
    fixed_arg <- spec_arg$fixed
    out <- file.path(study_arg,if(fixture_arg) "tests" else "fits",task_arg)
    .ivi_evaluator <<- new_evaluator(ctx,out,fixed_arg,"DEoptim",deadline_arg)
    ivi_objective <<- function(x) .ivi_evaluator(x)
    list(pid=Sys.getpid(),R_version=as.character(getRversion()),source_manifest=hash_file(file.path(study_arg,"manifests/input_manifest.tsv")),dll=o2simps_cpp_dll_info()$path,libpaths=.libPaths(),probe=score(ctx,ctx$endpoint)$L)
  },study,entry,task,spec,ctx$provenance,fixture,deadline)
  stopifnot(length(initialized)==8L,length(unique(vapply(initialized,`[[`,integer(1),"pid")))==8L,
    max(abs(vapply(initialized,`[[`,numeric(1),"probe")-score(ctx,ctx$endpoint)$L))<1e-8)
  atomic_rds(initialized,file.path(out,"worker_initialization.rds"))
  parallel::clusterExport(cl,"ivi_objective",envir=.GlobalEnv)
  t <- proc.time()[[3]]; limited <- FALSE
  result <- tryCatch(DEoptim::DEoptim(ivi_objective,lower=ctx$lower[free],upper=ctx$upper[free],
    control=list(initialpop=pop,NP=nrow(pop),itermax=spec$block_generations,strategy=2L,reltol=0,steptol=spec$block_generations+1L,trace=TRUE,cluster=cl)),
    error=function(e) {
      if(grepl("IVI_SOFT_TIME_LIMIT",conditionMessage(e),fixed=TRUE) || grepl("IVI_CALL_LIMIT",conditionMessage(e),fixed=TRUE)) {limited <<- TRUE; return(NULL)}
      stop(e)
    })
  worker_counts <- parallel::clusterCall(cl,function() c(as.list(.ivi_state),list(runtime_meter=meter_values())))
  master_counts <- c(as.list(.ivi_state),list(runtime_meter=meter_values()))
  actual_pop <- parallel::clusterCall(cl,function() list(free_names=attr(.ivi_evaluator,"free"),pid=Sys.getpid()))
  parallel::stopCluster(cl); cl <- NULL
  active <- function(z) {
    f <- file.path("/proc",as.character(z$pid),"stat")
    line <- tryCatch(readLines(f,warn=FALSE),error=function(e) character())
    length(line)>0L && !sub("^.*\\) ([A-Z]).*$","\\1",line[1])%in%c("Z","X")
  }
  for(i in 1:100) { if(!any(vapply(initialized,active,logical(1)))) break; Sys.sleep(.1) }
  workers_gone <- !any(vapply(initialized,active,logical(1)))
  stopifnot(workers_gone,all(vapply(actual_pop,function(z) identical(z$free_names,free),logical(1))))
  write_kv(list(pass=TRUE,workers_started=8L,workers_stopped=TRUE,free_coordinates=length(free),fixed_coordinates=length(fixed)),file.path(out,"worker_shutdown.tsv"))
  received <- unlist(lapply(c(worker_counts,list(master_counts)),`[[`,"received"),recursive=FALSE)
  observed <- min(length(received),nrow(pop))
  first <- received[head(order(vapply(received,`[[`,numeric(1),"started")),observed)]
  rowhash <- function(x) digest::digest(as.numeric(x),algo="sha256")
  requested_hashes <- vapply(seq_len(nrow(pop)),function(i) rowhash(pop[i,]),character(1))
  received_hashes <- vapply(first,function(z) rowhash(z$free),character(1))
  requested_counts<-table(requested_hashes);received_counts<-table(received_hashes)
  partial_subset<-all(names(received_counts)%in%names(requested_counts)) && all(received_counts<=requested_counts[names(received_counts)])
  complete_population<-observed==nrow(pop)
  stopifnot(partial_subset,is.null(result)||complete_population)
  if(complete_population)stopifnot(identical(sort(requested_hashes),sort(received_hashes)))
  atomic_rds(first,file.path(out,"actual_native_initial_evaluations.rds"))
  write_kv(list(pass=TRUE,requested_population=nrow(pop),native_initial_evaluations=length(first),initial_population_complete=complete_population,
    exact_population_multiset=complete_population,partial_population_subset=partial_subset),file.path(out,"initial_population_validation.tsv"))
  paths<-list.files(out,pattern='^improvement_.*[.]rds$',recursive=TRUE,full.names=TRUE)
  # Projecting a free/neighbor start to another density invalidates its old
  # record as an incumbent at this anchor. Retain only an identical own point.
  own_record<-if(!is.null(start$record)&&identical(start$record$parameters,start$parameters))start$record else NULL
  best<-if(length(paths))incumbents(out,ctx) else if(!is.null(own_record))own_record else stop('No complete feasible incumbent')
  if(!is.null(own_record)&&own_record$components$L<best$components$L)best<-own_record
  checked <- score(ctx,best$parameters)
  stopifnot(abs(checked$L-best$components$L)<1e-8,all(best$parameters[names(fixed)]==fixed))
  atomic_rds(result,file.path(out,"native_result.rds"))
  atomic_rds(list(record=best,checked=checked,fixed=fixed,policy=policy,budget_limited=limited,fixture=fixture,provenance=ctx$provenance),file.path(out,"accepted_candidate.rds"))
  atomic_rds(list(master=.Random.seed,workers=worker_counts,master_evaluator=master_counts),file.path(out,"final_rng_and_counters.rds"))
  atomic_rds(lapply(worker_counts,`[[`,"runtime_meter"),file.path(out,"worker_runtime_meters.rds"))
  write_kv(list(pass=TRUE,fixture=fixture,policy=policy,objective=checked$L,wall_seconds=proc.time()[[3]]-t,evaluations=sum(vapply(worker_counts,`[[`,integer(1),"n")),objective_seconds=sum(vapply(worker_counts,`[[`,numeric(1),"seconds")),iterations=if(is.null(result)) NA_integer_ else as.integer(result$optim$iter),target_iterations=spec$block_generations,budget_limited=limited,fixed_log10_density=if(length(fixed)) unname(fixed) else NA_real_,fixed_residual=if(length(fixed)) max(abs(best$parameters[names(fixed)]-fixed)) else 0,workers=8L,native_continuation_supported=FALSE),file.path(out,"result.tsv"))
  finalize_payload(out)
}
