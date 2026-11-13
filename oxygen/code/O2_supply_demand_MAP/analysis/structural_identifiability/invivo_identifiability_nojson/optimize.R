# Tumor optimizer orchestration using native DEoptim and production L-BFGS-B.
# A timing chain is not an unrestricted baseline or a density profile point.
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
    e$n <- e$n+1L; t <- proc.time()[[3]]
    if(e$n<=256L) e$received[[e$n]] <- list(started=as.numeric(Sys.time()),free=as.numeric(x))
    z <- score(ctx,p); e$seconds <- e$seconds+proc.time()[[3]]-t
    model_pid <- if(exists(".ivi_last_model_pid",envir=.GlobalEnv,inherits=FALSE)) get(".ivi_last_model_pid",envir=.GlobalEnv) else Sys.getpid()
    if(z$L<e$best) {
      pred <- predictions(ctx,p,TRUE)
      e$improvements <- e$improvements+1L
      record <- list(parameters=p,components=z,predictions=pred,method=method,evaluation=e$n,
        policy=ctx$policy,evaluation_mode=ctx$evaluation_mode,model_pid=model_pid,created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
        contract_sha256=hash_file(file.path(ctx$study,"manifests/contract.rds")),fixed=fixed,
        parameter_names=names(ctx$endpoint),free_names=free,source_manifest_sha256=hash_file(file.path(ctx$study,"manifests/input_manifest.tsv")))
      path <- file.path(pidout,sprintf("improvement_%08d.rds",e$improvements))
      if(file.exists(path)) stop("Immutable improvement already exists")
      atomic_rds(record,path)
      write_kv(list(path=basename(path),sha256=hash_file(path),objective=z$L),file.path(pidout,"BEST.tsv"))
      e$best <- z$L
    }
    if(e$n%%16L==0L || e$n==1L) write_kv(list(evaluations=e$n,best=e$best,objective_seconds=e$seconds,improvements=e$improvements),file.path(pidout,"evaluations.tsv"))
    z$L
  }
  attr(fn,"free") <- free
  fn
}
incumbents <- function(out) {
  paths <- list.files(file.path(out,"evaluators"),pattern="^improvement_.*\\.rds$",recursive=TRUE,full.names=TRUE)
  if(!length(paths)) stop("No complete feasible incumbent")
  records <- lapply(paths,readRDS)
  stopifnot(all(vapply(records,function(z) identical(z$source_manifest_sha256,hash_file(file.path(study,"manifests/input_manifest.tsv"))) && identical(z$contract_sha256,hash_file(file.path(study,"manifests/contract.rds"))),logical(1))))
  records[[which.min(vapply(records,function(z) z$components$L,numeric(1)))]]
}
run_global <- function(study,entry,task,policy,fixed,fixture=FALSE) {
  out <- file.path(study,if(fixture) "tests" else "fits",task); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  ctx <- load_contract(study,TRUE,policy)
  stopifnot(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))==8L)
  free <- setdiff(names(ctx$endpoint),names(fixed)); seed <- 17291L
  set.seed(seed); initial_rng <- .Random.seed
  center <- ctx$endpoint; if(length(fixed)) center[names(fixed)] <- fixed
  plan <- ctx$backend$.build_de_initialpop(if(fixture) 40L else 256L,ctx$lower[free],ctx$upper[free],center[free],"hybrid",.35,.08)
  pop <- plan$pop
  if(is.null(pop)) stop("Native population builder contract changed")
  stopifnot(ncol(pop)==length(free),nrow(pop)==if(fixture) 40L else 256L,all(is.finite(pop)))
  atomic_rds(list(free=pop,fixed=fixed,free_names=free,plan=plan,initial_rng=initial_rng,seed=seed),file.path(out,"actual_initial_population.rds"))
  role <- if(fixture) "validate_parallel" else "global"
  wall <- read_tsv(file.path(study,"manifests/budget.tsv"))
  deadline <- min(getOption("ivi.job.deadline",Inf),as.numeric(Sys.time())+wall$wall_seconds[wall$role==role]-300)
  .ivi_evaluator <<- new_evaluator(ctx,out,fixed,"DEoptim",deadline)
  ivi_objective <<- function(x) .ivi_evaluator(x)
  cl <- parallel::makePSOCKcluster(8L,outfile=file.path(out,"workers.log"))
  on.exit(try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
  initialized <- parallel::clusterCall(cl,function(study_arg,entry_arg,task_arg,policy_arg,fixed_arg,fixture_arg,deadline_arg) {
    study <<- study_arg
    source(file.path(entry_arg,"adapter.R"),local=.GlobalEnv)
    source(file.path(entry_arg,"optimize.R"),local=.GlobalEnv)
    ctx <- load_contract(study_arg,TRUE,policy_arg)
    out <- file.path(study_arg,if(fixture_arg) "tests" else "fits",task_arg)
    .ivi_evaluator <<- new_evaluator(ctx,out,fixed_arg,"DEoptim",deadline_arg)
    ivi_objective <<- function(x) .ivi_evaluator(x)
    list(pid=Sys.getpid(),R_version=as.character(getRversion()),source_manifest=hash_file(file.path(study_arg,"manifests/input_manifest.tsv")),dll=o2simps_cpp_dll_info()$path,libpaths=.libPaths(),probe=score(ctx,ctx$endpoint)$L)
  },study,entry,task,policy,fixed,fixture,deadline)
  stopifnot(length(initialized)==8L,length(unique(vapply(initialized,`[[`,integer(1),"pid")))==8L,
    max(abs(vapply(initialized,`[[`,numeric(1),"probe")-score(ctx,ctx$endpoint)$L))<1e-8)
  atomic_rds(initialized,file.path(out,"worker_initialization.rds"))
  parallel::clusterExport(cl,"ivi_objective",envir=.GlobalEnv)
  t <- proc.time()[[3]]; limited <- FALSE
  result <- tryCatch(DEoptim::DEoptim(ivi_objective,lower=ctx$lower[free],upper=ctx$upper[free],
    control=list(initialpop=pop,NP=nrow(pop),itermax=if(fixture) 3L else 1000L,strategy=2L,reltol=1e-4,steptol=25L,trace=TRUE,cluster=cl)),
    error=function(e) {
      if(grepl("IVI_SOFT_TIME_LIMIT",conditionMessage(e),fixed=TRUE)) {limited <<- TRUE; return(NULL)}
      stop(e)
    })
  worker_counts <- parallel::clusterCall(cl,function() as.list(.ivi_state))
  master_counts <- as.list(.ivi_state)
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
  first <- received[order(vapply(received,`[[`,numeric(1),"started"))[seq_len(nrow(pop))]]
  rowhash <- function(x) digest::digest(as.numeric(x),algo="sha256")
  requested_hashes <- vapply(seq_len(nrow(pop)),function(i) rowhash(pop[i,]),character(1))
  received_hashes <- vapply(first,function(z) rowhash(z$free),character(1))
  stopifnot(identical(sort(requested_hashes),sort(received_hashes)))
  atomic_rds(first,file.path(out,"actual_native_initial_evaluations.rds"))
  write_kv(list(pass=TRUE,requested_population=nrow(pop),native_initial_evaluations=length(first),exact_population_multiset=TRUE),file.path(out,"initial_population_validation.tsv"))
  best <- incumbents(out); checked <- score(ctx,best$parameters)
  stopifnot(abs(checked$L-best$components$L)<1e-8,all(best$parameters[names(fixed)]==fixed))
  atomic_rds(result,file.path(out,"native_result.rds"))
  atomic_rds(list(record=best,checked=checked,fixed=fixed,policy=policy,budget_limited=limited,fixture=fixture),file.path(out,"accepted_candidate.rds"))
  atomic_rds(list(master=.Random.seed,workers=worker_counts,master_evaluator=master_counts),file.path(out,"final_rng_and_counters.rds"))
  write_kv(list(pass=TRUE,fixture=fixture,policy=policy,objective=checked$L,wall_seconds=proc.time()[[3]]-t,evaluations=sum(vapply(worker_counts,`[[`,integer(1),"n")),objective_seconds=sum(vapply(worker_counts,`[[`,numeric(1),"seconds")),iterations=if(is.null(result)) NA_integer_ else as.integer(result$optim$iter),target_iterations=if(fixture) 3L else 1000L,budget_limited=limited,fixed_log10_density=if(length(fixed)) unname(fixed) else NA_real_,fixed_residual=if(length(fixed)) max(abs(best$parameters[names(fixed)]-fixed)) else 0,workers=8L,native_continuation_supported=FALSE),file.path(out,"result.tsv"))
  finalize_payload(out)
}
run_local <- function(study,task,policy,fixture=FALSE) {
  out <- file.path(study,if(fixture) "tests" else "fits",task); dir.create(out,recursive=TRUE,showWarnings=FALSE)
  stopifnot(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))==1L)
  ctx <- load_contract(study,TRUE,policy)
  fixed <- setNames(5.198970004336019,"log10_rho_2N")
  p <- if(fixture) ctx$endpoint else readRDS(file.path(study,"fits",sub("local_","global_",task),"accepted_candidate.rds"))$record$parameters
  p[names(fixed)] <- fixed; free <- setdiff(names(p),names(fixed))
  atomic_rds(p,file.path(out,"actual_local_start.rds"))
  role <- if(fixture) "validate_serial" else "local"
  wall <- read_tsv(file.path(study,"manifests/budget.tsv"))
  .ivi_evaluator <<- new_evaluator(ctx,out,fixed,"L-BFGS-B",min(getOption("ivi.job.deadline",Inf),as.numeric(Sys.time())+wall$wall_seconds[wall$role==role]-300))
  t <- proc.time()[[3]]; start <- .ivi_evaluator(p[free]); limited <- FALSE
  res <- tryCatch(optim(p[free],.ivi_evaluator,method="L-BFGS-B",lower=ctx$lower[free],upper=ctx$upper[free],control=list(maxit=if(fixture) 2L else 15000L)),error=function(e) {
    if(grepl("IVI_SOFT_TIME_LIMIT",conditionMessage(e),fixed=TRUE)) {limited <<- TRUE; return(NULL)}
    stop(e)
  })
  best <- incumbents(out); z <- score(ctx,best$parameters)
  pred <- predictions(ctx,best$parameters,TRUE)
  ind <- independent_components(ctx,best$parameters,pred)
  stopifnot(abs(z$L-best$components$L)<1e-8,max(abs(ind-unlist(z[names(ind)])))<1e-8,z$L<=start+1e-8,best$parameters[[names(fixed)]]==unname(fixed))
  atomic_rds(res,file.path(out,"native_result.rds"))
  atomic_rds(list(record=best,checked=z,predictions_full=pred,fixed=fixed,policy=policy,budget_limited=limited,fixture=fixture),file.path(out,"accepted_candidate.rds"))
  for(k in names(pred$tables)) write_tsv(pred$tables[[k]],file.path(out,paste0(k,"_predictions.tsv")))
  write_kv(list(pass=TRUE,fixture=fixture,policy=policy,start_objective=start,objective=z$L,wall_seconds=proc.time()[[3]]-t,evaluations=.ivi_state$n,objective_seconds=.ivi_state$seconds,improvements=.ivi_state$improvements,native_convergence=if(is.null(res)) NA_integer_ else as.integer(res$convergence),budget_limited=limited,fixed_log10_density=unname(fixed),fixed_residual=0,workers=1L),file.path(out,"result.tsv"))
  finalize_payload(out)
}
