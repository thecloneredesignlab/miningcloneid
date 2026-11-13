args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==3L)
study <- normalizePath(args[1],mustWork=TRUE); mode <- args[2]; task <- args[3]
entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"))
role <- switch(mode,replay="replay",validate_parallel="validate_parallel",validate_serial="validate_serial",global="global",local="local",check="check",stop("Unsupported allocated mode"))
budget <- read_tsv(file.path(study,"manifests/budget.tsv"))
options(ivi.job.deadline=as.numeric(Sys.time())+budget$wall_seconds[budget$role==role]-300)
out <- file.path(study,if(mode %in% c("replay","validate_parallel","validate_serial","check")) "tests" else "fits",task)
if(file.exists(file.path(out,"COMPLETE"))) stop("Task already complete; refusing overwrite")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(mode=="validate_parallel") Sys.setenv(MININGCLONEID_RCPP_REBUILD="TRUE")
ctx <- load_contract(study,TRUE,"production")
Sys.setenv(MININGCLONEID_RCPP_REBUILD="FALSE")
baseline_check <- score(ctx,ctx$endpoint)
stopifnot(abs(baseline_check$L-ctx$expected_score)<1e-8)
write_kv(list(pass=TRUE,objective=baseline_check$L,difference=baseline_check$L-ctx$expected_score),file.path(out,"baseline_recheck.tsv"))
if(mode=="replay") {
  t <- proc.time()[[3]]; z <- score(ctx,ctx$endpoint)
  stopifnot(abs(z$L-ctx$expected_score)<1e-8)
  pred <- predictions(ctx,ctx$endpoint,FALSE)
  ind <- independent_components(ctx,ctx$endpoint,pred)
  stopifnot(max(abs(ind-unlist(z[names(ind)])))<1e-8)
  for(k in c("burden","ploidy")) {
    old <- read_tsv(file.path(study,"source/seed25",paste0(if(k=="ploidy") "terminal_ploidy" else "burden","_fit.tsv")))
    new <- pred$tables[[k]]
    stopifnot(identical(as.character(old$harvest),as.character(new$harvest)),nrow(old)==nrow(new))
    for(n in intersect(names(old),names(new))) if(is.numeric(old[[n]]) && all(is.finite(old[[n]]))) {
      stopifnot(all(is.finite(new[[n]])),all(abs(old[[n]]-new[[n]])<=1e-8+1e-8*abs(old[[n]])))
    }
  }
  prior <- ctx$backend$compute_soft_prior_penalty(ctx$endpoint,ctx$cfg)
  write_tsv(data.frame(parameter=names(prior$terms),raw=unname(prior$terms),weight=ctx$cfg$lambda_prior,weighted=unname(prior$terms)*ctx$cfg$lambda_prior),file.path(out,"prior_terms.tsv"))
  for(k in names(pred$tables)) write_tsv(pred$tables[[k]],file.path(out,paste0(k,"_predictions.tsv")))
  atomic_rds(list(parameters=ctx$endpoint,components=z,predictions=pred),file.path(out,"replay.rds"))
  write_kv(c(list(pass=TRUE,expected=ctx$expected_score,difference=z$L-ctx$expected_score,wall_seconds=proc.time()[[3]]-t),component_values(z)),file.path(out,"result.tsv"))
  write_tsv(data.frame(component=names(ind),independent=ind,native=unlist(z[names(ind)]),difference=ind-unlist(z[names(ind)])),file.path(out,"independent_scores.tsv"))
  finalize_payload(out)
} else if(mode=="check") {
  source_task <- sub("^check_","",task)
  policy <- if(grepl("data_only",source_task,fixed=TRUE)) "data_only" else "production"
  ctx <- load_contract(study,TRUE,policy)
  candidate <- readRDS(file.path(study,"fits",source_task,"accepted_candidate.rds"))
  stopifnot(!isTRUE(candidate$fixture),identical(candidate$policy,policy),
    identical(candidate$record$source_manifest_sha256,hash_file(file.path(study,"manifests/input_manifest.tsv"))),
    identical(candidate$record$contract_sha256,hash_file(file.path(study,"manifests/contract.rds"))))
  p <- candidate$record$parameters; z <- score(ctx,p)
  stopifnot(all(p[names(candidate$fixed)]==candidate$fixed),abs(z$L-candidate$checked$L)<1e-8,abs(z$L-candidate$record$components$L)<1e-8)
  pred <- predictions(ctx,p,TRUE); ind <- independent_components(ctx,p,pred)
  stopifnot(max(abs(ind-unlist(z[names(ind)])))<1e-8,
    isTRUE(all.equal(pred$tables,candidate$record$predictions$tables,tolerance=1e-8)))
  atomic_rds(list(parameters=p,components=z,predictions=pred,policy=policy,source_task=source_task),file.path(out,"fresh_acceptance.rds"))
  write_kv(c(list(pass=TRUE,source_task=source_task,policy=policy,difference=z$L-candidate$checked$L,fixed_residual=max(abs(p[names(candidate$fixed)]-candidate$fixed))),component_values(z)),file.path(out,"result.tsv"))
  finalize_payload(out)
} else if(mode=="validate_parallel") {
  source(file.path(entry,"optimize.R"))
  source(file.path(entry,"test_runtime.R"))
  for(policy in c("production","data_only")) for(constrained in c(FALSE,TRUE)) {
    fixed <- if(constrained) setNames(5.198970004336019,"log10_rho_2N") else NULL
    run_global(study,entry,paste0(task,"_",policy,"_",if(constrained) "fixed" else "free"),policy,fixed,TRUE)
  }
  write_kv(list(pass=TRUE,fresh_private_compilation=TRUE,actual_eight_worker_cases=4L),file.path(out,"result.tsv"))
  finalize_payload(out)
} else if(mode=="validate_serial") {
  source(file.path(entry,"optimize.R")); policy <- if(grepl("data_only",task,fixed=TRUE)) "data_only" else "production"
  run_local(study,task,policy,TRUE)
} else if(mode=="global") {
  source(file.path(entry,"optimize.R")); policy <- if(grepl("data_only",task,fixed=TRUE)) "data_only" else "production"
  run_global(study,entry,task,policy,setNames(5.198970004336019,"log10_rho_2N"),FALSE)
} else if(mode=="local") {
  source(file.path(entry,"optimize.R")); policy <- if(grepl("data_only",task,fixed=TRUE)) "data_only" else "production"
  run_local(study,task,policy,FALSE)
} else stop("Unsupported task mode")
