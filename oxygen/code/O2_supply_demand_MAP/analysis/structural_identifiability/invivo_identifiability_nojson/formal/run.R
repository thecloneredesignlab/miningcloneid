args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==3L)
study <- normalizePath(args[1],mustWork=TRUE);mode <- args[2];task <- args[3]
entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"));source(file.path(entry,"helpers.R"));source(file.path(entry,"optimize.R"))
specs <- readRDS(file.path(study,"manifests/task_specs.rds"));spec <- specs[[task]]
stopifnot(!is.null(spec),identical(mode,spec$mode),identical(task,spec$task))
dispatch <- file.path(study,"dispatch",task)
m <- read_tsv(file.path(dispatch,"manifest.tsv"));stopifnot(all(vapply(file.path(dispatch,m$path),hash_file,character(1))==m$sha256))
context <- read_kv(file.path(dispatch,"context.tsv"))
stopifnot(identical(context$task,task),identical(context$mode,mode),identical(context$policy,spec$policy),
 identical(context$input_manifest_sha256,hash_file(file.path(study,"manifests/input_manifest.tsv"))),
 identical(context$task_specs_sha256,hash_file(file.path(study,"manifests/task_specs.rds"))),
 identical(context$release_manifest_sha256,hash_file(file.path(entry,"manifest.tsv"))))
generation <- readLines(file.path(study,"state/CURRENT"),warn=FALSE)
ledger <- read.delim(file.path(study,"state",generation,"jobs.tsv"),colClasses="character")
receipt <- ledger[ledger$task==task,,drop=FALSE]
stopifnot(nrow(receipt)==1L,receipt$dispatch_sha256==hash_file(file.path(dispatch,"manifest.tsv")),
 !nzchar(receipt$job_id) || receipt$job_id==Sys.getenv("SLURM_JOB_ID"))
budget <- read_tsv(file.path(study,"manifests/budget.tsv"));b <- budget[budget$role==spec$role,,drop=FALSE]
stopifnot(nrow(b)==1L,as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))==b$cores)
margin <- if(b$wall_seconds<=900L) 120L else if(b$wall_seconds<=1800L) 180L else 300L
deadline <- as.numeric(Sys.getenv("IVI_ALLOCATION_START"))+b$wall_seconds-margin
stopifnot(is.finite(deadline),as.numeric(Sys.time())<deadline)
options(ivi.job.deadline=deadline,ivi.root.task=task,ivi.provenance=c(context,list(dispatch_sha256=hash_file(file.path(dispatch,"manifest.tsv")))))
out <- file.path(study,if(mode%in%c("global","local")) "fits" else "tests",task)
if(file.exists(out)) stop("Task output already exists; refusing overwrite")
dir.create(out,recursive=TRUE)
if(mode=="engineering_parallel") Sys.setenv(MININGCLONEID_RCPP_REBUILD="TRUE")
ctx <- formal_context(study,spec);Sys.setenv(MININGCLONEID_RCPP_REBUILD="FALSE")
baseline <- validate_point(ctx,ctx$endpoint)
old <- readRDS(file.path(study,"reference/production_replay.rds"))
expected <- if(spec$policy=="production") old$components$L else old$components$L_data
stopifnot(abs(baseline$components$L-expected)<1e-8)
write_kv(list(pass=TRUE,objective=baseline$components$L,difference=baseline$components$L-expected),file.path(out,"baseline_recheck.tsv"))
atomic_rds(spec,file.path(out,"dispatch_spec.rds"))
if(mode=="engineering_parallel") {
 source(file.path(entry,"test_formal.R"))
 k <- 0L
 for(policy in c("production","data_only")) for(init in c("hybrid","uniform")) for(fixed in c(FALSE,TRUE)) {
  k <- k+1L;s <- spec;s$policy <- policy;s$init_mode <- init;s$family <- init;s$seed <- spec$seed+k
  s$fixed <- if(fixed) setNames(ctx$rho_grid[4],"log10_rho_2N") else NULL
  s$task <- paste0(task,"_",policy,"_",init,"_",if(fixed) "fixed" else "free")
  run_global(study,entry,s$task,s,TRUE)
 }
 write_kv(list(pass=TRUE,eight_worker_cases=8L,strict_control_reltol=spec$reltol,strict_control_patience=spec$steptol),file.path(out,"result.tsv"));finalize_payload(out)
} else if(mode=="engineering_serial") {
 for(fixed in c(FALSE,TRUE)) {
  s <- spec;s$fixed <- if(fixed) setNames(ctx$rho_grid[4],"log10_rho_2N") else NULL
  s$task <- paste0(task,"_",if(fixed) "fixed" else "free")
  run_local(study,s$task,s,TRUE)
 }
 # The direction loop is exercised on the real objective in each reference task
 # and its box/step bookkeeping has a separate analytic metadata oracle.
 write_kv(list(pass=TRUE,serial_cases=2L,policy=spec$policy),file.path(out,"result.tsv"));finalize_payload(out)
} else if(mode=="reference") {
 run_reference(study,task,spec)
} else if(mode=="global") {
 run_global(study,entry,task,spec,FALSE)
} else if(mode=="local") {
 run_local(study,task,spec,FALSE)
} else if(mode%in%c("fresh","directions","profile_fresh")) {
 run_acceptance(study,task,spec,mode!="fresh")
} else stop("Unsupported formal task mode")
