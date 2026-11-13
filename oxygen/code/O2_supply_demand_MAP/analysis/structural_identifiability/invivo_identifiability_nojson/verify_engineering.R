# Read-only fixture acceptance, followed by a release-specific engineering gate.
# This script performs no model evaluations.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1],mustWork=TRUE)
entry <- file.path(study,"releases/fresh_pipeline01")
source(file.path(entry,"adapter.R"))
contract <- readRDS(file.path(study,"manifests/contract.rds"))
jobs <- read_tsv(file.path(study,"tables/jobs.tsv"))
tasks <- c("fresh_parallel","fresh_serial_production","fresh_serial_data_only")
for(task in tasks) {
 j <- jobs[jobs$task==task,,drop=FALSE]
 stopifnot(nrow(j)==1L,j$state=="COMPLETED",j$exit_code=="0:0",j$release=="fresh_pipeline01")
 stopifnot(isTRUE(read_kv(file.path(study,"tests",task,"result.tsv"))$pass))
}
stopifnot(isTRUE(read_kv(file.path(study,"tests/fresh_parallel_unit/result.tsv"))$pass))
cases <- c(as.vector(outer(c("production","data_only"),c("free","fixed"),
 function(p,f) paste0("fresh_parallel_",p,"_",f))),tasks[-1])
records <- 0L
max_difference <- 0
for(task in cases) {
 out <- file.path(study,"tests",task)
 result <- read_kv(file.path(out,"result.tsv"))
 stopifnot(isTRUE(result$pass),isTRUE(result$fixture),!isTRUE(result$budget_limited))
 accepted <- readRDS(file.path(out,"accepted_candidate.rds"))
 stopifnot(isTRUE(accepted$fixture),accepted$record$evaluation_mode=="fresh_fork")
 difference <- abs(accepted$record$components$L-accepted$checked$L)
 stopifnot(difference<1e-8)
 max_difference <- max(max_difference,difference)
 if(startsWith(task,"fresh_parallel_")) {
  pop <- readRDS(file.path(out,"actual_initial_population.rds"))
  received <- readRDS(file.path(out,"actual_native_initial_evaluations.rds"))
  nfree <- if(endsWith(task,"_fixed")) 19L else 20L
  stopifnot(nrow(pop$free)==40L,ncol(pop$free)==nfree,length(received)==40L,
    isTRUE(read_kv(file.path(out,"initial_population_validation.tsv"))$exact_population_multiset),
    isTRUE(read_kv(file.path(out,"worker_shutdown.tsv"))$workers_stopped))
  workers <- readRDS(file.path(out,"worker_initialization.rds"))
  stopifnot(length(workers)==8L,length(unique(vapply(workers,`[[`,integer(1),"pid")))==8L)
 }
 paths <- list.files(file.path(out,"evaluators"),pattern="^improvement_.*[.]rds$",recursive=TRUE,full.names=TRUE)
 for(path in paths) {
  r <- readRDS(path)
  stopifnot(r$evaluation_mode=="fresh_fork",
   r$model_pid!=as.integer(sub("^pid_","",basename(dirname(path)))),
   identical(r$source_manifest_sha256,hash_file(file.path(study,"manifests/input_manifest.tsv"))),
   identical(r$contract_sha256,hash_file(file.path(study,"manifests/contract.rds"))),
   all(r$parameters[names(r$fixed)]==r$fixed))
  records <- records+1L
 }
}
gate <- file.path(study,"tables/engineering_gate.tsv")
history <- file.path(study,"tables/gate_history/engineering_gate_before_amendment02.tsv")
stopifnot(!file.exists(history));stopifnot(file.copy(gate,history))
write_kv(list(pass=TRUE,release="fresh_pipeline01",release_manifest_sha256=hash_file(file.path(entry,"manifest.tsv")),
 evaluation_mode="fresh_fork",optimizer_fixture_cases=6L,eight_worker_cases=4L,
 serial_cases=2L,fixture_improvement_records=records,score_reevaluation_max_difference=max_difference,
 model_evaluations_in_gate_script=0L,source_model_unchanged=TRUE,
 scope="Engineering acceptance only; not unrestricted baseline recovery or statistical calibration"),gate)
cat("Fresh-fork engineering gate passed;",records,"fixture improvement records; max score gap",max_difference,"\n")
