# Metadata/archived-state acceptance only: no simulator is loaded or evaluated.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2L)
study <- normalizePath(args[1]);action <- args[2]
entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"));source(file.path(entry,"helpers.R"))
contract <- readRDS(file.path(study,"manifests/contract.rds"))
specs <- readRDS(file.path(study,"manifests/task_specs.rds"))
protocol <- readRDS(file.path(study,"manifests/protocol.rds"))
backend <- load_backend(study,FALSE)
audit_record <- function(record) {
 ctx <- contract;ctx$backend <- backend
 ctx$cfg$parameter_table <- file.path(study,"source/seed25/parameter_table_input.csv")
 if(record$policy=="data_only") ctx$cfg$use_soft_prior <- FALSE
 independent <- independent_components(ctx,record$parameters,record$predictions)
 stopifnot(max(abs(independent-unlist(record$components[names(independent)])))<1e-8)
 prior <- backend$compute_soft_prior_penalty(record$parameters,ctx$cfg)
 stopifnot(abs(prior$L_prior_raw-record$components$L_prior_raw)<1e-8,
  abs(record$components$L-record$components$L_data-record$components$L_prior)<1e-8)
 invisible(TRUE)
}
jobs <- read.delim(file.path(study,"tables/jobs.tsv"),colClasses="character")
release <- "formal01"
gate_values <- function(pass,more=list()) c(list(pass=pass,release=release,
 input_manifest_sha256=hash_file(file.path(study,"manifests/input_manifest.tsv")),
 release_manifest_sha256=hash_file(file.path(study,"releases",release,"manifest.tsv"))),more)
scope <- function(spec) if(spec$mode%in%c("global","local")) "fits" else "tests"
closed_payload <- function(task) {
 j <- jobs[jobs$task==task,,drop=FALSE];stopifnot(nrow(j)==1L,j$state=="COMPLETED",j$exit_code=="0:0",j$release%in%c("formal01","formal02"))
 if(j$release=="formal02") {
  bridge <- read_tsv(file.path(study,"tables/release_bridge_formal02.tsv"))
  stopifnot(all(vapply(file.path(study,bridge$path),hash_file,character(1))==bridge$sha256))
  g <- read_kv(file.path(study,"tables/engineering_gate_formal02.tsv"))
  stopifnot(isTRUE(g$pass),identical(g$parent_engineering_gate_sha256,hash_file(file.path(study,"tables/engineering_gate.tsv"))))
 }
 root <- file.path(study,scope(specs[[task]]),task);verify_payload_R(root);root
}
candidate <- function(task) {
 root <- closed_payload(task)
 a <- readRDS(file.path(root,"accepted_candidate.rds"))
 stopifnot(identical(a$record$source_manifest_sha256,hash_file(file.path(study,"manifests/input_manifest.tsv"))),
  identical(a$record$contract_sha256,hash_file(file.path(study,"manifests/contract.rds"))),
  abs(a$record$components$L-a$checked$L)<1e-8)
 audit_record(a$record)
 list(record=a$record,policy=a$policy,source_task=task,source_scope=scope(specs[[task]]),
  payload_sha256=hash_file(file.path(root,"payload_manifest.tsv")))
}
evidence_manifest <- function(tasks,path) {
 rows <- lapply(tasks,function(task) {
  root <- closed_payload(task);relative <- paste(scope(specs[[task]]),task,"payload_manifest.tsv",sep="/")
  data.frame(path=relative,sha256=hash_file(file.path(root,"payload_manifest.tsv")))
 })
 write_tsv(do.call(rbind,rows),path)
}
if(action=="engineering") {
 parents <- c("engineering_parallel","engineering_serial_production","engineering_serial_data_only","reference_production","reference_data_only")
 for(task in parents) {
  root <- closed_payload(task);stopifnot(isTRUE(read_kv(file.path(root,"result.tsv"))$pass))
 }
 stopifnot(isTRUE(read_kv(file.path(study,"tests/engineering_parallel_unit/result.tsv"))$pass))
 verify_payload_R(file.path(study,"tests/engineering_parallel_unit"))
 for(policy in c("production","data_only")) for(init in c("hybrid","uniform")) for(fixed in c(FALSE,TRUE)) {
  task <- paste0("engineering_parallel_",policy,"_",init,"_",if(fixed) "fixed" else "free")
  root <- file.path(study,"tests",task);verify_payload_R(root)
  pop <- readRDS(file.path(root,"actual_initial_population.rds"));actual <- readRDS(file.path(root,"actual_spec.rds"))
  nfree <- if(fixed) 19L else 20L
  stopifnot(nrow(pop$free)==40L,ncol(pop$free)==nfree,actual$reltol==1e-6,actual$steptol==100L,
   identical(pop$plan$warm_start_used,init=="hybrid"),
   isTRUE(read_kv(file.path(root,"initial_population_validation.tsv"))$exact_population_multiset),
   isTRUE(read_kv(file.path(root,"worker_shutdown.tsv"))$workers_stopped))
  accepted <- readRDS(file.path(root,"accepted_candidate.rds"))
  stopifnot(!isTRUE(accepted$budget_limited),abs(accepted$record$components$L-accepted$checked$L)<1e-8,all(accepted$record$parameters[names(accepted$fixed)]==accepted$fixed))
  audit_record(accepted$record)
 }
 for(policy in c("production","data_only")) for(fixed in c(FALSE,TRUE)) {
  root <- file.path(study,"tests",paste0("engineering_serial_",policy,"_",if(fixed) "fixed" else "free"));verify_payload_R(root)
  a <- readRDS(file.path(root,"accepted_candidate.rds"));nfree <- if(fixed) 19L else 20L
  stopifnot(!isTRUE(a$budget_limited),length(a$record$free_names)==nfree,abs(a$record$components$L-a$checked$L)<1e-8)
  audit_record(a$record)
  root <- file.path(study,"tests",paste0("reference_",policy,"_direction_fixture_",if(fixed) "fixed" else "free"));verify_payload_R(root)
  a <- readRDS(file.path(root,"accepted_candidate.rds"));stopifnot(length(a$record$free_names)==nfree)
  audit_record(a$record)
  probes <- readRDS(file.path(root,"directional_probes.rds"));stopifnot(length(probes)>0L)
  for(p in probes) stopifnot(all(p$parameters>=contract$lower),all(p$parameters<=contract$upper),
   all(p$parameters[names(a$fixed)]==a$fixed),abs(p$components$L-p$objective)<1e-8)
 }
 write_kv(gate_values(TRUE,list(eight_worker_cases=8L,serial_local_cases=4L,reference_points=10L,
  full_direction_fixtures=4L,model_evaluations_in_gate=0L)),file.path(study,"tables/engineering_gate.tsv"))
 cat("Formal engineering accepted: eight DE cases, four serial cases, ten references and four real direction fixtures.\n")
} else if(action=="baseline") {
 summaries <- list();selected <- list();payload_tasks <- character();both <- TRUE;unresolved <- list()
 for(policy in c("production","data_only")) {
  local_tasks <- names(specs)[vapply(specs,function(s) s$mode=="local" && s$policy==policy && s$role%in%c("baseline_local","reserve_baseline_local"),logical(1))]
  endpoints <- list();directions <- list();rows <- list()
  for(task in local_tasks) {
   fresh <- sub("_local$","_fresh",task);direction <- sub("_local$","_directions",task)
   if(!all(c(task,fresh,direction)%in%jobs$task)) next
   j <- jobs[match(c(task,fresh,direction),jobs$task),,drop=FALSE]
   if(!all(j$state=="COMPLETED" & j$exit_code=="0:0")) next
   endpoint <- candidate(task);fr <- closed_payload(fresh);dr <- closed_payload(direction)
   f <- read_kv(file.path(fr,"result.tsv"));d <- read_kv(file.path(dr,"result.tsv"))
   stopifnot(isTRUE(f$pass),abs(f$difference)<1e-8,isTRUE(d$pass))
   endpoints[[task]] <- endpoint;directions[[task]] <- candidate(direction)
   payload_tasks <- c(payload_tasks,task,fresh,direction)
   independent <- specs[[task]]$role=="baseline_local";family <- specs[[task]]$family
   if(specs[[task]]$role=="reserve_baseline_local" && nzchar(specs[[task]]$source_task)) {
    source_root <- closed_payload(specs[[task]]$source_task)
    actual <- readRDS(file.path(source_root,"actual_spec.rds"))
    independent <- actual$selector=="seed" && (actual$init_mode=="uniform" || actual$start_seed%in%c(25L,74L,366L,292L,392L))
    family <- actual$init_mode
   }
   rows[[task]] <- data.frame(policy=policy,task=task,family=family,
    objective=endpoint$record$components$L,probe_best=directions[[task]]$record$components$L,
    native_code=read_kv(file.path(study,"fits",task,"result.tsv"))$native_convergence,
    budget_limited=read_kv(file.path(study,"fits",task,"result.tsv"))$budget_limited,
    independent=independent)
  }
  stopifnot(length(rows)>=4L,all(paste0("b_",policy,"_",c("h25","h392","u1","u2"),"_local")%in%names(rows)))
  table <- do.call(rbind,rows);best_index <- which.min(table$objective);best <- table$objective[best_index]
  table$gap <- table$objective-best
  recovered <- table$independent & table$gap<=protocol$baseline_recovery
  family_pass <- all(c("hybrid","uniform")%in%table$family[recovered])
  direction_best <- min(table$probe_best);direction_pass <- best-direction_best<=protocol$direction_improvement
  known_bound <- min(vapply(c(25L,74L,366L,292L,392L),function(seed)
   readRDS(file.path(study,"tests",paste0("reference_",policy),paste0("seed",seed,"_reference.rds")))$components$L,numeric(1)))
  external <- read_tsv(file.path(study,"reference_external/manifest.tsv"))
  row <- external[external$policy==policy,,drop=FALSE];stopifnot(nrow(row)==1L)
  stopifnot(hash_file(file.path(study,row$path))==row$sha256,
   row$image_sha256=="dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372")
  old_contract <- readRDS(file.path(study,"reference/pilot_contract.rds"))
  stopifnot(identical(old_contract$cfg,contract$cfg),identical(old_contract$scenarios,contract$scenarios))
  old_manifest <- read_tsv(file.path(study,"prior/pilot_input_manifest.tsv"))
  current_manifest <- read_tsv(file.path(study,"manifests/input_manifest.tsv"))
  core <- startsWith(old_manifest$path,"source/numerical/") | startsWith(old_manifest$path,"source/inputs/")
  stopifnot(all(old_manifest$sha256[core]==current_manifest$sha256[match(old_manifest$path[core],current_manifest$path)]))
  external_point <- readRDS(file.path(study,row$path))
  stopifnot(external_point$policy==policy,all(external_point$parameters>=contract$lower),all(external_point$parameters<=contract$upper))
  audit_record(external_point)
  known_bound <- min(known_bound,external_point$components$L)
  known_pass <- best-known_bound<=protocol$baseline_recovery
  pass <- family_pass && direction_pass && known_pass;both <- both && pass
  selected[[policy]] <- endpoints[[table$task[best_index]]]
  atomic_rds(selected[[policy]],file.path(study,"tables",paste0("baseline_",policy,".rds")))
  probe_index <- which.min(table$probe_best)
  atomic_rds(directions[[table$task[probe_index]]],file.path(study,"tables",paste0("baseline_",policy,"_probe_best.rds")))
  summaries[[policy]] <- table
  unresolved[[policy]] <- data.frame(policy=policy,pass=pass,family_recovery=family_pass,directional_pass=direction_pass,
   known_feasible_pass=known_pass,known_feasible_bound=known_bound,known_feasible_gap=best-known_bound,
   best_objective=best,best_probe=direction_best,directional_improvement=best-direction_best,
   recovered_independent_chains=sum(recovered),
   uniform_gap=min(table$gap[table$family=="uniform" & table$independent]),best_task=table$task[best_index])
 }
 write_tsv(do.call(rbind,summaries),file.path(study,"tables/baseline_chain_summary.tsv"))
 write_tsv(do.call(rbind,unresolved),file.path(study,"tables/baseline_policy_gates.tsv"))
 evidence_manifest(unique(payload_tasks),file.path(study,"tables/baseline_gate_payloads.tsv"))
 m <- read_tsv(file.path(study,"tables/baseline_gate_payloads.tsv"))
 for(policy in c("production","data_only")) {
  path <- paste0("tables/baseline_",policy,".rds")
  m <- rbind(m,data.frame(path=path,sha256=hash_file(file.path(study,path))))
 }
 write_tsv(m,file.path(study,"tables/baseline_gate_payloads.tsv"))
 external_paths <- c("reference_external/manifest.tsv",read_tsv(file.path(study,"reference_external/manifest.tsv"))$path)
 m <- rbind(m,data.frame(path=external_paths,sha256=vapply(file.path(study,external_paths),hash_file,character(1))))
 write_tsv(m,file.path(study,"tables/baseline_gate_payloads.tsv"))
 write_kv(gate_values(both,list(data_only_pass=unresolved$data_only$pass,production_pass=unresolved$production$pass,
  model_evaluations_in_gate=0L)),file.path(study,"tables/baseline_gate.tsv"))
 print(do.call(rbind,unresolved),row.names=FALSE)
} else if(action=="profile_select") {
 for(policy in c("production","data_only")) for(i in 1:7) {
  prefix <- paste0("d_",policy,"_",i)
  primary <- paste0(prefix,c("_local","_asc","_desc"))
  if(!all(primary%in%jobs$task)) next
  j <- jobs[match(primary,jobs$task),,drop=FALSE]
  if(!all(j$state=="COMPLETED" & j$exit_code=="0:0")) next
  choices <- primary
  for(task in paste0(prefix,c("_reserve_global","_reserve_local","_fresh","_reserve_fresh"))) {
   j <- jobs[jobs$task==task,,drop=FALSE]
   if(nrow(j)==1L && j$state=="COMPLETED" && j$exit_code=="0:0") choices <- c(choices,task)
  }
  values <- lapply(choices,candidate);k <- which.min(vapply(values,function(a)a$record$components$L,numeric(1)))
  atomic_rds(values[[k]],file.path(study,"tables",paste0("density_",policy,"_",i,"_selected.rds")))
 }
 cat("Selections published only for closed three-branch density anchors.\n")
} else if(action=="profile_review") {
 rows <- list();all_pass <- TRUE;baseline_consistent <- TRUE;counterexamples <- list()
 for(policy in c("production","data_only")) for(i in 1:7) {
  prefix <- paste0("d_",policy,"_",i);primary <- paste0(prefix,c("_local","_asc","_desc"))
  if(!all(primary%in%jobs$task)) next
  j <- jobs[match(primary,jobs$task),,drop=FALSE]
  if(!all(j$state=="COMPLETED" & j$exit_code=="0:0")) next
  fresh_tasks <- paste0(prefix,c("_fresh","_reserve_fresh"))
  fresh_tasks <- fresh_tasks[vapply(fresh_tasks,function(task) {
   j <- jobs[jobs$task==task,,drop=FALSE];nrow(j)==1L && j$state=="COMPLETED" && j$exit_code=="0:0"
  },logical(1))]
  if(!length(fresh_tasks)) next
  choices <- lapply(fresh_tasks,candidate)
  k <- which.min(vapply(choices,function(a)a$record$components$L,numeric(1)));best <- choices[[k]]
  fresh_root <- closed_payload(fresh_tasks[k]);fresh_result <- read_kv(file.path(fresh_root,"result.tsv"))
  branches <- lapply(primary,candidate);scores <- vapply(branches,function(a)a$record$components$L,numeric(1))
  objective <- best$record$components$L
  baseline <- readRDS(file.path(study,"tables",paste0("baseline_",policy,".rds")))
  baseline_value <- baseline$record$components$L;delta <- objective-baseline_value
  selected <- readRDS(file.path(study,"tables",paste0("density_",policy,"_",i,"_selected.rds")))
  selected_verified <- selected$record$components$L>=objective-1e-8
  recovered <- sum(scores-objective<=protocol$profile_recovery)
  direction_pass <- fresh_result$directional_improvement<=protocol$direction_improvement
  consistency <- delta>=-protocol$baseline_recovery
  pass <- selected_verified && recovered>=2L && direction_pass && consistency
  all_pass <- all_pass && pass;baseline_consistent <- baseline_consistent && consistency
  if(!consistency) {
   old <- counterexamples[[policy]]
   if(is.null(old) || objective<old$record$components$L) counterexamples[[policy]] <- best
  }
  rows[[prefix]] <- data.frame(policy=policy,density_index=i,log10_density=contract$rho_grid[i],density=10^contract$rho_grid[i],
   objective=objective,data=best$record$components$L_data,prior=best$record$components$L_prior,
   baseline_objective=baseline_value,raw_delta=delta,recovered_branches=recovered,
   global_branch=scores[1],ascending_branch=scores[2],descending_branch=scores[3],
   directional_improvement=fresh_result$directional_improvement,selected_verified=selected_verified,
   baseline_consistent=consistency,pass=pass,fresh_task=fresh_tasks[k])
  atomic_rds(best,file.path(study,"tables",paste0("profile_",policy,"_",i,"_verified.rds")))
 }
 if(length(rows)) write_tsv(do.call(rbind,rows),file.path(study,"tables/profile_point_review.tsv"))
 complete <- length(rows)==14L
 write_kv(gate_values(complete && all_pass,list(reviewed_anchors=length(rows),complete=complete,
  baseline_consistent=baseline_consistent,model_evaluations_in_gate=0L)),file.path(study,"tables/profile_gate.tsv"))
 if(!baseline_consistent) {
  history <- file.path(study,"tables/gate_history");dir.create(history,showWarnings=FALSE)
  gate <- file.path(study,"tables/baseline_gate.tsv")
  saved <- file.path(history,paste0("baseline_before_profile_counterexample_",hash_file(gate),".tsv"))
  if(!file.exists(saved)) stopifnot(file.copy(gate,saved))
  for(policy in names(counterexamples)) atomic_rds(counterexamples[[policy]],file.path(study,"tables",paste0("baseline_",policy,"_probe_best.rds")))
  write_kv(gate_values(FALSE,list(reason="A constrained candidate beats the matching baseline beyond the declared optimization tolerance")),gate)
 }
 if(length(rows)) print(do.call(rbind,rows)[,c("policy","density_index","raw_delta","recovered_branches","directional_improvement","pass")],row.names=FALSE)
} else stop("Unsupported metadata postprocessing action")
assert_namespaces()
