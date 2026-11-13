# Audit each collected completed fit once per immutable payload; no model calls.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=1L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"));source(file.path(entry,"helpers.R"))
contract <- readRDS(file.path(study,"manifests/contract.rds"));backend <- load_backend(study,FALSE)
tasks <- if(length(args)>1L) args[-1] else basename(list.dirs(file.path(study,"fits"),recursive=FALSE))
out <- file.path(study,"tables/record_audits");dir.create(out,recursive=TRUE,showWarnings=FALSE)
near <- function(a,b) all(is.finite(a)) && all(is.finite(b)) && all(abs(a-b)<=1e-8+1e-8*abs(b))
for(task in tasks) {
 root <- file.path(study,"fits",task)
 if(!dir.exists(root)) root <- file.path(study,"tests",task)
 verify_payload_R(root)
 summary_path <- file.path(out,paste0(task,"_summary.tsv"));payload <- hash_file(file.path(root,"payload_manifest.tsv"))
 if(file.exists(summary_path)) {
  prior <- read_kv(summary_path)
  stopifnot(isTRUE(prior$pass),identical(prior$payload_sha256,payload));next
 }
 rows <- list();maxdiff <- 0;all_scores <- numeric()
 for(directory in list.dirs(file.path(root,"evaluators"),recursive=FALSE,full.names=TRUE)) {
  files <- sort(list.files(directory,pattern="^improvement_[0-9]+[.]rds$",full.names=TRUE));if(!length(files)) next
  last <- Inf
  for(path in files) {
   r <- readRDS(path);p <- r$parameters;c <- r$components
   stopifnot(identical(names(p),names(contract$endpoint)),all(p>=contract$lower-1e-12),all(p<=contract$upper+1e-12),
    r$evaluation_mode=="fresh_fork",r$model_pid!=as.integer(sub("^pid_","",basename(directory))),
    identical(r$contract_sha256,hash_file(file.path(study,"manifests/contract.rds"))),
    identical(r$source_manifest_sha256,hash_file(file.path(study,"manifests/input_manifest.tsv"))),
    all(p[names(r$fixed)]==r$fixed),c$L<last,
    c$n_burden_obs_total==92L,c$n_ploidy_obs_total==6502L,c$n_necrosis_obs_total==6L)
   ctx <- contract;ctx$backend <- backend
   ctx$cfg$parameter_table <- file.path(study,"source/seed25/parameter_table_input.csv")
   if(r$policy=="data_only") ctx$cfg$use_soft_prior <- FALSE
   independent <- independent_components(ctx,p,r$predictions)
   diff <- max(abs(independent-unlist(c[names(independent)])))
   prior <- backend$compute_soft_prior_penalty(p,ctx$cfg)
   stopifnot(diff<1e-8,abs(prior$L_prior_raw-c$L_prior_raw)<1e-8,
    abs(c$L_data-c$L_b-c$L_p-c$L_n)<1e-8,abs(c$L-c$L_data-c$L_prior)<1e-8,
    abs(c$L_prior-if(ctx$cfg$use_soft_prior) ctx$cfg$lambda_prior*prior$L_prior_raw else 0)<1e-8)
   rho <- 10^p[['log10_rho_2N']]
   for(z in r$predictions$native) {
    n <- z$result
    stopifnot(all(is.finite(n$live_state_obs)),all(n$live_state_obs>=0),
     all(is.finite(n$dead_hypoxia_state_obs)),all(n$dead_hypoxia_state_obs>=0),
     all(is.finite(n$dead_buffer_state_obs)),all(n$dead_buffer_state_obs>=0),
     near(rowSums(n$live_state_obs),n$Ntot_live_obs),near(rowSums(n$dead_hypoxia_state_obs),n$Ntot_dead_hypoxia_obs),
     near(rowSums(n$dead_buffer_state_obs),n$Ntot_dead_buffer_obs),
     near(n$Ntot_live_obs+n$Ntot_dead_hypoxia_obs+n$Ntot_dead_buffer_obs,n$Ntot_total_obs),
     near(n$Vmm3_total_obs,n$Ntot_total_obs/rho),near(n$Vmm3_total_terminal,n$Ntot_total_terminal/rho),
     all(is.finite(n$O2_eff_obs)),all(is.finite(n$O2_target_obs)),near(sum(n$frac_N_live),1))
   }
   rows[[length(rows)+1L]] <- data.frame(task=task,evaluator=basename(directory),record=basename(path),
    method=r$method,evaluation=r$evaluation,objective=c$L,data=c$L_data,prior=c$L_prior,
    independent_component_max_gap=diff,sha256=hash_file(path))
   maxdiff <- max(maxdiff,diff);all_scores <- c(all_scores,c$L);last <- c$L
  }
  pointer <- read_kv(file.path(directory,"BEST.tsv"))
  stopifnot(pointer$path==basename(tail(files,1)),pointer$sha256==hash_file(tail(files,1)),pointer$objective==last)
 }
 accepted <- readRDS(file.path(root,"accepted_candidate.rds"))
 stopifnot(length(all_scores)>0L,abs(min(all_scores)-accepted$checked$L)<1e-8,
  abs(accepted$checked$L-accepted$record$components$L)<1e-8)
 write_tsv(do.call(rbind,rows),file.path(out,paste0(task,"_records.tsv")))
 write_kv(list(pass=TRUE,task=task,payload_sha256=payload,records=length(rows),independent_component_max_gap=maxdiff,
  model_evaluations=0L),summary_path)
 cat(task,":",length(rows),"records audited; maximum component gap",maxdiff,"\n")
}
assert_namespaces()
