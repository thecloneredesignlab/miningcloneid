# Validate immutable full-precision records without additional model evaluations.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=2L)
study <- normalizePath(args[1]);pkg <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(pkg,"adapter.R"));e <- load_backend(study,FALSE)
contract <- readRDS(file.path(study,"manifests/contract.rds"))
near <- function(a,b) all(is.finite(a)) && all(is.finite(b)) && all(abs(a-b)<=1e-8+1e-8*abs(b))
score_near <- function(a,b) all(is.finite(a)) && all(is.finite(b)) && all(abs(a-b)<=1e-8)
rows <- list();counter <- 0L
for(task in args[-1]) {
 root <- file.path(study,"fits",task)
 stopifnot(file.exists(file.path(root,"COMPLETE")))
 all_records <- list()
 for(d in list.dirs(file.path(root,"evaluators"),recursive=FALSE,full.names=TRUE)) {
  files <- sort(list.files(d,pattern="^improvement_[0-9]+[.]rds$",full.names=TRUE))
  if(!length(files)) next
  last <- Inf
  for(p in files) {
   r <- readRDS(p);v <- r$parameters;c <- r$components
   stopifnot(identical(names(v),names(contract$endpoint)),all(v>=contract$lower-1e-12),all(v<=contract$upper+1e-12),
     c$L<last,all(v[names(r$fixed)]==r$fixed),
     identical(r$contract_sha256,hash_file(file.path(study,"manifests/contract.rds"))),
     identical(r$source_manifest_sha256,hash_file(file.path(study,"manifests/input_manifest.tsv"))),
     score_near(c$L_data,c$L_b+c$L_p+c$L_n),score_near(c$L,c$L_data+c$L_prior),
     c$n_burden_obs_total==contract$expected_counts$burden,c$n_ploidy_obs_total==contract$expected_counts$ploidy,c$n_necrosis_obs_total==contract$expected_counts$necrosis)
   if(startsWith(task,"fresh_")) stopifnot(identical(r$evaluation_mode,"fresh_fork"),
     is.integer(r$model_pid),r$model_pid!=as.integer(sub("^pid_","",basename(d))))
   cfg <- contract$cfg;if(r$policy=="data_only") cfg$use_soft_prior <- FALSE
   prior <- e$compute_soft_prior_penalty(v,cfg)
   stopifnot(score_near(c$L_prior_raw,prior$L_prior_raw),score_near(c$L_prior,if(cfg$use_soft_prior) cfg$lambda_prior*prior$L_prior_raw else 0))
   ctx <- contract;ctx$cfg <- cfg;ctx$backend <- e
   independent <- independent_components(ctx,v,r$predictions)
   component_difference <- max(abs(independent-unlist(c[names(independent)])))
   stopifnot(component_difference<=1e-8)
   rho <- 10^v[['log10_rho_2N']]
   min_live <- Inf
   for(z in r$predictions$native) {
    n <- z$result
    stopifnot(all(is.finite(n$live_state_obs)),all(is.finite(n$dead_hypoxia_state_obs)),all(is.finite(n$dead_buffer_state_obs)),
      all(n$live_state_obs>=0),all(n$dead_hypoxia_state_obs>=0),all(n$dead_buffer_state_obs>=0),
      near(rowSums(n$live_state_obs),n$Ntot_live_obs),near(rowSums(n$dead_hypoxia_state_obs),n$Ntot_dead_hypoxia_obs),near(rowSums(n$dead_buffer_state_obs),n$Ntot_dead_buffer_obs),
      near(n$Ntot_live_obs+n$Ntot_dead_hypoxia_obs+n$Ntot_dead_buffer_obs,n$Ntot_total_obs),
      near(n$Vmm3_total_obs,n$Ntot_total_obs/rho),near(n$Vmm3_total_terminal,n$Ntot_total_terminal/rho),
      all(is.finite(n$O2_eff_obs)),all(is.finite(n$O2_target_obs)),near(sum(n$frac_N_live),1))
    min_live <- min(min_live,n$Ntot_live_terminal)
   }
   counter <- counter+1L
   rows[[counter]] <- data.frame(task=task,evaluator=basename(d),record=basename(p),method=r$method,worker_evaluation=r$evaluation,objective=c$L,data=c$L_data,prior=c$L_prior,independent_component_max_difference=component_difference,min_terminal_live=min_live,sha256=hash_file(p))
   all_records[[length(all_records)+1L]] <- r
   last <- c$L
  }
  ptr <- read_kv(file.path(d,"BEST.tsv"))
  stopifnot(identical(ptr$path,basename(tail(files,1))),identical(ptr$sha256,hash_file(tail(files,1))),ptr$objective==last)
 }
 accepted <- readRDS(file.path(root,"accepted_candidate.rds"))
 stopifnot(score_near(min(vapply(all_records,function(r)r$components$L,numeric(1))),accepted$checked$L),
   score_near(accepted$record$components$L,accepted$checked$L))
}
write_tsv(do.call(rbind,rows),file.path(study,"tables/improvement_archive_audit.tsv"))
write_kv(list(pass=TRUE,tasks=length(args)-1L,records=counter,
  independent_component_max_difference=max(vapply(rows,function(r) r$independent_component_max_difference,numeric(1))),
  model_evaluations=0L),file.path(study,"tables/improvement_archive_audit_summary.tsv"))
cat(counter,"immutable improvement records passed independent score arithmetic, source, constraint and native-state accounting checks.\n")
