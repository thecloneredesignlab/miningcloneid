# Controlled history test. All numerical calls use the unchanged frozen model.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=1L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"))
label <- if(length(args)>=3L) args[3] else "history_diagnostic_local"
out <- file.path(study,"tests",label);dir.create(out,recursive=TRUE,showWarnings=FALSE)
# This server only loads source/inputs. It never evaluates the model itself;
# each numerical call executes in a collected child fork of its clean memory.
cl <- parallel::makePSOCKcluster(1L,outfile=file.path(out,"clean_server.log"))
on.exit(try(parallel::stopCluster(cl),silent=TRUE))
parallel::clusterCall(cl,function(study_arg,entry_arg) {
 source(file.path(entry_arg,"adapter.R"),local=.GlobalEnv)
 assign(".clean_ctx",load_contract(study_arg,TRUE,"data_only"),envir=.GlobalEnv)
 TRUE
},study,entry)
fork_eval <- function(p,policy="data_only",with_predictions=FALSE) {
 parallel::clusterCall(cl,function(p,policy,with_predictions) {
  job <- parallel::mcparallel({
   ctx <- .clean_ctx;ctx$policy <- policy;ctx$cfg$use_soft_prior <- policy=="production"
   z <- score(ctx,p)
   list(components=z,predictions=if(with_predictions) predictions(ctx,p,TRUE) else NULL,pid=Sys.getpid())
  },mc.set.seed=FALSE,silent=TRUE)
  z <- parallel::mccollect(job,wait=TRUE)[[1]]
  if(inherits(z,"try-error")) stop(z)
  z
 },p,policy,with_predictions)[[1]]
}
ctx <- load_contract(study,TRUE,"data_only")
candidate <- readRDS(file.path(study,"fits/local_data_only/accepted_candidate.rds"))
p <- candidate$record$parameters
canonical <- fork_eval(p,"data_only",TRUE)
baseline <- fork_eval(ctx$endpoint,"production",TRUE)
stopifnot(abs(baseline$components$L-ctx$expected_score)<1e-8)
raw_first <- score(ctx,p)
rows <- list(data.frame(case="cold_parent_A",coordinate="none",raw=raw_first$L,canonical=canonical$components$L,difference=raw_first$L-canonical$components$L))
for(n in c("log10_o2_S0","log10_kappa_O","log10_eta_o2")) {
 b <- p;b[n] <- min(ctx$upper[n],p[n]+.002)
 if(b[n]==p[n]) b[n] <- max(ctx$lower[n],p[n]-.002)
 warm_b <- score(ctx,b);warm_a <- score(ctx,p)
 clean_b <- fork_eval(b);clean_a <- fork_eval(p)
 rows[[length(rows)+1L]] <- data.frame(case="B_then_A",coordinate=n,raw=warm_a$L,canonical=clean_a$components$L,difference=warm_a$L-clean_a$components$L)
 rows[[length(rows)+1L]] <- data.frame(case="A_then_B",coordinate=n,raw=warm_b$L,canonical=clean_b$components$L,difference=warm_b$L-clean_b$components$L)
 # A canonical A evaluated first seeds all of A's bins and can mask the bug.
 # Change a true matrix coefficient first, then seed the resource-perturbed B
 # before returning to A with the same matrix signature.
 cpoint <- p;cpoint['log10_lam_max'] <- min(ctx$upper['log10_lam_max'],p['log10_lam_max']+.01)
 if(cpoint['log10_lam_max']==p['log10_lam_max']) cpoint['log10_lam_max'] <- p['log10_lam_max']-.01
 invisible(score(ctx,cpoint));cold_b <- score(ctx,b);after_b <- score(ctx,p)
 rows[[length(rows)+1L]] <- data.frame(case="C_then_B",coordinate=n,raw=cold_b$L,canonical=clean_b$components$L,difference=cold_b$L-clean_b$components$L)
 rows[[length(rows)+1L]] <- data.frame(case="C_then_B_then_A",coordinate=n,raw=after_b$L,canonical=clean_a$components$L,difference=after_b$L-clean_a$components$L)
 stopifnot(abs(clean_a$components$L-canonical$components$L)<1e-8)
}
again <- fork_eval(p,"data_only",TRUE)
stopifnot(identical(component_values(canonical$components),component_values(again$components)),
 identical(canonical$predictions$tables,again$predictions$tables),
 max(abs(independent_components(ctx,p,canonical$predictions)-unlist(canonical$components[c("L_b","L_p","L_n")])))<1e-8)
write_tsv(do.call(rbind,rows),file.path(out,"history_probe.tsv"))
atomic_rds(list(candidate=candidate,canonical=canonical,repeat_canonical=again,baseline=baseline,raw_first=raw_first),file.path(out,"diagnostic.rds"))
write_kv(list(pass=TRUE,production_reference_reproduced=TRUE,isolated_repeats_identical=TRUE,
 stored_candidate_score=candidate$record$components$L,candidate_final_reevaluation=candidate$checked$L,
 canonical_candidate_score=canonical$components$L,
 checkpoint_mismatch=candidate$checked$L-candidate$record$components$L,
 maximum_history_difference=max(abs(do.call(rbind,rows)$difference)),
 allocation_cores=as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")),model_optimizer_calls=0L),file.path(out,"result.tsv"))
parallel::stopCluster(cl);cl <- NULL
finalize_payload(out)
