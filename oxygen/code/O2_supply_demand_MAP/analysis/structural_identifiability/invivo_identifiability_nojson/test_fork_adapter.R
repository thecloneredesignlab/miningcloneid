# Integration of fresh-equivalent evaluation; no optimizer is called.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=3L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"))
out <- file.path(study,"tests",args[3]);dir.create(out,recursive=TRUE,showWarnings=FALSE)
ctx <- load_contract(study,TRUE,"production")
stopifnot(identical(ctx$evaluation_mode,"fresh_fork"))
reference <- score(ctx,ctx$endpoint)
stopifnot(abs(reference$L-ctx$expected_score)<1e-8)
diagnostic <- readRDS(file.path(study,"tests/history_diagnostic/diagnostic.rds"))
p <- diagnostic$candidate$record$parameters
data <- load_contract(study,TRUE,"data_only")
rows <- list();pids <- integer();counter <- 0L
for(policy in c("production","data_only")) {
 cctx <- if(policy=="production") ctx else data
 first <- score(cctx,p);pids <- c(pids,get(".ivi_last_model_pid",.GlobalEnv))
 for(n in c("log10_o2_S0","log10_kappa_O","log10_eta_o2")) {
  b <- p;b[n] <- min(cctx$upper[n],p[n]+.002)
  if(b[n]==p[n]) b[n] <- p[n]-.002
  other <- score(cctx,b);pids <- c(pids,get(".ivi_last_model_pid",.GlobalEnv))
  again <- score(cctx,p);pids <- c(pids,get(".ivi_last_model_pid",.GlobalEnv))
  stopifnot(identical(component_values(first),component_values(again)))
  counter <- counter+1L
  rows[[counter]] <- data.frame(policy=policy,coordinate=n,first=first$L,after_other=again$L,difference=again$L-first$L)
 }
 pred <- predictions(cctx,p,TRUE)
 independent <- independent_components(cctx,p,pred)
 stopifnot(max(abs(independent-unlist(first[names(independent)])))<1e-8)
 if(policy=="data_only") stopifnot(abs(first$L-diagnostic$canonical$components$L)<1e-8)
 atomic_rds(list(parameters=p,components=first,predictions=pred),file.path(out,paste0(policy,"_candidate.rds")))
}
stopifnot(!any(pids==Sys.getpid()),!anyDuplicated(pids))
write_tsv(do.call(rbind,rows),file.path(out,"history_independence.tsv"))
write_kv(list(pass=TRUE,production_reference=reference$L,independent_evaluation_processes=length(pids),
 source_unchanged=TRUE,history_independence=TRUE,independent_prediction_score=TRUE,
 legacy_local_candidate_rejected=abs(diagnostic$candidate$record$components$L-diagnostic$canonical$components$L)>1e-8,
 optimizer_calls=0L,parent_pid=Sys.getpid()),file.path(out,"result.tsv"))
finalize_payload(out)
