# Inspect archived starting-point metadata only; no model loading or evaluation.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);source(file.path(study,"releases/fresh_pipeline01/adapter.R"))
contract <- readRDS(file.path(study,"manifests/contract.rds"))
seeds <- c(25L,74L,366L,292L,392L)
base <- readRDS(file.path(study,"planning_sources/seed25/fit_config.rds"))
parameters <- list();summary <- list();differences <- list()
for(i in seq_along(seeds)) {
 seed <- seeds[i];root <- file.path(study,"planning_sources",paste0("seed",seed))
 saved <- readRDS(file.path(root,"deoptim_result.rds"))
 p <- saved$optim$bestmem
 stopifnot(identical(names(p),names(contract$endpoint)),all(is.finite(p)),
   all(p>=contract$lower-1e-12),all(p<=contract$upper+1e-12))
 cfg <- readRDS(file.path(root,"fit_config.rds"))
 changed <- union(names(base),names(cfg));changed <- changed[!vapply(changed,function(n) identical(base[[n]],cfg[[n]]),logical(1))]
 differences[[i]] <- data.frame(seed=seed,field=if(length(changed)) changed else "none")
 parameters[[i]] <- data.frame(seed=seed,parameter=names(p),transformed=unname(p))
 delta <- (p-contract$endpoint)/(contract$upper-contract$lower)
 summary[[i]] <- data.frame(seed=seed,archived_objective=as.numeric(saved$optim$bestval),
   density_cells_per_mm3=10^p[["log10_rho_2N"]],scaled_rms_distance_from_seed25=sqrt(mean(delta^2)),
   scaled_max_distance_from_seed25=max(abs(delta)),fresh_evaluated_here=FALSE,
   parameter_input_sha256=hash_file(file.path(root,"parameter_table_input.csv")))
}
summary <- do.call(rbind,summary)
stopifnot(length(unique(summary$parameter_input_sha256))==1L)
write_tsv(summary,file.path(study,"tables/planning_candidate_summary.tsv"))
write_tsv(do.call(rbind,parameters),file.path(study,"tables/planning_candidate_parameters.tsv"))
write_tsv(do.call(rbind,differences),file.path(study,"tables/planning_candidate_config_differences.tsv"))
print(summary[,1:6],row.names=FALSE)
cat("Candidate metadata inspected; zero model evaluations. Archived scores are not fresh acceptance.\n")
