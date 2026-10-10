#!/usr/bin/env Rscript
# Recompute only invalid points on their unchanged canonical fixed-O2 matrices.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[[1L]]); part <- normalizePath(args[[2L]])
script <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
here <- dirname(normalizePath(script))
source(file.path(here,"../../simulation/o2/fixed_o2/run_fixed_o2_simulation.R"),local=TRUE)
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(jsonlite))
sourceCpp(file.path(here,"perron_high_precision.cpp"),cacheDir=Sys.getenv("EFAST_PERRON_CACHE"),rebuild=FALSE)
plan <- fromJSON(file.path(root,"slurm/submission_plan.json"))
metadata <- fromJSON(file.path(part,"metadata.json"))
points <- fromJSON(file.path(part,"recovery/failed_points.json"),simplifyVector=FALSE)
samples <- read.delim(gzfile(file.path(part,"samples.tsv.gz")),check.names=FALSE)
fit <- file.path(plan$fit_root,metadata$fit_seed)
params <- read.delim(file.path(fit,"best_params.tsv"))
base <- as.list(setNames(as.numeric(params$value),params$parameter))
config <- readRDS(file.path(fit,"fit_config.rds"))
cfg <- prepare_sim_cfg(config,list(),0,prepare_run_params(base,"invivo",config,0))
hash_file <- function(path) strsplit(system2("sha256sum",shQuote(path),stdout=TRUE)," ")[[1L]][[1L]]
reports <- list()
for (i in seq_along(points)) {
  row <- points[[i]]
  stopifnot(row$status=="ok",is.finite(as.numeric(row$dominant_mean_ploidy)),
    is.finite(as.numeric(row$dominant_growth_rate)))
  id <- as.integer(row$sample_id); oxygen <- as.numeric(row$O2_pct)
  values <- base
  for (name in setdiff(names(samples),"sample_id")) values[[name]] <- samples[[name]][[id]]
  rp <- prepare_run_params(values,"invivo",cfg,0)
  fm <- fixo2_fixed_matrix(globalenv(),cfg,rp,oxygen)
  p50 <- perron_high_precision(fm$M,50L)
  p100 <- perron_high_precision(fm$M,100L)
  ploidy50 <- sum(fm$ngrid*p50$vector)/cfg$N_UNIT
  ploidy100 <- sum(fm$ngrid*p100$vector)/cfg$N_UNIT
  double_residual <- max(abs(fm$M%*%p50$vector-p50$lambda*p50$vector))
  stopifnot(all(p50$vector>=0),all(p100$vector>=0),abs(sum(p50$vector)-1)<1e-12,
    abs(ploidy50-ploidy100)<1e-12,abs(p50$lambda-p100$lambda)<1e-12,
    sum(abs(p50$vector-p100$vector))<1e-10,double_residual<1e-12,
    abs(p50$lambda-as.numeric(row$dominant_growth_rate))<1e-8)
  matrix_path <- file.path(part,"recovery",sprintf("point%04d_matrix.rds",i))
  saveRDS(fm$M,matrix_path)
  reports[[i]] <- list(sample_id=id,O2_pct=oxygen,original=row,
    dominant_mean_ploidy=ploidy50,dominant_growth_rate=p50$lambda,
    eigenvector_nonnegative=TRUE,status="ok",matrix_sha256=hash_file(matrix_path),
    matrix_file=basename(matrix_path),vector50=p50$vector,vector100=p100$vector,
    double_residual=double_residual,high_precision_residual=p50$residual,
    leading_root_bound_width=p50$upper_minus_lambda,iterations50=p50$iterations,
    iterations100=p100$iterations,ploidy_precision_difference=abs(ploidy50-ploidy100),
    growth_precision_difference=abs(p50$lambda-p100$lambda),
    vector_precision_difference=sum(abs(p50$vector-p100$vector)),
    ploidy_change=ploidy50-as.numeric(row$dominant_mean_ploidy),
    growth_change=p50$lambda-as.numeric(row$dominant_growth_rate),
    spectral_gap_source="unchanged canonical double-precision spectrum diagnostic")
}
write_json(list(status="passed",n_points=length(reports),reports=reports),
  file.path(part,"recovery/corrections.json"),pretty=TRUE,auto_unbox=TRUE,digits=17)
cat("Certified",length(reports),"points at 50 and 100 decimal digits\n")
