#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
input <- args[[1L]]; out <- args[[2L]]
script <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
here <- dirname(normalizePath(script))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(jsonlite))
sourceCpp(file.path(here,"perron_high_precision_safe.cpp"), cacheDir=Sys.getenv("EFAST_PERRON_CACHE"), rebuild=FALSE)
manifest <- fromJSON(file.path(input,"matrix_manifest.json"),simplifyVector=FALSE)
reports <- list()
for (entry in manifest) {
  M <- readRDS(file.path(input,paste0(entry$tag,".matrix.rds")))
  ngrid <- read.delim(gzfile(file.path(input,paste0(entry$tag,".vector.tsv.gz"))))$N
  begin <- proc.time()[["elapsed"]]
  p50 <- perron_high_precision(M,50L)
  p100 <- perron_high_precision(M,100L)
  ploidy50 <- sum(ngrid*p50$vector)/22
  ploidy100 <- sum(ngrid*p100$vector)/22
  stopifnot(all(p50$vector>=0), all(p100$vector>=0),
    abs(ploidy50-ploidy100)<1e-12, abs(p50$lambda-p100$lambda)<1e-12,
    max(abs(M%*%p50$vector-p50$lambda*p50$vector))<1e-12)
  reports[[length(reports)+1L]] <- list(tag=entry$tag, original=entry$original,
    corrected_ploidy=ploidy50, corrected_growth=p50$lambda,
    ploidy_precision_difference=abs(ploidy50-ploidy100),
    growth_precision_difference=abs(p50$lambda-p100$lambda),
    original_ploidy_difference=abs(ploidy50-entry$original$dominant_mean_ploidy),
    original_growth_difference=abs(p50$lambda-entry$original$dominant_growth_rate),
    high_precision_residual=p50$residual, leading_root_bound_width=p50$upper_minus_lambda,
    iterations50=p50$iterations, iterations100=p100$iterations,
    elapsed_seconds=proc.time()[["elapsed"]]-begin)
  cat(entry$tag,"iterations",p50$iterations,"ploidy delta",ploidy50-entry$original$dominant_mean_ploidy,"\n")
}
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
write_json(list(status="passed",n_matrices=length(reports),reports=reports),out,auto_unbox=TRUE,pretty=TRUE,digits=17)
