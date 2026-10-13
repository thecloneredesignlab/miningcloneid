#!/usr/bin/env Rscript
# Validate the immutable recovery adapter against one actual sampled matrix.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[[1L]]); part <- normalizePath(args[[2L]])
script <- sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[[1L]])
here <- dirname(normalizePath(script))
suppressPackageStartupMessages(library(jsonlite))
con <- gzfile(file.path(part,"outputs.tsv.gz"),"rt")
lines <- readLines(con,n=2L); close(con)
row <- read.delim(text=lines,check.names=FALSE,colClasses="character")
stopifnot(nrow(row)==1L)
dir.create(file.path(part,"recovery"),showWarnings=FALSE)
write_json(list(as.list(row[1L,,drop=FALSE])),file.path(part,"recovery/failed_points.json"),
  pretty=TRUE,auto_unbox=TRUE)
status <- system2("Rscript",c("--vanilla",shQuote(file.path(here,"repair_eigen_points_safe.R")),shQuote(root),shQuote(part)))
stopifnot(status==0L)
proof <- fromJSON(file.path(part,"recovery/corrections.json"))
stopifnot(proof$status=="passed",proof$n_points==1L,
  proof$reports$double_residual<1e-12,proof$reports$ploidy_precision_difference<1e-12,
  proof$reports$growth_precision_difference<1e-12,proof$reports$vector_precision_difference<1e-10,
  abs(proof$reports$growth_change)<1e-8)
write_json(list(status="passed",n_matrices=1L,decimal_digits=c(50L,100L),
  double_residual=proof$reports$double_residual,
  leading_root_bound_width=proof$reports$leading_root_bound_width,
  ploidy_change=proof$reports$ploidy_change,growth_change=proof$reports$growth_change,
  matrix_sha256=proof$reports$matrix_sha256),file.path(root,"smoke/perron_validation.json"),
  auto_unbox=TRUE,pretty=TRUE,digits=17)
cat("Actual sampled matrix passed 50/100-digit recovery interface validation\n")
