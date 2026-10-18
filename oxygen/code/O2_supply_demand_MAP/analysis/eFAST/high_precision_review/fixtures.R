#!/usr/bin/env Rscript
args <- commandArgs(TRUE);stopifnot(length(args)==1L)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[[1L]])
here <- dirname(normalizePath(script))
suppressPackageStartupMessages(library(Rcpp));suppressPackageStartupMessages(library(jsonlite))
sourceCpp(file.path(here,'perron_review.cpp'),cacheDir=Sys.getenv('EFAST_REVIEW_CACHE'),rebuild=FALSE)
tests <- list()
for(digits in c(50L,100L)){
  symmetric <- matrix(c(-1,1,1,-1),2)
  z <- review_perron(symmetric,c(22,44),22,digits)
  stopifnot(z$converged,abs(as.numeric(z$lambda))<1e-20,abs(as.numeric(z$ploidy)-1.5)<1e-14)
  # Exact binary-double matrix; a vector component must survive decimal serialization.
  weak <- matrix(c(-1,1e-100,1,-2),2)
  v <- review_perron(weak,c(22,44),22,digits)
  stopifnot(v$converged,as.numeric(v$vector_decimal[[2L]])>0,as.numeric(v$vector_decimal[[2L]])<1e-90)
  bad <- tryCatch(review_perron(matrix(c(-1,-1,1,-2),2),c(22,44),22,digits),error=function(e)e)
  stopifnot(inherits(bad,'error'))
  tests[[as.character(digits)]] <- list(symmetric=z,weak=v,negative_off_diagonal_rejected=TRUE)
}
write_json(list(status='passed',tests=tests),file.path(args[[1L]],'smoke_fixtures.json'),auto_unbox=TRUE,pretty=TRUE,digits=17)
cat('Known-root, tiny-component and invalid-matrix checks passed\n')
