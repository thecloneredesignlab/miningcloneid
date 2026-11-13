# Compare immutable first-evaluation snapshots with dispatched starts; no model calls.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
study <- normalizePath(args[1]);entry <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,"adapter.R"))
rows <- list()
for(root in list.dirs(file.path(study,"tests/start_identity_snapshots"),recursive=FALSE,full.names=TRUE)) {
 start <- readRDS(file.path(root,"start.rds"));first <- readRDS(file.path(root,"first.rds"))
 origin <- read_tsv(file.path(root,"origins.tsv"))
 stopifnot(hash_file(file.path(root,"start.rds"))==origin$sha256[origin$kind=="start"],
  hash_file(file.path(root,"first.rds"))==origin$sha256[origin$kind=="first"],
  identical(first$parameters,start$parameters),first$evaluation==1L,
  identical(first$provenance$start_sha256,hash_file(file.path(root,"start.rds"))))
 previous <- start$resolution$record
 same <- !is.null(previous) && max(abs(previous$parameters-first$parameters))<=1e-12
 expected <- start$resolution$expected
 if(is.null(expected) && same) expected<-previous$components$L
 gap <- if(same && !is.null(expected)) first$components$L-expected else NA_real_
 if(same)stopifnot(abs(gap)<1e-8)
 rows[[length(rows)+1L]] <- data.frame(task=basename(root),same_source_point=same,source_score=if(same) expected else NA_real_,
  fresh_start_score=first$components$L,score_gap=gap,pass=TRUE)
}
stopifnot(length(rows)>0L)
write_tsv(do.call(rbind,rows),file.path(study,"tables/local_start_identity.tsv"))
print(do.call(rbind,rows),row.names=FALSE)
