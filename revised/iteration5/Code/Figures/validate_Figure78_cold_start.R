#!/usr/bin/env Rscript
# A real 1000-day, reduced-size forward run with forbidden-input guards.
script <- normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1]),mustWork=TRUE)
code <- dirname(script);workspace <- normalizePath(file.path(code,"../.."),mustWork=TRUE)
for(module in c("figure78_pipeline.R","figure78_model.R","figure7_simulation.R","figure8_empirical.R","figure8_simulation.R"))
  source(file.path(code,"util/analysis",module),local=globalenv())
args <- fp_args();args$profile <- "representative"
if(is.null(args$data_root))stop("Specify an isolated --data-root for cold-start acceptance")
args$workspace <- workspace
pre <- fp_init(7L,file.path(code,"run_Figure7.R"),args);args$data_root <- pre$data_root
if(any(file.exists(file.path(args$data_root,paste0("Figure",c(7L,8L)),"prepared/analysis_bundle.rds"))))
  stop("Cold-start acceptance requires empty scientific outputs, not a resume")
# Protection-only byte checks precede the guarded calculation. No reference
# values are decoded or passed to any numerical worker. This makes the public
# two-command acceptance workflow self-contained at any isolated audit root.
protected <- unlist(lapply(file.path(workspace,c("data/Figures/Figure7","data/Figures/Figure8")),
  list.files,recursive=TRUE,full.names=TRUE))
for(directory in file.path(workspace,c("Figures","manuscript/Figures"))) {
  files <- list.files(directory,pattern="\\.(pdf|png|svg)$",full.names=TRUE)
  protected <- c(protected,files[grepl("(fig(ure)?[78]|supp_fig[78])",basename(files),ignore.case=TRUE)])
}
protected <- sort(unique(protected[basename(protected)!=".DS_Store"]))
baseline <- data.frame(path=protected,sha256_before=unname(vapply(protected,digest::digest,"",file=TRUE,algo="sha256")))
fp_write(baseline,file.path(dirname(args$data_root),"protected_before.tsv"))
.cold_access <- list()
.cold_record <- function(path,reader) {
  if(!is.character(path)||length(path)!=1L||is.na(path)||!file.exists(path))return(invisible(NULL))
  path <- normalizePath(path,mustWork=TRUE)
  if(grepl("(^|/)(arxives|archives|frozen_sources)(/|$)",path)||
    startsWith(path,file.path(workspace,"data/Figures/Figure7/"))||
    startsWith(path,file.path(workspace,"data/Figures/Figure8/")))
    stop("Forbidden old scientific/code input during cold start: ",path)
  .cold_access[[length(.cold_access)+1L]] <<- data.frame(reader=reader,path=path)
  invisible(NULL)
}
trace("readRDS",where=baseenv(),tracer=quote(.GlobalEnv$.cold_record(file,"readRDS")),print=FALSE)
trace("source",where=baseenv(),tracer=quote(.GlobalEnv$.cold_record(file,"source")),print=FALSE)
trace("sys.source",where=baseenv(),tracer=quote(.GlobalEnv$.cold_record(file,"sys.source")),print=FALSE)
trace("fread",where=asNamespace("data.table"),tracer=quote(.GlobalEnv$.cold_record(
  if(!is.null(file))file else input,"fread")),print=FALSE)
for(figure in c(7L,8L)) {
  a <- fp_init(figure,file.path(code,paste0("run_Figure",figure,".R")),args)
  if(file.exists(file.path(a$run,"prepared/analysis_bundle.rds")))
    stop("Cold-start acceptance requires empty scientific outputs, not a resume: ",a$run)
  fp_prepare(a);fp_simulate(a);fp_aggregate(a);fp_validate(a)
}
untrace("readRDS",where=baseenv());untrace("source",where=baseenv());untrace("sys.source",where=baseenv())
untrace("fread",where=asNamespace("data.table"))
log <- unique(data.table::rbindlist(.cold_access))
fp_write(log,file.path(dirname(args$data_root),"cold_start_access_log.tsv"))
fp_write(data.frame(status="representative_cold_start_passed",forbidden_result_or_archive_reads=0L,
  traced_readers="R readRDS/fread/source/sys.source",traced_files=nrow(log),full_hpc_jobs=0L),
  file.path(dirname(args$data_root),"cold_start_acceptance.tsv"))
message("PASS: guarded cold starts; no retained Figure7/8 result or archived-code input")
