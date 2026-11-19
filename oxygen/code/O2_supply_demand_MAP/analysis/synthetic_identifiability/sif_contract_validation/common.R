# Isolated validation utilities. No fitting or synthetic-observation mode exists.
args <- commandArgs(TRUE)
stopifnot(length(args) >= 2L)
project <- normalizePath(args[1]); study <- normalizePath(args[2])
archive <- file.path(project,"oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson")
wp0 <- file.path(project,"oxygen/results/structural_identifiability/20261109_seed144_wp0")
contract <- file.path(study,"inputs/local_contract")
expected_image <- "dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372"
tab <- function(path) read.delim(path,stringsAsFactors=FALSE,check.names=FALSE)
emit <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  old <- options(digits=17);on.exit(options(old))
  tmp <- tempfile(".pending_",dirname(path));on.exit(unlink(tmp),add=TRUE)
  write.table(x,tmp,sep="\t",quote=TRUE,row.names=FALSE,na="NA",qmethod="double")
  stopifnot(file.rename(tmp,path))
}
kv <- function(x,path) emit(data.frame(key=names(x),type=vapply(x,function(v)if(is.logical(v))"boolean" else if(is.integer(v))"integer" else if(is.numeric(v))"double" else "string",character(1)),value=vapply(x,function(v)if(is.logical(v))tolower(as.character(v)) else as.character(v),character(1))),path)
guard <- function(path,preflight=FALSE) {
  forbidden <- intersect(loadedNamespaces(),c("rJava","V8","jsonlite","rjson","RJSONIO"))
  stopifnot(!length(forbidden))
  if(!preflight) stopifnot(identical(Sys.getenv("SIF_VALIDATION_IMAGE_SHA256"),expected_image),getRversion()=="4.4.2")
  rows <- lapply(sort(loadedNamespaces()),function(n)data.frame(namespace=n,version=as.character(utils::packageVersion(n)),loaded_path=if(n=="base")file.path(R.home("library"),"base") else getNamespaceInfo(asNamespace(n),"path")))
  emit(do.call(rbind,rows),file.path(path,"loaded_packages.tsv"))
}
load_defs <- function() {
  source(file.path(archive,"implementation/helpers.R"),local=.GlobalEnv)
  source(file.path(archive,"implementation/runtime.R"),local=.GlobalEnv)
}
record_environment <- function(path) {
  guard(path)
  writeLines(capture.output(sessionInfo()),file.path(path,"sessionInfo.txt"))
  kv(list(pass=TRUE,image_sha256=Sys.getenv("SIF_VALIDATION_IMAGE_SHA256"),R_version=as.character(getRversion()),library_paths=paste(.libPaths(),collapse=";"),loaded_BLAS=as.character(extSoftVersion()["BLAS"]),TZ=Sys.getenv("TZ"),locale=Sys.getlocale(),job_id=Sys.getenv("SLURM_JOB_ID"),allocated_cpus=Sys.getenv("SLURM_CPUS_PER_TASK")),file.path(path,"environment_fingerprint.tsv"))
}
