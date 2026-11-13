# Exact scalar/RDS round trips, including known cross-platform decimal edges.
args <- commandArgs(trailingOnly=TRUE);stopifnot(length(args)>=1L)
study <- normalizePath(args[1]);pkg <- dirname(normalizePath(sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])))
source(file.path(pkg,"adapter.R"))
label <- if(length(args)>=3L) args[3] else "metadata_roundtrip_local"
out <- file.path(study,"tests",label);dir.create(out,recursive=TRUE,showWarnings=FALSE)
values <- c(readRDS(file.path(study,"manifests/contract.rds"))$endpoint,
  rounding_edge=as.numeric("0x1.ef6ab1ac5bf8fp+3"),negative_edge=-as.numeric("0x1.ef6ab1ac5bf8fp+3"),
  small=.Machine$double.xmin,large=.Machine$double.xmax,zero=0,negative_zero=-0,
  infinity=Inf,negative_infinity=-Inf,missing=NA_real_,not_a_number=NaN)
set.seed(52717);values <- c(values,setNames(rnorm(1000),paste0("random_",1:1000)))
for(task in c("global_production","global_data_only")) {
 fs <- list.files(file.path(study,"fits",task,"evaluators"),pattern="^improvement_.*[.]rds$",recursive=TRUE,full.names=TRUE)
 for(f in fs) values <- c(values,setNames(readRDS(f)$components$L,paste0(task,"_",basename(dirname(f)),"_",basename(f))))
}
stopifnot(!anyDuplicated(names(values)))
write_kv(as.list(values),file.path(out,"roundtrip.tsv"))
actual <- unlist(read_kv(file.path(out,"roundtrip.tsv")))
stopifnot(identical(values,actual))
atomic_rds(values,file.path(out,"roundtrip.rds"))
restored <- readRDS(file.path(out,"roundtrip.rds"))
stopifnot(identical(values,restored))
misc <- list(literal_na="NA",empty="",missing=NA_character_,true=TRUE,false=FALSE,missing_boolean=NA,missing_integer=NA_integer_,quoted="A tab\tquote\" and newline\n")
write_kv(misc,file.path(out,"misc.tsv"))
stopifnot(identical(misc,read_kv(file.path(out,"misc.tsv"))))
write_kv(list(pass=TRUE,values=length(values),model_evaluations=0L,R_version=as.character(getRversion()),exact_ieee_roundtrip=TRUE),file.path(out,"result.tsv"))
finalize_payload(out)
