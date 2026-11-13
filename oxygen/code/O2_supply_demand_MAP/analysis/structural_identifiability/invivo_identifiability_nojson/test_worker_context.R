# Regression for serialized worker evaluator metadata. No tumor evaluation.
args <- commandArgs(trailingOnly=TRUE); stopifnot(length(args)==2L)
pkg <- normalizePath(args[1]); out <- normalizePath(args[2],mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
cl <- parallel::makePSOCKcluster(2L)
z <- tryCatch(parallel::clusterCall(cl,function(pkg,out) {
 source(file.path(pkg,"adapter.R"),local=.GlobalEnv)
 source(file.path(pkg,"optimize.R"),local=.GlobalEnv)
 assign("score",function(ctx,p) list(L=sum(p*p)),.GlobalEnv)
 assign("predictions",function(...) list(mock=TRUE),.GlobalEnv)
 assign("hash_file",function(path) {
  if(file.exists(path)) return(digest::digest(file=path,algo="sha256"))
  stopifnot(startsWith(path,"/unit-study/")); "unit-source"
 },.GlobalEnv)
 ctx <- list(endpoint=c(a=1,b=2),policy="unit",study="/unit-study")
 ev <- new_evaluator(ctx,out,NULL,"unit",Inf)
 stopifnot(ev(c(1,2))==5,ev(c(.5,1))==1.25)
 list(pass=TRUE,pid=Sys.getpid())
},pkg,out),finally=parallel::stopCluster(cl))
stopifnot(all(vapply(z,`[[`,logical(1),"pass")))
source(file.path(pkg,"adapter.R"))
write_kv(list(pass=TRUE,workers=2L,tumor_evaluations=0L,global_study_variable_required=FALSE),file.path(out,"result.tsv"))
