# Metadata-only negative-access tests. No data values, simulation or RNG draws.
a<-commandArgs(TRUE);study<-normalizePath(a[1]);out<-file.path(study,'verification/C0');paths<-readLines(file.path(out,'alias_canary_paths.txt'))
checks<-list()
for(i in seq_along(paths)){
 stopifnot(file.exists(paths[i]))
 denied<-tryCatch({readBin(paths[i],'raw',n=1);FALSE},error=function(e)TRUE)
 checks[[i]]<-data.frame(case=i,path=paths[i],denied=denied)
}
write.table(do.call(rbind,checks),file.path(out,'alias_access_checks.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
stopifnot(all(vapply(checks,function(x)isTRUE(x$denied),logical(1))))
