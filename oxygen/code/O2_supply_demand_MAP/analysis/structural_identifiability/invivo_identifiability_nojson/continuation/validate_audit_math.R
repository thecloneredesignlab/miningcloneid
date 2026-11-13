# Original scalar vs vectorized archived-output arithmetic; zero model calls.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'audit_math.R'));contract<-readRDS(file.path(study,'manifests/contract.rds'));backend<-load_backend(study,FALSE)
paths<-list.files(file.path(study,'reference/imports'),pattern='[.]rds$',full.names=TRUE)
points<-lapply(paths,function(path)readRDS(path)$original_record)
engineering<-file.path(study,'tests/e_parallel_cont06')
for(directory in list.dirs(file.path(engineering,'cases'),recursive=FALSE,full.names=TRUE))points[[length(points)+1L]]<-readRDS(file.path(directory,'accepted_candidate.rds'))$record
for(name in c('partial_first','short_partial_first','projected_global','tiny_population'))points[[length(points)+1L]]<-readRDS(file.path(engineering,name,'accepted_candidate.rds'))$record
rows<-list()
for(i in seq_along(points)){
 r<-points[[i]];ctx<-contract;ctx$backend<-backend;ctx$cfg$parameter_table<-file.path(study,'source/seed25/parameter_table_input.csv');if(r$policy=='data_only')ctx$cfg$use_soft_prior<-FALSE
 t<-proc.time()[[3]];old<-independent_components(ctx,r$parameters,r$predictions);slow<-proc.time()[[3]]-t
 t<-proc.time()[[3]];fast<-independent_components_vectorized(ctx,r$parameters,r$predictions);quick<-proc.time()[[3]]-t
 gap<-max(abs(old-fast));canonical_gap<-max(abs(fast-unlist(r$components[names(fast)])))
 stopifnot(gap<1e-10,canonical_gap<1e-8)
 rows[[i]]<-data.frame(case=i,policy=r$policy,fixed_coordinates=length(r$fixed),scalar_seconds=slow,vector_seconds=quick,arithmetic_gap=gap,canonical_gap=canonical_gap)
}
# Sparse/tiny masses exercise stable log-sum-exp without changing the data.
for(fractions in list(c(1,0,0),c(1e-300,.2,.8)))for(obs in list(c(-2,0,1,4),c(.5,1.1,2.9))){
 fractions<-fractions/sum(fractions);means<-c(0,1,3);sigma<-.12
 old<-vapply(obs,function(o){v<-log(pmax(fractions,1e-300))+dnorm(o,means,sigma,log=TRUE);m<-max(v);-(m+log(sum(exp(v-m))))},numeric(1))
 stopifnot(max(abs(old-mixture_nll_vectorized(obs,means,fractions,sigma)))<1e-10)
}
tab<-do.call(rbind,rows);write_tsv(tab,file.path(study,'tables/audit_math_validation_cases.tsv'))
write_kv(list(pass=TRUE,archived_cases=nrow(tab),analytic_cases=4L,maximum_arithmetic_gap=max(tab$arithmetic_gap),maximum_canonical_gap=max(tab$canonical_gap),scalar_seconds=sum(tab$scalar_seconds),vector_seconds=sum(tab$vector_seconds),input_manifest_sha256=hash_file(file.path(study,'manifests/input_manifest.tsv')),helper_sha256=hash_file(file.path(entry,'audit_math.R')),original_adapter_sha256=hash_file(file.path(entry,'adapter.R')),model_evaluations=0L),file.path(study,'tables/audit_math_validation.tsv'))
assert_namespaces();cat(nrow(tab),'archived native cases and4analytic cases passed; maximum arithmetic gap',max(tab$arithmetic_gap),'; zero model calls.\n')
