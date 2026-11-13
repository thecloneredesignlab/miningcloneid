# Archived optimizer-RNG/matrix replay only. No tumor trajectory/objective call.
# Use the specified R4.4.2 SIF for exact RNG/normal-generator reproduction.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);tasks<-args[-1]
entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'))
stopifnot(as.character(getRversion())=='4.4.2');contract<-readRDS(file.path(study,'manifests/contract.rds'));backend<-load_backend(study,FALSE)
jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character')
if(!length(tasks))tasks<-jobs$task[jobs$mode=='global'&jobs$state=='COMPLETED'&jobs$exit_code=='0:0']
rows<-list()
for(task in tasks){
 root<-file.path(study,'fits',task);if(!dir.exists(root))root<-file.path(study,'tests',task);verify_payload_R(root)
 paths<-list.files(root,pattern='^actual_initial_population[.]rds$',recursive=TRUE,full.names=TRUE)
 for(path in paths){
  block<-dirname(path);pop_record<-readRDS(path);spec<-readRDS(file.path(block,'actual_spec.rds'));start<-readRDS(file.path(block,'actual_start.rds'));pop<-pop_record$free;free<-pop_record$free_names
  lower<-contract$lower[free];upper<-contract$upper[free];clipped_entries<-clipped_rows<-0L
  if(is.null(spec$initial_population)){
   original<-backend$.sample_truncnorm_box
   backend$.sample_truncnorm_box<-function(n,center,lower,upper,sigma_frac=.1){
    before<-get('.Random.seed',.GlobalEnv);out<-original(n,center,lower,upper,sigma_frac);after<-get('.Random.seed',.GlobalEnv)
    assign('.Random.seed',before,.GlobalEnv)
    raw<-matrix(stats::rnorm(n*length(lower)),nrow=n,ncol=length(lower));raw<-sweep(raw,2,pmax(as.numeric(upper-lower)*sigma_frac,1e-12),'*');raw<-sweep(raw,2,as.numeric(center),'+')
    outside<-sweep(raw,2,as.numeric(lower),'<')|sweep(raw,2,as.numeric(upper),'>')
    clipped_entries<<-sum(outside);clipped_rows<<-sum(rowSums(outside)>0L);assign('.Random.seed',after,.GlobalEnv);out
   }
   assign('.Random.seed',pop_record$initial_rng,.GlobalEnv)
   rebuilt<-backend$.build_de_initialpop(nrow(pop),lower,upper,if(spec$init_mode=='uniform')NULL else start$parameters[free],spec$init_mode,.35,spec$sigma_frac)$pop
   backend$.sample_truncnorm_box<-original;stopifnot(identical(rebuilt,pop))
   semantics<-'exact_archived_initial_population_replay'
  }else{
   stopifnot(identical(spec$initial_population,pop));clipped_entries<-clipped_rows<-NA_integer_;semantics<-'inherited_population_no_new_gaussian_clipping'
  }
  at_lower<-sweep(pop,2,as.numeric(lower),'==');at_upper<-sweep(pop,2,as.numeric(upper),'==')
  stopifnot(all(sweep(pop,2,as.numeric(lower),'>=')),all(sweep(pop,2,as.numeric(upper),'<=')),ncol(pop)==20L-length(spec$fixed))
  rows[[length(rows)+1L]]<-data.frame(task=task,block=substring(block,nchar(root)+2L),mode=spec$init_mode,semantics=semantics,population=nrow(pop),free_coordinates=ncol(pop),clipped_entries=clipped_entries,clipped_rows=clipped_rows,boundary_entries=sum(at_lower|at_upper),boundary_rows=sum(rowSums(at_lower|at_upper)>0L),duplicate_rows=sum(duplicated(pop,MARGIN=1)),population_sha256=hash_file(path),pass=TRUE,model_evaluations=0L)
 }
}
if(length(rows))write_tsv(do.call(rbind,rows),file.path(study,'tables/population_geometry_audit.tsv'))
assert_namespaces();cat(length(rows),'archived population matrices validated, including raw clipping and boundary duplicates; zero model calls.\n')
