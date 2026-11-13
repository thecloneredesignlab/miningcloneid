# Full archived native predictive-output identity; zero model calls.
# Cache hit/build statistics are retained as diagnostics, not predictive quantities.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);tasks<-args[-1]
entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'));jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character')
if(!length(tasks))tasks<-jobs$task[jobs$state=='COMPLETED'&jobs$exit_code=='0:0'&jobs$mode%in%c('fresh','directions','profile_fresh','replay')]
all_rows<-list()
for(task in tasks){
 j<-jobs[jobs$task==task,,drop=FALSE];stopifnot(nrow(j)==1L,j$state=='COMPLETED',j$exit_code=='0:0');root<-file.path(study,'tests',task);verify_payload_R(root)
 payload<-hash_file(file.path(root,'payload_manifest.tsv'));rows<-list()
 if(j$mode=='replay'){
  paths<-list.files(root,pattern='^original_import[.]rds$',recursive=TRUE,full.names=TRUE)
  points<-lapply(paths,function(path){p<-dirname(path);start<-readRDS(path);a<-readRDS(file.path(p,'accepted_candidate.rds'));list(point=basename(p),old=start$record,predictions=a$record$predictions,components=a$record$components,parameters=a$record$parameters)})
 }else{
  start<-readRDS(file.path(root,'actual_start.rds'));fresh<-readRDS(file.path(root,'fresh_acceptance.rds'))
  points<-list(list(point='tested_center',old=start$record,predictions=fresh$predictions,components=fresh$components,parameters=fresh$parameters))
 }
 stopifnot(length(points)>0L)
 for(point in points){
  old<-point$old;now<-point$predictions;before<-old$predictions
  stopifnot(identical(point$parameters,old$parameters),abs(point$components$L-old$components$L)<1e-8,
   isTRUE(all.equal(now$tables,before$tables,tolerance=1e-8)),length(now$native)==length(before$native))
  maximum<-0;cache_identity<-TRUE;fields_total<-0L
  for(i in seq_along(now$native)){
   a<-now$native[[i]]$result;b<-before$native[[i]]$result;stopifnot(identical(names(a),names(b)))
   fields<-names(a)[!startsWith(names(a),'cache_')];diagnostics<-setdiff(names(a),fields)
   stopifnot(length(fields)>0L,isTRUE(all.equal(a[fields],b[fields],tolerance=1e-8)))
   maximum<-max(maximum,max(abs(unlist(a[fields],use.names=FALSE)-unlist(b[fields],use.names=FALSE))))
   cache_identity<-cache_identity&&isTRUE(all.equal(a[diagnostics],b[diagnostics],tolerance=0));fields_total<-fields_total+length(fields)
  }
  rows[[length(rows)+1L]]<-data.frame(task=task,point=point$point,policy=j$policy,scenarios=length(now$native),native_predictive_fields=fields_total,score_gap=point$components$L-old$components$L,native_max_absolute_gap=maximum,cache_diagnostic_identity=cache_identity,pass=TRUE,payload_sha256=payload,model_evaluations=0L)
 }
 tab<-do.call(rbind,rows);all_rows[[task]]<-tab;write_tsv(tab,file.path(study,'tables/prediction_audits',paste0(task,'_records.tsv')))
 write_kv(list(pass=TRUE,task=task,payload_sha256=payload,points=nrow(tab),native_max_absolute_gap=max(tab$native_max_absolute_gap),cache_diagnostic_identity=all(tab$cache_diagnostic_identity),model_evaluations=0L),file.path(study,'tables/prediction_audits',paste0(task,'_summary.tsv')))
 cat(task,':',nrow(tab),'full native output identities; maximum absolute gap',max(tab$native_max_absolute_gap),'\n')
}
if(length(all_rows))write_tsv(do.call(rbind,all_rows),file.path(study,'tables/prediction_identity_summary.tsv'))
assert_namespaces()
