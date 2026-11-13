# Profile arithmetic/lineage gates over closed native records; zero model calls.
args<-commandArgs(trailingOnly=TRUE);study<-normalizePath(args[1]);action<-args[2]
entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));source(file.path(entry,'helpers.R'));source(file.path(entry,'profile_rules.R'))
source(file.path(entry,'audit_math.R'));fast_audit<-install_validated_audit_math(study,entry)
verify_inputs(study)
jobs<-read.delim(file.path(study,'tables/jobs.tsv'),colClasses='character')
contract<-readRDS(file.path(study,'manifests/contract.rds'));tables<-file.path(study,'tables')
backend<-load_backend(study,FALSE)
audit_record<-function(r){
 stopifnot(r$source_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')),r$contract_sha256==hash_file(file.path(study,'manifests/contract.rds')),
  all(r$parameters>=contract$lower),all(r$parameters<=contract$upper),r$components$n_burden_obs_total==92L,r$components$n_ploidy_obs_total==6502L,r$components$n_necrosis_obs_total==6L)
 ctx<-contract;ctx$backend<-backend;ctx$cfg$parameter_table<-file.path(study,'source/seed25/parameter_table_input.csv');if(r$policy=='data_only')ctx$cfg$use_soft_prior<-FALSE
 stopifnot(max(abs(independent_components(ctx,r$parameters,r$predictions)-unlist(r$components[c('L_b','L_p','L_n')])))<1e-8)
}
closed<-function(j)sub('[+ ].*$','',j$state)%in%c('COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED')
successful<-function(j)j$state=='COMPLETED'&&j$exit_code=='0:0'
verify_evidence<-function(path){m<-read_tsv(path);stopifnot(all(vapply(file.path(study,m$path),hash_file,character(1))==m$sha256))}
epfile<-file.path(tables,'profile_epoch.tsv')
if(action=='freeze'){
 stopifnot(all(closed(jobs[jobs$kind%in%c('baseline','profile'),,drop=FALSE])))
 gate<-read_kv(file.path(tables,'baseline_gate.tsv'));stopifnot(gate$pass,gate$input_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')))
 verify_evidence(file.path(study,gate$evidence_manifest))
 epoch<-'e0001'
 if(file.exists(epfile)){
  previous<-read_kv(epfile);pg<-read_kv(file.path(tables,'profile_gate.tsv'))
  stopifnot(pg$free_baseline_invalidated,!pg$pass)
  # Epoch identity is unbounded. A new accepted free reference requires a new
  # frozen epoch; the previous baseline and every raw delta remain unchanged.
  epoch<-sprintf('e%04d',as.integer(sub('^e','',previous$epoch))+1L)
 }
 directory<-file.path(study,'profile_epochs',epoch);stopifnot(!dir.exists(directory));dir.create(directory,recursive=TRUE)
 for(policy in c('production','data_only')){
  a<-readRDS(file.path(tables,paste0('baseline_',policy,'.rds')))
  stopifnot(a$policy==policy,length(a$record$parameters)==20L)
  # Preserve the complete native ancestry. This is the free reference for all
  # raw deltas in this epoch, even if another feasible point is found later.
  atomic_rds(a,file.path(directory,paste0('baseline_',policy,'.rds')))
 }
 file.copy(file.path(tables,'baseline_gate.tsv'),file.path(directory,'baseline_gate.tsv'))
 file.copy(file.path(study,gate$evidence_manifest),file.path(directory,'baseline_gate_evidence.tsv'))
 file.copy(file.path(study,'manifests/density_grid.tsv'),file.path(directory,'density_grid.tsv'))
 paths<-list.files(directory,full.names=TRUE);write_tsv(data.frame(path=substring(paths,nchar(study)+2L),sha256=vapply(paths,hash_file,character(1))),file.path(directory,'manifest.tsv'))
 write_kv(list(epoch=epoch,manifest=paste0('profile_epochs/',epoch,'/manifest.tsv'),manifest_sha256=hash_file(file.path(directory,'manifest.tsv')),input_manifest_sha256=hash_file(file.path(study,'manifests/input_manifest.tsv'))),epfile)
 cat('Frozen accepted free baselines; no profile allocation submitted.\n')
} else if(action=='review'){
 ep<-read_kv(epfile);stopifnot(ep$manifest_sha256==hash_file(file.path(study,ep$manifest)),ep$input_manifest_sha256==hash_file(file.path(study,'manifests/input_manifest.tsv')))
 verify_evidence(file.path(study,ep$manifest));directory<-file.path(study,'profile_epochs',ep$epoch)
 epoch_jobs<-jobs$kind=='profile'&startsWith(jobs$task,paste0('p_',ep$epoch,'_'))
 stopifnot(system2('python3',c(shQuote(file.path(entry,'audit_profile_lineage.py')),shQuote(study)))==0L)
 rows<-list();evidence<-c('tables/profile_epoch.tsv',ep$manifest,'tables/profile_lineage_audit.tsv');counterexamples<-list()
 # A closed canonical global/unchecked local point is already a feasible
 # counterexample; do not wait for profile branch acceptance to protect it.
 for(i in which(epoch_jobs&jobs$state=='COMPLETED'&jobs$exit_code=='0:0')){
  j<-jobs[i,];scope<-if(j$mode%in%c('global','local'))'fits' else 'tests';p<-file.path(study,scope,j$task)
  verify_payload_R(p);a<-readRDS(file.path(p,'accepted_candidate.rds'));audit_record(a$record);free<-readRDS(file.path(directory,paste0('baseline_',j$policy,'.rds')))
  stopifnot(abs(a$checked$L-a$record$components$L)<1e-8,a$record$source_manifest_sha256==ep$input_manifest_sha256)
  evidence<-c(evidence,paste0(scope,'/',j$task,'/payload_manifest.tsv'))
  if(a$checked$L-free$checked$L < -0.002){
   a$policy<-j$policy;a$root_id<-j$root_id;a$root_family<-'profile';a$source_task<-j$task;a$source_scope<-scope
   a$free_ancestor_root<-free$root_id;a$free_ancestor_family<-free$root_family;a$profile_epoch<-ep$epoch
   name<-paste0('counterexample_',ep$epoch,'_',j$policy,'_',j$task);atomic_rds(a,file.path(tables,paste0(name,'.rds')))
   counterexamples[[j$task]]<-data.frame(policy=j$policy,epoch=ep$epoch,profile_task=j$task,selection=name,objective=a$checked$L,free_reference=free$checked$L,raw_delta=a$checked$L-free$checked$L,free_ancestor_root=free$root_id)
   evidence<-c(evidence,paste0('tables/',name,'.rds'))
  }
 }
 for(i in which(epoch_jobs&jobs$mode=='local')){
  j<-jobs[i,];check<-jobs[jobs$task==paste0(j$task,'_check'),,drop=FALSE]
  if(!successful(j)||nrow(check)!=1L||!successful(check))next
  p<-file.path(study,'fits',j$task);d<-file.path(study,'tests',check$task);verify_payload_R(p);verify_payload_R(d)
  prediction_path<-paste0('tables/prediction_audits/',check$task,'_summary.tsv');prediction_proof<-read_kv(file.path(study,prediction_path));stopifnot(prediction_proof$pass,prediction_proof$payload_sha256==hash_file(file.path(d,'payload_manifest.tsv')));evidence<-c(evidence,prediction_path)
  a<-readRDS(file.path(p,'accepted_candidate.rds'));b<-readRDS(file.path(d,'accepted_candidate.rds'));z<-read_kv(file.path(p,'result.tsv'));dr<-read_kv(file.path(d,'result.tsv'))
  audit_record(a$record);audit_record(b$record)
  index<-as.integer(j$density_index);rho<-contract$rho_grid[index];free<-readRDS(file.path(directory,paste0('baseline_',j$policy,'.rds')))
  stopifnot(dr$pass,abs(dr$difference)<1e-8,dr$fixed_residual==0,z$free_coordinates==19L,
   identical(a$record$provenance$root_id,j$root_id),identical(b$record$provenance$root_id,j$root_id),
   a$record$parameters['log10_rho_2N']==rho,b$record$parameters['log10_rho_2N']==rho,
   all(a$record$parameters>=contract$lower),all(a$record$parameters<=contract$upper),
   b$record$source_manifest_sha256==ep$input_manifest_sha256)
  best<-if(b$checked$L<a$checked$L)b else a;best$policy<-j$policy;best$root_id<-j$root_id;best$root_family<-'profile'
  best$source_task<-if(b$checked$L<a$checked$L)check$task else j$task;best$source_scope<-if(b$checked$L<a$checked$L)'tests' else 'fits'
  branch<-sub('^density_','',j$role);stopifnot(branch%in%c('global_local','ascending','descending'))
  # Different tasks or continuations within one branch do not create recovery.
  rows[[j$task]]<-data.frame(epoch=ep$epoch,policy=j$policy,index=index,density=10^rho,branch=branch,root_id=j$root_id,task=j$task,source_task=best$source_task,objective=best$checked$L,center_objective=a$checked$L,directional_gain=dr$directional_improvement,raw_delta=best$checked$L-free$checked$L,termination=z$termination,window_gain=z$window_gain)
  name<-paste0('profile_',ep$epoch,'_',j$policy,'_r',sprintf('%02d',index),'_',branch,'_selected')
  if(!file.exists(file.path(tables,paste0(name,'.rds')))||readRDS(file.path(tables,paste0(name,'.rds')))$checked$L>best$checked$L)atomic_rds(best,file.path(tables,paste0(name,'.rds')))
  evidence<-c(evidence,substring(file.path(p,'payload_manifest.tsv'),nchar(study)+2L),substring(file.path(d,'payload_manifest.tsv'),nchar(study)+2L),paste0('tables/',name,'.rds'))
  if(best$checked$L-free$checked$L < -0.002){
   best$free_ancestor_root<-free$root_id;best$free_ancestor_family<-free$root_family;best$profile_epoch<-ep$epoch
   name<-paste0('counterexample_',ep$epoch,'_',j$policy,'_',j$task);atomic_rds(best,file.path(tables,paste0(name,'.rds')))
   counterexamples[[j$task]]<-data.frame(policy=j$policy,epoch=ep$epoch,profile_task=j$task,selection=name,objective=best$checked$L,free_reference=free$checked$L,raw_delta=best$checked$L-free$checked$L,free_ancestor_root=free$root_id)
   evidence<-c(evidence,paste0('tables/',name,'.rds'))
  }
 }
 summary<-list();tab<-if(length(rows))do.call(rbind,rows) else NULL
 if(!is.null(tab))write_tsv(tab,file.path(tables,'profile_branch_history.tsv'))
 for(policy in c('production','data_only'))for(index in seq_along(contract$rho_grid)){
  x<-if(is.null(tab))NULL else tab[tab$policy==policy&tab$index==index,,drop=FALSE]
  summary[[paste(policy,index)]]<-cbind(data.frame(policy=policy,index=index,density=10^contract$rho_grid[index]),as.data.frame(summarize_profile_anchor(x)))
 }
 summary<-do.call(rbind,summary);write_tsv(summary,file.path(tables,'profile_anchor_gates.tsv'))
 invalid<-length(counterexamples)>0L
 if(invalid){
  write_tsv(do.call(rbind,counterexamples),file.path(tables,'profile_counterexamples.tsv'))
  # This gate blocks all later profile dispatch immediately; historical epoch
  # acceptance remains immutable and is never silently reinterpreted.
  write_kv(list(pass=FALSE,release=read_kv(file.path(tables,'execution_release.tsv'))$release,input_manifest_sha256=ep$input_manifest_sha256,evidence_manifest='tables/profile_review_evidence.tsv',reason='Density-constrained feasible counterexample below frozen free baseline'),file.path(tables,'baseline_gate.tsv'))
 }
 evidence<-unique(evidence);write_tsv(data.frame(path=evidence,sha256=vapply(file.path(study,evidence),hash_file,character(1))),file.path(tables,'profile_review_evidence.tsv'))
 pending<-any(!closed(jobs[epoch_jobs,,drop=FALSE]))
 write_kv(list(pass=all(summary$pass)&&!invalid&&!pending,epoch=ep$epoch,free_baseline_invalidated=invalid,pending_parents=pending,input_manifest_sha256=ep$input_manifest_sha256,evidence_manifest='tables/profile_review_evidence.tsv',model_evaluations=0L),file.path(tables,'profile_gate.tsv'))
 print(summary,row.names=FALSE)
} else stop('Unsupported profile metadata action')
assert_namespaces()
