# Metadata-only preparation: source bridge, original records, precise task registry.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2L)
study<-normalizePath(args[1]);old<-normalizePath(args[2]);entry<-dirname(normalizePath(sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])))
source(file.path(entry,'adapter.R'));e<-load_backend(study,FALSE)
c<-readRDS(file.path(study,'manifests/contract.rds'));pc<-readRDS(file.path(study,'reference/pilot_contract.rds'))
stopifnot(identical(c$cfg,pc$cfg),identical(c$scenarios,pc$scenarios))
for(label in c('pilot','formal')) {
 m<-read_tsv(file.path(study,'prior',label,'input_manifest.tsv'));take<-startsWith(m$path,'source/numerical/')|startsWith(m$path,'source/inputs/')
 stopifnot(sum(take)==8L,all(vapply(file.path(study,m$path[take]),hash_file,character(1))==m$sha256[take]))
}
cfg<-c$cfg;cfg$parameter_table<-file.path(study,'source/seed25/parameter_table_input.csv');cfg$use_soft_prior<-TRUE
old_specs<-readRDS(file.path(old,'manifests/task_specs.rds'))
original_root<-function(path) {
 if(startsWith(path,'reference_external/'))return(list(family='archive',id='pilot_archive_h25'))
 task<-strsplit(path,'/',fixed=TRUE)[[1]][2];seen<-character()
 repeat {
  stopifnot(!task%in%seen,!is.null(old_specs[[task]]));seen<-c(seen,task);sp<-old_specs[[task]]
  if(sp$mode=='global')return(list(family=if(sp$init_mode=='uniform')'uniform' else 'archive',id=task))
  if(nzchar(sp$source_task)){task<-sp$source_task;next}
  stopifnot(sp$mode=='local',sp$selector=='baseline_probe')
  task<-readRDS(file.path(old,'fits',task,'actual_local_start.rds'))$resolution$source
 }
}
target_value<-function(r,policy) r$components$L_data+if(policy=='production') cfg$lambda_prior*e$compute_soft_prior_penalty(r$parameters,cfg)$L_prior_raw else 0
pilot_d<-readRDS(file.path(old,'reference_external/data_only_pilot_constrained.rds'))
pilot_p<-readRDS(file.path(old,'reference_external/production_pilot_constrained.rds'))
best<-list(production=list(record=pilot_p,path='reference_external/production_pilot_constrained.rds',value=target_value(pilot_p,'production')),
 data_only=list(record=pilot_d,path='reference_external/data_only_pilot_constrained.rds',value=target_value(pilot_d,'data_only')))
# Complete retained improvement index: cross-policy comparisons use unchanged prior arithmetic.
index<-do.call(rbind,lapply(list.files(file.path(old,'tables/record_audits'),pattern='_records[.]tsv$',full.names=TRUE),read_tsv))
for(policy in names(best)) {
 eligible<-if(policy=='data_only') which(index$data<best[[policy]]$value) else which(index$data<best[[policy]]$value)
 for(k in eligible) {
  row<-index[k,];scope<-if(grepl('_directions$',row$task)) 'tests' else 'fits'
  rel<-file.path(scope,row$task,'evaluators',row$evaluator,row$record)
  r<-readRDS(file.path(old,rel));v<-target_value(r,policy)
  if(v<best[[policy]]$value)best[[policy]]<-list(record=r,path=rel,value=v)
 }
}
imports<-list();rows<-list()
add_import<-function(id,r,path,policy,root,family) {
 stopifnot(all(r$parameters>=c$lower),all(r$parameters<=c$upper),identical(names(r$parameters),names(c$endpoint)))
 z<-list(original_record=r,original_relative_path=path,original_file_sha256=hash_file(file.path(old,path)),
  prior_study_manifest_sha256=hash_file(file.path(old,'manifests/input_manifest.tsv')),original_root_id=original_root(path)$id,target_policy=policy,
  expected_objective=target_value(r,policy),root_id=root,root_family=family)
 atomic_rds(z,file.path(study,'reference/imports',paste0(id,'.rds')))
 rows[[id]]<<-data.frame(import_id=id,target_policy=policy,expected_objective=z$expected_objective,
  root_id=root,root_family=family,original_path=path,original_file_sha256=z$original_file_sha256,
  import_sha256=hash_file(file.path(study,'reference/imports',paste0(id,'.rds'))))
 imports[[id]]<<-z
}
for(policy in names(best)) {
 b<-best[[policy]]
 stopifnot(original_root(b$path)$family=='archive')
 cat('Selected',policy,'anchor:',b$path,'target Q',format(b$value,digits=17),'root',original_root(b$path)$id,'\n')
 add_import(paste0(policy,'_anchor'),b$record,b$path,policy,paste0('archive_h25_',policy),'archive')
 suffix<-if(policy=='production') 'u1' else 'u2'
 path<-paste0('tests/b_',policy,'_',suffix,'_directions/accepted_candidate.rds');r<-readRDS(file.path(old,path))$record
 add_import(paste0(policy,'_uniform_old'),r,path,policy,paste0('formal_',policy,'_',suffix),'uniform')
}
write_tsv(do.call(rbind,rows),file.path(study,'reference/imports.tsv'))
atomic_rds(list(prior_lifetime_core_seconds=222988L,prior_formal_core_seconds=192497L,repair_cap_core_seconds=129600L,
 profile_cap_core_seconds=316800L,formal_cap_core_seconds=648000L,lifetime_cap_core_seconds=691200L,
 score_tolerance=1e-8,prediction_absolute_tolerance=1e-8,prediction_relative_tolerance=1e-8,
 baseline_recovery=.002,direction_improvement=.002,profile_recovery=.01,steps=c(.001,.0001)),file.path(study,'manifests/protocol.rds'))
specs<-list()
add<-function(task,role,mode,policy='production',root_id='',root_family='reference',selector='reference',import_id='',source_task='',required_check='',density_index=0L,seed=1L,init_mode='uniform',sigma_frac=.08,stage='engineering') {
 specs[[task]]<<-list(task=task,role=role,mode=mode,policy=policy,root_id=root_id,root_family=root_family,family=root_family,
  selector=selector,import_id=import_id,source_task=source_task,required_check=required_check,density_index=as.integer(density_index),
  fixed=if(density_index>0L)setNames(c$rho_grid[density_index],'log10_rho_2N') else NULL,
  seed=as.integer(seed),start_seed=0L,init_mode=init_mode,sigma_frac=sigma_frac,NP=256L,itermax=1000L,reltol=1e-7,steptol=1000L,stage=stage)
}
add('engineering_parallel','engineering_parallel','engineering_parallel',seed=39100L)
for(policy in c('production','data_only')) {
 root<-imports[[paste0(policy,'_anchor')]]$root_id
 for(fixed in c(FALSE,TRUE))add(paste0('engineering_',policy,'_',if(fixed)'fixed' else 'free'),'engineering_serial','engineering_serial',policy,root,'archive','import',paste0(policy,'_anchor'),density_index=if(fixed)4L else 0L)
 for(kind in c('anchor','uniform_old')) {
  id<-paste0(policy,'_',kind);im<-imports[[id]]
  add(paste0('reference_',id),'source_replay','reference',policy,im$root_id,im$root_family,'import',id)
 }
 for(kind in c('anchor','uniform_old','uniform_new')) {
  prefix<-paste0('r_',policy,'_',kind)
  if(kind=='uniform_new') {
   root<-paste0('new_uniform_',policy);family<-'uniform';stage<-'extended'
   add(paste0(prefix,'_global'),'uniform_global','global',policy,root,family,seed=if(policy=='production')40101L else 40201L,stage=stage)
   add(paste0(prefix,'_local'),'uniform_local','local',policy,root,family,'task',source_task=paste0(prefix,'_global'),stage=stage)
  } else {
   im<-imports[[paste0(policy,'_',kind)]];root<-im$root_id;family<-im$root_family;stage<-'cheap'
   add(paste0(prefix,'_local'),if(kind=='anchor')'incumbent_local' else 'uniform_local','local',policy,root,family,'import',paste0(policy,'_',kind),stage=stage)
  }
  add(paste0(prefix,'_fresh'),'fresh','fresh',policy,root,family,'task',source_task=paste0(prefix,'_local'),stage=stage)
  add(paste0(prefix,'_directions'),'directions','directions',policy,root,family,'task',source_task=paste0(prefix,'_local'),required_check=paste0(prefix,'_fresh'),stage=stage)
 }
 # Conditional alternative identities share finite role caps; selectors are frozen at dispatch.
 for(kind in c('uniform_old','uniform_new')) {
  root<-if(kind=='uniform_old')imports[[paste0(policy,'_uniform_old')]]$root_id else paste0('new_uniform_',policy)
  add(paste0('r_',policy,'_reserve_',kind,'_global'),'reserve_global','global',policy,root,'uniform','root_selected',seed=if(policy=='production')41101L else 41201L,init_mode='hybrid',sigma_frac=.15,stage='reserve')
 }
 for(k in 1:2)for(family in c('archive','uniform')) {
  prefix<-paste0('r_',policy,'_reserve_',family,'_',k)
  add(paste0(prefix,'_local'),'reserve_local','local',policy,'dispatch',family,'root_selected',stage='reserve')
  add(paste0(prefix,'_fresh'),'reserve_fresh','fresh',policy,'dispatch',family,'task',source_task=paste0(prefix,'_local'),stage='reserve')
  add(paste0(prefix,'_directions'),'reserve_directions','directions',policy,'dispatch',family,'task',source_task=paste0(prefix,'_local'),required_check=paste0(prefix,'_fresh'),stage='reserve')
 }
 for(i in 1:7) {
  prefix<-paste0('d_',policy,'_',i);root<-paste0('density_',policy,'_',i)
  add(paste0(prefix,'_global'),'density_global','global',policy,root,'profile','baseline_best',density_index=i,seed=if(policy=='production')32100L+i else 32200L+i,init_mode='hybrid',stage='profile')
  add(paste0(prefix,'_local'),'density_global_local','local',policy,root,'profile','task',source_task=paste0(prefix,'_global'),density_index=i,stage='profile')
  for(direction in c('asc','desc')) {
   first<-if(direction=='asc')i==1L else i==7L;previous<-if(direction=='asc')i-1L else i+1L
   add(paste0(prefix,'_',direction),'density_continuation_local','local',policy,paste0('density_',policy,'_',direction),'profile',if(first)'baseline_best' else 'task',source_task=if(first)'' else paste0('d_',policy,'_',previous,'_',direction),density_index=i,stage='profile')
  }
  add(paste0(prefix,'_fresh'),'density_fresh','profile_fresh',policy,root,'profile','density_selected',density_index=i,stage='profile')
  add(paste0(prefix,'_reserve_global'),'reserve_density_global','global',policy,root,'profile','density_selected',density_index=i,seed=if(policy=='production')34100L+i else 34200L+i,init_mode='hybrid',stage='profile_reserve')
  add(paste0(prefix,'_reserve_local'),'reserve_density_local','local',policy,root,'profile','density_selected',density_index=i,stage='profile_reserve')
  add(paste0(prefix,'_reserve_fresh'),'reserve_density_fresh','profile_fresh',policy,root,'profile','density_selected',density_index=i,stage='profile_reserve')
  add(paste0(prefix,'_reserve_directions'),'reserve_density_directional','directions',policy,root,'profile','density_selected',density_index=i,stage='profile_reserve')
 }
}
atomic_rds(specs,file.path(study,'manifests/task_specs.rds'))
fields<-c('task','role','mode','policy','root_id','root_family','family','selector','import_id','source_task','required_check','density_index','seed','init_mode','sigma_frac','stage')
write_tsv(do.call(rbind,lapply(specs,function(s)as.data.frame(s[fields],stringsAsFactors=FALSE))),file.path(study,'manifests/task_plan.tsv'))
write_kv(list(pass=TRUE,model_evaluations=0L,prior_lifetime_core_seconds=222988L,prior_formal_core_seconds=192497L,registry_alternatives=length(specs)),file.path(study,'tables/metadata_preparation.tsv'))
assert_namespaces();cat('Prepared precise contract and',length(specs),'bounded task alternatives; zero model calls.\n')
