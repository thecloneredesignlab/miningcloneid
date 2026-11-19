# Privileged contract export only; no propagation and no response export.
f_setup<-function(study,out,env,check){
 b<-readRDS(file.path(study,'inputs/private/source_bundle.rds'));d<-b$design
 check('public_design_schema',identical(sort(names(b)),sort(c('design','observations','scope'))))
 check('complete20_bounds',nrow(d$bounds)==20&&setequal(d$bounds$parameter,cg_roles))
 check('two_cohorts',setequal(names(d$adapters),c('2N','4N')))
 neutral<-cg_neutral(env,d$bounds)
 for(n in cg_roles){if(n%in%names(d$cfg))check(paste0('neutral_',n),identical(as.numeric(d$cfg[[n]]),as.numeric(neutral[[n]])));k<-paste0(n,'_init');if(k%in%names(d$cfg))check(paste0('neutral_',k),identical(as.numeric(d$cfg[[k]]),as.numeric(neutral[[n]])))}
 check('no_active_fixed_parameter',!any(cg_roles%in%names(d$fixed)))
 check('v2_death',identical(d$cfg$ploidy_O2_death,'ploidy_related'))
 edges<-list();segments<-list()
 for(cohort in names(d$adapters)){
  a<-d$adapters[[cohort]];check(paste0('adapter_',cohort),a$cohort==cohort)
  ids<-vapply(a$segments,`[[`,character(1),'segment_id');check(paste0('unique_segments_',cohort),!anyDuplicated(ids))
  for(i in seq_along(a$segments)){
   s<-a$segments[[i]];check(paste0('no_observed_',cohort,'_',i),!'observed'%in%names(s))
   parent<-s$parent_segment_id;has<-!is.null(parent)&&length(parent)==1&&!is.na(parent)&&nzchar(parent)
   if(has)check(paste0('topological_parent_',cohort,'_',i),parent%in%ids[seq_len(i-1L)])
   kind<-if(has)'explicit_parent'else if(i>1)'production_previous_fraction_fallback'else'root'
   raw_initial<-as.numeric(s$initial_cells);initial_valid<-length(raw_initial)==1L&&is.finite(raw_initial)&&raw_initial>0
   initial_effective<-if(initial_valid)raw_initial else as.numeric(d$cfg$init_total_size)
   raw_final<-as.numeric(s$final_cells);final_valid<-length(raw_final)==1L&&is.finite(raw_final)&&raw_final>0
   next_initial<-if(i<length(a$segments))as.numeric(a$segments[[i+1L]]$initial_cells)else NA_real_
   next_valid<-length(next_initial)==1L&&is.finite(next_initial)&&next_initial>0
   target_effective<-if(final_valid)raw_final else if(next_valid)next_initial else NA_real_
   edges[[length(edges)+1L]]<-data.frame(cohort=cohort,ordinal=i,parent=if(has)parent else if(i>1)ids[i-1]else NA_character_,child=s$segment_id,kind=kind,initial_cells=raw_initial,initial_cells_effective=initial_effective,initial_source=if(initial_valid)'recorded_initial_cells'else'production_cfg_init_total_size_fallback',final_cells=raw_final,target_effective=target_effective,target_source=if(final_valid)'recorded_final_cells'else if(next_valid)'next_array_segment_recorded_initial'else'production_last_positive_day_fallback',oxygen_pct=as.numeric(s$oxygen_pct),duration_days=as.numeric(s$duration_days))
   segments[[length(segments)+1L]]<-data.frame(cohort=cohort,segment=s$segment_id,index=seq_along(s$obs_days_local),day=as.numeric(s$obs_days_local))
  }
 }
 edges<-do.call(rbind,edges);check('finite_effective_registered_counts',all(is.finite(edges$initial_cells_effective)&edges$initial_cells_effective>0))
 # These are explicitly enumerated public fields; no historical response payload.
 contract<-list(design=d,edges=edges,scope='Observation-free feasibility under frozen known conditioning;14 mechanisms plus6 nuisance registry')
 sf_rds(contract,file.path(out,'feasibility_contract.rds'));sf_emit(edges,file.path(out,'protocol_edges.tsv'));sf_emit(do.call(rbind,segments),file.path(out,'observation_days.tsv'));sf_emit(d$bounds,file.path(out,'exact_bound_exports.tsv'))
 raw<-sf_tab(file.path(study,'inputs/diagnostic/raw_bank.tsv'));u<-as.matrix(raw[,paste0('u',1:19)]);storage.mode(u)<-'double'
 check('raw1000_dimensions',identical(dim(u),c(1000L,19L)))
 check('all_raw_unit_hashes',all(vapply(seq_len(nrow(u)),function(i)sf_hash(as.numeric(u[i,]))==raw$candidate_sha[i],logical(1))))
 sf_rds(u,file.path(out,'old_units.rds'))
 RNGkind("L'Ecuyer-CMRG","Inversion","Rejection");set.seed(20261118);seed<-.Random.seed;rows<-list()
 for(origin in paste0('F_',LETTERS[1:4])){sf_rds(list(origin=origin,state=seed),file.path(out,'rng',paste0(origin,'.rds')));rows[[length(rows)+1L]]<-data.frame(origin=origin,RNG_sha=sf_hash(seed));seed<-parallel::nextRNGStream(seed)}
 sf_emit(do.call(rbind,rows),file.path(out,'rng_registry.tsv'))
 (function(){status<-system2('gcc',c('-shared','-fPIC','-O2',shQuote(file.path(study,'release/read_guard.c')),'-ldl','-o',shQuote(file.path(out,'read_guard.so'))));check('native_read_guard_compiled',status==0)})()
}
