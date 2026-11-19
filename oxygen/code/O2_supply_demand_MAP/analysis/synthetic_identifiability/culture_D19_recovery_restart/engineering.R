# Privileged engineering is isolated from blind origin workers.
d_setup<-function(study,out,env,check){
 contract<-readRDS(file.path(study,'inputs/public/contract.rds'));check('contract_digest',sf_hash(contract)=='3d765c1cd021deff5ebe32a48140200b3786cdcf05eea9c4943afa2c188ab2d4')
 # Freeze every bank/order and all streams before opening D19 responses/scale.
 locks<-list();populations<-list()
 for(letter in LETTERS[1:4]){
  origin<-paste0('D_',letter);path<-file.path(study,'inputs/banks',paste0(origin,'.rds'));bank<-readRDS(path);u<-do.call(rbind,lapply(bank$members,`[[`,'unit'))
  check(paste0(origin,'_bank_identity'),bank$origin==paste0('F_',letter)&&bank$contract_sha==sf_hash(contract)&&identical(dim(u),c(200L,19L)))
  check(paste0(origin,'_native_member_keys'),all(vapply(bank$members,function(m)identical(sf_hash(m$unit),m$key),logical(1))))
  sf_rds(list(origin=origin,population=u,keys=vapply(bank$members,`[[`,character(1),'key')),file.path(out,'populations',paste0(origin,'.rds')))
  populations[[origin]]<-u;locks[[origin]]<-list(bank_sha=digest::digest(file=path,algo='sha256'),population_sha=sf_hash(u))
 }
 RNGkind("L'Ecuyer-CMRG","Inversion","Rejection");set.seed(20261119);seed<-.Random.seed;rngs<-list();rows<-list()
 for(origin in c(paste0('D_',LETTERS[1:6]),'F_E','F_F','isolated_tests')){
  rngs[[origin]]<-seed;sf_rds(list(origin=origin,state=seed,RNG_sha=sf_hash(seed)),file.path(out,'rng',paste0(origin,'.rds')))
  rows[[length(rows)+1L]]<-data.frame(origin=origin,RNG_sha=sf_hash(seed));seed<-parallel::nextRNGStream(seed)
 }
 sf_emit(do.call(rbind,rows),file.path(out,'rng_registry.tsv'));sf_rds(locks,file.path(out,'bank_locks_before_case.rds'))
 case<-readRDS(file.path(study,'inputs/public/D19.rds'));truth<-readRDS(file.path(study,'inputs/private/truth.rds'))
 check('case_expected_data_and_design',case$kind=='D'&&identical(case$design,contract$design)&&case$data_sha=='63d58a0f0c8fa76cbd8eeacbbfec85611c71eeedc519a26c13ddfd81d8254595'&&identical(sf_hash(case$data),case$data_sha))
 check('native_known_scale',identical(case$known_sigma_growth,truth$sigma_growth)&&is.finite(case$known_sigma_growth)&&case$known_sigma_growth>0)
 source_sha<-digest::digest(file=file.path(study,'source_manifest.tsv'),algo='sha256')
 for(origin in names(rngs)){
  id<-list(origin=origin,contract_sha=sf_hash(contract),rng_initial_sha=sf_hash(rngs[[origin]]),source_sha=source_sha)
  if(startsWith(origin,'D_'))id<-c(id,list(data_sha=case$data_sha,case_file_sha=digest::digest(file=file.path(study,'inputs/public/D19.rds'),algo='sha256'),scale=case$known_sigma_growth,scale_sha=sf_hash(case$known_sigma_growth),bank_lock=locks[[origin]]))
  sf_rds(id,file.path(out,'identities',paste0(origin,'.rds')))
 }
 t<-contract$design$bounds;t$scientific_role<-ifelse(t$parameter%in%c('init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N','sigma_growth','sigma_kary'),'nuisance','mechanism');t$estimation_role<-ifelse(t$parameter=='sigma_growth','known_generating_control','unknown_active');t$native_unit_order<-match(t$parameter,f_free(contract$design)$parameter)
 sf_emit(t,file.path(out,'parameter_roles.tsv'));sf_emit(data.frame(parameter='sigma_growth',known_value=case$known_sigma_growth,native_scalar_sha=sf_hash(case$known_sigma_growth),source='original sealed D19 RDS scalar; direct fixed-field injection'),file.path(out,'scale_registry.tsv'))
 status<-system2('gcc',c('-shared','-fPIC','-O2',shQuote(file.path(study,'release/read_guard.c')),'-ldl','-o',shQuote(file.path(out,'read_guard.so'))));check('native_role_guard_compiled',status==0)
}
d_tests<-function(study,out,env,check){
 x<-d_read(study,'D_A');c<-d_context(env,x$case,x$contract,x$id);rng<-readRDS(file.path(study,'verification/T00/rng/isolated_tests.rds'));assign('.Random.seed',rng$state,.GlobalEnv)
 failed<-function(f)isTRUE(tryCatch({f();FALSE},error=function(e)TRUE))
 check('active19_fixed1',sum(c$map$table$status=='active')==19&&sum(c$map$table$status=='fixed')==1&&c$map$table$parameter[c$map$table$status=='fixed']=='sigma_growth')
 for(bound in c(0,1)){
  p<-d_decode(c,rep(bound,19));check(paste0('endpoints_',bound),all(vapply(1:19,function(i)identical(p[[c$free$parameter[i]]],if(bound==0)c$free$lower[i]else c$free$upper[i])||abs(p[[c$free$parameter[i]]]-if(bound==0)c$free$lower[i]else c$free$upper[i])<=64*.Machine$double.eps*max(1,abs(c$free$upper[i])),logical(1))))
 }
 zero<-d_decode(c,rep(0,19));check('exact_mechanism_nulls',identical(zero$p_misseg,0)&&identical(zero$buffer_beta,0)&&identical(zero$buffer_smax,0));check('inactive_labels_preserved',all(c('k_o_mis','buffer_beta','buffer_n_exp')%in%r_activity(zero)))
 check('wrong_dimension_rejected',failed(function()d_decode(c,rep(.5,18))))
 # Nonzero-error fixture is required; D=0 does not test the normalization.
 synthetic_pred<-list(growth=c$case$data$growth+seq_len(114)*1e-5,kary=c$case$data$prob);p<-d_decode(c,rep(.5,19));a<-d_loss(c,synthetic_pred,p)
 c2<-list2env(as.list(c),parent=emptyenv());c2$case$known_sigma_growth<-c$case$known_sigma_growth*2;c2$id$scale_sha<-sf_hash(c2$case$known_sigma_growth);p2<-p;p2$sigma_growth<-c2$case$known_sigma_growth;b<-d_loss(c2,synthetic_pred,p2)
 check('nonzero_scale_DG_quarter',a$growth>0&&abs(b$growth-a$growth/4)<=1e-14);check('nonzero_scale_errors_half',abs(b$rms_growth-a$rms_growth/2)<=1e-14&&abs(b$max_growth-a$max_growth/2)<=1e-14);check('scale_DK_unchanged',identical(a$kary,b$kary))
 p2<-p;p2$sigma_growth<-c$case$known_sigma_growth*.5;check('wrong_postdecode_scale_rejected',failed(function()d_loss(c,synthetic_pred,p2)))
 # All bucket caps, identity fields, interrupted/duplicate charging and RNG integrity.
 fixture<-file.path(out,'counter_fixture');dir.create(fixture)
 dummy<-function(u)list(feasible=FALSE)
 ev<-d_evaluator(c,fixture,'D_A','import',backend=dummy);ev$fn(rep(.5,19));ev$fn(rep(.5,19));s<-ev$snapshot();ev$close();check('duplicate_attempts_charged',s$counts['import']==2L)
 for(field in c('identity','scale_sha','rng_sha','origin')){bad<-s;bad[[field]]<-'stale';check(paste0('stale_',field,'_rejected'),failed(function()d_validate_state(c,bad,'D_A')))}
 bad<-s;bad$rng[2]<-bad$rng[2]+1L;check('stale_rng_payload_rejected',failed(function()d_validate_state(c,bad,'D_A')))
 for(bucket in names(d_limits)){
  q<-s;q$counts[bucket]<-d_limits[bucket];nested<-file.path(out,paste0('cap_',bucket));dir.create(nested);e<-d_evaluator(c,nested,'D_A','import',q,backend=dummy);n<-e$state$n
  check(paste0('cap_',bucket),tryCatch({e$fn(rep(.4,19),bucket);FALSE},recovery_budget_stop=function(e)TRUE));check(paste0('cap_not_extra_attempt_',bucket),e$state$n==n);e$close()
 }
 nested<-file.path(out,'resume');dir.create(nested);ev<-d_evaluator(c,nested,'D_A','import',s,backend=dummy);ev$fn(rep(.4,19));check('restored_counts_preserved',ev$state$counts['import']==3L);ev$close()
 # Native API fixtures deliberately visit actual +Inf regions.
 suppressPackageStartupMessages(library(DEoptim));suppressPackageStartupMessages(library(nloptr))
 check('pinned_optimizers',as.character(packageVersion('DEoptim'))=='2.2.8'&&as.character(packageVersion('nloptr'))=='2.2.1')
 pop<-matrix(runif(3800),200,19);pop[1,1]<-.75;seen<-list();calls<-0L;infs<-0L
 fn<-function(u){calls<<-calls+1L;if(calls<=200L)seen[[calls]]<<-as.numeric(u);if(u[1]>.5){infs<<-infs+1L;return(Inf)};sum((u-.2)^2)}
 de<-DEoptim::DEoptim(fn,rep(0,19),rep(1,19),control=list(NP=200L,itermax=1L,initialpop=pop,trace=FALSE,strategy=2L,F=.8,CR=.9,reltol=0,steptol=2L,VTR=-Inf,parallelType=0))
 check('native_population_exact_order',identical(unname(do.call(rbind,seen)),unname(pop)));check('native400_callbacks_with_Inf',calls==400&&de$optim$nfeval==400&&infs>0&&is.finite(de$optim$bestval));sf_rds(list(pop=de$member$pop,rng=.Random.seed),file.path(out,'toy_native_checkpoint.rds'))
 ni<-0L;nf<-0L;nf_fn<-function(u){nf<<-nf+1L;if(u[1]>=.5){ni<<-ni+1L;Inf}else (u[1]-.9)^2}
 nl<-nloptr::nloptr(.49,nf_fn,lb=0,ub=1,opts=list(algorithm='NLOPT_LN_SBPLX',maxeval=100L,maxtime=0,xtol_abs=1e-5,ftol_abs=1e-8,print_level=0))
 check('native_Subplex_actual_Inf_visits',ni>0&&nf>0&&is.finite(nl$objective))
 boom<-function(u)stop('UNRECOGNIZED_FIXTURE_ERROR')
 check('native_DE_unknown_error_propagates',failed(function()DEoptim::DEoptim(boom,rep(0,19),rep(1,19),control=list(NP=200L,itermax=1L,initialpop=pop,trace=FALSE))))
 check('native_Subplex_unknown_error_propagates',failed(function()nloptr::nloptr(.2,boom,lb=0,ub=1,opts=list(algorithm='NLOPT_LN_SBPLX',maxeval=10L,maxtime=0))))
 # Original passage rules, including day0 exclusion and earliest tie.
 sim<-function(v)list(Ntot_live_obs=v,live_state_obs=cbind(v*.4,v*.6));choose<-function(v,required,target)c$extract(sim(v),required,c(44,88),target,seq_along(v)-1)
 z<-choose(c(200,90,80),100,90);check('day0_cannot_rescue',!z$passage_executed&&z$available_cells==90)
 z<-choose(c(50,110,110),100,110);check('earliest_tie',z$selected_index==2)
 z<-choose(c(100,90,130,115),110,90);check('smallest_surplus_fallback',z$selected_index==4)
 check('finite_negative_growth_permitted',env$ivt_log_growth_rate(100,50,2)<0)
 # Full physical checks: truth is confined here; preselected bank row1 is not ranked.
 log<-file(file.path(out,'forward_attempts.tsv'),'a');writeLines('attempt\tevent\tlabel\tseconds\tfeasible',log);flush(log);n<-0L;times<-list();on.exit(close(log),add=TRUE)
 forward<-function(label,f){if(n>=64L)stop('Engineering external subcap');n<<-n+1L;writeLines(paste(n,'START',label,NA,NA,sep='\t'),log);flush(log);t<-proc.time()[3];z<-f();dt<-proc.time()[3]-t;times[[n]]<<-data.frame(label=label,seconds=dt,feasible=z$feasible);writeLines(paste(n,'END',label,dt,z$feasible,sep='\t'),log);flush(log);z}
 truth<-readRDS(file.path(study,'inputs/private/truth.rds'));reference<-readRDS(file.path(study,'inputs/private/reference.rds'))
 z<-forward('isolated_truth_positive',function()d_physics(c,truth));check('positive_control_feasible',z$feasible);check('positive_D_tolerance',z$loss$value<=1e-12);f_compare(z$pred,reference$pred,check,'sealed_reference')
 for(letter in LETTERS[1:4]){
  bank<-d_bank(study,paste0('D_',letter));m<-bank$members[[1]];p<-d_decode(c,m$unit)
  check(paste0(letter,'_all19_natural_bit_identity'),identical(p[c$free$parameter],m$params[c$free$parameter]));check(paste0(letter,'_placeholder_replaced'),identical(p$sigma_growth,c$case$known_sigma_growth))
  tampered<-m;tampered$params$sigma_growth<-NaN;check(paste0(letter,'_placeholder_sentinel_unused'),identical(d_decode(c,tampered$unit),p))
  z<-forward(paste0(letter,'_observed'),function()d_physics(c,p));check(paste0(letter,'_feasible'),z$feasible);f_compare(z$pred,m$pred,check,paste0(letter,'_scale_injection_physical'))
  r<-forward(paste0(letter,'_independent'),function()d_physics(c,p,TRUE));f_compare(z$pred,r$pred,check,paste0(letter,'_unmodified'))
  v<-d_independent_loss(c,r$pred);check(paste0(letter,'_independent_components'),max(abs(v-c(growth=z$loss$growth,kary=z$loss$kary,value=z$loss$value)))<=1e-8)
 }
 infeasible<-forward('low_growth_protocol_context',function()d_physics(c,d_decode(c,rep(0,19))))
 check('recognized_infeasible_context',!infeasible$feasible&&all(vapply(infeasible$lineages,function(l)l$feasible||(!is.null(l$error$parent)&&!is.null(l$error$child)&&l$error$selection$available_cells<l$error$selection$required_cells),logical(1))))
 sf_emit(f_ledger(infeasible,'predeclared_low_box_fixture',c$contract),file.path(out,'infeasibility_edges.tsv'))
 sf_emit(do.call(rbind,times),file.path(out,'complete_cost_samples.tsv'));sf_emit(data.frame(full_forward_attempts=n,toy_attempts_label='analytic only; separately counted'),file.path(out,'attempt_receipt.tsv'))
 pdf(file.path(out,'headless_fixture.pdf'));plot(0,0,xlab='Engineering fixture',ylab='Value');dev.off()
}
d_canaries<-function(study,out,env,check,mode){
 readbyte<-function(p)tryCatch({readBin(p,'raw',1);TRUE},error=function(e)FALSE)
 contract<-readRDS(file.path(study,'inputs/public/contract.rds'));check('own_public_contract',length(f_free(contract$design)$parameter)==19)
 paths<-c('inputs/private/truth.rds','inputs/private/confirmed_plan.md','inputs/private/private_fixture.txt','inputs/private/engineering_fixture.txt','inputs/public/alias_private.txt','inputs/public/alias_foreign.rds','verification/D_B_fixture/state.rds',paste0('inputs/private/',c('fit_data.Rds','jobs_2N.Rds','jobs_4N.Rds')))
 if(mode=='D'){
  check('D_allowed_case',readbyte(file.path(study,'inputs/public/D19.rds')));check('D_allowed_bank',readbyte(file.path(study,'inputs/banks/D_A.rds')))
  paths<-c(paths,'inputs/banks/D_B.rds','verification/T00/identities/D_B.rds','verification/T00/rng/D_B.rds')
 }else paths<-c(paths,'inputs/public/D19.rds','inputs/public/alias_target.rds','verification/T00/scale_registry.tsv','verification/T00/identities/D_A.rds','inputs/banks/D_A.rds')
 for(i in seq_along(paths)){p<-file.path(study,paths[i]);check(paste0('existing_canary_',i),file.exists(p));check(paste0('denied_canary_',i),!readbyte(p))}
 # Historical payload exists and is outside the study's canonical root.
 old<-file.path(dirname(study),'20261117_culture_D19_feasibility_01/inputs/private/truth.rds');check('historical_canary_exists',file.exists(old));check('historical_canary_denied',!readbyte(old))
 sf_emit(data.frame(path=paths,classification='EXISTING_FIXTURE_DENIED',mode=mode),file.path(out,'negative_read_canaries.tsv'))
}
