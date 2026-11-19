# New instrumentation, scientific-contract and API checks before blind work.
f_tests<-function(study,out,env,check){
 contract<-readRDS(file.path(study,'verification/P0/feasibility_contract.rds'));c<-f_context(env,contract)
 failed<-function(f)isTRUE(tryCatch({f();FALSE},error=function(e)TRUE))
 # Analytic trajectories exercise actual original selection functions.
 sim<-function(totals)list(Ntot_live_obs=totals,live_state_obs=cbind(totals*.4,totals*.6))
 choose<-function(totals,required,target,days=0:(length(totals)-1))c$extract(sim(totals),required,c(44,88),target,days)
 a<-choose(c(100,120,150),110,120);check('sufficient_target',a$selected_index==2&&a$passage_executed&&a$boundary_scale<=1)
 a<-choose(c(100,90,130,115),110,90);check('smallest_surplus_fallback',a$selected_index==4&&a$available_cells==115)
 a<-choose(c(200,90,80),100,90);check('day0_does_not_rescue',!a$passage_executed&&a$available_cells==90)
 a<-choose(c(50,100,110),100,100);check('exact_inoculum_equality',a$passage_executed&&a$boundary_scale==1)
 a<-choose(c(50,110,110),100,110);check('earliest_tie',a$selected_index==2)
 a<-choose(c(50,110,120),100,NA_real_);check('missing_target_last_positive',a$selected_index==3)
 check('invalid_nonfinite_simulation',failed(function()f_sim_valid(list(live_state_obs=matrix(NA_real_,2,133),Ntot_live_obs=c(1,2)))))
 check('negative_growth_permitted',env$ivt_log_growth_rate(100,50,2)<0)
 for(bound in c(0,1)){p<-f_decode(env,c$d,rep(bound,19));tab<-f_free(c$d);check(paste0('endpoints_',bound),all(vapply(seq_len(nrow(tab)),function(i)abs(p[[tab$parameter[i]]]-if(bound==0)tab$lower[i]else tab$upper[i])<=64*.Machine$double.eps*max(1,abs(tab$upper[i])),logical(1))))}
 check('wrong_dimension_rejected',failed(function()f_decode(env,c$d,rep(.5,18))))
 # Counts and RNG/checkpoint identity are tested with a deterministic toy backend.
 dummy<-function(u)list(F=sum(u*u),feasible=FALSE)
 ev<-f_evaluator(c,out,'tests','fixture',backend=dummy);ev$fn(rep(.4,19));ev$fn(rep(.3,19));saved<-ev$snapshot();ev$close()
 sf_rds(saved,file.path(out,'counter_fixture.rds'));restored<-readRDS(file.path(out,'counter_fixture.rds'))
 check('counter_roundtrip',restored$count==2&&identical(saved$contract_sha,sf_hash(contract)))
 nested<-file.path(out,'counter_resume');dir.create(nested)
 ev2<-f_evaluator(c,nested,'tests','fixture',restored,backend=dummy);ev2$fn(rep(.2,19));check('counts_retained_on_resume',ev2$state$counts==3);ev2$close()
 validate_checkpoint<-function(x){v<-f_evaluator(c,nested,'tests','fixture',x,backend=dummy);v$close()}
 wrong<-restored;wrong$contract_sha<-'stale';check('stale_contract_rejected',failed(function()validate_checkpoint(wrong)))
 wrong<-restored;wrong$origin<-'foreign';check('foreign_origin_rejected',failed(function()validate_checkpoint(wrong)))
 # Toy-only DE API checks do not consume biological forward attempts.
 suppressPackageStartupMessages(library(DEoptim));check('pinned_DEoptim',as.character(packageVersion('DEoptim'))=='2.2.8')
 set.seed(2517);pop<-matrix(runif(200*19),200,19);seen<-list();n<-0L
 toy<-function(u){n<<-n+1L;if(n<=200)seen[[n]]<<-u;sum(u*u)}
 de<-DEoptim::DEoptim(toy,rep(0,19),rep(1,19),control=list(NP=200L,itermax=1L,initialpop=pop,trace=FALSE,strategy=2L,F=.8,CR=.9,reltol=0,steptol=2L,VTR=-Inf,parallelType=0))
 check('native_population_order_and_count',identical(unname(do.call(rbind,seen)),unname(pop))&&n==400&&de$optim$nfeval==400)
 checkpoint<-list(pop=de$member$pop,rng=.Random.seed);sf_rds(checkpoint,file.path(out,'toy_native_checkpoint.rds'));check('native_pop_checkpoint',identical(readRDS(file.path(out,'toy_native_checkpoint.rds')),checkpoint))
 bank_exception<-function(u)stop(structure(list(message='bank target fixture',call=NULL),class=c('feasibility_bank_complete','error','condition')))
 check('native_bank_target_class_retained',tryCatch({DEoptim::DEoptim(bank_exception,rep(0,19),rep(1,19),control=list(NP=200L,itermax=1L,initialpop=pop,trace=FALSE));FALSE},feasibility_bank_complete=function(e)TRUE))
 # All real forwards below remain isolated positive/reference diagnostics.
 receipt<-file(file.path(out,'forward_attempts.tsv'),'a');writeLines('attempt\tevent\tlabel\tseconds\tfeasible',receipt);flush(receipt);count<-0L;times<-list()
 on.exit(close(receipt),add=TRUE)
 forward<-function(label,f){
  if(count>=125L)f_fail('Tests plus three toy counter attempts reached128');count<<-count+1L
  writeLines(paste(count,'START',label,NA,NA,sep='\t'),receipt);flush(receipt);t<-proc.time()[3];z<-f();secs<-proc.time()[3]-t
  times[[count]]<<-data.frame(label=label,seconds=secs,feasible=z$feasible);writeLines(paste(count,'END',label,sprintf('%.9g',secs),z$feasible,sep='\t'),receipt);flush(receipt);z
 }
 p<-readRDS(file.path(study,'inputs/private/truth.rds'));ref<-readRDS(file.path(study,'inputs/private/reference.rds'))
 original<-forward('original_positive',function()f_reference(c,p));check('original_positive_feasible',original$feasible);f_compare(original$pred,ref$pred,check,'sealed_production')
 observed<-forward('observed_positive',function()f_physics(c,p));check('observed_positive_feasible',observed$feasible);f_compare(observed$pred,original$pred,check,'observer_equivalence')
 for(name in c('sigma_growth','sigma_kary'))for(side in c('lower','upper')){
  q<-p;tab<-c$d$bounds;q[[name]]<-tab[[side]][match(name,tab$parameter)];z<-forward(paste(name,side,sep='_'),function()f_physics(c,q))
  check(paste0(name,'_',side,'_physical_feasibility'),z$feasible)
  check(paste0(name,'_',side,'_physical_states'),identical(z$pred$states,observed$pred$states)&&identical(z$pred$selection,observed$pred$selection))
 }
 for(i in 1:2){z<-forward(paste0('complete_cost_',i),function()f_physics(c,p));check(paste0('cost_positive_',i),z$feasible)}
 # Recognized errors preserve both cohorts and agree with original order.
 units<-readRDS(file.path(study,'verification/P0/old_units.rds'));q<-f_decode(env,c$d,units[1,]);z<-forward('old_first_observed',function()f_physics(c,q,capture=TRUE));orig<-forward('old_first_original',function()f_reference(c,q))
 first<-z$lineages[[which(!vapply(z$lineages,`[[`,logical(1),'feasible'))[1]]]$error
 check('recognized_first_failure_fields',!z$feasible&&!orig$feasible&&identical(first,orig$error))
 check('both_cohorts_observed',setequal(names(z$lineages),c('2N','4N')))
 # Algebraic merit and declared boundary-activity tests.
 pp<-p;pp$p_misseg<-0;qq<-pp;qq$k_o_mis<-c$d$bounds$upper[match('k_o_mis',c$d$bounds$parameter)];check('zero_induction_k_inactive',identical(env$.pmisseg_of_O2(1,pp,N=22:154),env$.pmisseg_of_O2(1,qq,N=22:154)))
 for(beta in c(0,.4)){
  a<-env$cpp_o2simps_pr_delta_vec(44L,.02,1e-8,buffer_smax=0,buffer_beta=beta,buffer_n_exp=.2,N_unit=22L);b<-env$cpp_o2simps_pr_delta_vec(44L,.02,1e-8,buffer_smax=0,buffer_beta=9,buffer_n_exp=9,N_unit=22L);check(paste0('smax_zero_',beta),identical(a$prob,b$prob))
 }
 a<-env$cpp_o2simps_pr_delta_vec(88L,.02,1e-8,buffer_smax=.8,buffer_beta=0,buffer_n_exp=.2,N_unit=22L);b<-env$cpp_o2simps_pr_delta_vec(88L,.02,1e-8,buffer_smax=.8,buffer_beta=0,buffer_n_exp=9,N_unit=22L);check('beta_zero_shape_inactive',identical(a$prob,b$prob))
 check('balanced_merit_zero_only_full',observed$F==0&&z$F>0&&is.finite(z$F))
 sf_emit(do.call(rbind,times),file.path(out,'complete_cost_samples.tsv'));sf_emit(data.frame(full_forward_attempts=count,counter_toy_attempts=3,all_candidate_attempts=count+3),file.path(out,'attempt_receipt.tsv'))
 pdf(file.path(out,'headless_fixture.pdf'));plot(0,0);dev.off()
 sf_emit(data.frame(package='DEoptim',version=as.character(packageVersion('DEoptim')),license=as.character(packageDescription('DEoptim')$License)),file.path(out,'optimizer_versions.tsv'))
 sf_rds(list(canary='existing engineering output, not a permitted construction input'),file.path(out,'state_canary.rds'))
}
f_canary<-function(study,out,env,check){
 allowed<-readRDS(file.path(study,'verification/P0/feasibility_contract.rds'));bank<-readRDS(file.path(study,'verification/P0/rng/F_A.rds'));check('permitted_design_and_own_rng',length(f_free(allowed$design)$parameter)==19&&bank$origin=='F_A')
 paths<-c(file.path(study,'inputs/private/truth.rds'),file.path(study,'inputs/private/source_bundle.rds'),file.path(study,'inputs/diagnostic/raw_bank.tsv'),file.path(study,'verification/P0/rng/F_B.rds'),file.path(study,'verification/P1/state_canary.rds'))
 for(i in seq_along(paths)){check(paste0('existing_fixture_',i),file.exists(paths[i]));check(paste0('forbidden_read_denied_',i),tryCatch({readBin(paths[i],'raw',n=1);FALSE},error=function(e)TRUE))}
 denied<-readLines(file.path(out,'read_denials.tsv'));check('canaries_only_denied',setequal(unique(denied),paths))
}
