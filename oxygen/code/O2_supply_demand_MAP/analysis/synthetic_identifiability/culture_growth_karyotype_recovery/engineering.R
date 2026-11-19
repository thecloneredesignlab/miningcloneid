# Bounded source, package/API, algebra and biological-interface checks.
r_engineering<-function(study,task,out,env,check){
 failed<-function(f)isTRUE(tryCatch({f();FALSE},error=function(e)TRUE))
 if(task=='I01'){
  suppressPackageStartupMessages(library(DEoptim));suppressPackageStartupMessages(library(nloptr))
  check('DEoptim_pinned_version',as.character(packageVersion('DEoptim'))=='2.2.8')
  sf_emit(data.frame(package=c('DEoptim','nloptr'),version=vapply(c('DEoptim','nloptr'),function(p)as.character(packageVersion(p)),character(1)),license=vapply(c('DEoptim','nloptr'),function(p)as.character(packageDescription(p)$License),character(1))),file.path(out,'optimizer_versions.tsv'))
  fingerprints<-list()
  for(pkg in c('DEoptim','nloptr'))for(f in list.files(system.file(package=pkg),recursive=TRUE,full.names=TRUE))if(grepl('[.](rdb|rdx|so)$',f))fingerprints[[length(fingerprints)+1L]]<-data.frame(package=pkg,file=f,sha256=digest::digest(file=f,algo='sha256'))
  sf_emit(do.call(rbind,fingerprints),file.path(out,'optimizer_files.tsv'))
  set.seed(1729);pop<-matrix(runif(40),20,2);inputs<-list();n<-0L
  toy<-function(u){n<<-n+1L;if(n<=20L)inputs[[n]]<<-u;if(u[1]<.1)Inf else (u[1]-.4)^2}
  de<-DEoptim::DEoptim(toy,c(0,0),c(1,1),control=list(NP=20L,itermax=2L,trace=FALSE,initialpop=pop,F=.8,CR=.9,strategy=2L,reltol=0,steptol=3L,VTR=-Inf,parallelType=0,storepopfrom=1L,storepopfreq=1L))
  check('native_DE_callback_count',n==60L&&de$optim$nfeval==n)
  actual<-do.call(rbind,inputs);keys<-function(a)sort(apply(a,1,function(v)paste(sprintf('%.17g',v),collapse=',')))
  check('explicit_native_population',identical(keys(actual),keys(pop)))
  check('native_population_checkpoint',is.matrix(de$member$pop)&&identical(dim(de$member$pop),c(20L,2L)))
  sf_rds(list(pop=de$member$pop,rng=.Random.seed),file.path(out,'toy_checkpoint.rds'));saved<-readRDS(file.path(out,'toy_checkpoint.rds'));assign('.Random.seed',saved$rng,.GlobalEnv)
  inputs<-list();n<-0L
  de2<-DEoptim::DEoptim(toy,c(0,0),c(1,1),control=list(NP=20L,itermax=1L,trace=FALSE,initialpop=saved$pop,F=.8,CR=.9,strategy=2L,reltol=0,steptol=2L,VTR=-Inf,parallelType=0))
  check('linked_population_first_callbacks',identical(keys(do.call(rbind,inputs)),keys(saved$pop))&&n==40L)
  nl_calls<-0L;nf<-function(u){nl_calls<<-nl_calls+1L;if(u[1]<.1)Inf else (u[1]-.4)^2}
  nl<-nloptr::nloptr(c(.7,.3),nf,lb=c(0,0),ub=c(1,1),opts=list(algorithm='NLOPT_LN_SBPLX',maxeval=20L,maxtime=0,xtol_abs=rep(1e-5,2),ftol_abs=1e-8,print_level=0))
  check('Subplex_finite_and_bounded',is.finite(nl$objective)&&all(nl$solution>=0&nl$solution<=1))
  sf_emit(data.frame(native_maxeval=20,actual_callbacks=nl_calls,native_status=nl$status,maxtime=0),file.path(out,'Subplex_api.tsv'))
  fn_calls<-0L;limited<-function(x){if(fn_calls>=3)stop(structure(list(message='fixture cap',call=NULL),class=c('recovery_budget_stop','error','condition')));fn_calls<<-fn_calls+1L;sum(x*x)}
  stopped<-tryCatch({nloptr::nloptr(c(.7,.3),limited,lb=c(0,0),ub=c(1,1),opts=list(algorithm='NLOPT_LN_SBPLX',maxeval=20L,maxtime=0));FALSE},recovery_budget_stop=function(e)TRUE)
  check('native_wrapper_cap_caught',stopped&&fn_calls==3L)
  poll<-r_poll(function(x)sum(x*x),c(.004,.4),diag(2),c(.001,.0001),1e-8,3L)
  check('partial_poll_not_complete',poll$budget_limited&&!all(poll$rows$complete))
  pdf(file.path(out,'headless_fixture.pdf'));plot(0,0);dev.off();check('headless_pdf',file.info(file.path(out,'headless_fixture.pdf'))$size>0)
  pub<-readRDS(file.path(study,'inputs/public/evaluation_bundle.rds'));reference<-readRDS(file.path(study,'inputs/private/reference_result.rds'));truth<-readRDS(file.path(study,'inputs/private/truth_parameters.rds'))
  expected<-list(kind='deterministic_expected_predictions_not_finite_sample',growth=reference$pred$growth,prob=lapply(reference$pred$kary,function(q)q/sum(q)))
  D<-list(kind='D',design=pub$design,data=expected,known_sigma_growth=truth$sigma_growth,data_sha=sf_hash(expected));N<-list(kind='N',design=pub$design,data=pub$observations,data_sha=sf_hash(pub$observations))
  check('frozen_noisy_payload',N$data_sha=='070c04be23adc7329f840ce83e8a2f0aae229e0bb8e32452c04c5a84444c1978')
  sf_rds(D,file.path(out,'cases/D19.rds'));sf_rds(N,file.path(out,'cases/N20.rds'))
  RNGkind("L'Ecuyer-CMRG","Inversion","Rejection");set.seed(20261117);state<-.Random.seed;rows<-list()
  for(origin in c(paste0('D_',LETTERS[1:6]),paste0('N_',LETTERS[1:6]))){sf_rds(list(origin=origin,state=state),file.path(out,'banks',paste0(origin,'.rds')));rows[[length(rows)+1L]]<-data.frame(origin=origin,RNG_sha=sf_hash(state),bank='unobserved proposals from registered independent stream');state<-parallel::nextRNGStream(state)}
  sf_emit(do.call(rbind,rows),file.path(out,'origin_registry.tsv'));sf_emit(pub$design$bounds,file.path(out,'public_bounds.tsv'))
  status<-system2('gcc',c('-shared','-fPIC','-O2',shQuote(file.path(study,'release/read_guard.c')),'-ldl','-o',shQuote(file.path(out,'read_guard.so'))));check('native_read_guard_compiled',status==0)
  return(invisible(NULL))
 }
 if(task=='I03'){
  case<-readRDS(file.path(study,'verification/I01/cases/N20.rds'));bank<-readRDS(file.path(study,'verification/I01/banks/N_A.rds'))
  check('public_N_no_known_truth_SD',is.null(case$known_sigma_growth)&&bank$origin=='N_A')
  paths<-c(file.path(study,'inputs/private/truth_parameters.rds'),file.path(study,'verification/I01/cases/D19.rds'),file.path(study,'inputs/private/reference_result.rds'),file.path(study,'verification/I02_R1/state.rds'))
  for(i in seq_along(paths))check(paste0('truth_or_D_read_denied_',i),failed(function()readRDS(paths[i])))
  denied<-readLines(file.path(out,'read_denials.tsv'));check('only_explicit_canary_denials',length(denied)>=3&&all(denied%in%paths))
  c<-r_context(env,case);neutral<-cg_neutral(env,case$design$bounds);u<-r_encode(c,neutral);p<-r_decode(c,u);check('N20_complete_injection',length(c$free$parameter)==20&&all(cg_roles%in%names(p)))
  return(invisible(NULL))
 }
 stopifnot(task=='I02')
 N<-readRDS(file.path(study,'verification/I01/cases/N20.rds'));D<-readRDS(file.path(study,'verification/I01/cases/D19.rds'))
 cn<-r_context(env,N);cd<-r_context(env,D);truth<-readRDS(file.path(study,'inputs/private/truth_parameters.rds'));ref<-readRDS(file.path(study,'inputs/private/reference_result.rds'))
 check('D19_and_N20_roles',nrow(cd$free)==19&&nrow(cn$free)==20&&!'sigma_growth'%in%cd$free$parameter)
 for(c in list(cd,cn))for(bound in c('lower','upper')){u<-rep(if(bound=='lower')0 else 1,nrow(c$free));p<-r_decode(c,u);check(paste0(c$case$kind,'_endpoint_',bound),all(vapply(seq_len(nrow(c$free)),function(i)abs(p[[c$free$parameter[i]]]-c$free[[bound]][i])<=64*.Machine$double.eps*max(1,abs(c$free$upper[i])),logical(1))))}
 check('missing_coordinate_rejected',failed(function()r_decode(cn,rep(.5,19))))
 prior<-NULL
 if(basename(out)=='I02_R1'){
  old<-sf_tab(file.path(study,'verification/I02/evaluations.tsv'));charged<-sum(old$event=='START')
  counters<-setNames(integer(length(r_limits)),names(r_limits));counters['external']<-charged
  prior<-list(counts=counters,best=list(score=Inf))
 }
 ev<-r_evaluator(cn,out,'engineering','external',prior);on.exit(ev$close(),add=TRUE)
 u<-r_encode(cn,truth);q<-ev$fn(u);check('decoded_truth_selection',identical(ev$state$last$pred$selection,ref$pred$selection))
 for(name in names(ref$pred$states))check(paste0('state_',name),max(abs(ev$state$last$pred$states[[name]]-ref$pred$states[[name]])/pmax(1,abs(ref$pred$states[[name]])))<=1e-10)
 dist<-r_loss(cd,ev$state$last$pred,truth);check('D_positive_zero',dist$value<=1e-12)
 full<-cg_score(env,N$design,ev$state$last$runs,r_decode(cn,u),ev$state$last$pred,N$data);check('noisy_fast_score_reference',abs(q-full$objective)<=1e-8)
 # Boundary activity is checked against actual exported source kernels.
 for(n in c(44L,88L)){
  z1<-env$cpp_o2simps_pr_delta_vec(n,.02,eps_tail=1e-8,buffer_smax=0,buffer_beta=.2,buffer_n_exp=.2,N_unit=22L)
  z2<-env$cpp_o2simps_pr_delta_vec(n,.02,eps_tail=1e-8,buffer_smax=0,buffer_beta=9,buffer_n_exp=9,N_unit=22L)
  check(paste0('smax_zero_kernel_',n),identical(z1$prob,z2$prob)&&identical(z1$mass_dropped,z2$mass_dropped))
  z1<-env$cpp_o2simps_pr_delta_vec(n,.02,eps_tail=1e-8,buffer_smax=.8,buffer_beta=0,buffer_n_exp=.2,N_unit=22L)
  z2<-env$cpp_o2simps_pr_delta_vec(n,.02,eps_tail=1e-8,buffer_smax=.8,buffer_beta=0,buffer_n_exp=9,N_unit=22L)
  check(paste0('beta_zero_kernel_',n),identical(z1$prob,z2$prob))
 }
 pp<-truth;pp$buffer_smax<-0;check('smax_classifier',setequal(r_activity(pp),c('buffer_beta','buffer_n_exp')))
 pp<-truth;pp$p_misseg<-0;pp2<-pp;pp2$k_o_mis<-.02
 check('zero_induction_function',identical(env$.pmisseg_of_O2(1,pp,N=22:154),env$.pmisseg_of_O2(1,pp2,N=22:154)))
 # Positive-control perturbations test wiring, not blind-origin feasibility.
 # Failed independent random fixtures remain in I02 and are never injected into
 # a population. Independent cost/feasibility admission still uses D_A's bank.
 feasible<-0L;fixture_names<-c('p_misseg','buffer_beta','init_sd_2N','sigma_kary')
 for(name in fixture_names){z<-u;i<-match(name,cn$free$parameter);z[i]<-z[i]+if(z[i]>.5)-1e-6 else 1e-6
  q<-ev$fn(z);check(paste0('perturbed_fixture_feasible_',name),ev$state$last$feasible)
  p<-r_decode(cn,z);full<-cg_score(env,N$design,ev$state$last$runs,p,ev$state$last$pred,N$data)
  check(paste0('nontruth_fast_full_',name),abs(full$objective-q)<=1e-8);feasible<-feasible+1L}
 check('representative_nontruth_scores',feasible>=4)
 # Source-shape, protocol infeasibility and candidate SD are never swallowed.
 p<-truth;p$sigma_growth<-truth$sigma_growth*2;q1<-r_loss(cn,ref$pred,truth)$value;q2<-r_loss(cn,ref$pred,p)$value;check('candidate_variance_not_fixed',abs(q1-q2)>1e-6)
 saved<-ev$snapshot(list(feasible_checks=feasible));sf_rds(saved,file.path(out,'state.rds'));sf_emit(data.frame(external_attempts=sum(saved$counts),positive_D=dist$value,representative_feasible=feasible),file.path(out,'engineering_receipt.tsv'))
 status<-system2('gcc',c('-shared','-fPIC','-O2',shQuote(file.path(study,'release/read_guard.c')),'-ldl','-o',shQuote(file.path(out,'read_guard.so'))));check('retry_namespace_read_guard_compiled',status==0)
}
