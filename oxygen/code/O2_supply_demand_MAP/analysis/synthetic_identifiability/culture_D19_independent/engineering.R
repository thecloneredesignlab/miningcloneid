blind_setup <- function(study, out, check) {
  case_path <- file.path(study, 'inputs/public/D19.rds')
  case <- readRDS(case_path)
  contract <- readRDS(file.path(study, 'inputs/public/contract.rds'))
  check('sealed_case_and_design', case$kind == 'D' && identical(case$design, contract$design) &&
        identical(case$data_sha, sf_hash(case$data)) && length(case$data$growth) == 114L &&
        length(case$data$prob) == 12L && all(lengths(case$data$prob) == 133L) &&
        is.double(case$known_sigma_growth) && is.finite(case$known_sigma_growth) && case$known_sigma_growth > 0)
  check('sealed_contract_digest', sf_hash(contract) == '3d765c1cd021deff5ebe32a48140200b3786cdcf05eea9c4943afa2c188ab2d4')
  source_sha <- digest::digest(file = file.path(study, 'source_manifest.tsv'), algo = 'sha256')
  registry <- sf_tab(file.path(study, 'inputs/start_registry.tsv'))
  physical <- list()
  for (origin in paste0('D_', c('A','B','D','F'))) {
    bank <- d_bank(study, origin)
    seed <- readRDS(file.path(study, 'inputs/rng', paste0(origin, '.rds')))
    row <- registry[registry$origin == origin, , drop = FALSE]
    check(paste0(origin, '_exact_release'), nrow(row) == 1L &&
          bank$origin == origin && bank$source_sha == row$source_sha &&
          identical(dim(bank$population), c(200L, 19L)) &&
          sf_hash(bank$population) == row$population_sha &&
          all(is.finite(bank$population)) && all(bank$population >= 0 & bank$population <= 1) &&
          bank$contract_sha == sf_hash(contract) &&
          seed$origin == origin && sf_hash(seed$state) == seed$RNG_sha && seed$RNG_sha == row$new_rng_sha &&
          paste(seed$state, collapse = ',') == row$new_rng_state && row$new_rng_sha != row$historical_rng_sha)
    physical[[origin]] <- bank$population[, f_free(contract$design)$parameter != 'sigma_kary', drop = FALSE]
    id <- list(origin = origin, contract_sha = sf_hash(contract), rng_initial_sha = sf_hash(seed$state),
               source_sha = source_sha, data_sha = case$data_sha,
               case_file_sha = digest::digest(file = case_path, algo = 'sha256'),
               scale = case$known_sigma_growth, scale_sha = sf_hash(case$known_sigma_growth),
               bank_lock = list(bank_sha = bank$source_sha, population_sha = sf_hash(bank$population)))
    sf_rds(id, file.path(out, 'identities', paste0(origin, '.rds')))
    sf_rds(seed, file.path(out, 'rng', paste0(origin, '.rds')))
  }
  check('all_800_physical_rows_distinct', !anyDuplicated(as.data.frame(do.call(rbind, physical))))
  status <- system2('gcc', c('-shared', '-fPIC', '-O2', '-Wall',
    shQuote(file.path(study, 'release/read_guard.c')), '-ldl', '-o', shQuote(file.path(out, 'read_guard.so'))))
  check('canonical_role_guard_compiled', status == 0L)
}

blind_canaries <- function(study, out, check) {
  readable <- function(path) {
    connection <- tryCatch(suppressWarnings(file(path, 'rb')), error = function(e) NULL)
    if (is.null(connection)) return(FALSE)
    on.exit(close(connection))
    length(readBin(connection, 'raw', 1L)) == 1L
  }
  allow <- file.path(study, c('inputs/public/contract.rds', 'inputs/public/D19.rds',
    'inputs/banks/D_A.rds', 'verification/T00/identities/D_A.rds', 'verification/T00/rng/D_A.rds'))
  deny <- file.path(study, c('inputs/private/truth.rds', 'inputs/private/reference.rds',
   'inputs/banks/D_B.rds', 'inputs/canaries/alias_truth.rds', 'inputs/canaries/alias_B.rds',
   'verification/D_B_fixture/state.rds', 'verification/T00/identities/D_B.rds',
   'verification/T00/rng/D_B.rds', 'inputs/rng/D_D.rds'))
  history <- file.path(dirname(study), '20261117_culture_D19_recovery_restart_01/verification/D_A_replay2b/state.rds')
  warm <- file.path(dirname(study), '20261118_culture_D19_C_feasibility_pilot_01/verification/V2/point.rds')
  deny <- c(deny, history, warm)
  for (i in seq_along(allow)) check(paste0('allowed_', i), readable(allow[i]), allow[i])
  for (i in seq_along(deny)) check(paste0('denied_', i), file.exists(deny[i]) && !readable(deny[i]), deny[i])
  sf_emit(data.frame(path = c(allow, deny), expected = rep(c('ALLOWED','DENIED'), c(length(allow),length(deny)))),
          file.path(out, 'read_canaries.tsv'))
}
d_tests<-function(study,out,env,check){
 x<-d_read(study,'D_A');c<-d_context(env,x$case,x$contract,x$id);set.seed(3361L)
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
 ev<-d_evaluator(c,fixture,'D_A','import',backend=dummy);ev$fn(rep(.5,19),'handoff');ev$fn(rep(.5,19),'handoff');s<-ev$snapshot();ev$close();check('duplicate_attempts_charged',s$counts['handoff']==2L)
 for(field in c('identity','scale_sha','rng_sha','origin')){bad<-s;bad[[field]]<-'stale';check(paste0('stale_',field,'_rejected'),failed(function()d_validate_state(c,bad,'D_A')))}
 bad<-s;bad$rng[2]<-bad$rng[2]+1L;check('stale_rng_payload_rejected',failed(function()d_validate_state(c,bad,'D_A')))
 for(bucket in names(d_limits)){
  q<-s;q$counts[bucket]<-d_limits[bucket];nested<-file.path(out,paste0('cap_',bucket));dir.create(nested);e<-d_evaluator(c,nested,'D_A','import',q,backend=dummy);n<-e$state$n
  check(paste0('cap_',bucket),tryCatch({e$fn(rep(.4,19),bucket);FALSE},recovery_budget_stop=function(e)TRUE));check(paste0('cap_not_extra_attempt_',bucket),e$state$n==n);e$close()
 }
 nested<-file.path(out,'resume');dir.create(nested);ev<-d_evaluator(c,nested,'D_A','import',s,backend=dummy);ev$fn(rep(.4,19),'handoff');check('restored_counts_preserved',ev$state$counts['handoff']==3L);ev$close()
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
 interface_fixture()
 check('bounded_restoration_interface_fixture',TRUE)
 a<-list(ok=TRUE,metrics=list(pass=TRUE,violation=0),candidate=list(score=8e-5))
 b<-list(ok=TRUE,metrics=list(pass=FALSE,violation=4.6e-12),candidate=list(score=7e-5))
 check('strict_pass_precedes_tiny_positive_violation',is_improved(a,b)&&!is_improved(b,a))
 toy<-matrix(0,1710,19);toy[1,1]<-1;res<-c(.2,rep(0,1709))
 proposal<-joint_step(rep(.5,19),res,toy,.03,1e-4)
 check('bounded_joint_step_fixture',proposal$predicted_gain>0&&proposal$proposal[1]<.5&&max(abs(proposal$proposal-.5))<=.03+1e-12)
 fake<-list(candidate=list(unit=rep(.5,19)),residual=res,cutoff_margin=1,branch='A')
 probe<-function(u,...)list(ok=TRUE,residual=res,cutoff_margin=1,branch='B')
 invalid<-joint_jacobian(fake,probe,function()1000L)
 check('derivative_rejects_wrong_branch',all(invalid$scheme=='blocked')&&all(invalid$matrix==0))
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
 z<-forward('isolated_truth_positive',function()d_physics(c,truth));check('positive_control_feasible',z$feasible);check('positive_scale_exact',identical(truth$sigma_growth,c$case$known_sigma_growth));check('positive_D_tolerance',z$loss$value<=1e-12);f_compare(z$pred,reference$pred,check,'sealed_reference')
 for(letter in c('A','B','D','F')){
  bank<-d_bank(study,paste0('D_',letter));unit<-as.numeric(bank$population[1,]);p<-d_decode(c,unit)
  check(paste0(letter,'_all19_natural_bit_identity'),identical(p[c$free$parameter],bank$natural[[1]]));check(paste0(letter,'_placeholder_replaced'),identical(p$sigma_growth,c$case$known_sigma_growth))
  z<-forward(paste0(letter,'_observed'),function()d_physics(c,p));check(paste0(letter,'_feasible'),z$feasible)
  r<-forward(paste0(letter,'_independent'),function()d_physics(c,p,TRUE));f_compare(z$pred,r$pred,check,paste0(letter,'_unmodified'))
  v<-d_independent_loss(c,r$pred);check(paste0(letter,'_independent_components'),max(abs(v-c(growth=z$loss$growth,kary=z$loss$kary,value=z$loss$value)))<=1e-8)
 }
 infeasible<-forward('low_growth_protocol_context',function()d_physics(c,d_decode(c,rep(0,19))))
 check('recognized_infeasible_context',!infeasible$feasible&&all(vapply(infeasible$lineages,function(l)l$feasible||(!is.null(l$error$parent)&&!is.null(l$error$child)&&l$error$selection$available_cells<l$error$selection$required_cells),logical(1))))
 sf_emit(f_ledger(infeasible,'predeclared_low_box_fixture',c$contract),file.path(out,'infeasibility_edges.tsv'))
 sf_emit(do.call(rbind,times),file.path(out,'complete_cost_samples.tsv'));sf_emit(data.frame(full_forward_attempts=n,toy_attempts_label='analytic only; separately counted'),file.path(out,'attempt_receipt.tsv'))
 pdf(file.path(out,'headless_fixture.pdf'));plot(0,0,xlab='Engineering fixture',ylab='Value');dev.off()
}
