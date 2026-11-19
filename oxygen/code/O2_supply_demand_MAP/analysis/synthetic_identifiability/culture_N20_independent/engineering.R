n_setup <- function(study,out,check) {
  case <- readRDS(file.path(study,'inputs/public/N20.rds'))
  contract <- readRDS(file.path(study,'inputs/public/contract.rds'))
  check('exact_noisy_public_case',identical(names(case),c('kind','design','data','data_sha')) &&
    case$kind=='N' && identical(case$data_sha,n_payload) && identical(sf_hash(case$data),n_payload) &&
    identical(case$design,contract$design))
  roles <- case$design$bounds
  nuisance <- c('init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N','sigma_growth','sigma_kary')
  roles$role <- ifelse(roles$parameter %in% nuisance,'nuisance','mechanism')
  roles$status <- 'active'
  check('all_roles_and_modalities',sum(roles$role=='mechanism')==14L&&sum(roles$role=='nuisance')==6L&&
    length(case$data$growth)==114L&&length(case$data$kary)==12L&&sum(lengths(case$data$kary))==220L)
  registry <- sf_tab(file.path(study,'inputs/start_registry.tsv'))
  source_sha <- digest::digest(file=file.path(study,'source_manifest.tsv'),algo='sha256')
  physical <- list()
  for (origin in n_origins) {
    bank <- n_bank(study,origin)
    seed <- readRDS(file.path(study,'inputs/rng',paste0(origin,'.rds')))
    row <- registry[registry$origin==origin,,drop=FALSE]
    check(paste0(origin,'_admission'),nrow(row)==1L&&bank$origin==origin&&seed$origin==origin&&
      identical(bank$roles,roles$parameter)&&identical(bank$bounds_sha,sf_hash(case$design$bounds))&&
      identical(dim(bank$population),c(200L,20L))&&sf_hash(bank$population)==row$population_sha&&
      row$source_sha==bank$source_sha&&row$optimizer_stream_sha==seed$RNG_sha&&sf_hash(seed$state)==seed$RNG_sha&&
      row$growth_stream_sha!=seed$RNG_sha)
    physical[[origin]] <- bank$population[,!(roles$parameter%in%c('sigma_kary','sigma_growth')),drop=FALSE]
    id <- list(origin=origin,data_sha=n_payload,contract_sha=sf_hash(contract),source_sha=source_sha,
               population_sha=row$population_sha,rng_initial_sha=seed$RNG_sha)
    sf_rds(id,file.path(out,'identities',paste0(origin,'.rds')))
    sf_rds(seed,file.path(out,'rng',paste0(origin,'.rds')))
  }
  check('1200_unchanged_physical_admissibility_rows',!anyDuplicated(as.data.frame(do.call(rbind,physical))))
  sf_emit(roles,file.path(out,'N20_parameter_roles.tsv'))
  sf_emit(data.frame(id=case$design$mask$passage_id,growth=case$design$mask$growth,
                    karyotype=case$design$mask$kary,karyotype_cells=case$design$mask$n_kary),
          file.path(out,'N20_observation_contract.tsv'))
  status <- system2('gcc',c('-shared','-fPIC','-O2','-Wall',
    shQuote(file.path(study,'release/read_guard.c')),'-ldl','-o',shQuote(file.path(out,'read_guard.so'))))
  check('canonical_role_guard_compiled',status==0L)
}
n_canaries <- function(study,out,check) {
  readable <- function(path) {
    f <- tryCatch(suppressWarnings(file(path,'rb')),error=function(e)NULL)
    if(is.null(f))return(FALSE)
    on.exit(close(f))
    length(readBin(f,'raw',1L))==1L
  }
  own <- file.path(study,c('inputs/public/N20.rds','inputs/public/contract.rds','inputs/banks/N_A.rds',
    'verification/T00/identities/N_A.rds','verification/T00/rng/N_A.rds'))
  denied <- c(file.path(study,c('inputs/private/truth.rds','inputs/private/reference.rds','inputs/canaries/alias_truth.rds',
    'inputs/banks/N_B.rds','inputs/canaries/alias_B.rds','verification/T00/identities/N_B.rds',
    'verification/T00/rng/N_B.rds','verification/N_B_fixture/state.rds')),
    file.path(dirname(study),'20261118_culture_D19_independent_01',c('inputs/public/D19.rds','verification/D_B_replay2b/state.rds')),
    file.path(dirname(study),'20261117_culture_N007_precision_01','REPORT.md'))
  # All core canaries are real files. A separate study's report may be absent.
  for(i in seq_along(own)) check(paste0('own_',i),readable(own[i]),own[i])
  for(i in seq_along(denied)) {
    if(i==length(denied)&&!file.exists(denied[i]))next
    check(paste0('deny_',i),file.exists(denied[i])&&!readable(denied[i]),denied[i])
  }
  sf_emit(data.frame(path=c(own,denied),expected=c(rep('ALLOW',length(own)),rep('DENY',length(denied)))),
          file.path(out,'read_canaries.tsv'))
}
n_tests <- function(study,out,env,check) {
  inp <- n_data(study,'N_A')
  ctx <- n_context(env,inp$case,inp$contract,inp$id)
  c <- ctx
  j <- match('sigma_growth',c$free$parameter)
  check('all_twenty_are_free',all(c$map$table$status=='active'))
  p0 <- n_decode(c,rep(0,20))
  check('exact_mechanism_nulls',p0$p_misseg==0&&p0$buffer_smax==0&&p0$buffer_beta==0&&
    all(c('k_o_mis','buffer_beta','buffer_n_exp') %in% r_activity(p0)))
  f <- file(file.path(out,'forward_attempts.tsv'),'w')
  on.exit(close(f),add=TRUE)
  writeLines('event\tlabel',f);flush(f)
  ev <- function(p,label,reference=FALSE) {
    writeLines(paste('START',label,sep='\t'),f);flush(f)
    x <- n_physics(c,p,reference)
    writeLines(paste('END',label,sep='\t'),f);flush(f)
    x
  }
  p <- n_decode(c,as.numeric(n_bank(study,'N_A')$population[1,]))
  a <- ev(p,'blind_first_row')
  check('blind_candidate_feasible',a$feasible)
  prod <- cg_score(env,c$d,lapply(a$lineages,`[[`,'run'),p,a$pred,c$case$data)
  check('production_component_agreement',max(abs(c(prod$growth,prod$kary,prod$objective)-
    c(a$loss$growth,a$loss$kary,a$loss$value)))<=1e-8)
  u <- as.numeric(n_bank(study,'N_A')$population[1,])
  u[j] <- if(u[j]<.5) .75 else .25
  p2 <- n_decode(c,u)
  b <- ev(p2,'independent_unknown_growth_scale',TRUE)
  check('candidate_growth_scale_changes_score',b$feasible && p2$sigma_growth!=p$sigma_growth &&
    abs(b$loss$growth-a$loss$growth)>1e-8 && abs(b$loss$kary-a$loss$kary)<=1e-10 &&
    max(abs(n_independent_loss(c,b$pred,p2)-unlist(b$loss[c('growth','kary','value')])))<=1e-8)
  f_compare(a$pred,b$pred,check,'growth_scale_does_not_propagate')
  v <- n_sigma_check(c,n_candidate(c,as.numeric(n_bank(study,'N_A')$population[1,]),a,'N_A','test'))
  check('variance_stationarity_fixture',v$score<=a$loss$value+1e-8)
  truth <- readRDS(file.path(study,'inputs/private/truth.rds'))
  t <- ev(truth,'isolated_generating_score')
  check('isolated_truth_score_valid',t$feasible&&is.finite(t$loss$value))
  z <- cg_score(env,c$d,lapply(t$lineages,`[[`,'run'),truth,t$pred,c$case$data)
  check('isolated_truth_score_agreement',abs(z$objective-t$loss$value)<=1e-8)
  sf_emit(data.frame(role='isolated generating reference without optimizer credit',Q=t$loss$value,
      growth_log_mean=t$loss$growth,karyotype_log_mean=t$loss$kary),
      file.path(out,'isolated_truth_score.tsv'))
  check('optimizer_package_versions',as.character(packageVersion('DEoptim'))=='2.2.8'&&
    as.character(packageVersion('nloptr'))=='2.2.1')
}
