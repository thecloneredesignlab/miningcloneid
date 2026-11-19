n_limits <- c(global=73400L,local=6000L,pattern=2000L,verify=1800L,handoff=800L)
n_origins <- paste0('N_',LETTERS[1:6])
n_payload <- '070c04be23adc7329f840ce83e8a2f0aae229e0bb8e32452c04c5a84444c1978'
n_data <- function(study,origin) {
  case <- readRDS(file.path(study,'inputs/public/N20.rds'))
  contract <- readRDS(file.path(study,'inputs/public/contract.rds'))
  id <- readRDS(file.path(study,'verification/T00/identities',paste0(origin,'.rds')))
  stopifnot(case$kind=='N',is.null(case$known_sigma_growth),
    identical(case$design,contract$design),identical(case$data_sha,n_payload),
    identical(case$data_sha,sf_hash(case$data)),identical(id$data_sha,case$data_sha),
    identical(id$contract_sha,sf_hash(contract)),id$origin==origin)
  list(case=case,contract=contract,id=id)
}
n_context <- function(env,case,contract,id) {
  c <- f_context(env,contract)
  c$case <- case
  c$id <- id
  c$free <- case$design$bounds
  stopifnot(nrow(c$free)==20L,identical(sort(c$free$parameter),sort(cg_roles)),
    identical(names(case$data$growth),as.character(case$design$mask$passage_id[case$design$mask$growth])),
    identical(names(case$data$kary),as.character(case$design$mask$passage_id[case$design$mask$kary])),
    length(case$data$growth)==114L,length(case$data$kary)==12L,
    all(lengths(case$data$kary)==case$design$mask$n_kary[case$design$mask$kary]),
    sum(lengths(case$data$kary))==220L,all(unlist(case$data$kary)%in%22:154),
    identical(unname(case$design$weights),c(1,1,0,0,0)))
  c$low <- env$o2prof_transform(c$free$lower,c$free)
  c$high <- env$o2prof_transform(c$free$upper,c$free)
  c$map <- list(table=transform(c$free,status='active',fixed_value=NA_real_),
                base_params=c(case$design$fixed,cg_neutral(env,c$free)))
  c$identity <- sf_hash(id)
  c$hist <- lapply(case$data$kary,function(x)tabulate(match(x,22:154),nbins=133L))
  c
}
n_decode <- function(c,u) {
  stopifnot(is.numeric(u),length(u)==20L,all(is.finite(u)),all(u>=0 & u<=1))
  free <- setNames(c$low+as.numeric(u)*(c$high-c$low),c$free$coordinate)
  p <- c$env$o2prof_from_free(c$map,free)
  stopifnot(is.finite(p$sigma_growth),p$sigma_growth>0)
  cg_candidate(c$d,p[cg_roles])
}
n_loss <- function(c,pred,p) {
  stopifnot(identical(names(pred$growth),names(c$case$data$growth)),identical(names(pred$kary),names(c$hist)),
            all(lengths(pred$kary)==133L))
  e <- c$case$data$growth-pred$growth
  growth <- mean(-log(p$sigma_growth)-log(2*pi)/2-(e/p$sigma_growth)^2/2)
  kary <- mean(vapply(names(c$hist),function(id) {
    q <- pred$kary[[id]]
    stopifnot(all(is.finite(q)),all(q>=0),abs(sum(q)-1)<=1e-12)
    h <- c$hist[[id]]
    sum(h*log(pmax(q,1e-12)))/sum(h)
  },numeric(1)))
  stopifnot(is.finite(growth+kary))
  list(value=-(growth+kary),growth=growth,kary=kary)
}
n_independent_loss <- function(c,pred,p) {
  g <- mean(dnorm(c$case$data$growth,pred$growth,p$sigma_growth,log=TRUE))
  k <- mean(vapply(names(c$case$data$kary),function(id){
    bins <- match(c$case$data$kary[[id]],22:154)
    mean(log(pmax(pred$kary[[id]][bins],1e-12)))
  },numeric(1)))
  c(growth=g,kary=k,value=-(g+k))
}
# Production propagation follows lineage order; the observed case follows the frozen mask order.
n_align_predictions <- function(case,pred) {
  growth <- names(case$data$growth)
  kary <- names(case$data$kary)
  order_matches <- function(x,ids) {
    length(x)==length(ids) && !anyDuplicated(names(x)) && !anyDuplicated(ids) && setequal(names(x),ids)
  }
  stopifnot(order_matches(pred$growth,growth),order_matches(pred$kary,kary),order_matches(pred$latent_kary,kary))
  pred$growth <- pred$growth[growth]
  pred$kary <- pred$kary[kary]
  pred$latent_kary <- pred$latent_kary[kary]
  pred
}
n_physics <- function(c,p,reference=FALSE) {
  z <- if(reference)f_reference(c,p)else f_physics(c,p)
  if(z$feasible){
    z$pred <- n_align_predictions(c$case,z$pred)
    z$loss <- n_loss(c,z$pred,p)
  }
  z
}
n_sigma_check <- function(c,anchor) {
  j <- match('sigma_growth',c$free$parameter)
  raw <- sqrt(mean((c$case$data$growth-anchor$pred$growth)^2))
  s <- min(c$free$upper[j],max(c$free$lower[j],raw))
  u <- anchor$unit
  u[j] <- min(1,max(0,(c$env$o2prof_transform(s,c$free[j,,drop=FALSE])-c$low[j])/(c$high[j]-c$low[j])))
  p <- n_decode(c,u)
  list(unit=u,unclipped=raw,scale=p$sigma_growth,score=n_independent_loss(c,anchor$pred,p)[['value']])
}
