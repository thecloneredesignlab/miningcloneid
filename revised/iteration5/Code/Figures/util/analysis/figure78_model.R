# Shared current model adapter and canonical stochastic-passage protocol.
fp_engine <- function(a) {
  m <- fp_model_manifest(a)
  relative <- substring(m$path,nchar(a$model_root)+2L)
  runtime <- file.path(a$run,"prepared/model_source")
  for(i in seq_along(relative)) {
    dest <- file.path(runtime,relative[i]);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if(!file.exists(dest)&&!file.copy(m$path[i],dest,overwrite=FALSE))stop("Could not stage canonical model source")
    if(!identical(fp_sha(dest),m$sha256[i]))stop("Current model staging mismatch: ",dest)
  }
  Sys.setenv(O2SIMPS_RCPP_CACHE_DIR=file.path(a$run,".cache/model"))
  source(file.path(runtime,"simulation/o2/fixed_o2/run_fixed_o2_simulation.R"),local=globalenv(),chdir=TRUE)
  env <- new.env(parent=baseenv())
  sys.source(file.path(runtime,"util/o2_supply_demand_map_invitro_lineage_simulation_utils.R"),envir=env)
  fun <- get("ivt_gaussian_initial_state",envir=env,inherits=FALSE);environment(fun) <- baseenv()
  assign("fp_gaussian_initial_state",fun,globalenv())
  kernels <- c("figure78_culture.cpp",if(a$figure==7L)paste0("figure7_",c("continuous","history","direct"),".cpp"))
  for(kernel in kernels)Rcpp::sourceCpp(file.path(a$code,"util/analysis",kernel),
    cacheDir=file.path(a$run,".cache/propagators"),showOutput=FALSE)
  invisible(runtime)
}
fp_ledger <- function(pr,o2,p=NA_real_) {
  rp <- pr$run_params;native <- pr$native_parameters
  if(length(p)!=1L)stop("One native p_misseg value is required")
  if(is.finite(p)) {
    if(p<.005-1e-12||p>.5+1e-12)stop("Native p_misseg outside joint-fit bounds")
    rp$p_misseg <- min(max(p,.005),.5)
  }
  rp$o2_S0 <- o2
  after <- vapply(names(native),function(n)as.numeric(rp[[n]]),0.0)
  allow <- names(native)%in%c("o2_S0",if(is.finite(p))"p_misseg")
  if(any(abs(after[!allow]-native[!allow])>1e-12))stop("Scenario changed an unapproved native parameter")
  cfg <- pr$config
  fx <- fixo2_fixed_matrix(globalenv(),cfg,rp,O2=o2)
  ng <- fx$ngrid;n <- length(ng)
  lam <- as.numeric(.lambda_eff_of_O2(o2,rp,N=ng));mu <- as.numeric(.mu_eff_of_O2(o2,rp,N=ng))
  pw <- as.numeric(rp$p_wgd %||% 0);W <- matrix(0,n,n)
  dst <- match(ng*2,ng);good <- !is.na(dst);W[cbind(dst[good],which(good))] <- 1
  B <- sweep(as.matrix(fx$G)+diag(lam,n)-sweep(W,2,lam*pw,"*"),2,lam*(1-pw),"/")
  zero <- lam*(1-pw)==0;if(any(zero))B[,zero] <- 0
  B[B<0&B> -1e-12] <- 0
  if(any(!is.finite(B))||min(B)< -1e-10||any(!is.finite(lam))||any(mu<0))stop("Invalid canonical production operator")
  list(M=as.matrix(fx$M),G=as.matrix(fx$G),B=B,W=W,lambda=lam,mu=mu,p_wgd=pw,
    ngrid=ng,N_unit=cfg$N_UNIT %||% 22L,p_misseg_used=rp$p_misseg,
    p_mis_effective=as.numeric(.pmisseg_of_O2(o2,rp,N=ng)))
}
fp_redirect <- function(l,threshold,beta) {
  if(!is.finite(beta)||beta<0||beta>1)stop("beta must be in [0,1]")
  low <- l$ngrid<threshold
  if(!any(low)||all(low))stop("Threshold must split the chromosome grid")
  old <- l$B
  for(j in which(!low)) {
    amount <- beta*sum(old[low,j]);l$B[low,j] <- (1-beta)*old[low,j]
    l$B[j,j] <- l$B[j,j]+amount
  }
  J <- sweep(l$B,2,l$lambda*(1-l$p_wgd),"*")+sweep(l$W,2,l$lambda*l$p_wgd,"*")
  l$G <- J-diag(l$lambda);l$M <- l$G-diag(l$mu)
  error <- max(abs(colSums(l$B)-colSums(old)))
  if(error>1e-12||min(l$B)< -1e-14)stop("Blocking changed viable daughter yield")
  l$intervention_qc <- data.frame(beta=beta,threshold_N=threshold,column_yield_error=error,passed=TRUE)
  l
}
fp_pure_initial <- function(l,initial) {
  index <- match(initial*l$N_unit,l$ngrid)
  if(is.na(index))stop("Pure initial cohort outside chromosome grid")
  x <- numeric(length(l$ngrid));x[index] <- 1;x
}
fp_initial <- function(pr,l,Q,fitted=TRUE) {
  pure <- function(k) {
    if(pr$model_context=="in vitro"&&fitted) {
      x <- fp_gaussian_initial_state(l$ngrid,pr$run_params[[paste0("init_mean_",k,"N")]],
        pr$run_params[[paste0("init_sd_",k,"N")]])
      if(is.null(x))stop("Invalid seed-specific fitted Gaussian initial distribution")
      x
    } else fp_pure_initial(l,k)
  }
  if(Q==0)return(pure(2))
  if(is.infinite(Q))return(pure(4))
  if(!is.finite(Q)||Q<0)stop("Invalid Q = 4N:2N")
  (pure(2)+Q*pure(4))/(1+Q)
}
fp_streams <- function(family,seed,o2,p,Q,n) {
  old_kind <- RNGkind();had <- exists(".Random.seed",globalenv(),inherits=FALSE)
  old_seed <- if(had)get(".Random.seed",globalenv()) else NULL
  on.exit({do.call(RNGkind,as.list(old_kind));if(had)assign(".Random.seed",old_seed,globalenv()) else
    if(exists(".Random.seed",globalenv(),inherits=FALSE))rm(".Random.seed",envir=globalenv())},add=TRUE)
  RNGkind("L'Ecuyer-CMRG","Inversion","Rejection")
  # Preserve the established stochastic-stream contract. This string is an RNG
  # domain separator, not a file/run version or a dependency on an old result.
  set.seed(digest::digest2int(paste("figure78_culture_repair_v1",family,seed,sprintf("%.6f",o2),
    sprintf("%.12f",p),Q,sep="|"),seed=20261007L))
  st <- .Random.seed;state <- label <- matrix(NA_integer_,7L,n)
  for(r in seq_len(n)) {
    state[,r] <- st;label[,r] <- parallel::nextRNGSubStream(st)
    st <- parallel::nextRNGSubStream(label[,r])
  }
  list(state=state,label=label)
}
fp_birth_coefficients <- function(l,threshold=66L) {
  low <- l$ngrid<threshold
  ordinary <- sweep(l$B,2,l$lambda*(1-l$p_wgd),"*")
  J <- ordinary+sweep(l$W,2,l$lambda*l$p_wgd,"*")
  KH <- ordinary;KH[!low,] <- 0;KH[,low] <- 0
  list(J=J,KH=KH,low=low,hl=colSums(J[low,,drop=FALSE])*(!low),
    ll=colSums(J[low,,drop=FALSE])*low)
}
fp_culture <- function(b,pr,l,x0,family,seed,o2,p,Q,threshold=66L,
    tagged=FALSE,daily=TRUE,events=FALSE,days=b$config$main_day %||% b$config$day,
    repeats=b$replicates,passage=TRUE) {
  c <- fp_birth_coefficients(l,threshold)
  bound <- max(pmax(colSums(l$G),0)-l$mu)
  no_cross <- passage&&days*bound<log(b$rule$target_live_cells[1]/b$rule$seed_cells[1])
  nr <- if(no_cross||!passage)1L else as.integer(repeats)
  streams <- fp_streams(family,seed,o2,p,Q,nr)
  o <- f78_culture_cpp(l$G,l$mu,c$KH,x0,l$ngrid/l$N_unit,c$hl,c$ll,c$low,
    as.integer(days),pr$config$DT,isTRUE(pr$config$Crowding),
    if(identical(pr$config$crowding,"logistic"))0L else 1L,pr$config$K,
    b$rule$seed_cells[1],b$rule$target_live_cells[1],streams$state,streams$label,
    passage,tagged,daily,events)
  if(no_cross&&any(o$passage_count!=0))stop("No-passage bound violated")
  o$deterministic_no_passage <- no_cross;o$executed_repeats <- nr;o$nominal_repeats <- as.integer(repeats)
  o
}
