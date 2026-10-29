# Figure7-only native-fit passage replay and explicitly assumed extrapolation.
# No Figure8/shared protocol modification; no input from an old prediction/audit.
.f7p_dir <- dirname(normalizePath(sys.frame(1)$ofile,mustWork=TRUE))
source(file.path(.f7p_dir,"figure78_paths.R"),local=globalenv())
source(file.path(.f7p_dir,"figure78_model.R"),local=globalenv())
f7p_scalar <- function(x,fallback=NA_real_) {
  if(length(x)&&is.finite(as.numeric(x)[1]))as.numeric(x)[1] else fallback
}
f7p_max <- function(x,empty=0)if(any(is.finite(x)))max(x[is.finite(x)]) else empty
f7p_distance <- function(x,y) {
  if(length(x)!=length(y)||any(xor(is.finite(x),is.finite(y))))return(Inf)
  f7p_max(abs(x-y))
}
f7p_write <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  data.table::fwrite(x,path,sep="\t",na="NA")
}
f7p_runtime <- function(workspace,run) {
  suppressPackageStartupMessages({library(data.table);library(Matrix);library(Rcpp);library(digest);library(jsonlite)})
  data.table::setDTthreads(2L)
  locked <- f78_locked_model_manifest(workspace)
  runtime <- file.path(run,"runtime/model_source")
  for(i in seq_len(nrow(locked))) {
    dest <- file.path(runtime,locked$relative_path[i]);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    if(!file.exists(dest))stopifnot(file.copy(locked$path[i],dest))
    stopifnot(identical(digest::digest(dest,file=TRUE,algo="sha256"),locked$sha256[i]))
  }
  Sys.setenv(O2SIMPS_RCPP_CACHE_DIR=file.path(run,"runtime/native_cache"))
  source(file.path(runtime,"simulation/o2/fixed_o2/run_fixed_o2_simulation.R"),local=globalenv(),chdir=TRUE)
  # The native selector is sourced unchanged. Its sole summary dependency below
  # is an exact copy of the upstream seven-line weighted-mean function.
  native <- new.env(parent=globalenv())
  native$ivt_weighted_mean_kary_N <- function(weights,grid_pre) {
    w <- as.numeric(weights)
    if(length(w)!=length(grid_pre))stop("weights and grid_pre length mismatch.")
    s <- sum(w)
    if(!is.finite(s)||s<=0)return(NA_real_)
    sum((w/s)*as.numeric(grid_pre))
  }
  sys.source(file.path(runtime,"util/o2_supply_demand_map_invitro_lineage_simulation_utils.R"),envir=native)
  Rcpp::sourceCpp(file.path(.f7p_dir,"figure7_passaged.cpp"),
    cacheDir=file.path(run,"runtime/instrumentation_cache"),showOutput=FALSE)
  f7p_write(locked,file.path(run,"results/locked_model_sources.tsv"))
  list(native=native,locked=locked)
}
f7p_load_fit <- function(workspace,family,fit_root=NULL) {
  shared <- jsonlite::fromJSON(file.path(workspace,"Code/config/figure78_inputs.json"))
  stopifnot(family%in%c("C01","C02"))
  info <- shared$families[[family]]
  root <- f78_fit_root(workspace,shared,fit_root)
  dir <- file.path(root,info$pair_id,paste0("seed",info$winner_seed))
  files <- file.path(dir,c("fit_result.rds","joint_best_params_long.tsv","fit_config.rds"))
  if(any(!file.exists(files)))stop("Complete original joint fit is required: ",dir)
  fit <- readRDS(files[1]);cfg <- fit$ctx$invitro_cfg
  stopifnot(is.list(cfg),identical(cfg$crowding,"logistic"),!isTRUE(cfg$o2_burden_feedback))
  native <- read_joint_scope_params(files[2],"invitro_effective")
  stopifnot(length(native)==33L,!anyDuplicated(names(native)),all(is.finite(native)))
  rp <- prepare_run_params(native,"invitro",cfg,5)
  fitted <- vapply(names(native),function(n)as.numeric(fit$invitro_run_params[[n]]),0.0)
  after <- vapply(names(native),function(n)as.numeric(rp[[n]]),0.0)
  stopifnot(f7p_distance(native,fitted)<1e-12,
    f7p_distance(after[names(native)!="o2_S0"],native[names(native)!="o2_S0"])<1e-12)
  list(family=family,seed=info$winner_seed,fit=fit,
    prepared=list(config=cfg,original_fit_config=cfg,run_params=rp,
      native_parameters=native,model_context="in vitro"),
    parity=data.table::data.table(family=family,parameter=names(native),fitted=fitted,
      loaded=as.numeric(native),passed=abs(native-fitted)<1e-12),
    manifest=data.table::data.table(path=files,
      sha256=vapply(files,digest::digest,"",file=TRUE,algo="sha256")))
}

f7p_protocol <- function(adapter,block=5L) {
  segments <- adapter$segments;ids <- vapply(segments,`[[`,"","segment_id")
  stopifnot(length(ids)>0L,!anyDuplicated(ids),block>0L)
  parent <- rep(NA_character_,length(ids));inventory <- vector("list",length(ids))
  for(j in seq_along(ids)) {
    s <- segments[[j]];p <- as.character(s$parent_segment_id)[1]
    if(!length(p)||is.na(p)||!nzchar(p)||!p%in%head(ids,j-1L))
      p <- if(j>1L)ids[j-1L] else NA_character_
    parent[j] <- p
    days <- as.numeric(s$obs_days_local);duration <- f7p_scalar(s$duration_days)
    stopifnot(is.finite(duration),duration>=1,identical(days,as.numeric(0:duration)))
    raw_i <- f7p_scalar(s$initial_cells);raw_f <- f7p_scalar(s$final_cells)
    target <- raw_f
    if((!is.finite(target)||target<=0)&&j<length(ids))target <- f7p_scalar(segments[[j+1L]]$initial_cells)
    inventory[[j]] <- data.table::data.table(index=j,segment_id=ids[j],parent=p,
      original_O2=f7p_scalar(s$oxygen_pct),initial_raw=raw_i,final_raw=raw_f,
      target=target,duration=duration,observation_days=paste(days,collapse=","),
      complete=is.finite(raw_i)&&raw_i>0&&is.finite(raw_f)&&raw_f>0)
  }
  inv <- data.table::rbindlist(inventory)
  leaves <- setdiff(ids,stats::na.omit(parent))
  group <- function(s)if(!grepl("::",s,fixed=TRUE))"control" else
    if(grepl("::O1$",s))"O1" else if(grepl("::O2$",s))"O2" else "other"
  trace <- function(terminal) {
    path <- character();id <- terminal
    repeat {
      if(id%in%path)stop("Cycle in native passage DAG")
      path <- c(id,path);id <- parent[match(id,ids)]
      if(is.na(id))break
    };path
  }
  routes <- templates <- list()
  for(g in c("control","O1","O2")) {
    candidates <- leaves[vapply(leaves,group,"")==g]
    stopifnot(length(candidates)>0)
    paths <- lapply(candidates,trace)
    durations <- vapply(paths,function(p)sum(inv$duration[match(p,ids)]),0.0)
    chosen <- order(-durations,candidates)[1]
    path <- paths[[chosen]];eligible <- path[inv$complete[match(path,ids)]]
    if(length(eligible)<block)stop("Not enough actual complete records for ",g)
    routes[[g]] <- path;templates[[g]] <- tail(eligible,block)
  }
  list(adapter=adapter,inventory=inv,routes=routes,templates=templates,block=block,
    short_terminal_ids=setdiff(leaves,vapply(routes,tail,"",1)))
}

f7p_loss <- function(l,pr,o2,pm) {
  rp <- pr$run_params;ng <- l$ngrid
  upper <- lower <- numeric(length(ng))
  for(j in seq_along(ng)) {
    z <- cpp_o2simps_pr_delta_vec(ng[j],l$p_mis_effective[j],1e-8,
      rp$buffer_smax,rp$buffer_beta,rp$buffer_n_exp,l$N_unit)
    dest <- ng[j]+z$ts
    upper[j] <- l$lambda[j]*((1-l$p_wgd)*sum(z$prob[dest>max(ng)])+
      l$p_wgd*(2*ng[j]>max(ng)))
    lower[j] <- l$lambda[j]*(1-l$p_wgd)*sum(z$prob[dest<min(ng)])
  }
  list(upper=upper,lower=lower)
}
f7p_operator <- function(input,o2,pm,top,beta=0,threshold=66L) {
  stopifnot(pm%in%c(.1,.3),top>=4,top==as.integer(top),o2>=0,o2<=20)
  pr <- input$prepared;pr$config$N_MAX <- as.integer(top*pr$config$N_UNIT)
  pr$run_params$p_misseg <- pm;pr$run_params$o2_S0 <- o2
  l <- fp_redirect(fp_ledger(pr,o2,pm),threshold,beta)
  c <- fp_birth_coefficients(l,threshold);loss <- f7p_loss(l,pr,o2,pm)
  sparse <- function(m)as(Matrix::Matrix(m,sparse=TRUE),"dgCMatrix")
  list(pr=pr,l=l,c=c,G=sparse(l$G),KH=sparse(c$KH),GU=sparse(l$G-c$KH),loss=loss,
    meta=data.table::data.table(family=input$family,seed_number=input$seed,
      O2_pct=o2,p_misseg=pm,maximum_ploidy=top,beta=beta,threshold_N=threshold))
}
f7p_sim_segment <- function(op,native,seg,inventory,x,d,u,elapsed,parent_cells=NA_real_,
    ancestry_known=TRUE) {
  cfg <- op$pr$config;l <- op$l;c <- op$c
  n0 <- f7p_scalar(seg$initial_cells,cfg$init_total_size)
  if(n0<=0)n0 <- cfg$init_total_size
  sim <- f7p_segment_cpp(op$G,op$KH,op$GU,l$lambda,l$mu,c$hl,c$ll,
    op$loss$upper,op$loss$lower,l$ngrid/l$N_unit,c$low,x,d,u,
    as.integer(seg$duration_days),cfg$DT,n0,isTRUE(cfg$Crowding),cfg$K)
  m <- data.table::as.data.table(sim$metrics);m[,local_day:=0:seg$duration_days]
  counts <- exp(m$log_cells)
  picked <- native$ivt_extract_passage_end_state(list(Ntot_live_obs=counts,
    live_state_obs=t(sweep(sim$X,2,counts,"*"))),n0,l$ngrid,inventory$target,seg$obs_days_local)
  at <- as.integer(picked$selected_day)+1L
  # This is native zero-mass fallback, not an inoculum feasibility stop. If it
  # ever occurs, ancestry is unidentified; do not invent a provenance label.
  fallback <- !is.finite(sum(sim$X[,at]*counts[at]))||sum(sim$X[,at]*counts[at])<=0
  before <- c(NA_real_,head(m$log_cells,-1L))
  m[,`:=`(log_start_live_cells=before,
    expected_HL_daughters=HL_per_start_cell*exp(before),
    expected_LL_daughters=LL_per_start_cell*exp(before),ancestry_identifiable=ancestry_known)]
  if(!ancestry_known) {
    cols <- grep("^(direct|other|historical|never)_",names(m),value=TRUE)
    m[,(cols):=NA_real_]
  }
  row <- data.table::data.table(segment_id=seg$segment_id,initial_cells=n0,
    raw_initial_cells=inventory$initial_raw,raw_final_cells=inventory$final_raw,
    original_O2=inventory$original_O2,target_cells=inventory$target,
    duration_days=seg$duration_days,observation_days=paste(seg$obs_days_local,collapse=","),
    parent_selected_cells=parent_cells,
    reseeding_scale=if(is.finite(parent_cells)&&parent_cells>0)n0/parent_cells else NA_real_,
    insufficient_stock_recorded_only=is.finite(parent_cells)&&parent_cells<n0,
    selected_day=picked$selected_day,selected_live_cells=picked$selected_live_cells,
    elapsed_start=elapsed,elapsed_selected_end=elapsed+picked$selected_day,
    uniform_fallback=fallback,ancestry_identifiable=ancestry_known,
    subset_error=sim$subset_error,roundoff_corrections=sim$roundoff_corrections)
  list(metrics=m,X=sim$X,D=sim$D[c$low,,drop=FALSE],U=sim$U[c$low,,drop=FALSE],
    next_x=picked$selected_frac,next_d=if(fallback)numeric(length(x)) else sim$D[,at],
    next_u=if(fallback)numeric(length(x)) else sim$U[,at],
    known=ancestry_known&&!fallback,row=row,selection=picked,
    initial=list(x=x,d=d,u=u,n0=n0),segment=seg)
}
f7p_prefix <- function(op,native,protocol,initial) {
  stopifnot(initial%in%c(2,4));out <- list();inv <- protocol$inventory
  for(j in seq_len(nrow(inv))) {
    item <- inv[j];parent <- item$parent;seg <- protocol$adapter$segments[[j]]
    seg$oxygen_pct <- op$meta$O2_pct
    if(is.na(parent)) {
      x <- fp_pure_initial(op$l,initial);d <- x*0;u <- x
      elapsed <- 0;count <- NA_real_;known <- TRUE
    } else {
      z <- out[[parent]];x <- z$next_x;d <- z$next_d;u <- z$next_u
      elapsed <- z$row$elapsed_selected_end;count <- z$selection$selected_live_cells;known <- z$known
    }
    out[[item$segment_id]] <- f7p_sim_segment(op,native,seg,item,x,d,u,elapsed,count,known)
  };out
}
f7p_route <- function(op,native,protocol,prefix,route,days=1000L,initial=2) {
  stopifnot(days>0,days==as.integer(days),route%in%c("control","O1","O2"))
  path <- protocol$routes[[route]];template <- protocol$templates[[route]]
  n <- length(op$l$ngrid);low <- op$c$low
  X <- matrix(NA_real_,n,days+1L);D <- U <- matrix(NA_real_,sum(low),days+1L)
  metrics <- events <- checkpoints <- list();elapsed <- 0;previous <- NULL;iteration <- 0L
  add <- function(z,phase,index,source_id) {
    remaining <- days-elapsed;used <- min(z$selection$selected_day,remaining)
    keep <- seq_len(used+1L);dest <- elapsed+keep
    # At a passage boundary report the pre-passage endpoint once. Next day's
    # birth window starts from the reseeded inoculum, without counting passage.
    if(elapsed>0){keep <- keep[-1L];dest <- dest[-1L]}
    X[,dest] <<- z$X[,keep,drop=FALSE];D[,dest] <<- z$D[,keep,drop=FALSE];U[,dest] <<- z$U[,keep,drop=FALSE]
    mm <- data.table::copy(z$metrics[keep]);mm[,`:=`(day=elapsed+local_day,
      route=route,phase=phase,route_segment=index,template_source_id=source_id)]
    metrics[[length(metrics)+1L]] <<- mm
    rr <- data.table::copy(z$row);rr[,`:=`(route=route,phase=phase,route_segment=index,
      template_source_id=source_id,retained_days=used,
      retained_end=elapsed+used,passage_after_segment_applied=elapsed+z$selection$selected_day<days,
      selected_window_truncated=used<z$selection$selected_day)]
    events[[length(events)+1L]] <<- rr
    if(index==1L||(phase=="extrapolation"&&elapsed>=500&&length(checkpoints)<2L))
      checkpoints[[length(checkpoints)+1L]] <<- list(initial=z$initial,segment=z$segment,
        row=rr,X=z$X,metrics=z$metrics,selection=z$selection)
    elapsed <<- elapsed+used;previous <<- z
  }
  for(id in path) {
    if(elapsed>=days)break
    iteration <- iteration+1L;add(prefix[[id]],"experimental_replay",iteration,id)
  }
  prefix_end <- elapsed
  k <- 0L
  while(elapsed<days) {
    k <- k+1L;which_template <- (k-1L)%%length(template)+1L;id <- template[which_template]
    item <- protocol$inventory[match(id,protocol$inventory$segment_id)]
    seg <- protocol$adapter$segments[[item$index]]
    seg$segment_id <- paste0(route,"::extension_",k,"::",id);seg$oxygen_pct <- op$meta$O2_pct
    z <- f7p_sim_segment(op,native,seg,item,previous$next_x,previous$next_d,previous$next_u,
      elapsed,previous$selection$selected_live_cells,previous$known)
    iteration <- iteration+1L;add(z,"extrapolation",iteration,id)
  }
  mm <- data.table::rbindlist(metrics);ev <- data.table::rbindlist(events)
  stopifnot(identical(mm$day,as.numeric(0:days))||identical(mm$day,0:days),
    !anyDuplicated(mm$day),all(is.finite(X)),max(abs(colSums(X)-1))<1e-10,
    f7p_max(abs(mm$A+mm$B-100))<1e-10,all(ev$retained_days>0),
    all(ev$selected_day>0),all(ev$selected_day<=ev$duration_days))
  list(meta=cbind(op$meta,initial_ploidy=initial,route=route,block=protocol$block),
    metrics=mm,events=ev,X=X,D=D,U=U,ngrid=op$l$ngrid,low=low,
    coefficients=list(lambda=op$l$lambda[low],ll=op$c$ll[low]),
    prefix_end=prefix_end,template=template,checkpoints=checkpoints,
    source_note="Conditional reseeding; no inventory feasibility gate. Protocol spread is not a CI.")
}
f7p_summary <- function(routes) {
  z <- data.table::rbindlist(lapply(routes,function(r)cbind(r$meta,
    r$metrics[,setdiff(names(r$metrics),names(r$meta)),with=FALSE])),fill=TRUE)
  # Do not pool daughter counts from distinct protocol scenarios. Summarize
  # paired route-specific shares with equal route weights.
  data.table::rbindlist(lapply(c(A="flux_percent",B="growth_percent"),function(metric) {
    column <- if(metric=="flux_percent")"A" else "B"
    z[,{
      v <- get(column);v <- v[is.finite(v)]
      q <- if(length(v))quantile(v,c(.25,.5,.75),type=7,names=FALSE) else rep(NA_real_,3)
      list(q25=q[1],median=q[2],q75=q[3],n_valid=length(v),n_protocol_routes=.N,
        interval_type="equal-route protocol IQR; not CI")
    },by=.(family,seed_number,O2_pct,p_misseg,maximum_ploidy,beta,threshold_N,initial_ploidy,block,day)][,
      `:=`(metric=metric,analysis="instantaneous",source_accounting="daily_viable_low_daughters")]
  }))
}
f7p_compare <- function(old,now,limits) {
  a <- old$metrics;z <- now$metrics;stopifnot(identical(a$day,z$day),nrow(now$X)>=nrow(old$X))
  cdf <- source_cdf <- 0
  for(j in seq_len(ncol(now$X))) {
    ax <- c(old$X[,j],rep(0,nrow(now$X)-nrow(old$X)))
    cdf <- max(cdf,max(abs(cumsum(ax)-cumsum(now$X[,j]))))
    for(g in c("D","U")) {
      aa <- old[[g]][,j];bb <- now[[g]][,j]
      if(sum(aa)>1e-12&&sum(bb)>1e-12)
        source_cdf <- max(source_cdf,max(abs(cumsum(aa/sum(aa))-cumsum(bb/sum(bb)))))
    }
  }
  calendars_same <- identical(old$events$template_source_id,now$events$template_source_id)&&
    identical(old$events$selected_day,now$events$selected_day)&&
    identical(old$events$retained_days,now$events$retained_days)
  A <- f7p_distance(a$A,z$A);B <- f7p_distance(a$B,z$B)
  mp <- f7p_distance(a$mean_ploidy,z$mean_ploidy);ln <- f7p_distance(a$log_cells,z$log_cells)
  upper <- f7p_max(z$upper_loss_fraction)
  cbind(now$meta,previous_maximum_ploidy=old$meta$maximum_ploidy,
    max_A_change_pp=A,max_B_change_pp=B,max_mean_ploidy_change=mp,
    max_log_cells_change=ln,max_state_cdf_change=cdf,max_upper_loss_fraction=upper,
    max_source_conditional_cdf_change=source_cdf,passage_calendar_identical=calendars_same,
    passed=A<=limits$share_pp&&B<=limits$share_pp&&mp<=limits$mean_ploidy_N&&
      ln<=limits$log_cells&&cdf<=limits$state_cdf&&upper<=limits$upper_loss_fraction&&
      calendars_same&&!any(now$events$uniform_fallback))
}
