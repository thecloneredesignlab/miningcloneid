# Local production runner, using the frozen Step 05 definitions and unmodified C++.
out <- "oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03"
bridge <- "oxygen/results/structural_identifiability/20261116_culture_step05_bridge_design_01"
step02 <- "oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step02"
source("oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/structural_first/culture_step05_bridge_design/api.R")

branch_context <- function() {
  env <- sf_backend(bridge, out, native=TRUE)
  b <- sf_context(bridge, env)
  keep <- c("lam_max", "p_mis_base", "p_misseg", "k_o_mis", "p_wgd",
            "buffer_smax", "buffer_beta", "buffer_n_exp", "alpha_o2", "gamma_growth",
            "mu_hp", "gamma_mu", "O2_crit", "n_O", "init_mean_2N", "init_sd_2N",
            "init_mean_4N", "init_sd_4N", "sigma_kary")
  bounds <- read.delim(file.path(step02,"anchor.tsv"))
  bounds <- bounds[match(keep, bounds$parameter), ]
  stopifnot(identical(bounds$parameter, keep), length(b$adapters) == 2L,
            identical(unname(vapply(b$adapters, function(a) length(a$segments), integer(1))), c(62L,60L)))
  logscale <- keep %in% c("lam_max","p_mis_base","p_misseg","k_o_mis","p_wgd",
                         "buffer_beta","mu_hp","O2_crit","sigma_kary")
  encode <- function(p) {
    x <- unlist(p[keep]); lo <- bounds$lower; hi <- bounds$upper
    x[logscale] <- log(x[logscale])
    lo[logscale] <- log(lo[logscale])
    hi[logscale] <- log(hi[logscale])
    as.numeric((x-lo)/(hi-lo))
  }
  decode <- function(u) {
    stopifnot(length(u)==19L,all(is.finite(u)),all(u>=0 & u<=1))
    lo <- bounds$lower; hi <- bounds$upper
    value <- lo+u*(hi-lo)
    value[logscale] <- exp(log(lo[logscale])+u[logscale]*(log(hi[logscale])-log(lo[logscale])))
    value[u==0] <- lo[u==0]
    value[u==1] <- hi[u==1]
    p <- b$p
    p[keep] <- as.list(value)
    p
  }
  obs <- lapply(b$objects$fit_data, env$ivt_observed_passage_summary)
  growth_ids <- names(obs)[vapply(obs,function(r)is.finite(r$observed_growth),FALSE)]
  kary_ids <- names(obs)[vapply(obs,function(r)length(r$observed_kary)>0,FALSE)]
  stopifnot(length(growth_ids)==114L,length(kary_ids)==12L)
  mask_inputs <- function(oxygens) {
    d <- expand.grid(N=22:154, O=unique(oxygens))
    d$n <- lapply(d$N,function(n)seq.int(0L,n))
    list(N=d$N,O=d$O,m=unlist(d$n),row=rep(seq_len(nrow(d)),lengths(d$n)))
  }
  standard_oxy <- sort(unique(unlist(lapply(b$adapters, function(a)vapply(a$segments,function(s)s$oxygen_pct,0)))))
  standard_m <- mask_inputs(standard_oxy)
  masks <- function(p, m=standard_m) {
    h <- p$O2_crit^p$n_O/(p$O2_crit^p$n_O + m$O^p$n_O)
    mu <- p$mu_hp * h * (m$N/44)^p$gamma_mu
    prob <- p$p_mis_base+p$p_misseg*mu/(p$k_o_mis+mu)
    dbinom(m$m,m$N[m$row],prob[m$row])>=1e-8
  }
  evaluate <- function(p, oxy=NULL, keep_details=FALSE, exact_kary=FALSE) {
    a <- b$adapters
    if (!is.null(oxy)) {
      a <- lapply(a,function(z) {
        z$segments <- lapply(z$segments,function(s){s$oxygen_pct<-oxy;s})
        z
      })
    }
    runs <- list()
    err <- NULL
    for (cohort in names(a)) {
      z <- tryCatch(env$ivt_run_lineage(a[[cohort]],b$cfg,p,b$core),
        invitro_protocol_infeasible=function(e){err<<-e;NULL})
      if (is.null(z)) return(list(ok=FALSE,failure=paste(err$cohort,err$parent_segment$segment_id,
         err$child_segment$segment_id,sep="|"),maximum=err$selection$available_cells,
         required=err$selection$required_cells))
      runs[[cohort]] <- z
    }
    if (!is.null(oxy)) m <- mask_inputs(oxy) else m <- standard_m
    y <- setNames(rep(NA_real_,114L),growth_ids)
    n <- b$core$grid_pre
    k <- matrix(NA_real_,nrow=length(n),ncol=length(kary_ids),
                dimnames=list(as.character(n),kary_ids))
    latent <- k
    days <- parent_days <- setNames(rep(NA_real_,122L),character(122L))
    fallback <- rep(FALSE,122L); edge_supply <- rep(NA_real_,122L)
    position <- 0L
    storage <- if (keep_details) vector("list",122L) else NULL
    for (cohort in names(runs)) for (seg in runs[[cohort]]$segment_results) {
      position <- position+1L
      label <- paste(cohort,seg$segment$segment_id,sep="|")
      names(days)[position] <- names(parent_days)[position] <- label
      days[position] <- seg$selection$selected_day
      if (!is.null(seg$parent_boundary_selection)) {
        par <- runs[[cohort]]$segment_results[[seg$segment$parent_segment_id]]
        q <- seg$parent_boundary_selection
        parent_days[position] <- q$selected_day
        fallback[position] <- q$selected_index != par$selection$selected_index
        edge_supply[position] <- sum(par$sim$live_state_obs[par$selection$selected_index,])/q$required_cells
      }
      gids <- intersect(seg$segment$data_ids,growth_ids)
      if (length(gids)) y[gids] <- env$ivt_log_growth_rate(seg$sim$Ntot_live_obs[1L],
                                        seg$selection$selected_live_cells,seg$selection$selected_day)
      kids <- intersect(seg$segment$data_ids,kary_ids)
      if (length(kids)) {
        z <- seg$selection$selected_frac
        if (!exact_kary && p$sigma_kary <= 0.015 && p$sigma_kary >= 0.01) {
          # Mathlib returns off-diagonal zeros here, including the Step 02 interval.
          kz <- z/sum(z)
        } else kz <- env$ivt_smooth_kary_distribution(n,z,p$sigma_kary)
        for (id in kids) {k[,id]<-kz;latent[,id]<-z}
      }
      if (keep_details) storage[[position]] <- list(label=label,states=seg$sim$live_state_obs,
         selection=seg$selection,boundary=seg$parent_boundary_selection)
    }
    stopifnot(position==122L,all(is.finite(y)),all(is.finite(k)),all(k>=0),
              max(abs(colSums(k)-1))<=1e-12,sum(is.finite(parent_days))==120L)
    list(ok=TRUE,growth=y,kary=k,latent=latent,days=days,parent_days=parent_days,
      fallback=fallback,supply=edge_supply,mask=masks(p,m),details=storage)
  }
  list(b=b,env=env,bounds=bounds,keep=keep,encode=encode,decode=decode,evaluate=evaluate,
       anchor=b$p,u=pmin(1,pmax(0,encode(b$p))),oxygens=standard_oxy)
}

growth_distance <- function(a,b) max(abs(a$growth-b$growth))
kary_distance <- function(a,b) max(colSums(abs(a$kary-b$kary))/2)
differences <- function(a,b) c(growth_max=growth_distance(a,b),
  growth_rms=sqrt(mean((a$growth-b$growth)^2)),kary_tv_max=kary_distance(a,b),
  kary_tv_rms=sqrt(mean((colSums(abs(a$kary-b$kary))/2)^2)),
  endpoint=sum(a$days!=b$days),mask=sum(xor(a$mask,b$mask)),
  reset_mode=sum(xor(a$fallback,b$fallback)),parent_day=sum(a$parent_days!=b$parent_days,na.rm=TRUE))

emit_parameters <- function(parameters,path,keep) {
  values <- do.call(rbind,lapply(names(parameters),function(id) {
    p <- parameters[[id]]
    data.frame(id=id,parameter=keep,binary64=vapply(p[keep],function(v)sprintf("%a",v),""))
  }))
  sf_emit(values,path)
}
