# Daily viable-daughter source partition and passive origin-population telemetry.
# Quantiles are across matched passage repeats, never across chromosomes/cells.
f7_count_summary <- function(d,meta,nominal) {
  fields <- c("HL_per_start_live_cell","LL_per_start_live_cell","log_start_live_cells")
  data.table::rbindlist(lapply(fields,function(n)cbind(meta,counter=n,fp_summaries(d[[n]],nominal))))
}
f7_population_summary <- function(d,meta,nominal) {
  data.table::rbindlist(lapply(seq_along(d$metrics),function(i) {
    g <- (i-1L)%/%length(d$metric_names)+1L;m <- (i-1L)%%length(d$metric_names)+1L
    cbind(meta,group=d$group_names[g],metric=d$metric_names[m],fp_summaries(d$metrics[[i]],nominal))
  }))
}
f7_distribution_summary <- function(d,meta,l,c) {
  # Conditional within-group cell distributions, averaged with EQUAL weight for
  # each nonempty repeat. Empty groups do not become fictitious zero-rate cells.
  k <- nrow(d$high_parent_percent);low <- which(c$low);out <- list()
  variables <- list(ploidy=l$ngrid[low]/l$N_unit,division_rate=l$lambda[low],low_birth_rate=c$ll[low])
  for(s in seq_along(d$snapshot_days)) {
    ix <- (s-1L)*k+seq_len(k);X <- d$snapshot_X[low,ix,drop=FALSE]
    U <- d$snapshot_U[low,ix,drop=FALSE];D <- d$snapshot_D[low,ix,drop=FALSE]
    groups <- list(all_low=X,direct_first_generation=D,other_low=X-D,historical_lineage=X-U,never_HL_lineage=U)
    for(g in names(groups))for(v in names(variables)) {
      weights <- pmax(groups[[g]],0);mass <- colSums(weights);valid <- mass>0
      # Rate CDFs include the actual contemporaneous crowding factor of each
      # repeat. Ploidy CDF uses the exact chromosome lattice; no KDE invention.
      values <- outer(variables[[v]],if(v=="ploidy")rep(1,k) else d$snapshot_crowding[ix],"*")
      grid <- if(v=="ploidy")sort(unique(variables[[v]])) else
        sort(unique(c(0,seq(0,max(variables[[v]]),length.out=101L))))
      cdfs <- matrix(NA_real_,k,length(grid))
      if(any(valid))for(r in which(valid)) {
        ord <- order(values[,r]);xx <- values[ord,r];cw <- cumsum(weights[ord,r]/mass[r])
        at <- findInterval(grid,xx);cdfs[r,] <- c(0,cw)[at+1L]
      }
      getq <- function(p)vapply(seq_along(grid),function(j)if(any(valid))
        unname(quantile(cdfs[valid,j],p,type=7)) else NA_real_,0.0)
      out[[length(out)+1L]] <- cbind(meta,group=g,day=d$snapshot_days[s],variable=v,
        data.table::data.table(value=grid,cdf=if(any(valid))colMeans(cdfs[valid,,drop=FALSE]) else NA_real_,
          q25=getq(.25),q75=getq(.75),n_valid_repeats=sum(valid),n_repeats=k))
    }
  }
  data.table::rbindlist(out)
}
