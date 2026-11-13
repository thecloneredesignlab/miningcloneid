# Static report from independently accepted full-data fits. No optimization.
d <- file.path(study, "reports"); dir.create(d, recursive = TRUE, showWarnings = FALSE)
cs <- read_table(file.path(study, "tables/candidates.tsv"))
quality <- read_table(file.path(study, "tables/point_quality.tsv"))
baseline <- read_kv(file.path(study, "tables/baseline_gate.tsv"))
cs <- cs[order(cs$objective, cs$task_id), ]
selected <- cs[!duplicated(cs$point), ]
selected <- selected[match(c("base", quality$point), selected$point), ]
fits <- lapply(selected$source_fit, function(path) readRDS(file.path(project, path)))
names(fits) <- selected$point
roles <- list(); components <- list(); growth <- list(); distributions <- list(); flows <- list(); days <- list(); curves <- list(); chromosomes <- list(); local_quality <- list()
for (i in seq_len(nrow(selected))) {
  row <- selected[i, ]; fit <- fits[[i]]; comp <- fit$best_components
  fixed <- if (row$point == "base") list() else setNames(list(quality$value[quality$point == row$point]), quality$target[quality$point == row$point])
  map <- original_map(ctx, fixed); tab <- map$table
  tab$estimate <- vapply(tab$parameter, function(n) fit$best_params[[n]], numeric(1))
  tab$point <- row$point;tab$task_id <- row$task_id;tab$normalized_coordinate <- NA_real_
  tab$normalized_coordinate[tab$status == "active"] <- o2review_unit_spec(map, fit$best_params)$unit
  tab$exact_lower <- tab$estimate == tab$lower;tab$exact_upper <- tab$estimate == tab$upper
  tab$near_lower <- tab$normalized_coordinate <= .001;tab$near_upper <- tab$normalized_coordinate >= .999
  roles[[i]] <- tab
  components[[i]] <- data.frame(point = row$point, task_id = row$task_id, t(score_components(comp)), check.names = FALSE)
  g <- comp$growth_df;g$point <- row$point;growth[[i]] <- g
  f <- comp$flow_overlay_df;f$point <- row$point;flows[[i]] <- f
  for (cohort in c("2N", "4N")) {
    run <- comp[[paste0("run_", cohort)]]
    for (j in seq_along(run$segment_results)) {
      seg <- run$segment_results[[j]]
      days[[length(days) + 1L]] <- data.frame(point = row$point, cohort = cohort, segment_order = j,
        segment_id = seg$segment$segment_id, selected_day = seg$selection$selected_day, selected_live_cells = seg$selection$selected_live_cells)
      distributions[[length(distributions) + 1L]] <- data.frame(point = row$point, cohort = cohort, segment_order = j,
        segment_id = seg$segment$segment_id, chromosome_number = run$grid_pre,
        live_fraction = seg$selection$selected_frac, propagated_cells = seg$selection$reseeded_state)
      for (pid in seg$segment$data_ids) {
        obs <- ctx$backend$ivt_observed_passage_summary(ctx$objects$fit_data[[pid]])
        if (!length(obs$observed_kary)) next
        grid <- sort(unique(c(run$grid_pre, obs$observed_kary)))
        K <- outer(run$grid_pre, run$grid_pre, function(a,b) dnorm(a,b,comp$sigma_kary))
        sums <- colSums(K); sums[!is.finite(sums) | sums <= 0] <- 1
        predicted <- as.numeric(sweep(K,2,sums,"/") %*% seg$selection$selected_frac)
        prob <- predicted[match(grid,run$grid_pre)];prob[is.na(prob)] <- 0
        chromosomes[[length(chromosomes)+1L]] <- data.frame(point=row$point,cohort=cohort,passage_id=pid,
          chromosome_number=grid,observed_fraction=as.numeric(table(factor(obs$observed_kary,levels=grid)))/length(obs$observed_kary),
          predicted_observation_fraction=prob,scoring_probability_with_floor=pmax(prob,1e-12))
      }
    }
  }
  for (oxygen in c(0,.5,2,5,20.5)) for (N in c(44L,88L)) for (m in 0:5) {
    p <- fit$best_params;survival <- p$buffer_smax * exp(-p$buffer_beta*(2*ctx$cfg$N_UNIT/N)^p$buffer_n_exp)
    probability <- ctx$backend$.pmisseg_of_O2(oxygen,p,N)
    z <- 1-probability+probability*survival
    curves[[length(curves)+1L]] <- data.frame(point=row$point,N=N,oxygen_pct=oxygen,error_size=m,
      p_per_chromosome=probability,per_error_survival=survival,error_size_survival=if(m==0) 1 else survival^m,
      conditional_surviving_error_probability=if(z>0) probability*survival/z else NA_real_,
      biological_viable_daughters_per_ordinary_division=2*z^N)
  }
}
for (i in seq_len(nrow(cs))) if(cs$mode[i]=="local") {
  path <- dirname(file.path(project,cs$source_fit[i]));l <- readRDS(file.path(path,"local_search.rds"));sb <- readRDS(file.path(path,"sbplx.rds"))
  local_quality[[length(local_quality)+1L]] <- data.frame(task_id=cs$task_id[i],point=cs$point[i],origin=cs$origin[i],family=cs$family[i],phase=cs$phase[i],
    objective=cs$objective[i],local_stable=l$tested_directional_stability,overall_budget_limited=l$budget_limited,
    native_sbplx_status=sb$status,native_sbplx_message=sb$message,evaluations=l$evaluations,
    improvement_records=length(list.files(file.path(path,"best_candidates/improvements"),pattern="\\.rds$")))
}
combine <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
save <- function(x,n) o2prof_write_tsv(x,file.path(d,n))
parameters <- combine(roles);comp <- combine(components);g <- combine(growth);dist <- combine(distributions);flow <- combine(flows);fun <- combine(curves);kary <- combine(chromosomes)
save(parameters,"nuisance_paths.tsv");save(comp,"winner_components.tsv");save(g,"growth_predictions.tsv")
save(dist,"full_chromosome_distributions.tsv");save(flow,"flow_overlays.tsv");save(kary,"observed_karyotype_distributions.tsv")
save(combine(days),"selected_days.tsv");save(fun,"mechanism_functions.tsv");save(combine(local_quality),"local_quality.tsv")
spread <- vapply(quality$point,function(point) {z <- cs$objective[cs$point==point & cs$mode=="local"];if(length(z)) diff(range(z)) else NA_real_},numeric(1))
quality$shape_below_recovery_resolution <- quality$delta<=.001 | quality$delta<=spread
save(quality,"profile_point_summary.tsv")
ledger <- list();diag_status <- list()
for (point in selected$point) {
  path <- file.path(study,"diagnostics",point)
  complete <- file.exists(file.path(path,"COMPLETE"))
  diag_status[[point]] <- data.frame(point=point,event_ledger_verified=complete)
  if(complete) {
    x <- read_table(file.path(path,"segment_event_ledger.tsv"));x$point <- point;ledger[[point]] <- x
  }
}
save(combine(diag_status),"event_diagnostic_status.tsv")
if(length(ledger)) {
  led <- combine(ledger);save(led,"segment_event_ledger.tsv");totals <- list()
  keys <- unique(led[c("point","cohort","recording_scope")])
  for(i in seq_len(nrow(keys))) {
    key <- keys[i,];x <- led[led$point==key$point & led$cohort==key$cohort & led$recording_scope==key$recording_scope,]
    ordinary <- sum(x$division)*(1-fits[[key$point]]$best_params$p_wgd)
    totals[[i]] <- data.frame(key,segments=nrow(x),ordinary_divisions=ordinary,
      biological_failed_daughters_per_ordinary_division=sum(x$ordinary_biological_nonviability)/ordinary,
      grid_lost_viable_daughters_per_ordinary_division=sum(x$ordinary_boundary_loss+x$ordinary_boundary_gain)/ordinary,
      numerical_pruning_per_ordinary_division=sum(x$numerical_tail_pruning)/ordinary,
      retained_ordinary_daughters_per_ordinary_division=sum(x$ordinary_retained_loss+x$ordinary_retained_gain+x$ordinary_unchanged)/ordinary)
  }
  save(combine(totals),"event_summary_by_recording_scope.tsv")
}
colors <- c("black","#D55E00","#0072B2","#009E73","#CC79A7","#E69F00","#56B4E9","#777777")
names(colors) <- selected$point
pdf(file.path(d,"culture_profile_review.pdf"),width=11,height=8.5,useDingbats=FALSE)
par(mfrow=c(1,2),mar=c(7,5,4,1))
for(target in c("p_misseg","buffer_beta")) {
  z <- quality[quality$target==target,];z <- z[order(z$value),];scale <- if(target=="p_misseg") 1e-4 else 1e-3
  x <- log1p(z$value/scale);endpoints <- cs[cs$point %in% z$point,]
  y <- endpoints$objective-baseline$objective
  plot(x,z$delta,type="n",xaxt="n",ylim=range(c(0,y),finite=TRUE),xlab="",ylab="Objective minus unrestricted reference",main=target)
  points(log1p(z$value[match(endpoints$point,z$point)]/scale),y,col="grey75",pch=16,cex=.55)
  points(x,z$delta,pch=ifelse(z$pass_gate,19,1),col="#194B7A",cex=1.3)
  axis(1,at=x,labels=format(z$value,digits=4),las=2);mtext("Fixed natural value; zero-capable display spacing",1,line=5.5,cex=.8)
  abline(h=0,lty=3,col="grey60")
  legend("topleft",c("Checks passed","Checks unresolved","All accepted endpoints"),pch=c(19,1,16),col=c("#194B7A","#194B7A","grey75"),cex=.7,bty="n",inset=.02)
}
par(mfrow=c(1,2),mar=c(10,5,4,1))
for(prefix in c("pm_","bb_")) {
  ids <- c("base",quality$point[startsWith(quality$point,prefix)])
  mat <- sapply(ids,function(id) parameters$normalized_coordinate[parameters$point==id])
  matplot(seq_len(nrow(mat)),mat,type="b",col=colors[ids],lty=1,pch=seq_along(ids),xaxt="n",yaxt="n",ylim=c(0,1.35),xlab="",ylab="Transformed nuisance coordinate",main=if(prefix=="pm_") "Induction compensation" else "Survival compensation")
  axis(1,at=seq_len(nrow(mat)),labels=parameters$parameter[parameters$point=="base"],las=2,cex.axis=.6);axis(2,at=seq(0,1,.2))
  legend("topright",ids,col=colors[ids],lty=1,pch=seq_along(ids),cex=.7,bty="n",inset=.02)
}
par(mfrow=c(2,2),mar=c(4,5,4,1))
for(N in c(44,88)) for(quantity in c("p_per_chromosome","per_error_survival")) {
  x <- fun[fun$N==N & fun$error_size==0,];ys <- lapply(selected$point,function(id) x[x$point==id,quantity]);oxygen <- c(0,.5,2,5,20.5)
  matplot(log1p(oxygen/.025),do.call(cbind,ys),type="l",col=colors,lty=1,xlab="Oxygen percent; log1p spacing",ylab=quantity,xaxt="n",main=paste("Parent N =",N))
  axis(1,at=log1p(oxygen/.025),labels=oxygen)
  legend("topright",selected$point,col=colors,lty=1,cex=.55,bty="n")
}
par(mfrow=c(1,1),mar=c(5,5,4,1))
plot(g$observed_growth,g$predicted_growth_rate,col=colors[g$point],pch=16,cex=.5,xlab="Observed growth rate",ylab="Predicted growth rate",main="Saved full-data growth predictions")
abline(0,1,lty=2);legend("bottomright",selected$point,col=colors,pch=16,cex=.7,bty="n",inset=.02)
par(mfrow=c(1,2),mar=c(5,5,4,1))
for(cohort in c("2N","4N")) {
  z <- kary[kary$cohort==cohort,];pid <- tail(unique(z$passage_id[z$point=="base"]),1)
  z <- z[z$passage_id==pid,];base_obs <- z[z$point=="base",]
  plot(base_obs$chromosome_number,base_obs$observed_fraction,type="h",col="grey65",xlab="Chromosome number",ylab="Observed/predicted fraction",ylim=range(c(z$observed_fraction,z$predicted_observation_fraction)),main=paste(cohort,"last scored karyotype"))
  for(id in selected$point){a <- z[z$point==id,];lines(a$chromosome_number,a$predicted_observation_fraction,col=colors[id])}
  legend("topright",selected$point,col=colors,lty=1,cex=.55,bty="n")
}
invisible(dev.off())
status <- if(all(quality$pass_gate)) "PASSED" else "UNRESOLVED"
lines <- c("# Seven point culture profile review","",paste("Operational seven-point status:",status),"",
  sprintf("The validated unrestricted reference is %.17g. All deltas use this reference; negative values are retained.",baseline$objective),"",
  "This study reviews four induction and three survival points in the original zero-extended domain. These are partial diagnostic profile-objective results, without calibrated confidence intervals, likelihood-ratio cutoffs or structural-identifiability claims.","",
  "## Results and optimization quality","",sprintf("Points passing all gates: %d of 7. See profile_point_summary.tsv and local_quality.tsv.",sum(quality$pass_gate)),"",
  "Independent recovery requires two new origins from A and B within 0.001 of the final minimum, each with its own completed global/local chain. Historical, transferred and neighbor origins remain separate. Native global and local caps, fine-mesh stability, final repeats and boundary hits are separate fields. Small score rises flagged below recovery resolution do not support precise local shape claims.","",
  "## Mechanisms and measurements","",
  "nuisance_paths.tsv and mechanism_functions.tsv show parameter compensation and model functions. Fixed/inactive coordinates are labeled. Growth, full chromosome distributions, flow overlays and selected days are exported with keys. Full-data component traces do not measure a stream's unique information contribution.","",
  "event_summary_by_recording_scope.tsv, when present, separates biological failed daughters, grid-lost viable daughters, numerical pruning and retained daughters per ordinary division. Only diagnostics with COMPLETE and a passed kernel/dense-event balance are included. Recorded-growth and other protocol segments are separate. Missing or failed ledgers are withheld; no flux claim is inferred from normalized distributions.","",
  "## Scope and resource closure","",
  "All new execution uses R, Python, shell and pinned C++ with typed TSV/RDS state; no Java, JavaScript or JSON execution/transport/reporting dependency. All allocations use xxlarge with no application concurrency cap. The new hard ceiling is 4000 core-hours, with a 120-core-hour non-fit cap, 36 global-origin attempts and 70 fitting-local attempts. Final accounting must be refreshed after this report and the controller close.","",
  "No additional thirteen-point review, expanded domain, stream deletion, COMBOS, tumor/joint fitting, statistical calibration, commit or push was performed.")
writeLines(lines,file.path(study,"REPORT.md"));assert_runtime(d);writeLines("Completed static English report",file.path(d,"COMPLETE"))
