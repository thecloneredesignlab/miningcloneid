#!/usr/bin/env Rscript
# Chart contract: stochastic variability and actual passage events, not a claim
# that the ensemble mean must differ from continuous culture. Existing Figure6
# palette and static headless PDF/PNG delivery; no manuscript text edits.
root <- normalizePath(Sys.getenv("FIGURE_WORKSPACE_ROOT", getwd()),mustWork=TRUE)
source(file.path(root,"Code","Figures","draw_Figure6.R"))
f6r_require_packages(c("ggplot2","patchwork","data.table"))
paths <- f6r_paths(root); run <- f6g_paths(paths)
trace_files <- list.files(run$run_root,pattern="^stochastic_trace_.*[.]rds$",full.names=TRUE)
stopifnot(length(trace_files)>0L)
trajectory_rows <- event_rows <- list()
for (path in trace_files) {
  object <- readRDS(path)
  for (trace in object$traces) {
    if (!trace$initial_ploidy %in% c(2,4)) next
    if (object$O2_pct==.5) {
      x <- trace$trajectories
      x <- x[!duplicated(x$day) & x$day%%10==0,,drop=FALSE]
      x$continuous <- trace$continuous[x$day+1L]
      x$pair_label <- object$endpoint$pair_label[[1]]
      x$initial <- paste0(trace$initial_ploidy,"N")
      x$p_misseg <- paste0("p_misseg = ",object$p_misseg)
      trajectory_rows[[length(trajectory_rows)+1L]] <- x
    }
    if (trace$initial_ploidy==4 && object$p_misseg==.005 && nrow(trace$events)) {
      x <- trace$events[trace$events$column==1,,drop=FALSE]
      x <- head(x,10L)
      x$pair_label <- object$endpoint$pair_label[[1]]
      x$oxygen <- paste0("Oxygen ",object$O2_pct,"%")
      event_rows[[length(event_rows)+1L]] <- x
    }
  }
}
trajectories <- as.data.frame(data.table::rbindlist(trajectory_rows))
events <- as.data.frame(data.table::rbindlist(event_rows))
precision <- as.data.frame(data.table::rbindlist(lapply(c("c01","c02"),function(family)
  f6r_read_tsv(file.path(run$run_root,paste0("stochastic_precision_",family,".tsv"))))))
precision <- precision[precision$initial_ploidy==4,]
palette <- c("C01"="#C99700","C02"="#6A3D9A")
base <- ggplot2::theme_bw(base_size=8,base_family="Helvetica") +
  ggplot2::theme(panel.grid.minor=ggplot2::element_blank(),legend.position="bottom",
    strip.background=ggplot2::element_rect(fill="#EEEEEE",colour="#BBBBBB"),
    plot.title=ggplot2::element_text(face="bold",size=10))
a <- ggplot2::ggplot(trajectories,ggplot2::aes(day,mean)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin=mean-sd,ymax=mean+sd),fill="#CC79A7",alpha=.22) +
  ggplot2::geom_line(ggplot2::aes(y=replicate_1),colour="#999999",linewidth=.25,alpha=.6) +
  ggplot2::geom_line(ggplot2::aes(y=replicate_2),colour="#999999",linewidth=.25,alpha=.6) +
  ggplot2::geom_line(ggplot2::aes(y=replicate_3),colour="#999999",linewidth=.25,alpha=.6) +
  ggplot2::geom_line(ggplot2::aes(y=continuous,colour="Continuous"),linewidth=.45,linetype="dashed") +
  ggplot2::geom_line(ggplot2::aes(colour="Stochastic mean"),linewidth=.45) +
  ggplot2::facet_grid(pair_label+initial~p_misseg) +
  ggplot2::scale_colour_manual(values=c("Continuous"="#222222","Stochastic mean"="#CC79A7"),name=NULL) +
  ggplot2::labs(title="A  Mean ploidy trajectories at fixed oxygen 0.5%",
    x="Experimental time (day)",y="Mean ploidy (N)") + base
if (nrow(events)) {
  b <- ggplot2::ggplot(events,ggplot2::aes(day,cells_before)) +
    ggplot2::geom_segment(ggplot2::aes(xend=day,yend=cells_after),colour="#999999",linewidth=.4) +
    ggplot2::geom_point(ggplot2::aes(colour="Before sampling"),size=1.4) +
    ggplot2::geom_point(ggplot2::aes(y=cells_after,colour="After sampling"),size=1.4,shape=1) +
    ggplot2::scale_colour_manual(values=c("Before sampling"="#0072B2","After sampling"="#CC79A7"),name=NULL) +
    ggplot2::scale_y_log10(labels=scales::label_number(scale=1e-6,suffix="M")) +
    ggplot2::facet_grid(pair_label~oxygen,scales="free_x") +
    ggplot2::labs(title="B  Actual integer-day passage populations",
      x="Experimental time (day)",y="Live cells") + base
} else {
  b <- ggplot2::ggplot() + ggplot2::theme_void()
}
c <- ggplot2::ggplot(precision,ggplot2::aes(O2_pct,within_endpoint_sd_final_N,colour=pair_label)) +
  ggplot2::geom_line(linewidth=.5) + ggplot2::facet_wrap(~p_misseg,nrow=1,
    labeller=ggplot2::label_both) + ggplot2::scale_colour_manual(values=palette,name="Cluster") +
  ggplot2::labs(title="C  Stochastic dispersion at day 10,000",
    x="Fixed oxygen (%)",y="Stochastic SD (N)") + base
combined <- (a / b / c) + patchwork::plot_layout(heights=c(2.6,1.3,1)) +
  patchwork::plot_annotation(title="Supplementary Figure 6-13. Stochastic passage diagnostics")
directory <- file.path(root,"data","Figures","Supp_Figure6_13")
dir.create(directory,recursive=TRUE,showWarnings=FALSE)
name <- "supp_fig6-13_stochastic_passage_diagnostics"
output <- f6r_save_plot(combined,file.path(directory,name),width=13.5,height=16.5,dpi=300)
published <- f6g_publish_plot(output,paths,name)
f6g_render_hash_validation(output,published,file.path(directory,paste0(name,"_render_validation.tsv")))
