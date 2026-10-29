# Same approved plot builders, data, dimensions, and style; no new calculations.
f78ro_table <- function(ctx,key,name,frame=FALSE) {
  cache <- paste(key,name,frame,sep="|")
  if (!exists(cache,envir=ctx$tables,inherits=FALSE)) {
    path <- f78ro_data(ctx,key,name)
    z <- data.table::fread(path)
    assign(cache,if(frame)as.data.frame(z) else z,envir=ctx$tables)
  }
  get(cache,envir=ctx$tables,inherits=FALSE)
}
f78ro_plot7 <- function(ctx,row) {
  stem <- row$stem;version <- row$analysis;a <- ctx$contexts[[row$key]]
  assign("f7_metric_label",f7_metric_labels(version,direct_growth=identical(row$key,"instantaneous7")),globalenv())
  if (stem %in% c("supp_fig7_1","supp_fig7_2")) return(f7_draw_sources(a,stem))
  if (stem %in% c("supp_fig7_1_ver2","supp_fig7_2_ver2")) return(f7_draw_lineage_diagnostics(a,stem))
  if (stem=="supp_fig7_8") return(f7_draw_origin_populations(ctx,row))
  z <- f78ro_table(ctx,row$key,if(version=="history")"history_time_summary.tsv" else "instant_time_summary.tsv",frame=TRUE)
  if(version=="instantaneous"&&(!"source_accounting"%in%names(z)||
    any(z$source_accounting!="daily_viable_low_daughters")))
    stop("Daily A/B figures require newly validated daily-source data; old B will not be relabeled")
  if (row$pages==8L) return(f78ro_matrix(a,z,if(grepl("supp_fig7_3",stem))"flux_percent" else "growth_percent",stem))
  spec <- data.table::fread(file.path(ctx$workspace,"Code/config/manifests/figure7_display_panels.tsv"))
  original <- if(stem=="figure7_ver2")"assembled_fig7" else sub("_ver2$","",stem)
  s <- spec[spec$stem==original,];stopifnot(nrow(s)==1L)
  ta <- if(version=="instantaneous")"A  High-parent source" else "A  Historical lineage share"
  tb <- if(version=="instantaneous")"B  Low-parent source" else "B  Historical-lineage proliferation"
  A <- f7_panel_block(z,s$O2_pct,s$p_misseg,paste0(ta,"; O2 = ",s$O2_pct,"%; p_misseg = ",s$p_misseg),"flux_percent")
  B <- f7_panel_block(z,s$O2_pct,s$p_misseg,paste0(tb,"; O2 = ",s$O2_pct,"%; p_misseg = ",s$p_misseg),"growth_percent")
  legend <- cowplot::get_legend(f7_line_plot(z[z$task_id==z$task_id[1],])+ggplot2::theme(legend.position="bottom"))
  main <- cowplot::plot_grid(cowplot::plot_grid(A,B,nrow=1),legend,ncol=1,rel_heights=c(1,.04))
  f78_publish(a,main,14,11,stem)
}
f78ro_plot8 <- function(ctx,row) {
  stem <- row$stem;a <- ctx$contexts[[row$key]]
  if (stem=="assembled_fig8") {
    A <- f8_median_iqr(a);B <- f8_labeled_author_figure5(a)
    C <- f78ro_table(ctx,"main8","Figure8C_curve_summary.tsv")
    D <- f78ro_table(ctx,"main8","Figure8D_curve_summary.tsv")
    limits <- range(C$q25,C$q75,D$q25,D$q75)
    pc <- f8_start_curve_plot(C,"C  In vivo",limits);pd <- f8_start_curve_plot(D,"D  In vitro",limits)
    bottom <- (pc+pd)+patchwork::plot_layout(guides="collect")
    bottom <- bottom&ggplot2::theme(legend.position="bottom")
    main <- (A$plot+B+patchwork::plot_layout(widths=c(1,1)))/bottom+patchwork::plot_layout(heights=c(1.35,1))
    return(f78_publish(a,main,13,12,stem))
  }
  if (stem %in% c("supp_fig8_1","supp_fig8_2")) {
    p <- if(stem=="supp_fig8_1").01 else .3
    empirical <- list(sites=f78ro_table(ctx,"main8","empirical_tissue_ploidy_summary.tsv"))
    summaries <- f78ro_table(ctx,"main8","curve_summary.tsv")
    plot <- f8_tissue_plot(empirical,summaries,"fixed_p_misseg",p,
      title=paste0("Tissue and model ploidy; p_misseg = ",p))+ggplot2::labs(subtitle=NULL)+f8_ab_theme()
    return(f78_publish(a,plot,8,8,stem))
  }
  if (stem %in% c("supp_fig8_3","supp_fig8_4")) {
    panels <- f8_mixture_panels(a)
    plot <- if(stem=="supp_fig8_3")panels$C+panels$D else panels$supp3+panels$supp4
    return(f78_publish(a,plot+patchwork::plot_layout(ncol=2),13,5.2,stem))
  }
  z <- f78ro_table(ctx,row$key,paste0(stem,"_curve_summary.tsv"))
  if (stem=="supp_fig8_5") return(f78_publish(a,f8_culture_family_plot(z),13,10,stem))
  context <- if(stem %in% c("supp_fig8_6","supp_fig8_8"))"in vitro" else "in vivo"
  if (stem %in% c("supp_fig8_6","supp_fig8_7")) return(f78_publish(a,f8_time_matrix_plot(z,context),20,9,stem))
  f78_publish(a,f8_family_time_plot(z,context),a$family_config$width_inches,a$family_config$height_inches,stem)
}
f78ro_plot <- function(ctx,row) {
  if (row$figure==7L) f78ro_plot7(ctx,row) else f78ro_plot8(ctx,row)
}
