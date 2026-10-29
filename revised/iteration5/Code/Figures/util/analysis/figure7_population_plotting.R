# Origin profiles are distinct from the complementary daily-parent source audit.
# CDFs are actual discrete cell-weighted distributions, not KDE violin artifacts.
f7_population_pages <- function(z,d) {
  groups <- list(direct=c("all_low","direct_first_generation","other_low"),
    historical=c("all_low","historical_lineage","never_HL_lineage"))
  captions <- list(direct=c("All low","First-generation H-to-L","Other low (includes descendants)"),
    historical=c("All low","H-to-L lineage (including descendants)","Never-H-to-L lineage"))
  metric_labels <- c(cell_share_percent="Current low-cell share (%)",low_birth_share_percent="Low-daughter production share (%)",
    mean_low_birth_rate="Low-daughter output (cell/day)",mean_division_rate="Division rate (1/day)",mean_ploidy="Mean ploidy (N)")
  colors <- c("#222222","#D55E00","#0072B2");plots <- list();index <- list()
  panels <- function(x) {
    x$panel <- factor(paste(x$model_context,x$pair_label,paste0(x$initial_ploidy,"N"),sep="\n"),
      levels=unlist(lapply(c("in vivo","in vitro"),function(e)unlist(lapply(c("C01","C02"),
        function(f)paste(e,f,paste0(c(4,2),"N"),sep="\n"))))))
    x
  }
  add <- function(p,kind,origin,o2,pm,beta,variable=NA_character_) {
    plots[[length(plots)+1L]] <<- p
    index[[length(index)+1L]] <<- data.frame(page=length(plots),kind=kind,origin=origin,
      O2_pct=o2,p_misseg=pm,beta=beta,variable=variable)
  }
  for(origin in names(groups))for(o2 in c(.5,5))for(pm in c(.1,.3)) {
    keep <- groups[[origin]];labels <- captions[[origin]]
    for(beta in c(0,.5,1)) {
      zz <- data.table::copy(z[O2_pct==o2&p_misseg==pm&q==beta&group%in%keep])
      zz <- panels(zz);zz$group <- factor(zz$group,levels=keep,labels=labels)
      zz$metric <- factor(zz$metric,levels=names(metric_labels),labels=metric_labels)
      data.table::setorderv(zz,c("panel","group","metric","day"))
      # Explicit separate paths prevent a ribbon or line from bridging an
      # invalid interval, including when only one protocol route is missing.
      zz[,validity_segment:=data.table::rleid(is.finite(median)&is.finite(q25)&is.finite(q75)),
        by=.(panel,group,metric)]
      p <- ggplot2::ggplot(zz,ggplot2::aes(day,median,color=group,fill=group))+
        ggplot2::geom_ribbon(ggplot2::aes(ymin=q25,ymax=q75,group=interaction(group,validity_segment)),alpha=.12,color=NA,na.rm=TRUE)+
        ggplot2::geom_line(ggplot2::aes(group=interaction(group,validity_segment)),linewidth=.4,na.rm=TRUE)+ggplot2::facet_grid(metric~panel,scales="free_y")+
        ggplot2::scale_color_manual(values=colors,drop=FALSE)+ggplot2::scale_fill_manual(values=colors,drop=FALSE)+
        ggplot2::scale_x_continuous(breaks=c(0,500,1000),limits=c(0,1000))+
        ggplot2::labs(title=paste0(if(beta==0)"Unblocked origin accounting" else "Blocked origin comparison",
          " | ",origin," | O2 = ",o2,"% | p_misseg = ",pm," | β = ",beta),x="Time (days)",y=NULL)+
        f78_plot_theme()+ggplot2::theme(strip.text=ggplot2::element_text(size=9),axis.text=ggplot2::element_text(size=8))
      add(p,"time",origin,o2,pm,beta)
    }
    for(variable_name in c("ploidy","division_rate","low_birth_rate")) {
      dd <- data.table::copy(d[O2_pct==o2&p_misseg==pm&q==0&group%in%keep&variable==variable_name])
      dd <- panels(dd);dd$group <- factor(dd$group,levels=keep,labels=labels)
      dd$snapshot <- factor(dd$day,levels=c(0,10,100,500,1000),labels=paste0(c(0,10,100,500,1000)," days"))
      xlab <- switch(variable_name,ploidy="Cell ploidy (N)",division_rate="Division rate (1/day)",
        low_birth_rate="Expected viable low-daughter output (cell/day)")
      p <- ggplot2::ggplot(dd,ggplot2::aes(value,cdf,color=group,fill=group))+
        ggplot2::geom_ribbon(ggplot2::aes(ymin=q25,ymax=q75),alpha=.12,color=NA,na.rm=TRUE)+
        ggplot2::geom_step(linewidth=.45,na.rm=TRUE)+ggplot2::facet_grid(snapshot~panel)+
        ggplot2::scale_color_manual(values=colors,drop=FALSE)+ggplot2::scale_fill_manual(values=colors,drop=FALSE)+
        ggplot2::scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
        ggplot2::labs(title=paste0("Unblocked within-group distributions | ",origin," | O2 = ",o2,
          "% | p_misseg = ",pm," | β = 0"),x=xlab,y="Cell-weighted cumulative fraction")+
        f78_plot_theme()+ggplot2::theme(strip.text=ggplot2::element_text(size=9),axis.text=ggplot2::element_text(size=8))
      add(p,"distribution",origin,o2,pm,0,variable_name)
    }
  }
  list(plots=plots,index=data.table::rbindlist(index))
}
f7_draw_origin_populations <- function(ctx,row) {
  a <- ctx$contexts[[row$key]]
  if(is.null(ctx$population_validity_table))stop("Conditional-mean reporting layer was not authenticated")
  pages <- f7_population_pages(data.table::fread(ctx$population_validity_table),
    f78ro_table(ctx,row$key,"low_population_distribution_summary.tsv"))
  stopifnot(length(pages$plots)==row$pages)
  f78ro_write(pages$index,file.path(ctx$audit,"supp_fig7_8_page_manifest.tsv"))
  if("pdf"%in%ctx$formats) {
    path <- file.path(ctx$stage,"Figures/supp_fig7_8.pdf")
    grDevices::cairo_pdf(path,width=24,height=13,onefile=TRUE)
    tryCatch(for(p in pages$plots)print(p),finally=grDevices::dev.off())
    file.copy(path,file.path(ctx$stage,"manuscript/Figures/supp_fig7_8.pdf"),overwrite=TRUE)
    ctx$exports[[basename(path)]] <- data.frame(stem=row$stem,format="pdf",width_inches=24,
      height_inches=13,dpi=NA_real_,sha256=f78ro_hash(path))
  }
  if("png"%in%ctx$formats) {
    # PNG is the first unblocked/direct page; the PDF contains all 48 labeled
    # pages. Do not compress 1,920 subplots into an unreadable contact sheet.
    old <- ctx$formats;ctx$formats <- "png"
    tryCatch(f78ro_publish(a,pages$plots[[1]],24,13,row$stem,dpi=180),finally={ctx$formats <- old})
  }
}
