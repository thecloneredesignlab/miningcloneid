# Same approved geometry/encodings; only fitted-p condition selection differs.
source(file.path(.f7fit_dir,"figure78_plot_common.R"),local=globalenv())
source(file.path(.f7fit_dir,"figure7_plotting.R"),local=globalenv())

f7fit_panel <- function(z,oxygen,title,metric) {
  canvas <- cowplot::ggdraw();cols <- c(C01="#C99700",C02="#6A3D9A")
  header <- function(txt,color,angle=0,size=9) {
    cowplot::ggdraw()+cowplot::draw_label(txt,color="white",angle=angle,
      fontfamily="Helvetica",fontface="bold",size=size)+
      ggplot2::theme(plot.background=ggplot2::element_rect(fill=color,colour="#BEBEBE",linewidth=.2))
  }
  xpos <- c(.18,.585);ww <- .395;hh <- .205;ytop <- c(.72,.505,.265,.05)
  canvas <- canvas+cowplot::draw_label(title,x=.18,y=.985,hjust=0,fontface="bold",size=13)
  for(ci in 1:2)canvas <- canvas+cowplot::draw_plot(header(names(cols)[ci],cols[ci]),xpos[ci],.938,ww,.035)
  context <- c("in vivo","in vivo","in vitro","in vitro");initial <- c(4,2,4,2)
  for(r in 1:4) {
    tag <- cowplot::ggdraw()+cowplot::draw_label(paste0(initial[r],"N"),angle=90,size=9)+
      ggplot2::theme(plot.background=ggplot2::element_rect(fill="#F2F2F2",colour="#BEBEBE",linewidth=.2))
    canvas <- canvas+cowplot::draw_plot(tag,.11,ytop[r],.028,hh)
    for(ci in 1:2) {
      zz <- z[z$model_context==context[r]&z$initial_ploidy==initial[r]&
        z$pair_label==names(cols)[ci]&z$O2_pct==oxygen&z$metric==metric,]
      stopifnot(nrow(zz)==3*1001,length(unique(zz$p_misseg))==1L)
      plot <- f7_line_plot(zz,metric,axes=FALSE,y_limits=c(0,100))+
        ggplot2::theme(legend.position="none",plot.margin=ggplot2::margin(3,3,3,3))
      if(r!=4)plot <- plot+ggplot2::theme(axis.text.x=ggplot2::element_blank(),axis.ticks.x=ggplot2::element_blank())
      if(ci!=1)plot <- plot+ggplot2::theme(axis.text.y=ggplot2::element_blank(),axis.ticks.y=ggplot2::element_blank())
      canvas <- canvas+cowplot::draw_plot(plot,xpos[ci],ytop[r],ww,hh)
    }
  }
  for(rr in c(1,3))canvas <- canvas+cowplot::draw_plot(header(context[rr],f78_context_colors[[context[rr]]],90),
    .07,ytop[rr+1],.028,2*hh+.01)
  canvas+cowplot::draw_label(f7_metric_label(metric),x=.025,y=.49,angle=90,size=10)+
    cowplot::draw_label("Time (days)",x=.59,y=.015,size=10)
}
f7fit_population_pages <- function(z,d) {
  groups <- list(direct=c("all_low","direct_first_generation","other_low"),
    historical=c("all_low","historical_lineage","never_HL_lineage"))
  captions <- list(direct=c("All low","First-generation H-to-L","Other low (includes descendants)"),
    historical=c("All low","H-to-L lineage (including descendants)","Never-H-to-L lineage"))
  ml <- c(cell_share_percent="Current low-cell share (%)",low_birth_share_percent="Low-daughter production share (%)",
    mean_low_birth_rate="Low-daughter output (cell/day)",mean_division_rate="Division rate (1/day)",mean_ploidy="Mean ploidy (N)")
  colors <- c("#222222","#D55E00","#0072B2");plots <- list();index <- list()
  panels <- function(x) {
    x$panel <- factor(paste(x$model_context,x$pair_label,paste0(x$initial_ploidy,"N"),sep="\n"),
      levels=unlist(lapply(c("in vivo","in vitro"),function(e)unlist(lapply(c("C01","C02"),
        function(f)paste(e,f,paste0(c(4,2),"N"),sep="\n"))))));x
  }
  add <- function(p,kind,origin,o2,beta,variable=NA_character_) {
    plots[[length(plots)+1L]] <<- p
    index[[length(index)+1L]] <<- data.table::data.table(page=length(plots),kind=kind,
      origin=origin,O2_pct=o2,condition_mode="fitted_p_misseg",beta=beta,variable=variable)
  }
  for(origin in names(groups))for(o2 in c(.5,5)) {
    keep <- groups[[origin]];labels <- captions[[origin]]
    for(beta in c(0,.5,1)) {
      zz <- data.table::copy(z[O2_pct==o2&q==beta&group%in%keep])
      zz <- panels(zz);zz$group <- factor(zz$group,levels=keep,labels=labels)
      zz$metric <- factor(zz$metric,levels=names(ml),labels=ml)
      data.table::setorderv(zz,c("panel","group","metric","day"))
      zz[,validity_segment:=data.table::rleid(is.finite(median)&is.finite(q25)&is.finite(q75)),by=.(panel,group,metric)]
      p <- ggplot2::ggplot(zz,ggplot2::aes(day,median,color=group,fill=group))+
        ggplot2::geom_ribbon(ggplot2::aes(ymin=q25,ymax=q75,group=interaction(group,validity_segment)),alpha=.12,color=NA,na.rm=TRUE)+
        ggplot2::geom_line(ggplot2::aes(group=interaction(group,validity_segment)),linewidth=.4,na.rm=TRUE)+
        ggplot2::facet_grid(metric~panel,scales="free_y")+
        ggplot2::scale_color_manual(values=colors,drop=FALSE)+ggplot2::scale_fill_manual(values=colors,drop=FALSE)+
        ggplot2::scale_x_continuous(breaks=c(0,500,1000),limits=c(0,1000))+
        ggplot2::labs(title=paste0(if(beta==0)"Unblocked origin accounting" else "Blocked origin comparison",
          " | ",origin," | O2 = ",o2,"% | fitted p_misseg | β = ",beta),x="Time (days)",y=NULL)+
        f78_plot_theme()+ggplot2::theme(strip.text=ggplot2::element_text(size=9),axis.text=ggplot2::element_text(size=8))
      add(p,"time",origin,o2,beta)
    }
    for(v in c("ploidy","division_rate","low_birth_rate")) {
      dd <- data.table::copy(d[O2_pct==o2&q==0&group%in%keep&variable==v])
      dd <- panels(dd);dd$group <- factor(dd$group,levels=keep,labels=labels)
      dd$snapshot <- factor(dd$day,levels=c(0,10,100,500,1000),labels=paste0(c(0,10,100,500,1000)," days"))
      lab <- switch(v,ploidy="Cell ploidy (N)",division_rate="Division rate (1/day)",low_birth_rate="Expected viable low-daughter output (cell/day)")
      p <- ggplot2::ggplot(dd,ggplot2::aes(value,cdf,color=group,fill=group))+
        ggplot2::geom_ribbon(ggplot2::aes(ymin=q25,ymax=q75),alpha=.12,color=NA,na.rm=TRUE)+
        ggplot2::geom_step(linewidth=.45,na.rm=TRUE)+ggplot2::facet_grid(snapshot~panel)+
        ggplot2::scale_color_manual(values=colors,drop=FALSE)+ggplot2::scale_fill_manual(values=colors,drop=FALSE)+
        ggplot2::scale_y_continuous(limits=c(0,1),breaks=c(0,.5,1))+
        ggplot2::labs(title=paste0("Unblocked within-group distributions | ",origin," | O2 = ",o2,
          "% | fitted p_misseg | β = 0"),x=lab,y="Cell-weighted cumulative fraction")+
        f78_plot_theme()+ggplot2::theme(strip.text=ggplot2::element_text(size=9),axis.text=ggplot2::element_text(size=8))
      add(p,"distribution",origin,o2,0,v)
    }
  }
  list(plots=plots,index=data.table::rbindlist(index))
}
f7fit_plot <- function(a) {
  suppressPackageStartupMessages({library(ggplot2);library(cowplot)})
  b <- f7fit_bundle(a)
  done <- jsonlite::fromJSON(file.path(a$run,"provenance/completion.json"))
  stopifnot(done$status=="validated_complete",done$fingerprint==b$fingerprint,
    done$bundle_sha256==fp_sha(file.path(a$run,"prepared/analysis_bundle.rds")),
    done$manifest_sha256==fp_sha(file.path(a$run,"provenance/validated_files.tsv")))
  m <- data.table::fread(file.path(a$run,"provenance/validated_files.tsv"))
  stopifnot(identical(fp_sha(file.path(a$run,m$relative_path)),m$sha256))
  stems <- c(a$policy$figure_stem,a$policy$supplement_stem)
  chosen <- if(is.null(a$args$products))stems else strsplit(a$args$products,",",fixed=TRUE)[[1]]
  stopifnot(length(chosen)>0,all(chosen%in%stems))
  canonical <- file.path(a$workspace,"data/Figures/Figure7/fitted_p_misseg")
  root <- if(is.null(a$args$output_root)&&identical(a$run,canonical))a$workspace else a$output_root
  dir.create(file.path(root,"Figures"),recursive=TRUE,showWarnings=FALSE)
  dir.create(file.path(root,"manuscript/Figures"),recursive=TRUE,showWarnings=FALSE)
  options(figure7.proliferation=TRUE,bitmapType="cairo")
  assign("f7_metric_label",f7_metric_labels("instantaneous"),globalenv())
  exports <- list()
  publish <- function(plots,stem,width,height,dpi) {
    pdf <- file.path(root,"Figures",paste0(stem,".pdf"))
    png <- file.path(root,"Figures",paste0(stem,".png"))
    grDevices::cairo_pdf(pdf,width=width,height=height,onefile=TRUE)
    tryCatch(for(p in plots)print(p),finally=grDevices::dev.off())
    ggplot2::ggsave(png,plots[[1]],width=width,height=height,units="in",dpi=dpi,bg="white",limitsize=FALSE)
    for(path in c(pdf,png))stopifnot(file.copy(path,file.path(root,"manuscript/Figures",basename(path)),overwrite=TRUE))
    exports[[length(exports)+1L]] <<- data.table::data.table(stem=stem,pages=length(plots),
      pdf_sha256=fp_sha(pdf),png_sha256=fp_sha(png),width=width,height=height)
  }
  if(stems[1]%in%chosen) {
    z <- as.data.frame(data.table::fread(file.path(a$run,"tables/instant_time_summary.tsv")))
    A <- f7fit_panel(z,.5,"A  High-parent source; O2 = 0.5%; fitted p_misseg","flux_percent")
    B <- f7fit_panel(z,.5,"B  Low-parent source; O2 = 0.5%; fitted p_misseg","growth_percent")
    legend <- cowplot::get_legend(f7_line_plot(z[z$task_id==z$task_id[1],])+ggplot2::theme(legend.position="bottom"))
    main <- cowplot::plot_grid(cowplot::plot_grid(A,B,nrow=1),legend,ncol=1,rel_heights=c(1,.04))
    publish(list(main),stems[1],14,11,300)
  }
  if(stems[2]%in%chosen) {
    pages <- f7fit_population_pages(data.table::fread(file.path(a$run,"postprocess/low_population_time_summary.tsv")),
      data.table::fread(file.path(a$run,"tables/low_population_distribution_summary.tsv")))
    stopifnot(length(pages$plots)==24L)
    publish(pages$plots,stems[2],24,13,180)
    fp_write(pages$index,file.path(root,"supp_fig7_8_fitted_p_misseg_page_manifest.tsv"))
  }
  fp_write(data.table::rbindlist(exports),file.path(root,"figure7_fitted_exports.tsv"))
  sources <- file.path(a$code,c("draw_Figure7_fitted_p_misseg.R",
    "util/analysis/figure7_fitted_plotting.R","util/analysis/figure7_plotting.R",
    "util/analysis/figure78_plot_common.R"))
  fp_write(data.table::data.table(path=sources,sha256=fp_sha(sources)),
    file.path(root,"figure7_fitted_render_sources.tsv"))
  jsonlite::write_json(list(status="rendered",fingerprint=b$fingerprint,products=chosen,
    input_manifest_sha256=done$manifest_sha256,render_sources_sha256=fp_sha(file.path(root,"figure7_fitted_render_sources.tsv")),
    exports_sha256=fp_sha(file.path(root,"figure7_fitted_exports.tsv")),simulation_calls=0L),
    file.path(root,"figure7_fitted_render_receipt.json"),auto_unbox=TRUE,pretty=TRUE)
  cat("FITTED_PM_RENDERED",paste(chosen,collapse=","),"\n")
}
