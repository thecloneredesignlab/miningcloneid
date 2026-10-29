# Unified Figure7/8 rendering; never prepares, simulates, or reads archives.
.f78ro_directory <- dirname(normalizePath(sys.frame(1)$ofile,mustWork=TRUE))
.f78ro_names <- c("figure78_paths.R","figure78_plot_common.R","figure7_plotting.R","figure7_population_validity.R","figure7_population_plotting.R","figure8_plotting.R",
  "figure78_render_io.R","figure78_provenance.R","figure78_render_plots.R")
for(.f78ro_name in .f78ro_names)
  source(file.path(.f78ro_directory,.f78ro_name),local=globalenv())
.f78ro_loaded_files <- c(file.path(.f78ro_directory,"figure78_render_only.R"),
  file.path(.f78ro_directory,.f78ro_names))
options(stringsAsFactors=FALSE,warn=1)

f78ro_registry <- function() {
  z <- data.frame(stem=c("assembled_fig7",paste0("supp_fig7_",1:7),"figure7_ver2",
    paste0("supp_fig7_",1:7,"_ver2"),"assembled_fig8",paste0("supp_fig8_",1:9)),
    figure=c(rep(7L,16),rep(8L,10)),
    analysis=c(rep("instantaneous",8),rep("history",8),rep("current",10)))
  z$key <- c("instantaneous7","lineage7","lineage7",rep("instantaneous7",5),rep("lineage7",8),
    rep("main8",3),rep("family_curves8",2),"main8",rep("time8",2),rep("family_time8",2))
  z$pages <- ifelse(z$stem%in%c("supp_fig7_3","supp_fig7_4","supp_fig7_3_ver2","supp_fig7_4_ver2"),8L,1L)
  rbind(z,data.frame(stem="supp_fig7_8",figure=7L,analysis="instantaneous",key="instantaneous7",pages=48L))
}
f78ro_contexts <- function(ctx) {
  make <- function(figure,key,group) {
    config <- jsonlite::fromJSON(file.path(ctx$workspace,"Code/config",paste0("figure",figure,".json")))
    list(figure=figure,key=key,group=group,workspace=ctx$workspace,
      code=file.path(ctx$workspace,"Code/Figures"),cfg=config,
      run=file.path(ctx$data_root,paste0("Figure",figure)),
      revision=list(low_fraction_tolerance=1e-14),family_config=config$family_time,render_ctx=ctx)
  }
  list(instantaneous7=make(7L,"instantaneous7","instantaneous"),
    lineage7=make(7L,"lineage7","lineage"),main8=make(8L,"main8","main"),
    family_curves8=make(8L,"family_curves8","family_curves"),
    time8=make(8L,"time8","time_snapshots"),family_time8=make(8L,"family_time8","family_time_snapshots"))
}
f78ro_required_tables <- function(row) {
  stem <- row$stem
  if(row$figure==7L) {
    if(stem=="supp_fig7_8")return(c("low_population_time_summary.tsv","low_population_distribution_summary.tsv"))
    if(stem%in%c("supp_fig7_1","supp_fig7_2"))
      return(c("source_daughter_flux_binned.tsv","winner_low_population_budget.tsv"))
    if(stem%in%c("supp_fig7_1_ver2","supp_fig7_2_ver2"))
      return(c("history_source_daughter_flux_binned.tsv","history_low_population_budget.tsv"))
    return(paste0(if(row$analysis=="history")"history" else "instant","_time_summary.tsv"))
  }
  switch(stem,assembled_fig8=c("empirical_tissue_ploidy_summary.tsv","curve_summary.tsv","author_figure5_points.tsv",
      "Figure8C_curve_summary.tsv","Figure8D_curve_summary.tsv"),
    supp_fig8_1=c("empirical_tissue_ploidy_summary.tsv","curve_summary.tsv"),
    supp_fig8_2=c("empirical_tissue_ploidy_summary.tsv","curve_summary.tsv"),
    supp_fig8_3=c("fitted_family_curve_summary.tsv","supp_fig8_3_curve_summary.tsv","supp_fig8_4_curve_summary.tsv"),
    supp_fig8_4=c("fitted_family_curve_summary.tsv","supp_fig8_3_curve_summary.tsv","supp_fig8_4_curve_summary.tsv"),
    supp_fig8_5="supp_fig8_5_curve_summary.tsv",paste0(stem,"_curve_summary.tsv"))
}
f78_render_cli <- function(script,default_scope="all",default_products=NULL,validate_only=FALSE) {
  args <- f78_args()
  known <- c("action","workspace","task_index","products","analysis","version","figure","stem","output_root","data_root",
    "audit_root","formats","list_products","preflight_only","help")
  if(length(setdiff(names(args),known)))stop("Unsupported render option: ",paste(setdiff(names(args),known),collapse=", "))
  if(!identical(args$action,"all"))stop("This is a retained-data entry point; computation actions are not supported")
  registry <- f78ro_registry()
  suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(data.table)})
  data.table::setDTthreads(2L)
  assign("JOINT_FAMILY_LEVELS",c("C01","C02"),globalenv())
  if(identical(args$help,"yes")) {
    cat("Unified retained-data renderer: draw_Figure7.R / draw_Figure8.R\n",
      "--products=stem,... --analysis=instantaneous|history|all --output-root=DIR --formats=pdf,png\n",
      "--workspace=DIR --data-root=DIR --preflight-only=yes --list-products=yes\n",
      "Numerical workers are not run. Legacy --version=instant|history is a compatibility alias only.\n")
    return(invisible(NULL))
  }
  if(identical(args$list_products,"yes")) {print(registry,row.names=FALSE);return(invisible(registry))}
  workspace <- normalizePath(if(nzchar(args$workspace))args$workspace else file.path(dirname(script),"../.."),mustWork=TRUE)
  f78ro_no_archive(workspace)
  scope <- default_scope
  if(scope=="repair") {
    figure <- args$figure %||% Sys.getenv("FIGURE78_REPAIR_FIGURE",unset="8")
    if(!figure%in%c("7","8"))stop("--figure must be 7 or 8")
    scope <- paste0("figure",figure)
  }
  if(scope=="versions7")scope <- "all7"
  analysis <- args$analysis
  if(!is.null(args$version)) {
    legacy <- switch(args$version,instant="instantaneous",history="history",stop("--version must be instant or history"))
    if(!is.null(analysis)&&analysis!=legacy)stop("Conflicting --analysis and --version options")
    analysis <- legacy
  }
  if(!is.null(analysis)&&!analysis%in%c("instantaneous","history","all"))stop("Unknown --analysis")
  products <- default_products
  if(!is.null(args$products))products <- strsplit(args$products,",",fixed=TRUE)[[1]]
  if(!is.null(args$stem))products <- args$stem
  if(is.null(products)) {
    products <- switch(scope,all=registry$stem,all7=registry$stem[registry$figure==7L],
      figure7=registry$stem[registry$figure==7L&registry$analysis=="instantaneous"],
      figure8=registry$stem[registry$figure==8L],stop("Unknown render scope: ",scope))
    if(scope=="figure7"&&!is.null(analysis))
      products <- registry$stem[registry$figure==7L&(analysis=="all"|registry$analysis==analysis)]
    if(scope=="all7"&&!is.null(analysis)&&analysis!="all")
      products <- intersect(products,registry$stem[registry$analysis==analysis])
  } else if(!is.null(analysis)&&analysis!="all") {
    if(any(!products%in%registry$stem[registry$analysis==analysis]))stop("Product selection conflicts with --analysis")
  }
  if(!length(products)||anyDuplicated(products)||any(!products%in%registry$stem))
    stop("Unknown, empty, or duplicated product selection: ",paste(products,collapse=", "))
  ctx <- new.env(parent=emptyenv());ctx$workspace <- workspace;ctx$registry <- registry;ctx$products <- products
  ctx$data_root <- f78_current_path(args$data_root %||% "data/Figures",workspace,mustWork=TRUE)
  ctx$model_sources <- f78_locked_model_manifest(workspace)
  f78ro_no_archive(ctx$data_root)
  ctx$output_root <- f78ro_absolute(args$output_root %||% workspace)
  f78ro_no_archive(ctx$output_root)
  forbidden <- c(ctx$data_root,file.path(workspace,c("data","Code","manuscript/figure_legends")))
  if(any(startsWith(paste0(ctx$output_root,"/"),paste0(forbidden,"/"))))
    stop("Output root must not be scientific data, code, or methods")
  ctx$formats <- strsplit(args$formats %||% "pdf,png",",",fixed=TRUE)[[1]]
  if(!length(ctx$formats)||anyDuplicated(ctx$formats)||any(!ctx$formats%in%c("pdf","png")))stop("Formats must be pdf and/or png")
  config <- file.path(workspace,"Code/config")
  manifest <- file.path(config,"manifests")
  selected_figures <- unique(registry$figure[match(products,registry$stem)])
  # Figure7-only Linux runs on headless compute nodes cannot use Xlib PNG.
  # Choose the already available Cairo device, without changing plot builders,
  # numerical inputs, local macOS output, or the Figure8/historical products.
  if(identical(selected_figures,7L)&&identical(unname(Sys.info()["sysname"]),"Linux")&&
    !nzchar(Sys.getenv("DISPLAY"))&&isTRUE(capabilities("cairo")))options(bitmapType="cairo")
  selected_rows <- registry[match(products,registry$stem),]
  retained_history7 <- any(selected_rows$figure==7L&selected_rows$analysis=="history")
  ctx$all_current <- !retained_history7&&
    all(file.exists(file.path(ctx$data_root,paste0("Figure",selected_figures),"prepared/analysis_bundle.rds")))
  if(ctx$all_current) {
    # A cold-start generated dataset needs none of the historical consolidation
    # manifests or frozen old workers. Explicitly selected historical products
    # still need their ORIGINAL bindings, even alongside a new current bundle.
    ctx$inputs <- data.table::data.table(figure=character(),group=character(),logical_name=character(),
      relative_path=character(),sha256=character())
    ctx$migration <- ctx$frozen_sources <- ctx$bindings <- data.table::data.table()
    ctx$layout <- list()
  } else {
    ctx$inputs <- data.table::fread(file.path(manifest,"figure78_current_plot_inputs.tsv"))
    ctx$migration <- data.table::fread(file.path(manifest,"figure78_data_migration.tsv"))
    ctx$frozen_sources <- data.table::fread(file.path(manifest,"figure78_frozen_sources.tsv"))
    ctx$bindings <- data.table::fread(file.path(manifest,"figure78_bundle_bindings.tsv"))
    ctx$layout <- jsonlite::fromJSON(file.path(manifest,"figure78_layout.json"))
  }
  f78ro_current_inputs(ctx)
  if(anyDuplicated(ctx$inputs$relative_path)||anyDuplicated(paste(ctx$inputs$figure,ctx$inputs$group,ctx$inputs$logical_name)))
    stop("Duplicate canonical input registration")
  ctx$checked <- data.frame(path=character(),relative_path=character(),sha256=character())
  ctx$authenticated <- ctx$model_sources[,c("path","sha256")]
  ctx$contexts <- f78ro_contexts(ctx);ctx$tables <- new.env(parent=emptyenv());ctx$exports <- list()
  ctx$allowed_stems <- unique(c(products,unlist(lapply(products,f78ro_component_stems))))
  ctx$allowed_files <- as.vector(outer(ctx$allowed_stems,paste0(".",ctx$formats),paste0))
  f78ro_preflight(ctx)
  audit_root <- f78ro_absolute(args$audit_root %||% file.path(ctx$output_root,"audit/figure78_render_only"))
  f78ro_no_archive(audit_root)
  if(any(startsWith(paste0(audit_root,"/"),paste0(forbidden,"/"))))stop("Audit root must not be scientific data, code, or methods")
  ctx$audit <- file.path(audit_root,paste0(format(Sys.time(),"%Y%m%d_%H%M%OS6"),"_",Sys.getpid()))
  if(dir.exists(ctx$audit))stop("Audit run already exists")
  dir.create(ctx$audit,recursive=TRUE,showWarnings=FALSE)
  ctx$stage <- file.path(ctx$audit,"staging")
  for(directory in c("Figures","manuscript/Figures"))dir.create(file.path(ctx$stage,directory),recursive=TRUE,showWarnings=FALSE)
  ctx$protected <- f78ro_image_snapshot(ctx)
  render_files <- unique(c(script,.f78ro_loaded_files,
    file.path(workspace,"Code/Figures/merge_Figure7_matrix_pdfs.py"),
    file.path(config,c("figure7.json","figure8.json")),
    file.path(manifest,c(if(!ctx$all_current)c("figure78_current_plot_inputs.tsv","figure78_data_migration.tsv",
      "figure78_frozen_sources.tsv","figure78_bundle_bindings.tsv","figure78_layout.json"),
      "figure7_display_panels.tsv","figure78_model_sources.tsv"))))
  ctx$render_sources <- data.frame(path=render_files,sha256=f78ro_hash(render_files))
  f78ro_write(ctx$protected,file.path(ctx$audit,"protected_before.tsv"))
  f78ro_write(ctx$render_sources,file.path(ctx$audit,"render_code_manifest.tsv"))
  f78ro_write(registry[match(products,registry$stem),],file.path(ctx$audit,"requested_products.tsv"))
  f78ro_write(ctx$checked,file.path(ctx$audit,"inputs_verified.tsv"))
  if(validate_only||identical(args$preflight_only,"yes")) {
    message("PASS: preflight only. Audit: ",ctx$audit);return(invisible(ctx$audit))
  }
  cache <- file.path(ctx$audit,"graphics_cache");dir.create(cache,showWarnings=FALSE)
  Sys.setenv(XDG_CACHE_HOME=cache)
  grDevices::pdf(NULL);null_device <- grDevices::dev.cur()
  on.exit(try(grDevices::dev.off(null_device),silent=TRUE),add=TRUE)
  options(figure7.proliferation=TRUE)
  for(stem in products)f78ro_plot(ctx,registry[match(stem,registry$stem),])
  f78ro_finish(ctx)
  invisible(ctx$audit)
}
