# Authenticate immutable numerical results after a path-only consolidation.
# Original receipts and RDS bundles are not edited or re-signed.
f78ro_path <- function(a,files) f78ro_data(a$render_ctx,a$key,files)

f78ro_data <- function(ctx,key,files) {
  a <- ctx$contexts[[key]]
  if(is.null(a)) stop("Unknown dataset group: ",key)
  if(paste0("Figure",a$figure)%in%names(ctx$current)) {
    historical <- a$figure==7L&&key=="lineage7"&&
      (isTRUE(ctx$reading_historical7)||all(startsWith(files,"history_")))
    rows <- ctx$inputs[ctx$inputs$figure==paste0("Figure",a$figure)&
      ctx$inputs$group==if(historical)"lineage" else "current",]
    ix <- match(files,rows$logical_name)
    if(anyNA(ix))stop("Current pipeline table is not registered: ",paste(files[is.na(ix)],collapse=", "))
    paths <- file.path(ctx$data_root,rows$relative_path[ix]);f78ro_input(ctx,paths);return(paths)
  }
  rows <- ctx$inputs[ctx$inputs$figure==paste0("Figure",a$figure)&ctx$inputs$group==a$group,]
  ix <- match(files,rows$logical_name)
  if(anyNA(ix)) stop("Unregistered dataset slot: ",key,"/",paste(files[is.na(ix)],collapse=", "))
  paths <- file.path(ctx$data_root,rows$relative_path[ix])
  f78ro_input(ctx,paths)
  paths
}

# Newly simulated data authenticate against their own current numerical bundle,
# never against a retired worker or a frozen copy of old figure code.
f78ro_current_inputs <- function(ctx) {
  ctx$current <- list()
  selected <- unique(ctx$registry$figure[match(ctx$products,ctx$registry$stem)])
  for(figure in selected) {
    tag <- paste0("Figure",figure);root <- file.path(ctx$data_root,tag)
    bundle <- file.path(root,"prepared/analysis_bundle.rds")
    receipt <- file.path(root,"provenance/completion.tsv")
    if(!file.exists(bundle))next
    f78ro_no_archive(c(bundle,receipt))
    if(!file.exists(receipt))stop("Current pipeline dataset has not completed validation: ",root)
    done <- data.table::fread(receipt)
    if(nrow(done)!=1L||done$status!="validated_complete"||done$profile!="full"||!isTRUE(done$production_ensemble))
      stop("Only validated full-ensemble pipeline results may be rendered as publication figures: ",root)
    table_manifest <- file.path(root,"provenance/plot_inputs.tsv")
    validation <- file.path(root,"provenance/numerical_validation.tsv")
    observed <- f78ro_hash(c(bundle,validation,table_manifest))
    if(!identical(observed,as.character(unlist(done[,.(bundle_sha256,validation_sha256,plot_manifest_sha256)],use.names=FALSE))))
      stop("Current pipeline receipt mismatch: ",root)
    tables <- data.table::fread(table_manifest)
    tables[,`:=`(figure=tag,group="current")]
    meta_paths <- c(bundle,receipt,table_manifest,validation,file.path(root,"provenance/current_code_manifest.tsv"))
    meta <- data.table::data.table(relative_path=substring(meta_paths,nchar(ctx$data_root)+2L),
      logical_name=basename(meta_paths),sha256=f78ro_hash(meta_paths),figure=tag,group="current")
    # A new current Figure7 generation does not re-sign or replace its preserved
    # historical tables/receipts. Keep their original bindings when present;
    # an isolated current-only campaign needs none of these legacy files.
    preserved <- NULL
    if(figure==7L) {
      preserved <- ctx$inputs[ctx$inputs$figure==tag&
        !ctx$inputs$relative_path%in%c(tables$relative_path,meta$relative_path),]
      preserved <- preserved[file.exists(file.path(ctx$data_root,preserved$relative_path)),]
    }
    ctx$inputs <- data.table::rbindlist(list(ctx$inputs[ctx$inputs$figure!=tag,],preserved,tables,meta),use.names=TRUE,fill=TRUE)
    ctx$current[[tag]] <- list(bundle=bundle,validation=validation,receipt=done)
  }
  invisible(ctx)
}
f78ro_current_preflight <- function(ctx,figure) {
  tag <- paste0("Figure",figure);record <- ctx$current[[tag]]
  b <- readRDS(record$bundle)
  if(!identical(b$fingerprint,record$receipt$fingerprint)||b$profile!="full"||b$figure!=figure)
    stop("Current pipeline bundle identity mismatch")
  checks <- data.table::fread(record$validation)
  if(!nrow(checks)||!all(checks$passed))stop("Current numerical validation is not complete")
  sources <- b$code_manifest
  # A verified HPC-generated bundle retains its original immutable bytes and
  # fingerprint. Resolve its Code-relative source slots in this checkout rather
  # than requiring the HPC absolute mount to exist on a laptop. Every current
  # source must still match the originally authenticated SHA256 exactly.
  source_slots <- sub("^.*(/Code/)","Code/",sources$path)
  if(any(!startsWith(source_slots,"Code/"))||any(grepl("(^|/)\\.\\.(/|$)",source_slots)))
    stop("Unrecognized current numerical source slot in bundle")
  sources$path <- file.path(ctx$workspace,source_slots)
  if(!identical(f78ro_hash(sources$path),sources$sha256))stop("Current numerical code no longer matches this generated dataset")
  ctx$authenticated <- unique(rbind(ctx$authenticated,as.data.frame(sources[,.(path,sha256)])))
  message("PASS: Figure",figure," current full simulation/data-generation provenance authenticated")
}

f78ro_receipt_paths <- function(ctx,paths,figure) {
  origin <- ctx$layout$original_workspace
  wanted_figure <- paste0("Figure",figure)
  source <- ctx$frozen_sources[ctx$frozen_sources$figure==wanted_figure,]
  vapply(paths,function(path) {
    relative <- if(startsWith(path,paste0(origin,"/"))) substring(path,nchar(origin)+2L) else path
    i <- match(relative,ctx$migration$old_relative_path)
    if(!is.na(i)) return(file.path(ctx$data_root,ctx$migration$relative_path[i]))
    absolute <- if(startsWith(path,"/"))path else file.path(origin,path)
    j <- match(absolute,source$original_path)
    if(is.na(j)) stop("Unmapped original source receipt: ",path)
    file.path(ctx$data_root,source$relative_path[j])
  },"")
}

f78ro_manifest_guard <- function(ctx,key,name) {
  m <- data.table::fread(f78ro_data(ctx,key,name))
  original <- if("relative_path"%in%names(m))m$relative_path else m$path
  paths <- unname(f78ro_receipt_paths(ctx,original,ctx$contexts[[key]]$figure))
  f78ro_input(ctx,paths)
  sha <- f78ro_hash(paths)
  if(!identical(sha,m$sha256)) stop("Frozen source/code manifest mismatch: ",key,"/",name)
  ctx$authenticated <- unique(rbind(ctx$authenticated,data.frame(path=paths,sha256=sha)))
}

f78ro_binding <- function(ctx,name) {
  b <- ctx$bindings[ctx$bindings$binding==name,]
  b <- b[order(b$order),]
  if(!nrow(b)||!identical(as.integer(b$order),seq.int(0L,nrow(b)-1L)))
    stop("Invalid frozen bundle binding: ",name)
  paths <- file.path(ctx$data_root,b$relative_path)
  f78ro_input(ctx,paths)
  observed <- f78ro_hash(paths)
  if(!identical(observed,b$sha256))stop("Frozen binding mismatch: ",name)
  ctx$authenticated <- unique(rbind(ctx$authenticated,data.frame(path=paths,sha256=observed)))
  paths
}

f78ro_preflight <- function(ctx) {
  selected <- ctx$registry[match(ctx$products,ctx$registry$stem),]
  for(i in seq_len(nrow(selected)))f78ro_data(ctx,selected$key[i],f78ro_required_tables(selected[i,]))
  figures <- paste0("Figure",unique(selected$figure))
  # Verify every retained scientific and frozen-source file belonging to a
  # selected figure, including records that are not read by that one panel.
  f78ro_input(ctx,file.path(ctx$data_root,ctx$inputs$relative_path[ctx$inputs$figure%in%figures]))
  for(figure in unique(selected$figure))if(paste0("Figure",figure)%in%names(ctx$current))f78ro_current_preflight(ctx,figure)
  if("supp_fig7_8"%in%ctx$products)f7pv_render_preflight(ctx)
  if("Figure7"%in%names(ctx$current)&&any(selected$figure==7L&selected$analysis=="history")) {
    # Authenticate historical products against their ORIGINAL numerical receipt,
    # never the new current campaign's fingerprint or grid certificates.
    ctx$reading_historical7 <- TRUE
    tryCatch({
      original <- ctx$inputs[ctx$inputs$figure=="Figure7"&ctx$inputs$group=="lineage"&
        ctx$inputs$logical_name=="completion.tsv",]
      if(nrow(original)) {
        f78ro_complete(ctx,"lineage7")
        f78ro_manifest_guard(ctx,"lineage7","code_manifest.tsv")
        old <- readRDS(f78ro_data(ctx,"lineage7","analysis_bundle.rds"))
        paths <- f78ro_binding(ctx,"figure7_code")
        if(!identical(old$code_sha,f78ro_hash(paths))||
          !identical(old$config,jsonlite::fromJSON(paths[1])))
          stop("Preserved Figure7 historical provenance mismatch")
      } else {
        # The previous deployed dataset explicitly retained these three older
        # summaries byte-for-byte, separately from its newly computed history.
        # Its obsolete replay directory is not required or restored. Authenticate
        # the preserved results against the ORIGINAL consolidation SHA manifest;
        # do not present them as newly simulated or as replay-code authenticated.
        wanted <- c("history_time_summary.tsv","history_source_daughter_flux_binned.tsv",
          "history_low_population_budget.tsv")
        retained <- ctx$inputs[ctx$inputs$figure=="Figure7"&ctx$inputs$group=="lineage"&
          ctx$inputs$logical_name%in%wanted,]
        if(nrow(retained)!=3L||!setequal(retained$logical_name,wanted))
          stop("Missing original retained historical SHA bindings")
        f78ro_data(ctx,"lineage7",wanted)
        certificate <- file.path(ctx$data_root,"Figure7/provenance/preserved_history.tsv")
        if(file.exists(certificate)) {
          f78ro_no_archive(certificate);preserved <- data.table::fread(certificate)
          i <- match(retained$logical_name,preserved$logical_name)
          if(nrow(preserved)!=3L||anyNA(i)||
            !identical(retained$sha256,as.character(preserved$preserved_sha256[i])))
            stop("Retained historical preservation certificate disagrees with original SHA bindings")
        }
        message("PASS: Figure7 historical summaries match original retained SHA256; not new simulations or replay-code authentication")
      }
    },finally={ctx$reading_historical7 <- FALSE})
  }
  if(7L%in%selected$figure&&!"Figure7"%in%names(ctx$current)) {
    f78ro_complete(ctx,"lineage7");f78ro_complete(ctx,"instantaneous7")
    f78ro_manifest_guard(ctx,"lineage7","code_manifest.tsv")
    f78ro_manifest_guard(ctx,"instantaneous7","code_manifest.tsv")
    f78ro_manifest_guard(ctx,"instantaneous7","source_manifest.tsv")
    b <- readRDS(f78ro_data(ctx,"lineage7","analysis_bundle.rds"))
    paths <- f78ro_binding(ctx,"figure7_code")
    if(!identical(b$code_sha,f78ro_hash(paths)))stop("Figure7 frozen replay-code authentication failed")
    if(!identical(b$config,jsonlite::fromJSON(paths[1])))stop("Figure7 original numerical configuration changed")
  }
  if(8L%in%selected$figure&&!"Figure8"%in%names(ctx$current)) {
    f78ro_complete(ctx,"main8");f78ro_complete(ctx,"time8")
    f78ro_manifest_guard(ctx,"main8","code_manifest.tsv")
    f78ro_manifest_guard(ctx,"time8","code_manifest.tsv")
    b <- readRDS(f78ro_data(ctx,"time8","analysis_bundle.rds"))
    paths <- f78ro_binding(ctx,"figure8_code")
    if(!identical(b$code_sha,f78ro_hash(paths)))stop("Figure8 frozen numerical-code authentication failed")
    if(!identical(b$config,jsonlite::fromJSON(paths[1])))stop("Figure8 original numerical configuration changed")
    engine <- f78ro_binding(ctx,"engine")
    if(!identical(b$engine_source_sha,f78ro_hash(engine)))stop("Figure8 original model-engine authentication failed")
    source <- f78ro_binding(ctx,"initializer")
    # Retained snapshots authenticate original receipts only. The initializer
    # definition is loaded from the private locked model, never a retired worker.
    private <- ctx$model_sources$path[ctx$model_sources$relative_path==
      "util/o2_supply_demand_map_invitro_lineage_simulation_utils.R"]
    if(length(private)!=1L||!identical(f78ro_hash(private),f78ro_hash(source)))
      stop("Private locked initializer differs from the retained numerical source")
    env <- new.env(parent=baseenv());sys.source(private,envir=env)
    fun <- get("ivt_gaussian_initial_state",envir=env,inherits=FALSE)
    environment(fun) <- baseenv()
    if(!identical(b$gaussian_initializer_sha,digest::digest(list(formals(fun),body(fun)),algo="sha256")))
      stop("Figure8 canonical Gaussian initializer authentication failed")
    if(any(selected$key=="family_time8")) {
      f78ro_complete(ctx,"family_time8","numerical_validation.tsv")
      f78ro_manifest_guard(ctx,"family_time8","source_manifest.tsv")
      f78ro_manifest_guard(ctx,"family_time8","code_manifest.tsv")
    }
  }
  message("PASS: canonical inputs, original numerical/code receipts, model and initializer authenticated")
}
