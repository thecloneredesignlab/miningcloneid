# Render-only I/O. Scientific inputs and historical receipts are never written.
f78ro_absolute <- function(path) {
  path <- path.expand(path)
  if (length(path)!=1L || is.na(path) || !nzchar(path)) stop("A nonempty directory path is required")
  if (!startsWith(path,"/")) path <- file.path(getwd(),path)
  parts <- character()
  for (part in strsplit(path,"/",fixed=TRUE)[[1]]) {
    if (!nzchar(part) || part==".") next
    if (part=="..") parts <- head(parts,-1L) else parts <- c(parts,part)
  }
  path <- paste0("/",paste(parts,collapse="/"))
  suffix <- character();ancestor <- path
  while (!file.exists(ancestor)) {suffix <- c(basename(ancestor),suffix);ancestor <- dirname(ancestor)}
  result <- normalizePath(ancestor,mustWork=TRUE)
  for (part in suffix) result <- file.path(result,part)
  result
}
f78ro_no_archive <- function(path) {
  path <- as.character(path)
  if (any(grepl("(^|[/\\\\])(arxives|archives)([/\\\\]|$)",path)))
    stop("Archive access is forbidden in render-only mode",call.=FALSE)
  for (p in path[file.exists(path)]) {
    real <- normalizePath(p,mustWork=TRUE)
    if (grepl("(^|/)(arxives|archives)(/|$)",real))
      stop("An input/output resolves into arxives: ",p,call.=FALSE)
  }
  invisible(path)
}
f78ro_hash <- function(paths) {
  f78ro_no_archive(paths)
  if (any(!file.exists(paths))) stop("Missing required file: ",paste(paths[!file.exists(paths)],collapse=", "),call.=FALSE)
  unname(vapply(paths,digest::digest,"",file=TRUE,algo="sha256"))
}
f78ro_write <- function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  data.table::fwrite(x,path,sep="\t",na="NA")
}
f78ro_input <- function(ctx,path) {
  f78ro_no_archive(path)
  relative <- substring(path,nchar(ctx$data_root)+2L)
  if (any(!startsWith(path,paste0(ctx$data_root,"/"))))
    stop("Scientific input outside the configured data root",call.=FALSE)
  i <- match(relative,ctx$inputs$relative_path)
  if (anyNA(i)) stop("Unregistered render input: ",paste(relative[is.na(i)],collapse=", "),call.=FALSE)
  for (j in seq_along(path)) {
    if (path[j] %in% ctx$checked$path) next
    if (!file.exists(path[j])) stop("Missing retained plotting input: ",path[j],
      "; no fallback, archive restoration, or simulation will be attempted",call.=FALSE)
    observed <- f78ro_hash(path[j])
    if (!identical(observed,ctx$inputs$sha256[i[j]]))
      stop("Retained input SHA-256 mismatch: ",path[j],call.=FALSE)
    ctx$checked <- rbind(ctx$checked,data.frame(path=path[j],relative_path=relative[j],sha256=observed))
  }
  invisible(path)
}
f78ro_complete <- function(ctx,key,validation="final_numerical_validation.tsv") {
  paths <- f78ro_data(ctx,key,c("completion.tsv",validation))
  done <- data.table::fread(paths[1]);checks <- data.table::fread(paths[2])
  if (!all(done$status %in% c("numerical_complete","validated_complete")) || !all(checks$passed))
    stop("Numerical completion/validation failed for ",key,call.=FALSE)
}
f78ro_component_stems <- function(stem) {
  if (!stem %in% c("supp_fig7_3","supp_fig7_4","supp_fig7_3_ver2","supp_fig7_4_ver2")) return(character())
  unlist(lapply(c("invivo","invitro"),function(context)
    unlist(lapply(c("c01","c02"),function(family)
      paste(stem,context,family,c("4N","2N"),sep="_")))))
}
f78ro_image_snapshot <- function(ctx) {
  files <- character()
  for (root in unique(c(ctx$workspace,ctx$output_root))) for (folder in c("Figures","manuscript/Figures")) {
    directory <- file.path(root,folder)
    if (!dir.exists(directory)) next
    paths <- list.files(directory,pattern="\\.(pdf|png|svg)$",full.names=TRUE)
    target <- root==ctx$output_root & basename(paths) %in% ctx$allowed_files
    files <- c(files,paths[!target])
  }
  files <- sort(unique(files))
  data.frame(path=files,sha256=f78ro_hash(files))
}
f78ro_guard_current <- function(ctx) {
  current <- f78ro_image_snapshot(ctx)
  if (!identical(ctx$protected,current)) stop("A non-target image was added, removed, or changed during rendering",call.=FALSE)
  if (!identical(ctx$checked$sha256,f78ro_hash(ctx$checked$path)))
    stop("Scientific inputs changed during rendering",call.=FALSE)
  if (nrow(ctx$authenticated) && !identical(ctx$authenticated$sha256,f78ro_hash(ctx$authenticated$path)))
    stop("Authenticated numerical source/code changed during rendering",call.=FALSE)
  if (!identical(ctx$render_sources$sha256,f78ro_hash(ctx$render_sources$path)))
    stop("Drawing code changed during rendering",call.=FALSE)
  invisible(TRUE)
}
f78ro_publish <- function(a,plot,width,height,stem,dpi=300) {
  ctx <- a$render_ctx
  if (!stem %in% ctx$allowed_stems) stop("Attempt to export a non-target product: ",stem,call.=FALSE)
  plot <- f78_normalize_annotation_theme(plot)
  for (format in ctx$formats) {
    path <- file.path(ctx$stage,"Figures",paste0(stem,".",format))
    if (format=="pdf") ggplot2::ggsave(path,plot,width=width,height=height,units="in",device=grDevices::cairo_pdf,bg="white",limitsize=FALSE)
    else ggplot2::ggsave(path,plot,width=width,height=height,units="in",dpi=dpi,bg="white",limitsize=FALSE)
    mirror <- file.path(ctx$stage,"manuscript/Figures",basename(path))
    if (!file.copy(path,mirror,overwrite=TRUE)) stop("Failed to stage the publication mirror")
    ctx$exports[[basename(path)]] <- data.frame(stem=stem,format=format,width_inches=width,
      height_inches=height,dpi=if(format=="png")dpi else NA_real_,sha256=f78ro_hash(path))
  }
  message("Rendered: ",stem)
  invisible(NULL)
}
f78ro_python <- function() {
  candidates <- c(Sys.getenv("FIGURE78_PYTHON",unset=""),
    "/Users/4482173/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3",Sys.which("python3"))
  candidates <- candidates[nzchar(candidates)&file.exists(candidates)]
  if (!length(candidates)) stop("Python with pypdf is required to package matrix PDFs")
  candidates[1]
}
f78ro_matrix <- function(a,z,metric,stem) {
  ctx <- a$render_ctx
  pages <- f7_matrix_pages(a,z,metric,stem,publish=FALSE)
  for (i in seq_along(pages$plots)) f78ro_publish(a,pages$plots[[i]],
    if(grepl("invitro",pages$labels[i]))17 else 15,10,pages$labels[i])
  if ("png" %in% ctx$formats) {
    temporary <- ctx$formats;ctx$formats <- "png"
    tryCatch(f78ro_publish(a,patchwork::wrap_plots(pages$plots,ncol=2),24,30,stem,dpi=120),
      finally={ctx$formats <- temporary})
  }
  if ("pdf" %in% ctx$formats) {
    f78ro_write(data.frame(page=seq_along(pages$labels),matrix_artifact=pages$labels),
      file.path(ctx$audit,paste0(stem,"_page_manifest.tsv")))
    status <- system2(f78ro_python(),c(shQuote(file.path(a$code,"merge_Figure7_matrix_pdfs.py")),
      "--workspace",shQuote(ctx$stage),"--stem",shQuote(stem),"--run",shQuote(ctx$audit),"--audit",shQuote(ctx$audit)))
    if (status!=0L) stop("Matrix PDF packaging failed: ",stem)
    path <- file.path(ctx$stage,"Figures",paste0(stem,".pdf"))
    ctx$exports[[basename(path)]] <- data.frame(stem=stem,format="pdf",width_inches=NA_real_,
      height_inches=10,dpi=NA_real_,sha256=f78ro_hash(path))
  }
}
f78ro_finish <- function(ctx) {
  expected <- as.vector(outer(ctx$allowed_stems,paste0(".",ctx$formats),paste0))
  if (!setequal(expected,names(ctx$exports))) stop("Requested export coverage is incomplete")
  for (name in expected) {
    paths <- file.path(ctx$stage,c("Figures","manuscript/Figures"),name)
    if (!all(file.exists(paths)) || length(unique(f78ro_hash(paths)))!=1L)
      stop("Staged export/mirror is incomplete: ",name)
  }
  f78ro_guard_current(ctx)
  # Back up only the selected targets; do not rewrite historical audit receipts.
  for (folder in c("Figures","manuscript/Figures")) for (name in expected) {
    dest <- file.path(ctx$output_root,folder,name)
    if (file.exists(dest)) {
      backup <- file.path(ctx$audit,"before",folder,name)
      dir.create(dirname(backup),recursive=TRUE,showWarnings=FALSE)
      if (!file.copy(dest,backup,overwrite=FALSE)) stop("Could not back up selected output: ",dest)
    }
  }
  for (folder in c("Figures","manuscript/Figures")) for (name in expected) {
    source <- file.path(ctx$stage,folder,name);dest <- file.path(ctx$output_root,folder,name)
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    temporary <- tempfile(pattern=paste0(".",name,"_"),tmpdir=dirname(dest))
    if (!file.copy(source,temporary,overwrite=FALSE) || !identical(f78ro_hash(source),f78ro_hash(temporary)) ||
        !file.rename(temporary,dest)) stop("Atomic output publication failed: ",dest)
  }
  f78ro_guard_current(ctx)
  f78ro_write(data.table::rbindlist(ctx$exports),file.path(ctx$audit,"artifact_manifest.tsv"))
  f78ro_write(ctx$checked,file.path(ctx$audit,"inputs_verified.tsv"))
  f78ro_write(ctx$authenticated,file.path(ctx$audit,"numerical_sources_authenticated.tsv"))
  f78ro_write(f78ro_image_snapshot(ctx),file.path(ctx$audit,"protected_after.tsv"))
  f78ro_write(data.frame(status="complete",products=length(ctx$products),files_per_mirror=length(expected),
    archive_reads=0L,new_simulations=0L,non_target_images_unchanged=TRUE),file.path(ctx$audit,"completion.tsv"))
  message("PASS: ",length(ctx$products)," requested products; no archive access or simulation. Audit: ",ctx$audit)
}
