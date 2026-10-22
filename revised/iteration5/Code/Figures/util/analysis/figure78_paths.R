# Figure7/8-only path policy. No source checkout or archived-result fallback.
f78_current_path <- function(path,workspace,mustWork=FALSE) {
  if(length(path)!=1L||is.na(path)||!nzchar(path))stop("A nonempty Figure7/8 path is required")
  path <- path.expand(path)
  if(!startsWith(path,"/"))path <- file.path(workspace,path)
  if(grepl("(^|/)(arxives|archives|frozen_sources)(/|$)|/figure[78]_[^/]*_v[0-9]+(/|$)",path))
    stop("Old-version/archive roots are not Figure7/8 runtime inputs: ",path)
  # Resolve existing ancestors as well as complete paths, including symlinks.
  tail <- character();parent <- path
  while(!file.exists(parent)) {
    if(identical(parent,dirname(parent)))stop("Cannot resolve path: ",path)
    tail <- c(basename(parent),tail);parent <- dirname(parent)
  }
  result <- normalizePath(parent,mustWork=TRUE)
  for(part in tail) {
    if(part=="..")result <- dirname(result) else if(part!=".")result <- file.path(result,part)
  }
  if(grepl("(^|/)(arxives|archives|frozen_sources)(/|$)|/figure[78]_[^/]*_v[0-9]+(/|$)",result))
    stop("Figure7/8 path resolves to an old-version/archive root: ",result)
  if(mustWork&&!file.exists(result))stop("Required Figure7/8 path is missing: ",result)
  result
}
f78_model_root <- function(workspace,requested=NULL) {
  expected <- f78_current_path("Code/vendor/figure78_model",workspace,mustWork=TRUE)
  override <- if(!is.null(requested))requested else Sys.getenv("FIGURE_MODEL_CODE_ROOT",unset="")
  if(nzchar(override)&&!identical(f78_current_path(override,workspace,mustWork=TRUE),expected))
    stop("Figure7/8 requires its private locked model: ",expected,
      "; an external/shared model override is not permitted")
  expected
}
f78_fit_root <- function(workspace,shared,requested=NULL) {
  value <- if(!is.null(requested))requested else Sys.getenv("FIGURE_JOINT_RESULT_ROOT",unset=shared$fit_root)
  f78_current_path(value,workspace,mustWork=TRUE)
}
f78_locked_model_manifest <- function(workspace,requested=NULL) {
  root <- f78_model_root(workspace,requested)
  lock <- utils::read.delim(file.path(workspace,"Code/config/manifests/figure78_model_sources.tsv"),
    stringsAsFactors=FALSE,check.names=FALSE)
  if(nrow(lock)!=12L||anyDuplicated(lock$relative_path)||
    any(grepl("(^/|(^|/)\\.\\.(/|$))",lock$relative_path)))stop("Invalid Figure7/8 model lock")
  paths <- file.path(root,lock$relative_path)
  if(any(!file.exists(paths)))stop("Private locked model source is missing")
  if(any(!startsWith(normalizePath(paths,mustWork=TRUE),paste0(root,"/"))))
    stop("Private model source resolves outside the dedicated model directory")
  observed <- unname(vapply(paths,digest::digest,"",file=TRUE,algo="sha256"))
  if(!identical(observed,lock$sha256))stop("Private Figure7/8 model SHA-256 differs from its lock")
  data.frame(path=paths,relative_path=lock$relative_path,sha256=observed,stringsAsFactors=FALSE)
}
