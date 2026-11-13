# Pilot bookkeeping only: the frozen production backend remains the fitter.

o2pilot_atomic_rds <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(".candidate_", tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(value, temporary)
  if (!file.rename(temporary, path)) stop("checkpoint_write_failed: atomic rename")
  invisible(path)
}

o2pilot_checkpoint <- function(params, objective, directory, fingerprint, phase = getOption("o2prof.phase", "unknown")) {
  if (!is.finite(objective)) stop("Nonfinite candidate cannot be checkpointed")
  path <- file.path(directory, paste0("pid_", Sys.getpid(), ".rds"))
  # .GlobalEnv belongs to the process and survives repeated DE closure
  # serialization. Disk is the fallback after cache reset/process restart.
  if (!exists(".o2pilot_best_cache", envir = .GlobalEnv, inherits = FALSE))
    assign(".o2pilot_best_cache", new.env(parent = emptyenv()), envir = .GlobalEnv)
  cache <- get(".o2pilot_best_cache", envir = .GlobalEnv)
  key <- paste(fingerprint, path, sep = "|")
  if (exists(key, envir = cache, inherits = FALSE)) {
    if (get(key, envir = cache) <= objective) return(invisible(FALSE))
  } else if (file.exists(path)) {
    previous <- readRDS(path)
    if (!identical(previous$fingerprint, fingerprint)) stop("Stale candidate checkpoint")
    assign(key, previous$objective, envir = cache)
    if (previous$objective <= objective) return(invisible(FALSE))
  }
  tryCatch(o2pilot_atomic_rds(list(params = params, objective = objective, phase = phase,
                                fingerprint = fingerprint, pid = Sys.getpid(),
                                recorded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), path),
           error = function(e) stop("checkpoint_write_failed: ", conditionMessage(e)))
  assign(key, objective, envir = cache)
  invisible(TRUE)
}

o2pilot_track_context <- function(ctx, directory, fingerprint) {
  raw <- ctx$raw_objective
  manifest <- ctx$manifest
  ctx$raw_objective <- function(run_params, fit_objects, cfg, ...) {
    comp <- raw(run_params = run_params, fit_objects = fit_objects, cfg = cfg, ...)
    o2prof_check_samples(comp, manifest)
    o2pilot_checkpoint(run_params, comp$objective, directory, fingerprint)
    comp
  }
  ctx
}

o2pilot_modules <- function(study) {
  file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP/util",
            c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R", "o2_supply_demand_map_pilot.R"))
}

o2pilot_verify_code <- function(project, study) {
  manifest <- read.delim(file.path(study, "manifests/code_manifest.tsv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(study, "implementation", manifest$path[i])
    if (o2prof_hash(path) != manifest$sha256[i]) stop("Frozen pilot implementation changed: ", path)
  }
  invisible(TRUE)
}

o2pilot_json <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
}
