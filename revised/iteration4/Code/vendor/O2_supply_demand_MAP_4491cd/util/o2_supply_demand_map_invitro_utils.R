# Canonical loader for the in-vitro fitting helper library.
#
# This file is intentionally a library loader, not an executable entrypoint.
# It resolves paths from its own source location so callers may source it from
# any working directory.

.o2sd_invitro_loader_dir <- local({
  frame_files <- Filter(
    nzchar,
    vapply(
      sys.frames(),
      function(env) {
        ofile <- env$ofile
        if (is.null(ofile)) return("")
        ofile <- as.character(ofile[[1]])
        if (!startsWith(ofile, "/")) {
          basename_candidate <- file.path(getwd(), basename(ofile))
          if (file.exists(basename_candidate)) {
            ofile <- basename_candidate
          }
        }
        normalizePath(ofile, mustWork = FALSE)
      },
      character(1)
    )
  )
  if (length(frame_files) == 0L) {
    stop("Cannot resolve o2_supply_demand_map_invitro_utils.R source directory.")
  }
  dirname(frame_files[[length(frame_files)]])
})

.o2sd_invitro_helper_files <- c(
  "o2_supply_demand_map_invitro_io_utils.R",
  "o2_supply_demand_map_invitro_lineage_utils.R",
  "o2_supply_demand_map_invitro_summary_utils.R",
  "o2_supply_demand_map_invitro_lineage_simulation_utils.R",
  "o2_supply_demand_map_invitro_objective_utils.R"
)

for (.o2sd_invitro_helper_name in .o2sd_invitro_helper_files) {
  .o2sd_invitro_helper_path <- file.path(
    .o2sd_invitro_loader_dir,
    .o2sd_invitro_helper_name
  )
  if (!file.exists(.o2sd_invitro_helper_path)) {
    stop("Missing canonical in-vitro helper file: ", .o2sd_invitro_helper_path)
  }
  sys.source(
    .o2sd_invitro_helper_path,
    envir = environment(),
    chdir = TRUE
  )
}

rm(
  .o2sd_invitro_helper_files,
  .o2sd_invitro_helper_name,
  .o2sd_invitro_helper_path,
  .o2sd_invitro_loader_dir
)
