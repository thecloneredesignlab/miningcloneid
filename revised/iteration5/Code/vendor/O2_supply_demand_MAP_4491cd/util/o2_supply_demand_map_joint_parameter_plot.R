# Backward-compatible data-helper loader.
#
# Plot constructors moved to:
#   vis/joint/o2_supply_demand_map_joint_parameter_plot_utils.R
# This historical path intentionally exposes only the table-building API.

.joint_parameter_legacy_dir <- local({
  frame_files <- Filter(nzchar, vapply(sys.frames(), function(env) {
    ofile <- env$ofile
    if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
  }, character(1)))
  own <- frame_files[
    basename(frame_files) == "o2_supply_demand_map_joint_parameter_plot.R"
  ]
  if (length(own)) return(dirname(own[[length(own)]]))
  normalizePath(getwd(), mustWork = FALSE)
})
source(
  file.path(
    .joint_parameter_legacy_dir,
    "o2_supply_demand_map_joint_parameter_utils.R"
  ),
  local = environment(),
  chdir = TRUE
)
rm(.joint_parameter_legacy_dir)
