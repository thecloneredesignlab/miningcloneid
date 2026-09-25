# Domain-module loader for in-vivo fitted-parameter simulation.
#
# Keep this file deliberately small: numerical implementations live under the
# functional subdirectories below, while the public generator only orchestrates
# reading fitted inputs and writing table contracts.

normalize_cfg_for_simulation <- function(cfg) {
  # Stored fitted configurations use the established non-fitting normalization
  # branch. Keep that semantic contract under a simulation-specific public name.
  normalize_sim_cfg_common(cfg, context = "viz")
}

.invivo_simulation_module_root <- local({
  frame_files <- Filter(
    nzchar,
    vapply(sys.frames(), function(env) {
      ofile <- env$ofile
      if (is.null(ofile)) "" else normalizePath(ofile, mustWork = FALSE)
    }, character(1))
  )
  if (length(frame_files)) dirname(frame_files[[length(frame_files)]]) else getwd()
})

for (.module_path in c(
  file.path("o2", "invivo_o2_simulation.R"),
  file.path("ploidy", "invivo_ploidy_simulation.R"),
  file.path("cin", "invivo_cin_simulation.R"),
  file.path("population", "invivo_population_simulation.R"),
  file.path("functional_response", "invivo_functional_response_simulation.R")
)) {
  source(file.path(.invivo_simulation_module_root, .module_path), local = environment())
}
rm(.module_path, .invivo_simulation_module_root)
