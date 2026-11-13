#!/usr/bin/env Rscript
# Verify actual out-of-original-domain values reach the unchanged R/C++ model.
args <- commandArgs(trailingOnly = TRUE)
project <- normalizePath(args[1], mustWork = TRUE); study <- normalizePath(args[2], mustWork = TRUE)
code <- file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP")
modules <- file.path(code, "util", c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R",
                                    "o2_supply_demand_map_pilot.R", "o2_supply_demand_map_review.R"))
for (path in modules) source(path)
o2pilot_verify_code(project, study)
settings <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
if (Sys.getenv("O2PROF_IMAGE_SHA256") != settings$image_sha256) stop("Wrong image fingerprint")
ctx <- o2prof_load_culture(file.path(project, settings$wp0_relative_path), modules)
helper <- file.path(study, "reporting_v01/culture_profile_diagnostics.R"); source(helper)
reference <- readRDS(file.path(project, settings$parent_relative_path, "fits/fit_090/accepted_fit.rds"))$best_params
records <- list()
for (domain in c("D_base", "D_sat")) {
  fixed <- if (domain == "D_base") list(p_misseg = 0) else list()
  map <- o2review_map(ctx, domain, fixed)
  p <- o2prof_project_start(map, reference)$params
  name <- if (domain == "D_base") "p_mis_base" else "k_o_mis"
  probe <- if (domain == "D_base") .015 else .1
  p[[name]] <- probe
  normalized <- ctx$backend$normalize_run_params_common(p, cfg = ctx$cfg)
  if (!identical(as.numeric(normalized[[name]]), probe)) stop("Expanded parameter was silently clipped")
  o2prof_to_free(map, normalized)
  for (oxygen in c(0, .5, 5, 20.5)) {
    # This checks all 133 chromosome states against the production C++ generator,
    # including daughter accounting, rather than only a symbolic parameter map.
    kernel <- pilotdiag_rates(ctx, normalized, ctx$cfg, oxygen)
    for (N in c(44L, 88L)) {
      mu <- ctx$backend$.mu_eff_of_O2(oxygen, normalized, N)
      expected <- normalized$p_mis_base + normalized$p_misseg * mu / (mu + normalized$k_o_mis)
      actual <- ctx$backend$.pmisseg_of_O2(oxygen, normalized, N)
      if (abs(actual - expected) > 1e-12) stop("Expanded-domain response does not follow the production function")
      if (abs(kernel$rates$p_per_chromosome[kernel$rates$N == N] - actual) > 1e-12) stop("Kernel/response probability mismatch")
      old <- normalized; old[[name]] <- if (domain == "D_base") .005 else .02
      old_probability <- ctx$backend$.pmisseg_of_O2(oxygen, old, N)
      if (abs(actual - old_probability) < 1e-8) stop("Expanded value did not change the intended mechanism")
      records[[length(records) + 1L]] <- data.frame(domain = domain, parameter = name, probe = probe,
          oxygen_pct = oxygen, N = N, actual = actual, expected = expected, original_limit_response = old_probability)
    }
  }
}
out <- file.path(study, "domain_validation_v01")
o2prof_write_tsv(do.call(rbind, records), file.path(out, "response_checks.tsv"))
o2pilot_json(list(pass = TRUE, tested_domains = c("D_base", "D_sat"), response_checks = length(records),
    generator_ledger_cases = 8L * (ctx$cfg$N_MAX - ctx$cfg$N_MIN + 1L),
    fitting_code_manifest_sha256 = o2prof_hash(file.path(study, "manifests/code_manifest.tsv")),
    helper_sha256 = o2prof_hash(helper), image_sha256 = settings$image_sha256,
    model_and_objective_changed = FALSE, optimization_performed = FALSE,
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(out, "validation_status.json"))
cat("Expanded values survive normalization and reach the production R/C++ mechanisms; no optimization performed.\n")
