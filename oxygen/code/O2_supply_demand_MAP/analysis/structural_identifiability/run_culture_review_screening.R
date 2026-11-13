#!/usr/bin/env Rscript
# Initialization-only stratification. No production bound or score changes.
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 4L) stop("Usage: run_culture_review_screening.R project study mode task")
  project <- normalizePath(args[1], mustWork = TRUE); study <- normalizePath(args[2], mustWork = TRUE)
  mode <- args[3]; tid <- args[4]
  code <- file.path(study, "implementation/oxygen/code/O2_supply_demand_MAP")
  modules <- file.path(code, "util", c("o2_supply_demand_map_profile.R", "o2_supply_demand_map_profile_culture.R",
                                      "o2_supply_demand_map_pilot.R", "o2_supply_demand_map_review.R"))
  for (path in modules) source(path)
  o2pilot_verify_code(project, study)
  settings <- jsonlite::read_json(file.path(study, "study.json"), simplifyVector = TRUE)
  validation <- jsonlite::read_json(file.path(study, "tests/validation_status.json"), simplifyVector = TRUE)
  if (!isTRUE(validation$integration_pass) || validation$code_manifest_sha256 != o2prof_hash(file.path(study, "manifests/code_manifest.tsv"))) stop("Original fitting implementation is not validated")
  if (Sys.getenv("O2PROF_IMAGE_SHA256") != settings$image_sha256) stop("Wrong runtime image")
  own_path <- file.path(study, "screening_v01/run_culture_review_screening.R")
  ctx <- o2prof_load_culture(file.path(project, settings$wp0_relative_path), modules)
  draw <- function(map, proposal, i) {
    p <- o2review_draw(map, if (proposal == "full_domain") "full_domain" else "growth_favoring", i)
    high <- "p_misseg" %in% names(map$fixed) && map$fixed$p_misseg >= .1
    if (high && proposal != "full_domain") {
      tab <- map$table
      ranges <- list(lam_max = c(1, 1.5), alpha_o2 = c(.5, .75), gamma_growth = c(1.5, 1.75),
          mu_hp = c(.005, .0065), gamma_mu = c(1.2, 1.4), buffer_beta = c(.001, .2),
          buffer_smax = c(.98, 1), buffer_n_exp = c(2, 8))
      k <- tab[tab$parameter == "k_o_mis", ]
      ranges$k_o_mis <- c(k$upper / 2, k$upper)
      for (n in names(ranges)) {
        row <- tab[tab$parameter == n, ]
        if (row$status != "active") next
        low <- max(row$lower, ranges[[n]][1]); upper <- min(row$upper, ranges[[n]][2])
        if (low > upper) stop("Initialization proposal outside declared domain")
        p[[n]] <- runif(1, low, upper)
      }
    }
    o2prof_project_start(map, p)$params
  }
  screen <- function(task, folder, target_count = 64L, max_draws = 5000L, seconds = 3500) {
    if (dir.exists(folder)) stop("Existing screening attempt cannot be overwritten")
    dir.create(folder, recursive = TRUE)
    fixed <- if (length(task$fixed)) as.list(task$fixed) else list()
    map <- o2review_map(ctx, task$domain, fixed)
    set.seed(task$seed); started <- as.numeric(Sys.time()); accepted <- list(); records <- list()
    for (i in seq_len(max_draws)) {
      if (as.numeric(Sys.time()) - started >= seconds || length(accepted) >= target_count) break
      proposal <- if (i %% 2L) "full_domain" else
        if ("p_misseg" %in% names(fixed) && fixed$p_misseg >= .1) "high_induction_viability_stratum" else "growth_favoring"
      p <- draw(map, proposal, i); reason <- ""; q <- NA_real_
      q <- tryCatch({
        comp <- ctx$raw_objective(p, ctx$objects, ctx$cfg); o2prof_check_samples(comp, ctx$manifest); comp$objective
      }, error = function(e) { reason <<- conditionMessage(e); NA_real_ })
      records[[i]] <- list(draw = i, proposal = proposal, params = p, objective = q, error = reason)
      if (is.finite(q)) accepted[[length(accepted) + 1L]] <- p
      if (i %% 25L == 0L) o2pilot_atomic_rds(records, file.path(folder, "screening_draws.rds"))
    }
    o2pilot_atomic_rds(records, file.path(folder, "screening_draws.rds"))
    if (length(accepted) < 8L) {
      o2pilot_json(list(pass = FALSE, feasible_members = length(accepted),
          reason = "Independent feasible-bank threshold not met"), file.path(folder, "screening_result.json"))
      stop("Independent feasible-bank screening unresolved")
    }
    # Preserve first-feasible ordering, without selecting low scores.
    feasible_draws <- which(vapply(records, function(r) is.finite(r$objective), logical(1)))
    for (i in seq_along(accepted)) if (!identical(accepted[[i]], records[[feasible_draws[i]]]$params)) stop("Feasibility-only order changed")
    anchor <- if (isTRUE(task$pure_independent)) accepted[[1]] else
      o2prof_project_start(map, readRDS(file.path(project, task$source_fit))$best_params)$params
    population <- o2prof_initial_population(map, 256L, list(anchor = anchor))$population
    for (i in seq_along(accepted)) population[i, ] <- o2prof_to_free(map, accepted[[i]])
    if (!isTRUE(task$pure_independent)) population[1, ] <- o2prof_to_free(map, anchor)
    saveRDS(population, file.path(folder, "population.rds")); saveRDS(anchor, file.path(folder, "anchor.rds"))
    origin <- data.frame(member = seq_len(256L), origin = ifelse(seq_len(256L) <= length(accepted), "independent_feasible_draw", "full_domain_uniform_transformed"),
        screening_draw = c(feasible_draws, rep(NA_integer_, 256L - length(feasible_draws))))
    if (!isTRUE(task$pure_independent)) { origin$origin[1] <- "transferred_accepted_vector"; origin$screening_draw[1] <- NA_integer_ }
    o2prof_write_tsv(origin, file.path(folder, "member_origins.tsv"))
    result <- list(pass = TRUE, feasible_members = length(accepted), screening_evaluations = length(records),
        seed = task$seed, pure_independent = task$pure_independent, source_fit = task$source_fit,
        proposal_selection = "First protocol-feasible draws, never objective ranking", initialization_only = TRUE,
        full_optimizer_domain_unchanged = TRUE, script_sha256 = o2prof_hash(own_path),
        member_origins_sha256 = o2prof_hash(file.path(folder, "member_origins.tsv")),
        population_sha256 = o2prof_hash(file.path(folder, "population.rds")),
        elapsed_seconds = as.numeric(Sys.time()) - started)
    o2pilot_json(result, file.path(folder, "screening_result.json")); result
  }
  if (mode == "validate") {
    tests <- list()
    for (i in 1:2) {
      task <- list(task_id = paste0("test_", i), fixed = list(p_misseg = .4), domain = "original",
                   pure_independent = TRUE, source_fit = NULL, seed = 710001L + i * 10003L)
      tests[[i]] <- screen(task, file.path(study, "screening_v01/tests", task$task_id), 8L, 500L, 500)
      map <- o2review_map(ctx, "original", task$fixed)
      folder <- file.path(study, "screening_v01/tests", task$task_id)
      population <- readRDS(file.path(folder, "population.rds"))
      for (member in seq_len(8L)) {
        p <- o2prof_from_free(map, setNames(population[member, ], colnames(population)))
        if (!identical(p$p_misseg, .4)) stop("High-induction target not exact")
        comp <- ctx$raw_objective(p, ctx$objects, ctx$cfg)
        if (abs(comp$objective - o2prof_independent_score(ctx, comp)$objective) > 1e-8) stop("Screened-member score verification failed")
      }
    }
    if (tests[[1]]$population_sha256 == tests[[2]]$population_sha256) stop("Independent population seeds produced identical populations")
    o2pilot_json(list(pass = TRUE, script_sha256 = o2prof_hash(own_path), populations = tests,
        production_objective_unchanged = TRUE, score_ranking_used = FALSE,
        completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)), file.path(study, "screening_v01/validation_status.json"))
    cat("Two independent high-induction initialization banks passed feasibility, exact-target and score checks. No optimization performed.\n")
  } else if (mode == "screen") {
    task <- jsonlite::read_json(file.path(study, "tasks", tid, "task.json"), simplifyVector = TRUE)
    validation <- jsonlite::read_json(file.path(study, "screening_v01/validation_status.json"), simplifyVector = TRUE)
    if (!isTRUE(validation$pass) || validation$script_sha256 != o2prof_hash(own_path) || task$screening_script_sha256 != o2prof_hash(own_path)) stop("Unvalidated screening source")
    screen(task, file.path(study, "populations", tid))
  } else stop("Unknown screening mode")
}
main()
