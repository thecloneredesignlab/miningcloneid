.load_invitro_observed_feasibility_refinement_helpers <- function() {
  env <- new.env(parent = baseenv())
  sys.source(
    file.path(
      repo_info$root,
      "oxygen",
      "scripts",
      "invitro_observed_feasibility_utils.R"
    ),
    envir = env
  )
  env
}

.load_invitro_observed_feasibility_refinement_driver <- function() {
  env <- new.env(parent = globalenv())
  sys.source(
    file.path(
      repo_info$root,
      "oxygen",
      "scripts",
      "scan_invitro_observed_feasibility_refinement.R"
    ),
    envir = env
  )
  env
}

.ivt_refinement_parameter_names <- c(
  "p_wgd", "p_mis_base", "gamma_growth", "mu_hp", "gamma_mu",
  "p_misseg", "k_o_mis", "buffer_smax", "buffer_beta",
  "buffer_n_exp", "O2_crit", "n_O"
)

.ivt_refinement_envelope <- function() {
  data.frame(
    p_wgd = c(1e-4, 5e-4, 2.5e-3),
    p_mis_base_upper = c(8e-4, 4e-4, 8e-5),
    eligible = TRUE,
    stringsAsFactors = FALSE
  )
}

.ivt_refinement_conservative_cap <- function(p_wgd, envelope) {
  vapply(as.numeric(p_wgd), function(value) {
    exact <- which(envelope$p_wgd == value)
    if (length(exact)) return(envelope$p_mis_base_upper[[exact[[1L]]]])
    right <- which(envelope$p_wgd > value)[[1L]]
    left <- right - 1L
    min(envelope$p_mis_base_upper[c(left, right)])
  }, numeric(1))
}

.ivt_refinement_anchor <- function(anchor_id = "seed10") {
  data.frame(
    anchor_id = anchor_id,
    p_wgd = 5e-4,
    p_mis_base = 1e-5,
    gamma_growth = 0.20,
    mu_hp = 0.002,
    gamma_mu = 1.5,
    p_misseg = 0.01,
    k_o_mis = 5e-4,
    buffer_smax = 0.95,
    buffer_beta = 1,
    buffer_n_exp = 6,
    O2_crit = 1,
    n_O = 2,
    stringsAsFactors = FALSE
  )
}

.ivt_refinement_good_target_row <- function() {
  data.frame(
    candidate_id = "candidate_good",
    status = "OK",
    protocol_feasibility_status = "PASS",
    n_scenarios = 6,
    n_summary_rows = 114,
    n_growth = 219,
    n_death = 90,
    n_flow_samples = 20,
    n_ploidy_samples = 12,
    n_low_o2_flow_groups = 6,
    n_insufficient_boundaries = 0,
    n_missing_predictions = 0,
    n_negative_predictions = 0,
    n_invalid_distributions = 0,
    control_flow_mass_ge3 = 0.05,
    control_flow_wgd_mass = 0.05,
    all_passage_boundaries_feasible = TRUE,
    all_passage_executed = TRUE,
    all_selected_at_or_after_last_observation = TRUE,
    all_selected_within_tolerance = TRUE,
    pooled_2N_A19_mean_N = 77.45,
    pooled_2N_A19_mass_ge70 = 0.625,
    pooled_4N_A17_mean_N = 82,
    decline_4N_to_A17_N = 6,
    flow_2N_low_A12_wgd_mass = 0.30,
    flow_2N_low_A18_wgd_mass = 0.35,
    flow_2N_low_A23_wgd_mass = 0.40,
    observed_flow_2N_low_A12_wgd_mass = 0.40,
    observed_flow_2N_low_A18_wgd_mass = 0.40,
    observed_flow_2N_low_A23_wgd_mass = 0.45,
    flow_2N_low_A12_bimodal_count = 2,
    flow_2N_low_A18_bimodal_count = 2,
    flow_2N_low_A23_bimodal_count = 2,
    stringsAsFactors = FALSE
  )
}

testthat::test_that("joint normoxic calibration produces a conditional p_mis_base envelope", {
  env <- .load_invitro_observed_feasibility_refinement_helpers()
  calibration <- data.frame(
    p_wgd = rep(c(1e-4, 5e-4, 2e-3, 2.5e-3), each = 3),
    p_mis_base = rep(c(1e-7, 1e-6, 1e-5), times = 4),
    high_mass_ge_lower = c(
      0.01, 0.03, 0.06,
      0.02, 0.06, 0.07,
      0.06, 0.07, 0.08,
      0.01, 0.03, 0.06
    ),
    status = "OK",
    protocol_feasibility_status = "PASS",
    stringsAsFactors = FALSE
  )

  envelope <- env$ivt_scan_select_conditional_envelope(calibration, control_limit = 0.05)

  testthat::expect_identical(envelope$p_wgd, c(1e-4, 5e-4, 2e-3, 2.5e-3))
  testthat::expect_equal(envelope$p_mis_base_upper, c(1e-6, 1e-7, NA_real_, 1e-6))
  testthat::expect_identical(envelope$n_prefix_pass, c(2L, 1L, 0L, 2L))
  testthat::expect_identical(envelope$eligible, c(TRUE, TRUE, FALSE, FALSE))
  testthat::expect_error(
    env$ivt_scan_select_conditional_envelope(
      calibration[, names(calibration) != "p_wgd"],
      control_limit = 0.05
    ),
    "p_wgd"
  )
})

testthat::test_that("focused and outer specifications lock the planned 12-dimensional ranges", {
  env <- .load_invitro_observed_feasibility_refinement_helpers()
  focused <- env$ivt_scan_extended_parameter_specs("focused")
  outer <- env$ivt_scan_extended_parameter_specs("outer")

  testthat::expect_identical(focused$param_symbol, .ivt_refinement_parameter_names)
  testthat::expect_identical(outer$param_symbol, .ivt_refinement_parameter_names)
  testthat::expect_equal(
    focused$lower,
    c(2e-4, 1e-7, 0.01, 2e-4, 0.3, 0.001, 1e-5, 0.90, 0.2, 4, 1e-4, 0.5)
  )
  testthat::expect_equal(
    focused$upper,
    c(2.5e-3, 5e-3, 0.35, 0.015, 3.5, 0.03, 0.005, 0.9999, 2, 10, 2.5, 5)
  )
  testthat::expect_equal(
    outer$lower,
    c(1e-4, 1e-7, 0.01, 2e-4, 0.1, 0.001, 1e-5, 0.80, 0.2, 4, 1e-4, 0.5)
  )
  testthat::expect_equal(
    outer$upper,
    c(2.5e-3, 5e-3, 0.50, 0.020, 3.5, 0.03, 0.020, 0.9999, 2, 10, 2.5, 5)
  )
  testthat::expect_identical(
    focused$transform,
    c(
      "log", "conditional_log", "log", "log", "linear", "log", "log",
      "complement_log", "log", "log", "log", "linear"
    )
  )
  testthat::expect_error(
    env$ivt_scan_extended_parameter_specs("unsupported"),
    "focused.*outer|outer.*focused"
  )
})

testthat::test_that("extended stage-1 design is deterministic, bounded, conditional, and anchor-safe", {
  env <- .load_invitro_observed_feasibility_refinement_helpers()
  envelope <- .ivt_refinement_envelope()
  anchors <- .ivt_refinement_anchor()
  anchors$p_misseg <- 0.04

  set.seed(177)
  before <- .Random.seed
  design_a <- env$ivt_scan_extended_design(
    n_focused = 24,
    n_outer = 8,
    seed = 101,
    envelope = envelope,
    anchors = anchors
  )
  testthat::expect_identical(.Random.seed, before)
  design_b <- env$ivt_scan_extended_design(
    n_focused = 24,
    n_outer = 8,
    seed = 101,
    envelope = envelope,
    anchors = anchors
  )

  testthat::expect_identical(design_a, design_b)
  testthat::expect_equal(nrow(design_a), 33L)
  testthat::expect_true(all(c(
    "candidate_id", "stage", "region", "parent_candidate_id", "is_anchor",
    "conditional_p_mis_base_upper", .ivt_refinement_parameter_names
  ) %in% names(design_a)))
  testthat::expect_identical(anyDuplicated(design_a$candidate_id), 0L)
  testthat::expect_equal(
    as.integer(table(design_a$region)),
    c(1L, 24L, 8L)
  )
  testthat::expect_identical(
    names(table(design_a$region)),
    c("anchor", "focused", "outer")
  )
  testthat::expect_true(all(design_a$stage == "stage1"))
  testthat::expect_true(all(is.na(design_a$parent_candidate_id)))

  generated <- design_a[!design_a$is_anchor, , drop = FALSE]
  focused <- generated$region == "focused"
  outer <- generated$region == "outer"
  focused_specs <- env$ivt_scan_extended_parameter_specs("focused")
  outer_specs <- env$ivt_scan_extended_parameter_specs("outer")
  for (name in .ivt_refinement_parameter_names) {
    focus_idx <- match(name, focused_specs$param_symbol)
    outer_idx <- match(name, outer_specs$param_symbol)
    testthat::expect_true(all(generated[[name]][focused] >= focused_specs$lower[[focus_idx]]))
    testthat::expect_true(all(generated[[name]][focused] <= focused_specs$upper[[focus_idx]]))
    testthat::expect_true(all(generated[[name]][outer] >= outer_specs$lower[[outer_idx]]))
    testthat::expect_true(all(generated[[name]][outer] <= outer_specs$upper[[outer_idx]]))
  }
  testthat::expect_true(all(is.finite(generated$conditional_p_mis_base_upper)))
  testthat::expect_true(all(
    generated$p_mis_base <= generated$conditional_p_mis_base_upper
  ))
  testthat::expect_equal(
    generated$conditional_p_mis_base_upper,
    .ivt_refinement_conservative_cap(generated$p_wgd, envelope)
  )
  focused_buffer <- generated$buffer_smax[focused]
  focused_buffer_u <- (
    log(1 - focused_buffer) - log(1 - 0.90)
  ) / (
    log(1 - 0.9999) - log(1 - 0.90)
  )
  testthat::expect_equal(sort(floor(focused_buffer_u * 24)), 0:23)

  anchor <- design_a[design_a$is_anchor, , drop = FALSE]
  testthat::expect_identical(anchor$candidate_id, "anchor_seed10")
  testthat::expect_equal(
    unlist(anchor[1, .ivt_refinement_parameter_names], use.names = FALSE),
    unlist(anchors[1, .ivt_refinement_parameter_names], use.names = FALSE)
  )
  testthat::expect_error(
    env$ivt_scan_extended_design(
      n_focused = 2,
      n_outer = 1,
      seed = 1,
      envelope = envelope,
      anchors = rbind(anchors, anchors)
    ),
    "anchor.*unique|unique.*anchor"
  )
})

testthat::test_that("target gate reports technical and phenotype layers independently", {
  env <- .load_invitro_observed_feasibility_refinement_helpers()
  good <- .ivt_refinement_good_target_row()
  technical_failure <- transform(good, candidate_id = "technical_failure", n_missing_predictions = 1)
  phenotype_failure <- transform(
    good,
    candidate_id = "phenotype_failure",
    flow_2N_low_A18_wgd_mass = observed_flow_2N_low_A18_wgd_mass - 0.1501
  )
  rows <- rbind(good, technical_failure, phenotype_failure)

  gated <- env$ivt_scan_target_gate(rows, control_limit = 0.05)

  expected_columns <- c(
    "technical_pass", "bimodality_pass", "wgd_mass_pass", "kary_2N_pass",
    "kary_4N_pass", "bimodality_gap", "wgd_mass_gap", "kary_2N_gap",
    "kary_4N_gap", "target_gap_score", "target_gate_failure_count",
    "phenotype_pass", "full_target_pass"
  )
  testthat::expect_true(all(expected_columns %in% names(gated)))
  testthat::expect_identical(gated$technical_pass, c(TRUE, FALSE, TRUE))
  testthat::expect_identical(gated$phenotype_pass, c(TRUE, TRUE, FALSE))
  testthat::expect_identical(gated$full_target_pass, c(TRUE, FALSE, FALSE))
  testthat::expect_equal(gated$target_gap_score[1:2], c(0, 0))
  testthat::expect_identical(gated$target_gate_failure_count, c(0L, 0L, 1L))
  testthat::expect_gt(gated$wgd_mass_gap[[3L]], 0)
  testthat::expect_gt(gated$target_gap_score[[3L]], 0)

  missing_flow_group <- transform(
    good,
    n_low_o2_flow_groups = 5,
    flow_2N_low_A18_bimodal_count = 1
  )
  missing_flow_gate <- env$ivt_scan_target_gate(missing_flow_group, control_limit = 0.05)
  testthat::expect_false(missing_flow_gate$technical_pass)
  testthat::expect_false(missing_flow_gate$bimodality_pass)
  testthat::expect_false(missing_flow_gate$full_target_pass)

  boundaries <- transform(
    good,
    pooled_2N_A19_mean_N = 77.45 + 8,
    pooled_2N_A19_mass_ge70 = 0.625 - 0.20,
    pooled_4N_A17_mean_N = 81.475 + 4,
    decline_4N_to_A17_N = 6,
    flow_2N_low_A12_wgd_mass = 0,
    observed_flow_2N_low_A12_wgd_mass = 0.15,
    flow_2N_low_A18_wgd_mass = 0.15,
    observed_flow_2N_low_A18_wgd_mass = 0,
    flow_2N_low_A23_wgd_mass = 0,
    observed_flow_2N_low_A23_wgd_mass = 0.15
  )
  boundary_gate <- env$ivt_scan_target_gate(boundaries, control_limit = 0.05)
  testthat::expect_true(boundary_gate$full_target_pass)
  testthat::expect_equal(boundary_gate$target_gap_score, 0)
  testthat::expect_false(
    env$ivt_scan_target_gate(
      transform(boundaries, decline_4N_to_A17_N = 6 + 1e-6),
      control_limit = 0.05
    )$kary_4N_pass
  )
})

testthat::test_that("stage-2 parent selection reserves coverage for nearest target gaps", {
  env <- .load_invitro_observed_feasibility_refinement_helpers()
  specs <- env$ivt_scan_extended_parameter_specs("focused")
  parents <- .ivt_refinement_anchor()[rep(1L, 4L), , drop = FALSE]
  parents$candidate_id <- c("z_nearest", "a_distant", "b_distant", "c_distant")
  parents$p_wgd <- c(5e-4, 2e-4, 1e-3, 2e-3)
  parents$hard_pass <- TRUE
  parents$pareto_front <- FALSE
  parents$full_target_pass <- FALSE
  parents$target_failure_class <- "four_n_close_two_n_fail"
  parents$target_gap_score <- c(0.01, 2, 1, 1.5)
  parents$objective <- c(50, 10, 20, 30)
  parents$bimodality_pass <- TRUE
  parents$wgd_mass_pass <- FALSE
  parents$kary_2N_pass <- TRUE
  parents$kary_4N_pass <- TRUE

  selected <- env$ivt_scan_select_refinement_parents(
    stage1_summary = parents,
    specs = specs,
    max_parents = 2
  )

  testthat::expect_equal(nrow(selected), 2L)
  testthat::expect_identical(selected$candidate_id[[1L]], "z_nearest")
  testthat::expect_true("z_nearest" %in% selected$candidate_id)
  testthat::expect_error(
    env$ivt_scan_select_refinement_parents(
      transform(parents, target_gap_score = NA_real_),
      specs = specs,
      max_parents = 2
    ),
    "finite target_gap_score"
  )
})

testthat::test_that("adaptive stage-2 design is local, deterministic, and records every parent", {
  env <- .load_invitro_observed_feasibility_refinement_helpers()
  specs <- env$ivt_scan_extended_parameter_specs("focused")
  envelope <- .ivt_refinement_envelope()
  parent_a <- .ivt_refinement_anchor("parent_a")
  parent_b <- .ivt_refinement_anchor("parent_b")
  parent_b[.ivt_refinement_parameter_names] <- list(
    1.5e-3, 2e-5, 0.08, 0.006, 2.5, 0.02, 0.002,
    0.985, 0.6, 8, 0.5, 3.5
  )
  stage1 <- rbind(parent_a, parent_b)
  stage1$candidate_id <- c("stage1_a", "stage1_b")
  stage1$hard_pass <- TRUE
  stage1$pareto_front <- TRUE
  stage1$full_target_pass <- FALSE
  stage1$target_failure_class <- c("two_n_close_four_n_fail", "four_n_close_two_n_fail")
  stage1$target_gap_score <- c(0.1, 0.2)
  stage1$objective <- c(10, 20)
  stage1$bimodality_pass <- c(TRUE, FALSE)
  stage1$wgd_mass_pass <- TRUE
  stage1$kary_2N_pass <- TRUE
  stage1$kary_4N_pass <- c(FALSE, TRUE)

  design_a <- env$ivt_scan_stage2_design(
    stage1_summary = stage1,
    specs = specs,
    n = 40,
    seed = 303,
    envelope = envelope,
    max_parents = 2,
    neighborhood_fraction = 0.05
  )
  design_b <- env$ivt_scan_stage2_design(
    stage1_summary = stage1,
    specs = specs,
    n = 40,
    seed = 303,
    envelope = envelope,
    max_parents = 2,
    neighborhood_fraction = 0.05
  )

  testthat::expect_identical(design_a, design_b)
  testthat::expect_equal(nrow(design_a), 40L)
  testthat::expect_identical(anyDuplicated(design_a$candidate_id), 0L)
  testthat::expect_true(all(design_a$stage == "stage2"))
  testthat::expect_true(all(design_a$region == "local_refinement"))
  testthat::expect_setequal(unique(design_a$parent_candidate_id), stage1$candidate_id)
  testthat::expect_true(all(is.finite(design_a$conditional_p_mis_base_upper)))
  testthat::expect_true(all(design_a$p_mis_base <= design_a$conditional_p_mis_base_upper))
  testthat::expect_equal(
    design_a$conditional_p_mis_base_upper,
    .ivt_refinement_conservative_cap(design_a$p_wgd, envelope)
  )

  for (name in .ivt_refinement_parameter_names) {
    spec_idx <- match(name, specs$param_symbol)
    testthat::expect_true(all(design_a[[name]] >= specs$lower[[spec_idx]]))
    testthat::expect_true(all(design_a[[name]] <= specs$upper[[spec_idx]]))
  }
  parent_idx <- match(design_a$parent_candidate_id, stage1$candidate_id)
  for (name in setdiff(.ivt_refinement_parameter_names, c("p_mis_base", "buffer_smax"))) {
    spec_idx <- match(name, specs$param_symbol)
    parent_value <- stage1[[name]][parent_idx]
    if (specs$transform[[spec_idx]] == "log") {
      distance <- abs(log(design_a[[name]]) - log(parent_value))
      allowed <- 0.05 * log(specs$upper[[spec_idx]] / specs$lower[[spec_idx]])
    } else {
      distance <- abs(design_a[[name]] - parent_value)
      allowed <- 0.05 * (specs$upper[[spec_idx]] - specs$lower[[spec_idx]])
    }
    testthat::expect_true(all(distance <= allowed + 1e-12))
  }

  restricted_envelope <- envelope
  restricted_envelope$eligible <- c(TRUE, TRUE, FALSE)
  restricted <- env$ivt_scan_stage2_design(
    stage1_summary = stage1[1, , drop = FALSE],
    specs = specs,
    n = 20,
    seed = 304,
    envelope = restricted_envelope,
    max_parents = 1,
    neighborhood_fraction = 0.05
  )
  testthat::expect_true(all(restricted$p_wgd <= 5e-4))
})

testthat::test_that("conditional robustness design evaluates 12 parameters in both directions", {
  env <- .load_invitro_observed_feasibility_refinement_driver()
  envelope <- .ivt_refinement_envelope()
  parent <- .ivt_refinement_anchor("parent")
  candidates <- parent[rep(1L, 11L), , drop = FALSE]
  candidates$candidate_id <- sprintf("target_%02d", seq_len(nrow(candidates)))
  candidates$full_target_pass <- TRUE
  candidates$pareto_front <- TRUE
  candidates$objective <- seq_len(nrow(candidates))
  candidates$is_anchor <- FALSE
  candidates$is_anchor[[1L]] <- TRUE

  design <- env$.ivt_ref_robustness_design(
    candidates = candidates,
    specs = env$ivt_scan_extended_parameter_specs("outer"),
    envelope = envelope,
    n_parents = 10,
    relative_change = 0.05
  )

  testthat::expect_equal(nrow(design), 240L)
  testthat::expect_identical(anyDuplicated(design$candidate_id), 0L)
  testthat::expect_false("target_01" %in% design$parent_candidate_id)
  testthat::expect_setequal(
    unique(design$perturbed_parameter),
    .ivt_refinement_parameter_names
  )
  testthat::expect_setequal(
    unique(design$perturbation_direction),
    c("minus", "plus")
  )
  testthat::expect_true(all(
    design$p_mis_base <= design$conditional_p_mis_base_upper
  ))
  testthat::expect_equal(
    design$effective_relative_change,
    design$requested_relative_change,
    tolerance = 1e-12
  )
  parent_idx <- match(design$parent_candidate_id, candidates$candidate_id)
  n_changed <- vapply(seq_len(nrow(design)), function(i) {
    sum(vapply(.ivt_refinement_parameter_names, function(name) {
      !isTRUE(all.equal(
        as.numeric(design[[name]][[i]]),
        as.numeric(candidates[[name]][[parent_idx[[i]]]]),
        tolerance = 1e-14
      ))
    }, logical(1)))
  }, integer(1))
  testthat::expect_true(all(n_changed == 1L))
})

testthat::test_that("checkpoint validation rejects stale stage-2 membership and parent metadata", {
  env <- .load_invitro_observed_feasibility_refinement_helpers()
  design <- data.frame(
    candidate_id = sprintf("stage2_%04d", 1:3),
    stage = "stage2",
    parent_candidate_id = c("stage1_a", "stage1_a", "stage1_b"),
    stringsAsFactors = FALSE
  )
  batches <- list(1:2, 3L)
  checkpoints <- list(
    list(
      schema_version = 2L,
      stage = "stage2",
      design_digest = "digest_stage2",
      batch_id = 1L,
      candidate_ids = design$candidate_id[1:2],
      summary = design[1:2, , drop = FALSE]
    ),
    list(
      schema_version = 2L,
      stage = "stage2",
      design_digest = "digest_stage2",
      batch_id = 2L,
      candidate_ids = design$candidate_id[3],
      summary = design[3, , drop = FALSE]
    )
  )

  testthat::expect_true(
    env$ivt_scan_validate_checkpoint_membership(
      checkpoints,
      design,
      batches,
      stage = "stage2",
      design_digest = "digest_stage2"
    )
  )
  stale_ids <- checkpoints
  stale_ids[[1]]$summary <- stale_ids[[1]]$summary[2:1, , drop = FALSE]
  testthat::expect_error(
    env$ivt_scan_validate_checkpoint_membership(
      stale_ids, design, batches, stage = "stage2", design_digest = "digest_stage2"
    ),
    "batch.*1|1.*batch"
  )
  stale_parent <- checkpoints
  stale_parent[[2]]$summary$parent_candidate_id <- "stage1_a"
  testthat::expect_error(
    env$ivt_scan_validate_checkpoint_membership(
      stale_parent, design, batches, stage = "stage2", design_digest = "digest_stage2"
    ),
    "parent|batch.*2|2.*batch"
  )
  testthat::expect_error(
    env$ivt_scan_validate_checkpoint_membership(
      checkpoints[1], design, batches, stage = "stage2", design_digest = "digest_stage2"
    ),
    "checkpoint|batch"
  )
  stale_stage <- checkpoints
  stale_stage[[1]]$stage <- "stage1"
  testthat::expect_error(
    env$ivt_scan_validate_checkpoint_membership(
      stale_stage, design, batches, stage = "stage2", design_digest = "digest_stage2"
    ),
    "stage"
  )
  stale_digest <- checkpoints
  stale_digest[[2]]$design_digest <- "old_design"
  testthat::expect_error(
    env$ivt_scan_validate_checkpoint_membership(
      stale_digest, design, batches, stage = "stage2", design_digest = "digest_stage2"
    ),
    "digest"
  )
  stale_schema <- checkpoints
  stale_schema[[1]]$schema_version <- 1L
  testthat::expect_error(
    env$ivt_scan_validate_checkpoint_membership(
      stale_schema, design, batches, stage = "stage2", design_digest = "digest_stage2"
    ),
    "schema"
  )
})
