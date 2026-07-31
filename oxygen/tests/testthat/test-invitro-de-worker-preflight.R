.load_invitro_de_preflight_api <- function() {
  backend_path <- file.path(
    repo_info$root,
    "oxygen",
    "code",
    "O2_supply_demand_MAP",
    "util",
    "o2_supply_demand_map_fit_invitro_backend.R"
  )
  targets <- c(
    "INVITRO_DE_PREFLIGHT_MAX_RETRIES",
    "INVITRO_DE_PENALTY_OBJECTIVE",
    "invitro_de_preflight_text",
    "write_invitro_de_preflight_audit",
    "initialize_invitro_deoptim_workers",
    "invitro_de_preflight_evaluate",
    "start_invitro_deoptim_cluster_with_preflight"
  )
  env <- new.env(parent = globalenv())
  expressions <- parse(backend_path)
  for (expr in expressions) {
    if (!is.call(expr) || !identical(expr[[1L]], as.name("<-"))) next
    target <- as.character(expr[[2L]])
    if (target %in% targets) eval(expr, envir = env)
  }
  missing <- targets[!vapply(targets, exists, logical(1), envir = env, inherits = FALSE)]
  if (length(missing)) {
    stop("Missing in-vitro DE preflight definitions: ", paste(missing, collapse = ", "))
  }
  env
}

.preflight_row <- function(ok,
                           status,
                           error = NA_character_,
                           objective = if (ok) 10 else 1e9,
                           penalty_reason = NA_character_) {
  data.frame(
    worker = 1L,
    role = "worker",
    ok = ok,
    status = status,
    objective = objective,
    penalty_reason = penalty_reason,
    error = error,
    stringsAsFactors = FALSE
  )
}

testthat::test_that("in-vitro DE preflight classifies objective failures", {
  env <- .load_invitro_de_preflight_api()

  pass <- env$invitro_de_preflight_evaluate(
    function(par) list(objective = sum(par)),
    c(1, 2)
  )
  testthat::expect_true(pass$ok)
  testthat::expect_identical(pass$status, "PASS")
  testthat::expect_equal(pass$objective, 3)

  penalty <- env$invitro_de_preflight_evaluate(
    function(par) list(
      objective = 1e9,
      penalty_reason = "simulation_error: stale worker backend"
    ),
    1
  )
  testthat::expect_false(penalty$ok)
  testthat::expect_identical(penalty$status, "PENALTY_OBJECTIVE")
  testthat::expect_match(penalty$penalty_reason, "stale worker backend")

  failure <- env$invitro_de_preflight_evaluate(
    function(par) stop("worker DLL unavailable"),
    1
  )
  testthat::expect_false(failure$ok)
  testthat::expect_identical(failure$status, "OBJECTIVE_ERROR")
  testthat::expect_match(failure$error, "worker DLL unavailable")
})

testthat::test_that("in-vitro DE preflight retries with fresh clusters until success", {
  env <- .load_invitro_de_preflight_api()
  audit_path <- tempfile("invitro_de_preflight_", fileext = ".tsv")
  on.exit(unlink(audit_path, force = TRUE), add = TRUE)
  starts <- 0L
  stopped <- integer(0)

  result <- env$start_invitro_deoptim_cluster_with_preflight(
    n_cores = 8L,
    objective_from_free = function(par) list(objective = 10),
    objective_value = function(par) 10,
    init_free = 1,
    cpp_info = list(wrapper_path = "wrapper", path = "dll"),
    audit_path = audit_path,
    max_retries = 5L,
    cluster_factory = function(n_cores) {
      starts <<- starts + 1L
      list(starts)
    },
    cluster_stopper = function(cluster) {
      stopped <<- c(stopped, cluster[[1L]])
    },
    worker_initializer = function(cluster, cpp_info) {
      if (cluster[[1L]] < 3L) {
        .preflight_row(
          ok = FALSE,
          status = "WORKER_CPP_INIT_ERROR",
          error = paste0("attempt ", cluster[[1L]], " could not load DLL")
        )
      } else {
        .preflight_row(ok = TRUE, status = "PASS_WRAPPER_DLL_LOADED")
      }
    },
    preflight_runner = function(cluster, objective_from_free, objective_value, init_free) {
      .preflight_row(ok = TRUE, status = "PASS")
    },
    sleep_fn = function(seconds) invisible(NULL)
  )

  testthat::expect_identical(starts, 3L)
  testthat::expect_identical(stopped, c(1L, 2L))
  testthat::expect_identical(result$attempts, 3L)
  testthat::expect_identical(result$retries, 2L)
  testthat::expect_identical(result$cluster[[1L]], 3L)
  audit <- utils::read.delim(audit_path, stringsAsFactors = FALSE)
  testthat::expect_identical(audit$retry_number, c(0L, 1L, 2L, 2L))
  testthat::expect_identical(audit$status, c(
    "WORKER_CPP_INIT_ERROR",
    "WORKER_CPP_INIT_ERROR",
    "PASS_WRAPPER_DLL_LOADED",
    "PASS"
  ))
})

testthat::test_that("in-vitro DE preflight stops and records errors after five retries", {
  env <- .load_invitro_de_preflight_api()
  audit_path <- tempfile("invitro_de_preflight_exhausted_", fileext = ".tsv")
  on.exit(unlink(audit_path, force = TRUE), add = TRUE)
  starts <- 0L
  stops <- 0L

  testthat::expect_error(
    env$start_invitro_deoptim_cluster_with_preflight(
      n_cores = 8L,
      objective_from_free = function(par) list(objective = 10),
      objective_value = function(par) 10,
      init_free = 1,
      cpp_info = list(wrapper_path = "wrapper", path = "dll"),
      audit_path = audit_path,
      max_retries = 5L,
      cluster_factory = function(n_cores) {
        starts <<- starts + 1L
        list(starts)
      },
      cluster_stopper = function(cluster) {
        stops <<- stops + 1L
      },
      worker_initializer = function(cluster, cpp_info) {
        .preflight_row(
          ok = FALSE,
          status = "WORKER_CPP_INIT_ERROR",
          error = paste0("worker DLL initialization failed on attempt ", cluster[[1L]])
        )
      },
      preflight_runner = function(cluster, objective_from_free, objective_value, init_free) {
        stop("preflight must not run when worker initialization fails")
      },
      sleep_fn = function(seconds) invisible(NULL)
    ),
    "failed after 6 attempts \\(5 retries\\)"
  )

  testthat::expect_identical(starts, 6L)
  testthat::expect_identical(stops, 6L)
  audit <- utils::read.delim(audit_path, stringsAsFactors = FALSE)
  testthat::expect_identical(audit$attempt, 1:6)
  testthat::expect_identical(audit$retry_number, 0:5)
  testthat::expect_true(all(audit$status == "WORKER_CPP_INIT_ERROR"))
  testthat::expect_match(tail(audit$error, 1L), "attempt 6")
})
