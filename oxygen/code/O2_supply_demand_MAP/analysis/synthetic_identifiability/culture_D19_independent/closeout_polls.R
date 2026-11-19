# Offline local-poll audit. A certificate stays with the exact anchor that was checked.
tsv <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

bound_step <- function(x, d, h) {
  crossing <- c((1 - x[d > 0]) / d[d > 0], -x[d < 0] / d[d < 0])
  step <- min(h, crossing)
  if (!is.finite(step) || step <= 0) return(NULL)
  u <- as.numeric(pmin(1, pmax(0, x + step*d)))
  if (identical(u, x)) NULL else u
}

poll_keys <- function(anchor, coordinates) {
  x <- anchor$unit
  stopifnot(length(x) == length(coordinates), length(x) == 19L, !anyDuplicated(coordinates))
  direction <- diag(length(x))
  pairs <- list(c('lam_max','alpha_o2',1), c('mu_hp','k_o_mis',1),
    c('p_misseg','p_mis_base',-1), c('buffer_smax','buffer_beta',1),
    c('init_sd_2N','sigma_kary',-1), c('init_sd_4N','sigma_kary',-1))
  for (pair in pairs) {
    d <- numeric(length(x))
    d[match(pair[1:2], coordinates)] <- c(1, as.numeric(pair[3L])) / sqrt(2)
    direction <- rbind(direction, d)
  }
  expected <- character()
  for (h in c(.001,.0001,.00002)) {
    for (j in seq_len(nrow(direction))) {
      for (sign in c(-1,1)) {
        u <- bound_step(x, sign*direction[j,], h)
        if (!is.null(u)) expected <- c(expected, digest::digest(u, algo = 'sha256'))
      }
    }
  }
  unique(expected)
}

audit_polls <- function(root, out) {
  case <- readRDS(file.path(root, 'inputs/public/D19.rds'))
  coordinates <- setdiff(case$design$bounds$parameter, 'sigma_growth')
  quality <- tsv(file.path(root, 'reports/verified_candidates.tsv'))
  selection <- tsv(file.path(root, 'reports/selection_lock.tsv'))
  rows <- list()
  for (origin in unique(quality$origin)) {
    final <- readRDS(file.path(root, 'verification', paste0(origin, '_replay2b'), 'state.rds'))
    stopifnot(final$origin == origin)
    for (phase in paste0(origin, c('_verify1','_verify2'))) {
      path <- file.path(root, 'verification', phase)
      anchor <- readRDS(file.path(path, 'checked_representative.rds'))
      polls <- tsv(file.path(path, 'verification_polls.tsv'))
      key <- anchor$key
      expected <- poll_keys(anchor, coordinates)
      cert <- final$local_certificates[[key]]
      same_final_poll <- identical(cert$task, phase)
      complete <- identical(as.character(polls$key), expected)
      stopifnot(anchor$origin == origin, all(polls$anchor_key == key), !anyDuplicated(polls$key),
        !anyNA(polls$score), identical(digest::digest(anchor$unit, algo = 'sha256'), key),
        all(polls$key %in% expected))
      gain <- max(0, anchor$score - polls$score)
      if (same_final_poll) {
        stopifnot(complete == cert$complete, abs(gain - cert$max_gain) <= 1e-12)
        rec <- quality[quality$origin == origin & quality$key == key, , drop = FALSE]
        stopifnot(nrow(rec) == 1L, rec$checked, rec$local_pass == (complete && gain <= 1e-5))
      }
      rows[[length(rows) + 1L]] <- data.frame(origin = origin, phase = phase, key = key,
        anchor_D = anchor$score, expected_directions = length(expected), evaluated_directions = nrow(polls),
        complete = complete, max_D_improvement = gain, local_pass = complete && gain <= 1e-5,
        associated_final_poll = same_final_poll,
        is_global_lowest_found = origin == selection$best_found_origin && key == selection$best_found_key)
    }
    if (origin == selection$best_found_origin) {
      stopifnot(identical(final$best$key, selection$best_found_key),
        abs(final$best$score - selection$best_found_D) <= 1e-12)
    }
  }
  result <- do.call(rbind, rows)
  matches <- result$is_global_lowest_found & result$associated_final_poll
  stopifnot(sum(matches) <= 1L,
    isTRUE(selection$best_found_local_pass) == any(result$local_pass[matches]))
  old <- options(digits = 17)
  on.exit(options(old))
  write.table(result, file.path(out, 'd19_local_poll_audit.tsv'), sep = '\t', row.names = FALSE, quote = TRUE)
  cat('OFFLINE_EXACT_VECTOR_LOCAL_POLL_AUDIT_PASS\n')
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  stopifnot(length(args) == 2L)
  audit_polls(args[1L], args[2L])
}
