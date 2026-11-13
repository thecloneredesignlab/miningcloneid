# Actual native algorithms, typed I/O, exact-zero roles and stale-checkpoint tests.
out <- file.path(study, "tests/unit"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
ctx <- o2prof_load_culture(file.path(project, settings$wp0_relative_path), modules)
m <- original_map(ctx); base <- o2prof_project_start(m, ctx$saved$best_params)$params
stopifnot(sum(m$table$status == "active") == 20)
positive <- original_map(ctx, list(p_misseg = .01)); stopifnot(sum(positive$table$status == "active") == 19)
for (fixed in list(list(p_misseg = 0), list(buffer_beta = 0))) {
  z <- original_map(ctx, fixed); p <- o2prof_project_start(z, base)$params
  stopifnot(sum(z$table$status == "active") == 18L)
  for (n in names(fixed)) stopifnot(p[[n]] == 0, o2prof_from_free(z, o2prof_to_free(z, p))[[n]] == 0)
}
pen <- o2prof_map(ctx$cfg$parameter_table, ctx$saved$best_params, "zero_extended", fixed = list(p_misseg = 0), penalty_active = "k_o_mis")
stopifnot(pen$table$status[pen$table$parameter == "k_o_mis"] == "active")
x <- list(zero = 0, negative = -1.2345678901234567, missing = NA_real_, infinite = Inf, false = FALSE, integer = 12L, label = "Quoted field with a tab\tand a quote\"")
write_kv(x, file.path(out, "typed_roundtrip.tsv")); stopifnot(identical(read_kv(file.path(out, "typed_roundtrip.tsv")), x))
saveRDS(base, file.path(out, "parameter_roundtrip.rds")); stopifnot(identical(base, readRDS(file.path(out, "parameter_roundtrip.rds"))))
free <- o2prof_to_free(m, base); o2prof_write_tsv(data.frame(coordinate = names(free), value = unname(free)), file.path(out, "free_roundtrip.tsv"))
stopifnot(identical(as.numeric(free), read_table(file.path(out, "free_roundtrip.tsv"))$value))
fn <- function(x) if (x[1] < .1) Inf else abs(x[1] - .4) + (x[2] - 1)^2
s <- nloptr::nloptr(c(.6, .6), fn, lb = c(0, 0), ub = c(1, 1), opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1000L, xtol_abs = c(1e-7, 1e-7), ftol_abs = 1e-10))
stopifnot(s$objective < 1e-5, all(s$solution >= 0), all(s$solution <= 1))
p <- o2review_pattern(fn, s$solution, diag(2), max_evaluations = 1000L)
stopifnot(p$tested_directional_stability, p$value <= s$objective + 1e-8)
# All-finite constant objective exercises native patience at 300 generations.
d <- DEoptim::DEoptim(function(x) 1, lower = c(0, 0), upper = c(1, 1), control = list(NP = 20L, itermax = 305L, reltol = 1e-7, steptol = 300L, trace = FALSE))
stopifnot(d$optim$iter >= 300L)
cp <- file.path(out, "checkpoint"); o2pilot_checkpoint(base, 4, cp, "valid")
o2pilot_checkpoint(base, 3, cp, "valid"); o2pilot_checkpoint(base, 2, cp, "valid")
stopifnot(!o2pilot_checkpoint(base, 2, cp, "valid"), length(list.files(file.path(cp, "improvements"))) == 3L)
old_records <- lapply(list.files(file.path(cp, "improvements"), full.names = TRUE), readRDS)
stopifnot(identical(vapply(old_records, function(z) z$objective, numeric(1)), c(4, 3, 2)))
# Real component/parameter/day vectors round-trip in an append-only history.
archive_out <- file.path(out, "real_improvement_archive")
p_bad <- base; p_bad$sigma_growth <- m$table$upper[m$table$parameter == "sigma_growth"]
c_good <- ctx$raw_objective(base, ctx$objects, ctx$cfg); c_bad <- ctx$raw_objective(p_bad, ctx$objects, ctx$cfg)
stopifnot(c_bad$objective > c_good$objective)
o2pilot_checkpoint(p_bad, c_bad$objective, file.path(archive_out, "best_candidates"), "archive_unit", "fixture", c_bad, 1L)
o2pilot_checkpoint(base, c_good$objective, file.path(archive_out, "best_candidates"), "archive_unit", "fixture", c_good, 2L)
archive_index <- audit_improvements(ctx, m, archive_out, "archive_unit", TRUE); stopifnot(nrow(archive_index) == 2L)
o2prof_install_culture_map(ctx, m, list(anchor = base))
fatal <- tryCatch(DEoptim::DEoptim(function(x) {
  tryCatch(stop("checkpoint_write_failed: injected unit failure"), error = function(e)
    ctx$backend$make_penalty_components(reason = paste0("simulation_error: ", conditionMessage(e))))
}, lower = c(0, 0), upper = c(1, 1), control = list(NP = 20L, itermax = 1L, trace = FALSE)), error = identity)
stopifnot(inherits(fatal, "error"), grepl("checkpoint_write_failed", conditionMessage(fatal)),
  ctx$backend$make_penalty_components(reason = "ordinary simulation infeasibility")$objective == 1e9)
stopifnot(inherits(tryCatch({o2pilot_checkpoint(base, 3, cp, "stale"); NULL}, error = identity), "error"))
# Preserve native production global feasibility guidance; do not introduce a new penalty.
stopifnot(ctx$backend$INVITRO_DE_PENALTY_OBJECTIVE == 1e9)
protocol_error <- NULL
set.seed(725991)
for (i in seq_len(32)) {
  p <- o2review_draw(m, "full_domain", i)
  e <- tryCatch({ctx$raw_objective(p, ctx$objects, ctx$cfg); NULL}, error = identity)
  if (inherits(e, "invitro_protocol_infeasible")) {protocol_error <- e; break}
}
stopifnot(inherits(protocol_error, "invitro_protocol_infeasible"),
  is.finite(ctx$backend$invitro_protocol_penalty_objective(protocol_error)))
ev <- o2review_evaluator(ctx, m, file.path(out, "production_infeasible_probe"), "unit_probe", "validation")
stopifnot(is.infinite(ev$fn(o2review_unit_spec(m, p)$unit)), ev$counts$failed == 1L)
write_kv(list(local_infeasible_value = Inf, global_guidance = "Unmodified production graded protocol penalty",
  protocol_error = conditionMessage(protocol_error)), file.path(out, "infeasibility_semantics.tsv"))
assert_runtime(out)
write_kv(list(append_only_parameter_component_archives = TRUE, immutable_records = 3L,
  real_component_records = nrow(archive_index), fatal_checkpoint_guard = TRUE), file.path(out, "improvement_archive_validation.tsv"))
write_kv(list(pass = TRUE, de_minimum_actual = as.integer(d$optim$iter), nloptr_version = as.character(packageVersion("nloptr")),
  sbplx_inf_compatible = TRUE, pattern_boundary_test = TRUE, exact_zero_mapping = TRUE, penalty_roles = TRUE), file.path(out, "result.tsv"))

# Every canonical positive target, exact null and fail-closed dispatcher.
registry <- read_table(file.path(study, "cases.tsv"))
for (i in seq_len(nrow(registry))) {
  t <- list(case = "fixed_target", point = registry$point[i], target = registry$target[i], value = registry$value[i])
  fixed <- fixed_for(t, study); map <- original_map(ctx, fixed)
  p <- o2prof_project_start(map, base)$params
  stopifnot(p[[t$target]] == t$value, sum(map$table$status == "active") == if(t$value==0) 18L else 19L)
  t$case <- "unknown"
  stopifnot(inherits(tryCatch(fixed_for(t,study),error=identity),"error"))
}
bad <- list(case="unrestricted",point="pm_05",target="none",value=NA_real_)
stopifnot(inherits(tryCatch(fixed_for(bad,study),error=identity),"error"))
for(target in c(1,2,3)) {
  fn <- function(u) (target+u-3)^2+(u-1)^2
  fit <- nloptr::nloptr(1,fn,lb=0,ub=3,opts=list(algorithm="NLOPT_LN_SBPLX",maxeval=1000L,xtol_abs=1e-8,ftol_abs=1e-12))
  stopifnot(abs(fit$solution-(2-target/2))<1e-4,abs(fit$objective-(target-2)^2/2)<1e-7)
}
imp <- read_table(file.path(study,"imports.tsv"))
path_for <- function(id) file.path(study,imp$relative_path[imp$import_id==id])
relative_for <- function(id) substring(path_for(id),nchar(project)+2L)
neighbor <- list(phase="neighbor",source_fit=relative_for("historical_pm_06"),fallback_source=relative_for("historical_pm_05"))
positive <- original_map(ctx,list(p_misseg=registry$value[registry$point=="pm_05"]))
actual <- prepare_local_start(ctx,positive,neighbor,project,file.path(out,"real_neighbor_transfer"))
stopifnot(!actual$skipped,actual$fit$best_params$p_misseg==positive$fixed$p_misseg)
source_raw <- ctx$raw_objective; probe_count <- 0L
infeasible <- function() stop(structure(list(message="Injected typed protocol infeasibility",call=NULL),class=c("invitro_protocol_infeasible","error","condition")))
neighbor$source_fit <- relative_for("historical_pm_05")
ctx$raw_objective <- function(...) {probe_count <<- probe_count+1L;if(probe_count==1L) infeasible();source_raw(...)}
mixed <- prepare_local_start(ctx,positive,neighbor,project,file.path(out,"mixed_neighbor_transfer"))
stopifnot(!mixed$skipped,tail(read_table(file.path(out,"mixed_neighbor_transfer/transfer_attempts.tsv"))$source_fraction,1)==.75)
ctx$raw_objective <- function(...) infeasible()
skipped <- prepare_local_start(ctx,positive,neighbor,project,file.path(out,"skipped_neighbor_transfer"))
stopifnot(skipped$skipped,nrow(read_table(file.path(out,"skipped_neighbor_transfer/transfer_attempts.tsv")))==4L)
ctx$raw_objective <- source_raw
source(file.path(code,"diagnostics.R"))
stopifnot(pilotdiag_offspring_identity(44,0,0)$biological_viable_daughters_per_ordinary_division==2)
invisible(pilotdiag_rates(ctx,base,ctx$cfg,2))
survival <- readRDS(path_for("historical_bb_00"))$best_params
invisible(pilotdiag_rates(ctx,survival,ctx$cfg,2))
for(r in imp$relative_path) audit_native_strings(readRDS(file.path(study,r)))
write_kv(list(pass=TRUE,canonical_targets=7L,unknown_case_rejected=TRUE,analytic_nuisance_refit=TRUE,
  real_neighbor=TRUE,mixed_neighbor=TRUE,infeasible_neighbor_skip=TRUE,production_kernel_ledgers=TRUE,
  clean_imports=12L),file.path(out,"profile_adapter_validation.tsv"))
