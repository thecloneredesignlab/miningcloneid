# Allow the one approved reserve global to target a missing initialization family.
# Numerical source, score, bounds, caps and native optimizer controls are immutable.
apply_reserve_override <- function(spec,context) {
 fields <- grep("^override_",names(context),value=TRUE)
 if(!length(fields)) return(spec)
 stopifnot(spec$role=="reserve_baseline_global",spec$mode=="global",is.null(spec$fixed),
  identical(sort(fields),sort(c("override_init_mode","override_start_seed"))))
 mode <- context$override_init_mode
 stopifnot(mode%in%c("hybrid","uniform"),grepl("^[0-9]+$",context$override_start_seed))
 seed <- as.integer(context$override_start_seed)
 stopifnot(is.finite(seed))
 if(mode=="uniform") stopifnot(seed==0L)
 else stopifnot(seed%in%c(25L,74L,366L,292L,392L))
 spec$init_mode <- mode;spec$family <- mode;spec$start_seed <- seed
 spec
}
