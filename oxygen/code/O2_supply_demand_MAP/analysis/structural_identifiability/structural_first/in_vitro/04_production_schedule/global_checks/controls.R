source("oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step03/backend.R")
cx <- branch_context()
stopifnot(identical(as.numeric(cx$oxygens), c(0,.1,.2,.3,.5,1,2,20.5)))
a <- cx$anchor
anchor <- cx$evaluate(a,keep_details=TRUE,exact_kary=TRUE)
stopifnot(anchor$ok)
prior <- readRDS(file.path(bridge,"verification/R01/result.rds"))
stopifnot(max(abs(anchor$growth-prior$pred$growth[names(anchor$growth)]))<1e-11,
  max(abs(anchor$kary-do.call(cbind,prior$pred$kary[colnames(anchor$kary)])))<1e-11,
  identical(unname(anchor$days),prior$pred$selection$day),
  identical(unname(anchor$parent_days),prior$pred$selection$parent_day),
  !any(anchor$fallback))
hill <- function(O,C,n) C^n/(C^n+O^n)
parameters <- list(anchor=a)
data <- list()
check <- function(id,p,q,O=NULL) {
  x <- cx$evaluate(p, O, keep_details=TRUE, exact_kary=TRUE)
  y <- cx$evaluate(q, O, keep_details=TRUE, exact_kary=TRUE)
  stopifnot(x$ok,y$ok)
  state <- max(vapply(seq_along(x$details),function(i) {
    v <- x$details[[i]]$states; w <- y$details[[i]]$states
    max(abs(v-w)/pmax(1,abs(v),abs(w)))
  },0))
  d <- differences(x,y)
  stopifnot(d["growth_max"]<1e-11,d["kary_tv_max"]<1e-11,
            d["endpoint"]==0,d["mask"]==0,d["parent_day"]==0,d["reset_mode"]==0)
  list(row=data.frame(control=id,O=if(is.null(O))"actual"else as.character(O),
                      growth_bitwise=identical(x$growth,y$growth),kary_bitwise=identical(x$kary,y$kary),
                      scaled_state_gap=state,as.list(d),check.names=FALSE),a=x,b=y)
}
o <- 20.5
alt <- a
alt$O2_crit <- a$O2_crit*0.99
alt$n_O <- a$n_O+0.01
q <- hill(o,a$O2_crit,a$n_O)/hill(o,alt$O2_crit,alt$n_O)
alt$alpha_o2 <- a$alpha_o2*q
alt$mu_hp <- a$mu_hp*q
stopifnot(all(cx$encode(alt)>=-1e-14&cx$encode(alt)<=1+1e-14))
single <- check("single_oxygen_fixed_hill_products",a,alt,o)
parameters$single_oxygen_analytic_member_rounded <- alt
data[[1L]] <- single$row
cat("Single oxygen structural control passed native replay\n")

wid <- a; wid$sigma_kary <- .0142
plateau <- check("multi_oxygen_binary64_width_plateau",a,wid)
stopifnot(plateau$row$growth_bitwise,plateau$row$kary_bitwise)
parameters$width_plateau <- wid
data[[2L]] <- plateau$row

s0 <- a; s0$buffer_smax <- 0
s0$p_misseg <- cx$bounds$lower[match("p_misseg",cx$keep)]
s0$k_o_mis <- cx$bounds$upper[match("k_o_mis",cx$keep)]
z <- s0; z$buffer_beta <- 5; z$buffer_n_exp <- 6
face <- check("multi_oxygen_zero_maximum_survival",s0,z)
stopifnot(face$row$growth_bitwise,face$row$kary_bitwise)
parameters$survival_face_a <- s0; parameters$survival_face_b <- z
data[[3L]] <- face$row

sf_emit(do.call(rbind,data),file.path(out,"controls.tsv"))
emit_parameters(parameters,file.path(out,"control_parameters.tsv"),cx$keep)
sf_rds(list(anchor=anchor,parameters=parameters),file.path(out,"controls.rds"))
cat("Controls passed. State and output replays used the original production C++.\n")
