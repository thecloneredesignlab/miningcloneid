# Archived-state and import regression. No propagation or optimization.
args <- commandArgs(TRUE)
root <- normalizePath(args[1L])
out <- normalizePath(args[2L])
for (f in c('production_bindings.R','interface.R','feasibility.R','search_helpers.R','n20_model.R','core.R','fit.R'))
  source(file.path(root,'release',f))
env <- sf_backend(root,out,native=FALSE)
a <- n_data(root,'N_A')
ctx <- n_context(env,a$case,a$contract,a$id)
x <- readRDS(file.path(root,'verification/N_A_bench/state.rds'))
n_validate_state(ctx,x,'N_A')
stopifnot(identical(dim(x$population),c(200L,20L)),x$counts[['global']]==400L,x$calls==400L,
          x$native_evaluations==400L,x$generations==1L,
          abs(n_independent_loss(ctx,x$best$pred,x$best$params)[['value']]-x$best$score)<=1e-8)
decode_error <- numeric(length(n_origins))
for (i in seq_along(n_origins)) {
  origin <- n_origins[i]
  input <- n_data(root,origin)
  bank <- n_bank(root,origin)
  rng <- readRDS(file.path(root,'verification/T00/rng',paste0(origin,'.rds')))
  stopifnot(identical(dim(bank$population),c(200L,20L)),
            identical(sf_hash(bank$population),input$id$population_sha),
            identical(sf_hash(rng$state),input$id$rng_initial_sha),
            identical(bank$roles,input$case$design$bounds$parameter))
  context <- n_context(sf_backend(root,out,native=FALSE),input$case,input$contract,input$id)
  errors <- vapply(seq_len(200L),function(j) {
    decoded <- n_decode(context,as.numeric(bank$population[j,]))
    expected <- unlist(bank$physical_native[[j]])
    values <- unlist(decoded[names(expected)])
    stopifnot(identical(names(values),names(expected)),decoded$sigma_growth>0)
    max(abs(values-expected)/pmax(1,abs(expected)))
  },numeric(1))
  # The offline host can differ from the production image by rounding in transforms.
  # Every actual worker import still checks exact native coordinates in the SIF.
  decode_error[i] <- max(errors)
  stopifnot(decode_error[i]<=64*.Machine$double.eps)
}
cat('Maximum offline transformed-coordinate decode discrepancy',max(decode_error),'\n')
fixture <- file.path(out,'negative_Q')
dir.create(fixture)
score <- n_evaluator(ctx,fixture,'N_A','local',backend=function(u)
  list(feasible=TRUE,loss=list(value=-.1,growth=.6,kary=-.5),pred=x$best$pred))
stopifnot(score$fn(x$best$unit)==-.1,score$state$counts['local']==1L)
score$close()
cat('N20_A400_HANDOFF_AND_SIX_IMPORTS_PASS; signed Q fixture passed; zero real forwards\n')
