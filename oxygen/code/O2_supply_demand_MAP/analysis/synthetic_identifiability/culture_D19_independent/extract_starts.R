# Offline start extraction. No observation case or generating vector is opened here.
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L)
archive <- args[1L]
out <- args[2L]
source(file.path(archive, 'release/production_bindings.R'))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
bank_paths <- c(A = 'inputs/banks/D_A.rds', B = 'inputs/banks/D_B.rds',
                D = 'inputs/banks/D_D.rds', F = 'verification/F_F_de4/bank.rds')
expected <- c(A = '8e3ebac616f384f4bc20b9c1ce3c3db800d81efc2f13078684330678cefa2115',
 B = 'e27665ffc68ae7a5a43dbed163e99fbba969b6d9be819d63189ea9e7582f1e44',
 D = 'f4e29c3532af47af3afc43045a6d3981e34a563e81d4452be259182d60b647b7',
 F = '4ca1eed31b60055222e50c4de9d1aec4f96ef0c4072304c4c9ef4161a8303eb1')
rows <- physical <- list()
roles <- readRDS(file.path(archive, 'inputs/public/contract.rds'))$design$bounds$parameter
free <- setdiff(roles, 'sigma_growth')
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
set.seed(2026111801L)
seed <- .Random.seed
for (letter in names(bank_paths)) {
  origin <- paste0('D_', letter)
  input <- file.path(archive, bank_paths[[letter]])
  hash <- digest::digest(file = input, algo = 'sha256')
  stopifnot(hash == expected[[letter]])
  bank <- readRDS(input)
  u <- do.call(rbind, lapply(bank$members, function(m) {
    stopifnot(identical(sf_hash(m$unit), m$key), all(is.finite(m$unit)), all(m$unit >= 0 & m$unit <= 1))
    m$unit
  }))
  requested <- readRDS(file.path(archive, 'verification', paste0(origin, '_bench/requested_population.rds')))
  actual <- readRDS(file.path(archive, 'verification', paste0(origin, '_bench/actual_initial_population.rds')))
  old_seed <- readRDS(file.path(archive, 'verification/T00/rng', paste0(origin, '.rds')))
  stopifnot(bank$origin == paste0('F_', letter), identical(dim(u), c(200L, 19L)),
            identical(u, requested$population), identical(unname(u), unname(actual)),
            !identical(seed, old_seed$state))
  record <- list(origin = origin, source_sha = hash, contract_sha = bank$contract_sha,
                 population = u, natural = lapply(bank$members, function(m) m$params[free]))
  sf_rds(record, file.path(out, 'banks', paste0(origin, '.rds')))
  sf_rds(list(origin = origin, state = seed, RNG_sha = sf_hash(seed)), file.path(out, 'rng', paste0(origin, '.rds')))
  rows[[origin]] <- data.frame(origin = origin, source_path = bank_paths[[letter]], source_sha = hash,
    historical_requested_order = TRUE, historical_actual_order = TRUE, members = nrow(u),
    coordinates = ncol(u), population_sha = sf_hash(u), new_rng_state = paste(seed, collapse = ','),
    new_rng_sha = sf_hash(seed), historical_rng_sha = sf_hash(old_seed$state))
  physical[[origin]] <- u[, free != 'sigma_kary', drop = FALSE]
  seed <- parallel::nextRNGStream(seed)
  rm(bank)
  invisible(gc())
}
stopifnot(!anyDuplicated(as.data.frame(do.call(rbind, physical))))
sf_emit(do.call(rbind, rows), file.path(out, 'start_registry.tsv'))
cat('START_EXTRACTION_PASS A/B/D/F, 800 physically distinct rows\n')
