# Construct fresh response-free starts before opening either scoring case.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
archive <- args[1L]
out <- args[2L]
source(file.path(archive,'release/production_bindings.R'))
bounds <- readRDS(file.path(archive,'inputs/public/contract.rds'))$design$bounds
roles <- bounds$parameter
free <- setdiff(roles,'sigma_growth')
stopifnot(length(roles)==20L, length(free)==19L, sum(roles=='sigma_growth')==1L)
letters <- LETTERS[1:6]
paths <- c(paste0('inputs/banks/D_',letters[1:4],'.rds'),paste0('verification/F_',letters[5:6],'_de4/bank.rds'))
hashes <- c('8e3ebac616f384f4bc20b9c1ce3c3db800d81efc2f13078684330678cefa2115',
 'e27665ffc68ae7a5a43dbed163e99fbba969b6d9be819d63189ea9e7582f1e44',
 '765a740842c7eab130ae7ace9305fafbc2fd03c2ce9033d0637bd6318e36c4ab',
 'f4e29c3532af47af3afc43045a6d3981e34a563e81d4452be259182d60b647b7',
 'ca56091aadad963c78590d0c5834388838ebfe4847572ee5d47867fa305f49e0',
 '4ca1eed31b60055222e50c4de9d1aec4f96ef0c4072304c4c9ef4161a8303eb1')
RNGkind("L'Ecuyer-CMRG","Inversion","Rejection")
set.seed(2026111802L)
stream <- .Random.seed
registry <- physical <- list()
for (i in seq_along(letters)) {
  path <- file.path(archive,paths[i])
  hash <- digest::digest(file=path,algo='sha256')
  stopifnot(hashes[i]==hash)
  bank <- readRDS(path)
  stopifnot(bank$origin==paste0('F_',letters[i]), length(bank$members)==200L)
  u <- do.call(rbind,lapply(bank$members,function(m) {
    stopifnot(identical(m$key,sf_hash(m$unit)),all(is.finite(m$unit)),all(m$unit>=0&m$unit<=1))
    m$unit
  }))
  stopifnot(identical(dim(u),c(200L,19L)))
  physical[[i]] <- u[, free!='sigma_kary', drop=FALSE]
  assign('.Random.seed',stream,.GlobalEnv)
  population <- matrix(NA_real_,nrow=200L,ncol=20L)
  population[,roles!='sigma_growth'] <- u
  population[,roles=='sigma_growth'] <- runif(200L)
  optimizer <- parallel::nextRNGStream(stream)
  origin <- paste0('N_',letters[i])
  sf_rds(list(origin=origin,source_sha=hash,bounds_sha=sf_hash(bounds),roles=roles,
              population=population,physical_native=lapply(bank$members,function(m)m$params[free])),
         file.path(out,'banks',paste0(origin,'.rds')))
  sf_rds(list(origin=origin,state=optimizer,RNG_sha=sf_hash(optimizer)),
         file.path(out,'rng',paste0(origin,'.rds')))
  registry[[i]] <- data.frame(origin=origin,source_path=paths[i],source_sha=hash,physical_order_unchanged=TRUE,
       population_sha=sf_hash(population),growth_stream_sha=sf_hash(stream),optimizer_stream_sha=sf_hash(optimizer),
       growth_unit_min=min(population[,roles=='sigma_growth']),growth_unit_max=max(population[,roles=='sigma_growth']))
  stream <- parallel::nextRNGStream(optimizer)
  rm(bank)
  invisible(gc())
}
stopifnot(!anyDuplicated(as.data.frame(do.call(rbind,physical))))
sf_emit(do.call(rbind,registry),file.path(out,'start_registry.tsv'))
cat('N20_STARTS_PASS six response-free populations with 1,200 distinct physical rows\n')
