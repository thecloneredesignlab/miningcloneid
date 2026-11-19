# New forward-only adapter. No optimizer, synthetic observation or legacy counter.
sf_tab <- function(p) read.delim(p,stringsAsFactors=FALSE,check.names=FALSE)
sf_emit <- function(x,p) {
  dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE)
  old<-options(digits=17);on.exit(options(old));q<-tempfile('.pending_',dirname(p))
  write.table(x,q,sep='\t',quote=TRUE,row.names=FALSE,na='NA');stopifnot(file.rename(q,p))
}
sf_rds <- function(x,p) {dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);q<-tempfile('.pending_',dirname(p));saveRDS(x,q);stopifnot(file.rename(q,p))}
sf_hash <- function(x) digest::digest(x,algo='sha256')
sf_inventory <- function(out) {
  ns<-sort(loadedNamespaces());stopifnot(!length(intersect(ns,c('rJava','V8','jsonlite','rjson','RJSONIO'))))
  sf_emit(data.frame(namespace=ns,version=vapply(ns,function(n)as.character(getNamespaceVersion(n)),character(1))),file.path(out,'loaded_namespaces.tsv'))
}
sf_backend <- function(study,out,native=FALSE) {
  stopifnot(getRversion()=='4.4.2',startsWith(normalizePath(tempdir()),paste0(normalizePath(file.path(out,'tmp')),'/')))
  suppressPackageStartupMessages(library(digest));suppressPackageStartupMessages(library(dplyr))
  env<-new.env(parent=.GlobalEnv);records<-list();src<-file.path(study,'inputs/production')
  files<-c('o2_supply_demand_map_shared.R','o2_supply_demand_map_common_semantics.R','model_O2_supply_demand_MAP.R',
    'o2_supply_demand_map_invitro_io_utils.R','o2_supply_demand_map_invitro_lineage_utils.R','o2_supply_demand_map_invitro_summary_utils.R',
    'o2_supply_demand_map_invitro_lineage_simulation_utils.R','o2_supply_demand_map_invitro_objective_utils.R','parameter_helpers.R')
  constants<-c('.chr_lengths_default_bp_1to22','.chr_lengths_default_ord_desc','.chr_lengths_default_denom','INVITRO_PASSAGE_IMPLEMENTATION','o2prof_symbols')
  for(file in files)for(e in parse(file.path(src,file)))if(is.call(e)&&identical(e[[1]],as.name('<-'))&&is.symbol(e[[2]])) {
    name<-as.character(e[[2]]);rhs<-e[[3]]
    func<-is.call(rhs)&&identical(rhs[[1]],as.name('function'))
    closure<-name%in%c('ivt_load_fit_objects','ivt_run_segment_fixed_o2')
    if(func||closure||name%in%constants) {
      eval(e,envir=env);records[[length(records)+1L]]<-data.frame(file=file,name=name,kind=if(func)'unchanged_function' else if(closure)'unchanged_pure_closure_definition' else 'unchanged_constant')
    }
  }
  env$.first_non_null_local<-env$o2sd_first_non_null
  stopifnot(identical(env$INVITRO_PASSAGE_IMPLEMENTATION,'v2'),!exists('.init_cpp_o2simps_backend',env,inherits=FALSE))
  sf_emit(do.call(rbind,records),file.path(out,'loaded_source_definitions.tsv'))
  if(native) {
    cache<-file.path(out,'tmp/native_cache');dir.create(cache,recursive=TRUE)
    Rcpp::sourceCpp(file.path(src,'model_O2_supply_demand_MAP.cpp'),env=env,cacheDir=cache,rebuild=FALSE,showOutput=FALSE,verbose=FALSE)
    stopifnot(exists('cpp_o2simps_simulate_one',env,inherits=FALSE));env$.USE_CPP_O2SIMPS_BACKEND<-TRUE
    sf_emit(data.frame(export=c('cpp_o2simps_simulate_one','cpp_o2simps_build_G_for_o2_triplet'),loading='direct sourceCpp into private cache; original timed bootstrap not evaluated'),file.path(out,'native_interface.tsv'))
  }
  env
}
