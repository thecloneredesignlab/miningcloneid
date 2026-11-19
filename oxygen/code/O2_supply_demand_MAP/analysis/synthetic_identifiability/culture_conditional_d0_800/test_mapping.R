# Metadata-only mapping checks; no simulator is sourced or evaluated.
root<-Sys.getenv('MININGCLONEID_REPO_ROOT');stopifnot(nzchar(root))
source(file.path(root,'oxygen/results/structural_identifiability/20261112_culture_v2_seven_point_review_nojson/implementation/helpers.R'))
wp0<-file.path(root,'oxygen/results/structural_identifiability/20261109_seed144_wp0/source/seed144')
p<-readRDS(file.path(wp0,'fit_result.rds'))$best_params;table<-file.path(wp0,'parameter_table_input.csv')
allmap<-o2prof_map(table,p,'zero_extended');p<-o2prof_project_start(allmap,p)$params
nuisance<-c('sigma_growth','sigma_kary','init_mean_2N','init_sd_2N','init_mean_4N','init_sd_4N');map<-o2prof_map(table,p,'zero_extended',p[nuisance])
stopifnot(sum(map$table$status=='active')==14L,sum(map$table$status=='fixed')==6L)
for(name in c('p_misseg','buffer_beta')) {
 q<-p;q[[name]]<-0;z<-o2prof_to_free(map,q);r<-o2prof_from_free(map,z)
 stopifnot(r[[name]]==0,identical(r[nuisance],p[nuisance]))
}
bad<-p[nuisance];bad$sigma_growth<-0
stopifnot(inherits(try(o2prof_map(table,p,'zero_extended',bad),silent=TRUE),'try-error'))
cat('Fourteen-free/six-fixed roles, exact zeros and positive-scale rejection PASS; no model calls\n')
