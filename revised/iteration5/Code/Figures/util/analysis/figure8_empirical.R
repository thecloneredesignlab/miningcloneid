# Rebuild the accepted empirical comparison from the author's raw local bundle.
f8_empirical <- function(a) {
  root <- file.path(a$workspace,"data/share_bundle_figures_3_5")
  paths <- file.path(root,c("data/Figure3_PCAWG_included_samples.csv",
    "data/Figure3_DepMap_included_lines.csv","originals/organ_lut.csv"))
  fp_safe(paths);raw <- lapply(paths,data.table::fread)
  stopifnot(nrow(raw[[1]])==2338L,nrow(raw[[2]])==1433L,
    !anyDuplicated(raw[[1]]$sample_id),!anyDuplicated(raw[[2]]$model_id),
    sum(is.finite(raw[[2]]$ploidy))==1185L,nrow(raw[[3]])==19L)
  sites <- data.table::rbindlist(lapply(seq_len(2L),function(i) {
    z <- raw[[i]][is.finite(ploidy)]
    z[,.(n=.N,mean=mean(ploidy),median=median(ploidy),q25=unname(quantile(ploidy,.25,type=7)),
      q75=unname(quantile(ploidy,.75,type=7))),by=.(site=figure3_label)][,cohort:=c("PCAWG","DepMap")[i]]
  }))
  common <- intersect(sites[cohort=="PCAWG"]$site,sites[cohort=="DepMap"]$site)
  sites[,shared_site:=site%in%common];stopifnot(length(common)==18L,nrow(sites)==36L)
  lut <- raw[[3]]
  points <- data.table::data.table(row_id=seq_len(nrow(lut)),organ=lut$organ,cancer_type=lut$type,
    normal_tissue=lut$Normal.tissue.type,wgd_fraction=lut$wgd,oxygen=lut$PtO2.normal.tissue,
    oxygen_sd=lut$SD,oxygen_low=lut$PtO2.normal.tissue-lut$SD,oxygen_high=lut$PtO2.normal.tissue+lut$SD,
    row_policy="original published-source row; no primary filtering or organ deduplication")
  fp_write(sites,file.path(a$run,"tables/empirical_tissue_ploidy_summary.tsv"))
  fp_write(points,file.path(a$run,"tables/author_figure5_points.tsv"))
  test <- stats::cor.test(points$wgd_fraction,points$oxygen)
  fp_write(data.frame(n=nrow(points),pearson_r=unname(test$estimate),p_value=test$p.value,
    fit="unweighted OLS; plot mean confidence band 95%; oxygen bars mean +/- SD"),
    file.path(a$run,"provenance/author_figure5_statistics.tsv"))
  list(sites=sites,points=points,manifest=data.table::data.table(path=paths,sha256=fp_sha(paths),role="author raw empirical input"))
}
