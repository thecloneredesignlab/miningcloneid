# Analytic metadata fixtures only: no real model, fits, scheduler or allocation.
args<-commandArgs(trailingOnly=TRUE);source(file.path(args[1],'metadata_stage.R'))
stopifnot(!metadata_profile_phase(FALSE,TRUE,FALSE),metadata_profile_phase(TRUE,TRUE,FALSE),!metadata_profile_phase(TRUE,FALSE,FALSE),!metadata_profile_phase(TRUE,TRUE,TRUE))
s<-tempfile('ivc_metadata_stage_fixture_');dir.create(s);dir.create(file.path(s,'tables'));dir.create(file.path(s,'fits'))
writeLines('original audit evidence',file.path(s,'tables/value.tsv'));writeLines('immutable closed fingerprint',file.path(s,'fits/fingerprint.tsv'))
frozen<-freeze_baseline_evidence(s,c('tables/value.tsv','fits/fingerprint.tsv'))
writeLines('later audit evidence',file.path(s,'tables/value.tsv'))
stopifnot(readLines(file.path(s,frozen[1]))=='original audit evidence',frozen[2]=='fits/fingerprint.tsv')
unlink(s,recursive=TRUE);cat('Metadata stage and immutable accepted-baseline evidence fixtures passed; zero model calls.\n')
