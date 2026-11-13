# Exercise the actual finite strategy body over synthetic metadata, zero model calls.
args<-commandArgs(trailingOnly=TRUE)
entry<-normalizePath(if(length(args)>1L)args[2] else args[1])
source(file.path(entry,'adapter.R'));source(file.path(entry,'resolution_rules.R'))
fixture<-function(after,eligible,initial_pass=FALSE){
 expected_eligible<-eligible
 study<-tempfile('completed_batch_fixture_');dir.create(study);dir.create(file.path(study,'tables'))
 on.exit(unlink(study,recursive=TRUE),add=TRUE)
 jobs<-data.frame(task=c('existing','global_h','global_u','polish_h','polish_u'),
  policy='production',state='COMPLETED',exit_code='0:0')
 slots<-data.frame(task=jobs$task,phase='baseline',
  slot_class=c('existing_local','baseline_global','baseline_global','baseline_polish','baseline_polish'),
  slot_key=c('existing_root','production:1:hybrid','production:1:uniform','production:1:hybrid','production:1:uniform'))
 decisions<-data.frame(task=jobs$task,decision='resolution_stop',complete=TRUE)
 policies<-data.frame(policy=c('production','data_only'),pass=c(initial_pass,TRUE),
  known_feasible_pass=TRUE,known_feasible_bound=c(13,12))
 chains<-data.frame(policy='production',family='uniform',root_id='uniform_fixture',verified_best=after)
 write_tsv(jobs,file.path(study,'tables/jobs.tsv'));write_tsv(slots,file.path(study,'tables/resolution_assignments.tsv'))
 write_tsv(decisions,file.path(study,'tables/resolution_local_decisions.tsv'))
 write_tsv(policies,file.path(study,'tables/baseline_policy_gates.tsv'))
 write_tsv(chains,file.path(study,'tables/baseline_chain_summary.tsv'))
 write_tsv(data.frame(policy='production',batch=1L,before_uniform=14,before_feasible=13),file.path(study,'tables/resolution_batches.tsv'))
 atomic_rds(list(fixture=TRUE,checked=list(L=13)),file.path(study,'tables/baseline_production.rds'))
 atomic_rds(list(fixture=TRUE,checked=list(L=after)),file.path(study,'tables/root_uniform_fixture_selected.rds'))
 collection<-file.path(study,'collections/c_a');dir.create(collection,recursive=TRUE)
 write_tsv(jobs,file.path(collection,'jobs.tsv'))
 write_tsv(data.frame(path='jobs.tsv',sha256=hash_file(file.path(collection,'jobs.tsv'))),file.path(collection,'manifest.tsv'))
 # Omit only the CLI/bootstrap lines; execute the real dispatch decision body.
 text<-readLines(file.path(entry,'strategy_review.R'))
 first<-which(startsWith(text,"jobs<-read_tsv("))[1L];stopifnot(!is.na(first))
 body<-parse(text=text[first:length(text)]);args<-c(study,'c_a')
 eval(body,envir=environment())
 review<-read_tsv(file.path(study,'tables/strategy_review.tsv'))
 stopifnot(review$decision[review$policy=='production']==if(initial_pass)'accepted' else if(expected_eligible)'next_global_batch' else 'finite_baseline_repair')
 fixture_outcome_path<-file.path(study,'tables/resolution_batch_outcomes.tsv')
 fixture_original_bytes<-readBin(fixture_outcome_path,'raw',file.info(fixture_outcome_path)$size)
 fixture_frozen<-read_tsv(fixture_outcome_path);stopifnot(identical(fixture_frozen$eligible,expected_eligible))
 if(!expected_eligible){
  # Two later repairs improve Uniform by far more than0.002. They cannot
  # retroactively qualify the already completed, non-qualifying first batch.
  jobs<-rbind(jobs,data.frame(task='repair1',policy='production',state='COMPLETED',exit_code='0:0'))
  slots<-rbind(slots,data.frame(task='repair1',phase='baseline',slot_class='baseline_repair',slot_key='production'))
  decisions<-rbind(decisions,data.frame(task='repair1',decision='resolution_stop',complete=TRUE))
  chains$verified_best<-13.99
  policies$pass<-c(FALSE,TRUE)
  write_tsv(jobs,file.path(study,'tables/jobs.tsv'));write_tsv(slots,file.path(study,'tables/resolution_assignments.tsv'))
  write_tsv(chains,file.path(study,'tables/baseline_chain_summary.tsv'))
  write_tsv(policies,file.path(study,'tables/baseline_policy_gates.tsv'))
  write_tsv(decisions,file.path(study,'tables/resolution_local_decisions.tsv'))
  eval(body,envir=environment())
  review<-read_tsv(file.path(study,'tables/strategy_review.tsv'))
  stopifnot(review$decision[review$policy=='production']=='resolution_unresolved')
  decisions$decision[decisions$task=='repair1']<-'continue'
  write_tsv(decisions,file.path(study,'tables/resolution_local_decisions.tsv'))
  eval(body,envir=environment())
  review<-read_tsv(file.path(study,'tables/strategy_review.tsv'))
  stopifnot(review$decision[review$policy=='production']=='finite_baseline_repair')
  jobs<-rbind(jobs,data.frame(task='repair2',policy='production',state='COMPLETED',exit_code='0:0'))
  slots<-rbind(slots,data.frame(task='repair2',phase='baseline',slot_class='baseline_repair',slot_key='production'))
  write_tsv(jobs,file.path(study,'tables/jobs.tsv'));write_tsv(slots,file.path(study,'tables/resolution_assignments.tsv'))
  eval(body,envir=environment())
  review<-read_tsv(file.path(study,'tables/strategy_review.tsv'))
  stopifnot(review$decision[review$policy=='production']=='resolution_unresolved',
   identical(fixture_original_bytes,readBin(fixture_outcome_path,'raw',file.info(fixture_outcome_path)$size)))
 }
 cat('Actual strategy completed-outcome fixture eligible=',expected_eligible,': PASS\n')
}
fixture(13.9981,FALSE);fixture(13.9979,TRUE);fixture(13.9981,FALSE,initial_pass=TRUE)
cat('Zero model calls.\n');assert_namespaces()
