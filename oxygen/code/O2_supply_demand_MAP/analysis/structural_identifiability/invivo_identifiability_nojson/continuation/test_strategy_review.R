# Labeled analytic controller metadata only; no model, scheduler or submission.
args<-commandArgs(trailingOnly=TRUE);entry<-normalizePath(args[1]);source(file.path(entry,'adapter.R'))
s<-tempfile('ivc_strategy_metadata_fixture_');dir.create(s);dir.create(file.path(s,'tables'))
job<-function(task,policy,mode='local',state='COMPLETED')data.frame(task=task,kind='baseline',policy=policy,mode=mode,state=state,exit_code='0:0')
jobs<-do.call(rbind,lapply(c('production','data_only'),function(policy)rbind(job(paste0('b_',policy,'_A_w10000'),policy),job(paste0('b_',policy,'_U_w10000'),policy))))
chains<-data.frame(task=jobs$task,policy=jobs$policy,family=rep(c('archive','uniform'),2),verified_best=rep(c(14,14.2),2),termination='call_limit')
policies<-data.frame(policy=c('production','data_only'),pass=FALSE,best_objective=14,known_feasible_bound=14)
run<-function(){write_tsv(jobs,file.path(s,'tables/jobs.tsv'));write_tsv(chains,file.path(s,'tables/baseline_chain_summary.tsv'));write_tsv(policies,file.path(s,'tables/baseline_policy_gates.tsv'));stopifnot(system2(file.path(R.home('bin'),'Rscript'),c('--vanilla',shQuote(file.path(entry,'strategy_review.R')),shQuote(s)),stdout=FALSE)==0L);read_tsv(file.path(s,'tables/strategy_review.tsv'))}
stopifnot(all(run()$decision=='continue_local'))
chains$termination<-'poll_stationary';stopifnot(all(run()$decision=='initial_global_batch'))
# A checked stationary root does not throttle another already-ready branch.
jobs<-rbind(jobs,job('b_production_other_w10001','production',state='RUNNING'));stopifnot(run()$decision[1]=='initial_global_batch');jobs<-head(jobs,-1)
batches<-data.frame(policy=c('production','data_only'),batch=1,own_root='fixture_uniform',before_uniform=14.2,before_feasible=14)
write_tsv(batches,file.path(s,'tables/search_batches.tsv'))
for(policy in c('production','data_only'))for(mode in c('hybrid','uniform')){
 task<-paste0('g_',policy,'_c0001_',mode,'_w10000');jobs<-rbind(jobs,job(task,policy,'global'));dir.create(file.path(s,'fits',task),recursive=TRUE)
 write_kv(list(window_gain=0,objective=14.2),file.path(s,'fits',task,'result.tsv'))
 polish<-paste0('b_',task,'_polish_w0001');jobs<-rbind(jobs,job(polish,policy));chains<-rbind(chains,data.frame(task=polish,policy=policy,family='uniform',verified_best=14.2,termination='poll_stationary'))
}
stopifnot(all(run()$decision=='strategy_stagnation'))
batches$before_uniform<-14.3;write_tsv(batches,file.path(s,'tables/search_batches.tsv'));z<-run();stopifnot(all(z$decision=='next_global_batch'),all(z$batch==2L))
write_kv(list(window_gain=1e-7,objective=14.2),file.path(s,'fits/g_production_c0001_hybrid_w10000/result.tsv'));stopifnot(run()$decision[1]=='continue_global')
unlink(s,recursive=TRUE);cat('Analytic strategy fixtures passed: incomplete continuation, concurrent readiness, true stagnation, protected gain and unbounded window IDs; zero model calls.\n')
