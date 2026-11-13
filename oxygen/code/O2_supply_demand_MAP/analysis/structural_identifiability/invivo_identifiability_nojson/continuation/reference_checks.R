current_record <- function(ctx,p,point,spec) {
 list(parameters=p,components=point$components,predictions=point$predictions,policy=spec$policy,
  free_names=setdiff(names(p),names(spec$fixed)),fixed=spec$fixed,evaluation_mode='fresh_fork',
  contract_sha256=hash_file(file.path(ctx$study,'manifests/contract.rds')),
  source_manifest_sha256=hash_file(file.path(ctx$study,'manifests/input_manifest.tsv')),provenance=ctx$provenance)
}
run_import_replay <- function(study,task,spec) {
 out<-file.path(study,'tests',task);dir.create(out,recursive=TRUE,showWarnings=FALSE)
 ctx<-formal_context(study,spec);start<-resolve_start(ctx,spec);p<-start$parameters
 point<-validate_point(ctx,p,expected=start$expected)
 stopifnot(isTRUE(all.equal(point$predictions$tables,start$record$predictions$tables,tolerance=1e-8)))
 atomic_rds(start,file.path(out,'original_import.rds'));record<-current_record(ctx,p,point,spec)
 proposal<-observation_proposal(ctx,record,NULL);check<-validate_point(ctx,proposal$parameters,expected=proposal$expected)
 stopifnot(abs(check$components$L_b-proposal$expected_burden)<1e-8,
  abs(check$components$L_p-point$components$L_p)<1e-8,abs(check$components$L_n-point$components$L_n)<1e-8)
 # Parameter-history probe remains inside the fresh-child contract.
 q<-p;q[['log10_o2_S0']]<-min(ctx$upper[['log10_o2_S0']],q[['log10_o2_S0']]+.002)
 invisible(score(ctx,q));again<-score(ctx,p);stopifnot(abs(again$L-point$components$L)<1e-8)
 atomic_rds(list(record=record,checked=point$components,policy=spec$policy,provenance=ctx$provenance),file.path(out,'accepted_candidate.rds'))
 atomic_rds(list(proposal=proposal,checked=check),file.path(out,'observation_validation.rds'))
 write_kv(list(pass=TRUE,policy=spec$policy,import_id=spec$import_id,objective=point$components$L,
  difference=point$components$L-start$expected,observation_difference=check$components$L-proposal$expected,
  original_record_retained=TRUE,history_repeat_difference=again$L-point$components$L,root_id=spec$root_id),file.path(out,'result.tsv'))
 write_kv(meter_values(),file.path(out,'runtime_meter.tsv'));finalize_payload(out)
}
