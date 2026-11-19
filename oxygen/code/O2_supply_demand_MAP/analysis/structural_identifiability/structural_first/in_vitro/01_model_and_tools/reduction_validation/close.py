"""Terminal accounting and read-only evidence closure; no scientific compute."""
from pathlib import Path
import argparse,fcntl,importlib.util
def module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def close(s,code):
 w=module('workflow',code/'workflow.py');collector=module('collect',code/'collect.py')
 if (s/'COMPLETE').exists():raise ValueError('Study is already closed')
 with (s/'state/JOB_LOCK').open('a') as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);jobs=w.reconcile(s)
  assert len(jobs)==76 and all(r['state']=='COMPLETED' and r['exit_code']=='0:0' for r in jobs),'Terminal allocation gate failed'
  assert all(w.payload_ok(s/'results'/r['task']) for r in jobs),'Unsealed or stale evidence'
  for release in sorted((s/'releases').iterdir()):w.verify(s,release)
  collector.collect(s,code)
  budget=w.position(jobs);actual=sum(r['actual'] for r in budget);assert sum(r['committed'] for r in budget)==0
  w.table(s/'reports/final_resource_accounting.tsv',budget,['phase','actual','committed','cap','remaining'])
  w.typed(s/'reports/final_status.tsv',dict(closed_utc=w.utc(),execution_complete=True,scientific_status='HOLD: no aggregate candidate passed every state/flux fixture',structural_classification_complete=False,step05_release=False,actual_allocated_core_seconds=actual,actual_allocated_core_hours=actual/3600,independent_core_hour_cap=96,live_reserved_core_hours=0,owned_parent_allocations=len(jobs),owned_completed_zero_exit=len(jobs),failed_allocations=0,retries=0,old_D0_scope='No controller or task modifications; read-only mounts used',java_used=False,javascript_used=False,json_used=False))
  reference={r['case_id']:r for r in w.read(s/'inputs/cases.tsv')};tails=w.read(s/'reports/all_tail_acceptance.tsv');comparisons=w.read(s/'reports/all_production_domain_contrasts.tsv');paired=w.read(s/'reports/aggregation_contrast_summary.tsv')
  lines=['# Interpretation, assumptions and registered holds','',
  'The numerical audit separated four questions; it did not classify structural identifiability.','',
  '## Initial lifting and dynamic closure','For the twenty-subclass candidate, the continuous full-original versus full-lifted contrast had maximum class TV 0.0288532, log-live gap 0.0530572 and mean chromosome gap 0.590175, but some per-cell flux errors still exceeded tolerance. Even after matching the lifted initial distribution, full-lifted versus coarse propagation reached class TV 0.209478 and mean chromosome gap 1.69450. This is evidence of dynamic closure loss under this fixed lifting, in addition to initial-information loss. These contrast maxima can occur at different times and fixtures and must not be added or used as causal percentages.','',
  'All four candidates failed the combined expected karyotype/flow gate for every paired scheme/fixture. Some narrow initial distributions are immediately changed by uniform lifting. The full expected readout should therefore be checked independently of the five-class state summary; there are no sampling-power or experimental-death claims here.','',
  '## Production time step','Production DT=0.05 and continuous-time propagation have a worst absolute log-live difference 0.288734 over these fourteen-day unpassaged fixtures. All 24 dt sequences converged within ten halvings under the registered adjacent-plus-reference test. This supports the validity of the continuous comparison reference but does not silently replace the production discrete model or prove equivalence over every experimental passage history.','',
  '## Finite boundaries versus extended biology','Upper-only and both-ended extension are distinct assumption cases. The lower extension permits chromosome states 1:21 to live and use the unchanged production rate functions; those states are excluded from the original production domain. This can create strong selection and large abundance differences and must not be described as correcting a known numerical bug.','']
  maximum=max(comparisons,key=lambda r:float(r['max_full_tv']));maximum_log=max(comparisons,key=lambda r:float(r['max_log_live']))
  r=reference[maximum['case_id']];lines.append(f'The largest finite-versus-extended full-distribution TV was {float(maximum["max_full_tv"]):.9g}, in {maximum["variant"]}, {r["parameter_reference"]}, {r["cohort"]}, {r["oxygen_pct"]} percent oxygen. The largest absolute log-live gap was {float(maximum_log["max_log_live"]):.9g}. These are extrapolation contrasts, not claims that the extended support is biologically correct.')
  lines+=['','The strict two-adjacent-domain gate remained unresolved for:','', '| Reference | Background | Oxygen percent | Upper-only | Both-ended |','|---|---|---:|---|---|']
  for r in tails:
   if r['upper_only_converged']!='true' or r['both_ends_converged']!='true':
    evidence=w.read(s/'results'/('tail_'+r['case_id'])/'boundary_flux_ledger.tsv')
    for branch,key in [('upper_only','upper_only_converged'),('both_ends','both_ends_converged')]:
     if r[key]!='true':
      latest=next(x for x in evidence if x['variant']==branch+'_1232');previous=next(x for x in evidence if x['variant']==branch+'_616');assert latest['upper_flux_gate']=='TRUE' and previous['upper_flux_gate']=='FALSE','Review unexpected domain hold before interpreting'
    lines.append(f'| {r["reference"]} | {r["cohort"]} | {r["oxygen_pct"]} | {r["upper_only_converged"]} | {r["both_ends_converged"]} |')
  lines+=['','In these held branches, the largest single domain (1232 chromosomes) met its own upper-flux/edge-occupancy criterion and its adjacent distribution/abundance comparison was small. The previous domain (616) did not meet the strict upper-flux criterion; two individually acceptable adjacent domains were required. No 2464-chromosome extension was authorized or run. This is a registered conservative hold, distinct from solver failure.','',
  '## Consequence for the structural-first plan','Do not promote the five-, ten- or twenty-class autonomous model to a surrogate for full production behavior. No paired full-history branch or symbolic batch was released. The next model decision should explicitly distinguish finite-support production semantics, positive-state extrapolation, and the response-function/within-class representation needed by a structural tool. This audit does not justify labeling any of the fourteen mechanisms identifiable or non-identifiable. Step 5 remains pending a separately reviewed structural result.','']
  (s/'reports/interpretation_and_holds.md').write_text('\n'.join(lines))
  text=(s/'REPORT.md').read_text();start=text.index('- Accounting snapshot:');end=text.index('\n',start)
  text=text[:start]+f'- Final terminal accounting: **{actual/3600:.9f} allocated core-hours** ({actual} allocated CPU-seconds), no live reservation, independent cap 96. All {len(jobs)} owned parent allocations COMPLETED with exit 0:0; all internal results and payload hashes validated. No reserve retry used.'+text[end:]
  text+='\nSee reports/interpretation_and_holds.md for initial-lifting versus closure contrasts, domain-assumption effects and the exact reason for the conservative adjacent-domain holds.\n';(s/'REPORT.md').write_text(text)
  forbidden=[];language=[]
  selected=[]
  for p in sorted(s.rglob('*')):
   if not p.is_file():continue
   if p.suffix.lower() in ('.json','.jsonl','.js','.java','.jar','.class'):forbidden.append(str(p))
   if p.relative_to(s).parts[0] in ('runtime','cache','state'):continue
   if p.name in ('collection_manifest.tsv','COMPLETE'):continue
   if p.suffix.lower() in ('.md','.tsv','.txt','.py','.r','.sh','.out','.err'):
    content=p.read_text(errors='replace')
    if any('\u4e00'<=c<='\u9fff' for c in content):language.append(str(p))
   selected.append(p)
  assert not forbidden and not language,(forbidden,language)
  w.typed(s/'reports/source_and_format_integrity.tsv',dict(source_files_verified=len(w.read(s/'inputs/external_sources.tsv')),releases_verified=len(list((s/'releases').iterdir())),input_manifest_verified=True,forbidden_artifact_files=0,Chinese_text_in_new_scientific_artifacts=0,collection_excludes='Private native caches, duplicate reference cache and transient locks; raw scientific evidence is included in results',source_guards='SIF SHA256 and frozen sources checked at allocation start and end'))
  selected.append(s/'reports/source_and_format_integrity.tsv')
  selected=sorted(set(selected));manifest=[dict(path=str(p.relative_to(s)),sha256=w.sha(p),bytes=p.stat().st_size) for p in selected];w.table(s/'collection_manifest.tsv',manifest,['path','sha256','bytes']);(s/'COMPLETE').write_text(w.sha(s/'collection_manifest.tsv'))
  print('CLOSED',len(jobs),'parents',actual,'allocated CPU-seconds',len(manifest),'files',sum(r['bytes'] for r in manifest),'bytes',flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();close(a.study.resolve(),Path(__file__).resolve().parent)
