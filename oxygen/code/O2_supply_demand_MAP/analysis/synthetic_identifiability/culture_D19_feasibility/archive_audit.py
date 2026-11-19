"""HPC archive and actual-read audit only. No model, RNG proposals or optimizer."""
from pathlib import Path
import csv,hashlib,os,re,sys
s=Path(sys.argv[1]).resolve();sys.path.insert(0,str(s/'repair_release'))
from common import read,table,sha,verify
out=s/'verification/C1';out.mkdir(exist_ok=True)
for name in ('source_manifest.tsv','control_manifest.tsv','closure_manifest.tsv','audit_manifest.tsv'):
 for r in read(s/name):
  p=s/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],str(p)
assert (s/'controller/FINISHED').read_text().strip()=='BANKS_READY_FOR_D19_REVIEW'
jobs=read(s/'accounting/jobs.tsv');phases=[];audit=[];bad=[];opened_total=0;project=Path(str(s).split('/oxygen/results/')[0])
for j in jobs:
 if j['task']=='C1':continue
 p=s/'verification'/j['task']
 for r in read(p/'phase_manifest.tsv'):
  q=p/r['path'];assert q.stat().st_size==int(r['bytes']) and sha(q)==r['sha256'],'Phase bytes/hash mismatch '+str(q)
 if j['task']=='P0':assert j['state']=='FAILED' and j['exit_code']=='1:0'
 else:
  assert j['state']=='COMPLETED' and j['exit_code']=='0:0'
  for name in ('gate.tsv','application_result.tsv','policy_gate.tsv'):assert read(p/name)[0]['status']=='PASS'
 phases.append(dict(task=j['task'],state=j['state'],exit_code=j['exit_code'],receipt='PASS',classification='retained_expected_original_failure' if j['task']=='P0' else 'validated'))
 if not j['task'].startswith('F_'):continue
 origin=j['task'][:3];opened=0;findings=[]
 for t in p.glob('trace.*'):
  with t.open(errors='replace') as f:
   for line in f:
    if not (re.match(r'(openat|open)\(',line) or re.match(r'<\.\.\. (openat|open) resumed>',line)):continue
    m=re.search(r'=\s*(\d+)<([^>]+)>\s*$',line)
    if not m:continue
    if 'O_WRONLY' in line or 'O_DIRECTORY' in line:continue
    target=m.group(2);opened+=1;target=target.removesuffix(' (deleted)');q=Path(target)
    if not target.startswith('/'):continue
    inside=q==s or s in q.parents
    if (q==project or project in q.parents) and not inside:findings.append('repository read outside own study: '+target)
    if q.name.lower() in ('fit_data.rds','jobs_2n.rds','jobs_4n.rds','truth.rds','truth_parameters.rds','source_bundle.rds','evaluation_bundle.rds','synthetic_dataset.rds'):findings.append('forbidden endpoint/response/fit object: '+target)
    if inside:
     rel=str(q.relative_to(s))
     if rel.startswith(('inputs/private/','inputs/diagnostic/','reports/','verification/P1/','verification/P1_')):findings.append('forbidden own-study input: '+rel)
     if rel.startswith('verification/F_') and not rel.startswith('verification/'+origin+'_'):findings.append('foreign-origin data: '+rel)
     if rel.startswith('verification/P0_R1/') and rel not in ('verification/P0_R1/feasibility_contract.rds','verification/P0_R1/read_guard.so','verification/P0_R1/rng/'+origin+'.rds'):findings.append('unapproved export: '+rel)
 audit.append(dict(task=j['task'],successful_payload_opens=opened,forbidden_reads=len(findings),status='FAIL' if findings else 'PASS'));opened_total+=opened;bad.extend(findings)
table(out/'phase_receipt_audit.tsv',phases);table(out/'actual_worker_read_audit.tsv',audit)
assert not bad,'Actual scientific read audit failure: '+'; '.join(bad[:10])
old=s.parent/'20261116_culture_growth_karyotype_recovery_01'
assert sha(old/'collection_manifest.tsv')=='16d38d3b5536d627926843f679620fb601d4ff508b1ca7f413601a012d58f170'
for r in read(old/'collection_manifest.tsv'):assert sha(old/r['path'])==r['sha256']
assert read(s/'verification/C0/alias_access_gate.tsv')[0]['status']=='PASS'
modules=sorted(sys.modules);assert not set(modules)&{'json','_json','simplejson','orjson','jpype','jnius'}
table(out/'imported_modules.tsv',[dict(module=n) for n in modules])
table(out/'audit_gate.tsv',[dict(status='PASS',scientific_phases=len(audit),successful_payload_opens=opened_total,forbidden_worker_reads=0,model_calls=0,optimizer_calls=0,old_collection='UNCHANGED',scope='Actual R workers/children; immutable inputs and code; no private response read; stricter closure barrier is prospective template')])
# Freeze a full authoritative HPC inventory. The currently writing C1 allocation
# and mutable root accounting/final reports are outside this intermediate set.
records=[]
for folder in ('inputs','runtime','release','repair_release','control_release','closure_release'):
 for p in sorted((s/folder).rglob('*')):
  if p.is_file():records.append(dict(path=str(p.relative_to(s)),bytes=p.stat().st_size,sha256=sha(p),authority='HPC'))
for p in sorted((s/'verification').rglob('*')):
 if p.is_file() and 'C1' not in p.relative_to(s/'verification').parts:records.append(dict(path=str(p.relative_to(s)),bytes=p.stat().st_size,sha256=sha(p),authority='HPC'))
table(s/'reports/authoritative_hpc_inventory.tsv',records)
table(out/'application_result.tsv',[dict(status='PASS',model_calls=0,optimizer_calls=0)]);table(out/'gate.tsv',[dict(status='PASS')])
table(out/'phase_manifest.tsv',[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='phase_manifest.tsv'])
print('ARCHIVE_AND_READ_AUDIT_PASS',len(records),'HPC files;',opened_total,'scientific payload opens; forbidden0; model_calls0')
