"""Post-selection archive/read audit. No response calculation or optimization."""
from pathlib import Path
import os,re,sys
s=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();sys.path.insert(0,str(s/'release'))
from common import *
assert (s/'controller/SELECTION_LOCKED').exists(),'Selection must be locked first'
verify(s)
for folder in ['control_release_v2','control_release_v3','assessment_release']:
 for r in read(s/folder/'manifest.tsv'):
  p=s/folder/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],'Control/postselection source conflict '+str(p)
jobs=read(s/'accounting/jobs.tsv');project=Path(str(s).split('/oxygen/results/')[0]);audits=[];phases=[];bad=[]
for j in jobs:
 t=j['task']
 if t.startswith('C'):continue
 p=s/'verification'/t
 assert j['state']=='COMPLETED' and j['exit_code']=='0:0','Retain/review failed records separately: '+t
 for r in read(p/'phase_manifest.tsv'):
  q=p/r['path'];assert q.stat().st_size==int(r['bytes']) and sha(q)==r['sha256'],'Phase receipt conflict '+str(q)
 for n in ['application_result.tsv','gate.tsv','policy_gate.tsv']:assert read(p/n)[0]['status']=='PASS',t+'/'+n
 phases.append(dict(task=t,slurm_state=j['state'],exit_code=j['exit_code'],receipt='PASS',application='PASS',policy='PASS'))
 if not t.startswith(('D_','F_')):continue
 origin=t[:3];count=0;findings=[]
 for trace in p.glob('trace.*'):
  with trace.open(errors='replace') as f:
   for line in f:
    if not (re.match(r'(openat|open)\(',line) or re.match(r'<\.\.\. (openat|open) resumed>',line)):continue
    m=re.search(r'=\s*(\d+)<([^>]+)>\s*$',line)
    if not m or 'O_WRONLY' in line or 'O_DIRECTORY' in line:continue
    target=m[2].removesuffix(' (deleted)');q=Path(target)
    if not target.startswith('/'):continue
    count+=1;inside=q==s or s in q.parents
    if project in q.parents and not inside:findings.append('Repository read outside new study: '+target)
    if q.name.lower() in {'fit_data.rds','jobs_2n.rds','jobs_4n.rds','truth.rds','truth_parameters.rds','source_bundle.rds','evaluation_bundle.rds','synthetic_dataset.rds'}:findings.append('Forbidden scientific payload: '+target)
    if inside:
     rel=str(q.relative_to(s));allowed=rel.startswith(('release/','runtime/','inputs/production/',f'verification/{origin}_')) or rel in {'inputs/public/contract.rds','inputs/public/parameter_roles.tsv','verification/T00/read_guard.so',f'verification/T00/identities/{origin}.rds',f'verification/T00/rng/{origin}.rds'}
     if t.startswith('D_'):allowed=allowed or rel in {'inputs/public/D19.rds',f'inputs/banks/{origin}.rds',f'controller/F_{origin[2]}_bank_path.txt'} or rel.startswith(f'verification/F_{origin[2]}_')
     else:allowed=allowed or rel==f'controller/{origin}_bank_path.txt'
     if not allowed:findings.append('Unapproved role input: '+rel)
 audits.append(dict(task=t,successful_payload_opens=count,forbidden_reads=len(findings),status='FAIL' if findings else'PASS'));bad+=findings
assert not bad,'Actual worker read conflict: '+'; '.join(bad[:10])
table(out/'actual_worker_read_audit.tsv',audits);table(out/'phase_receipt_audit.tsv',phases)
# Historical manifest identities are reused; original scientific payloads are untouched.
old=s.parent/'20261117_culture_D19_feasibility_01';assert sha(old/'local_collection_manifest.tsv')=='2199c8ed94bc93289fe0cfd3fe8d657f49ddfaf7ce2cfedeb6b27b7ec6fefb9a'
older=s.parent/'20261116_culture_growth_karyotype_recovery_01';assert sha(older/'collection_manifest.tsv')=='16d38d3b5536d627926843f679620fb601d4ff508b1ca7f413601a012d58f170'
for r in read(s/'accounting/HPC_preservation.tsv'):
 p=project/r['path'];assert (sha(p) if p.is_file() else 'NOT_PRESENT')==r['HPC_sha256'],'Pre-existing HPC content changed '+r['path']
records=[]
for folder in ['inputs','runtime','release','control_release_v2','control_release_v3','assessment_release']:
 for p in sorted((s/folder).rglob('*')):
  if p.is_file() and not p.is_symlink():records.append(dict(path=str(p.relative_to(s)),bytes=p.stat().st_size,sha256=sha(p),authority='HPC'))
for p in sorted((s/'verification').rglob('*')):
 if p.is_file() and out not in p.parents:records.append(dict(path=str(p.relative_to(s)),bytes=p.stat().st_size,sha256=sha(p),authority='HPC'))
table(s/'reports/authoritative_hpc_inventory.tsv',records)
mods=sorted(sys.modules);assert not set(mods)&{'json','_json','simplejson','orjson','jpype','jnius'};table(out/'imported_modules.tsv',[dict(module=m) for m in mods])
table(out/'audit_gate.tsv',[dict(status='PASS',scientific_phases=len(audits),successful_payload_opens=sum(r['successful_payload_opens'] for r in audits),forbidden_worker_reads=0,full_forward_calls=0,old_source_manifests='UNCHANGED',preexisting_dirty_content='UNCHANGED',HPC_inventory_files=len(records))])
table(out/'application_result.tsv',[dict(status='PASS',full_forward_calls=0,optimizer_calls=0)]);table(out/'gate.tsv',[dict(status='PASS')])
table(out/'policy_gate.tsv',[dict(status='PASS',coverage='postselection standard-library metadata inventory; scientific worker/child traces independently audited',excluded_imports=0)])
print('ARCHIVE_AND_ACTUAL_READ_AUDIT_PASS',len(records))
