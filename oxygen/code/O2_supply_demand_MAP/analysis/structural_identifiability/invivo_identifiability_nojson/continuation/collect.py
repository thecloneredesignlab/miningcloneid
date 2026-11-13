"""Read and collect closed owned payloads plus a coherent controller snapshot."""
from pathlib import Path
import argparse,subprocess,shlex,re,shutil
import workflow as w

def main():
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();s=a.study.resolve()
 if s.name!='20261112_invivo_continuation_141206Z':raise ValueError('Unexpected study')
 remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+s.name
 lean=(s/'tables/local_mirror_mode.tsv').exists() and w.kv(s/'tables/local_mirror_mode.tsv').get('mode')=='lean_metadata'
 ssh=['ssh','-o','BatchMode=yes','-o','ControlMaster=auto','-o','ControlPath=/private/tmp/ivc_ssh_%C','-o','ControlPersist=600']
 script="""from pathlib import Path
import sys,uuid,shutil
s=Path(REMOTE);sys.path.insert(0,str(s/'implementation'));import workflow as w
state=w.State(s)
with state.locked():
 rows=state.load();g=(s/'state/CURRENT').read_text().strip();collection_id='c_'+uuid.uuid4().hex;directory=s/'collections'/collection_id;directory.mkdir(parents=True)
 (directory/'jobs.tsv').write_bytes((s/'state'/g/'jobs.tsv').read_bytes())
 (directory/'raw_accounting.txt').write_bytes((s/'tables/raw_accounting.txt').read_bytes())
 controls=directory/'controller_snapshot';controls.mkdir()
 for name in ('execution_release.tsv','reviewed_failures.tsv','cache_source_manifest.tsv','metadata_release.tsv','metadata_execution_release.tsv','profile_epoch.tsv','search_batches.tsv','informed_roots.tsv','resolution_policy.tsv','resolution_assignments.tsv','resolution_phase_accounting.tsv','resolution_runtime.tsv','resolution_batches.tsv','resolution_terminal.tsv','resolution_budget_stops.tsv','resolution_baseline_terminal.tsv','resolution_final_phase_accounting.tsv','closure_audit.tsv','final_accounting.tsv','final_role_accounting.tsv','delivery_manifest.tsv','delivery_verified.tsv'):
  if (s/'tables'/name).exists():shutil.copy2(str(s/'tables'/name),str(controls/name))
 final_reports=directory/'final_reports';final_reports.mkdir()
 for name in ('final_accounting.md','resolution_execution_report.md'):
  if (s/'reports'/name).exists():shutil.copy2(str(s/'reports'/name),str(final_reports/name))
 w.table(directory/'origin.tsv',[dict(key='state_generation',value=g),dict(key='captured_utc',value=w.utc())],['key','value'])
 w.table(directory/'manifest.tsv',[dict(path=str(p.relative_to(directory)),sha256=w.sha(p)) for p in sorted(directory.rglob('*')) if p.is_file()],['path','sha256'])
files=['collections/'+collection_id+'/','state/'+g+'/']
for name in ('policies','policy_releases'):
 if (s/name).exists():files.append(name+'/')
inventory=[]
for j in rows:
 files.append('dispatch/'+j['task']+'/')
 if j['mode']=='metadata_audit':
  spec=w.kv(s/'dispatch'/j['task']/'spec.tsv');cid=spec['metadata_collection'];assert __import__('re').fullmatch('c_[0-9a-f]+',cid);files.append('collections/'+cid+'/')
 if not w.closed(j):continue
 scope='fits' if j['mode'] in ('local','global') else 'tests';root=s/scope/j['task']
 if j['state']=='COMPLETED' and j['exit_code']=='0:0':
  w.verify_payload(root)
  inventory.append(dict(task=j['task'],job_id=j['job_id'],mode=j['mode'],scope=scope,payload_sha256=w.sha(root/'payload_manifest.tsv'),full_payload_verified='true'))
  if not LEAN or j['mode']=='metadata_audit':files.append(scope+'/'+j['task']+'/')
  else:
   proof=directory/'remote_payload_manifests'/j['task'];proof.mkdir(parents=True)
   for filename in ('payload_manifest.tsv','COMPLETE','result.tsv','runtime_meter_master.tsv','loaded_namespaces.tsv','start_identity.tsv','reference_recheck.tsv'):
    if (root/filename).exists():shutil.copy2(str(root/filename),str(proof/filename))
 elif (root/'FAILED.tsv').exists():
  files.append(scope+'/'+j['task']+'/FAILED.tsv')
  if j['mode']=='metadata_audit' and (root/'logs').exists():files.append(scope+'/'+j['task']+'/logs/')
 for name in (j['task']+'_submission.txt',j['task']+'_receipt.txt',j['task']+'_'+j['job_id']+'.log'):
  if (s/'logs'/name).exists():files.append('logs/'+name)
 for name in ('image_sha256.txt','launch_exit_code.txt'):
  if (s/'runtime'/j['task']/name).exists():files.append('runtime/'+j['task']+'/'+name)
w.table(directory/'payload_inventory.tsv',inventory,['task','job_id','mode','scope','payload_sha256','full_payload_verified'])
w.table(directory/'manifest.tsv',[dict(path=str(p.relative_to(directory)),sha256=w.sha(p)) for p in sorted(directory.rglob('*')) if p.is_file() and p!=directory/'manifest.tsv'],['path','sha256'])
print('SNAP\\t'+collection_id+'\\t'+g);print('\\n'.join(files))
""".replace('REMOTE',repr(remote)).replace('LEAN',repr(lean))
 result=subprocess.run(ssh+['4482173@red.moffitt.org','python3','-'],input=script,universal_newlines=True,stdout=subprocess.PIPE,check=True)
 lines=result.stdout.splitlines();marker=lines.pop(0).split('\t');assert marker[0]=='SNAP' and re.fullmatch('c_[0-9a-f]+',marker[1]) and re.fullmatch('g_[0-9a-f]+',marker[2])
 assert all(not x.startswith('/') and '..' not in Path(x).parts for x in lines)
 listing=s/'tables/collection_files.txt';listing.write_text('\n'.join(lines)+'\n')
 copy_ssh=['ssh','-o','BatchMode=yes','-o','ControlMaster=no','-o','ControlPath=none','-o','ServerAliveInterval=30','-o','ServerAliveCountMax=6','-o','ConnectTimeout=30']
 subprocess.run(['rsync','-a','--recursive','--timeout=180','--files-from='+str(listing),'-e',' '.join(shlex.quote(x) for x in copy_ssh),'4482173@red.moffitt.org:'+remote+'/',str(s)+'/'],check=True)
 collection=s/'collections'/marker[1];generation=s/'state'/marker[2]
 w.verify_manifest(collection,collection/'manifest.tsv');w.verify_manifest(generation,generation/'manifest.tsv')
 assert w.read_table(generation/'jobs.tsv')==w.read_table(collection/'jobs.tsv')
 assert w.kv(collection/'origin.tsv')['state_generation']==marker[2]
 w.atomic_text(s/'tables/jobs.tsv',(collection/'jobs.tsv').read_text());w.atomic_text(s/'tables/raw_accounting.txt',(collection/'raw_accounting.txt').read_text())
 w.atomic_text(s/'state/CURRENT',marker[2]+'\n')
 rows=w.State(s).load();actual=sum(int(j['actual_core_seconds']) for j in rows if w.closed(j) and j['actual_core_seconds']);live=sum(int(j['reserved_core_seconds']) for j in rows if not w.closed(j))
 w.atomic_text(s/'STATUS.md','# Improvement-driven in vivo continuation\n\nCoherent collected generation '+marker[2]+'; collected '+w.utc()+'.\n\nNew closed actual '+str(actual)+' core-seconds; live reservations '+str(live)+'. Historical lifetime301846/formal271355remain charged. Resolution01 baseline220/profile360core-hour caps and finite non-resetting role allowances. Culture protected; pinned SIF. Profiles require both baselines.\n\n'+'\n'.join(j['task']+' '+j['job_id']+' '+j['state'] for j in rows if not w.closed(j))+'\n')
 inventory={r['task']:r for r in w.read_table(collection/'payload_inventory.tsv')}
 for j in w.State(s).load():
  if j['state']!='COMPLETED' or j['exit_code']!='0:0':continue
  if not lean or j['mode']=='metadata_audit':w.verify_payload(s/('fits' if j['mode'] in ('local','global') else 'tests')/j['task'])
  else:
   proof=collection/'remote_payload_manifests'/j['task'];assert inventory[j['task']]['full_payload_verified']=='true' and w.sha(proof/'payload_manifest.tsv')==inventory[j['task']]['payload_sha256']
   assert (proof/'COMPLETE').read_text().strip()==inventory[j['task']]['payload_sha256']
   local_proof=s/inventory[j['task']]['scope']/j['task'];local_proof.mkdir(parents=True,exist_ok=True)
   shutil.copy2(proof/'payload_manifest.tsv',local_proof/'payload_manifest.tsv')
   w.atomic_text(local_proof/'REMOTE_PAYLOAD_ONLY.txt','Lean mirror: full native/checkpoint data are authoritative on HPC; collection inventory verifies the complete remote payload.\n')
 w.table(s/'tables/remote_payload_inventory.tsv',list(inventory.values()),['task','job_id','mode','scope','payload_sha256','full_payload_verified'])
 packets=[j for j in rows if j['mode']=='metadata_audit' and j['state']=='COMPLETED' and j['exit_code']=='0:0']
 if lean and packets:
  packet=s/'tests'/packets[-1]['task']/'metadata_snapshot'
  excluded={'jobs.tsv','raw_accounting.txt','execution_release.tsv','reviewed_failures.tsv','cache_source_manifest.tsv','metadata_release.tsv','metadata_execution_release.tsv','collection_files.txt','profile_epoch.tsv','search_batches.tsv','informed_roots.tsv','local_mirror_mode.tsv','local_replica_eviction.tsv','resolution_policy.tsv','resolution_assignments.tsv','resolution_phase_accounting.tsv','resolution_runtime.tsv','resolution_batches.tsv','resolution_terminal.tsv','resolution_budget_stops.tsv'}
  for scope in ('tables','reports','plots'):
   for source in (packet/scope).rglob('*'):
    if not source.is_file() or (scope=='tables' and source.name in excluded):continue
    target=s/scope/source.relative_to(packet/scope);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
  w.typed(s/'tables/metadata_packet_source.tsv',dict(task=packets[-1]['task'],payload_sha256=inventory[packets[-1]['task']]['payload_sha256'],authority='SIF-pinned HPC archived-output validation; local lean mirror',model_evaluations=0))
 for source in (collection/'controller_snapshot').glob('*.tsv'):
  w.atomic_text(s/'tables'/source.name,source.read_text())
 for source in (collection/'final_reports').glob('*.md'):
  w.atomic_text(s/'reports'/source.name,source.read_text())
 print(('Lean metadata with authoritative HPC full-payload proofs' if lean else 'Closed full payloads')+' and immutable accounting collection '+marker[1]+' verified; no active fitting payload copied.')
if __name__=='__main__':main()
