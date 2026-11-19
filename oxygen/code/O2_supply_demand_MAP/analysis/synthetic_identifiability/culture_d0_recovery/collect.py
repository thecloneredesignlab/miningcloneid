"""Read-only collection after terminal accounting; full HPC proof, lean fit mirror."""
from pathlib import Path
import argparse,csv,hashlib,shlex,subprocess,sys
sys.dont_write_bytecode=True
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def read(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def kv(p):return {r['key']:r['value'] for r in read(p)}
def verify(root,rows):
    for r in rows:
        if sha(root/r['path'])!=r['sha256']:raise ValueError('Collection mismatch: '+r['path'])
def main():
    ap=argparse.ArgumentParser();ap.add_argument('study',type=Path)
    ap.add_argument('--ssh-control-path',type=Path)
    args=ap.parse_args();s=args.study.resolve()
    if s.name!='20261113_culture_d0_recovery_01':raise ValueError('Unexpected owned study')
    remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/synthetic_identifiability/'+s.name
    script=r'''from pathlib import Path
import csv,hashlib,uuid,shutil
s=Path(REMOTE)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def kv(p):return {r['key']:r['value'] for r in read(p)}
def table(p,rows,fields):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
proof=kv(s/'tables/delivery_verified.tsv');manifest=s/'tables/delivery_manifest.tsv'
assert proof['pass_gate']=='true' and proof['delivery_manifest_sha256']==sha(manifest)
assert kv(s/'tables/final_accounting.tsv')['active_parents']=='0'
listing=read(manifest)
for row in listing:assert sha(s/row['path'])==row['sha256'],row['path']
# The numerical finalizer omits delivery-proof basenames recursively to avoid
# self-reference. Bind the preserved earlier closure's two proof files here.
history=s/'closure_history/stop_01'
oldproof=kv(history/'tables/delivery_verified.tsv')
oldmanifest=history/'tables/delivery_manifest.tsv'
assert oldproof['pass_gate']=='true' and oldproof['delivery_manifest_sha256']==sha(oldmanifest)
assert sha(oldmanifest)=='2b954f43a463e666a86eaa1ef97467256e79ab284679bad290b4a970fdcc2edd'
for row in read(oldmanifest):assert sha(history/row['path'])==row['sha256'],row['path']
cid='c_'+uuid.uuid4().hex;d=s/'collection'/cid;d.mkdir(parents=True)
for name in ['delivery_manifest.tsv','delivery_verified.tsv','final_accounting.tsv']:
 shutil.copy2(str(s/'tables'/name),str(d/name))
files=[];fits=set()
for row in listing:
 p=Path(row['path'])
 if p.parts[0]=='fits':
  fits.add(p.parts[1])
  if p.name=='COMPLETE' or 'improvements' in p.parts or 'initial_inputs' in p.parts:continue
  if p.suffix=='.rds' and p.name not in ('accepted_candidate.rds','requested_initial_population.rds','actual_initial_population.rds','rng_initial.rds','rng_final.rds','native_deoptim.rds','native_sbplx.rds','native_pattern.rds'):continue
  if p.suffix=='.rds' and (s/p).stat().st_size>16*1024*1024:continue
 files.append(row['path'])
files+=['tables/delivery_manifest.tsv','tables/delivery_verified.tsv']
files+=['closure_history/stop_01/tables/delivery_manifest.tsv','closure_history/stop_01/tables/delivery_verified.tsv']
table(d/'fit_mirror_scope.tsv',[dict(task=t,complete_local_payload=False,authority='HPC full payload checked against final delivery manifest') for t in sorted(fits)],['task','complete_local_payload','authority'])
table(d/'captured_files.tsv',[dict(path=p,sha256=sha(s/p)) for p in sorted(set(files))],['path','sha256'])
table(d/'manifest.tsv',[dict(path=p.name,sha256=sha(p)) for p in sorted(d.iterdir()) if p.is_file()],['path','sha256'])
assert proof['delivery_manifest_sha256']==sha(manifest)
print(cid);print('\n'.join(sorted(set(files))+['collection/'+cid+'/']))
'''.replace('REMOTE',repr(remote))
    connection=['ssh','-o','LogLevel=ERROR','-o','BatchMode=yes','-o','ConnectTimeout=15','-o','ServerAliveInterval=15','-o','ServerAliveCountMax=2']
    if args.ssh_control_path:
        connection+=['-o','ControlPath='+str(args.ssh_control_path),'-o','ControlMaster=auto']
    ssh=connection+['4482173@red.moffitt.org','python3','-B','-']
    result=subprocess.run(ssh,input=script,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,check=True)
    lines=result.stdout.splitlines();cid=lines[0];assert cid.startswith('c_') and all(not x.startswith('/') and '..' not in Path(x).parts for x in lines[1:])
    (s/'collection').mkdir(exist_ok=True);f=s/'collection/files.txt';f.write_text('\n'.join(lines[1:])+'\n')
    subprocess.run(['rsync','-a','--recursive','--timeout=180','--files-from='+str(f),'-e',shlex.join(connection),'4482173@red.moffitt.org:'+remote+'/',str(s)+'/'],check=True)
    d=s/'collection'/cid;verify(d,read(d/'manifest.tsv'));verify(s,read(d/'captured_files.tsv'))
    for row in read(d/'fit_mirror_scope.tsv'):
        root=s/'fits'/row['task'];root.mkdir(parents=True,exist_ok=True)
        if (root/'COMPLETE').exists():raise ValueError('Partial fit mirror cannot carry COMPLETE')
        (root/'REMOTE_PAYLOAD_ONLY.txt').write_text('Full native/improvement payload remains authoritative on HPC; this is a verified lean mirror.\n')
    with (s/'collection/local_verification.tsv').open('w',newline='') as out:
        w=csv.DictWriter(out,['key','value'],delimiter='\t');w.writeheader()
        for k,v in dict(pass_gate=True,collection=cid,delivery_manifest_sha256=sha(s/'tables/delivery_manifest.tsv'),files_collected=len(read(d/'captured_files.tsv')),all_owned_parents_closed=True,mode='Full HPC validation and lean local fit mirror').items():w.writerow(dict(key=k,value=v))
    print('Verified collection',cid,'with',len(read(d/'captured_files.tsv')),'files; full fit payloads remain on HPC')
if __name__=='__main__':main()
