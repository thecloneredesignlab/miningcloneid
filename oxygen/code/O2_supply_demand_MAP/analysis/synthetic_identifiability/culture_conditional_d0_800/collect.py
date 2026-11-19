"""Verify a selected local mirror after canonical HPC completion; no fit calls."""
import argparse,re,shlex,shutil,subprocess,sys,os
from pathlib import Path,PurePosixPath
import workflow as w
REMOTE='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/synthetic_identifiability/20261113_culture_conditional_d0_seed144_02'
REMOTE=os.environ.get('D0_COLLECT_REMOTE_ROOT',REMOTE)
SSH=['ssh','-o','BatchMode=yes','-o','IdentitiesOnly=yes','-o','ConnectTimeout=15','-i','/Users/4482173/.ssh/id_ed25519']
HOST='4482173@red.moffitt.org'
RESERVE=2*1024**3
def eligible(row):
    p=PurePosixPath(row['path'])
    if p.is_absolute() or '..' in p.parts:raise ValueError('Unsafe manifest path')
    if p.suffix.lower() in {'.json','.jsonl','.ndjson','.js','.java','.jar','.class'}:raise ValueError('Forbidden artifact')
    return not (any(x in {'runtime','state','__pycache__','improvements'} for x in p.parts) or p.parts[0]=='collection' or re.fullmatch(r'best_\d+[.]rds',p.name))
def admitted(rows,free):
    selected=[r for r in rows if eligible(r)]
    total=sum(int(r['bytes']) for r in selected)
    if total+RESERVE>free:raise ValueError('Selected payload would violate local2GiB reserve')
    return selected,total
def remote_text(relative):
    return subprocess.run(SSH+[HOST,'cat '+shlex.quote(REMOTE+'/'+relative)],check=True,stdout=subprocess.PIPE,universal_newlines=True).stdout
def collect(project):
    project=Path(project).resolve();s=project/'oxygen/results/synthetic_identifiability'/Path(REMOTE).name
    out=s/'collection';out.mkdir(parents=True,exist_ok=True)
    marker=remote_text('COMPLETE').strip()
    if not re.fullmatch('[a-f0-9]{64}',marker):raise ValueError('HPC completion marker invalid')
    manifest=remote_text('tables/delivery_manifest.tsv');w.text(out/'hpc_full_delivery_manifest.tsv',manifest)
    if w.sha(out/'hpc_full_delivery_manifest.tsv')!=marker:raise ValueError('HPC full manifest does not match completion marker')
    rows=w.read(out/'hpc_full_delivery_manifest.tsv');selected,total=admitted(rows,shutil.disk_usage(str(project)).free)
    receipt=out/'LOCAL_COLLECTION_COMPLETE.tsv'
    if receipt.exists():
        if w.kv(receipt)['hpc_manifest_sha256']!=marker:raise ValueError('Closed local collection points at another HPC package')
        w.verify(s,out/'selected_manifest.tsv');print('Selected local mirror already verified');return
    w.table(out/'selected_manifest.tsv',selected,['path','sha256','bytes'])
    w.text(out/'selected_files.txt','\n'.join(r['path'] for r in selected)+'\n')
    subprocess.run(['rsync','-az','--partial','--files-from='+str(out/'selected_files.txt'),'-e',' '.join(shlex.quote(x) for x in SSH),HOST+':'+REMOTE+'/',str(s)+'/'],check=True)
    w.verify(s,out/'selected_manifest.tsv')
    texts=0
    for r in selected:
        p=s/r['path']
        if p.suffix in {'.R','.py','.sh','.md','.tsv','.csv','.txt'}:
            content=p.read_text(encoding='utf-8');texts+=1
            if re.search('[\u3400-\u4dbf\u4e00-\u9fff]',content):raise ValueError('Non-English CJK text in new artifact: '+r['path'])
    free=shutil.disk_usage(str(project)).free
    if free<RESERVE:raise ValueError('Post-transfer local reserve violated')
    w.typed(receipt,dict(pass_gate=True,hpc_manifest_sha256=marker,selected_manifest_sha256=w.sha(out/'selected_manifest.tsv'),selected_files=len(selected),selected_bytes=total,english_text_artifacts=texts,local_free_bytes=free,local_scope='Selected scientific payload plus full canonical HPC inventory',all_HPC_files_copied=False,large_improvement_histories='Canonical on HPC; hash inventory delivered',collector_sha256=w.sha(Path(__file__)),utc=w.utc()))
    print('Verified selected local delivery:',len(selected),'files,',total,'bytes; full checkpoints remain canonical on HPC')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('project',type=Path);collect(p.parse_args().project)
