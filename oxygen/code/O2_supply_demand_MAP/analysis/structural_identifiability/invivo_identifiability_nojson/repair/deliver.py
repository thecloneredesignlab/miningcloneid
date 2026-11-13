#!/usr/bin/env python3
"""Seal and mirror derived results after closure; never submit or copy active fits."""
import argparse
import shlex
import shutil
import subprocess
import tempfile
from pathlib import Path
import workflow as w


def main():
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();s=a.study.resolve()
    if s.name!='20261111_invivo_baseline_repair_045318Z': raise ValueError('Unexpected study')
    proof={r['key']:r['value'] for r in w.read_table(s/'tables/closure_audit.tsv')}
    rows=w.read_table(s/'tables/jobs.tsv')
    if proof.get('pass')!='true' or not all(w.closed(j) for j in rows): raise ValueError('Closure required')
    if w.State(s).load()!=rows: raise ValueError('Control snapshot differs')
    for path in Path(__file__).resolve().parent.iterdir():
        if path.is_file() and path.suffix in ('.py','.R','.sh'): shutil.copy2(str(path),str(s/'implementation'/path.name))
    manifest=s/'tables/final_delivery_manifest.tsv'
    receipt=s/'tables/delivery_transport_verification.tsv'
    derived=set()
    for name in ('tables','reports','figures','implementation'):
        derived.update(path for path in (s/name).rglob('*') if path.is_file() and path not in (manifest,receipt)
                       and '.DS_Store'!=path.name and '__pycache__' not in path.parts)
    evidence=set(derived)
    generation=(s/'state/CURRENT').read_text().strip()
    evidence.update((s/'manifests/input_manifest.tsv',s/'state/CURRENT',s/'state'/generation/'manifest.tsv',s/'state'/generation/'jobs.tsv'))
    evidence.update((s/'releases').glob('*/manifest.tsv'))
    for directory in ('fits','tests'):
        for path in (s/directory).glob('*/payload_manifest.tsv'):
            w.verify_payload(path.parent);evidence.update((path,path.parent/'COMPLETE'))
    for j in rows:
        evidence.update((s/'dispatch'/j['task']).iterdir())
        for name in (j['task']+'_submission.txt',j['task']+'_receipt.txt',j['task']+'_'+j['job_id']+'.log'):
            evidence.add(s/'logs'/name)
        for name in ('image_sha256.txt','launch_exit_code.txt'): evidence.add(s/'runtime'/j['task']/name)
    entries=[dict(path=str(path.relative_to(s)),bytes=path.stat().st_size,sha256=w.sha(path)) for path in sorted(evidence)]
    w.table(manifest,entries,['path','bytes','sha256'])
    protected={s/'tables/jobs.tsv',s/'tables/raw_accounting.txt'}
    mirror=sorted(str(path.relative_to(s)) for path in derived-protected)+['tables/final_delivery_manifest.tsv']
    ssh=['ssh','-o','BatchMode=yes','-o','IdentitiesOnly=yes','-i','/Users/4482173/.ssh/id_ed25519',
         '-o','ControlMaster=auto','-o','ControlPath=/private/tmp/ivi_repair_ssh_%C','-o','ControlPersist=600']
    target='4482173@red.moffitt.org'
    remote='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/structural_identifiability/'+s.name
    with tempfile.NamedTemporaryFile(mode='w',prefix='ivi_final_',suffix='.txt') as listing:
        listing.write('\n'.join(mirror)+'\n');listing.flush()
        subprocess.run(['rsync','-az','--files-from='+listing.name,'-e',' '.join(shlex.quote(x) for x in ssh),str(s)+'/',target+':'+remote+'/'],check=True)
    script="""from pathlib import Path
import sys
s=Path(REMOTE);sys.path.insert(0,str(s/'implementation'));import workflow as w
w.verify_manifest(s,s/'tables/final_delivery_manifest.tsv');w.verify_manifest(s,s/'manifests/input_manifest.tsv')
rows=w.State(s).load();assert rows==w.read_table(s/'tables/jobs.tsv') and all(w.closed(j) for j in rows)
for release in set(j['release'] for j in rows): w.verify_manifest(s/'releases'/release,s/'releases'/release/'manifest.tsv')
for scope in ('fits','tests'):
 for manifest in (s/scope).glob('*/payload_manifest.tsv'): w.verify_payload(manifest.parent)
print('Remote final bundle, source, releases and complete payload fingerprints verified.')
""".replace('REMOTE',repr(remote))
    result=subprocess.run(ssh+[target,'python3','-'],input=script,universal_newlines=True,stdout=subprocess.PIPE,check=True)
    w.table(receipt,[dict(key='pass',value='true'),dict(key='manifest_sha256',value=w.sha(manifest)),
                    dict(key='bound_metadata_and_payload_pointers',value=str(len(entries))),dict(key='created_utc',value=w.utc())],['key','value'])
    subprocess.run(['rsync','-az','-e',' '.join(shlex.quote(x) for x in ssh),str(receipt),target+':'+remote+'/tables/'],check=True)
    print(result.stdout.strip());print('Mirrored',len(mirror),'derived files; sealed',len(entries),'metadata and transitive payload/source pointers. No state or numerical-payload write.')


if __name__=='__main__': main()
