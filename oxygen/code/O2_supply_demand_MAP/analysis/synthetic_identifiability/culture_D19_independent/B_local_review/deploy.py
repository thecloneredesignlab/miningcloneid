"""Create a new, read-only scientific input set from the closed D19 archive."""
from pathlib import Path
import csv
import hashlib
import shutil
import sys

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for data in iter(lambda:f.read(1024*1024), b''):
            h.update(data)
    return h.hexdigest()

def table(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(p, rows):
    with p.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

KEY = '3691ca93d3c4d115dadf9afe0a699e76406149978680eb8056e646835f6bf716'
SRC_NAME = '20261118_culture_D19_independent_01'
DEST_NAME = '20261118_culture_D19_B_lowest_local_review_01'
SUPPORT = ('production_bindings.R','interface.R','feasibility.R','search_helpers.R','core.R',
           'joint_numerics.R','constraints.R','policy.py','verify_source.py','seal_phase.py','common.py')

def main(parent, supplied):
    parent = Path(parent).resolve()
    old = parent/SRC_NAME
    dest = parent/DEST_NAME
    prev = {r['path']:r for r in table(old/'source_manifest.tsv')}
    lock, = table(old/'reports/selection_lock.tsv')
    budget, = table(old/'accounting/final_budget.tsv')
    if lock['best_found_key']!=KEY or lock['best_found_origin']!='D_B' or lock['best_found_replay_verified']!='True':
        raise ValueError('The closed lowest point changed')
    if budget['active_owned_jobs']!='0' or int(budget['settled_core_seconds'])!=68747:
        raise ValueError('The original closed settlement changed')
    if old.joinpath('controller/TERMINAL').read_text().strip()!=budget['status']:
        raise ValueError('Inconsistent closed receipt')
    dest.mkdir()
    for rel in ('release','inputs/public','verification/D_B_anchor','verification/T00/identities','accounting','logs'):
        (dest/rel).mkdir(parents=True)
    provenance=[]
    def copy(rel, out=None, receipt=None):
        src=old/rel
        out=out or rel
        expected=receipt or prev[rel]
        h=sha(src)
        if h!=expected['sha256'] or src.stat().st_size!=int(expected['bytes']):
            raise ValueError('Closed source digest mismatch: '+rel)
        (dest/out).parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(src,dest/out)
        provenance.append(dict(archive_path=rel,supplement_path=out,sha256=h,bytes=src.stat().st_size))

    for rel in sorted(prev):
        if rel.startswith(('inputs/production/','runtime/')):
            copy(rel)
    for p in (old/'runtime').rglob('*'):
        if p.is_symlink():
            out=dest/p.relative_to(old)
            if p.resolve().parent!=p.parent.resolve() or not (dest/p.resolve().relative_to(old)).exists():
                raise ValueError('Unexpected runtime symlink')
            out.symlink_to(p.readlink())
    for name in SUPPORT:
        copy('release/'+name)
    for name in ('D19.rds','contract.rds'):
        copy('inputs/public/'+name)
    for phase, name in (('D_B_replay2a','replay_a.rds'),('D_B_replay2b','replay_b.rds'),('T00','read_guard.so'),('T00','identities/D_B.rds')):
        file='point_'+KEY+'.rds' if phase.startswith('D_B') else name
        path='verification/'+phase+'/'+file
        expected,=[r for r in table(old/'verification'/phase/'phase_manifest.tsv') if r['path']==file]
        output='verification/D_B_anchor/'+name if phase.startswith('D_B') else path
        copy(path,output,expected)
    for p in Path(supplied).iterdir():
        if p.is_file():
            target=dest/'release'/p.name
            if target.exists(): raise ValueError('Cannot replace production release file')
            shutil.copy2(p,target)
    write(dest/'release/source_and_anchor_manifest.tsv',provenance)
    manifest=[]
    for sub in ('inputs','release','runtime','verification'):
        for p in sorted((dest/sub).rglob('*')):
            if p.is_file() and not p.is_symlink():
                manifest.append(dict(path=p.relative_to(dest).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    write(dest/'source_manifest.tsv',manifest)
    print(dest)

if __name__=='__main__':
    main(*sys.argv[1:])
