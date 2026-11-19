"""Assemble the immutable single-case release from verified production artifacts."""
from pathlib import Path
import csv
import hashlib
import shutil
import sys

NAME = '20261118_culture_N20_independent_01'
CASE_SHA = 'd555fbe4ad52395d10dff169ae613035d8dfdd05da4869c6b9b28af12df263ed'
NOISY_SHA = 'cfe0803a09baf31c5c3244a752ca899395185301f57e939a148ba38e30ac9205'
SUPPORT = ('production_bindings.R','interface.R','feasibility.R','search_helpers.R',
           'verify_source.py','seal_phase.py','common.py')

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def read(path):
    with path.open(newline='') as stream:return list(csv.DictReader(stream,delimiter='\t'))

def table(path,rows):
    with path.open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def main(results,bundle):
    parent=Path(results).resolve()
    bundle=Path(bundle).resolve()
    old=parent/'20261118_culture_D19_independent_01'
    archive={r['path']:r for r in read(old/'source_manifest.tsv')}
    case=parent/'20261116_culture_growth_karyotype_recovery_01'/'verification/I01/cases/N20.rds'
    noisy=parent/'20261116_culture_growth_karyotype_conditional_01'/'public/synthetic_dataset.rds'
    if sha(case)!=CASE_SHA or sha(noisy)!=NOISY_SHA:
        raise RuntimeError('Frozen noisy input hash mismatch')
    gate,=read(bundle/'N20_review_gate.tsv')
    if gate['status']!='ELIGIBLE_FOR_SEPARATE_N20_REVIEW_AT_REGISTERED_D_RESOLUTION' or not (
       gate['anchor_poll_complete']=='TRUE' and int(gate['prediction_origins'])==3 and gate['N20_released']=='FALSE'):
        raise RuntimeError('Upstream resolution-aware review gate is not eligible')
    dest=parent/NAME
    dest.mkdir()
    for sub in ('verification','controller','accounting','logs','reports',
                'inputs/public','inputs/private','inputs/canaries'):
        (dest/sub).mkdir(parents=True)
    shutil.copytree(bundle/'release',dest/'release')
    shutil.copy2(bundle/'N20_review_gate.tsv',dest/'release/')
    receipts=[]
    def copy(rel):
        src=old/rel
        t=dest/rel
        row=archive[rel]
        if sha(src)!=row['sha256'] or src.stat().st_size!=int(row['bytes']):
            raise RuntimeError('Production-source archive mismatch: '+rel)
        t.parent.mkdir(parents=True,exist_ok=True)
        if t.exists():raise RuntimeError('Could not preserve inherited production source: '+rel)
        shutil.copy2(src,t)
        receipts.append(dict(path=rel,source='closed D19 immutable source archive',sha256=row['sha256']))
    for rel in sorted(archive):
        if rel.startswith(('inputs/production/','runtime/')):
            copy(rel)
    for p in (old/'runtime').rglob('*'):
        if p.is_symlink():
            out=dest/p.relative_to(old)
            if p.resolve().parent!=p.parent.resolve() or not (dest/p.resolve().relative_to(old)).exists():
                raise RuntimeError('Unexpected runtime symlink')
            out.symlink_to(p.readlink())
    for name in SUPPORT:
        copy('release/'+name)
    for rel in ('inputs/public/contract.rds','inputs/private/truth.rds','inputs/private/reference.rds'):
        copy(rel)
    shutil.copy2(case,dest/'inputs/public/N20.rds')
    receipts.append(dict(path='inputs/public/N20.rds',source='frozen N20 case',sha256=CASE_SHA))
    receipts.append(dict(path='not mounted: frozen synthetic_dataset.rds',source='upstream archived draw',sha256=NOISY_SHA))
    table(dest/'release/production_source_receipts.tsv',receipts)
    for name in ('banks','rng'):
        shutil.copytree(bundle/'starts'/name,dest/'inputs'/name)
    shutil.copy2(bundle/'starts/start_registry.tsv',dest/'inputs/start_registry.tsv')
    (dest/'inputs/canaries/alias_truth.rds').symlink_to('../private/truth.rds')
    (dest/'inputs/canaries/alias_B.rds').symlink_to('../banks/N_B.rds')
    (dest/'verification/N_B_fixture').mkdir()
    (dest/'verification/N_B_fixture/state.rds').write_text('Foreign origin path only.\n')
    rows=[]
    for sub in ('inputs','runtime','release'):
        for p in sorted((dest/sub).rglob('*')):
            if p.is_file() and not p.is_symlink():
                rows.append(dict(path=p.relative_to(dest).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    table(dest/'source_manifest.tsv',rows)
    table(dest/'accounting/authority.tsv',[dict(scope='N20 single frozen noisy dataset only',core_seconds=662400,
        max_full_forwards=504400,technical_retries=1,reserve_rounds=1,profiles=0,tumor=0,repeated_noisy=0)])
    print(dest)

if __name__=='__main__':
    main(*sys.argv[1:])
