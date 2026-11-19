"""Bind the one explicitly approved N20 repair validation and bounded A benchmark."""
from pathlib import Path
import csv
import hashlib
import shutil
import sys

ORIGINAL = '20261118_culture_N20_independent_01'
NAME = '20261118_culture_N20_score_and_A400_01'
SUPPORT = ('production_bindings.R','interface.R','feasibility.R','search_helpers.R','core.R',
           'fit.R','n20_model.R','policy.py','common.py','verify_source.py','seal_phase.py','runner.sh')
ORIGINAL_CASE_SHA = 'd555fbe4ad52395d10dff169ae613035d8dfdd05da4869c6b9b28af12df263ed'

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def read(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def table(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def main(parent, patch):
    parent=Path(parent).resolve()
    patch=Path(patch).resolve()
    source=parent/ORIGINAL
    target=parent/NAME
    previous={r['path']:r for r in read(source/'source_manifest.tsv')}
    ledger=read(source/'accounting/jobs.tsv')
    budget,=read(source/'accounting/final_budget.tsv')
    if source.joinpath('controller/TERMINAL').read_text().strip()!='VALIDATION_HOLD' or budget['status']!='VALIDATION_HOLD':
        raise RuntimeError('Original N20 study is not a closed validation hold')
    if {j['job_id'] for j in ledger}!={'22058898','22058906','22058909','22058943'}:
        raise RuntimeError('Original allocation set differs from the approved archive')
    if (budget['active_owned_jobs']!='0' or budget['settled_core_seconds']!='119' or
        budget['attempts']!='2' or budget['technical_retries']!='1' or
        sum(int(j['settled_seconds']) for j in ledger)!=119):
        raise RuntimeError('Original settlement differs from the approved history')
    if sha(source/'inputs/public/N20.rds')!=ORIGINAL_CASE_SHA:
        raise RuntimeError('Frozen noisy case changed')
    target.mkdir()
    for rel in ('release','inputs/public','inputs/private','verification/T00/identities',
                'verification/T00/rng','accounting','logs','continuity'):
        (target/rel).mkdir(parents=True)
    provenance=[]
    def copy(rel,rec=None,dest=None):
        source_path=source/rel
        path=dest or rel
        r=rec or previous[rel]
        d=sha(source_path)
        n=source_path.stat().st_size
        if d!=r['sha256'] or n!=int(r['bytes']):
            raise RuntimeError('Source digest mismatch: '+rel)
        (target/path).parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source_path,target/path)
        provenance.append(dict(archive_path=rel,supplement_path=path,bytes=n,sha256=d))

    for rel in sorted(previous):
        if rel.startswith(('inputs/production/','runtime/')):
            copy(rel)
    for p in (source/'runtime').rglob('*'):
        if p.is_symlink():
            q=target/p.relative_to(source)
            if p.resolve().parent!=p.parent.resolve() or not (target/p.resolve().relative_to(source)).exists():
                raise RuntimeError('Unexpected runtime symlink')
            q.symlink_to(p.readlink())
    for rel in ('inputs/public/N20.rds','inputs/public/contract.rds','inputs/private/truth.rds',
                'inputs/banks/N_A.rds','inputs/start_registry.tsv'):
        copy(rel)
    for name in SUPPORT:
        copy('release/'+name)
    for name in ('read_guard.so','identities/N_A.rds','rng/N_A.rds'):
        r,=[r for r in read(source/'verification/T00/phase_manifest.tsv') if r['path']==name]
        copy('verification/T00/'+name,r)
    for phase in ('T00','G'):
        base='verification/'+phase
        directory=source/base
        marks={r['path']:r for r in read(directory/'phase_manifest.tsv')}
        for f in ('checks.tsv','policy_gate.tsv','application_result.tsv'):
            if f=='checks.tsv' and not all(r['passed']=='TRUE' for r in read(directory/f)):
                raise RuntimeError('Prior source/access gate is incomplete: '+phase)
            if f!='checks.tsv' and read(directory/f)[0]['status']!='PASS':
                raise RuntimeError('Prior source/access gate is incomplete: '+phase)
            copy(base+'/'+f,marks[f],f'continuity/{phase}/{f}')
    past=('controller/TERMINAL','accounting/final_budget.tsv','accounting/jobs.tsv',
          'accounting/technical_retry.tsv','accounting/stop_record.tsv','accounting/candidate_attempts.tsv',
          'history/order_hold_22058909/final_budget.tsv','history/order_hold_22058909/stop_record.tsv')
    for name in past:
        f=source/name
        copy(name,dict(bytes=f.stat().st_size,sha256=sha(f)),'continuity/'+name)
    for phase in ('E','E__attempt1_22058909'):
        path=source/'verification'/phase
        marks={r['path']:r for r in read(path/'phase_manifest.tsv')}
        for name in ('checks.tsv','console.log','phase_exit_codes.tsv','forward_attempts.tsv'):
            copy('verification/'+phase+'/'+name,marks[name],f'continuity/{phase}/{name}')
    for name in ('main.R','production_score.R','engineering.R','launch_scope.sh','schedule.py',
                 'local_production_score_regression.tsv'):
        dst=target/'release'/name
        shutil.copy2(patch/name,dst)
    table(target/'release/source_and_history_manifest.tsv',provenance)
    manifest=[]
    for sub in ('release','inputs','runtime','verification','continuity'):
        for p in sorted((target/sub).rglob('*')):
            if p.is_file() and not p.is_symlink():
                manifest.append(dict(path=p.relative_to(target).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    table(target/'source_manifest.tsv',manifest)
    table(target/'accounting/authority.tsv',[dict(scope='One engineering allocation, then one A benchmark allocation after acceptance',
        previous_settled_core_seconds=119,previous_charged_forwards=2,total_core_second_cap=662400,
        engineering_allocations=1,benchmark_allocations=1,benchmark_forwards=400,
        other_blind_chains=0,profiles=0,tumor=0)])
    print(target)

if __name__=='__main__':
    main(*sys.argv[1:])
