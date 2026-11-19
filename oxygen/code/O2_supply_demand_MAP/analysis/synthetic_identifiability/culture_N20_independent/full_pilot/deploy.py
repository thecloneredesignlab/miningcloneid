"""Assemble the authorized pilot with A's exact archived checkpoint and ledgers."""
from pathlib import Path
import csv
import hashlib
import shutil
import sys

NAME = '20261118_culture_N20_full_01'
OLD_NAME = '20261118_culture_N20_independent_01'
BENCH_NAME = '20261118_culture_N20_score_and_A400_01'
SUPPORT = ('production_bindings.R','interface.R','feasibility.R','search_helpers.R',
           'n20_model.R','fit.R','production_score.R','policy.py','common.py',
           'verify_source.py','seal_phase.py','runner.sh')
CASE_SHA = 'd555fbe4ad52395d10dff169ae613035d8dfdd05da4869c6b9b28af12df263ed'
PHASES = {
    'E': ('checks.tsv','scoring_agreement.tsv','sample_alignment.tsv','forward_attempts.tsv'),
    'N_A_import': ('state.rds','requested_population.rds','checks.tsv','evaluations.tsv'),
    'N_A_bench': ('state.rds','requested_population.rds','actual_initial_population.rds',
                  'evaluations.tsv','callback_timings.tsv','checks.tsv','result.tsv'),
}

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for data in iter(lambda:f.read(1024*1024),b''):h.update(data)
    return h.hexdigest()

def read(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def table(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def manifest(path):
    return {r['path']:r for r in read(path)}

def main(parent,release):
    parent=Path(parent).resolve()
    release=Path(release).resolve()
    old=parent/OLD_NAME
    bench=parent/BENCH_NAME
    root=parent/NAME
    original=manifest(old/'source_manifest.tsv')
    admitted=manifest(bench/'source_manifest.tsv')
    settled,=read(bench/'accounting/cumulative_budget.tsv')
    b,=read(bench/'accounting/benchmark_N.tsv')
    jobs=read(bench/'accounting/jobs.tsv')
    history=read(old/'accounting/jobs.tsv')
    if (settled['cumulative_core_seconds']!='369' or settled['live_reserved_core_seconds']!='0'
        or settled['A_benchmark_forwards']!='400' or settled['engineering_forwards']!='3'
        or settled['historical_forwards']!='2' or b['status']!='PASS'
        or b['job_id']!='22059685' or b['continuation_released']!='False'
        or [j['job_id'] for j in jobs]!=['22059657','22059685']
        or any((j['state'],j['exit_code'])!=('COMPLETED','0:0') for j in jobs)
        or len(history)!=4 or sum(int(j['settled_seconds']) for j in history)!=119
        or sum(int(j['settled_core_seconds']) for j in jobs)!=250):
        raise RuntimeError('Archived cost gate or settled history changed')
    if sha(bench/'inputs/public/N20.rds')!=CASE_SHA or sha(old/'inputs/public/N20.rds')!=CASE_SHA:
        raise RuntimeError('Frozen public case changed')
    for rel in ('inputs/public/contract.rds','inputs/banks/N_A.rds','inputs/start_registry.tsv'):
        if admitted[rel]['sha256']!=original[rel]['sha256']:
            raise RuntimeError('Handoff identity conflicts with the original: '+rel)
    root.mkdir()
    for sub in ('release','inputs','runtime','verification','controller','reports','accounting','logs','continuity'):
        (root/sub).mkdir()
    provenance=[]
    def put(base,rel,receipt=None,destination=None):
        q=base/rel
        digest=sha(q)
        size=q.stat().st_size
        if receipt and (digest!=receipt['sha256'] or size!=int(receipt['bytes'])):
            raise RuntimeError('Sealed source mismatch: '+rel)
        dest=destination or rel
        target=root/dest
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():raise RuntimeError('Duplicate staged path: '+dest)
        shutil.copy2(q,target)
        provenance.append(dict(archive=base.name,archive_path=rel,continuation_path=dest,
                               bytes=size,sha256=digest))
        return digest
    for rel in sorted(admitted):
        if rel.startswith(('inputs/production/','runtime/')):
            put(bench,rel,admitted[rel])
    for p in (bench/'runtime').rglob('*'):
        if p.is_symlink():
            q=root/p.relative_to(bench)
            if p.resolve().parent!=p.parent.resolve() or not (root/p.resolve().relative_to(bench)).exists():
                raise RuntimeError('Unexpected runtime symlink')
            q.symlink_to(p.readlink())
    for name in SUPPORT:
        r='release/'+name
        put(bench,r,admitted[r])
    for r in ('inputs/public/N20.rds','inputs/public/contract.rds','inputs/start_registry.tsv','inputs/banks/N_A.rds'):
        put(bench,r,admitted[r])
    for r in ('inputs/private/reference.rds','inputs/private/truth.rds',*('inputs/banks/N_'+x+'.rds' for x in 'BCDEF')):
        put(old,r,original[r])
    t00=manifest(old/'verification/T00/phase_manifest.tsv')
    for r in ('read_guard.so',*(f'{d}/N_{x}.rds' for d in ('identities','rng') for x in 'ABCDEF')):
        if r in ('read_guard.so','identities/N_A.rds','rng/N_A.rds') and t00[r]['sha256']!=admitted['verification/T00/'+r]['sha256']:
            raise RuntimeError('A checkpoint parameter mapping, RNG, or guard changed')
        put(old,'verification/T00/'+r,t00[r])
    for phase,files in PHASES.items():
        base='verification/'+phase
        records=manifest(bench/base/'phase_manifest.tsv')
        expected='22059657' if phase=='E' else '22059685'
        for name in ('application_result.tsv','policy_gate.tsv','phase_exit_codes.tsv'):
            r,=read(bench/base/name)
            if (name=='phase_exit_codes.tsv' and any(int(v) for v in r.values())) or (
               name!='phase_exit_codes.tsv' and (r['status']!='PASS' or (name=='application_result.tsv' and r['job']!=expected))):
                raise RuntimeError('Handoff gate failed: '+phase+'/'+name)
            put(bench,base+'/'+name,records[name],'continuity/'+phase+'/'+name)
        put(bench,base+'/phase_manifest.tsv',destination='continuity/'+phase+'/phase_manifest.tsv')
        for name in files:
            dest=base+'/'+name if phase!='E' else 'continuity/'+phase+'/'+name
            put(bench,base+'/'+name,records[name],dest)
        checks=read(bench/base/'checks.tsv')
        if not checks or not all(r['passed']=='TRUE' for r in checks):
            raise RuntimeError('Handoff checks failed: '+phase)
    result={r['key']:r['value'] for r in read(root/'verification/N_A_bench/result.tsv')}
    if (result.get('status')!='READY' or result.get('origin')!='N_A' or result.get('stage')!='bench'
        or result.get('calls')!='400' or result.get('native_evaluations')!='400'
        or result.get('generations')!='1'):
        raise RuntimeError('A400 state receipt changed')
    for r in ('accounting/jobs.tsv','accounting/final_budget.tsv','accounting/technical_retry.tsv',
              'controller/TERMINAL','verification/E/console.log','verification/E__attempt1_22058909/console.log'):
        put(old,r,destination='continuity/original/'+r)
    for r in ('benchmark_N.tsv','cumulative_budget.tsv','jobs.tsv'):
        put(bench,'accounting/'+r,destination='continuity/'+r)
    for p in release.iterdir():
        if p.is_file():
            target=root/'release'/p.name
            if target.exists():raise RuntimeError('Cannot replace frozen implementation: '+p.name)
            shutil.copy2(p,target)
    table(root/'continuity/handoff_manifest.tsv',provenance)
    rows=[]
    for subtree in ('inputs','runtime','release','verification','continuity'):
        for p in sorted((root/subtree).rglob('*')):
            if p.is_file() and not p.is_symlink():
                rows.append(dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    table(root/'source_manifest.tsv',rows)
    table(root/'accounting/authority.tsv',[dict(
        scope='Complete one frozen N20 pilot with four primary origins and conditional E/F',
        historical_core_seconds=369,historical_engineering_forwards=5,A_checkpoint_forwards=400,
        core_hour_cap='NONE',max_forwards=504400,reserve_rounds=1,
        profiles=0,tumor=0,repeated_datasets=0)])
    print(root)

if __name__=='__main__':
    main(*sys.argv[1:])
