"""Reconcile the N20 engineering and A400 receipts with native Slurm records."""
from pathlib import Path
import shlex
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'culture_D19_independent'))
import collect_closeout as t

t.REMOTE='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/synthetic_identifiability/20261118_culture_N20_score_and_A400_01'
COMMON=('phase_manifest.tsv','checks.tsv','application_result.tsv','policy_gate.tsv',
        'phase_exit_codes.tsv','image_sha256.txt','loaded_namespaces.tsv')
PHASES={
    'E':('forward_attempts.tsv','sample_alignment.tsv','scoring_agreement.tsv'),
    'N_A_import':('requested_population.rds','evaluations.tsv'),
    'N_A_bench':('requested_population.rds','actual_initial_population.rds','evaluations.tsv','callback_timings.tsv'),
}
BUDGET=('accounting/authority.tsv','accounting/jobs.tsv','accounting/cumulative_budget.tsv',
        'accounting/benchmark_N.tsv','accounting/slurm_snapshot.txt')
HISTORY=('continuity/accounting/jobs.tsv','continuity/accounting/final_budget.tsv',
         'continuity/accounting/technical_retry.tsv','continuity/accounting/stop_record.tsv',
         'continuity/controller/TERMINAL','continuity/E/console.log',
         'continuity/E__attempt1_22058909/console.log')

def main(path):
    out=Path(path).resolve()
    out.mkdir(parents=True,exist_ok=True)
    base=('source_manifest.tsv','release/source_and_history_manifest.tsv','release/local_production_score_regression.tsv')
    t.collect(base+BUDGET+HISTORY,out)
    src=t.table(out/'source_manifest.tsv')
    for x in base[1:]+HISTORY:t.verified_file(out/x,out,src)
    jobs=t.table(out/'accounting/jobs.tsv')
    if [j['action'] for j in jobs]!=['engineering','benchmark']:
        raise RuntimeError('Unexpected scope')
    if any((r['state'],r['exit_code'])!=('COMPLETED','0:0') for r in jobs):
        raise RuntimeError('Allocation is incomplete')
    old=t.table(out/'continuity/accounting/jobs.tsv')
    ids=[j['job_id'] for j in jobs+old]
    if any(not x.isdigit() for x in ids) or len(set(ids))!=6:
        raise RuntimeError('Unexpected native allocation identities')
    command=('sacct -nP -j '+','.join(ids)+' --format=JobIDRaw,State%40,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS && '
       'printf "\\n--- LIVE ---\\n" && squeue -h -u 4482173 -o "%i"')
    native=t.ssh('bash -lc '+shlex.quote(command)).decode()
    (out/'accounting/independent_slurm_closeout.txt').write_text(native)
    settled,live=native.split('--- LIVE ---',1)
    if any(x.strip() in ids for x in live.splitlines()):
        raise RuntimeError('Historical or current owned allocation remains live')
    records={v[0]:v for line in settled.splitlines() if len(v:=line.split('|'))>=7}
    for j in old+jobs:
        id=j['job_id']
        z=records[id]
        cost=max(int(r[3] or 0)*int(r[4] or 0) for k,r in records.items() if k==id or k.startswith(id+'.'))
        expected=j.get('settled_seconds',j.get('settled_core_seconds'))
        if z[1].split()[0]!=j['state'] or z[2]!=j['exit_code'] or z[6]!='xxlarge' or cost!=int(expected):
            raise RuntimeError('Native allocation ledger mismatch: '+id)
    for phase,names in PHASES.items():
        folder='verification/'+phase
        t.collect([folder+'/'+n for n in COMMON+names],out)
        directory=out/folder
        rec=t.table(directory/'phase_manifest.tsv')
        for f in COMMON+names:
            if f!='phase_manifest.tsv':t.verified_file(directory/f,directory,rec)
        app,=t.table(directory/'application_result.tsv')
        role,=t.table(directory/'policy_gate.tsv')
        returned,=t.table(directory/'phase_exit_codes.tsv')
        expected=jobs[0 if phase=='E' else 1]['job_id']
        if app['status']!='PASS' or app['job']!=expected or role['status']!='PASS' or any(int(c) for c in returned.values()):
            raise RuntimeError('Application or policy failed: '+phase)
        if (directory/'image_sha256.txt').read_text().strip()!=t.IMAGE_SHA:
            raise RuntimeError('Container identity changed')
        if not all(r['passed']=='TRUE' for r in t.table(directory/'checks.tsv')):
            raise RuntimeError('Scientific fixture failed: '+phase)
    pay,=t.table(out/'accounting/cumulative_budget.tsv')
    if int(pay['cumulative_core_seconds'])!=119+sum(int(r['settled_core_seconds']) for r in jobs):
        raise RuntimeError('Cumulative settlement mismatch')
    if pay['live_reserved_core_seconds']!='0' or pay['historical_forwards']!='2' or pay['engineering_forwards']!='3' or pay['A_benchmark_forwards']!='400':
        raise RuntimeError('Forward attempt ledger mismatch')
    evaluation=t.table(out/'verification/N_A_bench/evaluations.tsv')
    start=[r for r in evaluation if r['event']=='START']
    finish=[r for r in evaluation if r['event']=='END']
    if len(evaluation)!=800 or len(start)!=400 or len(finish)!=400 or [
        (r['call'],r['key']) for r in start]!=[(r['call'],r['key']) for r in finish]:
        raise RuntimeError('A400 independent journal pairing failed')
    if any(r['feasible']!='TRUE' for r in finish[:200]):
        raise RuntimeError('An imported feasible start did not propagate in N20')
    print('VERIFIED_ENGINEERING_AND_A400',ids[:2],pay['cumulative_core_seconds'])

if __name__=='__main__':
    main(sys.argv[1])
