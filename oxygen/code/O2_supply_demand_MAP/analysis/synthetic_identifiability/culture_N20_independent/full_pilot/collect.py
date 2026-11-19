"""Verify and collect the native N20 pilot and its archived accounting continuity."""
from pathlib import Path
import re
import shlex
import sys
import tarfile
import time
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'culture_D19_independent'))
import collect_closeout as transfer

transfer.REMOTE='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/synthetic_identifiability/20261118_culture_N20_full_01'
SUCCESS={'N20_BOUNDED_OPTIMIZER_QUALIFIED','NOISY_ESTIMATE_PROVISIONAL_OPTIMIZER_UNRESOLVED'}
CONTINUITY=('continuity/handoff_manifest.tsv','continuity/benchmark_N.tsv','continuity/cumulative_budget.tsv',
            'continuity/jobs.tsv','continuity/original/accounting/jobs.tsv',
            'continuity/original/accounting/final_budget.tsv','continuity/original/accounting/technical_retry.tsv',
            'continuity/original/verification/E/console.log',
            'continuity/original/verification/E__attempt1_22058909/console.log')
ASESS=('application_result.tsv','all20_parameter_recovery.tsv','prediction_recovery.tsv',
       'response_recovery.tsv','sigma_kary_and_growth_diagnostics.tsv','root_distribution_comparison.tsv')
PHASE_RECEIPTS=('application_result.tsv','gate.tsv','phase_exit_codes.tsv','policy_gate.tsv')

def copy(paths, root):
    for attempt in range(3):
        try:
            transfer.collect(paths, root)
            return
        except (tarfile.ReadError, OSError, RuntimeError):
            if attempt==2:raise
            time.sleep(2)

def verified_phase(root,folder,names):
    copy([f'{folder}/{f}' for f in ('phase_manifest.tsv',)+tuple(names)],root)
    phase=root/folder
    manifest=transfer.table(phase/'phase_manifest.tsv')
    for f in names:transfer.verified_file(phase/f,phase,manifest)

def main(out):
    root=Path(out).resolve()
    root.mkdir(parents=True,exist_ok=True)
    base=('REPORT.md','controller/TERMINAL','source_manifest.tsv','release/search_helpers.R',
          'inputs/public/N20.rds','verification/N_A_bench/evaluations.tsv')
    paths=transfer.ssh(transfer.in_study('find accounting reports -maxdepth 1 -type f -printf "%p\\n"')).decode().splitlines()
    copy(base+CONTINUITY+tuple(f for f in paths if f.endswith(('.tsv','.txt'))),root)
    terminal=(root/'controller/TERMINAL').read_text().strip()
    for f in CONTINUITY+base[3:]:
        transfer.verified_file(root/f,root,transfer.table(root/'source_manifest.tsv'))
    fresh=transfer.ssh(transfer.in_study('/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S release/verify_source.py .')).decode()
    if 'IMMUTABLE_SOURCE_PASS' not in fresh:
        raise RuntimeError('Continuation source failed remote integrity')
    jobs=transfer.table(root/'accounting/jobs.tsv')
    history=transfer.table(root/'continuity/original/accounting/jobs.tsv')
    benchmark=transfer.table(root/'continuity/jobs.tsv')
    budget,=transfer.table(root/'accounting/final_budget.tsv')
    all_jobs=history+benchmark+jobs
    ids=[r['job_id'] for r in all_jobs]
    if any(not x.isdigit() for x in ids) or len(ids)!=len(set(ids)):
        raise RuntimeError('Unreconciled allocation ID')
    cmd=('sacct -nP -j '+','.join(ids)+' --format=JobIDRaw,State%40,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,QOS && '
         'printf "\\n--- LIVE ---\\n" && squeue -h -u 4482173 -o "%i"')
    native=transfer.ssh('bash -lc '+shlex.quote(cmd)).decode()
    (root/'accounting/independent_slurm_closeout.txt').write_text(native)
    rows,live=native.split('--- LIVE ---',1)
    if any(x.strip() in ids for x in live.splitlines()):
        raise RuntimeError('An owned allocation is still live')
    actual={z[0]:z for line in rows.splitlines() if len(z:=line.split('|'))>=7}
    settled=0
    for j in all_jobs:
        id=j['job_id']
        r=actual[id]
        s=max(int(t[3] or 0)*int(t[4] or 0) for k,t in actual.items() if k==id or k.startswith(id+'.'))
        if r[1].split()[0]!=j['state'] or r[2]!=j['exit_code'] or r[6]!='xxlarge' or (
            s!=int(j.get('settled_seconds',j.get('settled_core_seconds')))):
            raise RuntimeError('Slurm settlement mismatch: '+id)
        settled+=s
    attempts=transfer.table(root/'accounting/candidate_attempts.tsv')
    if settled!=int(budget['settled_core_seconds']) or budget['core_hour_cap']!='NONE' or int(budget['attempts'])!=sum(int(x['attempts']) for x in attempts):
        raise RuntimeError('Pilot accounting mismatch')
    if int(budget['attempts'])>504400 or int(budget['historical_core_seconds'])!=369 or budget['status']!=terminal:
        raise RuntimeError('Attempt scope or historic carryforward mismatch')
    if sum(r['event']=='START' for r in transfer.table(root/'verification/N_A_bench/evaluations.tsv'))!=400:
        raise RuntimeError('A400 charged handoff missing')
    retry_path=root/'accounting/infrastructure_retry.tsv'
    if retry_path.is_file():
        retry,=transfer.table(retry_path)
        prior=retry['prior_job']
        earlier=[j for j in jobs if j['job_id']==prior and j['state']=='FAILED' and j['settled_seconds']=='28']
        recovered=[j for j in jobs if j['purpose']=='infrastructure_retry' and j['task']=='N_A_local' and j['attempt']=='2']
        if len(earlier)!=1 or len(recovered)!=1 or budget['infrastructure_retries']!='1':
            raise RuntimeError('NSS recovery does not match the owner ledger')
        archived=retry['archived_phase']
        hold=retry['archived_hold']
        patch='recovery_nss_22065165'
        copy([f'{hold}/{n}' for n in ('TERMINAL','stop_record.tsv','final_budget.tsv','jobs.tsv','REPORT.md')]
              +[f'{patch}/{n}' for n in ('manifest.tsv','launcher.sh','controller.py','resume.py')],root)
        held_jobs=transfer.table(root/hold/'jobs.tsv')
        for j in held_jobs:
            actual_job=next(k for k in jobs if k['job_id']==j['job_id'])
            if any(actual_job[c]!=j[c] for c in ('task','attempt','state','exit_code','settled_seconds')):
                raise RuntimeError('Original hold allocation changed: '+j['job_id'])
        held_ids={j['job_id'] for j in held_jobs}
        if any(j['job_id'] not in held_ids and re.fullmatch(r'N_[A-D]_(import|bench|de[1-6]b?)',j['task']) for j in jobs):
            raise RuntimeError('A recovered primary origin repeated a global phase')
        if any(sum(int(r['attempts']) for r in attempts if r['origin']=='N_'+letter and r['bucket']=='global')!=73400
               for letter in 'ABCD'):
            raise RuntimeError('Primary global evaluation counts differ from the hold')
        patch_rows=transfer.table(root/patch/'manifest.tsv')
        for n in ('launcher.sh','controller.py','resume.py'):
            transfer.verified_file(root/patch/n,root/patch,patch_rows)
        verified_phase(root,archived,('console.log','phase_exit_codes.tsv','policy_gate.tsv'))
        if recovered[0]['state']=='COMPLETED':
            verified_phase(root,'verification/N_A_local',('node_preflight.tsv','evaluations.tsv')+PHASE_RECEIPTS)
            verified_phase(root,'verification/N_A_de6',('best_identity.tsv',))
            node,=transfer.table(root/'verification/N_A_local/node_preflight.tsv')
            if (node['status']!='PASS' or node['job_id']!=recovered[0]['job_id'] or
                node['node']==retry['prior_node'] or node['uid']!=node['container_uid'] or node['host_user']!='4482173'):
                raise RuntimeError('Unverified replacement node')
            checkpoint,=transfer.table(root/'verification/N_A_de6/best_identity.tsv')
            first,second=transfer.table(root/'verification/N_A_local/evaluations.tsv')[:2]
            if ((first['event'],second['event'])!=('START','END') or
                any(x['call']!='1' or x['bucket']!='handoff' or x['key']!=checkpoint['key'] for x in (first,second)) or
                abs(float(second['score'])-float(checkpoint['Q']))>1e-8):
                raise RuntimeError('Replacement local worker did not reproduce its archived checkpoint')
    if terminal in SUCCESS:
        q=transfer.table(root/'reports/optimizer_quality.tsv')
        groups={}
        origins=set()
        for row in q:
            key=row['key'];origin=row['origin'];origins.add(origin)
            if not re.fullmatch('[0-9a-f]{64}',key) or not re.fullmatch('N_[A-F]',origin):
                raise RuntimeError('Unrecognized optimized vector')
            for tag in ('first_process_task','second_process_task'):
                task=row[tag]
                if not re.fullmatch(re.escape(origin)+r'_replay[12][ab]',task):
                    raise RuntimeError('Unrecognized replay activity')
                groups.setdefault('verification/'+task,set()).add('point_'+key+'.rds')
        for origin in origins:
            for i in (1,2):
                groups[f'verification/{origin}_verify{i}']={'checked_representative.rds','verification_polls.tsv','variance_stationarity.tsv'}
        for f,names in sorted(groups.items()):
            verified_phase(root,f,sorted(names|set(PHASE_RECEIPTS)))
        verified_phase(root,'verification/P',ASESS)
        p,=transfer.table(root/'verification/P/application_result.tsv')
        if p['status']!='PASS' or p['full_forward_calls']!='0' or p['selection_sha']!=transfer.sha(root/'reports/selection_lock.tsv'):
            raise RuntimeError('Assessment did not consume the sealed selection')
    print('N20_TERMINAL_ARCHIVE_VERIFIED',terminal,settled,budget['attempts'])

if __name__=='__main__':
    main(sys.argv[1])
