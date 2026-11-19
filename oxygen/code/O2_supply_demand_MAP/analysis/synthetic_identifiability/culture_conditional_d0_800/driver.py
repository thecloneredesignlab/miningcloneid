"""One finite orchestration pass; no scheduler concurrency throttle."""
from pathlib import Path
import argparse,fcntl,shutil,subprocess
import workflow as w

def root(s,j):
    if j['mode'] in ('prepare','replay','replay_final','parallel','engineering_check','local_api_check'):return s/'engineering'/j['task']
    if j['mode']=='report':return s/'reports/science'
    if j['mode']=='report_resume01':return s/'reports/science_resume01'
    return s/('checks' if j['mode'].startswith(('globalcheck_','check_')) or j['mode']=='positive_check' else 'fits')/j['task']
def accepted(s,j):
    if not w.success(j):raise ValueError('Allocation not successful: '+j['task'])
    p=root(s,j);w.payload(p)
    if w.kv(p/'result.tsv')['pass']!='true':raise ValueError('Internal result failed: '+j['task'])
    if (s/'runtime'/j['task']/'image_sha256.txt').read_text().strip()!=w.IMAGE_SHA:raise ValueError('Image mismatch')
    if (s/'runtime'/j['task']/'launch_exit_code.txt').read_text().strip()!='0':raise ValueError('Launcher failed')
def finish(s,rows,endpoint):
    if any(not w.closed(j) for j in rows):return
    current_report='report_resume01' if (s/'provenance_parent').exists() else 'report'
    report=next((j for j in rows if j['mode']==current_report),None)
    if report is None:w.submit(s,current_report);return
    if not w.closed(report):return
    if w.success(report):accepted(s,report)
    else:endpoint='report_failed_after_'+endpoint
    w.verify_inputs(s);pos=w.position(rows)
    ids=','.join(j['job_id'] for j in rows)
    resource=subprocess.run([str(w.SLURM/'sacct'),'-n','-P','-j',ids,'--format=JobIDRaw,State%30,ExitCode,AllocCPUS,ElapsedRaw,MaxRSS,MaxVMSize,TotalCPU,ReqMem'],check=True,universal_newlines=True,stdout=subprocess.PIPE)
    fields=['job_id','state','exit_code','allocated_cpus','elapsed_seconds','max_rss','max_vm_size','total_cpu','requested_memory']
    w.table(s/'reports/slurm_resource_accounting.tsv',[dict(zip(fields,line.split('|'))) for line in resource.stdout.splitlines() if line.strip()],fields)
    # Preserve private native/source evidence without shipping binary caches.
    names={'model_source_manifest.tsv','runtime_packages.tsv','loaded_native_libraries.tsv','loaded_namespaces.tsv','image_sha256.txt','launch_exit_code.txt'}
    for p in (s/'runtime').rglob('*'):
        if p.is_file() and p.name in names:
            target=s/'reports/runtime_provenance'/p.relative_to(s/'runtime');target.parent.mkdir(parents=True,exist_ok=True)
            if target.exists():
                if w.sha(p)!=w.sha(target):raise ValueError('Prior runtime provenance differs; preserve original')
            else:shutil.copy2(str(p),str(target))
    actual=sum(int(j['actual_core_seconds']) for j in rows)
    w.typed(s/'tables/final_accounting.tsv',dict(endpoint=endpoint,owned_parents=len(rows),active_parents=0,actual_core_seconds=actual,authorized_core_seconds=800*3600,parameter_uniqueness_claim=False,physical_sampling_release=False,noisy_datasets=0,profiles=0,tumor_fits=0))
    lines=['# Conditional D0 recovery execution','', 'Endpoint: **'+endpoint+'**. Updated '+w.utc()+'.','',
      'One prespecified seed144 generating truth; fourteen free mechanisms; six generating nuisance fixed. This is production-count-conditioned deterministic expected-observation recovery, not finite noisy sampling, biological model validation, or a structural proof.','',
      'Physical growth/stopping, flow preprocessing/event denominators, and live/dead mapping remain held. No profiles, noisy replication, or tumor fitting was released. The closed earlier D0 startup failures remain in their own immutable study; the current independent ceiling is800core-hours.','',
      '| Phase | Actual allocated core-hours | Ceiling |','|---|---:|---:|']
    for name,z in pos.items():lines.append('| {} | {:.8f} | {:.0f} |'.format(name,z['actual']/3600,z['cap']/3600))
    lines+=['','Actual total: {:.8f}core-hours. Every failed allocation is counted once using parent AllocCPUS times ElapsedRaw.'.format(actual/3600),'',
      '| Task | Parent | State | Exit | Actual core-hours |','|---|---|---|---|---:|']
    for j in rows:lines.append('| {} | {} | {} | {} | {:.8f} |'.format(j['task'],j['job_id'],j['state'],j['exit_code'],int(j['actual_core_seconds'])/3600))
    science='science_resume01' if current_report=='report_resume01' else 'science'
    q=s/'reports'/science/'prediction_recovery.tsv'
    if q.exists():
        lines+=['','| Origin | Independent | Objective gap | Growth RMS / sigma | Karyotype TV | Flow-mass TV | Prediction target | Local stability |','|---|---|---:|---:|---:|---:|---|---|']
        for z in w.read(q):lines.append('| {origin} | {independent} | {delta} | {growth_rms} | {kary_tv} | {flow_tv} | {prediction_pass} | {local_stability} |'.format(**z))
    else:lines+=['','No verified blind-recovery result is available.']
    lines+=['','Parameter errors and one-percent transformed-domain resolution are in reports/science/parameter_recovery.tsv where available. Continuous errors, bounds, endpoint inactivity and function differences are kept separate from optimizer sufficiency. A finite optimizer failure is labelled unresolved, not nonidentifiability. Truth-initialized positive control never counts as independent. Estimates are selected by verified objective, never distance to truth.']
    if (s/'tables/STOP.tsv').exists():lines+=['','Stop receipt:']+['- '+k+': '+v for k,v in w.kv(s/'tables/STOP.tsv').items()]
    w.text(s/'REPORT.md','\n'.join(lines)+'\n');w.text(s/'STATUS.md','# Conditional D0 closed\n\n'+endpoint+'; all owned parents terminal.\n')
    listing=[]
    for p in sorted(s.rglob('*')):
        rel=p.relative_to(s)
        if not p.is_file() or 'runtime' in rel.parts or 'state' in rel.parts or '__pycache__' in rel.parts or rel.parts[0]=='collection' or p.name in ('delivery_manifest.tsv','delivery_verified.tsv','COMPLETE'):continue
        if p.suffix.lower() in ('.json','.jsonl','.ndjson','.js','.java','.jar','.class'):raise ValueError('Forbidden artifact: '+str(rel))
        listing.append(dict(path=str(rel),sha256=w.sha(p),bytes=p.stat().st_size))
    w.table(s/'tables/delivery_manifest.tsv',listing,['path','sha256','bytes']);w.verify(s,s/'tables/delivery_manifest.tsv')
    w.typed(s/'tables/delivery_verified.tsv',dict(pass_gate=True,files=len(listing),bytes=sum(x['bytes'] for x in listing),manifest_sha256=w.sha(s/'tables/delivery_manifest.tsv'),all_owned_parents_terminal=True))
    w.text(s/'COMPLETE',w.sha(s/'tables/delivery_manifest.tsv')+'\n');print('D0_TERMINAL',endpoint,actual,flush=True)
def advance(s):
    if (s/'COMPLETE').exists():print('Study already closed');return
    rows=w.reconcile(s);by={j['mode']:j for j in rows}
    if (s/'tables/STOP.tsv').exists():finish(s,rows,w.kv(s/'tables/STOP.tsv')['reason']);return
    unreviewed=[j for j in rows if w.closed(j) and not w.success(j)]
    reviewed={r['task'] for r in w.read(s/'tables/reviewed_failures.tsv')} if (s/'tables/reviewed_failures.tsv').exists() else set()
    if any(j['task'] not in reviewed for j in unreviewed):raise ValueError('HOLD: review failed allocations; no automatic retries')
    for mode in ('prepare','replay','parallel','engineering_check','replay_final','local_api_check'):
        j=by.get(mode)
        if j is None:w.submit(s,mode);return
        if not w.closed(j):return
        accepted(s,j)
    evidence=['tables/engineering_gate.tsv','tables/search_controls.tsv','truth/manifest.tsv','data/expected_observations.rds','design/blind_bundle.rds']
    evidence+=['engineering/'+by[m]['task']+'/payload_manifest.tsv' for m in ('prepare','replay','parallel','engineering_check','replay_final')]
    if not (s/'tables/engineering_evidence_r02.tsv').exists():w.table(s/'tables/engineering_evidence_r02.tsv',[dict(path=p,sha256=w.sha(pth)) for p in evidence for pth in [s/p]],['path','sha256'])
    gate=w.kv(s/'tables/engineering_gate.tsv')
    if gate['pass']!='true' or gate['cost_fits']!='true':w.typed(s/'tables/STOP.tsv',dict(reason='engineering_or_cost_hold',utc=w.utc()));return
    for mode in ('positive_local','positive_check'):
        j=by.get(mode)
        if j is None:w.submit(s,mode);return
        if not w.closed(j):return
        accepted(s,j)
    positive=w.kv(s/'checks'/by['positive_check']['task']/'result.tsv')
    if positive['prediction_recovery']!='true' or positive['local_directional_stability']!='true':w.typed(s/'tables/STOP.tsv',dict(reason='positive_control_hold',utc=w.utc()));return
    def chain(origin):
        for step in ('screen','global','globalcheck','local','check'):
            mode=step+'_'+origin;j=by.get(mode)
            if j is None:w.submit(s,mode);return False
            if not w.closed(j):return False
            accepted(s,j)
        return True
    if not chain('A'):return
    if not (s/'tables/benchmark_gate.tsv').exists():
        used=sum(int(j['actual_core_seconds']) for j in rows if j['mode'].endswith('_A'))
        w.typed(s/'tables/benchmark_gate.tsv',dict(pass_gate=used<=105*3600,**{'pass':used<=105*3600},measured_complete_chain_core_seconds=used,remaining_three_max_commitment=315*3600))
    # All remaining primary origins are eligible together. No running-chain cap.
    complete=[chain(origin) for origin in 'BCD']
    if not all(complete):return
    checks=[w.kv(s/'checks'/by['check_'+o]['task']/'result.tsv') for o in 'ABCD']
    good=lambda z:z['prediction_recovery']=='true' and z['local_directional_stability']=='true'
    if sum(good(z) for z in checks)<2:
        # One reserve round, at most two complete independent origin packages.
        if not all([chain(origin) for origin in 'EF']):return
        checks += [w.kv(s/'checks'/by['check_'+o]['task']/'result.tsv') for o in 'EF']
    w.typed(s/'tables/fitting_complete.tsv',dict(pass_gate=sum(good(z) for z in checks)>=2,independent_accepted=sum(good(z) for z in checks),endpoint='conditional_recovery_pass' if sum(good(z) for z in checks)>=2 else 'optimizer_unresolved',utc=w.utc()))
    finish(s,w.reconcile(s),w.kv(s/'tables/fitting_complete.tsv')['endpoint'])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);s=p.parse_args().study.resolve()
    with (s/'state/DRIVER_LOCK').open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB);advance(s)
