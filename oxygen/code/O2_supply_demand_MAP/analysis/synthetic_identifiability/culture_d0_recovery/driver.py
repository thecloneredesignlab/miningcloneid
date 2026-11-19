"""Finite host passes. Never dispatch against an older culture or tumor run."""
from pathlib import Path
import argparse,fcntl,subprocess,sys
import workflow as w

def root(s,j):return s/('engineering' if j['phase'] in ('engineering','buffer') else 'checks' if j['phase']=='verification' else 'fits')/j['task']
def accepted(s,j):
    if not w.success(j):raise ValueError('Allocation did not succeed: '+j['task'])
    p=root(s,j);w.payload(p)
    if w.kv(p/'result.tsv')['pass']!='true':raise ValueError('Internal result failed')
    if (s/'runtime'/j['task']/'image_sha256.txt').read_text().strip()!=w.IMAGE_SHA:raise ValueError('Image identity failed')
    if (s/'runtime'/j['task']/'launch_exit_code.txt').read_text().strip()!='0':raise ValueError('Launcher status failed')
def numeric_check_evidence(s,j):
    # Preserve the failed end-to-end gate. Only document assertions already
    # completed before the known A X11 presentation failure; never rerun it.
    if (j['task'],j['job_id'],j['state'],j['exit_code'],j['actual_core_seconds'],j['release'])!=('check_A','21904146','FAILED','1:0','22','r05'):
        raise ValueError('Not the reviewed A presentation failure')
    p=root(s,j);z=w.kv(p/'result.tsv');failure=w.kv(p/'FAILED.tsv')
    if failure['message']!='unable to start device PNG' or z['pass']!='true' or z['prediction_recovery']!='false':
        raise ValueError('Only known failed-target A numerical evidence may be retained')
    if z['origin']!='A' or z['model_forward_calls']!='2' or z['data_sha256']!=w.kv(s/'tables/engineering_gate.tsv')['data_sha256']:
        raise ValueError('A numerical evidence identity mismatch')
    proof=w.read(p/'candidate_checks.tsv')
    if {r['source_task'] for r in proof}!={'global_A','local_A_resume_01'} or len(proof)!=2:raise ValueError('A fresh candidate proof incomplete')
    if any(abs(float(r['difference']))>1e-8 or abs(float(r['fresh'])-float(r['independent']))>1e-8 for r in proof):raise ValueError('A numerical replay failed')
    if abs(float(z['objective'])-min(float(r['fresh']) for r in proof))>1e-8:raise ValueError('A score selection mismatch')
    parameters=w.read(p/'parameter_recovery.tsv')
    if len(parameters)!=20 or sum(r['estimated']=='TRUE' for r in parameters)!=14:raise ValueError('A parameter outputs incomplete')
    if not w.read(p/'function_recovery.tsv') or not (p/'verified_candidate.rds').exists():raise ValueError('A native/function outputs incomplete')
    spec=s/'dispatch/check_A/spec.tsv';binding=w.kv(spec);rel=s/'releases/r05'
    if w.sha(spec)!=j['dispatch_sha256'] or binding['release_manifest_sha256']!=w.sha(rel/'manifest.tsv'):raise ValueError('A source binding failed')
    w.verify(rel,rel/'manifest.tsv')
    if (s/'runtime/check_A/image_sha256.txt').read_text().strip()!=w.IMAGE_SHA or (s/'runtime/check_A/launch_exit_code.txt').read_text().strip()!='1':raise ValueError('A launch evidence failed')
    paths=[x for x in p.rglob('*') if x.is_file()]+[spec,s/'logs/check_A_21904146.log',s/'runtime/check_A/image_sha256.txt',s/'runtime/check_A/launch_exit_code.txt']
    manifest=s/'tables/check_A_numerical_payload.tsv'
    w.table(manifest,[dict(path=str(x.relative_to(s)),sha256=w.sha(x)) for x in sorted(paths)],['path','sha256'])
    w.typed(s/'tables/check_A_numerical_evidence.tsv',dict(numerical_replay_evidence_pass=True,allocation_success=False,
        end_environment_guard_completed=False,prediction_target_pass=False,model_forward_calls=2,
        plot_failure='X11 PNG device',source_release='r05',data_sha256=z['data_sha256'],payload_manifest_sha256=w.sha(manifest)))
def finalize(s,rows,endpoint):
    if any(not w.closed(j) for j in rows):raise ValueError('An owned allocation remains active')
    reviewed={r['task'] for r in w.read(s/'tables/reviewed_failures.tsv')} if (s/'tables/reviewed_failures.tsv').exists() else set()
    for j in rows:
        if w.success(j):accepted(s,j)
        elif j['task'] not in reviewed:raise ValueError('Unreviewed unsuccessful allocation')
        elif j['task']=='check_A':numeric_check_evidence(s,j)
    w.verify(Path(w.kv(s/'study.tsv')['project']),s/'manifests/protected_scope_manifest.tsv')
    pos=w.position(rows);actual=sum(int(j['actual_core_seconds']) for j in rows)
    checks={j['task'][-1]:w.kv(root(s,j)/'result.tsv') for j in rows if j['task'].startswith('check_')}
    both=set(checks)=={'A','B'} and all(c['prediction_recovery']=='true' for c in checks.values())
    w.typed(s/'tables/final_accounting.tsv',dict(complete=True,endpoint=endpoint,owned_parents=len(rows),active_parents=0,
        actual_core_seconds=actual,cap_core_seconds=128*3600,prediction_recovery_verified=both,
        parameter_uniqueness_claim=False,realistic_observation_release=False,new_invivo_core_hours_released=0,
        all_check_allocations_succeeded=all(w.success(j) for j in rows if j['task'].startswith('check_'))))
    lines=['# Culture D0 conditional expected-observation recovery','',
      'Updated '+w.utc()+'. Endpoint: **'+endpoint+'**.','',
      'This is a known-six-nuisance engineering study of fourteen free mechanisms. It is not the primary twenty-unknown noisy experiment. Physical-experiment holds remain. Existing culture runs/budgets and in-vivo closure are unchanged; tumor24core-hours remain unreleased.','',
      '| Phase | Allocated actual core-hours | Authorized cap |','|---|---:|---:|']
    for name,v in pos.items():lines.append('| {} | {:.8f} | {:.0f} |'.format(name,v['actual']/3600,v['cap']/3600))
    lines+=['','Total actual: {}core-seconds ({:.8f}core-hours), including all owned allocations. Parent AllocCPUS times ElapsedRaw is charged once; nested timers are not additional fees.'.format(actual,actual/3600),'',
        '| Origin | Verified objective | Max growth / true sigma | Max karyotype TV | Max flow TV | Engineering prediction target |','|---|---:|---:|---:|---:|---|']
    for branch,z in sorted(checks.items()):lines.append('| {} | {} | {} | {} | {} | {} |'.format(branch,z['objective'],z['max_growth_error'],z['max_kary_tv'],z['max_flow_tv'],z['prediction_recovery']))
    if not checks:lines.append('| Not released | | | | | Engineering or cost hold |')
    if any(j['task']=='check_A' and not w.success(j) for j in rows):
        lines+=['','Check A remains FAILED1:0: its two fresh numerical comparisons and parameter/function tables were completed before an X11 PNG-device error. The numerical evidence is separately bound in tables/check_A_numerical_evidence.tsv; the final environment guard and end-to-end check did not pass. Derived local figures use the saved tables. No additional A model call or allocation was made for this presentation repair.']
    if not checks:lines+=['','No blind start bank was generated or executed. The registered future policy forbids truth/old-fit injection.']
    lines+=['','Truth is a feasible validated seed144 reference selected before synthetic results; it is not claimed to be the real-data global optimum. Uniform start banks contain no truth/old-fit injection. Estimates are chosen by verified score; truth comparisons occur afterward. Parameter and function errors are continuous outputs. Different vectors with matching expected observations indicate candidate ambiguity, not a finite proof of structural nonidentifiability.','',
      'No realistic noisy dataset, repeated-data batch, additional truth, ablation, profile, or new tumor task is released by this package. All numerical workers use the pinned production SIF.','',
      '| Task | Slurm parent | State | Exit | Cores | Parent seconds | Actual core-hours |','|---|---|---|---|---:|---:|---:|']
    for j in rows:lines.append('| {} | {} | {} | {} | {} | {} | {:.8f} |'.format(j['task'],j['job_id'],j['state'],j['exit_code'],j['cores'],int(j['actual_core_seconds'])//int(j['cores']),int(j['actual_core_seconds'])/3600))
    w.text(s/'REPORT.md','\n'.join(lines)+'\n')
    if endpoint.startswith('engineering_failed'):
        slots=w.read(s/'accounting/engineering_forward_slots.tsv') if (s/'accounting/engineering_forward_slots.tsv').exists() else []
        with (s/'REPORT.md').open('a') as f:
            f.write('\nEngineering did not pass. Consumed full-forward slots: {}. No blind fit was released. '
              'Failed costs and attempt history remain attached; another preparation retry requires explicit authorization. '
              'The128core-hour ceiling, first8core-hour window and64forward-call ceiling have not increased.\n'.format(len(slots)))
    w.text(s/'STATUS.md','# Culture D0 pilot closed\n\n'+endpoint+';0active parents. Actual{:.8f}/128core-hours. Existing culture and tumor scopes protected.\n'.format(actual/3600))
    files=[]
    for d in ('inputs','manifests','releases','dispatch','truth','design','data','engineering','fits','checks','tables','accounting','logs'):
        files.extend(p for p in (s/d).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    for pattern in ('image_sha256.txt','launch_exit_code.txt','model_source_manifest.tsv'):
        files.extend((s/'runtime').rglob(pattern))
    files.extend(p for p in (s/'reports').rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    for name in ('authorizations','closure_history'):
        if (s/name).exists():files.extend(p for p in (s/name).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    files.extend([s/'study.tsv',s/'authorization.tsv',s/'APPROVED_PLAN.md',s/'REPORT.md',s/'STATUS.md'])
    listing=[dict(path=str(p.relative_to(s)),sha256=w.sha(p)) for p in sorted(set(files)) if p.name not in ('delivery_manifest.tsv','delivery_verified.tsv')]
    w.table(s/'tables/delivery_manifest.tsv',listing,['path','sha256']);w.verify(s,s/'tables/delivery_manifest.tsv')
    w.typed(s/'tables/delivery_verified.tsv',dict(pass_gate=True,files=len(listing),delivery_manifest_sha256=w.sha(s/'tables/delivery_manifest.tsv'),owned_parents=len(rows),all_parents_closed=True,authority='HPC full scientific payload and source verification'))
    print('D0_TERMINAL',endpoint,'actual_core_seconds',actual,'files',len(listing),flush=True)

def advance(s):
    if (s/'tables/final_accounting.tsv').exists():print('D0_TERMINAL already closed; no state mutation');return
    w.verify(Path(w.kv(s/'study.tsv')['project']),s/'manifests/protected_scope_manifest.tsv')
    rows=w.reconcile(s);by={j['task']:j for j in rows}
    reviewed={r['task'] for r in w.read(s/'tables/reviewed_failures.tsv')} if (s/'tables/reviewed_failures.tsv').exists() else set()
    failed=[j for j in rows if w.closed(j) and not w.success(j) and j['task'] not in reviewed]
    if failed:raise ValueError('CONTROLLER_HOLD: review failed allocation, preserve costs and64-slot counter: '+','.join(j['task'] for j in failed))
    engineering=[]
    for task in ('prepare','parallel','engineering_check'):
        matches=[j for j in rows if j['mode']==task]
        j=matches[-1] if matches else None
        if j is None:w.submit(s,task);return
        if not w.closed(j):return
        accepted(s,j)
        engineering.append(dict(mode=task,task=j['task']))
        w.table(s/'tables/engineering_sources.tsv',engineering,['mode','task'])
    gate=w.kv(s/'tables/engineering_gate.tsv')
    if gate['pass']!='true':raise ValueError('Engineering gate failed')
    evidence=['tables/engineering_gate.tsv','tables/search_controls.tsv','truth/manifest.tsv','data/expected_observations.rds','design/blind_bundle.rds']
    evidence+=['engineering/'+r['task']+'/payload_manifest.tsv' for r in engineering]
    w.table(s/'tables/engineering_evidence.tsv',[dict(path=p,sha256=w.sha(s/p)) for p in evidence],['path','sha256'])
    if gate['cost_fits']!='true':finalize(s,rows,'engineering_pass_cost_hold');return
    for branch in ('A','B'):
        task='global_'+branch;j=by.get(task)
        if j is None:w.submit(s,task);continue
        if not w.closed(j):continue
        accepted(s,j)
        task='local_'+branch;j=by.get(task)
        if j is not None and w.closed(j) and not w.success(j):
            if task!='local_A' or task not in reviewed:raise ValueError('Unreviewed local hold')
            # A reviewed failed parent never silently creates another attempt.
            # Continue B's already authorized work while waiting for A consent.
            if not (s/'authorizations/local_A_resume_01/approval.tsv').exists():continue
            task='local_A_resume_01';j=by.get(task)
        if j is None:w.submit(s,task,mode='local_'+branch);continue
        if not w.closed(j):continue
        accepted(s,j)
        if task=='local_A_resume_01':
            if int(w.kv(root(s,j)/'result.tsv')['evaluations'])+2>11183:raise ValueError('Cumulative local callback cap')
        aliases={r['origin']:r['task'] for r in w.read(s/'tables/local_sources.tsv')} if (s/'tables/local_sources.tsv').exists() else {}
        aliases[branch]=task;w.table(s/'tables/local_sources.tsv',[dict(origin=b,task=t) for b,t in sorted(aliases.items())],['origin','task'])
        task='check_'+branch;j=by.get(task)
        if j is None:w.submit(s,task);continue
        if w.closed(j):
            if w.success(j):accepted(s,j)
            elif task=='check_A' and task in reviewed:numeric_check_evidence(s,j)
            else:raise ValueError('Independent check allocation failed')
    rows=w.reconcile(s);by={j['task']:j for j in rows}
    ready_A='check_A' in by and w.closed(by['check_A']) and (w.success(by['check_A']) or 'check_A' in reviewed)
    if ready_A and 'check_B' in by and w.success(by['check_B']):
        if w.success(by['check_A']):accepted(s,by['check_A'])
        else:numeric_check_evidence(s,by['check_A'])
        accepted(s,by['check_B'])
        both=all(w.kv(root(s,by['check_'+b])/'result.tsv')['prediction_recovery']=='true' for b in ('A','B'))
        endpoint='conditional_prediction_recovery_pass' if both else 'conditional_recovery_optimizer_unresolved'
        if not w.success(by['check_A']):endpoint+='_check_A_artifact_failure'
        finalize(s,rows,endpoint)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('study',type=Path);s=a.parse_args().study.resolve()
    with (s/'state/DRIVER_LOCK').open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB);advance(s)
