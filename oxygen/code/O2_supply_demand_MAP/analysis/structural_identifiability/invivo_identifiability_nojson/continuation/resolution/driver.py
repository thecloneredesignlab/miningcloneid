"""Bounded single-owner passes for the approved resolution01 protocol."""
from pathlib import Path
import argparse,re,uuid,subprocess,fcntl
import workflow as w
import policy as p


def success(j):return j['state']=='COMPLETED' and j['exit_code']=='0:0'
def root(s,j):return s/('fits' if j['mode'] in ('local','global') else 'tests')/j['task']

def assert_reviewed(s,rows):
    file=s/'tables/reviewed_failures.tsv'
    reviewed={r['task'] for r in w.read_table(file)} if file.exists() else set()
    bad=[j['task'] for j in rows if w.closed(j) and not success(j) and j['task'] not in reviewed]
    if bad:raise ValueError('Unreviewed closed allocation holds the controller: '+','.join(bad))

def defer_old_representative(s,a,ledger,by,checks_by_fit,old_best):
    """Prefer current evidence; a later checked lower point supersedes old C02."""
    j=by[a['task']];seen=False
    for b in ledger:
        if b['task']==a['task']:seen=True;continue
        if not seen or b['slot_class']!=a['slot_class'] or b['slot_key']!=a['slot_key']:continue
        later=by[b['task']]
        if (later['root_id'],later['policy'],later['density_index'])!=(j['root_id'],j['policy'],j['density_index']):continue
        if not w.closed(later):return True
        if not success(later):continue
        proofs=checks_by_fit.get(later['task'],[])
        if not proofs or any(not w.closed(c) for c in proofs):return True
        for check in proofs:
            if not success(check):continue
            r=w.kv(root(s,check)/'result.tsv')
            if r.get('pass')=='true' and r.get('scope')=='fixed_center_all_four_scales_axis_coupled_observation' and r['complete']=='true' and float(r['representative_gap'])<=1e-4 and float(r['directional_improvement'])<=.002 and float(r['best_objective'])<=old_best:
                return True
    return False

def submit(s,task,**options):
    try:return w.submit(s,task,**options)
    except p.PhaseBudgetStop as error:
        resource_hold(s,w.snapshot(s),options.get('policy','all'),str(error),options['phase'],task)


def validation(s,rows):
    path=s/'tables/resolution_validation_gate.tsv'
    if path.exists():w.gate(s,path.name);return True
    cases=[j for j in rows if j['mode']=='resolution_validation']
    if len(cases)!=4 or not all(success(j) for j in cases):return False
    evidence=['releases/res01/manifest.tsv','releases/cont06/manifest.tsv']
    for file in ('adapter.R','helpers.R','control.R','machine.R','optimize.R','global.R','observation.R'):
        if w.sha(s/'releases/res01'/file)!=w.sha(s/'releases/cont06'/file):raise ValueError('Numerical method changed')
    for j in cases:
        a=root(s,j);w.verify_payload(a);z=w.kv(a/'result.tsv')
        if z['pass']!='true' or z['complete']!='true' or float(z['difference'])!=0 or float(z['fixed_residual'])!=0:
            raise ValueError('Incomplete real-map full-check validation')
        if int(z['free_coordinates'])!=(19 if int(j['density_index']) else 20):raise ValueError('Wrong map')
        evidence.append(str((a/'payload_manifest.tsv').relative_to(s)))
    ep=s/'tables/resolution_validation_evidence.tsv'
    w.table(ep,[dict(path=f,sha256=w.sha(s/f)) for f in evidence],['path','sha256'])
    w.typed(path,{'pass':True,'release':'res01','input_manifest_sha256':w.sha(s/'manifests/input_manifest.tsv'),
        'evidence_manifest':str(ep.relative_to(s)),'real_map_cases':4,'science_release':'cont06'})
    print('VALIDATED resolution01 free20/fixed19 both policies; cont06 numerical source bridge unchanged.')
    return True


def checks(s,rows):
    ledger=p.assignments(s);by={j['task']:j for j in rows}
    checks_by_fit={}
    for c in ledger:
        if c['slot_class']!='check':continue
        spec=w.kv(s/'dispatch'/c['task']/'spec.tsv')
        checks_by_fit.setdefault(spec['checked_fit'],[]).append(by[c['task']])
    for a in ledger:
        if a['slot_class'] not in p.LIMITS and not (a['slot_class']=='reference' and by[a['task']]['kind']=='baseline'):continue
        j=by[a['task']]
        if not success(j):continue
        existing=checks_by_fit.get(j['task'],[])
        if existing:
            latest=existing[-1]
            if not success(latest):continue
            z=w.kv(root(s,latest)/'result.tsv')
            # One additional representative check is evidence completion,
            # never another disguised 20000-callback fitting component.
            if len(existing)>=2 or float(z['representative_gap'])<=1e-4:continue
            if defer_old_representative(s,a,ledger,by,checks_by_fit,float(z['best_objective'])):continue
            fit_result=w.kv(root(s,j)/'result.tsv')
            count=sum(x['slot_class']==a['slot_class'] and x['slot_key']==a['slot_key'] for x in ledger)
            gain=fit_result.get('window_gain','0')
            if gain!='NA' and float(gain)>=.001 and count<p.LIMITS.get(a['slot_class'],0):continue
            if j['mode']=='global':continue
            source=latest['task']
        else:source=j['task']
        task='r_resolution01_'+j['task']+'_c'+str(len(existing)+1).zfill(2)
        submit(s,task,mode='resolution_check',role='resolution_check' if a['phase']=='baseline' else 'density_check',
          kind='baseline' if a['phase']=='baseline' else 'profile',release='res01',phase=a['phase'],slot_class='check',
          slot_key=j['task']+':full',checked_fit=j['task'],policy=j['policy'],root_id=j['root_id'],
          root_family=j['root_family'],selector='task',source_task=source,density_index=int(j['density_index']),
          cores=1,wall_seconds=900,margin_seconds=120,max_calls=500)


def capture(s,rows):
    state=w.State(s)
    with state.locked():
        current=state.load()
        if [(j['task'],j['state']) for j in current]!=[(j['task'],j['state']) for j in rows]:
            raise ValueError('Controller changed before immutable audit capture')
        generation=(s/'state/CURRENT').read_text().strip();identity='c_'+uuid.uuid4().hex
        directory=s/'collections'/identity;directory.mkdir()
        (directory/'jobs.tsv').write_bytes((s/'state'/generation/'jobs.tsv').read_bytes())
        (directory/'raw_accounting.txt').write_bytes((s/'tables/raw_accounting.txt').read_bytes())
        w.table(directory/'origin.tsv',[dict(key='state_generation',value=generation),dict(key='captured_utc',value=w.utc())],['key','value'])
        w.table(directory/'manifest.tsv',[dict(path=f.name,sha256=w.sha(f)) for f in sorted(directory.iterdir())],['path','sha256'])
    return identity,w.sha(directory/'manifest.tsv')


def audit(s,rows):
    audits=[j for j in rows if j['task'].startswith('a_resolution01_') and j['mode']=='metadata_audit']
    if any(not w.closed(j) for j in audits):return False
    # Baseline fitting also waits for outstanding full checks. Coalesce their
    # nearly simultaneous completions into one useful archived-output audit.
    if any(j['mode']=='resolution_check' and not w.closed(j) for j in rows):return False
    checks_ready=set()
    for j in rows:
        if j['mode']=='resolution_check' and success(j):checks_ready.add(j['task'])
    bridge=w.read_table(s/'policies/resolution01/bridge.tsv')
    if not audits:
        # One coherent initial audit covers all inherited endpoints and all
        # four validation maps; no per-root redundant audit wave.
        checked={w.kv(s/'dispatch'/t/'spec.tsv').get('checked_fit') for t in checks_ready}
        if not {j['task'] for j in bridge}<=checked:return False
    previous=[j for j in audits if success(j)]
    if previous:
        spec=w.kv(s/'dispatch'/previous[-1]['task']/'spec.tsv')
        old=w.read_table(s/'collections'/spec['metadata_collection']/'jobs.tsv')
        old_checks={j['task'] for j in old if j['mode']=='resolution_check' and success(j)}
        if checks_ready<=old_checks and (s/'tables/resolution_native_identity_gate.tsv').exists():return True
    identity,digest=capture(s,rows)
    phase='profile' if (s/'tables/profile_epoch.tsv').exists() else 'baseline'
    task='a_resolution01_w'+str(len(audits)+1).zfill(4)
    options=dict(mode='metadata_audit',role='archived_output_validation',kind='engineering',
        release=w.kv(s/'tables/resolution_runtime.tsv')['audit_release'],phase=phase,slot_class='audit',slot_key=task,
        root_id='metadata_reference',root_family='reference',selector='reference',
        cores=1,wall_seconds=3600,margin_seconds=300,max_calls=1,
        metadata_collection=identity,metadata_collection_sha256=digest)
    try:w.submit(s,task,**options)
    except p.PhaseBudgetStop as error:
        resource_hold(s,w.snapshot(s),'all',str(error),phase,task)
        options['slot_class']='closure'
        w.submit(s,task,**options)
    return False


def audited(s,j):
    f=s/'tables/record_audits'/(j['task']+'_summary.tsv')
    return f.exists() and w.kv(f)['pass']=='true' and w.kv(f)['payload_sha256']==w.sha(root(s,j)/'payload_manifest.tsv')


def baseline(s,rows):
    # A prior audit can still describe the previous strategy while a freshly
    # closed fit's complete check is running. Wait for that evidence and the
    # ensuing audit before using any stale strategy to issue another fit.
    if any(j['mode']=='resolution_check' and not w.closed(j) for j in rows):return
    ledger=p.assignments(s);by={j['task']:j for j in rows}
    policies={r['policy']:r for r in w.read_table(s/'tables/baseline_policy_gates.tsv')}
    d={r['task']:r for r in w.read_table(s/'tables/resolution_local_decisions.tsv')}
    families={}
    for row in w.read_table(s/'tables/baseline_chain_summary.tsv'):
        before=families.setdefault(row['root_id'],(row['policy'],row['family']))
        if before!=(row['policy'],row['family']):raise ValueError('One root cannot supply two independent lineage families')
    replay_path=s/'tables/resolution_target_replays.tsv'
    if replay_path.exists():
        pending=False
        for r in w.read_table(replay_path):
            key='target:'+r['policy']+':'+r['record_sha256']
            if any(a['slot_class']=='reference' and a['slot_key']==key for a in ledger):
                pending=True;continue
            task='q_resolution01_'+r['policy']+'_'+r['record_sha256'][:12]
            submit(s,task,mode='reference_reconstruct',role='target_policy_reference',kind='baseline',
                release='res03',phase='baseline',slot_class='reference',slot_key=key,
                reference_purpose='target_replay',policy=r['policy'],root_id=r['root_id'],root_family='informed',
                selector='selection',selection_name=r['selection'],cores=1,wall_seconds=900,margin_seconds=120,max_calls=1)
            pending=True
        if pending:return
    schedule=[];latest={}
    for a in ledger:
        if a['slot_class']=='existing_local':latest[a['slot_key']]=a
    for key,a in latest.items():
        j=by[a['task']];decision=d.get(j['task']);count=sum(b['slot_class']=='existing_local' and b['slot_key']==key for b in ledger)
        if not success(j) or not decision:continue
        selected=decision['decision']
        if selected=='continue' and count>=3:selected='window_cap_stop'
        schedule.append(dict(task=j['task'],root_id=key,used=count,maximum=3,native_termination=w.kv(root(s,j)/'result.tsv')['termination'],decision=selected,policy_version=p.VERSION))
        if selected!='continue' or policies[j['policy']]['pass'].lower()=='true':continue
        match=re.fullmatch('(.*)_w([0-9]+)',j['task']);task=match[1]+'_w'+str(int(match[2])+1).zfill(4)
        if task in by:continue
        submit(s,task,phase='baseline',slot_class='existing_local',slot_key=key,decision_fit=j['task'],
            release='cont06',policy=j['policy'],root_id=key,root_family=j['root_family'],
            selector=decision['selector'],source_task=decision['continuation_source'])
    if schedule:w.table(s/'tables/resolution_schedule.tsv',schedule,['task','root_id','used','maximum','native_termination','decision','policy_version'])
    # Each successful global gets exactly its own one polish, retaining an
    # improving check point only inside that same root/branch.
    for a in ledger:
        if a['slot_class']!='baseline_global':continue
        j=by[a['task']]
        if not success(j) or not audited(s,j) or policies[j['policy']]['pass'].lower()=='true':continue
        if any(b['slot_class']=='baseline_polish' and b['slot_key']==a['slot_key'] for b in ledger):continue
        dec=d.get(j['task'])
        if not dec:continue
        source=dec['check_task']
        submit(s,'b_'+j['task']+'_polish_w0001',phase='baseline',slot_class='baseline_polish',slot_key=a['slot_key'],
            polish_parent=j['task'],release='cont06',policy=j['policy'],root_id=j['root_id'],root_family='uniform',
            selector='task',source_task=source)
    strategies=w.read_table(s/'tables/strategy_review.tsv')
    batches=w.read_table(s/'tables/resolution_batches.tsv') if (s/'tables/resolution_batches.tsv').exists() else []
    for decision in strategies:
        policy=decision['policy'];mode=decision['decision']
        if mode in ('initial_global_batch','next_global_batch'):
            batch=int(decision['batch'])
            if any(r['policy']==policy and int(r['batch'])==batch for r in batches):continue
            # Reserve the complete two-branch plan including each polish and
            # both complete-check allocations before dispatching either branch.
            current=w.State(s).load();position=p.positions(current,p.assignments(s))['baseline']
            if position['committed']+2*(57600+8100+900)>position['cap']:
                current=w.snapshot(s);position=p.positions(current,p.assignments(s))['baseline']
                if position['committed']+2*(57600+8100+900)>position['cap']:
                    resource_hold(s,current,policy,'Insufficient complete two-branch batch commitment','baseline','batch_'+str(batch));continue
            chains=w.read_table(s/'tables/baseline_chain_summary.tsv')
            own=min((r for r in chains if r['policy']==policy and r['family']=='uniform'),key=lambda r:float(r['verified_best']))
            entry=dict(policy=policy,batch=str(batch),root_id=own['root_id'],before_uniform=own['verified_best'],
                before_feasible=policies[policy]['known_feasible_bound'],reason=decision['reason'],recorded_utc=w.utc())
            batches.append(entry);w.table(s/'tables/resolution_batches.tsv',batches,list(entry))
            for index,branch in enumerate(('hybrid','uniform')):
                task='g_'+policy+'_c'+str(batch).zfill(4)+'_'+branch+'_w0001'
                opts=dict(mode='global',role='baseline_global',phase='baseline',slot_class='baseline_global',
                    slot_key='{}:{}:{}'.format(policy,batch,branch),batch=batch,branch=branch,release='cont06',
                    policy=policy,root_family='uniform',init_mode=branch,seed=61000+100*batch+(10 if policy=='data_only' else 0)+index,
                    cores=8,wall_seconds=7200,margin_seconds=600,max_calls=200000)
                if branch=='hybrid':opts.update(root_id=own['root_id'],selector='selection',selection_name='root_'+own['root_id']+'_selected')
                else:opts.update(root_id='full_uniform_'+policy+'_c'+str(batch).zfill(4),selector='reference')
                submit(s,task,**opts)
        elif mode=='finite_baseline_repair':
            repairs=[a for a in ledger if a['slot_class']=='baseline_repair' and a['slot_key']==policy]
            if len(repairs)>=2 or any(not w.closed(by[a['task']]) for a in repairs):continue
            chains=[r for r in w.read_table(s/'tables/baseline_chain_summary.tsv') if r['policy']==policy]
            floor=float(policies[policy]['known_feasible_bound'])
            family=None
            for f in ('archive','uniform'):
                if not any(r['family']==f and r['representative_pass'].lower()=='true' and float(r['verified_best'])<=floor+.002 for r in chains):family=f;break
            choices=[r for r in chains if family is None or r['family']==family]
            own=min(choices,key=lambda r:float(r['verified_best']))
            task='b_resolution01_'+policy+'_repair_w'+str(len(repairs)+1).zfill(4)
            submit(s,task,phase='baseline',slot_class='baseline_repair',slot_key=policy,release='cont06',
                role='baseline_repair',policy=policy,root_id=own['root_id'],root_family=own['family'],
                selector='selection',selection_name='root_'+own['root_id']+'_selected',
                repair_reason='independent_recovery' if family else 'representative_check')


def record_budget_stop(s,policy,reason,phase='baseline',task=''):
    path=s/'tables/resolution_budget_stops.tsv';old=w.read_table(path) if path.exists() else []
    if not any(r['policy']==policy and r.get('phase','baseline')==phase for r in old):
        old.append(dict(policy=policy,phase=phase,task=task,decision='phase_budget_stop',reason=reason,recorded_utc=w.utc()))
        w.table(path,old,['policy','phase','task','decision','reason','recorded_utc'])

def resource_hold(s,rows,policy,reason,phase,task):
    # Full commitments can temporarily block a new allocation. Settle active
    # parents before calling that a finite resource endpoint; no cap is relaxed.
    assigned={a['task'] for a in p.assignments(s) if a['phase']==phase}
    if any(j['task'] in assigned and not w.closed(j) for j in rows):
        print('PHASE_BUDGET_WAIT',phase,task,reason,flush=True)
        return 'wait'
    record_budget_stop(s,policy,reason,phase,task)
    print('PHASE_BUDGET_STOP',phase,task,reason,flush=True)
    return 'stop'

def baseline_terminal(s,rows):
    if any(not w.closed(j) for j in rows):return False
    strategies={r['policy']:r['decision'] for r in w.read_table(s/'tables/strategy_review.tsv')}
    path=s/'tables/resolution_budget_stops.tsv'
    for stop in w.read_table(path) if path.exists() else []:
        if stop.get('phase','baseline')!='baseline':continue
        for policy in strategies:
            if stop['policy'] in (policy,'all') and strategies[policy]!='accepted':strategies[policy]='phase_budget_stop'
    if set(strategies)!={'production','data_only'}:raise ValueError('Missing score-policy outcome')
    if any(v not in ('accepted','resolution_unresolved','phase_budget_stop') for v in strategies.values()):return False
    both=all(v=='accepted' for v in strategies.values())
    w.typed(s/'tables/resolution_baseline_terminal.tsv',dict(complete=True,both_baselines_accepted=both,
        production=strategies['production'],data_only=strategies['data_only'],
        next_stage='conditional_profile_required' if both else 'profiles_not_released',model_evaluations=0))
    print('BASELINE_TERMINAL',strategies,'PROFILE_REQUIRED' if both else 'UNRESOLVED_NO_PROFILE',flush=True)
    return True


def advance(s):
    rows=w.snapshot(s)
    assert_reviewed(s,rows)
    if not validation(s,rows):print('Await four real full-check validation cases.');return
    checks(s,rows);rows=w.snapshot(s)
    if not audit(s,rows):print('Await complete endpoint checks/allocated archived-output validation.');return
    print('REVIEW_READY')
    if not (s/'tables/profile_epoch.tsv').exists():
        baseline(s,rows)
        baseline_terminal(s,w.snapshot(s))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('study',type=Path)
    study=parser.parse_args().study.resolve()
    with (study/'state/resolution_driver.lock').open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        advance(study)
