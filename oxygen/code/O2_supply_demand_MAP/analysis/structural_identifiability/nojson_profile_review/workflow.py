#!/usr/bin/env python3
"""Bounded seven-point culture study; transactional TSV control, no throttle."""
import argparse
from pathlib import Path
import subprocess
import time
import uuid
import ledger as L
from ledger import (State, read_table, write_table, read_kv, write_kv, atomic_text,
                    sha, now, is_closed, job_ok, get_job, accepted, committed)

def task(study, tid, mode, phase, point='base', origin='', family='', source=None,
         alternative=False, fallback=None, fixture=False):
    cases={r['point']:r for r in read_table(study/'cases.tsv')}
    if point!='base' and point not in cases: raise ValueError('Unknown fixed point')
    case=cases.get(point)
    seed=next((int(r['seed']) for r in read_table(study/'seeds.tsv') if r['task_id']==tid),None)
    if seed is None: raise ValueError('Undeclared task seed '+tid)
    return dict(task_id=tid,mode=mode,phase=phase,origin=origin or tid,family=family,
                case='fixed_target' if case else 'unrestricted',point=point,
                target=case['target'] if case else 'none',value=float(case['value']) if case else None,
                source_fit=source,fallback_source=fallback,seed=seed,alternative=alternative,
                min_generations=1 if fixture else 300,max_generations=1 if fixture else 1000)

def add(state,data,t):
    L.add_task(state,data,t)
    return t

def source(study,project,tid):
    return str((study/'fits'/tid/'accepted_fit.rds').relative_to(project))

def counts(data,phase=None):
    jobs=[j for j in data['jobs'] if phase is None or j['phase']==phase]
    reservations=[r for r in data['reservations'] if phase is None or r['phase']==phase]
    return len({j['task'] for j in jobs if j['role'] in {'screen','global'}}),sum(j['role']=='local' for j in jobs)+len(reservations)

def submit(state,data,project,tid,mode,role,phase):
    if get_job(data,tid,mode): return get_job(data,tid,mode)
    ng,nl=counts(data)
    rg,rl=counts(data,phase)
    if role=='screen' and (ng>=36 or (phase.startswith('reserve_') and rg>=4)):
        raise ValueError('Global-origin attempt cap reached')
    if role=='local':
        data['reservations'][:]=[r for r in data['reservations'] if r['task']!=tid]
        ng,nl=counts(data);rg,rl=counts(data,phase)
        if nl>=70 or (phase.startswith('reserve_') and rl>=8): raise ValueError('Local attempt cap reached')
    return L.submit(state,data,project,tid,mode,role,phase)

def pair_reservation(state,data,t):
    tid=t['task_id'].replace('_global','_local')
    if get_job(data,tid,'local') or any(r['task']==tid for r in data['reservations']): return
    if counts(data)[1]>=70 or (t['phase'].startswith('reserve_') and counts(data,t['phase'])[1]>=8):
        raise ValueError('Dependent local attempt cannot be reserved')
    if committed(data)+102>4000: raise ValueError('Global and dependent local cannot both be funded')
    data['reservations'].append(dict(task=tid,role='local',phase=t['phase'],reserved_ch='6'))
    state.event(data,'dependent_local_reserved',tid);state.save(data)

def release_pair(state,data,t):
    tid=t['task_id'].replace('_global','_local')
    data['reservations'][:]=[r for r in data['reservations'] if r['task']!=tid]
    state.save(data)

def ready(state,data,project,t):
    tid=t['task_id'];j=get_job(data,tid,t['mode'])
    if j is None:
        if t['mode']=='global': pair_reservation(state,data,t)
        submit(state,data,project,tid,t['mode'],t['mode'],t['phase']);return False
    if not is_closed(j): return False
    if not job_ok(j):
        if t['mode']=='global':release_pair(state,data,t)
        return True
    c=get_job(data,tid,'check')
    if c is None:submit(state,data,project,tid,'check','check',t['phase']);return False
    if not is_closed(c):return False
    if not job_ok(c):
        # Only a demonstrated infrastructure failure receives the one allowed retry.
        retried=any(x['task']==tid and x['mode']=='check_failed_1' for x in data['jobs'])
        if not retried and c['state'] in {'NODE_FAIL','BOOT_FAIL','PREEMPTED'}:
            c['mode']='check_failed_1';state.event(data,'artifact_check_infrastructure_retry',tid);state.save(data)
            submit(state,data,project,tid,'check','check',t['phase']);return False
        raise ValueError('Independent artifact check failed: '+tid)
    if accepted(data,state.study,tid) is None:raise ValueError('Missing validated payload '+tid)
    return True

def chain(state,data,project,t):
    add(state,data,t);screen=get_job(data,t['task_id'],'screen')
    if screen is None:submit(state,data,project,t['task_id'],'screen','screen',t['phase']);return False
    if not is_closed(screen):return False
    if not job_ok(screen):return True
    if not ready(state,data,project,t):return False
    result=accepted(data,state.study,t['task_id'])
    if result is None:return True
    local=task(state.study,t['task_id'].replace('_global','_local'),'local',t['phase'],t['point'],t['origin'],t['family'],source(state.study,project,t['task_id']))
    add(state,data,local);return ready(state,data,project,local)

def imports(study,project):
    rows=[]
    for r in read_table(study/'imports.tsv'):
        path=study/r['relative_path']
        rows.append(dict(task_id=r['import_id'],point=r['point'],origin=r['origin'],family=r['family'],
                         objective=float(r['objective']),local_stable=r['point']=='base',mode='import',phase=r['kind'],alternative=False,
                         source_fit=str(path.relative_to(project)),start_objective=float(r['objective'])))
    return rows

def candidates(data,study,project):
    rows=imports(study,project)
    for t in data['tasks']:
        if t['phase']=='fixture' or t['mode'] not in {'local','global'}:continue
        r=accepted(data,study,t['task_id'])
        if r and not r.get('skipped',False):
            r.update(phase=t['phase'],source_fit=source(study,project,t['task_id']));rows.append(r)
    return sorted(rows,key=lambda r:(r['objective'],r['task_id']))

def gates(data,study,project):
    cs=candidates(data,study,project);base=[r for r in cs if r['point']=='base'];bq=base[0]['objective']
    baseline_roots={}
    for r in base:
        if r['family'] in {'A','B'} and r['mode'] in {'local','import'} and r['origin'] not in baseline_roots:
            baseline_roots[r['origin']]=r
    near=[r for r in baseline_roots.values() if r['objective']-bq<=.001 and r['local_stable']]
    repeats=[r for r in base if r['mode']=='local' and r.get('alternative',False) and r['local_stable'] and r['objective']-bq<=1e-8 and r['start_objective']-r['objective']<1e-4]
    bg=dict(objective=bq,task_id=base[0]['task_id'],families=';'.join(sorted({r['family'] for r in near})),
            recovered_origins=len(near),repeat_pass=bool(repeats))
    bg['pass']=len(near)>=2 and {r['family'] for r in near}=={'A','B'} and bool(repeats)
    rows=[]
    for c in read_table(study/'cases.tsv'):
        point=c['point'];rs=[r for r in cs if r['point']==point];best=rs[0];q=best['objective'];roots={}
        for r in rs:
            if r['mode']!='local' or r['family'] not in {'A','B'} or r['origin'] in roots:continue
            gt=r['origin']+'_global'
            if accepted(data,study,gt) is not None:roots[r['origin']]=r
        raw=[r for r in roots.values() if r['objective']-q<=.001]
        good=[r for r in raw if r['local_stable']]
        repeat=[r for r in rs if r['mode']=='local' and r.get('alternative',False) and (r['phase']=='repeat' or r['phase'].startswith('reserve_')) and r['local_stable'] and r['objective']-q<=1e-8 and r['start_objective']-r['objective']<1e-4]
        local_pass=any(r['mode']=='local' and r['local_stable'] and r['objective']-q<=1e-8 for r in rs)
        independent=len(good)>=2 and {r['family'] for r in good}=={'A','B'}
        rows.append(dict(point=point,target=c['target'],value=float(c['value']),objective=q,delta=q-bq,
                         best_task=best['task_id'],independent_pass=independent,
                         recovered_origins=len(good),families=';'.join(sorted({r['family'] for r in good})),
                         raw_near_families=';'.join(sorted({r['family'] for r in raw})),local_pass=local_pass,
                         repeat_pass=bool(repeat),reference_pass=bg['pass'] and q>=bq-1e-8,
                         pass_gate=independent and local_pass and bool(repeat) and bg['pass'] and q>=bq-1e-8))
    return bg,rows,cs

def emit(state,data,project,status):
    bg,points,cs=gates(data,state.study,project)
    write_table(state.study/'tables/point_quality.tsv',points,list(points[0]))
    write_kv(state.study/'tables/baseline_gate.tsv',bg)
    cols=['task_id','point','origin','family','mode','phase','objective','local_stable','source_fit','start_objective']
    write_table(state.study/'tables/candidates.tsv',[{k:r[k] for k in cols} for r in cs],cols)
    write_table(state.study/'tables/jobs.tsv',data['jobs'],L.JOB_FIELDS)
    text='# Seven point culture review status\n\nStatus: %s. Updated: %s.\n\nBaseline Q: %.17g; baseline gate: %s. Point gates passed: %d/7.\n\nAllocated use or commitments: %.6f / 4000 core-hours. Non-fit: %.6f / 120. Global-origin attempts: %d/36; fitting-local attempts including standing reservations: %d/70.\n\nNo application concurrency throttle. No tumor fitting, additional grid points or expanded domains are authorized.\n'%(status,now(),bg['objective'],bg['pass'],sum(r['pass_gate'] for r in points),committed(data),committed(data,'nonfit'),*counts(data))
    atomic_text(state.study/'STATUS.md',text)
    return bg,points,cs

def best(cs,point):return next(r for r in cs if r['point']==point)

def round_number(study):
    for n in (1,2):
        d=study/'reserve'/str(n)
        if not (d/'CLOSED').exists():return n
    return None

def open_round(state,data,n,reason):
    path=state.study/'reserve'/str(n)/'decision.tsv'
    if not path.exists():
        write_kv(path,dict(round=n,reason=reason,utc=now(),baseline_q=gates(data,state.study,state.study.parents[3])[0]['objective']))
        state.event(data,'reserve_round_opened',str(n)+' '+reason);state.save(data)

def next_local_id(data,phase):
    used=[r for r in data['tasks'] if r['phase']==phase and r['task_id'].startswith(phase+'_repair_')]
    return phase+'_repair_'+str(len(used)+1)

def repair_baseline(state,data,project):
    bg,pg,cs=gates(data,state.study,project)
    deficit=not bg['pass'] or min(r['delta'] for r in pg)<-1e-8
    if not deficit:return True
    n=round_number(state.study)
    if n is None:return False
    phase='reserve_'+str(n);open_round(state,data,n,'baseline_priority')
    active=[read_kv(state.study/'tasks'/(r['task_id']+'.tsv')) for r in data['tasks'] if r['phase']==phase and r['point']=='base']
    pending=False
    for t in active:
        ok=chain(state,data,project,t) if t['mode']=='global' else ready(state,data,project,t)
        pending=not ok or pending
    if pending:return None
    bg,pg,cs=gates(data,state.study,project)
    if bg['pass'] and min(r['delta'] for r in pg)>=-1e-8:
        # An early intervention consumes this reserve round, even if few slots were used.
        if not (state.study/'reserve'/str(n)/'global_choices.tsv').exists():atomic_text(state.study/'reserve'/str(n)/'CLOSED',now()+' baseline repaired\n')
        return True
    rg,rl=counts(data,phase)
    if rl>=8:
        atomic_text(state.study/'reserve'/str(n)/'CLOSED',now()+' baseline unresolved at round cap\n');return None
    base_best=best(cs,'base');lowest=cs[0]
    need_local=lowest['objective']<bg['objective']-1e-8 or not bg['repeat_pass']
    if need_local:
        tid=next_local_id(data,phase)
        t=task(state.study,tid,'local',phase,'base',lowest['origin'],lowest['family'],lowest['source_fit'],True)
        add(state,data,t);ready(state,data,project,t);return None
    # Missing cross-family baseline recovery: at most one new A and B origin per round.
    missing={'A','B'}-set(bg['families'].split(';'))
    targets=sorted(missing) if missing else ['A','B']
    launched=False
    for family in targets:
        tid='base_R'+str(n)+'_'+family+'_global'
        if any(r['task_id']==tid for r in data['tasks']):continue
        if counts(data,phase)[0]>=4:break
        t=task(state.study,tid,'global',phase,'base','base_R'+str(n)+'_'+family,family)
        add(state,data,t);chain(state,data,project,t);launched=True
    if not launched:
        atomic_text(state.study/'reserve'/str(n)/'CLOSED',now()+' baseline independent recovery unresolved\n')
    return None

def reserve_points(state,data,project,n):
    phase='reserve_'+str(n);open_round(state,data,n,'point_review')
    path=state.study/'reserve'/str(n)/'global_choices.tsv'
    if not path.exists():
        bg,pg,cs=gates(data,state.study,project);choices=[]
        for p in pg:
            if p['independent_pass']:continue
            missing={'A','B'}-set(p['raw_near_families'].split(';'))
            if not missing:continue # Repair local support before spending another global.
            if len(choices)>=4:break
            completed={fam:sum(r['point']==p['point'] and r['family']==fam and r['mode']=='global' for r in cs) for fam in ('A','B')}
            fam=sorted(missing,key=lambda f:(completed[f],f))[0]
            choices.append(dict(point=p['point'],family=fam,origin=p['point']+'_R'+str(n)+'_'+fam))
        write_table(path,choices,['point','family','origin'])
        state.event(data,'reserve_global_choices',str(n)+' '+sha(path));state.save(data)
    done=[]
    for r in read_table(path):
        t=task(state.study,r['origin']+'_global','global',phase,r['point'],r['origin'],r['family'])
        done.append(chain(state,data,project,t))
    if not all(done):return False
    # General local repairs are chosen serially because they depend on updated minima.
    while counts(data,phase)[1]<8:
        bg,pg,cs=gates(data,state.study,project)
        if all(p['pass_gate'] for p in pg):break
        pending=[read_kv(state.study/'tasks'/(r['task_id']+'.tsv')) for r in data['tasks'] if r['phase']==phase and r['task_id'].startswith(phase+'_repair_')]
        for t in pending:
            if not ready(state,data,project,t):return False
        chosen=None
        for p in pg:
            if p['pass_gate']:continue
            eligible=[r for r in cs if r['point']==p['point'] and r['objective']-p['objective']<=.001 and not r['local_stable'] and r['family'] in {'A','B'}]
            if eligible:chosen=eligible[0];break
            if not p['local_pass'] or not p['repeat_pass']:chosen=best(cs,p['point']);break
        if chosen is None:break
        tid=next_local_id(data,phase);t=task(state.study,tid,'local',phase,chosen['point'],chosen['origin'],chosen['family'],chosen['source_fit'],True)
        add(state,data,t);ready(state,data,project,t);return False
    # Close all already released general repairs even at the slot ceiling.
    for r in data['tasks']:
        if r['phase']==phase and r['mode']=='local' and r['task_id'].startswith(phase+'_repair_'):
            if not ready(state,data,project,read_kv(state.study/'tasks'/(r['task_id']+'.tsv'))):return False
    atomic_text(state.study/'reserve'/str(n)/'CLOSED',now()+' reserve decision closed\n');return True

def drive(state,data,project):
    L.refresh(state,data);L.verify_code(state.study)
    if (state.study/'tables/resource_stop.tsv').exists():return 'complete_unresolved_resource'
    release=read_kv(state.study/'settings.tsv').get('validation_release',1)
    suffix='' if release==1 else '_release_2'
    for tid,mode,role in [('validation'+suffix,'validate','validate'),('validation_local'+suffix,'validate_local','check')]:
        j=get_job(data,tid,mode)
        if j is None:submit(state,data,project,tid,mode,role,'validation');return 'validation'
        if not is_closed(j):return 'validation'
        if not job_ok(j):raise ValueError('Validation release failed; scientific fitting held')
    for name in ['global_validation.tsv','local_validation.tsv']:
        v=read_kv(state.study/'tests'/name)
        if not v['pass'] or v['code_sha256']!=sha(state.study/'manifests/code_manifest.tsv'):raise ValueError('Stale validation gate')
    if not (state.study/'tests/validation.tsv').exists():
        write_kv(state.study/'tests/validation.tsv',dict(pass_gate=True,code_sha256=sha(state.study/'manifests/code_manifest.tsv'),utc=now()))
        write_kv(state.study/'manifests/method_freeze.tsv',dict(code_sha256=sha(state.study/'manifests/code_manifest.tsv'),utc=now()))
    imported={r['import_id']:r for r in read_table(state.study/'imports.tsv')}
    guards=[]
    for i,key in enumerate(['base_reference','base_A2'],1):
        r=imported[key];t=task(state.study,'guard_'+str(i),'local','guard','base',r['origin'],r['family'],str((state.study/r['relative_path']).relative_to(project)),True)
        add(state,data,t);guards.append(ready(state,data,project,t))
    if not all(guards):return 'unrestricted_guard'
    repair=repair_baseline(state,data,project)
    if repair is False:return 'complete_unresolved_reference'
    if repair is None:return 'baseline_reserve_repair'
    outstanding=[]
    for row in data['tasks']:
        if row['mode']=='global' and row['phase'].startswith('reserve_') and get_job(data,row['task_id'],'global') is not None:
            outstanding.append(chain(state,data,project,read_kv(state.study/'tasks'/(row['task_id']+'.tsv'))))
    if not all(outstanding):return 'completing_committed_reserve_chains'
    done=[]
    for c in read_table(state.study/'cases.tsv'):
        p=c['point'];r=imported['historical_'+p]
        t=task(state.study,p+'_historical','local','historical',p,r['origin'],'',str((state.study/r['relative_path']).relative_to(project)))
        add(state,data,t);done.append(ready(state,data,project,t))
        for fam in ('A','B'):
            for k in (1,2):
                origin=p+'_'+fam+str(k);t=task(state.study,origin+'_global','global','primary',p,origin,fam)
                done.append(chain(state,data,project,t))
    if not all(done):return 'primary_chains'
    snap=state.study/'manifests/neighbor_sources.tsv'
    if not snap.exists():
        cs=candidates(data,state.study,project)
        write_table(snap,[dict(point=c['point'],source_fit=best(cs,c['point'])['source_fit']) for c in read_table(state.study/'cases.tsv')],['point','source_fit'])
    frozen={r['point']:r['source_fit'] for r in read_table(snap)};done=[]
    for target,seq in [('pm',['pm_00','pm_05','pm_06','pm_09']),('bb',['bb_00','bb_07','bb_08'])]:
        for direction,order in [('forward',seq),('reverse',list(reversed(seq)))]:
            previous=frozen[order[0]]
            for point in order[1:]:
                tid=point+'_neighbor_'+direction
                t=task(state.study,tid,'local','neighbor',point,tid,'',previous,False,frozen[point])
                add(state,data,t)
                if not ready(state,data,project,t):done.append(False);break
                r=accepted(data,state.study,tid)
                previous=frozen[point] if r is None or r.get('skipped') else source(state.study,project,tid)
            else:done.append(True)
    if not all(done):return 'bidirectional_neighbor_continuation'
    cs=candidates(data,state.study,project);done=[]
    for c in read_table(state.study/'cases.tsv'):
        tid=c['point']+'_repeat';existing=next((r for r in data['tasks'] if r['task_id']==tid),None)
        if existing:t=read_kv(state.study/'tasks'/(tid+'.tsv'))
        else:
            r=best(cs,c['point']);t=task(state.study,tid,'local','repeat',c['point'],r['origin'],r['family'],r['source_fit'],True);add(state,data,t)
        done.append(ready(state,data,project,t))
    if not all(done):return 'final_point_repeats'
    for n in (1,2):
        bg,pg,cs=gates(data,state.study,project)
        if all(p['pass_gate'] for p in pg):return 'complete_pass'
        if (state.study/'reserve'/str(n)/'CLOSED').exists():continue
        if not reserve_points(state,data,project,n):return 'point_reserve_'+str(n)
    return 'complete_pass' if all(p['pass_gate'] for p in gates(data,state.study,project)[1]) else 'complete_unresolved_points'

def validation(project,study,local=False):
    audit_locked(study)
    run=lambda mode,tid:subprocess.run(['bash',str(study/'implementation/run.sh'),str(project),str(study),mode,tid],check=True)
    if not local:
        subprocess.run([__import__('sys').executable,str(study/'implementation/test_workflow.py')],check=True)
        run('unit','unit')
        for r in read_table(study/'imports.tsv'):
            reps=1 if r['import_id'] in {'base_A1','base_B4'} else 3
            for k in range(1,reps+1):run('replay',r['import_id']+'_rep'+str(k))
    cases=['base','pm_00','bb_00','pm_05','bb_07']
    for point in cases:
        tid='fixture_'+point+('_local' if local else '_global')
        run('local' if local else 'global',tid);run('check',tid)
    audit_locked(study)
    write_kv(study/'tests'/('local_validation.tsv' if local else 'global_validation.tsv'),
             {'pass':True,'code_sha256':sha(study/'manifests/code_manifest.tsv'),'utc':now(),'fixtures':5,'fresh_process_replays':0 if local else 32})

def audit_locked(study):
    # A recursive scan must not race the rename of an unpublished TSV generation.
    with State(study).locked():L.verify_code(study);L.audit(study)

def controller(project,study,once=False):
    state=State(study);start=time.monotonic()
    while True:
        with state.locked():
            data=state.load()
            try:
                status=drive(state,data,project);emit(state,data,project,status)
                if status.startswith('complete_'):
                    drained=[]
                    for r in data['tasks']:
                        if r['mode'] not in {'global','local'} or r['phase']=='fixture':continue
                        j=get_job(data,r['task_id'],r['mode'])
                        if j is not None:
                            try:drained.append(ready(state,data,project,read_kv(study/'tasks'/(r['task_id']+'.tsv'))))
                            except ValueError as exc:
                                if not (study/'tables/resource_stop.tsv').exists() or not resource_limit(exc):raise
                                # Unfunded acceptance stays unaccepted; raw output is preserved.
                                drained.append(is_closed(j))
                    if not all(drained):status='closing_submitted_fits_and_checks'
                    elif data['reservations']:
                        state.event(data,'unlaunched_dependent_slots_released_at_stop',';'.join(r['task'] for r in data['reservations']))
                        data['reservations']=[];state.save(data)
                if status.startswith('complete_') and not (study/'tables/resource_stop.tsv').exists():
                    snapshot=study/'manifests/diagnostic_sources.tsv'
                    if not snapshot.exists():
                        cs=candidates(data,study,project)
                        write_table(snapshot,[{k:best(cs,point)[k] for k in ['point','origin','family','source_fit']} for point in ['base']+[r['point'] for r in read_table(study/'cases.tsv')]],['point','origin','family','source_fit'])
                    diagnostic_done=[]
                    for r in read_table(snapshot):
                        tid='diagnostic_'+r['point'];t=task(study,tid,'diagnostic','diagnostic',r['point'],r['origin'],r['family'],r['source_fit']);add(state,data,t)
                        j=get_job(data,tid,'diagnostic')
                        if j is None:submit(state,data,project,tid,'diagnostic','check','diagnostic');diagnostic_done.append(False)
                        else:diagnostic_done.append(is_closed(j))
                    if not all(diagnostic_done):status='final_event_diagnostics'
                if status.startswith('complete_'):
                    # Valid work already submitted closes even when a scientific gate is unresolved.
                    live=[j for j in data['jobs'] if not is_closed(j) and j['role']!='controller']
                    if live:status='closing_submitted_work'
                    else:
                        report=get_job(data,'final_report','report')
                        if report is None:submit(state,data,project,'final_report','report','report','report');status='final_report'
                        elif not is_closed(report):status='final_report'
                        elif not job_ok(report) or not (study/'reports/COMPLETE').exists():raise ValueError('Final report failed')
                if status.startswith('complete_'):
                    write_kv(study/'tables/final_status.tsv',dict(status=status,utc=now()))
                    atomic_text(study/'CLOSURE.md',(study/'STATUS.md').read_text()+'\nFinal accounting must be refreshed after this controller closes.\n');return
                if once:return
                if time.monotonic()-start>42300:
                    submit(state,data,project,'controller_'+uuid.uuid4().hex[:8],'controller','controller','control');return
            except BaseException as exc:
                if resource_limit(exc):
                    state.event(data,'resource_stop',str(exc));state.save(data)
                    write_kv(study/'tables/resource_stop.tsv',dict(utc=now(),reason=str(exc)))
                    emit(state,data,project,'resource_stop_closing')
                    if once:return
                    continue
                state.event(data,'dependent_submissions_held',str(exc));state.save(data)
                atomic_text(study/'BLOCKED.txt',now()+' '+str(exc)+'\n');emit(state,data,project,'blocked');raise
        time.sleep(30)

def resource_limit(exc):
    return isinstance(exc,ValueError) and any(text in str(exc) for text in ['allocation ceiling','attempt cap','cannot both be funded','cannot be reserved'])

def main():
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['start','controller','validate','validate_local','audit']);p.add_argument('project',type=Path);p.add_argument('study',type=Path);p.add_argument('--once',action='store_true');a=p.parse_args()
    if a.mode=='audit':audit_locked(a.study)
    elif a.mode=='validate':validation(a.project,a.study)
    elif a.mode=='validate_local':validation(a.project,a.study,True)
    elif a.mode=='controller':controller(a.project,a.study,a.once)
    elif a.mode=='start':
        state=State(a.study)
        with state.locked():
            release=read_kv(a.study/'settings.tsv').get('validation_release',1)
            if release not in (1,2):raise ValueError('Validation release allowance exhausted')
            suffix='' if release==1 else '_release_2'
            data=state.load();submit(state,data,a.project,'validation'+suffix,'validate','validate','validation')
            if not any(j['role']=='controller' and not is_closed(j) for j in data['jobs']):submit(state,data,a.project,'controller_start'+suffix,'controller','controller','control')

if __name__=='__main__':main()
