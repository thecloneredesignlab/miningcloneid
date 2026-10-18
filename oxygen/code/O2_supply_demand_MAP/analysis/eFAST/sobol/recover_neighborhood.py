#!/usr/bin/env python3
"""Resume the original N4096 analysis without recalculating certified cells.

Original model inputs/snapshots are immutable. Legacy statistical identities are
accepted only after matching the original B-block identity and input manifest.
Independent shards never write the original statistical checkpoint.
"""
import argparse
from contextlib import contextmanager
import fcntl
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil

import numpy as np
import core as c
import full_core as f
import full_report as report

EXPECTED_MISSING = [25318, 27551, 27552, 33636, 33637, 41044, 45450,
                    48408, 53810, 63791, 68684, 69826, 71559, 74956]
CHUNK = 16


@contextmanager
def locked(path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def blank():
    shape = (len(f.NS), 5, len(c.OXYGEN), 2, len(f.METRICS))
    return (np.full(shape, np.nan), np.full(shape, np.nan),
            np.full(shape[:-1], np.nan), np.zeros(shape[1:4], bool))


def validate_arrays(values, widths, variance, done):
    shape = (len(f.NS), 5, len(c.OXYGEN), 2, len(f.METRICS))
    if (values.shape != shape or widths.shape != shape or variance.shape != shape[:-1]
            or done.shape != shape[1:4] or done.dtype != np.dtype(bool)
            or any(a.dtype != np.dtype('float64') for a in (values, widths, variance))):
        raise ValueError('Statistical checkpoint axes/dtypes differ')
    if (not np.isfinite(variance[:, done]).all() or np.any(variance[:, done] < 0)
            or np.isinf(values[:, done]).any() or np.isinf(widths[:, done]).any()):
        raise ValueError('Invalid certified statistical cells')
    if any(np.isfinite(a[:, ~done]).any() for a in (values, widths, variance)):
        raise ValueError('Uncertified checkpoint cells contain results')


def legacy_keys(p, case, units, bounds):
    return {f.identity(p, case, r, c.STEP-1,
                       f.scale_unit(units[r][c.STEP-1::c.STEP], bounds)) for r in range(1, 6)}


def load_legacy(source, allowed):
    """Recover a valid temporary checkpoint only if it extends the durable one."""
    candidates = []; rejected = []
    for name in ('indices_checkpoint.npz', 'indices_checkpoint.npz.incomplete'):
        path = Path(source)/name
        if not path.exists(): continue
        try:
            with np.load(path, allow_pickle=False) as z:
                if str(z['identity']) not in allowed:
                    raise ValueError('Legacy identity does not match original B block')
                arrays = tuple(z[k] for k in ('values', 'widths', 'variance', 'done'))
                identity = str(z['identity'])
            validate_arrays(*arrays)
            candidates.append((arrays, dict(path=str(path), sha256=c.sha(path),
                                            legacy_identity=identity, n_done=int(arrays[3].sum()))))
        except Exception as error:
            if name == 'indices_checkpoint.npz': raise
            rejected.append(dict(path=str(path), error=str(error)))
    if not candidates:
        return blank(), dict(path=None, n_done=0, rejected_temporary=rejected)
    selected = max(candidates, key=lambda x: x[1]['n_done'])
    for arrays, record in candidates:
        overlap = arrays[3]
        if not selected[0][3][overlap].all(): raise ValueError('Checkpoints have conflicting progress')
        for old, new in zip(arrays[:3], selected[0][:3]):
            np.testing.assert_array_equal(old[:, overlap], new[:, overlap])
    selected[1]['rejected_temporary'] = rejected
    return selected


def verify_summary(target):
    receipt = c.load_json(target/'completion.json')
    if receipt['status'] != 'complete': raise ValueError('Incomplete case receipt')
    fields = {'indices_summary.npz': 'summary_sha256',
              'indices_repeats.npz': 'repeat_indices_sha256'}
    for filename, field in fields.items():
        if c.sha(target/filename) != receipt[field]: raise ValueError('Changed completed case: '+str(target))
    for name in ('indices_summary.tsv.gz', 'numerical_qc.tsv'):
        if not (target/name).is_file(): raise ValueError('Missing completed source table')
    return receipt


def check_runtime(p):
    if version('SALib') != '1.5.2': raise ValueError('Require original SALib 1.5.2')
    sif = Path(p['sif']['path'])
    if (sif.stat().st_size != p['sif']['bytes'] or sif.stat().st_mtime_ns != p['sif']['mtime_ns']
            or p['sif']['sha256'] != c.SIF_SHA256): raise ValueError('Original SIF changed')


def prepare(root, recovery, revision):
    plan_path = root/'slurm/submission_plan.json'; p = c.load_json(plan_path); c.configure(p)
    check_runtime(p)
    if (p['N'], p['n_repeats'], len(c.OXYGEN), p['n_parameters'], p['n_pairs']) != (4096, 5, 201, 14, 91):
        raise ValueError('Original full analysis contract differs')
    if c.OXYGEN != [i*.025 for i in range(201)] and not np.allclose(c.OXYGEN, np.arange(201)*.025, rtol=0, atol=1e-14):
        raise ValueError('Oxygen grid differs')
    if (recovery/'plan.json').exists(): raise ValueError('Recovery already prepared; reuse its plan')
    for name, digest in p['code_sha256'].items():
        if c.sha(Path(p['code_snapshot_root'])/name) != digest: raise ValueError('Original frozen source changed')
    if c.sha(f.HERE/'evaluate.R') != p['engine_sha256']: raise ValueError('Recovery model engine differs')
    if c.sha(p['sif']['path']) != p['sif']['sha256']: raise ValueError('SIF SHA256 differs')
    units = {}
    for repeat in range(1, 6):
        path = root/'inputs'/('unit_R%02d.npz'%repeat)
        if c.sha(path) != p['unit_sha256'][str(repeat)]: raise ValueError('Unit samples changed')
        with np.load(path) as z: units[repeat] = z['unit']
    ready = []; waiting = []; missing = []; cases = []; reused = []; source_records = []
    total_done = 0; cache_bytes = 0
    for case in p['cases']:
        source = f.case_dir(root, case); bounds = c.read(source/'parameter_bounds.tsv')
        if c.sha(source/'parameter_bounds.tsv') != case['bounds_sha256']:
            raise ValueError('Case bounds changed')
        if c.load_json(source/'case.json') != case: raise ValueError('Case metadata changed')
        absent = {}
        for repeat in range(1, 6):
            for bi, block in enumerate(c.BLOCKS):
                path = source/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz'
                receipt_path = Path(str(path)+'.receipt.json')
                task_id = (case['case_index']-1)*150+(repeat-1)*30+bi+1
                if not path.is_file() or not receipt_path.is_file():
                    missing.append(dict(task_id=task_id, case_id=case['case_id'], repeat=repeat, block=block,
                                        model_checkpoint=(path.parent/'checkpoint.npz').is_file()))
                    absent.setdefault(repeat, []).append(task_id)
                else:
                    rr = c.load_json(receipt_path)
                    if rr['status'] != 'complete' or rr['n_points'] != f.N*len(c.OXYGEN):
                        raise ValueError('Malformed raw completion receipt')
                    source_records.append(dict(path=str(path.relative_to(root)), sha256=rr['sha256'],
                                               receipt_sha256=c.sha(receipt_path)))
        allowed = legacy_keys(p, case, units, bounds)
        if (source/'completion.json').exists():
            rr = verify_summary(source)
            if rr['analysis_identity'] not in allowed: raise ValueError('Unexpected legacy completed identity')
            if absent: raise ValueError('Completed seed is missing model outputs')
            reused.append(dict(case_id=case['case_id'], scope=case['scope'],
                               completion_sha256=c.sha(source/'completion.json'),
                               summary_sha256=rr['summary_sha256'], repeat_sha256=rr['repeat_indices_sha256'],
                               table_sha256=c.sha(source/'indices_summary.tsv.gz'),
                               qc_sha256=c.sha(source/'numerical_qc.tsv')))
            continue
        if case['scope'] != 'neighborhood10pct': raise ValueError('Global must already be complete')
        arrays, migration = load_legacy(source, allowed)
        total_done += int(arrays[3].sum())
        checkpoint = recovery/'checkpoints'/(case['case_id']+'.npz')
        f.atomic_npz(checkpoint, values=arrays[0], widths=arrays[1], variance=arrays[2], done=arrays[3])
        item = dict(case=case, legacy=migration, checkpoint_sha256=c.sha(checkpoint), task_ids=[])
        for repeat in range(1, 6):
            oi = np.flatnonzero(~arrays[3][repeat-1].all(axis=1)).tolist()
            cache_bytes += len(oi)*f.N*c.STEP*2*8
            for start in range(0, len(oi), CHUNK):
                oxygen = oi[start:start+CHUNK]
                task = dict(case_id=case['case_id'], case_index=case['case_index'], repeat=repeat,
                            oxygen_indices=oxygen,
                            needed_outputs=[[k for k in range(2) if not arrays[3][repeat-1, j, k]] for j in oxygen],
                            pending_model_tasks=absent.get(repeat, []))
                (waiting if task['pending_model_tasks'] else ready).append(task)
        cases.append(item)
        if len(cases)%25 == 0: print('Audited pending seeds', len(cases), flush=True)
    if sorted(r['task_id'] for r in missing) != EXPECTED_MISSING:
        raise ValueError('Missing model set differs from authorized 14; inspect before submission')
    if len(cases) != 355 or len(reused) != 146: raise ValueError('Original 355/145/global partition changed')
    if shutil.disk_usage(root).free < cache_bytes+300*1024**3: raise ValueError('Insufficient cache and reserve space')
    tasks = ready+waiting; by_case = {r['case']['case_id']: r for r in cases}
    for tid, task in enumerate(tasks, 1):
        task['task_id'] = tid; by_case[task['case_id']]['task_ids'].append(tid)
    snapshot = f.REPO
    code_hashes = {str(path.relative_to(snapshot)): c.sha(path) for path in snapshot.rglob('*')
                   if path.is_file() and '__pycache__' not in path.parts}
    imports=[f.HERE/name for name in ('recover_neighborhood.py','core.py','full_core.py','full_report.py','report.py','run.py','evaluate.R')]
    imports += [f.HERE.parent/name for name in ('efast.py','neighborhood.py')]
    runtime_hashes={str(path.relative_to(snapshot)):c.sha(path) for path in imports}
    c.table(recovery/'source_outputs.tsv', source_records)
    c.table(recovery/'missing_model_blocks.tsv', missing)
    c.table(recovery/'reused_cases.tsv', reused)
    for row in cases: c.save_json(recovery/'migrations'/(row['case']['case_id']+'.json'), row)
    rp = dict(source_root=str(root), source_plan_sha256=c.sha(plan_path), source_git_commit=p['git_commit'],
              recovery_git_commit=revision, code_snapshot_root=str(snapshot), code_sha256=code_hashes,runtime_sha256=runtime_hashes,
              cases=cases, reused_cases=reused, model_tasks=missing, tasks=tasks,
              n_ready_tasks=len(ready), n_waiting_tasks=len(waiting), chunk_oxygen=CHUNK,
              n_saved_cells=total_done, n_remaining_cells=len(cases)*5*len(c.OXYGEN)*2-total_done,
              cache_bytes=cache_bytes, sif=p['sif'], created_utc=c.now(),
              convergence_policy='diagnostic_only', qos='xxlarge', time='12:00:00',
              specified_node=None, array_concurrency_limit=None)
    c.save_json(recovery/'plan.json', rp)
    c.save_json(recovery/'status.json', dict(status='prepared', n_pending_seeds=len(cases),
        n_model_tasks=len(missing), n_shards=len(tasks), n_saved_cells=total_done,
        n_remaining_cells=rp['n_remaining_cells'], updated_utc=c.now()))
    print(json.dumps({k:rp[k] for k in ['n_ready_tasks','n_waiting_tasks','n_saved_cells','n_remaining_cells','cache_bytes']}),flush=True)


def context(recovery):
    rp = c.load_json(recovery/'plan.json')
    expected = os.environ.get('RECOVERY_PLAN_SHA256')
    if expected and c.sha(recovery/'plan.json') != expected: raise ValueError('Recovery plan changed')
    root = Path(rp['source_root'])
    if c.sha(root/'slurm/submission_plan.json') != rp['source_plan_sha256']: raise ValueError('Original plan changed')
    p = c.load_json(root/'slurm/submission_plan.json'); c.configure(p); check_runtime(p)
    if str(f.REPO) != rp['code_snapshot_root']: raise ValueError('Use immutable recovery snapshot')
    for name, digest in rp['runtime_sha256'].items():
        # R sources/model math are checked too; no mutable checkout code is imported.
        if c.sha(f.REPO/name) != digest: raise ValueError('Frozen recovery source changed: '+name)
    return root, p, rp


def source_cache(root, recovery, p, rp, task):
    case = p['cases'][task['case_index']-1]; repeat = task['repeat']
    directory = recovery/'cache'/case['case_id']/('R%02d'%repeat)
    directory.mkdir(parents=True, exist_ok=True)
    data = directory/'outputs.npy'; receipt = directory/'completion.json'
    oxygen = sorted({oi for t in rp['tasks'] if t['case_id']==case['case_id'] and t['repeat']==repeat
                     for oi in t['oxygen_indices']})
    with locked(directory/'cache.lock'):
        if receipt.exists():
            rr = c.load_json(receipt)
            if (rr['plan_sha256'] != c.sha(recovery/'plan.json') or rr['oxygen_indices'] != oxygen
                    or data.stat().st_size != rr['bytes']): raise ValueError('Cache receipt changed')
            return data, oxygen, rr
        source = f.case_dir(root, case); bounds = c.read(source/'parameter_bounds.tsv')
        if c.sha(source/'parameter_bounds.tsv') != case['bounds_sha256']: raise ValueError('Bounds changed')
        unit_path=root/'inputs'/('unit_R%02d.npz'%repeat)
        if c.sha(unit_path)!=p['unit_sha256'][str(repeat)]: raise ValueError('Unit sample changed')
        with np.load(unit_path) as z: unit=z['unit']
        tmp = data.with_suffix('.npy.incomplete')
        mm = np.lib.format.open_memmap(tmp, mode='w+', dtype=np.float64,
                                      shape=(len(oxygen), f.N*c.STEP, 2))
        records=[]
        for bi, block in enumerate(c.BLOCKS):
            path = source/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz'
            identity=f.identity(p,case,repeat,bi,f.scale_unit(unit[bi::c.STEP],bounds))
            if not f.complete(path,identity): raise ValueError('Missing source block: '+str(path))
            rr=c.load_json(str(path)+'.receipt.json')
            with np.load(path,allow_pickle=False) as z:
                outputs=z['outputs']
                if outputs.shape!=(f.N,len(c.OXYGEN),2) or outputs.dtype!=np.float64 or not np.array_equal(z['oxygen'],c.OXYGEN):
                    raise ValueError('Source axes or precision changed')
                if not np.isfinite(outputs).all(): raise ValueError('Nonfinite raw model output')
                mm[:,bi::c.STEP,:]=outputs[:,oxygen,:].transpose(1,0,2)
            records.append(dict(path=str(path.relative_to(root)),sha256=rr['sha256']))
        mm.flush(); del mm
        with tmp.open('rb') as h: os.fsync(h.fileno())
        tmp.replace(data)
        # Per oxygen hashes allow consumers to check only the slices they read.
        mm=np.load(data,mmap_mode='r'); import hashlib
        hashes=[hashlib.sha256(np.asarray(mm[i]).tobytes()).hexdigest() for i in range(len(oxygen))]
        del mm
        rr=dict(status='complete',plan_sha256=c.sha(recovery/'plan.json'),oxygen_indices=oxygen,
                bytes=data.stat().st_size,slice_sha256=hashes,source_outputs=records,created_utc=c.now())
        c.save_json(receipt,rr)
        data.chmod(0o444)
    return data,oxygen,rr


def shard_identity(recovery, task):
    import hashlib
    return hashlib.sha256((c.sha(recovery/'plan.json')+json.dumps(task,sort_keys=True)).encode()).hexdigest()


def calculate_shard(recovery, task, ordered, oxygen_lookup, analyze):
    """Each analyzer result is persisted, including results at the smaller N."""
    directory=recovery/'shards'/('%05d'%task['task_id']); directory.mkdir(parents=True,exist_ok=True)
    path=directory/'checkpoint.npz'; receipt=directory/'completion.json'
    identity=shard_identity(recovery,task)
    shape=(len(f.NS),len(task['oxygen_indices']),2,len(f.METRICS))
    requested=np.zeros(shape[:-1],bool)
    for j,ks in enumerate(task['needed_outputs']): requested[:,j,ks]=True
    with locked(directory/'worker.lock'):
        if receipt.exists():
            rr=c.load_json(receipt)
            if rr['identity']!=identity or rr['sha256']!=c.sha(path): raise ValueError('Shard receipt changed')
            return
        if path.exists():
            with np.load(path) as z:
                if str(z['identity'])!=identity: raise ValueError('Shard identity changed')
                values,widths,variance,done=[z[k] for k in ('values','widths','variance','done')]
            if values.shape!=shape or widths.shape!=shape or variance.shape!=shape[:-1] or done.shape!=shape[:-1]:
                raise ValueError('Shard axes changed')
        else:
            values=np.full(shape,np.nan); widths=np.full(shape,np.nan)
            variance=np.full(shape[:-1],np.nan); done=~requested
        for j,oi in enumerate(task['oxygen_indices']):
            for k in task['needed_outputs'][j]:
                for ni,n in enumerate(f.NS):
                    if done[ni,j,k]: continue
                    seed=2026101300+task['repeat']+10000+n+oi*10+k
                    si,var=analyze(ordered[oxygen_lookup[oi],:n*c.STEP,k],seed)
                    values[ni,j,k],widths[ni,j,k]=f.metric_arrays(si)
                    variance[ni,j,k]=var; done[ni,j,k]=True
                    f.atomic_npz(path,identity=np.array(identity),values=values,widths=widths,variance=variance,done=done)
            print('Shard',task['task_id'],task['case_id'],'repeat',task['repeat'],'O2 index',oi,'saved',flush=True)
        if not done.all(): raise ValueError('Unfinished shard')
        c.save_json(receipt,dict(status='complete',identity=identity,sha256=c.sha(path),
                               n_analyzer_cells=int(requested.sum()),completed_utc=c.now()))


def shard(root,recovery,p,rp,task_id):
    task=rp['tasks'][task_id-1]
    if task['task_id']!=task_id: raise ValueError('Task mapping differs')
    # A completed shard does not need to open or construct the source cache again.
    directory=recovery/'shards'/('%05d'%task_id)
    if (directory/'completion.json').exists():
        rr=c.load_json(directory/'completion.json')
        if rr['identity']!=shard_identity(recovery,task) or rr['sha256']!=c.sha(directory/'checkpoint.npz'):
            raise ValueError('Completed shard changed')
        return
    path,oxygen,receipt=source_cache(root,recovery,p,rp,task)
    ordered=np.load(path,mmap_mode='r'); lookup={oi:i for i,oi in enumerate(oxygen)}
    import hashlib
    for oi in task['oxygen_indices']:
        if hashlib.sha256(np.asarray(ordered[lookup[oi]]).tobytes()).hexdigest()!=receipt['slice_sha256'][lookup[oi]]:
            raise ValueError('Source cache slice changed')
    bounds=c.read(f.case_dir(root,p['cases'][task['case_index']-1])/'parameter_bounds.tsv')
    calculate_shard(recovery,task,ordered,lookup,lambda y,seed:c.analyze(bounds,y,seed=seed))


def merge_case(recovery,rp,item):
    case=item['case']; path=recovery/'checkpoints'/(case['case_id']+'.npz')
    if c.sha(path)!=item['checkpoint_sha256']: raise ValueError('Migrated checkpoint changed')
    with np.load(path) as z: values,widths,variance,done=[z[k] for k in ('values','widths','variance','done')]
    original_done=done.copy(); originals=[a[:,done].copy() for a in (values,widths,variance)]
    for tid in item['task_ids']:
        task=rp['tasks'][tid-1]; directory=recovery/'shards'/('%05d'%tid)
        rr=c.load_json(directory/'completion.json'); file=directory/'checkpoint.npz'
        if rr['identity']!=shard_identity(recovery,task) or rr['sha256']!=c.sha(file):
            raise ValueError('Missing or changed completed shard')
        with np.load(file) as z:
            if not z['done'].all(): raise ValueError('Partial shard')
            vv,ww,var=z['values'],z['widths'],z['variance']
        for j,oi in enumerate(task['oxygen_indices']):
            for k in task['needed_outputs'][j]:
                if done[task['repeat']-1,oi,k]: raise ValueError('Overlapping computation claims')
                values[:,task['repeat']-1,oi,k]=vv[:,j,k]
                widths[:,task['repeat']-1,oi,k]=ww[:,j,k]
                variance[:,task['repeat']-1,oi,k]=var[:,j,k]
                done[task['repeat']-1,oi,k]=True
    if not done.all(): raise ValueError('Missing five-repeat/O2/output cells')
    for a,old in zip((values,widths,variance),originals): np.testing.assert_array_equal(a[:,original_done],old)
    validate_arrays(values,widths,variance,done)
    return values,widths,variance


def finalize_case(root,recovery,p,rp,case_id):
    item=rp['cases'][case_id-1]; case=item['case']; target=f.case_dir(root,case)
    with locked(target/'summary_recovery.lock'):
        if (target/'completion.json').exists():
            rr=verify_summary(target)
            if rr.get('recovery_plan_sha256')!=c.sha(recovery/'plan.json'): raise ValueError('Different case completion exists')
            return
        values,widths,variance=merge_case(recovery,rp,item)
        f.atomic_npz(target/'indices_repeats.npz',values=values,bootstrap_half_width=widths,variance=variance,
            resolutions=np.array(f.NS),oxygen=np.asarray(c.OXYGEN),metrics=np.array(f.METRICS),outputs=np.array(c.OUTPUTS))
        mean,spread,delta,status=f.diagnostics(values)
        outvar=np.mean(variance[-1],axis=0); width=np.median(widths[-1],axis=0)
        f.atomic_npz(target/'indices_summary.npz',mean=mean,repeat_range=spread,resolution_delta=delta,
                     status=status,variance=outvar,bootstrap_half_width_median=width,oxygen=np.asarray(c.OXYGEN))
        from itertools import chain
        rows=report.summary_rows(case,c.OXYGEN,mean,spread,delta,status,outvar,width); first=next(rows)
        c.efast.write_table(target/'indices_summary.tsv.gz',list(first),chain([first],rows))
        numerical=[]; receipts=[]
        for repeat in range(1,6):
            for block in c.BLOCKS:
                path=target/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz'
                rr=c.load_json(str(path)+'.receipt.json')
                if not path.exists() or rr['status']!='complete': raise ValueError('Missing model output')
                receipts.append(dict(path=str(path.relative_to(root)),sha256=rr['sha256']))
                numerical.append(dict(repeat=repeat,block=block,n_points=rr['n_points'],n_recovered=rr['n_recovered'],
                    minimum_spectral_gap=rr['minimum_spectral_gap'],maximum_residual=rr['maximum_residual']))
        c.table(target/'numerical_qc.tsv',numerical)
        report.plot_indices(mean,c.OXYGEN,target/'figures',case['scope']+'/'+case['case_id']+'; mean of five repeats')
        c.save_json(target/'completion.json',dict(status='complete',analysis_identity=c.sha(recovery/'plan.json')+':'+case['case_id'],
            case=case,N=f.N,repeats=5,oxygen_points=len(c.OXYGEN),summary_sha256=c.sha(target/'indices_summary.npz'),
            repeat_indices_sha256=c.sha(target/'indices_repeats.npz'),summary_table_sha256=c.sha(target/'indices_summary.tsv.gz'),
            numerical_qc_sha256=c.sha(target/'numerical_qc.tsv'),output_receipts=receipts,n_recovered=sum(r['n_recovered'] for r in numerical),
            n_undefined_indices=int((status==0).sum()),n_not_passed_indices=int((status==2).sum()),
            recovery_plan_sha256=c.sha(recovery/'plan.json'),reused_statistical_cells=item['legacy']['n_done'],
            all_original_checkpoint_cells_preserved=True,completed_utc=c.now(),included_in_full_analysis=True,
            convergence_policy='diagnostic_only'))


def audit(root,recovery,p,rp):
    for record in rp['reused_cases']:
        case=next(x for x in p['cases'] if x['case_id']==record['case_id'])
        target=f.case_dir(root,case); verify_summary(target)
        if c.sha(target/'completion.json')!=record['completion_sha256']: raise ValueError('Reused completed case changed')
        if c.sha(target/'indices_summary.tsv.gz')!=record['table_sha256'] or c.sha(target/'numerical_qc.tsv')!=record['qc_sha256']:
            raise ValueError('Reused completed tables changed')
    for item in rp['cases']:
        rr=verify_summary(f.case_dir(root,item['case']))
        if rr.get('recovery_plan_sha256')!=c.sha(recovery/'plan.json'): raise ValueError('Unexpected restored completion')
    c.save_json(recovery/'all_cases_validated.json',dict(status='passed',n_cases=501,n_neighborhoods=500,
        original_completed_cases_unchanged=True,n_saved_cells=rp['n_saved_cells'],N=4096,repeats=5,oxygen_points=201,
        convergence_policy='diagnostic_only',validated_utc=c.now()))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('role',choices=['prepare','shard','case','audit','aggregate','finish'])
    parser.add_argument('--recovery-dir',required=True); parser.add_argument('--root')
    parser.add_argument('--revision'); parser.add_argument('--task-id',type=int)
    args=parser.parse_args(); recovery=Path(args.recovery_dir).resolve()
    if args.role=='prepare': prepare(Path(args.root).resolve(),recovery,args.revision); return
    root,p,rp=context(recovery)
    try:
        if args.role=='shard': shard(root,recovery,p,rp,args.task_id)
        elif args.role=='case': finalize_case(root,recovery,p,rp,args.task_id)
        elif args.role=='audit': audit(root,recovery,p,rp)
        elif args.role=='aggregate':
            if c.load_json(recovery/'all_cases_validated.json')['status']!='passed': raise ValueError('Cases not validated')
            report.aggregate_cases(root,p)
        else:
            report.finish_report(root,p)
            c.save_json(recovery/'status.json',dict(status='complete',n_restored_seeds=len(rp['cases']),
                n_model_blocks=len(rp['model_tasks']),n_shards=len(rp['tasks']),completed_utc=c.now(),convergence_policy='diagnostic_only'))
            rows=c.read(root/'collection_manifest.tsv'); relative=str((recovery/'status.json').relative_to(root))
            rows=[r for r in rows if r['path']!=relative]
            rows.append(dict(path=relative,bytes=(recovery/'status.json').stat().st_size,sha256=c.sha(recovery/'status.json')))
            c.table(root/'collection_manifest.tsv',sorted(rows,key=lambda r:r['path']))
    except Exception as error:
        c.save_json(recovery/'failures'/('%s_%s.json'%(args.role,args.task_id or 'main')),
                    dict(status='failed',error=str(error),recorded_utc=c.now()))
        raise


if __name__=='__main__': main()
