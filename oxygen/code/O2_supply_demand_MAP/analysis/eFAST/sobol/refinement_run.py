#!/usr/bin/env python3
"""Append certified Sobol tails; original N4096 model outputs are read-only.

Selection/manifest code uses only the standard library (RED Python 3.6+).
Numerical commands run in the original SALib 1.5.2 SIF and import the original
frozen model/sampling helpers. No command evaluates an old parameter sample.
"""
import argparse
import contextlib
import csv
import datetime
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

SOURCE_N = 4096
GLOBAL_NS = (8192, 16384, 32768)
EXPECTED_PAIRS = (
    ('p_mis_base', 'buffer_smax'), ('p_mis_base', 'buffer_beta'),
    ('buffer_smax', 'buffer_beta'), ('buffer_smax', 'O2_crit'),
    ('buffer_beta', 'mu_hp'), ('buffer_beta', 'O2_crit'))
OUTPUTS = ('dominant_mean_ploidy', 'dominant_growth_rate')


def now():
    return datetime.datetime.utcnow().isoformat() + 'Z'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for b in iter(lambda: handle.read(8 * 1024 * 1024), b''): h.update(b)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(path, obj):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    os.replace(str(tmp), str(path))


def read_table(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else Path(path).open()) as h:
        return list(csv.DictReader(h, delimiter='\t'))


def write_table(path, rows, fields):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    opener = gzip.open(tmp, 'wt') if str(path).endswith('.gz') else tmp.open('w')
    with opener as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    os.replace(str(tmp), str(path))


def selected_rows(rows, pairs):
    """Union failures across pairs AND outputs, without filtering by estimate."""
    result = {}
    for r in rows:
        pair = (r['parameter_i'], r['parameter_j'])
        if r['index'] != 'S2' or pair not in pairs: continue
        if r['convergence_status'] not in ('passed', 'not_passed', 'undefined'):
            raise ValueError('Unknown convergence status')
        if r['convergence_status'] == 'passed': continue
        oi = int(round(float(r['O2_pct']) / .025))
        if not 0 <= oi <= 200 or abs(oi * .025 - float(r['O2_pct'])) > 1e-10:
            raise ValueError('Unexpected oxygen point')
        result.setdefault(oi, []).append('__'.join(pair) + '/' + r['output'])
    return result


def select(source, output):
    source, output = Path(source), Path(output)
    if (output / 'selection.json').exists(): raise ValueError('Selection already frozen')
    p = load(source / 'slurm/submission_plan.json')
    if load(source / 'status.json')['status'] != 'complete' or p['N'] != SOURCE_N:
        raise ValueError('Require completed N4096 source')
    peaks = {}
    summary = source / 'neighborhood10pct/summaries/indices_scope_summary.tsv.gz'
    for r in read_table(summary):
        if r['index'] == 'S2':
            key = (r['parameter_i'], r['parameter_j'])
            peaks[key] = max(peaks.get(key, -float('inf')), float(r['value_mean']))
    pairs = {key for key, value in peaks.items() if value >= .01}
    if pairs != set(EXPECTED_PAIRS): raise ValueError('Approved six-pair selection changed')
    targets = []; cases = []; hashes = {}
    for case in p['cases']:
        directory = source / case['scope'] / 'full' / case['case_id']
        table = directory / 'indices_summary.tsv.gz'
        hashes[str(table.relative_to(source))] = sha(table)
        if case['scope'] == 'global': continue
        cells = selected_rows(read_table(table), pairs)
        if cells: cases.append(dict(case=case, oxygen_indices=sorted(cells)))
        for oi in sorted(cells):
            targets.append(dict(fit_seed=case['fit_seed'], oxygen_index=oi,
                                O2_pct=format(oi * .025, '.10g'),
                                failed_metrics=';'.join(sorted(cells[oi]))))
    if len(cases) != 497 or len(targets) != 55089:
        raise ValueError('Approved 497/55089 target population changed')
    write_table(output / 'targets.tsv.gz', targets,
                ['fit_seed', 'oxygen_index', 'O2_pct', 'failed_metrics'])
    save(output / 'selection.json', dict(status='frozen', created_utc=now(),
         source_plan_sha256=sha(source / 'slurm/submission_plan.json'),
         source_status_sha256=sha(source / 'status.json'),
         scope_summary_sha256=sha(summary), source_table_sha256=hashes,
         pairs=[list(x) for x in EXPECTED_PAIRS], cases=cases,
         n_seeds=497, n_seed_oxygen=55089, n_repeats=5,
         global_resolutions=list(GLOBAL_NS), neighborhood_resolution=8192,
         targets_sha256=sha(output / 'targets.tsv.gz'),
         selection_rule='scope median S2 peak >=0.01; any nonpassing key pair/output; union seed/O2',
         convergence_policy='diagnostic_only; original 500 endpoints retained'))
    print('Frozen 497 seeds / 55089 seed-oxygen targets', flush=True)


def task_layout(selection):
    cases = selection['cases']
    model_counts = {'g8192':150, 'g16384':300, 'g32768':600, 'n8192':len(cases)*150}
    analysis_counts = {'g8192':130, 'g16384':130, 'g32768':130, 'n8192':len(cases)*5}
    return model_counts, analysis_counts


def decode_model(scope, n, task, cases):
    if scope == 'global':
        segments = {8192:1, 16384:2, 32768:4}[n]
        if not 1 <= task <= 150 * segments: raise ValueError('Invalid model task')
        segment, within = divmod(task - 1, 150)
        repeat, block = divmod(within, 30)
        begin = n // 2 + segment * SOURCE_N
        return 'global', repeat + 1, block, begin, list(range(201))
    if n != 8192 or not 1 <= task <= len(cases)*150: raise ValueError('Invalid neighborhood task')
    ci, within = divmod(task - 1, 150); repeat, block = divmod(within, 30)
    row = cases[ci]
    return row['case']['case_id'], repeat + 1, block, SOURCE_N, row['oxygen_indices']


@contextlib.contextmanager
def locked(path):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as h:
        fcntl.flock(h, fcntl.LOCK_EX); yield


def runtime(root):
    global np, c, f
    import numpy as np
    p = load(Path(root) / 'plan.json'); old = Path(p['source_root'])
    if sha(Path(root)/'selection.json') != p['selection_sha256']:
        raise ValueError('Frozen target selection changed')
    if sha(Path(__file__)) != p['controller_sha256']['refinement_run.py']:
        raise ValueError('Frozen numerical controller changed')
    prepared=Path(root)/'preparation.json'
    if prepared.exists() and load(prepared)['plan_sha256']!=sha(Path(root)/'plan.json'):
        raise ValueError('Prepared execution plan changed')
    original = load(old / 'slurm/submission_plan.json')
    sys.path.insert(0, str(Path(original['code_snapshot_root']) /
                    'oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol'))
    import core as c
    import full_core as f
    c.configure(original)
    if sha(old / 'slurm/submission_plan.json') != p['source_plan_sha256']:
        raise ValueError('Source plan changed')
    if c.sha(f.HERE / 'evaluate.R') != original['engine_sha256']:
        raise ValueError('Frozen evaluator changed')
    sif = Path(original['sif']['path'])
    if sif.stat().st_size != original['sif']['bytes'] or sif.stat().st_mtime_ns != original['sif']['mtime_ns']:
        raise ValueError('SIF identity changed')
    return p, original, old


def case_for(original, cid):
    return next(x for x in original['cases'] if x['case_id'] == cid)


def scope_root(p, scope):
    return Path(p['global_root'] if scope == 'global' else p['neighborhood_root'])


def block_path(p, case, repeat, block, begin):
    return scope_root(p, case['scope']) / 'increments' / case['case_id'] / (
        'base%05d_%05d' % (begin+1, begin+SOURCE_N)) / ('R%02d' % repeat) / c.BLOCKS[block] / 'outputs.npz'


def unit_path(root, repeat, block, begin):
    return Path(root) / 'units' / ('R%02d_B%02d_from%05d.npz' % (repeat, block, begin))


def prepare(root):
    p, original, old = runtime(root); root = Path(root)
    if (root / 'preparation.json').exists():
        if load(root/'preparation.json')['plan_sha256'] != sha(root/'plan.json'):
            raise ValueError('Prepared plan changed')
        return
    if shutil.disk_usage(root).free < 650 * 1024**3:
        raise ValueError('Require 650 GiB free for new-only outputs, flags, checkpoints and reserve')
    subprocess.run([sys.executable,str(root/'slurm/code_snapshot/test_refinement.py')],check=True)
    if sha(Path(original['sif']['path'])) != c.SIF_SHA256: raise ValueError('SIF digest changed')
    # Audit every canonical R/model helper using the original immutable snapshot.
    for name, expected in original['code_sha256'].items():
        if name.endswith('.R') or '/model/' in name or name.endswith(('efast.py','neighborhood.py','perron_high_precision_safe.cpp')):
            if sha(Path(original['code_snapshot_root'])/name) != expected: raise ValueError('Frozen model/helper changed: '+name)
    unit_ranges = [dict(parameter=q, encoded_lower=0., encoded_upper=1.,
                       natural_lower=0., natural_upper=1., transform='identity') for q in c.ACTIVE]
    records = {}
    for repeat in range(1,6):
        prefix = old/'inputs'/('unit_R%02d.npz'%repeat)
        if sha(prefix) != original['unit_sha256'][str(repeat)]: raise ValueError('Old unit design changed')
        unit = c.samples(unit_ranges, 32768, 2026101300+repeat)
        with np.load(prefix, allow_pickle=False) as z:
            if not np.array_equal(unit[:SOURCE_N*c.STEP], z['unit']): raise ValueError('N4096 prefix differs')
        for begin in range(SOURCE_N, 32768, SOURCE_N):
            for bi in range(30):
                path = unit_path(root, repeat, bi, begin)
                data = unit[begin*c.STEP:(begin+SOURCE_N)*c.STEP][bi::c.STEP]
                if path.exists():
                    with np.load(path) as z:
                        if not np.array_equal(z['unit'],data): raise ValueError('Prepared tail changed')
                else: f.atomic_npz(path, unit=data)
                records[str(path.relative_to(root))] = sha(path)
        print('Exact nested prefix and seven new-only tails verified, repeat', repeat, flush=True)
    save(root/'preparation.json',dict(status='complete', plan_sha256=sha(root/'plan.json'),
         unit_sha256=records, n_verified_prefixes=5, new_sample_min=4097,
         model_evaluations_repeated=0, numerical_tests='passed in original SIF',
         free_bytes=shutil.disk_usage(root).free, completed_utc=now()))


def audit(root, scope, task):
    p, original, old = runtime(root); root=Path(root)
    case = original['cases'][0] if scope=='global' else original['cases'][task]
    cid=case['case_id']; dest=root/'audits'/(cid+'.json')
    with locked(str(dest)+'.lock'):
        if dest.exists():
            if load(dest)['plan_sha256'] != sha(root/'plan.json'): raise ValueError('Audit plan changed')
            return
        directory=f.case_dir(old,case); receipt=load(directory/'completion.json')
        if receipt['status']!='complete': raise ValueError('Incomplete original case')
        selection=load(root/'selection.json')
        table=directory/'indices_summary.tsv.gz'
        expected=selection['source_table_sha256'][str(table.relative_to(old))]
        if sha(table)!=expected or receipt['summary_table_sha256']!=expected:
            raise ValueError('Frozen source statistics changed')
        if scope!='global':
            cells=selected_rows(read_table(table),set(EXPECTED_PAIRS))
            target=next((x['oxygen_indices'] for x in selection['cases'] if x['case']['case_id']==cid),[])
            if sorted(cells)!=target: raise ValueError('Target mask differs from frozen diagnostics')
        bounds=c.read(directory/'parameter_bounds.tsv')
        if sha(directory/'parameter_bounds.tsv')!=case['bounds_sha256']: raise ValueError('Bounds changed')
        for name,field in [('best_params.tsv','best_params_sha256'),('fit_config.rds','fit_config_sha256')]:
            if sha(Path(original['fit_root'])/case['fit_seed']/name)!=case[field]: raise ValueError('Fit input changed')
        outputs=[]
        selected=(scope=='global' or any(x['case']['case_id']==cid for x in selection['cases']))
        if selected:
            expected_outputs={r['path']:r['sha256'] for r in receipt['output_receipts']}
            for repeat in range(1,6):
                with np.load(old/'inputs'/('unit_R%02d.npz'%repeat)) as z: unit=z['unit']
                for bi, block in enumerate(c.BLOCKS):
                    path=directory/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz'
                    key=f.identity(original,case,repeat,bi,f.scale_unit(unit[bi::30],bounds))
                    if not f.complete(path,key): raise ValueError('Old model block incomplete')
                    rr=load(str(path)+'.receipt.json')
                    if expected_outputs[str(path.relative_to(old))]!=rr['sha256']:
                        raise ValueError('Case/model receipt differs')
                    st=path.stat(); outputs.append(dict(path=str(path.relative_to(old)),sha256=rr['sha256'],bytes=st.st_size,mtime_ns=st.st_mtime_ns))
            for name,field in [('indices_repeats.npz','repeat_indices_sha256'),('indices_summary.npz','summary_sha256')]:
                if sha(directory/name)!=receipt[field]: raise ValueError('Old index cache changed')
        save(dest,dict(status='complete',case_id=cid,plan_sha256=sha(root/'plan.json'),
             outputs=outputs,source_completion_sha256=sha(directory/'completion.json'),
             source_indices_sha256=receipt['repeat_indices_sha256'],
             n_model_blocks_reused=len(outputs),completed_utc=now()))


def checked_old_outputs(root,p,old,case,repeat,oxygen_indices):
    audit=load(Path(root)/'audits'/(case['case_id']+'.json'))
    records={x['path']:x for x in audit['outputs']}; blocks=[]
    for block in c.BLOCKS:
        path=f.case_dir(old,case)/'runs'/('R%02d'%repeat)/'blocks'/block/'outputs.npz'
        record=records[str(path.relative_to(old))]; st=path.stat()
        if (st.st_size,st.st_mtime_ns)!=(record['bytes'],record['mtime_ns']): raise ValueError('Audited old block changed')
        with np.load(path,allow_pickle=False) as z:
            if z['outputs'].shape!=(4096,201,2) or not np.array_equal(z['oxygen'],c.OXYGEN): raise ValueError('Old model axes changed')
            blocks.append(z['outputs'][:,oxygen_indices,:])
    return c.assemble(blocks)


def evaluate(root,scope,n,task,max_chunks=None):
    p,original,old=runtime(root); root=Path(root)
    cid,repeat,bi,begin,ois=decode_model(scope,n,task,p['cases'])
    case=case_for(original,cid); audit_path=root/'audits'/(cid+'.json')
    if load(audit_path)['status']!='complete': raise ValueError('Case not audited')
    bounds=c.read(f.case_dir(old,case)/'parameter_bounds.tsv')
    unitfile=unit_path(root,repeat,bi,begin); prep=load(root/'preparation.json')
    if sha(unitfile)!=prep['unit_sha256'][str(unitfile.relative_to(root))]: raise ValueError('Tail unit changed')
    with np.load(unitfile) as z: x=f.scale_unit(z['unit'],bounds)
    path=block_path(p,case,repeat,bi,begin); directory=path.parent
    key=hashlib.sha256(json.dumps(dict(plan=sha(root/'plan.json'),case=case,repeat=repeat,
        block=bi,begin=begin,oxygen_indices=ois,unit=sha(unitfile)),sort_keys=True).encode()).hexdigest()
    with locked(directory/'worker.lock'):
        if f.complete(path,key): return
        oxygen=np.array(c.OXYGEN)[ois]; checkpoint=directory/'checkpoint.npz'
        values,flags,done,gaps,residuals,elapsed=f.checkpoint(checkpoint,key,4096,oxygen)
        scratch=Path(os.environ['SOBOL_SCRATCH'])/'evaluation'; scratch.mkdir(parents=True,exist_ok=True)
        c.table(scratch/'samples.tsv.gz',[dict(sample_id=i+1,**dict(zip(c.ACTIVE,map(float,row)))) for i,row in enumerate(x)])
        chunks=0
        while not done.all():
            if shutil.disk_usage(root).free<64*1024**3: raise ValueError('Storage reserve below 64 GiB; resume from checkpoint')
            positions=np.flatnonzero(~done)[:8]; current=oxygen[positions].tolist()
            meta=dict(method='Sobol_second_order',scope=case['scope'],case_id=cid,fit_seed=case['fit_seed'],
                n_samples=4096,parameters=list(c.ACTIVE),oxygen_pct=current,
                samples_sha256=sha(scratch/'samples.tsv.gz'),best_params_sha256=case['best_params_sha256'],
                fit_config_sha256=case['fit_config_sha256'],validate_reference=False,
                sample_id_offset=begin,source_prefix_reused=4096)
            save(scratch/'metadata.json',meta)
            subprocess.run(['Rscript','--vanilla',str(f.HERE/'evaluate.R'),str(scratch),original['fit_root'],
                str(old/'inputs/fixed_o2_dominant_ploidy_201grid.tsv')],check=True)
            v,s,g,r=f.read_model(scratch/'outputs.tsv.gz',4096,current)
            rr=load(str(scratch/'outputs.tsv.gz')+'.receipt.json')
            if rr['samples_sha256']!=meta['samples_sha256'] or rr['metadata_sha256']!=sha(scratch/'metadata.json') or rr['outputs_sha256']!=sha(scratch/'outputs.tsv.gz'):
                raise ValueError('New model provenance differs')
            values[:,positions]=v; flags[:,positions]=s; gaps[positions]=g; residuals[positions]=r
            done[positions]=True; elapsed+=rr['elapsed_seconds']
            if (scratch/'recovery').exists():
                shutil.copytree(scratch/'recovery',directory/'recovery',dirs_exist_ok=True); shutil.rmtree(scratch/'recovery')
            f.atomic_npz(checkpoint,identity=np.array(key),oxygen=oxygen,outputs=values,solver_flags=flags,
                done=done,minimum_gap=gaps,maximum_residual=residuals,elapsed_seconds=np.array(elapsed))
            print(cid,repeat,c.BLOCKS[bi],begin,'oxygen',int(done.sum()),'/',len(ois),flush=True)
            chunks+=1
            if max_chunks is not None and chunks>=max_chunks and not done.all():
                save(root/'smoke.json',dict(status='complete',checkpoint=str(checkpoint),
                     checkpoint_sha256=sha(checkpoint),certified_new_model_points=int(done.sum())*4096,
                     original_model_points_evaluated=0,production_will_resume_same_checkpoint=True,
                     completed_utc=now()))
                return
        f.atomic_npz(path,identity=np.array(key),oxygen=oxygen,oxygen_indices=np.array(ois),
            outputs=values,solver_flags=flags,minimum_gap=gaps,maximum_residual=residuals)
        save(str(path)+'.receipt.json',dict(status='complete',identity=key,sha256=sha(path),
            sample_begin=begin,sample_end=begin+4096,n_points=4096*len(ois),
            oxygen_indices=ois,elapsed_seconds=elapsed,n_recovered=int((flags==2).sum()),
            completed_utc=now(),job_id=os.environ.get('SLURM_JOB_ID')))
        checkpoint.unlink()


def analysis_task(p,scope,task):
    if scope=='global':
        if not 1<=task<=130: raise ValueError('Invalid global analysis task')
        rep,chunk=divmod(task-1,26); return 'global',rep+1,list(range(chunk*8,min(201,chunk*8+8)))
    if not 1<=task<=len(p['cases'])*5: raise ValueError('Invalid analysis task')
    ci,rep=divmod(task-1,5); row=p['cases'][ci]
    return row['case']['case_id'],rep+1,row['oxygen_indices']


def analysis_path(p,scope,n,cid,repeat,ois):
    return scope_root(p,scope)/('N%d'%n)/'shards'/cid/('R%02d_O%03d.npz'%(repeat,ois[0]))


def analyze(root,scope,n,task):
    p,original,old=runtime(root); root=Path(root)
    cid,repeat,ois=analysis_task(p,scope,task); case=case_for(original,cid)
    path=analysis_path(p,scope,n,cid,repeat,ois); key=sha(root/'plan.json')+'/%s/%d/%d'%(cid,n,repeat)+'/'+str(ois)
    with locked(str(path)+'.lock'):
        if f.complete(path,key): return
        if path.exists():
            with np.load(path,allow_pickle=False) as z:
                if str(z['identity'])!=key or not np.array_equal(z['oxygen_indices'],ois): raise ValueError('Analysis checkpoint changed')
                values,widths,variance,done=[z[k] for k in ('values','widths','variance','done')]
        else:
            shape=(2,len(ois),2,len(f.METRICS)); values=np.full(shape,np.nan); widths=np.full(shape,np.nan)
            variance=np.full(shape[:-1],np.nan); done=np.zeros((len(ois),2),bool)
            if n==8192:
                baseline=f.case_dir(old,case)/'indices_repeats.npz'
                if sha(baseline)!=load(root/'audits'/(cid+'.json'))['source_indices_sha256']:
                    raise ValueError('Audited N4096 statistical baseline changed')
                with np.load(baseline,allow_pickle=False) as z:
                    values[0]=z['values'][-1,repeat-1,ois]; widths[0]=z['bootstrap_half_width'][-1,repeat-1,ois]
                    variance[0]=z['variance'][-1,repeat-1,ois]
            else:
                lower=analysis_path(p,scope,n//2,cid,repeat,ois)
                if not f.complete(lower): raise ValueError('Lower resolution statistics missing')
                with np.load(lower) as z:
                    values[0]=z['values'][-1]; widths[0]=z['widths'][-1]; variance[0]=z['variance'][-1]
        # Load each model block once per repeat shard; never rerun the model here.
        ordered=[checked_old_outputs(root,p,old,case,repeat,ois)]
        for begin in range(4096,n,4096):
            blocks=[]
            for bi in range(30):
                source=block_path(p,case,repeat,bi,begin)
                if not f.complete(source): raise ValueError('New model block missing/corrupt: '+str(source))
                with np.load(source,allow_pickle=False) as z:
                    index={int(o):i for i,o in enumerate(z['oxygen_indices'])}
                    blocks.append(z['outputs'][:,[index[o] for o in ois],:])
            ordered.append(c.assemble(blocks)); del blocks
        ordered=np.concatenate(ordered,axis=0)
        if ordered.shape!=(n*30,len(ois),2) or not np.isfinite(ordered).all(): raise ValueError('Incomplete ordered design')
        bounds=c.read(f.case_dir(old,case)/'parameter_bounds.tsv')
        for j,oi in enumerate(ois):
            for k in range(2):
                if done[j,k]: continue
                si,var=c.analyze(bounds,ordered[:,j,k],seed=2026101300+repeat+10000+n+oi*10+k)
                values[1,j,k],widths[1,j,k]=f.metric_arrays(si); variance[1,j,k]=var; done[j,k]=True
                f.atomic_npz(path,identity=np.array(key),oxygen_indices=np.array(ois),values=values,widths=widths,variance=variance,done=done)
                print(cid,n,repeat,oi,k,'index cell certified',flush=True)
        save(str(path)+'.receipt.json',dict(status='complete',identity=key,sha256=sha(path),
            completed_utc=now(),model_evaluations_repeated=0))


def plot_summary(mean,oxygen,directory,title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap,TwoSlopeNorm
    directory.mkdir(parents=True,exist_ok=True)
    for index in ('S1','ST','S2'):
        chosen=[i for i,m in enumerate(f.METRICS) if m[0]==index]
        labels=[f.METRICS[i][1]+(' × '+f.METRICS[i][2] if index=='S2' else ' ['+c.efast.GROUP[f.METRICS[i][1]]+']') for i in chosen]
        fig,axes=plt.subplots(1,2,figsize=(16,23 if index=='S2' else 6),layout='constrained')
        for k,ax in enumerate(axes):
            data=mean[:,k,chosen].T; finite=data[np.isfinite(data)]
            lo=min(0.,float(finite.min())) if len(finite) else 0.; hi=max(.001,float(finite.max())) if len(finite) else .001
            color='#e66101' if index=='ST' else '#54278f'
            cmap=LinearSegmentedColormap.from_list('sobol',['#2166ac','white',color] if lo<0 else ['white',color]); cmap.set_bad('#d9d9d9')
            norm=TwoSlopeNorm(vmin=lo,vcenter=0,vmax=hi) if lo<0 else None
            im=ax.imshow(np.ma.masked_invalid(data),aspect='auto',cmap=cmap,norm=norm,**({} if norm else dict(vmin=0,vmax=hi)))
            ax.set_xticks(range(0,len(oxygen),20),[format(oxygen[i],'.2g') for i in range(0,len(oxygen),20)])
            ax.set_yticks(range(len(labels)),labels if k==0 else ['']*len(labels),fontsize=7)
            ax.set_xlabel('Oxygen (%)'); ax.set_title(OUTPUTS[k])
            bar=fig.colorbar(im,ax=ax,orientation='horizontal',pad=.025,fraction=.025,label=index+'; raw estimates, blue = negative sampling estimates')
            ticks=[lo,0,hi] if lo<0 else [0,hi/2,hi]; bar.set_ticks(ticks); bar.set_ticklabels([format(x,'.3g') for x in ticks])
        fig.suptitle(title+'\nAll estimates retained; convergence is diagnostic only',fontsize=11)
        name='s2_pair_by_oxygen' if index=='S2' else index.lower()+'_parameter_by_oxygen'
        fig.savefig(directory/(name+'.pdf')); fig.savefig(directory/(name+'.png'),dpi=120); plt.close(fig)


def comparison_curves(scopeout,scope,series,coverage=None):
    """Compare the same population at each oxygen; never mix target eligibility."""
    import matplotlib.pyplot as plt
    pairs=[('lam_max','p_misseg'),('lam_max','p_mis_base'),('lam_max','O2_crit'),
           ('lam_max','mu_hp'),('p_misseg','mu_hp')] if scope=='global' else EXPECTED_PAIRS
    directory=scopeout/'resolution_comparison'; directory.mkdir(parents=True,exist_ok=True)
    records=[]
    for kind,items in [('S2',pairs),('ST',[(q,'') for q in c.ACTIVE])]:
        fig,axes=plt.subplots(len(items),2,figsize=(14,3*len(items)),layout='constrained',squeeze=False)
        for row,(pi,pj) in enumerate(items):
            m=f.METRICS.index((kind,pi,pj))
            for k,out in enumerate(OUTPUTS):
                ax=axes[row,k]
                for label,data in series:
                    ax.plot(c.OXYGEN,data[:,k,m],label=label,lw=1)
                    for oi,o in enumerate(c.OXYGEN):
                        records.append(dict(index=kind,parameter_i=pi,parameter_j=pj,output=out,
                            resolution=label,O2_pct=o,value=data[oi,k,m],
                            n_matched_endpoints=1 if coverage is None else int(coverage[oi])))
                ax.axhline(0,color='gray',lw=.5); ax.set_xlabel('Oxygen (%)'); ax.set_ylabel(kind)
                ax.set_title(pi+(' × '+pj if pj else '')+'; '+out)
                if row==0: ax.legend(fontsize=8)
        fig.suptitle('Global resolution comparison' if scope=='global' else 'Matched targeted endpoints only; coverage varies with oxygen')
        fig.savefig(directory/(kind.lower()+'_curves.pdf')); fig.savefig(directory/(kind.lower()+'_curves.png'),dpi=110); plt.close(fig)
    write_table(directory/'curves_source.tsv.gz',records,list(records[0]))


def report(root,scope):
    p,original,old=runtime(root); root=Path(root); scopeout=scope_root(p,scope)
    stages=GLOBAL_NS if scope=='global' else (8192,)
    rows_diag=[]; series=[]
    if scope=='global':
        with np.load(f.case_dir(old,original['cases'][0])/'indices_summary.npz') as z:
            series.append(('N4096',z['mean']))
    for n in stages:
        case_rows=[dict(case=original['cases'][0],oxygen_indices=list(range(201)))] if scope=='global' else p['cases']
        cube=[]; masks=[]; statuses=[]
        if scope!='global':
            for case in original['cases'][1:]:
                with np.load(f.case_dir(old,case)/'indices_summary.npz') as z:
                    cube.append(z['mean']); statuses.append(z['status'])
            cube=np.stack(cube); baseline=cube.copy(); statuses=np.stack(statuses); used=np.full((500,201),4096,int)
        for row in case_rows:
            case=row['case']; cid=case['case_id']; ois=row['oxygen_indices']; vv=[]; ww=[]; va=[]
            for repeat in range(1,6):
                parts=[ois[i:i+8] for i in range(0,len(ois),8)] if scope=='global' else [ois]
                vparts=[]; wparts=[]; aparts=[]
                for part in parts:
                    source=analysis_path(p,scope,n,cid,repeat,part)
                    if not f.complete(source): raise ValueError('Statistical shard incomplete')
                    with np.load(source) as z:
                        if not z['done'].all(): raise ValueError('Partial statistical shard')
                        vparts.append(z['values']); wparts.append(z['widths']); aparts.append(z['variance'])
                vv.append(np.concatenate(vparts,axis=1)); ww.append(np.concatenate(wparts,axis=1)); va.append(np.concatenate(aparts,axis=1))
            values=np.stack(vv,axis=1); widths=np.stack(ww,axis=1); variance=np.stack(va,axis=1)
            mean,spread,delta,status=f.diagnostics(values)
            target=scopeout/('N%d'%n)/'cases'/cid; target.mkdir(parents=True,exist_ok=True)
            f.atomic_npz(target/'indices_summary.npz',mean=mean,repeat_range=spread,resolution_delta=delta,status=status,oxygen_indices=np.array(ois))
            def records():
                for j,oi in enumerate(ois):
                    for k,out in enumerate(OUTPUTS):
                        for m,(idx,pi,pj) in enumerate(f.METRICS):
                            yield dict(scope=scope,fit_seed=case['fit_seed'],N=n,O2_pct=c.OXYGEN[oi],output=out,index=idx,parameter_i=pi,parameter_j=pj,
                                value_mean=mean[j,k,m],repeat_range=spread[j,k,m],resolution_delta=delta[j,k,m],resolution_low=n//2,
                                bootstrap_half_width_median=np.median(widths[1,:,j,k,m]),convergence_status=['undefined','passed','not_passed'][status[j,k,m]])
            records_iter=records(); first=next(records_iter)
            from itertools import chain
            write_table(target/'indices_summary.tsv.gz',chain([first],records_iter),list(first))
            rows_diag.append(dict(scope=scope,case_id=cid,N=n,n_oxygen=len(ois),n_passed=int((status==1).sum()),n_not_passed=int((status==2).sum()),n_undefined=int((status==0).sum())))
            if scope=='global':
                cube=mean[None]; statuses=status[None]
            else:
                si=int(cid[4:])-1; cube[si,ois]=mean; statuses[si,ois]=status; used[si,ois]=n
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            median=np.nanmedian(cube,axis=0); q10=np.nanquantile(cube,.1,axis=0); q90=np.nanquantile(cube,.9,axis=0)
        target=scopeout/('N%d'%n); title='Global N%d; mean of five repeats'%n if scope=='global' else '500-endpoint median; MIXED N8192 targeted / N4096 elsewhere'
        plot_summary(median,c.OXYGEN,target/'figures',title)
        if scope=='global': series.append(('N%d'%n,median))
        def scope_records():
            for oi,o in enumerate(c.OXYGEN):
                for k,out in enumerate(OUTPUTS):
                    for m,(idx,pi,pj) in enumerate(f.METRICS):
                        yield dict(scope=scope,N=str(n) if scope=='global' else 'mixed_4096_8192',O2_pct=o,output=out,index=idx,parameter_i=pi,parameter_j=pj,
                            value_mean=median[oi,k,m],endpoint_q10=q10[oi,k,m],endpoint_q90=q90[oi,k,m],n_endpoints=len(cube),
                            n_defined=int(np.isfinite(cube[:,oi,k,m]).sum()),n_passed=int((statuses[:,oi,k,m]==1).sum()),
                            n_at_new_N=1 if scope=='global' else int((used[:,oi]==n).sum()),convergence_exclusion='none')
        rows=scope_records(); first=next(rows); write_table(target/'summaries/indices_scope_summary.tsv.gz',chain([first],rows),list(first))
        if scope!='global':
            write_table(target/'N_used.tsv.gz',(dict(fit_seed='seed%d'%(si+1),O2_pct=o,N_used=int(used[si,oi])) for si in range(500) for oi,o in enumerate(c.OXYGEN)),['fit_seed','O2_pct','N_used'])
            # Matched targets only: unchanged N4096 cells cannot enter an N8192 median.
            targeted=np.where((used==8192)[:,:,None,None],cube,np.nan)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning); only=np.nanmedian(targeted,axis=0)
                old_only=np.nanmedian(np.where((used==8192)[:,:,None,None],baseline,np.nan),axis=0)
            plot_summary(only,c.OXYGEN,target/'targeted_only_figures','Targeted N8192 only; endpoint coverage varies with oxygen')
            series=[('N4096 matched targets',old_only),('N8192 matched targets',only)]
            coverage=(used==8192).sum(axis=0)
            comparison_curves(scopeout,scope,series,coverage)
        save(target/'completion.json',dict(status='complete',scope=scope,N=n if scope=='global' else 'mixed_4096_8192',n_repeats=5,
            original_endpoints_retained=1 if scope=='global' else 500,convergence_policy='diagnostic_only',completed_utc=now()))
    write_table(scopeout/'diagnostics.tsv',rows_diag,list(rows_diag[0]))
    if scope=='global': comparison_curves(scopeout,scope,series)
    (scopeout/'analysis_notes.md').write_text(
        '# Incremental Sobol refinement\n\n'
        'Original N4096 model evaluations are reused without recalculation. Five repeats are retained.\n'
        'All estimates, including diagnostic failures and undefined indices, remain in the output.\n'
        'Repeat range <=0.05 and adjacent-resolution mean delta <=0.02 are operational checks.\n'
        'Negative estimates are sampling error, not direction or antagonism. Bootstrap half widths\n'
        'are retained in case source tables. Key curves compare identical target populations.\n'
        'Neighborhood mixed summaries retain all 500 original endpoints; N_used identifies each\n'
        'cell resolution. Targeted-only figures have varying endpoint coverage and cannot be\n'
        'described as uniform N8192 sampling of all 500 neighborhoods. Existing conditional\n'
        'response surfaces are unchanged and remain conditional illustrations, not marginalized S2.\n')
    save(scopeout/'status.json',dict(status='complete',resolutions=list(stages),completed_utc=now(),original_N4096_read_only=True))


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('mode',choices=['select','prepare','audit','smoke','model','analyze','report'])
    a.add_argument('--root'); a.add_argument('--source'); a.add_argument('--out')
    a.add_argument('--scope',choices=['global','neighborhood10pct']); a.add_argument('--n',type=int)
    a.add_argument('--task',type=int)
    args=a.parse_args()
    if args.mode=='select': select(args.source,args.out); return
    try:
        if args.mode=='prepare': prepare(args.root)
        elif args.mode=='audit': audit(args.root,args.scope,args.task)
        elif args.mode=='smoke': evaluate(args.root,'global',8192,1,max_chunks=1)
        elif args.mode=='model': evaluate(args.root,args.scope,args.n,args.task)
        elif args.mode=='analyze': analyze(args.root,args.scope,args.n,args.task)
        else: report(args.root,args.scope)
    except Exception as exc:
        save(Path(args.root)/'failures'/('%s_%s_%s_%s.json'%(args.mode,args.scope,args.n,args.task)),
             dict(error=str(exc),recorded_utc=now(),job_id=os.environ.get('SLURM_JOB_ID')))
        raise


if __name__=='__main__': main()
