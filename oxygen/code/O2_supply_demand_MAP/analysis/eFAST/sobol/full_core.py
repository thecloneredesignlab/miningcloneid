"""Full 501-domain Sobol contracts, compact outputs and resumable checkpoints."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import os

import numpy as np
import core as c

N=4096
NS=[2048,4096]
N_CASES=501
N_TASKS=N_CASES*5*c.STEP
METRICS=[('S1',p,'') for p in c.ACTIVE]+[('ST',p,'') for p in c.ACTIVE]+[
    ('S2',c.ACTIVE[i],c.ACTIVE[j]) for i,j in c.PAIRS]
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[5]


def decode(task_id):
    if not 1<=task_id<=N_TASKS: raise ValueError('Task ID outside full manifest')
    case_index,within=divmod(task_id-1,5*c.STEP)
    repeat,block=divmod(within,c.STEP)
    return case_index+1,repeat+1,block


def case_dir(root,case):
    return Path(root)/case['scope']/'full'/case['case_id']


def plan(root):
    p=c.load_json(Path(root)/'slurm/submission_plan.json')
    engine=HERE/'evaluate.R'
    if c.sha(engine)!=p['engine_sha256']: raise ValueError('Canonical engine changed')
    if str(REPO)!=p['code_snapshot_root']: raise ValueError('Use the frozen code snapshot')
    sif=Path(p['sif']['path'])
    if sif.stat().st_size!=p['sif']['bytes'] or sif.stat().st_mtime_ns!=p['sif']['mtime_ns']:
        raise ValueError('Verified runtime changed')
    c.configure(p)
    return p


def samples(root,case,repeat,block=None):
    bounds=c.read(case_dir(root,case)/'parameter_bounds.tsv')
    p=c.load_json(Path(root)/'slurm/submission_plan.json')
    unit_file=Path(root)/'inputs'/('unit_R%02d.npz'%repeat)
    if c.sha(unit_file)!=p['unit_sha256'][str(repeat)]: raise ValueError('Unit design changed')
    with np.load(unit_file) as z: x=z['unit']
    if block is not None: x=x[block::c.STEP]
    return scale_unit(x,bounds)


def scale_unit(unit,bounds):
    lo=np.array([float(r['encoded_lower']) for r in bounds])
    hi=np.array([float(r['encoded_upper']) for r in bounds])
    x=unit.copy(); np.multiply(x,hi-lo,out=x); np.add(x,lo,out=x)
    for j,r in enumerate(bounds):
        if r['transform']=='log10': x[:,j]=10.**x[:,j]
    return x


def identity(p,case,repeat,block,x):
    obj=dict(code_git_commit=p['git_commit'],engine=p['engine_sha256'],sif=p['sif']['sha256'],
             scope=case['scope'],case_id=case['case_id'],repeat=repeat,block=block,
             samples_sha256=hashlib.sha256(np.asarray(x,dtype='<f8').tobytes()).hexdigest(),
             bounds_sha256=case['bounds_sha256'],best_params_sha256=case['best_params_sha256'],
             fit_config_sha256=case['fit_config_sha256'],oxygen=p['oxygen_pct'])
    return hashlib.sha256(json.dumps(obj,sort_keys=True).encode()).hexdigest()


def atomic_npz(path,**values):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.incomplete')
    with open(tmp,'wb') as handle:
        np.savez_compressed(handle,**values); handle.flush(); os.fsync(handle.fileno())
    tmp.replace(path)


def checkpoint(path,key,n,oxygen):
    path=Path(path)
    if not path.exists():
        return (np.full((n,len(oxygen),2),np.nan),np.zeros((n,len(oxygen)),np.uint8),
                np.zeros(len(oxygen),bool),np.full(len(oxygen),np.nan),
                np.zeros(len(oxygen)),0.)
    with np.load(path) as z:
        if str(z['identity'])!=key or not np.array_equal(z['oxygen'],oxygen):
            raise ValueError('Checkpoint identity changed')
        values,flags,done=z['outputs'],z['solver_flags'],z['done']
        gaps,residuals=z['minimum_gap'],z['maximum_residual']
        elapsed=float(z['elapsed_seconds'])
    if values.shape!=(n,len(oxygen),2) or flags.shape!=(n,len(oxygen)) or done.shape!=(len(oxygen),):
        raise ValueError('Checkpoint axes differ')
    if not np.isfinite(values[:,done]).all() or not np.isin(flags[:,done],[1,2]).all():
        raise ValueError('Invalid completed checkpoint points')
    if np.any(flags[:,~done]) or np.isfinite(values[:,~done]).any():
        raise ValueError('Uncertified checkpoint points')
    return values,flags,done,gaps,residuals,elapsed


def read_model(path,n,oxygen):
    values=np.full((n,len(oxygen),2),np.nan)
    flags=np.zeros((n,len(oxygen)),np.uint8)
    gaps=np.full(len(oxygen),np.inf); residuals=np.zeros(len(oxygen))
    lookup={round(float(o),10):j for j,o in enumerate(oxygen)}
    with gzip.open(path,'rt',newline='') as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            i=int(row['sample_id'])-1; j=lookup.get(round(float(row['O2_pct']),10))
            if not 0<=i<n or j is None or flags[i,j]: raise ValueError('Duplicate/unexpected model point')
            if row['status']!='ok' or row['eigenvector_nonnegative']!='TRUE':
                raise ValueError('Invalid model point; never delete rows')
            values[i,j]=[float(row[k]) for k in c.OUTPUTS]
            if not np.isfinite(values[i,j]).all(): raise ValueError('Nonfinite model output')
            if row['solver'] not in ('double','decimal50_100'): raise ValueError('Unknown solver')
            flags[i,j]=1 if row['solver']=='double' else 2
            gap=float(row['spectral_gap']); residual=float(row['eigen_residual'])
            if not np.isfinite([gap,residual]).all(): raise ValueError('Nonfinite eigenpair diagnostics')
            gaps[j]=min(gaps[j],gap); residuals[j]=max(residuals[j],residual)
    if not flags.all(): raise ValueError('Incomplete model rectangle')
    return values,flags,gaps,residuals


def complete(path,key=None):
    path=Path(path); receipt=Path(str(path)+'.receipt.json')
    if not path.exists() or not receipt.exists(): return False
    r=c.load_json(receipt)
    if r.get('status')!='complete' or (key is not None and r['identity']!=key):
        raise ValueError('Completion identity differs')
    if c.sha(path)!=r['sha256']: raise ValueError('Completed output checksum differs')
    return True


def metric_arrays(si):
    values=[]; widths=[]
    for index,pi,pj in METRICS:
        i=c.ACTIVE.index(pi)
        if index=='S2':
            j=c.ACTIVE.index(pj); values.append(si[index][i,j]); widths.append(si[index+'_conf'][i,j])
        else: values.append(si[index][i]); widths.append(si[index+'_conf'][i])
    return np.asarray(values),np.asarray(widths)


def diagnostics(values):
    """values = N2048/N4096 x five repeats x oxygen x output x metric."""
    valid=np.isfinite(values[-1]).all(axis=0)
    mean=np.mean(values,axis=1)
    spread=np.ptp(values[-1],axis=0)
    delta=np.abs(mean[-1]-mean[-2])
    status=np.where(valid&np.isfinite(delta),np.where((spread<=.05)&(delta<=.02),1,2),0).astype(np.uint8)
    return mean[-1],spread,delta,status
