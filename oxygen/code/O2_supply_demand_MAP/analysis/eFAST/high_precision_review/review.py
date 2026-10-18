#!/usr/bin/env python3
"""Frozen point selection and independent decimal audit. Never resample FAST."""
import argparse
import collections
import csv
import datetime
from decimal import Decimal,localcontext
import gzip
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
PARAMETERS=('lam_max','p_mis_base','p_misseg','k_o_mis','buffer_smax','buffer_beta','buffer_n_exp','p_wgd','alpha_o2','gamma_growth','mu_hp','gamma_mu','O2_crit','n_O')


def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(8*1024**2),b''):h.update(block)
    return h.hexdigest()


def read(p):
    with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def table(p,rows):
    rows=list(rows);assert rows
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


def dump(p,value):
    p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def prepare(root):
    import numpy as np
    parent=root.parent;repo=HERE.parents[5]
    assert not (root/'selection.json').exists(),'Refuse to replace a frozen selection'
    root.mkdir(parents=True,exist_ok=True)
    rows=read(parent/'influential_response_points.tsv');trajectories={r['trajectory_id']:r for r in read(parent/'trajectory_manifest.tsv')}
    chosen=[r for r in rows if (r['fit_seed'] in ('seed138','seed101') and 6.5<=float(r['ploidy'])<=7.00000001)
            or (r['fit_seed']=='seed174' and r['M']=='8' and r['replicate']=='5' and float(r['ploidy'])>1.01)]
    assert len(chosen)==106
    selected={};relations=[]
    def add(tid,sid,j,role):
        key=(tid,int(sid),int(j));selected.setdefault(key,set()).add(role);return key
    for r in chosen:
        tid=r['trajectory_id'];sid=int(r['sample_id']);j=int(round(float(r['O2_pct'])/.025))
        with np.load(parent/'processed_cache'/('trajectory_%03d.npz'%int(tid))) as z:
            ci=int(np.argmin(abs(z['case_oxygen']-j*.025)));assert abs(z['case_oxygen'][ci]-j*.025)<1e-10
            y=z['case_y'][:,ci,0];u=z['u'];assert abs(y[sid-1]-float(r['ploidy']))<1e-11
            lo,hi=np.quantile(y,[.4,.6]);ordinary=(y>=lo)&(y<=hi)&(y<2)
            if r['fit_seed']=='seed174':ordinary=(abs(y-1)<.01)
            assert ordinary.any()
            distance=np.sum((u-u[sid-1])**2,axis=1);distance[~ordinary]=np.inf;control=int(np.argmin(distance))+1
        abnormal=add(tid,sid,j,'abnormal');normal=add(tid,control,j,'ordinary_control')
        relations.append(dict(abnormal_key=':'.join(map(str,abnormal)),related_key=':'.join(map(str,normal)),relationship='ordinary_control',normalized_parameter_distance=float(np.sqrt(distance[control-1]))))
        for neighbor in (j-1,j+1):
            if 0<=neighbor<=200:
                key=add(tid,sid,neighbor,'oxygen_neighbor')
                relations.append(dict(abnormal_key=':'.join(map(str,abnormal)),related_key=':'.join(map(str,key)),relationship='oxygen_neighbor',normalized_parameter_distance=0.))
    points=[];sources={};plans={};fit_proofs={}
    for tid in sorted({key[0] for key in selected},key=int):
        tr=trajectories[tid];path=parent/'raw_cache'/tr['part_path'];n=int(tr['N'])
        for name,field in [('samples.tsv.gz','samples_sha256'),('outputs.tsv.gz','output_sha256'),('metadata.json','metadata_sha256')]:
            assert sha(path/name)==tr[field]
            sources[str(Path(tr['part_path'])/name)]=dict(path=str(Path(tr['part_path'])/name),sha256=tr[field])
        samples=np.loadtxt(gzip.open(path/'samples.tsv.gz','rt'),skiprows=1,delimiter='\t')
        data=np.loadtxt(gzip.open(path/'outputs.tsv.gz','rt'),skiprows=1,delimiter='\t',usecols=(0,1,2,3,4)).reshape(n,201,5)
        np.testing.assert_array_equal(data[:,:,0],np.broadcast_to(np.arange(1,n+1)[:,None],(n,201)))
        np.testing.assert_allclose(data[:,:,1],np.broadcast_to(np.arange(201)*.025,(n,201)),atol=1e-12,rtol=0)
        for key in sorted(k for k in selected if k[0]==tid):
            _,sid,j=key;values=data[sid-1,j];assert int(samples[sid-1,0])==sid
            rec=dict(point_id=len(points)+1,point_key=':'.join(map(str,key)),trajectory_id=tid,fit_seed=tr['fit_seed'],N=tr['N'],M=tr['M'],
                replicate=tr['replicate'],parameter=tr['parameter'],sample_id=sid,O2_index=j,O2_pct=format(j*.025,'.17g'),
                roles=','.join(sorted(selected[key])),original_ploidy=format(values[2],'.17g'),original_growth=format(values[3],'.17g'),original_gap=format(values[4],'.17g'),
                source_part=tr['part_path'],samples_sha256=tr['samples_sha256'],outputs_sha256=tr['output_sha256'],metadata_sha256=tr['metadata_sha256'])
            rec.update({name:format(samples[sid-1,k+1],'.17g') for k,name in enumerate(PARAMETERS)})
            points.append(rec)
        print('Selected trajectory',tid,tr['fit_seed'],'points',sum(k[0]==tid for k in selected),flush=True)
    for config in ('N8193_M4','N4097_M8'):
        plan_path=parent.parent/config/'slurm/submission_plan.json';plan=json.loads(plan_path.read_text())
        plans[config]=dict(path=str(plan_path.relative_to(repo)),sha256=sha(plan_path),
                           fit_root=plan['fit_root'],sif=plan['sif'],sif_sha256=plan['sif_sha256'],archives=plan['archives'],code_sha256=plan['code_sha256'])
        for r in read(parent.parent/config/'seed_selection.tsv'):
            if r['fit_seed'] in ('seed138','seed101','seed174'):
                fit_proofs[r['fit_seed']]={k:v for k,v in r.items() if k.endswith('_sha256')}
    table(root/'point_manifest.tsv',points);table(root/'relationships.tsv',relations);table(root/'source_files.tsv',sources.values())
    bykey={p['point_key']:p for p in points};smoke=[]
    for seed,parameter in [('seed138','p_mis_base'),('seed101','p_mis_base'),('seed174','O2_crit')]:
        candidate=next(p for p in points if p['fit_seed']==seed and p['parameter']==parameter and 'abnormal' in p['roles']
                       and (p['N']=='8193' if seed=='seed138' else p['M']=='8'))
        control=next(r['related_key'] for r in relations if r['abnormal_key']==candidate['point_key'] and r['relationship']=='ordinary_control')
        smoke.extend((candidate['point_id'],bykey[control]['point_id']))
    dump(root/'smoke_points.json',smoke)
    dump(root/'selection.json',dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        n_points=len(points),n_abnormal_candidates=len(chosen),abnormal_counts=dict(collections.Counter(r['fit_seed'] for r in chosen)),
        role_counts=dict(collections.Counter(role for roles in selected.values() for role in roles)),
        selection_scope='All qualifying top-three influential points at previously disclosed representative O2 cases; not a census across all O2',
        near_seven_threshold=6.5,seed174_scope='M8 phase 5, all five focused parameters, response >1.01',
        ordinary_control='Closest normalized 14D vector among middle 40-60% ploidy <2; seed174 |ploidy-1|<.01',
        oxygen_neighbors='Same exact parameter vector at existing O2 grid index +/-1',
        point_manifest_sha256=sha(root/'point_manifest.tsv'),relationships_sha256=sha(root/'relationships.tsv'),
        source_figure_table_sha256=sha(parent/'influential_response_points.tsv'),source_trajectory_manifest_sha256=sha(parent/'trajectory_manifest.tsv'),
        source_plans=plans,fit_proofs=fit_proofs,code_sha256={p.name:sha(p) for p in sorted(HERE.iterdir()) if p.is_file() and p.suffix in ('.py','.R','.cpp','.sbatch')},
        qos='xxlarge',time='12:00:00',specified_node=None,array_concurrency_limit=None,new_FAST_samples=0))
    print('Frozen',len(points),'points',flush=True)


def audit_point(path):
    report=json.loads((path/'result.json').read_text());assert report['status']=='evaluated'
    for name,digest in report['files_sha256'].items():assert sha(path/name)==digest,name
    result={k:report['identity'][k] for k in ('point_id','fit_seed','N','M','replicate','parameter','sample_id','O2_pct','roles')}
    result.update(original_ploidy=float(report['identity']['original_ploidy']),original_growth=float(report['identity']['original_growth']),
        reproduced_ploidy=report['double']['normalized_ploidy'],reproduced_growth=report['double']['lambda_real'],
        raw_negative_mass=report['double']['negative_l1_mass'],normalization_ploidy_change=report['double']['clipped_ploidy']-report['double']['signed_ploidy'],
        raw_relative_residual=report['double']['raw_relative_residual'],normalized_relative_residual=report['double']['normalized_relative_residual'],
        strongly_connected=report['matrix']['strongly_connected'],matrix_sha256=report['files_sha256']['matrix.rds'])
    all_solved=True
    with localcontext() as ctx:
        ctx.prec=140
        solved={}
        for digits in (50,100):
            filename=path/('precision_%d.json'%digits);p=json.loads(filename.read_text());solved[digits]=p
            for side in ('right','left'):
                v=p[side]
                if v.get('error') or not v.get('converged',False):all_solved=False;continue
                x=[Decimal(a) for a in v['vector_decimal']];assert len(x)==report['matrix']['n'] and min(x)>0
                assert abs(sum(x)-1)<Decimal('1e-35' if digits==50 else '1e-75')
                assert Decimal(v['relative_residual_inf'])<Decimal(v['relative_tolerance'])
                assert abs(Decimal(v['relative_upper_minus_lambda']))<Decimal(v['relative_tolerance'])
                assert Decimal(v['vector_change_l1'])<Decimal(v['vector_tolerance'])
                result['residual_%s_%d'%(side,digits)]=v['relative_residual_inf']
        if all_solved:
            a,b=solved[50],solved[100]
            for side in ('right','left'):
                l1=sum(abs(Decimal(x)-Decimal(y)) for x,y in zip(a[side]['vector_decimal'],b[side]['vector_decimal']))
                result[side+'_precision_L1_difference']=str(l1)
                if l1>=Decimal('1e-10'):all_solved=False
            delta=abs(Decimal(a['right']['ploidy'])-Decimal(b['right']['ploidy']))
            lam_delta=abs(Decimal(a['right']['lambda'])-Decimal(b['right']['lambda']))
            lr=abs(Decimal(b['right']['lambda'])-Decimal(b['left']['lambda']))
            result.update(ploidy_50=a['right']['ploidy'],ploidy_100=b['right']['ploidy'],growth_100=b['right']['lambda'],
                ploidy_precision_difference=str(delta),growth_precision_difference=str(lam_delta),left_right_eigenvalue_difference=str(lr))
            x=[Decimal(s) for s in b['right']['vector_decimal']];w=[Decimal(s) for s in b['left']['vector_decimal']]
            product=sum(u*v for u,v in zip(x,w));result['left_right_overlap']=str(product)
            result['eigenvalue_condition_2norm']=str((sum(u*u for u in x)*sum(v*v for v in w)).sqrt()/abs(product))
            result['ploidy_change']=float(Decimal(b['right']['ploidy'])-Decimal(report['identity']['original_ploidy']))
            result['growth_change']=float(Decimal(b['right']['lambda'])-Decimal(report['identity']['original_growth']))
            if delta>=Decimal('1e-10') or lam_delta>=Decimal('1e-12') or lr>=Decimal('1e-12'):all_solved=False
    result['precision_verified']=all_solved
    result['classification']='unresolved'
    if all_solved:
        if abs(result['ploidy_change'])>1e-6:result['classification']='double_response_changed'
        elif result['original_ploidy']>=6.5:result['classification']='high_ploidy_retained'
        else:result['classification']='response_retained'
    # A reducible matrix requires separate uniqueness adjudication.
    result['unique_vector_certified']=all_solved and result['strongly_connected']
    return result


def collect(root,smoke=False):
    selection=json.loads((root/'selection.json').read_text());assert sha(root/'point_manifest.tsv')==selection['point_manifest_sha256']
    rows=read(root/'point_manifest.tsv');ids=json.loads((root/'smoke_points.json').read_text()) if smoke else [int(r['point_id']) for r in rows]
    results=[];missing=[]
    for i in ids:
        path=root/'points'/('point_%04d'%i)
        if not (path/'result.json').exists():missing.append(i);continue
        results.append(audit_point(path))
    if missing:raise RuntimeError('Missing point results: '+str(missing))
    fields=list(dict.fromkeys(k for r in results for k in r))
    table(root/('smoke_summary.tsv' if smoke else 'point_summary.tsv'),[{k:r.get(k,'') for k in fields} for r in results])
    record=dict(status='passed' if smoke else 'complete',n_points=len(results),n_precision_verified=sum(r['precision_verified'] for r in results),
        classifications=dict(collections.Counter(r['classification'] for r in results)),selection_sha256=sha(root/'selection.json'),
        point_manifest_sha256=sha(root/'point_manifest.tsv'),new_FAST_samples=0,
        note='Workflow completeness is separate from per-point precision convergence and original-double agreement')
    if smoke:
        assert all(r['precision_verified'] for r in results),'Smoke solver did not certify all demonstration cases'
        fixtures=json.loads((root/'smoke_fixtures.json').read_text());assert fixtures['status']=='passed'
    dump(root/('smoke_receipt.json' if smoke else 'completion.json'),record);print(json.dumps(record,indent=2))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('prepare','collect','smoke-audit'))
    p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root.resolve()
    if a.command=='prepare':prepare(root)
    else:collect(root,a.command=='smoke-audit')


if __name__=='__main__':main()
