#!/usr/bin/env python3
"""Exercise cross-batch reports with explicitly synthetic indices; no model/estimator calls."""
import argparse
from pathlib import Path
import shutil
import tempfile
import numpy as np
from common import (load,save,read,table,sha,now,resolutions,previous_resolution,
                    efast_directory,OUTPUTS,METRICS)
import plot_helpers as gp
from report import final_report,scientific_tables
from verify_exports import validate
from finish_report import finish


def check(reference,receipt):
    reference=Path(reference)
    with tempfile.TemporaryDirectory(prefix='efast_refinement_SYNTHETIC_') as temp:
        root=Path(temp)
        shutil.copytree(reference/'inputs',root/'inputs',ignore=shutil.ignore_patterns('designs','__pycache__'))
        shutil.copytree(reference/'sobol',root/'sobol')
        (root/'validation').mkdir();shutil.copyfile(reference/'validation/method_alignment.json',root/'validation/method_alignment.json')
        shutil.copyfile(reference/'preparation.json',root/'preparation.json')
        p=load(reference/'plan.json');p.update(reference_resolutions=p['resolutions'],resolutions=[32769,65537,131073],
            previous_root=str(reference),model_tasks=15890,model_points=3227362530,synthetic_test_only=True)
        save(root/'plan.json',p)
        for n in p['reference_resolutions']:
            directory=root/'inputs/previous_efast'/('N%d'%n);directory.mkdir(parents=True)
            for path in (reference/'efast'/('N%d'%n)).iterdir():
                if path.is_file():shutil.copyfile(path,directory/path.name)
        with np.load(reference/'efast/N16385/indices_repeats.npz') as z:
            base={k:z[k].copy() for k in z.files}
        for i,n in enumerate(p['resolutions'],1):
            d=efast_directory(root,p,n);(d/'figures').mkdir(parents=True)
            values=base['values']+i*.000123
            mean=values.mean(axis=0);spread=np.ptp(values,axis=0);low=previous_resolution(p,n)
            with np.load(efast_directory(root,p,low)/'indices_repeats.npz') as z:delta=np.abs(mean-z['values'].mean(axis=0))
            rows=[]
            for oi,o2 in enumerate(base['oxygen']):
                for k,out in enumerate(OUTPUTS):
                    for j,(index,param,_) in enumerate(METRICS):
                        rows.append(dict(method='eFAST',scope='global',N=n,M=4,O2_pct=o2,output=out,index=index,
                            parameter_i=param,parameter_j='',value_mean=mean[oi,k,j],repeat_min=values[:,oi,k,j].min(),
                            repeat_max=values[:,oi,k,j].max(),repeat_range=spread[oi,k,j],resolution_low=low,
                            resolution_delta=delta[oi,k,j],convergence_status='passed' if spread[oi,k,j]<=.05 and delta[oi,k,j]<=.02 else 'not_passed',included=True))
            table(d/'indices_summary.tsv.gz',rows)
            np.savez_compressed(d/'indices_repeats.npz',**dict(base,values=values))
            np.savez_compressed(d/'indices_summary.npz',mean=mean,repeat_range=spread,resolution_delta=delta,oxygen=base['oxygen'])
            save(d/'completion.json',dict(status='complete',N=n,M=4,n_repeats=5,n_oxygen=201,n_indices=11256,
                                         summary_sha256=sha(d/'indices_summary.tsv.gz'),synthetic_test_only=True))
            classes=gp.classify(values,METRICS,base['oxygen'])
            for row in classes:row['bootstrap_unit']='whole independent phase-repeat curve'
            table(d/'classification.tsv',classes)
            scientific_tables(root,d,values,base['oxygen'],METRICS,'eFAST',n)
        final_report(root)
        # Every source checksum in this synthetic fixture is checked by its own
        # completion receipts. Production additionally verifies the transfer manifest.
        save(root/'collection_manifest.json',dict(files=[],synthetic_test_only=True))
        validate(root)
        save(root/'publication/manifest.json',dict(files=[]))
        finish(root)
        rows=read(root/'comparison/method_pointwise_comparison.tsv.gz')
        expected=(1+len(resolutions(p)))*11256
        assert len(rows)==expected
        bridged=[r for r in rows if int(r['efast_N'])==32769]
        assert len(bridged)==11256 and all(r['diagnostic_status']!='no_adjacent_new_design' for r in bridged)
        assert max(abs(float(r['resolution_delta'])-.000123) for r in bridged)<1e-12
        budgets=read(root/'comparison/design_budget.tsv')
        assert sum(int(r['parameter_oxygen_design_positions']) for r in budgets
                   if r['method']=='eFAST' and r['evaluated_in_this_batch']=='True')==3227362530
        assert next(int(r['parameter_oxygen_design_positions']) for r in budgets if r['method']=='Sobol')==3951820800
        text=(root/'interpretation_zh.md').read_text();assert '131073' in text and str(p['model_points']) in text
        assert 'gamma_mu的ST氧偏好在两种方法间相反' not in text
        save(receipt,dict(status='passed',synthetic_test_only=True,reference_root=str(reference),
             certified_reference_resolutions=p['reference_resolutions'],new_resolutions=p['resolutions'],
             comparison_rows_verified=expected,boundary_delta_verified=.000123,
             arbitrary_series_count=8,new_result_synthesis_exercised=True,model_calls=0,estimator_calls=0,
             temporary_outputs_discarded=True,completed_utc=now()))
    print('Synthetic refinement report integration passed; no scientific results generated')


if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__);a.add_argument('--reference',type=Path,required=True)
    a.add_argument('--receipt',type=Path,required=True);args=a.parse_args();check(args.reference,args.receipt)
