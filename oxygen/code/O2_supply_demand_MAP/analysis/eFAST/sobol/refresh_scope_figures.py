#!/usr/bin/env python3
"""Redraw main scope figures from saved arrays; never call a sensitivity analyzer."""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True)
    args=parser.parse_args(); root=Path(args.root).resolve()
    retry=root/'recovery_20261018/retry5734'
    if json.loads((root/'status.json').read_text())['status']!='complete':
        raise ValueError('Original final report must be complete first')
    frozen=root/'recovery_20261018/slurm/code_snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol'
    sys.path.insert(0,str(frozen))
    import core as c
    import numpy as np
    c.configure(json.loads((root/'slurm/submission_plan.json').read_text()))
    module_path=Path(__file__).parent/'full_report.py'
    spec=importlib.util.spec_from_file_location('scope_figure_layout_patch',module_path)
    report=importlib.util.module_from_spec(spec); spec.loader.exec_module(report)
    sources={scope:root/scope/'summaries/scope_summary.npz' for scope in ('global','neighborhood10pct')}
    before={scope:sha(path) for scope,path in sources.items()}
    rendered=[]
    for scope,path in sources.items():
        with np.load(path,allow_pickle=False) as arrays:
            report.plot_indices(arrays['median'],arrays['oxygen'],root/scope/'figures',
                'Global full-bound sensitivity' if scope=='global' else 'Median of 500 endpoint five-repeat means',
                None if scope=='global' else arrays['pass_fraction_all_endpoints'])
        names=['s1_parameter_by_oxygen','st_parameter_by_oxygen','s2_pair_by_oxygen']
        if scope=='neighborhood10pct': names.append('s2_diagnostic_fraction')
        rendered.extend(root/scope/'figures'/(name+suffix) for name in names for suffix in ('.pdf','.png'))
    after={scope:sha(path) for scope,path in sources.items()}
    if before!=after: raise ValueError('Saved sensitivity arrays changed during rendering')
    receipt=dict(status='complete',job_id=os.environ['SLURM_JOB_ID'],
                 source_revision=os.environ['LAYOUT_REVISION'],plot_source_sha256=sha(module_path),
                 input_arrays_sha256=before,input_arrays_unchanged=True,
                 sensitivity_calculations_repeated=0,n_redrawn_figures=7,
                 files=[dict(path=str(path.relative_to(root)),bytes=path.stat().st_size,sha256=sha(path)) for path in rendered],
                 completed_utc=c.now())
    receipt_path=retry/'figure_layout_patch.json'; c.save_json(receipt_path,receipt)
    with (root/'collection_manifest.tsv').open() as handle:
        rows={r['path']:r for r in csv.DictReader(handle,delimiter='\t')}
    for item in receipt['files']:
        if item['path'] not in rows: raise ValueError('Main figure absent from original inventory')
        rows[item['path']]=item
    relative=str(receipt_path.relative_to(root))
    rows[relative]=dict(path=relative,bytes=receipt_path.stat().st_size,sha256=sha(receipt_path))
    temp=root/'collection_manifest.tsv.layout_tmp'
    with temp.open('w') as handle:
        writer=csv.DictWriter(handle,fieldnames=['path','bytes','sha256'],delimiter='\t')
        writer.writeheader(); writer.writerows(rows[name] for name in sorted(rows))
    os.replace(str(temp),str(root/'collection_manifest.tsv'))
    print(json.dumps(dict(job_id=receipt['job_id'],redrawn_figures=7,arrays_unchanged=True,analyses_repeated=0)),flush=True)


if __name__=='__main__': main()
