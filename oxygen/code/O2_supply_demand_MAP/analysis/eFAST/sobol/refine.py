"""Checked reuse of completed nested Sobol prefixes, without changing the model."""
import csv
import gzip
from pathlib import Path
import subprocess

import numpy as np
import core as c


def audit_source(root, current_hashes, cases, oxygen, resolutions):
    source=c.load_json(root/'slurm/submission_plan.json')
    status=c.load_json(root/'status.json')
    if status.get('status')!='complete' or source['cases']!=cases:
        raise ValueError('Reuse requires a complete matching representative pilot')
    if source['oxygen_pct']!=oxygen or source['n_repeats']!=5:
        raise ValueError('Reuse oxygen/repeats differ')
    if source['max_N']!=min(resolutions) or source['max_N']>=max(resolutions):
        raise ValueError('Source resolution must match the lowest refinement resolution')
    if source['sif']['sha256']!=c.SIF_SHA256:
        raise ValueError('Reuse SIF differs')
    # Sampling prefixes are compared exactly below. Every file executed by the
    # R model and its input helpers must still match the submitted source.
    required=[p for p in source['code_sha256'] if p.endswith('.R') or '/model/' in p or
              p.endswith(('efast.py','neighborhood.py','perron_high_precision_safe.cpp'))]
    if not required:
        raise ValueError('Missing model implementation hashes')
    for name in required:
        if current_hashes.get(name)!=source['code_sha256'][name]:
            raise ValueError('Cannot reuse changed model/helper: '+name)
    for name,expected in source['input_sha256'].items():
        if c.sha(root/name)!=expected:
            raise ValueError('Reuse source input changed: '+name)
    summaries={}
    for case in cases:
        target=c.case_directory(root,case)
        receipt=c.load_json(target/'completion.json')
        if receipt.get('status')!='complete' or c.sha(target/'indices_summary.tsv')!=receipt['summary_sha256'] or c.sha(target/'indices_repeats.tsv.gz')!=receipt['indices_sha256']:
            raise ValueError('Reuse case summary changed/incomplete: '+case['case_id'])
        summaries[str(target/'indices_summary.tsv')]=receipt['summary_sha256']
    return dict(root=str(root),source_N=source['max_N'],
                source_plan_sha256=c.sha(root/'slurm/submission_plan.json'),
                source_status_sha256=c.sha(root/'status.json'),
                source_git_commit=source['git_commit'],verified_model_files=required,summary_sha256=summaries,
                policy='exact nested parameter prefix; immutable completed output receipts; unchanged model and SIF')


def check_design_prefix(root,case,repeat,new_samples,bounds,scramble_seed):
    target=c.case_directory(root,case)
    design=target/'runs'/('R%02d'%repeat)
    meta=c.load_json(design/'design.json')
    if meta['scramble_seed']!=scramble_seed or meta['parameters']!=list(c.ACTIVE):
        raise ValueError('Reuse scramble/parameter order differs')
    if c.sha(design/'ordered_samples.npz')!=meta['ordered_samples_sha256']:
        raise ValueError('Reuse ordered samples changed')
    if c.sha(target/'parameter_bounds.tsv')!=meta['bounds_sha256']:
        raise ValueError('Reuse bounds changed')
    old_bounds=c.read(target/'parameter_bounds.tsv')
    fields=('encoded_lower','encoded_upper','natural_lower','natural_upper')
    for old,new in zip(old_bounds,bounds):
        if old['parameter']!=new['parameter'] or old['transform']!=new['transform'] or any(
                float(old[k])!=float(new[k]) for k in fields):
            raise ValueError('Reuse parameter distribution differs')
    with np.load(design/'ordered_samples.npz') as old:
        if not np.array_equal(new_samples[:meta['N']*c.STEP],old['natural']):
            raise ValueError('Reuse prefix is not exactly the old design')


def source_block(root,relative,n):
    part=root/relative
    receipt_path=Path(str(part/'outputs.tsv.gz')+'.receipt.json')
    receipt=c.load_json(receipt_path)
    meta=c.load_json(part/'metadata.json')
    if meta['N']!=n or receipt.get('status')!='complete':
        raise ValueError('Incomplete/mismatched source block')
    for name,field in [('samples.tsv.gz','samples_sha256'),('metadata.json','metadata_sha256')]:
        if c.sha(part/name)!=receipt[field]:
            raise ValueError('Reuse source block changed: '+str(part/name))
    return dict(part=str(part),N=n,receipt_sha256=c.sha(receipt_path),
                outputs_sha256=receipt['outputs_sha256'],samples_sha256=receipt['samples_sha256'],
                metadata_sha256=receipt['metadata_sha256'])


def fuse_outputs(prefix,tail,target,n_old,n_new,oxygen):
    """Restore global sample IDs and require a complete nonduplicate rectangle."""
    # Parse both sources before promoting the new file; incomplete rows must
    # never be silently removed or passed to the Sobol analyzer.
    c.output_array(prefix,n_old,oxygen)
    c.output_array(tail,n_new,oxygen)
    fields=None
    temporary=target.with_name(target.name+'.incomplete.gz')
    with gzip.open(temporary,'wt',newline='') as out:
        for path,offset in [(prefix,0),(tail,n_old)]:
            with gzip.open(path,'rt',newline='') as handle:
                reader=csv.DictReader(handle,delimiter='\t')
                if fields is None:
                    fields=reader.fieldnames
                    writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t',lineterminator='\n')
                    writer.writeheader()
                elif reader.fieldnames!=fields:
                    raise ValueError('Prefix/increment output columns differ')
                for row in reader:
                    row['sample_id']=int(row['sample_id'])+offset
                    writer.writerow(row)
    c.output_array(temporary,n_old+n_new,oxygen)
    temporary.replace(target)


def evaluate_increment(root,part,plan,evaluator):
    meta=c.load_json(part/'metadata.json'); reuse=meta['reuse']
    source=Path(reuse['part']); prefix=source/'outputs.tsv.gz'
    source_receipt=Path(str(prefix)+'.receipt.json')
    if c.sha(source_receipt)!=reuse['receipt_sha256'] or c.sha(prefix)!=reuse['outputs_sha256']:
        raise ValueError('Reuse output/receipt changed')
    for name,field in [('metadata.json','metadata_sha256'),('samples.tsv.gz','samples_sha256')]:
        if c.sha(source/name)!=reuse[field]: raise ValueError('Reuse source changed: '+name)
    previous=c.load_json(source_receipt)
    rows=c.read(part/'samples.tsv.gz'); old=c.read(source/'samples.tsv.gz'); n=reuse['N']
    if len(old)!=n or len(rows)!=meta['n_samples'] or n>=len(rows):
        raise ValueError('Invalid prefix/increment sizes')
    for i in range(n):
        if int(old[i]['sample_id'])!=i+1 or int(rows[i]['sample_id'])!=i+1 or any(
                float(old[i][p])!=float(rows[i][p]) for p in c.ACTIVE):
            raise ValueError('Prefix parameter vector changed')
    incremental=part/'incremental'; incremental.mkdir(exist_ok=True)
    tail=[dict(row,sample_id=i+1) for i,row in enumerate(rows[n:])]
    if (incremental/'samples.tsv.gz').exists():
        saved=c.read(incremental/'samples.tsv.gz')
        if len(saved)!=len(tail) or any(int(a['sample_id'])!=int(b['sample_id']) or any(
                float(a[p])!=float(b[p]) for p in c.ACTIVE) for a,b in zip(saved,tail)):
            raise ValueError('Existing increment samples differ')
    else:
        c.table(incremental/'samples.tsv.gz',tail)
    tail_meta=dict(meta); tail_meta.pop('reuse')
    tail_meta.update(n_samples=len(tail),sample_id_offset=n,validate_reference=False,
                     samples_sha256=c.sha(incremental/'samples.tsv.gz'),
                     parent_metadata_sha256=c.sha(part/'metadata.json'))
    c.save_json(incremental/'metadata.json',tail_meta)
    out=incremental/'outputs.tsv.gz'
    if not c.complete(out):
        subprocess.run(['Rscript','--vanilla',str(evaluator),str(incremental),plan['fit_root'],
                        str(root/'inputs/fixed_o2_dominant_ploidy_201grid.tsv')],check=True)
    receipt=c.load_json(str(out)+'.receipt.json')
    if receipt['metadata_sha256']!=c.sha(incremental/'metadata.json') or receipt['samples_sha256']!=c.sha(incremental/'samples.tsv.gz'):
        raise ValueError('Increment provenance differs')
    target=part/'outputs.tsv.gz'
    fuse_outputs(prefix,out,target,n,len(tail),meta['oxygen_pct'])
    c.save_json(str(target)+'.receipt.json',dict(status='complete',n_rows=len(rows)*len(meta['oxygen_pct']),
        metadata_sha256=c.sha(part/'metadata.json'),samples_sha256=c.sha(part/'samples.tsv.gz'),
        outputs_sha256=c.sha(target),reference_checks=previous['reference_checks'],
        common_reference=previous['common_reference'],reuse=reuse,
        incremental_receipt_sha256=c.sha(str(out)+'.receipt.json'),
        incremental_sample_id_offset=n,elapsed_new_seconds=receipt['elapsed_seconds'],
        elapsed_reused_seconds=previous['elapsed_seconds'],
        elapsed_seconds=previous['elapsed_seconds']+receipt['elapsed_seconds'],
        host=receipt['host'],job_id=receipt['job_id']))
