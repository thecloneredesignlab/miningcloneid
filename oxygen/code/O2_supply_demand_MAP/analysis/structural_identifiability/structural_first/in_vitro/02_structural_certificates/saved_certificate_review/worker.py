"""Separate saved verification and native execution entry points for finite controls."""
import datetime
import importlib.util
import os
from pathlib import Path
import random
import sys
import shutil
from artifact import atomic_table,start_call,reconstruct_C03,verify_matrix,select_matrix
from base_contract import read,table,sha,static_contract,select_native_matrix

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);return module

def flat(v,sp,np):
    if isinstance(v,sp.MatrixBase):v=list(v)
    if isinstance(v,np.ndarray):v=v.tolist()
    if isinstance(v,(list,tuple)):return [a for p in v for a in flat(p,sp,np)]
    return [v]

def interfaces(study,sp,np,out):
    rows=[]
    for case in ('C04','C05'):
        folder=study/'fixtures'/case
        static_contract(case,(folder/'model.py').read_text(),(folder/'options.py').read_text())
        model=load('checked_fixture_'+case,folder/'model.py')
        namespace=dict(model=model,np=np,zeros=sp.zeros,w1vector=[],w1vector_dot=[])
        exec(compile((study/'inputs/validated_packing_span.py').read_text(),'verified_pinned_packing','exec'),namespace)
        assert namespace['xaug'].shape==namespace['faug'].shape==(3,1)
        assert sp.Matrix(namespace['xaug'])==sp.Matrix([model.x1,model.a,model.b])
        assert sp.Matrix(namespace['faug'])==sp.Matrix([model.f[0],0,0])
        assert all(v.is_positive is True for v in (model.x1,model.a,model.b))
        rows.append(dict(case=case,static_contract=True,native_array_contract=True,solver_calls=0))
    atomic_table(out/'future_interface_gate.tsv',rows)

def collect_native(study,case,sp,out):
    source=study/'controls'/case
    metadata=read(source/'native_completion.tsv')[0]
    assert metadata['native_returned']=='True' and metadata['case']==case
    assert metadata['input_sha256']==sha(study/'fixtures'/case/'model.py')
    for row in read(source/'native_artifact_inventory.tsv'):
        assert sha(source/row['path'])==row['sha256'],'Native artifact changed'
    results=source/'native_app/results'
    final=select_matrix(results,int(metadata['nd']))
    summary=list(results.glob('id_results_*.txt'));assert len(summary)==1
    text=summary[0].read_text();assert 'Total execution time:' in text
    if case=='C05':
        from artifact import summary_sets
        ids,nonids=summary_sets(text);assert not ids and nonids=={'a','b'}
    else:assert 'All its parameters are locally structurally identifiable.' in text
    result=verify_matrix(case,metadata,final,sp,out,source/'native_final_frame.txt')
    result.update(run_id=study.name,input_sha256=metadata['input_sha256'],native_calls_during_collection=0)
    atomic_table(out/'result.tsv',[result]);return result

def native(study,case,sp,np,out):
    if case not in ('C04','C05'):raise ValueError('Native case forbidden')
    assert not (out/'call_started.tsv').exists() and not (out/'native_completion.tsv').exists()
    app=out/'native_app';assert not app.exists();(app/'models').mkdir(parents=True)
    (app/'models/__init__.py').write_text('')
    root=study/'upstream/StrikePy/StrikePy'
    shutil.copy2(root/'strike_goldd.py',app/'strike_goldd.py');shutil.copytree(root/'functions',app/'functions')
    shutil.copy2(study/'fixtures'/case/'model.py',app/'models/case_model.py');shutil.copy2(study/'fixtures'/case/'options.py',app/'options.py')
    static_contract(case,(app/'models/case_model.py').read_text(),(app/'options.py').read_text())
    sys.path.insert(0,str(app));os.chdir(app)
    options=load('options',app/'options.py');load('models',app/'models/__init__.py')
    model=load('models.case_model',app/'models/case_model.py');core=load('native_strike',app/'strike_goldd.py')
    assert sys.modules['options'] is options
    captures={};calls=[];subs=[]
    def profile(frame,event,arg):
        if event!='return':return
        file=Path(frame.f_code.co_filename)
        if frame.f_code.co_name=='strike_goldd' and file==app/'strike_goldd.py':captures.update(frame.f_locals)
        if frame.f_code.co_name=='generic_rank' and file==study/'runtime/private_prefix/symbtools/core.py':
            values=frame.f_locals;idx=len(calls)
            calls.append(dict(call=idx,seed=values.get('seed'),rank=str(arg),precisions='native 100,200,300; conditional 400'))
            for symbol,value in values.get('rnst',[]):subs.append(dict(call=idx,symbol=str(symbol),value=str(value)))
    seed=20261116+int(case[1:]);random.seed(seed);np.random.seed(seed)
    start_call(out/'call_started.tsv',case,study.name)
    returned=False;error=''
    try:
        sys.setprofile(profile);core.strike_goldd();returned=True
    except BaseException as exc:
        error=type(exc).__name__;raise
    finally:
        sys.setprofile(None)
        metadata=dict(case=case,run_id=study.name,native_returned=returned,exception_type=error,seed=seed,input_sha256=sha(study/'fixtures'/case/'model.py'),native_core_sha256=sha(app/'strike_goldd.py'),nd=captures.get('nd',''),rank=captures.get('rango',''),identifiable=','.join(sorted(str(v) for v in flat(captures.get('p_id',[]),sp,np))),nonidentifiable=','.join(sorted(str(v) for v in flat(captures.get('p_un',[]),sp,np))),measured=','.join(str(v) for v in flat(captures.get('meas_x',[]),sp,np)),observable=','.join(str(v) for v in flat(captures.get('obs_states',[]),sp,np)),column_order=','.join(str(v) for v in flat(captures.get('xaug',[]),sp,np)),positive_real_symbols=all(v.is_positive is True and v.is_real is True for v in (model.x1,model.a,model.b)),order_provenance='persisted_native_return_frame',completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        if returned and 'onx' in captures:
            frame_text='onx = '+str(sp.Matrix(captures['onx']).tolist())
            p=out/'native_final_frame.txt'
            with p.open('x') as stream:stream.write(frame_text);stream.flush();os.fsync(stream.fileno())
        table(out/'rank_calls.tsv',calls,['call','seed','rank','precisions'])
        table(out/'random_substitutions.tsv',subs,['call','symbol','value'])
        inventory=[]
        for p in sorted(app.rglob('*')):
            if p.is_file():inventory.append(dict(path=str(p.relative_to(out)),sha256=sha(p),bytes=p.stat().st_size))
        for name in ('native_final_frame.txt','rank_calls.tsv','random_substitutions.tsv'):
            p=out/name
            if p.is_file():inventory.append(dict(path=name,sha256=sha(p),bytes=p.stat().st_size))
        atomic_table(out/'native_artifact_inventory.tsv',inventory)
        atomic_table(out/'native_completion.tsv',[metadata])
    assert returned
    collect_native(study,case,sp,out)

def main():
    study=Path(sys.argv[1]);mode=sys.argv[2]
    if mode not in ('review_C03','native_C04','native_C05','collect_C04','collect_C05'):raise ValueError('Entry point rejected before scientific imports')
    is_native=mode.startswith('native_');case=mode[-3:]
    out=study/('controls' if is_native else 'review' if case=='C03' else 'validation_retry')/case
    out.mkdir(parents=True,exist_ok=True)
    if (out/'result.tsv').exists():raise ValueError('Immutable completed result exists')
    for row in read(study/'source_manifest.tsv'):assert sha(study/row['path'])==row['sha256'],'Frozen source drift'
    def guard(event,args):
        root=args[0].split('.')[0] if event=='import' else ''
        if root=='ipydex':raise ImportError('Approved debugger exclusion')
        forbidden={'json','_json','simplejson','orjson','ujson','rapidjson','jpype','jnius','execjs','py_mini_racer'}
        if not is_native:forbidden|={'symbtools','strike_goldd','native_strike','functions'}
        bad=root in forbidden
        if not is_native and event=='exec':
            name=str(getattr(args[0],'co_filename',''));bad=bad or name.endswith('strike_goldd.py') or '/symbtools/' in name
        if bad:
            with (out/'blocked_imports.tsv').open('a') as stream:stream.write('BLOCKED\t'+event+'\t'+str(args[0])+'\n')
            raise ImportError('Forbidden execution path')
    sys.addaudithook(guard)
    sys.path[:0]=[str(study/'runtime/private_prefix'),'/home/4482173/.local/lib/python3.11/site-packages',*[p for p in os.environ.get('PYTHONPATH','').split(':') if p]]
    try:
        import numpy as np
        import sympy as sp
        assert np.__version__=='2.4.6' and sp.__version__=='1.12'
        if mode=='review_C03':
            source=study/'saved/C03'
            core_text=(study/'upstream/StrikePy/StrikePy/strike_goldd.py').read_text()
            metadata=reconstruct_C03(source,core_text)
            atomic_table(out/'final_order_evidence.tsv',[dict(case='C03',nd=3,origin=metadata['order_provenance'],console_sha256=sha(source/'console.log'),core_sha256=sha(study/'upstream/StrikePy/StrikePy/strike_goldd.py'),branch='nd=3 >= augmented_dimension=3; branch does not decrement nd',source_job_id=21966795,original_slurm_state='FAILED_collector_exit_preserved')])
            result=verify_matrix('C03',metadata,metadata['selected_matrix'],sp,out)
            lower=dict(metadata);lower['nd']=2
            verify_matrix('C03',lower,select_native_matrix(list((source/'native_app/results').glob('obs_ident_matrix_*.txt')),2),sp,out/'intermediate_order2')
            interfaces(study,sp,np,out)
            result.update(run_id=study.name,new_native_C03_calls=0,source_job_id=21966795,input_sha256=sha(source/'native_app/models/case_model.py'))
            atomic_table(out/'result.tsv',[result])
            atomic_table(out/'call_ledger.tsv',[dict(case='C03',new_native_calls=0,randomized_rank_calls=0,scope='saved_artifact_verification_only')])
        elif is_native:
            import symbtools
            assert symbtools.__version__=='0.4.1'
            native(study,case,sp,np,out)
        else:collect_native(study,case,sp,out)
        print('VERIFIED_RESULT_READY',mode,flush=True)
    finally:
        with (out/'imported_modules.tsv').open('w') as stream:
            stream.write('module\tpath\n')
            for name,module in sorted(sys.modules.items()):stream.write(name+'\t'+str(getattr(module,'__file__','builtin'))+'\n')

if __name__=='__main__':main()
