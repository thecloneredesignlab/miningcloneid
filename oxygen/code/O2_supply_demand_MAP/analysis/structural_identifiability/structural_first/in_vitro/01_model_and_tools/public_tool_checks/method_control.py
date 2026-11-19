"""Small native StrikePy controls only; no culture model or recovery calculation."""
from pathlib import Path
import csv
import importlib.util
import os
import random
import runpy
import shutil
import sys

study=Path(sys.argv[1]);case=sys.argv[2]
assert case in ('C01','C02','C03','C04','C05')
out=study/'controls'/case
os.environ['CULTURE_PYTHON_RUN']='../controls/'+case
boot=runpy.run_path(str(study/'releases/r02/runtime_python.py'),run_name='runtime_control_bootstrap')
sym=boot['sympy'];np=boot['numpy']
seed=20261116+int(case[1:]);random.seed(seed);np.random.seed(seed)
app=out/'native_app';assert not app.exists(),'No stale result reuse'
(app/'models').mkdir(parents=True)
(app/'models/__init__.py').write_text('')
native=study/'sources/StrikePy/StrikePy'
shutil.copy2(native/'strike_goldd.py',app/'strike_goldd.py')
models={
 'C02':"x=[x1];p=[[a]];h=[x1];u=[];w=[];f=[-a*x1]",
 'C03':"x=[x1];p=[[a],[b]];h=[x1];u=[];w=[];f=[-(a+b)*x1]",
 'C04':"x=[x1];p=[[a],[b]];h=[x1];u=[];w=[];f=[-a*x1**b]",
 'C05':"x=[x1];p=[[a],[b]];h=[x1];u=[];w=[];f=[-sym.exp(a+b)*x1]"}
if case=='C01':
    model=(study/'controls/published_C2M_source.py').read_text()
else:
    model="import sympy as sym\nx1,a,b=sym.symbols('x1 a b',positive=True)\n"+models[case]+"\nvariables_locales=locals().copy()\n"
(app/'models/case_model.py').write_text(model)
(app/'options.py').write_text("modelname='case_model'\ncheckObser=1\nmaxLietime=60\nnnzDerU="+('[2]' if case=='C01' else '[]')+"\nnnzDerW=[]\nprev_ident_pars=[]\n")
sys.path.insert(0,str(app));os.chdir(app)
spec=importlib.util.spec_from_file_location('native_control_'+case,app/'strike_goldd.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
try:
    module.strike_goldd()
    summaries=list((app/'results').glob('id_results_*.txt'));assert len(summaries)==1
    text=summaries[0].read_text()
    positive=case in ('C01','C02','C04')
    all_id=('All its parameters are locally structurally identifiable.' in text or
            'All its parameters are structurally identifiable.' in text)
    negative=('These parameters are unidentifiable:' in text and '[a, b]' in text)
    match=all_id if positive else negative and not all_id
    with (out/'control_result.tsv').open('w',newline='') as stream:
        writer=csv.writer(stream,delimiter='\t');writer.writerow(['case','seed','expected','observed_match','scope'])
        writer.writerow([case,seed,'all_parameters_locally_identifiable' if positive else 'individual_a_b_nonidentifiable',match,'native small-control result only; no project classification'])
    if not match:
        raise SystemExit('Native result does not match the preregistered analytic control')
    print('NATIVE_METHOD_CONTROL_COMPLETED',case,flush=True)
finally:
    with (out/'final_imported_modules.tsv').open('w') as stream:
        stream.write('module\tpath\n')
        for name,module in sorted(sys.modules.items()):
            stream.write(name+'\t'+str(getattr(module,'__file__','builtin'))+'\n')
