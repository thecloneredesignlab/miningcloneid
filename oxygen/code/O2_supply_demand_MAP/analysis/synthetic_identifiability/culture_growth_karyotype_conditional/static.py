"""Pre-allocation static and metadata regressions; no local model or CAS."""
import ast
import importlib.util
from pathlib import Path
import shutil
import subprocess
import sys
from common import *
from policy import audit
project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();code=Path(__file__).parent
checks=[]
def check(n,v=True):assert v,n;checks.append(dict(check=n,passed=True))
roles_valid(read(study/'inputs/parameter_roles.tsv'));check('original14_plus6')
for p in code.glob('*.py'):
 tree=ast.parse(p.read_text());check('Python_syntax_'+p.name)
 for node in ast.walk(tree):
  if isinstance(node,(ast.Import,ast.ImportFrom)):
   names=[a.name.split('.')[0] for a in node.names] if isinstance(node,ast.Import) else [(node.module or '').split('.')[0]]
   check('excluded_import_'+p.name+'_'+str(node.lineno),not set(names)&{'json','_json','jpype','jnius','execjs'})
  if isinstance(node,ast.Call):
   name=getattr(node.func,'attr',getattr(node.func,'id',''));check('finite_no_timeout_'+p.name+'_'+str(node.lineno),name not in ('alarm','setitimer','minimize','DEoptim') and not any(k.arg=='timeout' for k in node.keywords))
for p in code.glob('*.R'):
 subprocess.run([str(project/'scripts/agentRrunner.sh'),'-e','invisible(parse(file='+repr(str(p))+'))'],check=True,stdout=subprocess.DEVNULL);check('R_parse_'+p.name)
subprocess.run(['bash','-n',str(code/'runner.sh')],check=True);check('shell_parse')
check('initial_budget_reservation',reservation(0,0,43200,0)==7200)
try:reservation(8000,0,43200,0)
except ValueError:check('reject_overcommit')
else:raise AssertionError('Budget rejection')
for name,excluded,outside,expected in [('complete',False,False,'PASS'),('bad_module',True,False,'FAIL'),('outside_write',False,True,'FAIL')]:
 d=study/'reports/static_fixtures'/name;d.mkdir(parents=True)
 (d/'trace.100').write_text('openat(AT_FDCWD, "/usr/lib/libc.so.6", O_RDONLY) = 3</usr/lib/libc.so.6>\n'+('openat(AT_FDCWD, "/tmp/foreign.tsv", O_WRONLY|O_CREAT) = 3</tmp/foreign.tsv>\n' if outside else ''))
 table(d/'imported_modules.tsv',[dict(module='json' if excluded else 'sys',path='fixture')]);check('trace_regression_'+name,audit(d,'Python')['status']==expected)
interface=(code/'interface.R').read_text();worker=(code/'worker.R').read_text()
check('no_flow_score_component','ivt_flow_loglik_df' not in interface and 'ivt_flow_loglik_df' not in worker)
check('forward_allowance','used<4' in interface)
check('observation_allowance','used<1027' in interface)
check('no_optimization_entry',all(x not in interface+worker for x in ('DEoptim(','nloptr(','optim(','d0_loss(')))
check('candidate_SD_not_fixed','p$sigma_growth' in interface and 'map$fixed$sigma_growth' not in interface)
check('cold_bootstrap_bypassed',"env$.init_cpp_o2simps_backend(" not in (code/'production_bindings.R').read_text())
check('RNG_registry_and_QA','20261116' in interface and '1824' in worker and '1:64' in worker and '65:1024' in worker)
check('closed_case_truth_not_rebuilt','anchor_bundle.rds' in worker)
table(study/'reports/static_tests.tsv',checks)
assert not (study/'release/common.py').exists()
for p in code.iterdir():
 if p.is_file():shutil.copy2(p,study/'release'/p.name)
table(study/'source_manifest.tsv',[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)) for root in ('inputs','runtime','release') for p in sorted((study/root).rglob('*')) if p.is_file()])
print('PREALLOCATION_STATIC_PASS',len(checks))
