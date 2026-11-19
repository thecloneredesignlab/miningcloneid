"""No model evaluation. Parse/API/decision/accounting fixtures before release."""
from pathlib import Path
import ast
import importlib.util
import subprocess
import sys
import tempfile
from common import *
from controller import cost_admission,result_quality
p=Path(__file__).resolve().parent;study=Path(sys.argv[1]);rows=[]
def check(n,ok):
 rows.append(dict(check=n,passed=bool(ok)))
 if not ok:raise AssertionError(n)
for f in p.glob('*.py'):ast.parse(f.read_text(),filename=str(f))
check('python_syntax',True)
subprocess.run(['bash','-n',str(p/'runner.sh')],check=True);check('shell_syntax',True)
subprocess.run([str(p.parents[5]/'scripts/agentRrunner.sh'),'-e','for(f in list.files("'+str(p)+'",pattern="[.]R$",full.names=TRUE))parse(f);cat("PARSE_ONLY_PASS\\n")'],check=True)
check('R_parse_only',True)
roles_valid(read(study/'inputs/parameter_roles.tsv'));check('original14_plus6_roles',True)
check('forbidden_artifact_extensions_absent',not any(f.suffix.lower() in ('.json','.jsonl','.js','.java','.jar') for f in p.rglob('*')))
a,b=cost_admission([dict(seconds=.5,feasible=True)]*400);check('cost_floor_and_total',a==1 and b==87600)
a,b=cost_admission([dict(seconds=5,feasible=True)]*400);check('expensive_search_rejected',b>100800)
base=dict(score='1',prediction_pass='TRUE',local_pass='TRUE',complete_poll='TRUE',best_found='1',key='a',best_found_key='a')
check('two_D_required',not result_quality([base],'D') and result_quality([base,base],'D'))
check('two_N_required',not result_quality([base],'N') and result_quality([base,base],'N'))
check('lower_unchecked_candidate_blocks_N',not result_quality([base,dict(base,best_found='.99',best_found_key='b')],'N'))
check('incomplete_poll_blocks_N',not result_quality([dict(base,complete_poll='FALSE')]*2,'N'))
check('one_origin_allocation_safe',reservation(0,0,43200,0)==43200)
try:reservation(CAP,0,43200,0);ok=False
except ValueError:ok=True
check('budget_exhaustion_stops',ok)
for r in read(study/'inputs/evidence_registry.tsv'):check('sealed_input_'+r['path'],sha(study/r['path'])==r['sha256'])
table(study/'reports/static_tests.tsv',rows);print('STATIC_PASS',len(rows))
