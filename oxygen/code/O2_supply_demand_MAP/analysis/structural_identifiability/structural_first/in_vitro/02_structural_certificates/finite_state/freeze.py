"""Static checks and a fresh immutable release; no scientific imports."""
import ast
from pathlib import Path
import shutil
import subprocess
import sys
from protocol import read,table,sha,roles_valid,successful_worker,reservation,CAP

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();release=sys.argv[3];code=Path(__file__).parent
dest=study/'releases'/release;assert not dest.exists()
checks=[]
def good(name,condition=True):assert condition,name;checks.append(dict(check=name,passed=True))
roles_valid(read(study/'inputs/parameter_roles.tsv'));good('original_twenty_roles')
good('semantic_SIF_gate',read(study/'verification/B0/semantic_gate.tsv')[0]['status']=='PASS')
good('semantic_SIF_policy',read(study/'verification/B0/policy_gate.tsv')[0]['status']=='PASS')
for r in read(study/'inputs/frozen_source_verification.tsv'):good('production:'+r['path'],sha(project/r['path'])==r['sha256'])
for p in code.glob('*.py'):
    tree=ast.parse(p.read_text());good('syntax_'+p.name)
    for n in ast.walk(tree):
        if isinstance(n,(ast.Import,ast.ImportFrom)):
            names=[a.name.split('.')[0] for a in n.names] if isinstance(n,ast.Import) else [(n.module or '').split('.')[0]]
            good('imports_'+p.name+'_'+str(n.lineno),not set(names)&{'json','jpype','jnius','symbtools','strike_goldd'})
        if isinstance(n,ast.Call):
            attr=n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id if isinstance(n.func,ast.Name) else ''
            good('no_timer_'+p.name+'_'+str(n.lineno),attr not in ('alarm','setitimer','generic_rank','strike_goldd') and not any(k.arg=='timeout' for k in n.keywords))
for p in code.glob('*.sh'):subprocess.run(['bash','-n',str(p)],check=True);good('shell_'+p.name)
for p in code.glob('*.R'):
    subprocess.run([str(project/'scripts/agentRrunner.sh'),'-e','invisible(parse(file='+repr(str(p))+'))'],check=True,stdout=subprocess.DEVNULL);good('R_parse_'+p.name)
for gate,app,policy,state,exitcode,expected in [([dict(status='PASS')],[],[dict(status='PASS')],'COMPLETED','0:0',False),
 ([dict(status='PASS')],[dict(status='PASS')],[],'COMPLETED','0:0',False),
 ([dict(status='PASS')],[dict(status='PASS')],[dict(status='PASS')],'COMPLETED','0:0',True)]:
    good('partial_collection_gate_'+str(len(checks)),successful_worker(gate,app,policy,state,exitcode)==expected)
good('site_and_budget_reservation',reservation(9,43200,43200)==43200)
dest.mkdir()
for p in code.iterdir():
    if p.is_file() and p.suffix in ('.py','.R','.sh'):shutil.copy2(p,dest/p.name)
table(study/'reports'/('static_'+release+'_tests.tsv'),checks)
manifest=[]
for parent in ('inputs','runtime',str(dest.relative_to(study))):
    for p in sorted((study/parent).rglob('*')):
        if p.is_file():manifest.append(dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)))
table(study/('source_manifest_'+release+'.tsv'),manifest)
print('FROZEN_RELEASE_STATIC_GATE_PASSED',release,len(checks),'checks;',len(manifest),'source files')
