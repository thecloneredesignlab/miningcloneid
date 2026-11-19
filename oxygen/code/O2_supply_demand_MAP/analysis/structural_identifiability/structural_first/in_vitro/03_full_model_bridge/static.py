"""Pre-allocation metadata and import/trace regressions; no local CAS."""
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
def check(name,condition=True):assert condition,name;checks.append(dict(check=name,passed=True))
roles_valid(read(study/'inputs/parameter_roles.tsv'));check('14_plus_6_inventory')
for name in ('observer_extension.py','geometry.py'):
    p=code/name;tree=ast.parse(p.read_text());main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    init=next(n.lineno for n in ast.walk(main) if isinstance(n,ast.Call) and isinstance(n.func,ast.Name) and n.func.id=='boot')
    check('imports_after_boot_'+name,all(n.lineno>init for n in ast.walk(main) if isinstance(n,(ast.Import,ast.ImportFrom))))
    spec=importlib.util.spec_from_file_location('static_'+p.stem,p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    class StopBeforeMath(Exception):pass
    def stub(*args):raise StopBeforeMath()
    m.boot=stub;old=sys.argv;sys.argv=[str(p),str(study),'O01']
    try:m.main()
    except StopBeforeMath:check('entry_stub_'+name,'flint' not in sys.modules)
    finally:sys.argv=old
for p in code.glob('*.py'):
    tree=ast.parse(p.read_text());check('syntax_'+p.name)
    for n in ast.walk(tree):
        if isinstance(n,(ast.Import,ast.ImportFrom)):
            names=[a.name.split('.')[0] for a in n.names] if isinstance(n,ast.Import) else [(n.module or '').split('.')[0]]
            check('excluded_import_'+p.name+'_'+str(n.lineno),not set(names)&{'json','_json','jpype','jnius','strike_goldd','symbtools','ipydex'})
        if isinstance(n,ast.Call):
            name=n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id if isinstance(n.func,ast.Name) else ''
            check('no_timeout_optimizer_'+p.name+'_'+str(n.lineno),name not in ('alarm','setitimer','generic_rank','minimize','DEoptim') and not any(k.arg=='timeout' for k in n.keywords))
for p in code.glob('*.sh'):subprocess.run(['bash','-n',str(p)],check=True);check('shell_'+p.name)
for p in code.glob('*.R'):
    subprocess.run([str(project/'scripts/agentRrunner.sh'),'-e','invisible(parse(file='+repr(str(p))+'))'],check=True,stdout=subprocess.DEVNULL);check('R_parse_'+p.name)
for name,which,expected in [('complete','normal','PASS'),('missing_inventory','missing','FAIL'),('excluded_inventory','excluded','FAIL'),('relative_compiler','relative','PASS'),('outside_write','outside','FAIL'),('unresolved_fd','fd','FAIL'),('missing_trace','none','FAIL')]:
    d=study/'reports/policy_fixtures'/name;d.mkdir(parents=True)
    if which!='none':
        line='openat(AT_FDCWD, "/usr/lib/libc.so.6", O_RDONLY) = 3</usr/lib/libc.so.6>\n'
        if which=='relative':
            line+='chdir("'+str(d/'tmp')+'") = 0\nclone(child_stack=NULL, flags=SIGCHLD) = 101\nopenat(AT_FDCWD, "obj.o", O_WRONLY|O_CREAT, 0600) = 4<'+str(d/'tmp/obj.o')+'>\n'
            (d/'trace.101').write_text('unlinkat(5<'+str(d/'tmp')+'>, "obj.o", 0) = 0\n')
        if which=='outside':line+='openat(AT_FDCWD, "/tmp/outside.tsv", O_WRONLY|O_CREAT) = 3</tmp/outside.tsv>\n'
        if which=='fd':line+='openat(7, "unknown.tsv", O_WRONLY|O_CREAT) = 3\n'
        (d/'trace.100').write_text(line)
    if which!='missing':table(d/'imported_modules.tsv',[dict(module='json' if which=='excluded' else 'sys',path='fixture')])
    check('trace_policy_'+name,audit(d,'Python')['status']==expected)
check('budget_initial_reservation',reservation(0,0,43200,0)==28800)
try:reservation(40000,0,43200,0)
except ValueError:check('budget_overcommit_rejected')
else:raise AssertionError('Budget admission')
check('bounded_forward_registration',"used<74L" in (code/'api.R').read_text())
check('original_timed_bootstrap_not_called',"exists('.init_cpp_o2simps_backend'" in (code/'api.R').read_text() and "env$.init_cpp_o2simps_backend(" not in (code/'api.R').read_text())
table(study/'reports/static_tests.tsv',checks)
assert not (study/'release/common.py').exists()
for p in code.iterdir():
    if p.is_file() and p.suffix in ('.py','.R','.sh'):shutil.copy2(p,study/'release'/p.name)
table(study/'source_manifest.tsv',[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)) for root in ('inputs','runtime','release') for p in sorted((study/root).rglob('*')) if p.is_file()])
print('STEP5_STATIC_GATE_PASS',len(checks),'checks')
