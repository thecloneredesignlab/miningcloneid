"""Adversarial collector/policy metadata fixtures, not scientific-worker traces."""
import ast
from pathlib import Path
import sys
from common import read,table,sha,verify,accepted,reservation
from policy import audit

study=Path(sys.argv[1]).resolve();root=study/'reports/adversarial_metadata_fixtures';root.mkdir(parents=True)
checks=[]
def check(name,condition):assert condition,name;checks.append(dict(check=name,passed=True))
for name,trace,inventory,expected in (
    ('complete',True,'normal','PASS'),('missing_trace',False,'normal','FAIL'),
    ('missing_inventory',True,'missing','FAIL'),('excluded_inventory_literal',True,'excluded','FAIL'),
    ('outside_write_literal',True,'normal','FAIL'),('excluded_path_literal',True,'normal','FAIL')):
    d=root/name;d.mkdir()
    if trace:
        target='/usr/lib/libc.so.6';flags='O_RDONLY'
        if name=='outside_write_literal':target='/tmp/another_scope.tsv';flags='O_WRONLY|O_CREAT'
        if name=='excluded_path_literal':target='/tmp/excluded.json'
        (d/'trace.fixture').write_text('openat(AT_FDCWD, "'+target+'", '+flags+') = 3\n')
    if inventory!='missing':table(d/'imported_modules.tsv',[dict(module='json' if inventory=='excluded' else 'sys',path='fixture')])
    check('policy_'+name,audit(d,'Python')['status']==expected)
fixture=root/'source_identity';fixture.mkdir();source=fixture/'source.tsv';source.write_text('identity\nB01_W1\n')
table(fixture/'source_manifest.tsv',[dict(path='source.tsv',bytes=source.stat().st_size,sha256=sha(source))])
verify(fixture);check('source_identity_accepted',True)
for changed in ('identity\nB01_W2\n','identity\nB04_W1\n',''):
    source.write_text(changed)
    try:verify(fixture)
    except ValueError:check('stale_swapped_partial_rejected_'+str(len(checks)),True)
    else:raise AssertionError('Source identity mismatch admitted')
for state,app,policy in (('OUT_OF_MEMORY','MISSING','MISSING'),('COMPLETED','PASS','MISSING'),('COMPLETED','MISSING','PASS')):
    check('incomplete_return_'+str(len(checks)),not accepted(state,'0:0',app,policy,'PASS'))
for settled,live in ((0,0),(10,0),(60,3600),(10000,0)):
    try:seconds=reservation(settled,live,43200,0)
    except ValueError:check('budget_stop_'+str(settled),settled+live+3600+1800>=14400)
    else:check('budget_admission_'+str(settled),settled+live+seconds+3600+1800<=14400)
manage=(Path(__file__).parent/'manage.py').read_text()
check('pending_submission_reconciliation_guard',"all(j['job_id'].isdigit() for j in jobs)" in manage)
check('duplicate_dispatch_guard',"task not in {j['task'] for j in jobs}" in manage)
for p in Path(__file__).parent.glob('*.py'):ast.parse(p.read_text())
table(study/'reports/adversarial_metadata_tests.tsv',checks)
print('ADVERSARIAL_METADATA_TESTS_PASS',len(checks))
