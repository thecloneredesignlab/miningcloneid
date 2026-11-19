"""Freeze one tested SIF harness repair; preserve all original phase returns."""
from pathlib import Path
import re
import shutil
import subprocess
import sys
from common import *
from policy import audit

project=Path(sys.argv[1]).resolve();study=Path(sys.argv[2]).resolve();code=Path(__file__).parent
assert kv(study/'accounting/retry_ledger.tsv')['used']=='0'
j=next(j for j in read(study/'accounting/jobs.tsv') if j['task']=='P03')
assert j['state']=='FAILED' and j['exit_code']=='77:0'
old=study/'verification/P03'
assert read(old/'phase_exit_codes.tsv')[0]==dict(application_exit_code='0',policy_exit_code='77',source_after_exit_code='0')
assert len(read(old/'passage_contract_tests.tsv'))==9 and all(r['passed']=='TRUE' for r in read(old/'passage_contract_tests.tsv'))
assert (old/'blocked_access.tsv').read_text().splitlines()==['WRITE_OUTSIDE_PHASE\topen\t/dev/tty']*5
table(study/'reports/failed_P03_manifest.tsv',[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(old.rglob('*')) if p.is_file()])
dest=study/'retry_release';assert not dest.exists();dest.mkdir()
for name in ('common.py','policy.py'):shutil.copy2(code/name,dest/name)
r=(study/'release/passage.R').read_text()
assert r.count('stopifnot(task == "P03")')==1
r=r.replace('stopifnot(task == "P03")','stopifnot(task == "P03_R1")')
r=r.replace('  write_tsv(data.frame(stage="before_allowlisted_expression_loading",task=task),"bootstrap_receipt.tsv")',
'''  write_tsv(data.frame(stage="before_allowlisted_expression_loading",task=task),"bootstrap_receipt.tsv")
  stopifnot(startsWith(normalizePath(tempdir()),paste0(normalizePath(file.path(out,"tmp")),"/")))
  write_tsv(data.frame(tempdir=tempdir(),scope="phase-owned writable temporary directory"),"temporary_scope.tsv")''')
(dest/'passage.R').write_text(r)
s=(study/'release/runner.sh').read_text()
s=s.replace('task=$2','task=$2\ntest "$task" = P03_R1')
s=s.replace('mkdir -p "$out/tmp/home"','mkdir -p "$out/tmp/home"\nchmod 1777 "$out/tmp"')
s=s.replace('if [ "$task" = P03 ]; then','if [ "$task" = P03_R1 ]; then')
s=s.replace('--bind "$study/inputs:$study/inputs:ro"','--bind "$out/tmp:$out/tmp:rw" --bind "$study/inputs:$study/inputs:ro"')
s=s.replace('"$study/release/passage.R"','"$study/retry_release/passage.R"').replace('"$study/release/policy.py"','"$study/retry_release/policy.py"')
s=s.replace('"$py" -B -S "$study/release/verify_source.py" "$study" > "$out/source_before.txt"',
            '"$py" -B -S "$study/retry_release/verify_retry.py" "$study" > "$out/source_before.txt"')
s=s.replace('"$py" -B -S "$study/release/verify_source.py" "$study" > "$out/source_after.txt" 2>&1',
            '"$py" -B -S "$study/retry_release/verify_retry.py" "$study" > "$out/source_after.txt" 2>&1')
(dest/'runner.sh').write_text(s)
(dest/'verify_retry.py').write_text('''from pathlib import Path
import sys
from common import verify,read,sha
study=Path(sys.argv[1]);verify(study)
for r in read(study/'retry_source_manifest.tsv'):
 p=study/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256']
print('ORIGINAL_AND_RETRY_SOURCE_MANIFESTS_PASS')
''')
# This is a narrow classifier change, not an allowance for terminal/file writes.
# The native guard remains unchanged and continues to reject /dev/tty opens.
checks=[]
for name,log,trace,expected in (
 ('denied_tty_only','WRITE_OUTSIDE_PHASE\topen\t/dev/tty\n','', 'PASS'),
 ('persistent_path_denial','WRITE_OUTSIDE_PHASE\topen\t/tmp/foreign.tsv\n','', 'FAIL'),
 ('excluded_path_denial','BLOCKED\topen\t/tmp/excluded.json\n','', 'FAIL'),
 ('successful_terminal_open','', 'openat(AT_FDCWD, "/dev/tty", O_RDWR) = 3\n','FAIL'),
 ('outside_temporary_directory','', 'mkdir("/tmp/foreign", 0700) = 0\n','FAIL')):
    folder=study/'reports/retry_policy_fixtures'/name;folder.mkdir(parents=True)
    (folder/'blocked_access.tsv').write_text(log)
    (folder/'trace.fixture').write_text('openat(AT_FDCWD, "/usr/lib/libc.so.6", O_RDONLY) = 3\n'+trace)
    table(folder/'loaded_namespaces.tsv',[dict(namespace='base',version='fixture')])
    result=audit(folder,'R')['status'];assert result==expected,name
    checks.append(dict(check=name,passed=True,expected=expected))
subprocess.run(['bash','-n',str(dest/'runner.sh')],check=True)
subprocess.run([str(project/'scripts/agentRrunner.sh'),'-e','invisible(parse(file='+repr(str(dest/'passage.R'))+'))'],check=True,stdout=subprocess.DEVNULL)
checks.append(dict(check='shell_and_R_parse',passed=True,expected='PASS'))
table(study/'reports/retry_static_tests.tsv',checks)
table(study/'retry_source_manifest.tsv',[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(dest.iterdir()) if p.is_file()])
(study/'reports/technical_retry.md').write_text('''# One authorized technical retry: SIF harness only

Original P03 application and all9 interface fixtures passed; policy failed.
The guard denied five shell /dev/tty probes. These are nonpersistent terminal
queries, not excluded language/library access. The unchanged native guard still
denies them; the repaired collector records these exact denied queries without
calling them persistent-file violations. Other denials and any successful
out-of-scope terminal/file mutation remain failures, tested by negative fixtures.

The original trace also showed access(TMPDIR,W_OK) failing and R falling back to
an empty /tmp/Rtmp directory. The retry explicitly binds a phase-owned writable
temporary root and gives only that new scratch directory mode1777. An early R
assertion requires its actual tempdir to lie there. Successful mkdir/unlink/
rename operations are included in the repaired policy audit, not silently ignored.

No numerical model, parameter, observer, fixture input/expectation, saved matrix,
certificate or source function changes. Only the P03 task/output label, temporary
root contract and narrow policy classification change. All failed originals are
preserved. The9 interface fixtures are revalidated for the repaired environment;
no project structural case, certificate or boundary calculation is repeated.
''')
table(study/'accounting/retry_ledger.tsv',[dict(key='used',value=1),dict(key='maximum',value=1),dict(key='failed_task',value='P03'),dict(key='failed_job',value=j['job_id']),dict(key='retry_task',value='P03_R1')])
print('ONE_TECHNICAL_RETRY_FROZEN_AND_STATIC_TESTED',len(checks))
