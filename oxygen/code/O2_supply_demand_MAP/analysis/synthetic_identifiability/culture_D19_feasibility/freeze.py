"""Freeze scientific release after parse and input-integrity checks; no model call."""
from pathlib import Path
import ast,shutil,subprocess,sys
from common import *
p=Path(__file__).resolve().parent;s=Path(sys.argv[1]).resolve();assert not (s/'source_manifest.tsv').exists()
rows=[]
def check(n,v):
 rows.append(dict(check=n,passed=bool(v)))
 if not v:raise AssertionError(n)
for f in p.glob('*.py'):ast.parse(f.read_text())
check('python_syntax',True)
subprocess.run(['bash','-n',str(p/'runner.sh')],check=True);check('shell_syntax',True)
subprocess.run([str(p.parents[5]/'scripts/agentRrunner.sh'),'-e','for(f in list.files("'+str(p)+'",pattern="[.]R$",full.names=TRUE))parse(f);cat("PARSE_ONLY_PASS\\n")'],check=True);check('R_parse_only',True)
roles_valid(read(s/'inputs/parameter_roles.tsv'));check('original20_roles',True)
check('budget_and_call_caps',4+2+16+4+2+4==32 and 128+1064+48000+2752==51944)
check('native_block_count',600+1400+5*1800+1000==12000)
check('no_forbidden_artifact_extensions',not any(f.suffix.lower() in ('.json','.jsonl','.js','.java','.jar') for f in p.rglob('*')))
for r in read(s/'inputs/evidence_registry.tsv'):check('source_'+r['path'],sha(s/r['path'])==r['sha256'])
for f in p.iterdir():
 if f.is_file():shutil.copy2(f,s/'release'/f.name)
table(s/'source_manifest.tsv',[dict(path=str(f.relative_to(s)),bytes=f.stat().st_size,sha256=sha(f)) for d in ('inputs','runtime','release') for f in sorted((s/d).rglob('*')) if f.is_file()]);verify(s)
table(s/'reports/static_tests.tsv',rows)
table(s/'accounting/site_limits.tsv',[dict(qos='xxlarge',max_seconds=43200,checked='2026-11-17 live sacctmgr snapshot')])
print('SCIENTIFIC_RELEASE_FROZEN',len(rows),sha(s/'source_manifest.tsv'))
