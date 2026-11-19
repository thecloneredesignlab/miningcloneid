"""One documented interface-only repair; no model/data calculation locally."""
from pathlib import Path
import ast
import shutil
import subprocess
import sys
from common import *
study=Path(sys.argv[1]).resolve();project=Path(sys.argv[2]).resolve();dest=study/'repair_release';assert not dest.exists()
assert kv(study/'accounting/retry_ledger.tsv')['used']=='0'
assert read(study/'accounting/jobs.tsv')[0]['task']=='I01' and read(study/'accounting/jobs.tsv')[0]['state']=='FAILED'
shutil.copytree(study/'release',dest)
p=dest/'interface.R';s=p.read_text();old='''  values<-setNames(lapply(seq_len(nrow(t)),function(i)env$o2prof_transform(z[[t$coordinate[i]]],t[i,,drop=FALSE],inverse=TRUE)),t$parameter)
  cg_candidate(design,values)'''
new='''  # The original helper owns endpoint-roundoff semantics; its base is independent neutral.
  t$status<-'active';t$fixed_value<-NA_real_
  map<-list(table=t,base_params=c(design$fixed,cg_neutral(env,t)))
  values<-env$o2prof_from_free(map,z[t$coordinate])
  cg_candidate(design,values[cg_roles])'''
assert old in s;s=s.replace(old,new);p.write_text(s)
p=dest/'worker.R';s=p.read_text().replace("task<-args[2];out<-file.path(study,'verification',task)","task_id<-args[2];task<-sub('_R1$','',task_id);out<-file.path(study,'verification',task_id)")
s=s.replace("'release/","'repair_release/").replace('verification/I01/','verification/I01_R1/')
needle="   neutral<-cg_neutral(env,design$bounds)\n   for(n in cg_roles)"
insert="""   for(which_bound in c('lower','upper')) {
     endpoint<-setNames(env$o2prof_transform(t[[which_bound]],t),t$coordinate)
     back<-cg_from_free(env,design,endpoint)
     check(paste0('all20_endpoint_roundtrip_',which_bound),all(vapply(seq_len(nrow(t)),function(i){v<-back[[t$parameter[i]]];v>=t$lower[i]&&v<=t$upper[i]&&abs(v-t[[which_bound]][i])<=64*.Machine$double.eps*max(1,abs(t$lower[i]),abs(t$upper[i]))},logical(1))))
   }
   neutral<-cg_neutral(env,design$bounds)
   for(n in cg_roles)"""
assert needle in s;s=s.replace(needle,insert);p.write_text(s)
p=dest/'read_guard.c';s=p.read_text().replace('"I01","P01"','"I01","I01_R1","P01"').replace('i<6','i<7');p.write_text(s)
p=dest/'runner.sh';s=p.read_text().replace('$study/release/','$study/repair_release/').replace('$study/verification/I01/','$study/verification/I01_R1/');s=s.replace('--bind "$study/inputs:', '--bind "$study/repair_release:$study/repair_release:ro" --bind "$study/inputs:');p.write_text(s)
p=dest/'verify_source.py';p.write_text('''from pathlib import Path
import sys
from common import verify,read,sha
s=Path(sys.argv[1]);verify(s)
for r in read(s/'repair_manifest.tsv'):
 p=s/r['path'];assert sha(p)==r['sha256'] and p.stat().st_size==int(r['bytes'])
print('BASE_AND_REPAIR_SOURCE_PASS')
''')
checks=[]
for p in dest.glob('*.py'):ast.parse(p.read_text());checks.append(dict(check='syntax_'+p.name,passed=True))
for p in dest.glob('*.R'):
 subprocess.run([str(project/'scripts/agentRrunner.sh'),'-e','invisible(parse(file='+repr(str(p))+'))'],check=True,stdout=subprocess.DEVNULL);checks.append(dict(check='R_parse_'+p.name,passed=True))
subprocess.run(['bash','-n',str(dest/'runner.sh')],check=True);checks.append(dict(check='shell_parse',passed=True))
assert 'env$o2prof_from_free(map,z[t$coordinate])' in (dest/'interface.R').read_text();checks.append(dict(check='original_helper_endpoint_rule',passed=True))
assert 'all20_endpoint_roundtrip_' in (dest/'worker.R').read_text();checks.append(dict(check='both_endpoints_regression_registered',passed=True))
table(study/'reports/repair_static_tests.tsv',checks)
table(study/'repair_manifest.tsv',[dict(path=str(p.relative_to(study)),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(dest.iterdir()) if p.is_file()])
table(study/'accounting/retry_ledger.tsv',[dict(key='used',value=1),dict(key='maximum',value=1)])
table(study/'accounting/technical_retry.tsv',[dict(ordinal=1,failed_task='I01',failed_job='21986841',repair_task='I01_R1',reason='Strict natural-bound check rejected log-transformed endpoint roundoff; new adapter had omitted original helper projection semantics',repair='Reuse source-identical o2prof_from_free with an independent neutral public base; register exact lower/upper regressions',scientific_question_changed=False,prior_forward_calls=0,prior_observations=0)])
(study/'reports/technical_repair.md').write_text('''# Single interface repair

I01 failed before propagation or generation. The new adapter independently inverted transforms then enforced exact natural bounds, unlike the existing validated helper. The repair reuses unchanged o2prof_from_free, including64 machine-epsilon endpoint-roundoff tolerance and projection onto exact declared endpoints. It does not enlarge the domain, replace zero with epsilon or project materially invalid proposals. All20 lower/upper endpoint regressions are added. The old immutable release, trace, checks and failed job remain. One global technical retry is consumed; any second technical failure stops dependent work. No scientific point or generating truth is changed.
''')
print('INTERFACE_REPAIR_STATIC_PASS',len(checks))
