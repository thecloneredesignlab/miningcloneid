"""Finite constructor/global windows; tool limitation is not nonidentifiability."""
from pathlib import Path
import csv,os,signal,subprocess,sys,time
s=Path(sys.argv[1]);release=Path(sys.argv[2]);task=sys.argv[3];out=s/'results'/task;out.mkdir(exist_ok=True)
kind='conditional_k' if task=='global_conditional_k' else 'full_ten'
cmd=[str(s/'toolchain/julia-1.10.13/bin/julia'),'--startup-file=no','--history-file=no','--project='+str(s/'toolchain/environment'),str(release/'julia_case.jl'),str(out),kind]
started=time.monotonic();reason='completed';status='unresolved';peak_rss_kb=0
with (out/'native_stdout.txt').open('w') as stdout,(out/'native_stderr.txt').open('w') as stderr:
 p=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,start_new_session=True)
 while p.poll() is None:
  elapsed=time.monotonic()-started
  try:
   rss=[line for line in Path('/proc/'+str(p.pid)+'/status').read_text().splitlines() if line.startswith('VmRSS:')]
   if rss:peak_rss_kb=max(peak_rss_kb,int(rss[0].split()[1]))
  except FileNotFoundError:pass
  if peak_rss_kb>24*1024*1024:reason='registered_24_GiB_constructor_or_solver_soft_limit';break
  if not (out/'ODE_READY').exists() and elapsed>900:reason='registered_900_second_construction_limit';break
  if elapsed>17000:reason='registered_global_time_limit';break
  time.sleep(1)
 if p.poll() is None:
  os.killpg(p.pid,signal.SIGTERM)
  try:p.wait(timeout=20)
  except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
 if (out/'GLOBAL_COMPLETE').exists() and p.returncode==0:status='tool_result_available'
 elif (out/'construction_limit.tsv').exists():reason='registered_polynomial_term_limit'
 elif p.returncode:reason='tool_process_limit_or_failure_preserved'
with (out/'result.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['key','value']);w.writerows([['execution_valid','true'],['structural_status',status],['reason',reason],['native_exit_code',p.returncode],['peak_sampled_rss_kb',peak_rss_kb],['seconds',time.monotonic()-started],['condition',kind],['synthetic_datasets',0],['fit_calls',0],['step05_release','false']])
print(task,status,reason,p.returncode,flush=True)
