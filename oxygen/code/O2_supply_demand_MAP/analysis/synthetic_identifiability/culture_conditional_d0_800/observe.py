"""Read-only local progress observer. HPC dispatch does not depend on this process."""
import csv,io,os,re,subprocess,time
ROOT='/share/lab_crd/taoli/Project/StructuralAndPracticalIdentifiability/oxygen/results/synthetic_identifiability/20261113_culture_conditional_d0_seed144_02'
ROOT=os.environ.get('D0_OBSERVE_ROOT',ROOT)
REMOTE='s='+ROOT+'; if test -f "$s/COMPLETE"; then echo D0_COMPLETE; fi; if test -f "$s/tables/STOP.tsv"; then echo D0_STOP; cat "$s/tables/STOP.tsv"; fi; cat "$s/accounting/jobs.tsv"; for f in "$s"/logs/global_[A-F]_[0-9]*.log; do test -f "$f" || continue; printf "%s " "${f##*/}"; tail -n 1 "$f"; done'
COMMAND=['ssh','-o','BatchMode=yes','-o','IdentitiesOnly=yes','-o','ControlMaster=auto','-o','ControlPersist=300','-o','ControlPath=/private/tmp/culture_d0_seed144_02_%C','-o','ConnectTimeout=15','-o','ServerAliveInterval=10','-o','ServerAliveCountMax=2','-i','/Users/4482173/.ssh/id_ed25519','4482173@red.moffitt.org',REMOTE]
def main():
    last_jobs=None;last_progress={};failures=0
    print('Read-only observer PID',os.getpid(),'; explicit SSH identity, reused connection',flush=True)
    for tick in range(4000):
        try:
            r=subprocess.run(COMMAND,stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,timeout=35)
            if r.returncode:raise RuntimeError(r.stderr.strip())
            lines=r.stdout.splitlines();header=next(i for i,l in enumerate(lines) if l.startswith('task\tphase\t'))
            data=[lines[header]]+[l for l in lines[header+1:] if l.count('\t')==len(lines[header].split('\t'))-1]
            rows=list(csv.DictReader(io.StringIO('\n'.join(data)),delimiter='\t'));key=tuple((j['task'],j['state'],j['exit_code']) for j in rows)
            if key!=last_jobs:
                active=[j for j in rows if j['state'].split('+')[0] not in {'COMPLETED','FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED'}]
                print(time.strftime('%Y-%m-%d %H:%M:%S %Z'),'Owned parents',len(rows),'; active:', '; '.join('{}={}({})'.format(j['task'],j['state'],j['job_id']) for j in active) or 'none','; terminal failures:',sum(j['state']=='FAILED' for j in rows),flush=True);last_jobs=key
            for line in lines:
                m=re.search(r'^(global_[A-F]_\d+[.]log) Iteration:\s*(\d+) bestvalit:\s*([0-9.eE+-]+)',line)
                if m:
                    name,generation,score=m.groups()
                    if last_progress.get(name)!=(generation,score):
                        print(time.strftime('%Y-%m-%d %H:%M:%S %Z'),name,'generation='+generation,'objective='+score,flush=True);last_progress[name]=(generation,score)
            markers=[x for x in lines if x in ('D0_COMPLETE','D0_STOP')]
            if markers:print('OBSERVER_TERMINAL_MARKER',*markers,flush=True);break
            failures=0
        except Exception as e:
            failures+=1;print('Read-only connectivity warning',failures,str(e),flush=True)
            if failures>=3:print('OBSERVER_HOLD; HPC controller remains independent',flush=True);break
        time.sleep(55)
if __name__=='__main__':main()
