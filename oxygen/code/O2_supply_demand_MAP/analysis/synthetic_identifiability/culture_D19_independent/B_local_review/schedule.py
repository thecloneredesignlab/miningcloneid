"""Finite Slurm reservations and application gates for the B-only supplement."""
import fcntl
from pathlib import Path
import subprocess
import sys

from common import read, sha, table, verify

LIMIT = 7200
TERMINAL = {'COMPLETED','FAILED','TIMEOUT','CANCELLED','NODE_FAIL','OUT_OF_MEMORY',
            'PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED','SPECIAL_EXIT'}
INFRASTRUCTURE = {'NODE_FAIL','PREEMPTED','BOOT_FAIL'}
STAGES = ('admit','poll','replayA','replayB')
RESERVATION = dict(admit=900,poll=2700,replayA=900,replayB=900,retry=1800)
FIELD = ('task','stage','attempt','job_id','state','exit_code','settled_seconds','live_reserved_seconds','charged_calls')

class Schedule:
    def __init__(self, root):
        self.root=Path(root)
        self.ledger=self.root/'accounting/jobs.tsv'
        self.rows=read(self.ledger) if self.ledger.exists() else []

    def save(self):
        table(self.ledger,self.rows,FIELD)

    def folder(self,row):
        return self.root/'verification'/row['task']

    def refresh(self):
        if not self.rows:
            return
        ids=[r['job_id'] for r in self.rows]
        if any(not i.isdigit() for i in ids):
            raise RuntimeError('Unreconciled sbatch reservation')
        output=subprocess.check_output(['sacct','-nP','-j',','.join(ids),
            '--format=JobIDRaw,State%40,ExitCode,AllocCPUS,ElapsedRaw,QOS'],text=True)
        (self.root/'accounting/slurm_snapshot.txt').write_text(output)
        records={v[0]:v for line in output.splitlines() if len(v:=line.split('|'))>=6}
        for row in self.rows:
            jid=row['job_id']
            if jid not in records:
                continue
            v=records[jid]
            if v[5] not in ('xxlarge',''): raise RuntimeError('Unexpected QoS')
            state=v[1].split()[0]
            row.update(state=state,exit_code=v[2])
            row['charged_calls']=str(len(list((self.folder(row)/'starts').glob('*.rds'))))
            if state in TERMINAL:
                actual=max(int(x[3] or 0)*int(x[4] or 0) for k,x in records.items()
                           if k==jid or k.startswith(jid+'.'))
                row.update(settled_seconds=str(actual),live_reserved_seconds='0')
        self.save()
        self.check_budget()

    def check_budget(self):
        settled=sum(int(r['settled_seconds']) for r in self.rows)
        reserved=sum(int(r['live_reserved_seconds']) for r in self.rows)
        calls=sum(int(r['charged_calls']) for r in self.rows)
        if settled+reserved>LIMIT or calls>155:
            raise RuntimeError('Supplement ceiling exceeded')
        return settled,reserved,calls

    def latest(self, stage):
        return next((r for r in reversed(self.rows) if r['stage']==stage),None)

    def valid(self, stage):
        row=self.latest(stage)
        if not row or (row['state'],row['exit_code'])!=('COMPLETED','0:0'):
            return False
        p=self.folder(row)
        manifest=p/'phase_manifest.tsv'
        if not manifest.is_file():
            return False
        for m in read(manifest):
            f=p/m['path']
            if f.stat().st_size!=int(m['bytes']) or sha(f)!=m['sha256']:
                raise RuntimeError('Phase receipt digest mismatch: '+str(f))
        app,=read(p/'application_result.tsv')
        policy,=read(p/'policy_gate.tsv')
        codes,=read(p/'phase_exit_codes.tsv')
        if app['job']!=row['job_id'] or app['status']!='PASS' or policy['status']!='PASS' or any(int(k) for k in codes.values()):
            raise RuntimeError('Independent application/source/policy gate failed: '+stage)
        denied=p/'read_denials.tsv'
        if denied.exists() and denied.read_text().strip():
            raise RuntimeError('Unexpected scientific read denial')
        return True

    def submit(self, stage, retry=False):
        self.refresh()
        verify(self.root)
        if stage not in STAGES:
            raise RuntimeError('Unregistered stage')
        previous=self.latest(stage)
        if retry:
            if not previous or previous['state'] not in INFRASTRUCTURE or any(r['attempt']=='2' for r in self.rows):
                raise RuntimeError('Retry requires a diagnosed first infrastructure loss')
            f=self.folder(previous)
            starts=list((f/'starts').glob('*.rds'))
            records=list((f/'points').glob('*.rds'))
            # Partial evaluations require sealed policy receipts to become eligible checkpoints.
            if records:
                policy=f/'policy_gate.tsv'
                manifest=f/'phase_manifest.tsv'
                if stage!='poll' or not policy.exists() or not manifest.exists() or read(policy)[0]['status']!='PASS':
                    raise RuntimeError('Partial checkpoint lacks a valid sealed policy audit')
                sealed={m['path']:m for m in read(manifest)}
                for info in records+starts+[policy]:
                    m=sealed[info.relative_to(f).as_posix()]
                    if info.stat().st_size!=int(m['bytes']) or sha(info)!=m['sha256']:
                        raise RuntimeError('Invalid partial checkpoint')
            if len(starts)-len(records)>2:
                raise RuntimeError('Unfinished attempts exceed retry allowance')
        elif previous:
            raise RuntimeError('Stage already booked')
        if stage!='admit' and not self.valid('admit'):
            raise RuntimeError('Admission gate has not passed')
        if stage.startswith('replay'):
            if not self.valid('poll'):
                raise RuntimeError('Poll has not passed')
            poll,=read(self.folder(self.latest('poll'))/'local_resolution_audit.tsv')
            if poll['lower_found']!='TRUE': raise RuntimeError('No lower point requires replay')
        cap=RESERVATION['retry' if retry else stage]
        used,reserved,calls=self.check_budget()
        if used+reserved+cap>LIMIT or calls>=155:
            raise RuntimeError('New work exceeds the independent cap')
        task='D_B_'+stage+('_retry' if retry else '')
        row=dict(task=task,stage=stage,attempt='2' if retry else '1',job_id='submission_pending',
                 state='RESERVED',exit_code='',settled_seconds='0',live_reserved_seconds=str(cap),charged_calls='0')
        self.rows.append(row)
        self.save()
        args=['sbatch','--parsable','--account=tao.li2','--qos=xxlarge','--ntasks=1','--cpus-per-task=1',
              '--hint=nomultithread','--mem=16G','--no-requeue','--time='+str(cap//60),
              '--job-name=d19B_local_'+stage,'--output='+str(self.root/'logs'/(task+'_%j.txt')),
              str(self.root/'release/launch.sh'),str(self.root),task]
        jid=subprocess.check_output(args,text=True).strip().split(';')[0]
        if not jid.isdigit():
            raise RuntimeError('Uncertain job ID; native reconciliation required')
        row.update(job_id=jid,state='SUBMITTED')
        self.save()
        print('SUBMITTED',task,jid)

    def status(self):
        self.refresh()
        for row in self.rows:
            stage=row['stage']
            science='PASS' if self.latest(stage) is row and self.valid(stage) else 'PENDING_OR_FAILED'
            print(*(row[k] for k in FIELD),science,sep='\t')
        print('SETTLED_RESERVED_AND_CALLS',*self.check_budget(),sep='\t')

def main(root,action,*args):
    study=Schedule(root)
    with (study.root/'accounting/owner.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if action=='status':study.status()
        elif action=='submit':study.submit(args[0])
        elif action=='retry':study.submit(args[0],retry=True)
        else:raise RuntimeError('Unknown scheduler action')

if __name__=='__main__':
    main(*sys.argv[1:])
