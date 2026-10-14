#!/usr/bin/env python3
"""Temporarily reprioritize only pending model tasks; restore after global allocation."""
import argparse
import datetime
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def command(args):
    return subprocess.check_output(args,universal_newlines=True,timeout=45).strip()


def indices(job):
    parent,part=job.split('_',1)
    part=part.strip('[]')
    result=set()
    for token in part.split(','):
        if not re.fullmatch(r'\d+(-\d+)?',token): raise ValueError('Unexpected array selector')
        pair=token.split('-'); lo=int(pair[0]); hi=int(pair[-1])
        result.update(range(lo,hi+1))
    return parent,result


def selector(parent,values):
    values=sorted(values); ranges=[]; start=previous=values[0]
    for value in values[1:]:
        if value==previous+1: previous=value; continue
        ranges.append(str(start) if start==previous else '%d-%d'%(start,previous))
        start=previous=value
    ranges.append(str(start) if start==previous else '%d-%d'%(start,previous))
    return parent+'_['+','.join(ranges)+']'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--global-job',required=True)
    parser.add_argument('--array-ids',required=True)
    parser.add_argument('--out-dir',required=True)
    parser.add_argument('--max-wait',type=int,default=600)
    args=parser.parse_args(); os.environ['SLURM_BITSTR_LEN']='0'
    parents=args.array_ids.split(',')
    if any(not p.isdigit() for p in parents): raise ValueError('Numeric array IDs required')
    out=Path(args.out_dir); out.mkdir(parents=True,exist_ok=True)
    lock=open(str(out/'priority_adjustment.lock'),'a'); fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    path=out/'priority_adjustment.json'
    if path.exists(): raise ValueError('Existing adjustment record; inspect before repeating')
    record=dict(created_utc=now(),global_job=args.global_job,array_ids=parents,
                controller_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                status='preparing',changes=[],restore_results=[],running_tasks_targeted=False)
    def save():
        tmp=path.with_suffix('.json.tmp'); tmp.write_text(json.dumps(record,indent=2)+'\n')
        os.replace(str(tmp),str(path))
    def stopped(signum,frame): raise RuntimeError('Controller received signal '+str(signum))
    signal.signal(signal.SIGTERM,stopped); signal.signal(signal.SIGINT,stopped)
    save()
    if not command(['squeue','-h','-j',args.global_job,'-o','%T']).startswith('PENDING'):
        record['status']='global_already_allocated_or_terminal'; save(); return
    try:
        lines=command(['squeue','-h','-t','PD','-j',args.array_ids,'-o','%i|%y|%j']).splitlines()
        for line in lines:
            job,nice,name=line.split('|',2); parent,tasks=indices(job)
            if parent not in parents or name not in ['sobol_full_evaluate','sobol_full_evaluate_retry']:
                raise ValueError('Unexpected pending job scope')
            if min(tasks)<=150 and name=='sobol_full_evaluate':
                raise ValueError('This adjustment must target neighborhood workers only')
            item=dict(selector=job,parent=parent,task_ids=sorted(tasks),original_nice=int(nice),
                      temporary_nice=max(int(nice),100000),captured_state='PENDING',requested_utc=now())
            record['changes'].append(item); save()
            command(['scontrol','update','JobId='+job,'Nice='+str(item['temporary_nice'])])
            item['applied_utc']=now(); save()
        record['status']='waiting_for_global_allocation'; save()
        deadline=time.monotonic()+args.max_wait
        while time.monotonic()<deadline:
            state=command(['squeue','-h','-j',args.global_job,'-o','%T|%R|%N'])
            record['last_global_state']=state; record['last_check_utc']=now(); save()
            if not state.startswith('PENDING'):
                record['global_allocation_or_terminal_utc']=now(); break
            time.sleep(10)
        else:
            record['wait_limit_reached']=True
    except Exception as error:
        record['error']=str(error)
        raise
    finally:
        record['status']='restoring'; save()
        restore_errors=[]
        for item in record['changes']:
            try:
                live=command(['squeue','-r','-h','-j',item['parent'],'-o','%i'])
                ids=set()
                for job in live.splitlines():
                    parent,tasks=indices(job); ids.update(tasks)
                remain=ids.intersection(item['task_ids'])
                if remain:
                    job=selector(item['parent'],remain)
                    command(['scontrol','update','JobId='+job,'Nice='+str(item['original_nice'])])
                record['restore_results'].append(dict(parent=item['parent'],restored_live_tasks=len(remain),
                                                       original_nice=item['original_nice'],restored_utc=now()))
            except Exception as error:
                restore_errors.append(dict(parent=item['parent'],error=str(error)))
            save()
        record['restore_errors']=restore_errors
        record['status']='restore_failed' if restore_errors else 'restored'
        record['finished_utc']=now(); save()
        print(json.dumps({k:v for k,v in record.items() if k!='changes'},indent=2),flush=True)
        if restore_errors: raise RuntimeError('Nice restoration requires attention')


if __name__=='__main__': main()
