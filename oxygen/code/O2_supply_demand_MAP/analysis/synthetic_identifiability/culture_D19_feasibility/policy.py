"""Trace/inventory audit with inherited cwd and annotated directory descriptors."""
from pathlib import Path
import os
import re
import sys
from common import read,table

def audit(folder,kind):
    folder=Path(folder);bad=[];events=0;writes=0;mutations=0;terminal=0
    traces={p.name.split('.')[-1]:p for p in folder.glob('trace.*')}
    if not traces:bad.append('Missing scientific trace')
    for p in folder.glob('blocked*'):
        for line in p.read_text().splitlines():
            if kind=='R' and line=='WRITE_OUTSIDE_PHASE\topen\t/dev/tty':terminal+=1
            else:bad.append('Unexpected denied access: '+line)
    children=set()
    for p in traces.values():
        for line in p.read_text(errors='replace').splitlines():
            m=re.search(r'(?:clone3?|fork|vfork).*?=\s*(\d+)\s*$',line)
            if m:children.add(m.group(1))
    visited=set()
    def resolve(path,cwd,fd=None):
        if path.startswith('/'):return Path(os.path.normpath(path))
        if fd and fd!='AT_FDCWD':
            match=re.search(r'<([^>]+)>',fd)
            if not match:raise ValueError('Unresolved relative dirfd: '+fd)
            cwd=match.group(1)
        return Path(os.path.normpath(os.path.join(str(cwd),path)))
    def inside(path):return path==Path('/dev/null') or path==folder or folder in path.parents
    def scan(pid,initial):
        nonlocal events,writes,mutations
        if pid in visited or pid not in traces:return
        visited.add(pid);cwd=Path(initial)
        for line in traces[pid].read_text(errors='replace').splitlines():
            result=re.search(r'=\s*(-?\d+)(?:<[^>]*>)?\s*$',line)
            if not result or int(result.group(1))<0:continue
            ret=int(result.group(1))
            if re.search(r'(?:clone3?|fork|vfork)',line) and str(ret) in traces:scan(str(ret),cwd)
            cd=re.match(r'chdir\("([^"\n]+)"\)',line)
            if cd:cwd=resolve(cd.group(1),cwd)
            fcd=re.match(r'fchdir\(\d+<([^>]+)>\)',line)
            if fcd:cwd=Path(fcd.group(1))
            m=re.match(r'(openat|open|execve)\(',line)
            if m:
                pairs=re.findall(r'(AT_FDCWD|\d+(?:<[^>]+>)?),\s*"([^"\n]+)"',line)
                quoted=re.findall(r'"([^"\n]+)"',line)
                if not quoted:continue
                fd,path=pairs[0] if m.group(1)=='openat' and pairs else (None,quoted[0])
                try:p=resolve(path,cwd,fd)
                except ValueError as e:bad.append(str(e));continue
                events+=1;lower=str(p).lower();base=p.name.lower()
                if '.json' in lower or any(v in lower for v in ('libjvm','libjava.','libv8','libmozjs','libjson')) or base in ('java','node','nodejs','js') or lower.endswith(('.js','.java','.jar')):bad.append('Excluded successful event: '+line)
                if any(flag in line for flag in ('O_WRONLY','O_RDWR','O_CREAT','O_TRUNC','O_APPEND')):
                    writes+=1
                    if not inside(p):bad.append('Successful write outside phase: '+str(p))
            if re.match(r'(mkdir|mkdirat|rmdir|unlink|unlinkat|rename|renameat|renameat2)\(',line):
                mutations+=1
                pairs=re.findall(r'(AT_FDCWD|\d+(?:<[^>]+>)?),\s*"([^"\n]+)"',line)
                targets=[(None,p) for p in re.findall(r'"([^"\n]+)"',line)] if not pairs else pairs
                for fd,path in targets:
                    try:p=resolve(path,cwd,fd)
                    except ValueError as e:bad.append(str(e));continue
                    if not inside(p):bad.append('Mutation outside phase: '+str(p))
    for pid in traces:
        if pid not in children:scan(pid,folder)
    for pid in traces:
        if pid not in visited:bad.append('Unresolved process ancestry: '+pid)
    inv=folder/('loaded_namespaces.tsv' if kind=='R' else 'imported_modules.tsv')
    if not inv.exists():bad.append('Missing completed inventory')
    else:
        for r in read(inv):
            n=r['namespace'] if kind=='R' else r['module'].split('.')[0]
            if n in {'json','_json','simplejson','orjson','ujson','jpype','jnius','execjs','py_mini_racer','symbtools','strike_goldd','rJava','V8','jsonlite','rjson','RJSONIO'}:bad.append('Excluded inventory: '+n)
    record=dict(status='FAIL' if bad else 'PASS',trace_files=len(traces),successful_events=events,write_events=writes,mutation_events=mutations,denied_terminal_probes=terminal,coverage='scientific worker/children; initial cwd=phase output; fork/chdir and annotated dirfd audited; launchers excluded',findings='; '.join(bad))
    table(folder/'policy_gate.tsv',[record]);return record

if __name__=='__main__':
    r=audit(sys.argv[1],sys.argv[2]);print('POLICY',r['status'])
    if r['status']!='PASS':raise SystemExit(77)
