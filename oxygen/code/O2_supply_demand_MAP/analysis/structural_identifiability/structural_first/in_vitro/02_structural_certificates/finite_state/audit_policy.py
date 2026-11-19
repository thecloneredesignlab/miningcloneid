"""Metadata inspection of the traced scientific worker; no CAS imports."""
from pathlib import Path
import re
import sys
from protocol import read,table

folder=Path(sys.argv[1]);kind=sys.argv[2];bad=[];events=0
traces=list(folder.glob('trace*'))
if not traces:bad.append('No scientific-worker process/file trace')
for p in folder.glob('blocked*'):
    if p.stat().st_size:bad.append('Blocked activity: '+p.name)
for p in traces:
    for line in p.read_text(errors='replace').splitlines():
        m=re.search(r'(?:openat\([^,]+,\s*"|open\("|execve\(")([^"\n]+)"',line)
        result=re.search(r'=\s*(-?\d+)\s*$',line)
        if m and result and int(result.group(1))>=0:
            events+=1;path=m.group(1).lower();base=Path(path).name
            if '.json' in path or '/sys/java/' in path or any(s in path for s in ('libjvm','libjava.','libv8','libmozjs')) or 'json' in base and '.so' in base or base in ('java','node','nodejs','js') or path.endswith(('.js','.java','.jar')):
                bad.append(line)
if kind=='R':
    inv=folder/'loaded_namespaces.tsv'
    if not inv.exists():bad.append('Missing completed R namespace inventory')
    else:
        for r in read(inv):
            if r['namespace'] in ('rJava','V8','jsonlite','rjson','RJSONIO'):bad.append(r['namespace'])
else:
    inv=folder/'imported_modules.tsv'
    if not inv.exists():bad.append('Missing completed Python module inventory')
    else:
        for r in read(inv):
            if r['module'].split('.')[0] in ('json','_json','simplejson','orjson','ujson','jpype','jnius','execjs','py_mini_racer','symbtools','strike_goldd'):bad.append(r['module'])
table(folder/'policy_gate.tsv',[dict(status='FAIL' if bad else 'PASS',trace_files=len(traces),successful_events=events,coverage='scientific worker and children; file/process trace and language module inventory; container launcher excluded',findings='; '.join(bad))])
print('RUNTIME_POLICY_GATE','FAIL' if bad else 'PASS',events)
if bad:raise SystemExit(77)
