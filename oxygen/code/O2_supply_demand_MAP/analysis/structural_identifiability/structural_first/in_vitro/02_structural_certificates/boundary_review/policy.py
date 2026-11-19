"""Trace/inventory acceptance, independently of scientific application exit."""
from pathlib import Path
import re
import sys
from common import read,table

def audit(folder,kind):
    folder=Path(folder);bad=[];events=0;writes=0;terminal_denials=0;mutations=0
    traces=list(folder.glob('trace.*'))
    if not traces:bad.append('Missing scientific file/process trace')
    for p in folder.glob('blocked*'):
        if p.stat().st_size:
            for line in p.read_text().splitlines():
                if kind=='R' and line=='WRITE_OUTSIDE_PHASE\topen\t/dev/tty':terminal_denials+=1
                else:bad.append('Blocked persistent/excluded activity: '+p.name+': '+line)
    for p in traces:
        for line in p.read_text(errors='replace').splitlines():
            result=re.search(r'=\s*(-?\d+)\s*$',line)
            if result and int(result.group(1))>=0 and re.match(r'(?:mkdir|mkdirat|rmdir|unlink|unlinkat|rename|renameat|renameat2)\(',line):
                mutations+=1
                for target in re.findall(r'"([^"\n]+)"',line):
                    pp=Path(target)
                    if not pp.is_absolute() or not (pp==folder or folder in pp.parents):bad.append('Mutation outside phase: '+line)
            m=re.search(r'(?:openat\([^,]+,\s*"|open\("|execve\(")([^"\n]+)"',line)
            if not m or not result or int(result.group(1))<0:continue
            events+=1;path=m.group(1);lower=path.lower();base=Path(lower).name
            if '.json' in lower or any(s in lower for s in ('libjvm','libjava.','libv8','libmozjs','libjson')) or base in ('java','node','nodejs','js') or lower.endswith(('.js','.java','.jar')):bad.append('Forbidden successful event: '+line)
            if any(f in line for f in ('O_WRONLY','O_RDWR','O_CREAT','O_TRUNC','O_APPEND')):
                writes+=1;pout=Path(path)
                if not pout.is_absolute() or not (pout==Path('/dev/null') or pout==folder or folder in pout.parents):bad.append('Successful write outside phase: '+line)
    inv=folder/('loaded_namespaces.tsv' if kind=='R' else 'imported_modules.tsv')
    if not inv.exists():bad.append('Missing completed module/namespace inventory')
    else:
        for r in read(inv):
            name=r['namespace'] if kind=='R' else r['module'].split('.')[0]
            if name in {'json','_json','simplejson','orjson','ujson','jpype','jnius','execjs','py_mini_racer','symbtools','strike_goldd','rJava','V8','jsonlite','rjson','RJSONIO'}:bad.append('Excluded inventory entry: '+name)
    record=dict(status='FAIL' if bad else 'PASS',trace_files=len(traces),successful_events=events,successful_write_events=writes,
        successful_mutation_events=mutations,denied_terminal_probes=terminal_denials,
        coverage='scientific worker/children only; container and scheduler launchers excluded',findings='; '.join(bad))
    table(folder/'policy_gate.tsv',[record]);return record

if __name__=='__main__':
    row=audit(sys.argv[1],sys.argv[2]);print('SCIENTIFIC_POLICY',row['status'],row['successful_events'])
    if row['status']!='PASS':raise SystemExit(77)
