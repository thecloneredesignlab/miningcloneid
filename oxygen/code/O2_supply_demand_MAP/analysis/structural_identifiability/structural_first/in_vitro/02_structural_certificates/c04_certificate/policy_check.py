from pathlib import Path
import re
import sys
from protocol import read,table
folder=Path(sys.argv[1]);bad=[]
for p in folder.glob('blocked*'):
    if p.stat().st_size:bad.append(p.name)
for p in folder.glob('trace*'):
    for line in p.read_text().splitlines():
        m=re.search(r'(?:openat\([^,]+,\s*"|open\("|execve\(")([^"\n]+)"',line);r=re.search(r'=\s*(-?\d+)\s*$',line)
        if m and r and int(r.group(1))>=0:
            s=m.group(1).lower();base=Path(s).name
            if '.json' in s or '/sys/java/' in s or any(v in s for v in ('libjvm','libjava.','libv8','libmozjs')) or ('json' in base and '.so' in base) or base in ('java','node','nodejs','js') or s.endswith('.js'):bad.append(line)
for row in read(folder/'imported_modules.tsv'):
    if row['module'].split('.')[0] in ('symbtools','strike_goldd','native_strike','json','_json','simplejson','orjson','ujson','jpype','jnius','execjs','py_mini_racer'):bad.append(row['module'])
table(folder/'policy_gate.tsv',[dict(status='FAIL' if bad else 'PASS',findings='; '.join(bad))])
if bad:raise SystemExit(77)
