"""Read-only progress: signed closed-phase summaries plus active journals only."""
from pathlib import Path
import csv,math,sys,time
s=Path(sys.argv[1]).resolve();last=None

def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
while True:
 jobs=read(s/'accounting/jobs.tsv');parts=[]
 for origin in ['D_'+l for l in 'ABCDEF']+['F_E','F_F']:
  phases=[j for j in jobs if j['task'].startswith(origin+'_')]
  if not phases:continue
  total=0;lowest=math.inf
  for j in phases:
   p=s/'verification'/j['task'];result=p/'result.tsv'
   if j['state']=='COMPLETED' and result.exists():
    kv={r['key']:r['value'] for r in read(result)};total+=int(kv.get('calls',0))
    identity=p/'best_identity.tsv'
    if identity.exists():lowest=min(lowest,float(read(identity)[0]['D']))
    if origin.startswith('F_') and kv.get('status')=='BANK_READY':lowest=0
    continue
   journal=p/'evaluations.tsv'
   if journal.exists():
    with journal.open() as f:
     for r in csv.DictReader(f,delimiter='\t'):
      if r.get('event')=='START':total+=1
      elif r.get('event')=='END':
       try:lowest=min(lowest,float(r['score' if origin.startswith('D_') else 'F']))
       except (KeyError,TypeError,ValueError):pass
  parts.append(origin+':'+str(total)+','+format(lowest,'.6g')+','+phases[-1]['task'].split('_')[-1]+','+phases[-1]['state'])
 for t in ['C0','C1']:
  j=next((j for j in jobs if j['task']==t),None)
  if j:parts.append(t+':'+j['state'])
 snapshot=' | '.join(parts)
 if snapshot!=last:print(time.strftime('%Y-%m-%d %H:%M:%S'),snapshot,flush=True);last=snapshot
 for marker in ['STOP','FINISHED']:
  p=s/'controller'/marker
  if p.exists():print(marker,p.read_text().strip(),flush=True);raise SystemExit(0)
 if '--once' in sys.argv[2:]:break
 time.sleep(30)
