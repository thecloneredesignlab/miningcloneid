"""Prepare/verify a thin English deliverable set after resident completion."""
from pathlib import Path
import argparse,os,sys
from common import read,table,sha,verify

def select(s):
 s=Path(s).resolve();assert (s/'controller/FINISHED').exists(),'Do not collect an incomplete study as final'
 assert read(s/'reports/final_gate.tsv')[0]['active_owned_jobs']=='0';verify(s)
 inventory=read(s/'reports/authoritative_hpc_inventory.tsv');allowed=[]
 def keep(rel):
  p=Path(rel);parts=p.parts
  if any(x in parts for x in ['tmp','members','improvements','failures']):return False
  if p.name.startswith('trace.') or p.name=='evaluations.tsv':return False
  if parts[0] in ['inputs','runtime','release','control_release_v2','control_release_v3','assessment_release','testing_release']:return True
  if parts[0]!='verification':return True
  if p.suffix in ['.tsv','.txt','.log','.pdf'] or p.name=='read_guard.so':return True
  if parts[1]=='T00' and p.suffix=='.rds':return True
  if parts[1].startswith('D_') and p.suffix=='.rds':return True
  if parts[1].startswith('F_') and p.name in ['bank.rds','requested_population.rds','actual_initial_population.rds','native_DE_result.rds']:return True
  return False
 for r in inventory:
  if keep(r['path']):allowed.append(dict(path=r['path'],bytes=r['bytes'],sha256=r['sha256']))
 # Final metadata/report files are produced after the scientific inventory is sealed.
 extras=[]
 # C1 seals the scientific inventory before its own receipts exist. Include
 # those finished receipts, scheduler logs and the collector's bound source.
 for folder in ['accounting','reports','controller','testing_release','verification/C1','logs','collect_release','collect_release_v2']:
  for p in sorted((s/folder).rglob('*')):
   if p.is_file() and p.name not in {'owner.lock'} and (p.suffix in {'.tsv','.md','.txt','.log','.py'} or p.name in {'FINISHED','SELECTION_LOCKED'}):extras.append(p)
 for p in s.iterdir():
  if p.is_file() and p.suffix in {'.tsv','.md','.txt'}:extras.append(p)
 known={r['path'] for r in allowed}
 for p in extras:
  rel=str(p.relative_to(s))
  if rel in known or p.name in ['local_collection_selection.tsv','artifact_locations.tsv','local_collection_manifest.tsv','local_collection_files.txt']:continue
  known.add(rel);allowed.append(dict(path=rel,bytes=p.stat().st_size,sha256=sha(p)))
 allowed.sort(key=lambda r:r['path']);table(s/'reports/local_collection_selection.tsv',allowed)
 (s/'reports/local_collection_files.txt').write_text('\n'.join(r['path'] for r in allowed)+'\n')
 selected={r['path'] for r in allowed};locations={r['path']:r for r in inventory}
 locations.update({r['path']:r for r in allowed if r['path'] not in locations})
 table(s/'reports/artifact_locations.tsv',[dict(path=r['path'],bytes=r['bytes'],sha256=r['sha256'],HPC_authority='FULL_SCIENTIFIC_ARCHIVE',local_scope='PLANNED_THIN_COPY' if r['path'] in selected else'HPC_ONLY') for r in sorted(locations.values(),key=lambda r:r['path'])])
 print('THIN_SELECTION_READY',len(allowed),'files',sum(int(r['bytes']) for r in allowed),'bytes; raw traces/failure snapshots/large constructor states remain on HPC')

def validate(s,repo):
 s=Path(s).resolve();repo=Path(repo).resolve();records=read(s/'reports/local_collection_selection.tsv')
 for r in records:
  p=s/r['path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],'Local subset mismatch '+r['path']
 for r in read(s/'accounting/preservation.tsv'):assert sha(repo/r['path'])==r['sha256'],'Local pre-existing content changed '+r['path']
 table(s/'local_collection_manifest.tsv',records)
 selected={r['path'] for r in records};locations=read(s/'reports/artifact_locations.tsv')
 for r in locations:r['local_scope']='VERIFIED_LOCAL_COPY' if r['path'] in selected else'HPC_ONLY'
 table(s/'reports/artifact_locations.tsv',locations)
 table(s/'reports/local_collection_gate.tsv',[dict(status='PASS',files=len(records),bytes=sum(int(r['bytes']) for r in records),raw_authority='HPC',preserved_local_tracked_files=len(read(s/'accounting/preservation.tsv')),scope='Verified thin subset; no claim that omitted raw files are locally present')])
 print('LOCAL_THIN_COLLECTION_PASS',len(records))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['select','validate']);p.add_argument('study');p.add_argument('--repo');a=p.parse_args()
 if a.mode=='select':select(a.study)
 else:assert a.repo;validate(a.study,a.repo)
