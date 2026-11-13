"""Audit immutable profile dispatch ancestry; no Slurm/model invocation."""
import argparse
from pathlib import Path
import workflow as w
from profile_windows import verify_branch_source

def main():
 p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();s=a.study.resolve()
 ep=w.kv(s/'tables/profile_epoch.tsv');jobs=w.read_table(s/'tables/jobs.tsv');by={j['task']:j for j in jobs};rows=[]
 for j in jobs:
  if j['kind']!='profile' or not j['task'].startswith('p_'+ep['epoch']+'_'):continue
  d=s/'dispatch'/j['task'];w.verify_manifest(d,d/'manifest.tsv')
  if w.sha(d/'manifest.tsv')!=j['dispatch_sha256'] or w.sha(d/'spec.tsv')!=j['spec_sha256']:raise ValueError('Profile dispatch identity changed')
  options=w.kv(d/'spec.tsv');verify_branch_source(s,ep['epoch'],j['task'],options,by)
  rows.append({'task':j['task'],'epoch':ep['epoch'],'pass':'true','dispatch_sha256':j['dispatch_sha256'],'model_evaluations':0})
 w.table(s/'tables/profile_lineage_audit.tsv',rows,['task','epoch','pass','dispatch_sha256','model_evaluations'])
 print('Verified',len(rows),'profile dispatch histories; zero model calls.')
if __name__=='__main__':main()
