#!/usr/bin/env python3
"""One bounded pass for the already recorded data-only reserve activation."""
import argparse
from pathlib import Path
import workflow as w
import advance

SELECTED=('b_data_only_reserve_global','b_data_only_reserve_probe_local',
          'b_data_only_reserve_g_local','b_data_only_reserve_probe_fresh',
          'b_data_only_reserve_g_fresh','b_data_only_reserve_probe_directions',
          'b_data_only_reserve_g_directions')

def main():
    p=argparse.ArgumentParser();p.add_argument('study',type=Path);a=p.parse_args();s=a.study.resolve()
    if s.name!='20261111_invivo_identifiability_formal_010554Z': raise ValueError('Unexpected study')
    activation=s/'reports/reserve_activation.md'
    if not activation.exists(): raise ValueError('No recorded reserve activation')
    rows=advance.snapshot(s);by={j['task']:j for j in rows}
    primary=[j for j in rows if not j['role'].startswith('reserve_')]
    if len(primary)!=37 or not all(j['state']=='COMPLETED' and j['exit_code']=='0:0' for j in primary): raise ValueError('Primary evidence not fully closed')
    if not all(j['task'] in SELECTED for j in rows if j['role'].startswith('reserve_')): raise ValueError('Different reserve activation already exists')
    if any(w.closed(j) and (j['state']!='COMPLETED' or j['exit_code']!='0:0') for j in rows): raise ValueError('Unexpected allocation failure')
    plans={r['task']:r for r in w.read_table(s/'manifests/task_plan.tsv')}
    for task in SELECTED:
        if task in by: continue
        plan=plans[task]
        parents=[plan[k] for k in ('source_task','required_check') if plan[k]]
        if all(t in by and by[t]['state']=='COMPLETED' and by[t]['exit_code']=='0:0' for t in parents):
            w.submit(s,task);by={j['task']:j for j in w.State(s).load()}
    rows=advance.snapshot(s);advance.capture_start_snapshots(s,rows);advance.publish_status(s,rows)
    print('Selected reserve pass complete; no loop or alternative activation.')
    print('Formal actual-plus-live:',sum(map(w.charge,rows)),'core-seconds; lifetime:',30491+sum(map(w.charge,rows)))
    for j in rows:
        if not w.closed(j): print('ACTIVE',j['task'],j['job_id'],j['state'])

if __name__=='__main__': main()
