"""Component contrasts of verified exported points; no model or Slurm calls."""
from pathlib import Path
import argparse,math
import workflow as w

def diagnose(s):
    s=Path(s);source=s/'tables/selected_point_components.tsv'
    points=w.read_table(source);jobs=w.read_table(s/'tables/jobs.tsv');families={j['root_id']:j['root_family'] for j in jobs}
    rows=[]
    for policy in ('production','data_only'):
        baseline=next(r for r in points if r['point']=='baseline_'+policy)
        uniform=min((r for r in points if r['policy']==policy and r['kind']=='root_candidate' and families.get(r['root_id'])=='uniform'),key=lambda r:float(r['L']))
        gap=float(uniform['L'])-float(baseline['L'])
        deltas={c:float(uniform[c])-float(baseline[c]) for c in ('L_b','L_p','L_n','L_prior')}
        if not math.isclose(sum(deltas.values()),gap,rel_tol=0,abs_tol=1e-8):raise ValueError('Component/objective contrast mismatch')
        rows.append(dict(policy=policy,baseline_root=baseline['root_id'],uniform_root=uniform['root_id'],
            baseline_objective=baseline['L'],uniform_objective=uniform['L'],actual_minimum_gap=gap,
            burden_difference=deltas['L_b'],ploidy_difference=deltas['L_p'],necrosis_difference=deltas['L_n'],
            weighted_prior_difference=deltas['L_prior'],source_components_sha256=w.sha(source),model_evaluations=0))
    w.table(s/'tables/recovery_component_differences.tsv',rows,list(rows[0]))
    return rows

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('study',type=Path)
    for row in diagnose(parser.parse_args().study.resolve()):print(row)
