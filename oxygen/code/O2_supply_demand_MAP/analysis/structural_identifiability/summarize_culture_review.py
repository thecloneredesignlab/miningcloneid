#!/usr/bin/env python3
"""Summarize verified review candidates, independent recovery and owned costs."""
import argparse
import collections
import csv
import importlib.util
import io
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone


def table(path, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:path.write_text('');return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w')as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader();writer.writerows(rows)


def markdown(rows, fields):
    return '| '+' | '.join(fields)+' |\n|'+'|'.join('---'for _ in fields)+'|\n'+'\n'.join(
        '| '+' | '.join(str(r.get(k,''))for k in fields)+' |'for r in rows)


def accounting(module, study, state):
    selected={r['job_id']:dict(kind=r['kind'],key=key,receipt=r['receipt'])
              for key,r in state['reservations'].items()if r.get('job_id')}
    fields='JobID,JobName,State,ExitCode,ElapsedRaw,AllocCPUS,MaxRSS,TotalCPU,Start,End,QOS,Timelimit'
    p=subprocess.run([module.SLURM+'sacct','-j',','.join(sorted(selected)),'-P','--format='+fields],
        stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True,timeout=60)
    if p.returncode:raise RuntimeError(p.stderr)
    raw=[{k:v for k,v in row.items()if k}for row in csv.DictReader(io.StringIO(p.stdout),delimiter='|')]
    parents={row['JobID']:row for row in raw if row['JobID']in selected}
    rows=[];totals=collections.Counter()
    for job,metadata in sorted(selected.items()):
        if job not in parents:continue
        row=dict(parents[job],**metadata)
        row['allocated_core_hours']=int(row['ElapsedRaw'])*int(row['AllocCPUS'])/3600
        totals[row['kind']]+=row['allocated_core_hours'];rows.append(row)
    table(study/'tables/scheduler_accounting.tsv',rows)
    module.write(study/'manifests/scheduler_accounting_raw.json',dict(records=raw,scope='Only jobs in this review state; parent allocations counted once'))
    result=dict(generated_utc=module.now(),selected_parent_jobs=len(selected),accounted_parent_jobs=len(parents),
        missing_job_ids=sorted(set(selected)-set(parents)),
        unfinished_job_ids=[r['JobID']for r in rows if r['State'].split()[0]in('PENDING','RUNNING','COMPLETING','CONFIGURING','SUSPENDED')],
        allocated_core_hours_by_role=dict(totals),total_allocated_core_hours=sum(totals.values()),
        QOS_counts=dict(collections.Counter(r['QOS']for r in rows)),
        scheduler_state_counts=dict(collections.Counter(r['State']for r in rows)),
        scope='This follow-up only; includes tests, screening, canceled metadata and failed attempts; live totals are interim')
    module.write(study/'reports/resource_summary.json',result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project',type=Path);parser.add_argument('study',type=Path)
    parser.add_argument('--snapshot');parser.add_argument('--accounting',action='store_true')
    args=parser.parse_args();project=args.project.resolve();study=args.study.resolve()
    active=json.loads((study/'manifests/orchestrator_versions.json').read_text())['active']
    spec=importlib.util.spec_from_file_location('review',study/active['path']);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    r=module.Review(project,study);r.verify();live=r.state
    resources=accounting(module,study,live)if args.accounting else module.read(study/'reports/resource_summary.json')if(study/'reports/resource_summary.json').exists()else None
    if args.snapshot:r.state=module.read(study/'report_inputs'/args.snapshot/'review_state_snapshot.json')
    endpoints=[];points=[];parameters=[];boundary=[]
    domains=sorted({p['domain']for p in r.state['points'].values()})
    baseline_rows=[]
    for domain in domains:
        bid='base_'+domain;values=r.values(bid)
        if not values:continue
        best=min(values,key=lambda x:x['objective']);gate=module.recovered(values,True)
        baseline_rows.append(dict(domain=domain,**gate))
        baseline_metadata=module.read((project/best['source_fit']).parent/'validated_result.json')
        cases=[(pid,p)for pid,p in r.state['points'].items()if p['domain']==domain and p['fixed']]
        if domain=='original':
            for pid,p in r.old['points'].items():
                if pid not in r.state['points']:
                    cases.append((pid,dict(domain='original',fixed={p['target']:p['value']},priority='not_released',task_ids=[])))
        for pid,point in cases:
            released=pid in r.state['points'];values=r.values(pid)if released else r.old_values(pid)
            quality=module.recovered(values);target,value=next(iter(point['fixed'].items()))
            if not values:
                points.append(dict(domain=domain,point_id=pid,target=target,value=value,released_in_review=released,reason='no_verified_candidate'));continue
            winner=min(values,key=lambda x:x['objective'])
            source=(project/winner['source_fit']).parent;metadata=module.read(source/'validated_result.json')
            original=min(r.old_values(pid),key=lambda x:x['objective'])if domain=='original'and pid in r.old['points']else None
            row=dict(domain=domain,point_id=pid,target=target,value=value,released_in_review=released,
                priority=point['priority'],accepted_candidates=len(values),baseline_task_id=best['task_id'],
                baseline_objective=best['objective'],delta=winner['objective']-best['objective'],
                old_point_objective=original['objective']if original else None,
                old_to_new_objective_improvement=original['objective']-winner['objective']if original else None,**quality)
            for modality in['growth','karyotype','flow']:
                row[modality+'_component_delta']=baseline_metadata['component_means'][modality]-metadata['component_means'][modality]
            if abs(sum(row[m+'_component_delta']for m in['growth','karyotype','flow'])-row['delta'])>1e-8:
                raise RuntimeError('Component traces do not reconstruct delta')
            points.append(row)
            for result in values:
                owner=(project/result['source_fit']).parent;detail=module.read(owner/'validated_result.json')
                endpoints.append(dict(domain=domain,point_id=pid,target=target,value=value,task_id=result['task_id'],
                    objective=result['objective'],delta=result['objective']-best['objective'],
                    new_pure_independent=result.get('pure_independent',False),mode=result.get('mode','parent'),
                    tested_directional_stability=result.get('tested_directional_stability',False),
                    budget_limited=detail.get('budget_limited'),actual_generations=detail.get('actual_generations'),
                    local_improvement=detail.get('local_improvement'),source_fit=result['source_fit']))
            with(source/'parameter_result.tsv').open()as stream:
                for p in csv.DictReader(stream,delimiter='\t'):
                    inactive=(p['parameter']=='k_o_mis'and winner['params']['p_misseg']==0)or(
                        p['parameter']=='buffer_n_exp'and winner['params']['buffer_beta']==0)or(
                        p['parameter']in('buffer_beta','buffer_n_exp')and winner['params']['buffer_smax']==0)
                    parameters.append(dict(p,domain=domain,point_id=pid,target=target,target_value=value,
                        winning_task_id=winner['task_id'],data_inactive_at_endpoint=inactive,
                        estimate_interpretable=p['status']=='active'and not inactive,
                        bound_class='physical_probability'if p['parameter']=='buffer_smax'else'declared_search_limit'))
    for tid,task in r.state['tasks'].items():
        f=study/'fits'/tid/'boundary_probes.tsv'
        if task['status']=='accepted'and f.exists():
            with f.open()as stream:
                boundary.extend(dict(row,task_id=tid,point_id=task['point_id'],domain=task['domain'])for row in csv.DictReader(stream,delimiter='\t'))
    table(study/'tables/review_point_quality.tsv',points);table(study/'tables/review_endpoint_quality.tsv',endpoints)
    table(study/'tables/review_parameter_paths.tsv',parameters);table(study/'tables/review_boundary_probes.tsv',boundary)
    table(study/'tables/review_baselines.tsv',baseline_rows)
    ledger=[]
    for tid,t in r.state['tasks'].items():
        row={k:t.get(k)for k in['task_id','point_id','domain','priority','mode','start_class','pure_independent','status','reserve','source_fit','job_id','execution_mode','reason']}
        file=study/'fits'/tid/'validated_result.json'
        if t['status']=='accepted'and file.exists():
            detail=module.read(file);row.update(objective=detail['objective'],actual_generations=detail['actual_generations'],
                budget_limited=detail['budget_limited'],tested_directional_stability=detail['tested_directional_stability'])
        ledger.append(row)
    table(study/'tables/review_execution_ledger.tsv',ledger)
    summary=dict(generated_utc=module.now(),execution_status=r.state['status'],snapshot=args.snapshot,
        baseline_gates=baseline_rows,point_checks_passed=sum(p.get('pass_gate',False)for p in points),
        registered_optimization_tasks=len(r.state['tasks']),global_attempts=r.state['global_attempts'],local_attempts=r.state['local_attempts'],
        task_status_counts=dict(collections.Counter(t['status']for t in r.state['tasks'].values())),
        reserve_global_attempts=r.state['reserve_global_attempts'],reserve_local_attempts=r.state['reserve_local_attempts'],
        domains=domains,scientific_hold=r.state.get('scientific_hold'),domain_release=r.state.get('domain_release'),
        total_allocated_or_committed_core_hours=module.commitment(live['reservations']),
        new_allocation_cap=10000,inference_calibration='absent',output_language='English',tumor_fitting=False)
    module.write(study/'reports/review_summary.json',summary)
    nulls=[p for p in points if float(p['value'])==0]
    report='''# Culture profile optimization review

Generated: {utc}. Execution status: **{status}**. Snapshot: {snapshot}.

This report separates accepted feasible candidates, independent-population recovery, tested local directional stability and sensitivity to search bounds. No formal confidence interval, likelihood-ratio threshold or p value is used. All simulations retain the production-v2 culture protocol and complete growth/chromosome/flow observations. Tumor fitting is outside this stage.

## Unrestricted references

{baselines}

The original-domain reducer retains parent and new accepted candidates. Expanded-domain sentinels use their own references. A lower restricted candidate requires unrestricted repair; negative deltas remain visible. A failed reference gate prevents a scientific claim of completed optimization.

## Exact mechanism-zero candidates

{nulls}

Component deltas are traces along complete-data fits, not separately refitted measurement-stream profiles. Unreleased original points remain labeled parent-pilot context. Two fresh independent populations within 0.001 and a completed local check are required for point recovery. Native global iteration caps, local tested directions and historical optimization failures retain their separate meanings.

## Execution and allocation

Registered optimizer tasks: {tasks}; global attempts: {global_count}; local attempts: {local_count}. Status counts: {counts}. Reserve global/local attempts: {reserve_global}/{reserve_local}. Total actual or conservatively committed allocation at this state: {committed:.6f} core-hours; new hard ceiling 10000. Receipt-selected consumed allocation: {actual}. Every job uses xxlarge, whose verified single-job limit is 12 hours. There is no application concurrency throttle. Failed screens and replaced metadata jobs remain in the ledger and cost scope. Failed screens are not independent optimized endpoints.

The high-induction initialization stratum changes proposals only, uses first-feasible draws without score ranking, and retains full-domain random population members. Its source, validation and member origins are fingerprinted. It does not change model bounds or create a prior. Physical probability limits remain; arbitrary search limits require the separately declared domain comparisons.

## Artifacts and interpretation limits

- tables/review_point_quality.tsv and review_endpoint_quality.tsv: minima, raw score differences, component traces, independent recovery and local support.
- tables/review_parameter_paths.tsv and review_boundary_probes.tsv: active/fixed/inactive coordinates, physical/search limits and inward probes.
- tables/review_execution_ledger.tsv and scheduler_accounting.tsv: all attempts and receipt-selected parent allocations.
- tables/review_evaluation_log_summary.tsv: finite evaluation fractions by phase, separate worker-call timing and local selected-day variation.
- diagnostics/<snapshot>/<domain>/: complete predictions, parameter/function paths and dense production-kernel event accounting.
- report_inputs/<snapshot>/sources.json: read-only source vectors and hashes; derived render state is separate from raw review state.

Numerical grid loss is distinct from biological failed daughters. Full-history event totals can be dominated by assumed prehistory; headline flux comparisons must use the same recorded-count/duration subset for every candidate. Growth-count conditioning survives future growth-score deletion. No stream deletion, grid expansion, COMBOS or tumor study is launched by this report.
'''.format(utc=module.now(),status=r.state['status'],snapshot=args.snapshot or'live metadata',
        baselines=markdown(baseline_rows,['domain','best_task','objective','pass_gate','fresh_populations_near_minimum','stable_local_checks']),
        nulls=markdown(nulls,['domain','target','objective','delta','pass_gate','growth_component_delta','karyotype_component_delta','flow_component_delta']),
        tasks=len(r.state['tasks']),global_count=r.state['global_attempts'],local_count=r.state['local_attempts'],
        counts=summary['task_status_counts'],reserve_global=r.state['reserve_global_attempts'],reserve_local=r.state['reserve_local_attempts'],
        committed=module.commitment(live['reservations']),actual=resources['total_allocated_core_hours']if resources else'not yet refreshed')
    (study/'reports/culture_review_report.md').write_text(report)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
