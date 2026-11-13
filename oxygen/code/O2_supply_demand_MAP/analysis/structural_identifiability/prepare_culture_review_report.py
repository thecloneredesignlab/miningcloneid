#!/usr/bin/env python3
"""Build disposable renderer inputs without changing accepted study artifacts."""
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import hashlib
from datetime import datetime, timezone


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project',type=Path);parser.add_argument('study',type=Path)
    parser.add_argument('snapshot',help='New, nonexisting report snapshot name')
    args=parser.parse_args();project=args.project.resolve();study=args.study.resolve()
    registry=json.loads((study/'manifests/orchestrator_versions.json').read_text())['active']
    spec=importlib.util.spec_from_file_location('review',study/registry['path'])
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    r=module.Review(project,study);r.verify()
    raw=(study/'state.json').read_bytes();r.state=json.loads(raw)
    root=study/'report_inputs'/args.snapshot
    if root.exists():raise RuntimeError('A report snapshot must not be overwritten')
    root.mkdir(parents=True)
    (root/'review_state_snapshot.json').write_bytes(raw)
    raw_hash=hashlib.sha256(raw).hexdigest()
    sources=[];domain_manifest=[];point_rows=[]
    for domain in sorted({p['domain']for p in r.state['points'].values()}):
        folder=root/domain
        for d in ['fits','scripts','manifests']: (folder/d).mkdir(parents=True,exist_ok=True)
        (folder/'implementation').symlink_to(study/'implementation',target_is_directory=True)
        (folder/'scripts/culture_profile_diagnostics.R').symlink_to(Path(__file__).resolve().parent/'culture_profile_diagnostics.R')
        settings=dict(wp0_relative_path=r.config['wp0_relative_path'],optimization_agreement_tolerance=.001,
                      domain=domain,review_state_sha256=raw_hash,output_language='English')
        module.write(folder/'study.json',settings)
        tasks={};points={};qualities={}
        bid='base_'+domain
        values=r.values(bid)
        if not values:continue
        baseline_gate=module.recovered(values,True)

        def expose(value,pid,stage):
            tid=value['task_id'];source=project/value['source_fit'];owner=source.parent
            if tid.startswith('parent_'):
                original=module.read(owner/'validated_result.json')
                if module.sha(source)!=original['files_sha256']['accepted_fit.rds']:
                    raise RuntimeError('Parent accepted vector changed')
                trace=original['optimizer_trace'];phase=original['winner_phase'];budget=original['budget_limited']
                params=value['params'];origin='parent_pilot'
            else:
                original=module.read(owner/'validated_result.json')
                guard=module.read(owner/'artifact_validation.json')
                if not guard['pass']or module.sha(owner/'artifact_validation.json')!=original['artifact_validation_sha256']:
                    raise RuntimeError('Review acceptance guard changed')
                for name,digest in original['files_sha256'].items():
                    if module.sha(owner/name)!=digest:raise RuntimeError('Accepted review payload changed')
                trace=original.get('native_optimizer_trace')or dict(iter_completed=0,local_convergence=None)
                phase=original['mode'];budget=original['budget_limited'];params=value['params'];origin=original['start_class']
            target=folder/'fits'/tid;target.mkdir(exist_ok=True)
            for name in ['accepted_fit.rds','parameter_result.tsv','parameter_map.rds']:
                (target/name).symlink_to(owner/name)
            compatible=dict(accepted=True,launch_completed=True,objective=value['objective'],
                component_means=original['component_means'],optimizer_trace=trace,winner_phase=phase,
                budget_limited=budget,params=params,tested_directional_stability=value.get('tested_directional_stability',False),
                files_sha256={'accepted_fit.rds':module.sha(source)},source_task=tid,source_class=origin,
                derived_renderer_metadata=True)
            module.write(target/'validated_result.json',compatible)
            tasks[tid]=dict(status='accepted',stage=stage,role=phase,primary_family=origin)
            sources.append(dict(domain=domain,point_id=pid,task_id=tid,source_fit=str(source.relative_to(project)),
                accepted_fit_sha256=module.sha(source),source_metadata_sha256=module.sha(owner/'validated_result.json')))
            return tid

        for value in values:expose(value,bid,'baseline')
        cases=[(pid,p)for pid,p in r.state['points'].items()if p['domain']==domain and p['fixed']]
        if domain=='original':
            # Unreleased original points remain visible as parent-pilot context.
            for pid,old in r.old['points'].items():
                if pid not in r.state['points']:
                    cases.append((pid,dict(domain='original',parent_point=pid,
                        fixed={old['target']:old['value']},priority='not_released',task_ids=[])))
        for pid,point in cases:
            released=pid in r.state['points']
            values=r.values(pid)if released else r.old_values(pid)
            gate=module.recovered(values)
            gate['released_in_review']=released
            qualities[pid]=gate
            target,value=next(iter(point['fixed'].items()))
            ids=[expose(v,pid,'profile')for v in values]
            points[pid]=dict(point_id=pid,target=target,value=value,task_ids=ids)
            point_rows.append(dict(domain=domain,point_id=pid,target=target,value=value,
                attempts=len(point['task_ids']),**gate))
        if domain=='original':
            with (folder/'old_profile_minima.tsv').open('w')as stream:
                writer=csv.DictWriter(stream,fieldnames=['point_id','target','value','objective'],delimiter='\t');writer.writeheader()
                for pid,p in r.old['points'].items():
                    best=min(r.old_values(pid),key=lambda x:x['objective'])
                    writer.writerow(dict(point_id=pid,target=p['target'],value=p['value'],objective=best['objective']))
        module.write(folder/'state.json',dict(tasks=tasks,points=points,status=r.state['status'],baseline_gate=baseline_gate))
        module.write(folder/'review_quality.json',qualities)
        domain_manifest.append(dict(domain=domain,input_directory=str(folder),
            artifact_directory=str(study/'diagnostics'/args.snapshot/domain),baseline_gate=baseline_gate))
    module.write(root/'sources.json',dict(generated_utc=datetime.now(timezone.utc).isoformat(),
        review_state_sha256=raw_hash,parent_state_sha256=r.config['parent_state_sha256'],
        description='Derived renderer metadata with read-only links to accepted vectors; this is not a fitting study',sources=sources))
    module.write(root/'domains.json',domain_manifest)
    module.write(root/'quality.json',point_rows)
    print(json.dumps(dict(snapshot=str(root),domains=len(domain_manifest),source_vectors=len(sources),state_sha256=raw_hash),indent=2))


if __name__=='__main__':main()
