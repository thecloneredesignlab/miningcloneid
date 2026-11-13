#!/usr/bin/env python3
"""Join rendered function traces and recording-scope fluxes for each domain."""
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import math


def rows(path):
    with path.open()as stream:return list(csv.DictReader(stream,delimiter='\t'))


def write_table(path, values):
    fields=list(dict.fromkeys(k for row in values for k in row))
    with path.open('w')as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(values)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('project',type=Path);parser.add_argument('study',type=Path);parser.add_argument('snapshot')
    args=parser.parse_args();project=args.project.resolve();study=args.study.resolve()
    helper=Path(__file__).resolve().parent/'summarize_culture_profile_findings.py'
    spec=importlib.util.spec_from_file_location('scope',helper);scope_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(scope_module)
    inputs=study/'report_inputs'/args.snapshot
    domains=json.loads((inputs/'domains.json').read_text());settings=json.loads((study/'study.json').read_text())
    manifest=project/settings['wp0_relative_path']/'manifests/observation_manifest.tsv'
    scope=scope_module.recording_scope(rows(manifest));summary=[]
    all_traces=[];all_totals=[]
    for domain in domains:
        folder=Path(domain['artifact_directory']);name=domain['domain']
        status=json.loads((folder/'reports/diagnostic_artifacts_status.json').read_text())
        baseline_id=status['baseline_task_id']
        baseline=json.loads((Path(domain['input_directory'])/'fits'/baseline_id/'validated_result.json').read_text())
        parameters={baseline_id:baseline['params']}
        nuisance=rows(folder/'tables/nuisance_paths.tsv');paths={};inactive={}
        for row in nuisance:
            paths.setdefault(row['point_id'],{})[row['parameter']]=float(row['estimate'])
            parameters.setdefault(row['winning_task_id'],{})[row['parameter']]=float(row['estimate'])
            if row['data_inactive_at_endpoint']=='TRUE':inactive.setdefault(row['point_id'],[]).append(row['parameter'])
        functions={(r['point_id'],int(r['N']),float(r['oxygen_pct']),int(r['error_size'])):r for r in rows(folder/'tables/state_function_values.tsv')}
        endpoint={(r['point_id'],r['task_id']):r for r in rows(folder/'tables/profile_endpoints.tsv')}
        traces=[]
        for point in rows(folder/'tables/profile_points_validated.tsv'):
            if point['objective']=='NA':continue
            pid=point['point_id'];trace=dict(point,domain=name,**paths[pid])
            trace['data_inactive_parameters']=';'.join(sorted(inactive.get(pid,[])))
            selected=endpoint[pid,point['winner']]
            for modality in['growth','karyotype','flow']:
                trace[modality+'_component_delta']=float(selected[modality+'_contribution'])+baseline['component_means'][modality]
            if abs(sum(trace[m+'_component_delta']for m in['growth','karyotype','flow'])-float(point['delta']))>1e-8:
                raise RuntimeError('Rendered component traces do not reconstruct the point delta')
            for N in[44,88]:
                for oxygen,label in[(0.,'O0'),(20.5,'O20p5')]:
                    f=functions[pid,N,oxygen,1]
                    p=float(f['p_per_chromosome']);s=float(f['survival']);z=1-p+p*s
                    trace['p_N%d_%s'%(N,label)]=p
                    trace['q_conditional_N%d_%s'%(N,label)]=p*s/z if z>0 else None
                    trace['Y_biological_N%d_%s'%(N,label)]=2*z**N
                    trace['division_N%d_%s'%(N,label)]=float(f['division_rate_per_day'])
                    trace['death_N%d_%s'%(N,label)]=float(f['death_rate_per_day'])
                trace['survival_N%d_m1'%N]=float(functions[pid,N,0.,1]['survival'])
            trace['mu_hp_over_k_o_mis']=trace['mu_hp']/trace['k_o_mis']if trace['p_misseg']>0 else None
            traces.append(trace);all_traces.append(trace)
        write_table(folder/'tables/compensation_and_component_traces.tsv',traces)
        labeled,totals=scope_module.flux_summaries(rows(folder/'tables/passage_flux_summaries.tsv'),scope,parameters)
        write_table(folder/'tables/segment_flux_recording_scope.tsv',list(scope.values()))
        write_table(folder/'tables/passage_flux_scope_audit.tsv',labeled)
        for row in totals:
            row['domain']=name
            row['ordinary_grid_loss_per_ordinary_division']=(row['ordinary_boundary_loss']+row['ordinary_boundary_gain'])/row['ordinary_divisions']
            row['viable_tail_pruning_per_ordinary_division']=row['numerical_tail_pruning']/row['ordinary_divisions']
        write_table(folder/'tables/flux_summary_by_recording_scope.tsv',totals);all_totals.extend(totals)
        baseline_flux=[row for row in totals if row['task_id']==baseline_id]
        unrecorded=sum(row['maternal_divisions']for row in baseline_flux if row['recording_scope']=='counts_or_duration_unrecorded')
        full=sum(row['maternal_divisions']for row in baseline_flux)
        summary.append(dict(domain=name,baseline_task_id=baseline_id,baseline_objective=status['baseline_objective'],
            recorded_segments=sum(v['recording_scope']=='counts_and_duration_recorded'for v in scope.values()),
            unrecorded_segments=sum(v['recording_scope']=='counts_or_duration_unrecorded'for v in scope.values()),
            baseline_unrecorded_division_share=unrecorded/full,flux_segment_ledgers=status['flux_segment_count']))
    write_table(study/'tables/diagnostic_compensation_traces.tsv',all_traces)
    write_table(study/'tables/diagnostic_flux_by_recording_scope.tsv',all_totals)
    (study/'reports/diagnostic_aggregation.json').write_text(json.dumps(dict(snapshot=args.snapshot,domains=summary,
        interpretation='Function traces, not derived-target profiles; numerical loss is separate from biological nonviability; recorded-scope totals remain conditional on production history'),indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
