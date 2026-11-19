"""Summarize verified terminal TSVs; never run a model or dispatch a job."""
from pathlib import Path
import argparse, csv, hashlib, math

def read(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def kv(p):
    return {r['key']: r['value'] for r in read(p)}

def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()

def number(x):
    return '{:.8g}'.format(float(x))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('study', type=Path)
    s = ap.parse_args().study.resolve()
    proof = kv(s/'collection/local_verification.tsv')
    assert proof['pass_gate'].lower() == 'true'
    assert proof['delivery_manifest_sha256'] == digest(s/'tables/delivery_manifest.tsv')
    d = s/'collection'/proof['collection']
    captured = read(d/'captured_files.tsv')
    for row in captured:
        assert digest(s/row['path']) == row['sha256'], row['path']
    final = kv(s/'tables/final_accounting.tsv')
    jobs = read(s/'accounting/jobs.tsv')
    assert final['active_parents'] == '0'
    actual = sum(int(j['actual_core_seconds']) for j in jobs)
    assert actual == int(final['actual_core_seconds']) <= 128*3600
    gates = kv(s/'tables/engineering_gate.tsv')
    source = {r['mode']:r['task'] for r in read(s/'tables/engineering_sources.tsv')}
    prep = kv(s/'engineering'/source['prepare']/'result.tsv')
    truth_q = float(prep['truth_expected_objective'])
    assert gates['engineering_full_forward_slots'] == '64'
    assert len(read(s/'accounting/engineering_forward_slots.tsv')) == 64
    assert gates['pass'] == gates['cost_fits'] == 'true'
    tracked = ['tables/final_accounting.tsv','accounting/jobs.tsv',
        'tables/engineering_gate.tsv','tables/engineering_sources.tsv',
        'accounting/engineering_forward_slots.tsv',
        'engineering/'+source['prepare']+'/result.tsv']
    out = ['# Culture D0 recovery: verified terminal analysis', '',
        'Endpoint: **'+final['endpoint']+'**. This report derives only from the verified terminal TSV/RDS delivery; no additional forward evaluations are performed.', '',
        'The experiment is deterministic expected-observation recovery conditional on the preserved culture protocol/count histories and six known nuisance values. Fourteen mechanisms are estimated. It does not resolve the physical-observation contract holds or establish uniqueness, confidence intervals, or realistic noisy-data recovery.', '',
        '## Engineering and cost', '',
        'All 64 engineering forward slots were consumed. Truth expected score: '+number(truth_q)+'. The score includes entropy terms and is not expected to be zero. Native reference identity, synthetic response consumption, eight-worker identity and the production DEoptim API fixture passed.', '',
        'Measured feasible single-call p95: '+number(gates['measured_feasible_p95_seconds'])+' seconds; parallel warm throughput: '+number(gates['parallel_warm_calls_per_second'])+' calls/second. Registered conservative forecast per complete chain: '+number(gates['forecast_complete_chain_core_hours'])+' core-hours, below 58.', '',
        '| Task | Parent | Slurm / exit | Cores | Parent elapsed seconds | Charged core-hours |',
        '|---|---|---|---:|---:|---:|']
    for j in jobs:
        cost=int(j['actual_core_seconds']); cores=int(j['cores'])
        assert cost % cores == 0
        out.append('| {} | {} | {} / {} | {} | {} | {:.8f} |'.format(j['task'],j['job_id'],j['state'],j['exit_code'],cores,cost//cores,cost/3600))
    engineering = sum(int(j['actual_core_seconds']) for j in jobs if j['phase'] in ('engineering','buffer'))
    out += ['', 'Engineering including both historical startup failures: {:.8f} core-hours. Total actual: {} core-seconds = {:.8f}/128 core-hours. All parent costs are charged exactly once; nested numerical timers do not add fees.'.format(engineering/3600,actual,actual/3600), '',
        '## Blind fitting and fresh validation', '',
        'Each origin starts from 256 independent full-box Uniform vectors in the registered transformed domain. Truth and old fitted vectors were not injected. Global and local accepted candidates were independently replayed in a fresh SIF process; final estimates were selected by score before truth comparison.', '',
        '| Origin / stage | Callbacks | Feasible callbacks | Termination | Objective | Numerical wall seconds |',
        '|---|---:|---:|---|---:|---:|']
    if final.get('all_check_allocations_succeeded')=='false':
        evidence=kv(s/'tables/check_A_numerical_evidence.tsv')
        assert evidence['allocation_success']=='false' and evidence['end_environment_guard_completed']=='false'
        i=out.index('## Blind fitting and fresh validation')+2
        out[i:i]=['Check A remains FAILED / 1:0 after its two numerical replay assertions and parameter/function exports completed. X11 PNG device creation failed; the final environment guard and complete end-to-end check did not pass. Its numerical evidence is separately hash-bound; local table-derived figures do not change the allocation state. No additional A model calls or HPC allocation were used for presentation repair.','']
        tracked+=['tables/check_A_numerical_evidence.tsv','tables/check_A_numerical_payload.tsv']
    checks={}; params={}; funcs={}
    aliases={r['origin']:r['task'] for r in read(s/'tables/local_sources.tsv')} if (s/'tables/local_sources.tsv').exists() else {}
    if aliases:tracked.append('tables/local_sources.tsv')
    for b in ('A','B'):
        for stage in ('global_','local_'):
            task=aliases.get(b,stage+b) if stage=='local_' else stage+b
            base='fits/'+task; z=kv(s/base/'result.tsv'); logs=list((s/base).glob('evaluations_*.tsv'))
            records=[r for p in logs for r in read(p)]
            assert len(records)==int(z['evaluations'])
            feasible=sum(r['feasible'].lower()=='true' for r in records)
            out.append('| {} / {} | {} | {} | {} | {} | {} |'.format(b,task,len(records),feasible,z['termination'],number(z['objective']),number(z['wall_seconds'])))
            tracked += [base+'/result.tsv']+[str(p.relative_to(s)) for p in logs]
        base='checks/check_'+b
        z=checks[b]=kv(s/base/'result.tsv')
        assert z['pass']=='true' and z['data_sha256']==gates['data_sha256']
        assert float(z['objective']) >= truth_q-1e-8
        for r in read(s/base/'candidate_checks.tsv'):
            assert abs(float(r['difference']))<=1e-8
            assert abs(float(r['fresh'])-float(r['independent']))<=1e-8
        params[b]={r['parameter']:r for r in read(s/base/'parameter_recovery.tsv')}
        funcs[b]=read(s/base/'function_recovery.tsv')
        assert sum(r['estimated'].lower()=='true' for r in params[b].values())==14
        for r in params[b].values():
            if r['estimated'].lower()!='true':assert float(r['natural_error'])==0
        tracked += [base+'/'+name for name in ('result.tsv','candidate_checks.tsv','parameter_recovery.tsv','function_recovery.tsv')]
    if aliases.get('A')=='local_A_resume_01':
        old=[r for p in (s/'fits/local_A').glob('evaluations_*.tsv') for r in read(p)]
        new=[r for p in (s/'fits/local_A_resume_01').glob('evaluations_*.tsv') for r in read(p)]
        assert len(old)==2 and len(old)+len(new)<=11183
        local_jobs=[j for j in jobs if j['task'] in ('local_A','local_A_resume_01')]
        assert sum(int(j['actual_core_seconds']) for j in local_jobs)<=7200
        out += ['', 'Origin A retained its failed 22-core-second / 2-callback local startup. Its explicitly authorized continuation used only the remaining original stage allowance. Cumulative A local callbacks: {}/11,183; cumulative A local core-seconds: {}/7,200.'.format(len(old)+len(new),sum(int(j['actual_core_seconds']) for j in local_jobs))]
    out += ['', '| Origin | Selected stage | Verified score | Score minus truth | Max growth / true sigma | Max karyotype TV | Max flow TV | All three <=0.01 |',
        '|---|---|---:|---:|---:|---:|---:|---|']
    for b,z in checks.items():
        values=[float(z[k]) for k in ('max_growth_error','max_kary_tv','max_flow_tv')]
        assert all(v<=.01 for v in values)==(z['prediction_recovery']=='true')
        out.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(b,z['source_task'],number(z['objective']),number(float(z['objective'])-truth_q),*(number(v) for v in values),z['prediction_recovery']))
    out += ['', '## Parameter recovery', '',
        'Signed range-scaled error is (estimate - truth)/(registered natural upper - lower). It is a descriptive error, not an uncertainty interval. The six known nuisance coordinates are unchanged by construction and are not parameter-recovery successes.', '',
        '| Free mechanism | Truth | Origin A estimate | A signed range error | Origin B estimate | B signed range error |',
        '|---|---:|---:|---:|---:|---:|']
    for n,a in params['A'].items():
        if a['estimated'].lower()!='true':continue
        b=params['B'][n];assert a['truth']==b['truth']
        out.append('| {} | {} | {} | {} | {} | {} |'.format(n,number(a['truth']),number(a['estimate']),number(a['scaled_error']),number(b['estimate']),number(b['scaled_error'])))
    out += ['', '## Function recovery', '',
        'These are descriptive maximum absolute errors across assigned oxygen levels and chromosome counts 44/88, computed by the audited rate and offspring operators. Different functions have different units; no shared numerical acceptance threshold is imposed.', '',
        '| Function | Origin A max absolute error | Origin B max absolute error |',
        '|---|---:|---:|']
    for name in sorted({r['function_name'] for r in funcs['A']}):
        vals=[max(abs(float(r['error'])) for r in funcs[b] if r['function_name']==name) for b in ('A','B')]
        out.append('| {} | {} | {} |'.format(name,*(number(v) for v in vals)))
    out += ['', '## Interpretation and scope', '',
        'A missed prediction target leaves optimization unresolved at the authorized call/time windows; it does not demonstrate structural nonidentifiability. Matching predictions with displaced parameters identify candidate practical ambiguity under this conditional observation design. Neither result authorizes profiles or additional fits.', '',
        'Existing culture budgets/tasks and in-vivo closure authorities were independently hash-checked by the terminal controller. The proposed in-vivo 24 core-hours remain unreleased. Remaining D0 authorization is not transferred to any other study.', '',
        'Full scientific fit/improvement payloads remain on HPC. This local delivery is a hash-verified lean fit mirror and deliberately carries no fit COMPLETE markers. The captured delivery manifest and collection proof bind the results used here.']
    report=d/'RECOVERY_ANALYSIS.md';report.write_text('\n'.join(out)+'\n')
    tracked += [str(report.relative_to(s)), 'collection/local_verification.tsv']
    with (d/'analysis_manifest.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,['path','sha256'],delimiter='\t',lineterminator='\n');w.writeheader()
        w.writerows(dict(path=p,sha256=digest(s/p)) for p in sorted(set(tracked)))
    print(report)

if __name__=='__main__':
    main()
