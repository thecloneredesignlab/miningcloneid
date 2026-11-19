#!/usr/bin/env python3
"""Standard-library validation only; no optimization or data-generation mode."""
from pathlib import Path
import csv
import hashlib
import sys
import datetime
sys.dont_write_bytecode = True


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def kv(path, values):
    write(path, [dict(key=k, type='boolean' if isinstance(v, bool) else 'integer' if isinstance(v, int) else 'double' if isinstance(v, float) else 'string', value=str(v).lower() if isinstance(v, bool) else str(v)) for k, v in values.items()])


def values(path):
    return {r['key']: r['value'] for r in read(path)}


def verify(project, study, label):
    rows = []
    for manifest, root in [('code_manifest.tsv', study / 'implementation'), ('input_manifest.tsv', study), ('existing_source_manifest.tsv', project)]:
        for r in read(study / 'manifests' / manifest):
            path = root / r['path']
            actual = sha(path)
            assert actual == r['sha256'], 'Source/input fingerprint mismatch: ' + str(path)
            rows.append(dict(group=manifest, path=str(path), sha256=actual, pass_gate=True))
    for r in read(study / 'cases.tsv'):
        assert sha(project / r['path']) == r['sha256'], 'Existing fit changed'
    write(study / 'results' / label / 'source_input_integrity.tsv', rows)
    print('Verified pinned source/input hashes:', len(rows), flush=True)


def summarize(project, study):
    verify(project, study, 'summary')
    meta = values(study / 'results/metadata/validation.tsv')
    assert meta['pass'] == 'true' and meta['container_runtime_verified'] == 'true'
    assert (study / 'results/metadata/COMPLETE').is_file()
    assert values(study / 'reports/independent_comparison.tsv')['pass'] == 'true'
    receipts = []
    propagation_count = 0
    for case in read(study / 'cases.tsv'):
        for rep in [1, 2]:
            folder = study / 'results/replays' / (case['case'] + '_rep' + str(rep))
            assert (folder / 'COMPLETE').is_file(), 'Incomplete fixed-fit slot: ' + str(folder)
            v = values(folder / 'validation.tsv')
            assert v['pass'] == 'true' and v['source_sha256'] == case['sha256']
            assert int(v['model_propagations']) == 1 and int(v['optimizer_calls']) == 0
            assert int(v['synthetic_datasets']) == 0 and abs(float(v['score_difference'])) <= 1e-8
            propagation_count += int(v['model_propagations'])
            receipts.append(dict(case=case['case'], replicate=rep, objective=v['objective'], archived_objective=v['archived_objective'], independent_objective=v['independent_objective'], score_difference=v['score_difference'], source_sha256=v['source_sha256'], pass_gate=True))
    assert propagation_count == 8
    slots = list((study / 'results/full_forward_slots').iterdir())
    assert len(slots) == 8
    wiring = read(study / 'results/replays/production_seed144_rep1/parameter_wiring_audit.tsv')
    assert len(wiring) == 20 and sum(r['scientific_role'] == 'mechanism' for r in wiring) == 14
    assert all(r['pass'] == 'TRUE' for r in wiring)
    write(study / 'reports/fixed_fit_replays.tsv', receipts)
    write(study / 'reports/parameter_wiring_audit.tsv', wiring)
    layers = [dict(layer='environment', status='PASS', interpretation='Pinned SIF, actual namespaces and compiled source verified'), dict(layer='mapping_and_numerical_interfaces', status='PASS', interpretation='Twenty coordinates and eight completed fixed-fit fresh-process replays'), dict(layer='physical_experiment_contract', status='HELD', interpretation='Original SOP, count dependence, exact flow processing and live/dead operator/time gaps remain'), dict(layer='synthetic_recovery_or_identifiability', status='NOT_PERFORMED', interpretation='No synthetic observations, fitting, confidence intervals or identifiability claim')]
    write(study / 'reports/validation_layers.tsv', layers)
    gaps = read(study / 'inputs/local_contract/unresolved_inputs.tsv')
    write(study / 'reports/unresolved_inputs.tsv', gaps)
    utc = datetime.datetime.now(datetime.timezone.utc).isoformat()
    kv(study / 'reports/validation.tsv', dict(pass_gate=True, utc=utc, environment_pass=True, mapping_and_replay_pass=True, physical_experiment_release=False, full_forward_evaluations=8, optimizer_calls=0, synthetic_datasets=0, forward_slots_consumed=8))
    lines = ['# Production-SIF validation of culture Steps 01-03', '', 'Environment and fixed-parameter interface validation passed. The faithful physical-experiment contract remains held. No optimization or synthetic observation generation was performed.', '', '## Runtime and input identity', '', 'Image SHA256: dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372. Actual loaded versions and paths, R/C++ compiler evidence, isolated cold compilation and source/input hashes are retained in results/. Production and existing controllers were not changed.', '', '## Parameters and existing-fit replay', '', 'Fourteen mechanisms plus four initial-distribution and two error nuisance coordinates are retained. Seven prior mapping cases, sixty deterministic coordinate/boundary roundtrips and two invalid error-scale checks passed in the SIF. All twenty coordinates reached their declared normalization/operator paths; numerical plateaus are not interpreted as inactivity or nonidentifiability.', '', 'Four existing fits were replayed in two fresh processes each, using real observations and unchanged parameter vectors. Existing frozen score/state checks remain in force. No new truth, starts, optimizer, profiles or synthetic dataset was created.', '', '| Case | Repeat | Replayed Q | Difference from archive |', '|---|---|---|---|']
    lines += ['| '+r['case']+' | '+str(r['replicate'])+' | '+r['objective']+' | '+r['score_difference']+' |' for r in receipts]
    lines += ['', '## Observation contract', '', 'All 114 growth summaries reproduce the verified multi-time-point OLS definition. Saved selection and child-boundary routines reproduce the retained real-data trajectories. Karyotype sample sizes, original flow histogram counts and processed-density units remain distinct. Live/dead scoring remains disabled; the existing 6.1293-day assay-versus-selected-time discrepancy is reproduced, not repaired.', '', 'The missing laboratory stopping/inoculum SOP, raw-count error/dependence, exact 500-bin/G0/G1 flow processing release, retained event sizes, classification/retention and physical assay-time/reassignment records remain open. The smax-zero special activity case remains a separately registered mapper limitation. Container agreement does not resolve these scientific inputs.', '', 'Tie and infeasibility branches without a matching retained input are source-audited rather than claimed as new empirical coverage. Neither the environment checks nor fixed-fit replays establish parameter recovery, biological correctness or statistical confidence.', '', '## Scope and resource closure', '', 'All allocations use red/xxlarge, no requeue and no application concurrency throttle. The independent ceiling is six allocated core-hours, including at most one infrastructure retry. Final terminal scheduler states, elapsed times and resource totals are recorded separately in accounting/final_parent_accounting.tsv and COMPLETE after all owned allocations close. Prior culture and in-vivo budgets are not reused.', '', 'All new files use English R/Python/shell/TSV/RDS/plain text; no Java, JavaScript or JSON execution/state/reporting dependency. Loaded-namespace checks apply to every R process. Source image/code/input hashes and model propagation counts are checked separately from scientific and collection status.']
    (study / 'REPORT.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    (study / 'reports/COMPLETE').write_text('Container numerical validation complete; physical design release remains held.\n')
    print('Summary complete: environment and interfaces PASS; physical design HELD; eight forward calls, zero optimizations.', flush=True)


if __name__ == '__main__':
    assert len(sys.argv) == 5
    mode, project, study, label = sys.argv[1:]
    project, study = Path(project), Path(study)
    if mode == 'verify':
        verify(project, study, label)
    elif mode == 'summary':
        summarize(project, study)
    else:
        raise ValueError('No fitting or generation mode exists')
