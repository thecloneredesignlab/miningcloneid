#!/usr/bin/env python3
"""Verify captured metadata and held budgets. No scientific computations."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
from pathlib import Path


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def kv(path):
    return {r['key']: r['value'] for r in rows(path)}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, fields, records):
    with path.open('w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', required=True, type=Path)
    args = parser.parse_args()
    r = args.run.resolve()
    if 'invivo_structural_' not in r.name:
        raise SystemExit('Not a new independent tumor audit namespace')
    t = r / 'tables'
    checks = []

    def check(name, condition, evidence):
        checks.append(dict(check=name, passed=str(bool(condition)).lower(), evidence=str(evidence)))

    inv = rows(r / 'inputs/source_inventory.tsv')
    for f in inv:
        if f['frozen_copy'] == 'not_copied_not_executed':
            continue
        p = r / f['frozen_copy']
        check('frozen_input:' + f['frozen_copy'], sha(p) == f['sha256'] and p.stat().st_size == int(f['bytes']), f['sha256'])
    local_by_role = {Path(f['source']).name: f for f in inv if f['role'] in ('old_numerical_source', 'old_contract')}
    hpc_hashes = rows(r / 'snapshots/hpc/source_hashes.tsv')
    sif = [f for f in hpc_hashes if f['path'].endswith('.sif')]
    check('specified_SIF_fresh_hash', len(sif) == 1 and sif[0]['sha256'] == 'dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372', sif)
    for f in hpc_hashes:
        if f['path'].endswith('.sif'):
            continue
        local = local_by_role[Path(f['path']).name]
        check('HPC_local_source:' + Path(f['path']).name, local['sha256'] == f['sha256'], f['sha256'])
    compare = []
    for f in inv:
        if f['role'] != 'current_checkout_comparison':
            continue
        old = local_by_role[Path(f['source']).name]
        compare.append(dict(current_path=f['source'], current_sha256=f['sha256'], old_sha256=old['sha256'], identical=str(f['sha256'] == old['sha256']).lower()))
    write(t / 'current_to_frozen_source_comparison.tsv', ['current_path', 'current_sha256', 'old_sha256', 'identical'], compare)
    check('current_checkout_comparison_complete', len(compare) == 5, 'comparison does not upgrade production authority')

    def scheduler(scope, ledger, id_key, cost_key):
        base = r / 'snapshots/hpc' / scope
        fields = ['job_id', 'state', 'exit_code', 'cpus', 'elapsed_seconds', 'start', 'end']
        with (base / 'sacct.tsv').open() as f:
            sacct = [dict(zip(fields, line.rstrip('\n').split('|')[:7])) for line in f if line.strip()]
        ids = {j[id_key] for j in ledger}
        sids = {j['job_id'] for j in sacct}
        check(scope + ':scoped_parent_set', ids == sids and len(sacct) == len(ids), len(ids))
        terminal = ('COMPLETED', 'FAILED', 'CANCELLED', 'TIMEOUT', 'OUT_OF_MEMORY', 'NODE_FAIL', 'PREEMPTED', 'BOOT_FAIL', 'DEADLINE')
        check(scope + ':all_terminal', all(j['state'].split()[0].rstrip('+') in terminal for j in sacct), Counter(j['state'] for j in sacct))
        check(scope + ':scoped_queue_empty', not (base / 'squeue.tsv').read_text().strip(), 'exact owned-parent IDs only')
        actual = sum(int(j['cpus']) * int(j['elapsed_seconds']) for j in sacct)
        check(scope + ':allocated_bill_matches_ledger', actual == sum(int(j[cost_key]) for j in ledger), actual)
        return sacct, actual

    old = r / 'snapshots/hpc/old_invivo'
    jobs = rows(old / 'jobs.tsv')
    sacct, old_actual = scheduler('old_invivo', jobs, 'job_id', 'actual_core_seconds')
    check('old_parent_count_125', len(jobs) == 125, len(jobs))
    check('old_states_113_completed_10_failed_2_cancelled', Counter(j['state'].split()[0] for j in sacct) == Counter(COMPLETED=113, FAILED=10, CANCELLED=2), Counter(j['state'] for j in sacct))
    check('old_density_profile_not_launched', not any('profile' in j['role'].lower() or j['task'].startswith('p_') for j in jobs), '125 owned task/role rows')
    mode = kv(old / 'tables/resolution_closure_mode.tsv')
    check('old_closure_only', mode['mode'] == 'closure_only' and mode['new_synthetic_authorized_core_seconds'] == '0' and mode['culture_changes_authorized'] == 'false', mode['mode'])
    phase = {p['phase']: p for p in rows(old / 'tables/resolution_final_phase_accounting.tsv')}
    term = kv(old / 'tables/resolution_terminal.tsv')
    audit = kv(old / 'tables/closure_audit.tsv')
    b, p = phase['baseline'], phase['profile']
    check('old_baseline_bill_and_cap', int(b['actual_core_seconds']) == 318729 and int(b['cap_core_seconds']) == 792000, b)
    check('old_profile_zero_and_unreleased_360', int(p['actual_core_seconds']) == 0 and int(p['cap_core_seconds']) == 1296000, p)
    check('old_current_bill', old_actual == 405407 == int(audit['new_actual_core_seconds']), old_actual)
    lifetime = int(term['revised_historical_core_seconds']) + int(b['actual_core_seconds']) + int(p['actual_core_seconds'])
    check('old_lifetime_accounting_identity', lifetime == 707253 == 301846 + old_actual == int(audit['lifetime_actual_core_seconds']), lifetime)
    check('old_delivery_proof_passed', kv(old / 'tables/delivery_verified.tsv')['pass'] == 'true' and audit['pass'] == 'true', 'existing full-payload proof read; 83196 files not rehashed by this audit')
    gates = rows(old / 'tables/baseline_policy_gates.tsv')
    check('old_both_recovery_gates_remain_false', len(gates) == 2 and all(g['pass'] == 'FALSE' and g['family_recovery'] == 'FALSE' for g in gates), [(g['policy'], g['uniform_gap']) for g in gates])

    culture = r / 'snapshots/hpc/culture_step04'
    cj = rows(culture / 'accounting/jobs.tsv')
    _, cactual = scheduler('culture_step04', cj, 'job_id', 'actual')
    check('culture_step04_six_parents_940_seconds', len(cj) == 6 and cactual == 940, cactual)
    ca = kv(culture / 'reports/step04_acceptance.tsv')
    check('upstream_hold_preserved', ca['reduction_adequacy_gate'] == 'false' and ca['main_structural_case_gate'] == 'false' and ca['step05_release'] == 'false', ca)
    check('upstream_local_HPC_gate_identical', sha(culture / 'reports/step04_acceptance.tsv') == sha(r / 'inputs/culture_gate/reports/step04_acceptance.tsv'), 'fresh source compared to captured local evidence')

    param = rows(t / 'production_parameter_table.tsv')
    roles = rows(t / 'parameter_roles.tsv')
    active = {p['param_name'] for p in param if p['estimate'] == 'TRUE'}
    check('all_20_active_coordinates_accounted', len(active) == 20 and {p['coordinate'] for p in roles if p['estimated'] == 'true'} == active, len(active))
    check('fixed_tau_not_fitted_coordinate', len(param) == len(roles) == 21 and [p['natural_parameter'] for p in roles if p['estimated'] == 'false'] == ['tau_O2'], '20 estimated + fixed tau')
    settings = {p['setting']: p['value'] for p in rows(t / 'production_scalar_settings.tsv')}
    expected = {'Crowding':'FALSE', 'fit_tau_O2':'FALSE', 'fit_treatment':'FALSE', 'harvest_init_multiplier':'FALSE', 'ploidy_O2_death':'ploidy_related', 'start_with':'chr_number'}
    check('production_switches_exact', all(settings[k] == v for k, v in expected.items()), expected)
    check('production_scale_conventions', float(settings['init_total_size']) == float(settings['o2_Nref']) == 1e6 and float(settings['tau_O2']) == .1, 'not experimental scale measurements')
    initial = rows(t / 'archived_initial_live_states.tsv')
    for cohort, state in [('2N',44),('4N',88)]:
        nz = [x for x in initial if float(x['initial_' + cohort]) != 0]
        check('pure_archived_IC_' + cohort, len(nz) == 1 and int(nz[0]['chromosome_state']) == state and float(nz[0]['initial_' + cohort]) == 1e6, nz)
    scenarios = rows(t / 'scenario_inventory.tsv')
    check('actual_observation_counts', len(scenarios) == 8 and sum(int(s['burden_raw']) for s in scenarios) == 120 and sum(int(s['burden_positive']) for s in scenarios) == 92 and sum(int(s['burden_zero']) for s in scenarios) == 28 and sum(int(s['chromosome_records']) for s in scenarios) == 6502 and sum(int(s['necrosis_record_count']) for s in scenarios) == 6, '8 scenarios;120 raw/92 positive;6502 chromosome estimates;6 terminal histology records')
    check('actual_untreated_and_histology_count', all(float(s['dose']) == 0 and int(s['necrosis_record_count']) in (0,1) for s in scenarios), 'no invented binomial denominator')

    cases = rows(t / 'symbolic_case_manifest.tsv')
    check('eight_finite_held_cases', len(cases) == len({c['case_id'] for c in cases}) == 8 and all(c['status'] == 'held_upstream' and c['structural_classification'] == 'not_run' and int(c['max_cores']) == 1 and int(c['max_wall_seconds']) == 7200 for c in cases), 'case equations remain formulation_pending')
    budget = rows(t / 'budget_ledger.tsv')
    primary = [b for b in budget if b['package'] == 'tumor_structural_32']
    check('32_hour_conditional_ceiling', sum(int(b['authorized_ceiling_core_seconds']) for b in primary) == 115200, '4 + 16 + 8 + 4 hours')
    check('no_compute_released_committed_or_spent', all(int(b[k]) == 0 for b in budget for k in ['released_core_seconds', 'committed_core_seconds', 'actual_core_seconds']), 'bridge synthetic old24 all unreleased')
    check('old24_and_bridge_no_new_authorization', all(int(b['authorized_ceiling_core_seconds']) == 0 for b in budget if b['package'] != 'tumor_structural_32'), '32 does not release proposed8 bridge or old24')
    releases = rows(t / 'release_gates.tsv')
    check('upstream_and_downstream_gate_holds', all(g['passed'] == 'false' for g in releases if g['gate_id'] >= 'G04'), 'no automatic release path exists')
    structural = rows(t / 'invivo_structural_results.tsv')
    check('no_fabricated_structural_conclusions', len(structural) == 21 and all(x['classification'] in ('not_run','fixed_convention') and x['proof_status'] == 'no_project_structural_analysis' for x in structural), 'includes probabilistic sigma coordinate as not_run')
    receipt = kv(r / 'snapshots/hpc/capture_receipt.tsv')
    check('read_only_capture_zero_model_calls', all(receipt[k] == '0' for k in ['remote_files_written','model_calls','controller_calls','slurm_submissions']) and receipt['old_pointer_stable'] == 'true', receipt['end_utc'])
    write(t / 'audit_verification.tsv', ['check', 'passed', 'evidence'], checks)
    summary = [dict(key=k, value=v) for k,v in [
        ('verified_utc',datetime.now(timezone.utc).isoformat()), ('metadata_checks',len(checks)),
        ('metadata_checks_passed',sum(c['passed']=='true' for c in checks)),
        ('old_owned_parents',len(jobs)), ('old_actual_core_seconds',old_actual),
        ('old_lifetime_core_seconds',lifetime), ('culture_step04_actual_core_seconds',cactual),
        ('new_HPC_actual_core_seconds',0), ('new_HPC_commitment_core_seconds',0),
        ('new_symbolic_calls',0), ('new_native_model_calls',0), ('new_fitting_calls',0),
        ('status','held_upstream_culture_and_formulation'),
        ('scientific_structural_gate_passed','false')]]
    write(t / 'audit_summary.tsv', ['key','value'], summary)
    failed = [c['check'] for c in checks if c['passed'] != 'true']
    if failed:
        raise SystemExit('Audit failed: ' + '; '.join(failed))
    print('%d/%d metadata checks passed; scientific gates remain held; new HPC cost 0.' % (len(checks),len(checks)))


if __name__ == '__main__':
    main()
