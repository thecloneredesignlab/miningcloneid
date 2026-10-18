#!/usr/bin/env python3
"""Freeze and audit the approved 166-vector chromosome-domain experiment."""
import argparse
import collections
import datetime
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[5]
sys.path.insert(0, str(HERE.parent / 'high_precision_review'))
from review import sha, read, table, dump, audit_point, PARAMETERS

LIMITS = (154, 220, 308, 440)


def selected(rows):
    return [r for r in rows if r['fit_seed'] in ('seed138', 'seed101')
            and set(r['roles'].split(',')) & {'abnormal', 'ordinary_control'}]


def point_dir(root, task_id):
    return root / 'points' / ('task_%04d' % int(task_id))


def prepare(root):
    assert not (root / 'selection.json').exists(), 'Frozen selection already exists'
    prior = root.parent / 'high_precision_review'
    previous = json.loads((prior / 'selection.json').read_text())
    assert sha(prior / 'point_manifest.tsv') == previous['point_manifest_sha256']
    cases = selected(read(prior / 'point_manifest.tsv'))
    assert len(cases) == 166
    counts = collections.Counter((r['fit_seed'], r['roles']) for r in cases)
    assert counts == {('seed138', 'abnormal'): 11, ('seed138', 'ordinary_control'): 7,
                      ('seed101', 'abnormal'): 80, ('seed101', 'ordinary_control'): 68}
    root.mkdir(parents=True, exist_ok=True)
    tasks, proofs = [], {}
    evidence = {r['point_id']: r for r in read(prior / 'point_evidence.tsv')}
    for case_id, row in enumerate(cases, 1):
        path = prior / 'points' / ('point_%04d' % int(row['point_id']))
        report = json.loads((path / 'result.json').read_text())
        assert report['identity'] == row
        assert evidence[row['point_id']]['precision_verified'] == 'True'
        files = {name: report['files_sha256'][name] for name in
                 ('matrix.rds', 'context.rds', 'precision_50.json', 'precision_100.json')}
        for name, digest in files.items():
            assert sha(path / name) == digest
        files['result.json'] = sha(path / 'result.json')
        proofs[row['point_id']] = files
        for limit in LIMITS:
            tasks.append(dict(task_id=len(tasks) + 1, case_id=case_id, N_MAX=limit, **row))
    # Six fixed cases x all four limits: each seed has boundary, interior-high,
    # and ordinary responses. Select interior examples by least boundary mass.
    smoke_cases = []
    for seed in ('seed138', 'seed101'):
        high = [(i + 1, r) for i, r in enumerate(cases)
                if r['fit_seed'] == seed and r['roles'] == 'abnormal']
        score = lambda z: float(evidence[z[1]['point_id']]['upper_boundary_mass_100'])
        smoke_cases += [max(high, key=score)[0], min(high, key=score)[0]]
        smoke_cases += [next(i + 1 for i, r in enumerate(cases)
                             if r['fit_seed'] == seed and r['roles'] == 'ordinary_control')]
    assert len(set(smoke_cases)) == 6
    smoke = [int(t['task_id']) for t in tasks if t['case_id'] in smoke_cases]
    table(root / 'task_manifest.tsv', tasks)
    table(root / 'case_manifest.tsv', [dict(case_id=i + 1, **r) for i, r in enumerate(cases)])
    keys = {r['point_key'] for r in cases}
    relations = [r for r in read(prior / 'relationships.tsv')
                 if r['relationship'] == 'ordinary_control'
                 and r['abnormal_key'] in keys and r['related_key'] in keys]
    assert len(relations) == 91
    table(root / 'relationships.tsv', relations)
    dump(root / 'smoke_tasks.json', smoke)
    sources = dict(previous['source_plans'])
    inherited = {str(p.relative_to(REPO)): sha(p) for p in
                 (HERE.parent / 'high_precision_review' / n for n in
                  ('review.py', 'fixtures.R', 'perron_review.cpp'))}
    utility = HERE.parent / '../../util/o2_supply_demand_map_postfit_probe_utils.R'
    inherited[str(utility.resolve().relative_to(REPO))] = sha(utility)
    dump(root / 'selection.json', dict(
        created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        n_cases=166, n_tasks=664, n_expanded_tasks=498, n_baseline_rebuilds=166,
        N_MAX=list(LIMITS), N_MIN=22, N_UNIT=22, smoke_tasks=smoke,
        role_counts={'seed138': {'abnormal': 11, 'ordinary_control': 7},
                     'seed101': {'abnormal': 80, 'ordinary_control': 68}},
        prior_root=str(prior.relative_to(REPO)), prior_selection_sha256=sha(prior / 'selection.json'),
        prior_manifest_sha256=sha(prior / 'point_manifest.tsv'), baseline_proofs=proofs,
        task_manifest_sha256=sha(root / 'task_manifest.tsv'),
        case_manifest_sha256=sha(root / 'case_manifest.tsv'),
        relationships_sha256=sha(root / 'relationships.tsv'),
        source_plans=sources, fit_proofs=previous['fit_proofs'],
        code_sha256={p.name: sha(p) for p in sorted(HERE.iterdir())
                     if p.is_file() and p.suffix in ('.py', '.R', '.sbatch')},
        inherited_code_sha256=inherited, new_FAST_samples=0, unchanged_FAST_N_M=True,
        boundary='drop', eps_tail=1e-8,
        precision_scope='50/100-digit eigensolves of canonical double matrices; baseline reuses certified precision outputs after exact matrix/context match',
        diagnostic_thresholds=dict(ploidy_abs=0.01, growth_abs=1e-6, growth_rel=1e-3,
                                   distribution_TV=1e-3, upper_22_mass=1e-4),
        diagnostic_scope='308 to 440 only; thresholds are numerical diagnostics, not biological validation or exclusion rules',
        qos='xxlarge', time='12:00:00', specified_node=None, array_concurrency_limit=None))
    print('Frozen 166 cases, 664 tasks, 24 smoke tasks', flush=True)


def preflight(root, image=False):
    s = json.loads((root / 'selection.json').read_text())
    for name in ('task_manifest', 'case_manifest', 'relationships'):
        assert sha(root / (name + '.tsv')) == s[name + '_sha256'], name
    for name, digest in s['code_sha256'].items():
        assert sha(HERE / name) == digest, name
    for name, digest in s['inherited_code_sha256'].items():
        assert sha(REPO / name) == digest, name
    prior = REPO / s['prior_root']
    assert sha(prior / 'selection.json') == s['prior_selection_sha256']
    assert sha(prior / 'point_manifest.tsv') == s['prior_manifest_sha256']
    checked = {}
    for plan in s['source_plans'].values():
        for name, digest in plan['code_sha256'].items():
            p = HERE.parent / name
            if str(p) not in checked:
                checked[str(p)] = sha(p)
            assert checked[str(p)] == digest, name
        if image:
            for p, digest in [(Path(plan['sif']), plan['sif_sha256']),
                              (Path(plan['archives']['rcpp_template.tar.gz']['path']),
                               plan['archives']['rcpp_template.tar.gz']['sha256'])]:
                if str(p) not in checked:
                    checked[str(p)] = sha(p)
                assert checked[str(p)] == digest, str(p)
    return s


def audit(root, row):
    path = point_dir(root, row['task_id'])
    report = json.loads((path / 'result.json').read_text())
    assert report['identity'] == row
    assert report['runtime']['selection_sha256'] == sha(root / 'selection.json')
    assert report['runtime']['code_sha256'] == json.loads((root / 'selection.json').read_text())['code_sha256']
    z = audit_point(path)
    z.update(task_id=int(row['task_id']), case_id=int(row['case_id']),
             source_point_id=int(row['point_id']), N_MAX=int(row['N_MAX']),
             baseline_exact=report['domain']['baseline_matrix_exact'],
             common_block_max_difference=report['domain']['common_block_max_difference'],
             rates_exact_on_shared_states=report['domain']['rates_exact_on_shared_states'],
             precision_reused=report['domain']['precision_reused'])
    z['classification'] = 'precision_verified' if z['precision_verified'] else 'unresolved'
    return z


def collect(root, smoke=False):
    s = preflight(root)
    rows = read(root / 'task_manifest.tsv')
    if smoke:
        rows = [r for r in rows if int(r['task_id']) in s['smoke_tasks']]
    results = [audit(root, row) for row in rows]
    fields = list(dict.fromkeys(k for r in results for k in r))
    table(root / ('smoke_summary.tsv' if smoke else 'task_summary.tsv'),
          [{k: r.get(k, '') for k in fields} for r in results])
    record = dict(status='passed' if smoke else 'complete', n_tasks=len(results),
                  n_precision_verified=sum(r['precision_verified'] for r in results),
                  n_unique_vectors_certified=sum(r['unique_vector_certified'] for r in results),
                  selection_sha256=sha(root / 'selection.json'), new_FAST_samples=0,
                  n_baseline_exact=sum(r['N_MAX'] == 154 and r['baseline_exact'] for r in results))
    if smoke:
        assert all(r['precision_verified'] and r['unique_vector_certified'] for r in results), 'Unresolved smoke precision'
        for name in ('smoke_fixtures.json', 'cache_portability.json'):
            assert json.loads((root / name).read_text())['status'] == 'passed'
    dump(root / ('smoke_receipt.json' if smoke else 'completion.json'), record)
    print(json.dumps(record, indent=2), flush=True)


def compare(a, b, thresholds):
    """Padding preserves escaped mass; never renormalize to the common domain."""
    with localcontext() as ctx:
        ctx.prec = 140
        x, y = [Decimal(v) for v in a['vector_decimal']], [Decimal(v) for v in b['vector_decimal']]
        x += [Decimal(0)] * (len(y) - len(x))
        assert len(x) == len(y)
        tv = sum(abs(u - v) for u, v in zip(x, y)) / 2
        dp = abs(Decimal(b['ploidy']) - Decimal(a['ploidy']))
        dg = abs(Decimal(b['lambda']) - Decimal(a['lambda']))
        tail = sum(y[-22:])
        stable = (dp < Decimal(str(thresholds['ploidy_abs']))
                  and dg < Decimal(str(thresholds['growth_abs']))
                  + Decimal(str(thresholds['growth_rel'])) * max(abs(Decimal(a['lambda'])), abs(Decimal(b['lambda'])))
                  and tv < Decimal(str(thresholds['distribution_TV']))
                  and tail < Decimal(str(thresholds['upper_22_mass'])))
        return dict(ploidy_abs_change=str(dp), growth_abs_change=str(dg),
                    distribution_TV=str(tv), upper_22_mass_larger=str(tail), diagnostic_stable=stable)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=('prepare', 'preflight', 'collect', 'smoke-audit'))
    p.add_argument('--root', type=Path, required=True)
    a = p.parse_args()
    root = a.root.resolve()
    if a.command == 'prepare': prepare(root)
    elif a.command == 'preflight': preflight(root)
    else: collect(root, a.command == 'smoke-audit')
