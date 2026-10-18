#!/usr/bin/env python3
"""Postprocess fixed outputs: diagonal competition and net-growth balance.

This supplementary script is not part of the frozen HPC evaluation code.
It consumes the unchanged matrices, decimal vectors and saved rate diagnostics.
"""
from decimal import Decimal, localcontext
import argparse
import collections
import gzip
import json
from pathlib import Path
import numpy as np
from boundary import read, table, dump, sha, point_dir


def main(root):
    evidence = read(root / 'task_evidence.tsv')
    output = []
    baseline, baseline_case = None, None
    for row in evidence:
        if row['precision_verified'] != 'True': continue
        path = point_dir(root, row['task_id'])
        rates = read(path / 'transition_rates.tsv')
        matrix = np.loadtxt(gzip.open(path / 'matrix.tsv.gz', 'rt'), delimiter='\t')
        if row['N_MAX'] == '154':
            baseline, baseline_case = matrix.copy(), row['case_id']
        assert baseline_case == row['case_id']
        shared_equal = np.array_equal(matrix[:133, :133], baseline)
        assert shared_equal, ('Changed common matrix block', row['task_id'])
        z = json.loads((path / 'precision_100.json').read_text())['right']
        states = [int(r['N']) for r in rates]
        diag = np.diag(matrix)
        diag_index = int(np.argmax(diag))
        mode_index = states.index(int(row['peak_N']))
        with localcontext() as ctx:
            ctx.prec = 140
            x = [Decimal(v) for v in z['vector_decimal']]
            mom = sum(Decimal(r['lambda_division']) * v for r, v in zip(rates, x))
            death = sum(Decimal(r['mu']) * v for r, v in zip(rates, x))
            colsum = [sum(Decimal.from_float(float(a)) for a in matrix[:, j]) for j in range(len(x))]
            net = sum(a * v for a, v in zip(colsum, x))
            error = abs(net - Decimal(z['lambda']))
            assert error < Decimal('1e-60'), (row['task_id'], str(error))
            # Algebraic balance, using retained canonical matrix columns.
            retained = net + mom + death
            rate_roundoff = Decimal('1e-13') * max(abs(mom), abs(death), abs(net), Decimal('1e-30'))
            assert retained >= -rate_roundoff, ('Negative inflow beyond rate roundoff', row['task_id'])
            rec = {k: row[k] for k in ('task_id', 'case_id', 'source_point_id', 'fit_seed',
                   'N_MAX', 'roles', 'N', 'M', 'replicate', 'parameter', 'sample_id', 'O2_pct',
                   'ploidy_100', 'growth_100', 'peak_N', 'upper_mass', 'upper_22_mass',
                   'old_domain_exceeded_mass')}
            rec.update(max_diagonal_N=states[diag_index], max_diagonal=float(diag[diag_index]),
                       common_block_equal_independent=shared_equal,
                       mode_diagonal=float(diag[mode_index]), mode_mu=rates[mode_index]['mu'],
                       mode_division=rates[mode_index]['lambda_division'],
                       leading_minus_max_diagonal=str(Decimal(z['lambda']) - Decimal.from_float(float(diag[diag_index]))),
                       weighted_division_mother_exit=str(mom), weighted_death=str(death),
                       retained_offspring_inflow_balance_estimate=str(retained),
                       inflow_balance_roundoff_threshold=str(rate_roundoff),
                       inflow_balance_sign_resolved=abs(retained) > rate_roundoff,
                       net_from_column_balance=str(net),
                       balance_residual=str(error))
            output.append(rec)
    table(root / 'mechanism_summary.tsv', output)
    comparisons = read(root / 'ceiling_comparisons.tsv')
    classification = []
    for case in read(root / 'case_manifest.tsv'):
        rs = [r for r in comparisons if r['case_id'] == case['case_id']]
        last = next(r for r in rs if (r['smaller_max'], r['larger_max']) == ('308', '440'))
        total = next(r for r in rs if (r['smaller_max'], r['larger_max']) == ('154', '440'))
        last_resolved = last['comparison_resolved'] == 'True'
        total_resolved = total['comparison_resolved'] == 'True'
        label = 'unresolved'
        if last_resolved and total_resolved:
            if last['diagnostic_stable'] == 'True':
                label = ('boundary_sensitive_then_stable' if Decimal(total['ploidy_abs_change']) >= Decimal('.01')
                         else 'little_ploidy_change_domain_stable')
            else:
                label = 'not_stable_by_domain_diagnostics'
        rec = {k: case[k] for k in ('case_id', 'point_id', 'fit_seed', 'roles', 'N', 'M', 'replicate',
                                    'parameter', 'sample_id', 'O2_pct')}
        rec['domain_classification'] = label
        for prefix, r in (('154_to_440', total), ('308_to_440', last)):
            for k in ('ploidy_abs_change', 'growth_abs_change', 'distribution_TV', 'upper_22_mass_larger', 'diagnostic_stable'):
                rec[prefix + '_' + k] = r.get(k, '')
        classification.append(rec)
    table(root / 'case_boundary_classification.tsv', classification)
    counts = collections.Counter((r['fit_seed'], r['roles'], r['domain_classification']) for r in classification)
    table(root / 'domain_classification_counts.tsv', [dict(fit_seed=k[0], role=k[1], classification=k[2], n_cases=v)
          for k, v in sorted(counts.items())])
    dump(root / 'mechanism_validation.json', dict(status='passed', n_tasks=len(output),
         independent_shared_blocks_equal=all(r['common_block_equal_independent'] for r in output),
         script_sha256=sha(Path(__file__)), task_evidence_sha256=sha(root / 'task_evidence.tsv'),
         explanation='Matrix diagonal is a local self-state rate, not the eigenvalue. Net balance is evaluated with the complete eigenvector.',
         rate_scope='Model builder and saved transition rates are double precision; vector arithmetic is Decimal. Retained inflow is a balance estimate from printed double rates. Cancellation can yield signed roundoff near zero; the sign-resolution flag and roundoff threshold identify these cases. It is not an independent biological estimate.'))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--root', type=Path, required=True)
    main(p.parse_args().root.resolve())
