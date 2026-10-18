#!/usr/bin/env python3
"""Independent saved-vector verification and chromosome-ceiling comparisons."""
import argparse
import collections
from decimal import Decimal, localcontext
import csv
import gzip
import json
from pathlib import Path
import numpy as np
from boundary import read, table, dump, sha, audit, preflight, point_dir, compare, LIMITS


def verify(root):
    selection = preflight(root)
    rows = read(root / 'task_manifest.tsv')
    evidence, independent, solutions = [], [], {}
    with gzip.open(root / 'complete_vectors.tsv.gz', 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=('task_id', 'case_id', 'N_MAX', 'N',
                                'right_100_decimal', 'left_100_decimal'), delimiter='\t')
        writer.writeheader()
        for row in rows:
            path = point_dir(root, row['task_id'])
            result = audit(root, row)
            matrix = np.loadtxt(gzip.open(path / 'matrix.tsv.gz', 'rt'), delimiter='\t')
            nmax = int(row['N_MAX'])
            states = np.arange(22, nmax + 1)
            assert matrix.shape == (len(states), len(states))
            rates = read(path / 'transition_rates.tsv')
            assert [int(r['N']) for r in rates] == list(states)
            with localcontext() as context:
                context.prec = 160
                entries = [(i, j, Decimal.from_float(float(matrix[i, j])))
                           for i, j in zip(*np.nonzero(matrix))]
                for digits in (50, 100):
                    precision = json.loads((path / ('precision_%d.json' % digits)).read_text())
                    for side in ('right', 'left'):
                        z = precision[side]
                        if z.get('error') or not z.get('converged', False):
                            continue
                        x = [Decimal(v) for v in z['vector_decimal']]
                        assert len(x) == len(states) and min(x) > 0
                        assert abs(sum(x) - 1) < Decimal('1e-40' if digits == 50 else '1e-90')
                        lam = Decimal(z['lambda'])
                        mx, rownorm = [Decimal(0)] * len(x), [Decimal(0)] * len(x)
                        for i, j, a in entries:
                            if side == 'left': i, j = j, i
                            mx[i] += a * x[j]
                            rownorm[i] += abs(a)
                        residue = max(abs(a - lam * b) for a, b in zip(mx, x))
                        relative = residue / ((max(rownorm) + abs(lam)) * max(x))
                        bound = abs(max(a / b for a, b in zip(mx, x)) - lam) / max(rownorm)
                        assert relative < Decimal(z['relative_tolerance']) * 10
                        assert bound < Decimal(z['relative_tolerance']) * 10
                        independent.append(dict(task_id=row['task_id'], digits=digits, side=side,
                                                relative_residual=str(relative), relative_upper_gap=str(bound)))
                        if side == 'right':
                            ploidy = sum(Decimal(int(n)) * v for n, v in zip(states, x)) / 22
                            assert abs(ploidy - Decimal(z['ploidy'])) < Decimal('1e-45' if digits == 50 else '1e-95')
                p100 = json.loads((path / 'precision_100.json').read_text())
                if result['precision_verified']:
                    z, left = p100['right'], p100['left']
                    solutions[int(row['task_id'])] = z
                    x = [Decimal(v) for v in z['vector_decimal']]
                    result.update(ploidy_100=z['ploidy'], growth_100=z['lambda'],
                                  upper_mass=str(x[-1]), upper_22_mass=str(sum(x[-22:])),
                                  old_domain_exceeded_mass=str(sum(x[133:])),
                                  peak_N=int(states[max(range(len(x)), key=x.__getitem__)]),
                                  mean_N=str(Decimal(z['ploidy']) * 22))
                    for key in ('upper_dropped_rate', 'upper_misseg_dropped_rate', 'upper_wgd_dropped_rate',
                                'lower_dropped_rate', 'total_boundary_dropped_rate', 'misseg_nonviable_rate'):
                        result['weighted_' + key] = str(sum(Decimal(r[key]) * v for r, v in zip(rates, x)))
                    for n, a, b in zip(states, z['vector_decimal'], left['vector_decimal']):
                        writer.writerow(dict(task_id=row['task_id'], case_id=row['case_id'], N_MAX=nmax,
                                             N=int(n), right_100_decimal=a, left_100_decimal=b))
            evidence.append(result)
            if len(evidence) % 40 == 0:
                print('Independently verified', len(evidence), '/', len(rows), flush=True)
    fields = list(dict.fromkeys(k for r in evidence for k in r))
    table(root / 'task_evidence.tsv', [{k: r.get(k, '') for k in fields} for r in evidence])
    table(root / 'independent_residuals.tsv', independent)
    comparisons = []
    for case in read(root / 'case_manifest.tsv'):
        tasks = [r for r in rows if r['case_id'] == case['case_id']]
        bymax = {int(r['N_MAX']): int(r['task_id']) for r in tasks}
        for small, large in ((154, 220), (220, 308), (308, 440), (154, 440)):
            rec = dict(case_id=case['case_id'], source_point_id=case['point_id'], fit_seed=case['fit_seed'],
                       roles=case['roles'], N=case['N'], M=case['M'], replicate=case['replicate'],
                       parameter=case['parameter'], sample_id=case['sample_id'], O2_pct=case['O2_pct'],
                       smaller_max=small, larger_max=large)
            a, b = solutions.get(bymax[small]), solutions.get(bymax[large])
            rec['comparison_resolved'] = a is not None and b is not None
            if rec['comparison_resolved']:
                rec.update(compare(a, b, selection['diagnostic_thresholds']))
                rec.update(ploidy_smaller=a['ploidy'], ploidy_larger=b['ploidy'],
                           growth_smaller=a['lambda'], growth_larger=b['lambda'])
                # For an unchanged principal block, enlarging a Metzler domain
                # cannot lower its spectral bound. Report/check roundoff scale.
                common_exact = all(r['common_block_max_difference'] == 0 for r in evidence
                                   if r['case_id'] == int(case['case_id']))
                rec['common_block_exact'] = common_exact
                if common_exact:
                    assert Decimal(b['lambda']) >= Decimal(a['lambda']) - Decimal('1e-12')
            comparisons.append(rec)
    fields = list(dict.fromkeys(k for r in comparisons for k in r))
    table(root / 'ceiling_comparisons.tsv', [{k: r.get(k, '') for k in fields} for r in comparisons])
    groups = []
    for seed in ('seed138', 'seed101'):
        for role in ('abnormal', 'ordinary_control'):
            for limit in LIMITS:
                subset = [r for r in evidence if r['fit_seed'] == seed and r['roles'] == role and r['N_MAX'] == limit]
                good = [r for r in subset if r['precision_verified']]
                gs = dict(fit_seed=seed, role=role, N_MAX=limit, n_tasks=len(subset), n_verified=len(good))
                if good:
                    for key in ('ploidy_100', 'growth_100', 'upper_mass', 'upper_22_mass', 'old_domain_exceeded_mass',
                                'weighted_upper_dropped_rate', 'weighted_upper_misseg_dropped_rate', 'weighted_upper_wgd_dropped_rate'):
                        vals = np.array([float(r[key]) for r in good])
                        for stat, value in (('min', vals.min()), ('median', np.median(vals)), ('max', vals.max())):
                            gs[key + '_' + stat] = float(value)
                    gs['n_upper_mass_ge99pct'] = sum(Decimal(r['upper_mass']) >= Decimal('.99') for r in good)
                    gs['n_negative_growth'] = sum(Decimal(r['growth_100']) < 0 for r in good)
                groups.append(gs)
    fields = list(dict.fromkeys(k for r in groups for k in r))
    table(root / 'group_summary.tsv', [{k: r.get(k, '') for k in fields} for r in groups])
    dump(root / 'independent_validation.json', dict(status='passed', n_tasks=len(rows),
         n_precision_verified=sum(r['precision_verified'] for r in evidence), n_residual_checks=len(independent),
         n_leading_root_bound_checks=len(independent), n_comparisons=len(comparisons), new_FAST_samples=0,
         matrix_source='17-digit round-trip matrix parsed to exact binary double; Decimal.from_float, 160 digits',
         state_domain='Only N_MAX changed. Wide transition kernels are flux diagnostics, not eigensolve domains.'))


def figures(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.size': 10, 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    evidence = read(root / 'task_evidence.tsv')
    rows = [r for r in evidence if r['precision_verified'] == 'True']
    colors = {'abnormal': '#D55E00', 'ordinary_control': '#0072B2'}
    out = root / 'figures'; out.mkdir(exist_ok=True)
    def save(fig, name):
        fig.savefig(out / (name + '.pdf'))
        fig.savefig(out / (name + '.png'), dpi=150)
        plt.close(fig)
    def lines(ax, seed, field, log=False):
        for case_id in sorted({r['case_id'] for r in rows if r['fit_seed'] == seed}, key=int):
            rs = sorted([r for r in rows if r['case_id'] == case_id], key=lambda r: int(r['N_MAX']))
            values = [float(r[field]) for r in rs]
            if log: values = np.maximum(values, 1e-120)
            ax.plot([int(r['N_MAX']) for r in rs], values, 'o-', markersize=2, lw=.6,
                    alpha=.35, color=colors[rs[0]['roles']])
        ax.set_xticks(LIMITS); ax.set_xlabel('Chromosome ceiling N_MAX')
        if log: ax.set_yscale('log')
    legend = [Line2D([0], [0], color=c, label=k.replace('_', ' ')) for k, c in colors.items()]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), squeeze=False)
    for j, seed in enumerate(('seed138', 'seed101')):
        lines(axes[0, j], seed, 'ploidy_100'); axes[0, j].set(title=seed, ylabel='Dominant mean ploidy', ylim=(.7, 20.3))
        lines(axes[1, j], seed, 'growth_100'); axes[1, j].set_ylabel('Asymptotic net live-cell growth rate')
    fig.suptitle('Same parameter vectors and O₂: response versus chromosome ceiling')
    fig.legend(handles=legend, loc='lower center', ncol=2, frameon=False)
    fig.tight_layout(rect=(0, .045, 1, .95)); save(fig, 'ceiling_responses')
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8), squeeze=False)
    for j, seed in enumerate(('seed138', 'seed101')):
        lines(axes[0, j], seed, 'upper_22_mass', True)
        axes[0, j].axhline(1e-4, color='0.5', ls=':', lw=1)
        axes[0, j].set(title=seed, ylabel='Mass in highest 22 states')
        lines(axes[1, j], seed, 'old_domain_exceeded_mass')
        axes[1, j].set(ylabel='Mass above original N=154', ylim=(-.03, 1.03))
    fig.suptitle('Boundary occupancy and mass beyond the original domain\nTop-panel zeros display at 10⁻¹²⁰; dotted line is the diagnostic threshold')
    fig.legend(handles=legend, loc='lower center', ncol=2, frameon=False)
    fig.tight_layout(rect=(0, .045, 1, .92)); save(fig, 'boundary_mass')
    selection = json.loads((root / 'selection.json').read_text())
    tasks = read(root / 'task_manifest.tsv')
    cases = sorted({r['case_id'] for r in tasks if int(r['task_id']) in selection['smoke_tasks']}, key=int)
    palette = {154: '#777777', 220: '#0072B2', 308: '#009E73', 440: '#CC79A7'}
    for log in (False, True):
        fig, axes = plt.subplots(2, 3, figsize=(14, 8), squeeze=False)
        for ax, case_id in zip(axes.flat, cases):
            chosen = sorted([r for r in rows if r['case_id'] == case_id], key=lambda r: int(r['N_MAX']))
            if not chosen: continue
            for r in chosen:
                z = json.loads((point_dir(root, r['task_id']) / 'precision_100.json').read_text())['right']
                x = np.arange(22, int(r['N_MAX']) + 1)
                if log:
                    y = [float(Decimal(v).log10()) for v in z['vector_decimal']]
                else: y = [float(v) for v in z['vector_decimal']]
                ax.plot(x, y, lw=1.2, color=palette[int(r['N_MAX'])], label='N_MAX=' + r['N_MAX'])
            r = chosen[0]
            ax.axvline(154, color='0.3', ls=':', lw=.8)
            ax.set(title='%s | point %s | %s\n%s | O₂ %.4g%%' %
                   (r['fit_seed'], r['source_point_id'], r['roles'].replace('_', ' '), r['parameter'], float(r['O2_pct'])),
                   xlabel='Chromosome state N', ylabel='log₁₀ normalized mass' if log else 'Normalized mass', xlim=(20, 445))
        axes.flat[0].legend(fontsize=8, frameon=False)
        fig.suptitle('Complete 100-digit dominant right eigenvectors: six fixed smoke cases')
        fig.tight_layout(rect=(0, 0, 1, .94), h_pad=2)
        save(fig, 'state_distributions_log' if log else 'state_distributions')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('command', choices=('verify', 'figures'))
    p.add_argument('--root', type=Path, required=True)
    a = p.parse_args(); globals()[a.command](a.root.resolve())
