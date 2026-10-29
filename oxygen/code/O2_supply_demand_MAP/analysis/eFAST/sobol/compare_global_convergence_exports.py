#!/usr/bin/env python3
"""Validate and compare downloaded N65536/N131072 global Sobol indices.

Only exported sensitivity arrays are read. No model, SALib estimator, cluster
submission, neighborhood, or original plotting module is imported or invoked.
"""
import argparse
import datetime
import hashlib
import importlib.metadata
import itertools
import json
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

NS = (65536, 131072)
OUTPUTS = ('dominant_mean_ploidy', 'dominant_growth_rate')
LABELS = dict(zip(OUTPUTS, ('Ploidy', 'Net live-cell growth rate')))
FOCAL = frozenset(('p_mis_base', 'mu_hp', 'gamma_mu'))
KEY = ['O2_pct', 'output', 'index', 'parameter_i', 'parameter_j']
COLORS = {65536: '#6b819b', 131072: '#5b247a'}
RANGE_LIMIT, DELTA_LIMIT = .05, .02


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 ** 2), b''):
            h.update(block)
    return h.hexdigest()


def save(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def table(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep='\t', index=False, float_format='%.17g')


def category(frame):
    focal = frame.parameter_i.isin(FOCAL) | frame.parameter_j.isin(FOCAL)
    return np.where(frame['index'] == 'S2', np.where(focal, 'S2_focal', 'S2_other'), frame['index'])


def selection(frame, kind):
    if kind == 'S2_all':
        return frame['index'] == 'S2'
    if kind == 'S2_focal':
        return (frame['index'] == 'S2') & (frame.parameter_i.isin(FOCAL) | frame.parameter_j.isin(FOCAL))
    return frame['index'] == kind


def band_name(oxygen):
    return np.where(oxygen <= 1, 'Low [0,1]', np.where(oxygen < 3, 'Mid (1,3)', 'High [3,5]'))


def metric_label(metric):
    return ' x '.join(x for x in metric[1:] if x)


def validate(root, collection, out):
    manifest = json.loads((collection / 'remote_manifest.json').read_text())
    receipt = json.loads((collection / 'collection_receipt.json').read_text())
    assert receipt['status'] == 'complete' and receipt['all_file_hashes_verified']
    assert sha(collection / 'remote_manifest.json') == receipt['manifest_sha256']
    verified = []
    for record in manifest['files']:
        base = root if record['path'].split('/')[0] in ['N%d' % n for n in NS] else collection / 'provenance'
        path = base / record['path']
        assert path.is_file() and path.stat().st_size == record['size']
        assert sha(path) == record['sha256'], str(path)
        verified.append(dict(path=str(path), sha256=record['sha256']))
    reference = collection / 'provenance/reference/parameter_bounds.tsv'
    bounds = pd.read_csv(reference, sep='\t', keep_default_na=False)
    external = manifest['external_evidence'][0]
    assert sha(reference) == external['sha256']
    values, summaries, frames, originals = {}, {}, {}, {}
    metrics, oxygen, active = None, None, None
    checks = dict(export_hashes_verified=len(verified), external_bounds_sha256_verified=True,
                  model_calls=0, neighborhood_files_read=0, resolutions={},
                  tolerance_for_export_float_roundtrip=1e-12)
    for n in NS:
        directory = root / ('N%d' % n)
        completion = json.loads((directory / 'completion.json').read_text())
        assert completion['status'] == 'complete' and completion['N'] == n
        assert completion['n_repeats'] == 5 and completion['n_oxygen'] == 201
        assert completion['n_indices'] == 47838
        assert sha(directory / 'indices_summary.tsv.gz') == completion['summary_sha256']
        df = pd.read_csv(directory / 'indices_summary.tsv.gz', sep='\t', keep_default_na=False)
        df['O2_pct'] = df.O2_pct.round(12)
        assert len(df) == 47838 and not df.duplicated(KEY).any()
        assert set(df.scope) == {'global'} and set(df.N) == {n}
        assert set(df.output) == set(OUTPUTS)
        assert set(df.resolution_low) == {n // 2} and set(df.resolution_high) == {n}
        current_metrics = list(df[(df.O2_pct == 0) & (df.output == OUTPUTS[0])][['index', 'parameter_i', 'parameter_j']].itertuples(index=False, name=None))
        assert len(current_metrics) == len(set(current_metrics)) == 119
        if metrics is None:
            metrics = current_metrics
            active = [p for idx, p, q in metrics if idx == 'ST']
            assert len(active) == len(set(active)) == 14 and FOCAL.issubset(set(active))
            assert len(bounds) == 14 and set(bounds['parameter']) == set(active)
            pairs = {frozenset((p, q)) for idx, p, q in metrics if idx == 'S2'}
            assert pairs == {frozenset(pair) for pair in itertools.combinations(active, 2)}
            assert sum(idx == 'S2' and bool({p, q} & FOCAL) for idx, p, q in metrics) == 36
        else:
            assert metrics == current_metrics
        with np.load(directory / 'indices_repeats.npz', allow_pickle=False) as archive:
            v = archive['values'].copy()
            current_oxygen = archive['oxygen'].copy()
            assert np.array_equal(archive['resolutions'], [n // 2, n])
            assert v.shape == (2, 5, 201, 2, 119)
            assert archive['bootstrap_half_width'].shape == v.shape
            assert archive['variance'].shape == (2, 5, 201, 2)
            assert np.isfinite(archive['variance']).all() and (archive['variance'] > 0).all()
        assert np.allclose(current_oxygen, np.arange(201) / 40, rtol=0, atol=1e-12)
        if oxygen is None:
            oxygen = current_oxygen
        else:
            assert np.array_equal(oxygen, current_oxygen)
        # Independent calculation: do not call production f.diagnostics.
        mean = np.add.reduce(v[1], axis=0) / 5
        previous_mean = np.add.reduce(v[0], axis=0) / 5
        spread = np.max(v[1], axis=0) - np.min(v[1], axis=0)
        delta = np.abs(mean - previous_mean)
        finite = np.isfinite(v[1]).all(axis=0) & np.isfinite(delta)
        passed = finite & (spread <= RANGE_LIMIT) & (delta <= DELTA_LIMIT)
        status = np.where(~finite, 0, np.where(passed, 1, 2)).astype(np.uint8)
        with np.load(directory / 'indices_summary.npz', allow_pickle=False) as archive:
            for name, expected in [('mean', mean), ('repeat_range', spread), ('resolution_delta', delta)]:
                assert np.allclose(archive[name], expected, rtol=0, atol=1e-12, equal_nan=True), name
            assert np.array_equal(archive['status'], status)
            assert np.array_equal(archive['oxygen'], oxygen)
        records = [(round(float(o), 12), output, idx, pi, pj) for o in oxygen for output in OUTPUTS for idx, pi, pj in metrics]
        calculated = pd.DataFrame(records, columns=KEY)
        calculated['N'] = n
        calculated['value_mean'] = mean.ravel()
        calculated['repeat_min'] = np.min(v[1], axis=0).ravel()
        calculated['repeat_max'] = np.max(v[1], axis=0).ravel()
        calculated['repeat_range'] = spread.ravel()
        calculated['resolution_delta'] = delta.ravel()
        calculated['status_code'] = status.ravel()
        calculated['convergence_status'] = np.asarray(['undefined', 'passed', 'not_passed'])[status.ravel()]
        calculated['O2_band'] = band_name(calculated.O2_pct)
        joined = calculated.merge(df, on=KEY + ['N'], validate='one_to_one', suffixes=('_calculated', '_export'))
        assert len(joined) == 47838
        for column in ['value_mean', 'repeat_range', 'resolution_delta']:
            assert np.allclose(joined[column + '_calculated'], joined[column + '_export'], rtol=0, atol=1e-12, equal_nan=True), column
        assert (joined.convergence_status_calculated == joined.convergence_status_export).all()
        failures = pd.read_csv(directory / 'not_passed.tsv.gz', sep='\t', keep_default_na=False)
        failures['O2_pct'] = failures.O2_pct.round(12)
        expected_failures = df[df.convergence_status != 'passed']
        assert len(failures) == len(expected_failures)
        assert set(map(tuple, failures[KEY].to_numpy())) == set(map(tuple, expected_failures[KEY].to_numpy()))
        classes = pd.read_csv(directory / 'classification.tsv', sep='\t', keep_default_na=False)
        assert len(classes) == 119 and not classes.duplicated(['index', 'parameter_i', 'parameter_j']).any()
        ordering = pd.read_csv(directory / 'row_order.tsv', sep='\t', keep_default_na=False)
        assert ordering.groupby('panel').size().to_dict() == {'S1': 14, 'S2': 91, 'S2_focal': 36, 'ST': 14}
        values[n] = v
        summaries[n] = dict(mean=mean, spread=spread, delta=delta, status=status)
        frames[n], originals[n] = calculated, df
        checks['resolutions'][str(n)] = dict(rows=47838, repeats=5, oxygen_points=201, parameters=14,
            all_pairs=91, focal_pairs=36, n_passed=int(passed.sum()), n_not_passed=int((status == 2).sum()),
            n_undefined=int((status == 0).sum()), npz_and_tsv_match_independent_diagnostics=True)
    assert np.array_equal(values[65536][1], values[131072][0], equal_nan=True)
    checks['N65536_reused_indices_bitwise_equal'] = True
    checks['reused_indices_shape'] = list(values[65536][1].shape)
    checks['input_sha256'] = {str(p.relative_to(root)): sha(p) for n in NS for p in [root / ('N%d/indices_repeats.npz' % n), root / ('N%d/indices_summary.npz' % n), root / ('N%d/indices_summary.tsv.gz' % n)]}
    save(out / 'validation_data.json', dict(verified_utc=now(), **checks))
    return dict(values=values, summaries=summaries, frames=frames, originals=originals,
                metrics=metrics, oxygen=oxygen, active=active, checks=checks)


def analysis(root, out, data):
    frames, oxygen, metrics = data['frames'], data['oxygen'], data['metrics']
    combined = pd.concat([frames[n] for n in NS], ignore_index=True)
    summary_rows, band_rows, metric_rows = [], [], []
    for n in NS:
        for output in OUTPUTS:
            for kind in ('S1', 'ST', 'S2_all', 'S2_focal'):
                sub = frames[n][(frames[n].output == output) & selection(frames[n], kind)]
                summary_rows.append(dict(N=n, output=output, index_group=kind, n_total=len(sub),
                    n_passed=int((sub.status_code == 1).sum()), n_not_passed=int((sub.status_code == 2).sum()),
                    n_undefined=int((sub.status_code == 0).sum()), pass_fraction=float((sub.status_code == 1).mean()),
                    n_range_failed=int((sub.repeat_range > RANGE_LIMIT).sum()),
                    n_delta_failed=int((sub.resolution_delta > DELTA_LIMIT).sum()),
                    repeat_range_median=float(sub.repeat_range.median()), repeat_range_p95=float(sub.repeat_range.quantile(.95)),
                    repeat_range_max=float(sub.repeat_range.max()), delta_median=float(sub.resolution_delta.median()),
                    delta_p95=float(sub.resolution_delta.quantile(.95)), delta_max=float(sub.resolution_delta.max())))
                for band, b in sub.groupby('O2_band', sort=False):
                    band_rows.append(dict(N=n, output=output, index_group=kind, O2_band=band,
                        n_total=len(b), n_passed=int((b.status_code == 1).sum()),
                        pass_fraction=float((b.status_code == 1).mean()), repeat_range_max=float(b.repeat_range.max()),
                        delta_max=float(b.resolution_delta.max())))
            focal = frames[n][(frames[n].output == output) & (selection(frames[n], 'ST') | selection(frames[n], 'S2_focal'))]
            for (idx, pi, pj), sub in focal.groupby(['index', 'parameter_i', 'parameter_j'], sort=False):
                sub = sub.sort_values('O2_pct')
                row = dict(N=n, output=output, index=idx, parameter_i=pi, parameter_j=pj,
                    n_oxygen=len(sub), n_passed=int((sub.status_code == 1).sum()),
                    repeat_range_max=float(sub.repeat_range.max()), resolution_delta_max=float(sub.resolution_delta.max()),
                    peak_value=float(sub.value_mean.max()), peak_O2_pct=float(sub.loc[sub.value_mean.idxmax(), 'O2_pct']))
                for lo, hi, label in [(0, 1, 'low'), (1, 3, 'mid'), (3, 5, 'high')]:
                    b = sub[(sub.O2_pct >= lo) & (sub.O2_pct <= hi)]
                    row[label + '_integral_mean'] = float(np.trapezoid(b.value_mean, b.O2_pct) / (hi - lo))
                metric_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    metric_summary = pd.DataFrame(metric_rows)
    table(summary, out / 'diagnostic_summary.tsv')
    table(pd.DataFrame(band_rows), out / 'diagnostic_summary_by_O2_band.tsv')
    table(metric_summary, out / 'parameter_and_pair_summary.tsv')
    focus = combined[selection(combined, 'ST') | selection(combined, 'S2_focal')].copy()
    table(focus, out / 'curves_source.tsv.gz')
    columns = KEY + ['value_mean', 'repeat_min', 'repeat_max', 'repeat_range', 'resolution_delta', 'convergence_status']
    comparison = frames[65536][columns].merge(frames[131072][columns], on=KEY, suffixes=('_N65536', '_N131072'), validate='one_to_one')
    comparison['mean_change'] = comparison.value_mean_N131072 - comparison.value_mean_N65536
    comparison['absolute_mean_change'] = comparison.mean_change.abs()
    comparison['repeat_range_change'] = comparison.repeat_range_N131072 - comparison.repeat_range_N65536
    comparison['newly_passed'] = (comparison.convergence_status_N65536 != 'passed') & (comparison.convergence_status_N131072 == 'passed')
    comparison['lost_pass'] = (comparison.convergence_status_N65536 == 'passed') & (comparison.convergence_status_N131072 != 'passed')
    table(comparison, out / 'all_index_comparison.tsv.gz')
    focus_comparison = comparison[selection(comparison, 'ST') | selection(comparison, 'S2_focal')]
    table(focus_comparison, out / 'ST_and_36_S2_comparison.tsv.gz')
    failures = comparison[comparison.convergence_status_N65536 != 'passed']
    table(failures, out / 'previously_not_passed_tracking.tsv')
    rows = []
    for n in NS:
        for rep in range(5):
            for oi, o in enumerate(oxygen):
                for k, output in enumerate(OUTPUTS):
                    for m, (idx, pi, pj) in enumerate(metrics):
                        if idx == 'ST' or (idx == 'S2' and bool({pi, pj} & FOCAL)):
                            rows.append((n, rep + 1, o, output, idx, pi, pj, data['values'][n][1, rep, oi, k, m]))
    table(pd.DataFrame(rows, columns=['N', 'repeat', 'O2_pct', 'output', 'index', 'parameter_i', 'parameter_j', 'value']), out / 'five_repeat_source.tsv.gz')
    historical = pd.read_csv(root / 'N131072/resolution_comparison/curves_source.tsv.gz', sep='\t', keep_default_na=False)
    historical['O2_pct'] = historical.O2_pct.round(12)
    assert not historical.duplicated(['N'] + KEY).any()
    assert sorted(historical.N.unique()) == [4096, 8192, 16384, 32768, 65536, 131072]
    hist_rows = []
    for (n, output), sub in historical.groupby(['N', 'output'], sort=True):
        for kind in ('ST', 'S2_focal'):
            selected = sub[selection(sub, kind)]
            assert len(selected) == 201 * (14 if kind == 'ST' else 36)
            hist_rows.append(dict(N=n, output=output, index_group=kind, n_total=len(selected),
                n_passed=int((selected.diagnostic_status == 1).sum()), pass_fraction=float((selected.diagnostic_status == 1).mean())))
    history = pd.DataFrame(hist_rows)
    for n in NS:
        h = historical[historical.N == n].merge(frames[n], on=['N'] + KEY, validate='one_to_one', suffixes=('_history', '_current'))
        assert len(h) == 20100 and (h.diagnostic_status == h.status_code).all()
        assert np.allclose(h.value_mean_history, h.value_mean_current, rtol=0, atol=1e-12)
    table(history, out / 'historical_ST_and_focal_S2_pass_rates.tsv')
    ranks = []
    for output in OUTPUTS:
        for kind in ('ST', 'S2_focal'):
            sub = focus_comparison[(focus_comparison.output == output) & selection(focus_comparison, kind)]
            for o, b in sub.groupby('O2_pct'):
                a = b.value_mean_N65536.rank().to_numpy()
                c = b.value_mean_N131072.rank().to_numpy()
                rho = float(np.corrcoef(a, c)[0, 1]) if np.std(a) and np.std(c) else np.nan
                ranks.append(dict(output=output, index_group=kind, O2_pct=o, spearman_rank_correlation=rho,
                    max_absolute_mean_change=float(b.absolute_mean_change.max())))
    ranking = pd.DataFrame(ranks)
    table(ranking, out / 'ranking_stability.tsv')
    data.update(summary=summary, metric_summary=metric_summary, comparison=comparison,
                focus=focus, failures=failures, history=history, ranking=ranking)


def oxygen_axis(ax):
    ax.set_xscale('symlog', base=10, linthresh=.025, linscale=1)
    ax.set_xlim(0, 5)
    ax.set_xticks([0, .025, .1, .5, 1, 3, 5], ['0', '.025', '.1', '.5', '1', '3', '5'])
    ax.set_xlabel('Fixed oxygen (%)')
    for x in (1, 3):
        ax.axvline(x, color='#777777', linestyle='--', linewidth=.7, zorder=0)
    transform = ax.xaxis.get_transform()
    for lo, hi, label in [(0, 1, 'Low'), (1, 3, 'Mid'), (3, 5, 'High')]:
        mid = transform.inverted().transform(np.mean(transform.transform(np.array([lo, hi]))))
        ax.text(mid, 1.015, label, transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=7)


def plot_all(out, data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullLocator
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.titlesize': 10, 'axes.labelsize': 9, 'xtick.labelsize': 8,
        'ytick.labelsize': 8, 'axes.spines.top': False, 'axes.spines.right': False,
        'pdf.fonttype': 42, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})
    figures = out / 'figures'
    figures.mkdir(exist_ok=True)
    meta = dict(Creator='compare_global_convergence_exports.py', CreationDate=None, ModDate=None)
    legend = [Line2D([], [], color=COLORS[n], linestyle='--' if n == 65536 else '-',
                     linewidth=1.5, label='N=%d' % n) for n in NS]
    def export(fig, name):
        fig.savefig(figures / (name + '.pdf'), bbox_inches='tight', metadata=meta)
        fig.savefig(figures / (name + '.png'), dpi=150, bbox_inches='tight')
        plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5))
    group_colors = {'ST': '#5b247a', 'S2_focal': '#c57629'}
    for k, output in enumerate(OUTPUTS):
        for kind in ('ST', 'S2_focal'):
            color = group_colors[kind]
            h = data['history'][(data['history'].output == output) & (data['history'].index_group == kind)].sort_values('N')
            axes[k, 0].plot(h.N, h.pass_fraction * 100, color=color, marker='o', label='ST' if kind == 'ST' else '36 focal S2 pairs')
            s = data['summary'][(data['summary'].output == output) & (data['summary'].index_group == kind)].sort_values('N')
            for j, prefix in [(1, 'repeat_range'), (2, 'delta')]:
                axes[k, j].plot(s.N, s[prefix + '_max'], color=color, marker='o', label=kind.replace('_focal', ' focal') + ' maximum')
                axes[k, j].plot(s.N, s[prefix + '_p95'], color=color, linestyle='--', marker='s', markersize=4, label=kind.replace('_focal', ' focal') + ' 95th percentile')
        axes[k, 0].set(xscale='log', ylim=(0, 105), ylabel='Passed cells (%)', title=LABELS[output] + ': diagnostic pass rate')
        ticks = [4096, 8192, 16384, 32768, 65536, 131072]
        axes[k, 0].set_xticks(ticks, [str(n) for n in ticks], rotation=30, ha='right')
        axes[k, 0].legend(fontsize=8, loc='lower right')
        for j, threshold, title, ylabel in [(1, RANGE_LIMIT, 'Five-repeat range', 'max(repeats) - min(repeats)'),
                                           (2, DELTA_LIMIT, 'Adjacent-resolution change', 'Absolute change in five-repeat mean')]:
            ax = axes[k, j]
            ax.axhline(threshold, color='#444444', linestyle=':', linewidth=1.1, label='Diagnostic limit %.2f' % threshold)
            ax.set(xscale='log', title=LABELS[output] + ': ' + title, ylabel=ylabel)
            ax.set_xticks(NS, [str(n) for n in NS])
            ax.set_xlim(60000, 145000)
            ax.set_ylim(0, max(ax.get_ylim()[1], threshold * 1.15))
            ax.legend(fontsize=6.8, loc='upper right')
        for j in range(3):
            axes[k, j].xaxis.set_minor_locator(NullLocator())
            axes[k, j].set_xlabel('Base sample count N')
            axes[k, j].grid(axis='y', alpha=.18)
            axes[k, j].text(-.13, 1.08, chr(65 + k * 3 + j), transform=axes[k, j].transAxes, weight='bold', fontsize=12)
    fig.suptitle('Global Sobol convergence: ST and 36 focal S2 pairs', fontsize=15, y=.985)
    fig.text(.5, .018, 'Five independent scramble repeats; 201 fixed oxygen levels; thresholds: repeat range <=0.05 and adjacent-N mean change <=0.02.\n'
             'N65536 compares with N32768; N131072 compares with N65536. Undefined cells count as not passed.', ha='center', fontsize=9)
    fig.subplots_adjust(left=.075, right=.985, bottom=.15, top=.90, wspace=.32, hspace=.55)
    export(fig, 'convergence_overview')

    def curve(ax, metric, output):
        m, k = data['metrics'].index(metric), OUTPUTS.index(output)
        for n in NS:
            v = data['values'][n][1, :, :, k, m]
            ax.fill_between(data['oxygen'], v.min(axis=0), v.max(axis=0), color=COLORS[n], alpha=.13, linewidth=0)
            ax.plot(data['oxygen'], v.mean(axis=0), color=COLORS[n], linestyle='--' if n == 65536 else '-', linewidth=1.3)
        ax.axhline(0, color='#777777', linewidth=.6, zorder=0)
        oxygen_axis(ax)
        ax.set_title(metric_label(metric), pad=24)
        ax.set_ylabel(metric[0])
        ax.grid(axis='y', alpha=.13)

    bad = data['failures'][['output', 'index', 'parameter_i', 'parameter_j']].drop_duplicates()
    if len(bad):
        with PdfPages(figures / 'previously_unstable_ST.pdf', metadata=meta) as pdf:
            for page_number, row in enumerate(bad.itertuples(index=False), 1):
                metric = (row.index, row.parameter_i, row.parameter_j)
                m, k = data['metrics'].index(metric), OUTPUTS.index(row.output)
                fig, axes = plt.subplots(3, 1, figsize=(10, 11))
                curve(axes[0], metric, row.output)
                for n in NS:
                    s = data['summaries'][n]
                    style = '--' if n == 65536 else '-'
                    axes[1].plot(data['oxygen'], s['spread'][:, k, m], color=COLORS[n], linestyle=style, label='N=%d' % n)
                    axes[2].plot(data['oxygen'], s['delta'][:, k, m], color=COLORS[n], linestyle=style, label='%d -> %d' % (n // 2, n))
                old_status = data['summaries'][65536]['status'][:, k, m]
                ids = old_status != 1
                axes[1].scatter(data['oxygen'][ids], data['summaries'][65536]['spread'][ids, k, m], s=14, marker='x', color=COLORS[65536], label='Previously not passed')
                for ax, threshold, ylabel in [(axes[1], RANGE_LIMIT, 'Five-repeat range'), (axes[2], DELTA_LIMIT, 'Adjacent-N mean change')]:
                    oxygen_axis(ax)
                    ax.axhline(threshold, color='#444444', linestyle=':', linewidth=1, label='Diagnostic limit %.2f' % threshold)
                    ax.set_ylabel(ylabel)
                    ax.set_ylim(0, max(ax.get_ylim()[1], threshold * 1.12))
                    ax.grid(axis='y', alpha=.15)
                    ax.legend(fontsize=8)
                axes[0].legend(handles=legend, fontsize=9)
                fig.suptitle('Previously unstable %s: %s / %s' % (row.index, row.parameter_i, LABELS[row.output]), fontsize=14)
                fig.text(.5, .025, 'Shaded bands are five-repeat min-max ranges, not confidence intervals.\nOxygen axis: linear from 0 to 0.025%; logarithmic above 0.025%. Dashed boundaries: 1% and 3%.', ha='center', fontsize=9)
                fig.subplots_adjust(top=.91, bottom=.10, left=.12, right=.98, hspace=.54)
                pdf.savefig(fig, bbox_inches='tight')
                fig.savefig(figures / ('previously_unstable_ST_%02d.png' % page_number), dpi=150, bbox_inches='tight')
                plt.close(fig)

    top_rows = {}
    for output in OUTPUTS:
        top = data['metric_summary'][(data['metric_summary'].N == 131072) & (data['metric_summary'].output == output) & (data['metric_summary']['index'] == 'S2')]
        top = top.sort_values(['peak_value', 'parameter_i', 'parameter_j'], ascending=[False, True, True]).head(6)
        top_rows[output] = list(top[['index', 'parameter_i', 'parameter_j']].itertuples(index=False, name=None))
    fig, axes = plt.subplots(6, 2, figsize=(13, 17))
    for k, output in enumerate(OUTPUTS):
        for row, metric in enumerate(top_rows[output]):
            curve(axes[row, k], metric, output)
    fig.legend(handles=legend, loc='upper center', bbox_to_anchor=(.5, .967), ncol=2)
    fig.suptitle('Focal S2: six pairs per output with the largest N131072 peak', y=.99, fontsize=14)
    fig.text(.5, .017, 'Selection: descending maximum five-repeat mean S2 over 201 oxygen levels, separately for each output.\n'
             'All 36 pairs are in the appendix. Bands: five-repeat min-max, not CIs. Negative estimates: sampling noise.\n'
             'Each panel has its own y scale; both N values share that scale. Oxygen boundaries: 1% and 3%.', ha='center', fontsize=9)
    fig.subplots_adjust(top=.87, bottom=.09, left=.08, right=.98, hspace=.95, wspace=.25)
    for k, output in enumerate(OUTPUTS):
        position = axes[0, k].get_position()
        fig.text((position.x0 + position.x1) / 2, .925, LABELS[output], ha='center', va='center', fontsize=11, weight='bold')
    export(fig, 'focal_S2_top_pairs')

    for kind, name in [('ST', 'ST_all_parameters'), ('S2', 'S2_all_36_focal_pairs')]:
        chosen = [m for m in data['metrics'] if m[0] == kind and (kind == 'ST' or bool(set(m[1:]) & FOCAL))]
        if kind == 'S2':
            peaks = data['metric_summary'][(data['metric_summary'].N == 131072) & (data['metric_summary'].output == OUTPUTS[0]) & (data['metric_summary']['index'] == 'S2')]
            chosen = list(peaks.sort_values(['peak_value', 'parameter_i', 'parameter_j'], ascending=[False, True, True])[['index', 'parameter_i', 'parameter_j']].itertuples(index=False, name=None))
        with PdfPages(figures / (name + '.pdf'), metadata=meta) as pdf:
            for page, begin in enumerate(range(0, len(chosen), 4), 1):
                current = chosen[begin:begin + 4]
                fig, axes = plt.subplots(4, 2, figsize=(13, 12))
                for row in range(4):
                    for k, output in enumerate(OUTPUTS):
                        if row < len(current):
                            curve(axes[row, k], current[row], output)
                        else:
                            axes[row, k].set_visible(False)
                fig.suptitle('%s: N65536 vs N131072 | page %d/%d' % (name.replace('_', ' '), page, (len(chosen) + 3) // 4), fontsize=14, y=.99)
                fig.legend(handles=legend, loc='upper center', bbox_to_anchor=(.5, .966), ncol=2)
                fig.text(.5, .02, 'Bands: min-max across five scramble repeats, not confidence intervals. Negative estimates are sampling noise.\n'
                         'Independent y scales by panel; shared y scale for both N values within each panel. O2 boundaries: 1% and 3%.', ha='center', fontsize=9)
                fig.subplots_adjust(top=.84, bottom=.10, left=.08, right=.98, hspace=.95, wspace=.25)
                for k, output in enumerate(OUTPUTS):
                    position = axes[0, k].get_position()
                    fig.text((position.x0 + position.x1) / 2, .915, LABELS[output], ha='center', va='center', fontsize=11, weight='bold')
                pdf.savefig(fig, bbox_inches='tight')
                fig.savefig(figures / ('%s_%02d.png' % (name, page)), dpi=140, bbox_inches='tight')
                plt.close(fig)


def validate_assets(root, out):
    from PIL import Image, ImageStat
    from pypdf import PdfReader
    rows = []
    for directory in [root / 'N65536', root / 'N131072', out / 'figures']:
        for path in sorted(directory.rglob('*')):
            if path.suffix not in ('.pdf', '.png'):
                continue
            record = dict(path=str(path), size=path.stat().st_size, sha256=sha(path))
            if path.suffix == '.pdf':
                reader = PdfReader(str(path), strict=True)
                assert len(reader.pages) > 0
                record['pages'] = len(reader.pages)
                record['text_characters'] = sum(len(page.extract_text() or '') for page in reader.pages)
                assert record['text_characters'] > 20
                assert all(page.get_contents() is not None for page in reader.pages)
                record['check'] = 'parsed_pages_content_and_text'
            else:
                with Image.open(path) as image:
                    image.verify()
                with Image.open(path) as image:
                    record['pixels'] = list(image.size)
                    small = image.convert('L')
                    small.thumbnail((128, 128))
                    record['thumbnail_stddev'] = ImageStat.Stat(small).stddev[0]
                    assert record['thumbnail_stddev'] > 1
                record['check'] = 'decoded_crc_and_nonblank_pixels'
            rows.append(record)
    save(out / 'validation_assets.json', dict(verified_utc=now(), n_assets=len(rows), files=rows,
        scope='All collected PDF/PNG and all new comparison PDF/PNG; structural checks, not a claim of visual inspection of every original page.'))
    return rows


def render_pdf_checks(root, out):
    """Render all collected/new PDFs; retain contact sheets and exact commands."""
    from concurrent.futures import ThreadPoolExecutor
    import tempfile
    import textwrap
    from PIL import Image, ImageDraw, ImageFont, ImageStat
    from pypdf import PdfReader
    executable = shutil.which('pdftoppm')
    if not executable:
        raise RuntimeError('pdftoppm is required for --render-pdfs')
    temporary = Path(tempfile.mkdtemp(prefix='global_sobol_pdf_qa_'))
    paths = sorted(p for directory in [root / 'N65536', root / 'N131072', out / 'figures'] for p in directory.rglob('*.pdf'))
    previous = out / 'pdf_render_validation.json'
    reusable = {r['pdf']: r for r in json.loads(previous.read_text())['records']} if previous.exists() else {}
    def render(path):
        old = reusable.get(str(path))
        if old and old['sha256'] == sha(path) and all(Path(p).is_file() for p in old['rendered_pages']):
            return dict(old, reused_previous_render=True)
        prefix = temporary / '__'.join(path.relative_to(root).with_suffix('').parts)
        command = [executable, '-scale-to', '2000' if out in path.parents else '1600', '-png', str(path), str(prefix)]
        p = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=180)
        assert p.returncode == 0, p.stderr
        pages = sorted(temporary.glob(prefix.name + '-*.png'))
        assert len(pages) == len(PdfReader(str(path)).pages)
        for page in pages:
            with Image.open(page) as image:
                image.load()
                small = image.convert('L')
                small.thumbnail((128, 128))
                assert ImageStat.Stat(small).stddev[0] > 1
        return dict(pdf=str(path), sha256=sha(path), pages=len(pages), command=command,
                    stderr=p.stderr, rendered_pages=[str(p) for p in pages], reused_previous_render=False)
    with ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(render, paths))
    qa = out / 'validation'
    qa.mkdir(exist_ok=True)
    for label, current in [('new_comparison', True), ('original_exports', False)]:
        entries = [(record, page, i) for record in records if (out in Path(record['pdf']).parents) == current
                   for i, page in enumerate(record['rendered_pages'], 1)]
        width, height, columns = 500, 440, 4
        sheet = Image.new('RGB', (width * columns, height * ((len(entries) + columns - 1) // columns)), 'white')
        draw = ImageDraw.Draw(sheet)
        font = ImageFont.load_default(size=13)
        for index, (record, page, number) in enumerate(entries):
            x, y = (index % columns) * width, (index // columns) * height
            with Image.open(page) as image:
                tile = image.convert('RGB')
                tile.thumbnail((width - 16, height - 48))
                sheet.paste(tile, (x + (width - tile.width) // 2, y + 4))
            title = str(Path(record['pdf']).relative_to(root)) + ' p%d' % number
            draw.multiline_text((x + 5, y + height - 42), '\n'.join(textwrap.wrap(title, 64)), fill='#222222', font=font)
        sheet.save(qa / (label + '_contact_sheet.png'))
    save(out / 'pdf_render_validation.json', dict(verified_utc=now(), n_pdfs=len(records),
        n_pages=sum(r['pages'] for r in records), renderer=executable, records=records,
        note='All PDFs rendered and every raster decoded/nonblank. Human visual inspection is documented separately in visual_review.md.'))
    return records


def intervals(values):
    values = sorted(values)
    groups = []
    for o in values:
        if not groups or o - groups[-1][-1] > .025000001:
            groups.append([o])
        else:
            groups[-1].append(o)
    return '、'.join(('%.3f%%' % g[0]) if len(g) == 1 else ('%.3f–%.3f%%' % (g[0], g[-1])) for g in groups)


def markdown_table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(str(value) for value in row) + ' |' for row in rows])


def write_report(root, collection, out, data):
    summary, detail = data['summary'], data['comparison']
    receipt = json.loads((collection / 'collection_receipt.json').read_text())
    lines = ['# 全域 Sobol：N65536 与 N131072 收敛比较', '',
        '## 结论', '',
        '本次只读取已完成的全域索引及其五次重复结果；未调用模型、未重算 Sobol 估计、未读取或恢复邻域分析。', '',
        'N65536 的未通过项共 **%d** 个；N131072 的未通过项为 **%d** 个。N131072 的 **%s/47,838** 个指标单元通过现有操作性收敛诊断。' %
        (len(data['failures']), int((data['frames'][131072].status_code != 1).sum()), format(int((data['frames'][131072].status_code == 1).sum()), ',')), '',
        '## 1. 同步与完整性', '',
        '- 从 HPC 同步两个分辨率的 %d 个导出文件，共 %.2f MiB；逐文件校验 SHA-256，传输前后远端清单一致。另校验并保存14参数上下界参考表。' % (receipt['n_files'], receipt['total_bytes'] / 1024 ** 2),
        '- 两档各有5次重复、201个固定氧点、两个输出、14个 S1、14个 ST 和91个不重复 S2 参数对，每档共47,838行。',
        '- 检查 gzip/NPZ 压缩容器、数组维度、氧点、参数集合、重复键、完成回执及汇总表校验值；从重复索引独立复算均值、重复范围、相邻分辨率差和诊断状态，与 NPZ/TSV 一致。',
        '- N131072 文件内的 N65536 五次重复索引与独立 N65536 文件逐元素完全相同。',
        '- 合并表时将 O₂ 标签规范到小数点后12位，仅消除文本浮点表示差异；不改变索引或统计计算。',
        '- 所有原始 PDF/PNG 进行了结构/解码检查；可视检查覆盖范围另见 `visual_review.md`。', '',
        '## 2. 判据与比较范围', '',
        '逐个“输出 × 指标/参数对 × O₂”单元检查：五次重复的最大值减最小值 ≤0.05，并且相邻分辨率五次重复均值的绝对差 ≤0.02；非有限值记为 undefined，不计通过。', '',
        '- N65536：相邻分辨率为 N32768→N65536。',
        '- N131072：相邻分辨率为 N65536→N131072。',
        '- ST：全部14个参数，每个输出2,814个单元。',
        '- 关键 S2：涉及 `p_mis_base`、`mu_hp` 或 `gamma_mu` 的36个不重复参数对，每个输出7,236个单元；三个焦点参数互相组成的参数对只算一次。',
        '- 氧区通过率统计使用互斥分组：低氧0–1%（含1%），中氧1–3%（不含两端），高氧3–5%（含3%）。带平均值采用区间内线性梯形积分除以区间宽度；不代表对氧浓度分布取期望。', '',
        '## 3. ST 与36对关键 S2 的改善', '']
    rows = []
    for output in OUTPUTS:
        for kind in ('ST', 'S2_focal'):
            a = summary[(summary.N == 65536) & (summary.output == output) & (summary.index_group == kind)].iloc[0]
            b = summary[(summary.N == 131072) & (summary.output == output) & (summary.index_group == kind)].iloc[0]
            rows.append(['倍性' if output == OUTPUTS[0] else '增长率', 'ST' if kind == 'ST' else '36对关键S2',
                '%d/%d（%.2f%%）' % (a.n_passed, a.n_total, a.pass_fraction * 100),
                '%d/%d（%.2f%%）' % (b.n_passed, b.n_total, b.pass_fraction * 100),
                '%.6f → %.6f' % (a.repeat_range_max, b.repeat_range_max), '%.6f → %.6f' % (a.delta_max, b.delta_max)])
    lines += [markdown_table(['输出', '指标', 'N65536通过数', 'N131072通过数', '最大重复范围', '最大相邻N均值差'], rows), '',
        '通过率以指标单元为分母，不是重复次数或模型分片数量。不同 N 的“相邻差”比较对象不同，见上节。所有均值差采用绝对量，避免对接近零的 S2 计算误导性的相对百分比。', '',
        '表中四组重复范围的中位数、95分位数和最大值在本轮均下降；单个参数对或氧点的波动不一定随 N 单调减小。', '']
    lines += ['### 此前未稳定的 ST', '']
    for (output, idx, pi, pj), sub in data['failures'].groupby(['output', 'index', 'parameter_i', 'parameter_j'], sort=False):
        m = data['metrics'].index((idx, pi, pj));k = OUTPUTS.index(output)
        old_range = data['summaries'][65536]['spread'][:, k, m]
        new_range = data['summaries'][131072]['spread'][:, k, m]
        lines += ['- **%s / `%s` / %s**：原有%d个未通过氧点，位于%s。N65536 的重复范围最大为 %.6f，N131072 降至 %.6f；这%d项在 N131072 均已通过。' %
            ('倍性' if output == OUTPUTS[0] else '增长率', pi, idx, len(sub), intervals(sub.O2_pct), old_range.max(), new_range.max(), int(sub.newly_passed.sum()))]
    failed655 = data['frames'][65536][data['frames'][65536].status_code != 1]
    lines += ['', '- 原有未通过项中，重复范围超门槛的有%d项，相邻分辨率均值差超门槛的有%d项。两种原因可同时出现，不能将两者简单相加。' %
        (int((failed655.repeat_range > RANGE_LIMIT).sum()), int((failed655.resolution_delta > DELTA_LIMIT).sum())), '',
        '### 关键 S2 的稳定程度', '']
    for output in OUTPUTS:
        a = summary[(summary.N == 65536) & (summary.output == output) & (summary.index_group == 'S2_focal')].iloc[0]
        b = summary[(summary.N == 131072) & (summary.output == output) & (summary.index_group == 'S2_focal')].iloc[0]
        sub = detail[(detail.output == output) & selection(detail, 'S2_focal')]
        rank = data['ranking'][(data['ranking'].output == output) & (data['ranking'].index_group == 'S2_focal')]
        lines += ['- **%s**：关键 S2 通过率 %.2f%% → %.2f%%；重复范围中位数 %.6f → %.6f，95分位数 %.6f → %.6f。两档均值的最大绝对变化为 %.6f；36对参数的逐氧点排序 Spearman 相关中位数为 %.4f、最小值为 %.4f。' %
            ('倍性' if output == OUTPUTS[0] else '增长率', a.pass_fraction * 100, b.pass_fraction * 100,
             a.repeat_range_median, b.repeat_range_median, a.repeat_range_p95, b.repeat_range_p95,
             sub.absolute_mean_change.max(), rank.spearman_rank_correlation.median(), rank.spearman_rank_correlation.min())]
    lines += ['', '关键 S2 在 N65536 已通过；本轮进一步量化其在更高采样量下的重复间波动和分辨率变化。接近零的参数对即使数值很稳定，其细微名次仍可能变动；排序指标不能替代效应大小。', '',
        '### N131072 中较强的关键参数对', '',
        '下表分别按各输出在201个氧点上的最大五次重复均值 S2 排序，列出前6对；完整36对均保留在源表和附图中，不对其设置效应量筛选门槛。', '']
    for output in OUTPUTS:
        top = data['metric_summary'][(data['metric_summary'].N == 131072) & (data['metric_summary'].output == output) & (data['metric_summary']['index'] == 'S2')]
        top = top.sort_values(['peak_value', 'parameter_i', 'parameter_j'], ascending=[False, True, True]).head(6)
        lines += ['**%s**' % ('倍性' if output == OUTPUTS[0] else '增长率'), '',
            markdown_table(['参数对', '峰值S2', '峰值O₂（%）', '最大重复范围'],
                [[r.parameter_i + ' × ' + r.parameter_j, '%.5f' % r.peak_value, '%.3f' % r.peak_O2_pct, '%.5f' % r.repeat_range_max]
                 for r in top.itertuples()]), '']
    lines += ['## 4. 图怎么读', '',
        '- [收敛总览](figures/convergence_overview.pdf)：左列比较历史通过率，中列比较重复范围的最大值和95分位数，右列比较相邻 N 均值绝对差；水平虚线为原判据。',
        '- [原未稳定ST](figures/previously_unstable_ST.pdf)：依次展示均值曲线及五次重复范围、重复范围、相邻 N 差。叉号标出 N65536 原未通过的氧点。',
        '- [较强关键S2](figures/focal_S2_top_pairs.pdf)：蓝灰虚线为 N65536，紫色实线为 N131072，阴影是各自五次重复的最小值到最大值，不是95%置信区间。',
        '- [全部14参数ST](figures/ST_all_parameters.pdf)及[全部36对关键S2](figures/S2_all_36_focal_pairs.pdf)：同一面板的两个 N 共用纵轴，不同面板可用不同纵轴；两类输出分别比较。',
        '- O₂保留零点：0–0.025%采用线性段，更高氧采用对数段；1%和3%为低、中、高氧分界。负的 S2 估计值保留，不解释为负交互。', '',
        '## 5. 结论的边界', '',
        'N131072 在既定阈值下通过全体单元的诊断，支持使用这档结果描述当前全域设计内的敏感性大小与两两交互。它不是无限采样极限的误差证明，也不等于交互显著性检验；0.05与0.02是既定操作性门槛。', '',
        '五次重复的范围不是概率覆盖率已校准的置信区间。S2不提供交互方向、协同/拮抗或补偿机制；这些问题仍需后续响应面或条件分析。此处没有重新进行该类分析。', '',
        '## 6. 文件与复现', '',
        '- 原始导出：`../N65536/`、`../N131072/`，包括原有图、源数据、行序/注释和 `resolution_comparison/`。',
        '- 同步与校验回执：`../%s/collection_receipt.json`、`verified_files.tsv`。' % collection.name,
        '- 独立数据验证：`validation_data.json`；PDF/PNG结构检查：`validation_assets.json`。',
        '- 统计汇总：`diagnostic_summary.tsv`、`diagnostic_summary_by_O2_band.tsv`。',
        '- 完整逐点比较：`all_index_comparison.tsv.gz`；重点比较：`ST_and_36_S2_comparison.tsv.gz`。',
        '- 五次重复源值：`five_repeat_source.tsv.gz`；曲线汇总源值：`curves_source.tsv.gz`。',
        '- 55项逐项追踪：`previously_not_passed_tracking.tsv`；参数/参数对汇总：`parameter_and_pair_summary.tsv`。',
        '- 历史通过率：`historical_ST_and_focal_S2_pass_rates.tsv`；排序稳定性：`ranking_stability.tsv`。',
        '- 复现方法：`README.md`；运行环境：`environment.json`；本次全部衍生文件哈希：`artifact_manifest.json`。', '']
    (out / 'interpretation_zh.md').write_text('\n'.join(lines))


def finalize(root, collection, out):
    distributions = ['numpy', 'pandas', 'matplotlib', 'Pillow', 'pypdf', 'contourpy', 'cycler', 'fonttools', 'kiwisolver', 'pyparsing', 'packaging', 'python-dateutil']
    versions = {name: importlib.metadata.version(name) for name in distributions}
    save(out / 'environment.json', dict(python=sys.version, executable=sys.executable, packages=versions,
         comparison_script=str(Path(__file__).resolve()), comparison_script_sha256=sha(__file__),
         collector_script_sha256=sha(Path(__file__).with_name('collect_global_convergence_exports.py')),
         completed_utc=now(), model_calls=0))
    (out / 'requirements.txt').write_text(''.join('%s==%s\n' % pair for pair in versions.items()))
    # Commands below use repository-relative paths, not a machine-specific env.
    script = 'oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/compare_global_convergence_exports.py'
    relroot = 'oxygen/results/eFAST/sobol/global_convergence_20261021'
    text = '''# Local global Sobol convergence comparison

This directory contains derivative analysis of completed global indices only.
No model evaluations or sensitivity re-estimation are performed. Neighborhoods
remain archived. The original N65536/N131072 exports are preserved unchanged.

## Reproduce from the repository root

Use Python 3.12 with the versions in `requirements.txt` (the pinned local plotting
environment, not a change to the HPC SIF):

```sh
python3 -m pip install -r {root}/comparison_N65536_N131072/requirements.txt
python3 {script} --root {root} --collection {root}/{collection}
```

Add `--render-pdfs` to render every original/comparison PDF with Poppler and
regenerate the validation contact sheets; this requires `pdftoppm` on PATH.
Manual visual-review notes are separate from the automated rendering record.

For transfer, use `collect_global_convergence_exports.py` in the same code folder
with `--remote-root` and `--local-root`. It selects only completed exports and
small provenance files, verifies source/destination hashes, and preserves
different existing files by refusing to overwrite them.

```sh
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/collect_global_convergence_exports.py --remote-root /share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/global_convergence_20261021 --local-root {root}
```

The exact collector executed for this download is archived as
`../{collection}/provenance/collector_used.py`; its SHA-256 matches the collection
receipt. The current collector also supports filling missing exports in a
partially populated result directory while preserving existing identical files.

## Scientific definitions

- Two output-specific analyses; all 14 ST factors and 36 distinct S2 pairs that
  involve p_mis_base, mu_hp, or gamma_mu.
- Five independent scramble repeats, 201 oxygen points, unchanged global bounds.
- Repeat range <=0.05 and adjacent-resolution mean absolute change <=0.02.
- N65536 uses N32768 as its previous resolution; N131072 uses N65536.
- Repeat min-max bands are not confidence intervals; S2 sign is not a direction
  or compensation measure. Negative estimates are preserved as sampling noise.
- The original exported arrays and TSVs are checked independently before plotting.
- Oxygen keys are rounded to 12 decimals only for joins; index values and
  statistical calculations are not rounded. Original files stay byte-identical.
- Whole-curve ranking correlations include near-zero pairs; interpret them with
  absolute effect sizes. Top-six display selection uses maximum mean S2 across
  oxygen for each output, with all 36 pairs retained in the appendix and sources.

See `interpretation_zh.md`, `validation_data.json`, `validation_assets.json`,
`visual_review.md`, and the source/manifest files for evidence and coverage.
'''.format(root=relroot, script=script, collection=collection.name)
    (out / 'README.md').write_text(text)
    records = [dict(path=str(p.relative_to(out)), size=p.stat().st_size, sha256=sha(p))
               for p in sorted(out.rglob('*')) if p.is_file() and p.name != 'artifact_manifest.json']
    save(out / 'artifact_manifest.json', dict(created_utc=now(), files=records,
         input_manifest=str(collection / 'remote_manifest.json'), input_manifest_sha256=sha(collection / 'remote_manifest.json'),
         script_sha256=sha(__file__), model_calls=0, neighborhood_files_read=0))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True)
    parser.add_argument('--collection', required=True)
    parser.add_argument('--output')
    parser.add_argument('--validate-only', action='store_true')
    parser.add_argument('--render-pdfs', action='store_true')
    args = parser.parse_args()
    root, collection = Path(args.root).resolve(), Path(args.collection).resolve()
    out = Path(args.output).resolve() if args.output else root / 'comparison_N65536_N131072'
    out.mkdir(parents=True, exist_ok=True)
    data = validate(root, collection, out)
    analysis(root, out, data)
    if not args.validate_only:
        plot_all(out, data)
        validate_assets(root, out)
        if args.render_pdfs:
            render_pdf_checks(root, out)
        write_report(root, collection, out, data)
        finalize(root, collection, out)
    print(json.dumps(dict(output=str(out), validation='passed',
        failures_N65536=len(data['failures']), failures_N131072=int((data['frames'][131072].status_code != 1).sum()),
        model_calls=0, figures_generated=not args.validate_only)))


if __name__ == '__main__':
    main()
