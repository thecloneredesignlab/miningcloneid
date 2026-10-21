"""Global-only plots and convergence comparisons from certified index arrays."""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm, to_rgb
from matplotlib.patches import Patch

COLORS = {'High O2': '#B2182B', 'Low O2': '#2166AC', 'No resolved O2 preference': '#8C8C8C',
          'Insensitive': '#8C8C8C', 'Not sensitive': '#8C8C8C'}
FOCAL = {'p_mis_base', 'mu_hp', 'gamma_mu'}


def weights(oxygen, lo, hi):
    ids = np.flatnonzero((oxygen >= lo) & (oxygen <= hi))
    w = np.zeros(len(oxygen))
    dx = np.diff(oxygen[ids]) / (2 * (hi - lo))
    w[ids[:-1]] += dx
    w[ids[1:]] += dx
    return w


def bh(p):
    p = np.asarray(p)
    ix = np.argsort(p)
    q = np.empty_like(p)
    q[ix] = np.minimum(1, np.minimum.accumulate((p[ix] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return q


def classify(values, metrics, oxygen):
    """Resample five whole repeat curves; linear integration commutes with mean."""
    loww, highw = weights(oxygen, 0, 1), weights(oxygen, 3, 5)
    low = np.einsum('rom,o->rm', values[:, :, 0, :], loww)
    high = np.einsum('rom,o->rm', values[:, :, 0, :], highw)
    draws = np.random.default_rng(5826).integers(0, 5, size=(5000, 5))
    boot = np.mean((low - high)[draws], axis=1)
    p = np.minimum(1, 2 * (np.minimum((boot <= 0).sum(axis=0), (boot >= 0).sum(axis=0)) + 1) / 5001)
    q = np.ones(len(metrics))
    for kind in ['S1', 'ST', 'S2']:
        ids = [i for i, m in enumerate(metrics) if m[0] == kind]
        q[ids] = bh(p[ids])
    rows = []
    for i, (kind, pi, pj) in enumerate(metrics):
        l, h = float(low[:, i].mean()), float(high[:, i].mean())
        group = ('High O2' if h > l else 'Low O2') if max(l, h) > 0 and q[i] < .05 and l != h else 'No resolved O2 preference'
        rows.append(dict(index=kind, parameter_i=pi, parameter_j=pj, low_score=l, high_score=h, sort_side='High' if h >= l else 'Low',
            sensitivity_group=group, bh_q=float(q[i]), bootstrap_sign_p=float(p[i]),
            bootstrap_delta_q025=float(np.quantile(boot[:, i], .025)), bootstrap_delta_q975=float(np.quantile(boot[:, i], .975)),
            classification_output='dominant_mean_ploidy', bootstrap_unit='whole scramble-repeat curve', bootstrap_reps=5000,
            peak_threshold='none', interpretation='exploratory O2 preference; not numerical convergence or calibrated significance'))
    return rows


def row_order(ids, classes):
    return sorted(ids, key=lambda i: (0, -classes[i]['high_score'], i) if classes[i]['high_score'] >= classes[i]['low_score']
                  else (1, classes[i]['low_score'], i))


def oxygen_axis(ax):
    ax.set_xscale('symlog', base=10, linthresh=.025, linscale=1)
    ax.set_xlim(0, 5)
    ax.set_xticks([0, .025, .1, .5, 1, 3, 5], ['0', '.025', '.1', '.5', '1', '3', '5'])
    ax.set_xlabel('Fixed oxygen (%)', fontsize=8)
    for x in [1, 3]:
        ax.axvline(x, ls='--', lw=.6, color='black')
    transform = ax.xaxis.get_transform()
    for lo, hi, label in [(0, 1, 'Low'), (1, 3, 'Mid'), (3, 5, 'High')]:
        mid = transform.inverted().transform(np.mean(transform.transform(np.array([lo, hi]))))
        ax.text(mid, 1.01, label, transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=7)


def savefig(fig, stem):
    fig.canvas.draw()
    fig.savefig(str(stem) + '.pdf', bbox_inches='tight')
    fig.savefig(str(stem) + '.png', dpi=140, bbox_inches='tight')
    plt.close(fig)


def panel(fig, cell, data, ids, kind, metrics, classes, oxygen, process, correlation, palette, title):
    # Give the colorbar its own row; stealing space from nested axes can place
    # it over the next panel or over the last S2 rows in tall combined figures.
    grid = cell.subgridspec(3, 2, width_ratios=[.16, 1], height_ratios=[1, .15, .04], wspace=.025, hspace=.15)
    ann, ax = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    cax = fig.add_subplot(grid[2, 1])
    rgb = []
    for i in ids:
        _, pi, pj = metrics[i]
        members = [pi, pj] if pj else [pi, pi]
        colors = [to_rgb(palette[process[p]]) for p in members] + [to_rgb(COLORS.get(correlation[p], '#8C8C8C')) for p in members]
        if kind != 'S2':
            colors += [to_rgb(COLORS[classes[i]['sensitivity_group']])] * 2
        rgb.append(colors)
    ann.imshow(np.asarray(rgb), aspect='auto', interpolation='nearest', extent=(-.5, len(rgb[0]) - .5, len(ids) - .5, -.5))
    ann.set_yticks(range(len(ids)), [metrics[i][1] + (' × ' + metrics[i][2] if metrics[i][2] else '') for i in ids], fontsize=6)
    ann.set_xticks([.5, 2.5] + ([4.5] if kind != 'S2' else []), ['P', 'C'] + (['S'] if kind != 'S2' else []))
    ann.xaxis.tick_top()
    ann.tick_params(axis='both', length=0)
    edges = np.r_[0, (oxygen[:-1] + oxygen[1:]) / 2, 5]
    matrix = data[:, ids].T
    finite = matrix[np.isfinite(matrix)]
    vmax = max(float(np.max(np.abs(finite))) if kind == 'S2' else float(np.max(finite)), 1e-10)
    if np.any(finite < 0):
        cmap = LinearSegmentedColormap.from_list('signed', ['#2166ac', 'white', '#54278f'])
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    else:
        cmap = LinearSegmentedColormap.from_list('positive', ['white', '#54278f'])
        norm = Normalize(0, vmax)
    im = ax.pcolormesh(edges, np.arange(len(ids) + 1) - .5, np.ma.masked_invalid(matrix), shading='flat', cmap=cmap, norm=norm, rasterized=True)
    ax.set_ylim(len(ids) - .5, -.5)
    ax.set_yticks([])
    oxygen_axis(ax)
    ax.set_title(title, pad=24, fontsize=10)
    fig.colorbar(im, cax=cax, orientation='horizontal', label=kind)


def render(root, n, p, c, f, r, values, mean, status):
    root = Path(root)
    dest = root / ('N%d' % n)
    figures = dest / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    oxygen = np.asarray(c.OXYGEN)
    classes = classify(values[-1], f.METRICS, oxygen)
    r.write_table(dest / 'classification.tsv', classes, list(classes[0]))
    process = {x['parameter']: x['parameter_group'] for x in r.read_table(root / 'inputs/figure4_parameter_groups_source.tsv')}
    correlation = {x['parameter']: x['o2_association_group'] for x in r.read_table(root / 'inputs/figure4b_o2_classification_source.tsv')}
    palette_rows = r.read_table(root / 'inputs/figure4_parameter_group_palette_source.tsv')
    palette = {x['parameter_group']: x['color'] for x in palette_rows}
    kinds = ['S1', 'ST', 'S2', 'S2_focal']
    chosen = {kind: row_order([i for i, m in enumerate(f.METRICS) if m[0] == kind.split('_')[0] and
              (kind != 'S2_focal' or bool(set(m[1:]) & FOCAL))], classes) for kind in kinds}
    orderrows = [dict(panel=kind, rank=rank + 1, **classes[i]) for kind in kinds for rank, i in enumerate(chosen[kind])]
    r.write_table(dest / 'row_order.tsv', orderrows, list(orderrows[0]))
    names = {'S1': 's1_parameter_by_oxygen', 'ST': 'st_parameter_by_oxygen', 'S2': 's2_pair_by_oxygen', 'S2_focal': 's2_focal_pairs_by_oxygen'}
    for kind in kinds:
        fig = plt.figure(figsize=(19, 25 if kind == 'S2' else 13 if kind == 'S2_focal' else 7))
        grid = fig.add_gridspec(1, 2, wspace=.65)
        for k, label in enumerate(['Dominant mean ploidy', 'Net live-cell growth rate']):
            panel(fig, grid[k], mean[:, k], chosen[kind], kind.split('_')[0], f.METRICS, classes, oxygen, process, correlation, palette, label)
        fig.suptitle('Global N%d; five repeats; %s' % (n, kind), y=.99)
        fig.text(.5, .015, 'P: process; C: original Figure 4B correlation; S: ploidy Sobol O2 preference (exploratory).\n'
                 'Rows: High >= Low, High descending; then Low > High, Low ascending. Negative estimates are sampling noise.\n'
                 'O2: log above 0.025%, linear below to retain zero; dashed boundaries at 1% and 3%.', ha='center', fontsize=8)
        fig.subplots_adjust(bottom=.15 if kind in ['S1', 'ST'] else .08, top=.9)
        savefig(fig, figures / names[kind])
    fig = plt.figure(figsize=(20, 22))
    grid = fig.add_gridspec(3, 2, wspace=.65, hspace=.35, height_ratios=[14, 14, 36])
    for row, kind in enumerate(['S1', 'ST', 'S2_focal']):
        for k, label in enumerate(['Ploidy', 'Growth']):
            panel(fig, grid[row, k], mean[:, k], chosen[kind], kind.split('_')[0], f.METRICS, classes, oxygen, process, correlation, palette, kind + ': ' + label)
    fig.suptitle('Global N%d; S1, ST and 36 focal S2 pairs; five repeats' % n)
    fig.subplots_adjust(top=.94, bottom=.1)
    savefig(fig, figures / 'global_s1_st_focal_s2')
    fig = plt.figure(figsize=(14, 3))
    handles = [Patch(facecolor=x['color'], label=x['display_label']) for x in palette_rows]
    handles += [Patch(facecolor=COLORS[x], label=x) for x in ['High O2', 'Low O2', 'No resolved O2 preference']]
    fig.legend(handles=handles, loc='center', ncol=3, frameon=False)
    fig.text(.5, .05, 'S annotation is exploratory, uses ploidy repeat curves, and does not determine row order or convergence.', ha='center', fontsize=9)
    savefig(fig, figures / 'annotation_legend')
    compare(root, n, p, c, f, r, mean, status, values)


def compare(root, n, p, c, f, r, mean, status, values):
    series = []
    with np.load(Path(p['source_initial_root']) / 'global/full/global/indices_summary.npz') as z:
        series.append((4096, z['mean'], z['status']))
    for size in [8192, 16384, 32768] + ([65536] if n > 65536 else []):
        path = (Path(p['source_global_root']) / ('N%d/cases/global/indices_summary.npz' % size)) if size < 65536 else root / 'N65536/indices_summary.npz'
        with np.load(path) as z:
            series.append((size, z['mean'], z['status']))
    series.append((n, mean, status))
    directory = root / ('N%d/resolution_comparison' % n)
    directory.mkdir(parents=True, exist_ok=True)
    rows, transitions = [], []
    baseline = series[0][2]
    for size, data, st in series:
        for k, out in enumerate(c.OUTPUTS):
            for kind in ['S1', 'ST', 'S2']:
                ids = [i for i, m in enumerate(f.METRICS) if m[0] == kind]
                selected = st[:, k, ids]
                rows.append(dict(N=size, output=out, index=kind, n_total=selected.size, n_passed=int((selected == 1).sum()),
                    n_undefined=int((selected == 0).sum()), pass_fraction=float((selected == 1).mean())))
                if size == n:
                    oldst = baseline[:, k, ids]
                    transitions.append(dict(N=n, output=out, index=kind, newly_passed=int(((oldst != 1) & (selected == 1)).sum()),
                        lost_pass=int(((oldst == 1) & (selected != 1)).sum()), still_not_passed=int(((oldst != 1) & (selected != 1)).sum())))
    r.write_table(directory / 'diagnostic_summary.tsv', rows, list(rows[0]))
    r.write_table(directory / 'diagnostic_transitions_from_N4096.tsv', transitions, list(transitions[0]))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
    for k, out in enumerate(c.OUTPUTS):
        for kind in ['S1', 'ST', 'S2']:
            rr = [x for x in rows if x['output'] == out and x['index'] == kind]
            axes[k].plot([x['N'] for x in rr], [x['pass_fraction'] for x in rr], marker='o', label=kind)
        axes[k].set(xscale='log', ylim=(0, 1.02), xlabel='Base sample count N', ylabel='Diagnostic pass fraction', title=out)
        axes[k].legend()
    savefig(fig, directory / 'diagnostic_pass_fraction')
    records = []
    for kind in ['ST', 'S2']:
        ids = [i for i, m in enumerate(f.METRICS) if m[0] == kind and (kind == 'ST' or bool(set(m[1:]) & FOCAL))]
        fig, axes = plt.subplots(len(ids), 2, figsize=(14, len(ids) * 2.3), squeeze=False, layout='constrained')
        for row, m in enumerate(ids):
            for k in range(2):
                ax = axes[row, k]
                for size, data, st in series:
                    ax.plot(c.OXYGEN, data[:, k, m], label='N%d' % size, lw=.8)
                    for oi, oxygen in enumerate(c.OXYGEN):
                        records.append(dict(N=size, output=c.OUTPUTS[k], index=kind, parameter_i=f.METRICS[m][1], parameter_j=f.METRICS[m][2],
                            O2_pct=oxygen, value_mean=float(data[oi, k, m]), diagnostic_status=int(st[oi, k, m])))
                ax.fill_between(c.OXYGEN, values[-1, :, :, k, m].min(axis=0), values[-1, :, :, k, m].max(axis=0), alpha=.13, color='black', label='Current five-repeat range')
                oxygen_axis(ax)
                ax.set_title(' × '.join(x for x in f.METRICS[m][1:] if x) + '; ' + c.OUTPUTS[k], fontsize=8, pad=24)
                if row == 0: ax.legend(fontsize=6)
        savefig(fig, directory / (kind.lower() + '_curves'))
    r.write_table(directory / 'curves_source.tsv.gz', records, list(records[0]))
