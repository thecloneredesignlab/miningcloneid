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


def panel(fig, cell, data, ids, kind, metrics, classes, oxygen, process, correlation, palette, title, vmax_override=None):
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
    if vmax_override is not None:
        vmax = max(float(vmax_override), 1e-10)
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
