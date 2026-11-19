"""Independent reverse chain-rule check of the saved real-kernel minor."""
import csv
import sys
from pathlib import Path

from certify import (arb, arb_mat, ctx, exact, positive, anchor, operators,
                     COLS, GRID, N, D, DT, CAPACITY, FALLBACK, BRIDGE,
                     CONTRACT, STEP1, count_or_none, read, table, model)
from algebra import decode, independent_minor_check, nonzero
from collections import defaultdict

OUT = Path(__file__).parent
BASE = OUT / 'anchor_384'
I = len(model.MECHANISMS)


def dot(x, y):
    return (x.transpose() * y)[0, 0]


def thaw(rows, columns):
    source = read(BASE / 'real_kernel_19_matrix.tsv')
    result = arb_mat(len(rows), len(columns))
    for r in source:
        result[int(r['row']), int(r['column'])] = decode(r)
    return result


def precompute(out):
    ctx.prec = 512
    p, _ = anchor('anchor')
    schedule = read(CONTRACT / 'experimental_schedule.tsv')
    branches = {(r['cohort'], r['segment']): r for r in read(BASE / 'passage_branches.tsv')}
    slots = {}
    for r in schedule:
        slots.update({pid: (r['cohort'], r['segment_id']) for pid in r['data_ids'].split(';')})
    operators_by_oxygen = operators(p, sorted({r['oxygen_pct'] for r in schedule}), out)
    dt = exact(DT)
    guide = {}
    groups = {}
    blocks = {}
    with Path(sys.argv[1]).open(newline='') as handle:
        native = {(r[0],r[1],int(r[2])): list(map(float, r[3:])) for r in list(csv.reader(handle, delimiter='\t'))[1:]}
    maximum = 0.
    array_count = 0
    for r in schedule:
        key = r['cohort'], r['segment_id']
        branch = branches[key]
        n = count_or_none(r['initial_cells'])
        if n is None:
            n = arb(FALLBACK) / CAPACITY
        parent_key = r['cohort'], r['parent_segment_id']
        if branch['parent_day'] == 'NA':
            mean, sd = (p['init_' + field + '_' + key[0]] for field in ('mean','sd'))
            # Different normalization expression from algebra.initial.
            mass = [(-((j-mean)/sd)**2 / 2).exp() for j in GRID]
            x = arb_mat([[n*z/sum(mass)] for z in mass])
        else:
            ancestor = groups[parent_key][int(exact(branch['parent_day']).mid().fmpq()) * 20]
            x = n / sum(ancestor.entries()) * ancestor
        g, dg, death, ddeath = operators_by_oxygen[r['oxygen_pct']]
        if r['oxygen_pct'] not in guide:
            # Matrix orientation is rebuilt separately for reverse parameter multiplication.
            dgt = [arb_mat([[dg[i*I+j,k] for k in range(N)] for i in range(N)]) for j in range(I)]
            guide[r['oxygen_pct']] = (g.transpose(), arb_mat([[z] for z in death]), dgt,
                                      arb_mat([[ddeath[j][i] for j in range(I)] for i in range(N)]))
        blocks[key] = (r['oxygen_pct'], n)
        steps = int(r['nominal_duration_days']) * 20
        hist = []
        for step in range(steps+1):
            hist.append(x)
            if step % 20 == 0:
                reference = native[key + (step//20,)]
                err = sum(abs(float(x[i,0].mid()) * CAPACITY - reference[i]) for i in range(N))
                maximum = max(maximum, err / max(1., sum(reference)))
                array_count += 1
            if step < steps:
                div = g*x * (1 - sum(x.entries()))
                x = x + dt * (div - arb_mat([[death[i] * x[i,0]] for i in range(N)]))
        groups[key] = hist
    assert array_count == len(native) and maximum < 5e-11, (array_count, len(native), maximum)
    return p, groups, guide, blocks, branches, slots, array_count, maximum


def reverse(label, p, groups, guide, blocks, branches, slots):
    key = slots[label['passage']]
    dt = exact(DT)
    step = int(label['local_day']) * 20
    x = groups[key][step]
    total = sum(x.entries())
    width_direct = arb(0)
    repeat = [arb(0) for _ in COLS]
    if label['stream'] == 'growth':
        lam = arb_mat([[1/(total*int(label['local_day']))] for _ in GRID])
    else:
        bin_id = int(label['chromosome'])
        width = p['sigma_kary']
        k_row, width_d = [], []
        for n in GRID:
            square = [(j-n)**2 for j in GRID]
            weights = [(-arb(z)/(2*width**2)).exp() for z in square]
            denom = sum(weights)
            k = weights[bin_id-GRID[0]] / denom
            weighted_square = sum(w*z for w,z in zip(weights,square)) / denom
            k_row.append(k)
            width_d.append(k*((bin_id-n)**2-weighted_square)/width**3)
        k_row = arb_mat([[v] for v in k_row])
        value = dot(k_row, x) / total
        width_direct = dot(arb_mat([[v] for v in width_d]), x) / total
        lam = arb_mat([[(v - value)/total] for v in k_row.entries()])
    while True:
        oxygen, n = blocks[key]
        gt, death, dg, ddeath = guide[oxygen]
        g = gt.transpose()
        for substep in range(step-1, -1, -1):
            x = groups[key][substep]
            crowd = 1 - sum(x.entries())
            lamx = arb_mat([[lam[i,0]*x[i,0]] for i in range(N)])
            death_rate_gradient = lamx.transpose() * ddeath
            for j in range(I):
                repeat[j] += dt*(crowd*dot(lam, dg[j]*x) - death_rate_gradient[0,j])
            # VJP from A(x)=x+dt[(1-1^Tx)Gx-Dx].
            density = dot(lam, g*x)
            lam = lam + dt*(gt*lam*crowd -
                            arb_mat([[density + death[i,0]*lam[i,0]] for i in range(N)]))
        x = groups[key][0]
        parent_id = branches[key]['parent_segment']
        if parent_id == 'NA':
            cohort = key[0]
            mean = p['init_mean_'+cohort]
            sd = p['init_sd_'+cohort]
            pl = x/n
            for name, scores in [('init_mean_', [(j-mean)/sd**2 for j in GRID]),
                                 ('init_sd_', [(j-mean)**2/sd**3 for j in GRID])]:
                center = sum(pl[i,0]*s for i,s in enumerate(scores))
                derivative = arb_mat([[x[i,0]*(s-center)] for i,s in enumerate(scores)])
                repeat[COLS.index(name+cohort)] += dot(lam, derivative)
            break
        parent_step = int(branches[key]['parent_day']) * 20
        key = key[0], parent_id
        step = parent_step
        parent_x = groups[key][step]
        parent_total = sum(parent_x.entries())
        expected = dot(lam, parent_x) / parent_total
        lam = (n / parent_total) * arb_mat([[value-expected] for value in lam.entries()])
    return repeat + [width_direct]


def main():
    out = OUT / 'independent_review'
    out.mkdir(exist_ok=True)
    rows = read(BASE / 'real_kernel_19_rows.tsv')
    cols = read(BASE / 'real_kernel_19_cols.tsv')
    p, groups, guide, blocks, branches, slots, array_count, native_gap = precompute(out)
    original = thaw(rows, cols)
    reverse_rows = []
    for i, label in enumerate(rows):
        x = reverse(label, p, groups, guide, blocks, branches, slots)
        scaled = [z * arb((1,int(c['rescale_power_of_2']))) for z,c in zip(x,cols)]
        for j, ball in enumerate(scaled):
            assert original[i,j].overlaps(ball), (i,j,original[i,j],ball)
        reverse_rows.append(scaled)
        print(f'Reverse certificate row {i+1}/{len(rows)}', flush=True)
    rebuilt = arb_mat(reverse_rows)
    determinant = rebuilt.det()
    verified, bound = independent_minor_check(rebuilt)
    assert nonzero(determinant) and verified
    from certificates import save_matrix
    save_matrix(out / 'reverse_real_kernel_19_matrix.tsv', rebuilt)
    table(out / 'review_summary.tsv', [dict(precision=ctx.prec, daily_native_state_arrays=array_count,
        max_daily_native_L1_error_over_live_total=native_gap, matrix_entries_overlapping=len(rows)*len(cols),
        reverse_determinant=str(determinant), reverse_neumann_bound=str(bound))])


if __name__ == '__main__':
    main()
