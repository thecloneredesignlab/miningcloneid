"""Outward-ball forward sensitivities for the sparse production endpoint design."""
import argparse
import csv
import hashlib
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import flint
import numpy as np
from flint import arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / (
    'oxygen/code/O2_supply_demand_MAP/analysis/structural_identifiability/'
    'structural_first/culture_stageB_finite'
)))
import model
from algebra import encode, exact, initial, independent_minor_check, nonzero, rates
from certificates import save_matrix
from protocol import INITIALS, MECHANISMS, read, table

COLS = MECHANISMS + INITIALS
OBS = COLS + ('sigma_kary',)
GRID = model.GRID
N = len(GRID)
D = len(COLS)
DT = '1/20'
CAPACITY = 10**12
FALLBACK = 10**6
BRIDGE = ROOT / 'oxygen/results/structural_identifiability/20261116_culture_step05_bridge_design_01'
ARCHIVE = ROOT / 'oxygen/results/structural_identifiability/20261109_seed144_wp0/source/seed144'
CONTRACT = ROOT / 'oxygen/results/synthetic_identifiability/20261113_culture_contract_steps01_03'
STEP1 = ROOT / 'oxygen/results/structural_identifiability/20261117_culture_actual_schedule_step01'


def bounds(value):
    return arb(value.lower()), arb(value.upper())


def positive(value, label):
    if not value.is_finite() or not value > 0:
        raise ValueError(f'{label}: strict positivity unresolved: {value}')
    return arb(value.lower())


def anchor(case):
    point = {}
    metadata = []
    for row in read(Path(__file__).parent / 'anchor.tsv'):
        decimal_bounds = tuple(Fraction(row[key]) for key in ('lower', 'upper'))
        value = Fraction(float.fromhex(row['binary64']))
        boundary = ''
        for name, b in zip(('lower', 'upper'), decimal_bounds):
            # The native bound is the binary64 rounding of the declared decimal endpoint.
            if float(value) == float(b):
                boundary = name
                value = b
                if case == 'interior':
                    value += (decimal_bounds[1] - decimal_bounds[0]) * Fraction(
                        1 if name == 'lower' else -1, 10**8
                    )
                break
        lo, hi = decimal_bounds
        assert lo <= value <= hi, (row['parameter'], value, lo, hi)
        metadata.append(dict(parameter=row['parameter'], native_anchor=row['binary64'],
                             lower=str(lo), upper=str(hi), anchor_boundary=boundary,
                             exact_value=str(value)))
        point[row['parameter']] = exact(value)
    assert set(point) == set(OBS + ('sigma_growth',))
    return point, metadata


def fixed_rates(p, u, oxygen):
    if not oxygen.is_zero():
        return rates(p, u, oxygen)
    # For C,n>0 and O=0, H=1 identically. C and n have zero partials.
    inv_u = 1 / u
    ug, um = u ** p['gamma_growth'], u ** p['gamma_mu']
    exponent = inv_u ** p['buffer_n_exp']
    den = 1 + p['alpha_o2'] * ug
    lam = p['lam_max'] / den
    mu = p['mu_hp'] * um
    half = p['k_o_mis'] + mu
    prob = p['p_mis_base'] + p['p_misseg'] * mu / half
    log_survival = p['buffer_smax'].log() - p['buffer_beta'] * exponent
    dl = dict.fromkeys(MECHANISMS, arb(0))
    dm, ds = dl.copy(), dl.copy()
    dl.update(lam_max=1 / den, alpha_o2=-lam * ug / den,
              gamma_growth=-lam * p['alpha_o2'] * ug * u.log() / den)
    dm.update(mu_hp=um, gamma_mu=mu * u.log())
    dp = {name: p['p_misseg'] * p['k_o_mis'] * dm[name] / half**2 for name in MECHANISMS}
    dp.update(p_mis_base=arb(1), p_misseg=mu / half, k_o_mis=-p['p_misseg'] * mu / half**2)
    ds.update(buffer_smax=1 / p['buffer_smax'], buffer_beta=-exponent,
              buffer_n_exp=-p['buffer_beta'] * exponent * inv_u.log())
    return (lam, mu, prob, log_survival), (dl, dm, dp, ds)


def operators(p, oxygen_values, out):
    # Retain the previously source-reviewed 133-state event router and manual gradients.
    model.rates = fixed_rates
    result = {}
    branch = []
    masks = []
    for oxygen in oxygen_values:
        op = model.operator(p, exact(oxygen))
        minimum_gap = min(r['minimum_margin'] for r in op['masks'])
        survival_upper = max(v['survival'].upper() for v in op['values'])
        margin = min(positive(1 - exact(DT) * (v['lam'] + v['mu']), 'Euler diagonal')
                     for v in op['values'])
        positive(1 - survival_upper, 'survival clamp')
        branch.append(dict(oxygen_pct=oxygen, kept=sum(r['mask'].count('1') for r in op['masks']),
                           minimum_probability_margin=str(minimum_gap),
                           minimum_euler_diagonal_margin=str(margin),
                           maximum_per_copy_survival_upper=str(survival_upper),
                           mask_sha256=op['mask_sha256']))
        masks.extend(dict(oxygen_pct=oxygen, chr_count=r['N'], kept_mask=r['mask'],
                          margin=str(r['minimum_margin'])) for r in op['masks'])
        dG = op['dG']
        packed = arb_mat([[dG[j][i, k] for k in range(N)]
                          for i in range(N) for j in range(len(MECHANISMS))])
        result[oxygen] = (op['G'], packed, op['mu'], op['dmu'])
        print('operator', oxygen, 'mask margin', minimum_gap, flush=True)
    table(out / 'operator_branches.tsv', branch)
    table(out / 'probability_masks.tsv', masks)
    return result


def colsum(matrix, j=0):
    return sum(matrix[i, j] for i in range(N))


def normalized_reseed(x, jac, total, initial_size):
    scale = initial_size / total
    dtotal = [colsum(jac, j) / total for j in range(D)]
    return x * scale, arb_mat([[(jac[i, j] - x[i, 0] * dtotal[j]) * scale for j in range(D)]
                               for i in range(N)])


def propagate(x, jac, stop_step, need_derivatives, days, operator, key):
    G, dG, mu, dmu = operator
    dt = exact(DT)
    one_minus_death = [1 - dt * v for v in mu]
    dt_dmu = [[dt * column[i] for column in dmu] for i in range(N)]
    wanted = dict(zip((int(Fraction(day) / Fraction(DT)) for day in days), days))
    states = {}
    minimum = arb(1)
    for step in range(stop_step + 1):
        total = colsum(x)
        gap = positive(1 - total, f'{key} step {step} logistic branch')
        minimum = min(minimum, gap)
        if step in wanted:
            states[wanted[step]] = (x, jac if step <= need_derivatives else None, total)
        if step == stop_step:
            break
        crowd = dt * (1 - total)
        gx = G * x
        if step < need_derivatives:
            gj, q = G * jac, dG * x
            dtotal = [dt * colsum(jac, j) for j in range(D)]
            entries = []
            m = len(MECHANISMS)
            for i in range(N):
                xi = x[i, 0]
                row = [one_minus_death[i] * jac[i, j] + crowd * gj[i, j] - gx[i, 0] * dtotal[j]
                       for j in range(D)]
                for j in range(m):
                    row[j] += crowd * q[m * i + j, 0] - xi * dt_dmu[i][j]
                entries.append(row)
            jac = arb_mat(entries)
        x = arb_mat([[one_minus_death[i] * x[i, 0] + crowd * gx[i, 0]] for i in range(N)])
    return states, minimum


def select(states, target, required, verify_day, label):
    days = sorted(states, key=Fraction)
    days = [day for day in days if Fraction(day) > 0] or days
    margins = []
    if target is None:
        score_day, method = days[-1], 'last_without_target'
    else:
        errors = {day: abs(states[day][2] - target) for day in days}
        score_day = min(days, key=lambda day: errors[day].mid())
        margins.extend(positive(errors[day] - errors[score_day], f'{label} nearest gap {day}')
                       for day in days if day != score_day)
        method = 'nearest_target'
    chosen = score_day
    if required is not None:
        gap = states[score_day][2] - required
        if gap > 0:
            margins.append(positive(gap, f'{label} sufficient target-day supply'))
            method = 'target_day_supply'
        elif gap < 0:
            margins.append(positive(-gap, f'{label} insufficient target-day supply'))
            eligible = []
            for day in days:
                diff = states[day][2] - required
                if diff > 0:
                    eligible.append(day)
                margins.append(positive(diff if diff > 0 else -diff, f'{label} eligibility {day}'))
            if not eligible:
                raise ValueError(f'{label} no eligible parent state')
            chosen = min(eligible, key=lambda day: states[day][2].mid())
            margins.extend(positive(states[day][2] - states[chosen][2], f'{label} fallback gap {day}')
                           for day in eligible if day != chosen)
            method = 'eligible_fallback'
        else:
            raise ValueError(f'{label} supply equality unresolved')
    if chosen != verify_day:
        raise ValueError(f'{label} expected {verify_day} but certified {chosen}')
    return chosen, method, str(min(margins)) if margins else 'fixed_by_design'


def count_or_none(value):
    return None if value in ('NA', '') or Fraction(value) <= 0 else exact(value) / CAPACITY


def growth_rows(jac, total, known_initial, day):
    # R's 1e-12 count floor is inactive when both counts are greater than it.
    positive(total * CAPACITY - exact('1e-12'), 'growth endpoint count floor')
    positive(known_initial * CAPACITY - exact('1e-12'), 'growth initial count floor')
    return [colsum(jac, j) / (total * exact(day)) for j in range(D)] + [arb(0)]


def kary_rows(x, jac, total, kernel):
    smoothed, derivative = kernel
    frac = x / total
    dtotal = [colsum(jac, j) for j in range(D)]
    df = arb_mat([[(jac[i, j] - frac[i, 0] * dtotal[j]) / total for j in range(D)] for i in range(N)])
    all_j, width_j = smoothed * df, derivative * frac
    return [[all_j[i, j] for j in range(D)] + [width_j[i, 0]] for i in range(N - 1)]


def run(p, case, out):
    schedule = read(CONTRACT / 'experimental_schedule.tsv')
    replay = {(r['cohort'], r['segment_id']): r for r in read(BRIDGE / 'verification/R01/selection.tsv')}
    passages = {r['passage_id']: r for r in read(STEP1 / 'passage_crosswalk.tsv')}
    assert len(schedule) == 122 and len(passages) == 131
    archive = defaultdict(dict)
    for row in read(ARCHIVE / 'invitro_daily_counts.tsv'):
        archive[(row['cohort'], row['segment_id'])][row['day']] = float(row['live_cells'])
    for row in schedule:
        key = (row['cohort'], row['segment_id'])
        assert replay[key]['day'] == row['selected_day']
        assert Fraction(key[1].split('::')[0].split('_')[-1]) == Fraction(row['oxygen_pct'])
        assert (Fraction(row['oxygen_pct']) / Fraction('0.01')).denominator == 1
    ops = operators(p, sorted({r['oxygen_pct'] for r in schedule}, key=Fraction), out)
    kernel = model.kary_kernel(p['sigma_kary'])
    keep_to = {key: int(Fraction(r['day']) / Fraction(DT)) for key, r in replay.items()}
    for row in schedule:
        key = (row['cohort'], row['segment_id'])
        if row['parent_segment_id'] != 'NA':
            parent = (row['cohort'], row['parent_segment_id'])
            keep_to[parent] = max(keep_to[parent], int(Fraction(replay[key]['parent_day']) / Fraction(DT)))

    results, observations, certified_branches, labels = {}, [], [], []
    max_native_error = 0.0
    identity_rows = []
    for i, row in enumerate(schedule):
        cohort, seg = key = row['cohort'], row['segment_id']
        size = count_or_none(row['initial_cells'])
        if size is None:
            size = arb(FALLBACK) / CAPACITY
        parent = (cohort, row['parent_segment_id'])
        reset_day, reset_method, boundary_margin = 'NA', 'gaussian_root', 'fixed_by_design'
        if parent in results:
            states, target = results[parent]
            reset_day, reset_method, boundary_margin = select(states, target, size, replay[key]['parent_day'], key)
            x, jac, total = states[reset_day]
            x, jac = normalized_reseed(x, jac, total, size)
        else:
            assert parent[1] == 'NA'
            x, jac = initial(p, cohort, COLS, N0=size)
        target = count_or_none(row['final_cells'])
        if target is None and i + 1 < len(schedule) and schedule[i + 1]['cohort'] == cohort:
            target = count_or_none(schedule[i + 1]['initial_cells'])
        expected_target = count_or_none(row['target_live_cells'])
        assert (target is None and expected_target is None) or (target - expected_target).contains(0)
        days = row['candidate_days'].split(';')
        steps = int(Fraction(row['nominal_duration_days']) / Fraction(DT))
        states, crowd_margin = propagate(x, jac, steps, keep_to[key], days, ops[row['oxygen_pct']], key)
        chosen, method, selection_margin = select(states, target, None, replay[key]['day'], key)
        results[key] = (states, target)
        certified_branches.append(dict(cohort=cohort, segment=seg, oxygen_pct=row['oxygen_pct'],
                                       end_day=row['nominal_duration_days'], selected_day=chosen,
                                       selection_rule=method, selection_margin_as_fraction_K=selection_margin,
                                       parent_segment=parent[1], parent_day=reset_day, reset_rule=reset_method,
                                       reset_margin_as_fraction_K=boundary_margin,
                                       min_one_minus_L_over_K=str(crowd_margin)))
        if case == 'anchor':
            for day, (_, _, total) in states.items():
                native = archive[key][day]
                error = abs(float(total.mid()) * CAPACITY - native) / max(1.0, abs(native))
                max_native_error = max(max_native_error, error)
        x, jac, total = states[chosen]
        if any(passages[pid]['recorded_ols_growth_per_day'] != 'NA' for pid in row['data_ids'].split(';')):
            gr = growth_rows(jac, total, size, chosen)
        if any(int(passages[pid]['karyotype_cell_count']) for pid in row['data_ids'].split(';')):
            kr = kary_rows(x, jac, total, kernel)
            identity_kr = kary_rows(x, jac, total, (arb_mat([[int(j == k) for k in range(N)] for j in range(N)]),
                                                        arb_mat(N, N)))
        for pid in row['data_ids'].split(';'):
            if passages[pid]['recorded_ols_growth_per_day'] != 'NA':
                labels.append(dict(stream='growth', passage=pid, oxygen_pct=row['oxygen_pct'],
                                   cohort=cohort, local_day=chosen, chromosome=''))
                observations.append(gr)
                identity_rows.append(gr)
            if int(passages[pid]['karyotype_cell_count']):
                for j, n in enumerate(GRID[:-1]):
                    labels.append(dict(stream='kary_probability', passage=pid, oxygen_pct=row['oxygen_pct'],
                                       cohort=cohort, local_day=chosen, chromosome=n))
                    observations.append(kr[j])
                    identity_rows.append(identity_kr[j])
        if (i + 1) % 10 == 0 or i == len(schedule) - 1:
            print(f'segments {i + 1}/{len(schedule)}', flush=True)
    assert len(observations) == 114 + 12 * 132
    if case == 'anchor':
        assert max_native_error < 1e-9, max_native_error
    print('Maximum scaled saved-native total gap', max_native_error, flush=True)
    table(out / 'passage_branches.tsv', certified_branches)
    return observations, identity_rows, labels, max_native_error


def emin(value):
    mid = value.mid()
    if mid.is_zero():
        return None
    m, e = mid.man_exp()
    return abs(int(m)).bit_length() + int(e)


def nominate(rows, col_indices):
    cscale = []
    for j in col_indices:
        nonzeros = (e for row in rows if (e := emin(row[j])) is not None)
        cscale.append(-max(nonzeros))
    scaled = [[row[j] * arb((1, e)) for j, e in zip(col_indices, cscale)] for row in rows]
    # Column-scale the width before converting midpoints to binary64. No rank
    # conclusion relies on this numerical pivot strategy.
    matrix = np.array([[float(x.mid()) for x in row] for row in scaled])
    U, singular, _ = np.linalg.svd(matrix, full_matrices=False)
    assert np.isfinite(matrix).all()
    pool = U.copy()
    picked = []
    for _ in col_indices:
        index = int(np.argmax(np.linalg.norm(pool, axis=1)))
        assert index not in picked
        picked.append(index)
        direction = pool[index] / np.linalg.norm(pool[index])
        for _ in range(2):
            pool -= (pool @ direction)[:, None] * direction[None, :]
        pool[picked] = 0
    return scaled, cscale, picked, singular


def certify(rows, labels, cols, tag, out, rank_records):
    scaled, cscale, picked, singular = nominate(rows, cols)
    M = arb_mat([scaled[i] for i in picked])
    determinant = M.det()
    passed, neumann = independent_minor_check(M)
    record = dict(case=tag, dimension=len(cols), determinant=str(determinant),
                  determinant_excludes_zero=nonzero(determinant), neumann_upper=str(neumann),
                  neumann_certified=passed, smallest_scaled_singular_midpoint=float(singular[-1]),
                  largest_scaled_singular_midpoint=float(singular[0]), precision=ctx.prec)
    rank_records.append(record)
    table(out / 'rank_summary.tsv', rank_records)
    table(out / f'{tag}_rows.tsv', [dict(matrix_row=i, **labels[j]) for i, j in enumerate(picked)])
    table(out / f'{tag}_cols.tsv', [dict(matrix_column=i, parameter=OBS[j], rescale_power_of_2=e)
                                   for i, (j, e) in enumerate(zip(cols, cscale))])
    save_matrix(out / f'{tag}_matrix.tsv', M)
    table(out / f'{tag}_determinant.tsv', [encode(determinant)])
    print(record, flush=True)
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=['anchor', 'interior'], required=True)
    parser.add_argument('--precision', type=int, default=384)
    args = parser.parse_args()
    ctx.prec = args.precision
    ctx.threads = 1
    out = Path(__file__).parent / f'{args.case}_{args.precision}'
    out.mkdir(exist_ok=True)
    p, metadata = anchor(args.case)
    table(out / 'exact_parameters.tsv', metadata)
    assert p['O2_crit'] > exact('1e-12') and p['n_O'] > 0 and p['sigma_kary'] > 0
    assert p['buffer_smax'] > 0 and (p['buffer_smax'] <= 1 or (p['buffer_smax'] - 1).is_zero())
    assert all(p[name] > 0 for name in MECHANISMS)
    h = (-(1 / (2 * p['sigma_kary']**2))).exp()
    print('First off-diagonal unnormalized weight', h, flush=True)
    rows, identity_rows, labels, max_error = run(p, args.case, out)
    records = []
    certify(rows, labels, range(len(OBS)), 'real_kernel_19', out, records)
    assert all(row[-1].is_zero() for row in identity_rows)
    certify(identity_rows, labels, range(D), 'identity_kernel_18', out, records)
    growth = [i for i, label in enumerate(labels) if label['stream'] == 'growth']
    certify([rows[i] for i in growth], [labels[i] for i in growth], range(D),
            'growth_only_18', out, records)
    width_lower, width_upper = exact('0.014'), exact('0.015')
    assert p['sigma_kary'] > width_lower and p['sigma_kary'] < width_upper
    # The closed interval's raw diagonal is greater than one. Its largest
    # off-diagonal is below 2^-1075 even if the omitted normalizer 1/sqrt(2pi)
    # is bounded above by one. The derivative proof retains that tiny mass.
    bound = (-(1 / (2 * width_upper**2))).exp() / width_lower
    assert bound < arb((1, -1075))
    table(out / 'observer_and_replay.tsv', [dict(
        first_unnormalized_off_diagonal_weight=str(h), binary64_width_lower=str(width_lower),
        binary64_width_upper=str(width_upper), off_diagonal_density_upper=str(bound.upper()),
        half_smallest_subnormal=str(arb((1, -1075))), max_scaled_distance_to_archived_totals=max_error)])
    table(out / 'runtime.tsv', [dict(python=sys.version.split()[0], flint=flint.__version__,
                                   numpy=np.__version__, precision=ctx.prec)])
    for path in [Path(__file__), Path(__file__).parent / 'anchor.tsv',
                 CONTRACT / 'experimental_schedule.tsv', BRIDGE / 'verification/R01/selection.tsv',
                 STEP1 / 'passage_crosswalk.tsv', ARCHIVE / 'invitro_daily_counts.tsv',
                 Path(model.__file__), Path(sys.modules['algebra'].__file__)]:
        print('sha256', path.relative_to(ROOT), hashlib.sha256(path.read_bytes()).hexdigest(), flush=True)


if __name__ == '__main__':
    main()
