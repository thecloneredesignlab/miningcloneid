"""Fixed viable inocula, with no visible dead stock at seeding.

Certify output rank using only admissible directions of the fourteen
independent restart karyotype distributions. Keep all mechanisms free.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
import sys

from flint import arb, arb_mat, ctx

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / "20261118_culture_experimental_OLS_rank_01"
sys.path.insert(0, str(PREVIOUS))
import certify_experimental as cert
from algebra import encode, exact
from protocol import MECHANISMS, read, table


def mean_rows(psg, samples, phase, k):
    if phase == "ols":
        # Day zero has zero derivative because the visible seed count is fixed.
        jac = [sum(exact(w) * samples[t][3][j] / samples[t][2]
                   for t, _, w in psg.samples[1:]) for j in range(k)]
        return [(psg.id, jac)]
    # Every modeled day-zero count is the fixed recorded count, with zero derivative.
    return [(record_id, samples[tick][3][:k]) for tick, record_id, _ in psg.samples[1:]]


def design_columns(passages, roots, phase):
    sizes = Counter()
    for p in passages.values():
        if p.samples:
            sizes[roots[p.id]] += 1 if phase == "ols" else len(p.samples) - 1
    active = 114 if phase == "ols" else 221
    cert.need(len(sizes) == 14 and sum(sizes.values()) == active, "missing nonconstant output rows")
    selected = {}
    for root, size in sizes.items():
        ref = 46 if passages[root].cohort == "2N" else 88
        choices = [i for i, n in enumerate(cert.model.GRID) if n != ref]
        # Coordinate changes transfer mass between a selected class and ref.
        selected[root] = [choices[j * (len(choices) - 1) // max(1, size - 1)] for j in range(size)]
    return active, selected


def summary(phase, visibility, bits, out, active, verified):
    total = 114 if phase == "ols" else 335
    all_verified = verified == active
    table(out / "rank.tsv", [dict(
        phase=phase, visibility=visibility, precision_bits=bits, registered_observations=total,
        constant_day_zero_rows=total - active, nuisance_rank_lower_bound=verified,
        nuisance_and_joint_rank_upper_bound=active, full_nuisance_rank_certified=all_verified,
        mechanisms_with_exact_local_compensation=14 if all_verified else "",
    )])
    print(f"{phase} {visibility} at {bits} bits: rank interval [{verified}, {active}]", flush=True)
    return all_verified


def run(phase, visibility, bits, out):
    ctx.prec = bits
    out.mkdir(parents=True, exist_ok=True)
    passages, roots, uses = cert.load_design(out)
    p = cert.parameter_point(out)
    ops = cert.operators(p, sorted({p.oxygen for p in passages.values()}, key=Fraction), out)
    decay = 1 - exact(cert.DT) * p["k_clear"]
    active, selected = design_columns(passages, roots, phase)

    samples_by_passage = {}
    children = defaultdict(list)
    for passage in passages.values():
        if passage.parent:
            children[(passage.parent, passage.donor_tick)].append(passage)
    matrices, row_ids = defaultdict(list), defaultdict(list)
    evidence, branches, matrix_entries, root_ranks = [], [], [], []
    verified = 0
    remaining = Counter(roots.values())

    for psg in passages.values():
        root_id = roots[psg.id]
        cols = selected[root_id]
        k = len(cols)
        eta = cert.etas(visibility, psg.cohort)
        if psg.parent:
            x, jac, _, _, L = samples_by_passage[psg.parent][psg.donor_tick]
            pi, pi_j = cert.normalized_state(x, jac, L, [arb(0)] * k)
        else:
            pi, ref = cert.profile(psg.cohort)
            pi_j = arb_mat([[int(i == j) - int(i == ref) for j in cols] for i in range(cert.N)])

        a = exact(Fraction(psg.seed_count, cert.CAPACITY))
        wanted = set(uses.get(psg.id, ())) | {r[0] for r in psg.samples}
        # a equals the initial visible count; hence the initial visible dead
        # contribution is zero. The local a column from the reused solver is
        # deliberately excluded from all rank matrices and all propagated resets.
        samples, crowd, visible = cert.propagate(pi, pi_j, a, a, wanted, ops[psg.oxygen], eta, decay)
        if psg.samples:
            # The first output is analytically constant, not an estimated rank.
            cert.need((samples[0][2] - a).contains(0), f"day-zero mismatch {psg.id}")
            cert.need(all(v.contains(0) for v in samples[0][3][:k]), f"day-zero nonconstant {psg.id}")
        for tick in uses.get(psg.id, ()):
            demand = sum(c.seed_count for c in children[(psg.id, tick)])
            available = samples[tick][4] * cert.CAPACITY
            margin = available - demand
            cert.need(margin > 0, f"donor cannot supply all fixed inocula: {psg.id} tick {tick}")
            branches.append(dict(passage_id=psg.id, donor_tick=tick, simultaneous_children=len(children[(psg.id, tick)]),
                                 fixed_demand_cells=demand, minimum_margin_cells=str(margin.lower())))
        if psg.samples:
            for label, values in mean_rows(psg, samples, phase, k):
                row_ids[root_id].append(label)
                matrices[root_id].append(values)
            for tick, label, _ in psg.samples:
                evidence.append(dict(record_id=label, passage_id=psg.id, euler_step=tick,
                                     witness_visible_cells=str(samples[tick][2] * cert.CAPACITY)))
        samples_by_passage[psg.id] = {t: samples[t] for t in uses.get(psg.id, ())}
        remaining[root_id] -= 1
        if not remaining[root_id]:
            matrix = arb_mat(matrices[root_id])
            cert.need(matrix.nrows() == matrix.ncols() == k, "wrong root-block shape")
            try:
                bound = cert.midpoint_preconditioned_bound(matrix)
                success = bool(bound < 1)
            except (ZeroDivisionError, ValueError):
                bound, success = arb("inf"), False
            if success:
                verified += k
            matrix_entries.extend(dict(root_id=root_id, output_id=label, chromosome=cert.model.GRID[col],
                                       **encode(matrix[i, j]))
                                  for i, label in enumerate(row_ids[root_id])
                                  for j, col in enumerate(cols))
            root_ranks.append(dict(root_id=root_id, dimension=k, certified=success, neumann_error_upper=str(bound)))
            print(root_id, k, success, str(bound)[:130], flush=True)

    table(out / "root_minors.tsv", matrix_entries)
    table(out / "root_ranks.tsv", root_ranks)
    table(out / "donor_supply.tsv", branches)
    table(out / "witness_counts.tsv", evidence)
    cert.need(len(evidence) == 335 and len(root_ranks) == 14, "incomplete rank run")
    return summary(phase, visibility, bits, out, active, verified)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", required=True, choices=("ols", "counts"))
    parser.add_argument("--visibility", required=True, choices=("shared", "split"))
    parser.add_argument("--bits", type=int, default=768)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if not run(args.phase, args.visibility, args.bits, args.out):
        sys.exit(2)
