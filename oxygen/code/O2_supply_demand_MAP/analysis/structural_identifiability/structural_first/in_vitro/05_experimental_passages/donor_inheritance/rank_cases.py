"""Apply OLS-first gating and certify Jacobian minors for eight nested designs."""
from collections import defaultdict

from flint import arb, arb_mat, ctx
import numpy as np

import inherit as study
from algebra import decode, encode, independent_minor_check
from protocol import INITIALS, MECHANISMS, read, table

HERE = study.HERE
COHORTS = ("2N", "4N")


def block(root, phase, visibility, transfer):
    by_cohort = {}
    for cohort in COHORTS:
        items = read(HERE / root / f"{cohort}_{phase}_derivatives.tsv")
        rows = defaultdict(dict)
        for edge in items:
            rows[edge["output_id"]][edge["column"]] = decode(edge)
        by_cohort[cohort] = rows
    columns = list(MECHANISMS)
    if visibility == "shared":
        columns += ["eta_H", "eta_M"]
    else:
        columns.extend(f"eta_{cause}_{cohort}" for cohort in COHORTS for cause in ("H", "M"))
    if root == "gaussian":
        columns += INITIALS
    if transfer == "tilt":
        for cohort in COHORTS:
            row = next(iter(by_cohort[cohort].values()))
            columns.extend(c for c in row if c.startswith("tilt_"))
    ids = []
    matrix = []
    for cohort, rows in by_cohort.items():
        for key, data in rows.items():
            ids.append(key)
            matrix.append([data.get(k if k not in ("eta_H", "eta_M") else f"{k}_{cohort}", arb(0))
                           for k in columns])
    assert len(matrix) == (114 if phase == "ols" else 221) and len(ids) == len(set(ids))
    return columns, ids, matrix


def row_selector(data):
    """Greedy row selection. Float values nominate rows, balls certify them."""
    orthogonal = np.array([[float(x.mid()) for x in row] for row in data])
    used = []
    for _ in range(len(data[0])):
        lengths = np.sum(orthogonal**2, axis=1)
        if used:
            lengths[used] = -1
        pivot = int(np.argmax(lengths))
        if lengths[pivot] <= 0:
            break
        used.append(pivot)
        u = orthogonal[pivot] / np.sqrt(lengths[pivot])
        for _ in range(2):
            orthogonal -= (orthogonal @ u)[:, None] * u[None, :]
    return used


def higher_precision_selection(values):
    rows = [[x.mid() for x in r] for r in values]
    col_ids = list(range(len(rows[0])))
    row_ids = list(range(len(rows)))
    n = len(col_ids)
    completed = 0
    for k in range(n):
        i_best, j_best, magnitude = k, k, arb(0)
        for i in range(k, len(rows)):
            for j in range(k, n):
                size = abs(rows[i][j])
                if size > magnitude:
                    i_best, j_best, magnitude = i, j, size
        if magnitude.is_zero():
            break
        rows[k], rows[i_best] = rows[i_best], rows[k]
        row_ids[k], row_ids[i_best] = row_ids[i_best], row_ids[k]
        if j_best != k:
            for row in rows:
                row[k], row[j_best] = row[j_best], row[k]
            col_ids[k], col_ids[j_best] = col_ids[j_best], col_ids[k]
        pivot = rows[k][k]
        for i in range(k+1, len(rows)):
            scale = (rows[i][k] / pivot).mid()
            for j in range(k+1, n):
                rows[i][j] = (rows[i][j] - scale * rows[k][j]).mid()
            rows[i][k] = arb(0)
        completed = k+1
    return row_ids[:completed], col_ids[:completed]


def certify(phase, root, visibility, transfer, receipts):
    cols, ids, original = block(root, phase, visibility, transfer)
    n = len(cols)
    powers = []
    for i in range(n):
        mx = max(abs(x[i].mid()) for x in original)
        powers.append(0 if mx.is_zero() else mx.man_exp()[1] + abs(int(mx.man_exp()[0])).bit_length())
    scaled = [[v / arb((1, power)) for v, power in zip(row, powers)] for row in original]
    selected = row_selector(scaled)
    chosen_columns = list(range(n))
    success, bound = False, arb("inf")
    if len(selected) == n:
        chosen = arb_mat([scaled[i] for i in selected])
        try:
            bound = study.cert.midpoint_preconditioned_bound(chosen)
            success = bool(bound < 1)
        except (ZeroDivisionError, ValueError):
            pass
    if not success:
        selected, chosen_columns = higher_precision_selection(scaled)
        while selected:
            chosen = arb_mat([[scaled[i][k] for k in chosen_columns] for i in selected])
            try:
                bound = study.cert.midpoint_preconditioned_bound(chosen)
                if bound < 1:
                    break
            except (ZeroDivisionError, ValueError):
                pass
            selected, chosen_columns = selected[:-1], chosen_columns[:-1]
        success = len(selected) == n
    certified = len(selected)
    independent_bound = arb("inf")
    if certified:
        independent, independent_bound = independent_minor_check(chosen)
        assert independent and independent_bound < 1
        receipts.extend(dict(phase=phase, root=root, visibility=visibility, transfer=transfer,
                             observation_id=ids[i], column=cols[k], scale_exponent=powers[k],
                             **encode(chosen[row_idx, col_idx]))
                        for row_idx, i in enumerate(selected)
                        for col_idx, k in enumerate(chosen_columns))
    tag = "all_14_locally_identifiable" if success else "full_rank_not_certified"
    record = dict(observation=phase, dp_root=root, visibility=visibility, transfer=transfer,
                  observation_rows=len(ids), mechanism_columns=14, nuisance_columns=n-14,
                  certified_rank_lower=certified, rank_upper=n,
                  independent_neumann_error_bound=str(independent_bound), verdict=tag)
    print(phase, root, visibility, transfer, certified, "/", n, tag, flush=True)
    return record


def run():
    ctx.prec = 512
    rows, receipts, params = [], [], []
    for root in ("hist", "gaussian"):
        for transfer in ("exact", "tilt"):
            for visibility in ("shared", "split"):
                ols = certify("ols", root, visibility, transfer, receipts)
                rows.append(ols)
                result = ols
                if ols["verdict"] != "all_14_locally_identifiable":
                    result = certify("counts", root, visibility, transfer, receipts)
                    rows.append(result)
                for mechanism in MECHANISMS:
                    params.append(dict(mechanism=mechanism, dp_root=root, visibility=visibility,
                                       transfer=transfer, evidence=result["observation"],
                                       verdict="certified_local" if result["verdict"] == "all_14_locally_identifiable"
                                       else "requires_individual_analysis"))
    table(HERE / "rank_summary.tsv", rows)
    table(HERE / "rank_minors.tsv", receipts)
    table(HERE / "mechanism_verdicts_all_cases.tsv", params)
    table(HERE / "parameter_summary.tsv", [
        dict(mechanism=name, **{
            f"{r['dp_root']}_{r['transfer']}_{r['visibility']}": f"{r['verdict']} ({r['evidence']})"
            for r in params if r["mechanism"] == name
        }) for name in MECHANISMS
    ])
    return rows


if __name__ == "__main__":
    run()
