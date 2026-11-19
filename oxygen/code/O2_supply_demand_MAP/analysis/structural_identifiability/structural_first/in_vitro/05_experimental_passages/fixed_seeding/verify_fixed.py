"""Check fixed-inoculum certificates and replay root-composition derivatives."""
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from flint import arb, arb_mat, ctx

import certify_fixed as fixed
from algebra import decode, independent_minor_check
from protocol import MECHANISMS, read, table

cert = fixed.cert
BASE = fixed.HERE


def grouped(rows, key):
    out = defaultdict(list)
    for row in rows:
        out[row[key]].append(row)
    return out


def decode_blocks(case_dir, expected_ids):
    blocks = grouped(read(case_dir / "root_minors.tsv"), "root_id")
    entries_seen = set()
    assert len(blocks) == 14
    worst_bound, dims = arb(0), 0
    first_columns, first_values = {}, {}
    for root, rows in blocks.items():
        outputs = tuple(dict.fromkeys(r["output_id"] for r in rows))
        cols = tuple(dict.fromkeys(int(r["chromosome"]) for r in rows))
        assert len(rows) == len(outputs) * len(cols) and len(outputs) == len(cols)
        m = {(r["output_id"], int(r["chromosome"])): decode(r) for r in rows}
        assert len(rows) == len(m)
        matrix = arb_mat([[m[out, col] for col in cols] for out in outputs])
        ok, bound = independent_minor_check(matrix)
        assert ok and bound < 1, (root, bound)
        worst_bound = max(worst_bound, bound)
        dims += len(outputs)
        entries_seen.update(outputs)
        first_columns[root] = cols[0]
        first_values.update({o: float(m[o, cols[0]].mid()) for o in outputs})
    assert entries_seen == expected_ids and len(entries_seen) == dims
    return dims, float(worst_bound.upper()), first_columns, first_values


def make_numeric_model():
    with TemporaryDirectory(prefix="fixed_inoculum_replay_") as tmp:
        tmp = Path(tmp)
        ctx.prec = 192
        p = cert.parameter_point(tmp)
        design, _, _ = cert.load_design(tmp)
        operators = cert.operators(p, sorted({d.oxygen for d in design.values()}), tmp)
        return float(p["k_clear"]), {
            o: (np.array([[float(G[i,j]) for j in range(cert.N)] for i in range(cert.N)]),
                np.array([float(H[0,j]) for j in range(cert.N)]),
                np.array([float(M[0,j]) for j in range(cert.N)]))
            for o, (G, H, M) in operators.items()
        }


def replay(case_dir, phase, visibility, numerical, perturbed_chr, epsilon):
    clearance, operators = numerical
    schedule = grouped(read(cert.STEP2 / "ols_calendar_map.tsv"), "passage_id")
    active = read(case_dir / "active_passages.tsv")
    children = grouped(active, "retained_parent_id")
    states = {}
    counts = {}
    outputs = {}
    for row in active:
        key = row["passage_id"]
        root = row["root_id"]
        if row["retained_parent_id"]:
            supply = states[row["retained_parent_id"]][int(row["parent_donor_step"])]
            pi = supply / supply.sum()
        else:
            center = 46 if row["cohort"] == "2N" else 88
            chr_values = np.array(cert.model.GRID)
            pi = 1 / (1 + ((chr_values - center) / 9)**4)
            pi /= pi.sum()
            direction = perturbed_chr[root]
            pi[cert.model.GRID.index(direction)] += epsilon
            pi[cert.model.GRID.index(center)] -= epsilon
            assert np.all(pi > 0)
        out = sorted(schedule.get(key, ()), key=lambda r: int(r["euler_step"]))
        wanted = set(int(s["euler_step"]) for s in out)
        donor_ticks = set(int(c["parent_donor_step"]) for c in children.get(key, ()))
        wanted.update(donor_ticks)
        x = int(row["seed_corrected_count"]) / cert.CAPACITY * pi
        H = M = 0.
        G, direct, biological = operators[row["oxygen_pct"]]
        eta = (0.4, 0.6) if visibility == "split" and row["cohort"] == "4N" else (1/3, 2/3)
        observed = {}
        donor = {}
        for step in range(max(wanted) + 1):
            live = x.sum()
            if step in wanted:
                observed[step] = live + eta[0] * H + eta[1] * M
            if step in donor_ticks:
                donor[step] = x.copy()
            s = 0.05 * (1-live)
            H += 0.05 * (direct @ x - clearance * H)
            M += s * (biological @ x) - 0.05 * clearance * M
            x = x + s * (G @ x) - 0.05 * (direct * x)
        states[key] = donor
        if out:
            counts.update({r["record_id"]: observed[int(r["euler_step"])] * cert.CAPACITY for r in out})
            if phase == "ols":
                outputs[key] = sum(np.log(observed[int(r["euler_step"])] * cert.CAPACITY) *
                                   int(r["ols_weight_numerator"]) / int(r["ols_weight_denominator"]) for r in out)
            else:
                outputs.update({r["record_id"]: observed[int(r["euler_step"])] for r in out[1:]})
    return counts, outputs


def main():
    numeric = make_numeric_model()
    schedule = read(cert.STEP2 / "ols_calendar_map.tsv")
    checks = []
    for phase in ("ols", "counts"):
        for visibility in ("shared", "split"):
            case_dir = BASE / f"{phase}_{visibility}"
            need_ids = {r["passage_id"] for r in schedule} if phase == "ols" else {
                r["record_id"] for r in schedule if int(r["euler_step"]) > 0
            }
            ctx.prec = 1024
            dims, bound, first_cols, first_vals = decode_blocks(case_dir, need_ids)
            baseline, _ = replay(case_dir, phase, visibility, numeric, first_cols, 0)
            _, plus = replay(case_dir, phase, visibility, numeric, first_cols, 1e-7)
            _, minus = replay(case_dir, phase, visibility, numeric, first_cols, -1e-7)
            certified = {r["record_id"]: float(arb(r["witness_visible_cells"]).mid())
                         for r in read(case_dir / "witness_counts.tsv")}
            assert len(baseline) == 335 and baseline.keys() == certified.keys()
            assert plus.keys() == minus.keys() == first_vals.keys()
            count_err = max(abs(baseline[k] / v - 1) for k, v in certified.items())
            diff_err = max(abs((plus[k] - minus[k]) / 2e-7 - first_vals[k]) for k in plus)
            scale = max(abs(v) for v in first_vals.values())
            rel_scale = diff_err / scale
            assert count_err < 5e-13 and rel_scale < 1e-4, (phase, visibility, count_err, rel_scale)
            checks.append(dict(phase=phase, visibility=visibility, dimension=dims,
                               max_exact_midpoint_neumann_bound=bound, relative_count_error=count_err,
                               root_derivative_abs_error_over_max_root_derivative=rel_scale))
            print(phase, visibility, checks[-1], flush=True)
    table(BASE / "verification_summary.tsv", checks)
    table(BASE / "parameter_verdicts.tsv", [dict(
        mechanism=name, ols_shared_nuisance_and_joint_rank="114/114", ols_split_nuisance_and_joint_rank="114/114",
        raw_shared_nuisance_and_joint_rank="221/221", raw_split_nuisance_and_joint_rank="221/221",
        ols_conclusion="exact local compensation through restart composition only",
        counts_conclusion="exact local compensation through restart composition only",
        initial_viable_cells="fixed to recorded correctedCount on every retained segment",
        day_zero_visible_dead_cells="zero; mandated by fixed live and observed visible counts",
    ) for name in MECHANISMS])


if __name__ == "__main__":
    main()
