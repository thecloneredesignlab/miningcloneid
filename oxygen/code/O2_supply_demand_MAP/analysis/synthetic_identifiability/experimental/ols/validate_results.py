"""Replay selected fits and check that every reported decision uses the frozen rule."""
from collections import defaultdict
from hashlib import sha256

import numpy as np

from fit_stage import model
from inference import SLACK, THRESHOLD, likelihood_ratio
from practical import Configuration, HERE, MECHANISMS, read, row_x, write
from recovery import CONVERGE_LR, KKT_TOL


def run():
    k = model()
    c = Configuration(k)
    assert k.n == 114
    for r in read(HERE / "source_provenance.tsv"):
        assert sha256((HERE.parents[2] / r["path"]).read_bytes()).hexdigest() == r["sha256"], r["path"]
    data = defaultdict(list)
    for r in read(HERE / "noisy_slope_inputs.tsv"):
        data[r["scenario"], int(r["dataset"])].append((r["passage_id"], float(r["ols_slope"])))
    truths = {r["scenario"]: np.array([float(r[q]) for q in c.names])
              for r in read(HERE / "synthetic_truth.tsv")}
    expected = {(case, i) for case in truths for i in range(1, 101)}
    assert set(data) == expected
    assert all(tuple(p for p, _ in v) == tuple(k.ids) for v in data.values())
    fits = read(HERE / "qualified_recovery.tsv")
    assert len(fits) == len(expected) and {(r["scenario"], int(r["dataset"])) for r in fits} == expected
    main = read(HERE / "main_constrained_best.tsv")[0]
    summary = {r["parameter"]: r for r in read(HERE / "parameter_practical_verdicts.tsv")}
    reported = {(r["scenario"], r["parameter"]): r for r in read(HERE / "synthetic_parameter_recovery.tsv")}
    assert set(summary) == set(MECHANISMS)
    assert set(reported) == {(case, name) for case in truths for name in MECHANISMS}
    profiles = defaultdict(list)
    for r in read(HERE / "main_profiles_envelope.tsv"):
        if r["feasible"] == "True":
            assert abs(float(r["lr"]) - likelihood_ratio(float(r["loss"]), float(main["loss"]), 114)) < 1e-10
            profiles[r["profile_parameter"]].append(r)
    assert set(profiles) == set(MECHANISMS)
    output = []
    for case, truth in truths.items():
        qualified = []
        errors, donors = [], []
        for r in fits:
            key = r["scenario"], int(r["dataset"])
            if key[0] != case or r["qualified"] != "True":
                continue
            x = row_x(c, r)
            target = np.array([y for _, y in data[key]])
            f, _, g = k.constrained(c.unmap(x)[0], target, c, gradient=False)
            assert g >= -SLACK and abs(g - float(r["min_aggregate_donor_margin_fraction"])) < 1e-9
            assert float(r["max_projected_gradient"]) <= KKT_TOL
            assert float(r["independent_LR_gap"]) <= CONVERGE_LR
            errors.append(abs(f - float(r["loss"])))
            donors.append(g)
            qualified.append(x)
        assert max(errors) < 1e-8
        counts = np.sum(np.abs(np.array(qualified) - truth) <= .10, axis=0)
        for j, name in enumerate(MECHANISMS):
            row = reported[case, name]
            assert int(row["assigned"]) == 100 and int(row["qualified"]) == len(qualified)
            assert int(row["recovered_count"]) == counts[j]
        output.append(dict(scenario=case, datasets=100, slope_rows=11400, qualified=len(qualified),
                           largest_half_RSS_discrepancy=max(errors), smallest_supply_margin=min(donors)))
    for name, row in summary.items():
        x = float(row["estimate_normalized"])
        outside = any(abs(float(r["fixed_normalized"]) - x) > .1000000001
                      and float(r["lr"]) <= THRESHOLD for r in profiles[name])
        lower, upper = (row[f"profile95_normalized_{edge}"] for edge in ("lower", "upper"))
        narrow = bool(lower and upper) and max(abs(float(lower)-x), abs(float(upper)-x)) <= .10
        recovery_passes = all(int(reported[case, name]["recovered_count"]) >= 90 for case in truths)
        assert (row["verdict"] == "precision_gate_passed_conditional") == (narrow and not outside and recovery_passes)
        assert (row["verdict"] == "insufficient_profile_precision") == outside
    write(HERE / "final_validation.tsv", output)
    for row in output:
        print(row)


if __name__ == "__main__":
    run()
