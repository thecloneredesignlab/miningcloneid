#!/usr/bin/env python3
"""Audit completed Sobol exports against the culture contract; never run a model."""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_table(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_table(path, rows):
    with path.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def key(row):
    return (round(float(row["O2_pct"]), 12), row["output"], row["index"],
            row["parameter_i"], row["parameter_j"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("project", type=Path)
    parser.add_argument("old_project", type=Path)
    parser.add_argument("study", type=Path)
    args = parser.parse_args()
    project, old, study = (p.resolve() for p in (args.project, args.old_project, args.study))
    sobol = old / "oxygen/results/eFAST/sobol"
    comparison = sobol / "global_convergence_20261021/comparison_N65536_N131072"
    wp0 = project / "oxygen/results/structural_identifiability/20261109_seed144_wp0"
    output = study / "sobol_bridge"
    output.mkdir(exist_ok=True)
    sources = []

    def source(path):
        sources.append(dict(path=str(path), sha256=sha(path), bytes=path.stat().st_size))
        return path

    def load_json(path):
        return json.loads(source(path).read_text())

    source(Path(__file__).resolve())
    source(project / "docs/sobol_reference_integration.md")
    rows = read_table(source(comparison / "all_index_comparison.tsv.gz"))
    bounds = read_table(source(sobol / "inputs/global_bounds.tsv"))
    parameters = [r["parameter"] for r in bounds]
    assert len(set(parameters)) == 14
    assert len(rows) == 47838 and len({key(r) for r in rows}) == len(rows)
    assert {round(float(r["O2_pct"]), 12) for r in rows} == {round(i * .025, 12) for i in range(201)}
    assert Counter((r["output"], r["index"]) for r in rows) == Counter({
        (out, index): count for out in ["dominant_mean_ploidy", "dominant_growth_rate"]
        for index, count in [("S1", 2814), ("ST", 2814), ("S2", 18291)]})
    expected_pairs = {frozenset(x) for x in itertools.combinations(parameters, 2)}
    cells = defaultdict(list)
    for r in rows:
        cells[(round(float(r["O2_pct"]), 12), r["output"])].append(r)
    assert len(cells) == 402
    for cell in cells.values():
        for index in ["S1", "ST"]:
            assert {r["parameter_i"] for r in cell if r["index"] == index} == set(parameters)
        pair_cell = [r for r in cell if r["index"] == "S2"]
        assert len(pair_cell) == 91
        assert {frozenset([r["parameter_i"], r["parameter_j"]]) for r in pair_cell} == expected_pairs
    for out in ["dominant_mean_ploidy", "dominant_growth_rate"]:
        assert {frozenset([r["parameter_i"], r["parameter_j"]]) for r in rows
                if r["index"] == "S2" and r["output"] == out} == expected_pairs
    assert all(r["convergence_status_N131072"] == "passed" for r in rows)
    assert all(float(r["repeat_range_N131072"]) <= .05 and
               float(r["resolution_delta_N131072"]) <= .02 for r in rows)

    # This repeat export covers ST and 36 focal S2 pairs, not S1 or all 91 S2.
    repeats = read_table(source(comparison / "five_repeat_source.tsv.gz"))
    grouped = defaultdict(list)
    for r in repeats:
        if r["N"] == "131072":
            grouped[key(r)].append((int(r["repeat"]), float(r["value"])))
    assert len(grouped) == 20100
    lookup = {key(r): r for r in rows}
    max_mean_error = max_range_error = 0.0
    for k, values in grouped.items():
        assert sorted(i for i, _ in values) == [1, 2, 3, 4, 5]
        v = [v for _, v in values]
        max_mean_error = max(max_mean_error, abs(sum(v) / 5 - float(lookup[k]["value_mean_N131072"])))
        max_range_error = max(max_range_error, abs(max(v) - min(v) - float(lookup[k]["repeat_range_N131072"])))
    assert max_mean_error < 1e-12 and max_range_error < 1e-12
    compact = read_table(source(comparison / "parameter_and_pair_summary.tsv"))
    for r in compact:
        if r["N"] != "131072" or r["index"] != "ST":
            continue
        values = sorted((float(x["O2_pct"]), float(x["value_mean_N131072"])) for x in rows
                        if x["output"] == r["output"] and x["index"] == "ST" and x["parameter_i"] == r["parameter_i"])
        for lo, hi, field in [(0, 1, "low_integral_mean"), (1, 3, "mid_integral_mean"), (3, 5, "high_integral_mean")]:
            band = [(x, y) for x, y in values if lo <= x <= hi]
            average = sum((b[0] - a[0]) * (a[1] + b[1]) / 2 for a, b in zip(band, band[1:])) / (hi - lo)
            assert abs(average - float(r[field])) < 1e-12
    for name in ["interpretation_zh.md", "validation_final.json"]:
        source(comparison / name)
    for name in ["model_equations.md", "parameter_equation_map.tsv", "source_manifest.json"]:
        source(sobol / "model_structure" / name)
    audit = load_json(sobol / "inputs/source_audit.json")
    plan = load_json(sobol / "global_convergence_20261021/plan.json")
    completion = load_json(sobol / "global_convergence_20261021/N131072/completion.json")
    global_bounds = source(sobol / "full_N4096/global/full/global/parameter_bounds.tsv")
    global_case = load_json(sobol / "full_N4096/global/full/global/case.json")
    assert plan["scope"] == global_case["scope"] == completion["scope"] == "global"
    assert completion["status"] == "complete" and completion["N"] == 131072
    assert completion["n_repeats"] == 5
    assert sha(global_bounds) == global_case["bounds_sha256"]
    assert read_table(global_bounds) == bounds
    assert sha(global_bounds) in plan["source_files"].values()
    contract = load_json(wp0 / "manifests/production_contract.json")
    configuration = load_json(wp0 / "manifests/runtime_model_configuration.json")
    cfg = configuration["runtime_cfg"]
    assert [cfg[n] for n in ["N_MIN", "N_MAX", "N_UNIT"]] == [22, 154, 22]
    assert cfg["O2_growth"] is True and cfg["ploidy_O2_death"] == "ploidy_related"
    source(study / "study.json")
    source(study / "manifests/code_manifest.tsv")
    production_table = source(wp0 / "source/seed144/parameter_table_input.csv")
    with production_table.open() as stream:
        culture = [r for r in csv.DictReader(stream) if r["use_invitro_fit"].lower() == "true"]
    assert len(culture) == 20
    linear = {"buffer_smax", "gamma_growth", "gamma_mu", "n_O", "init_mean_2N", "init_mean_4N"}
    old_map = {r["parameter"]: r for r in bounds}
    domain_rows = []
    for r in culture:
        name = r["param_symbol"]
        old_row = old_map.get(name, {})
        domain_rows.append(dict(parameter=name, old_lower=old_row.get("natural_lower", ""),
            old_upper=old_row.get("natural_upper", ""), old_sampling_distribution=old_row.get("sampling_distribution", "not sampled"),
            culture_production_lower=r["lower_bound"], culture_upper=r["upper_bound"],
            culture_extended_lower=0 if name in ["p_misseg", "buffer_beta"] else r["lower_bound"],
            culture_transform="log1p" if name in ["p_misseg", "buffer_beta"] else "identity" if name in linear else "log10",
            role="mechanism" if name in old_map else "initial state" if name.startswith("init_") else "observation error"))
    write_table(output / "parameter_domain_comparison.tsv", domain_rows)
    code_rows = []
    model_base = Path("oxygen/code/O2_supply_demand_MAP")
    for relative in ["model/model_O2_supply_demand_MAP.cpp", "model/model_O2_supply_demand_MAP.R",
                     "util/o2_supply_demand_map_common_semantics.R", "util/o2_supply_demand_map_invitro_objective_utils.R",
                     "util/o2_supply_demand_map_fit_invitro_backend.R"]:
        row = dict(file=relative)
        for label, root in [("old", old), ("current", project), ("production_v2", wp0 / "source/code")]:
            p = root / model_base / relative
            row[label + "_sha256"] = sha(source(p)) if p.exists() else "absent"
        row["old_matches_production"] = row["old_sha256"] == row["production_v2_sha256"]
        code_rows.append(row)
    assert all(r["old_matches_production"] for r in code_rows[:2])
    write_table(output / "code_comparison.tsv", code_rows)
    for relative in ["analysis/eFAST/sobol/evaluate.R", "simulation/o2/fixed_o2/run_fixed_o2_simulation.R",
                     "analysis/eFAST/sobol/collect_global_convergence_exports.py",
                     "analysis/eFAST/sobol/compare_global_convergence_exports.py"]:
        source(old / model_base / relative)

    target_st = [r for r in compact if r["N"] == "131072" and r["index"] == "ST" and
                 r["parameter_i"] in ["p_misseg", "buffer_beta"]]
    write_table(output / "target_ST_reference.tsv", target_st)
    pairs = [("p_misseg", "mu_hp"), ("p_misseg", "gamma_mu"), ("p_misseg", "k_o_mis"),
             ("p_misseg", "p_mis_base"), ("buffer_beta", "buffer_smax"),
             ("buffer_beta", "buffer_n_exp"), ("buffer_smax", "buffer_n_exp")]
    pair_rows = []
    for pair in pairs:
        for out in ["dominant_mean_ploidy", "dominant_growth_rate"]:
            subset = [r for r in rows if r["output"] == out and r["index"] == "S2" and
                      frozenset([r["parameter_i"], r["parameter_j"]]) == frozenset(pair)]
            peak = max(subset, key=lambda r: float(r["value_mean_N131072"]))
            pair_rows.append(dict(output=out, parameter_i=pair[0], parameter_j=pair[1],
                peak_S2=float(peak["value_mean_N131072"]), peak_O2_pct=float(peak["O2_pct"]),
                interpretation="nondirectional interaction; compensation untested"))
    write_table(output / "candidate_pair_reference.tsv", pair_rows)
    hypotheses = [dict(target="p_misseg", track="mu_hp;gamma_mu;k_o_mis;p_mis_base;O2_crit;n_O",
        diagnostic="Compare nuisance paths, effective death rates, saturation mu/(mu+k_o_mis), total and induced per-copy error, and daughter flux at matched states and actual passage histories.",
        test="Does refitting preserve the observations and effective output as induction changes?", consequence="Candidate compensation only after refit; no direction inferred from S2."),
        dict(target="buffer_beta", track="buffer_smax;buffer_beta;buffer_n_exp",
        diagnostic="Compare survival at matched parental chromosome count and error size, viable daughter kernels and integrated output; track growth/death and observation-error compensation too.",
        test="Is the survival function or viable output more stable than individual coordinates?", consequence="Small old S2 does not rule out culture compensation."),
        dict(target="both", track="all 20 fitted coordinates;selected passage days;growth;karyotype;flow",
        diagnostic="Retain all nuisance parameters and boundary/search flags. Compare component traces and complete distribution predictions.",
        test="Which ambiguity remains for a later targeted stream-deletion refit or finite-time experimental comparison?", consequence="No stream attribution without its own unrestricted and restricted refits.")]
    write_table(output / "profile_diagnostic_hypotheses.tsv", hypotheses)
    record = dict(generated_utc=datetime.now(timezone.utc).isoformat(), method="completed global Sobol N131072",
        status="passed_read_only_bridge_audit", model_calls=0, jobs_submitted=0, archived_neighborhood_resumed=False,
        original_fit_configuration_changed=False, global_index_cells=len(rows), all_91_pairs_verified=True,
        all_cells_pass_operational_criteria=True, focal_repeat_cells_verified=len(grouped),
        max_focal_repeat_mean_error=max_mean_error, max_focal_repeat_range_error=max_range_error,
        source_audit_caveat="source_audit.json describes neighborhood-derived preparation; global case/bounds and frozen global plan are the authority for this result.",
        old_source_root=audit["source_root"], culture_contract_sha256=sha(wp0 / "manifests/production_contract.json"),
        shared_configuration=dict(N_MIN=22, N_MAX=154, N_UNIT=22, O2_growth=True, ploidy_O2_death="ploidy_related"),
        culture_runtime_configuration=configuration["runtime_cfg"],
        source_files=sources)
    (output / "audit.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({k:record[k] for k in ["status", "global_index_cells", "focal_repeat_cells_verified", "model_calls", "jobs_submitted"]}, indent=2))


if __name__ == "__main__":
    main()
