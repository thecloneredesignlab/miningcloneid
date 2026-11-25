"""Rebuild the specimen, observer replay and case tables without optimization."""

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime
from hashlib import sha256
from math import exp, isclose, log
from pathlib import Path
import os
import subprocess
import sys

from paths import (
    ALIGNMENT as ALIGN, FIT, PROFILES as PROFILE, ROOT,
    STRUCTURAL as STRUCT, SYNTHETIC as SYN, output_directory,
)

FREEZE = SYN / "20261122_experimental_joint_observation_freeze_01"
REVIEW = SYN / "20261122_experimental_joint_identifiability_limited_review_01"
FIXED = STRUCT / "20261118_culture_experimental_fixed_inoculum_01" / "counts_shared"
INHERIT = STRUCT / "20261118_culture_experimental_inheritance_01" / "hist"
CONTRACT = SYN / "20261113_culture_contract_steps01_03"
FCS = Path.home() / "Documents/Projects/FCS/outputs/01a1cc52-329e-7b21-b548-3c0ee7890a7c/gate_counts.csv"
STEP = 0.05
N_STATES = 133
EVIDENCE = {}
CHECKS = []


def digest(path):
    with path.open("rb") as stream:
        return sha256(stream.read()).hexdigest()


def read(path):
    EVIDENCE[path] = digest(path)
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write(path, records):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def index(rows, *fields):
    keyed = {tuple(row[key] for key in fields): row for row in rows}
    assert len(keyed) == len(rows), f"Duplicate keys in {fields}"
    return keyed


def check(name, result, observed, requirement):
    CHECKS.append(dict(check=name, observed=observed, requirement=requirement, status="PASS" if result else "FAIL"))
    assert result, f"{name}: observed {observed}; required {requirement}"


def canon(value):
    return value.replace("SUM159_NLS", "SUM-159_NLS")


def verify_sources():
    for directory in (FREEZE, REVIEW):
        receipt = directory / "SHA256SUMS"
        EVIDENCE[receipt] = digest(receipt)
        entries = [line.split(maxsplit=1) for line in receipt.read_text().splitlines()]
        mismatch = [name for expected, name in entries if digest(directory / name) != expected]
        check(directory.name + "_archive_hashes", not mismatch, len(entries) - len(mismatch), len(entries))
    source_rows = read(FREEZE / "source_manifest.tsv")
    mismatches = []
    for row in source_rows:
        locator, expected = row["locator"], row["sha256"]
        if " at " in locator:
            repo, revision = locator.split(" at ", 1)
            data = subprocess.check_output(["git", "-C", repo, "show", revision])
            actual = sha256(data).hexdigest()
        else:
            path = Path(locator)
            actual = digest(path)
            EVIDENCE[path] = actual
        if actual != expected:
            mismatches.append(locator)
    check("original_source_receipts", not mismatches, len(source_rows) - len(mismatches), len(source_rows))
    for row in read(FIT / "joint_source_hashes.tsv"):
        path = Path(row["path"])
        EVIDENCE[path] = digest(path)
        assert EVIDENCE[path] == row["sha256"], path


def specimen_tables(out):
    specimens = read(FREEZE / "observation_manifest.tsv")
    index(specimens, "observation_key")
    saved = [r for r in read(FIT / "readout_predictions.tsv") if r["mode"] == "primary"]
    prediction_index = index([dict(r, event_id=canon(r["event_id"])) for r in saved], "readout", "event_id")
    seed = set()
    for branch in read(FIT / "a7k_event_contract.tsv"):
        seed.add(branch["seed"])
    samples = []
    for row in specimens:
        readout, role = row["readout"], row["analysis_role"]
        event = row["confirmed_event_id"]
        key = (readout, canon(event))
        prediction = prediction_index.get(key)
        if prediction:
            assert row["timestamp_as_recorded"] and row["passage_id"] == prediction["passage_id"]
            assert isclose(float(row["nominal_o2_pct"]), float(prediction["oxygen_pct"]))
            assert int(row["calendar_day_since_seed"]) == float(prediction["model_age_days"])
            if readout == "death_fraction":
                assert row["numerator"] == prediction["numerator"] and row["denominator"] == prediction["denominator"]
                assert int(row["denominator"]) > 0 and 0 <= int(row["numerator"]) <= int(row["denominator"])
            if readout == "corrected_count":
                assert float(row["measurement_value"]) == float(prediction["observed"])
            if readout == "karyotype":
                assert row["sample_size"] == prediction["sample_size"]
            decision = "joint_output"
        elif role == "seeding_conditioning" or (role == "branch_seeding_conditioning" and event in seed):
            decision = "conditioned_seed_record"
        elif role == "initialization_proxy_only":
            decision = "conditioned_root_proxy_in_B_fixed"
        elif role.startswith("held_") or role == "partial_classification_sensitivity":
            decision = "held_pending_evidence"
        else:
            decision = "outside_registered_graph"
        samples.append(dict(row, decision=decision,
                            model_day=prediction["model_age_days"] if prediction else "",
                            admission_evidence="Existing source and experimental observer cross-check" if prediction else row["evidence_note"]))
    outputs = [r for r in samples if r["decision"] == "joint_output"]
    n_assays = Counter(r["readout"] for r in outputs)
    check("unique_joint_observations", len(outputs) == len(saved) == 426, len(outputs), 426)
    check("joint_output_counts", n_assays == dict(corrected_count=225, death_fraction=195, karyotype=6),
          "/".join(str(n_assays[k]) for k in ("corrected_count", "death_fraction", "karyotype")), "225/195/6")
    check("known_seed_records", sum(r["decision"] == "conditioned_seed_record" for r in samples) == 118, 118, 118)
    event_counts = Counter(canon(r["confirmed_event_id"]) for r in outputs)
    overlap = dict(sorted(Counter(event_counts.values()).items()))
    check("event_overlap", overlap == {1: 30, 2: 189, 3: 6}, overlap, "30 one-readout, 189 two-readout and 6 three-readout events")
    death_rows = [r for r in specimens if r["readout"] == "death_fraction"]
    definitions = sum(r["count_definition_evidence"].startswith("whole_flask_corrected_count_including") for r in death_rows)
    partitions = Counter(r["classifier_object_partition_verified"] for r in death_rows)
    check("image_definition_coverage", definitions == 98 and partitions["yes"] == 46,
          f"{definitions} count definitions and {dict(partitions)} partitions", "98 definitions and 46 verified partitions")
    matched = defaultdict(Counter)
    for r in samples:
        matched[r["passage_id"]][r["decision"], r["readout"]] += 1
    coverage = []
    for passage, counts in sorted(matched.items()):
        rows = [r for r in samples if r["passage_id"] == passage]
        detected_events = {canon(r["confirmed_event_id"]) for r in rows if r["decision"] == "joint_output"}
        coverage.append(dict(passage_id=passage, cohort=rows[0]["cohort"], lineage=rows[0]["lineage"],
                             observed_count_events=counts["joint_output", "corrected_count"],
                             observed_death_events=counts["joint_output", "death_fraction"],
                             observed_karyotypes=counts["joint_output", "karyotype"],
                             karyotyped_cells=sum(int(r["sample_size"]) for r in rows
                                                  if r["readout"] == "karyotype" and r["decision"] == "joint_output"),
                             same_event_triplets=sum(event_counts[event] == 3 for event in detected_events),
                             held_flow=counts["held_pending_evidence", "flow_density"],
                             held_karyotypes=counts["held_pending_evidence", "karyotype"],
                             conditioned_seeds=counts["conditioned_seed_record", "corrected_count"]))
    if FCS.is_file():
        EVIDENCE[FCS] = digest(FCS)
    write(out / "specimen_decisions.tsv", samples)
    write(out / "coverage_by_passage.tsv", coverage)
    write(out / "open_evidence.tsv", [dict(r, disposition="open") for r in read(FREEZE / "missing_information.tsv")])
    write(out / "event_mask.tsv", [dict(observation_key=r["observation_key"], modality=r["readout"],
                                 event=canon(r["confirmed_event_id"]), timestamp=r["timestamp_as_recorded"],
                                 passage=r["passage_id"], cohort=r["cohort"], oxygen_pct=r["nominal_o2_pct"],
                                 model_day=r["model_day"], sample_size=r["sample_size"], denominator=r["denominator"])
                             for r in outputs])
    return samples, outputs, saved


def graph_table(samples, out):
    graph = read(INHERIT / "passage_design.tsv")
    old_graph = index(read(FIXED / "active_passages.tsv"), "passage_id")
    crosswalk = index(read(STRUCT / "20261118_culture_experimental_OLS_step01/count_event_mapping.tsv"), "record_id")
    cut = {key[0] for key, row in old_graph.items() if not row["retained_parent_id"]}
    certified_cuts = {r["child_passage_id"] for r in read(FIXED / "cut_edges.tsv")}
    assert cut == certified_cuts and len(cut) == 14
    rows = []
    for p in graph:
        key, parent = p["passage_id"], p["parent_id"]
        old_parent = old_graph[key,]["retained_parent_id"]
        assert old_parent == ("" if key in cut else parent)
        if old_parent:
            assert int(old_graph[key,]["parent_donor_step"]) == int(p["donor_step"])
        assert float(old_graph[key,]["seed_corrected_count"]) == float(p["fixed_live_cells"])
        donor_record = ""
        donor_timestamp = ""
        if parent:
            source = crosswalk[key,]
            assert source["physical_parent_passage_id"] == parent
            donor_record = source["passage_seed_donor_record_id"]
            donor_timestamp = source["passage_seed_donor_timestamp_as_recorded"]
            assert donor_record and donor_timestamp
            if (parent,) in crosswalk:
                parent_seed = crosswalk[parent,]["seed_timestamp_as_recorded"]
                donor_day = (datetime.fromisoformat(donor_timestamp).date() - datetime.fromisoformat(parent_seed).date()).days
                assert donor_day * 20 == int(p["donor_step"])
        rows.append(dict(passage_id=key, cohort=p["cohort"], oxygen_pct=p["oxygen_pct"],
                         fixed_seed_live_objects=p["fixed_live_cells"], donor_passage=parent,
                         donor_record=donor_record, donor_timestamp=donor_timestamp,
                         donor_model_day=int(p["donor_step"]) * STEP if parent else "",
                         mixture_rule_free_mixtures="free_simplex" if key in cut else "inherit_donor",
                         mixture_rule_inheritance="inherit_donor" if parent else "condition_on_proxy_or_use_unknown_gaussian",
                         historical_supply="original main-chain bookkeeping"))
    counts = index([r for r in samples if r["readout"] == "corrected_count"], "confirmed_event_id")
    for b in read(FIT / "a7k_event_contract.tsv"):
        donor = counts[b["donor"],]
        seed = counts[b["seed"],]
        harvest = counts[b["harvest"],]
        delta = (datetime.fromisoformat(b["harvest_timestamp"]) - datetime.fromisoformat(b["seed_timestamp"]))
        assert isclose(delta.total_seconds() / 86400, float(b["elapsed_age_days"]), abs_tol=5.1e-9)
        assert seed["timestamp_as_recorded"] == b["seed_timestamp"]
        assert harvest["timestamp_as_recorded"] == b["harvest_timestamp"]
        assert donor["timestamp_as_recorded"] == b["donor_timestamp"]
        assert float(donor["nominal_o2_pct"]) == 0 and int(donor["calendar_day_since_seed"]) == 5
        assert float(seed["measurement_value"]) == float(b["seed_count"])
        rows.append(dict(passage_id=seed["passage_id"], cohort=seed["cohort"], oxygen_pct=b["nominal_oxygen_pct"],
                         fixed_seed_live_objects=b["seed_count"], donor_passage=donor["passage_id"],
                         donor_record=b["donor"], donor_timestamp=b["donor_timestamp"],
                         donor_model_day=5, mixture_rule_free_mixtures="inherit_donor",
                         mixture_rule_inheritance="inherit_donor",
                         historical_supply="stock transfer with separate reset and no added historical demand"))
    check("graph_segments_and_free_roots", len(rows) == 120, f"{len(rows)}/{len(cut)}", "120/14")
    indices = {r["passage_id"]: i for i, r in enumerate(rows)}
    assert all(not r["donor_passage"] or indices[r["donor_passage"]] < i for i, r in enumerate(rows))
    write(out / "passage_graph.tsv", rows)
    return rows


def replay(saved, outputs, out):
    os.environ["JOINT_KERNEL"] = "1"
    sys.path.insert(0, str(FIT))
    from common import load_x
    from joint import JointConfig, JointModel

    endpoint = read(FIT / "selected_primary.tsv")[0]
    EVIDENCE[FIT / "joint.dylib"] = digest(FIT / "joint.dylib")
    for name in ("common.py", "joint.py"):
        EVIDENCE[FIT / name] = digest(FIT / name)
    model = JointModel()
    config = JointConfig(model)
    p, _ = config.unmap(load_x(config, endpoint))
    assert p[14] == p[16] and p[15] == p[17]
    channels, live_states = model.readouts(p)
    margin, _ = model.supply(p)
    stock = read(FIT / "seed_composition.tsv")
    saved_stock = {r["branch"]: r for r in stock if r["mode"] == "primary"}
    tv, _, _, seed_sizes = model.composition(p, details=True)
    saved_supply = [float(r["value"]) for r in read(FIT / "constraint_margins.tsv") if r["case"] == "selected_primary"]
    assert len(margin) == len(saved_supply) == 110
    max_gap_error = max(abs(float(a) - b) for a, b in zip(margin, saved_supply))
    check("original_supply_replay", min(margin) >= 0 and max_gap_error < 1e-10,
          f"minimum {min(margin):.9g}; max gap error {max_gap_error:.3g}", "110 nonnegative margins and max gap error below 1e-10")
    seed_error = max(abs(float(saved_stock[b["branch"]]["seed_count"]) - n)
                     for b, n in zip(read(FIT / "a7k_event_contract.tsv"), seed_sizes))
    check("a7k_reseed", len(tv) == 4 and max(tv) < 1e-12 and seed_error < 1e-6,
          f"max TV {max(tv):.3g}; maximum seed error {seed_error:.3g}", "Four seed totals and strict donor fractions")
    by_prediction = index([dict(r, event_id=canon(r["event_id"])) for r in saved], "readout", "event_id")
    c_errors, d_errors, signal = [], [], []
    for slot in model.indices:
        event = canon(model.slot_ids[slot])
        row = by_prediction["corrected_count", event]
        l, h, m = map(float, channels[:, slot])
        expected = (l + h + m) * 1e12
        c_errors.append(abs(expected / float(row["predicted"]) - 1))
        signal.append(dict(readout="visible_count", event=canon(event), predicted=expected, recorded=row["observed"]))
    a, loading, q = float(p[55]), float(p[56]), float(p[57])
    eps = 1e-12
    for slot in model.death_rows:
        event = canon(model.slot_ids[slot])
        row = by_prediction["death_fraction", event]
        l, h, m = map(float, channels[:, slot])
        base = a + log((h + q*m + eps) / (l + (1-q)*m + eps))
        residual = log(float(by_prediction["corrected_count", event]["observed"])) - log(max(1, (l+h+m)*1e12))
        probability = (1-2*eps) / (1+exp(-(base + loading*residual))) + eps
        d_errors.append(abs(probability - float(row["predicted"])))
        signal.append(dict(readout="death_probability_zero_residual", event=event,
                           predicted=(1-2*eps) / (1+exp(-base)) + eps, recorded=row["observed"]))
    check("count_observer_replay", max(c_errors) < 1e-10, max(c_errors), "maximum relative discrepancy below 1e-10")
    check("death_observer_replay", max(d_errors) < 1e-10, max(d_errors), "maximum absolute discrepancy below 1e-10")
    chromosomes = read(ALIGN / "cell_alignment.tsv")
    index(chromosomes, "assay_origin_id", "cell_index")
    by_cell = Counter((r["assay_origin_id"], int(r["observed_N_22_autosomes"])) for r in chromosomes)
    by_bin = index([r for r in read(FIT / "karyotype_distributions.tsv") if r["mode"] == "primary"], "event_id", "autosomes")
    k_errors = []
    k_total = 0
    for event, state in zip(model.all_kary, live_states):
        assert by_prediction["karyotype", canon(event)]["event_id"] == canon(event)
        normalization = sum(state)
        assert normalization > 0
        total = 0
        for n, value in enumerate(state, start=22):
            record = by_bin[event, str(n)]
            assert int(record["observed_cells"]) == by_cell[event, n]
            total += by_cell[event, n]
            k_errors.append(abs((1-eps) * value/normalization + eps/N_STATES - float(record["fitted_probability"])))
        assert total == int(by_prediction["karyotype", canon(event)]["sample_size"])
        k_total += total
        prob_sum = sum((1-eps) * value/normalization + eps/N_STATES for value in state)
        assert abs(prob_sum - 1.) < 1e-12
    check("chromosome_histogram_replay", len(by_bin) == 798 and k_total == 104 and max(k_errors) < 1e-12,
          f"798 bins, {k_total} cells, max probability error {max(k_errors):.3g}", "Six complete normalized histograms")
    event_modalities = defaultdict(set)
    for row in outputs:
        event_modalities[canon(row["confirmed_event_id"])].add(row["readout"])
    branches = read(FIT / "a7k_event_contract.tsv")
    n_joint = sum(len(event_modalities[canon(b["harvest"])]) == 3 for b in branches)
    check("a7k_three_readout_joins", n_joint == 4, n_joint, 4)
    write(out / "saved_endpoint_signal_replay.tsv", signal)


def cases(graph, out):
    mech_rows = [r for r in read(PROFILE / "parameter_and_nuisance_map.tsv")
                 if r["contract"] == "J" and r["condition"] == "unrestricted"]
    history = index(read(CONTRACT / "parameter_roles.tsv"), "parameter")
    mechanistic = [r for r in mech_rows if (r["parameter"],) in history
                   and history[r["parameter"],]["scientific_role"] == "mechanism"]
    assert len(mechanistic) == 14
    # One chromosome coordinate per sample is redundant because the probabilities sum to one.
    masks = {"C": (225, 0, 0), "CD": (225, 195, 0), "CK": (225, 0, 792),
             "DK": (0, 195, 792), "CDK": (225, 195, 792)}
    variants = [
        ("free_mixtures", "primary", 14 * 132, 0, "fixed", "calendar_day_4"),
        ("fixed_proxy_inheritance", "co_primary", 0, 0, "fixed", "calendar_day_4"),
        ("gaussian_root_inheritance", "root_sensitivity", 4, 0, "fixed", "calendar_day_4"),
        ("a7k_tv_inheritance", "composition_sensitivity", 0, 4, "fixed", "calendar_day_4"),
        ("elapsed_inheritance", "timing_sensitivity", 0, 0, "fixed", "elapsed_A7K_time"),
        ("elapsed_free_mixtures", "timing_sensitivity", 14 * 132, 0, "fixed", "elapsed_A7K_time"),
        ("free_mixtures_free_viable_seed", "seed_sensitivity", 14 * 132, 0, "118_vary_with_C0_conditioned", "calendar_day_4"),
    ]
    registered, roles, contrasts = [], [], []
    law_coordinates = ["count_sigma", "rho", "shared_count_loading", "log10_concentration"]
    old_scales = ["sigma_growth", "sigma_kary"]
    initials = ["init_mean_2N", "init_sd_2N", "init_mean_4N", "init_sd_4N"]
    for variant, priority, root_dims, tilt_dims, seeds, timing in variants:
        chosen_masks = list(masks) if priority in ("primary", "co_primary", "root_sensitivity") else ["CDK"]
        for mask in chosen_masks:
            c, d, k = masks[mask]
            label = f"{variant}_{mask}"
            columns = 14 + 2 + 2 * bool(d) + root_dims + tilt_dims + (118 if seeds != "fixed" else 0)
            registered.append(dict(case=label, priority=priority, mask=mask, count_rows=c, death_rows=d,
                                   independent_karyotype_rows=k, total_output_rows=c+d+k, mechanism_coordinates=14,
                                   root_coordinates=root_dims, a7k_tilt_coordinates=tilt_dims,
                                   free_viable_seed_coordinates=118 if seeds != "fixed" else 0,
                                   deterministic_signal_coordinates=columns, transfer_calendar=timing,
                                   proof_status="registered_pending_validation"))
            def role(name, status, dimension, bound, binding):
                roles.append(dict(case=label, coordinate=name, role=status, dimensions=dimension, domain=bound,
                                  source_or_binding=binding))
            for row in mechanistic:
                name = row["parameter"]
                role(name, "mechanism", 1, f"[{row['lower']}, {row['upper']}]", f"native index {row['native_index']}")
            role("eta_H_2N", "signal", 1, "[0, 1]", "native indices 14 and 16; shared eta_H")
            role("eta_M_2N", "signal", 1, "[0, 1]", "native indices 15 and 17; shared eta_M")
            role("class_log_odds", "signal" if d else "inactive", int(bool(d)), "[-4, 4]", "death observer, native index 55")
            role("misseg_dead_share", "signal" if d else "inactive", int(bool(d)), "[0, 1]", "death observer, native index 57")
            for name in initials:
                original = history[name,]
                status = "initial_state_unknown" if root_dims == 4 else "replaced_by_simplex" if root_dims > 4 else "conditioned_proxy"
                role(name, status, int(root_dims == 4), f"[{original['lower']}, {original['upper']}]", "historical root Gaussian family")
            for name in old_scales:
                original = history[name,]
                role(name, "historical_observer_extension", 0, f"[{original['lower']}, {original['upper']}]", "historical culture observer")
            for name in law_coordinates:
                applicable = (bool(c) if name in ("count_sigma", "rho") else bool(d))
                if name == "shared_count_loading":
                    applicable = bool(c and d)
                role(name, "law_only_extension" if applicable else "fixed_or_absent", 0,
                     {"count_sigma": "(0, infinity)", "rho": "[0, 0.99]", "shared_count_loading": "[-4, 4]",
                      "log10_concentration": "[0, 6]"}[name], "count/conditional classification law")
            if root_dims > 4:
                for r in graph:
                    if r["mixture_rule_free_mixtures"] == "free_simplex":
                        role(f"pi:{r['passage_id']}", "initial_state_unknown", 132,
                             "positive probability simplex; omit n=154", "frozen passage graph")
            if tilt_dims:
                for cohort in ("2N", "4N"):
                    for axis in ("mean", "spread"):
                        role(f"A7K_{cohort}_{axis}_tilt", "initial_state_unknown", 1, "[-32, 32], branch TV <=0.05",
                             "Step 3 TV cohort-shared tilt")
            for r in graph:
                if not r["donor_passage"]:
                    # The two auxiliary upstream seed totals are fixed in all registered comparisons.
                    continue
                if seeds != "fixed":
                    role(f"seed_live:{r['passage_id']}", "initial_state_unknown", 1, "0<L0<C0",
                         "Conditioned C0; allocate initial visible remainder equally between visible pools")
            dims = sum(r["dimensions"] for r in roles if r["case"] == label)
            assert dims == columns, (label, dims, columns)
        if len(chosen_masks) == 5:
            for smaller, larger, stream in (("C", "CD", "D"), ("C", "CK", "K"), ("CD", "CDK", "K"),
                                             ("CK", "CDK", "D"), ("DK", "CDK", "C")):
                contrasts.append(dict(smaller=f"{variant}_{smaller}", larger=f"{variant}_{larger}", changed_readout=stream,
                                      signal_rule="Use the same graph and common coordinates; inactive observer columns are zero.",
                                      law_rule="For count deletion, zero loading on both sides. Full-loading extension needs separate proof."))
        elif variant != "fixed_proxy_inheritance":
            reference = "free_mixtures" if "free" in variant else "fixed_proxy_inheritance"
            contrasts.append(dict(smaller=f"{reference}_CDK", larger=f"{variant}_CDK", changed_readout="none",
                                  signal_rule="Paired initial-state or timing sensitivity; no readout effect inferred.",
                                  law_rule="Retain the same readout and conditioning event IDs."))
    write(out / "case_register.tsv", registered)
    write(out / "parameter_roles_by_case.tsv", roles)
    write(out / "contrast_register.tsv", contrasts)


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    verify_sources()
    specimens, outputs, saved_predictions = specimen_tables(out)
    passage_graph = graph_table(specimens, out)
    replay(saved_predictions, outputs, out)
    cases(passage_graph, out)
    write(out / "interface_replay.tsv", CHECKS)
    write(out / "source_receipts.tsv", [dict(path=str(path), sha256=hash_value) for path, hash_value in sorted(EVIDENCE.items())])
    import numpy
    import scipy
    write(out / "environment.tsv", [dict(component=name, value=value) for name, value in
                              (("python", sys.version.split()[0]), ("python_executable", sys.executable),
                               ("numpy", numpy.__version__), ("scipy", scipy.__version__),
                               ("git_head", subprocess.check_output(["git", "-C", ROOT, "rev-parse", "HEAD"], text=True).strip()))])
    files = sorted(out.glob("*.tsv"))
    (out / "SHA256SUMS").write_text("".join(f"{digest(f)}  {f.name}\n" for f in files))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=output_directory, required=True)
    main(parser.parse_args().out)
