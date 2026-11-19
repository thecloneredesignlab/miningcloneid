#!/usr/bin/env python3
"""Audit the experimental count-to-segment crosswalk from existing snapshots."""

import argparse
import csv
from collections import Counter, defaultdict
from datetime import datetime
from decimal import Decimal
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
RESULTS = ROOT / "oxygen/results"
CONTRACT = RESULTS / "synthetic_identifiability/20261113_culture_contract_steps01_03"
SCHEDULE = RESULTS / "structural_identifiability/20261117_culture_actual_schedule_step01"
ARCHIVE = RESULTS / "structural_identifiability/20261109_seed144_wp0/source/seed144"
DEFAULT_SNAPSHOT = ROOT.parent / (
    "soft_couping_org/data/InVitroData/cloneid_passaging_sum159_snapshot_20260731.tsv"
)
DEFAULT_DEAD_ANNOTATIONS = DEFAULT_SNAPSHOT.with_name("sum159_dead_cell_db_enriched_20260731.tsv")
DEFAULT_PREPARED_COUNTS = ROOT.parents[1] / (
    "CloneID/Deployments/sum159-red-pilot-diagnostic-001/operator_project/"
    "calibration/prepared_inputs/count_growth_observations.tsv"
)


def read(path, delimiter="\t"):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter=delimiter))


def write(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def require(condition, context, detail):
    if not condition:
        raise ValueError(f"{context}: {detail}")


def key_by(rows, *columns):
    lookup = {}
    for row in rows:
        key = tuple(row[column] for column in columns)
        require(key not in lookup, key, f"duplicate key for {columns}")
        lookup[key] = row
    return lookup


def grouped(rows, column):
    groups = defaultdict(list)
    for row in rows:
        groups[row[column]].append(row)
    return groups


def value(raw):
    if not raw or raw in {"NA", "NaN", "NULL"}:
        return None
    return Decimal(raw)


def segment_key(row):
    return (row["cohort"], row["segment_id"])


def numbers(raw):
    return [value(item) for item in raw.split(";")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--prepared-counts", type=Path, default=DEFAULT_PREPARED_COUNTS)
    parser.add_argument("--dead-annotations", type=Path, default=DEFAULT_DEAD_ANNOTATIONS)
    parser.add_argument("--out-dir", type=Path, default=HERE)
    args = parser.parse_args()

    manifest = read(CONTRACT / "growth_raw_count_manifest.tsv")
    manifest_ids = key_by(manifest, "record_id")
    manifest_groups = grouped(manifest, "passage_id")
    metadata = key_by(read(ROOT / "oxygen/data/metadata.csv", ","), "id")
    records = read(args.snapshot)
    snapshot = key_by(records, "id")
    raw_children = grouped(records, "passaged_from_id1")
    all_child_seeds = defaultdict(list)
    for seed in records:
        if seed["event"] == "seeding":
            for parent_id in {seed["passaged_from_id1"], seed["passaged_from_id2"]} - {""}:
                all_child_seeds[parent_id].append(seed)

    passages = read(SCHEDULE / "passage_crosswalk.tsv")
    passage_by_id = key_by(passages, "passage_id")
    passage_children = grouped(passages, "direct_parent_record_id")
    ols_passages = {
        p["passage_id"]: p for p in passages if p["recorded_ols_growth_per_day"] != "NA"
    }
    require(ols_passages.keys() == manifest_groups.keys(), "scope", "passage sets differ")
    observations = key_by(
        read(ROOT / "data/InVitroData_LTEE/invitro_passage_observations.tsv"),
        "passage_id",
    )
    schedule = key_by(read(CONTRACT / "experimental_schedule.tsv"), "cohort", "segment_id")
    lineage = key_by(read(ARCHIVE / "invitro_lineage_summary.tsv"), "passage_id")
    archive_rows = read(ARCHIVE / "invitro_daily_counts.tsv")
    key_by(archive_rows, "cohort", "segment_id", "day")
    daily = {(*segment_key(row), value(row["day"])): row for row in archive_rows}
    require(len(daily) == len(archive_rows), "archive", "numeric day keys are ambiguous")
    dt = value(key_by(read(ARCHIVE / "fit_summary.tsv"), "metric")[("dt",)]["value"])
    dt_seconds = dt * 86400
    prepared = key_by(
        [r for r in read(args.prepared_counts) if (r["event_id"],) in manifest_ids],
        "event_id",
    )
    dead_annotations = key_by(
        [r for r in read(args.dead_annotations) if (r["db_source_id"],) in manifest_ids],
        "db_source_id",
    )

    rows = []
    segments = set()
    for passage_id, counts in manifest_groups.items():
        passage = ols_passages[passage_id]
        k = (passage["cohort"], passage["production_segment_id"])
        sched = schedule[k]
        saved = lineage[(passage_id,)]
        obs = observations[(passage_id,)]
        require(k not in segments, passage_id, "OLS groups share a production segment")
        segments.add(k)
        require(
            passage_id in sched["data_ids"].split(";")
            and segment_key(saved) == k
            and segment_key(obs) == (k[0], k[1].partition("::")[0])
            and passage["production_parent_segment_id"] == sched["parent_segment_id"]
            == saved["parent_segment_id"]
            and sched["parent_segment_id"].partition("::")[0] == (obs["parent_segment_id"] or "NA"),
            passage_id, "production lineage does not match",
        )
        duration = value(sched["nominal_duration_days"])
        oxygen = value(sched["oxygen_pct"])
        selected_day = value(sched["selected_day"])
        require(
            duration == value(saved["duration_days"])
            == value(passage["production_nominal_duration_days"])
            == value(obs["passage_duration"])
            and oxygen == value(saved["oxygen_pct"])
            == value(obs["oxygen_pct"]) == value(passage["assigned_oxygen_pct"])
            == value(passage["oxygen_pct_from_count_metadata_media_join"])
            and selected_day == value(saved["selected_day"])
            == value(passage["seed144_selected_local_day"])
            and passage["production_uses_fallback_duration"] == "FALSE",
            passage_id, "duration, oxygen, or selected day mismatch",
        )
        count_ids = {row["record_id"] for row in counts}
        seed = snapshot[(passage_id,)]
        seed_stamp = datetime.fromisoformat(seed["date"])
        seed_count = value(seed["correctedCount"])
        donor = snapshot[(seed["passaged_from_id1"],)]
        donor_passage = passage_by_id[(passage["physical_parent_passage_id"],)]
        transfer_seconds = int((seed_stamp - datetime.fromisoformat(donor["date"])).total_seconds())
        require(
            seed["event"] == "seeding" and seed["date"] == passage["seeding_time"]
            and seed["passaged_from_id1"] == passage["direct_parent_record_id"]
            and donor["event"] == "harvest" and donor["date"] == passage["direct_parent_record_time"]
            and donor["passaged_from_id1"] == passage["physical_parent_passage_id"]
            and donor_passage["production_segment_id"] == sched["parent_segment_id"]
            and donor_passage["cohort"] == k[0] and transfer_seconds >= 0
            and seed_count == value(passage["seed_corrected_count"])
            == value(obs["initial_cells"]) == value(sched["initial_cells"])
            == value(saved["initial_cells"])
            and abs(seed_count - value(daily[(*k, Decimal(0))]["live_cells"])) < Decimal("0.5"),
            passage_id, "seeding record or day-zero model total mismatch",
        )
        source_ids = {passage_id}
        source_ids.update(r["id"] for r in raw_children[passage_id]
                          if value(r["correctedCount"]) is not None)
        require(count_ids == source_ids, passage_id, "source count membership differs")
        counts.sort(key=lambda row: (value(row["relative_day"]), row["record_id"]))
        require(
            [c["record_id"] for c in counts] == passage["ols_count_record_ids"].split(";")
            and [value(c["relative_day"]) for c in counts]
            == numbers(passage["ols_count_calendar_days"])
            and len({c["relative_day"] for c in counts}) == len(counts),
            passage_id, "prior calendar-day crosswalk or within-passage order differs",
        )
        last_id = counts[-1]["record_id"]
        last_count = value(counts[-1]["corrected_count"])
        require(
            last_count == value(obs["final_cells"]) == value(sched["final_cells"])
            == value(saved["final_cells"]) == value(saved["target_live_cells"])
            == value(sched["target_live_cells"])
            and value(counts[-1]["relative_day"]) == duration,
            passage_id, "last count does not match model endpoint target input",
        )
        schedule_days = set(numbers(sched["candidate_days"]))
        for i, count in enumerate(counts):
            record_id = count["record_id"]
            raw = snapshot[(record_id,)]
            stamp = datetime.fromisoformat(raw["date"])
            raw_age_seconds = int((stamp - seed_stamp).total_seconds())
            calendar_day = (stamp.date() - seed_stamp.date()).days
            corrected = value(raw["correctedCount"])
            archive = daily[(*k, Decimal(calendar_day))]
            source_is_metadata = count["source"] == "oxygen/data/metadata.csv"
            require(
                stamp.date().isoformat() == count["date"] and raw_age_seconds >= 0
                and calendar_day == value(count["relative_day"])
                and raw["date"] == passage["ols_count_timestamps"].split(";")[i]
                and calendar_day in schedule_days and calendar_day <= duration
                and value(archive["selected_day"]) == selected_day
                and value(archive["target_live_cells"]) == last_count,
                record_id, "date, calendar day, or archived model day/target mismatch",
            )
            require(
                corrected is not None and corrected > 0
                and corrected == value(count["corrected_count"])
                and raw["passaged_from_id1"] == count["physical_parent"]
                and raw["event"] == ("seeding" if record_id == passage_id else "harvest"),
                record_id, "corrected count, event, or parent mismatch",
            )
            if source_is_metadata:
                meta = metadata[(record_id,)]
                require(
                    meta["passage_id"] == passage_id and meta["date"] == count["date"]
                    and meta["passaged_from_id1"] == raw["passaged_from_id1"]
                    and value(meta["correctedCount"]) == corrected
                    and value(meta["num_date"]) == calendar_day and value(meta["oxygen_pct"]) == oxygen,
                    record_id, "metadata export does not agree",
                )
            else:
                require(
                    Path(count["source"]).name == args.snapshot.name
                    and (record_id,) not in metadata,
                    record_id, "snapshot-only provenance does not agree",
                )
            semantics = ""
            if (record_id,) in prepared:
                annotation = prepared[(record_id,)]
                require(
                    annotation["event_timestamp_utc"][:19].replace("T", " ") == raw["date"]
                    and value(annotation["corrected_cell_count"]) == corrected,
                    record_id, "independent prepared-input timestamp/count differs",
                )
                semantics = annotation["count_semantics"]

            dead = dead_annotations.get((record_id,))
            if dead is not None:
                dead_count = value(dead["dead_count"])
                denominator = value(dead["eligible_denominator"])
                require(
                    dead["db_event_datetime"] == raw["date"]
                    and dead["db_event"] == raw["event"]
                    and dead["cohort"] == k[0]
                    and dead["cell_state"] == "dead_cell"
                    and denominator > 0 and 0 <= dead_count <= denominator
                    and abs(value(dead["observed_dead_fraction"]) - dead_count / denominator)
                    <= Decimal("0.00000001"),
                    record_id, "dead-object annotation or denominator differs",
                )

            children = sorted(passage_children[record_id], key=lambda p: p["seeding_time"])
            for child in children:
                child_seed = snapshot[(child["passage_id"],)]
                require(
                    raw["event"] == "harvest" and child_seed["event"] == "seeding"
                    and child_seed["passaged_from_id1"] == record_id
                    and child_seed["date"] == child["seeding_time"]
                    and child["direct_parent_record_time"] == raw["date"]
                    and datetime.fromisoformat(child_seed["date"]) >= stamp,
                    record_id, "documented harvest to child-seed link differs",
                )
            raw_seeds = sorted(all_child_seeds[record_id], key=lambda child: (child["date"], child["id"]))
            require(
                {p["passage_id"] for p in children}
                == {child["id"] for child in raw_seeds if (child["id"],) in passage_by_id}
                and all(child["date"] >= raw["date"] for child in raw_seeds),
                record_id, "full snapshot and crosswalk child-seed links differ",
            )
            count_role = "seed" if i == 0 else "last_count" if record_id == last_id else "intermediate_count"
            relative_to_end = (
                "after" if calendar_day > selected_day
                else "at" if calendar_day == selected_day else "before"
            )
            rows.append({
                "record_id": record_id,
                "passage_id": passage_id,
                "cohort": k[0],
                "oxygen_pct": str(oxygen),
                "model_segment_id": k[1],
                "model_parent_segment_id": sched["parent_segment_id"],
                "input_segment_id_before_branch_split": obs["segment_id"],
                "physical_parent_passage_id": passage["physical_parent_passage_id"],
                "first_parent_record_id": raw["passaged_from_id1"],
                "second_parent_record_id": raw["passaged_from_id2"],
                "passage_seed_donor_record_id": donor["id"],
                "passage_seed_donor_timestamp_as_recorded": donor["date"],
                "donor_to_seed_elapsed_seconds_as_recorded": transfer_seconds,
                "snapshot_event": raw["event"],
                "within_passage_count_role": count_role,
                "seed_timestamp_as_recorded": seed["date"],
                "count_timestamp_as_recorded": raw["date"],
                "ols_calendar_day": calendar_day,
                "archived_model_day": archive["day"],
                "age_seconds_from_recorded_timestamps": raw_age_seconds,
                "age_days_from_recorded_timestamps": f"{Decimal(raw_age_seconds) / 86400:.9f}",
                "timestamp_age_on_model_dt_grid": int(raw_age_seconds % dt_seconds == 0),
                "model_nominal_duration_days": str(duration),
                "timestamp_age_exceeds_model_duration": int(raw_age_seconds > duration * 86400),
                "seed144_selected_model_day": str(selected_day),
                "ols_day_vs_seed144_selected_day": relative_to_end,
                "seed144_incoming_parent_reseed_day": passage["seed144_parent_reseed_local_day"],
                "linked_child_seed_ids_in_131": ";".join(p["passage_id"] for p in children),
                "linked_child_seed_timestamps_in_131": ";".join(p["seeding_time"] for p in children),
                "linked_child_seed_ids_in_snapshot": ";".join(child["id"] for child in raw_seeds),
                "linked_child_seed_timestamps_in_snapshot": ";".join(child["date"] for child in raw_seeds),
                "db_cellCount": raw["cellCount"],
                "db_correctedCount": raw["correctedCount"],
                "db_Countess": raw["Countess"],
                "db_growthType": raw["growthType"],
                "db_comment": raw["comment"],
                "passage_seed_correctedCount": seed["correctedCount"],
                "production_initial_live_cells": sched["initial_cells"],
                "production_target_live_cells": sched["target_live_cells"],
                "manifest_source": count["source"],
                "sum159_pilot_count_semantics_if_matched": semantics,
                "morphology_dead_objects_if_matched": dead["dead_count"] if dead else "",
                "morphology_eligible_objects_if_matched": dead["eligible_denominator"] if dead else "",
                "morphology_classification_status_if_matched": dead["classification_status"] if dead else "",
            })

    categories = Counter((r["cohort"], r["snapshot_event"]) for r in rows)
    end_relation = Counter(r["ols_day_vs_seed144_selected_day"] for r in rows)
    count_source = Counter(r["manifest_source"] == "oxygen/data/metadata.csv" for r in rows)
    after_end = [r for r in rows if r["ols_day_vs_seed144_selected_day"] == "after"]
    beyond_duration = [r for r in rows if r["timestamp_age_exceeds_model_duration"]]
    linked_counts = [r for r in rows if r["linked_child_seed_ids_in_131"]]
    full_snapshot_links = [r for r in rows if r["linked_child_seed_ids_in_snapshot"]]
    positive_dead = [r for r in rows if (value(r["morphology_dead_objects_if_matched"]) or 0) > 0]
    last_counts = [r for r in rows if r["within_passage_count_role"] == "last_count"]
    seeds = [r for r in rows if r["within_passage_count_role"] == "seed"]
    maximum_transfer_seconds = max(r["donor_to_seed_elapsed_seconds_as_recorded"] for r in seeds)
    require(all(r["within_passage_count_role"] == "last_count" for r in full_snapshot_links),
            "outgoing links", "an OLS donor is not its group's last count")
    require((len(rows), len(segments), len(passages)) == (335, 114, 131),
            "scope", "expected 335 counts, 114 OLS segments and 131 recorded passages")

    checks = {
        "count_rows_with_snapshot_and_archived_calendar_day_match": len(rows),
        "unique_ols_passages_and_model_segments": len(segments),
        "2N_passages": categories["2N", "seeding"],
        "2N_count_rows": categories["2N", "seeding"] + categories["2N", "harvest"],
        "4N_passages": categories["4N", "seeding"],
        "4N_count_rows": categories["4N", "seeding"] + categories["4N", "harvest"],
        "seeding_count_rows": len(seeds),
        "harvest_count_rows": sum(r["snapshot_event"] == "harvest" for r in rows),
        "intermediate_count_rows": sum(r["within_passage_count_role"] == "intermediate_count" for r in rows),
        "last_count_rows": len(last_counts),
        "metadata_and_snapshot_match_rows": count_source[True],
        "snapshot_only_rows": count_source[False],
        "seed_correctedCount_matches_production_initial_live_cells_passages": len(segments),
        "last_correctedCount_matches_production_target_live_cells_passages": len(segments),
        "minimum_ols_passage_calendar_duration_days": min(r["ols_calendar_day"] for r in last_counts),
        "maximum_ols_passage_calendar_duration_days": max(r["ols_calendar_day"] for r in last_counts),
        "calendar_day_before_selected": end_relation["before"],
        "calendar_day_at_selected": end_relation["at"],
        "calendar_day_after_selected": end_relation["after"],
        "passages_with_calendar_day_after_selected": len({r["passage_id"] for r in after_end}),
        "last_count_day_equals_selected_day_passages": sum(
            r["ols_day_vs_seed144_selected_day"] == "at" for r in last_counts
        ),
        "timestamp_age_exceeds_nominal_duration_rows": len(beyond_duration),
        "timestamp_age_exceeds_nominal_duration_and_calendar_day_after_selected_rows": sum(
            r["timestamp_age_exceeds_model_duration"] for r in after_end
        ),
        "model_integration_dt_days": dt,
        "timestamp_age_off_model_dt_grid_rows": sum(not r["timestamp_age_on_model_dt_grid"] for r in rows),
        "timestamp_age_minus_calendar_day_min_days": min(
            Decimal(r["age_days_from_recorded_timestamps"]) - r["ols_calendar_day"] for r in rows
        ),
        "timestamp_age_minus_calendar_day_max_days": max(
            Decimal(r["age_days_from_recorded_timestamps"]) - r["ols_calendar_day"] for r in rows
        ),
        "snapshot_count_rows_with_child_seed_in_131": len(linked_counts),
        "child_seed_links_in_131_from_ols_count_rows": sum(
            len(r["linked_child_seed_ids_in_131"].split(";")) for r in linked_counts
        ),
        "count_rows_with_child_seed_in_full_snapshot": len(full_snapshot_links),
        "child_seed_links_in_full_snapshot_from_ols_count_rows": sum(
            len(r["linked_child_seed_ids_in_snapshot"].split(";")) for r in full_snapshot_links
        ),
        "incoming_physical_parent_matches_model_parent_for_ols_passages": len(seeds),
        "incoming_donor_to_seed_elapsed_over_one_day_passages": sum(
            r["donor_to_seed_elapsed_seconds_as_recorded"] > 86400 for r in seeds
        ),
        "incoming_donor_to_seed_elapsed_over_seven_days_passages": sum(
            r["donor_to_seed_elapsed_seconds_as_recorded"] > 7 * 86400 for r in seeds
        ),
        "maximum_incoming_donor_to_seed_elapsed_days": f"{Decimal(maximum_transfer_seconds) / 86400:.6f}",
        "count_rows_with_Countess_value": sum(bool(r["db_Countess"]) for r in rows),
        "count_rows_with_growthType_value": sum(bool(r["db_growthType"]) for r in rows),
        "count_rows_with_second_parent_id": sum(bool(r["second_parent_record_id"]) for r in rows),
        "rows_with_numeric_cellCount_different_from_correctedCount": sum(
            value(r["db_cellCount"]) != value(r["db_correctedCount"]) for r in rows
        ),
        "rows_noted_images_unavailable": sum("images unavailable" in r["db_comment"] for r in rows),
        "sum159_pilot_rows_matching_event_timestamp_count": len(prepared),
        "sum159_pilot_matching_rows_annotated_visible_objects_including_dead": sum(
            r["sum159_pilot_count_semantics_if_matched"]
            == "whole_flask_corrected_count_including_live_MN_and_dead_objects" for r in rows
        ),
        "morphology_annotation_matches_by_event_timestamp": len(dead_annotations),
        "morphology_matches_with_dead_objects_classified": len(positive_dead),
        "morphology_matches_with_dead_objects_classified_seeds": sum(
            r["within_passage_count_role"] == "seed" for r in positive_dead
        ),
        "morphology_matches_with_dead_objects_classified_harvests": sum(
            r["snapshot_event"] == "harvest" for r in positive_dead
        ),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write(args.out_dir / "count_event_mapping.tsv", rows)
    write(args.out_dir / "audit_checks.tsv", [{"check": k, "value": v} for k, v in checks.items()])
    print(f"Audited {len(rows)} count events in {len(segments)} production segments.")


if __name__ == "__main__":
    main()
