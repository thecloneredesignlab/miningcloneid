#!/usr/bin/env python3
"""Write the WP0 provenance contract/report from verified replay artifacts."""
import csv
import hashlib
import json
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()


def read_tsv(relative):
    with (root / relative).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


status = json.loads((root / "replay/replay_status.json").read_text())
settings = json.loads((root / "manifests/replay_settings.json").read_text())
execution = json.loads((root / "manifests/execution.json").read_text())
configuration = json.loads((root / "manifests/runtime_model_configuration.json").read_text())
checks = read_tsv("audit/replay_checks.tsv")
frozen = read_tsv("manifests/input_manifest.tsv")
hash_changes = [row["path"] for row in frozen if sha(root / row["path"]) != row["sha256"]]
assert status["pass"] and all(row["pass"] == "TRUE" for row in checks)
assert not hash_changes, hash_changes
assert status["forward_evaluations"] == 1 and not status["optimization"]
package = read_tsv("audit/production_package_md5.tsv")
assert len(package) == 25 and all(row["matches"] == "True" for row in package)
baseline = {row["metric"]: row for row in read_tsv("audit/baseline_replay.tsv")}
ledger = read_tsv("audit/objective_components.tsv")
observations = read_tsv("manifests/observation_manifest.tsv")
selections = read_tsv("audit/selected_passage_days.tsv")
boundaries = read_tsv("audit/conditioning_audit.tsv")
package_versions = read_tsv("audit/runtime_package_comparison.tsv")
version_mismatches = [row for row in package_versions if row["version_match"] != "True"]
nominal_changes = sum(float(row["nominal_duration_days"]) != float(row["selected_day"]) for row in selections)
boundary_changes = sum(row["differs_from_observation_day"] == "TRUE" for row in boundaries)
missing_targets = sum(row["target_live_cells"] == "NA" for row in selections)
image_hash = (root / "source/provenance/replay_image_sha256.txt").read_text().split()[0]
script_hash = sha(root / "scripts/replay_culture_baseline.R")
fingerprint_fields = {
    "input_manifest_sha256": sha(root / "manifests/input_manifest.tsv"),
    "replay_script_sha256": script_hash,
    "settings_sha256": sha(root / "manifests/replay_settings.json"),
    "replay_image_sha256": image_hash,
}
contract = {
    "study_id": root.name,
    "scope": "WP0 production-source audit and one no-optimization culture baseline replay",
    "status": "NUMERICAL_REPLAY_PASS_WITH_HISTORICAL_CONTAINER_IDENTITY_CAVEAT",
    "production_reference": {
        "seed": 144, "commit": settings["production_commit"], "branch": "soft_couping_org",
        "task_git_dirty_status": "clean",
        "task_provenance": "source/provenance/task144/run_provenance.tsv",
        "source_archive": "source/code/production_1a5d861.tar",
        "archive_sha256": sha(root / "source/code/production_1a5d861.tar"),
        "remote_fit_directory": "/share/lab_crd/taoli/Project/soft_couping_org/oxygen/results/"
        "fit_invitro_unified_np256_500seed_all_xxlarge_r442_exact_20260828_145253/seed144",
        "iteration5_manifest_files_verified": len(package),
        "parameter_authority": "source/seed144/fit_result.rds$best_params (full precision)",
        "configuration_authority": "source/seed144/fit_result.rds$cfg",
        "allowed_configuration_changes": ["parameter_table path relocation only"],
        "optimization_status": {"DE_iterations": 415, "DE_iteration_budget": 1000,
                                "NP_requested": 80, "NP_used": 256,
                                "local_accepted": True, "local_convergence_code": 1,
                                "certified_global_optimum": False},
    },
    "objective": {
        "definition": "Q = -mean(growth log scores) - mean(sample-mean karyotype log scores) - mean(flow cross-entropies)",
        "weights": settings["weights"], "normalization_counts": {"growth": 114, "karyotype": 12, "flow": 20},
        "karyotyped_cells": 220, "unique_passage_observation_ids": len(observations),
        "growth_model": "Gaussian score for interval growth; fitted sigma_growth",
        "karyotype_model": "Gaussian smoothing of selected chromosome distribution; fitted sigma_kary; probability floor 1e-12",
        "flow_model": "Normalized density on observed grid, observed grid masses normalized within sample; density floor 1e-12",
        "flow_kernel_sd": float(baseline["sigma_flow_ploidy"]["replay"]),
        "flow_kernel_sd_policy": "max(2 * median spacing of pooled observed flow grids, 1/22); not a fitted coordinate",
        "prior": "No prior contribution; prior-related cfg fields exist but this objective does not use them",
        "soft_coupling": "absent", "phenotype_regularization": "absent",
        "death_score": "disabled in production loader, no scored death observations",
        "failure_policy": {
            "generic_simulation_failure_objective": 1e9,
            "protocol_infeasible_objective": "min(1e6 + 1e4*remaining_segments + 1e3*min(relative_shortfall,10), 1e9-1)",
            "baseline_failure_penalty": 0,
            "retained_sample_rule_for_pilot": "Freeze the observed sample set; reject/penalize missing predictions instead of filtering them from averages",
        },
        "formal_confidence_intervals_or_LR_thresholds": False,
    },
    "numerics": {"DT_days": 0.05, "chromosome_grid": [22, 154], "N_UNIT": 22,
                 "fallback_max_passage_days": 14, "initial_total_size_default": 1e6,
                 "fixed_assigned_oxygen": True, "ploidy_O2_death": "ploidy_related",
                 "full_saved_configuration": "manifests/runtime_model_configuration.json"},
    "passage_contract": {
        "implementation": "v2", "segments": len(selections), "parent_child_boundaries": len(boundaries),
        "observation_selection": "Among positive-day candidates, closest live count to segment final-cell target; ties use first index; if no positive target, last candidate",
        "target_fallback": "If segment final count is unavailable, try next segment initial count; otherwise last candidate",
        "parent_to_child": "Use parent v2 extraction for requested child inoculum; require sufficient live cells; nearest eligible surplus then first index if necessary; downsample, never create cells",
        "nominal_duration_vs_selected_day_differences": nominal_changes,
        "missing_positive_targets": missing_targets,
        "observation_vs_propagation_day_differences_at_this_endpoint": boundary_changes,
        "initial_and_final_counts_remain_conditioning_inputs_after_growth_score_deletion": True,
    },
    "replay": {**status, "absolute_tolerance": settings["absolute_tolerance"],
               "relative_tolerance": settings["relative_tolerance"], "selected_day_tolerance": 0,
               "objective_difference": float(baseline["objective"]["difference"]),
               "frozen_input_files_verified_unchanged": len(frozen)},
    "runtime": {"image": execution["image_path"], "sha256": image_hash,
                "entrypoint": execution["rscript_entrypoint"], "R": "4.4.2", "platform": "x86_64-linux",
                "OS": "Rocky Linux 8.10", "BLAS": "FlexiBLAS OPENBLAS; LAPACK 3.12.0",
                "compiler": "GCC 13.3.0; R CXX: g++ -std=gnu++17",
                "compiler_evidence": "replay/container_compiler_and_environment.txt",
                "code_archive_HPC_snapshot_package_comparison": package_versions,
                "code_archive_HPC_snapshot_package_mismatches": version_mismatches,
                "package_resolution_correction_evidence": "audit/package_resolution_correction.json",
                "historical_launch_digest_proven": False,
                "historical_identity_caveat": execution["historical_image_identity"]},
    "fingerprint_fields": fingerprint_fields,
    "run_fingerprint_sha256": hashlib.sha256(json.dumps(fingerprint_fields, sort_keys=True).encode()).hexdigest(),
    "authorized_work_completed": True,
    "pilot_engineering_or_fits_started": False,
}
(root / "manifests/production_contract.json").write_text(json.dumps(contract, indent=2, ensure_ascii=False) + "\n")
table = "\n".join(
    f"| {row['term']} | {row['denominator']} | {row['weight']} | {float(row['objective_contribution']):.15g} |"
    for row in ledger
)
report = f"""# Culture seed144: production provenance and baseline replay without optimization

The replay passed {status['checks']:,} checks with zero failures. WP0 performed one forward evaluation and no optimization or parameter scan. These results establish reproducibility of this production endpoint; they do not establish identifiability or mechanism necessity.

## Production provenance

- The complete seed144 bundle came from the HPC directory recorded in production_contract.json. The original directory was preserved.
- All 25 seed144 files in the iteration5 production manifest matched their MD5 checksums, including fit_result.rds.
- Task144 recorded clean production commit `{settings['production_commit']}` on 2026-08-28 at 21:59:54 EDT. The frozen source archive was exported from that commit.
- All {len(frozen)} frozen source/input SHA256 hashes remained unchanged after replay. The four HPC data files and parameter table matched their local counterparts.
- Full-precision parameters and cfg came from fit_result.rds. Only the parameter-table path was relocated; rounded figure parameters were not used to reconstruct the endpoint.
- The checkout's R/C++ model core matched production, but its culture passage, scoring and fitting helpers differed. See source_comparison.tsv.

## Numerical results and independent scoring

Archived Q: `{status['archived_objective']:.15g}`; replayed Q: `{status['replay_objective']:.15g}`; independently recomputed Q: `{status['independent_objective']:.15g}`. Replay minus archive: `{float(baseline['objective']['difference']):.17g}`.

The objective combines within-modality average scores. It is not a likelihood to which ordinary likelihood-ratio cutoffs can be applied directly. Independent arithmetic recomputed Gaussian growth scores, smoothed chromosome probabilities and flow-density cross-entropies from raw observations and selected predictions.

| Term | Averaging denominator | Weight | Contribution to Q |
|---|---:|---:|---:|
{table}

The retained set contains 131 unique passage/sample IDs: 66 in 2N and 65 in 4N. Scoring uses 114 growth records, 12 chromosome samples containing 220 cells, and 20 flow distributions. Modalities overlap; these counts are not independent biological replicates. Flow-grid points and simulated daily values are not additional observations.

Checks compared observation identities, scores, growth/flow predictions, complete live/dead chromosome-state trajectories, distributions, oxygen trajectories, selected states and propagation states. Numerical tolerance was fixed before execution as `abs(diff) <= 1e-8 + 1e-8 * max(abs(reference), abs(replay))`; IDs and observation/propagation days required exact agreement. All 7,766 checks passed. Nonzero floating-point differences remain in replay_checks.tsv.

## Passage and conditioning contract

Production v2 contains {len(selections)} simulated segments and {len(boundaries)} parent-child boundaries. Observation states are selected by proximity to the target final live-cell count; ties choose the first candidate. For {missing_targets} segments without a usable target, the final candidate day was used. Selected days differed from nominal durations in {nominal_changes} segments, so nominal durations cannot replace the production selection rule.

Propagation separately requires sufficient cells for the child inoculum and downsamples an eligible state. If no state supplies the inoculum, the protocol is infeasible. All 120 baseline propagations used downsample_to_observed_inoculum. Observation and propagation days differed at {boundary_changes} boundaries in this endpoint; this does not guarantee their equality elsewhere in parameter space.

Deleting the growth score can leave the same counts involved in state selection and propagation. Such comparisons are information-contribution audits unless a separate leakage-free validation protocol is designed.

## Runtime and remaining provenance limitation

- Successful Slurm job: 21797874, COMPLETED, exit 0; one core, 4 GB, 39 seconds total scheduler wall time.
- Forward evaluation and scoring: {status['forward_seconds']:.3f} seconds. Audit script including compilation: {status['total_seconds']:.2f} seconds. Slurm MaxRSS: 806,472 KiB. A single forward evaluation does not measure optimization cost.
- SIF: `{execution['image_path']}`; SHA256: `{image_hash}`.
- Runtime entry: /usr/local/bin/o2-hpc-exact-rscript through scripts/agentRrunner.sh. R 4.4.2; Rocky Linux 8.10; FlexiBLAS OPENBLAS/LAPACK 3.12.0; GCC 13.3.0 and R CXX g++ -std=gnu++17.
- Package metadata correction on 2026-11-09: the earlier table selected lower-priority installed copies instead of loaded versions. Original sessionInfo records Rcpp 1.1.1-1.1 and jsonlite 2.0.0. All seven checked packages match the archived HPC snapshot after correction. Raw replay and frozen inputs were unchanged; see package_resolution_correction.json.
- Initial job 21797852 failed while loading Matrix, before any forward evaluation. Generic Rscript had bypassed the image-specific library-path setup. Both attempts' logs are retained.

The historical task confirms container execution but does not record its SIF path/digest, and its Slurm environment was not retained. The old image_runtime_lock.tsv cannot establish identity with the current SIF. This replay fingerprints the current image and reproduces the endpoint numerically; exact historical launch-image identity remains unproven.

## Contract for subsequent work

1. Adapt the verified production v2 culture implementation rather than the legacy checkout passage path.
2. The culture objective contains no prior, soft-coupling or phenotype regularization contribution; death scoring is disabled. Prior-related cfg fields do not imply an active penalty, so a duplicate data-only/penalized study is unnecessary.
3. Preserve production failure barriers: generic simulation errors score 1e9; protocol infeasibility has progress/shortfall penalties. Baseline failure penalty is zero.
4. Freeze retained observations and reject or penalize missing predictions. The replay had no missing scored growth prediction or discarded required sample.
5. Production v2 already injects the parameter-table initial vector into DE member 1. Extend that mechanism to actual endpoint and neighbor vectors, verifying the population received by DEoptim.
6. Seed144's local convergence code is 1; DE stopped at 415/1000 iterations. It is a reproducible reference and starting candidate, not a certified global optimum. The zero-extended domain needs its own unrestricted refit.
7. At WP0 completion, mapping adaptation, restricted optimization, profiles, COMBOS and manuscript rebuilding had not begun. Subsequent work has separate phase manifests.

## Evidence

- manifests/production_contract.json: model, score, conditioning, runtime and tolerances.
- manifests/input_manifest.tsv and audit/production_package_md5.tsv: frozen hashes and iteration5 comparison.
- manifests/observation_manifest.tsv: IDs, lineage, oxygen, times, modalities and conditioning roles.
- audit/objective_components.tsv, baseline_replay.tsv and replay_checks.tsv: score ledger and checks.
- audit/selected_passage_days.tsv and conditioning_audit.tsv: observation selection and propagation.
- replay/selected_distributions.tsv, replay_components.rds and flow_overlay_df.tsv: full predictions.
- replay/sessionInfo.txt, installed_packages.tsv, source/provenance and logs: runtime and raw evidence.

Production fits and manuscript assets were preserved. This English report replaces the earlier Chinese report without changing the numerical evidence or historical-image caveat.
"""
(root / "reports/baseline_replay_report.md").write_text(report)
with (root / "audit/frozen_input_integrity.tsv").open("w") as stream:
    writer = csv.writer(stream, delimiter="\t")
    writer.writerow(["files_checked", "changed_files", "pass"])
    writer.writerow([len(frozen), len(hash_changes), not hash_changes])
print(f"Report and production contract written: {root}")
