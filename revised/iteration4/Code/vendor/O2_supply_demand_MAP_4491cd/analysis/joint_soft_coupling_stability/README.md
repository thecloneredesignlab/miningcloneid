# Joint soft-coupling stability analysis

This folder answers which of the 14 joint-fit parameters are consistently
lower in vivo (ClassA), approximately coupled (ClassB), or higher in vivo
(ClassC), first within each pair and then across pairs. The legacy/default
classification uses the symmetric `1.1` threshold:

- ClassA: `0 < in_vivo/in_vitro < 1/1.1`
- ClassB: `1/1.1 <= in_vivo/in_vitro <= 1.1`
- ClassC: `in_vivo/in_vitro > 1.1`

Boundary values belong to ClassB. Non-positive or non-finite ratios are
retained as `Invalid` and excluded from class denominators.

An explicit asymmetric interval is also supported. For example,
`--class_lower_bound=0.8 --class_upper_bound=1.2
--class_boundary_rule=outer_inclusive` gives ClassA `ratio <= 0.8`, ClassB
`0.8 < ratio < 1.2`, and ClassC `ratio >= 1.2`. The exact bounds and boundary
ownership are written into every master row and `analysis_config.tsv`.

## Files

| File | Exact responsibility | Main outputs |
|---|---|---|
| `build_joint_soft_coupling_master_table.R` | Discovers completed pair directories, joins soft-coupling rows to seed objectives, independently recomputes each ratio, validates the 14-parameter/key contract, attaches pair metadata and the audited biological-process map. | `soft_coupling_master_long.tsv`, excluded rows, input-quality/config/pair/process tables, manifest. |
| `reclassify_joint_soft_coupling_master_table.R` | Reclassifies a validated materialized master under a new interval without changing ratios or fit diagnostics; records the source path and MD5 and carries forward validated metadata tables. | A new independent master/config/manifest set for sensitivity reruns. |
| `compare_joint_coupling_classification_versions.R` | Verifies identical keys/ratios, then compares source versus current calls at seed-parameter, pair-parameter, cross-pair, and pair-balanced class-share levels. | `tables/classification_comparison/*.tsv` and manifest. |
| `analyze_within_pair_soft_coupling.R` | Treats seeds within a pair as repeated fitted solutions; calculates A/B/C proportions, union, strict intersection, 80/90/95% cores, dominant class, entropy, raw/log2 moments, ClassB-only direction around 1, seed-set Jaccard, threshold/objective sensitivity, and feasibility/boundary rates. | `within_pair_*.tsv`, `threshold_sensitivity_pair_balanced.tsv`, `objective_sensitivity_pair_balanced.tsv`, and manifest. |
| `analyze_between_pair_soft_coupling.R` | Treats each pair equally; calculates cross-pair dominant-class agreement, graded/strict support, pair-bootstrap intervals, leave-one-pair-out stability, pair Jensen-Shannon similarity, within/between variance and ICC. It also materializes union/stable/strict pair-set patterns for UpSet-style displays. Pooled seeds are explicitly secondary. | `between_pair_*.tsv`, `class_intersection_pattern_summary.tsv`, `class_pair_set_membership.tsv`, and manifest. |
| `summarize_soft_coupling_processes.R` | Joins parameter stability to the audited process mapping and prepares process-level and report-ready summaries. | `parameter_class_process_summary.tsv`, `parameter_process_consensus.tsv`, `process_class_summary.tsv`, report summary. |

## Input contract

Each pair must contain:

```text
extra_results/joint_soft_coupling_all.tsv
extra_results/seed_summary.tsv
```

The key is `pair_id + seed + parameter`. Every included seed must have exactly
the 14 parameters listed in `o2jca_parameter_levels()`. The reported ratio is
checked against `vivo_natural / vitro_natural` with relative tolerance `1e-10`.

## Statistical scope

The six target pairs use different in-vivo warm-up representatives but the
same in-vitro `vt_seed10` anchor. Between-pair stability therefore means
stability across those in-vivo families conditional on that shared anchor.
It is not six independent in-vitro experiments.

## Direct execution

Use the feature runner rather than launching these scripts separately:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/runner/joint_coupling/run_joint_coupling_pipeline.R \
  --result_root=/absolute/path/to/fit_joint_multi_warmup_result \
  --output_root=/absolute/path/to/oxygen/results/analysis/joint_coupling/run_name \
  --class_threshold=1.1
```

To build a complete sensitivity version from an existing materialized result:

```bash
Rscript oxygen/code/O2_supply_demand_MAP/runner/joint_coupling/run_joint_coupling_pipeline.R \
  --source_analysis_root=/absolute/path/to/existing_joint_coupling_result \
  --output_root=/absolute/path/to/new_independent_result \
  --class_lower_bound=0.8 \
  --class_upper_bound=1.2 \
  --class_boundary_rule=outer_inclusive
```

`result_root` is a read-only input. If `--output_root` is omitted, the runner
uses `oxygen/results/analysis/joint_coupling/<basename(result_root)>`. Any
output path inside the fitting tree is rejected.

Tests: `test-joint-soft-coupling-stability.R` and
`test-joint-coupling-stage-split.R`.
