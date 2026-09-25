# HPC joint coupling analysis

| File | Responsibility |
|---|---|
| `submit_joint_coupling_analysis.sh` | Validates `RESULT_ROOT`/`OUTPUT_ROOT`, rejects output inside the fitting tree, creates `<output_root>/logs`, exports the exact ratio-class contract and optional `FIXED_O2_SOURCE_ROOT`, and submits one resource-configured Slurm job; `--dry_run=TRUE` prints the exact command. |
| `run_joint_coupling_analysis.sub` | Requires both fit/output roots, loads R on the compute node, and invokes the canonical runner with the classification and fixed-O2 source contracts. |

Default production input:

```text
/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540
```

Default production output:

```text
/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/joint_coupling/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540
```

The worker contains no statistical, classification, plotting, or reporting
logic. Logs are written below `<output_root>/logs/`.
The default wall time is 12 hours, matching the `xxlarge` QOS maximum.

The default remains the symmetric 1.1-fold rule. An outer-inclusive 0.8/1.2
run is requested with `--class_lower_bound=0.8 --class_upper_bound=1.2
--class_boundary_rule=outer_inclusive` and a distinct `--output_root`.
Supply `--fixed_o2_source_root=PATH` when automatic discovery under the
project's `oxygen/results/analysis/warm_up_joint_fitting_results_extra/`
directory is ambiguous or unavailable.
