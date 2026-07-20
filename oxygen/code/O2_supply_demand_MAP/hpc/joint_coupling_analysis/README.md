# HPC joint coupling analysis

| File | Responsibility |
|---|---|
| `submit_joint_coupling_analysis.sh` | Validates `RESULT_ROOT`/`OUTPUT_ROOT`, rejects output inside the fitting tree, creates `<output_root>/logs`, and submits one resource-configured Slurm job; `--dry_run=TRUE` prints the exact command. |
| `run_joint_coupling_analysis.sub` | Requires both roots, loads R on the compute node, and invokes the canonical runner with `--result_root` and `--output_root`. |

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
