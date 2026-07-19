# Joint warm-start runners

These entrypoints transform completed independent in-vivo and in-vitro fits
into parameter tables consumed by a later joint fit. They do not analyze a
post-fit simulation, so they live in `runner/` rather than `analysis/`.

| File | Function | Inputs | Outputs |
|---|---|---|---|
| `make_joint_soft_coupling_parameters_table.R` | Build one labelled center/delta joint soft-coupling start table. | One in-vivo and one in-vitro seed directory, each with `best_params.tsv`. | Optimizer-scale CSV with center and context-delta parameters. |
| `make_warm_start_from_separate_fit_results.py` | Build the four paired/best-seed warm-start variants used by the separate-fit workflow. | Completed in-vivo and in-vitro run directories and optional selection settings. | CSV/XLSX warm-start tables and selection metadata. |

Both scripts resolve the repository from their own canonical path and can be
called from any working directory.
