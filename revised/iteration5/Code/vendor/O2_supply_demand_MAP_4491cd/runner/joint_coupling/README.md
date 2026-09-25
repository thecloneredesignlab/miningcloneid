# Joint coupling pipeline runner

`run_joint_coupling_pipeline.R` is the canonical local and HPC entrypoint. It
sequences, without reimplementing, these stages:

```text
ratio master -> within pair -> between pair -> process map
ploidy classification -> within pair -> between pair -> Cat x Class
fixed-O2 seed mapping -> steady-state curve-class summaries
visualization -> report
```

For a raw fit tree, use `--result_root=PATH`. For a threshold-only sensitivity
rerun, use `--source_analysis_root=PATH` to reuse a validated ratio master and
unchanged ploidy-classification tables. Optional `--output_root=PATH` defaults to
`oxygen/results/analysis/joint_coupling/<basename(result_root)>`. Useful controls include
`--class_threshold=1.1`, or the explicit trio
`--class_lower_bound=0.8 --class_upper_bound=1.2
--class_boundary_rule=outer_inclusive`, plus `--max_pairs=N`, `--max_seeds=N`,
`--n_boot=2000`, `--permutations=999`, `--run_vis=TRUE|FALSE`,
`--run_report=TRUE|FALSE`, and `--dry_run=TRUE|FALSE`.
The fixed-O2 source may be supplied explicitly with
`--fixed_o2_source_root=PATH`; otherwise the runner requires exactly one
complete matching source under
`oxygen/results/analysis/warm_up_joint_fitting_results_extra/`.

The fitting `result_root` is read-only. The runner rejects any `output_root`
inside it. The output tree is:

```text
<output_root>/
  tables/{soft_coupling,ploidy_coupling,fixed_o2_ploidy_classification}/
  figures/{overview,soft_coupling,ploidy_categories,fixed_o2_ploidy_classification,category_association,appendix}/
  report/
  manifests/
  logs/
```

There are no per-layer output overrides: one root prevents split provenance.
`manifests/pipeline_config.tsv` records the input/output roots, source-analysis
root when used, exact class bounds, and boundary ownership;
the fixed-O2 source root is also recorded so the 3,000 synthetic seeds remain
traceable to their joint-fit pair and seed;
`manifests/manifest_index.tsv` inventories stage manifests.
