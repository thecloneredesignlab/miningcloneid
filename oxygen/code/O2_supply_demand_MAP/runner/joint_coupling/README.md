# Joint coupling pipeline runner

`run_joint_coupling_pipeline.R` is the canonical local and HPC entrypoint. It
sequences, without reimplementing, these stages:

```text
ratio master -> within pair -> between pair -> process map
ploidy classification -> within pair -> between pair -> Cat x Class
visualization -> report
```

Required: `--result_root=PATH`. Optional `--output_root=PATH` defaults to
`oxygen/results/analysis/joint_coupling/<basename(result_root)>`. Useful controls include
`--class_threshold=1.1`, `--max_pairs=N`, `--max_seeds=N`,
`--n_boot=2000`, `--permutations=999`, `--run_vis=TRUE|FALSE`,
`--run_report=TRUE|FALSE`, and `--dry_run=TRUE|FALSE`.

The fitting `result_root` is read-only. The runner rejects any `output_root`
inside it. The output tree is:

```text
<output_root>/
  tables/{soft_coupling,ploidy_coupling}/
  figures/{overview,soft_coupling,ploidy_categories,category_association,appendix}/
  report/
  manifests/
  logs/
```

There are no per-layer output overrides: one root prevents split provenance.
`manifests/pipeline_config.tsv` records the input/output roots and threshold;
`manifests/manifest_index.tsv` inventories stage manifests.
