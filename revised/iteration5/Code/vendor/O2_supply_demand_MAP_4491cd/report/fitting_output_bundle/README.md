# Fitting-output bundle exporter

`build_fitting_output_bundle.py` builds an auditable delivery package from an
existing multi-warmup joint fit, its exact in-vivo/in-vitro separate fits, and
the requested fixed-O2 supporting analyses. It is consume-only: it copies
recorded numerical/configuration/provenance outputs and validates recorded
objectives and source relationships without rerunning a model or reading the
contents of RDS files.

## Current requested export

Run from the repository root on the HPC checkout:

```bash
python3 oxygen/code/O2_supply_demand_MAP/report/fitting_output_bundle/build_fitting_output_bundle.py \
  --joint-root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_joint_multi_warmup_tsne_sigmaN0p1216_objmin_500seed_20260714_021540 \
  --invivo-root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invivo_O2_buffering_500seed \
  --invitro-root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed \
  --pooled-embedding-root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/best_fit_parameter_feature/04_combine_parameter_landscape/pooled_embedding_curve_class \
  --fixo2-eigen-embedding-root=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/best_fit_parameter_feature/06_combine_FixO2_eigen_attractor/pooled_embedding_curve_class \
  --fixed-o2-classification-tables=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification/dense-grid_monotonicity_classification/tables \
  --fixed-o2-regression-classification-tables=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/best_fit_parameter_feature/03_dense-grid_monotonicity_classification/monotonicity_classification/dense-grid_monotonicity_regression_classification/tables \
  --fixed-o2-attractor-tables=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/analysis/best_fit_parameter_feature/01_fixed_o2/FixO2_invivo_500seed/attractors/tables \
  --output-dir=/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fitting_output_bundle_20260722
```

The exporter writes the structured directory, a sibling `.tar.gz`, and a
sibling `.tar.gz.sha256`. It fails rather than replacing an existing result;
use `--force` only when intentionally replacing that exact export path.

## Output contract

- `README.md`: source roots, selection semantics, and bundle layout.
- `selected_results.tsv`: pair-level joint best seeds and global separate-fit
  best seeds with objective values.
- `source_relationships.tsv`: exact separate-fit anchors, joint run paths, and
  selected joint seed relationships.
- `supporting_analysis_sources.tsv`: exact source/bundle roots, report
  relationship, direct-input status, file count, and byte total for each
  supporting-analysis directory.
- `VALIDATION.tsv`: objective-component, best-seed, report-reference, and
  supporting-directory copy-count checks.
- `FILE_MANIFEST.tsv`: destination/source mapping, file sizes, and SHA-256.
- `SHA256SUMS.txt`: checksums for all bundle contents.
- `joint_fit/`: multi-warmup overview plus all six selected joint-seed outputs.
- `separate_fits/`: all-seed summaries plus exact anchor/global-best outputs.
- `supporting_analysis/pooled_embedding_curve_class/`: complete directory
  underlying `pooled_embedding_curve_class_report.html`.
- `supporting_analysis/fixo2_eigen_attractor_embedding_curve_class/`: complete
  directory underlying
  `fixo2_eigen_attractor_embedding_curve_class_report.html`.
- `supporting_analysis/fixed_o2_analysis_tables/`: both directly referenced
  dense-grid fixed-O2 classification table directories plus foundational
  analytical attractor tables.

Non-selected joint seed directories and their raster/PDF visualizations are not
duplicated. The two complete supporting report directories retain their HTML,
PDF, PNG, and tabular artifacts. No `SHARING_NOTE.md` is generated.
