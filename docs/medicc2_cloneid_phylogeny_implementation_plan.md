# MEDICC2 CloneID Copy-Number Phylogeny Implementation Plan

## Goal

Build a portable workflow to rerun a prior MEDICC2 copy-number phylogeny analysis on a new, partially overlapping set of CloneID cells. The immediate biological use case is to test whether cells from nominal 4N-derived tumors are phylogenetically closer to the 4N cell-line population or to a possible 2N contaminant/outgrowth population. The workflow must be reusable outside the original repository and must not depend on hard-coded SUM159 origins, report paths, Dockerfile paths, or branch names.

The prior workflow is useful only as reference context. It established a practical data shape and runtime choice:

- extract CloneID GenomePerspective total-copy-number profiles;
- write a wide copy-number matrix with `origin`, `clone_label`, and chromosome/segment columns such as `1:1-249250621`;
- convert the wide matrix to MEDICC2 long TSV with `sample_id`, `chrom`, `start`, `end`, `cn_total`;
- run MEDICC2 1.1.2 with total copy numbers;
- optionally run pooled, per-origin/per-replicate, and event-calling/heatmap analyses;
- postprocess final Newick trees and event tables into WGD-annotated tree figures.

## Analysis Questions

Primary question:

- Do 4N-injected tumor cells share copy-number ancestry with the 4N cell-line reference, or are they closer to the 2N cell-line reference?

Secondary questions:

- Are 4N-injected tumors homogeneous, or do they contain multiple ancestry-like groups?
- Are there shared private copy-number events that support 4N-lineage origin versus 2N contamination/outgrowth?
- Does adding partially overlapping new cells change the placement of previously analyzed cells?
- Are WGD events inferred on internal branches, singleton leaves, or absent after including the new reference/tumor cells?

## Proposed Portable Directory Layout

Each analysis should write to a new analysis root and never overwrite old outputs.

```text
<analysis_root>/
  config/
    analysis_config.yaml
  inputs/
    selected_cells.tsv
    previous_metadata.json              # optional
    preexported_cloneid_profiles.csv    # optional if not querying CloneID live
  intermediate/
    wide_cn_matrix.csv
    cell_selection_audit.csv
    wide_cn_matrix.heatmap.png
    wide_cn_matrix.heatmap.pdf
  medicc2/
    pooled/
      input.tsv
      output/
      output_events/
    by_group/
      <group_slug>/
        input.tsv
        output/
        output_events/
  figures/
    pooled_wgd_annotated_tree.png
    <group_slug>_wgd_annotated_tree.png
  qc/
    validation_report.json
    overlap_with_previous_analysis.csv
    sample_id_map.csv
  manifests/
    run_manifest.json
    commands.sh
```

## Configuration Schema

Use one config file to drive all steps.

```yaml
analysis_slug: cloneid_4n_tumor_ancestry_medicc2
analysis_root: /path/to/output/cloneid_4n_tumor_ancestry_medicc2

profile_source:
  mode: cloneid_db          # cloneid_db | preexported_wide_csv
  perspective: GenomePerspective
  cloneid_connection: default
  preexported_wide_csv: null

selected_cells_file: /path/to/selected_cells.tsv
previous_metadata_file: /path/to/previous_metadata.json

copy_number_matrix:
  include_autosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
  exclude_sample_id_prefixes: ["SP_1.0_"]
  allow_sex_chromosomes: false
  segment_label_regex: "^[0-9XY]+:[0-9]+-[0-9]+$"
  sample_id_template: "{origin}__{clone_label}"

grouping:
  group_field: group
  run_pooled: true
  run_by_group: true
  groups_to_run: null

medicc2:
  image: cloneid-medicc2:1.1.2
  platform: linux/amd64
  executable_mode: docker       # docker | singularity | local
  total_copy_numbers: true
  input_allele_columns: cn_total
  run_events: true
  plot_mode: heatmap

postprocessing:
  plot_wgd_annotated_trees: true
  artificial_normal_names: ["diploid"]
```

## Input Requirements

### Cell Selection File

The user-provided selection file is the central portable input. It should be CSV or TSV with at least:

```text
origin
clone_label
include
group
```

Recommended columns:

```text
origin
clone_label
cell_id
sample_id
include
exclude_reason
group
replicate
source_role
expected_ploidy_label
notes
```

Example `source_role` values for this project:

- `cell_line_2N`
- `cell_line_4N`
- `tumor_from_4N_injection`
- `tumor_from_2N_injection`
- `normal_or_diploid_reference`

The workflow should treat `origin + clone_label` as the default CloneID lookup key. If `sample_id` is absent, construct it deterministically from `sample_id_template`. Do not use bare `clone_label` as the MEDICC2 sample ID unless uniqueness has been validated across all origins.

### CloneID Profile Source

Support two modes:

1. Live CloneID extraction:
   - connect to CloneID through a configured access method;
   - request GenomePerspective or equivalent total-copy-number profiles;
   - extract only requested origins/cells.

2. Pre-exported table:
   - read a wide matrix already containing the selected cell profiles;
   - still run the same schema, provenance, and validation checks.

This keeps the workflow usable in a fresh repo or workflow directory where CloneID credentials are unavailable.

### Optional Previous Metadata

If a previous-analysis metadata JSON is provided, use it only for overlap/audit reporting. The workflow must not require the old report path.

The metadata should provide or be transformed into:

- prior analysis slug;
- prior origins;
- prior row count;
- prior sample IDs or `origin + clone_label` pairs;
- prior chromosome/segment labels.

## Extraction Layer

Implement `extract_cloneid_cn_profiles`.

Responsibilities:

- read selected cells;
- filter `include == true`;
- query or read copy-number profiles;
- retain only configured chromosomes/segments;
- exclude accidental control/spike-in profiles such as `SP_1.0_` unless explicitly requested;
- write a wide matrix independent of repository paths.

Expected wide matrix schema:

```text
origin,clone_label,sample_id,group,replicate,source_role,<segment_1>,...,<segment_n>
```

Segment columns may be whole chromosomes or finer bins, but must use parseable labels:

```text
chrom:start-end
```

For whole-chromosome data, labels may look like:

```text
1:1-249250621
2:1-243199373
```

The extraction layer should also write:

- `intermediate/cell_selection_audit.csv`
- `intermediate/wide_cn_matrix.csv`
- `intermediate/wide_cn_matrix.heatmap.png`
- `intermediate/wide_cn_matrix.heatmap.pdf`
- `qc/sample_id_map.csv`
- `manifests/run_manifest.json` with extraction metadata.

The audit file should include:

```text
origin
clone_label
requested_sample_id
resolved_sample_id
include_requested
profile_found
included_in_matrix
missing_reason
group
replicate
source_role
```

## Overlap Audit

Implement `audit_cell_selection_overlap`.

If previous metadata is available, compute:

- cells in both prior and new runs;
- cells only in the prior run;
- cells only in the new run;
- origins retained, dropped, and newly added;
- segment/chromosome label compatibility;
- row counts by origin, group, replicate, and source role.

Write:

```text
qc/overlap_with_previous_analysis.csv
qc/validation_report.json
```

Overlap should be based on stable keys:

```text
origin + clone_label
```

and, when present:

```text
sample_id
cell_id
```

Do not assume `sample_id` is stable across runs unless it is explicitly supplied in both.

## MEDICC2 Conversion

Implement `convert_wide_cn_to_medicc2`.

Input:

- `intermediate/wide_cn_matrix.csv`

Output:

- `medicc2/pooled/input.tsv`
- optional `medicc2/by_group/<group_slug>/input.tsv`

Long TSV schema:

```text
sample_id
chrom
start
end
cn_total
```

Coordinate conversion:

- parse `chrom:start-end`;
- convert input `start=1` to BED-like `start=0`;
- keep `end` as supplied;
- reject malformed labels.

Validation:

- every included row has nonempty `sample_id`;
- `sample_id` values are unique;
- every sample has the same ordered segment set;
- chromosomes are configured autosomes 1-22 unless sex chromosomes are explicitly enabled;
- copy-number values are numeric and nonnegative;
- no duplicate `sample_id + chrom + start + end` rows are emitted;
- group-specific TSVs contain enough samples to justify a tree.

The converter should support:

- `--include-group <group>`;
- `--exclude-group <group>`;
- `--sample-id-column <col>`;
- `--group-column <col>`;
- `--segments autosomes|all|file:<path>`;
- `--dry-run` to validate and report planned outputs without writing MEDICC2 inputs.

## Runtime Environment

Use a reproducible MEDICC2 runtime. Preferred:

```text
MEDICC2 1.1.2
linux/amd64
gcc 11
openfst 1.8.2
```

The prior workflow used a Docker image equivalent to:

```text
cloneid-medicc2:1.1.2
```

The portable workflow should expose:

- image tag;
- platform;
- container engine: Docker or Singularity;
- host analysis root;
- container mount point, e.g. `/work`;
- command log path.

Do not assume a Dockerfile path exists. If building is needed, provide a separate environment-build instruction or use a prebuilt image reference.

## MEDICC2 Runs

Implement `run_medicc2`.

### Pooled Run

Canonical command shape:

```bash
medicc2 <input.tsv> <output_dir> \
  --total-copy-numbers \
  --input-allele-columns cn_total \
  --prefix <prefix>
```

Containerized Docker shape:

```bash
docker run --rm --platform linux/amd64 \
  -v "<analysis_root>:/work" \
  <medicc2_image> \
  /work/medicc2/pooled/input.tsv \
  /work/medicc2/pooled/output \
  --total-copy-numbers \
  --input-allele-columns cn_total \
  --prefix <analysis_slug>_pooled
```

### Event Calling / Heatmap Run

If configured:

```bash
docker run --rm --platform linux/amd64 \
  -v "<analysis_root>:/work" \
  <medicc2_image> \
  /work/medicc2/pooled/input.tsv \
  /work/medicc2/pooled/output_events \
  --total-copy-numbers \
  --input-allele-columns cn_total \
  --events \
  --plot heatmap \
  --prefix <analysis_slug>_pooled_events
```

### Per-Group / Per-Replicate Runs

For each configured group:

```text
medicc2/by_group/<group_slug>/input.tsv
medicc2/by_group/<group_slug>/output/
medicc2/by_group/<group_slug>/output_events/
```

Run group-specific analyses when:

- there are enough cells per group;
- the question benefits from within-group topology;
- the pooled tree may be dominated by between-group differences.

For the 4N contamination question, recommended run modes are:

1. pooled: `2N_cell_line + 4N_cell_line + 4N_tumors + optional_2N_tumors`;
2. pooled event-calling run;
3. optional group-specific trees for each tumor/replicate if cell counts are sufficient;
4. optional reference-only tree with `2N_cell_line + 4N_cell_line` to characterize baseline separation.

## Postprocessing

Implement `collect_medicc2_outputs`.

Collect expected files from each run:

- `<prefix>_summary.tsv`
- `<prefix>_pairwise_distances.tsv`
- `<prefix>_branch_lengths.tsv`
- `<prefix>_final_cn_profiles.tsv`
- `<prefix>_final_tree.new`
- `<prefix>_final_tree.xml`
- `<prefix>_final_tree.png`
- `<prefix>_cn_profiles_heatmap.pdf` when heatmap plotting is requested;
- `<prefix>_copynumber_events_df.tsv` when `--events` is requested.

Implement `plot_wgd_annotated_tree`.

Inputs:

- final Newick tree: `<prefix>_final_tree.new`;
- event table: `<prefix>_copynumber_events_df.tsv`.

Output:

- `figures/<prefix>_wgd_annotated_tree.png`

The plotting implementation should be portable. Either:

- vendor a small script that parses Newick and highlights event-table rows where `type == WGD`; or
- define a dependency-managed plotting module using standard packages.

The plot should summarize:

- internal WGD-bearing branches;
- WGD singleton leaves;
- artificial normal/diploid leaf if MEDICC2 adds one;
- sample labels colored by group/source role if metadata is available.

## Validation And QC

### Before MEDICC2

Required checks:

- all requested included cells were found or explicitly reported missing;
- excluded cells are listed with reasons;
- no duplicate sample IDs remain;
- chromosome/segment labels parse correctly;
- segment order is deterministic;
- autosomes 1-22 are retained unless configured otherwise;
- no accidental control/spike-in profiles such as `SP_1.0_` are included unless explicitly requested;
- copy-number values are numeric and nonnegative;
- each sample has complete profiles across all selected segments;
- overlap with prior cell set is summarized when prior metadata is provided;
- group/replicate/source-role counts match expectations;
- generated filenames and prefixes are deterministic.

### After MEDICC2

Required checks:

- expected output files exist for each configured run;
- final Newick tree contains all expected samples plus any artificial normal/diploid sample;
- summary, branch-length, pairwise-distance, and final-CN tables parse cleanly;
- event table exists and parses if `--events` was requested;
- WGD annotations distinguish internal branch events from leaf/singleton events;
- command lines, image tag, platform, input hashes, timestamps, and output paths are recorded in the run manifest.

## Testing Plan

### Unit Tests With Synthetic Data

Create a tiny wide matrix:

- 4 samples;
- 3 chromosomes or segments;
- two groups;
- one duplicated `clone_label` across origins to test deterministic `sample_id` construction.

Test:

- segment parsing and coordinate conversion;
- sample ID uniqueness validation;
- group filtering;
- missing value rejection;
- autosome filtering;
- output TSV row count equals `n_samples * n_segments`.

### Dry-Run Tests On Real Metadata

Use a real selected-cells file without querying profiles or running MEDICC2.

Dry-run should produce:

- planned origins;
- planned included/excluded counts;
- expected output paths;
- overlap summary if previous metadata is supplied;
- validation failures for missing required fields.

### Runtime Smoke Test

If the MEDICC2 image is available, run the synthetic converted TSV through MEDICC2 and verify:

- summary exists;
- final tree exists;
- tree has expected sample labels plus artificial normal/diploid if inserted.

## Suggested Scripts Or Functions To Implement

The workflow can be one CLI with subcommands or separate scripts. Required functional boundaries:

```text
extract_cloneid_cn_profiles
audit_cell_selection_overlap
convert_wide_cn_to_medicc2
run_medicc2
collect_medicc2_outputs
plot_wgd_annotated_tree
write_run_manifest
```

Recommended CLI:

```bash
cloneid-medicc2-workflow extract --config config/analysis_config.yaml
cloneid-medicc2-workflow convert --config config/analysis_config.yaml
cloneid-medicc2-workflow validate --config config/analysis_config.yaml
cloneid-medicc2-workflow run --config config/analysis_config.yaml
cloneid-medicc2-workflow collect --config config/analysis_config.yaml
cloneid-medicc2-workflow plot-wgd --config config/analysis_config.yaml
cloneid-medicc2-workflow all --config config/analysis_config.yaml
```

Every command should support:

```text
--dry-run
--force false by default
--log-level
```

`--force` should be required to overwrite an existing analysis root.

## Deliverables

Each completed run should report:

- config file;
- selected cells file;
- wide copy-number matrix CSV;
- cell selection audit CSV;
- metadata JSON;
- optional QC heatmap PNG/PDF;
- pooled MEDICC2 TSV;
- per-group MEDICC2 TSVs when configured;
- MEDICC2 summaries;
- pairwise distances;
- branch lengths;
- final CN profiles;
- Newick/XML/PNG trees;
- event tables and heatmaps when `--events --plot heatmap` is used;
- WGD-annotated tree PNGs;
- overlap summary with previous run when prior metadata is supplied;
- run manifest with commands, image, platform, input hashes, timestamps, and software versions.

## Interpretation Guidance For The 4N Contamination Question

The manuscript-relevant readout should not rely on total ploidy alone. It should emphasize lineage-informative copy-number event patterns.

Support for 4N lineage origin:

- 4N-injected tumor cells cluster with 4N cell-line reference cells;
- 4N tumors share CN events private to the 4N reference and absent from the 2N reference;
- 4N tumor cells remain closer to 4N reference than 2N reference in pairwise MEDICC2 distances.

Support for 2N contamination/outgrowth:

- 4N-injected tumor cells cluster with the 2N cell-line reference;
- 4N tumors share CN events private to the 2N reference and absent from the 4N reference;
- only a small or absent subset of 4N-injected tumor cells falls near the 4N reference.

Ambiguous outcomes:

- 4N tumor cells fall between references;
- 2N and 4N references are not cleanly separated by copy-number profiles;
- tumor cells split across multiple clades;
- WGD placement is unstable across pooled versus group-specific analyses.

Ambiguous outcomes should trigger additional analyses rather than being forced into a binary contamination/no-contamination conclusion.
