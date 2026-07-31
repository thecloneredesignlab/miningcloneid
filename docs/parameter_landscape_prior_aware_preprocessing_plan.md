# Parameter Landscape Prior-Aware Preprocessing Plan

## Goal

Update the parameter-landscape dimensionality-reduction workflow so parameter distances can be interpreted relative to the optimizer search space, while keeping the existing empirical z-score UMAP workflow as the reference version.

The scientific question is whether fitted solutions occupy nearby regions of the allowed parameter volume and whether fitted solutions form visually or linearly separable regimes. For that question, each parameter should be measured after applying the same transform used for optimization and, in the rescaled sensitivity version, as a fraction of its allowed transformed range.

## Accepted Scope

This implementation scope is:

- keep the current empirical z-score parameter workflow unchanged as the reference version
- add prior-aware `[0, 1]` parameter rescaling as an additional major version
- add first-class PCA and t-SNE outputs beside UMAP
- retain one PCA-to-UMAP sensitivity workflow under `UMAPs/`
- keep pooled in vivo/in vitro objective coloring on the original objective scales
- preserve existing plot colors, themes, sizes, symbols, and legend style

Derived-only biological phenotype embeddings remain deferred.

## Directory Layout

Each dataset directory has parallel reduction folders:

```text
parameter_landscape/
  invivo/
    UMAPs/
    PCAs/
    TSNEs/
  invitro/
    UMAPs/
    PCAs/
    TSNEs/
  pooled_invivo_invitro/
    UMAPs/
    PCAs/
    TSNEs/
  FixO2Modes/
```

`UMAPs/` contains direct UMAP outputs and PCA-to-UMAP outputs. `PCAs/` contains true PCA coordinate plots with `PCA1/PCA2`. `TSNEs/` contains true t-SNE coordinate plots with `tSNE1/tSNE2` and axes labelled `t-SNE 1` and `t-SNE 2`.

## Preprocessing Modes

### zscore

This is the existing behavior and remains unchanged:

1. Apply the same log10 transform list already used by the parameter UMAP workflow.
2. Standardize each selected feature empirically with `scale()`.
3. Use the standardized matrix for UMAP, PCA, t-SNE, clustering from input features, and PCA-to-UMAP.

Existing z-score UMAP filenames are preserved.

### prior_unit

This is the new single-dataset rescaled version:

1. Read optimizer parameter bounds from the seed/input parameter table.
2. Use the optimizer transform encoded by the parameter table.
3. Transform the fitted or initial value on that optimizer scale.
4. Rescale by transformed lower and upper bounds:

```text
x_prior_unit = (g(x) - lower_g) / (upper_g - lower_g)
```

No empirical z-score scaling is applied after prior-unit scaling. PCA centers the prior-unit matrix internally before computing PCs.

## Pooled Scaling Rules

Both pooled prior-unit rules are implemented.

### context_prior_unit

In vivo rows are rescaled by in vivo optimizer bounds. In vitro rows are rescaled by in vitro optimizer bounds.

This answers: where does each result sit within its own model-fitting prior volume?

### common_prior_unit

For each shared parameter, the common transformed lower bound is the minimum of the in vivo and in vitro transformed lower bounds, and the common transformed upper bound is the maximum of the two transformed upper bounds.

Both in vivo and in vitro rows are then rescaled using this shared transformed interval.

This answers: how far apart are in vivo and in vitro rows on one common transformed parameter scale?

## Reductions

The workflow supports:

- `umap`: direct UMAP, written under `UMAPs/`
- `pca`: true PCA coordinate plots, written under `PCAs/`
- `tsne`: true t-SNE coordinate plots, written under `TSNEs/`
- `pca_umap`: PCA-to-UMAP sensitivity output, written under `UMAPs/`

PCA-to-UMAP is retained as a robustness check, not as the replacement for true PCA plots.

## PCA Interpretation

PCA-space plots are useful because they show a linear projection whose axes and variance fractions can be reported directly. For the prior-unit version, PCA is run after prior-range rescaling and centering, so the loadings describe variation relative to the allowed transformed prior volume.

UMAP and t-SNE remain nonlinear visualization methods. Cluster labels estimated from their final two-dimensional coordinates should be interpreted as visual annotations. Clustering from the input feature matrix remains available as a companion analysis.

## t-SNE Guard

Full initial-plus-best t-SNE can be computationally expensive for large pooled tables. The default workflow skips full t-SNE unless `run_full_tsne=TRUE`.

Sampled and best-only t-SNE outputs are still generated when `tsne` is requested.

## Metadata

Each generated coordinate table receives a companion preprocessing metadata table.

For z-score outputs, the table records empirical center and scale.

For prior-unit outputs, the table also records:

- dataset context
- natural parameter name
- transformed optimizer parameter name
- transform
- transformed lower and upper bounds
- natural lower and upper bounds

## Report Behavior

The report generator remains compatible with existing UMAP outputs and now understands reduction-specific clustered-output folders:

- `ByUMAPCoordinates`
- `ByPCACoordinates`
- `ByTSNECoordinates`
- `ByInputFeatures`

When corresponding clustered outputs exist, the report can include best-only UMAP, PCA, t-SNE, prior-unit, and PCA-to-UMAP parameter panels.

## Validation

Validation should include:

1. R parse/source checks for all parameter-landscape scripts.
2. Bash syntax check for the Slurm submitter.
3. Smoke test of z-score and prior-unit PCA/t-SNE.
4. Smoke test of PCA-to-UMAP.
5. Smoke test of pooled `context_prior_unit` and `common_prior_unit`.
6. Header checks for PCA/t-SNE coordinate tables and pooled distance tables.
7. Confirmation that pooled objective columns remain raw objective values, not 0-1 normalized values.

## Deferred Work

Derived-only biological phenotype embeddings should be implemented separately from raw-parameter embeddings. If raw parameters and derived features are combined later, the feature blocks should be weighted explicitly so the result is not determined only by the number of columns in each block.
