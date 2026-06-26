# Parameter Landscape Prior-Aware Preprocessing Implementation Plan

## Goal

Update the parameter-landscape PCA/UMAP workflow so distances are interpreted relative to the optimizer search space, not relative to the empirical spread of whichever rows happen to be plotted.

The main scientific question is whether fitted solutions occupy the same parameter or biological regime, whether they sit near one another within the allowed prior volume, and whether the optimizer moved from random initial populations into constrained regions of that volume. For that question, each parameter should be measured as displacement across its allowed prior range after applying the same transform used for optimization.

## Current Behavior

The current `parameter_landscape_clustering` workflow already handles raw unit differences partially:

- positive-scale parameters listed in `umap_log10_parameter_set()` are transformed as `log10(x)`
- all selected columns are then standardized as empirical z-scores with `scale()`
- PCA is run on the already-standardized matrix with no additional centering or scaling
- UMAP is run either directly on the standardized features or on the retained PCA scores

This is better than running UMAP on raw natural-scale parameter values. However, it does not explicitly use parameter prior bounds. As a result, a parameter can appear important because it varies strongly within the plotted rows even if that variation is small relative to its allowed search range. Conversely, a parameter can appear unimportant if the plotted rows have broad empirical variance, even when the fitted values occupy a biologically meaningful part of the prior.

This issue is most important for best-only embeddings, where empirical z-scoring is based only on the fitted solutions. It is also relevant for initial-plus-best embeddings because the result depends on how many initial points are included and whether those points truly represent the prior distribution.

## Recommended Preprocessing

Use prior-aware transformed coordinates as the primary feature matrix.

For each parameter `j`:

1. Select the same parameter set currently used for the corresponding analysis.
2. Apply the appropriate transform `g_j()` to the parameter value.
3. Apply the same transform to the lower and upper prior bounds.
4. Normalize the transformed value by the transformed prior width.

Recommended formula:

```text
x_prior_unit_ij =
  (g_j(x_ij) - g_j(lower_j)) /
  (g_j(upper_j) - g_j(lower_j))
```

Optionally center this value for distance calculations:

```text
x_prior_centered_ij = x_prior_unit_ij - 0.5
```

This makes each feature dimension comparable as a fraction of its allowed search interval. A distance of `0.1` in one parameter then has the same broad interpretation as a distance of `0.1` in another: the solution moved about 10% of that parameter's allowed transformed range.

## Transform Rules

Use transforms that match the optimizer scale or the biologically meaningful distance scale.

Recommended default rules:

- use `log10(x)` for positive rate, scale, and threshold parameters spanning orders of magnitude
- use `log10(x)` for very small probabilities where multiplicative change is the meaningful comparison
- use `logit(x)` for probabilities or fractions where both low and high values are plausible and distance near 0 and 1 should be treated symmetrically
- use identity transform for bounded shape parameters where additive differences are meaningful on the natural scale

The implementation should avoid maintaining a separate ad hoc transform list for the UMAP workflow if the optimization backends already encode parameter transforms. Prefer reading or reusing the optimizer parameter table so plotted coordinates and optimizer coordinates remain synchronized.

## Primary Feature Spaces

Implement and report at least three feature-space modes.

### 1. Raw Parameter Prior Space

This should become the primary UMAP/PCA input for raw parameter analyses.

Features:

- selected fitted parameters only
- transformed by parameter-specific optimizer/biology transform
- normalized by transformed prior range
- no empirical z-scoring as the primary scaling step

This is the best answer to: "Where are the fitted solutions within the allowed prior volume?"

### 2. Derived Biological Phenotype Space

Use model-implied biological features as a separate embedding, not only appended to raw parameters.

Candidate features:

- effective proliferation rate at representative O2 levels for 2N and 4N
- effective death rate at representative O2 levels for 2N and 4N
- net growth rate
- turnover rate
- 4N relative fitness under hypoxia
- O2 sensitivity thresholds
- effective missegregation burden under low O2
- viability after missegregation for 2N and 4N daughters
- oxygen supply regime at representative tumor burdens

Scale these features using biologically interpretable units where possible. If a feature does not have a clear prior range, use robust scaling such as median and MAD, and document that this feature space is phenotype-relative rather than prior-relative.

This is the best answer to: "Do different parameterizations imply the same biological behavior?"

### 3. Combined Parameter Plus Phenotype Space

Keep this as a sensitivity analysis, not as the primary embedding.

If raw parameters and derived features are combined, assign explicit group weights so the result is not determined merely by whichever group contributes more columns.

Example:

```text
combined_distance =
  0.5 * parameter_prior_space_distance +
  0.5 * phenotype_space_distance
```

Operationally this can be implemented by scaling each block so the total block contribution is controlled before concatenation.

## PCA Recommendation

Do not make PCA-to-UMAP the primary visualization for the current 14-18 dimensional parameter spaces. Direct UMAP on prior-normalized transformed parameters should be the main result because the input dimensionality is modest and each coordinate has a clear interpretation.

Use PCA-to-UMAP as a sensitivity check.

If PCA is used:

- run PCA on the prior-normalized transformed matrix
- avoid a fixed `pca_n = 10` as the primary rule
- retain enough PCs to explain a high variance threshold, for example 95-99%
- also allow retaining all PCs when the goal is only orthogonalization/noise smoothing
- report the number of retained PCs and cumulative variance in the output table and figure caption

For best-only fits, be cautious about discarding low-variance PCs. A low-variance parameter axis can still distinguish local optima if all fitted values are tightly clustered except for a few biologically meaningful directions.

## Implementation Tasks

### A. Add Prior-Aware Feature Builder

Add a helper in `parameter_landscape_utils.R` that builds a feature matrix from:

- a data frame of parameter values
- a parameter set
- a parameter metadata table containing lower bound, upper bound, and transform

The helper should:

- validate that every selected parameter has a finite lower and upper bound
- validate that transformed lower and upper bounds are finite
- validate that transformed upper is greater than transformed lower
- transform values and bounds consistently
- return prior-unit or centered-prior-unit coordinates
- return an attributes table describing transform, lower, upper, transformed lower, transformed upper, and transformed width

Suggested function names:

- `parameter_prior_transform_metadata(dataset, seed_dir = NULL)`
- `transform_umap_features_prior_unit(df, params, metadata, center = TRUE)`

### B. Preserve Current Z-Score Workflow As Legacy/Sensitivity Mode

Do not remove the existing empirical z-score workflow immediately. Add a preprocessing option with explicit names:

- `preprocess = "prior_unit"` for the new primary mode
- `preprocess = "zscore"` for the current behavior
- optionally `preprocess = "prior_unit_zscore"` only as a diagnostic, not as the default

The output filenames and reduction labels should include the preprocessing mode so old and new figures cannot be confused.

### C. Apply To Direct UMAP, PCA-UMAP, And Clustering

Use the same feature matrix for:

- direct UMAP
- PCA input
- input-feature clustering

Do not use one preprocessing mode for embedding and another for cluster assignment unless the output name and caption state that explicitly.

UMAP-coordinate clustering should remain secondary and should be described as a visual annotation, not as evidence of real parameter-space clusters.

### D. Add Best-Only And Initial-Plus-Best Variants

Generate both:

- initial-plus-best prior-space UMAP
- best-only prior-space UMAP

The initial-plus-best version answers whether fitted solutions move from random prior samples toward a constrained region. The best-only version answers whether fitted solutions are mutually close or form separated regimes.

For initial-plus-best plots, keep sampling behavior explicit. Record the number of initial rows actually used, whether sampling was per seed, and the number of best-fit rows.

### E. Add Derived-Only Biological Embeddings

Add a mode that uses only derived biological features. This should not include the raw parameters unless explicitly requested.

For derived features, write a companion metadata table that records:

- feature definition
- units
- O2 level or tumor burden used, if applicable
- scaling rule
- whether the feature was computed from in vivo, in vitro, or pooled context

This metadata is necessary because derived feature distances are otherwise harder to interpret than prior-unit parameter distances.

### F. Update Reports And Captions

Update generated report text so every UMAP/PCA panel states:

- input feature set
- preprocessing mode
- whether lower/upper prior bounds were used
- transform applied to each parameter class
- whether PCA was used before UMAP
- if PCA was used, number of PCs retained and cumulative variance
- whether clusters were estimated from input features or from UMAP coordinates

The methods language should emphasize that prior-unit preprocessing makes distances comparable across parameters as fractions of the allowed transformed prior range.

## Validation Plan

Add lightweight validation before running large plots:

1. Parse all scripts in `parameter_landscape_clustering`.
2. Unit-test transform behavior on a small synthetic metadata table.
3. Verify that a value at the lower bound maps to `0`, midpoint maps near `0.5` on identity-scale parameters, and upper bound maps to `1`.
4. Verify that log-scale parameters are normalized after log transformation, not before.
5. Verify that non-positive values fail for log-scale parameters.
6. Verify that invalid bounds fail loudly.
7. Verify that PCA variance output is produced from the prior-unit matrix.
8. Run a small smoke test using a tiny fixture table, generating direct UMAP and PCA-UMAP outputs.

For real-result runs, write the preprocessing metadata table next to the UMAP coordinate table. The metadata table should make the figure reproducible without reading the code.

## Acceptance Criteria

The implementation is complete when:

- prior-aware preprocessing is available and is the default for parameter-space UMAPs
- legacy z-score preprocessing remains available as an explicit sensitivity option
- direct UMAP and PCA-UMAP both use the selected preprocessing consistently
- PCA retention can be variance-threshold based
- best-only and initial-plus-best figures include preprocessing mode in filenames or labels
- derived-only biological phenotype UMAPs can be generated separately from raw-parameter UMAPs
- every output includes a metadata table listing parameter transforms and prior bounds
- validation tests cover transform and bound handling

## Interpretation Guidance

Use direct UMAP on prior-normalized transformed parameters as the primary parameter-space visualization.

Use PCA-to-UMAP only as a robustness check. If direct UMAP and PCA-to-UMAP agree, that supports the stability of the visual structure. If they disagree, inspect which PCs were discarded and whether low-variance directions separate fitted regimes.

Use derived-only phenotype UMAPs to answer whether apparently different parameter solutions imply materially different biology. Raw parameter UMAPs and biological phenotype UMAPs should be interpreted together: separated parameter regimes that collapse in phenotype space may represent equivalent biological solutions, whereas separation in both spaces is stronger evidence for distinct biological optima.
