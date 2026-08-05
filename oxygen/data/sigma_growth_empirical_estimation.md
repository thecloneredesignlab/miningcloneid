# Empirical estimation of `sigma_growth`

## Purpose

`sigma_growth` is the shared per-day standard deviation of the in vitro
log-growth observation error. In the current absolute viable-cell-count
likelihood, the per-day uncertainty is propagated separately to every
experimental post-seeding count as

\[
\log C_{jk}^{\mathrm{obs}} \sim
\mathcal{N}\left(
  \log C_{jk}^{\mathrm{model}},
  (\sigma_{\mathrm{growth}} t_{jk})^2
\right),
\]

where \(t_{jk}>0\) is the measured day within passage \(j\). Day 0 defines the
model's initial condition and is not entered again as a likelihood
observation. Thus, `sigma_growth` retains units of per-day log growth rather
than being a direct cell-count standard deviation.

The timepoint log-likelihoods are first averaged within each passage. The
existing lineage, cohort, and global averaging is then applied to those
passage means, so passages with more intermediate measurements do not receive
more objective weight than passages with only an endpoint count. The stored
`g` value is derived from the same cell counts and is not added as a second
likelihood term.

The likelihood timepoints are loaded from `oxygen/data/metadata.csv` using
`passage_id`, `num_date`, and `correctedCount`. The current table supplies 217
post-seeding counts for 112 formal passages. The O1 and O2 A23 passages are
absent from that table, so their existing `fit_data.Rds` endpoint counts are
used explicitly as fallbacks, giving 219 observations across 114 passages.

## Data source

The estimate uses:

```text
oxygen/ploidyOxygen/data/fit_objects/fit_data.Rds
```

The stored `g` field is the experimental per-day log growth rate. O1 and O2
provide paired biological lineages at the same cohort and lineage passage:

- 23 paired passages for the 2N cohort;
- 22 paired passages for the 4N cohort;
- 45 paired passages in total.

Control lineages are not included because the current data contain no matched
control replicate at each passage.

## Estimator

For each matched cohort and passage, define

\[
d_j = g_{j,\mathrm{O1}} - g_{j,\mathrm{O2}}.
\]

If the two lineages have independent errors with a common standard deviation
\(\sigma_{\mathrm{growth}}\), then

\[
\operatorname{Var}(d_j) = 2\sigma_{\mathrm{growth}}^2,
\]

and the paired estimator is

\[
\widehat{\sigma}_{\mathrm{growth}} =
\frac{\operatorname{SD}(d_j)}{\sqrt{2}}.
\]

## Results

| Dataset | Number of pairs | Estimated `sigma_growth` (day\(^{-1}\)) |
|---|---:|---:|
| 2N | 23 | 0.031753 |
| 4N | 22 | 0.033106 |
| Combined | 45 | 0.032100 |

Additional diagnostics from the combined pairs:

- mean O1 minus O2 difference: -0.000180 day\(^{-1}\);
- robust MAD-based estimate: 0.021384 day\(^{-1}\);
- approximate normal-theory 95% interval: 0.026574 to 0.040549 day\(^{-1}\);
- O1/O2 growth-rate correlation: 0.921808.

The interval is approximate because sequential passage pairs within a lineage
are not independent biological replicates.

## Parameter-table values

The fitted parameter uses a rounded empirical center and deliberately wider
optimization bounds:

```text
init_value  = 0.032
lower_bound = 0.015
upper_bound = 0.080
```

The lower bound extends below the robust estimate. The upper bound extends
above the approximate confidence interval to allow biological variation and
serial dependence, while preventing `sigma_growth` from expanding enough to
absorb large structural model errors. These values replace the previous broad
range of 0.001 to 1.

For interpretation, at a five-day observation:

- `sigma_growth = 0.032` gives a log-count standard deviation of 0.16 and a
  one-standard-deviation multiplicative factor of `exp(0.16) = 1.17`;
- the previous seed10 estimate of 0.5069 gives a log-count standard deviation
  of approximately 2.53 and a multiplicative factor of approximately 12.6.

## Reproducible calculation

Run from the repository root:

```r
fit_data <- readRDS(
  "oxygen/ploidyOxygen/data/fit_objects/fit_data.Rds"
)

ids <- names(fit_data)
keep <- grepl("_(2N|4N)_(O1|O2)_A[0-9]+_seed$", ids)
ids <- ids[keep]
parts <- strsplit(ids, "_", fixed = TRUE)

growth <- data.frame(
  sample_id = ids,
  g = vapply(fit_data[ids], function(x) as.numeric(x$g), numeric(1)),
  cohort = vapply(parts, `[`, character(1), 3),
  lineage = vapply(parts, `[`, character(1), 4),
  passage = as.integer(sub(
    "A", "", vapply(parts, `[`, character(1), 5), fixed = TRUE
  )),
  stringsAsFactors = FALSE
)

o1 <- growth[growth$lineage == "O1", c("cohort", "passage", "g")]
o2 <- growth[growth$lineage == "O2", c("cohort", "passage", "g")]
paired <- merge(
  o1,
  o2,
  by = c("cohort", "passage"),
  suffixes = c("_O1", "_O2")
)
paired$difference <- paired$g_O1 - paired$g_O2

estimate_sigma <- function(difference) {
  stats::sd(difference) / sqrt(2)
}

estimate_sigma(paired$difference)
estimate_sigma(paired$difference[paired$cohort == "2N"])
estimate_sigma(paired$difference[paired$cohort == "4N"])

robust_sigma <- stats::mad(paired$difference) / sqrt(2)
n_pairs <- nrow(paired)
sigma_hat <- estimate_sigma(paired$difference)
approximate_ci <- sqrt(
  (n_pairs - 1) * sigma_hat^2 /
    stats::qchisq(c(0.975, 0.025), df = n_pairs - 1)
)
```

## Limitations

This estimate combines measurement variation and true O1/O2 biological-lineage
variation. It is therefore an effective observation-error estimate for the
current deterministic shared-trajectory model, not a pure technical counting
error estimate. Repeated seed and harvest counts from the same culture would be
required to estimate technical counting error separately.
