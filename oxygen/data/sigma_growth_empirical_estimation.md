# Empirical context for `sigma_growth`

## Purpose

`sigma_growth` is the shared standard deviation of the in vitro
passage-average net growth-rate likelihood, in day\(^{-1}\). For passage \(j\),
all valid post-seeding count times are first reduced to one rate by fitting a
zero-intercept line to log live-cell fold change:

\[
y_{jt}=\log\left(\frac{C_{jt}}{C_{j0}}\right), \qquad
\widehat g_j=\frac{\sum_{t>0}t\,y_{jt}}{\sum_{t>0}t^2}.
\]

The growth observation model is

\[
\widehat g_j^{\mathrm{obs}} \sim
\mathcal{N}\left(
  \widehat g_j^{\mathrm{model}},
  \sigma_{\mathrm{growth}}^2
\right).
\]

Day 0 supplies the denominator and is not an additional likelihood unit. Every
positive measurement time, including the final experimental count, contributes
to the slope. Each passage then contributes exactly one likelihood value before
the existing lineage, cohort, and global averaging. Therefore, passages with
two and five positive count times have equal objective weight.

The individual measured counts are still exported as diagnostics, but they are
not independently scored as absolute-count likelihood units.

## Empirical scale check

The scale check uses:

```text
oxygen/data/metadata.csv
```

For each formal O1/O2 passage, a passage-average rate is computed with the same
estimator used by the objective. O1 and O2 rates are paired by cohort and
passage number. This gives 88 lineage-passage rates and 44 matched O1/O2 pairs.

For matched pair \(j\), define

\[
d_j=\widehat g_{j,\mathrm{O1}}-\widehat g_{j,\mathrm{O2}}.
\]

If the two lineages have independent rate errors with the same standard
deviation, a lower-bound scale estimate is

\[
\widehat\sigma_{\mathrm{growth}}=
\frac{\operatorname{SD}(d_j)}{\sqrt{2}}.
\]

The resulting estimates are:

| Dataset | Matched passage pairs | Estimated `sigma_growth` (day\(^{-1}\)) |
|---|---:|---:|
| 2N | 22 | 0.033420 |
| 4N | 22 | 0.037161 |
| Combined | 44 | 0.035284 |

Additional combined diagnostics are:

- mean O1 minus O2 growth-rate difference: -0.000838 day\(^{-1}\);
- robust MAD-based estimate: 0.026422 day\(^{-1}\);
- passage-pair bootstrap 95% interval: 0.026830 to 0.042443 day\(^{-1}\);
- bootstrap median: 0.034672 day\(^{-1}\).

This paired estimate is a lower-bound reference rather than the fitted value:
pairing removes cohort-passage effects shared by O1 and O2, while the
deterministic model must also absorb unmodeled passage-specific biological and
protocol variation.

## Parameter-table values

The fitted parameter retains a broad range so the model can estimate the
effective passage-level discrepancy:

```text
init_value  = 0.200 day^-1
lower_bound = 0.050 day^-1
upper_bound = 0.500 day^-1
```

The lower bound is already above the paired-lineage estimate. The upper bound
allows substantial passage-to-passage variation without making the growth
component unboundedly weak.

## Reproducible calculation

Run from the repository root:

```r
metadata <- read.csv(
  "oxygen/data/metadata.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
metadata <- metadata[
  grepl(
    "_NLS_(2N|4N)_(O1|O2)_A[0-9]+_seed$",
    metadata$passage_id
  ),
  c("passage_id", "num_date", "correctedCount")
]

parts <- regmatches(
  metadata$passage_id,
  regexec(
    "^SUM-159_NLS_(2N|4N)_(O1|O2)_A([0-9]+)_seed$",
    metadata$passage_id
  )
)
metadata$cohort <- vapply(parts, `[`, character(1), 2)
metadata$lineage <- vapply(parts, `[`, character(1), 3)
metadata$passage <- as.integer(vapply(parts, `[`, character(1), 4))
metadata$day <- as.numeric(metadata$num_date)
metadata$count <- as.numeric(metadata$correctedCount)

passage_key <- interaction(
  metadata$cohort,
  metadata$lineage,
  metadata$passage,
  drop = TRUE
)
rates <- do.call(rbind, lapply(split(metadata, passage_key), function(d) {
  day_zero <- d[d$day == 0 & is.finite(d$count) & d$count > 0, , drop = FALSE]
  positive <- d[
    d$day > 0 & is.finite(d$day) & is.finite(d$count) & d$count > 0,
    ,
    drop = FALSE
  ]
  if (nrow(day_zero) != 1L || nrow(positive) == 0L) return(NULL)
  data.frame(
    cohort = d$cohort[[1]],
    lineage = d$lineage[[1]],
    passage = d$passage[[1]],
    rate = sum(
      positive$day * log(positive$count / day_zero$count)
    ) / sum(positive$day^2)
  )
}))

o1 <- rates[rates$lineage == "O1", c("cohort", "passage", "rate")]
o2 <- rates[rates$lineage == "O2", c("cohort", "passage", "rate")]
names(o1)[[3]] <- "rate_O1"
names(o2)[[3]] <- "rate_O2"
paired <- merge(o1, o2, by = c("cohort", "passage"))
paired$difference <- paired$rate_O1 - paired$rate_O2

sigma_hat <- stats::sd(paired$difference) / sqrt(2)
sigma_robust <- stats::mad(paired$difference) / sqrt(2)
sigma_by_cohort <- vapply(
  split(paired$difference, paired$cohort),
  function(x) stats::sd(x) / sqrt(2),
  numeric(1)
)

set.seed(5826)
pair_rows <- seq_len(nrow(paired))
bootstrap_sigma <- replicate(10000, {
  sampled_rows <- sample(pair_rows, length(pair_rows), replace = TRUE)
  stats::sd(paired$difference[sampled_rows]) / sqrt(2)
})
bootstrap_interval <- stats::quantile(
  bootstrap_sigma,
  c(0.025, 0.5, 0.975)
)

sigma_hat
sigma_robust
sigma_by_cohort
bootstrap_interval
```

## Limitations

The fitted uncertainty combines count measurement error, day-0 count
uncertainty, true lineage variation, and deterministic-model discrepancy. It is
not a pure technical counting-error estimate. Repeated independent counts from
the same culture and timepoint would be required to identify technical error
separately.
