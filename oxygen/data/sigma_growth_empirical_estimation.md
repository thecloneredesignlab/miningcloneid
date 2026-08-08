# Empirical estimation of `sigma_growth`

## Purpose

`sigma_growth` is the shared, constant standard deviation of the in vitro
measurement-day log viable-cell-count error. The absolute viable-cell-count
likelihood is

\[
\log C_{jk}^{\mathrm{obs}} \sim
\mathcal{N}\left(
  \log C_{jk}^{\mathrm{model}},
  \sigma_{\mathrm{growth}}^2
\right),
\]

where every experimental post-seeding count receives the same weight for the
same log-count residual, irrespective of its measurement day. Equivalently,
the same observed-to-predicted fold error has the same likelihood penalty at
an early or late timepoint.

Day 0 defines the model initial condition and is not entered again as a
likelihood observation. Timepoint log-likelihoods are averaged within each
passage before the existing lineage, cohort, and global averaging, so passages
with more intermediate counts do not receive more objective weight solely
because they were sampled more often. The stored endpoint growth-rate value
`g` is derived from the same counts and is not added as a separate likelihood.

This replaces the previous model

\[
\operatorname{SD}\{\log C(t)\}=\sigma_{\mathrm{growth}}t,
\]

which gave the same log-count residual less weight at later measurement days.
Consequently, `sigma_growth` is now dimensionless rather than day\(^{-1}\).

## Data source and pairing

The estimate uses:

```text
oxygen/data/metadata.csv
```

For each formal O1/O2 passage, the day-0 count is used to compute the observed
log fold change at every post-seeding measurement:

\[
y_{j\ell t}=\log\left(\frac{C_{j\ell t}}{C_{j\ell 0}}\right),
\]

where \(j\) indexes cohort and passage, \(\ell\) is O1 or O2, and \(t\) is a
measurement day. O1 and O2 are paired only when cohort, passage number, and
measurement day all match. This produces:

- 95 matched post-seeding timepoints;
- 44 matched cohort-passage pairs;
- measurement days 1, 2, 3, 4, 5, 6, 7, 9, 10, 12, and 14.

## Estimator

For each matched timepoint, define

\[
d_{jt}=y_{j,\mathrm{O1},t}-y_{j,\mathrm{O2},t}.
\]

If O1 and O2 have independent errors with a shared constant log-count standard
deviation, then

\[
\operatorname{Var}(d_{jt})=2\sigma_{\mathrm{growth}}^2,
\]

and the pooled estimator is

\[
\widehat{\sigma}_{\mathrm{growth}}=
\frac{\operatorname{SD}(d_{jt})}{\sqrt{2}}.
\]

The uncertainty interval is obtained by resampling the 44 cohort-passage pairs
as clusters, preserving repeated timepoints from the same passage.

## Results

| Dataset | Matched timepoints | Estimated `sigma_growth` |
|---|---:|---:|
| 2N | 48 | 0.153560 |
| 4N | 47 | 0.239822 |
| Combined | 95 | 0.201082 |

Additional combined diagnostics:

- mean O1 minus O2 log-fold-change difference: 0.015602;
- robust MAD-based estimate: 0.128525;
- passage-cluster bootstrap 95% interval: 0.144948 to 0.252659;
- bootstrap median: 0.196174.

The day-stratified estimates are not monotone in time. A Gaussian likelihood
comparison using the same 95 paired differences gives:

| Error-scale model | Fitted scale | AIC |
|---|---:|---:|
| Constant SD | 0.200014 | 33.674 |
| SD proportional to \(\sqrt{t}\) | 0.104772 | 37.026 |
| SD proportional to \(t\) | 0.063088 | 66.867 |

The constant-SD model has the lowest AIC, while the former linear-time model is
strongly disfavored by these paired data.

## Parameter-table values

The fitted parameter uses a rounded empirical center and deliberately wider
bounds:

```text
init_value  = 0.200
lower_bound = 0.080
upper_bound = 0.400
```

The lower bound lies below the robust estimate. The upper bound exceeds the
cluster-bootstrap interval while preventing `sigma_growth` from expanding
without limit to absorb structural model error.

For interpretation:

- `sigma_growth = 0.20` corresponds to a one-standard-deviation multiplicative
  factor of `exp(0.20) = 1.22` at every measurement day;
- `sigma_growth = 0.08` corresponds to a factor of 1.08;
- `sigma_growth = 0.40` corresponds to a factor of 1.49.

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

day_zero <- metadata[
  metadata$day == 0,
  c("cohort", "lineage", "passage", "count")
]
names(day_zero)[[4]] <- "count0"
metadata <- merge(
  metadata,
  day_zero,
  by = c("cohort", "lineage", "passage")
)
metadata <- metadata[metadata$day > 0, ]
metadata$log_fold_change <- log(metadata$count / metadata$count0)

o1 <- metadata[
  metadata$lineage == "O1",
  c("cohort", "passage", "day", "log_fold_change")
]
o2 <- metadata[
  metadata$lineage == "O2",
  c("cohort", "passage", "day", "log_fold_change")
]
names(o1)[[4]] <- "log_fold_change_O1"
names(o2)[[4]] <- "log_fold_change_O2"
paired <- merge(o1, o2, by = c("cohort", "passage", "day"))
paired$difference <-
  paired$log_fold_change_O1 - paired$log_fold_change_O2
paired$pair_id <- paste(paired$cohort, paired$passage, sep = "_A")

sigma_hat <- stats::sd(paired$difference) / sqrt(2)
sigma_robust <- stats::mad(paired$difference) / sqrt(2)
sigma_by_cohort <- vapply(
  split(paired$difference, paired$cohort),
  function(x) stats::sd(x) / sqrt(2),
  numeric(1)
)

set.seed(5826)
pair_ids <- unique(paired$pair_id)
bootstrap_sigma <- replicate(10000, {
  sampled_ids <- sample(pair_ids, length(pair_ids), replace = TRUE)
  sampled_difference <- unlist(lapply(sampled_ids, function(pair_id) {
    paired$difference[paired$pair_id == pair_id]
  }), use.names = FALSE)
  stats::sd(sampled_difference) / sqrt(2)
})
bootstrap_interval <- stats::quantile(
  bootstrap_sigma,
  c(0.025, 0.5, 0.975)
)

fit_scale_model <- function(scale) {
  scale <- as.numeric(scale)
  fit <- stats::optim(
    c(mean(paired$difference), log(0.2)),
    function(theta) {
      -sum(stats::dnorm(
        paired$difference,
        mean = theta[[1]],
        sd = sqrt(2) * exp(theta[[2]]) * scale,
        log = TRUE
      ))
    }
  )
  c(
    mean_difference = fit$par[[1]],
    sigma = exp(fit$par[[2]]),
    negative_log_likelihood = fit$value,
    AIC = 2 * fit$value + 2 * length(fit$par)
  )
}

scale_comparison <- rbind(
  constant = fit_scale_model(rep(1, nrow(paired))),
  sqrt_time = fit_scale_model(sqrt(paired$day)),
  linear_time = fit_scale_model(paired$day)
)

sigma_hat
sigma_robust
sigma_by_cohort
bootstrap_interval
scale_comparison
```

## Limitations

This estimate combines cell-count measurement variation, day-0 count
uncertainty, and true O1/O2 biological-lineage variation. It is therefore an
effective observation-error estimate for the current deterministic
shared-trajectory model, not a pure technical counting-error estimate.
Repeated independent counts from the same culture at each measurement day
would be required to estimate technical error separately.
