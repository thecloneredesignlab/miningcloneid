# Figure 5F posterior-target audit

## Scope

This audit fixes the evidence boundary for the requested three-distribution
comparison in a single Figure 5F: prior, posterior, and optimizer endpoints.
It uses the joint-fit implementation at miningcloneid commit
`83953a874401e42cd176432786f889a896adc959` and treats the 500-seed joint-fit
result tree as read-only.

## Verified optimizer target

The saved joint objective is replayable from the historical code and saved
inputs. It combines:

1. an in-vivo empirical loss;
2. an in-vitro empirical loss;
3. the configured in-vivo soft-prior penalty;
4. the Welsch soft-coupling penalty; and
5. any enabled constraint penalty (disabled in the fitted runs used here).

The saved configuration uses unit in-vivo and in-vitro weights, unit weights
for the three in-vitro data modalities, `lambda_prior = 0.03`, and a Welsch
coupling constant of 0.4. The selected runs use 14 paired parameters and hard
parameter bounds.

## Why the exact optimizer target is not a strict likelihood posterior

The in-vivo C++ objective computes valid pointwise negative log-likelihood
terms, but the optimized burden loss is the mean of per-tumor mean losses.
Terminal ploidy is averaged first within each tumor and then balanced equally
between the 2N and 4N cohorts. Necrosis uses an averaged squared standardized
logit residual. These normalized losses are not the sum of the observation-level
log likelihoods.

The in-vitro objective likewise computes growth and karyotype log-likelihood
contributions, but averages them separately by passage before assigning equal
modality weights. The flow term is the cross-entropy of a smoothed empirical
G0/G1 density against the model density. The saved flow-density table does not
define an event-level sampling model or an effective event count for that
cross-entropy term.

Therefore, sampling a density proportional to `exp(-objective)` is a valid
bounded generalized/Gibbs posterior for the fitted composite loss, but it is
not the strict Bayesian posterior from a product of the recorded observation
likelihoods. Optimizer endpoints also remain optimizer endpoints and must not
be renamed posterior draws.

## Scientifically valid release choices

### A. Exact-target generalized posterior (recommended for direct comparison)

Use the exact replayed joint objective as a Gibbs loss. Plot `Prior`,
`Generalized posterior`, and `Optimizer endpoints` in every parameter-family
facet. This preserves direct comparability to all 500 optimizer endpoints and
does not change the fitted target.

### B. Strict posterior with flow held out of conditioning

Build a conventional likelihood from the raw in-vivo burden, terminal
karyotype, necrosis, in-vitro growth, and in-vitro karyotype contributions.
Use the flow-density data only for posterior-predictive validation. This is a
strict posterior, but it is not conditioned on exactly the same data target as
the optimizer endpoints, so differences cannot be attributed solely to
uncertainty or identifiability.

### C. Strict posterior with a new flow observation model

Specify and justify an event- or count-level observation model for calibrated
flow data, including its effective sample size, smoothing, gating, and
overdispersion. This would permit a strict full-data posterior, but it is a new
methodological assumption not present in the frozen joint-fitting code or
saved density-table contract.

## Frozen plotting contract after target selection

Release choice: **A. Exact-target generalized posterior**. The panel and all
supporting text must use the label `Generalized posterior`; the shorter label
`Bayesian posterior` is not permitted for this composite-loss target.

The final panel remains one Figure 5F. For each of the 14 paired parameters,
C01, C02, and C03 are displayed in separate aligned facets. Each family uses
the one warm-start pair selected by the lowest separate in-vivo objective:

- C01: `tsne_vi_seed366_C01Sc01_vt_seed10`;
- C02: `tsne_vi_seed25_C02Sc01_vt_seed10`;
- C03: `tsne_vi_seed311_C03Sc02_vt_seed10`.

The x axis is `log2(in vivo / in vitro)`. Prior is a neutral gray reference,
the posterior is a filled distribution, and the 500 optimizer endpoints for
the selected pair are shown as a distinct dot-dash outline with their median
and 5th--95th percentile span. That span is descriptive rather than an
inferential interval. Distribution type is encoded by line/fill and point
shape as well as color, and cluster identities remain C01, C02, and C03
throughout Figure 5.
