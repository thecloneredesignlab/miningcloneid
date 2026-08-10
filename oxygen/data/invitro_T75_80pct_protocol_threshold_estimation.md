# In-vitro T75 80% confluence protocol thresholds

## Purpose

The in-vitro lineage simulator needs a passage-feasibility threshold that is
independent of the observed final viable-cell count from each individual
culture cycle. The cultures used T75 flasks and were passaged at approximately
80% confluence. This file records the fixed cohort-level live-cell thresholds
used to approximate that protocol.

These thresholds are fixed preprocessing constants. They are not fitted model
parameters, likelihood observations, or passage-specific targets.

## Source data and inclusion rule

The calibration used the 114 formal SUM-159 LTEE culture cycles reconstructed
from `oxygen/data/metadata.csv`. The same rows are materialized as
`data/InVitroData_LTEE/invitro_passage_observations.tsv` in the companion
`HypoxiaLTEEFigures` analysis repository. Included passage identifiers matched
`SUM-159_NLS_[24]N_(C|O1|O2)_A[0-9]+_seed` and had a finite positive
`final_cells` value.

The observed endpoint count is treated as a noisy measurement of the cell
number near the approximate 80% confluence operation point. Robust cohort
medians were used so that a single unusual passage cannot define the protocol
threshold.

| Cohort | Passages | Median final cells | Q25 | Q75 | Rounded protocol threshold |
|---|---:|---:|---:|---:|---:|
| 2N | 58 | 6,996,355.5 | 5,833,908.5 | 7,845,877.0 | 7,000,000 |
| 4N | 56 | 5,636,641.0 | 4,744,963.0 | 6,731,379.8 | 5,600,000 |

The corresponding full-confluence density estimate is

\[
d_{100\%} = \frac{N_{80\%}}{0.8 \times 75\ \mathrm{cm^2}}.
\]

This gives approximately 116,667 cells/cm2 for 2N and 93,333 cells/cm2 for
4N. The lower cell count for 4N is compatible with a larger average adherent
cell footprint, but the thresholds remain empirical approximations and should
be revisited if independent image-based confluence calibration becomes
available.

## Model contract

The T75 estimates remain cohort-level diagnostics shared across C, O1, and O2:

```text
2N protocol_threshold_cells = 7.0e6
4N protocol_threshold_cells = 5.6e6
```

They do not control passage execution. For a nonterminal culture cycle, the
effective feasibility threshold is

```text
max(observed_final_cells, next_passage_initial_cells)
```

Within the interval from the last experimental observation day through the
configured passage-time tolerance, the selected state is the state whose viable
cell count is closest to this threshold from above. A state below the threshold
is never eligible. If the threshold is not reached, the candidate is
protocol-infeasible and the downstream lineage is not reseeded. Reseeding may
downsample to the observed next inoculum but may never upsample.

The individual observed `final_cells` remains among the measured counts used to
estimate the passage-average growth rate. It is not scored as a separate
absolute-count likelihood unit. Independently, it is also the passage-selection
target, so passage feasibility remains a hard condition distinct from the
growth likelihood.

The machine-readable contract is
`oxygen/data/invitro_T75_80pct_protocol_thresholds.tsv`.
