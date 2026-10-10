# In-vivo neighborhood eFAST

Each of 500 optimizer endpoints defines an independent rectangular neighborhood: best value +/- 10% of the full natural parameter span, clipped to the original bounds. The inherited log10 or identity sampling transform is retained.

Indices are computed separately for every neighborhood and phase. Five phase indices are averaged within each fit seed; the heatmaps show medians across fit seeds. Source tables also report means, IQR and valid-seed counts. Raw outputs from different neighborhoods are never pooled into one FAST analysis.

The bootstrap resamples complete fitted-seed sensitivity curves with their five-phase means retained. Optimizer endpoints are not independent posterior draws; bootstrap classifications describe repeatability across these endpoints and are not posterior significance claims.

Mechanism comparisons use group means of parameter indices. Total effects overlap through interactions; these means are not additive group variance fractions. The original Figure 4B correlation annotation supplies directional context.

All 500 endpoints are retained regardless of convergence diagnostics. summaries/seed_convergence_status.tsv and seed_convergence_diagnostics.tsv report repeat and resolution checks separately. Resolution stability is tested only for the multi-N pilot endpoints; a stable five-phase range alone does not establish resolution convergence.

## Fraction of evaluable fit seeds supporting the proposed shift

- dominant_mean_ploidy, S1: 0.076 (38/500); requires low-O2 death dominance and higher-O2 increases in both missegregation and buffering.
- dominant_growth_rate, S1: 0.052 (26/500); requires low-O2 death dominance and higher-O2 increases in both missegregation and buffering.
- dominant_mean_ploidy, ST: 0.066 (33/500); requires low-O2 death dominance and higher-O2 increases in both missegregation and buffering.
- dominant_growth_rate, ST: 0.056 (28/500); requires low-O2 death dominance and higher-O2 increases in both missegregation and buffering.

Interpret these fractions with the pilot resolution and phase diagnostics. Constant-output trajectories have undefined indices and remain recorded as NaN; reported denominators identify the evaluable subset.
