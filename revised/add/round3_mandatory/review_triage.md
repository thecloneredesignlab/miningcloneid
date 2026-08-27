# Round-3 review triage and scientific decision record

## Decision standard

An item is **must change** when it affects causal interpretation, identifiability, predictive adequacy, numerical validity, or auditability. An item is **reasonable but bounded** when the concern is scientifically sound but the requested implementation would require new data or an out-of-scope model redesign. An item is **editorial** when it improves clarity without changing the evidence.

## Adjudication

| ID | Review issue | Decision | Required response |
|---|---|---|---|
| R1 | Causal language exceeds the observational/model-conditioned comparison | Must change | Use model-based, association-level language in title, abstract, Results, and Discussion; distinguish generated latent trajectories from directly measured endpoints. |
| R2 | Optimizer seeds and solution clusters are treated as uncertainty | Must change | Call them numerical multimodality/robustness analyses, not formal statistical uncertainty; retain no unsupported confidence interpretation. |
| R3 | A single in-vitro anchor and saturated soft-coupling form may determine the joint contrast | Must change, partly refit-dependent | Audit same-region and other-region near-optimal anchors from existing fits; define a refit matrix for alternative anchors and coupling/regularization forms. Do not claim completed robustness until refits exist. |
| R4 | Predictive adequacy is shown mainly for representative fits | Must change | Report all-sample in-vitro likelihood/residual diagnostics, all-tumor in-vivo burden and terminal-ploidy diagnostics, and objective-near envelopes. Label these as in-sample checks; held-out refits remain a separate task. |
| R5 | Fixed-O2 eigenvector equilibria lack finite-time validation | Must change | Audit the existing matrix-exponential finite-time trajectories, spectral gaps, convergence flags, and class labels; restrict equilibrium claims where validation fails. |
| R6 | Boundary loss, CIN-associated nonviability, WGD implementation, and necrosis are not fully auditable | Must change | Quantify boundary versus within-grid nonviable loss on the evaluated grid; state the WGD kernel implementation; expose the necrosis export defect; reconstruct predictions from stored states and verify the objective while distinguishing this audit repair from an upstream exporter fix. |
| R7 | O1/O2 terminal mismatch needs quantitative rather than qualitative treatment | Must change | Report lineage-specific log scores, total-variation distance, Wasserstein-1 distance, Jensen-Shannon distance, and the high-state predictive check against the shared fitted distribution. |
| R8 | Longer figure legends and explicit data/model distinctions | Reasonable/editorial | Adopt continuation legends and mark directly observed, fitted, and simulated quantities consistently. |
| R9 | Full external validation or new wet-lab experiments | Reasonable but not achievable from current data | State as a limitation and future validation requirement; do not fabricate evidence or substitute an in-sample check. |
| R10 | An IF threshold can be guaranteed by revision alone | Not scientifically decidable | Target rigor, transparency, and falsifiability; do not claim that any analysis guarantees a journal impact factor. |

## Package reconciliation

The two review packages agree on the central scientific risks but differ in strictness. The stricter package is used as the controlling standard because it directly addresses evidence strength and model auditability. The deep-research report is used as a secondary synthesis: recommendations are accepted only when traceable to the manuscript, frozen inputs, or executable model outputs.

## Non-negotiable manuscript stance

1. The cross-context comparison is model-conditioned and association-level.
2. Seed ensembles characterize numerical solution multiplicity, not sampling uncertainty.
3. Representative curves are illustrations; conclusions must be tied to all-sample or objective-near summaries.
4. Fixed-O2 equilibrium summaries are conditional on finite-time convergence and spectral-gap diagnostics.
5. Missing necrosis predictions are disclosed as an export/audit defect, not silently interpreted.
6. No claim is added unless its numerical support is generated into this directory.
