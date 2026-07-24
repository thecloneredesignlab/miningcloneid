# LTEE Hypoxia Manuscript Figure-Set Ideas

## Gate A Closed — Approved July 24, 2026

The user approved the revised five-main-figure architecture for drafting:

1. matched experimental design plus the chromosome-number observations that
   motivate the study;
2. the mechanistic model;
3. separate in-vitro fit adequacy and inferred reshaping mechanism;
4. separate in-vivo fit adequacy, latent resource trajectories, fixed-O2
   associations, and solution landscape; and
5. joint-fit adequacy followed by cross-context parameter and function
   comparisons.

Historical Figure 6 is parked. The rejected Figure 3E ablation is removed from
the active proposal. Compact direct observed–predicted blocks in Figures 3–5
are approved, with objective distributions and optimizer diagnostics
supplementary. No optimizer refit, production manuscript edit, necrosis-export
resolution, Figure 6 work, or new analytical grid is authorized.

## Recorded Decisions

- **Figure 1 approved:** retain the matched design and add the chromosome
  trajectories/distribution changes used to motivate the Results claim.
- **Figure 2 retained:** use one integrated model schematic and distinguish
  observations, inputs, latent variables, state transitions, and predictions.
- **Figure 3A revised:** show growth rate over passage only, not the historical
  multi-row composite.
- **Figure 3 controls corrected:** the left sides of historical panels A and B
  are the control comparators. Preserve those control-versus-deprived contrasts.
- **Figure 3E rejected:** do not propose a standalone negative-control or
  restricted-model panel and do not initiate the associated refit.
- **Figures 3–5 fit quality required:** readers should see observed–predicted
  agreement before fitted mechanisms are interpreted.
- **Figure 4 prior simplification rejected:** keep the low/intermediate/high O2
  feature result and solution-landscape evidence in the main figure.
- **Figure 5 source basis approved:** use the six July joint-pair winners as one
  declared sensitivity universe.
- **Figure 5 reordered:** fit quality first, the full same-universe parameter
  ratios second, then proliferation, missegregation, and survival functions.
- **Figure 5 compact-ratio replacement rejected:** retain the full ratio view in
  the main figure rather than replacing it with a compact generic robustness
  panel.
- **Figure 6 and journal constraints deferred:** revisit Figure 6 after Figures
  1–5 are drafted; journal dimensions and limits are required before polishing,
  not before Gate A drafting.
- **Figure 4 solution regions confirmed:** retain `vi_C01`, `vi_C02`, and
  `vi_C03`, all of which are already displayed in the pooled embedding. No
  two-regime remapping is needed.
- **Fit-quality density approved:** use compact direct observed–predicted blocks
  in Figures 3–5 and place objective distributions and optimizer diagnostics in
  supplementary material.
- **Five-figure architecture approved:** Gate A is closed and figure drafting is
  authorized subject to the explicit scope exclusions above.

## Evidence Corrections That Affect The Plan

### Fit quality is not statistical confidence

The available artifacts support direct in-sample observed–predicted comparisons,
objective ranks across optimizer starts, and sensitivity across multiple fitted
solutions. They do not include held-out validation, profile likelihoods,
bootstraps, formal confidence intervals, or a saved Hessian/Jacobian suitable for
parameter uncertainty.

The figures should therefore communicate two things:

1. **fit adequacy:** how closely the selected fits reproduce the observations
   used for fitting; and
2. **optimization/solution robustness:** whether the interpretation persists
   across competitive starts or across the six approved joint-fit winners.

Bands across fitted solutions are sensitivity envelopes, not confidence
intervals. “Best weighted-MAP fit among 500 starts” is supportable where
applicable; “confident parameter estimate” is not.

### Historical Figure 4's pooled embedding is intentional

Historical Figure 4 panels A–E are in-vivo analyses, while historical panel F
deliberately places the in-vivo circles and in-vitro triangles in the same
pooled embedding. That shared reference landscape is part of the panel's
scientific role and must remain unchanged. Do not filter either context,
recompute a context-specific embedding, or alter the saved geometry. If the
panel later requires technical regeneration for manuscript packaging, it should
be reproduced from the same saved pooled coordinates with the same point
universe and overlays.

### Exact Figure 4 feature interpretation

At O2 = 0, the strongest reproduced one-feature separator is `mu_hp`, the
hypoxia-associated high-ploidy death-strength term (AUC approximately 0.849). At
O2 = 5, the strongest is `p_mis_base`, baseline per-chromosome missegregation
probability (AUC approximately 0.903). The explicitly buffering-related
`buffer_beta` ranks lower at high O2 (AUC approximately 0.672).

These AUCs discriminate model-defined lower- versus higher-ploidy fixed-O2
attractor modes. They do not estimate causal influence on endpoint ploidy.

### All three in-vivo solution regions are retained

The formal in-vivo landscape contains three primary solution regions:
`vi_C01` (99/500), `vi_C02` (385/500), and `vi_C03` (16/500). All three are
already displayed in the pooled embedding and will remain there. No additional
two-regime mapping or cluster-count decision is required.

## Working Five-Figure Architecture

### Figure 1 — Matched systems and motivating chromosome-number changes

- **A:** simplified matched-lineage and sampling timeline.
- **B:** observed in-vitro chromosome-number changes across the informative
  passages, retaining 2N/4N and control/deprived context.
- **C:** observed in-vivo starting-state versus terminal chromosome-number
  distributions with compact burden/harvest context.

The in-vivo observations are start/end comparisons, not longitudinal tumor
karyotyping. The figure motivates environment-dependent modeling without
assigning a fitted mechanism.

### Figure 2 — Mechanistic model overview

Use one integrated schematic with five linked callouts:

- context-specific resource input and latent effective resource state;
- chromosome-number-dependent proliferation and death;
- stress-linked missegregation and WGD generation;
- ploidy-dependent survival following missegregation; and
- the resulting chromosome-state distribution.

Every element should be marked as an observed input, modeled state, fitted
function, state transition, or predicted output.

### Figure 3 — In-vitro fit adequacy and inferred reshaping mechanism

- **A:** observed and predicted growth rate over passage only, faceted by 2N/4N
  and control/deprived. This is the first fit-adequacy view.
- **B:** predicted chromosome-state distributions across passage with observed
  karyotype markers overlaid; retain control on the left and deprived on the
  right. This supplies both the second fit check and the WGD/high-ploidy then
  lower-ploidy reshaping result.
- **C:** fitted post-missegregation survival function.
- **D:** exact-fit nonviable-daughter fraction versus missegregation rate across
  reference ploidies, retained as an integrated model consequence.

There is no Figure 3E. Panels C–D show the mechanism used by the selected fit;
without an ablation they do not prove that the mechanism is necessary or exclude
all strong-elimination alternatives.

The horizontal coordinate in the current fit table is lineage passage, not
literal elapsed time; a calendar-time axis would require an explicit join to the
experiment dates. The normoxic controls were included in fitting and therefore
are fitted baseline comparators, not held-out validation. They also end earlier
than the deprived branches, so late deprived behavior has no contemporaneous
control.

Recommended supplementary diagnostics are the 500-start objective distribution,
objective-component breakdown, and the optimizer-population identifiability
proxy. Seed 10 ranks first, but 12 starts lie within 0.01 objective units; the
local optimizer ended with convergence code 1, five active parameters are at
bounds, and the DEoptim convergence flag is 0/500. These diagnostics must not be
presented as formal uncertainty analysis.

### Figure 4 — In-vivo fit adequacy, resource regimes, and solution landscape

- **A:** compact observed–predicted fit block: tumor burden plus terminal
  chromosome-number distributions for selected seed 25.
- **B:** target versus model-inferred effective-O2 trajectories, explicitly
  labeled as latent.
- **C:** compact low/intermediate/high fixed-O2 feature triptych at O2 = 0, 1,
  and 5, retaining the scientific result in the main figure.
- **D:** the existing pooled in-vivo/in-vitro parameter embedding, retained
  unchanged with both contexts, the saved geometry, and the existing solution
  overlays.
- **E:** focused `p_mis_base` and `n_O` distributions across the displayed
  solution regions.

The selected seed is the best total weighted-MAP objective among 500 starts, but
it ranks 52nd for the ploidy component and 189th for the burden component.
Twenty-nine fits are within 1% and 241 within 5% of its total objective. The
solution multiplicity is therefore an uncertainty result, not clutter to hide.

Predicted necrosis values are `NA` in the exported `necrosis_fit.tsv`, so
necrosis must not appear in the fit-quality block until that export inconsistency
is resolved. Full parameter screens and objective diagnostics can be
supplementary, but the O2 triptych and an interpretable landscape remain main.

### Figure 5 — Joint-fit adequacy and context-specific functions

- **A:** in-sample observed–predicted summaries for both contexts, using the six
  approved July pair winners. Show in-vitro growth/karyotype and in-vivo
  burden/terminal-ploidy evidence; do not use the incomplete necrosis export.
- **B:** the full all-six in-vivo/in-vitro parameter-ratio view on a log-ratio
  scale, with all six values plus a median/range summary for every coupled
  parameter.
- **C:** in-vivo versus in-vitro proliferation functions.
- **D:** in-vivo versus in-vitro stress-linked missegregation functions.
- **E:** in-vivo versus in-vitro post-missegregation survival functions.

Panels B–E must be regenerated from the same six-winner universe. Across-solution
traces or envelopes show sensitivity, not replicate uncertainty. The existing
generators operate one winner at a time, so drafting will require an
all-six aggregation layer but no optimizer refit.

The six winners span total objectives 18.852–19.978. Their stored validation
checks confirm selection provenance and objective arithmetic, not predictive or
held-out validation. All pair runs report 500/500 DEoptim convergence flags, but
no local L-BFGS-B refinement was accepted.

## Approved Fit-Quality Encoding

The approved main-text encoding is direct and compact:

- observations as points, empirical distributions, or thin reference marks;
- fitted values as lines, densities, or paired predicted distributions;
- one small seed/rank annotation rather than an objective bar dominating the
  panel; and
- competitive-solution variation in a supplementary diagnostic or as a light
  sensitivity envelope where it materially changes interpretation.

An objective-component bar alone is not sufficient because it does not show
where or how the fit reproduces the data. Conversely, packing every diagnostic
into the first panel would recreate the density problem in historical Figure
3A.

## Gate A Resolution

Gate A is closed. The approved drafting basis is:

1. the five-figure architecture described above;
2. all three existing in-vivo solution regions retained in the unchanged pooled
   embedding; and
3. compact direct observed–predicted blocks in Figures 3–5, with objective
   distributions and optimizer diagnostics supplementary.

## Drafting Scope Boundaries

- final figure assembly and manuscript caption edits;
- production-code changes beyond bounded figure generators;
- optimizer refits or restricted-model comparisons;
- the Figure 6 fixed-O2 analytical grid and model-selection framing;
- formal resolution of the necrosis prediction export;
- journal-specific dimensions, color rules, and main-figure limits; and
- statistical uncertainty claims not supported by the current artifacts.
