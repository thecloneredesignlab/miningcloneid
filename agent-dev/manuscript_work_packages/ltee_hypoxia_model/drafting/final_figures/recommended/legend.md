# Recommended figure legends

## Figure 1

**Figure 1. Matched experimental systems and motivating chromosome-number
changes.** **(A)** Simplified design for matched SUM159 near-diploid (2N) and
near-tetraploid (4N) lineages. Culture strips show the longest tracked control
and oxygen-deprived paths, with known target O2 encoded from light to dark;
triangles mark flow cytometry and circles mark chromosome counts. Tumor strips
show burden-measurement times and terminal chromosome sampling. Tumor resource
state was not measured and is latent in the model. **(B)** Observed cell-level
chromosome numbers in culture. Lines and ribbons show medians and interquartile
ranges; the shared starting reference begins both summary trajectories. Control
series end earlier than oxygen-deprived series. **(C)** Starting cell-line
reference karyotypes and terminal tumor-cell chromosome distributions. Boxes
show distribution summaries, violins show terminal cell distributions, and
diamonds show individual tumor medians (four tumors per starting-ploidy
cohort). Starting and terminal measurements are distinct samples rather than
longitudinal karyotyping of the same tumor.

## Figure 2

**Figure 2. Model overview: resource state links growth, chromosome-state
generation, and survival filtering.** **(A)** In culture, the prescribed O2
schedule is an external input; the starting chromosome-number distribution is
an initial state. In tumors, effective O2 is a latent supply-demand state whose
demand term is computed from the simulated live population. Measured tumor
burden is a fit target, not the resource-state driver. **(B)** Fitted
proliferation and death functions depend on chromosome number and effective O2.
**(C)** Stress-linked missegregation generates a modeled distribution of
multi-chromosome shifts, `N -> N +/- delta N`, whereas WGD generates doubled
states, `N -> 2N`. **(D)** A fitted ploidy-dependent survival filter determines
which missegregated daughters persist. **(E)** These processes generate the
predicted chromosome-state distribution and predicted growth, live-burden, and
ploidy summaries, which are compared with in-vitro growth/karyotype/flow and
in-vivo burden/terminal-karyotype observations. The displayed output histogram
is schematic.

## Figure 3

**Figure 3. Direct fit quality and fitted chromosome-segregation functions.**
**A,** Observed growth rates (points) and fitted growth rates from selected seed
10 (lineage-connected lines) across lineage passage, shown separately for
control and oxygen-deprived 2N and 4N branches. **B,** Fitted chromosome-state
fractions (heatmaps) with individual observed karyotypes (open points) overlaid
at sampled lineage passages; all saved chromosome states and control branches
are retained. **C,** Fitted post-missegregation survival as a function of
chromosome number. Dashed lines and labeled points mark the 2N (`N = 44`) and
4N (`N = 88`) reference states. **D,** Model-implied nonviable-daughter fraction
versus per-chromosome missegregation probability for 1.5N–5N reference states,
using the saved oxygen-sweep curves. Controls were included in fitting and are
not held-out validation data. Lineage passage denotes branch depth rather than
calendar time. Panels C and D are fitted/model-implied quantities, and no
parameter-uncertainty intervals are available.

## Figure 4

**Figure 4. In-vivo fit adequacy, latent oxygen dynamics, and the pooled
fitted-solution landscape.** **A,** Observed versus fitted tumor burden for
selected seed 25 (points; dashed line, identity) and weighted observed versus
fitted terminal chromosome-number distributions across the eight in-vivo
endpoints. Seed 25 is the lowest-total-objective weighted-MAP fit among 500
starts; these are in-sample comparisons. Necrosis is omitted because its
exported predictions are unavailable. **B,** Prescribed oxygen target (orange
dashed line, drawn above the blue trace) and latent model-implied effective
oxygen (blue line) for the eight exposure schedules. Effective oxygen is a
fitted state, not a tumor oxygen measurement. **C,** Strongest univariate
separator of the two model-defined fixed-oxygen attractor modes at 0%, 1%, and
5% O2. At 0% O2, `mu_hp` is strongest (AUC = 0.849; larger in the lower-ploidy
mode); at 1% and 5% O2, `p_mis_base` is strongest (AUC = 0.741 and 0.903,
respectively; larger in the higher-ploidy mode). AUC quantifies discrimination,
not causal influence. **D,** Preserved pooled in-vivo/in-vitro t-SNE point
universe with the saved coordinates, initial-sample context, objective-value
overlays, geometry, and labeled in-vivo (`vi_C01`–`vi_C03`) and in-vitro
(`vt_C01`–`vt_C02`) solution regions. No embedding or model fit was recomputed.
**E,** Distributions of fitted `log10(p_mis_base)` and `n_O` across the three
formal in-vivo regions (`vi_C01`, n = 99; `vi_C02`, n = 385; `vi_C03`, n = 16);
white boxes show medians and interquartile ranges. Regions denote fitted
solution families rather than biological tumor subtypes, and proximity in the
embedding is descriptive.

## Figure 5

**Figure 5. Joint-fit adequacy and context-specific fitted functions.**
**A,** Direct in-sample observed-versus-predicted comparisons for within-series
normalized tumor burden, terminal mean chromosome number, passage growth rate,
and mean karyotype. Pale points show predictions from all six approved July
joint-pair winners; filled points show their median, and vertical bars show
their range where legible. Dashed diagonals denote equality. **B,** In-vivo to
in-vitro ratios for all 14 soft-coupled parameters on a common log2 scale.
Colored points are the six winners (color denotes warm-start region), gray
segments span their range, black diamonds mark their median, and the vertical
dashed line marks a ratio of one. **C,** Fitted proliferation rate versus
oxygen for 2N and 4N reference states. **D,** Fitted per-chromosome
missegregation rate versus oxygen for the same states. **E,** Fitted survival
after missegregation versus chromosome number. In C–E, orange solid curves are
in vivo, blue dashed curves are in vitro, thin curves are the six selected
solutions, and thick curves are pointwise medians. C–D use linear display
interpolation on a common 0–5% oxygen grid. All fit comparisons are in-sample,
and across-winner ranges and traces represent solution sensitivity rather than
confidence intervals. Necrosis is omitted because predicted values are
unavailable in the saved joint-fit exports.

## Supplementary fit-quality and optimizer diagnostics

**Supplementary figure. Optimizer-start distributions and solution
multiplicity.** **(A)** Objective difference from the best result across 500
starts for the separate in-vitro and in-vivo fits; highlighted points are seed
10 and seed 25, respectively. **(B)** Objective difference within each of the
six approved 500-start joint-pair runs; open points mark the selected winner of
each run. All six are rank-1/zero-difference winners and are shown in separate
warm-start-pair facets. In A and B, the y-axis is log10; exact zero differences
are displayed at `1e-4`. **(C)** Stored run diagnostics,
competitive-solution counts, and active-parameter boundary counts. These
quantities describe optimization and solution sensitivity. They do not
constitute held-out validation, posterior uncertainty, confidence intervals,
or identifiability analysis.
