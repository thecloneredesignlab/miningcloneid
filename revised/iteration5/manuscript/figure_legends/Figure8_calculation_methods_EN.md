# Figure8 and supplementary figures: parameters, finite-time simulation, passage, and statistical summaries

Execution update (2026-10-17): the complete current pipeline is `run_Figure8.R` (prepare → simulate → aggregate → validate → plot), with one unversioned `data/Figures/Figure8` directory. See [complete pipeline guide](../../Code/Figures/Figure78_pipeline.md). The scientific definitions below remain unchanged. Versioned run names and old helper names elsewhere in this document identify historical published receipts only; they are not runtime dependencies of the new pipeline. Existing published numerical examples remain historical reference values until a full regeneration is separately authorized.

Documentation and implementation-review date: 2026-10-15 (America/New_York).

Chinese counterpart: [Figure8 calculation methods](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/manuscript/figure_legends/Figure8_calculation_methods.md). Both versions use corresponding sections, formulas, numerical examples, and recalculation code. This document records the current implementation; it does not introduce a new model or modify existing figures.

The key statistical hierarchy is: population mean ploidy in one simulation → average over culture-passage repeats within one fitted seed → median and IQR across fitted seeds. Displaying an ensemble median does not change the model's internal population-mean ploidy into a median over individual cells.

## 1. Current figure scope and version selection

### 1.1 Main and supplementary figure map

“Native fitted p” means each seed's own final `p_misseg`; “fixed p” overrides only that native parameter, not the effective missegregation probability. Model O₂ is imposed in percent, with 0.1% steps for oxygen-response curves.

| Figure | Content and environment | p mode | Initial populations | Time | Seed aggregation |
| --- | --- | --- | --- | --- | --- |
| Figure8A | PCAWG/DepMap organ ploidy and three model anchors | Native fitted p for models | Model 2N start; fitted culture distribution | Model day 1000 | Models pool C01 500+C02 500; empirical groups remain organ/cohort-specific |
| Figure8B | Author Figure5: WGD fraction versus organ oxygen pressure | Not applicable | Not applicable | Not applicable | 19 original author-source rows |
| Figure8C | in vivo, O₂=0–5% | Native fitted p | 2N, 4N, 1:1 mixture | Day 1000 | Pool 1000 |
| Figure8D | in vitro, O₂=0–5% | Native fitted p | 2N, 4N, 1:1 mixture | Day 1000 | Pool 1000 after within-seed passage-repeat averaging |
| supp_fig8_1 | Fixed-p counterpart of A | 0.01 | As in A | Day 1000 | Models pool 1000 |
| supp_fig8_2 | Fixed-p counterpart of A | 0.3 | As in A | Day 1000 | Models pool 1000 |
| supp_fig8_3 | in vivo, O₂=0–5%, C01/C02 separate | Native fitted p | Seven Q values | Day 1000 | 500 per family |
| supp_fig8_4 | in vivo, O₂=0–5%, two fixed-p panels | 0.01, 0.3 | Seven Q values | Day 1000 | Pool 1000 in each panel |
| supp_fig8_5 | in vitro, O₂=0–5%, fixed-p rows and C01/C02 columns | 0.01, 0.3 | Seven Q values | Day 1000 | 500 per family after repeat averaging |
| supp_fig8_6 | in vitro, O₂=0–20%, 2×5 time matrix | Native fitted p | 2N, 4N, 1:1 | Days 100–1000, every 100 days | Pool 1000 |
| supp_fig8_7 | in vivo, O₂=0–5%, 2×5 time matrix | Native fitted p | 2N, 4N, 1:1 | Same times | Pool 1000 |
| supp_fig8_8 | Family split of supp_fig8_6, C01 above C02 | Native fitted p | As in supp_fig8_6 | Same times | 500 in each family block |
| supp_fig8_9 | Family split of supp_fig8_7, C01 above C02 | Native fitted p | As in supp_fig8_7 | Same times | 500 in each family block |

### 1.2 Current complete pipeline

`run_Figure8.R` reads every original full seed/context input and the author raw empirical files, independently generates finite-time predictions, and writes all pooled/family plotting tables into one current Figure8 directory. The fitted culture configuration and stochastic passage are applied in the same current worker; snapshots and terminal results share the same trajectory. `draw_Figure8.R` renders these newly generated tables without old runs, fallbacks, or archived sources.

The retained published data remain supported by a separate read-only compatibility path. Earlier violin plots, C/D representing separate families, fixed-p main C/D, and a single-time supp_fig8_6 remain superseded. Historical version names below document previous receipts, not required computational stages.

## 2. Empirical data and source material

### 2.1 Local author-supplied bundle

The read-only input bundle is:

```text
/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/data/share_bundle_figures_3_5
```

| File relative to the bundle root | Role |
| --- | --- |
| `data/Figure3_PCAWG_included_samples.csv` | The author's mapped/included 2338 tumor records for Figure3 |
| `data/Figure3_DepMap_included_lines.csv` | The author's mapped/included 1433 cell-line records; 1185 have finite ploidy |
| `originals/organ_lut.csv` | Saved 19-row Figure5 input controlling current B points and SDs |
| `supplementary/Supplementary_Table_3.csv` | Cross-check of Figure5 labels, WGD fractions, oxygen means, and SDs |
| `scripts/oxygen_plots.Rmd` | Author's plotting and ordinary-least-squares analysis source |
| `originals/nowicka.txt` | Earlier 20-row oxygen/WGD data, used only for source-version comparison |
| `documentation/METHODS.md` | Bundle provenance, mappings, missingness, and interpretation limitations |

Saved means, SDs, or quartiles are not resampled as if they were individual oxygen measurements. No nonexistent observations are reconstructed from a figure image.

### 2.2 Article and scope of the reproduction

The empirical comparison relates to Beck et al., [Energetics of whole genome doubling and genomic instability, Cancer Letters, 2025, DOI:10.1016/j.canlet.2025.217878](https://pubmed.ncbi.nlm.nih.gov/40517938/). Its Figure3 compares ploidy medians and quartiles; Figure5 relates tumor WGD fractions to organ oxygen pressure, with oxygen mean±SD bars.

Current A reaggregates the author's mapped records with finite ploidy. Current B preserves the 19-row author source rather than constructing a cleaned independent-organ cohort. Other preliminary mean-ploidy–oxygen analyses and primary-only/donor-level summaries retained in the bundle or outputs are not the current B input.

## 3. Fitted parameters, configurations, and model provenance

### 3.1 One thousand complete parameter sets

Figure8 uses 500 joint-fitting seeds per family, C01 and C02, for 1000 total. Every seed retains its own complete final natural-scale parameters. This is neither a one-best-seed-per-family analysis nor a prediction using only 14 soft-coupling parameters.

The fit root is:

```text
/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/results/fit_joint_unified_global_invitro_500seed_all_xxlarge_r442_exact_20260828_145253
```

| Family | Family directory | Seeds |
| --- | --- | --- |
| C01 | `fit_joint_tsne_vi_seed366_C01_vt_seed144` | 1–500 |
| C02 | `fit_joint_tsne_vi_seed25_C02_vt_seed144` | 1–500 |

In each seed's `joint_best_params_long.tsv`, `shared_invivo` supplies 26 in vivo parameters and `invitro_effective` supplies 33 in vitro parameters. `complete_input_signature`, input manifests, and parameter-level passthrough tables retain identities and values.

No fitted-endpoint deduplication is performed: identical parameters or predictions from two optimization seeds still occupy two ensemble positions. Parameters are not shuffled between seeds. No additional selection is imposed according to convergence, survival direction, or agreement with the expected trend. Parameters inactive in a fixed-O₂ experiment retain their original values rather than being replaced with common values.

### 3.2 Authoritative in vitro configuration

The repaired culture configuration is read strictly from each seed's `fit_result.rds$ctx$invitro_cfg`, without fallback to in vivo `fit_config.rds`. All 33 parameters are checked against the final parameter table, including fitted Gaussian initialization parameters.

The current 1000 culture configurations record `DT=0.05` days, `Crowding=TRUE`, `crowding="logistic"`, and `K=10^12`. These are verified configuration values; their equality does not mean individual configurations were omitted from prediction.

### 3.3 Native p versus effective probability

Native-p figures retain each seed's fitted `p_misseg`. Fixed-p supplements use `figure6_force_p_misseg()` to change only that natural-scale scalar to 0.01 or 0.3, preserving seed-specific `p_mis_base`, `k_o_mis`, and other parameters.

The effective probability used by the model is:

\[
p_{\mathrm{mis,eff}}(n,O_2)
=\operatorname{clip}_{[0,1]}\left[
p_{\mathrm{mis,base}}
+p_{\mathrm{misseg}}\frac{\mu(n,O_2)}{\mu(n,O_2)+k_{o,\mathrm{mis}}}
\right].
\]

Fixing the native scalar therefore does not fix the effective probability in every chromosome state. Rates remain model-derived across oxygen, environment, family, and seed. Figure8 does not apply the Figure7 β blockade.

### 3.4 Actual model implementation

Frozen inputs point to:

```text
/Users/4482173/Documents/GitHub/soft_couping_org/oxygen/code/O2_supply_demand_MAP
```

The closure includes the fixed-oxygen entry point, `model/model_O2_supply_demand_MAP.R/.cpp`, and associated utilities. `f78_engine()` copies required files byte-for-byte into the run's `.model_snapshot` and verifies SHA256 before loading. A similarly named oxygen directory in the figure repository does not silently replace that model.

`FIGURE_MODEL_CODE_ROOT` and `FIGURE_JOINT_RESULT_ROOT` can select sources during new preparation. Published provenance instead follows frozen manifests, operators, and fingerprints. Documentation does not switch sources or recompile the model.

## 4. Population ploidy and initial conditions

### 4.1 Units of a single population's mean ploidy

The current chromosome grid is 22–154, with 133 states and `N_UNIT=22`. Let \(X_{sri}(t)\) be living cells for seed s, culture repeat r, and state i, with chromosome count \(n_i\):

\[
P_{sr}(t)=\frac{\sum_i(n_i/22)X_{sri}(t)}{\sum_iX_{sri}(t)}.
\]

The unit is N. This is the cell-number-weighted mean ploidy of the living population at that time, not the initial population's mean, not a statistic including dead cells, and not a median over cells. In vivo has no random-passage repeat index r.

The implementation stores normalized states \(x_i=X_i/\sum X\) and tracks total population separately. The same quantity is \(\sum_i(n_i/22)x_i\). This avoids growth-related floating-point overflow without changing the weighting.

### 4.2 Pure-state in vivo starts

Pure in vivo 2N and 4N starts place unit mass at 44 and 88 chromosomes. They are not sampled from fitted culture distributions. Each seed evaluates these starts with its own model at the imposed O₂.

### 4.3 Fitted Gaussian in vitro starts

Culture “2N” and “4N” identify initial cohorts, not a requirement that every cell has exactly 44 or 88 chromosomes. Each seed supplies its own `init_mean_2N`, `init_sd_2N`, or corresponding 4N parameters:

\[
g_{s,k,i}=\frac{\phi(n_i;m_{s,k},\sigma_{s,k})}
{\sum_j\phi(n_j;m_{s,k},\sigma_{s,k})},\qquad k\in\{2N,4N\}.
\]

\(\phi\) is the normal density evaluated on the chromosome-count grid; its parameter units match that grid. The actual helper, `ivt_gaussian_initial_state()`, requires a finite mean and positive SD, evaluates `dnorm` on the existing finite grid, and normalizes the weights. Density is not treated as an absolute cell count, nor integrated over continuous bins.

### 4.4 Mixture ratio Q

\[
Q=\frac{\text{initial 4N-cohort cell number}}{\text{initial 2N-cohort cell number}},\qquad
f_2=\frac{1}{1+Q},\quad f_4=\frac{Q}{1+Q}.
\]

Pure 2N is Q=0; pure 4N is Q=∞; 2N:4N ratios 1:1–1:5 correspond to Q=1–5. Combining two individually normalized starting vectors gives:

\[
x_s(0;Q)=f_2x_{s,2}(0)+f_4x_{s,4}(0).
\]

Main curves use Q=0, 1, ∞; seven-ratio supplements additionally use Q=2, 3, 4, 5. All mixtures and pure starts have the same total inoculum; 1:1 does not double the inoculum. Q is an initial cohort-number ratio, not a terminal lineage ratio, and a 1:1 mixture is not a pure-3N population.

## 5. in vivo: fixed-oxygen finite-time propagation

### 5.1 Continuous-growth operator

The model supplies state-specific division rates λ, direct-death rates μ, viable ordinary daughter weights B, and WGD weights W. Rows index daughter states and columns maternal states:

\[
J=B\operatorname{diag}[\lambda(1-p_w)]
+W\operatorname{diag}(\lambda p_w),\qquad
G=J-\operatorname{diag}(\lambda),\qquad M=G-\operatorname{diag}(\mu).
\]

B/W already contain daughter multiplicity and viability; daughter counts must not be multiplied again. With fixed O₂, parameters, and environment:

\[
\dot X=MX,\qquad X(t)=\exp(tM)X(0).
\]

Main panels use t=1000 days; time matrices use the stated finite t. No starting state is forced to a dominant eigenvector, and convergence is not a prerequisite.

### 5.2 Spectral shifting is numerical stabilization only

For \(a=\max\operatorname{Re}\operatorname{eig}(M)\), compute \(\exp[t(M-aI)]\) and retain the common factor \(\exp(at)\) in log mass. That factor cancels in normalized ploidy but must be retained when combining cell numbers from different initial cohorts.

Using an eigenvalue to stabilize the exponential is not using an equilibrium eigenvector as the prediction. Time snapshots apply consecutive 100-day exponentials and retain state-normalization factors and log mass between intervals.

### 5.3 Mixtures require terminal population weighting

The continuous system is linear in the initial vector. Let \(m_2(t),m_4(t)\) be pure-start population expansion factors and \(P_2(t),P_4(t)\) their population means. Then:

\[
P_Q(t)=\frac{f_2m_2P_2+f_4m_4P_4}{f_2m_2+f_4m_4}
=w_2P_2+(1-w_2)P_4,
\]

\[
w_2=\operatorname{logistic}[\ell_2-\ell_4-\log Q],\qquad
\ell_k=\log m_k,\quad 0<Q<\infty.
\]

The pure limits use w₂=1 or 0. Initial fractions must not directly average terminal ploidies unless terminal expansion factors are equal.

Main C and seven-Q in vivo sources use this exact superposition. Time matrices directly propagate the pure-2N, 1:1, and pure-4N vectors together. Their finite-time equivalence was checked. A mixture close to the 2N curve can reflect terminal numerical dominance of the 2N-origin cohort without any incorrect initialization as pure 2N.

## 6. in vitro: finite-time dynamics and random passage

### 6.1 Culture is not a constant linear matrix plus dilution

Every seed retains its authoritative culture configuration and complete parameters. In the current fixed-O₂ production path, the actual living population is \(\mathcal N(t)\), with crowding:

\[
c(t)=\max[0,1-\mathcal N(t)/K_{\mathrm{cap}}],\qquad K_{\mathrm{cap}}=10^{12}.
\]

The absolute-count Euler update on 0.05-day steps is:

\[
X^+=X+\Delta t\left[c(t)GX-\operatorname{diag}(\mu)X\right],
\qquad \Delta t=0.05\ \mathrm{day}.
\]

Ordinary/WGD production and maternal removal are crowding-scaled; direct death is not. Normalized composition and total count are updated together, and crowding is recalculated after passage as well. The earlier erroneous constant-linear culture propagator is not substituted for this process.

Imposed O₂ and the following random-threshold passage protocol are the approved prediction design. Preserving fitting dynamics and configuration does not mean reproducing all variable-oxygen conditions or every original experimental passage date.

### 6.2 Inoculum and threshold provenance

The saved canonical rule uses existing equal-lineage geometric summaries of six lineages and 114 maintenance passages. Before integerization, the inoculum is 569799.21840352, the expansion ratio approximately 10.84729935146, and the target 6180782.69225091.

Specifically, within lineage ℓ, geometric means of valid passage inocula S_ℓj and expansion ratios T_ℓj/S_ℓj give S_ℓ and R_ℓ. The six lineage means then receive equal geometric weight:

\[
S_* = \exp\left[\frac{1}{6}\sum_{\ell=1}^{6}\log S_\ell\right],\qquad
R_* = \exp\left[\frac{1}{6}\sum_{\ell=1}^{6}\log R_\ell\right],\qquad
N_{\mathrm{seed}}=\operatorname{round}(S_*),\quad
N_{\mathrm{threshold}}=\lceil S_*R_*\rceil.
\]

Thus lineages are equally weighted rather than pooling 114 records and automatically giving more influence to lineages with more passages. The target uses a ceiling and the inoculum rounding. The representative approximately 5 days is the median of lineage-specific passage-duration medians.

Production reseeds 569799 living cells and uses a threshold of 6180783. These are existing, data-derived, approved common protocol values, not parameters fitted independently for each of 1000 seeds or an arbitrary single observed cell count.

### 6.3 Passage-event sequence

1. Continue from the current post-inoculation state, with 20 Euler substeps per actual day.
2. After growth on each actual integer day, check total living cells. Passage on the first day meeting threshold; otherwise keep waiting.
3. Use that day's actual population, not an assumed exact-threshold population. Systematically randomized integerization converts expected chromosome-state abundances to integer counts, with the total rounded to the nearest integer.
4. Sequential hypergeometric draws implement multivariate hypergeometric sampling without replacement, retaining exactly 569799 cells.
5. Update composition, reset total count to the inoculum, and continue along the same elapsed-time axis.

Approximately 5 days in the rule is a representative experimental duration, not mandatory passage every five days. Time is neither reset nor rewound. Without threshold crossing, no artificial passage is inserted. Discarded passage cells are not direct biological deaths.

### 6.4 Repeats, RNG, and the deterministic exception

Each seed×O₂×p mode×initial population has 200 nominal passage repeats. Family, fitted seed, O₂, p, and Q/initial condition identify reproducible RNG streams; repeats have separate substreams. An optimizer seed identifier and a passage-repeat index are different entities.

The implementation first computes a conservative rate bound:

\[
r_{\mathrm{ub}}=\max_j\left\{\max\left[\sum_iG_{ij},0\right]-\mu_j\right\}.
\]

If \(1000r_{\mathrm{ub}}<\log(6180783/569799)\), threshold crossing is impossible for c∈[0,1]. All nominal repeats would coincide, so one deterministic Euler trajectory is evaluated, recording `executed_repeats=1` and `n_stochastic_repeats=200` separately. A no-passage culture condition is not changed into an in vivo calculation.

### 6.5 Culture mixtures are simulated directly

Each mixed culture starts from its seed-specific mixed fitted Gaussian vector and directly undergoes density feedback and threshold sampling. Composition affects growth, threshold times, and samples; mixtures cannot be reconstructed by averaging pure endpoint ploidies or applying in vivo linear superposition.

The v11 repair regenerated culture predictions and actual passage events, without importing old v4/v9/v10 culture endpoints or passages. v12 uses the same authenticated engine to record snapshots from daily values along uninterrupted trajectories; it does not restart at each snapshot.

## 7. Medians and IQRs: three statistical levels

### 7.1 Average passage repeats within seed first

Within a condition, each culture repeat yields population mean \(P_{sr}(t)\). The seed-level endpoint or snapshot is its arithmetic repeat mean:

\[
Y_s(t)=\frac{1}{R_s}\sum_{r=1}^{R_s}P_{sr}(t).
\]

Repeats have equal weight, without reweighting by terminal repeat population. Cell-number weighting already occurs inside each repeat's population mean. The deterministic no-passage exception uses R_s=1. In vivo Y_s is that seed's deterministic population mean.

### 7.2 Then calculate medians and quartiles across seeds

For seed set S in the specified condition:

\[
\mathrm{center}=Q_{0.5}\{Y_s:s\in S\},\qquad
\mathrm{lower}=Q_{0.25}\{Y_s:s\in S\},\qquad
\mathrm{upper}=Q_{0.75}\{Y_s:s\in S\}.
\]

Pooled S contains 1000 `(pair_label, seed_number)` identities; a within-family S contains 500. The calculation does not pool all 200000 culture repeats before taking quantiles or treat cells, passage repeats, and optimizer seeds as the same sampling unit.

R uses type-7 quantiles: for sorted y, \(h=1+(n-1)p\), with linear interpolation between \(y_{\lfloor h\rfloor}\) and the next value. Bars/ribbons span Q25 to Q75; IQR width is Q75−Q25.

### 7.3 Correct order for pooling families

Concatenate 500 C01 Y_s values and 500 C02 Y_s values before computing quantiles. Each seed has equal weight and each family occupies 500 positions. Family medians or quartiles are not averaged; objective value, terminal population, and repeated optimization endpoints do not determine additional weights.

### 7.4 What these intervals are not

Model IQR is descriptive dispersion of predictions from the optimizer-fit ensemble, not a 90%/95% confidence interval, posterior credible interval, or biological-replicate variability. Monte Carlo SE of the repeat mean is stored separately as `stochastic_mcse` and is not the plotted seed IQR.

Identical or nearly identical predictions can produce thin or zero-width IQRs. Their actual widths are shown without jitter or artificial expansion. Repeated optimization endpoints are not removed.

## 8. Figure8A: empirical and model ploidy medians/IQRs

### 8.1 Empirical sample aggregation

Within the author's `figure3_label` groups, PCAWG and DepMap are summarized separately over finite `ploidy` values. The cohorts share 18 mapped organ categories, giving 36 empirical groups and three model groups, or 39 displayed rows.

PCAWG input contains 2338 unique `sample_id` records. DepMap contains 1433 unique `model_id` records, of which 248 lack finite ploidy, leaving 1185 contributing records. Missing values are not zero, and included records without ploidy are not counted as valid ploidy n.

This reproduction adds no primary-only, purity, WGD-uncertainty, or new QC filter and does not first average by donor. The author's primary-tumor caption wording does not fully match the included local specimen types; that limitation is retained explicitly. DepMap 24Q2 is source attribution, not an independently fingerprinted archive release.

### 8.2 Three model anchors

Each model anchor uses day 1000, a 2N-cohort start, and 1000 pooled seeds: in vivo O₂=0.5%, in vitro O₂=0.5%, and in vitro O₂=20%. Culture starts from the fitted 2N Gaussian and includes actual threshold passage.

Both empirical and model rows show N-unit medians and Q25–Q75, but their distribution units differ: tumor/cell-line records versus seed predictions after repeat averaging. Matching units and graphics do not imply identical uncertainty or paired observations.

## 9. Figure8B: reproduction of the author's Figure5

### 9.1 Points, coordinates, and error bars

The horizontal coordinate is the supplied tumor WGD fraction, not average ploidy. The vertical coordinate is reported organ oxygen pressure, with the displayed axis labelled mmHg following the original figure context. The 19 `organ_lut.csv` rows match Supplementary Table3 in labels, WGD, oxygen means, and SDs. One row lacks SD and therefore has no error bar.

Every available vertical SD bar is:

\[
[\bar O-\mathrm{SD},\ \bar O+\mathrm{SD}].
\]

This is neither IQR nor a measurement-count-based confidence interval. Negative lower endpoints are retained as supplied mean±SD instead of clipping to a new statistic; they are not evidence of negative measured oxygen.

### 9.2 Regression and confidence band

The unweighted ordinary-least-squares model is:

\[
O_k=\alpha+\beta_{\mathrm{OLS}}w_k+\varepsilon_k,
\]

where w is WGD fraction. `geom_smooth(method="lm", se=TRUE, level=.95)` draws its default 95% regression-mean confidence band, not a new-observation prediction interval or organ-measurement CI. SD bars do not weight the regression. The coefficient β_OLS is unrelated to Figure7 blockade β.

### 9.3 Row dependence, mappings, and statistical-version discrepancies

Exact source reproduction retains repeated kidney/marrow rows, different brain/lung oxygen summaries linked to shared broad tumor groups, and the author's apparent suprarenal→thyroid mapping problem. These are not 19 independent organ cohorts. No silent organ deduplication, anatomical correction, subtype allocation, or primary-only filtering is applied.

The local 19 rows yield Pearson R=0.616415699813 and two-sided p=0.00494393936871. The article caption reports R=0.62 and p=0.0013; the earlier 20-row data give approximately R=0.6670465 and p=0.0013154. This establishes a source-version/statistic inconsistency without proving its cause. Caption p is not hard-coded as the recomputed result.

The local oxygen tables lack individual measurements, measurement n, complete SD derivation, and explicit unit metadata. mmHg follows the original figure's pressure context rather than a new calibration from model oxygen percent. The regression reproduces the source descriptively and does not establish a new independent cross-organ causal result.

## 10. Figure8C/D: main curves for the two environments

C and D both use day 1000, O₂=0–5%, 51 oxygen coordinates, three starts, and 1000 seeds per coordinate: 153 median/IQR points each. All seeds retain native fitted p. C is continuous in vivo growth without culture passage; D is in vitro with fitted configuration and threshold passage.

All curves are solid: 2N blue `#0072B2`, 4N pink `#CC79A7`, and 1:1 orange `#E69F00`. Q annotations are omitted from the main panels. Colors and shared axes do not alter the statistics.

A's low-oxygen in vivo 2N anchor matches C; its low-oxygen culture anchor matches D; its 20% culture anchor belongs to the extended supplement. D is the exact 0–5% numerical subset of the day-1000 supp_fig8_6 panel.

## 11. supp_fig8_1–5: fixed p and multiple mixture ratios

supp_fig8_1/2 repeat A's empirical comparison and three model conditions, overriding native p to 0.01 or 0.3. Empirical records do not change when model p changes. Each model anchor still pools 1000 seeds, with culture repeats averaged first.

supp_fig8_3 shows native-p in vivo seven-Q curves separately for C01/C02, 500 seeds per coordinate. supp_fig8_4 instead has two matched fixed-p=0.01/0.3 in vivo panels, each pooling 1000; it is not family-separated in the same way.

The seven Q values are 0, 1, 2, 3, 4, 5, ∞. Fixed-0.01 results come from v7 mixture-ratio outputs. Fixed-0.3 results come from the v8 expansion of authenticated complete-parameter pure endpoints. Both use the terminal-count weights in section 5.

supp_fig8_5 is a 2×2 fixed-p culture matrix: C01/C02 columns, 0.01/0.3 rows, seven Q values, and 500 seeds per curve point. All culture mixtures are directly simulated. Its 1428 summary points equal two families×two p values×seven Q values×51 oxygen coordinates. Additional fixed-p pure-2N anchors at 20% belong to the A-style comparison in supplements 1/2, not the supp_fig8_5 curves.

## 12. supp_fig8_6/7: finite-time matrices

Both are chronological row-major 2×5 matrices: days 100, 200, 300, 400, 500 above days 600, 700, 800, 900, 1000. Each element shows oxygen–ploidy median/IQR curves for three starts.

supp_fig8_6 covers in vitro 0–20% at 201 oxygen coordinates, with 6030000 seed×time×condition records and 6030 summary points. supp_fig8_7 covers in vivo 0–5% at 51 coordinates, with 1530000 records and 1530 summary points.

In vivo propagates three initial vectors through consecutive shifted 100-day exponentials. Each culture repeat is one uninterrupted 1000-day process, retaining required times. A snapshot is recorded after any passage on that integer day, not after restarting or reseeding specifically at day 100, 200, or another panel time.

The day-1000 panels recover the preceding 603 culture summary points and 153 main-C points. The previous single-panel supp_fig8_6 is preserved in the time-snapshot audit's `before/`, not used as the current supplement.

## 13. supp_fig8_8/9: family-resolved matrices

Supplement 8 groups supplement 6's existing seed snapshots by family; supplement 9 does the same for supplement 7. C01 is above C02, and each family has a complete 2×5 chronological block with common axes, starts, and colors.

Each family×time×start×oxygen slot computes type-7 quantiles from that family's 500 Y_s values. Supplement 8 therefore contains 12060 summary points and supplement 9 contains 3060. There are no new simulations, p overrides, passage samples, or seed deduplication.

Repooling the same raw snapshots recovers supplements 6/7. Averaging the two family medians or IQR endpoints in supplements 8/9 does not generally recover the pooled summary.

## 14. Numerical examples and independent checks

### 14.1 Native p, day 1000, 1000 pooled seeds

The following N-unit values were independently reaggregated from saved seed endpoints.

| Environment | O₂ (%) | Start/Q | Q25 | Median | Q75 |
| --- | ---: | --- | ---: | ---: | ---: |
| in vivo | 0.5 | 2N / 0 | 1.02247325325 | 1.94841629820 | 2.77625918248 |
| in vivo | 0.5 | 1:1 / 1 | 1.08828221615 | 2.09871161189 | 2.77625918248 |
| in vivo | 0.5 | 4N / ∞ | 2.38555689654 | 2.77376137893 | 2.77625918249 |
| in vitro | 0.5 | 2N / 0 | 1.05054978945 | 1.05113082191 | 1.05177089750 |
| in vitro | 0.5 | 1:1 / 1 | 1.06759770428 | 1.06844011929 | 1.06928354859 |
| in vitro | 0.5 | 4N / ∞ | 3.58926924619 | 3.58950572552 | 3.58960961060 |
| in vitro | 20 | 2N / 0 | 5.99349544150 | 5.99360786483 | 5.99370919149 |
| in vitro | 20 | 1:1 / 1 | 5.99349607836 | 5.99360367889 | 5.99370376412 |
| in vitro | 20 | 4N / ∞ | 5.99350403120 | 5.99360572860 | 5.99372045777 |

A's three model anchors are the in vivo 0.5%/2N, in vitro 0.5%/2N, and in vitro 20%/2N conditions above, matching `Figure8A_display_summary.tsv`.

### 14.2 Scope of this documentation review

Independent quartile calculations cover C, D, day-1000 supplements 6/7, day-1000 family-specific supplements 8/9, and all supplement-5 points: 4002 summary coordinates. Maximum discrepancy from saved summaries is approximately 6.1×10⁻¹⁴, consistent with floating-point/text serialization precision, below the 10⁻¹⁰ verification tolerance.

The 36 empirical A groups were also recomputed from sample-level inputs, with maximum discrepancy approximately 4.9×10⁻¹⁵. The 19-row B Pearson statistics were recalculated. Terminal-population weighting was checked for every native-p in vivo Q=1–5 mixture endpoint.

These are independent calculations on existing results, not new model integrations or passage runs. Within-seed averaging and original per-repeat numerical checks are documented from implemented workers and existing validation receipts; not every production task's full repeat arrays was reaudited locally in this documentation pass.

## 15. Read-only recalculation code

### 15.1 Model seed summaries and temporal/family subsets

This R code reads existing endpoint and summary tables only: no simulation, Rcpp compilation, RNG, or file writes. The main culture table has approximately 1319000 rows and requires corresponding reading memory, but no complete production-task RDS collection.

```r
library(data.table)
workspace <- "/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5"
base <- file.path(workspace, "data/Figures/Figure8")
vivo_run <- file.path(base, "figure8_fitted_family_20261006_v8")
culture_run <- file.path(base, "figure8_invitro_repair_20261007_v11")
time_run <- file.path(base, "figure8_time_snapshots_20261007_v12")
family_run <- file.path(base, "figure8_family_time_snapshots_20261008_v13")
vivo <- fread(file.path(vivo_run, "fitted_mixture_seed_endpoints.tsv"),
  select = c("pair_label", "seed_number", "O2_pct", "Q", "mean_ploidy"))
culture <- fread(file.path(culture_run, "culture_seed_endpoints.tsv"),
  select = c("pair_label", "seed_number", "O2_pct", "Q", "p_setting", "mean_ploidy"))
summarize <- function(z, keys) {
  z[, .(q25 = unname(quantile(mean_ploidy, .25, type = 7)),
        median = median(mean_ploidy),
        q75 = unname(quantile(mean_ploidy, .75, type = 7)), n = .N), by = keys]
}
check <- function(calculated, saved, keys, label) {
  z <- merge(calculated, saved, by = keys, suffixes = c(".calc", ".saved"))
  error <- max(vapply(c("q25", "median", "q75"), function(k) {
    max(abs(z[[paste0(k, ".calc")]] - z[[paste0(k, ".saved")]]))
  }, 0.0))
  stopifnot(nrow(z) == nrow(calculated), nrow(z) == nrow(saved), error < 1e-10)
  cat("PASS", label, nrow(z), "points; max error", error, "\n")
}
keys <- c("O2_pct", "Q")
C <- summarize(vivo[Q %in% c(0, 1, Inf)], keys)
D <- summarize(culture[p_setting == "fitted" & O2_pct <= 5], keys)
extended <- summarize(culture[p_setting == "fitted"], keys)
check(C, fread(file.path(culture_run, "Figure8C_curve_summary.tsv")), keys, "Figure8C")
check(D, fread(file.path(culture_run, "Figure8D_curve_summary.tsv")), keys, "Figure8D")
check(extended, fread(file.path(time_run, "supp_fig8_6_curve_summary.tsv"))[day == 1000],
      keys, "supp_fig8_6/day1000")
check(C, fread(file.path(time_run, "supp_fig8_7_curve_summary.tsv"))[day == 1000],
      keys, "supp_fig8_7/day1000")
family_keys <- c("pair_label", keys)
check(summarize(culture[p_setting == "fitted"], family_keys),
      fread(file.path(family_run, "supp_fig8_8_curve_summary.tsv"))[day == 1000],
      family_keys, "supp_fig8_8/day1000")
check(summarize(vivo[Q %in% c(0, 1, Inf)], family_keys),
      fread(file.path(family_run, "supp_fig8_9_curve_summary.tsv"))[day == 1000],
      family_keys, "supp_fig8_9/day1000")
fixed_keys <- c(family_keys, "p_setting")
fixed_saved <- fread(file.path(culture_run, "supp_fig8_5_curve_summary.tsv"))
fixed_saved[, p_setting := as.character(p_setting)]
check(summarize(culture[p_setting != "fitted" & O2_pct <= 5], fixed_keys),
      fixed_saved, fixed_keys, "supp_fig8_5")
print(rbindlist(list(cbind(context = "in vivo", C[abs(O2_pct - .5) < 1e-10]),
                    cbind(context = "in vitro", extended[O2_pct %in% c(.5, 20)]))),
      digits = 12)
```

### 15.2 Empirical quantiles and B's Pearson statistics

```r
library(data.table)
workspace <- "/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5"
src <- file.path(workspace, "data/share_bundle_figures_3_5")
run <- file.path(workspace, "data/Figures/Figure8/figure8_invitro_repair_20261007_v11")
pc <- fread(file.path(src, "data/Figure3_PCAWG_included_samples.csv"))
dm <- fread(file.path(src, "data/Figure3_DepMap_included_lines.csv"))
stopifnot(nrow(pc) == 2338L, nrow(dm) == 1433L,
          sum(is.finite(dm$ploidy)) == 1185L)
raw <- rbindlist(list(
  pc[is.finite(ploidy), .(site = figure3_label, ploidy, cohort = "PCAWG")],
  dm[is.finite(ploidy), .(site = figure3_label, ploidy, cohort = "DepMap")]))
calculated <- raw[, .(median = median(ploidy),
  q25 = unname(quantile(ploidy, .25, type = 7)),
  q75 = unname(quantile(ploidy, .75, type = 7))), by = .(site, cohort)]
saved <- fread(file.path(run, "empirical_tissue_ploidy_summary.tsv"))
z <- merge(calculated, saved, by = c("site", "cohort"),
           suffixes = c(".calc", ".saved"))
stopifnot(nrow(z) == 36L)
for (k in c("median", "q25", "q75")) {
  stopifnot(max(abs(z[[paste0(k, ".calc")]] - z[[paste0(k, ".saved")]])) < 1e-12)
}
author <- fread(file.path(src, "originals/organ_lut.csv"))
stopifnot(nrow(author) == 19L, sum(is.na(author$SD)) == 1L)
test <- cor.test(author$wgd, author$PtO2.normal.tissue, method = "pearson")
print(c(n = nrow(author), R = unname(test$estimate), p = test$p.value), digits = 12)
fit <- lm(PtO2.normal.tissue ~ wgd, data = author)
print(coef(fit), digits = 12)
```

## 16. Numerical provenance, code, and interpretation limits

### 16.1 Numerical-source index

These runs are under `data/Figures/Figure8/` in the current workspace. Historical internal names do not redefine current figure numbering.

| Run | Current role and principal files |
| --- | --- |
| `figure8_empirical_full1000_20261001_v4` | Initial complete-parameter freeze, empirical tables, and preserved in vivo anchors; its old culture configuration is not current authority |
| `figure8_mixture_ratios_20261006_v7` | Fixed-0.01 seven-Q in vivo `mixture_ratio_seed_endpoints.tsv` and summaries |
| `figure8_fitted_family_20261006_v8` | Native-p `fitted_pure_seed_endpoints.tsv`, `fitted_mixture_seed_endpoints.tsv`, family summaries, and fixed-0.3 expansion `supp_fig8_4_seed_endpoints.tsv` |
| `figure8_invitro_repair_20261007_v11` | Authoritative repaired `culture_seed_endpoints.tsv`, main C/D, A anchors, and supplement 5; 1319000 culture seed/condition endpoints in 6100 tasks |
| `figure8_time_snapshots_20261007_v12` | `supp_fig8_6/7_seed_snapshots.tsv.gz` and corresponding `*_curve_summary.tsv`, from uninterrupted trajectories |
| `figure8_family_time_snapshots_20261008_v13` | `supp_fig8_8/9_curve_summary.tsv`, regrouped from v12 without new simulation |

The 1319000 v11 endpoints comprise 603000 native-p three-start results over 0–20%, 714000 fixed-p seven-start results over 0–5%, and 2000 fixed-p pure-2N anchors at 20%. v12 has 5040 tasks. Production receipts distinguish seed endpoints, time snapshots, and plotted summary coordinates.

Each run's `source_manifest.tsv`, `code_manifest.tsv`, and fingerprint bind inputs and code; `seed_input_manifest.tsv` binds family/seed identities. The culture repair's `authoritative_fit_inputs.tsv` records per-seed configuration sources and hashes, and `full_native_parameter_passthrough.tsv` records parameter-level value checks. These manifests and actual summaries form the evidence chain; directory dates or a `completion.tsv` alone do not establish a calculation's definition.

### 16.2 Implementation index

| Purpose | Code |
| --- | --- |
| Main entry point and repair precedence | [draw_Figure8.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/draw_Figure8.R), [draw_Figure78_invitro_repair.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/draw_Figure78_invitro_repair.R) |
| Complete parameters and initial distributions | [figure78_full_parameters.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_full_parameters.R): `f78v4_freeze`, `f78v4_initial` |
| Native-p finite-time in vivo and mixtures | [figure8_fitted_family_20261006.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_fitted_family_20261006.R): `f8f_pure_endpoint`, `f8f_expand`, `f8f_summary` |
| Fixed-0.01 multiple Q | [figure8_mixture_ratios_20261006.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_mixture_ratios_20261006.R): `f8q_compute` |
| Culture initial mixtures and seed summaries | [figure8_pooled_passage_20261007.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_pooled_passage_20261007.R): `f8p_initial`, `f8p_summary`; its early propagation function is not the current repaired worker |
| Authoritative configuration and current culture worker | [figure78_invitro_repair.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_invitro_repair.R): `f78i_correct_prepared`, `f78i_propagate`, `f78i_worker8` |
| Euler and random threshold sampling | [figure78_invitro_propagator.cpp](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_invitro_propagator.cpp): `f78i_euler_cpp` |
| Repaired summaries and A anchors | [figure78_invitro_repair_postprocessing.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_invitro_repair_postprocessing.R): `f78i_aggregate` |
| Empirical A reconstruction | [figure8_empirical_comparison.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_empirical_comparison.R): `f8e_build`; other older exploratory B tables are not current B |
| Current author-Figure5 input | [figure78_revision_20261005.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_revision_20261005.R): `f78r_empirical5` |
| A/B drawing | [figure78_full_drawing.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure78_full_drawing.R): `f78v4_tissue_plot`, `f78v5_author_figure5`; [figure8_ab_display_20261007.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_ab_display_20261007.R) |
| Finite-time snapshots | [figure8_time_snapshots_20261007.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_time_snapshots_20261007.R): `f8t_vivo`, `f8t_vitro`, `f8t_aggregate` |
| Family regrouping and type-7 quantiles | [figure8_family_time_snapshots_20261008.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure8_family_time_snapshots_20261008.R): `f8tf_summary`, `f8tf_prepare` |
| Native-p override and paths | [figure6_robustness.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure6_robustness.R): `figure6_force_p_misseg`, `f6r_paths` |
| Data-derived common passage protocol | [figure6_full_range_q10.R](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/Figures/util/analysis/figure6_full_range_q10.R): `f6g_canonical_passage_rule` |

Key configurations follow. Historical panel names remain subject to the current entry-point and version precedence in section 1.

| Configuration | Role |
| --- | --- |
| [figure8_empirical_comparison.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_empirical_comparison.json) | Initial complete-parameter task settings, not authority for superseded culture dynamics |
| [figure8_mixture_ratios_20261006.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_mixture_ratios_20261006.json) | Fixed-0.01 seven-Q in vivo extension |
| [figure8_three_start_20261006.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_three_start_20261006.json) | Fixed-0.3 pure-endpoint three-start parent source; later v8 supplies the seven-Q expansion |
| [figure8_fitted_family_20261006.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_fitted_family_20261006.json) | Native-p in vivo source, families, and mixtures |
| [figure78_invitro_repair_20261007.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure78_invitro_repair_20261007.json) | Current culture-configuration repair, oxygen ranges, passage repeats, and summaries |
| [figure8_time_snapshots_20261007.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_time_snapshots_20261007.json) | Times, oxygen ranges, and snapshot definitions for supplements 6/7 |
| [figure8_family_time_snapshots_20261008.json](/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/Code/config/figure8_family_time_snapshots_20261008.json) | Family regrouping and type-7 quantiles for supplements 8/9 |

### 16.3 Validation records and locally reproducible scope

Existing production receipts are under these audit roots and are not rewritten by documentation:

```text
/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/audit/figure78_invitro_repair_20261007
/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/audit/figure8_time_snapshots_20261007
/Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/iteration5/audit/figure8_family_time_snapshots_20261008
```

They cover authoritative parameters/configurations, reference replay of experimental lineages, timestep checks, snapshot-storage invariance of endpoints/passages, complete seed identities/counts/parameters, endpoint–time-matrix consistency, charts, and publication hashes. Existing family-split receipts contain 23 numerical, 44 chart, 52 PDF/PNG, and two visual-review checks: 121 total. These are original production checks, not fresh model runs during documentation.

Complete seed endpoints, compressed snapshots, and summary tables are local, but complete v11/v12 production-task RDS collections are not. An unfinished checkpoint is not a completed production task. Section 15 therefore recalculates from available complete tables without claiming new inspection of all repeat arrays. Auditing that level requires retrieving the original production archive rather than silently rerunning simulations. This documentation pass uses no SSH/HPC submission, refitting, resampling, or model simulation.

### 16.4 Biological and numerical limitations

- Repairing culture configuration does not guarantee reproduction of an observed hypoxic 2N ploidy increase. Native-p, 0.5%-oxygen, day-1000 2N and 1:1 culture predictions remain close to the lower boundary and are not adjusted to enforce a trend.
- The current 22–154 grid represents 1–7N. Original repair diagnostics show substantial low-oxygen 2N/1:1 terminal mass in states 22–24; some high-oxygen results approach higher ploidy. No new expanded-grid sensitivity experiment establishes that plateaus are grid-independent.
- Saved `minimum_log_population` is the minimum across repeats' terminal values, not the minimum over the full trajectory. Terminal boundary mass does not describe boundary mass at every time.
- Fractional expected cell abundances are not observations of that many single cells. Population-size, denominator, and boundary diagnostics remain relevant to curve interpretation.
- Optimizer seeds are neither biological replicates nor posterior samples. Passage repeats are not additional fitted seeds. Model IQR, empirical IQR, oxygen SD, and OLS confidence bands are not interchangeable.
- Imposed model O₂ percent is not calibrated to organ mmHg. Original mapping and row-dependence limitations remain, and claims of model–empirical agreement or causality require additional evidence.
