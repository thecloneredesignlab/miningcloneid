# Independent in-vivo fixed-oxygen eFAST interpretation

The independent in-vivo design varied 14 operator parameters over the fitted
transformed bounds (uniform in transformed space; log-uniform on natural scale
for `log10_*` parameters). At each of the 201 Figure 4 oxygen values (0–5%),
the fixed-oxygen operator supplied dominant mean ploidy and its leading
eigenvalue, the asymptotic net live-cell growth rate. SALib 1.5.2 used
`M=4`, `N=129,257,513`, and two independent phase seeds at each resolution.
Figure 4B fitted-solution Spearman correlations are retained as a separate
[directional panel](figures/figure4b_spearman_direction.png). The four FAST
heatmaps share [one aligned panel](figures/efast_four_panel.png).

**The proposed low-to-high oxygen mechanism shift is not supported as a
general statement.** At low oxygen (0–1%), death-specific parameters
(`mu_hp`, `gamma_mu`) did not dominate ploidy: their mean per-parameter total
effect was 0.379, compared with 0.472 for growth, 0.465 for missegregation,
and 0.404 for buffering. The largest low-oxygen ploidy first-order index was
`buffer_beta` (0.189); `mu_hp` and `gamma_mu` were 0.021 and 0.022.
At high oxygen (3–5%), the ploidy missegregation total-effect mean was 0.493
and death fell to 0.298. Buffering was 0.362, below its low-oxygen value;
individual buffering parameters changed in different directions. Shared
oxygen-stress parameters (`O2_crit`, `n_O`) are not death-specific.

Growth showed a narrower partial shift. Death had the largest group-mean
first-order index at low oxygen (0.088; `mu_hp` alone 0.165), but fell to
0.006 at high oxygen. High-oxygen growth instead favored growth parameters
(group mean first-order 0.105; `lam_max` 0.316) and baseline missegregation
(`p_mis_base` 0.137). Buffering's group-mean first-order index stayed at or
below 0.002 in both bands.

**Numerical limit.** Ploidy first-order indices were relatively stable:
at `N=513`, the 90th percentile of the absolute phase-repeat difference was
0.046, and the `N=257` to `N=513` mean-index difference was 0.035.
Ploidy total effects remained less stable (corresponding 90th percentiles
0.371 and 0.312), so fine parameter ordering in the ploidy ST heatmap is
exploratory. Growth total effects were more stable (0.073 and 0.062).
At `N=513`, about 30–32% of high-oxygen operator evaluations had a leading
spectral gap below `1e-4`, versus about 5.7–5.9% at low oxygen. Near-degenerate
leading modes may contribute to unstable asymptotic ploidy estimates; this
is an inference, not a demonstrated cause of the FAST variability.

The [index table](indices.tsv), [convergence table](convergence.tsv),
[mechanism band table](mechanism_band_summary.tsv),
[replicate band table](mechanism_band_replicates.tsv), and
[spectral-gap audit](spectral_gap_qc.tsv) support these statements. Group
values are **means of parameter indices**, not additive fractions of variance:
total effects overlap through interactions. Sensitivities are conditional on
the stated independent input ranges and distributions, and have no sign;
the Figure 4B correlation panel provides direction across fitted solutions.
