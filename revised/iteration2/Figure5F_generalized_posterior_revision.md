# Figure 5F family-conditioned generalized-posterior audit

> **SUPERSEDED FOR THE ACTIVE MANUSCRIPT AND FIGURE (2026-08-17).** The active
> Figure 5F no longer uses generalized-posterior or HPC sampling products. This
> file is retained only as historical provenance. The current contract and
> revision record are `data/Figures/Figure5/figure5f_chart_contract.md` and
> `Figure5F_prior_optimizer_revision.md`.

Status: the previous common-target production chain was stopped. A redesigned
family-conditioned workflow is staged for production. All writable outputs,
checkpoints, logs, figures, tables, and manuscript changes remain under
`revised/iteration2`; the external joint-fit result tree is read-only.

## Scientific correction

C01, C02, and C03 are treated as three distinct generalized-posterior targets,
not as dispersed initializations of one common target. For family (f), the
target at learning temperature (T) is

\[
\pi_{T,f}(\theta\mid D)\propto
\exp[-L_{\mathrm{data}}(\theta)/T]
\exp[-L_{\mathrm{regularization}}(\theta)]
I(\theta\in\mathcal B\cap\mathcal V_f).
\]

The three frozen supports \(\mathcal V_f\) are the selected-MAP Voronoi cells
defined by squared Euclidean distance in the exact bound-normalized
18-dimensional in-vivo mechanistic parameter space. The coordinates are the
four ordinary in-vivo oxygen coordinates plus the 14 in-vivo members of the
paired mechanistic parameters. Ties use the fixed order C01, C02, C03. The
partition excludes t-SNE axes, optimizer-population labels, in-vitro nuisance
coordinates, and fitted endpoint frequency.

The selected MAP for each family remains the pair chosen by the original
separate-in-vivo objective:

- C01: `tsne_vi_seed366_C01Sc01_vt_seed10`, selected joint seed 472.
- C02: `tsne_vi_seed25_C02Sc01_vt_seed10`, selected joint seed 497.
- C03: `tsne_vi_seed311_C03Sc02_vt_seed10`, selected joint seed 18.

These optimizer endpoints initialize and define the frozen numerical basins;
they are not posterior samples or biological replicates.

## Objective and prior

`L_data` is the weighted in-vivo fitted-data loss plus the weighted in-vitro
objective. `L_regularization` contains the configured in-vivo Gaussian
soft-prior terms, Welsch cross-context coupling penalty, and active feasibility
penalty. At `T=1`, their sum reproduces the audited historical joint objective.
Only `L_data` is tempered.

The objective replay against miningcloneid commit
`83953a874401e42cd176432786f889a896adc959` agrees with the three saved selected
MAP objectives to a maximum absolute difference of `8.53e-14` (required
tolerance `1e-8`). This verifies those audited inputs; it does not claim global
identity with an unrecoverable dirty working tree.

The displayed prior is also conditioned separately on C01, C02, and C03. It is
generated from the exact transformed bounds, Welsch coupling term, and
configured in-vivo Gaussian soft-prior terms, then rejected unless its in-vivo
mechanistic coordinates lie in the requested frozen basin. Production retains
100,000 prior draws per family and parameter. An independent 20,000-draw
replicate per family is used only for quantile-stability QA and is not plotted.

## Independent dynamic sampling

Each family has two independent replica-exchange ladders, R1 and R2, with
temperatures `0.5,1,2,4,8`. `T=1` is primary, `T=0.5/2` are sensitivity
targets, and `T=4/8` support within-family exchange. Cross-family MAP jumps are
not valid for these conditional targets and are disabled. A proposal crossing
the frozen family boundary is rejected before the expensive objective is
evaluated.

The initial target is 20,000 iterations with 3,000 warm-up iterations and no
thinning. Diagnostics are evaluated separately for each family. A family passes
one checkpoint only when all of the following hold:

- all 14 `T=1` ratio diagnostics have rank-normalized split R-hat at most 1.01
  and bulk/tail ESS at least 400;
- all 40 `T=1` sampling coordinates have R-hat at most 1.05 and bulk ESS at
  least 100;
- all retained states at `T=0.5,1,2` remain inside that family's frozen basin;
- every local acceptance rate is finite;
- every adjacent replica-swap acceptance rate is at least 0.05.

A family stops only after two consecutive passing checkpoints. If it does not
pass, only that family's R1/R2 pair is extended by 10,000 iterations from its
two exact atomic checkpoints. Prior states, retained draws, adaptive scales,
acceptance counters, swap counters, and random-number states are reused; no
completed objective evaluation is repeated. The frozen safety maximum is
200,000 iterations per family. Reaching it without two consecutive passes
halts that family for review rather than continuing indefinitely.

The previous common-target files and checkpoints are not reused by this
workflow. The family-conditioned namespace is
`data/Figures/Figure5/generalized_posterior_family_conditioned/`, with
checkpoints under
`tmp/figure5f_generalized_posterior_family_conditioned/`.

## HPC execution contract

The production launcher is
`Code/Figures/Figure5F/hpc/submit_Figure5F_family_conditioned.sh`.
It uses `--qos=small`, `--time=2-00:00:00`, loads `ml R/4.4`, and invokes
`/app/eb/software/R/4.4.2-gfbf-2024a/bin/Rscript` directly. Each active family
sampling job uses two Slurm tasks times five CPUs, one task for each independent
ladder and five workers for its temperature proposals. Aggregation and control
are file-backed and family-specific.

After all three families pass twice, one locked combine job validates family
identities, iteration counts, target versions, configuration signatures,
dimensions, and finite values. It then runs the frozen input audit, builds the
Figure 5F products and Supplementary Table, renders Figure 5, and writes only
under `revised/iteration2`.

## Interpretation boundary

Optimizer endpoint densities and 5th--95th percentile spans remain descriptive
numerical-search summaries. They are not generalized-posterior draws,
confidence intervals, bootstrap samples, or biological replication. A Figure
5F contrast can be called concentrated only relative to its matching
family-conditioned prior and after all three family-specific convergence and
boundary checks pass. This is an operational assessment of practical
identifiability conditional on the model, loss scale, bounds, regularization,
family conditioning, and data; it is not a proof of structural identifiability.

## Output contract

- Root primary draws: `data/Figures/Figure5/figure5f_posterior_draws.tsv`.
- Family sampler outputs:
  `data/Figures/Figure5/generalized_posterior_family_conditioned/`.
- Standalone panel:
  `Figures/figure5_subpanels/figure5f_prior_generalized_posterior.{png,pdf}`.
- Assembled Figure 5: `Figures/assembled_fig5.{png,pdf}`.
- Manuscript copy: `manuscript/Figures/assembled_fig5.{png,pdf}`.
- Supplementary values:
  `manuscript/tables/data/supp_figure5f_generalized_posterior_values.tsv`.
- Supplementary TeX:
  `manuscript/tables/supp_figure5f_generalized_posterior.tex`.

Final numerical results must remain unpopulated until all three distinct
family-conditioned posteriors pass the frozen release contract.
