# Independent in-vivo eFAST companion to Figure 4B

`efast.py` generates SALib 1.5.2 FAST trajectories and analyzes the two outputs
of the **fixed-oxygen** operator. `evaluate_fixed_o2.R` sources the canonical
`simulation/o2/fixed_o2/run_fixed_o2_simulation.R` and calls its
`fixo2_dominant_attractor_one` for every parameter vector and oxygen value.
It verifies the unchanged seed25 operator against three Figure 4 reference
points before evaluating a design. The leading real eigenvalue is the
asymptotic net live-cell growth rate in day^-1; the normalized leading
eigenvector determines dominant mean ploidy.

## Input and sampling contract

- Independent **in-vivo** fitting run: `parameter_table.csv` supplies exact
  transformed optimizer bounds; `seed25/best_params.tsv` and
  `seed25/fit_config.rds` supply fixed model settings and reference validation.
- Figure 4 data: `fixed_o2_dominant_ploidy_201grid.tsv` verifies the 201-point
  oxygen grid from 0 to 5% in 0.025% increments and the seed25 reference;
  `continuous_ploidy_spearman_by_o2.tsv` supplies the existing directional
  correlation panel, copied to the result folder with a recorded SHA-256.
- Fourteen parameters present in the fixed-oxygen operator are sampled
  independently. Transformed `log10_*` bounds are sampled uniformly in log10
  space (log-uniform on natural scale); identity bounds are sampled uniformly.
  These priors cover documented fit ranges. The 500 optimizer solutions do not
  constitute a FAST design and are not resampled as one.
- `o2_S0`, `kappa_O`, `eta_o2`, and `k_clear` are Figure 4B parameters whose
  direct effects are structurally absent when oxygen is externally fixed;
  their sensitivity indices are undefined for this 14-factor design. The
  heatmaps display these as grey `N/A` rows to retain Figure 4B alignment.
- `O2_crit` and `n_O` control the shared oxygen-stress response used by growth,
  death, and missegregation functions. They are reported separately from the
  death-specific `mu_hp` and `gamma_mu` when comparing mechanisms.
- The 201 oxygen points use the same FAST parameter vectors within each
  resolution/replicate. FAST trajectory order is preserved, including if
  numerical evaluation fails: failures stop analysis rather than silently
  removing rows.

## Execution

On `hpctpa3pc0028`, with the SALib 1.5.2 SIF and direct shell access:

```bash
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/run_efast.sh pilot
bash oxygen/code/O2_supply_demand_MAP/analysis/eFAST/run_efast.sh full
```

The pilot uses N=129, M=4, one phase seed, four Figure 4 oxygen values. The
full analysis uses N=129 and 257, M=4, two independent phase seeds at each
resolution, all 201 oxygen values. Each design contains `14*N` model vectors.
`EFAST_WORKERS` controls local fork workers (default 4). The runner enforces
the named node and does not submit a Slurm job. The full run is restartable:
an existing nonempty `outputs.tsv.gz` is retained. For a fresh rerun of one
design, remove only that design's output first.
The runner clears inherited R settings so that `Matrix` and other packages
come from the validated SIF rather than an incompatible HPC home library.

## Outputs

`oxygen/results/eFAST/` contains `parameter_ranges.tsv`, the copied Figure 4B
correlation source, full FAST samples and operator output tables under `runs/`,
`indices.tsv` with S1/ST per replicate, `convergence.tsv` with replicate ranges
and change from the highest resolution, parameter and mechanism summaries for
0–1% and 3–5% oxygen, four index heatmaps, a directional
correlation panel, and `interpretation.md`. The sample metadata records seeds,
distributions, input hashes, oxygen points, and SALib version.
The top-level `environment.tsv` records the evaluator image and Git commit;
`analysis_manifest.tsv` records the postprocessing commit and summary hashes.
S1 and ST are variance contributions without a sign; the correlation panel retains
direction. S1/ST are conditional on these chosen input ranges/distributions,
not a posterior uncertainty decomposition. Summing ST across parameters
double counts shared interactions.
