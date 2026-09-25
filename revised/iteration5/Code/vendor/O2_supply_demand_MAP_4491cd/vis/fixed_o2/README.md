# Fixed-O2 visualization

This directory renders figures from existing numerical and analysis tables. It
does not fit parameters, load model objects, estimate trajectories, or generate
missing simulations.

## Files

- `render_fixed_o2_figures.R`: figure CLI. It requires `--simulation_dir` and
  `--analysis_dir`, writes PDFs under `--vis_dir`, and records only a figure
  provenance manifest. `--parts` supports `attractors`, `counterfactual`,
  `validation`, and `agreement`.
- `fixed_o2_plot_utils.R`: plotting functions extracted from the historical
  fixed-O2 monolith. It contains attractor/reliability panels, counterfactual
  trajectory panels, expm/eigen/Euler validation panels, phase-plane panels,
  and analytical-versus-simulation scatter/residual/Bland-Altman panels.

The visualization CLI treats TSV inputs as immutable and writes no numerical
result tables.
