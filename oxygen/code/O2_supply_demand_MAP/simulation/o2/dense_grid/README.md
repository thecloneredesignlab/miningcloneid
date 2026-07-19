# Dense fixed-O2 simulation

`generate_dense_grid_simulation_tables.R` evaluates fitted parameter sets on a
dense O2 grid. `part=monotonicity` materializes attractor features;
`part=initial_ploidy` materializes 2N/4N finite-time trajectories. HPC modes are
`build_tasks`, `run_tasks`, `merge_daily_seed`, and `merge`; task lists and chunks
remain below the selected output directory's `hpc/` and `simulation_tables/`.

This producer never classifies curves or draws figures.
