# Local Docker execution

These wrappers execute the canonical O2 fitting and analysis code in
`zafiro/o2_supply_demand_map:r44`. They do not copy source into the image: the
current repository checkout is mounted read/write at its existing absolute
path, so outputs and logs remain on the host.

The exact OCI index, amd64 manifest, and corresponding HPC SIF hashes are
recorded in `../image_runtime_lock.tsv`.

## Verify the image

```bash
docker pull zafiro/o2_supply_demand_map:r44
bash oxygen/code/O2_supply_demand_MAP/Docker/local/verify_environment.sh
```

The verifier checks all locked R and repository-wide Python dependency
versions, then reports R, Python, Git, and aria2 versions from the image.

## Full fitting

The production-like wrapper defaults to `--fitting_mode=all`: it runs 500
in-vivo seeds, then 500 in-vitro seeds, automatically passes both generated
result directories into primary-cluster selection, and fits 500 joint seeds
per selected pair.

```bash
bash oxygen/code/O2_supply_demand_MAP/Docker/local/run_full_fit_docker.sh \
  --dry_run=TRUE
```

After reviewing the dry run, remove `--dry_run=TRUE` to execute. Seed counts can
be set explicitly:

```bash
bash oxygen/code/O2_supply_demand_MAP/Docker/local/run_full_fit_docker.sh \
  --invivo_total_seeds=500 \
  --invitro_total_seeds=500 \
  --joint_total_seeds=500 \
  --n_cores=12
```

For individual modes or any canonical fit option:

```bash
bash oxygen/code/O2_supply_demand_MAP/Docker/local/run_o2_fit_docker.sh \
  --fitting_mode=invivo \
  --total_seeds=10 \
  --n_cores=8
```

## Full analysis

Supply the selected or completed seed directories. Each provided scope runs
simulation, analysis, visualization, and report generation. The full
best-fit-parameter feature workflow runs by default.

```bash
bash oxygen/code/O2_supply_demand_MAP/Docker/local/run_full_analysis_docker.sh \
  --invivo_fit_dir=/absolute/path/to/invivo/seed1 \
  --invitro_fit_dir=/absolute/path/to/invitro/seed1 \
  --joint_fit_dir=/absolute/path/to/joint/seed1 \
  --dry_run=TRUE
```

Run a single post-fit scope with:

```bash
bash oxygen/code/O2_supply_demand_MAP/Docker/local/run_postfit_docker.sh \
  --fit_dir=/absolute/path/to/seed1 \
  --scope=joint \
  --dry_run=TRUE
```

For multi-warmup joint coupling analysis, append both
`--joint_result_root=...` and `--joint_output_root=...` to the full-analysis
command.

## Runtime overrides

- `O2SD_DOCKER_IMAGE`: image tag; defaults to
  `zafiro/o2_supply_demand_map:r44`.
- `O2SD_DOCKER_PLATFORM`: defaults to `linux/amd64`.
- `O2SD_DOCKER_RCPP_CACHE`: host temporary directory mounted over the protected
  model's Rcpp cache; defaults outside the repository.
- `O2SD_DOCKER_BINDS`: comma-separated additional Docker volume
  specifications, such as `/data:/data:ro,/results:/results`.
- `O2SD_INVIVO_TOTAL_SEEDS`, `O2SD_INVITRO_TOTAL_SEEDS`, and
  `O2SD_JOINT_TOTAL_SEEDS`: separate-fit seed counts and joint seeds per
  selected pair for the full-fit wrapper.

Use additional binds for any input or output path outside this repository.
These scripts do not mount certificates, SSH keys, Docker credentials, or the
host home directory. Only project-prefixed analytical controls and
thread-count environment variables are forwarded; host certificate variables
are not forwarded.
