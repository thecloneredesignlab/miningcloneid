# O2 container execution

This directory provides execution wrappers for the environment defined and
locked at repository root under `container/`.

- `hpc/` mirrors every file under `../hpc/` and runs R/Python through the
  verified Apptainer SIF while leaving Slurm orchestration on the host.
- `local/` runs canonical fitting and analysis through the published Docker
  image.
- `image_runtime_lock.tsv` records the exact artifacts verified by these
  wrappers.

No certificate, private key, Docker credential, SSH directory, or host home is
copied or mounted by these wrappers.
