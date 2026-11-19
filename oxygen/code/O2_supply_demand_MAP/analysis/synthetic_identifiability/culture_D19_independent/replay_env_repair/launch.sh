#!/bin/bash -l
# Supply the native allocation identity through the pinned clean environment.
set -euo pipefail
study=$1
phase=$2
case "$phase" in D_*|P) ;; *) printf '%s\n' 'Unregistered recovery phase' >&2; exit 78;; esac
export APPTAINERENV_SLURM_JOB_ID="${SLURM_JOB_ID:?Expected a native allocation ID}"
if [ "$phase" = D_B_replay1b ]; then
  image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
  received=$(apptainer exec --cleanenv "$image" /bin/sh -c 'printf "%s" "${SLURM_JOB_ID:-}"')
  test "$received" = "$SLURM_JOB_ID"
  printf 'D19_REPAIR_CLEANENV_JOB_ID_PASS %s\n' "$received"
fi
exec "$study/release/runner.sh" "$study" "$phase"
