#!/bin/bash -l
# Validate the new node and the pinned SIF before any scientific worker starts.
set -euo pipefail
study=$1
task=$2
attempt=$3
if [[ "$task" != P && ! "$task" =~ ^N_[A-F]_(import|bench|de1b|de[2-6]|local|pattern|verify[12]|replay[12][ab])$ ]]; then
  exit 78
fi
uid=$(id -u)
[[ "$uid" == 107865 ]]
entry=$(getent passwd "$uid")
IFS=: read -r account password mapped_uid remaining <<< "$entry"
[[ "$account" == 4482173 && "$mapped_uid" == "$uid" ]]
node=$(hostname -s)
out="$study/verification/$task"
mkdir -p "$out"
if [[ "$attempt" == 2 ]]; then
  [[ "$task" == N_A_local && "$node" != hpctpa3pc0027 ]]
  image=/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact.sif
  digest=$(sha256sum "$image")
  [[ "${digest%% *}" == dfbd7ebe1209814ca25358e56ae4d5553c6dc19852f5567ce369ab23bebe3372 ]]
  container_uid=$(apptainer exec --cleanenv "$image" /bin/sh -c 'id -u')
  [[ "$container_uid" == "$uid" ]]
  temp="$out/.pending_node_preflight_${SLURM_JOB_ID:?}"
  printf 'status\tjob_id\tnode\tuid\thost_user\tcontainer_uid\nPASS\t%s\t%s\t%s\t%s\t%s\n' \
    "$SLURM_JOB_ID" "$node" "$uid" "$account" "$container_uid" > "$temp"
  mv "$temp" "$out/node_preflight.tsv"
else
  [[ "$attempt" == 1 ]]
fi
exec "$study/release/runner.sh" "$study" "$task"
