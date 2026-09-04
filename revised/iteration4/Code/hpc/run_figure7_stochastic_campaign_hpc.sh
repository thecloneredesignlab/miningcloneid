#!/usr/bin/env bash
# Background campaign: validate/pilot -> lock R -> complete grid -> draw -> seal.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
ITERATION_ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd -P)"
cd "${ITERATION_ROOT}"
RUN_ID="${1:?Supply a new Figure7 stochastic campaign run id}"
BASELINE="${2:?Supply the validated continuous baseline directory}"
[[ "${RUN_ID}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]{0,64}$ ]] || exit 2
[[ "$(hostname -s)" == hpctpa3pc0028 ]] || { echo "Wrong compute node" >&2; exit 2; }
DATA_ROOT="${ITERATION_ROOT}/data/Figures/Figure7/fixed_pmisseg_v1/finite_time_full_q10_runs"
PILOT_ID="${RUN_ID}_pilot"
CAMPAIGN_DIR="${ITERATION_ROOT}/audit/hpc_figure7_stochastic/${RUN_ID}"
mkdir -p "${CAMPAIGN_DIR}"
CAMPAIGN_STAGE=INITIALIZING
campaign_status() {
  printf 'run_id\tstage\texit_code\tupdated_at\n%s\t%s\t%s\t%s\n' \
    "${RUN_ID}" "$1" "$2" "$(date -Iseconds)" > "${CAMPAIGN_DIR}/status.tsv"
}
trap 'rc=$?; if (( rc != 0 )); then campaign_status "FAILED_${CAMPAIGN_STAGE}" "$rc"; fi' EXIT
CAMPAIGN_STAGE=PILOT
campaign_status "${CAMPAIGN_STAGE}" 0
bash "${SCRIPT_DIR}/run_figure7_full_range_hpc.sh" --n-core=60 --pilot-only \
  "--run-id=${PILOT_ID}" "--reuse-continuous-run=${BASELINE}" --replicates=20
BUDGET="${DATA_ROOT}/${PILOT_ID}/stochastic_pilot_budget.tsv"
[[ -f "${BUDGET}" ]] || { echo "Pilot budget absent" >&2; exit 1; }
REPLICATES="$(awk -F '\t' 'NR==2 {print $1}' "${BUDGET}")"
[[ "${REPLICATES}" =~ ^(20|50|100)$ ]] || { echo "Unapproved pilot repeat count" >&2; exit 1; }
CAMPAIGN_STAGE=FULL_GRID_AND_RENDER
campaign_status "${CAMPAIGN_STAGE}" 0
printf 'Locked stochastic replicate count: %s\n' "${REPLICATES}"
bash "${SCRIPT_DIR}/run_figure7_full_range_hpc.sh" --n-core=60 --o2-chunk-size=10 \
  "--run-id=${RUN_ID}" "--reuse-continuous-run=${BASELINE}" "--replicates=${REPLICATES}"
for name in stochastic_pilot_budget.tsv stochastic_pilot_results.tsv stochastic_pilot_selection.tsv; do
  cp -- "${DATA_ROOT}/${PILOT_ID}/${name}" "${DATA_ROOT}/${RUN_ID}/${name}"
done
CAMPAIGN_STAGE=SEAL_RETURN
campaign_status "${CAMPAIGN_STAGE}" 0
python3 "${SCRIPT_DIR}/figure7_return_manifest.py" --root "${ITERATION_ROOT}" --run-id "${RUN_ID}"
campaign_status COMPLETE 0
printf 'Campaign complete. Return manifest: audit/hpc_figure7_full_range/%s/return_files.txt\n' "${RUN_ID}"
