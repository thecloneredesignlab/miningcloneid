#!/usr/bin/env bash
# Launch this driver with nohup on RED; it survives laptop/SSH disconnection.
set -euo pipefail
repo=${SOBOL_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
out=${SOBOL_OUT_DIR:-$repo/oxygen/results/eFAST/sobol/full_N4096}
rec=${SOBOL_RECOVERY_DIR:-$out/recovery_20261018}
sif=${SOBOL_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
mkdir -p "$rec/slurm/logs"
exec 8>"$rec/slurm/submission.lock"
flock -n 8 || { echo 'Recovery submission already active' >&2; exit 1; }
[[ ! -e "$rec/slurm/submission_jobs.tsv" ]] || { echo 'Existing recovery ledger; inspect before resubmitting' >&2; exit 1; }
snapshot="$rec/slurm/code_snapshot"
[[ ! -e "$snapshot" ]] || { echo 'Recovery snapshot exists; inspect interrupted preparation' >&2; exit 1; }
mkdir -p "$snapshot"
revision=$(git -C "$repo" rev-parse HEAD)
git -C "$repo" archive "$revision" oxygen/code/O2_supply_demand_MAP | tar -xf - -C "$snapshot"
python3 - "$snapshot" <<'PY'
import pathlib,sys
p=pathlib.Path(sys.argv[1])
for q in p.rglob('*'): q.chmod(0o555 if q.is_dir() else 0o444)
p.chmod(0o555)
PY
code="$snapshot/oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol"
printf 'role\tjob_id\tarray_spec\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$rec/slurm/submission_jobs.tsv"
common=(--parsable --qos=xxlarge --time=12:00:00 --cpus-per-task=1)
exports="ALL,SOBOL_OUT_DIR=$out,SOBOL_RECOVERY_DIR=$rec,SOBOL_SIF=$sif,RECOVERY_REVISION=$revision"
submit() {
  local role=$1 array=$2 mem=$3 dependency=$4 job
  local options=(--mem="$mem" --job-name="sobol_resume_$role" --output="$rec/slurm/logs/${role}_%A_%a.out" --error="$rec/slurm/logs/${role}_%A_%a.err" --export="$exports,SOBOL_ROLE=$role")
  [[ -z "$array" ]] || options+=(--array="$array")
  [[ -z "$dependency" ]] || options+=(--dependency="$dependency")
  job=$(sbatch "${common[@]}" "${options[@]}" "$code/run_neighborhood_recovery.sbatch")
  job=${job%%;*}
  printf '%s\t%s\t%s\txxlarge\t12:00:00\t1\t%s\t%s\t%s\n' "$role" "$job" "$array" "$mem" "$dependency" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$rec/slurm/submission_jobs.tsv"
  echo "Submitted $role $job" >&2
  printf '%s' "$job"
}
prepare=$(submit prepare '' 8G '')
while :; do
  state=$(sacct -X -n -P -j "$prepare" --format=State | head -1 | cut -d '|' -f 1)
  case "$state" in
    COMPLETED) break;;
    FAILED|TIMEOUT|CANCELLED*|NODE_FAIL|OUT_OF_MEMORY) echo "Preparation ended $state; no recovery computations submitted" >&2; exit 1;;
  esac
  sleep 30
done
[[ -s "$rec/plan.json" ]] || { echo 'Preparation did not create a plan' >&2; exit 1; }
RECOVERY_PLAN_SHA256=$(sha256sum "$rec/plan.json"); RECOVERY_PLAN_SHA256=${RECOVERY_PLAN_SHA256%% *}
exports="$exports,RECOVERY_PLAN_SHA256=$RECOVERY_PLAN_SHA256"
read -r ready total seeds model_ids < <(python3 - "$rec/plan.json" <<'PY'
import json,sys
p=json.load(open(sys.argv[1]))
print(p['n_ready_tasks'],len(p['tasks']),len(p['cases']),','.join(str(x['task_id']) for x in p['model_tasks']))
PY
)
model=$(submit model "$model_ids" 4G '')
ready_job=$(submit shard "1-$ready" 8G '')
waiting_job=$(submit shard "$((ready+1))-$total" 8G "afterok:$model")
cases=$(submit case "1-$seeds" 8G "afterany:$ready_job:$waiting_job")
audit=$(submit audit '' 8G "afterany:$cases")
aggregate=$(submit aggregate '' 16G "afterok:$audit")
responses=$(submit response 1-54 4G "afterok:$aggregate")
finish=$(submit finish '' 8G "afterok:$responses")
# Only these two original, pending, irrecoverably blocked workflow stages are replaced.
python3 - <<'PY'
import subprocess
text=subprocess.check_output(['squeue','-h','-j','20839426,20839427','-o','%T'],universal_newlines=True)
if any(s.strip()!='PENDING' for s in text.splitlines()): raise SystemExit('Old downstream stage is no longer pending; inspect before cancelling')
PY
scancel 20839426 20839427
python3 - "$rec" "$revision" <<'PY'
import csv,datetime,json,pathlib,sys
r=pathlib.Path(sys.argv[1]); p=json.load(open(str(r/'plan.json')))
with open(str(r/'slurm/submission_jobs.tsv')) as h: jobs=list(csv.DictReader(h,delimiter='\t'))
s=dict(status='submitted',git_commit=sys.argv[2],n_pending_seeds=len(p['cases']),n_shards=len(p['tasks']),
       n_model_blocks=len(p['model_tasks']),jobs=jobs,replaced_pending_jobs=['20839426','20839427'],
       convergence_policy='diagnostic_only',submitted_utc=datetime.datetime.utcnow().isoformat()+'Z')
(r/'slurm/launch_receipt.json').write_text(json.dumps(s,indent=2)+'\n')
(r/'status.json').write_text(json.dumps(s,indent=2)+'\n')
print(json.dumps(s,indent=2))
PY
