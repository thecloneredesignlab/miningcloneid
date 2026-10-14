#!/usr/bin/env bash
# Queue the full pilot while smoke waits for the user's xxlarge CPU allowance.
# Analysis remains gated by Slurm afterok; the frozen numerical code is unchanged.
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
baseline=${EFAST_BASELINE_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
out=${EFAST_OUT_DIR:-$baseline/pilot/N4097}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
exec 8>"$out/slurm/submit.lock"
flock -n 8 || { echo 'Another N4097 submission is in progress' >&2; exit 1; }
[[ ! -s "$out/slurm/trajectory_job_id.txt" ]] || { echo 'Trajectory array already submitted' >&2; exit 1; }
smoke=$(cat "$out/slurm/smoke_job_id.txt")
[[ "$smoke" =~ ^[0-9]+$ ]] || exit 2
smoke_state=$(scontrol show job -o "$smoke")
apptainer exec --cleanenv "$sif" python3 - "$out" "$code" "$sif" "$smoke" "$smoke_state" <<'PY'
import hashlib,json,pathlib,subprocess,sys
root,code=map(pathlib.Path,sys.argv[1:3]); sys.path.insert(0,str(code))
import resolution_4097 as run
plan=run.context(root,verify_sources=True)
assert pathlib.Path(sys.argv[3]).resolve()==pathlib.Path(plan['sif'])
fields=dict(x.split('=',1) for x in sys.argv[5].split() if '=' in x)
state=fields['JobState']
assert fields['JobName']=='efast_4097_smoke' and fields['QOS']=='xxlarge'
assert fields['UserId'].startswith('4482173(')
assert fields['Command']==str(code/'run_resolution_4097.sbatch')
assert state in ('PENDING','RUNNING','CONFIGURING','COMPLETING','COMPLETED'), state
receipt=root/'smoke/smoke_receipt.json'
plan_hash=hashlib.sha256((root/'slurm/submission_plan.json').read_bytes()).hexdigest()
if receipt.exists():
    r=json.loads(receipt.read_text())
    assert r['status']=='passed' and r['submission_plan_sha256']==plan_hash
if state=='COMPLETED': assert fields['ExitCode']=='0:0' and receipt.exists()
submitter=code/'queue_resolution_4097_after_smoke.sh'
request=dict(submission_plan_sha256=plan_hash,smoke_job_id=sys.argv[4],smoke_state_at_queue=state,
    smoke_receipt_present_at_queue=receipt.exists(),main_dependency='afterok:'+sys.argv[4],
    reason='Queue while preflight waits for xxlarge capacity; evaluate only after successful smoke',
    submitter_sha256=hashlib.sha256(submitter.read_bytes()).hexdigest(),
    submission_git_commit=subprocess.check_output(['git','-C',str(code),'rev-parse','HEAD'],universal_newlines=True).strip())
(root/'slurm/deferred_submission.json').write_text(json.dumps(request,indent=2)+'\n')
PY
common=(--parsable --qos=xxlarge --time=12:00:00)
exports="ALL,EFAST_REPO_ROOT=$repo,EFAST_OUT_DIR=$out,EFAST_BASELINE_DIR=$baseline,EFAST_SIF=$sif"
record() {
  printf '%s\t%s\t%s\txxlarge\t12:00:00\t%s\t%s\t%s\t%s\n' \
    "$1" "$2" "$3" "$4" "$5" "$6" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$out/slurm/submission_jobs.tsv"
}
[[ -s "$out/slurm/submission_jobs.tsv" ]] || { echo 'Missing smoke submission record' >&2; exit 1; }
trajectory=$(sbatch "${common[@]}" --array=1-420 --cpus-per-task=1 --mem=4G \
  --dependency="afterok:$smoke" --job-name=efast_4097_trajectory \
  --output="$out/slurm/logs/trajectory_%A_%a.out" --error="$out/slurm/logs/trajectory_%A_%a.err" \
  --export="$exports,EFAST_ARRAY_ROLE=trajectory" "$code/run_resolution_4097.sbatch")
trajectory=${trajectory%%;*}
printf '%s\n' "$trajectory" > "$out/slurm/trajectory_job_id.txt"
record trajectory "$trajectory" 1-420 1 4G "afterok:$smoke"
summary=$(sbatch "${common[@]}" --array=1-6 --cpus-per-task=1 --mem=8G \
  --dependency="afterany:$trajectory" --job-name=efast_4097_summary \
  --output="$out/slurm/logs/summary_%A_%a.out" --error="$out/slurm/logs/summary_%A_%a.err" \
  --export="$exports,EFAST_ARRAY_ROLE=seed_summary" "$code/run_resolution_4097.sbatch")
summary=${summary%%;*}
record seed_summary "$summary" 1-6 1 8G "afterany:$trajectory"
final=$(sbatch "${common[@]}" --cpus-per-task=4 --mem=32G --dependency="afterany:$summary" \
  --job-name=efast_4097_compare --output="$out/slurm/logs/finalize_%j.out" \
  --error="$out/slurm/logs/finalize_%j.err" \
  --export="$exports,EFAST_ARRAY_ROLE=finalize" "$code/run_resolution_4097.sbatch")
final=${final%%;*}
record finalize "$final" '' 4 32G "afterany:$summary"
cat "$out/slurm/submission_jobs.tsv"
