#!/usr/bin/env bash
set -euo pipefail
repo=${EFAST_REPO_ROOT:-/share/lab_crd/taoli/Project/soft_couping_eFAST}
out=${EFAST_OUT_DIR:-$repo/oxygen/results/eFAST/neighborhood10pct}
sif=${EFAST_SIF:-/share/lab_crd/taoli/Docker/o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif}
code="$repo/oxygen/code/O2_supply_demand_MAP/analysis/eFAST"
exec 8>"$out/slurm/eigen_safe_submit.lock"
flock -n 8 || exit 1
[[ ! -s "$out/slurm/numerical_recovery_safe_jobs.tsv" ]] || { echo "Safe-method jobs already submitted" >&2; exit 1; }
ids=$(python3 - "$out/slurm/numerical_recovery_jobs.tsv" <<'PY'
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
print(next(r['job_id'] for r in reversed(rows) if r['role']=='repair'))
print(next(r['job_id'] for r in reversed(rows) if r['role']=='seed_summary'))
PY
)
parent=${ids%%$'\n'*}; summary=${ids#*$'\n'}
summary_state=$(squeue -h -j "$summary" -o '%T' | sort -u)
[[ "$summary_state" == PENDING ]] || { echo "Existing summary array must still be pending for a dependency update" >&2; exit 1; }
python3 - "$out" "$parent" <<'PY'
import pathlib,subprocess,sys
root=pathlib.Path(sys.argv[1])
data=subprocess.check_output(['sacct','-n','-X','-j',sys.argv[2],'--state=FAILED','--format=JobID%32','-P'],universal_newlines=True)
ids=sorted(int(x.strip().split('|')[0].split('_')[1]) for x in data.splitlines() if '_' in x)
assert ids,'No failed targets'
(root/'slurm/eigen_recovery/safe_task_ids.txt').write_text(','.join(map(str,ids))+'\n')
PY
apptainer exec --cleanenv "$sif" env OPENBLAS_NUM_THREADS=1 python3 "$code/recover_eigen_safe.py" prepare --out-dir "$out" > "$out/slurm/eigen_recovery/safe_prepare.log"
spec=$(python3 - "$out/slurm/numerical_recovery_safe_plan.json" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))['array_spec'])
PY
)
job=$(sbatch --parsable --qos=xxlarge --time=12:00:00 --cpus-per-task=1 --mem=4G --array="$spec" \
  --job-name=efast_safe_recovery --output="$out/slurm/logs/safe_recovery_%A_%a.out" \
  --error="$out/slurm/logs/safe_recovery_%A_%a.err" \
  --export="ALL,EFAST_REPO_ROOT=$repo,EFAST_OUT_DIR=$out,EFAST_SIF=$sif,EFAST_RECOVERY_ROLE=repair" "$code/run_eigen_safe.sbatch")
job=${job%%;*}
printf 'role\tjob_id\tarray_spec\tqos\ttime\tcpus\tmem\tdependency\tsubmitted_utc\n' > "$out/slurm/numerical_recovery_safe_jobs.tsv"
printf 'repair_safe\t%s\t%s\txxlarge\t12:00:00\t1\t4G\t\t%s\n' "$job" "$spec" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$out/slurm/numerical_recovery_safe_jobs.tsv"
scontrol update "JobId=$summary" "Dependency=afterany:$parent:$job"
python3 - "$out" "$parent" "$summary" "$job" <<'PY'
import csv,datetime,json,pathlib,sys
root=pathlib.Path(sys.argv[1]); parent,summary,safe=sys.argv[2:]
path=root/'slurm/numerical_recovery_jobs.tsv'
rows=list(csv.DictReader(path.open(),delimiter='\t'))
for row in rows:
    if row['job_id']==summary: row['dependency']='afterany:'+parent+':'+safe
with path.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
record=dict(recorded_utc=datetime.datetime.utcnow().isoformat()+'Z',seed_summary_job_id=summary,
    original_dependency='afterany:'+parent,updated_dependency='afterany:'+parent+':'+safe,
    new_recovery_job_id=safe,reason='Include the safeguarded M-matrix recovery before the existing seed summaries.')
(root/'slurm/numerical_recovery_dependency_update.json').write_text(json.dumps(record,indent=2)+'\n')
PY
cat "$out/slurm/numerical_recovery_safe_jobs.tsv"
