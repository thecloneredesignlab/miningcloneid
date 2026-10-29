#!/usr/bin/env bash
# Submission is explicit. There is deliberately no array concurrency throttle.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
source "${script_dir}/../config/figure78_hpc_paths.env"
campaign=;phase=;submit=no
for arg in "$@"; do
  case "$arg" in
    --campaign-root=*) campaign="${arg#*=}" ;;
    --phase=small|--phase=full) phase="${arg#*=}" ;;
    --submit=yes) submit=yes ;;
    *) echo "Use --campaign-root=ABSOLUTE --phase=small|full --submit=yes" >&2; exit 2 ;;
  esac
done
[[ "$campaign" == /* && -n "$phase" && "$submit" == yes ]] || exit 2
array="$(python3 - "$campaign" "$phase" <<'PY'
import csv,json,sys
from pathlib import Path
p=Path(sys.argv[1]);phase=sys.argv[2]
meta=json.loads((p/'prepared/preparation.json').read_text())
if phase=='full':
    gate=json.loads((p/'provenance/small_acceptance.json').read_text())
    assert gate['status']=='passed' and gate['fingerprint']==meta['fingerprint']
with (p/'prepared/task_manifest.tsv').open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
ids=[]
for r in rows:
    if phase=='small' and r['small_acceptance']!='TRUE':continue
    i=int(r['task_id'])
    if (p/'receipts'/('task_%05d.tsv'%i)).exists():continue
    ids.append(str(i))
assert ids,'No pending task; collect existing results instead'
print(','.join(ids))
PY
)"
mkdir -p "${campaign}/logs"
exec sbatch --parsable --partition=red --qos=xxlarge --cpus-per-task=1 --mem=16G --time=12:00:00 \
  --array="$array" --export=ALL --job-name="figure7_fitted_${phase}" \
  --output="${campaign}/logs/${phase}_%A_%a.out" \
  "${script_dir}/run_Figure7_fitted_worker.sh" "$FIGURE78_WORKSPACE" "$campaign" "$phase"
