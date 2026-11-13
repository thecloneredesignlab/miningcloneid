#!/usr/bin/env python3
"""Validate the bounded error-category export against the archived v03 metrics."""
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

study=Path(sys.argv[1]); snapshot=sys.argv[2]; code=Path(__file__).resolve().parent
manifest=json.loads((code/'source_manifest.json').read_text())
for name,digest in manifest['sources'].items():
    if hashlib.sha256((code/name).read_bytes()).hexdigest()!=digest:
        raise RuntimeError('Search audit source changed: '+name)
subprocess.check_call(['python3',str(code/'audit_culture_review_search.py'),str(study),'--snapshot',snapshot])
csv.field_size_limit(100*1024**2)
def read(path):
    with path.open() as stream:return list(csv.DictReader(stream,delimiter='\t'))
old=read(study/'tables/archive_v03/review_evaluation_log_summary.tsv')
new=read(study/'tables/review_evaluation_log_summary.tsv')
indexed={(r['task_id'],r['phase']):r for r in new}
if len(old)!=len(new):raise RuntimeError('Evaluation summary row count changed')
fields=['logged_evaluations','finite_objectives','finite_fraction','evaluations_with_error',
        'evaluation_seconds_sum','log_files','distinct_selected_day_vectors',
        'changed_selected_day_vectors_vs_first']
for row in old:
    candidate=indexed[(row['task_id'],row['phase'])]
    if any(row[k]!=candidate[k]for k in fields):raise RuntimeError('Existing evaluation metrics changed')
errors=read(study/'tables/review_evaluation_error_categories.tsv'); counts={}
for row in errors:
    key=(row['task_id'],row['phase']);counts[key]=counts.get(key,0)+int(row['count'])
for key,row in indexed.items():
    if counts.get(key,0)!=int(row['evaluations_with_error']):raise RuntimeError('Error category totals do not reconcile')
maximum=max(len(row['error_counts_json'])for row in new)
if maximum>32767:raise RuntimeError('Inline error summary exceeds a common spreadsheet-cell limit')
csv.field_size_limit(131072)
read(study/'tables/review_evaluation_log_summary.tsv')
record=dict(generated_utc=datetime.now(timezone.utc).isoformat(),pass_format_check=True,
    snapshot=snapshot,evaluation_metrics_unchanged=True,summary_rows=len(new),
    error_category_rows=len(errors),maximum_inline_error_characters=maximum,
    previous_maximum_inline_error_characters=max(len(row['error_counts_json'])for row in old),
    source_manifest=manifest,optimization_performed=False,
    explanation='Numeric candidate values made raw error strings nearly unique and created oversized JSON cells. Bounded normalized categories now have a separate count table. The v03 export and every raw log are preserved; scores, parameters, predictions and evaluation metrics are unchanged.')
(study/'reports/search_audit_format_correction.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
