#!/usr/bin/env python3
"""Summarize completed evaluation logs without rerunning the culture model."""
import argparse
import collections
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re


def finite_number(value):
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def error_category(message):
    """Group protocol reasons while keeping numeric candidate details in raw logs."""
    cohort = re.search(r'(?:^|[ ;])cohort=([^;]+)', message)
    cohort = cohort.group(1).strip() if cohort else ''
    if message.startswith('protocol_infeasible:'):
        parts = [p.strip() for p in message.split(':', 1)[1].split(';')]
        reason = '; '.join(p for p in parts if p and '=' not in p)
        return ('protocol_infeasible', cohort, reason or 'unspecified_protocol_constraint')
    pattern = r'[-+]?(?:\d*\.\d+|\d+\.?\d*)(?:[eE][-+]?\d+)?'
    family = re.sub(pattern, '<number>', message.split(':', 1)[0])[:100]
    reason = re.sub(pattern, '<number>', message)[:240]
    return (family, cohort, reason)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('study', type=Path)
parser.add_argument('--snapshot')
args = parser.parse_args(); study = args.study.resolve()
state_path = study/'report_inputs'/args.snapshot/'review_state_snapshot.json' if args.snapshot else study/'state.json'
raw_state = state_path.read_bytes(); state = json.loads(raw_state); output = []; errors = []
for tid, task in state['tasks'].items():
    if task['status'] != 'accepted':
        continue
    owner = study/'fits'/tid
    files = sorted((owner/'evaluation_logs').glob('*.tsv'))
    if (owner/'local_evaluations.tsv').exists():
        files.append(owner/'local_evaluations.tsv')
    if not files:
        raise RuntimeError('Accepted task has no evaluation log: ' + tid)
    groups = {}
    for file in files:
        with file.open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                phase = row['phase']
                group = groups.setdefault(phase, dict(count=0, finite=0, errors=collections.Counter(),
                    examples={}, elapsed_sum=0., days=set(), first_days=None, changed_days=0, files=set()))
                group['count'] += 1; group['files'].add(file.name)
                group['finite'] += finite_number(row['objective']) is not None
                group['elapsed_sum'] += finite_number(row.get('seconds')) or 0.
                if row.get('error'):
                    category = error_category(row['error'])
                    group['errors'][category] += 1
                    group['examples'].setdefault(category, row['error'][:500])
                days = row.get('selected_days')
                if days and finite_number(row['objective']) is not None:
                    if group['first_days'] is None: group['first_days'] = days
                    group['days'].add(days)
                    group['changed_days'] += days != group['first_days']
    for phase, group in groups.items():
        output.append(dict(task_id=tid, point_id=task['point_id'], domain=task['domain'],
            mode=task['mode'], pure_independent=task['pure_independent'], phase=phase,
            logged_evaluations=group['count'], finite_objectives=group['finite'],
            finite_fraction=group['finite']/group['count'],
            evaluations_with_error=sum(group['errors'].values()),
            evaluation_seconds_sum=group['elapsed_sum'], log_files=len(group['files']),
            distinct_selected_day_vectors=len(group['days']) if group['days'] else None,
            changed_selected_day_vectors_vs_first=group['changed_days'] if group['days'] else None,
            error_categories=len(group['errors']),
            inline_error_categories_truncated=len(group['errors']) > 10,
            error_counts_json=json.dumps({' | '.join(k):v for k,v in group['errors'].most_common(10)}, sort_keys=True)))
        for (family, cohort, reason), count in group['errors'].items():
            errors.append(dict(task_id=tid, phase=phase, error_family=family, cohort=cohort,
                reason=reason, count=count, example_error=group['examples'][(family,cohort,reason)]))
path = study/'tables/review_evaluation_log_summary.tsv'
with path.open('w') as stream:
    fields=list(output[0]) if output else ['task_id']
    writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(output)
with (study/'tables/review_evaluation_error_categories.tsv').open('w') as stream:
    fields=['task_id','phase','error_family','cohort','reason','count','example_error']
    writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(errors)
record=dict(generated_utc=datetime.now(timezone.utc).isoformat(), snapshot=args.snapshot,
    review_state_sha256=hashlib.sha256(raw_state).hexdigest(),
    accepted_tasks=sum(t['status']=='accepted' for t in state['tasks'].values()), rows=len(output),
    error_category_rows=len(errors),
    scope='Completed accepted tasks only; global DE, preflight and post-DE calls remain separate. Inline error counts show at most ten normalized categories; the error-category table retains all category counts and raw logs retain every full message. Evaluation seconds sum worker-call wall times and are not allocated core-hours. Selected-day changes accompany parameter probes and do not establish a causal source of objective kinks.')
(study/'reports/search_log_audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
