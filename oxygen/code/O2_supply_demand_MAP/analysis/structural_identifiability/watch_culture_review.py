#!/usr/bin/env python3
"""Read-only state watcher; never submits or changes a job."""
import argparse
import collections
from datetime import datetime, timezone
import json
from pathlib import Path
import time

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('study', type=Path)
parser.add_argument('--duration', type=int, default=14400)
args = parser.parse_args()
started = time.monotonic(); last = None; printed = 0
while time.monotonic() - started < args.duration:
    try:
        state = json.loads((args.study / 'state.json').read_text())
        result = dict(status=state['status'], tasks=len(state['tasks']),
            counts=dict(collections.Counter(t['status'] for t in state['tasks'].values())),
            baseline_gate=state.get('baseline_gate'),
            reserve_global=state['reserve_global_attempts'],
            reserve_local=state['reserve_local_attempts'],
            controller_error=state.get('controller_error'),
            scientific_hold=state.get('scientific_hold'),
            last_event=state['events'][-1:])
        key = json.dumps(result, sort_keys=True)
        if key != last or time.monotonic() - printed >= 900:
            print(datetime.now(timezone.utc).isoformat(), key, flush=True)
            last = key; printed = time.monotonic()
        if state['tasks'] and all(t['status'] in ('accepted', 'failed', 'rejected', 'screen_failed')
                for t in state['tasks'].values()) and state['status'] not in (
                    'review_running', 'validation_pending', 'controller_validation_pending'):
            break
    except Exception as error:
        print('READ_ONLY_WATCH_ERROR', repr(error), flush=True)
    time.sleep(60)
print('READ_ONLY_WATCH_WINDOW_CLOSED', flush=True)
