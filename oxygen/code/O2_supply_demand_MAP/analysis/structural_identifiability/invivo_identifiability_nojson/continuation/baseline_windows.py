"""One baseline-local coordination pass, gated by completed SIF record audits.

Global exploration, profile release and scientific closure remain explicit
decisions. No local-window, sweep, call-lifetime or aggregate-time ceiling.
"""
from pathlib import Path
import argparse
import metadata_windows as metadata
import next_windows as continuation
import workflow as w


def covered(current_rows, audited_rows):
    return metadata.checked_endpoints(current_rows) <= metadata.checked_endpoints(audited_rows)


def advance(study):
    s = Path(study).resolve()
    w.advance(s, 'checks')
    metadata.advance(s)
    rows = w.State(s).load()
    audits = [j for j in rows if j['mode'] == 'metadata_audit']
    if any(not w.closed(j) for j in audits):
        print('COORDINATOR await_allocated_metadata_audit'); return
    audits = [j for j in audits if metadata.successful(j)]
    if not audits:
        raise ValueError('No completed allocated metadata proof')
    latest = audits[-1]
    payload = s / 'tests' / latest['task']
    w.verify_payload(payload)
    if w.kv(payload / 'result.tsv')['pass'] != 'true':
        raise ValueError('Internal metadata validation failed')
    spec = w.kv(s / 'dispatch' / latest['task'] / 'spec.tsv')
    collection = s / 'collections' / spec['metadata_collection']
    if w.sha(collection / 'manifest.tsv') != spec['metadata_collection_sha256']:
        raise ValueError('Metadata validation collection changed')
    w.verify_manifest(collection, collection / 'manifest.tsv')
    audited = w.read_table(collection / 'jobs.tsv')
    if not covered(rows, audited):
        print('COORDINATOR await_current_endpoint_audit'); return
    for task in metadata.checked_endpoints(rows):
        j = next(j for j in rows if j['task'] == task)
        targets = [task]
        if j['kind'] == 'baseline' and j['mode'] == 'local':
            targets += [task + '_directions']
        for target in targets:
            parent = next(j for j in rows if j['task'] == target)
            root = s / ('fits' if parent['mode'] in ('local', 'global') else 'tests') / target
            proof = w.kv(s / 'tables/record_audits' / (target + '_summary.tsv'))
            if proof['pass'] != 'true' or proof['payload_sha256'] != w.sha(root / 'payload_manifest.tsv'):
                raise ValueError('Missing current native-record proof: ' + target)
    gate = w.kv(s / 'tables/baseline_gate.tsv')
    if gate['input_manifest_sha256'] != w.sha(s / 'manifests/input_manifest.tsv'):
        raise ValueError('Baseline review belongs to different frozen input')
    w.verify_manifest(s, s / gate['evidence_manifest'])
    continuation.advance(s)
    for row in w.read_table(s / 'tables/strategy_review.tsv'):
        print('SCIENTIFIC_DECISION', row['policy'], row['decision'], row['batch'])


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('study', type=Path)
    advance(p.parse_args().study)
