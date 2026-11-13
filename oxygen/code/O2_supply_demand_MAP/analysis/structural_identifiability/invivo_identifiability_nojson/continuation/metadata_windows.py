"""One bounded archived-output validation dispatch; no model or search calls."""
from pathlib import Path
import argparse
import re
import uuid
import workflow as w


def successful(j):
    return j['state'] == 'COMPLETED' and j['exit_code'] == '0:0'


def checked_endpoints(rows):
    by = {j['task']: j for j in rows}
    ready = set()
    for j in rows:
        if not successful(j) or j['kind'] == 'engineering':
            continue
        if j['mode'] == 'local' and j['kind'] == 'baseline':
            names = [j['task'] + suffix for suffix in ('_fresh', '_directions')]
            if all(name in by and successful(by[name]) for name in names):
                ready.add(j['task'])
        elif j['mode'] in ('global', 'profile_fresh'):
            ready.add(j['task'])
    return ready


def advance(study):
    s = Path(study)
    rows = w.snapshot(s)
    audits = [j for j in rows if j['mode'] == 'metadata_audit']
    if any(not w.closed(j) for j in audits):
        print('Await current allocated metadata validation.'); return
    reviewed = {r['task'] for r in w.read_table(s / 'tables/reviewed_failures.tsv')}
    if any(w.closed(j) and not successful(j) and j['task'] not in reviewed for j in rows):
        raise ValueError('Review execution failure before archived-output validation')
    previous = [j for j in audits if successful(j)]
    old = []
    if previous:
        spec = w.kv(s / 'dispatch' / previous[-1]['task'] / 'spec.tsv')
        directory = s / 'collections' / spec['metadata_collection']
        if w.sha(directory / 'manifest.tsv') != spec['metadata_collection_sha256']:
            raise ValueError('Last validation input collection changed')
        w.verify_manifest(directory, directory / 'manifest.tsv')
        old = w.read_table(directory / 'jobs.tsv')
    if checked_endpoints(rows) <= checked_endpoints(old):
        print('No newly closed checked endpoint/population; no redundant audit.'); return
    state = w.State(s)
    with state.locked():
        current = state.load()
        if [(j['task'], j['state']) for j in current] != [(j['task'], j['state']) for j in rows]:
            raise ValueError('Controller generation changed; reconcile before audit dispatch')
        generation = (s / 'state/CURRENT').read_text().strip()
        identity = 'c_' + uuid.uuid4().hex
        directory = s / 'collections' / identity
        directory.mkdir()
        (directory / 'jobs.tsv').write_bytes((s / 'state' / generation / 'jobs.tsv').read_bytes())
        (directory / 'raw_accounting.txt').write_bytes((s / 'tables/raw_accounting.txt').read_bytes())
        w.table(directory / 'origin.tsv', [dict(key='state_generation', value=generation),
                dict(key='captured_utc', value=w.utc())], ['key', 'value'])
        w.table(directory / 'manifest.tsv', [dict(path=p.name, sha256=w.sha(p))
                for p in sorted(directory.iterdir())], ['path', 'sha256'])
    numbers = [int(re.fullmatch('a_metadata_w([0-9]+)', j['task']).group(1))
               for j in audits]
    task = 'a_metadata_w' + str(max(numbers or [0]) + 1).zfill(4)
    release = w.kv(s / 'tables/metadata_execution_release.tsv')['release']
    w.submit(s, task, mode='metadata_audit', role='archived_output_validation',
             release=release, kind='engineering', root_id='metadata_reference',
             root_family='reference', selector='reference', cores=1,
             wall_seconds=3600, margin_seconds=300, max_calls=1,
             metadata_collection=identity,
             metadata_collection_sha256=w.sha(directory / 'manifest.tsv'))
    w.snapshot(s)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('study', type=Path)
    advance(p.parse_args().study.resolve())
