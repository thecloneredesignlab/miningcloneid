"""Post-allocation billing, including the last archived-output audit itself.

Pure integer/archived metadata work. This does not grant scientific closure.
"""
from pathlib import Path
import argparse
import workflow as w


def summarize(jobs, raw):
    if not jobs or any(not w.closed(j) or not j['job_id'] for j in jobs):
        raise ValueError('Final accounting requires every owned parent closed')
    if len({j['job_id'] for j in jobs}) != len(jobs):
        raise ValueError('Duplicate parent identity')
    roles = {}
    for j in jobs:
        r = raw[j['job_id']]
        seconds = int(r[3]) * int(r[4])
        if (r[1], r[2], int(r[3]), seconds) != (
                j['state'], j['exit_code'], int(j['cores']),
                int(j['actual_core_seconds'])):
            raise ValueError('Incoherent parent accounting: ' + j['task'])
        successful = j['state'] == 'COMPLETED' and j['exit_code'] == '0:0'
        role = roles.setdefault(j['role'], dict(
            role=j['role'], parents=0, successful_parents=0,
            unsuccessful_parents=0, actual_core_seconds=0))
        role['parents'] += 1
        role['successful_parents'] += int(successful)
        role['unsuccessful_parents'] += int(not successful)
        role['actual_core_seconds'] += seconds
    return list(roles.values())


def finalize(study):
    s = Path(study)
    jobs = w.read_table(s / 'tables/jobs.tsv')
    raw = {r.split('|')[0]: r.split('|') for r in
           (s / 'tables/raw_accounting.txt').read_text().splitlines() if r}
    roles = summarize(jobs, raw)
    actual = sum(r['actual_core_seconds'] for r in roles)
    w.table(s / 'tables/final_role_accounting.tsv', roles, [
        'role', 'parents', 'successful_parents', 'unsuccessful_parents',
        'actual_core_seconds'])
    w.typed(s / 'tables/final_accounting.tsv', dict(
        all_owned_parents_closed=True, owned_parents=len(jobs),
        new_actual_core_seconds=actual,
        historical_lifetime_core_seconds=301846,
        lifetime_actual_core_seconds=301846 + actual,
        historical_formal_plus_new_core_seconds=271355 + actual,
        active_parents=0, aggregate_core_hour_ceiling='none',
        accounting_model_evaluations=0))
    lines = [
        '# Final owned allocation accounting', '', 'Updated ' + w.utc() + '.', '',
        'All owned parent allocations are closed. This billing snapshot includes '
        'the final metadata validation allocation, which was still active when '
        'its scientific report was generated. Scientific acceptance is recorded '
        'separately in the baseline, profile and closure gates.', '',
        'New actual: {} core-seconds ({:.8f} core-hours).'.format(actual, actual / 3600), '',
        'Historical lifetime: 301846 core-seconds (83.84611111 core-hours). '
        'Lifetime actual: {} core-seconds ({:.8f} core-hours).'.format(
            301846 + actual, (301846 + actual) / 3600), '',
        'Historical formal plus all new work: {} core-seconds ({:.8f} core-hours). '
        'This is an accounting view, with no aggregate ceiling.'.format(
            271355 + actual, (271355 + actual) / 3600), '',
        '| Role | Parents | Successful | Unsuccessful | Actual core-hours |',
        '|---|---:|---:|---:|---:|']
    for r in roles:
        lines.append('| {role} | {parents} | {successful_parents} | '
                     '{unsuccessful_parents} | {hours:.8f} |'.format(
                         hours=r['actual_core_seconds'] / 3600, **r))
    lines += ['', 'Parent AllocCPUS multiplied by ElapsedRaw is the charging '
              'authority. Failed, cancelled and partial work remains charged. '
              'Child timers and inherited optimizer counters are not added again. '
              'Culture accounting is outside this owned ledger and remains unchanged.', '']
    w.atomic_text(s / 'reports/final_accounting.md', '\n'.join(lines))
    return actual


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('study', type=Path)
    print('Final new actual core-seconds:', finalize(p.parse_args().study.resolve()))
