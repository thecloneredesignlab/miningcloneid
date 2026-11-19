import tempfile
from pathlib import Path

from controller import ATTEMPT_CAPS, CAP, ENVELOPES, ORIGINS, PHASES, benchmark, phase_buckets, site_seconds

assert sum(ENVELOPES.values()) == CAP
assert 4 * sum(ATTEMPT_CAPS.values()) + 400 == 336400
assert sum(200 * (g + 1) for g in (1, 59, 60, 60, 60, 60, 60)) == ATTEMPT_CAPS['global']
assert ORIGINS == ('D_A','D_B','D_D','D_F') and 'joint' in PHASES and 'restore' in PHASES
assert site_seconds('12:00:00|') == 43200
assert site_seconds('1-02:03:04') == 93784
assert benchmark([dict(seconds='0.8', feasible='TRUE')] * 200 +
                 [dict(seconds='0.2', feasible='FALSE')] * 200) == 1.0
assert benchmark([dict(seconds='2', feasible='TRUE')] * 400) == 2.5
with tempfile.TemporaryDirectory() as d:
    f = Path(d) / 'calls.tsv'
    f.write_text('call\tevent\tbucket\n1\tSTART\tglobal\n1\tEND\tglobal\n2\tSTART\tlocal\n')
    sums = phase_buckets(f)
    assert sums['global'] == sums['local'] == 1 and sum(sums.values()) == 2
print('Allocation and interrupted-attempt fixtures passed')
