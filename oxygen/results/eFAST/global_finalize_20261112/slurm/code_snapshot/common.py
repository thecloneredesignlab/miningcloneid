"""Global-only eFAST contracts; host-side code supports Python 3.6+."""
import contextlib
import csv
import datetime
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path

INITIAL_NS = (2049, 4097, 8193, 16385)
ACTIVE = ('lam_max', 'p_mis_base', 'p_misseg', 'k_o_mis', 'buffer_smax',
          'buffer_beta', 'buffer_n_exp', 'p_wgd', 'alpha_o2', 'gamma_growth',
          'mu_hp', 'gamma_mu', 'O2_crit', 'n_O')
OUTPUTS = ('dominant_mean_ploidy', 'dominant_growth_rate')
METRICS = [('S1', p, '') for p in ACTIVE] + [('ST', p, '') for p in ACTIVE]
CHUNK = 1024
REPEATS = 5
M = 4


def now():
    return datetime.datetime.utcnow().isoformat() + 'Z'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for part in iter(lambda: handle.read(8 * 1024 ** 2), b''):
            digest.update(part)
    return digest.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    os.replace(str(tmp), str(path))


def read(path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(str(path), 'rt', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def table(path, rows, fields=None):
    path = Path(path)
    rows = list(rows)
    fields = fields or list(rows[0])
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == '.gz' else open
    tmp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    with opener(str(tmp), 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(str(tmp), str(path))


@contextlib.contextmanager
def locked(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def valid_n(n, m=M):
    # This is the SALib design requirement, not a maximum N.
    if n <= 4 * m * m or n % 2 != 1:
        raise ValueError('Use an odd eFAST N > 4*M^2')
    return n


def model_count(n):
    valid_n(n)
    return ((n + CHUNK - 1) // CHUNK) * REPEATS * len(ACTIVE)


def decode(n, task):
    if not 1 <= task <= model_count(n):
        raise ValueError('Task outside this global eFAST design')
    chunks = (n + CHUNK - 1) // CHUNK
    design, segment = divmod(task - 1, chunks)
    repeat, parameter = divmod(design, len(ACTIVE))
    begin = segment * CHUNK
    return repeat + 1, parameter, begin, min(n, begin + CHUNK)


def decode_analysis(task):
    if not 1 <= task <= REPEATS * len(ACTIVE):
        raise ValueError('Analysis task outside the 70 phase/trajectory combinations')
    repeat, parameter = divmod(task - 1, len(ACTIVE))
    return repeat + 1, parameter


def model_path(root, n, repeat, parameter, begin):
    return Path(root) / 'efast' / ('N%d' % n) / 'blocks' / ('R%02d' % repeat) / ACTIVE[parameter] / ('from%06d' % begin) / 'outputs.npz'


def next_n(n):
    return 2 * valid_n(n) - 1
