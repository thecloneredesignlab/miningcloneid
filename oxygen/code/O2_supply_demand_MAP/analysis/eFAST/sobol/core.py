"""Independent-input second-order Sobol designs and immutable result helpers."""
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import efast
import neighborhood

ACTIVE = efast.ACTIVE
OUTPUTS = efast.OUTPUTS
PAIRS = [(i, j) for i in range(len(ACTIVE)) for j in range(i + 1, len(ACTIVE))]
PILOT_SEEDS = ['seed25', 'seed464', 'seed165', 'seed2', 'seed214', 'seed138', 'seed66', 'seed395']
OXYGEN = [0., .5, 2.5, 5.]
RESOLUTIONS = [256, 512, 1024]
REPEATS = 5
SIF_SHA256 = '0b60f6cdab8a91f6bbeea1ab6cf02dd79660a295f24bc6e8ad6c5fed04d982ad'
STEP = 2 * len(ACTIVE) + 2
BLOCKS = ['A'] + ['AB_' + p for p in ACTIVE] + ['BA_' + p for p in ACTIVE] + ['B']


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    return efast.sha256(path)


def read(path):
    return efast.read_table(path)


def table(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        if not rows:
            raise ValueError('Explicit columns required for an empty table')
        fields = list(rows[0])
    efast.write_table(path, fields, rows)


def save_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.incomplete')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def load_json(path):
    return json.loads(Path(path).read_text())


def bounds_for(ranges, best=None):
    if best is None:
        return [dict(r) for r in ranges]
    return neighborhood.neighborhood_bounds(ranges, best)


def problem(ranges):
    if [r['parameter'] for r in ranges] != list(ACTIVE):
        raise ValueError('Parameter order must match the canonical fixed-O2 design')
    return dict(num_vars=len(ACTIVE), names=list(ACTIVE),
                bounds=[[float(r['encoded_lower']), float(r['encoded_upper'])] for r in ranges])


def samples(ranges, n, scramble_seed):
    from SALib.sample import sobol
    if n < 2 or n & (n - 1):
        raise ValueError('Sobol N must be a power of two')
    encoded = sobol.sample(problem(ranges), n, calc_second_order=True,
                           scramble=True, seed=scramble_seed)
    natural = encoded.copy()
    for j, r in enumerate(ranges):
        if r['transform'] == 'log10':
            natural[:, j] = 10. ** encoded[:, j]
        if np.any(natural[:, j] < float(r['natural_lower']) - 1e-12) or np.any(
                natural[:, j] > float(r['natural_upper']) + 1e-12):
            raise ValueError('Sample outside natural bounds')
    return natural


def analyze(ranges, y, seed=1):
    from SALib.analyze import sobol
    y = np.asarray(y, dtype=float)
    if y.ndim != 1 or len(y) % STEP or not np.isfinite(y).all():
        raise ValueError('Incomplete/nonfinite ordered second-order design')
    variance = float(np.var(y))
    tolerance = (1e-12 * max(1., float(np.max(np.abs(y))))) ** 2
    if variance <= tolerance:
        return {k: np.full(len(ACTIVE) if k.startswith(('S1', 'ST')) else (len(ACTIVE), len(ACTIVE)),
                           np.nan) for k in ('S1', 'ST', 'S2', 'S1_conf', 'ST_conf', 'S2_conf')}, variance
    result = sobol.analyze(problem(ranges), y, calc_second_order=True, num_resamples=200,
                           conf_level=.95, seed=seed)
    return result, variance


def output_array(path, n, oxygen):
    """Require every ordered sample/O2 pair, retaining both outputs together."""
    values = np.full((n, len(oxygen), len(OUTPUTS)), np.nan)
    seen = np.zeros((n, len(oxygen)), dtype=bool)
    lookup = {round(float(o), 10): j for j, o in enumerate(oxygen)}
    qc = dict(n_points=0, n_recovered=0, minimum_spectral_gap=None, maximum_residual=0.)
    for row in read(path):
        i = int(row['sample_id']) - 1
        j = lookup.get(round(float(row['O2_pct']), 10))
        if not 0 <= i < n or j is None or seen[i, j]:
            raise ValueError('Duplicate/unexpected sample row')
        if row['status'] != 'ok' or row['eigenvector_nonnegative'] != 'TRUE':
            raise ValueError('Invalid model evaluation; rows must not be deleted')
        values[i, j] = [float(row[k]) for k in OUTPUTS]
        if not np.isfinite(values[i, j]).all():
            raise ValueError('Nonfinite output')
        seen[i, j] = True
        gap = float(row['spectral_gap'])
        if np.isfinite(gap):
            qc['minimum_spectral_gap'] = gap if qc['minimum_spectral_gap'] is None else min(gap, qc['minimum_spectral_gap'])
        qc['maximum_residual'] = max(qc['maximum_residual'], float(row.get('eigen_residual', 0)))
        qc['n_recovered'] += row.get('solver', '') == 'decimal50_100'
        qc['n_points'] += 1
    if not seen.all():
        raise ValueError('Incomplete ordered model output')
    return values, qc


def assemble(block_arrays):
    if len(block_arrays) != STEP or any(a.shape != block_arrays[0].shape for a in block_arrays):
        raise ValueError('Require all A/B/AB/BA blocks with identical axes')
    # SALib interleaves all 30 matrix rows for each base Sobol point.
    return np.stack(block_arrays, axis=1).reshape((-1,) + block_arrays[0].shape[1:])


def case_directory(root, case):
    return Path(root) / case['scope'] / 'pilot' / case['case_id']


def prefix(natural, n):
    return natural[:n * STEP]


def complete(path):
    receipt = Path(str(path) + '.receipt.json')
    return Path(path).is_file() and receipt.is_file() and load_json(receipt)['outputs_sha256'] == sha(path)


def configure(plan):
    """Use the submitted resolution ladder rather than a plotting default."""
    global RESOLUTIONS, OXYGEN
    ns = list(map(int, plan['resolutions']))
    if len(ns) < 2 or ns != sorted(set(ns)) or any(n < 2 or n & (n-1) for n in ns):
        raise ValueError('Require at least two increasing power-of-two resolutions')
    RESOLUTIONS = ns
    OXYGEN = list(map(float, plan['oxygen_pct']))
