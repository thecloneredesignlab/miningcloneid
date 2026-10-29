#!/usr/bin/env python3
"""Global-only continuation of certified N32768 outputs; no old model calls.

Host-side helpers use Python 3.6+. Numerical work runs in the original SIF.
The original frozen engine and SALib estimator are imported read-only.
"""
import argparse
import contextlib
import datetime
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

BASE_N = 32768
CHUNK = 4096
NS = (65536, 131072)
REPEATS = 5
BLOCK_COUNT = 30
OXYGEN_COUNT = 201
FOCAL = ('p_mis_base', 'mu_hp', 'gamma_mu')
OUTPUTS = ('dominant_mean_ploidy', 'dominant_growth_rate')


def now():
    return datetime.datetime.utcnow().isoformat() + 'Z'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for part in iter(lambda: handle.read(8 * 1024**2), b''):
            h.update(part)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def save(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp.' + str(os.getpid()))
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    os.replace(str(tmp), str(path))


@contextlib.contextmanager
def locked(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        yield


def model_count(n):
    if n not in NS:
        raise ValueError('Only the two approved global resolutions are allowed')
    return (n // 2 // CHUNK) * REPEATS * BLOCK_COUNT


def decode(n, task):
    if not 1 <= task <= model_count(n):
        raise ValueError('Model task outside approved tail')
    segment, within = divmod(task - 1, REPEATS * BLOCK_COUNT)
    repeat, block = divmod(within, BLOCK_COUNT)
    return repeat + 1, block, n // 2 + segment * CHUNK


def analysis_task(task):
    if not 1 <= task <= 130:
        raise ValueError('Analysis task outside five repeats / 201 oxygen points')
    repeat, chunk = divmod(task - 1, 26)
    return repeat + 1, list(range(chunk * 8, min(201, chunk * 8 + 8)))


def unit_path(root, repeat, bi, begin):
    return Path(root) / 'units' / ('R%02d_B%02d_from%06d.npz' % (repeat, bi, begin))


def model_path(root, repeat, bi, begin):
    return Path(root) / 'increments' / ('base%06d_%06d' % (begin + 1, begin + CHUNK)) / ('R%02d' % repeat) / c.BLOCKS[bi] / 'outputs.npz'


def model_identity(root, repeat, bi, begin, unit_sha256, bounds_sha256):
    return hashlib.sha256(json.dumps(dict(plan=sha(Path(root) / 'plan.json'), repeat=repeat, block=bi, begin=begin,
        unit=unit_sha256, bounds=bounds_sha256), sort_keys=True).encode()).hexdigest()


def shard_path(root, n, repeat, ois):
    return Path(root) / ('N%d' % n) / 'shards' / 'global' / ('R%02d_O%03d.npz' % (repeat, ois[0]))


def runtime(root):
    global np, r, c, f
    import numpy as np
    root = Path(root)
    p = load(root / 'plan.json')
    if p['scope'] != 'global' or p['resolutions'] != list(NS) or p['n_repeats'] != 5:
        raise ValueError('Unapproved execution scope')
    for name, expected in p['controller_sha256'].items():
        if sha(root / 'slurm/code_snapshot' / name) != expected:
            raise ValueError('Frozen continuation code changed: ' + name)
    if sha(__file__) != p['controller_sha256']['global_convergence.py']:
        raise ValueError('Use the frozen continuation controller')
    oldref = Path(p['source_refinement_root'])
    if sha(oldref / 'plan.json') != p['source_refinement_plan_sha256']:
        raise ValueError('Source refinement identity changed')
    spec = importlib.util.spec_from_file_location('certified_source_refinement', oldref / 'slurm/code_snapshot/refinement_run.py')
    r = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(r)
    source, original, old = r.runtime(oldref)
    c, f = r.c, r.f
    if len(c.ACTIVE) != 14 or len(c.OXYGEN) != 201 or not np.allclose(c.OXYGEN, np.linspace(0, 5, 201), rtol=0, atol=1e-12):
        raise ValueError('Original model axes differ')
    if original['cases'][0]['scope'] != 'global':
        raise ValueError('Source first case is not global')
    if tuple(c.OUTPUTS) != OUTPUTS:
        raise ValueError('Original output definitions differ')
    for path, expected in p['source_files'].items():
        if sha(path) != expected:
            raise ValueError('Certified global source changed: ' + path)
    if (root / 'preparation.json').exists() and load(root / 'preparation.json')['plan_sha256'] != sha(root / 'plan.json'):
        raise ValueError('Preparation identity changed')
    return p, source, original, old


def file_record(path):
    st = Path(path).stat()
    return dict(path=str(path), sha256=sha(path), bytes=st.st_size, mtime_ns=st.st_mtime_ns)


def checked_record(path, records):
    rec = records[str(path)]
    st = Path(path).stat()
    if (st.st_size, st.st_mtime_ns) != (rec['bytes'], rec['mtime_ns']):
        raise ValueError('Audited source file changed: ' + str(path))


def prepare(root):
    p, source, original, old = runtime(root)
    root = Path(root)
    if (root / 'preparation.json').exists():
        return
    if shutil.disk_usage(root).free < 128 * 1024**3:
        raise ValueError('Require 128 GiB free for global tails, recovery records and reserve')
    subprocess.run([sys.executable, str(root / 'slurm/code_snapshot/test_global_convergence.py')], check=True)
    if sha(original['sif']['path']) != c.SIF_SHA256:
        raise ValueError('Original SIF digest changed')
    for name, expected in original['code_sha256'].items():
        if sha(Path(original['code_snapshot_root']) / name) != expected:
            raise ValueError('Frozen original helper changed: ' + name)
    case = original['cases'][0]
    directory = f.case_dir(old, case)
    bounds = c.read(directory / 'parameter_bounds.tsv')
    if sha(directory / 'parameter_bounds.tsv') != case['bounds_sha256']:
        raise ValueError('Global bounds changed')
    for name, field in [('best_params.tsv', 'best_params_sha256'), ('fit_config.rds', 'fit_config_sha256')]:
        if sha(Path(original['fit_root']) / case['fit_seed'] / name) != case[field]:
            raise ValueError('Original global model input changed')
    unit_ranges = [dict(parameter=q, encoded_lower=0., encoded_upper=1., natural_lower=0., natural_upper=1., transform='identity') for q in c.ACTIVE]
    units, outputs, audit_rows = {}, {}, []
    sourceprep = load(Path(p['source_refinement_root']) / 'preparation.json')
    expected_base = {x['path']: x['sha256'] for x in load(directory / 'completion.json')['output_receipts']}
    for repeat in range(1, 6):
        high = c.samples(unit_ranges, NS[-1], 2026101300 + repeat)
        low = c.samples(unit_ranges, NS[0], 2026101300 + repeat)
        if not np.array_equal(low, high[:NS[0] * 30]):
            raise ValueError('N65536/N131072 prefixes differ')
        del low
        prefix = old / 'inputs' / ('unit_R%02d.npz' % repeat)
        if sha(prefix) != original['unit_sha256'][str(repeat)]:
            raise ValueError('N4096 unit prefix changed')
        with np.load(prefix) as z:
            if not np.array_equal(z['unit'], high[:4096 * 30]):
                raise ValueError('Original scrambled prefix differs')
        for begin in range(0, BASE_N, CHUNK):
            for bi in range(30):
                u = high[begin * 30:(begin + CHUNK) * 30][bi::30]
                if begin == 0:
                    path = directory / 'runs' / ('R%02d' % repeat) / 'blocks' / c.BLOCKS[bi] / 'outputs.npz'
                    identity = f.identity(original, case, repeat, bi, f.scale_unit(u, bounds))
                else:
                    oldunit = r.unit_path(p['source_refinement_root'], repeat, bi, begin)
                    if sha(oldunit) != sourceprep['unit_sha256'][str(oldunit.relative_to(p['source_refinement_root']))]:
                        raise ValueError('Source increment sample file changed')
                    with np.load(oldunit) as z:
                        if not np.array_equal(z['unit'], u):
                            raise ValueError('N32768 old increment differs from new prefix')
                    path = r.block_path(source, case, repeat, bi, begin)
                    identity = hashlib.sha256(json.dumps(dict(plan=p['source_refinement_plan_sha256'], case=case, repeat=repeat,
                        block=bi, begin=begin, oxygen_indices=list(range(201)), unit=sha(oldunit)), sort_keys=True).encode()).hexdigest()
                if not f.complete(path, identity):
                    raise ValueError('Required old model block missing: ' + str(path))
                rec = file_record(path)
                if begin == 0 and expected_base[str(path.relative_to(old))] != rec['sha256']:
                    raise ValueError('Original model completion manifest differs')
                outputs[str(path)] = rec
                with np.load(path) as z:
                    if z['outputs'].shape != (4096, 201, 2) or not np.array_equal(z['oxygen'], c.OXYGEN):
                        raise ValueError('Original output axes differ')
                    flags = z['solver_flags']
                    if not np.isfinite(z['outputs']).all() or not np.isin(flags, [1, 2]).all():
                        raise ValueError('Invalid certified original output')
                    audit_rows.append(dict(repeat=repeat, block=c.BLOCKS[bi], sample_begin=begin, n_recovered=int((flags == 2).sum()),
                        minimum_gap=float(np.min(z['minimum_gap'])), maximum_residual=float(np.max(z['maximum_residual']))))
        for begin in range(BASE_N, NS[-1], CHUNK):
            for bi in range(30):
                path = unit_path(root, repeat, bi, begin)
                data = high[begin * 30:(begin + CHUNK) * 30][bi::30]
                if path.exists():
                    with np.load(path) as z:
                        if not np.array_equal(data, z['unit']):
                            raise ValueError('Prepared new unit file differs')
                else:
                    f.atomic_npz(path, unit=data)
                units[str(path.relative_to(root))] = sha(path)
        print('Certified N32768 prefix and new-only tails, repeat', repeat, flush=True)
        del high
    # Freeze N32768 lower-resolution statistical shards as well as model blocks.
    for task in range(1, 131):
        repeat, ois = analysis_task(task)
        path = Path(p['source_global_root']) / 'N32768/shards/global' / ('R%02d_O%03d.npz' % (repeat, ois[0]))
        identity = p['source_refinement_plan_sha256'] + '/global/32768/%d/' % repeat + str(ois)
        if not f.complete(path, identity):
            raise ValueError('N32768 statistics incomplete')
        with np.load(path) as z:
            if not z['done'].all() or not np.array_equal(z['oxygen_indices'], ois):
                raise ValueError('Source statistical axes differ')
        outputs[str(path)] = file_record(path)
    r.write_table(root / 'source_solver_audit.tsv', audit_rows, list(audit_rows[0]))
    save(root / 'source_audit.json', dict(status='complete', files=outputs, n_model_blocks=1200, n_statistical_shards=130, completed_utc=now()))
    save(root / 'preparation.json', dict(status='complete', plan_sha256=sha(root / 'plan.json'), unit_sha256=units,
        source_audit_sha256=sha(root / 'source_audit.json'), n_verified_prefixes=5, new_sample_min=32769,
        model_evaluations_repeated=0, completed_utc=now(), free_bytes=shutil.disk_usage(root).free))


def evaluate(root, n, task, smoke=False):
    p, source, original, old = runtime(root)
    root = Path(root)
    prep = load(root / 'preparation.json')
    repeat, bi, begin = decode(n, task)
    case = original['cases'][0]
    unitfile = unit_path(root, repeat, bi, begin)
    if sha(unitfile) != prep['unit_sha256'][str(unitfile.relative_to(root))]:
        raise ValueError('New unit design changed')
    bounds = c.read(f.case_dir(old, case) / 'parameter_bounds.tsv')
    with np.load(unitfile) as z:
        x = f.scale_unit(z['unit'], bounds)
    path = model_path(root, repeat, bi, begin)
    key = model_identity(root, repeat, bi, begin, sha(unitfile), case['bounds_sha256'])
    with locked(path.parent / 'worker.lock'):
        if f.complete(path, key):
            return
        checkpoint = path.parent / 'checkpoint.npz'
        oxygen = np.asarray(c.OXYGEN)
        values, flags, done, gaps, residuals, elapsed = f.checkpoint(checkpoint, key, CHUNK, oxygen)
        scratch = Path(os.environ['SOBOL_SCRATCH']) / 'evaluation'
        scratch.mkdir(parents=True, exist_ok=True)
        c.table(scratch / 'samples.tsv.gz', [dict(sample_id=i + 1, **dict(zip(c.ACTIVE, map(float, row)))) for i, row in enumerate(x)])
        while not done.all():
            if shutil.disk_usage(root).free < 64 * 1024**3:
                raise ValueError('Storage reserve below 64 GiB; resume saved checkpoint')
            positions = np.flatnonzero(~done)[:8]
            current = oxygen[positions].tolist()
            meta = dict(method='Sobol_second_order', scope='global', case_id='global', fit_seed=case['fit_seed'], n_samples=CHUNK,
                parameters=list(c.ACTIVE), oxygen_pct=current, samples_sha256=sha(scratch / 'samples.tsv.gz'),
                best_params_sha256=case['best_params_sha256'], fit_config_sha256=case['fit_config_sha256'],
                validate_reference=False, sample_id_offset=begin, source_prefix_reused=BASE_N)
            save(scratch / 'metadata.json', meta)
            subprocess.run(['Rscript', '--vanilla', str(f.HERE / 'evaluate.R'), str(scratch), original['fit_root'],
                str(old / 'inputs/fixed_o2_dominant_ploidy_201grid.tsv')], check=True)
            v, s, g, residual = f.read_model(scratch / 'outputs.tsv.gz', CHUNK, current)
            receipt = load(str(scratch / 'outputs.tsv.gz') + '.receipt.json')
            if receipt['samples_sha256'] != meta['samples_sha256'] or receipt['metadata_sha256'] != sha(scratch / 'metadata.json') or receipt['outputs_sha256'] != sha(scratch / 'outputs.tsv.gz'):
                raise ValueError('Model output provenance differs')
            values[:, positions], flags[:, positions] = v, s
            gaps[positions], residuals[positions], done[positions] = g, residual, True
            elapsed += receipt['elapsed_seconds']
            if (scratch / 'recovery').exists():
                shutil.copytree(scratch / 'recovery', path.parent / 'recovery', dirs_exist_ok=True)
                shutil.rmtree(scratch / 'recovery')
            f.atomic_npz(checkpoint, identity=np.array(key), oxygen=oxygen, outputs=values, solver_flags=flags,
                done=done, minimum_gap=gaps, maximum_residual=residuals, elapsed_seconds=np.array(elapsed))
            print(n, repeat, c.BLOCKS[bi], begin, 'oxygen', int(done.sum()), '/201', flush=True)
            if smoke:
                save(root / 'smoke.json', dict(status='complete', checkpoint=str(checkpoint), checkpoint_sha256=sha(checkpoint),
                    sample_begin=begin, certified_new_model_points=int(done.sum()) * CHUNK, original_model_points_evaluated=0,
                    production_resumes_same_checkpoint=True, elapsed_seconds=elapsed, completed_utc=now()))
                return
        f.atomic_npz(path, identity=np.array(key), oxygen=oxygen, oxygen_indices=np.arange(201), outputs=values,
            solver_flags=flags, minimum_gap=gaps, maximum_residual=residuals)
        save(str(path) + '.receipt.json', dict(status='complete', identity=key, sha256=sha(path), sample_begin=begin,
            sample_end=begin + CHUNK, elapsed_seconds=elapsed, n_recovered=int((flags == 2).sum()),
            n_points=CHUNK * 201, model_evaluations_repeated=0, completed_utc=now(), job_id=os.environ.get('SLURM_JOB_ID')))
        checkpoint.unlink()


def assemble(root, p, source, original, old, n, repeat, ois):
    audit = load(Path(root) / 'source_audit.json')
    prep = load(Path(root) / 'preparation.json')
    if sha(Path(root) / 'source_audit.json') != prep['source_audit_sha256']:
        raise ValueError('Source audit changed')
    case = original['cases'][0]
    parts = []
    for begin in range(0, n, CHUNK):
        blocks = []
        for bi, block in enumerate(c.BLOCKS):
            if begin == 0:
                path = f.case_dir(old, case) / 'runs' / ('R%02d' % repeat) / 'blocks' / block / 'outputs.npz'
            elif begin < BASE_N:
                path = r.block_path(source, case, repeat, bi, begin)
            else:
                path = model_path(root, repeat, bi, begin)
            if begin < BASE_N:
                checked_record(path, audit['files'])
            else:
                unitsha = prep['unit_sha256'][str(unit_path(root, repeat, bi, begin).relative_to(root))]
                if not f.complete(path, model_identity(root, repeat, bi, begin, unitsha, case['bounds_sha256'])):
                    raise ValueError('New tail block incomplete: ' + str(path))
            with np.load(path) as z:
                if z['outputs'].shape != (CHUNK, 201, 2) or not np.array_equal(z['oxygen'], c.OXYGEN):
                    raise ValueError('Model axes differ')
                blocks.append(z['outputs'][:, ois, :])
        parts.append(c.assemble(blocks))
    ordered = np.concatenate(parts, axis=0)
    if ordered.shape != (n * 30, len(ois), 2) or not np.isfinite(ordered).all():
        raise ValueError('Incomplete ordered Sobol design')
    return ordered


def analyze(root, n, task):
    if n not in NS:
        raise ValueError('Unapproved resolution')
    p, source, original, old = runtime(root)
    root = Path(root)
    repeat, ois = analysis_task(task)
    path = shard_path(root, n, repeat, ois)
    key = sha(root / 'plan.json') + '/%d/%d/%s' % (n, repeat, ois)
    with locked(str(path) + '.lock'):
        if f.complete(path, key):
            return
        if path.exists():
            with np.load(path) as z:
                if str(z['identity']) != key or not np.array_equal(z['oxygen_indices'], ois):
                    raise ValueError('Statistical checkpoint identity changed')
                values, widths, variance, done = [z[k] for k in ('values', 'widths', 'variance', 'done')]
        else:
            shape = (2, len(ois), 2, len(f.METRICS))
            values, widths = np.full(shape, np.nan), np.full(shape, np.nan)
            variance, done = np.full(shape[:-1], np.nan), np.zeros((len(ois), 2), bool)
            lower = shard_path(p['source_global_root'] if n == NS[0] else root, n // 2, repeat, ois)
            if n == NS[0]:
                checked_record(lower, load(root / 'source_audit.json')['files'])
            elif not f.complete(lower, sha(root / 'plan.json') + '/%d/%d/%s' % (n // 2, repeat, ois)):
                raise ValueError('Lower resolution statistical shard incomplete')
            with np.load(lower) as z:
                if not z['done'].all() or not np.array_equal(z['oxygen_indices'], ois):
                    raise ValueError('Lower resolution axes differ')
                values[0], widths[0], variance[0] = z['values'][-1], z['widths'][-1], z['variance'][-1]
        ordered = assemble(root, p, source, original, old, n, repeat, ois)
        bounds = c.read(f.case_dir(old, original['cases'][0]) / 'parameter_bounds.tsv')
        for j, oi in enumerate(ois):
            for k in range(2):
                if done[j, k]:
                    continue
                si, var = c.analyze(bounds, ordered[:, j, k], seed=2026101300 + repeat + 10000 + n + oi * 10 + k)
                values[1, j, k], widths[1, j, k] = f.metric_arrays(si)
                variance[1, j, k], done[j, k] = var, True
                f.atomic_npz(path, identity=np.array(key), oxygen_indices=np.array(ois), values=values, widths=widths, variance=variance, done=done)
                print(n, repeat, oi, k, 'index cell certified', flush=True)
        save(str(path) + '.receipt.json', dict(status='complete', identity=key, sha256=sha(path), completed_utc=now(), model_evaluations_repeated=0))


def report(root, n):
    if n not in NS:
        raise ValueError('Unapproved resolution')
    p, source, original, old = runtime(root)
    root = Path(root)
    dest = root / ('N%d' % n)
    if (dest / 'completion.json').exists():
        return
    vv, ww, va = [], [], []
    for repeat in range(1, 6):
        parts = []
        for oi in range(0, 201, 8):
            path = shard_path(root, n, repeat, list(range(oi, min(201, oi + 8))))
            key = sha(root / 'plan.json') + '/%d/%d/%s' % (n, repeat, list(range(oi, min(201, oi + 8))))
            if not f.complete(path, key):
                raise ValueError('Statistical shard missing')
            with np.load(path) as z:
                if not z['done'].all():
                    raise ValueError('Incomplete statistical shard')
                parts.append([z[k] for k in ('values', 'widths', 'variance')])
        vv.append(np.concatenate([a[0] for a in parts], axis=1))
        ww.append(np.concatenate([a[1] for a in parts], axis=1))
        va.append(np.concatenate([a[2] for a in parts], axis=1))
    values, widths, variance = np.stack(vv, axis=1), np.stack(ww, axis=1), np.stack(va, axis=1)
    mean, spread, delta, status = f.diagnostics(values)
    f.atomic_npz(dest / 'indices_repeats.npz', values=values, bootstrap_half_width=widths, variance=variance, oxygen=np.asarray(c.OXYGEN), resolutions=np.array([n // 2, n]))
    f.atomic_npz(dest / 'indices_summary.npz', mean=mean, repeat_range=spread, resolution_delta=delta, status=status, oxygen=np.asarray(c.OXYGEN))
    rows = []
    for oi, oxygen in enumerate(c.OXYGEN):
        for k, output in enumerate(OUTPUTS):
            for m, (idx, pi, pj) in enumerate(f.METRICS):
                rows.append(dict(scope='global', N=n, O2_pct=oxygen, output=output, index=idx, parameter_i=pi, parameter_j=pj,
                    value_mean=mean[oi, k, m], repeat_range=spread[oi, k, m], resolution_delta=delta[oi, k, m],
                    resolution_low=n // 2, resolution_high=n, bootstrap_half_width_median=float(np.median(widths[1, :, oi, k, m])),
                    output_variance_mean=float(np.mean(variance[1, :, oi, k])), convergence_status=['undefined', 'passed', 'not_passed'][status[oi, k, m]]))
    r.write_table(dest / 'indices_summary.tsv.gz', rows, list(rows[0]))
    failed = [row for row in rows if row['convergence_status'] != 'passed']
    r.write_table(dest / 'not_passed.tsv.gz', failed, list(rows[0]))
    from global_convergence_plot import render
    render(root, n, p, c, f, r, values, mean, status)
    save(dest / 'completion.json', dict(status='complete', scope='global', N=n, n_repeats=5, n_oxygen=201,
        n_indices=len(rows), summary_sha256=sha(dest / 'indices_summary.tsv.gz'), convergence_policy='diagnostic_only', completed_utc=now()))
    save(root / 'status.json', dict(status='complete' if n == NS[-1] else 'N65536_complete', completed_N=n,
        resolutions=list(NS), scope='global', model_evaluations_repeated=0, completed_utc=now()))


def main():
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument('mode', choices=['prepare', 'smoke', 'model', 'analyze', 'report'])
    a.add_argument('--root', required=True)
    a.add_argument('--n', type=int, default=65536)
    a.add_argument('--task', type=int, default=1)
    args = a.parse_args()
    try:
        if args.mode == 'prepare': prepare(args.root)
        elif args.mode == 'smoke': evaluate(args.root, 65536, 1, smoke=True)
        elif args.mode == 'model': evaluate(args.root, args.n, args.task)
        elif args.mode == 'analyze': analyze(args.root, args.n, args.task)
        else: report(args.root, args.n)
    except Exception as exc:
        save(Path(args.root) / 'failures' / ('%s_N%d_T%d.json' % (args.mode, args.n, args.task)),
            dict(error=str(exc), utc=now(), job_id=os.environ.get('SLURM_JOB_ID')))
        raise


if __name__ == '__main__':
    main()
