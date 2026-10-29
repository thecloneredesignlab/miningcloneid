#!/usr/bin/env python3
"""Contracts for global-only submission, exact tails, retries and estimator parity."""
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import global_convergence as g
import submit_global_convergence as s


class PolicyTests(unittest.TestCase):
    def test_all_new_sample_blocks_are_covered_once(self):
        seen = set()
        for n, count in [(65536, 1200), (131072, 2400)]:
            self.assertEqual(g.model_count(n), count)
            for task in range(1, count + 1):
                repeat, block, begin = g.decode(n, task)
                self.assertNotIn((repeat, block, begin), seen)
                seen.add((repeat, block, begin))
                self.assertGreaterEqual(begin, 32768)
                self.assertLess(begin, n)
        for rep in range(1, 6):
            for block in range(30):
                self.assertEqual(sorted(x[2] for x in seen if x[:2] == (rep, block)), list(range(32768, 131072, 4096)))
        self.assertEqual(len(seen), 3600)
        for n, task in [(32768, 1), (65536, 0), (131072, 2401)]:
            with self.assertRaises(ValueError): g.decode(n, task)

    def test_statistics_cover_five_complete_oxygen_grids(self):
        seen = set()
        for task in range(1, 131):
            rep, ois = g.analysis_task(task)
            for oi in ois:
                self.assertNotIn((rep, oi), seen)
                seen.add((rep, oi))
        self.assertEqual(seen, {(r, o) for r in range(1, 6) for o in range(201)})

    def test_submission_has_only_two_global_levels_and_no_diagnostic_gate(self):
        rows = s.stages()
        self.assertEqual(len(rows), 8)
        self.assertEqual(sum(x[3] for x in rows if x[1] == 'model'), 3600)
        self.assertEqual(sum(x[3] for x in rows if x[1] == 'analyze'), 260)
        self.assertEqual(next(x[-1] for x in rows if x[0] == 'model_131072'), ['report_65536'])
        self.assertFalse(any('neighborhood' in str(x) for x in rows))
        for row in rows:
            args = s.sbatch_args(Path('/tmp/global_only'), row, ['123'])
            self.assertIn('--qos=xxlarge', args)
            self.assertIn('--time=12:00:00', args)
            self.assertIn('--hold', args)
            self.assertFalse(any(x.startswith(('--nodelist', '--nodes=')) for x in args))
            self.assertFalse(any('%' in x for x in args if x.startswith('--array=')))

    def test_colored_and_plain_startup_errors_retry_but_payload_errors_do_not(self):
        retry = Path(__file__).with_name('retry_global_startup.sh')
        for colored, identity, expected_calls, expected_code in [(True, True, 2, 0), (False, True, 2, 0), (False, False, 1, 7)]:
            with tempfile.TemporaryDirectory() as tmp:
                tmp = Path(tmp)
                count, child = tmp / 'count', tmp / 'child.sh'
                count.write_text('0')
                message = ("printf '\\033[31mFATAL: \\033[0m user: unknown userid 107865\\n'" if colored else
                           "printf 'FATAL: user: unknown userid 107865\\n'") if identity else "printf 'ValueError: actual model failure\\n'"
                child.write_text('#!/usr/bin/env bash\nn=$(cat "$1")\nn=$((n+1))\necho "$n" > "$1"\nif (( n == 1 )); then\n' + message + '\nexit ' + ('127' if identity else '7') + '\nfi\nexit 0\n')
                cmd = ['bash', '-c', 'set -euo pipefail; scratch=$1; source "$2"; export SOBOL_STARTUP_DELAY_SECONDS=0; launch bash "$3" "$4"',
                       'test', str(tmp), str(retry), str(child), str(count)]
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                self.assertEqual(result.returncode, expected_code, result.stderr.decode())
                self.assertEqual(int(count.read_text()), expected_calls)


@unittest.skipUnless(importlib.util.find_spec('numpy') and importlib.util.find_spec('SALib'), 'Numerical tests require SALib 1.5.2 environment')
class NumericalTests(unittest.TestCase):
    def test_all_five_high_resolution_prefixes_and_tail_interleaving(self):
        import numpy as np
        from SALib.sample import sobol
        problem = dict(num_vars=14, names=['p%d' % i for i in range(14)], bounds=[[0., 1.]] * 14)
        for rep in range(1, 6):
            high = sobol.sample(problem, 131072, calc_second_order=True, scramble=True, seed=2026101300 + rep)
            for n in [4096, 32768, 65536]:
                low = sobol.sample(problem, n, calc_second_order=True, scramble=True, seed=2026101300 + rep)
                self.assertTrue(np.array_equal(low, high[:n * 30]))
                del low
            for begin in [32768, 65536, 126976]:
                tail = high[begin * 30:(begin + 4096) * 30]
                joined = np.stack([tail[b::30] for b in range(30)], axis=1).reshape((-1, 14))
                self.assertTrue(np.array_equal(tail, joined))
            del high

    def test_old_plus_new_outputs_equal_one_shot_estimator(self):
        import numpy as np
        from SALib.sample import sobol
        from SALib.analyze import sobol as analyzer
        problem = dict(num_vars=3, names=['a', 'b', 'c'], bounds=[[0., 1.]] * 3)
        x = sobol.sample(problem, 1024, calc_second_order=True, scramble=True, seed=42)
        model = lambda a: a[:, 0] + a[:, 1] * a[:, 2] + a[:, 0] * a[:, 1]
        prefix = model(x[:256 * 8])
        calls = []
        for lo, hi in [(256, 512), (512, 1024)]:
            calls.append((lo, hi))
            prefix = np.concatenate([prefix, model(x[lo * 8:hi * 8])])
            a = analyzer.analyze(problem, prefix, calc_second_order=True, seed=77, num_resamples=20)
            b = analyzer.analyze(problem, model(x[:hi * 8]), calc_second_order=True, seed=77, num_resamples=20)
            for key in ['S1', 'ST', 'S2', 'S1_conf', 'ST_conf', 'S2_conf']:
                self.assertTrue(np.allclose(a[key], b[key], equal_nan=True))
        self.assertEqual(calls, [(256, 512), (512, 1024)])

    def test_production_assembler_preserves_sample_block_and_oxygen_order(self):
        import numpy as np
        from types import SimpleNamespace
        blocks = ['B%d' % i for i in range(30)]
        def fake_load(path):
            text = str(path)
            if '/runs/' in text:
                begin, bi = 0, int(Path(path).parent.name[1:])
            elif '/old/' in text:
                begin, bi = 2, int(Path(path).parent.name[1:])
            else:
                begin, bi = int(Path(path).parents[2].name.split('_')[0][4:]) - 1, int(Path(path).parent.name[1:])
            output = np.empty((2, 201, 2))
            for row in range(2):
                output[row, :, 0] = (begin + row) * 30 + bi
                output[row, :, 1] = np.arange(201)
            class Z(dict):
                def __enter__(self): return self
                def __exit__(self, *args): pass
            return Z(outputs=output, oxygen=np.linspace(0, 5, 201))
        c = SimpleNamespace(BLOCKS=blocks, OXYGEN=np.linspace(0, 5, 201), assemble=lambda bs: np.stack(bs, axis=1).reshape((-1, len(ois), 2)))
        f = SimpleNamespace(case_dir=lambda old, case: Path('/source'), complete=lambda *args: True)
        r = SimpleNamespace(block_path=lambda *args: Path('/old') / blocks[args[-2]] / 'outputs.npz')
        ois = [0, 9, 200]
        unitsha = {str(g.unit_path(Path('/new'), 1, bi, begin).relative_to('/new')): 'ok' for begin in [4, 6] for bi in range(30)}
        with patch.object(g, 'CHUNK', 2), patch.object(g, 'BASE_N', 4), patch.object(g, 'np', np, create=True), patch.object(g, 'c', c, create=True), patch.object(g, 'f', f, create=True), patch.object(g, 'r', r, create=True), patch.object(g, 'load', lambda path: {'files': {}, 'source_audit_sha256': 'ok', 'unit_sha256': unitsha}), patch.object(g, 'sha', lambda path: 'ok'), patch.object(g, 'checked_record'), patch.object(np, 'load', fake_load):
            actual = g.assemble(Path('/new'), {}, {}, {'cases': [dict(bounds_sha256='ok')]}, Path('/source'), 8, 1, ois)
        self.assertTrue(np.array_equal(actual[:, 0, 0], np.arange(8 * 30)))
        self.assertTrue(np.array_equal(actual[0, :, 1], ois))

    def test_o2_sort_and_bootstrap_match_whole_curve_resampling(self):
        import numpy as np
        import global_convergence_plot as plot
        rng = np.random.default_rng(123)
        data = rng.random((5, 201, 2, 3))
        oxygen = np.linspace(0, 5, 201)
        metrics = [('S1', 'a', ''), ('S1', 'b', ''), ('ST', 'c', '')]
        rows = plot.classify(data, metrics, oxygen)
        draws = np.random.default_rng(5826).integers(0, 5, size=(5000, 5))
        lw, hw = plot.weights(oxygen, 0, 1), plot.weights(oxygen, 3, 5)
        for i, row in enumerate(rows):
            delta = (data[draws, :, 0, i].mean(axis=1) @ (lw - hw))
            self.assertAlmostEqual(row['bootstrap_delta_q025'], float(np.quantile(delta, .025)))
        classes = [dict(low_score=.2, high_score=.6), dict(low_score=.7, high_score=.1), dict(low_score=.3, high_score=.1), dict(low_score=.2, high_score=.8)]
        self.assertEqual(plot.row_order(range(4), classes), [3, 0, 2, 1])


if __name__ == '__main__':
    unittest.main(verbosity=2)
