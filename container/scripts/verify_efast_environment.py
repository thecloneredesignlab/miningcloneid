#!/usr/bin/env python3
"""Check the pinned SALib runtime and its eFAST sampling/analysis path."""

from __future__ import annotations

import importlib.metadata
import sys

import numpy as np
from SALib.analyze import fast
from SALib.sample import fast_sampler
from SALib.test_functions import Ishigami


def main() -> int:
    assert sys.version_info[:2] == (3, 9), sys.version
    assert importlib.metadata.version("SALib") == "1.5.2"

    problem = {
        "num_vars": 3,
        "names": ["x1", "x2", "x3"],
        "bounds": [[-np.pi, np.pi]] * 3,
    }
    n = 1025
    samples = fast_sampler.sample(problem, n, M=4, seed=5826)
    assert samples.shape == (n * problem["num_vars"], problem["num_vars"])
    assert np.all(np.isfinite(samples))
    assert np.all(samples >= -np.pi) and np.all(samples <= np.pi)

    outputs = Ishigami.evaluate(samples)
    indices = fast.analyze(problem, outputs, M=4, num_resamples=20, seed=5826)
    expected = {
        "S1": np.array([0.314, 0.442, 0.0]),
        "ST": np.array([0.558, 0.442, 0.244]),
    }
    for name, reference in expected.items():
        actual = np.asarray(indices[name], dtype=float)
        assert actual.shape == reference.shape
        assert np.all(np.isfinite(actual))
        assert np.all(np.abs(actual - reference) < 0.12), (name, actual)
        print(f"{name}\t" + "\t".join(f"{value:.6f}" for value in actual))
    print("eFAST environment verification: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
