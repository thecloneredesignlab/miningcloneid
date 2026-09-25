#!/usr/bin/env python3
"""Deprecated compatibility wrapper for report/fit_results collection."""

from pathlib import Path
import runpy


SCRIPT = Path(__file__).resolve().parents[2] / "report" / "fit_results" / "collect_best_separate_fit_reports.py"
runpy.run_path(str(SCRIPT), run_name="__main__")
