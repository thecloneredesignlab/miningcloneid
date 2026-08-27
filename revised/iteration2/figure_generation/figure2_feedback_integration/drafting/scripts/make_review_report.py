#!/usr/bin/env python3
from __future__ import annotations

import base64
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / "initial_subpanels" / "prior_figure2" / "assembled_fig2.png"
REVISED = ROOT / "final_figures" / "recommended" / "assembled_fig2.png"
OUT = ROOT / "review_report.html"


def data_url(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def main() -> None:
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Figure 2 feedback integration</title>
<style>
body {{ font-family: Arial, Helvetica, sans-serif; color: #222; max-width: 1200px;
  margin: 32px auto; padding: 0 24px; line-height: 1.45; }}
h1, h2 {{ margin-bottom: .35rem; }}
.status {{ border-left: 5px solid #009E73; background: #edf7f3;
  padding: 12px 16px; margin-bottom: 20px; }}
.grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; align-items: start; }}
.card {{ border: 1px solid #ddd; padding: 14px; background: white; }}
img {{ width: 100%; height: auto; display: block; }}
code {{ background: #f4f4f4; padding: 1px 4px; }}
@media (max-width: 900px) {{ .grid {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<h1>Figure 2 feedback integration</h1>
<div class="status"><strong>Recommended revision implemented.</strong> The
four-panel model overview is replaced by a three-panel vector figure. No fit,
simulation, or raster-derived panel was used.</div>
<div class="grid">
  <section class="card"><h2>Prior figure</h2><img alt="Prior four-panel Figure 2"
    src="{data_url(PRIOR)}"></section>
  <section class="card"><h2>Revised figure</h2><img alt="Revised three-panel Figure 2"
    src="{data_url(REVISED)}"></section>
</div>
<h2>Directive coverage</h2>
<ul>
  <li>Panel A now makes death hazard, rather than oxygen itself, the direct
  mediator of the effective chromosome-missegregation probability.</li>
  <li>The feedback returns adapted ploidy composition to population-average
  death even at fixed O2.</li>
  <li>WGD remains a separate constant-probability-per-division branch.</li>
  <li>Panel B is restricted to direct chromosome-burden fitness effects.</li>
  <li>Panel C compares correlated N-m and N+m daughters and isolates
  ploidy-dependent survival buffering.</li>
  <li>The prior fixed-O2 operator panel is deferred to the fixed-O2 analysis.</li>
</ul>
<h2>Visual QC</h2>
<p>The final 7.1-inch figure was inspected at full resolution and on the
rendered manuscript page. Labels and arrows do not overlap, no content is
clipped, panel order matches the legend, and every element is regenerated from
the package-local R source.</p>
</body>
</html>
"""
    OUT.write_text(html, encoding="utf-8")


if __name__ == "__main__":
    main()
