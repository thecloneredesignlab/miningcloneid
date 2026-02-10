#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import yaml

ROOT = Path(__file__).resolve().parents[1]
ODE_CFG = ROOT / "config" / "ode_params.yml"
MCTS_CFG = ROOT / "config" / "mcts_params.yml"
OUT_DIR = ROOT / "generated"
OUT_DIR.mkdir(exist_ok=True)

def _require(d: Dict[str, Any], key: str) -> Any:
    if key not in d:
        raise KeyError(f"Missing required key: {key}")
    return d[key]

def _latex_table_header(caption: str, label: str, footnotesize: bool = False) -> str:
    size = "\\footnotesize\n" if footnotesize else ""
    return (
        "\\begin{table}[H]\n"
        "\\centering\n"
        f"{size}"
        f"\\caption{{{caption}}}\n"
        f"\\label{{{label}}}\n"
    )

def render_ode_table(spec: Dict[str, Any]) -> str:
    caption = _require(spec, "caption")
    label = _require(spec, "label")
    columns = _require(spec, "columns")
    groups = _require(spec, "groups")

    tex = _latex_table_header(caption=caption, label=label, footnotesize=False)
    tex += "\\begin{tabularx}{\\textwidth}{l l X l}\n"
    tex += "\\toprule\n"
    tex += " & ".join([c["header"] for c in columns]) + " \\\\\n"
    tex += "\\midrule\n"

    for g in groups:
        title = _require(g, "title")
        tex += f"\\multicolumn{{4}}{{l}}{{{title}}} \\\\\n"
        for row in _require(g, "rows"):
            tex += (
                f"{row['symbol']} & {row['code_variable']} & "
                f"{row['description']} & {row['value_source']} \\\\\n"
            )
        tex += "\\midrule\n"

    tex += "\\bottomrule\n"
    tex += "\\end{tabularx}\n"
    tex += "\\end{table}\n"
    return tex

def render_mcts_table(spec: Dict[str, Any]) -> str:
    caption = _require(spec, "caption")
    label = _require(spec, "label")
    footnotesize = bool(spec.get("footnotesize", False))
    columns = _require(spec, "columns")
    rows = _require(spec, "rows")

    tex = _latex_table_header(caption=caption, label=label, footnotesize=footnotesize)
    tex += "\\begin{tabular}{l p{7cm} l}\n"
    tex += "\\toprule\n"
    tex += " & ".join([c["header"] for c in columns]) + " \\\\\n"
    tex += "\\midrule\n"

    for r in rows:
        tex += f"{r['symbol']} & {r['description']} & {r['value_unit']} \\\\\n"

    tex += "\\bottomrule\n"
    tex += "\\end{tabular}\n"
    tex += "\\end{table}\n"
    return tex

def load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def main() -> None:
    ode_cfg = load_yaml(ODE_CFG)
    mcts_cfg = load_yaml(MCTS_CFG)

    ode_tex = render_ode_table(ode_cfg)
    mcts_tex = render_mcts_table(mcts_cfg)

    (OUT_DIR / "table_model_parameters.tex").write_text(ode_tex, encoding="utf-8")
    (OUT_DIR / "table_mcts_parameters.tex").write_text(mcts_tex, encoding="utf-8")

    print("Wrote:")
    print(f" - {OUT_DIR / 'table_model_parameters.tex'}")
    print(f" - {OUT_DIR / 'table_mcts_parameters.tex'}")

if __name__ == "__main__":
    main()

