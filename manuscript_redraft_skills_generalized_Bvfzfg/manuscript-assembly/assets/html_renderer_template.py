#!/usr/bin/env python3
"""Template renderer for a self-contained manuscript HTML review draft.

Copy this file into an assembly package, adapt `renderer_config.json` or the
defaults below, and run it from any location:

    python rebuild/manuscript/build_manuscript_html.py --assembly-root .

Expected package inputs, relative to the assembly root:

- source/manuscript_sections/*.md
- source/manuscript_sections/results/*.md
- source/integrated_figure_legends.md
- source/figure_set_manifest.csv
- assets/figures/*.png

The script intentionally performs mechanical rendering only. It validates
presence and consistency, embeds figures as base64 data URIs, and keeps audit
material out of journal-facing prose.
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import json
import re
from pathlib import Path
from typing import Any


DEFAULT_CONFIG: dict[str, Any] = {
    "title": "Manuscript Draft",
    "output": "draft/manuscript_draft.html",
    "manifest": "source/figure_set_manifest.csv",
    "legends": "source/integrated_figure_legends.md",
    "sections_dir": "source/manuscript_sections",
    "results_dir": "source/manuscript_sections/results",
    "figures_dir": "assets/figures",
    "legend_validation_report": "validation/draft_render_legend_validation_report.md",
    "section_order": ["abstract", "introduction", "results", "discussion", "methods"],
    "results_order": [],
    "figure_placements": {},
    "audit_sources": [],
    "embed_images": True,
}


def read_config(root: Path, explicit_config: Path | None) -> dict[str, Any]:
    config = dict(DEFAULT_CONFIG)
    config_path = explicit_config or root / "rebuild" / "manuscript" / "renderer_config.json"
    if config_path.exists():
        loaded = json.loads(config_path.read_text())
        if not isinstance(loaded, dict):
            raise ValueError(f"Renderer config must be a JSON object: {config_path}")
        config.update(loaded)
    return config


def e(text: object) -> str:
    return html.escape(str(text), quote=True)


def slug(text: str) -> str:
    value = re.sub(r"[^a-zA-Z0-9]+", "-", text.strip().lower()).strip("-")
    return value or "section"


def display_label(figure_id: str) -> str:
    if re.match(r"^S\d+", figure_id):
        suffix = figure_id.removeprefix("S")
        if suffix.endswith("_continued"):
            return f"Figure S{suffix.removesuffix('_continued')} continued"
        return f"Figure S{suffix}"
    if re.match(r"^F\d+", figure_id):
        return f"Figure {figure_id.removeprefix('F')}"
    if figure_id.startswith("Figure_S"):
        suffix = figure_id.removeprefix("Figure_S")
        if suffix.endswith("_continued"):
            return f"Figure S{suffix.removesuffix('_continued')} continued"
        return f"Figure S{suffix}"
    if figure_id.startswith("Figure_"):
        return f"Figure {figure_id.removeprefix('Figure_')}"
    return figure_id.replace("_", " ")


def label_to_figure_id(label: str) -> str:
    label = " ".join(label.strip().split())
    match = re.match(r"^Figure\s+(S?)(\d+)([A-Z]?)(?:\s+(continued))?$", label)
    if not match:
        raise ValueError(f"Could not parse figure label: {label}")
    is_supp, number, letter, continued = match.groups()
    figure_id = f"S{number}{letter}" if is_supp else f"F{number}"
    if continued:
        figure_id += "_continued"
    return figure_id


def figure_anchor(figure_id: str) -> str:
    return f"fig-{slug(figure_id)}"


def figure_sort_key(figure_id: str) -> tuple[int, int, str, int]:
    normalized = figure_id
    if normalized.startswith("Figure_"):
        normalized = normalized.removeprefix("Figure_")
    if normalized.startswith("S"):
        suffix = normalized.removeprefix("S")
        continued = 1 if suffix.endswith("_continued") else 0
        suffix = suffix.removesuffix("_continued")
        match = re.match(r"^(\d+)([A-Z]?)$", suffix)
        if match:
            return (1, int(match.group(1)), match.group(2), continued)
    if normalized.startswith("F"):
        number = normalized.removeprefix("F")
        if number.isdigit():
            return (0, int(number), "", 0)
    return (2, 10_000, normalized, 0)


def read_markdown(path: Path) -> tuple[dict[str, str], str]:
    text = path.read_text().strip()
    if not text:
        raise ValueError(f"Empty Markdown file: {path}")
    if not text.startswith("---"):
        return {}, text

    lines = text.splitlines()
    end = None
    for index, line in enumerate(lines[1:], start=1):
        if line.strip() == "---":
            end = index
            break
    if end is None:
        raise ValueError(f"Unclosed front matter in {path}")

    metadata: dict[str, str] = {}
    for line in lines[1:end]:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if ":" not in line:
            raise ValueError(f"Malformed front matter line in {path}: {line}")
        key, value = line.split(":", 1)
        metadata[key.strip()] = value.strip().strip('"')
    return metadata, "\n".join(lines[end + 1 :]).strip()


def extract_named_body(body: str, preferred_headings: list[str]) -> str:
    lines = body.splitlines()
    for heading in preferred_headings:
        pattern = re.compile(rf"^##\s+{re.escape(heading)}\s*$")
        start = None
        for index, line in enumerate(lines):
            if pattern.match(line):
                start = index + 1
                break
        if start is None:
            continue
        end = len(lines)
        for index in range(start, len(lines)):
            if lines[index].startswith("## "):
                end = index
                break
        text = "\n".join(lines[start:end]).strip()
        if text:
            return text
    return body.strip()


def markdown_inline(text: str, known_figures: set[str]) -> str:
    code_spans: list[str] = []

    def stash_code(match: re.Match[str]) -> str:
        code_spans.append(f"<code>{e(match.group(1))}</code>")
        return f"@@CODE{len(code_spans) - 1}@@"

    text = re.sub(r"`([^`]+)`", stash_code, text)
    escaped = e(text)
    escaped = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", escaped)
    escaped = link_figure_refs(escaped, known_figures)
    for index, code_html in enumerate(code_spans):
        escaped = escaped.replace(f"@@CODE{index}@@", code_html)
    return escaped


def markdown_blocks(text: str, known_figures: set[str]) -> str:
    chunks: list[str] = []
    blocks = [block.strip() for block in re.split(r"\n\s*\n", text.strip()) if block.strip()]
    for block in blocks:
        heading = re.match(r"^(#{1,4})\s+(.+)$", block)
        if heading:
            level = max(2, min(4, len(heading.group(1))))
            chunks.append(f"<h{level}>{markdown_inline(heading.group(2), known_figures)}</h{level}>")
        elif all(line.startswith("- ") for line in block.splitlines() if line.strip()):
            items = "\n".join(
                f"<li>{markdown_inline(line.strip()[2:], known_figures)}</li>"
                for line in block.splitlines()
                if line.strip()
            )
            chunks.append(f'<ul class="manuscript-list">\n{items}\n</ul>')
        else:
            chunks.append(f"<p>{markdown_inline(block, known_figures)}</p>")
    return "\n".join(chunks)


def callout_to_figure_id(label: str, known_figures: set[str]) -> str:
    figure_id = label_to_figure_id(label)
    if figure_id in known_figures:
        return figure_id
    panel_parent = re.match(r"^(F|S)(\d+)[A-Z]$", figure_id)
    if panel_parent:
        parent = f"{panel_parent.group(1)}{panel_parent.group(2)}"
        if parent in known_figures:
            return parent
    return figure_id


def link_figure_refs(html_text: str, known_figures: set[str]) -> str:
    def replace(match: re.Match[str]) -> str:
        figure_label = match.group(1)
        figure_id = callout_to_figure_id(figure_label, known_figures)
        return f'<a href="#{figure_anchor(figure_id)}">{figure_label}</a>'

    return re.sub(r"\b(Figure\s+S?\d+[A-Z]?(?:\s+continued)?)\b", replace, html_text)


def choose_manifest_column(fieldnames: list[str], candidates: list[str], label: str) -> str:
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(f"Manifest lacks a {label} column; tried {candidates}")


def load_figures(root: Path, config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    manifest_path = root / config["manifest"]
    figures_dir = root / config["figures_dir"]
    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Manifest has no header: {manifest_path}")
        figure_col = choose_manifest_column(
            reader.fieldnames,
            ["current_figure_name", "figure_id", "figure", "id"],
            "figure id",
        )
        panel_col = next((name for name in ["panel_id", "panel", "subpanel"] if name in reader.fieldnames), None)
        asset_col = next((name for name in ["asset_path", "final_image", "image_path"] if name in reader.fieldnames), None)
        rows = list(reader)

    figures: dict[str, dict[str, Any]] = {}
    for row in rows:
        figure_id = row.get(figure_col, "").strip()
        if not figure_id or figure_id in {"NA", "no_output"}:
            continue
        fig = figures.setdefault(
            figure_id,
            {
                "representative": row,
                "panels": [],
                "main_or_supplement": "main" if figure_id.startswith("F") else "supplement",
            },
        )
        panel_id = row.get(panel_col, "") if panel_col else ""
        if panel_id not in {"", "whole_figure", "no_output", "NA"}:
            fig["panels"].append(row)
        fig["representative"] = row

    for figure_id, fig in figures.items():
        asset_path = Path(fig["representative"].get(asset_col, "")) if asset_col else Path()
        image_path = root / asset_path if asset_path and not asset_path.is_absolute() else asset_path
        if not image_path or not image_path.exists():
            image_path = figures_dir / f"{figure_id}.png"
        if not image_path.exists():
            raise FileNotFoundError(f"Missing image for {figure_id}: {image_path}")
        if config.get("embed_images", True):
            encoded = base64.b64encode(image_path.read_bytes()).decode("ascii")
            fig["html_image_path"] = f"data:image/png;base64,{encoded}"
        else:
            fig["html_image_path"] = str(image_path.relative_to(root))
    return figures


def load_legends(root: Path, config: dict[str, Any]) -> dict[str, dict[str, str]]:
    legends_path = root / config["legends"]
    legends: dict[str, dict[str, str]] = {}
    current: str | None = None
    body_lines: list[str] = []

    def finish_current() -> None:
        if current is None:
            return
        body = "\n".join(body_lines).strip()
        if not body:
            raise ValueError(f"Empty legend body for {current}")
        legends[current]["body"] = body

    for line in legends_path.read_text().splitlines():
        header = re.match(r"^##\s+(Figure\s+S?\d+[A-Z]?(?:\s+continued)?)\.?\s*(.*)$", line)
        if header:
            finish_current()
            current = label_to_figure_id(header.group(1))
            if current in legends:
                raise ValueError(f"Duplicate legend block for {current}")
            title = header.group(2).strip().lstrip(".").strip() or display_label(current)
            legends[current] = {"title": title, "body": ""}
            body_lines = []
            continue
        if current:
            body_lines.append(line)
    finish_current()
    return legends


def write_legend_validation_report(root: Path, config: dict[str, Any], figures: dict[str, Any], legends: dict[str, Any]) -> None:
    expected = set(figures)
    observed = set(legends)
    missing = sorted(expected - observed, key=figure_sort_key)
    extra = sorted(observed - expected, key=figure_sort_key)
    status = "PASS" if not missing and not extra else "FAIL"
    lines = [
        "# Draft Render Legend Validation",
        "",
        f"Status: {status}",
        "",
        f"Expected rendered figure legends: {len(expected)}",
        f"Observed legend blocks: {len(observed)}",
        "",
    ]
    if missing:
        lines.extend(["Missing legend blocks:", *[f"- {display_label(item)}" for item in missing], ""])
    if extra:
        lines.extend(["Extra legend blocks:", *[f"- {display_label(item)}" for item in extra], ""])
    if status == "PASS":
        lines.append("Every rendered figure has exactly one legend block and no extra figure legend blocks were found.")
    report = root / config["legend_validation_report"]
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines) + "\n")
    if status != "PASS":
        raise ValueError(f"Legend validation failed; see {report}")


def load_section(root: Path, section_id: str, config: dict[str, Any]) -> dict[str, str]:
    path = root / config["sections_dir"] / f"{section_id}.md"
    metadata, body = read_markdown(path)
    text = extract_named_body(body, ["Manuscript Text", "Text"])
    return {
        "id": section_id,
        "title": metadata.get("title", section_id.replace("_", " ").title()),
        "text": text,
        "path": str(path.relative_to(root)),
    }


def load_result_sections(root: Path, config: dict[str, Any]) -> list[dict[str, str]]:
    results_dir = root / config["results_dir"]
    configured = [str(item) for item in config.get("results_order", [])]
    if configured:
        paths = [results_dir / (item if item.endswith(".md") else f"{item}.md") for item in configured]
    else:
        paths = sorted(results_dir.glob("*.md"))
    if not paths:
        raise FileNotFoundError(f"No Results sidecars found in {results_dir}")

    sections = []
    for path in paths:
        metadata, body = read_markdown(path)
        section_id = metadata.get("section_id", path.stem)
        sections.append(
            {
                "id": section_id,
                "title": metadata.get("title", section_id.replace("_", " ").title()),
                "text": extract_named_body(body, ["Results Text", "Manuscript Text", "Text"]),
                "path": str(path.relative_to(root)),
            }
        )
    return sections


def figure_block(figure_id: str, figures: dict[str, Any], legends: dict[str, Any], *, compact: bool = False) -> str:
    legend = legends[figure_id]
    css_class = "figure-card compact" if compact else "figure-card"
    return f"""
<figure class="{css_class}" id="{figure_anchor(figure_id)}">
  <img src="{e(figures[figure_id]["html_image_path"])}" alt="{e(display_label(figure_id) + '. ' + legend["title"])}" loading="lazy">
  <figcaption>
    <strong>{e(display_label(figure_id))}. {e(legend["title"])}</strong>
    <div class="legend-body">{markdown_blocks(legend["body"], set(figures))}</div>
  </figcaption>
</figure>
""".strip()


def render_toc(sections: list[dict[str, str]]) -> str:
    links = "\n".join(f'  <a href="#{slug(section["id"])}">{e(section["title"])}</a>' for section in sections)
    return f"""
<nav class="toc" aria-label="Table of contents">
  <strong>Contents</strong>
{links}
  <a href="#figures">Figure Index</a>
  <a href="#audit-trail">Audit Trail</a>
</nav>
""".strip()


def render_figure_index(figures: dict[str, Any], legends: dict[str, Any]) -> str:
    blocks = "\n".join(
        figure_block(figure_id, figures, legends, compact=True)
        for figure_id in sorted(figures, key=figure_sort_key)
    )
    return f"""
<section id="figures">
  <h2>Figure Index</h2>
  <div class="supplement-list">{blocks}</div>
</section>
""".strip()


def render_audit_trail(config: dict[str, Any], source_paths: list[str]) -> str:
    sources = [
        config["manifest"],
        config["legends"],
        config["sections_dir"],
        config["results_dir"],
        *[str(item) for item in config.get("audit_sources", [])],
        *source_paths,
    ]
    unique_sources = list(dict.fromkeys(sources))
    items = "\n".join(f"<li><code>{e(path)}</code></li>" for path in unique_sources)
    return f"""
<section id="audit-trail">
  <h2>Audit Trail</h2>
  <p>This review copy was generated mechanically from package-local source, figure, and legend inputs. Provenance, feedback, and validation notes are kept out of the journal-facing body and retained in package audit files.</p>
  <ul class="provenance-list">{items}</ul>
</section>
""".strip()


def render_html(root: Path, config: dict[str, Any]) -> str:
    figures = load_figures(root, config)
    legends = load_legends(root, config)
    write_legend_validation_report(root, config, figures, legends)
    known_figures = set(figures)

    sections: list[dict[str, str]] = []
    result_sections = load_result_sections(root, config)
    figure_placements = {str(key): str(value) for key, value in config.get("figure_placements", {}).items()}
    for section_id in config["section_order"]:
        if section_id == "results":
            sections.append({"id": "results", "title": "Results", "text": "", "path": ""})
            sections.extend(result_sections)
        else:
            sections.append(load_section(root, section_id, config))

    body_chunks = []
    for section in sections:
        if section["id"] == "results" and not section["text"]:
            body_chunks.append('<section id="results"><h2>Results</h2></section>')
            continue
        figure_html = ""
        placed_figure = figure_placements.get(section["id"])
        if placed_figure:
            if placed_figure not in figures:
                raise ValueError(f"Configured figure placement not found: {placed_figure}")
            figure_html = "\n" + figure_block(placed_figure, figures, legends)
        body_chunks.append(
            f"""
<section id="{slug(section["id"])}">
  <h2>{e(section["title"])}</h2>
  {markdown_blocks(section["text"], known_figures)}
  {figure_html}
</section>
""".strip()
        )

    source_paths = [section["path"] for section in sections if section.get("path")]
    return HTML_TEMPLATE.format(
        title=e(config["title"]),
        toc=render_toc(sections),
        body="\n\n".join(body_chunks),
        figure_index=render_figure_index(figures, legends),
        audit_trail=render_audit_trail(config, source_paths),
    )


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{title}</title>
  <style>
    :root {{
      --ink: #1b1f23;
      --muted: #5d6874;
      --line: #d7dde3;
      --paper: #ffffff;
      --wash: #f5f7f9;
      --accent: #1f6f78;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background: var(--paper);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      font-size: 17px;
      line-height: 1.58;
    }}
    a {{ color: var(--accent); text-decoration-thickness: 0.08em; }}
    .page {{
      display: grid;
      grid-template-columns: minmax(190px, 250px) minmax(0, 1fr);
      gap: 32px;
      max-width: 1320px;
      margin: 0 auto;
      padding: 32px 28px 72px;
    }}
    header {{
      grid-column: 1 / -1;
      border-bottom: 1px solid var(--line);
      padding: 36px 0 28px;
    }}
    .kicker {{
      color: var(--accent);
      font-size: 0.82rem;
      font-weight: 700;
      letter-spacing: 0.08em;
      text-transform: uppercase;
      margin: 0 0 12px;
    }}
    h1 {{
      margin: 0;
      max-width: 980px;
      font-size: clamp(2.1rem, 4.6vw, 4.2rem);
      line-height: 1.04;
      font-weight: 760;
    }}
    h2 {{ margin: 0 0 18px; font-size: 1.75rem; line-height: 1.18; }}
    h3 {{ margin: 0 0 14px; font-size: 1.28rem; line-height: 1.25; }}
    p {{ max-width: 900px; margin: 0 0 16px; }}
    .toc {{
      position: sticky;
      top: 20px;
      align-self: start;
      padding: 16px 0;
      font-size: 0.92rem;
    }}
    .toc strong {{
      display: block;
      color: var(--muted);
      margin: 0 0 8px;
      font-size: 0.78rem;
      text-transform: uppercase;
      letter-spacing: 0.08em;
    }}
    .toc a {{
      display: block;
      padding: 6px 0;
      color: var(--ink);
      text-decoration: none;
    }}
    main {{ min-width: 0; }}
    section {{ padding: 34px 0; border-bottom: 1px solid var(--line); }}
    .figure-card {{
      margin: 28px 0 8px;
      border: 1px solid var(--line);
      border-radius: 6px;
      overflow: hidden;
      background: #fff;
    }}
    .figure-card img {{
      display: block;
      width: 100%;
      max-height: 86vh;
      object-fit: contain;
      background: white;
      border-bottom: 1px solid var(--line);
    }}
    .figure-card figcaption {{
      padding: 14px 16px 16px;
      color: var(--muted);
      font-size: 0.9rem;
      line-height: 1.45;
    }}
    .figure-card figcaption strong {{
      display: block;
      color: var(--ink);
      margin-bottom: 8px;
      font-size: 0.98rem;
    }}
    .legend-body p {{ max-width: none; margin: 0 0 10px; }}
    .supplement-list {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(310px, 1fr));
      gap: 18px;
      margin-top: 18px;
    }}
    .figure-card.compact {{ margin: 0; }}
    .figure-card.compact img {{ max-height: 620px; }}
    code {{
      font-family: "SFMono-Regular", Consolas, "Liberation Mono", monospace;
      font-size: 0.88em;
      background: var(--wash);
      border: 1px solid var(--line);
      border-radius: 4px;
      padding: 0.08rem 0.26rem;
    }}
    .manuscript-list {{ max-width: 900px; margin: 0 0 16px 22px; padding: 0; }}
    .manuscript-list li {{ margin: 0 0 8px; }}
    .provenance-list {{
      margin: 0 0 0 20px;
      padding: 0;
      color: var(--muted);
      font-size: 0.95rem;
    }}
    @media (max-width: 900px) {{
      .page {{ display: block; padding: 22px 18px 54px; }}
      .toc {{ position: static; border-bottom: 1px solid var(--line); margin-bottom: 16px; }}
      h1 {{ font-size: 2.2rem; }}
      .supplement-list {{ grid-template-columns: 1fr; }}
    }}
  </style>
</head>
<body>
  <div class="page">
    <header>
      <p class="kicker">Manuscript draft</p>
      <h1>{title}</h1>
    </header>
    {toc}
    <main>
      {body}
      {figure_index}
      {audit_trail}
    </main>
  </div>
</body>
</html>
"""


def validate_output(html_text: str, figures: dict[str, Any], embed_images: bool) -> None:
    srcs = re.findall(r'<img src="([^"]+)"', html_text)
    if len(srcs) != len(figures):
        raise ValueError(f"Rendered {len(srcs)} images, expected {len(figures)}")
    if embed_images:
        non_embedded = [src for src in srcs if not src.startswith("data:image/png;base64,")]
        if non_embedded:
            raise ValueError("Non-embedded image sources: " + ", ".join(non_embedded))
        for src in srcs:
            base64.b64decode(src.removeprefix("data:image/png;base64,"), validate=True)
    anchors = re.findall(r'id="([^"]+)"', html_text)
    duplicates = sorted({anchor for anchor in anchors if anchors.count(anchor) > 1})
    if duplicates:
        raise ValueError("Duplicate HTML anchors: " + ", ".join(duplicates))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-root", default=".", help="Assembly package root.")
    parser.add_argument("--config", type=Path, help="Optional renderer_config.json path.")
    args = parser.parse_args()

    root = Path(args.assembly_root).resolve()
    config = read_config(root, args.config)
    figures = load_figures(root, config)
    html_text = render_html(root, config)
    validate_output(html_text, figures, bool(config.get("embed_images", True)))

    output = root / config["output"]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html_text)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
