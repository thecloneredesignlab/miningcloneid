#!/usr/bin/env python3
"""Convert a text/LaTeX manuscript into a navigable HTML report.

The converter is intentionally dependency-free. It handles the manuscript
features needed for quick review reports: LaTeX section headings, paragraphs,
MathJax-compatible equations, item lists, simple tabular blocks, and figures.
It is not a replacement for a full LaTeX compiler.
"""

from __future__ import annotations

import argparse
import base64
import html
import mimetypes
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


SECTION_LEVELS = {
    "part": 1,
    "chapter": 1,
    "section": 2,
    "subsection": 3,
    "subsubsection": 4,
    "paragraph": 5,
    "subparagraph": 6,
}

MATH_ENVS = {
    "equation",
    "equation*",
    "align",
    "align*",
    "gather",
    "gather*",
    "multline",
    "multline*",
    "eqnarray",
    "eqnarray*",
}

COMMON_IMAGE_EXTENSIONS = (".png", ".jpg", ".jpeg", ".svg", ".webp", ".gif", ".pdf")
LABEL_DISPLAY_TEXT: dict[str, str] = {}


@dataclass
class Section:
    id: str
    title: str
    level: int
    parent_id: str


@dataclass
class Figure:
    id: str
    label: str
    caption: str
    number_label: str


@dataclass
class FigureBlock:
    html: str
    figure: Figure


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a plain text or LaTeX manuscript to a self-contained HTML report."
    )
    parser.add_argument("--input", required=True, help="Input .tex/.txt/.md file.")
    parser.add_argument("--output", required=True, help="Output HTML file.")
    parser.add_argument("--title", default=None, help="Override report title.")
    parser.add_argument(
        "--hide-source",
        action="store_true",
        help="Omit the input file path from the visible report header.",
    )
    parser.add_argument(
        "--asset-mode",
        choices=("embed", "link"),
        default="embed",
        help="Embed image assets as data URIs or link them from the HTML.",
    )
    parser.add_argument(
        "--figure-placement",
        choices=("first-reference", "source"),
        default="first-reference",
        help="Place figures at their first text reference or preserve their source order.",
    )
    parser.add_argument(
        "--mathjax-url",
        default="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js",
        help="MathJax script URL used by the generated HTML.",
    )
    return parser.parse_args()


def strip_tex_comments(text: str) -> str:
    out_lines = []
    for line in text.splitlines():
        escaped = False
        kept = []
        for ch in line:
            if ch == "%" and not escaped:
                break
            kept.append(ch)
            escaped = ch == "\\" and not escaped
            if ch != "\\":
                escaped = False
        out_lines.append("".join(kept).rstrip())
    return "\n".join(out_lines)


def find_matching_brace(text: str, open_index: int) -> int:
    depth = 0
    escaped = False
    for i in range(open_index, len(text)):
        ch = text[i]
        if ch == "\\" and not escaped:
            escaped = True
            continue
        if ch == "{" and not escaped:
            depth += 1
        elif ch == "}" and not escaped:
            depth -= 1
            if depth == 0:
                return i
        escaped = False
    return -1


def extract_braced_after(text: str, start_index: int) -> tuple[str, int] | None:
    open_index = text.find("{", start_index)
    if open_index < 0:
        return None
    close_index = find_matching_brace(text, open_index)
    if close_index < 0:
        return None
    return text[open_index + 1 : close_index], close_index + 1


def extract_command_arg(text: str, command: str) -> str | None:
    match = re.search(r"\\" + re.escape(command) + r"\s*(?:\[[^\]]*\])?\s*{", text)
    if not match:
        return None
    extracted = extract_braced_after(text, match.end() - 1)
    return extracted[0] if extracted else None


def extract_reference_labels(text: str) -> list[str]:
    labels: list[str] = []
    for command in ("ref", "eqref", "autoref", "cref", "Cref"):
        pattern = "\\" + command
        start = text.find(pattern)
        while start >= 0:
            extracted = extract_braced_after(text, start + len(pattern))
            if not extracted:
                break
            content, end = extracted
            labels.extend(label.strip() for label in content.split(",") if label.strip())
            start = text.find(pattern, end)
    return labels


def document_body(text: str) -> str:
    begin = re.search(r"\\begin\s*{\s*document\s*}", text)
    end = re.search(r"\\end\s*{\s*document\s*}", text)
    if begin and end and begin.end() < end.start():
        return text[begin.end() : end.start()]
    return text


def resolve_tex_include_path(path_text: str, base_dir: Path, search_dirs: Iterable[Path]) -> Path | None:
    raw = Path(path_text.strip())
    names = [raw]
    if raw.suffix == "":
        names.append(Path(str(raw) + ".tex"))
    for name in names:
        if name.is_absolute() and name.exists():
            return name.resolve()
    for base in [base_dir, *search_dirs]:
        for name in names:
            if name.is_absolute():
                continue
            candidate = base / name
            if candidate.exists():
                return candidate.resolve()
    return None


def expand_tex_inputs(text: str, base_dir: Path, search_dirs: Iterable[Path], seen: set[Path] | None = None) -> str:
    seen = set() if seen is None else seen
    out: list[str] = []
    i = 0
    while i < len(text):
        match = re.search(r"\\(input|include)\s*{", text[i:])
        if not match:
            out.append(text[i:])
            break
        start = i + match.start()
        command = match.group(1)
        out.append(text[i:start])
        extracted = extract_braced_after(text, i + match.end() - 1)
        if not extracted:
            out.append(text[start : i + match.end()])
            i = i + match.end()
            continue
        path_text, end = extracted
        include_path = resolve_tex_include_path(path_text, base_dir, search_dirs)
        if include_path is None:
            out.append(f"\n\\begin{{quote}}Missing included TeX file: {path_text}\\end{{quote}}\n")
        elif include_path in seen:
            out.append(f"\n\\begin{{quote}}Skipped recursive TeX include: {include_path}\\end{{quote}}\n")
        else:
            included = include_path.read_text(encoding="utf-8")
            expanded = expand_tex_inputs(included, include_path.parent, search_dirs, seen | {include_path})
            out.append(f"\n% BEGIN expanded input: {include_path}\n{expanded}\n% END expanded input: {include_path}\n")
        i = end
    return "".join(out)


def slugify(text: str, fallback: str = "section") -> str:
    value = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?", " ", text)
    value = re.sub(r"[{}$\\]", " ", value)
    value = re.sub(r"[^A-Za-z0-9]+", "-", value).strip("-").lower()
    return value or fallback


def unique_id(base: str, used: set[str]) -> str:
    base = slugify(base)
    candidate = base
    i = 2
    while candidate in used:
        candidate = f"{base}-{i}"
        i += 1
    used.add(candidate)
    return candidate


def label_id(label: str) -> str:
    return "label-" + slugify(label, fallback="ref")


def reference_display_text(label: str, command: str = "ref") -> str:
    if label in LABEL_DISPLAY_TEXT:
        if command == "autoref":
            if label.startswith("fig:"):
                return f"Figure {LABEL_DISPLAY_TEXT[label]}"
            if label.startswith("tab:") or label.startswith("table:"):
                return f"Table {LABEL_DISPLAY_TEXT[label]}"
        return LABEL_DISPLAY_TEXT[label]
    return label


def reference_link_html(label: str, command: str = "ref") -> str:
    target_id = label_id(label)
    text = reference_display_text(label, command)
    attrs = [f'href="#{html.escape(target_id, quote=True)}"']
    if label.startswith("fig:") or label.startswith("tab:") or label.startswith("table:"):
        attrs.append(f'data-ref-preview-target="{html.escape(target_id, quote=True)}"')
        attrs.append(f'data-ref-preview-label="{html.escape(text, quote=True)}"')
    return f'<a {" ".join(attrs)}>{html.escape(text)}</a>'


def protect_html(value: str, replacements: list[str]) -> str:
    token = f"\x00HTML{len(replacements)}\x00"
    replacements.append(value)
    return token


def restore_html(value: str, replacements: list[str]) -> str:
    for _ in range(len(replacements) + 1):
        before = value
        for i, replacement in enumerate(replacements):
            value = value.replace(f"\x00HTML{i}\x00", replacement)
        if value == before:
            break
    return value


def convert_text_macros(text: str, replacements: list[str]) -> str:
    macros = {
        "textbf": ("<strong>", "</strong>"),
        "bfseries": ("<strong>", "</strong>"),
        "textit": ("<em>", "</em>"),
        "emph": ("<em>", "</em>"),
        "texttt": ("<code>", "</code>"),
        "underline": ("<u>", "</u>"),
        "textsuperscript": ("<sup>", "</sup>"),
        "textsubscript": ("<sub>", "</sub>"),
    }
    changed = True
    while changed:
        changed = False
        for macro, tags in macros.items():
            pattern = "\\" + macro
            start = text.find(pattern)
            while start >= 0:
                after = start + len(pattern)
                extracted = extract_braced_after(text, after)
                if not extracted:
                    break
                content, end = extracted
                inner = render_inline_tex(content)
                token = protect_html(f"{tags[0]}{inner}{tags[1]}", replacements)
                text = text[:start] + token + text[end:]
                changed = True
                start = text.find(pattern, start + len(token))
    return text


def convert_reference_macros(text: str, replacements: list[str]) -> str:
    def replace_command(command: str, builder) -> None:
        nonlocal text
        pattern = "\\" + command
        start = text.find(pattern)
        while start >= 0:
            extracted = extract_braced_after(text, start + len(pattern))
            if not extracted:
                break
            content, end = extracted
            token = protect_html(builder(content), replacements)
            text = text[:start] + token + text[end:]
            start = text.find(pattern, start + len(token))

    replace_command("ref", lambda x: reference_link_html(x, "ref"))
    replace_command("autoref", lambda x: reference_link_html(x, "autoref"))
    replace_command("eqref", lambda x: f'({reference_link_html(x, "eqref")})')
    replace_command("cite", lambda x: f'<span class="citation">[{html.escape(x)}]</span>')
    replace_command("citep", lambda x: f'<span class="citation">[{html.escape(x)}]</span>')
    replace_command("citet", lambda x: f'<span class="citation">{html.escape(x)}</span>')
    replace_command("url", lambda x: f'<a href="{html.escape(x, quote=True)}">{html.escape(x)}</a>')
    return text


def split_math_segments(text: str) -> list[tuple[str, bool]]:
    segments: list[tuple[str, bool]] = []
    i = 0
    buf = []
    while i < len(text):
        delimiter = None
        if text.startswith("\\(", i):
            delimiter = ("\\(", "\\)")
        elif text.startswith("\\[", i):
            delimiter = ("\\[", "\\]")
        elif text.startswith("$$", i):
            delimiter = ("$$", "$$")
        elif text[i] == "$" and (i == 0 or text[i - 1] != "\\"):
            delimiter = ("$", "$")
        if delimiter:
            if buf:
                segments.append(("".join(buf), False))
                buf = []
            start_delim, end_delim = delimiter
            end = text.find(end_delim, i + len(start_delim))
            if end < 0:
                buf.append(text[i])
                i += 1
            else:
                segments.append((text[i : end + len(end_delim)], True))
                i = end + len(end_delim)
        else:
            buf.append(text[i])
            i += 1
    if buf:
        segments.append(("".join(buf), False))
    return segments


def render_inline_tex(text: str) -> str:
    rendered = []
    for segment, is_math in split_math_segments(text):
        if is_math:
            rendered.append(html.escape(segment, quote=False))
            continue
        replacements: list[str] = []
        segment = segment.replace("~", " ")
        segment = segment.replace(r"\&", "&").replace(r"\%", "%").replace(r"\_", "_")
        segment = segment.replace(r"\$", "$").replace(r"\#", "#").replace(r"\{", "{").replace(r"\}", "}")
        segment = convert_text_macros(segment, replacements)
        segment = convert_reference_macros(segment, replacements)
        segment = re.sub(r"\\(textit|textbf|emph|texttt|underline)\s+", "", segment)
        segment = re.sub(r"\\(noindent|centering|small|footnotesize|normalsize|large|Large)\b", "", segment)
        segment = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?(?:{([^{}]*)})?", lambda m: m.group(1) or "", segment)
        segment = html.escape(segment, quote=False)
        segment = restore_html(segment, replacements)
        rendered.append(segment)
    return "".join(rendered)


def plain_text_from_tex(text: str) -> str:
    value = re.sub(r"\$([^$]*)\$", r"\1", text)
    value = re.sub(r"\\\((.*?)\\\)", r"\1", value)
    value = re.sub(r"\\\[(.*?)\\\]", r"\1", value, flags=re.S)
    value = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^\]]*\])?{([^{}]*)}", r"\1", value)
    value = re.sub(r"[{}]", "", value)
    return re.sub(r"\s+", " ", value).strip()


def remove_tex_command_with_arg(text: str, command: str) -> str:
    pattern = "\\" + command
    start = text.find(pattern)
    while start >= 0:
        extracted = extract_braced_after(text, start + len(pattern))
        if not extracted:
            break
        _content, end = extracted
        while end < len(text) and text[end] in " \t\r\n":
            end += 1
        if end + 1 < len(text) and text[end : end + 2] == r"\\":
            end += 2
        text = text[:start] + text[end:]
        start = text.find(pattern, start)
    return text


def parse_heading_line(line: str) -> tuple[str, int, str, str] | None:
    match = re.match(
        r"^\s*\\(part|chapter|section|subsection|subsubsection|paragraph|subparagraph)\*?"
        r"\s*(?:\[[^\]]*\])?\s*{",
        line,
    )
    if not match:
        return None
    extracted = extract_braced_after(line, match.end() - 1)
    if not extracted:
        return None
    title, end = extracted
    command = match.group(1)
    return command, SECTION_LEVELS[command], title.strip(), line[end:].strip()


def begin_environment(line: str) -> str | None:
    match = re.match(r"^\s*\\begin\s*{\s*([^}]+)\s*}", line)
    return match.group(1).strip() if match else None


def collect_environment(lines: list[str], start: int, env: str) -> tuple[str, int]:
    collected = [lines[start]]
    i = start + 1
    end_pattern = re.compile(r"\\end\s*{\s*" + re.escape(env) + r"\s*}")
    while i < len(lines):
        collected.append(lines[i])
        if end_pattern.search(lines[i]):
            return "\n".join(collected), i + 1
        i += 1
    return "\n".join(collected), i


def parse_graphic_paths(raw_text: str, input_dir: Path) -> list[Path]:
    out: list[Path] = []
    start = raw_text.find(r"\graphicspath")
    while start >= 0:
        extracted = extract_braced_after(raw_text, start + len(r"\graphicspath"))
        if not extracted:
            break
        content, end = extracted
        for entry in re.findall(r"{([^{}]+)}", content):
            path = Path(entry)
            out.append(path if path.is_absolute() else (input_dir / path))
        start = raw_text.find(r"\graphicspath", end)
    return out


def find_repo_root(start: Path) -> Path:
    current = start.resolve()
    for parent in (current, *current.parents):
        if (parent / ".git").exists():
            return parent
    return current


def candidate_image_paths(path_text: str, search_dirs: Iterable[Path]) -> list[Path]:
    raw = Path(path_text)
    candidates: list[Path] = []
    names = [raw]
    if raw.suffix == "":
        names = [Path(str(raw) + ext) for ext in COMMON_IMAGE_EXTENSIONS]
    for base in search_dirs:
        for name in names:
            candidates.append(name if name.is_absolute() else base / name)
    return candidates


def resolve_image_path(path_text: str, search_dirs: Iterable[Path]) -> Path | None:
    for candidate in candidate_image_paths(path_text, search_dirs):
        if candidate.exists():
            return candidate.resolve()
    return None


def asset_src(path: Path, output_dir: Path, mode: str) -> str:
    if mode == "link":
        return os.path.relpath(path, output_dir)
    mime = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{data}"


def render_figure(
    env_text: str,
    search_dirs: list[Path],
    output_dir: Path,
    mode: str,
    number_label: str,
    float_name: str = "Figure",
    continued: bool = False,
) -> tuple[str, Figure | None]:
    label = extract_command_arg(env_text, "label") or ""
    caption = extract_command_arg(env_text, "caption") or ""
    include_match = re.search(r"\\includegraphics\s*(?:\[[^\]]*\])?\s*{([^{}]+)}", env_text)
    if not include_match:
        return "", None
    image_text = include_match.group(1).strip()
    image_path = resolve_image_path(image_text, search_dirs)
    figure_id = label_id(label) if label else "figure-" + slugify(image_text, fallback="image")
    caption_html = render_inline_tex(caption) if caption else html.escape(image_text)
    if image_path and image_path.suffix.lower() != ".pdf":
        src = asset_src(image_path, output_dir, mode)
        original_src = asset_src(image_path, output_dir, "link")
        image_html = (
            f'<img src="{html.escape(src, quote=True)}" '
            f'alt="{html.escape(plain_text_from_tex(caption) or image_text, quote=True)}" '
            f'data-source-path="{html.escape(str(image_path), quote=True)}" '
            f'data-lightbox-image tabindex="0">'
        )
        source_link = (
            f'<p class="figure-source"><a href="{html.escape(original_src, quote=True)}" '
            f'target="_blank" rel="noopener" data-source-viewer>'
            f"Open source image</a></p>"
        )
    elif image_path:
        src = asset_src(image_path, output_dir, "link")
        image_html = f'<p class="asset-link"><a href="{html.escape(src, quote=True)}">Open figure file: {html.escape(image_path.name)}</a></p>'
        source_link = ""
    else:
        image_html = f'<p class="missing-asset">Missing image: {html.escape(image_text)}</p>'
        source_link = ""
    float_label = f"{float_name} {number_label}"
    if continued:
        float_label += " (continued)"
    block = (
        f'<figure id="{html.escape(figure_id, quote=True)}" class="tex-figure">'
        f"{image_html}"
        f'<figcaption><span class="float-label">{html.escape(float_label)}.</span> {caption_html}</figcaption>'
        f"{source_link}"
        f"</figure>"
    )
    return block, Figure(
        id=figure_id,
        label=label,
        caption=plain_text_from_tex(caption) or image_text,
        number_label=number_label,
    )


def parse_figure_counter_reset(text: str) -> int | None:
    match = re.search(r"\\setcounter\s*\{\s*figure\s*\}\s*\{\s*(-?\d+)\s*\}", text)
    return int(match.group(1)) if match else None


def parse_figure_name(text: str) -> str | None:
    match = re.search(r"\\renewcommand\s*\{\s*\\figurename\s*\}\s*\{([^{}]+)\}", text)
    if not match:
        return None
    name = plain_text_from_tex(match.group(1)).strip()
    return name or None


def parse_table_counter_reset(text: str) -> int | None:
    match = re.search(r"\\setcounter\s*\{\s*table\s*\}\s*\{\s*(-?\d+)\s*\}", text)
    return int(match.group(1)) if match else None


def parse_table_name(text: str) -> str | None:
    match = re.search(r"\\renewcommand\s*\{\s*\\tablename\s*\}\s*\{([^{}]+)\}", text)
    if not match:
        return None
    name = plain_text_from_tex(match.group(1)).strip()
    return name or None


def is_numbered_table(env: str, env_text: str) -> bool:
    return env in {"table", "longtable"} and bool(extract_command_arg(env_text, "caption"))


def collect_figure_blocks(lines: list[str], search_dirs: list[Path], output_dir: Path, mode: str) -> dict[str, list[FigureBlock]]:
    figure_blocks: dict[str, list[FigureBlock]] = {}
    current_label = ""
    current_number = 0
    current_name = "Figure"
    i = 0
    while i < len(lines):
        counter_reset = parse_figure_counter_reset(lines[i])
        if counter_reset is not None:
            current_number = counter_reset
            i += 1
            continue
        figure_name = parse_figure_name(lines[i])
        if figure_name is not None:
            current_name = figure_name
            i += 1
            continue
        env = begin_environment(lines[i])
        if env == "figure":
            env_text, next_i = collect_environment(lines, i, env)
            is_continued = r"\ContinuedFloat" in env_text and current_number > 0
            if not is_continued:
                current_number += 1
            number_label = str(current_number)
            block, figure = render_figure(
                env_text,
                search_dirs,
                output_dir,
                mode,
                number_label=number_label,
                float_name=current_name,
                continued=is_continued,
            )
            if block and figure and figure.label:
                current_label = figure.label
                figure_blocks.setdefault(current_label, []).append(FigureBlock(html=block, figure=figure))
            elif block and figure and current_label and r"\ContinuedFloat" in env_text:
                figure_blocks.setdefault(current_label, []).append(FigureBlock(html=block, figure=figure))
            i = next_i
        else:
            i += 1
    return figure_blocks


def collect_figure_labels(lines: list[str]) -> dict[str, str]:
    figure_labels: dict[str, str] = {}
    current_number = 0
    i = 0
    while i < len(lines):
        counter_reset = parse_figure_counter_reset(lines[i])
        if counter_reset is not None:
            current_number = counter_reset
            i += 1
            continue
        env = begin_environment(lines[i])
        if env == "figure":
            env_text, next_i = collect_environment(lines, i, env)
            is_continued = r"\ContinuedFloat" in env_text and current_number > 0
            if not is_continued:
                current_number += 1
            label = extract_command_arg(env_text, "label") or ""
            if label:
                figure_labels[label] = str(current_number)
            i = next_i
        else:
            i += 1
    return figure_labels


def collect_references_outside_figures(lines: list[str]) -> set[str]:
    labels: set[str] = set()
    i = 0
    while i < len(lines):
        env = begin_environment(lines[i])
        if env == "figure":
            _env_text, i = collect_environment(lines, i, env)
            continue
        labels.update(extract_reference_labels(lines[i]))
        i += 1
    return labels


def collect_table_labels(lines: list[str]) -> dict[str, str]:
    table_labels: dict[str, str] = {}
    table_number = 0
    i = 0
    while i < len(lines):
        counter_reset = parse_table_counter_reset(lines[i])
        if counter_reset is not None:
            table_number = counter_reset
            i += 1
            continue
        env = begin_environment(lines[i])
        if env in {"tabular", "longtable", "center", "table"}:
            env_text, next_i = collect_environment(lines, i, env)
            if r"\begin{tabular" in env_text or r"\begin{longtable" in env_text or env in {"tabular", "longtable"}:
                if is_numbered_table(env, env_text):
                    table_number += 1
                    label = extract_command_arg(env_text, "label") or ""
                    if label:
                        table_labels[label] = str(table_number)
            i = next_i
        else:
            i += 1
    return table_labels


def render_math_environment(env_text: str) -> str:
    return f'<div class="math-block">{html.escape(env_text, quote=False)}</div>'


def render_tabular(env_text: str) -> str:
    begin_match = re.search(r"\\begin\s*{\s*(tabular|longtable)\s*}", env_text)
    if not begin_match:
        return f"<pre>{html.escape(env_text)}</pre>"
    table_env = begin_match.group(1)
    spec_start = env_text.find("{", begin_match.end())
    if spec_start < 0:
        return f"<pre>{html.escape(env_text)}</pre>"
    spec_end = find_matching_brace(env_text, spec_start)
    if spec_end < 0:
        return f"<pre>{html.escape(env_text)}</pre>"
    end_match = re.search(r"\\end\s*{\s*" + re.escape(table_env) + r"\s*}", env_text[spec_end + 1 :])
    if not end_match:
        return f"<pre>{html.escape(env_text)}</pre>"
    content = env_text[spec_end + 1 : spec_end + 1 + end_match.start()]
    content = remove_tex_command_with_arg(content, "caption")
    content = remove_tex_command_with_arg(content, "label")
    content = re.sub(r"\\(toprule|midrule|bottomrule|hline|endfirsthead|endhead|endfoot|endlastfoot)\b", "", content)
    rows = []
    for row in re.split(r"\\\\", content):
        row = row.strip()
        if not row:
            continue
        cells = [render_inline_tex(cell.strip()) for cell in row.split("&")]
        rows.append(cells)
    if not rows:
        return ""
    max_cols = max(len(row) for row in rows)
    rendered_rows = []
    for row_index, row in enumerate(rows):
        tag = "th" if row_index == 0 else "td"
        cells = row + [""] * (max_cols - len(row))
        rendered_rows.append("<tr>" + "".join(f"<{tag}>{cell}</{tag}>" for cell in cells) + "</tr>")
    return '<div class="table-wrap"><table class="tex-table">' + "".join(rendered_rows) + "</table></div>"


def render_table_block(
    env_text: str,
    table_number: int | None,
    table_name: str = "Table",
    unnumbered_index: int = 0,
) -> tuple[str, str]:
    caption = extract_command_arg(env_text, "caption") or ""
    label = extract_command_arg(env_text, "label") or ""
    if label:
        table_id = label_id(label)
    elif table_number is not None:
        table_id = f"table-{table_number}"
    else:
        table_id = f"table-unnumbered-{unnumbered_index}"
    if label and table_number is not None:
        LABEL_DISPLAY_TEXT[label] = str(table_number)
    caption_html = render_inline_tex(caption) if caption else ""
    title_html = ""
    if caption_html:
        title_html = '<p class="table-caption">'
        if table_number is not None:
            title_html += (
                f'<span class="float-label">{html.escape(table_name)} '
                f'{table_number}.</span> '
            )
        title_html += f"{caption_html}</p>"
    return (
        f'<div id="{html.escape(table_id, quote=True)}" class="table-block">'
        f"{title_html}{render_tabular(env_text)}</div>",
        label,
    )


def render_list(env_text: str, ordered: bool) -> str:
    content = re.sub(r"^.*?\\begin\s*{[^}]+}", "", env_text, flags=re.S)
    content = re.sub(r"\\end\s*{[^}]+}.*$", "", content, flags=re.S)
    items = [item.strip() for item in re.split(r"\\item(?:\[[^\]]*\])?", content) if item.strip()]
    tag = "ol" if ordered else "ul"
    return f"<{tag}>" + "".join(f"<li>{render_inline_tex(item)}</li>" for item in items) + f"</{tag}>"


def render_paragraph(lines: list[str]) -> str:
    text = " ".join(line.strip() for line in lines).strip()
    if not text:
        return ""
    text = re.sub(r"\\\\\s*", "<br>", text)
    return f"<p>{render_inline_tex(text)}</p>"


def nav_html(sections: list[Section]) -> str:
    children: dict[str, list[Section]] = {"": []}
    for section in sections:
        children.setdefault(section.parent_id, []).append(section)
        children.setdefault(section.id, [])

    def render_node(section: Section) -> str:
        child_html = "".join(render_node(child) for child in children.get(section.id, []))
        link = (
            f'<a href="#{html.escape(section.id, quote=True)}" '
            f'data-target="{html.escape(section.id, quote=True)}">{render_inline_tex(section.title)}</a>'
        )
        if child_html:
            return (
                f'<details class="nav-branch nav-depth-{min(section.level, 6)}">'
                f"<summary>{link}</summary>{child_html}</details>"
            )
        return f'<div class="nav-leaf nav-depth-{min(section.level, 6)}">{link}</div>'

    if not sections:
        return '<p class="nav-empty">No sections found.</p>'
    return '<div class="toc-tree">' + "".join(render_node(section) for section in children[""]) + "</div>"


def render_document(
    raw_text: str,
    input_path: Path,
    output_path: Path,
    title_override: str | None,
    asset_mode: str,
    mathjax_url: str,
    hide_source: bool = False,
    figure_placement: str = "first-reference",
) -> str:
    input_dir = input_path.parent.resolve()
    output_dir = output_path.parent.resolve()
    repo_root = find_repo_root(input_dir)
    initial_search_dirs = [input_dir, Path.cwd().resolve(), repo_root]
    expanded_text = expand_tex_inputs(raw_text, input_dir, initial_search_dirs, seen={input_path.resolve()})
    graphic_dirs = [p.resolve() for p in parse_graphic_paths(expanded_text, input_dir)]
    search_dirs = [input_dir, *graphic_dirs, Path.cwd().resolve(), repo_root]
    title = title_override or extract_command_arg(expanded_text, "title") or input_path.stem
    body = document_body(strip_tex_comments(expanded_text))
    lines = body.splitlines()
    LABEL_DISPLAY_TEXT.clear()
    LABEL_DISPLAY_TEXT.update(collect_figure_labels(lines))
    LABEL_DISPLAY_TEXT.update(collect_table_labels(lines))
    figure_blocks = collect_figure_blocks(lines, search_dirs, output_dir, asset_mode)
    referenced_figure_labels = (
        collect_references_outside_figures(lines).intersection(figure_blocks)
        if figure_placement == "first-reference"
        else set()
    )

    used_ids: set[str] = set()
    sections: list[Section] = []
    figures: list[Figure] = []
    html_parts: list[str] = []
    paragraph: list[str] = []
    placed_figure_labels: set[str] = set()
    current_figure_label = ""
    current_figure_number = 0
    current_figure_name = "Figure"
    table_number = 0
    current_table_name = "Table"
    unnumbered_table_count = 0
    open_levels: list[int] = []
    parent_by_level: dict[int, str] = {}

    def append_referenced_figures(tex_text: str) -> None:
        for label in extract_reference_labels(tex_text):
            if label not in referenced_figure_labels or label in placed_figure_labels:
                continue
            for figure_block in figure_blocks[label]:
                html_parts.append(figure_block.html)
                figures.append(figure_block.figure)
            placed_figure_labels.add(label)

    def flush_paragraph() -> None:
        nonlocal paragraph
        tex_text = " ".join(line.strip() for line in paragraph).strip()
        block = render_paragraph(paragraph)
        if block:
            html_parts.append(block)
            append_referenced_figures(tex_text)
        paragraph = []

    def close_to_level(level: int) -> None:
        while open_levels and open_levels[-1] >= level:
            html_parts.append("</section>")
            open_levels.pop()

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped:
            flush_paragraph()
            i += 1
            continue
        if re.match(r"^\\(maketitle|clearpage|newpage|tableofcontents)\b", stripped):
            flush_paragraph()
            i += 1
            continue
        figure_counter_reset = parse_figure_counter_reset(stripped)
        if figure_counter_reset is not None:
            flush_paragraph()
            current_figure_number = figure_counter_reset
            i += 1
            continue
        figure_name = parse_figure_name(stripped)
        if figure_name is not None:
            flush_paragraph()
            current_figure_name = figure_name
            i += 1
            continue
        table_counter_reset = parse_table_counter_reset(stripped)
        if table_counter_reset is not None:
            flush_paragraph()
            table_number = table_counter_reset
            i += 1
            continue
        table_name = parse_table_name(stripped)
        if table_name is not None:
            flush_paragraph()
            current_table_name = table_name
            i += 1
            continue
        if re.match(r"^\\(begingroup|endgroup|scriptsize|footnotesize|small|normalsize|setlength|setcounter|newlength|renewcommand)\b", stripped):
            flush_paragraph()
            i += 1
            continue
        heading = parse_heading_line(line)
        if heading:
            flush_paragraph()
            _command, level, heading_title, leftover = heading
            close_to_level(level)
            parent_level = max([l for l in parent_by_level if l < level], default=0)
            parent_id = parent_by_level.get(parent_level, "")
            section_id = unique_id(heading_title, used_ids)
            sections.append(Section(id=section_id, title=heading_title, level=level, parent_id=parent_id))
            parent_by_level[level] = section_id
            for deeper in [l for l in list(parent_by_level) if l > level]:
                del parent_by_level[deeper]
            h_level = min(max(level, 1), 6)
            html_parts.append(
                f'<section id="{html.escape(section_id, quote=True)}" '
                f'class="tex-section tex-level-{level}" data-nav-target>'
                f"<h{h_level}>{render_inline_tex(heading_title)}</h{h_level}>"
            )
            open_levels.append(level)
            if leftover:
                paragraph.append(leftover)
            i += 1
            continue
        env = begin_environment(line)
        if env:
            flush_paragraph()
            env_text, next_i = collect_environment(lines, i, env)
            if env == "figure":
                is_continued = r"\ContinuedFloat" in env_text and current_figure_number > 0
                if not is_continued:
                    current_figure_number += 1
                block, figure = render_figure(
                    env_text,
                    search_dirs,
                    output_dir,
                    asset_mode,
                    number_label=str(current_figure_number),
                    float_name=current_figure_name,
                    continued=is_continued,
                )
                if figure and figure.label:
                    current_figure_label = figure.label
                    LABEL_DISPLAY_TEXT[figure.label] = figure.number_label
                figure_group_label = figure.label if figure and figure.label else ""
                if not figure_group_label and r"\ContinuedFloat" in env_text:
                    figure_group_label = current_figure_label
                should_keep_at_original_position = figure_group_label not in referenced_figure_labels
                if block and should_keep_at_original_position:
                    html_parts.append(block)
                if block and figure and should_keep_at_original_position:
                    figures.append(figure)
                    if figure.label:
                        placed_figure_labels.add(figure.label)
            elif env in MATH_ENVS:
                html_parts.append(render_math_environment(env_text))
            elif env in {"itemize", "enumerate"}:
                html_parts.append(render_list(env_text, ordered=(env == "enumerate")))
                append_referenced_figures(env_text)
            elif env in {"tabular", "longtable", "center", "table"}:
                if r"\begin{tabular" in env_text or r"\begin{longtable" in env_text or env in {"tabular", "longtable"}:
                    if is_numbered_table(env, env_text):
                        table_number += 1
                        rendered_table_number: int | None = table_number
                    else:
                        unnumbered_table_count += 1
                        rendered_table_number = None
                    table_html, _label = render_table_block(
                        env_text,
                        rendered_table_number,
                        table_name=current_table_name,
                        unnumbered_index=unnumbered_table_count,
                    )
                    html_parts.append(table_html)
                else:
                    paragraph.append(env_text)
                append_referenced_figures(env_text)
            else:
                paragraph.append(env_text)
            i = next_i
            continue
        label = extract_command_arg(stripped, "label")
        if label and stripped.startswith(r"\label"):
            flush_paragraph()
            html_parts.append(f'<span id="{html.escape(label_id(label), quote=True)}" class="label-anchor"></span>')
            i += 1
            continue
        paragraph.append(line)
        i += 1

    flush_paragraph()
    for label in figure_blocks:
        if label in referenced_figure_labels and label not in placed_figure_labels:
            for figure_block in figure_blocks[label]:
                html_parts.append(figure_block.html)
                figures.append(figure_block.figure)
            placed_figure_labels.add(label)
    close_to_level(0)
    figure_count = len(figures)
    subtitle = f"{len(sections)} sections"
    if figure_count:
        subtitle += f" · {figure_count} figures"
    content = "\n".join(html_parts)
    return build_html_shell(
        title=plain_text_from_tex(title),
        nav=nav_html(sections),
        content=content,
        subtitle=subtitle,
        source_path=input_path,
        mathjax_url=mathjax_url,
        hide_source=hide_source,
    )


def build_html_shell(
    title: str,
    nav: str,
    content: str,
    subtitle: str,
    source_path: Path,
    mathjax_url: str,
    hide_source: bool = False,
) -> str:
    source_meta = "" if hide_source else (
        f'\n        <p class="report-meta">Source: {html.escape(str(source_path))}</p>'
    )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <script>
    window.MathJax = {{
      tex: {{
        inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
        displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
        processEscapes: true,
        tags: 'ams'
      }},
      options: {{ skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre'] }}
    }};
  </script>
  <script async src="{html.escape(mathjax_url, quote=True)}"></script>
  <style>
{report_css()}
  </style>
</head>
<body>
  <div class="shell">
    <aside class="sidebar">
      <div class="side-head">
        <div class="kicker">Text to HTML</div>
        <div class="side-title">{html.escape(title)}</div>
        <div class="side-subtitle">{html.escape(subtitle)}</div>
      </div>
      <details class="nav-group" open>
        <summary>Contents</summary>
        {nav}
      </details>
    </aside>
    <main class="main">
      <header class="report-card">
        <h1>{html.escape(title)}</h1>{source_meta}
      </header>
      {content}
    </main>
  </div>
  <div class="image-lightbox" id="image-lightbox" aria-hidden="true">
    <div class="image-lightbox-toolbar" aria-label="Image viewer controls">
      <button type="button" data-lightbox-action="fit">Fit</button>
      <button type="button" data-lightbox-action="actual">1:1</button>
      <button type="button" data-lightbox-action="close">Close</button>
    </div>
    <div class="image-lightbox-stage">
      <img class="image-lightbox-img" alt="">
    </div>
    <div class="image-lightbox-caption"></div>
  </div>
  <div class="ref-popover" id="ref-popover" aria-hidden="true">
    <div class="ref-popover-panel" role="dialog" aria-label="Reference preview">
      <div class="ref-popover-head">
        <strong class="ref-popover-title">Reference preview</strong>
        <a class="ref-popover-jump" href="#">Jump to location</a>
      </div>
      <div class="ref-popover-body"></div>
    </div>
  </div>
  <script>
{report_script()}
  </script>
</body>
</html>
"""


def report_css() -> str:
    return """
html{scroll-behavior:smooth}body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif;background:#f5f6f8;color:#1d2730;line-height:1.58}
.shell{display:flex;gap:26px;max-width:1740px;margin:0 auto;padding:22px}.sidebar{position:sticky;top:22px;align-self:flex-start;width:330px;max-height:calc(100vh - 44px);overflow:auto;border:1px solid #d7dde4;border-radius:8px;background:#fff;box-shadow:0 8px 22px rgba(20,32,45,.08);scrollbar-gutter:stable}.side-head{padding:16px;background:#24384c;color:#fff}.kicker{font-size:11px;letter-spacing:.08em;text-transform:uppercase;opacity:.78;font-weight:700}.side-title{font-size:17px;font-weight:700;line-height:1.25;margin-top:4px}.side-subtitle{font-size:12px;line-height:1.35;margin-top:4px;opacity:.84}
.nav-group{border-bottom:1px solid #e4e9ef}.nav-group>summary{cursor:pointer;list-style:none;padding:10px 14px;font-size:12px;font-weight:700;letter-spacing:.04em;text-transform:uppercase;color:#334e68;background:#f7fafc}.nav-group summary::-webkit-details-marker,.nav-branch summary::-webkit-details-marker{display:none}.nav-group>summary:before,.nav-branch>summary:before{content:'+';display:inline-block;width:16px;color:#6b7c8f}.nav-group[open]>summary:before,.nav-branch[open]>summary:before{content:'-'}.toc-tree{padding:8px 8px 12px}.nav-branch,.nav-leaf{margin:3px 0;border-radius:7px}.nav-depth-1{margin-left:0}.nav-depth-2{margin-left:8px}.nav-depth-3{margin-left:16px}.nav-depth-4{margin-left:24px}.nav-depth-5{margin-left:32px}.nav-depth-6{margin-left:40px}.nav-branch>summary{cursor:pointer;list-style:none;padding:7px 8px;font-size:12px;line-height:1.3;color:#21384d;background:#fff;border-radius:7px}.nav-branch>summary:hover,.nav-leaf a:hover{background:#eef4fa}.nav-branch summary a,.nav-leaf a{text-decoration:none;color:inherit}.nav-leaf a{display:block;padding:7px 10px;border-radius:7px;font-size:12px;line-height:1.3;color:#21384d}.nav-leaf a.active,.nav-branch summary a.active{background:#dcecf8;color:#102a43;font-weight:700;box-shadow:inset 3px 0 0 #2f80c0}.nav-branch summary a.active{display:inline-block;border-radius:6px;padding:3px 6px}.nav-empty{font-size:12px;color:#6b7c8f;margin:8px 14px 12px}
.main{flex:1;min-width:0;max-width:1280px}.report-card,.tex-section{margin-bottom:24px;padding:22px;border:1px solid #d6dde6;border-radius:8px;background:#fff;box-shadow:0 8px 22px rgba(20,32,45,.05);scroll-margin-top:24px}.report-card h1{margin:0 0 8px;font-size:28px;line-height:1.2}.report-meta{margin:0;color:#516274;font-size:13px}.tex-section h1,.tex-section h2,.tex-section h3,.tex-section h4,.tex-section h5,.tex-section h6{line-height:1.25;margin:0 0 14px;color:#142536}.tex-section h2{font-size:24px}.tex-section h3{font-size:21px}.tex-section h4{font-size:18px}.tex-section h5,.tex-section h6{font-size:15px}p{margin:0 0 14px}.citation{color:#526477}.label-anchor{position:relative;top:-24px}.math-block{overflow:auto;margin:14px 0;padding:12px 14px;border:1px solid #e1e8f0;border-radius:8px;background:#fbfdff}.tex-figure,.table-block{margin:18px 0 24px;padding:12px;border:1px solid #d9e0e8;border-radius:8px;background:#fbfcfe;scroll-margin-top:24px}.tex-figure img{display:block;max-width:100%;height:auto;margin:0 auto;border-radius:4px}.tex-figure img[data-lightbox-image]{cursor:zoom-in}.tex-figure img[data-lightbox-image]:focus{outline:3px solid #2f80c0;outline-offset:3px}.tex-figure figcaption{font-size:13px;line-height:1.45;color:#33475c;margin-top:10px}.figure-source{margin:6px 0 0;font-size:11px}.figure-source a{color:#5f7488;text-decoration:none}.figure-source a:hover{text-decoration:underline}.float-label{font-weight:800;color:#172b3f}.missing-asset,.asset-link{font-size:13px;color:#6b3440;background:#fff4f4;border:1px solid #f1caca;padding:10px;border-radius:6px}.table-caption{font-size:13px;line-height:1.45;color:#33475c;font-weight:700;margin:0 0 8px}.table-wrap{overflow:auto;margin:12px 0 6px}.tex-table{border-collapse:collapse;width:100%;font-size:13px}.tex-table th,.tex-table td{padding:8px 10px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top}.tex-table th{background:#f7f9fb}code{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:#f1f4f7;border-radius:4px;padding:1px 4px}a{color:#1f6aa5}
.image-lightbox{position:fixed;inset:0;z-index:9999;display:none;background:rgba(8,14,22,.92);color:#fff;overflow:hidden}.image-lightbox.is-open{display:block}.image-lightbox-stage{position:absolute;inset:0;display:flex;align-items:center;justify-content:center;overflow:hidden;touch-action:none;cursor:grab}.image-lightbox-stage.is-dragging{cursor:grabbing}.image-lightbox-img{position:absolute;left:50%;top:50%;max-width:none;max-height:none;transform-origin:center center;user-select:none;-webkit-user-drag:none;will-change:transform}.image-lightbox-toolbar{position:absolute;z-index:2;top:14px;right:14px;display:flex;gap:8px}.image-lightbox-toolbar button{appearance:none;border:1px solid rgba(255,255,255,.28);border-radius:6px;background:rgba(17,28,42,.85);color:#fff;font-size:12px;font-weight:700;padding:8px 10px;cursor:pointer}.image-lightbox-toolbar button:hover{background:rgba(38,58,82,.95)}.image-lightbox-caption{position:absolute;z-index:2;left:18px;right:18px;bottom:14px;max-height:18vh;overflow:auto;padding:10px 12px;border-radius:7px;background:rgba(17,28,42,.82);font-size:13px;line-height:1.4;color:#eef4fa}
.ref-popover{position:fixed;inset:0;z-index:8000;display:none;background:rgba(8,14,22,.34);padding:6vh 24px;box-sizing:border-box}.ref-popover.is-open{display:flex;align-items:flex-start;justify-content:center}.ref-popover-panel{width:min(1040px,96vw);max-height:88vh;display:flex;flex-direction:column;border:1px solid #cbd5df;border-radius:8px;background:#fff;box-shadow:0 22px 72px rgba(10,20,32,.34);overflow:hidden}.ref-popover-head{display:flex;align-items:center;justify-content:space-between;gap:16px;padding:12px 14px;border-bottom:1px solid #dde5ed;background:#f7fafc}.ref-popover-title{font-size:14px;color:#172b3f}.ref-popover-jump{font-size:13px;font-weight:700;text-decoration:none;color:#1f6aa5;white-space:nowrap}.ref-popover-jump:hover{text-decoration:underline}.ref-popover-body{overflow:auto;padding:14px;max-height:calc(88vh - 50px)}.ref-popover-body .tex-figure,.ref-popover-body .table-block{margin:0;padding:10px;box-shadow:none}.ref-popover-body .figure-source{display:none}.ref-popover-body .table-wrap{max-height:62vh}.ref-popover-body [id]{scroll-margin-top:0}.ref-preview-zoom-tools{display:flex;justify-content:flex-end;gap:8px;margin:0 0 8px}.ref-preview-zoom-tools button{appearance:none;border:1px solid #c6d3df;border-radius:6px;background:#f7fafc;color:#24384c;font-size:12px;font-weight:700;padding:6px 9px;cursor:pointer}.ref-preview-zoom-tools button:hover{background:#eef4fa}.ref-preview-zoom-stage{position:relative;height:min(58vh,680px);min-height:280px;overflow:hidden;touch-action:none;cursor:grab;border:1px solid #d9e0e8;border-radius:7px;background:#f8fafc}.ref-preview-zoom-stage.is-dragging{cursor:grabbing}.ref-preview-zoom-img{position:absolute;left:50%;top:50%;display:block;max-width:none!important;max-height:none!important;width:auto;height:auto;transform-origin:center center;user-select:none;-webkit-user-drag:none;will-change:transform}
@media (max-width:1100px){.shell{display:block}.sidebar{position:static;width:auto;margin-bottom:20px}}@media print{.shell{display:block}.sidebar{display:none}.report-card,.tex-section{box-shadow:none}}
"""


def report_script() -> str:
    return """
(function(){
  var links = Array.prototype.slice.call(document.querySelectorAll('.sidebar a[data-target]'));
  var byId = {};
  links.forEach(function(a){ byId[a.getAttribute('data-target')] = a; });
  var targets = Array.prototype.slice.call(document.querySelectorAll('[data-nav-target]')).filter(function(el){ return !!byId[el.id]; });
  function ancestorDetails(active){
    var keep = [];
    var group = active ? active.closest('details') : null;
    while(group){ keep.push(group); group = group.parentElement ? group.parentElement.closest('details') : null; }
    return keep;
  }
  function collapseInactive(active){
    var keep = ancestorDetails(active);
    document.querySelectorAll('.toc-tree details.nav-branch').forEach(function(group){
      if(keep.indexOf(group) === -1){ group.open = false; }
    });
    keep.forEach(function(group){ group.open = true; });
  }
  function setActive(id){
    links.forEach(function(a){ a.classList.toggle('active', a.getAttribute('data-target') === id); });
    var active = byId[id];
    if(active){
      collapseInactive(active);
      var item = active.closest('.nav-leaf,details');
      if(item && item.scrollIntoView){ item.scrollIntoView({block:'nearest'}); }
    }
  }
  links.forEach(function(a){ a.addEventListener('click', function(){ setActive(a.getAttribute('data-target')); }); });
  if('IntersectionObserver' in window && targets.length){
    var visible = {};
    var observer = new IntersectionObserver(function(entries){
      entries.forEach(function(entry){
        if(entry.isIntersecting){ visible[entry.target.id] = entry.boundingClientRect.top; }
        else{ delete visible[entry.target.id]; }
      });
      var best = null;
      Object.keys(visible).forEach(function(id){
        if(best === null || Math.abs(visible[id]) < Math.abs(visible[best])){ best = id; }
      });
      if(best !== null){ setActive(best); }
    }, {rootMargin:'-18% 0px -68% 0px', threshold:[0,0.05,0.25]});
    targets.forEach(function(el){ observer.observe(el); });
  } else if(targets.length) {
    setActive(targets[0].id);
  }
})();
(function(){
  var popover = document.getElementById('ref-popover');
  if(!popover){ return; }
  var panel = popover.querySelector('.ref-popover-panel');
  var title = popover.querySelector('.ref-popover-title');
  var body = popover.querySelector('.ref-popover-body');
  var jump = popover.querySelector('.ref-popover-jump');
  var activeTargetId = null;
  function scrubCloneIds(node){
    if(node.nodeType !== 1){ return; }
    if(node.hasAttribute('id')){ node.removeAttribute('id'); }
    Array.prototype.forEach.call(node.querySelectorAll('[id]'), function(child){ child.removeAttribute('id'); });
    Array.prototype.forEach.call(node.querySelectorAll('img[data-lightbox-image]'), function(img){
      img.removeAttribute('data-lightbox-image');
      img.removeAttribute('tabindex');
    });
  }
  function setupPreviewZoom(container){
    Array.prototype.forEach.call(container.querySelectorAll('.tex-figure img'), function(img){
      var tools = document.createElement('div');
      tools.className = 'ref-preview-zoom-tools';
      tools.innerHTML = '<button type="button" data-preview-action="fit">Fit</button><button type="button" data-preview-action="actual">1:1</button>';
      var stage = document.createElement('div');
      stage.className = 'ref-preview-zoom-stage';
      img.parentNode.insertBefore(tools, img);
      img.parentNode.insertBefore(stage, img);
      stage.appendChild(img);
      img.classList.add('ref-preview-zoom-img');
      var scale = 1;
      var x = 0;
      var y = 0;
      var pointers = new Map();
      var dragging = false;
      var lastPoint = null;
      var lastPinchDistance = null;
      function clampScale(value){ return Math.max(0.05, Math.min(14, value)); }
      function apply(){ img.style.transform = 'translate(calc(-50% + ' + x + 'px), calc(-50% + ' + y + 'px)) scale(' + scale + ')'; }
      function fitScale(){
        if(!img.naturalWidth || !img.naturalHeight){ return 1; }
        return clampScale(Math.min((stage.clientWidth - 24) / img.naturalWidth, (stage.clientHeight - 24) / img.naturalHeight, 1));
      }
      function fit(){ scale = fitScale(); x = 0; y = 0; apply(); }
      function actual(){ scale = 1; x = 0; y = 0; apply(); }
      function setScaleAround(nextScale, clientX, clientY){
        nextScale = clampScale(nextScale);
        var rect = stage.getBoundingClientRect();
        var cx = clientX == null ? rect.left + rect.width / 2 : clientX;
        var cy = clientY == null ? rect.top + rect.height / 2 : clientY;
        var dx = cx - (rect.left + rect.width / 2) - x;
        var dy = cy - (rect.top + rect.height / 2) - y;
        var ratio = nextScale / scale;
        x -= dx * (ratio - 1);
        y -= dy * (ratio - 1);
        scale = nextScale;
        apply();
      }
      function pinchDistance(){
        var pts = Array.from(pointers.values());
        if(pts.length < 2){ return null; }
        return Math.hypot(pts[0].clientX - pts[1].clientX, pts[0].clientY - pts[1].clientY);
      }
      function pinchCenter(){
        var pts = Array.from(pointers.values());
        if(pts.length < 2){ return null; }
        return {clientX:(pts[0].clientX + pts[1].clientX) / 2, clientY:(pts[0].clientY + pts[1].clientY) / 2};
      }
      tools.querySelector('[data-preview-action="fit"]').addEventListener('click', fit);
      tools.querySelector('[data-preview-action="actual"]').addEventListener('click', actual);
      stage.addEventListener('wheel', function(event){
        event.preventDefault();
        var factor = Math.exp(-event.deltaY * 0.0014);
        setScaleAround(scale * factor, event.clientX, event.clientY);
      }, {passive:false});
      stage.addEventListener('dblclick', function(event){
        event.preventDefault();
        if(Math.abs(scale - fitScale()) < 0.02){ actual(); } else { fit(); }
      });
      stage.addEventListener('pointerdown', function(event){
        event.preventDefault();
        stage.setPointerCapture(event.pointerId);
        pointers.set(event.pointerId, {clientX:event.clientX, clientY:event.clientY});
        if(pointers.size === 1){
          dragging = true;
          lastPoint = {clientX:event.clientX, clientY:event.clientY};
          stage.classList.add('is-dragging');
        } else if(pointers.size === 2) {
          dragging = false;
          stage.classList.remove('is-dragging');
          lastPinchDistance = pinchDistance();
        }
      });
      stage.addEventListener('pointermove', function(event){
        if(!pointers.has(event.pointerId)){ return; }
        event.preventDefault();
        pointers.set(event.pointerId, {clientX:event.clientX, clientY:event.clientY});
        if(pointers.size >= 2){
          var dist = pinchDistance();
          var center = pinchCenter();
          if(dist && lastPinchDistance && center){ setScaleAround(scale * (dist / lastPinchDistance), center.clientX, center.clientY); }
          lastPinchDistance = dist;
          return;
        }
        if(dragging && lastPoint){
          x += event.clientX - lastPoint.clientX;
          y += event.clientY - lastPoint.clientY;
          lastPoint = {clientX:event.clientX, clientY:event.clientY};
          apply();
        }
      });
      function releasePointer(event){
        if(pointers.has(event.pointerId)){ pointers.delete(event.pointerId); }
        if(pointers.size < 2){ lastPinchDistance = null; }
        if(pointers.size === 0){ dragging = false; lastPoint = null; stage.classList.remove('is-dragging'); }
      }
      stage.addEventListener('pointerup', releasePointer);
      stage.addEventListener('pointercancel', releasePointer);
      if(img.complete){ fit(); } else { img.addEventListener('load', fit, {once:true}); }
    });
  }
  function closePreview(){
    popover.classList.remove('is-open');
    popover.setAttribute('aria-hidden', 'true');
    body.innerHTML = '';
    activeTargetId = null;
  }
  function openPreview(link){
    var targetId = link.getAttribute('data-ref-preview-target');
    var target = targetId ? document.getElementById(targetId) : null;
    if(!target){ return; }
    activeTargetId = targetId;
    title.textContent = link.getAttribute('data-ref-preview-label') || link.textContent.trim() || 'Reference preview';
    jump.setAttribute('href', '#' + targetId);
    body.innerHTML = '';
    var clone = target.cloneNode(true);
    scrubCloneIds(clone);
    body.appendChild(clone);
    popover.classList.add('is-open');
    popover.setAttribute('aria-hidden', 'false');
    setupPreviewZoom(body);
  }
  document.querySelectorAll('a[data-ref-preview-target]').forEach(function(link){
    link.addEventListener('click', function(event){
      event.preventDefault();
      event.stopPropagation();
      openPreview(link);
    });
  });
  jump.addEventListener('click', function(event){
    if(!activeTargetId){ return; }
    event.preventDefault();
    var target = document.getElementById(activeTargetId);
    closePreview();
    if(target && target.scrollIntoView){ target.scrollIntoView({behavior:'smooth', block:'start'}); }
    if(history.pushState){ history.pushState(null, '', '#' + activeTargetId); }
    else { window.location.hash = activeTargetId; }
  });
  popover.addEventListener('click', function(event){
    if(event.target === popover){ closePreview(); }
  });
  panel.addEventListener('click', function(event){
    event.stopPropagation();
  });
  document.addEventListener('keydown', function(event){
    if(popover.classList.contains('is-open') && event.key === 'Escape'){ closePreview(); }
  });
})();
(function(){
  function textContentFromFigure(link){
    var figure = link.closest('figure');
    var caption = figure ? figure.querySelector('figcaption') : null;
    return caption ? caption.textContent.trim() : 'Source image';
  }
  function sourceImageFromFigure(link){
    var figure = link.closest('figure');
    var img = figure ? figure.querySelector('img') : null;
    return img ? (img.currentSrc || img.src) : link.href;
  }
  function writeViewer(win, src, titleText){
    win.document.open();
    win.document.write('<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title></title><style>html,body{margin:0;height:100%;overflow:hidden;background:#08101a;color:#fff;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Arial,sans-serif}.stage{position:fixed;inset:0;overflow:hidden;touch-action:none;cursor:grab}.stage.dragging{cursor:grabbing}.viewer{position:absolute;left:50%;top:50%;max-width:none;max-height:none;transform-origin:center center;user-select:none;-webkit-user-drag:none;will-change:transform}.toolbar{position:fixed;z-index:2;top:14px;right:14px;display:flex;gap:8px}.toolbar button{appearance:none;border:1px solid rgba(255,255,255,.28);border-radius:6px;background:rgba(17,28,42,.88);color:#fff;font-size:12px;font-weight:700;padding:8px 10px;cursor:pointer}.caption{position:fixed;z-index:2;left:14px;right:14px;bottom:14px;max-height:18vh;overflow:auto;padding:10px 12px;border-radius:7px;background:rgba(17,28,42,.82);font-size:13px;line-height:1.4;color:#eef4fa}</style></head><body><div class="toolbar"><button type="button" data-action="fit">Fit</button><button type="button" data-action="actual">1:1</button></div><div class="stage"><img class="viewer" alt=""></div><div class="caption"></div></body></html>');
    win.document.close();
    win.document.title = titleText || 'Source image';
    var doc = win.document;
    var stage = doc.querySelector('.stage');
    var viewer = doc.querySelector('.viewer');
    var caption = doc.querySelector('.caption');
    var scale = 1;
    var x = 0;
    var y = 0;
    var pointers = new Map();
    var dragging = false;
    var lastPoint = null;
    var lastPinchDistance = null;
    function clampScale(value){ return Math.max(0.04, Math.min(16, value)); }
    function apply(){ viewer.style.transform = 'translate(calc(-50% + ' + x + 'px), calc(-50% + ' + y + 'px)) scale(' + scale + ')'; }
    function fitScale(){
      if(!viewer.naturalWidth || !viewer.naturalHeight){ return 1; }
      return clampScale(Math.min((win.innerWidth - 72) / viewer.naturalWidth, (win.innerHeight - 118) / viewer.naturalHeight, 1));
    }
    function fit(){ scale = fitScale(); x = 0; y = 0; apply(); }
    function actual(){ scale = 1; x = 0; y = 0; apply(); }
    function setScaleAround(nextScale, clientX, clientY){
      nextScale = clampScale(nextScale);
      var rect = stage.getBoundingClientRect();
      var cx = clientX == null ? rect.left + rect.width / 2 : clientX;
      var cy = clientY == null ? rect.top + rect.height / 2 : clientY;
      var dx = cx - (rect.left + rect.width / 2) - x;
      var dy = cy - (rect.top + rect.height / 2) - y;
      var ratio = nextScale / scale;
      x -= dx * (ratio - 1);
      y -= dy * (ratio - 1);
      scale = nextScale;
      apply();
    }
    function pinchDistance(){
      var pts = Array.from(pointers.values());
      if(pts.length < 2){ return null; }
      return Math.hypot(pts[0].clientX - pts[1].clientX, pts[0].clientY - pts[1].clientY);
    }
    function pinchCenter(){
      var pts = Array.from(pointers.values());
      if(pts.length < 2){ return null; }
      return {clientX:(pts[0].clientX + pts[1].clientX) / 2, clientY:(pts[0].clientY + pts[1].clientY) / 2};
    }
    caption.textContent = titleText || '';
    viewer.addEventListener('load', fit);
    viewer.src = src;
    viewer.alt = titleText || '';
    doc.querySelector('[data-action="fit"]').addEventListener('click', fit);
    doc.querySelector('[data-action="actual"]').addEventListener('click', actual);
    stage.addEventListener('wheel', function(event){
      event.preventDefault();
      var factor = Math.exp(-event.deltaY * 0.0014);
      setScaleAround(scale * factor, event.clientX, event.clientY);
    }, {passive:false});
    stage.addEventListener('dblclick', function(event){
      event.preventDefault();
      if(Math.abs(scale - fitScale()) < 0.02){ actual(); } else { fit(); }
    });
    stage.addEventListener('pointerdown', function(event){
      event.preventDefault();
      stage.setPointerCapture(event.pointerId);
      pointers.set(event.pointerId, {clientX:event.clientX, clientY:event.clientY});
      if(pointers.size === 1){
        dragging = true;
        lastPoint = {clientX:event.clientX, clientY:event.clientY};
        stage.classList.add('dragging');
      } else if(pointers.size === 2) {
        dragging = false;
        stage.classList.remove('dragging');
        lastPinchDistance = pinchDistance();
      }
    });
    stage.addEventListener('pointermove', function(event){
      if(!pointers.has(event.pointerId)){ return; }
      event.preventDefault();
      pointers.set(event.pointerId, {clientX:event.clientX, clientY:event.clientY});
      if(pointers.size >= 2){
        var dist = pinchDistance();
        var center = pinchCenter();
        if(dist && lastPinchDistance && center){ setScaleAround(scale * (dist / lastPinchDistance), center.clientX, center.clientY); }
        lastPinchDistance = dist;
        return;
      }
      if(dragging && lastPoint){
        x += event.clientX - lastPoint.clientX;
        y += event.clientY - lastPoint.clientY;
        lastPoint = {clientX:event.clientX, clientY:event.clientY};
        apply();
      }
    });
    function releasePointer(event){
      if(pointers.has(event.pointerId)){ pointers.delete(event.pointerId); }
      if(pointers.size < 2){ lastPinchDistance = null; }
      if(pointers.size === 0){ dragging = false; lastPoint = null; stage.classList.remove('dragging'); }
    }
    stage.addEventListener('pointerup', releasePointer);
    stage.addEventListener('pointercancel', releasePointer);
    win.addEventListener('resize', fit);
    win.addEventListener('keydown', function(event){
      if(event.key === '0'){ fit(); }
      if(event.key === '1'){ actual(); }
    });
  }
  document.querySelectorAll('a[data-source-viewer]').forEach(function(link){
    link.addEventListener('click', function(event){
      var win = window.open('', '_blank');
      if(!win){ return; }
      event.preventDefault();
      writeViewer(win, sourceImageFromFigure(link), textContentFromFigure(link));
    });
  });
})();
(function(){
  var overlay = document.getElementById('image-lightbox');
  if(!overlay){ return; }
  var stage = overlay.querySelector('.image-lightbox-stage');
  var viewer = overlay.querySelector('.image-lightbox-img');
  var caption = overlay.querySelector('.image-lightbox-caption');
  var scale = 1;
  var x = 0;
  var y = 0;
  var pointers = new Map();
  var lastPinchDistance = null;
  var dragging = false;
  var lastPoint = null;
  var pointerMoved = false;
  var backgroundPointerDown = false;
  var suppressNextStageClick = false;
  function clampScale(value){
    return Math.max(0.08, Math.min(12, value));
  }
  function applyTransform(){
    viewer.style.transform = 'translate(calc(-50% + ' + x + 'px), calc(-50% + ' + y + 'px)) scale(' + scale + ')';
  }
  function fitScale(){
    if(!viewer.naturalWidth || !viewer.naturalHeight){ return 1; }
    var maxW = Math.max(120, window.innerWidth - 72);
    var maxH = Math.max(120, window.innerHeight - 120);
    return clampScale(Math.min(maxW / viewer.naturalWidth, maxH / viewer.naturalHeight, 1));
  }
  function setScaleAround(nextScale, clientX, clientY){
    nextScale = clampScale(nextScale);
    var rect = stage.getBoundingClientRect();
    var cx = clientX == null ? rect.left + rect.width / 2 : clientX;
    var cy = clientY == null ? rect.top + rect.height / 2 : clientY;
    var dx = cx - (rect.left + rect.width / 2) - x;
    var dy = cy - (rect.top + rect.height / 2) - y;
    var ratio = nextScale / scale;
    x -= dx * (ratio - 1);
    y -= dy * (ratio - 1);
    scale = nextScale;
    applyTransform();
  }
  function fit(){
    scale = fitScale();
    x = 0;
    y = 0;
    applyTransform();
  }
  function actual(){
    scale = 1;
    x = 0;
    y = 0;
    applyTransform();
  }
  function openFrom(img){
    viewer.src = img.currentSrc || img.src;
    viewer.alt = img.alt || '';
    var figcaption = img.closest('figure') ? img.closest('figure').querySelector('figcaption') : null;
    caption.textContent = figcaption ? figcaption.textContent.trim() : (img.alt || '');
    overlay.classList.add('is-open');
    overlay.setAttribute('aria-hidden', 'false');
    document.body.style.overflow = 'hidden';
    pointers.clear();
    lastPinchDistance = null;
    dragging = false;
    pointerMoved = false;
    backgroundPointerDown = false;
    suppressNextStageClick = false;
    stage.classList.remove('is-dragging');
    if(viewer.complete){ fit(); }
  }
  function close(){
    overlay.classList.remove('is-open');
    overlay.setAttribute('aria-hidden', 'true');
    document.body.style.overflow = '';
    viewer.removeAttribute('src');
    pointers.clear();
    lastPinchDistance = null;
    dragging = false;
    pointerMoved = false;
    backgroundPointerDown = false;
    suppressNextStageClick = false;
    stage.classList.remove('is-dragging');
  }
  function pinchDistance(){
    var pts = Array.from(pointers.values());
    if(pts.length < 2){ return null; }
    var dx = pts[0].clientX - pts[1].clientX;
    var dy = pts[0].clientY - pts[1].clientY;
    return Math.hypot(dx, dy);
  }
  function pinchCenter(){
    var pts = Array.from(pointers.values());
    if(pts.length < 2){ return null; }
    return {
      clientX: (pts[0].clientX + pts[1].clientX) / 2,
      clientY: (pts[0].clientY + pts[1].clientY) / 2
    };
  }
  document.querySelectorAll('img[data-lightbox-image]').forEach(function(img){
    img.addEventListener('click', function(event){
      event.preventDefault();
      openFrom(img);
    });
    img.addEventListener('keydown', function(event){
      if(event.key === 'Enter' || event.key === ' '){
        event.preventDefault();
        openFrom(img);
      }
    });
  });
  viewer.addEventListener('load', fit);
  overlay.addEventListener('click', function(event){
    if(suppressNextStageClick){
      suppressNextStageClick = false;
      return;
    }
    if((event.target === overlay || event.target === stage) && backgroundPointerDown && !pointerMoved){ close(); }
    backgroundPointerDown = false;
    pointerMoved = false;
  });
  overlay.querySelectorAll('[data-lightbox-action]').forEach(function(button){
    button.addEventListener('click', function(event){
      event.stopPropagation();
      var action = button.getAttribute('data-lightbox-action');
      if(action === 'fit'){ fit(); }
      if(action === 'actual'){ actual(); }
      if(action === 'close'){ close(); }
    });
  });
  stage.addEventListener('wheel', function(event){
    if(!overlay.classList.contains('is-open')){ return; }
    event.preventDefault();
    var factor = Math.exp(-event.deltaY * 0.0014);
    setScaleAround(scale * factor, event.clientX, event.clientY);
  }, {passive:false});
  stage.addEventListener('dblclick', function(event){
    event.preventDefault();
    if(Math.abs(scale - fitScale()) < 0.02){ actual(); } else { fit(); }
  });
  stage.addEventListener('pointerdown', function(event){
    if(event.target !== viewer && event.target !== stage){ return; }
    event.preventDefault();
    stage.setPointerCapture(event.pointerId);
    pointers.set(event.pointerId, {clientX:event.clientX, clientY:event.clientY});
    pointerMoved = false;
    backgroundPointerDown = event.target === stage;
    if(pointers.size === 1){
      dragging = event.target === viewer;
      lastPoint = {clientX:event.clientX, clientY:event.clientY};
      if(dragging){ stage.classList.add('is-dragging'); }
    } else if(pointers.size === 2) {
      dragging = false;
      stage.classList.remove('is-dragging');
      lastPinchDistance = pinchDistance();
    }
  });
  stage.addEventListener('pointermove', function(event){
    if(!pointers.has(event.pointerId)){ return; }
    event.preventDefault();
    var previous = pointers.get(event.pointerId);
    if(previous && Math.hypot(event.clientX - previous.clientX, event.clientY - previous.clientY) > 3){
      pointerMoved = true;
    }
    pointers.set(event.pointerId, {clientX:event.clientX, clientY:event.clientY});
    if(pointers.size >= 2){
      var dist = pinchDistance();
      var center = pinchCenter();
      if(dist && lastPinchDistance && center){
        setScaleAround(scale * (dist / lastPinchDistance), center.clientX, center.clientY);
      }
      lastPinchDistance = dist;
      return;
    }
    if(dragging && lastPoint){
      x += event.clientX - lastPoint.clientX;
      y += event.clientY - lastPoint.clientY;
      lastPoint = {clientX:event.clientX, clientY:event.clientY};
      applyTransform();
    }
  });
  function releasePointer(event){
    if(pointers.has(event.pointerId)){ pointers.delete(event.pointerId); }
    if(pointers.size < 2){ lastPinchDistance = null; }
    if(pointers.size === 0){
      if(pointerMoved || !backgroundPointerDown){ suppressNextStageClick = true; }
      dragging = false;
      lastPoint = null;
      stage.classList.remove('is-dragging');
    }
  }
  stage.addEventListener('pointerup', releasePointer);
  stage.addEventListener('pointercancel', releasePointer);
  window.addEventListener('resize', function(){
    if(overlay.classList.contains('is-open')){ fit(); }
  });
  document.addEventListener('keydown', function(event){
    if(!overlay.classList.contains('is-open')){ return; }
    if(event.key === 'Escape'){ close(); }
    if(event.key === '0'){ fit(); }
    if(event.key === '1'){ actual(); }
  });
})();
"""


def main() -> int:
    args = parse_args()
    input_path = Path(args.input).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()
    if not input_path.exists():
        print(f"Missing input file: {input_path}", file=sys.stderr)
        return 2
    output_path.parent.mkdir(parents=True, exist_ok=True)
    raw_text = input_path.read_text(encoding="utf-8")
    html_text = render_document(
        raw_text=raw_text,
        input_path=input_path,
        output_path=output_path,
        title_override=args.title,
        asset_mode=args.asset_mode,
        mathjax_url=args.mathjax_url,
        hide_source=args.hide_source,
        figure_placement=args.figure_placement,
    )
    output_path.write_text(html_text, encoding="utf-8")
    print(f"Wrote HTML report: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
