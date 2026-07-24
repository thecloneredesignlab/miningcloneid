# Baseline Manuscript Build

## Authoritative Input

- Source: `ltee_hypoxia_model.tex`
- Lines: 954
- Bytes: 80,037
- SHA-256: `2fcc918e6d095a5b2962eac9bc8d2b581ba0d15a724fb6f975c61bd21e3ff3b0`
- Git state: ignored and untracked because `.gitignore` contains two global
  `*.tex` rules

## Toolchain

- `latexmk` 4.77
- TeX Live 2022 pdfTeX, XeTeX, and LuaHBTeX
- BibTeX 0.99d
- Biber 2.17

## Baseline Result

Command:

```bash
latexmk -pdf -g -interaction=nonstopmode -halt-on-error \
  -file-line-error ltee_hypoxia_model.tex
```

Result: exit 12; no PDF. The first fatal error is the missing Figure 1 image at
`ltee_hypoxia_model.tex:138`.

A draft-graphics diagnostic build reached page 15 before stopping on unescaped
ampersands at lines 235 and 237.

## Missing Manuscript Inputs

| TeX line | Missing path |
|---:|---|
| 138 | `oxygen/figures/assembled_fig1.png` |
| 157 | `oxygen/figures/assembled_fig3.png` |
| 169 | `oxygen/figures/assembled_fig4.png` |
| 179 | `oxygen/figures/assembled_fig5.png` |
| 191 | `oxygen/figures/assembled_fig6.png` |
| 257 | `oxygen/figures/FixedO2_vs_ploidy.png` |
| 944 | `figures/tao_model_parameters.tex` |
| 947 | `figures/tao_fixed_parameters.tex` |
| 951 | `references_Zotero_TaoLi.bib` |

Figure 2 has no `\includegraphics` call. The manuscript also contains 13 `XX`
placeholders, three explicit panel-to-be-generated statements, an empty
supplementary in-vitro experiment section, and an author-facing Figure 6B design
instruction.

These non-figure dependencies do not block ideation. They do block a
manuscript-ready compiled handoff.
