# Multi-warmup report

## File registry

### `render_multi_warmup_results_report.R`

- Reads existing manifests, summaries, figures, integrated extra-results HTML,
  and provenance files.
- Assembles `multi-warm-up_results.html` without computing analysis or figures.
- Called by `runner/multi_warmup/render_multi_warmup_results.R`.
