# Documentation runner

| File | Responsibility | Inputs | Outputs |
|---|---|---|---|
| `generate_code_file_registry.R` | Inspect every source file in the workflow and regenerate the per-file documentation registry. | `--workflow_root` and optional output paths. | `docs/CODE_FILE_REGISTRY.md` and `docs/code_file_registry.tsv`. |

Run it after structural changes and before the full regression gate.
