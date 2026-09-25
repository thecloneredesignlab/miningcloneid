# Deprecated 06 combine entrypoints

The three R files in this directory are thin compatibility entrypoints. Canonical code now lives in:

- `analysis/combined_fixo2_eigen/` for slope and annotated classification tables;
- `vis/combined_fixo2_eigen/` for figures;
- `report/combined_fixo2_eigen/` for HTML assembly;
- `runner/combined_fixo2_eigen/` for ordered execution.

Existing result paths and command-line options remain supported. New workflows should call the canonical runner.
