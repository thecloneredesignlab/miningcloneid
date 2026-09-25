# Stage 0 regression baseline

Baseline date: 2026-07-16
Branch: `soft_coupling`
Commit: `873c929d38c0d7f0c8e37724fed95b0d52b349d9`

The baseline was produced before moving any in-vitro helper. Scientific
regression outputs were written under `/tmp`, not into `oxygen/results`.

## Static and unit baseline

- R parse: 117/117 files passed.
- Shell/Slurm syntax: 24/24 files passed `bash -n`.
- Unit suite: 57 tests, 187 expectations, zero failures, warnings, or skips.
- The full suite ran against an isolated HEAD archive so its Rcpp cache was
  created under `/tmp`, not in the protected working-tree model directory.

## In-vitro numerical baseline

- Public `ivt_*` functions: 58, including exact formals.
- Optimizer parameters: 20.
- `fit_data` entries: 131.
- 2N jobs: 38.
- 4N jobs: 37.
- seed10 objective: `3.8525352626059366`.
- total log-likelihood: `-3.8525352626059366`.
- growth log-likelihood: `0.17059602541136337`.
- ploidy log-likelihood: `-3.0962217415551456`.
- flow log-likelihood: `-0.92690954646215451`.

Derived-table row counts:

| Table | Rows |
|---|---:|
| lineage summary | 131 |
| growth likelihood | 114 |
| ploidy likelihood | 12 |
| flow likelihood | 20 |
| flow overlay | 8000 |
| distribution summary | 9975 |
| distribution quantiles | 3750 |
| daily counts | 462 |
| observed karyotype | 220 |
| observed flow | 4000 |

The Stage 1 comparison tolerance is `1e-12` for the objective and deterministic
numeric tables.

## In-vitro visualization baseline

Using the existing seed10 table fixture and a new empty `/tmp` output directory:

- exit status: 0
- outputs: 8 PDF, 8 PNG, and one manifest
- successful panels: selected-O2 live, daily counts, lineage growth, lineage
  ploidy, burden decomposition, growth/ploidy/burden composite, flow density,
  and distribution heatmap
- unavailable in this fixture: identifiability and functional-response panels

The comparison ignores the dynamic `generated_at` value and normalizes the
recorded fit directory. Figure files are checked by name, page/dimensions, and
rendered pixels rather than PDF metadata bytes.

## Protected-tree baseline

The five tracked protected files have aggregate SHA-256:

```text
baba702decd469c42e1412d22f1d549b54190905eb18fcf938fc1af5bf384e81
```

The full protected tree contains 18 regular files, including the existing
ignored Rcpp cache and `.DS_Store` files. Its baseline aggregate is:

```text
0cd49a43164f1ff6771a3a4bea281820292a69c6c3131879ee33c4a14c644d59
```

Every stage compares SHA-256, size, and high-resolution modification time for
all 18 files by running:

```bash
Rscript oxygen/tests/check_immutable_o2_core.R --full
```
