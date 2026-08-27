# In-vitro structure fix: two-seed audit

## Scope

This directory contains only the audit code, two-seed parameter-replay results,
and assessment for the in-vitro structural correction. The model
implementation remains in:

`/Users/4482173/Documents/GitHub/soft_coupling/oxygen/code/O2_supply_demand_MAP`

No manuscript file, existing 500-seed result, git commit, remote branch, or HPC
job was changed or created.

The requested read-only HPC baseline was not mounted on this machine:

`/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed`

The validation therefore does not substitute or compare a similarly named local
result tree as though it were that baseline. It performs two explicit parameter
replays only, reading the parameter vectors from:

- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed/seed10/best_params.tsv`
- `/Users/4482173/Documents/GitHub/soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed/seed340/best_params.tsv`

No optimization is run.

## Validation command

The original completed run used:

```bash
Rscript /Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/add/invitro_structure_fix/code/run_invitro_structure_fix_validation.R
```

The script writes through a process-specific temporary directory and renames it
to `results/` only after all assertions pass. It refuses to run if `results/`
already exists, so rerunning cannot overwrite this audit. The review-fix replay
used a new directory and left the original results untouched:

```bash
INVITRO_STRUCTURE_FIX_RESULTS_DIR=results_review_fixes \
  Rscript /Users/4482173/Documents/GitHub/HypoxiaLTEEFigures/revised/add/invitro_structure_fix/code/run_invitro_structure_fix_validation.R
```

`results_review_fixes/` is the current result set. It includes the corrected
source-state assignment for non-root karyotype landmarks.

## Result map

- `results_review_fixes/validation_overview.tsv`: objective components, row
  counts, boundary mode counts, and top-level pass/fail flags for each seed.
- `results_review_fixes/scenario_cumulative_time.tsv`: recorded, endpoint,
  selected-alias, diagnostic-closest, and expected cumulative days for all six
  scenarios.
- `results_review_fixes/likelihood_uniqueness.tsv`: passage-ID uniqueness for growth,
  karyotype, and flow likelihoods.
- `results_review_fixes/shared_parameter_summary.tsv`: confirms that all six
  scenarios use the same replayed parameter value for every parameter.
- `results_review_fixes/input_and_code_provenance.tsv`: exact input/code paths
  and MD5 checksums used by the replay.
- `results_review_fixes/run_context.tsv`: branch, HEAD, baseline accessibility,
  numerical configuration, and non-overwrite policy.
- `results_review_fixes/seed*/invitro_passage_audit.tsv`: requested per-passage
  audit fields.
- `results_review_fixes/seed*/passage_boundary_checks.tsv`: boundary mode and
  full-state rule checks.
- `results_review_fixes/seed*/passage_boundary_state_vectors.tsv`: endpoint and
  reseeded values for every one of the 133 chromosome-state components at every
  passage.
- `results_review_fixes/seed*/invitro_*_loglik.tsv` and
  `invitro_objective_hierarchy.tsv`: likelihood units and their
  unit-to-lineage-to-cohort aggregation.

## Replay summary

| Seed | Objective | Scenarios | Passage rows | Insufficient boundaries | Downsample boundaries | Terminal passages |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 4.00749381259844 | 6 | 114 | 0 | 108 | 6 |
| 340 | 5.88644083078414 | 6 | 114 | 33 | 75 | 6 |

For both seeds:

- cumulative endpoint time is exactly 29, 29, 122, 122, 125, and 121 days for
  2N-C, 4N-C, 2N-O1, 2N-O2, 4N-O1, and 4N-O2, respectively;
- every `selected_day` compatibility value equals the fixed endpoint;
- all boundary scales are at most 1;
- every sufficient boundary is an exact proportional downsample;
- every insufficient seed340 boundary carries the complete state vector
  unchanged (maximum absolute component difference 0);
- growth, karyotype, and flow have 114/12/20 unique likelihood passage IDs and
  zero duplicates;
- `SUM-159_NLS_2N_A6M_seed` and `SUM-159_NLS_4N_A4M_seed` are evaluated at
  `2N-C-A1` and `4N-C-A2`, while only the true depth-zero anchors remain in
  `INITIAL`;
- all scenarios use the shared parameter vector.

For seed340, the diagnostic closest-day totals for 2N-O1 and 2N-O2 are each 100
days, but those values are not used for growth, likelihood, endpoint extraction,
or state propagation. Both model trajectories advance through the full 122
recorded days.

## Full refit boundary

These results validate structure and deterministic execution with two existing
parameter vectors; they are not replacement fitted parameters. A future
500-seed refit, if separately authorized, should write to a new root such as:

`/share/lab_crd/lab_crd/taoli/Project/miningcloneid_soft_coupling/oxygen/results/fit_invitro_O2_buffering_500seed_time_lineage_fix`

No such run was submitted here.
