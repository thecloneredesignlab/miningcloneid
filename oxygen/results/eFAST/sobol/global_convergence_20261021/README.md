# Global Sobol continuation approved 2026-10-21

Only the separate in-vivo global fixed-oxygen model is evaluated. Neighborhood
work remains archived. The original 14 factors, 91 pairs, 201 oxygen points,
five scramble seeds, bounds, transforms, model and SALib 1.5.2 are retained.

## Certified reuse

- N65536 appends base samples 32769–65536 to certified N32768 outputs.
- N131072 appends 65537–131072 to certified N65536 outputs.
- Preparation proves all five prefixes bit-for-bit, audits 1,200 existing
  model blocks and 130 N32768 statistical shards, and records solver flags,
  minimum spectral gaps and maximum residuals. Old model outputs are read-only.
- The model smoke task certifies eight new oxygen slices into the production
  checkpoint. Production resumes that checkpoint without recalculating them.
- Model tasks checkpoint every eight oxygen slices; statistics checkpoint each
  repeat/oxygen/output cell. Completed shards are reused on recovery.
- Original and appended outputs are interleaved in base-sample/block order.
  No model calls are made during statistical processing or plotting.

## Submission and outputs

The current branch is committed/pushed and synchronized before submission:

```sh
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/submit_global_convergence.py --revision COMMITTED_SHA
```

Eight parent stages: prepare, smoke, model/analyze/report N65536, then
model/analyze/report N131072. There are 1,200 + 2,400 model tasks and 130 + 130
statistical tasks. Every scientific task uses xxlarge, an initial 12 hours,
one CPU, no requested node and no array concurrency cap. Every parent has an
independent monitor: add six hours when fewer than six hours remain, checking
every 300 seconds until the target leaves the queue. Monitor exits alone do
not establish scientific success. Submission is idempotent; resume the ledger
instead of submitting a duplicate stage.

`slurm/jobs.json`, `launch_receipt.json`, and
`time_extension/submission_jobs.tsv` identify exact submitted jobs. New model
tails live in `increments/`; `N65536/` and `N131072/` contain indices, source
tables, annotations, row orders, figures, and resolution comparisons. The
existing N4096–N32768 results remain unchanged.

## Diagnostics and plots

Both levels run in full. Convergence is diagnostic only: five-repeat range
<=0.05 and adjacent-resolution mean difference <=0.02. Undefined and
nonpassing cells remain in all source tables. Reports retain bootstrap half
widths and output variances and compare newly passed/lost/still failing cells.
ST curves cover all 14 factors; focal S2 covers the 36 pairs involving
p_mis_base, mu_hp, or gamma_mu. Other S2 pairs remain in the full heatmap/table.

S1 and ST use the same white-to-purple palette with separate output scales.
Oxygen uses real cell coordinates, symlog with a 0–0.025% linear segment to
retain zero, and dashed 1%/3% boundaries. Rows use ploidy scores: High>=Low,
High descending, followed by Low>High, Low ascending; growth follows that
ploidy row order. S1/ST have process, original Figure 4B correlation and current
Sobol O2-preference annotations. S2 has only the first two, split by pair member.
The third annotation uses whole-repeat bootstrap curves, no 0.3 threshold,
and never determines row order. It is exploratory, not calibrated significance
or a numerical convergence test. Negative S1/S2 estimates are sampling noise.

The approved upper bound for this continuation is N131072. A further increase
is not submitted automatically, and passing every diagnostic is not promised.
