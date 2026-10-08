# Per-endpoint in-vivo eFAST neighborhoods

This folder contains a separate sensitivity analysis around each of the 500
independent **in-vivo** optimizer endpoints. It does not use in-vitro or joint-fit bounds.

For parameter j and endpoint s, with the original natural bounds L and U:

```
lower(s,j) = max(L(j), best(s,j) - 0.10 * (U(j) - L(j)))
upper(s,j) = min(U(j), best(s,j) + 0.10 * (U(j) - L(j)))
```

Thus the untruncated total width is 20% of the original natural range.
Log10 parameters are log-uniform within these natural bounds; other parameters
are uniform. Parameters are varied independently and no local refitting is done.
Fourteen parameters affect the fixed-oxygen operator. The oxygen grid is the
same 201 points, 0–5% by 0.025%, as Figure 4. Each endpoint uses N=513, M=4
and **five independent phase repetitions** if the pilot convergence gate passes.

## Execution and gates

`run_neighborhood.sh pipeline` audits 500 endpoints and validates their original
operators at 0, 2.5 and 5% oxygen. It runs three representative endpoints at
N=129 and four oxygen points, five phases each, and records runtime projections.
It then extends the global designs to five phases and checks ten representative
endpoints at N=129,257,513 on the full oxygen grid, five phases each.
The pilot covers the best/worst fitting endpoints, clipped neighborhoods,
seed25, and farthest points in normalized transformed parameter space.

Full execution requires every pilot output/index/endpoint to have all cells
defined, a p90 five-phase index range <=0.10, and a p90 absolute change from
the N=257 mean to N=513 mean <=0.05. These are declared practical tolerances
in variance-fraction units, not mathematical guarantees. A failed gate saves
`pilot/convergence_gate.json` with `needs_review`; it does not launch all 500.
Increasing N changes computational cost and requires reviewing that diagnostic.

No Slurm submission is used. The runner enforces hpctpa3pc0028 and the validated
SALib 1.5.2 SIF checksum. Only one neighborhood controller may run per root;
`EFAST_WORKERS` controls fork workers and defaults to 16.

## Outputs and interpretation

- `seed_manifest.tsv`, `neighborhood_bounds.tsv`, `operator_input_audit.json`,
  and `input_manifest.json`: all endpoints, clipping, references and input hashes.
- `bounds/seed*.tsv`: per-endpoint bounds on HPC.
- `pilot/runtime_projection.json`, `pilot/convergence_diagnostics.tsv`,
  `pilot/convergence_gate.json`: measured costs and convergence checks.
- `runs/seed*/runs/N513_R*/`: ordered samples, raw outputs, indices and completion
  receipts on HPC. A complete design is reused only after its hashes match.
  Partial designs are recomputed; completed designs and completed endpoint
  summaries are retained. Evaluation is streamed in chunks to bound memory.
- `summaries/phase_indices_*.npz`: exact phase indices in 100-endpoint shards with labeled axes
  `(fit_seed, phase, index, parameter, O2, output)`; five phases stay together.
- `summaries/seed_mean_indices_*.tsv.gz`: per-endpoint five-phase means, SD,
  ranges, and valid-phase counts. No pooling of raw outputs across endpoints.
- `convergence.tsv`: across-endpoint mean, median, Q25, Q75 and valid-seed counts.
- `summaries/seed_parameter_bands.tsv.gz`, `rank_consistency.tsv`,
  `mechanism_band_summary.tsv`, `mechanism_support.tsv`: window scores, ranks,
  group means and fraction of evaluable endpoints supporting each mechanism claim.
- `figures/efast_four_panel.pdf/png`: medians across endpoints, separate S1/ST
  panels, separate ploidy/growth scales, three row annotations. The original
  Figure 4B association group remains; the new sensitivity group follows ploidy
  S1 or ST. Complete endpoint curves are bootstrap units after phase averaging.
- `raw_output_inventory.tsv`, `collection_manifest.tsv`, `interpretation.md`:
  raw-output hashes, collected artifact hashes and interpretation.

Raw neighborhood samples/outputs and per-run caches stay on HPC; compact source
tables, labeled index arrays, figures and provenance are tracked in Git.
Projected raw compressed outputs are about 90.7 GiB before pilot data. The full
design requires 3,608,955,000 operator evaluations. Runtime projections are
empirical estimates and can change with node load and parameter values.

Indices from numerically constant trajectories are undefined and recorded as
NaN with variance and valid-count diagnostics. No failed sample row is removed.
Optimizer endpoints are not independent posterior draws. Bootstrap groups and
support fractions describe these fitted endpoints, not posterior probabilities.
Total-effect group means overlap through interactions and are not additive
group variance fractions. The existing Spearman panel retains direction.
