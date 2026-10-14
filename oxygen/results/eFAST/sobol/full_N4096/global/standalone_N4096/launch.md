# Independent global Sobol post-processing

Global job: **20862092**. Submitted 2026-10-14; started at 12:06 EDT.

The original 150 global model blocks are complete. This job computes S1/ST
for 14 parameters and S2 for all 91 pairs, for both dominant mean ploidy and
net live-cell growth, at N4096 with five independent repeats and all 201
fixed oxygen conditions. The existing N2048 prefix supplies resolution
diagnostics. Convergence remains diagnostic only.

Model source code: 471f6a84afbf2363078a9cefce12a75cd0ee691f.
Frozen post-processing code: 01e7ffc4e9b1dec7921e039f75ac9bcd5928f0ef.
SIF: o2_supply_demand_map_r442_hpc_exact_salib152_20260923.sif.
Resources: xxlarge, 12h, 1 CPU, 16G; no chosen node; no neighborhood dependency.
Eighteen local tests and three tests in the actual job SIF passed.

The user approved temporary reprioritization of pending neighborhood workers
to let this one-CPU job start. All original Nice values were restored as soon
as the global job obtained its allocation, without interrupting running model
workers. The exact affected task IDs and restoration results are in
slurm/priority_adjustment.json. Existing model code snapshots, model inputs,
array concurrency settings and downstream dependencies are preserved.

This is a launch record, not a claim that figures are already complete.
Completion requires this directory's status.json to say complete and job
20862092 to be COMPLETED with exit 0:0. Expected deliverables are three PDF/PNG
pairs in figures/, 47,838 summary rows, 478,380 per-repeat rows, diagnostic
tables, interpretation and a checksum inventory.

One-off background collection can wait for successful completion:

```bash
python3 oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/collect_global.py --wait
```

Collection copies only manifest-listed artifacts and verifies their hashes.
Raw global model arrays remain on HPC. Collection requires the local computer
and its network connection; HPC computation continues independently.
