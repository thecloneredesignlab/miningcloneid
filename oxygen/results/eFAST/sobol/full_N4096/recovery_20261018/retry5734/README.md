# Restore Sobol shard 5734 and seed283

Controller commit: `5c598f7221f90085ae4d8e364f28ec659bd3dfc2`. Frozen statistical code: `f64ebfe689030b132de476593356d61c35b5b0fb`. Original model/sampling commit: `471f6a84afbf2363078a9cefce12a75cd0ee691f`.

The preflight certified 12,058 completed statistical shards, 499 completed neighborhoods, one completed global domain, and all 14 repaired model blocks. Only shard 5734 (seed283, repeat 2, oxygen indices 160–175) and case array index 173 required restoration. Existing source completion hashes are recorded in preflight.json.

All analysis stages use xxlarge, initial 12 hours, one CPU, no requested/excluded node and no array concurrency cap. Every downstream dependency requires success. Convergence is diagnostic only.

| Role | Array | Job | Initial monitor |
|---|---|---|---|
| shard | 5734 | 21139597 | 21139603 |
| case | 173 | 21139598 | 21139604 |
| audit | scalar | 21139599 | 21139605 |
| aggregate | scalar | 21139600 | 21139606 |
| response | 1-54 | 21139601 | 21139607 |
| finish | scalar | 21139602 | 21139608 |

Monitors use small/1 CPU/1G, check every 300 seconds and add six hours when a running task has less than six hours remaining. They verify internal task ownership and TimeLimit readback, and hand over before their own 30-day-12-hour limit. Their targets.json registry authorizes only this replacement chain.

The obsolete pending chain 21129866/21129867/21129868 and its monitors 21138286/21138287/21138288 were cancelled after replacement registration. The first new shard was released only after all six new monitors started.

HPC Python 3.6.8 passed all 20 operational tests, including the observed short and long container-identity failures and Slurm script spooling. Scientific calculations use the unchanged original SALib 1.5.2 SIF and immutable recovery plan.

Scientific results retain canonical paths under full_N4096. Final completion/accounting validation and collection receipts are recorded separately. To collect completed results, run oxygen/code/O2_supply_demand_MAP/analysis/eFAST/sobol/collect_full.py; it requires the final job and 501-domain audit to have passed and verifies all manifest checksums.


## Completion and collection

All 60 analysis/report/layout tasks finished COMPLETED with exit 0:0. Shard 5734
and seed283 are complete; all 355 originally unfinished neighborhoods have now
been restored. The audit passed for all 500 neighborhoods plus the global domain.
All 54 response surfaces and the final report were generated. The seven main
scope figures were redrawn only to make colorbar ticks legible; the saved scope
arrays were checksum-identical before and after rendering.

The local artifact_collection_receipt.json verifies 19,828 files (1,619,105,148
bytes) against the frozen HPC inventory. All 499 prior neighborhood completion
records and the global completion record remain byte-identical. Convergence
remains diagnostic only. final_validation.json contains Slurm accounting and
reuse checks; figure_rendering_validation.json records the representative PDF
rendering review.

The primary scope figures, summary source tables, model structure, inventory and
recovery receipts are versioned in Git. The complete detailed endpoint collection
is available in this local results tree and on HPC; collect_full.py retrieves and
verifies it from the inventory. Large model arrays and statistical checkpoints
remain on HPC.

All 54 conditional response figures, their source tables and the nine reference parameter sets are also versioned in Git.
