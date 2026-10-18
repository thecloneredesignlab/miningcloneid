# N4096 neighborhood recovery submission

Submission is recorded here; this is not a certificate that all analysis jobs finished.

- Frozen recovery code: `f64ebfe689030b132de476593356d61c35b5b0fb`.
- Original model code: `471f6a84afbf2363078a9cefce12a75cd0ee691f`; original SALib 1.5.2 SIF is unchanged.
- Restore 355 neighborhoods, retaining all 145 completed neighborhoods and the global results.
- Reuse 345,516 certified repeat/oxygen/output cells; calculate only 368,034 missing cells.
- Submit exactly 14 model blocks that previously failed before model execution.
- Statistical tasks: 11,903 ready shards and 156 shards waiting for model repair.
- Each shard covers one seed/repeat and at most 16 unfinished oxygen points.
- N4096 and the nested N2048 diagnostic, five repeats, 201 oxygen points, 14 factors/91 pairs,
  original bootstrap seeds and 200 resamples are retained. Convergence is diagnostic only.
- Every analyzer result is atomically saved; original checkpoints remain untouched.
- Raw outputs, source caches and statistical NPZ files stay on HPC.
- All stages request xxlarge and 12 hours. No node or array concurrency cap is specified.
- The actual original SIF passed 21 tests in preparation job 21126820 (COMPLETED, exit 0:0).
- Only the obsolete pending response/report jobs 20839426 and 20839427 are replaced.

| Stage | Job ID | Array | Dependency |
|---|---|---|---|
| prepare | 21126820 | none | none |
| model | 21126822 | 25318,27551,27552,33636,33637,41044,45450,48408,53810,63791,68684,69826,71559,74956 | none |
| shard | 21126823 | 1-11903 | none |
| shard | 21126824 | 11904-12059 | afterok:21126822 |
| case | 21129864 | 1-355 | afterany:21126823:21126824 |
| audit | 21129865 | none | afterany:21129864 |
| aggregate | 21129866 | none | afterok:21129865 |
| response | 21129867 | 1-54 | afterok:21129866 |
| finish | 21129868 | none | afterok:21129867 |

HPC recovery directory: `/share/lab_crd/taoli/Project/soft_couping_eFAST/oxygen/results/eFAST/sobol/full_N4096/recovery_20261018`.
Submission UTC: 2026-10-18T03:49:03.462293Z.
Metadata collection SHA256 checks passed for 362 source files; see
`launch_collection_manifest.tsv` and `execution_validation.json`.
The local plan and original-output inventory are losslessly archived as
`plan.json.gz` and `source_outputs.tsv.gz`; their decompressed SHA256 values
match the working HPC originals.
