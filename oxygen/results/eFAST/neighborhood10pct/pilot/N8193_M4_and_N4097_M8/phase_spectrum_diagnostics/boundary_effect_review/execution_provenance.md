# Execution and provenance

- Initial evaluation-code commit: `944c4c20673f91c3a39abde510ecbb1044e1c6a8`.
- Smoke job: `21140535`, COMPLETED, exit 0:0, 24 domain tasks, 7m41s.
- Main array: `21140564`, 640/640 tasks COMPLETED, exit 0:0.
- All 664 task outputs were collected; all 166 baseline matrices matched exactly.
- Main submission HEAD: `0f7f3a488b99722a98afec4717b7f50ab95409be`.
  This descendant adds separate Sobol work. Every frozen evaluation/model file
  and input hash still matched the approved boundary experiment. Git SHA alone
  is not used as runtime-equivalence evidence.
- All jobs use xxlarge, 12h, 1 CPU/4GB. No node constraint or array cap.
- The SIF SHA256 is
  `0b60f6cdab8a91f6bbeea1ab6cf02dd79660a295f24bc6e8ad6c5fed04d982ad`.
- The frozen selection SHA256 is
  `82dcd39a4d4b1b0b8ff8e69840b860d93a693691168a541f7fb10ed232ac5b57`.

Before fast-forwarding the HPC checkout, 4,475 existing untracked files from
the previous precision review were compared byte-for-byte with the incoming
tracked versions. All matched; an independent backup was retained at
`/share/lab_crd/taoli/Project/soft_couping_eFAST_sync_backups/boundary_20261018_173629`.

The main-submission gate checked terminal Slurm success and the smoke receipt.
One early receipt lookup was temporarily unavailable after terminal accounting
became visible; the main array was submitted only after the receipt and logs
were confirmed. No smoke or scientific task was resubmitted.

The original high-precision and FAST outputs remain separate and unchanged.
`remote_collection_manifest.tsv` records the frozen HPC outputs before local
postprocessing. `artifact_manifest.tsv` adds the independent checks, figures,
interpretation and supplementary postprocessing script.
