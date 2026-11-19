# Independent tumor Step16 audit

This directory owns only the new tumor audit. It never imports a legacy
controller or model and contains no Slurm submission, cancellation, optimizer,
simulator, symbolic solver, or automatic gate-release path.

capture_authority.py creates a NEW invivo_structural run, freezes selected source
and metadata, and reads coherent old-owner state, scoped sacct/squeue and source/
SIF hashes over key-authenticated SSH. It refuses an existing destination. The
remote program writes only a stdout tar stream; no HPC file is changed.

export_contract.R reads the copied RDS and exports scalars, archived initial
vectors and observation metadata. It does not source or construct a model.
verify_audit.py verifies provenance, accounting, coordinate coverage and held
release conditions from those files. It is a metadata checker, not a structural
test or numerical replay.

Run from the repository root, with a fresh approved audit namespace:

    python3 -B .../invivo_step16_audit/capture_authority.py --run <new-run>
    Rscript --vanilla .../invivo_step16_audit/export_contract.R <new-run>
    python3 -B .../invivo_step16_audit/verify_audit.py --run <new-run>

The 20261113_invivo_structural_01 package is held at the upstream culture gate.
Future computation requires reviewed culture structural/bridge/workflow evidence
and an explicit tumor model/case freeze within the already conditional 32-hour
authorization. Full-model bridge and recovery have no released budget here.
