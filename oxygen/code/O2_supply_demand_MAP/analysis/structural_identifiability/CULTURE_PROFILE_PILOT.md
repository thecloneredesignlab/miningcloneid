# Approved culture v2 pilot

The user approved this study on 2026-11-09: five unrestricted searches, about
100 restricted searches, and up to 20 reserve searches. All jobs use `xxlarge`.
There is no application concurrency limit or array `%N` throttle. Each fit is
submitted as an independent Slurm job as soon as its scientific inputs are
available. Short, single-writer reducers serialize metadata only, not fits.

The immutable study is `oxygen/results/structural_identifiability/20261109_culture_v2_pilot`.
Its `implementation` snapshot, code manifest, source/input provenance, complete
start bank, task manifests, actual native DE populations, scheduler accounting,
and independent score checks establish the execution contract. The baseline
must pass the independent-start and stronger-search gates before profiles.

The production simulator, passage-v2 rules, objective, bounds except the two
registered zero extensions, and modality normalization are unchanged. The old
WP0/WP1 outputs and their implementation snapshots remain historical evidence.

Each fit uses eight workers, 16 GB, 256 DE members, up to 500 generations and
50 local iterations. Reserve attempts use up to 1,000/100 and 12 hours. The
fit reservation ceiling is 8,640 core-hours; validation/controllers/reporting
have another 60. New attempts use new IDs and seeds, never overwritten output
directories. DE's internal state is not resumed.

The small additional R module checkpoints complete, improved feasible vectors
per PID, using a process-global cache that survives serialized DE closures and
an atomic on-disk fallback. All saved minima and the native endpoint are replayed
before selecting the lowest verified result. Native termination remains a
separate diagnostic. This is not a confidence-interval analysis.

Run the frozen `culture_profile_pilot.py submit PROJECT STUDY` after the external
validation job has completed successfully. Reducers then advance eligible work.
`advance` can reconcile an accounting-delay event, and `status` is read-only.
Submission intents without returned IDs require scheduler audit before reuse.

Positive profile points have 19 active coordinates; exact nulls have 18. The
grid contains zero, the unrestricted optimum and eight log-spaced positive
probes. Each point has five different projected primary starts. Missing or
duplicated neighbors are explicitly replaced by independent starts. Reverse
null checks, unresolved-point searches and at most one curvature refinement
per target consume the common 20-fit reserve. The curvature threshold is only
a numerical grid-allocation rule.

Every profile uses its own current unrestricted baseline; improved constrained
fits trigger baseline repair. Curves contain verified minima and search flags,
not averages, likelihood-ratio thresholds, p values or formal confidence limits.
