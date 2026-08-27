# Anchor and joint-function sensitivity

Seed 132 is a same-cluster near-optimal in-vitro anchor (objective delta 0.0007 from seed 10), whereas seed 157 is the best anchor in the other in-vitro cluster (objective delta 0.419).
Across the six frozen joint winners, 11/14 split parameters have a common in-vivo-versus-in-vitro direction in all six solutions, but a median 100.0% of active parameter penalties are in the saturating Welsch region.

The existing winner ensemble therefore supports a numerical direction-stability statement conditional on the seed-10 anchor and current Welsch penalty. It does not establish anchor or coupling-form invariance. The P0/P1 refit configurations are written to `hpc/joint_anchor_function_sensitivity_task_matrix.tsv` and remain unsubmitted.
