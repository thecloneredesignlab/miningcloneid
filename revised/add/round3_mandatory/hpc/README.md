# HPC sensitivity preflight

No job in this directory has been submitted.

`joint_anchor_function_sensitivity_task_matrix.tsv` defines:

- P0: six in-vivo warm starts crossed with the same-cluster seed-132 and other-cluster seed-157 in-vitro anchors under the current coupling settings.
- P1: the six in-vivo warm starts paired with seed 10 under less-saturating Welsch settings (`c=1` and `c=10`).

`heldout_validation_design_matrix.tsv` records leave-one-tumor-out and leave-one-lineage-out folds. These folds require a code extension that removes the held-out unit from the objective while retaining it for prediction; they are not represented as runnable or completed.

Submission requires a separate user confirmation after the exact HPC checkout, command, resource request, output root, and job-array scope have been previewed.

