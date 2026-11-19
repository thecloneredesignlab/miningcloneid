#!/bin/bash -l
set -eu
study=$1
release=$2
bash "$study/releases/$release/job_python.sh" "$study" "$release" B2_prepare prepare_b2.py
bash "$study/releases/$release/job_sif.sh" "$study" "$release" B2_native native_kernel.R
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S "$study/releases/$release/audit_policy.py" "$study/verification/B2_native" R
bash "$study/releases/$release/job_python.sh" "$study" "$release" B2 compare_b2.py
