#!/bin/bash -l
set -eu
study=$1
module load sympy/1.12-gfbf-2023a IPython/8.14.0-GCCcore-12.3.0
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0
/app/eb/software/Python/3.11.3-GCCcore-12.3.0/bin/python3 -B -S "$study/audit_release/archive_audit.py" "$study"
