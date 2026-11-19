#!/bin/bash -l
# The allocation can execute the engineering gate or the first A benchmark only.
set -euo pipefail
study=$1
action=$2
case "$action" in
 engineering)
   exec "$study/release/runner.sh" "$study" E
   ;;
 benchmark)
   "$study/release/runner.sh" "$study" N_A_import
   exec "$study/release/runner.sh" "$study" N_A_bench
   ;;
 *) exit 78 ;;
esac
